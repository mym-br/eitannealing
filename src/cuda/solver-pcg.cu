#include "solver-pcg.h"

#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

#include "settings.h"
#include "utils.h"

//#include "preconditioner.h"

using namespace cgl;

#define USE_CONSOLIDATED_KERNELS

// multiplicacao matriz vetor e subtracao (r = b - A * x)
// solver triangular inferior e superior (usando apenas o primeiro bloco)
__global__ void cpcg_mult_subtr_solver(int size, numType * aData, numType * precondData, int * aIndices, int * aRowLength,
	int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset,
	numType * b, numType * x, numType * r, numType * z, numType * partial, int blocks);

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_tot_esc_add_mmv_inner(int size, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * partial, int blocks);

// totalizacao (intra-bloco, gamma = pt.(A*p))
// escalar e vetor soma (x += (rmod/gamma) * p)
// escalar e vetor subtracao (r -= (rmod/gamma) * q)
// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
// solver triangular superior (z = inv(M) * r)
// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
__global__ void cpcg_tot_esc_add_sub_solver(int size, numType * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * x, numType * r, numType * z, numType * p, numType * q,
	numType * rmod, numType * gamma, numType * partial, int blocks);

// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_inner(int size, numType * r, numType * z, numType * partial, int blocks);

PCGSolverCPJDS::PCGSolverCPJDS(MatrixCPJDSManager * mgr, MatrixCPJDS *M, Vector * b) {
	this->mgr = mgr;
	this->A = M;
	this->b = b;
	this->x = new Vector(M->matrixData.n);

	r = new Vector(M->matrixData.n);
	z = new Vector(M->matrixData.n);
	p = new Vector(M->matrixData.n);
	q = new Vector(M->matrixData.n);
	u = new Vector(M->matrixData.n);
	partial = new Vector(M->matrixData.n);

	rmod = new Number(0);
	rmod_prev = new Number(1);
	rmod_aux = new Number(1);
	alpha = new Number(0);
	beta = new Number(0);
	gamma = new Number(0);

	size = A->matrixData.n;
	blocks = ceil((double)size / BLOCKSIZE);

	streamInit();
}

#ifdef USE_CONSOLIDATED_KERNELS

// Setup and calculate the 1st iteration	
void PCGSolverCPJDS::init() {
	//m_preconditioner(A, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	x->reset(stream);
	r->reset(stream);
	z->reset(stream);
	p->reset(stream);
	q->reset(stream);
	u->reset(stream);
	partial->reset(stream);

	numType * aData = (*A).matrixData.data;
	numType * precond = (*A).preconditionedData;
	int * aIndices = (*A).matrixData.indices;
	int * aRowLength = (*A).matrixData.rowLength;
	int * aRowSize = (*A).matrixData.rowSize;
	int * aColOffset = (*A).matrixData.colOffset;

	int colorCount = (*A).matrixColors.colorCount;
	int * colors = (*A).matrixColors.colors_d;
	int * colorsColOffset = (*A).matrixColors.colorsColOffsetSize_d;
	numType * zData = z->getData();
	numType * rData = r->getData();
	numType * xData = x->getData();
	numType * bData = b->getData();

	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	numType * partialData = partial->getData();

	it = 0;

	cpcg_mult_subtr_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, precond, aIndices, aRowLength,
		aRowSize, aColOffset, colorCount, colors, colorsColOffset,
		bData, xData, rData, zData, partialData, blocks);
	//LOGV(z->getData(),1, "z");
	//p = z;
	//z->copyTo(p);

	//cudaError_t cudaStatus = cudaGetLastError();
	//if (cudaStatus != cudaSuccess) {
	//	std::ostringstream msg;
	//	msg << "cpcg_mult_subtr_solver kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
	//	msg.flush(); LOG(&msg.str()[0]);
	//}

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);
	//cudaStatus = cudaGetLastError();
	//if (cudaStatus != cudaSuccess) {
	//	std::ostringstream msg;
	//	msg << "cpcg_inner kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
	//	msg.flush(); LOG(&msg.str()[0]);
	//}

	//LOGV(partialData, blocks, "partialData ==============================================");
}

void PCGSolverCPJDS::doIteration(int iteration) {
	it++;

	numType * aData = (*A).matrixData.data;
	numType * precond = (*A).preconditionedData;
	int * aIndices = (*A).matrixData.indices;
	int * aRowLength = (*A).matrixData.rowLength;
	int * aRowSize = (*A).matrixData.rowSize;
	int * aColOffset = (*A).matrixData.colOffset;

	int colorCount = (*A).matrixColors.colorCount;
	int * colors = (*A).matrixColors.colors_d;
	int * colorsColOffset = (*A).matrixColors.colorsColOffsetSize_d;
	numType * zData = z->getData();
	numType * rData = r->getData();
	numType * xData = x->getData();
	numType * pData = p->getData();
	numType * qData = q->getData();

	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	numType * partialData = partial->getData();

	//LOGV(partialData, blocks, "partialData ==============================================");

	//LOGV(pData, p->getSize(), "p before==============================================");
	//LOGV(zData, z->getSize(), "z ==============================================");
	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	cpcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
		zData, pData, qData, rmodData, rmod_prevData, partialData, blocks);

	//LOGV(rmodData, 1, "rmod");
	//LOGV(rmod_prevData, 1, "rmodprev");
	//LOGV(pData, p->getSize(), "p after ==============================================");
	//LOGV(qData, q->getSize(), "q ==============================================");


	//// Check for any errors launching the kernel
	//cudaError_t cudaStatus = cudaGetLastError();
	//if (cudaStatus != cudaSuccess) {
	//	std::ostringstream msg;
	//	msg << "cpcg_tot_esc_add_mmv_inner kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
	//	msg.flush(); LOG(&msg.str()[0]);
	//}

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
		xData, rData, zData, pData, qData, rmodData, gammaData, partialData, blocks);

	//LOGV(gammaData, 1, "gamma");
	//LOGV(xData, x->getSize(), "x ==============================================");
	//LOGV(rData, r->getSize(), "r ==============================================");


	//// Check for any errors launching the kernel
	//cudaStatus = cudaGetLastError();
	//if (cudaStatus != cudaSuccess) {
	//	std::ostringstream msg;
	//	msg << "cpcg_tot_esc_add_sub_solver kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
	//	msg.flush(); LOG(&msg.str()[0]);
	//}

	rmod_aux = rmod_prev; // variavel auxiliar para rotacao de ponteiros
	rmod_prev = rmod; // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo
	rmod = rmod_aux; // prepara novo rmod

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	//LOGV(partialData, partial->getSize(), "partialData ==============================================");

	//// Check for any errors launching the kernel
	//cudaStatus = cudaGetLastError();
	//if (cudaStatus != cudaSuccess) {
	//	std::ostringstream msg;
	//	msg << "cpcg_inner kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
	//	msg.flush(); LOG(&msg.str()[0]);
	//}
}

#else

// Setup and calculate the 1st iteration	
void PCGSolverCPJDS::init() {
	//m_preconditioner(A, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	x->reset();
	r->reset();
	z->reset();
	p->reset();
	q->reset();
	u->reset();
	partial->reset();

	numType * zData = z->getData();
	numType * rData = r->getData();
	numType * xData = x->getData();
	numType * bData = b->getData();

	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	numType * partialData = partial->getData();

	it = 0;
	//r = b - A * x;
	mgr->mult(*A, x, u); // u = A * x

	//numType *fileData = new numType[A.matrixData.n*A.matrixData.n];
	//for (int row = 0; row < A.matrixData.n; row++) {
	//	for (int col = 0; col < A.matrixData.n; col++) {
	//		int idx = mgr->coordinates2Index(row, col);
	//		numType data = 0;
	//		if (idx > -1) {
	//			int dcol = A.cpuData.indices[idx];
	//			data = A.cpuData.data[idx];
	//			//data = stiffness.cpuData.dataAccept[idx];
	//		}
	//		fileData[row * A.matrixData.n + col] = data;
	//	}
	//}
	//LOGM2(fileData, A.matrixData.n, A.matrixData.n, "Matrix A", LOGCPU);

	//numType *fileDataPrecon = new numType[A.matrixData.n*A.matrixData.n];
	//for (int row = 0; row < A.matrixData.n; row++) {
	//	for (int col = 0; col < A.matrixData.n; col++) {
	//		int idx = mgr->coordinates2Index(row, col);
	//		numType data = 0;
	//		if (idx > -1) {
	//			int dcol = A.cpuData.indices[idx];
	//			data = A.cpuData.precond[idx];
	//			//data = stiffness.cpuData.dataAccept[idx];
	//		}
	//		fileDataPrecon[row * A.matrixData.n + col] = data;
	//	}
	//}
	//LOGM2(fileDataPrecon, A.matrixData.n, A.matrixData.n, "Preconditioned Matrix A", LOGCPU);


	//LOGV(x->getData(), x->getSize(), "x");
	//LOGV(u->getData(), u->getSize(), "A*x");

	b->subtr(u, r); // r = b - u = b - A * x
	//LOGV(r->getData(), r->getSize(), "r = b - u = b - A * x");

	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	// solving L * u = r
	mgr->solve(*A, r, u);
	//LOGV(u->getData(), u->getSize(), "u");
	// now solve Lt * z = u
	mgr->solve_t(*A, u, z);
	//LOGV(z->getData(), z->getSize(), "z");

	//p = z;
	z->copyTo(p);
}

void PCGSolverCPJDS::doIteration(int iteration) {
	it++;

	// rmod = z.dot(r);
	r->inner(z, rmod);//isto nao deveria ser recalculado
	//LOGV(rmod->getData(), 1, "rmod");

	if (iteration != -1) {
		char filename[100];
		sprintf(filename, "LOG-MatrixA-%d.txt", iteration);
		mgr->saveToFile(filename, *A, A->matrixData.data, false);
	}

	// q = A * p;
	mgr->mult(*A, p, q);
	//LOGV(q->getData(), q->getSize(), "q");
	// gamma = q.dot(p);
	p->inner(q, gamma);
	//mgr->multInner(A, p, q, PCGSolverCPJDS::streams);
	//mgr->multInner2(A, p, q, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);

	//LOGV(gamma->getData(), 1, "gamma");
	//LOGV(p->getData(), p->getSize(), "p");
	// x += alpha * p; alpha = rmod / gamma
	x->scalarAdd(rmod, gamma, p, NULL);
	//LOGV(x->getData(), x->getSize(), "x");

	//r -= alpha * q; alpha = rmod / gamma
	r->scalarSubtr(rmod, gamma, q, NULL);
	//LOGV(q->getData(), q->getSize(), "q ================================================================");
	//LOGV(r->getData(), r->getSize(), "r ================================================================");

	// substitutes both above
	// totalizes partial to calculate gamma
	// x += alpha * p; alpha = rmod / gamma
	// r -= alpha * q; alpha = rmod / gamma
	//p->scalarAddSubtrTotalization(rmod, gamma, x, r, q, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);

	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	// solving L * u = r
	mgr->solve(*A, r, u);
	// now solve Lt * z = u
	mgr->solve_t(*A, u, z);
	//mgr->solve_complete(A, r, z, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	//mgr->solve_and_inner(A, r, z, u, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	//LOGV(z->getData(), z->getSize(), "z ================================================================");

	//rmod_1 = rmod;
	rmod_aux = rmod_prev;
	rmod_prev = rmod;
	rmod = rmod_aux;

	// rmod = z.dot(r);
	z->inner(r, rmod);
	//LOGV(rmod->getData(), 1, "rmod ================================================================");

	// p = z + beta * p; beta = rmod / rmod_prev
	//p->scalarAddTotalization(rmod, rmod_prev, z, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	p->scalar(rmod, rmod_prev, u);
	u->sum(z, p);

	//LOGV(rmod->getData(), 1, "rmod");
	//LOGV(rmod_prev->getData(), 1, "rmodprev");
	//LOGV(p->getData(), p->getSize(), "pfinal ================================================================");
}

#endif

void PCGSolverCPJDS::streamInit() {
	cudaStreamCreate(&this->stream);
}

void PCGSolverCPJDS::streamDestroy() {
	cudaStreamDestroy(this->stream);
}

// multiplicacao matriz vetor e subtracao (r = b - A * x)
// solver triangular inferior e superior - usando apenas o primeiro bloco (cor precisa caber nessas threads)
// z = inv(M) * r
__global__ void cpcg_mult_subtr_solver(int size, numType * aData, numType * precondData, int * aIndices, int * aRowLength,
	int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, 
	numType * b, numType * x, numType * r, numType * z, numType * partial, int blocks) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// matrix-vector multiplication
	if (row < size) {
		int colorIdx = 0;

		for (; colorIdx < colorCount; colorIdx++) {
			if (row < colors[colorIdx]) {
				break;
			}
		}
		colorIdx--; // must decrease 1 due to for adding a ++ even at the break

		int colorStart = colors[colorIdx];
		int colorColOffset = colorsColOffset[colorIdx];

		// row size (length + padding zeros)
		int rowSize = aRowSize[row];

		numType sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * x[idx]; // NOT coalesced!

			__syncthreads(); // synchronization so all threads load from memory
		}
		r[row] = b[row] - sum; // coalesced, resultado da multiplicacao matriz-vetor, ja subtraindo b
	}

	//if (row == 0) {
	//	printf("b[%d] = %g\n", b[row]);
	//	printf("r[%d] = %g\n", r[row]);
	//}

	__syncthreads();

	// lower triangular solver
	partial[row] = 0.0;
	for (int k = 0; k < colorCount; k++) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		for (row = tidx + colorStart; row < colorEnd; row += BLOCKSIZE) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				numType rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * partial[idx];
				}
				//__syncthreads();
			}
			partial[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
		}
		__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
	}

	__syncthreads();
	//obs: most threads will just pass by, but that's ok, as we need only the number of threads equal to the
	//largest color group size - as long as it fits in a single block (we need it to use __syncthreads)

	// upper triangular kernel
	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		for (row = tidx + colorStart; row < colorEnd; row += BLOCKSIZE) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				numType rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx > row && idx > -1) { // main diagonal can be skiped
					sum += rowData * z[idx];
				}
				//__syncthreads();
			}
			// using partial result from previous (lower triangular) solver
			z[row] = (partial[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
		}
		__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
	}

	__syncthreads();

	//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
}

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_tot_esc_add_mmv_inner(int size, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * z, numType * p, numType * q, 
	numType * rmod, numType * rmod_prev, numType * partial, int blocks) {
	//printf("blocks = %d\n", blocks);

	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// partial totalization
	__shared__ numType total;
	if (tidx == 0) {
		total = 0;
		for (int i = 0; i < blocks; i++) {
			total += partial[i];//partial must be initialized with all zeroes
		}
		if (row == 0) {
			rmod[0] = total;
		}
	}

	__syncthreads();

	//total = rmod[0];
	// scalar add
	const numType beta = total / rmod_prev[0]; //rmod_prev must be initialized with non-zero
	if (row < size) {
		// p = z + beta * p
		p[row] = z[row] + beta * p[row];
	}

	__syncthreads();

	// matrix-vector multiplication
	cache[tidx] = 0.0; //cache inicializado para realizar produto interno ao mesmo tempo da MMV
	if (row < size) {
		int colorIdx = 0;
		for (; colorIdx < colorCount; colorIdx++) {
			if (row < colors[colorIdx]) {
				break;
			}
		}
		colorIdx--; // must decrease 1 due to for adding a ++ even at the break

		int colorStart = colors[colorIdx];
		int colorColOffset = colorsColOffset[colorIdx];

		// row size (length + padding zeros)
		int rowSize = aRowSize[row];

		numType sum = 0;
		//__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * p[idx]; // NOT coalesced!

			//__syncthreads(); // synchronization so all threads load from memory
		}
		q[row] = sum; // coalesced, resultado da multiplicacao matriz-vetor
		cache[tidx] = sum * p[row]; //produto interno: multiplicacao ponto a ponto
	}

	__syncthreads();

	//produto interno: totalizacao dentro do bloco
	partial[row] = 0;
	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (tidx < half) {
			cache[tidx] += cache[tidx + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (tidx == 0) {
		partial[blockIdx.x] = cache[0];
	}

	__syncthreads();
}

// totalizacao (intra-bloco, gamma = pt.(A*p))
// escalar e vetor soma (x += (rmod/gamma) * p)
// escalar e vetor subtracao (r -= (rmod/gamma) * q)
// solver triangular inferior (z = inv(M) * r) - utiliza somente um bloco (cor precisa caber nessas threads)
// solver triangular superior (z = inv(M) * r) - utiliza somente um bloco (cor precisa caber nessas threads)
// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
__global__ void cpcg_tot_esc_add_sub_solver(int size, numType * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * x, numType * r, numType * z, numType * p, numType * q,
	numType * rmod, numType * gamma, numType * partial, int blocks) {

	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// partial totalization
	__shared__ numType total;
	if (tidx == 0) {
		total = 0;
		for (int i = 0; i < blocks; i++) {
			total += partial[i];//partial must be initialized with all zeroes
		}
		if (row == 0) {
			gamma[0] = total;
		}
	}

	__syncthreads();

	// scalar add and subtraction
	const numType alpha = rmod[0] / total; //rmod_prev must be initialized with non-zero
	if (row < size) {
		// x += alpha * p
		x[row] += alpha * p[row];
		// r -= alpha * q
		r[row] -= alpha * q[row];
	}

	__syncthreads();

	// lower triangular solver
	partial[row] = 0.0;
	if (blockIdx.x == 0) {
		for (int k = 0; k < colorCount; k++) {

			int colorStart = colors[k];
			int colorEnd = colors[k + 1];
			int colorColOffset = colorsColOffset[k];

			for (row = tidx + colorStart; row < colorEnd; row += BLOCKSIZE) {
				int rowStep = (row - colorStart);
				int rowSize = aRowSize[row];

				numType sum = 0;
				//__syncthreads();

				for (int j = 1; j < rowSize; j++) { // first element is main diagonal
					// colorColOffset already includes colorOffset (thus, color's first row)
					int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

					numType rowData = precondData[offset]; // coalesced
					int idx = aIndices[offset]; // coalesced
					if (idx < row) { // main diagonal can be skiped
						sum += rowData * partial[idx];
					}
					//__syncthreads();
				}
				partial[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
			}
			__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
		}

		//__syncthreads();
		//obs: most threads will just pass by, but that's ok, as we need only the number of threads equal to the
		//largest color group size - as long as it fits in a single block (we need it to use __syncthreads)

		// upper triangular kernel
		for (int k = colorCount - 1; k >= 0; k--) {

			int colorStart = colors[k];
			int colorEnd = colors[k + 1];
			int colorColOffset = colorsColOffset[k];

			for (row = tidx + colorStart; row < colorEnd; row += BLOCKSIZE) {
				int rowStep = (row - colorStart);
				int rowSize = aRowSize[row];

				numType sum = 0;
				//__syncthreads();

				for (int j = 1; j < rowSize; j++) { // first element is main diagonal
					// colorColOffset already includes colorOffset (thus, color's first row)
					int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

					numType rowData = precondData[offset]; // coalesced
					int idx = aIndices[offset]; // coalesced
					if (idx > row && idx > -1) { // main diagonal can be skiped
						sum += rowData * z[idx];
					}
					//__syncthreads();
				}
				// using partial result from previous (lower triangular) solver
				z[row] = (partial[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
			}
			__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
		}

		//__syncthreads();
		//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here

	}
}

// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_inner(int size, numType * r, numType * z, numType * partial, int blocks) {

	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// produto interno
	cache[tidx] = 0.0; //cache inicializado para realizar produto interno ao mesmo tempo do solver
	if (row < size) {
		cache[tidx] = z[row] * r[row];//produto interno: multiplicacao ponto a ponto
	}

	__syncthreads();

	//produto interno: totalizacao dentro do bloco
	partial[row] = 0;
	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (tidx < half) {
			cache[tidx] += cache[tidx + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (tidx == 0) {
		partial[blockIdx.x] = cache[0];
	}
	// totalizacao intra blocos sera realizada no proximo kernel

	__syncthreads();
}