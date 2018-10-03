#include "solver-pcg.h"

#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <iomanip>

#include "settings.h"
#include "utils.h"
//#include "preconditioner.h"

void checkedCudaEventRecord(cudaEvent_t &event) {
#ifdef CGTIMING 
	cudaEventRecord(event);
#endif
}

using namespace cgl;

// multiplicacao matriz vetor e subtracao (r = b - A * x)
// solver triangular inferior e superior (usando apenas o primeiro bloco)
__global__ void cpcg_mult_subtr_solver(int size, numType * aData, numType * precondData, int * aIndices, int * aRowLength,
	int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset,
	numType * b, numType * x, numType * r, numType * z, numType * partial, int blocks);

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_tot_esc_add_mmv_inner0(int size, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * partialrmod, numType * partialgamma, int blocks);

// totalizacao (intra-bloco, gamma = pt.(A*p))
// escalar e vetor soma (x += (rmod/gamma) * p)
// escalar e vetor subtracao (r -= (rmod/gamma) * q)
// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
// solver triangular superior (z = inv(M) * r)
// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
__global__ void cpcg_tot_esc_add_sub_solver0(int size, numType * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * x, numType * r, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * gamma, numType * partial, int blocks);

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_tot_esc_add_mmv_inner(int size, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * partialrmod, numType * partialgamma, int blocks);

// totalizacao (intra-bloco, gamma = pt.(A*p))
// escalar e vetor soma (x += (rmod/gamma) * p)
// escalar e vetor subtracao (r -= (rmod/gamma) * q)
// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
// solver triangular superior (z = inv(M) * r)
// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
__global__ void cpcg_tot_esc_add_sub_solver(int size, numType * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * x, numType * r, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * gamma, numType * partial, int blocks);

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
	gamma = new Number(0);

	size = A->matrixData.n;
	blocks = ceil((double)size / BLOCKSIZE);
	partial2 = new Vector(this->blocks);

	totalItTime = totalTriangularTime = totalSpmvTime = 0;
	streamInit();
}

void PCGSolverCPJDS::init(double res) {
	x->reset(stream);
	this->init(x, res);
}

// Setup and calculate the 1st iteration	
void PCGSolverCPJDS::init(Vector *x0, double res) {
	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
#endif // CGTIMING

	numType *data_h = new numType[1];
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
	numType * partialData = partial->getData();

	it = 0;

	// Initial x value
	x0->copyTo(x);

	//r = b - A * x;
	mgr->mult(*A, x, u); // u = A * x
	b->subtr(u, r); // r = b - u = b - A * x
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	z->copyTo(p); //p = z;
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	mgr->mult(*A, p, q); // q = A * p;
	p->inner(q, gamma); // gamma = q.dot(p);
	#ifdef CALCULATE_ERRORS
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2_1 = rmod2 = *data_h;
	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;
	if (rmod2 < res) { it = 0; return; }

	// Error calculations
	r0norm2 = rmod2;
	r0norm = sqrt(r0norm2);
	beta = 0;
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2;
	#endif

	// Now do the first 3 iterations
	checkedCudaEventRecord(startTotal);
	x->scalarAdd(rmod, gamma, p, NULL); // 1...
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = alpha;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	c = rt1 / r1;
	s = eta_p1 / r1;
	wt[0] = 1 / rt1;
	w[0] = 1 / r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;
	if (rmod2 < res) { it = 1; return; }

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING

	checkedCudaEventRecord(startTotal);
	x->scalarAdd(rmod, gamma, p, NULL); // 2...
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = c * alpha - s * eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c * eta + s * alpha;	// r_2,2 = c_1*eta_2
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[1] = wt[1] = -r2 * w[0];
	wt[1] /= rt1;
	w[1] /= r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;
	if (rmod2 < res) { it = 2; return; }

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING

	checkedCudaEventRecord(startTotal);
	x->scalarAdd(rmod, gamma, p, NULL); // 3...
	if (rmod2 < res) { rmod2_1 = rmod2; it = 2; return; }
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = c * alpha - s * c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c_1 * c*eta + s * alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1 * eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[2] = wt[2] = -(r3*w[0] + r2 * w[1]);
	wt[2] /= rt1;
	w[2] /= r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif //CGTIMING

#ifdef CALCULATE_ERRORS
	err[0] = wt[0] * wt[0];
	err[1] = w[0] * w[0] + wt[1] * wt[1];
	err[2] = w[0] * w[0] + w[1] * w[1] + wt[2] * wt[2];
#endif
	it = 3;

	delete data_h;
}

void PCGSolverCPJDS::doIteration(int iteration) {
	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
#endif // CGTIMING

	checkedCudaEventRecord(startTotal);
	numType *data_h = new numType[1];
	it++;
	x->scalarAdd(rmod, gamma, p, NULL);
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c * alpha - s * c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c_1 * c*eta + s * alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1 * eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// FIXME: GET RID OF THOSE STUPID BUFFERS!!!!!!!!!!!!!!!!!!!!!
	if (it<360) {
		w[it - 1] = wt[it - 1] = -(r3*w[it - 3] + r2 * w[it - 2]);
		wt[it - 1] /= rt1;
		w[it - 1] /= r1;
	}
	else w[0] = it;
	err[it - 1] = w[it - 2] * w[it - 2] + wt[it - 1] * wt[it - 1] - wt[it - 2] * wt[it - 2];

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	delete data_h;
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING
}

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

			//__syncthreads(); // synchronization so all threads load from memory
		}
		r[row] = b[row] - sum; // coalesced, resultado da multiplicacao matriz-vetor, ja subtraindo b
	}

	// lower triangular solver
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
						sum += rowData * z[idx];
					}
					//__syncthreads();
				}
				z[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
			}
			__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
		}
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
				z[row] = (z[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
			}
			__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
		}
		//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
	}
}

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_tot_esc_add_mmv_inner0(int size, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * partialrmod, numType * partialgamma, int blocks) {
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
			total += partialrmod[i];//partial must be initialized with all zeroes
		}
		if (row == 0) {
			rmod[0] = total;
		}
	}

	__syncthreads();

	//total = rmod[0];
	// p = z for the first iteration
	if (row < size) {
		p[row] = z[row];
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
			sum += rowData * z[idx]; // NOT coalesced!

									 //__syncthreads(); // synchronization so all threads load from memory
		}
		q[row] = sum; // coalesced, resultado da multiplicacao matriz-vetor
		cache[tidx] = sum * p[row]; //produto interno: multiplicacao ponto a ponto
	}

	__syncthreads();

	//produto interno: totalizacao dentro do bloco
	int half = BLOCKSIZE >> 1;
	while (half != 0) {
		if (tidx < half) {
			cache[tidx] += cache[tidx + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (tidx == 0) {
		partialgamma[blockIdx.x] = cache[0];
	}

	__syncthreads();
}

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcg_tot_esc_add_mmv_inner(int size, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * z, numType * p, numType * q, 
	numType * rmod, numType * rmod_prev, numType * partialrmod, numType * partialgamma, int blocks) {
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
			total += partialrmod[i];//partial must be initialized with all zeroes
		}
		if (row == 0) {
			rmod[0] = total;
		}
	}
	__syncthreads();
	
	//total = rmod[0];
	// scalar add
	const numType beta = total / rmod_prev[0]; //rmod_prev must be initialized with non-zero
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
		int rowSize = aRowLength[row];

		numType sum = 0;
		//__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * (z[idx] + beta * p[idx]); // NOT coalesced!
			//__syncthreads(); // synchronization so all threads load from memory
		}
		q[row] = sum; // coalesced, resultado da multiplicacao matriz-vetor
		cache[tidx] = sum * (z[row] + beta * p[row]); //produto interno: multiplicacao ponto a ponto
	}

	__syncthreads();

	//produto interno: totalizacao dentro do bloco
	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (tidx < half) {
			cache[tidx] += cache[tidx + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (tidx == 0) {
		partialgamma[blockIdx.x] = cache[0];
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
	numType * rmod, numType * rmod_prev, numType * gamma, numType * partial, int blocks) {
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
	const numType beta = rmod[0] / rmod_prev[0]; //rmod_prev must be initialized with non-zero
	if (row < size) {
		// p = z + beta * p
		p[row] = z[row] + beta * p[row];
		// x += alpha * p
		x[row] += alpha * p[row];
		// r -= alpha * q
		r[row] -= alpha * q[row];

		// lower triangular solver
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
							sum += rowData * z[idx];
						}
						//__syncthreads();
					}
					z[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
				}
				__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
			}

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
					z[row] = (z[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
				}
				__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
			}
			//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
		}
	}
}

// totalizacao (intra-bloco, gamma = pt.(A*p))
// escalar e vetor soma (x += (rmod/gamma) * p)
// escalar e vetor subtracao (r -= (rmod/gamma) * q)
// solver triangular inferior (z = inv(M) * r) - utiliza somente um bloco (cor precisa caber nessas threads)
// solver triangular superior (z = inv(M) * r) - utiliza somente um bloco (cor precisa caber nessas threads)
// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
__global__ void cpcg_tot_esc_add_sub_solver0(int size, numType * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * x, numType * r, numType * z, numType * p, numType * q,
	numType * rmod, numType * rmod_prev, numType * gamma, numType * partial, int blocks) {
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
	const numType beta = total / rmod_prev[0]; //rmod_prev must be initialized with non-zero
	if (row < size) {
		// p = z
		p[row] = z[row];
		// x += alpha * p
		x[row] += alpha * p[row];
		// r -= alpha * q
		r[row] -= alpha * q[row];

		// lower triangular solver
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
							sum += rowData * z[idx];
						}
						//__syncthreads();
					}
					z[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
				}
				__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
			}

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
					z[row] = (z[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
				}
				__syncthreads();// make all threads in the first block stop and wait - the others can be discarded
			}
			//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
		}
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

void PCGSolverCPJDS2::init(Vector *x0, double res) {
	//m_preconditioner(A, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	x->reset(stream); x_1->reset(stream);
	r->reset(stream);
	z->reset(stream);
	p->reset(stream);
	q->reset(stream);
	u->reset(stream);
	partial->reset(stream);
	partial2->reset(stream);

	//numType * aData = (*A).matrixData.data;
	//numType * precond = (*A).preconditionedData;
	//int * aIndices = (*A).matrixData.indices;
	//int * aRowLength = (*A).matrixData.rowLength;
	//int * aRowSize = (*A).matrixData.rowSize;
	//int * aColOffset = (*A).matrixData.colOffset;

	//int colorCount = (*A).matrixColors.colorCount;
	//int * colors = (*A).matrixColors.colors_d;
	//int * colorsColOffset = (*A).matrixColors.colorsColOffsetSize_d;
	numType * zData = z->getData();
	numType * rData = r->getData();
	numType * xData = x->getData();
	numType * bData = b->getData();
	numType * pData = p->getData();
	numType * qData = q->getData();
	numType * partialData = partial->getData();
	numType * partialData2 = partial2->getData();

	it = 0;

	// Initial x value
	x0->copyTo(x);

	// multiplicacao matriz vetor e subtracao (r = b - A * x)
	// solver triangular inferior e superior - usando apenas o primeiro bloco (cor precisa caber nessas threads)
	// z = inv(M) * r
	cpcg_mult_subtr_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, (*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(),
		(*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(),
			bData, xData, rData, zData, partialData, blocks);

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	doIteration0((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);
	doIteration1((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);
	if (rmod2 < res) { it = 1; return; }
	doIteration2((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);
	if (rmod2 < res) { it = 2; return; }
	doIteration3((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);

#ifdef CALCULATE_ERRORS
	err[0] = wt[0] * wt[0];
	err[1] = w[0] * w[0] + wt[1] * wt[1];
	err[2] = w[0] * w[0] + w[1] * w[1] + wt[2] * wt[2];
#endif
}

void PCGSolverCPJDS2::doIteration0(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData, numType * partialData2) {
	numType *data_h = new numType[1];
	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	cpcg_tot_esc_add_mmv_inner0 << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_tot_esc_add_sub_solver0 << <blocks, BLOCKSIZE, 0, stream >> >
		(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;
	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Error calculations
	r0norm2 = rmod2;
	r0norm = sqrt(r0norm2);
	beta = 0;
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	delete data_h;
}

void PCGSolverCPJDS2::doIteration1(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData, numType * partialData2) {
	numType *data_h = new numType[1];
	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING
	it++;
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	checkedCudaEventRecord(startTri);
	cpcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = alpha;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	c = rt1 / r1;
	s = eta_p1 / r1;
	wt[0] = 1 / rt1;
	w[0] = 1 / r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING
	delete data_h;
}

void PCGSolverCPJDS2::doIteration2(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData, numType * partialData2) {
	numType *data_h = new numType[1];
	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING
	it++;
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	checkedCudaEventRecord(startTri);
	cpcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = c * alpha - s * eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c * eta + s * alpha;	// r_2,2 = c_1*eta_2
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[1] = wt[1] = -r2 * w[0];
	wt[1] /= rt1;
	w[1] /= r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING
	delete data_h;
}

void PCGSolverCPJDS2::doIteration3(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData, numType * partialData2) {
	numType *data_h = new numType[1];
	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING
	it++;
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	checkedCudaEventRecord(startTri);
	cpcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = c * alpha - s * c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c_1 * c*eta + s * alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1 * eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[2] = wt[2] = -(r3*w[0] + r2 * w[1]);
	wt[2] /= rt1;
	w[2] /= r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING
	delete data_h;
}

void PCGSolverCPJDS2::doIteration(int iteration) {
	numType *data_h = new numType[1];

	//numType * aData = (*A).matrixData.data;
	//numType * precond = (*A).preconditionedData;
	//int * aIndices = (*A).matrixData.indices;
	//int * aRowLength = (*A).matrixData.rowLength;
	//int * aRowSize = (*A).matrixData.rowSize;
	//int * aColOffset = (*A).matrixData.colOffset;

	//int colorCount = (*A).matrixColors.colorCount;
	//int * colors = (*A).matrixColors.colors_d;
	//int * colorsColOffset = (*A).matrixColors.colorsColOffsetSize_d;
	numType * zData = z->getData();
	numType * rData = r->getData();
	numType * xData = x->getData();
	numType * x2Data = x_1->getData();
	numType * pData = p->getData();
	numType * qData = q->getData();

	numType * rmodData = rmod->getData();
	numType * rmod_prevData = rmod_prev->getData();
	numType * gammaData = gamma->getData();

	numType * partialData = partial->getData();
	numType * partialData2 = partial2->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING

	it++;
	cudaMemcpy(x2Data, xData, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, (*A).matrixData.data.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(),
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	checkedCudaEventRecord(startTri);
	cpcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
		(size, (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(),
			xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c * alpha - s * c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c_1 * c*eta + s * alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1 * eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(numType), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// FIXME: GET RID OF THOSE STUPID BUFFERS!!!!!!!!!!!!!!!!!!!!!
	if (it<360) {
		w[it - 1] = wt[it - 1] = -(r3*w[it - 3] + r2 * w[it - 2]);
		wt[it - 1] /= rt1;
		w[it - 1] /= r1;
	}
	else {
		w[0] = it;
	}

	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[it] = lalpha;
	eta_[it] = leta;*/

	err[it - 1] = w[it - 2] * w[it - 2] + wt[it - 1] * wt[it - 1] - wt[it - 2] * wt[it - 2];

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);
	checkedCudaEventRecord(stopTotal);

#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING

	delete data_h;
}