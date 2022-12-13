#ifdef CGROUPS
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
#include <cooperative_groups.h>

#include "settings.h"
#include "utils.h"
//#include "preconditioner.h"

using namespace cgl;
namespace cg = cooperative_groups;

// multiplicacao matriz vetor e subtracao (r = b - A * x)
// solver triangular inferior e superior - usando apenas o primeiro bloco (cor precisa caber nessas threads)
// z = inv(M) * r
__global__ void cpcgcg_mult_subtr_solver(int size, double * aData, double * precondData, int * aIndices, int * aRowLength,
	int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset,
	double * b, double * x, double * r, double * z, double * partial, int blocks) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;
	cg::grid_group grid = cg::this_grid();

	// matrix-vector multiplication
	int colorIdx = 0;

	for (; colorIdx < colorCount; colorIdx++) {
		if (row < colors[colorIdx]) {
			break;
		}
	}
	colorIdx--; // must decrease 1 due to for adding a ++ even at the break

	int colorStart = colors[colorIdx];
	int colorColOffset = colorsColOffset[colorIdx];

	double sum = 0;
	__syncthreads();
	if (row < size) {
		for (int j = 0; j < aRowSize[row]; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			double rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * x[idx]; // NOT coalesced!
		}
		r[row] = b[row] - sum; // coalesced, resultado da multiplicacao matriz-vetor, ja subtraindo b
	}

	cg::sync(grid);
	tidx = blockDim.x * blockIdx.x + threadIdx.x;
	// lower triangular solver
	for (int k = 0; k < colorCount; k++) {
		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		row = tidx + colorStart;
		if (row < colorEnd) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];
			double sum = 0;
			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				double rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) sum += rowData * z[idx]; // main diagonal can be skiped
			}
			z[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
		}
		cg::sync(grid);// make all threads stop and wait
	}

	// upper triangular kernel
	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		row = tidx + colorStart;
		if (row < colorEnd) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];
			double sum = 0;
			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				double rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx > row && idx > -1) sum += rowData * z[idx]; // main diagonal can be skiped
			}
			// using partial result from previous (lower triangular) solver
			z[row] = (z[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
		}
		cg::sync(grid);// make all threads stop and wait
	}
	//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
}

// totalizacao (intra-blocos, rmod = zt.r)
// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
// matriz-vetor (q = A * p)
// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcgcg_tot_esc_add_mmv_inner0(int size, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * z, double * p, double * q,
	double * rmod, double * rmod_prev, double * partialrmod, double * partialgamma, int blocks) {
	//printf("blocks = %d\n", blocks);

	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// partial totalization
	__shared__ double total;
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

		double sum = 0;
		//__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			double rowData = aData[offset]; // coalesced
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
__global__ void cpcgcg_tot_esc_add_mmv_inner(int size, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * z, double * p, double * q, 
	double * rmod, double * rmod_prev, double * partialrmod, double * partialgamma, int blocks) {
	//printf("blocks = %d\n", blocks);

	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// partial totalization
	__shared__ double total;
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
	const double beta = total / rmod_prev[0]; //rmod_prev must be initialized with non-zero
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

		double sum = 0;
		//__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			double rowData = aData[offset]; // coalesced
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
__global__ void cpcgcg_tot_esc_add_sub_solver(int size, double * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * x, double * r, double * z, double * p, double * q,
	double * rmod, double * rmod_prev, double * gamma, double * partial, int blocks) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;
	cg::grid_group grid = cg::this_grid();

	// partial totalization
	__shared__ double total;
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
	const double alpha = rmod[0] / total; //rmod_prev must be initialized with non-zero
	const double beta = rmod[0] / rmod_prev[0]; //rmod_prev must be initialized with non-zero
	if (row < size) {
		// p = z + beta * p
		p[row] = z[row] + beta * p[row];
		// x += alpha * p
		x[row] += alpha * p[row];
		// r -= alpha * q
		r[row] -= alpha * q[row];
	}

	cg::sync(grid);
	tidx = blockDim.x * blockIdx.x + threadIdx.x;
	// lower triangular solver
	for (int k = 0; k < colorCount; k++) {
		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		row = tidx + colorStart;
		if (row < colorEnd) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];
			double sum = 0;
			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				double rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) sum += rowData * z[idx]; // main diagonal can be skiped
			}
			z[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
		}
		cg::sync(grid);// make all threads stop and wait
	}

	// upper triangular kernel
	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		row = tidx + colorStart;
		if (row < colorEnd) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];
			double sum = 0;
			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				double rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx > row && idx > -1) sum += rowData * z[idx]; // main diagonal can be skiped
			}
			// using partial result from previous (lower triangular) solver
			z[row] = (z[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
		}
		cg::sync(grid);// make all threads stop and wait
	}
	//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
}

// totalizacao (intra-bloco, gamma = pt.(A*p))
// escalar e vetor soma (x += (rmod/gamma) * p)
// escalar e vetor subtracao (r -= (rmod/gamma) * q)
// solver triangular inferior (z = inv(M) * r) - utiliza somente um bloco (cor precisa caber nessas threads)
// solver triangular superior (z = inv(M) * r) - utiliza somente um bloco (cor precisa caber nessas threads)
// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
__global__ void cpcgcg_tot_esc_add_sub_solver0(int size, double * precondData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * x, double * r, double * z, double * p, double * q,
	double * rmod, double * rmod_prev, double * gamma, double * partial, int blocks) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;
	cg::grid_group grid = cg::this_grid();

	// partial totalization
	__shared__ double total;
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
	const double alpha = rmod[0] / total; //rmod_prev must be initialized with non-zero
	const double beta = total / rmod_prev[0]; //rmod_prev must be initialized with non-zero
	if (row < size) {
		// p = z
		p[row] = z[row];
		// x += alpha * p
		x[row] += alpha * p[row];
		// r -= alpha * q
		r[row] -= alpha * q[row];
	}

	cg::sync(grid);
	tidx = blockDim.x * blockIdx.x + threadIdx.x;
	// lower triangular solver
	for (int k = 0; k < colorCount; k++) {
		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		row = tidx + colorStart;
		if (row < colorEnd) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];
			double sum = 0;
			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				double rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) sum += rowData * z[idx]; // main diagonal can be skiped
			}
			z[row] = (r[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];
		}
		cg::sync(grid);// make all threads stop and wait
	}

	// upper triangular kernel
	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		row = tidx + colorStart;
		if (row < colorEnd) {
			int rowStep = (row - colorStart);
			int rowSize = aRowSize[row];
			double sum = 0;
			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + rowStep; // coalesced?

				double rowData = precondData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx > row && idx > -1) sum += rowData * z[idx]; // main diagonal can be skiped
			}
			// using partial result from previous (lower triangular) solver
			z[row] = (z[row] - sum) / precondData[aColOffset[colorColOffset] + rowStep];// resultado do solver linear
		}
		cg::sync(grid);// make all threads stop and wait
	}
	//obs: it is a shame we lost all threads not in the first blocks, as because of that we cannot copy z to p here
}

// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
__global__ void cpcgcg_inner(int size, double * r, double * z, double * partial, int blocks) {

	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
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

void PCGSolverConsolidatedCPJDSCG::init(Vector *x0, double res) {
	//m_preconditioner(A, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	x->reset(stream); x_1->reset(stream);
	r->reset(stream);
	z->reset(stream);
	p->reset(stream);
	q->reset(stream);
	u->reset(stream);
	partial->reset(stream);
	partial2->reset(stream);

	//double * aData = (*A).matrixData.data;
	//double * precond = (*A).preconditionedData;
	//int * aIndices = (*A).matrixData.indices;
	//int * aRowLength = (*A).matrixData.rowLength;
	//int * aRowSize = (*A).matrixData.rowSize;
	//int * aColOffset = (*A).matrixData.colOffset;

	//int colorCount = (*A).matrixColors.colorCount;
	//int * colors = (*A).matrixColors.colors_d;
	//int * colorsColOffset = (*A).matrixColors.colorsColOffsetSize_d;
	double * zData = z->getData();
	double * rData = r->getData();
	double * xData = x->getData();
	double * bData = b->getData();
	double * pData = p->getData();
	double * qData = q->getData();
	double * partialData = partial->getData();
	double * partialData2 = partial2->getData();

	double *matrixData_d = (*A).matrixData.data.get();
	double *preconditionedData_d = (*A).preconditionedData.get();
	int *indices_d = (*A).matrixData.indices.get();
	int *rowLength_d = (*A).matrixData.rowLength.get();
	int *rowSize_d = (*A).matrixData.rowSize.get();
	int *colOffset_d = (*A).matrixData.colOffset.get();
	int colorCount = (*A).matrixColors.colorCount;
	int *colors_d = (*A).matrixColors.colors_d.get();
	int *colorsColOffsetSize_d = (*A).matrixColors.colorsColOffsetSize_d.get();

	it = 0;

	// Initial x value
	x0->copyTo(x);

	// multiplicacao matriz vetor e subtracao (r = b - A * x)
	// solver triangular inferior e superior - usando apenas o primeiro bloco (cor precisa caber nessas threads)
	// z = inv(M) * r
	int sMemSize = sizeof(double) * blocks;
	dim3 dimGrid(blocks, 1, 1), dimBlock(BLOCKSIZE, 1, 1);
	void *kernelArgs[] = {(void*)&size, &matrixData_d, &preconditionedData_d, &indices_d, &rowLength_d,
		&rowSize_d, &colOffset_d, (void*)&colorCount, &colors_d, &colorsColOffsetSize_d,
		&bData, &xData, &rData, &zData, &partialData, (void*)&blocks};
	cudaLaunchCooperativeKernel((void *)cpcgcg_mult_subtr_solver, dimGrid, dimBlock, kernelArgs, sMemSize);

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcgcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	doIteration0((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);
	doIteration1((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);
	#ifdef CALCULATE_ERRORS
	if (rmod2 < res) { it = 1; return; }
	#endif
	doIteration2((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);
	#ifdef CALCULATE_ERRORS
	if (rmod2 < res) { it = 2; return; }
	#endif
	doIteration3((*A).matrixData.data.get(), (*A).preconditionedData.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(), zData, rData, xData, pData, qData, partialData, partialData2);

#ifdef CALCULATE_ERRORS
	err[0] = wt[0] * wt[0];
	err[1] = w[0] * w[0] + wt[1] * wt[1];
	err[2] = w[0] * w[0] + w[1] * w[1] + wt[2] * wt[2];
#endif
}

void PCGSolverConsolidatedCPJDSCG::doIteration0(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2) {
	double *data_h = new double[1];
	double * rmodData = rmod->getData();
	double * rmod_prevData = rmod_prev->getData();
	double * gammaData = gamma->getData();
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(double), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	cpcgcg_tot_esc_add_mmv_inner0 << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	double *preconditionedData_d = (*A).preconditionedData.get();
	int *indices_d = (*A).matrixData.indices.get();
	int *rowLength_d = (*A).matrixData.rowLength.get();
	int *rowSize_d = (*A).matrixData.rowSize.get();
	int *colOffset_d = (*A).matrixData.colOffset.get();
	int *colors_d = (*A).matrixColors.colors_d.get();
	int *colorsColOffsetSize_d = (*A).matrixColors.colorsColOffsetSize_d.get();
	void *kernelArgs[] = { (void*)&size, &preconditionedData_d, &indices_d, &rowLength_d, &rowSize_d, &colOffset_d, (void*)&colorCount, &colors_d, &colorsColOffsetSize_d,
		&xData, &rData, &zData, &pData, &qData, &rmodData, &rmod_prevData, &gammaData, &partialData2, (void*)&blocks};
	int sMemSize = sizeof(double) * blocks;
	dim3 dimGrid(blocks, 1, 1), dimBlock(BLOCKSIZE, 1, 1);
	cudaLaunchCooperativeKernel((void *)cpcgcg_tot_esc_add_sub_solver0, dimGrid, dimBlock, kernelArgs, sMemSize);
	//cpcgcg_tot_esc_add_sub_solver0 << <blocks, BLOCKSIZE, 0, stream >> >
	//	(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
	//		xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;
	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
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
	cpcgcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

	delete data_h;
}

void PCGSolverConsolidatedCPJDSCG::doIteration1(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2) {
	double *data_h = new double[1];
	double * rmodData = rmod->getData();
	double * rmod_prevData = rmod_prev->getData();
	double * gammaData = gamma->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING
	it++;
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(double), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcgcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	double *preconditionedData_d = (*A).preconditionedData.get();
	int *indices_d = (*A).matrixData.indices.get();
	int *rowLength_d = (*A).matrixData.rowLength.get();
	int *rowSize_d = (*A).matrixData.rowSize.get();
	int *colOffset_d = (*A).matrixData.colOffset.get();
	int *colors_d = (*A).matrixColors.colors_d.get();
	int *colorsColOffsetSize_d = (*A).matrixColors.colorsColOffsetSize_d.get();
	void *kernelArgs[] = { (void*)&size, &preconditionedData_d, &indices_d, &rowLength_d, &rowSize_d, &colOffset_d, (void*)&colorCount, &colors_d, &colorsColOffsetSize_d,
		&xData, &rData, &zData, &pData, &qData, &rmodData, &rmod_prevData, &gammaData, &partialData2, (void*)&blocks };
	int sMemSize = sizeof(double) * blocks;
	dim3 dimGrid(blocks, 1, 1), dimBlock(BLOCKSIZE, 1, 1);
	checkedCudaEventRecord(startTri);
	cudaLaunchCooperativeKernel((void *)cpcgcg_tot_esc_add_sub_solver, dimGrid, dimBlock, kernelArgs, sMemSize);
	//cpcgcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
	//	(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
	//		xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
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
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcgcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

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

void PCGSolverConsolidatedCPJDSCG::doIteration2(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2) {
	double *data_h = new double[1];
	double * rmodData = rmod->getData();
	double * rmod_prevData = rmod_prev->getData();
	double * gammaData = gamma->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING
	it++;
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(double), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcgcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	double *preconditionedData_d = (*A).preconditionedData.get();
	int *indices_d = (*A).matrixData.indices.get();
	int *rowLength_d = (*A).matrixData.rowLength.get();
	int *rowSize_d = (*A).matrixData.rowSize.get();
	int *colOffset_d = (*A).matrixData.colOffset.get();
	int *colors_d = (*A).matrixColors.colors_d.get();
	int *colorsColOffsetSize_d = (*A).matrixColors.colorsColOffsetSize_d.get();
	void *kernelArgs[] = { (void*)&size, &preconditionedData_d, &indices_d, &rowLength_d, &rowSize_d, &colOffset_d, (void*)&colorCount, &colors_d, &colorsColOffsetSize_d,
		&xData, &rData, &zData, &pData, &qData, &rmodData, &rmod_prevData, &gammaData, &partialData2, (void*)&blocks };
	int sMemSize = sizeof(double) * blocks;
	dim3 dimGrid(blocks, 1, 1), dimBlock(BLOCKSIZE, 1, 1);
	checkedCudaEventRecord(startTri);
	cudaLaunchCooperativeKernel((void *)cpcgcg_tot_esc_add_sub_solver, dimGrid, dimBlock, kernelArgs, sMemSize);
	//cpcgcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
	//	(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
	//		xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
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
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcgcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

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

void PCGSolverConsolidatedCPJDSCG::doIteration3(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2) {
	double *data_h = new double[1];
	double * rmodData = rmod->getData();
	double * rmod_prevData = rmod_prev->getData();
	double * gammaData = gamma->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING
	it++;
	cudaMemcpy(x_1->getData(), xData, (size_t)size * sizeof(double), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcgcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, aData, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	double *preconditionedData_d = (*A).preconditionedData.get();
	int *indices_d = (*A).matrixData.indices.get();
	int *rowLength_d = (*A).matrixData.rowLength.get();
	int *rowSize_d = (*A).matrixData.rowSize.get();
	int *colOffset_d = (*A).matrixData.colOffset.get();
	int *colors_d = (*A).matrixColors.colors_d.get();
	int *colorsColOffsetSize_d = (*A).matrixColors.colorsColOffsetSize_d.get();
	void *kernelArgs[] = { (void*)&size, &preconditionedData_d, &indices_d, &rowLength_d, &rowSize_d, &colOffset_d, (void*)&colorCount, &colors_d, &colorsColOffsetSize_d,
		&xData, &rData, &zData, &pData, &qData, &rmodData, &rmod_prevData, &gammaData, &partialData2, (void*)&blocks };
	int sMemSize = sizeof(double) * blocks;
	dim3 dimGrid(blocks, 1, 1), dimBlock(BLOCKSIZE, 1, 1);
	checkedCudaEventRecord(startTri);
	cudaLaunchCooperativeKernel((void *)cpcgcg_tot_esc_add_sub_solver, dimGrid, dimBlock, kernelArgs, sMemSize);
	//cpcgcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
	//	(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
	//		xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
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
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
	#endif

	std::swap(rmod_prev, rmod); // rmod ja foi calculado e usado (para calcular alpha, no kernel anterior), pode-se setar como antigo

	// produto interno (rmod = zt . r, somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	cpcgcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);

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

void PCGSolverConsolidatedCPJDSCG::doIteration(int iteration) {
	double *data_h = new double[1];

	//double * aData = (*A).matrixData.data;
	//double * precond = (*A).preconditionedData;
	//int * aIndices = (*A).matrixData.indices;
	//int * aRowLength = (*A).matrixData.rowLength;
	//int * aRowSize = (*A).matrixData.rowSize;
	//int * aColOffset = (*A).matrixData.colOffset;

	//int colorCount = (*A).matrixColors.colorCount;
	//int * colors = (*A).matrixColors.colors_d;
	//int * colorsColOffset = (*A).matrixColors.colorsColOffsetSize_d;
	double * zData = z->getData();
	double * rData = r->getData();
	double * xData = x->getData();
	double * x2Data = x_1->getData();
	double * pData = p->getData();
	double * qData = q->getData();

	double * rmodData = rmod->getData();
	double * rmod_prevData = rmod_prev->getData();
	double * gammaData = gamma->getData();

	double * partialData = partial->getData();
	double * partialData2 = partial2->getData();

	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
	checkedCudaEventRecord(startTotal);
#endif // CGTIMING

	it++;
	cudaMemcpy(x2Data, xData, (size_t)size * sizeof(double), cudaMemcpyDeviceToDevice);

	// totalizacao (intra-blocos, rmod = zt.r)
	// escalar e vetor soma (p = z + (rmod/rmod_prev) * p)
	// matriz-vetor (q = A * p)
	// produto interno (gamma = pt . (A * p), somente ate totalizacao inter-blocos), precisa ser sincronizado (feito entre kernels)
	checkedCudaEventRecord(startSpmv);
	cpcgcg_tot_esc_add_mmv_inner << <blocks, BLOCKSIZE, 0, stream >> >
		(size, (*A).matrixData.data.get(), (*A).matrixData.indices.get(), (*A).matrixData.rowLength.get(), (*A).matrixData.rowSize.get(), (*A).matrixData.colOffset.get(), (*A).matrixColors.colorCount, (*A).matrixColors.colors_d.get(), (*A).matrixColors.colorsColOffsetSize_d.get(),
			zData, pData, qData, rmodData, rmod_prevData, partialData, partialData2, blocks);
	checkedCudaEventRecord(stopSpmv);

	// totalizacao (intra-bloco, gamma = pt.(A*p))
	// escalar e vetor soma (x += (rmod/gamma) * p)
	// escalar e vetor subtracao (r -= (rmod/gamma) * q)
	// solver triangular inferior (z = inv(M) * r), precisa ser sincronizado (feito entre kernels)
	// solver triangular superior (z = inv(M) * r)
	// precondData: vetor de dados do precondicionador (estrutura identica a matriz completa)
	double *preconditionedData_d = (*A).preconditionedData.get();
	int *indices_d = (*A).matrixData.indices.get();
	int *rowLength_d = (*A).matrixData.rowLength.get();
	int *rowSize_d = (*A).matrixData.rowSize.get();
	int *colOffset_d = (*A).matrixData.colOffset.get();
	int colorCount = (*A).matrixColors.colorCount;
	int *colors_d = (*A).matrixColors.colors_d.get();
	int *colorsColOffsetSize_d = (*A).matrixColors.colorsColOffsetSize_d.get();
	void *kernelArgs[] = { (void*)&size, &preconditionedData_d, &indices_d, &rowLength_d, &rowSize_d, &colOffset_d, (void*)&colorCount, &colors_d, &colorsColOffsetSize_d,
		&xData, &rData, &zData, &pData, &qData, &rmodData, &rmod_prevData, &gammaData, &partialData2, (void*)&blocks };
	int sMemSize = sizeof(double) * blocks;
	dim3 dimGrid(blocks, 1, 1), dimBlock(BLOCKSIZE, 1, 1);
	checkedCudaEventRecord(startTri);
	cudaLaunchCooperativeKernel((void *)cpcgcg_tot_esc_add_sub_solver, dimGrid, dimBlock, kernelArgs, sMemSize);
	//cpcgcg_tot_esc_add_sub_solver << <blocks, BLOCKSIZE, 0, stream >> >
	//	(size, precond, aIndices, aRowLength, aRowSize, aColOffset, colorCount, colors, colorsColOffset,
	//		xData, rData, zData, pData, qData, rmodData, rmod_prevData, gammaData, partialData2, blocks);
	checkedCudaEventRecord(stopTri);

	#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmodData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
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
	cudaMemcpy(data_h, gammaData, (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
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
	cpcgcg_inner << <blocks, BLOCKSIZE, 0, stream >> >(size, rData, zData, partialData, blocks);
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

#endif