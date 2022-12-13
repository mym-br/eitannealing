#include "solvercublas.h"
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

using namespace Cublas;

CG_Solver::~CG_Solver() {
	/* Free device memory */
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_omega);
	cudaFree(d_zm1);
	cudaFree(d_zm2);
	cudaFree(d_rm2);
}

CG_Solver::CG_Solver(Matrix *_A, float *_b, Precond *_precond) : A(_A), precond(_precond), buffer(NULL) {
	// Initialize x
	int N = A->N;
	x = (float *)malloc(sizeof(float)*N);
	for (int i = 0; i < N; i++) x[i] = 0.0;

	/* Allocate required memory */
	cudaMalloc((void **)&d_x, N * sizeof(float));
	cudaMalloc((void **)&d_r, N * sizeof(float));
	cudaMemcpy(d_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, _b, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc((void **)&d_y, N * sizeof(float));
	cudaMalloc((void **)&d_p, N * sizeof(float));
	cudaMalloc((void **)&d_omega, N * sizeof(float));
	cudaMalloc((void **)&d_zm1, (N) * sizeof(float));
	cudaMalloc((void **)&d_zm2, (N) * sizeof(float));
	cudaMalloc((void **)&d_rm2, (N) * sizeof(float));

	/* Wrap raw data into cuSPARSE generic API objects */
	matA = NULL;
	cusparseCreateCsr(&matA, N, N, A->nz, A->d_row, A->d_col, A->d_val,
									CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
									CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
	vecp = NULL;
	cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_32F);
	vecomega = NULL;
	cusparseCreateDnVec(&vecomega, N, d_omega, CUDA_R_32F);
	
	/* Allocate workspace for cuSPARSE */
	const float floatone = 1.0;
	const float floatzero = 0.0;
	size_t bufferSize = 0;
	size_t tmp = 0;
	int stmp = 0;
	cusparseSpMV_bufferSize(
		CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA, vecp,
		&floatzero, vecomega, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &tmp);
	if (tmp > bufferSize) {
		bufferSize = stmp;
	}
	cusparseScsrilu02_bufferSize(
		CusparseHandle::Instance().getHandle(), N, A->nz, A->descr, A->d_val, A->d_row, A->d_col, precond->infoILU, &stmp);
	if (stmp > bufferSize) {
		bufferSize = stmp;
	}
	cusparseScsrsv2_bufferSize(
		CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, A->nz, precond->descrL, A->d_val,
		A->d_row, A->d_col, precond->infoL, &stmp);
	if (stmp > bufferSize) {
	bufferSize = stmp;
	}
	cusparseScsrsv2_bufferSize(
		CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, A->nz, precond->descrU, A->d_val,
		A->d_row, A->d_col, precond->infoU, &stmp);
	if (stmp > bufferSize) {
	bufferSize = stmp;
	}
	cudaMalloc(&buffer, bufferSize);

	/* Perform analysis for ILU(0) */
	cusparseScsrilu02_analysis(
		CusparseHandle::Instance().getHandle(), A->N, A->nz, A->descr, A->d_val, A->d_row, A->d_col, precond->infoILU,
		CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	/* Copy A data to ILU(0) vals as input*/
	cudaMemcpy(precond->d_valsILU0, A->d_val, A->nz * sizeof(float),
							cudaMemcpyDeviceToDevice);

	/* generate the ILU(0) factors */
	cusparseScsrilu02(CusparseHandle::Instance().getHandle(), A->N, A->nz, A->descr, precond->d_valsILU0,
									A->d_row, A->d_col, precond->infoILU,
									CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	/* perform triangular solve analysis */
	cusparseScsrsv2_analysis(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
							 A->N, A->nz, precond->descrL, precond->d_valsILU0, A->d_row, A->d_col, precond->infoL,
							 CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	cusparseScsrsv2_analysis(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
							 A->N, A->nz, precond->descrU, precond->d_valsILU0, A->d_row, A->d_col, precond->infoU,
							 CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	k = 0;
	cublasSdot(CublasHandle::Instance().getHandle(), N, d_r, 1, d_r, 1, &r1);
}

void CG_Solver::doIteration() {
	const float tol = 1e-12f;
	const int max_iter = 1000;
	const float floatone = 1.0;
	const float floatzero = 0.0;
	float alpha, beta;
	float numerator, denominator, nalpha;
	cusparseStatus_t cusparseStatus;
	int N = A->N;

	// preconditioner application: d_zm1 = U^-1 L^-1 d_r
    cusparseStatus = cusparseScsrsv2_solve(
        CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, A->nz, &floatone,
        precond->descrL, precond->d_valsILU0, A->d_row, A->d_col, precond->infoL, d_r, d_y,
        CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);
	//checkCudaErrors(cusparseStatus);
	cusparseStatus = cusparseScsrsv2_solve(
        CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, A->nz, &floatone,
        precond->descrU, precond->d_valsILU0, A->d_row, A->d_col, precond->infoU, d_y, d_zm1,
        CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);
	//checkCudaErrors(cusparseStatus);

	k++;

	if (k == 1)
	{
		cublasScopy(CublasHandle::Instance().getHandle(), N, d_zm1, 1, d_p, 1);
	}
	else
	{
		cublasSdot(CublasHandle::Instance().getHandle(), N, d_r, 1, d_zm1, 1, &numerator);
		cublasSdot(CublasHandle::Instance().getHandle(), N, d_rm2, 1, d_zm2, 1, &denominator);
		beta = numerator / denominator;
		cublasSscal(CublasHandle::Instance().getHandle(), N, &beta, d_p, 1);
		cublasSaxpy(CublasHandle::Instance().getHandle(), N, &floatone, d_zm1, 1, d_p, 1);
	}

	cusparseSpMV(
        CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA, vecp,
        &floatzero, vecomega, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, buffer);
	cublasSdot(CublasHandle::Instance().getHandle(), N, d_r, 1, d_zm1, 1, &numerator);
	cublasSdot(CublasHandle::Instance().getHandle(), N, d_p, 1, d_omega, 1, &denominator);
	alpha = numerator / denominator;
	cublasSaxpy(CublasHandle::Instance().getHandle(), N, &alpha, d_p, 1, d_x, 1);
	cublasScopy(CublasHandle::Instance().getHandle(), N, d_r, 1, d_rm2, 1);
	cublasScopy(CublasHandle::Instance().getHandle(), N, d_zm1, 1, d_zm2, 1);
	nalpha = -alpha;
	cublasSaxpy(CublasHandle::Instance().getHandle(), N, &nalpha, d_omega, 1, d_r, 1);
	cublasSdot(CublasHandle::Instance().getHandle(), N, d_r, 1, d_r, 1, &r1);
}

float *CG_Solver::getX() {
	cudaMemcpy(x, d_x, A->N * sizeof(float), cudaMemcpyDeviceToHost);
	return x;
}

Precond::~Precond() {
	/* Destroy parameters */
	cusparseDestroyCsrilu02Info(infoILU);
	cusparseDestroyMatDescr(descrL);
	cusparseDestroyMatDescr(descrU);
}

Precond *Precond::createPrecond(Matrix *A) {
	Precond *precond = new Precond;
	cusparseStatus_t cusparseStatus;
	
	/* Create ILU(0) info object */
	precond->infoILU = NULL;
	cusparseCreateCsrilu02Info(&precond->infoILU);
	  
	/* Create L factor descriptor and triangular solve info */
	precond->descrL = NULL;
	cusparseStatus = cusparseCreateMatDescr(&precond->descrL);
	cusparseSetMatType(precond->descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(precond->descrL, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(precond->descrL, CUSPARSE_FILL_MODE_LOWER);
	cusparseSetMatDiagType(precond->descrL, CUSPARSE_DIAG_TYPE_UNIT);
	precond->infoL = NULL;
	cusparseCreateCsrsv2Info(&precond->infoL);

	precond->descrU = NULL;
	cusparseStatus = cusparseCreateMatDescr(&precond->descrU);
	cusparseSetMatType(precond->descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(precond->descrU, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(precond->descrU, CUSPARSE_FILL_MODE_UPPER);
	cusparseSetMatDiagType(precond->descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	precond->infoU = NULL;
	cusparseCreateCsrsv2Info(&precond->infoU);

	return precond;
}

Matrix::~Matrix() {
	cudaFree(d_col);
	cudaFree(d_row);
	cudaFree(d_val);
}

void Matrix::cudaMemcpyCublasMatrix(Matrix *A) {
	A->descr = 0;
	cusparseCreateMatDescr(&A->descr);

	/* Define the properties of the matrix */
	cusparseSetMatType(A->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(A->descr, CUSPARSE_INDEX_BASE_ZERO);

	/* Allocate required memory */
	int N = A->N, nz = A->nz;
	cudaMalloc((void **)&A->d_col, nz * sizeof(int));
	cudaMalloc((void **)&A->d_row, (N + 1) * sizeof(int));
	cudaMalloc((void **)&A->d_val, nz * sizeof(float));
	cudaMemcpy(A->d_col, A->J, nz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(A->d_row, A->I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(A->d_val, A->val, nz * sizeof(float), cudaMemcpyHostToDevice);

}

CusparseHandle::CusparseHandle() {
	hdl = 0;
	cusparseStatus_t cusparseStatus = cusparseCreate(&hdl);
}

CusparseHandle::~CusparseHandle() {
	cusparseDestroy(hdl);
}

CublasHandle::CublasHandle() {
	hdl = 0;
	cublasStatus_t cublasStatus = cublasCreate(&hdl);
}

CublasHandle::~CublasHandle() {
	cublasDestroy(hdl);
}
