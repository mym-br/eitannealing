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

CG_Solver::CG_Solver(Matrix *_A, float *_b, Precond *_precond) : A(_A), precond(_precond) {
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

	// Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
	cusparseStatus = cusparseScsrsv_solve(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, precond->descrL,
		precond->d_valsILU0, A->d_row, A->d_col, precond->infoA, d_r, d_y);
	//checkCudaErrors(cusparseStatus);

	// Back Substitution
	cusparseStatus = cusparseScsrsv_solve(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, precond->descrU,
		precond->d_valsILU0, A->d_row, A->d_col, precond->info_u, d_y, d_zm1);
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

	cusparseScsrmv(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, precond->nzILU0, &floatone, precond->descrU, A->d_val, A->d_row, A->d_col, d_p, &floatzero, d_omega);
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
	cusparseDestroySolveAnalysisInfo(infoA);
	cusparseDestroySolveAnalysisInfo(info_u);
	cudaFree(d_valsILU0);
}

Precond *Precond::createPrecond(Matrix *A) {
	Precond *precond = new Precond;

	cusparseStatus_t cusparseStatus;
	/* create the analysis info object for the A matrix */
	precond->infoA = 0;
	cusparseStatus = cusparseCreateSolveAnalysisInfo(&precond->infoA);

	/* Perform the analysis for the Non-Transpose case */
	cusparseStatus = cusparseScsrsv_analysis(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
		A->N, A->nz, A->descr, A->d_val, A->d_row, A->d_col, precond->infoA);

	/* Copy A data to ILU0 vals as input*/
	cudaMalloc((void **)&precond->d_valsILU0, A->nz * sizeof(float));
	cudaMemcpy(precond->d_valsILU0, A->d_val, A->nz * sizeof(float), cudaMemcpyDeviceToDevice);

	/* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
	cusparseStatus = cusparseScsrilu0(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, A->N, A->descr, precond->d_valsILU0, A->d_row, A->d_col, precond->infoA);

	///* Create info objects for the ILU0 preconditioner */
	cusparseCreateSolveAnalysisInfo(&precond->info_u);

	precond->descrL = 0;
	cusparseStatus = cusparseCreateMatDescr(&precond->descrL);
	cusparseSetMatType(precond->descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(precond->descrL, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(precond->descrL, CUSPARSE_FILL_MODE_LOWER);
	cusparseSetMatDiagType(precond->descrL, CUSPARSE_DIAG_TYPE_UNIT);

	precond->descrU = 0;
	cusparseStatus = cusparseCreateMatDescr(&precond->descrU);
	cusparseSetMatType(precond->descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(precond->descrU, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(precond->descrU, CUSPARSE_FILL_MODE_UPPER);
	cusparseSetMatDiagType(precond->descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	cusparseStatus = cusparseScsrsv_analysis(CusparseHandle::Instance().getHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, A->N, A->nz, precond->descrU, A->d_val, A->d_row, A->d_col, precond->info_u);

	precond->nzILU0 = 2 * A->N - 1;

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
