#include "solvercublas.h"
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

CGCUBLAS_Solver::~CGCUBLAS_Solver() {
	/* Destroy parameters */
	cusparseDestroySolveAnalysisInfo(infoA);
	cusparseDestroySolveAnalysisInfo(info_u);

	/* Destroy contexts */
	cusparseDestroy(cusparseHandle);
	cublasDestroy(cublasHandle);

	/* Free device memory */
	cudaFree(d_col);
	cudaFree(d_row);
	cudaFree(d_val);
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_omega);
	cudaFree(d_valsILU0);
	cudaFree(d_zm1);
	cudaFree(d_zm2);
	cudaFree(d_rm2);
}

void CGCUBLAS_Solver::cudaInitialize(float *val, int M, int N, int nz, int *I, int *J, float *rhs) {
	this->M = M; this->N = N; this->nz = nz;
	x = (float *)malloc(sizeof(float)*N);
	for (int i = 0; i < N; i++) x[i] = 0.0;

	/* Create CUBLAS context */
	cublasHandle = 0;
	cublasStatus_t cublasStatus;
	cublasStatus = cublasCreate(&cublasHandle);

	/* Create CUSPARSE context */
	cusparseHandle = 0;
	cusparseStatus_t cusparseStatus;
	cusparseStatus = cusparseCreate(&cusparseHandle);

	/* Description of the A matrix*/
	descr = 0;
	cusparseStatus = cusparseCreateMatDescr(&descr);

	/* Define the properties of the matrix */
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	/* Allocate required memory */
	cudaMalloc((void **)&d_col, nz * sizeof(int));
	cudaMalloc((void **)&d_row, (N + 1) * sizeof(int));
	cudaMalloc((void **)&d_val, nz * sizeof(float));
	cudaMalloc((void **)&d_x, N * sizeof(float));
	cudaMalloc((void **)&d_y, N * sizeof(float));
	cudaMalloc((void **)&d_r, N * sizeof(float));
	cudaMalloc((void **)&d_p, N * sizeof(float));
	cudaMalloc((void **)&d_omega, N * sizeof(float));

	cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice);

	nzILU0 = 2 * N - 1;
	cudaMalloc((void **)&d_valsILU0, nz * sizeof(float));
	cudaMalloc((void **)&d_zm1, (N) * sizeof(float));
	cudaMalloc((void **)&d_zm2, (N) * sizeof(float));
	cudaMalloc((void **)&d_rm2, (N) * sizeof(float));
}

void CGCUBLAS_Solver::calculatePrecond() {
	cusparseStatus_t cusparseStatus;
	/* create the analysis info object for the A matrix */
	infoA = 0;
	cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);

	//checkCudaErrors(cusparseStatus);

	/* Perform the analysis for the Non-Transpose case */
	cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		N, nz, descr, d_val, d_row, d_col, infoA);

	//checkCudaErrors(cusparseStatus);

	/* Copy A data to ILU0 vals as input*/
	cudaMemcpy(d_valsILU0, d_val, nz * sizeof(float), cudaMemcpyDeviceToDevice);

	/* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
	cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsILU0, d_row, d_col, infoA);

	//checkCudaErrors(cusparseStatus);

	/* Create info objects for the ILU0 preconditioner */
	cusparseCreateSolveAnalysisInfo(&info_u);

	descrL = 0;
	cusparseStatus = cusparseCreateMatDescr(&descrL);
	cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
	cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

	descrU = 0;
	cusparseStatus = cusparseCreateMatDescr(&descrU);
	cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
	cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val, d_row, d_col, info_u);

	k = 0;
	cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
}

void CGCUBLAS_Solver::doIteration() {
	const float tol = 1e-12f;
	const int max_iter = 1000;
	const float floatone = 1.0;
	const float floatzero = 0.0;
	float alpha, beta;
	float numerator, denominator, nalpha;
	cusparseStatus_t cusparseStatus;

	// Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
	cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, descrL,
		d_valsILU0, d_row, d_col, infoA, d_r, d_y);
	//checkCudaErrors(cusparseStatus);

	// Back Substitution
	cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, descrU,
		d_valsILU0, d_row, d_col, info_u, d_y, d_zm1);
	//checkCudaErrors(cusparseStatus);

	k++;

	if (k == 1)
	{
		cublasScopy(cublasHandle, N, d_zm1, 1, d_p, 1);
	}
	else
	{
		cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
		cublasSdot(cublasHandle, N, d_rm2, 1, d_zm2, 1, &denominator);
		beta = numerator / denominator;
		cublasSscal(cublasHandle, N, &beta, d_p, 1);
		cublasSaxpy(cublasHandle, N, &floatone, d_zm1, 1, d_p, 1);
	}

	cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nzILU0, &floatone, descrU, d_val, d_row, d_col, d_p, &floatzero, d_omega);
	cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
	cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
	alpha = numerator / denominator;
	cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
	cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
	cublasScopy(cublasHandle, N, d_zm1, 1, d_zm2, 1);
	nalpha = -alpha;
	cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
	cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
}

float *CGCUBLAS_Solver::getX() {
	cudaMemcpy(x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);
	return x;
}


CGCUBLAS_Precond *CGCUBLAS_Precond::createPrecond(Eigen::SparseMatrix<float, 0, int> *A) {
	CGCUBLAS_Precond *precond = new CGCUBLAS_Precond;

	//cusparseStatus_t cusparseStatus;
	///* create the analysis info object for the A matrix */
	//cusparseSolveAnalysisInfo_t infoA = 0;
	//cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);

	////checkCudaErrors(cusparseStatus);

	///* Perform the analysis for the Non-Transpose case */
	//cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
	//	N, nz, descr, d_val, d_row, d_col, infoA);

	////checkCudaErrors(cusparseStatus);

	///* Copy A data to ILU0 vals as input*/
	//cudaMemcpy(d_valsILU0, d_val, nz * sizeof(float), cudaMemcpyDeviceToDevice);

	///* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
	//cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsILU0, d_row, d_col, infoA);

	////checkCudaErrors(cusparseStatus);

	///* Create info objects for the ILU0 preconditioner */
	//cusparseCreateSolveAnalysisInfo(&info_u);

	//descrL = 0;
	//cusparseStatus = cusparseCreateMatDescr(&descrL);
	//cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
	//cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
	//cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
	//cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

	//descrU = 0;
	//cusparseStatus = cusparseCreateMatDescr(&descrU);
	//cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
	//cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
	//cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
	//cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	//cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val, d_row, d_col, info_u);

	return precond;
}

void CGCUBLAS_Matrix::cudaMemcpyCublasMatrix(CGCUBLAS_Matrix *A) {
	A->descr = 0;
	cusparseCreateMatDescr(&A->descr);

	/* Define the properties of the matrix */
	cusparseSetMatType(A->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(A->descr, CUSPARSE_INDEX_BASE_ZERO);

	/* Allocate required memory */
	int N = A->N, nz = A->nz;
	cudaMalloc((void **)&A->d_row, (N + 1) * sizeof(int));
	cudaMalloc((void **)&A->d_val, nz * sizeof(float));
	cudaMemcpy(A->d_col, A->J, nz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(A->d_row, A->I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(A->d_val, A->val, nz * sizeof(float), cudaMemcpyHostToDevice);

}