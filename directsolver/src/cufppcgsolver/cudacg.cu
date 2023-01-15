#include "cudacg.h"
// CUDA Runtime
#include <cuda_runtime.h>

// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cusparse.h>
#include <cublas_v2.h>

// Utilities and system includes
//#include "cudasamples/helper_functions.h"
//#include "cudasamples/helper_cuda.h"
#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
// Market Matrix I/O
extern "C" {
#include "mm/mmio.h"
}

std::tuple<long, long, int> runCusparseCublasCG(std::vector<int> &I, std::vector<int> &J, std::vector<double> &val, std::vector<double> &rhs, std::vector<double> &x, int M, int N, int nz, double tol, int max_iter) {
	int *d_col, *d_row;
	double r0, r1, alpha, beta;
	double *d_val, *d_x;
	double *d_zm1, *d_zm2, *d_rm2;
	double *d_r, *d_p, *d_omega, *d_y;
	double *d_valsILU0;
	void *buffer = NULL;
	double dot, numerator, denominator, nalpha;
	const double floatone = 1.0;
	const double floatzero = 0.0;

	int nErrors = 0;

	auto t1 = std::chrono::high_resolution_clock::now();
	/* Create CUBLAS context */
	cublasHandle_t cublasHandle = 0;
	cublasStatus_t cublasStatus;
	cublasStatus = cublasCreate(&cublasHandle);

	/* Create CUSPARSE context */
	cusparseHandle_t cusparseHandle = 0;
	cusparseStatus_t cusparseStatus;
	cusparseStatus = cusparseCreate(&cusparseHandle);

	/* Description of the A matrix*/
	cusparseMatDescr_t descr = 0;
	cusparseStatus = cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	/* Allocate required memory */
	cudaMalloc((void **)&d_col, nz * sizeof(int));
	cudaMalloc((void **)&d_row, (N + 1) * sizeof(int));
	cudaMalloc((void **)&d_val, nz * sizeof(double));
	cudaMalloc((void **)&d_x, N * sizeof(double));
	cudaMalloc((void **)&d_y, N * sizeof(double));
	cudaMalloc((void **)&d_r, N * sizeof(double));
	cudaMalloc((void **)&d_p, N * sizeof(double));
	cudaMalloc((void **)&d_omega, N * sizeof(double));
	cudaMalloc((void **)&d_valsILU0, nz * sizeof(double));
	cudaMalloc((void **)&d_zm1, (N) * sizeof(double));
	cudaMalloc((void **)&d_zm2, (N) * sizeof(double));
	cudaMalloc((void **)&d_rm2, (N) * sizeof(double));

	/* Wrap raw data into cuSPARSE generic API objects */
	cusparseSpMatDescr_t matA = NULL;
	cusparseStatus = cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val,
									CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
									CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
	cusparseDnVecDescr_t vecp = NULL;
	cusparseStatus = cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_64F);
	cusparseDnVecDescr_t vecomega = NULL;
	cusparseStatus = cusparseCreateDnVec(&vecomega, N, d_omega, CUDA_R_64F);

	/* Initialize problem data */
	cudaMemcpy(d_col, J.data(), nz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, I.data(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_val, val.data(), nz * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x.data(), N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, rhs.data(), N * sizeof(double), cudaMemcpyHostToDevice);

	/* Create ILU(0) info object */
	csrilu02Info_t infoILU = NULL;
	cusparseStatus = cusparseCreateCsrilu02Info(&infoILU);

	/* Create L factor descriptor and triangular solve info */
	cusparseMatDescr_t descrL = NULL;
	cusparseStatus = cusparseCreateMatDescr(&descrL);
	cusparseStatus = cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseStatus = cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
	cusparseStatus = cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
	cusparseStatus = cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);
	csrsv2Info_t infoL = NULL;
	cusparseStatus = cusparseCreateCsrsv2Info(&infoL);
  
	/* Create U factor descriptor and triangular solve info */
	cusparseMatDescr_t descrU = NULL;
	cusparseStatus = cusparseCreateMatDescr(&descrU);
	cusparseStatus = cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseStatus = cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
	cusparseStatus = cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
	cusparseStatus = cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	csrsv2Info_t infoU = NULL;
	cusparseStatus = cusparseCreateCsrsv2Info(&infoU);
  
	/* Allocate workspace for cuSPARSE */
	size_t bufferSize = 0;
	size_t tmp = 0;
	int stmp = 0;
	cusparseStatus = cusparseSpMV_bufferSize(
		cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA, vecp,
		&floatzero, vecomega, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &tmp);
	if (tmp > bufferSize) {
	  bufferSize = stmp;
	}
	cusparseStatus = cusparseDcsrilu02_bufferSize(
		cusparseHandle, N, nz, descr, d_val, d_row, d_col, infoILU, &stmp);
	if (stmp > bufferSize) {
	  bufferSize = stmp;
	}
	cusparseStatus = cusparseDcsrsv2_bufferSize(
		cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrL, d_val,
		d_row, d_col, infoL, &stmp);
	if (stmp > bufferSize) {
	  bufferSize = stmp;
	}
	cusparseStatus = cusparseDcsrsv2_bufferSize(
		cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val,
		d_row, d_col, infoU, &stmp);
	if (stmp > bufferSize) {
	  bufferSize = stmp;
	}
	cudaMalloc(&buffer, bufferSize);

	/* Preconditioned Conjugate Gradient using ILU.
	--------------------------------------------
	Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Algorithm 10.3.1  */

	/* Perform analysis for ILU(0) */
	cusparseStatus = cusparseDcsrilu02_analysis(
		cusparseHandle, N, nz, descr, d_val, d_row, d_col, infoILU,
		CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);
  
	/* Copy A data to ILU(0) vals as input*/
	cudaMemcpy(d_valsILU0, d_val, nz * sizeof(double),
							   cudaMemcpyDeviceToDevice);
  
	/* generate the ILU(0) factors */
	cusparseStatus = cusparseDcsrilu02(cusparseHandle, N, nz, descr, d_valsILU0,
									  d_row, d_col, infoILU,
									  CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_analyser = t2 - t1;
	std::cout << "Cublas analyser on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count() << " us." << std::endl;
  
	/* perform triangular solve analysis */
	cusparseStatus = cusparseDcsrsv2_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
								N, nz, descrL, d_valsILU0, d_row, d_col, infoL,
								CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	cusparseStatus = cusparseDcsrsv2_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
								N, nz, descrU, d_valsILU0, d_row, d_col, infoU,
								CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

	/* reset the initial guess of the solution to zero */
	for (int i = 0; i < N; i++) x[i] = 0.0;

	cudaMemcpy(d_r, rhs.data(), N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x.data(), N * sizeof(double), cudaMemcpyHostToDevice);

	t1 = std::chrono::high_resolution_clock::now();
	int k = 0;
	cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

#ifdef CGTIMING
	double totalItTime, totalTriangularTime, totalSpmvTime;
	totalItTime = totalTriangularTime = totalSpmvTime = 0;
	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
#endif // CGTIMING

	double tolsqr = tol > 0 ? tol*tol : -1;
	while (r1 > tolsqr && k <= max_iter)
	{
#ifdef CGTIMING
		cudaEventRecord(startTotal);
		cudaEventRecord(startTri);
#endif // CGTIMING
		// preconditioner application: d_zm1 = U^-1 L^-1 d_r
		cusparseStatus = cusparseDcsrsv2_solve(
			cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, &floatone,
			descrL, d_valsILU0, d_row, d_col, infoL, d_r, d_y,
			CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);
		cusparseStatus = cusparseDcsrsv2_solve(
			cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, &floatone,
			descrU, d_valsILU0, d_row, d_col, infoU, d_y, d_zm1,
			CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);
#ifdef CGTIMING
		cudaEventRecord(stopTri);
#endif // CGTIMING

		k++;

		if (k == 1)
		{
			cublasDcopy(cublasHandle, N, d_zm1, 1, d_p, 1);
		}
		else
		{
			cublasDdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
			cublasDdot(cublasHandle, N, d_rm2, 1, d_zm2, 1, &denominator);
			beta = numerator / denominator;
			cublasDscal(cublasHandle, N, &beta, d_p, 1);
			cublasDaxpy(cublasHandle, N, &floatone, d_zm1, 1, d_p, 1);
		}
#ifdef CGTIMING
		cudaEventRecord(startSpmv);
#endif // CGTIMING
		cusparseSpMV(
			cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA, vecp,
			&floatzero, vecomega, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, buffer);
#ifdef CGTIMING
		cudaEventRecord(stopSpmv);
#endif // CGTIMING
		cublasDdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
		cublasDdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
		alpha = numerator / denominator;
		cublasDaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		cublasDcopy(cublasHandle, N, d_r, 1, d_rm2, 1);
		cublasDcopy(cublasHandle, N, d_zm1, 1, d_zm2, 1);
		nalpha = -alpha;
		cublasDaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
#ifdef CGTIMING
		cudaEventRecord(stopTotal);
#endif // CGTIMING
#ifdef CGTIMING
		cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
		float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
		cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
		totalItTime += (double)(1e3 * msTotal);
		totalTriangularTime += (double)(1e3 * msTri);
		totalSpmvTime += (double)(1e3 *  msSpmv);
#endif // CGTIMING
	}
	cudaMemcpy(x.data(), d_x, N * sizeof(double), cudaMemcpyDeviceToHost); 
	t2 = std::chrono::high_resolution_clock::now();

	///************************/
	///* now write out result */
	///************************/
	std::chrono::duration<double> time_executor = t2 - t1;
	std::cout << "Cublas executor on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count() << " us. Final residual is " << sqrt(r1) << " after " << k << " iterations." << std::endl;

#ifdef CGTIMING
	totalItTime /= (double)k; totalTriangularTime /= (double)k; totalSpmvTime /= (double)k;
	std::cout << "Average cublas/cusparse iteration time breakdown: " << totalTriangularTime << " (triangular solver) " << totalSpmvTime << " (spmv) " << totalItTime - totalTriangularTime - totalSpmvTime << " (remaining) " << totalItTime << " (total)." << std::endl;
#endif // CGTIMING

	/* Destroy descriptors */
	cusparseDestroyCsrsv2Info(infoU);
	cusparseDestroyCsrsv2Info(infoL);
	cusparseDestroyCsrilu02Info(infoILU);
	cusparseDestroyMatDescr(descrL);
	cusparseDestroyMatDescr(descrU);
	cusparseDestroyMatDescr(descr);
	cusparseDestroySpMat(matA);
	cusparseDestroyDnVec(vecp);
	cusparseDestroyDnVec(vecomega);

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

	return std::make_tuple(std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count(), std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count(), k);
}