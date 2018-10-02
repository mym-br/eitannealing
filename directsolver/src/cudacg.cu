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

std::tuple<long, long> runCusparseCublasCG(std::vector<int> &I, std::vector<int> &J, std::vector<float> &val, std::vector<float> &rhs, std::vector<float> &x, int M, int N, int nz, float tol, int max_iter) {
	int *d_col, *d_row;
	float r0, r1, alpha, beta;
	float *d_val, *d_x;
	float *d_zm1, *d_zm2, *d_rm2;
	float *d_r, *d_p, *d_omega, *d_y;
	float *d_valsILU0;
	float *valsILU0;
	float dot, numerator, denominator, nalpha;
	const float floatone = 1.0;
	const float floatzero = 0.0;

	int nErrors = 0;

	auto t1 = std::chrono::high_resolution_clock::now();
	/* Create CUBLAS context */
	cublasHandle_t cublasHandle = 0;
	cublasStatus_t cublasStatus;
	cublasStatus = cublasCreate(&cublasHandle);

	cublasStatus;

	/* Create CUSPARSE context */
	cusparseHandle_t cusparseHandle = 0;
	cusparseStatus_t cusparseStatus;
	cusparseStatus = cusparseCreate(&cusparseHandle);

	cusparseStatus;

	/* Description of the A matrix*/
	cusparseMatDescr_t descr = 0;
	cusparseStatus = cusparseCreateMatDescr(&descr);

	cusparseStatus;

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

	cudaMemcpy(d_col, J.data(), nz * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, I.data(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_val, val.data(), nz * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x.data(), N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, rhs.data(), N * sizeof(float), cudaMemcpyHostToDevice);

	/* Preconditioned Conjugate Gradient using ILU.
	--------------------------------------------
	Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Algorithm 10.3.1  */

	int nzILU0 = 2 * N - 1;
	valsILU0 = (float *)malloc(nz * sizeof(float));

	cudaMalloc((void **)&d_valsILU0, nz * sizeof(float));
	cudaMalloc((void **)&d_zm1, (N) * sizeof(float));
	cudaMalloc((void **)&d_zm2, (N) * sizeof(float));
	cudaMalloc((void **)&d_rm2, (N) * sizeof(float));

	/* create the analysis info object for the A matrix */
	cusparseSolveAnalysisInfo_t infoA = 0;
	cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);

	/* Perform the analysis for the Non-Transpose case */
	cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		N, nz, descr, d_val, d_row, d_col, infoA);

	/* Copy A data to ILU0 vals as input*/
	cudaMemcpy(d_valsILU0, d_val, nz * sizeof(float), cudaMemcpyDeviceToDevice);

	/* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
	cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsILU0, d_row, d_col, infoA);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_analyser = t2 - t1;
	std::cout << "Cublas analyser on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count() << " us." << std::endl;

	/* Create info objects for the ILU0 preconditioner */
	cusparseSolveAnalysisInfo_t info_u;
	cusparseCreateSolveAnalysisInfo(&info_u);

	cusparseMatDescr_t descrL = 0;
	cusparseStatus = cusparseCreateMatDescr(&descrL);
	cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
	cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

	cusparseMatDescr_t descrU = 0;
	cusparseStatus = cusparseCreateMatDescr(&descrU);
	cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
	cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val, d_row, d_col, info_u);

	/* reset the initial guess of the solution to zero */
	for (int i = 0; i < N; i++) x[i] = 0.0;

	cudaMemcpy(d_r, rhs.data(), N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x.data(), N * sizeof(float), cudaMemcpyHostToDevice);

	t1 = std::chrono::high_resolution_clock::now();
	int k = 0;
	cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

#ifdef CGTIMING
	double totalItTime, totalTriangularTime, totalSpmvTime;
	totalItTime = totalTriangularTime = totalSpmvTime = 0;
	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
#endif // CGTIMING

	while (r1 > tol*tol && k <= max_iter)
	{
#ifdef CGTIMING
		cudaEventRecord(startTotal);
		cudaEventRecord(startTri);
#endif // CGTIMING
		// Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
		cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, descrL,
			d_valsILU0, d_row, d_col, infoA, d_r, d_y);

		// Back Substitution
		cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, descrU,
			d_valsILU0, d_row, d_col, info_u, d_y, d_zm1);
#ifdef CGTIMING
		cudaEventRecord(stopTri);
#endif // CGTIMING

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
#ifdef CGTIMING
		cudaEventRecord(startSpmv);
#endif // CGTIMING
		cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nzILU0, &floatone, descrU, d_val, d_row, d_col, d_p, &floatzero, d_omega);
#ifdef CGTIMING
		cudaEventRecord(stopSpmv);
#endif // CGTIMING
		cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
		cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
		alpha = numerator / denominator;
		cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
		cublasScopy(cublasHandle, N, d_zm1, 1, d_zm2, 1);
		nalpha = -alpha;
		cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
#ifdef CGTIMING
		cudaEventRecord(stopTotal);
#endif // CGTIMING
#ifdef CGTIMING
		cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
		float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
		cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
		totalItTime += (float)(1e3 * msTotal);
		totalTriangularTime += (float)(1e3 * msTri);
		totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING
	}
	cudaMemcpy(x.data(), d_x, N * sizeof(float), cudaMemcpyDeviceToHost); 
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

	/* Destroy parameters */
	cusparseDestroySolveAnalysisInfo(infoA);
	cusparseDestroySolveAnalysisInfo(info_u);

	/* Destroy contexts */
	cusparseDestroy(cusparseHandle);
	cublasDestroy(cublasHandle);

	/* Free device memory */
	free(valsILU0);
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

	return std::make_tuple(std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count(), std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count());
}