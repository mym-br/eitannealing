#ifndef SOLVERCUBLAS_H_
#define SOLVERCUBLAS_H_

// CUDA Runtime
#include <cuda_runtime.h>

// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cusparse.h>
#include <cublas_v2.h>

namespace Eigen {
	template<typename _Scalar, int _Options, typename _Index> class SparseMatrix;
	//template<typename _Scalar, int _Options, typename _Index> class Matrix;
}

struct CGCUBLAS_Matrix {
	static CGCUBLAS_Matrix *createCublasMatrix(Eigen::SparseMatrix<float, 0, int> *A);
	cusparseMatDescr_t descr;
	float *val, *d_val;
	int *I, *J, *d_col, *d_row;
	int N, nz;

private:
	static void cudaMemcpyCublasMatrix(CGCUBLAS_Matrix *A);
};

struct CGCUBLAS_Precond {
	static CGCUBLAS_Precond *createPrecond(Eigen::SparseMatrix<float, 0, int> *A);
	float *d_valsILU0;
	cusparseSolveAnalysisInfo_t infoA;
	cusparseSolveAnalysisInfo_t info_u;
	cusparseMatDescr_t descrL;
	cusparseMatDescr_t descrU;
};

class CGCUBLAS_Solver {
public:
	CGCUBLAS_Solver(Eigen::SparseMatrix<float, 0, int> *A, float *b);
	~CGCUBLAS_Solver();
	void calculatePrecond();
	void doIteration();
	float *getX();

private:
	void cudaInitialize(float *val, int M, int N, int nz, int *I, int *J, float *rhs);
	float *x;

	//int k, M, N, nz, *I, *J;
	//float *val;
	//const float tol = 1e-12f;
	//float *x, *rhs;

	//int *d_col, *d_row;
	//float *d_val, *d_x;
	//float *d_r, *d_p, *d_omega, *d_y;

	//int nzILU0;
	//float *d_valsILU0;
	//float *valsILU0;
	//float *d_zm1, *d_zm2, *d_rm2;
	//cusparseSolveAnalysisInfo_t infoA = 0;
	//cusparseSolveAnalysisInfo_t info_u;
	//cusparseMatDescr_t descrL = 0;
	//cusparseMatDescr_t descrU = 0;

	//float r0, r1, alpha, beta;
	//float dot, numerator, denominator, nalpha;
	//const float floatone = 1.0;
	//const float floatzero = 0.0;

	///* Create CUBLAS context */
	//cublasHandle_t cublasHandle = 0;
	//cusparseHandle_t cusparseHandle = 0;
	//cusparseMatDescr_t descr = 0;
	cusparseHandle_t cusparseHandle;
	cusparseMatDescr_t descr;
	int M, N, nz;
	float *d_val, *d_x;
	int *d_col, *d_row;
	float *d_valsILU0;
	int k;
	cublasHandle_t cublasHandle;
	float r1;
	float *d_r, *d_y;
	float *d_zm1, *d_zm2, *d_rm2;
	int nzILU0;
	float *d_p, *d_omega;
	cusparseSolveAnalysisInfo_t infoA;
	cusparseSolveAnalysisInfo_t info_u;
	cusparseMatDescr_t descrL;
	cusparseMatDescr_t descrU;
};

#endif // SOLVERCUBLAS_H_
