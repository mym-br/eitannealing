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

namespace Cublas {

	class CusparseHandle {
	public:
		static CusparseHandle& Instance() {
			// Since it's a static variable, if the class has already been created,
			// it won't be created again.
			// And it **is** thread-safe in C++11.
			static CusparseHandle myInstance;

			// Return a reference to our instance.
			return myInstance;
		}

		// delete copy and move constructors and assign operators
		CusparseHandle(CusparseHandle const&) = delete;             // Copy construct
		CusparseHandle(CusparseHandle&&) = delete;                  // Move construct
		CusparseHandle& operator=(CusparseHandle const&) = delete;  // Copy assign
		CusparseHandle& operator=(CusparseHandle &&) = delete;      // Move assign
		cusparseHandle_t getHandle() { return hdl; }

	protected:
		cusparseHandle_t hdl;
		CusparseHandle();
		~CusparseHandle();
	};

	class CublasHandle {
	public:
		static CublasHandle& Instance() {
			// Since it's a static variable, if the class has already been created,
			// it won't be created again.
			// And it **is** thread-safe in C++11.
			static CublasHandle myInstance;

			// Return a reference to our instance.
			return myInstance;
		}

		// delete copy and move constructors and assign operators
		CublasHandle(CublasHandle const&) = delete;             // Copy construct
		CublasHandle(CublasHandle&&) = delete;                  // Move construct
		CublasHandle& operator=(CublasHandle const&) = delete;  // Copy assign
		CublasHandle& operator=(CublasHandle &&) = delete;      // Move assign
		cublasHandle_t getHandle() { return hdl; }

	protected:
		cublasHandle_t hdl;
		CublasHandle();
		~CublasHandle();
	};

	struct Matrix {
		~Matrix();
		static Matrix *createCublasMatrix(Eigen::SparseMatrix<float, 0, int> *A);
		cusparseMatDescr_t descr;
		float *val, *d_val;
		int *I, *J, *d_col, *d_row;
		int N, nz;

	private:
		static void cudaMemcpyCublasMatrix(Matrix *A);
	};

	struct Precond {
		~Precond();
		static Precond *createPrecond(Matrix *A);
		float *d_valsILU0;
		int nzILU0;
		cusparseSolveAnalysisInfo_t infoA;
		cusparseSolveAnalysisInfo_t info_u;
		cusparseMatDescr_t descrL;
		cusparseMatDescr_t descrU;
	};

	class CG_Solver {
	public:
		CG_Solver(Matrix *_A, float *_b, Precond *_precond);
		~CG_Solver();
		void doIteration();
		float *getX();

	private:
		Matrix *A;
		Precond *precond;
		float *x;
		int k;
		float r1;
		float *d_x;
		float *d_r, *d_y;
		float *d_zm1, *d_zm2, *d_rm2;
		float *d_p, *d_omega;
	};
}
#endif // SOLVERCUBLAS_H_
