#ifndef SOLVERCUBLAS_H_
#define SOLVERCUBLAS_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>

// CUDA Runtime
#include <cuda_runtime.h>

// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cusparse.h>
#include <cublas_v2.h>

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
		static Matrix *createCublasMatrix(Eigen::SparseMatrix<double> *A);
		cusparseMatDescr_t descr;
		double *val, *d_val;
		int *I, *J, *d_col, *d_row;
		int N, nz;

	private:
		static void cudaMemcpyCublasMatrix(Matrix *A);
	};

	struct Precond {
		~Precond();
		static Precond *createPrecond(Matrix *A);
		double *d_valsILU0;
		csrilu02Info_t infoILU;
		cusparseMatDescr_t descrL;
		cusparseMatDescr_t descrU;
		csrsv2Info_t infoU;
		csrsv2Info_t infoL;
	};

	class CG_Solver {
	public:
		CG_Solver(Matrix *_A, double *_b, Precond *_precond);
		~CG_Solver();
		void doIteration();
		double *getX();

	private:
		Matrix *A;
		cusparseSpMatDescr_t matA;
		cusparseDnVecDescr_t vecp, vecomega;
		Precond *precond;
		double *x;
		int k;
		double r1;
		double *d_x;
		double *d_r, *d_y;
		double *d_zm1, *d_zm2, *d_rm2;
		double *d_p, *d_omega;
		void *buffer;
	};
}
#endif // SOLVERCUBLAS_H_
