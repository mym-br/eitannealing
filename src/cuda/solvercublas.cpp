#include "solvercublas.h"

#include <Eigen/Core>
#include "../basematrix.h"

CGCUBLAS_Solver::CGCUBLAS_Solver(Eigen::SparseMatrix<float, Eigen::ColMajor> *A, float *b) {
	cudaInitialize(A->valuePtr(), A->outerSize(), A->outerSize(), A->nonZeros(), A->outerIndexPtr(), A->innerIndexPtr(), b);
}


CGCUBLAS_Matrix *CGCUBLAS_Matrix::createCublasMatrix(Eigen::SparseMatrix<float, 0, int> *A) {
	CGCUBLAS_Matrix *mat = new CGCUBLAS_Matrix;
	mat->N = A->outerSize();
	mat->nz = A->nonZeros();
	mat->val = A->valuePtr();
	mat->I = A->outerIndexPtr();
	mat->J = A->innerIndexPtr();
	CGCUBLAS_Matrix::cudaMemcpyCublasMatrix(mat);
	return mat;
}