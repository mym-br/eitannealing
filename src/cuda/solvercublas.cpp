#include "solvercublas.h"

#include <Eigen/Core>
#include "../basematrix.h"

using namespace Cublas;

CG_Solver::CG_Solver(Eigen::SparseMatrix<float, Eigen::ColMajor> *A, float *b) {
	cudaInitialize(A->valuePtr(), A->outerSize(), A->outerSize(), A->nonZeros(), A->outerIndexPtr(), A->innerIndexPtr(), b);
}


Matrix *Matrix::createCublasMatrix(Eigen::SparseMatrix<float, 0, int> *A) {
	Matrix *mat = new Matrix;
	mat->N = A->outerSize();
	mat->nz = A->nonZeros();
	mat->val = A->valuePtr();
	mat->I = A->outerIndexPtr();
	mat->J = A->innerIndexPtr();
	Matrix::cudaMemcpyCublasMatrix(mat);
	return mat;
}