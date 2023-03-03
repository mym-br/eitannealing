#include "solvercublas.h"

#include <Eigen/Core>
#include "../basematrix.h"

using namespace Cublas;

Matrix *Matrix::createCublasMatrix(Eigen::SparseMatrix<double> *A) {
	Matrix *mat = new Matrix;
	mat->N = A->outerSize();
	mat->nz = A->nonZeros();
	mat->val = A->valuePtr();
	mat->I = A->outerIndexPtr();
	mat->J = A->innerIndexPtr();
	Matrix::cudaMemcpyCublasMatrix(mat);
	return mat;
}