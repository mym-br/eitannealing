/*
 * partial_incomplete_cholesky.cpp
 *
 *  Created on: Aug 20, 2010
 *      Author: thiago
 */

#include "partial_incomplete_cholesky.h"

extern "C" {
#include <cblas.h>
}

#include "solver.h"

Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >
	preparePartialMatrix(const matrix &base, int ncols) {

	// Prepare a copy of the original matrix, stripped of the non-diagonal elements on the ncols first columns
	Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >
		result(base.cols(), base.cols());

	result.startFill(base.nonZeros());
	int i;
	for(i=0;i<ncols;i++) {
		result.fill(i,i) = base.coeff(i,i);
	}

	for(;i<base.cols();i++) {
		matrix::InnerIterator it(base, i);
		while(it) {
			result.fill(it.row(), it.col()) = it.value();
			++it;
		}
	}
	result.endFill();
	return result;
}

PartialIncompleteLLT::PartialIncompleteLLT(
		const Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >& matrix,
		int diagonalRows): SparseIncompleteLLT(preparePartialMatrix(matrix, diagonalRows))
{
}
