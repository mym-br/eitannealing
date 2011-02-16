/*
 * partial_incomplete_choleski.h
 *
 *  Created on: Aug 20, 2010
 *      Author: thiago
 */

#ifndef PARTIAL_INCOMPLETE_CHOLESKY_H_
#define PARTIAL_INCOMPLETE_CHOLESKY_H_

#include "incomplete_cholesky.h"

class PartialIncompleteLLT : public SparseIncompleteLLT
{
  public:
	PartialIncompleteLLT(const Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >& matrix, int diagonalRows=0);
};

#endif /* PARTIAL_INCOMPLETE_CHOLESKY_H_ */
