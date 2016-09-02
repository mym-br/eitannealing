/*
 * basematrix.h
 *
 *  Created on: Aug 29, 2016
 *      Author: thiago
 */

#ifndef BASEMATRIX_H_
#define BASEMATRIX_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>

typedef  double Scalar;
// Only the lower triangular coefficients are present!
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> matrix;

#endif	// BASEMATRIX_H_