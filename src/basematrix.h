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
typedef  std::complex<double> Complex;
// Only the lower triangular coefficients are present!
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> matrix;
typedef Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> matrixcomplex;

// Dense vector
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vectorx;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> vectorxcomplex;

#endif	// BASEMATRIX_H_