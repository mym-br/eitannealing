/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_QR_COMPLEX2_HPP_
#define INCOMPLETE_QR_COMPLEX2_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <unordered_map>
#include <iostream>

#include "basematrix.h"
#include "incomplete_qr_builder.h"
#include "sparse_incomplete_qr.h"

typedef SparseIncompleteQR<std::complex<double> > SparseIncompleteQRComplex2;

#endif /* INCOMPLETE_QR_COMPLEX2_HPP_ */
