/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_QR_COMPLEX_HPP_
#define INCOMPLETE_QR_COMPLEX_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>

#include "basematrix.h"
#include "incomplete_qr_builder.h"

template<class scalar> class SparseIncompleteQR
{
  protected:
    std::vector<double> idiagonal;
    std::vector<std::vector<double> > rowsR;
    std::vector<std::vector<double> > rowsI;
    std::vector<std::vector<double> > colsR;
    std::vector<std::vector<double> > colsI;

    typedef SparseIncompleteQRBuilder<scalar> QRMatrixType;
    QRMatrixType m_matrix;


  public:

    SparseIncompleteQR(const QRMatrixType& matrix);

    bool solveInPlaceC(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const;
    bool solveInPlaceCT(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const;
};



#endif /* INCOMPLETE_QR_COMPLEX_HPP_ */
