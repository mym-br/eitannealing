/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_LDLT_HPP_
#define INCOMPLETE_LDLT_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "basematrix.h"


class SparseIncompleteLDLT
{
  protected:

    typedef matrix CholMatrixType;
    CholMatrixType m_matrix;
    Eigen::DiagonalMatrix<Scalar> Di;

    mutable bool lInfNormCalc;
    mutable double lInfNorm;

  public:

    SparseIncompleteLDLT(const CholMatrixType& matrix);

    /** \returns the lower triangular matrix L */
    inline const CholMatrixType& matrixL(void) const { return m_matrix; }

    bool solveInPlace(Eigen::VectorXd &b) const;
    void multInPlace(Eigen::VectorXd &b) const;

    double getLINFinityNorm() const;

};



#endif /* INCOMPLETE_LDLT_HPP_ */
