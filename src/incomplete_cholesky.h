/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_CHOLESKY_HPP_
#define INCOMPLETE_CHOLESKY_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "basematrix.h"

class SparseIncompleteLLT
{
  protected:

    typedef matrix CholMatrixType;
    CholMatrixType m_matrix;

    mutable bool lInfNormCalc;
    mutable double lInfNorm;

  public:

    SparseIncompleteLLT(const CholMatrixType& matrix);

    /** \returns the lower triangular matrix L */
    inline const CholMatrixType& matrixL(void) const { return m_matrix; }

    bool solveInPlace(Eigen::VectorXd &b) const;
    bool solveInPlaceT(Eigen::VectorXd &b) const {
      return solveInPlace(b);	// symmetric preconditioner!
    }
    bool halfSolveInPlace(Eigen::VectorXd &b) const;
    bool halfSolveInPlaceT(Eigen::VectorXd &b) const;

    void halfMultInPlace(Eigen::VectorXd &b) const;
    void multInPlace(Eigen::VectorXd &b) const;

    double getLINFinityNorm() const;

};



#endif /* INCOMPLETE_CHOLESKY_HPP_ */
