/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_CHOLESKYCOMPLEX_HPP_
#define INCOMPLETE_CHOLESKYXOMPLEX_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "basematrix.h"

// FIXME: Switch to LDmLT!


class SparseIncompleteLLTComplex
{
  protected:

    typedef matrixcomplex CholMatrixTypeComplex;
	CholMatrixTypeComplex m_matrix;

    mutable bool lInfNormCalc;
    mutable double lInfNorm;

  public:

	  SparseIncompleteLLTComplex(const CholMatrixTypeComplex& matrix);

    /** \returns the lower triangular matrix L */
	  inline const CholMatrixTypeComplex& matrixL(void) const { return m_matrix; }

    bool solveInPlace(Eigen::VectorXcd &b) const;
    bool solveInPlaceT(Eigen::VectorXcd &b) const {
      return solveInPlace(b);	// symmetric preconditioner!
    }
    bool halfSolveInPlace(Eigen::VectorXcd &b) const;
    bool halfSolveInPlaceT(Eigen::VectorXcd &b) const;

    void halfMultInPlace(Eigen::VectorXcd &b) const;
    void multInPlace(Eigen::VectorXcd &b) const;

    double getLINFinityNorm() const;

};



#endif /* INCOMPLETE_CHOLESKYXOMPLEX_HPP_ */
