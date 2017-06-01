/*
 * incomplete_cholesky.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: thiago
 */

#include "incomplete_choleskycomplex.h"

// TODO: Now EIGEN3 has a built-in sparse incomplete Cholesky

SparseIncompleteLLTComplex::SparseIncompleteLLTComplex(const matrixcomplex& matrix) :
		m_matrix(matrix), lInfNormCalc(false)
{
	m_matrix = Eigen::SimplicialLLT<matrixcomplex, Eigen::Lower, Eigen::NaturalOrdering< matrixcomplex::StorageIndex >>(matrix).matrixL();
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLTComplex::solveInPlace(Eigen::VectorXcd &b)  const
{
  const int size = m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().solveInPlace(b);
  m_matrix.triangularView<Eigen::Lower>().transpose().solveInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLTComplex::halfSolveInPlace(Eigen::VectorXcd &b)  const
{
  const int size = m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().solveInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLTComplex::halfSolveInPlaceT(Eigen::VectorXcd &b)  const
{
  const int size = m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().transpose().solveInPlace(b);

  return true;
}

/*
 * Actually this is the L^(-1) infinity norm!
 * Here's the rationale:
 *  If every non-diagonal element of the stiffness matrix
 *  is negative, then it follows that every off-diagonal
 *  element of L is negative (and its diagonal is positive). 
 *  Then every element of L^(-1) is positive and it's
 *  L infinity norm is the max w_i where L w = ones.
 * The problem is, for some problematic meshes (with obtuse
 *  angles) there are a few non-negative off-diagonal elements.
 * Then there might be some non-negative off-diagonal L elements.
 * In this case, there *MIGHT* (verify it!) be negative coefficients
 *  of L^(-1) and we could underestimate the L infinity norm here.
 */
double SparseIncompleteLLTComplex::getLINFinityNorm() const
{
	if(this->lInfNormCalc) {
		return this->lInfNorm;
	}

	Eigen::VectorXcd w(Eigen::VectorXd::Ones(this->m_matrix.rows()));

	this->m_matrix.triangularView<Eigen::Lower>().solveInPlace(w);

	double max = 0;
	for(int i = 0; i<w.rows();i++)
		if(max < abs(w[i]))
			max = abs(w[i]);
	this->lInfNorm = max;
	this->lInfNormCalc = true;
	return max;
}

void SparseIncompleteLLTComplex::multInPlace(Eigen::VectorXcd &b) const
{
	b = this->m_matrix*b;
	b = this->m_matrix*b;
}

void SparseIncompleteLLTComplex::halfMultInPlace(Eigen::VectorXcd &b) const
{
	b = this->m_matrix.transpose()*b;
}
