/*
 * incomplete_cholesky.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: thiago
 */

#include "incomplete_cholesky.h"

// TODO: Now EIGEN3 has a built-in sparse incomplete Cholesky

SparseIncompleteLLT::SparseIncompleteLLT(const matrix& matrix):
		m_matrix(matrix), lInfNormCalc(false)
{
	int col=0;
	// just the diagonal for Electrode nodes!
	/*for(;col<31;col++) {
		Eigen::SparseMatrix<double, Eigen::LowerTriangular>::InnerIterator it(m_matrix, col);
		// just the diagonal, thanks!
		it.valueRef() = sqrt(it.value());
		while(it) {
			++it;
			it.valueRef() = 0;
		}
	}*/
	for(;col<m_matrix.cols();col++) {
		matrix::InnerIterator it(m_matrix, col);
		int *outer = m_matrix.outerIndexPtr();
		Scalar *data = m_matrix.valuePtr();
		// 1st element is the diagonal one
		if(it.value()<0)
			throw std::exception();
		double isqrtDiagonal = 1/std::sqrt(it.value());
		// Multiply the whole column
		Eigen::Map<Eigen::VectorXd>(data+outer[col], outer[col+1]-outer[col]) *= isqrtDiagonal;
		// This is not unlike a sparse vector-vector multiplication
		while(++it) {
			matrix::InnerIterator source(it);
			matrix::InnerIterator target(m_matrix, it.row());
			while(target && source) {
				// Sweep and subtract on coincident rows
				//	This should be relatively quick, as both target and source have very few
				//		non-zero entries
				if(target.row()==source.row()) {
					target.valueRef() -= source.value()*it.value();
					++target; ++source;
				}
				while(target && target.row()<source.row()) ++target;
				while(source && source.row()<target.row()) ++source;
			}
		}
	}

}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLT::solveInPlace(Eigen::VectorXd &b)  const
{
  const matrix::Index size = m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().solveInPlace(b);
  m_matrix.triangularView<Eigen::Lower>().transpose().solveInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLT::halfSolveInPlace(Eigen::VectorXd &b)  const
{
  const matrix::Index size = m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().solveInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLT::halfSolveInPlaceT(Eigen::VectorXd &b)  const
{
  const matrix::Index size = m_matrix.rows();
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
double SparseIncompleteLLT::getLINFinityNorm() const
{
	if(this->lInfNormCalc) {
		return this->lInfNorm;
	}

	Eigen::VectorXd w(Eigen::VectorXd::Ones(this->m_matrix.rows()));

	this->m_matrix.triangularView<Eigen::Lower>().solveInPlace(w);

	double max = 0;
	for(int i = 0; i<w.rows();i++)
		if(max < w[i])
			max = w[i];
	this->lInfNorm = max;
	this->lInfNormCalc = true;
	return max;
}

void SparseIncompleteLLT::multInPlace(Eigen::VectorXd &b) const
{
	b = this->m_matrix*b;
	b = this->m_matrix*b;
}

void SparseIncompleteLLT::halfMultInPlace(Eigen::VectorXd &b) const
{
	b = this->m_matrix.transpose()*b;
}
