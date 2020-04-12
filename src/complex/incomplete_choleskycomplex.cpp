/*
 * incomplete_cholesky.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: thiago
 */

#include "incomplete_choleskycomplex.h"
#include <iostream>

// TODO: Now EIGEN3 has a built-in sparse incomplete Cholesky

SparseIncompleteLLTComplex::SparseIncompleteLLTComplex(const matrixcomplex& matrix, bool init) :
		m_matrix(matrix), lInfNormCalc(false)
{
	if(init) m_matrix = Eigen::SimplicialLLT<matrixcomplex, Eigen::Lower, Eigen::NaturalOrdering< matrixcomplex::StorageIndex >>(matrix).matrixL();
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLTComplex::solveInPlace(Eigen::VectorXcd &b)  const
{
  const unsigned int size = (unsigned int)m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().solveInPlace(b);
  m_matrix.triangularView<Eigen::Lower>().transpose().conjugate().solveInPlace(b);

  return true;
}

bool SparseIncompleteLLTComplex::solveInPlace2(Eigen::VectorXcd &b) const
{
	const unsigned int size = (unsigned int)m_matrix.rows();
	assert(size == b.rows());
	CholMatrixTypeComplex m_matrixh = m_matrix.conjugate();
	m_matrixh.triangularView<Eigen::Lower>().solveInPlace(b); // b <- w
	m_matrixh.triangularView<Eigen::Lower>().transpose().conjugate().solveInPlace(b); // b <- z
	m_matrix.triangularView<Eigen::Lower>().solveInPlace(b); // b <- y
	m_matrix.triangularView<Eigen::Lower>().transpose().conjugate().solveInPlace(b); // b <- x

	return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLTComplex::halfSolveInPlace(Eigen::VectorXcd &b)  const
{
  const unsigned int size = (unsigned int)m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower>().solveInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLTComplex::halfSolveInPlaceT(Eigen::VectorXcd &b)  const
{
  const unsigned int size = (unsigned int)m_matrix.rows();
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

SparseIncompleteLLTComplex2::SparseIncompleteLLTComplex2(const CholMatrixTypeComplex& matrix) : SparseIncompleteLLTComplex(matrix, false) {
	int col = 0;
	
	for (; col<m_matrix.cols(); col++) {
	//for (; col<2; col++) {
		matrixcomplex::InnerIterator it(m_matrix, col);
		int *outer = m_matrix.outerIndexPtr();
		Complex *data = m_matrix.valuePtr();
		// 1st element is the diagonal one. FIXME: complex exception condition?
		//if (it.value()<0)
			//throw std::exception();
		Complex isqrtDiagonal = Complex(1,0) / std::sqrt(it.value());
		Complex itValueBackup = it.value();
		// Multiply the whole column
		Eigen::Map<Eigen::VectorXcd>(data + outer[col], outer[col + 1] - outer[col]) *= isqrtDiagonal;
		// This is not unlike a sparse vector-vector multiplication
		while (++it) {
			matrixcomplex::InnerIterator source(it);
			matrixcomplex::InnerIterator target(m_matrix, it.row());
			while (target && source) {
				// Sweep and subtract on coincident rows
				//	This should be relatively quick, as both target and source have very few
				//		non-zero entries
				if (target.row() == source.row()) {
					target.valueRef() -= source.value()*std::conj(it.value());
					//std::cout << "(" << target.row() + 1 << ", " << target.col() + 1 << ") -= (" << source.row() + 1 << ", " << source.col() + 1 << ")*(" << it.row() + 1 << ", " << it.col() + 1 << ") = " << target.value() << std::endl;
					++target; ++source;
				}
				while (target && target.row()<source.row()) ++target;
				while (source && source.row()<target.row()) ++source;
			}
		}
	}
}

bool SparseIncompleteLLTComplex2::solveInPlace2(Eigen::VectorXcd &b) const
{
	const unsigned int size = (unsigned int)m_matrix.rows();
	assert(size == b.rows());
	CholMatrixTypeComplex m_matrixh = m_matrix.conjugate();
	m_matrixh.triangularView<Eigen::Lower>().solveInPlace(b); // b <- w
	m_matrixh.triangularView<Eigen::Lower>().transpose().solveInPlace(b); // b <- z
	m_matrix.triangularView<Eigen::Lower>().solveInPlace(b); // b <- y
	m_matrix.triangularView<Eigen::Lower>().transpose().solveInPlace(b); // b <- x

	return true;
}