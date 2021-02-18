/*
 * incomplete_cholesky.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: thiago
 */

#include "incomplete_LDLT.h"

// TODO: Diagonal is now redundant!

SparseIncompleteLDLT::SparseIncompleteLDLT(const matrix& matrix):
		m_matrix(matrix), invD(matrix.rows()), lInfNormCalc(false)
{	
	for(int col=0;col<m_matrix.cols();col++) {
		matrix::InnerIterator it(m_matrix, col);
		int *outer = m_matrix.outerIndexPtr();
		Scalar *data = m_matrix.valuePtr();
		// 1st element is the diagonal one
		if(it.value()<0)
			throw std::exception();
		Scalar D = it.value();
		Scalar invDi = 1/D;
		// Multiply the whole column
		Eigen::Map<Eigen::VectorXd>(data+outer[col], outer[col+1]-outer[col]) *= invDi;
		// This is not unlike a sparse vector-vector multiplication
		while(++it) {
			matrix::InnerIterator source(it);
			matrix::InnerIterator target(m_matrix, it.row());
			while(target && source) {
				// Sweep and subtract on coincident rows
				//	While n^2 on the number of non-zeros, this should be relatively quick,
				//	as both target and source have very few
				//		non-zero entries
				if(target.row()==source.row()) {
					target.valueRef() -= source.value()*it.value()*D;
					++target; ++source;
				}
				while(target && target.row()<source.row()) ++target;
				while(source && source.row()<target.row()) ++source;
			}
		}
		invD[col] = invDi;
	}

}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLDLT::solveInPlace(Eigen::VectorXd &b)  const
{
  const matrix::Index size = m_matrix.rows();
  assert(size==b.rows());
  m_matrix.triangularView<Eigen::Lower|Eigen::UnitDiag>().solveInPlace(b);
  // Array in Eigen3 means coefficient-wise multiplication
  b.array() *= invD;
  m_matrix.triangularView<Eigen::Lower|Eigen::UnitDiag>().transpose().solveInPlace(b);

  return true;
}


double SparseIncompleteLDLT::getLINFinityNorm() const
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

/*void SparseIncompleteLDLT::multInPlace(Eigen::VectorXd &b) const
{
	b = this->m_matrix*b;
	b = this->m_matrix*b;
}*/
