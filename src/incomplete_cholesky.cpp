/*
 * incomplete_cholesky.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: thiago
 */

#include "incomplete_cholesky.h"

#ifdef USE_CBLAS
extern "C" {
#include <cblas.h>
}
#else

void cblas_dscal(int n, double alpha, double *x, int inc)
{
	int i;
	for(i=0;i<n;i++) {
		*x *= alpha;
		x += inc;
	}
}
#endif


SparseIncompleteLLT::SparseIncompleteLLT(const Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >& matrix):
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
		Eigen::SparseMatrix<double, Eigen::LowerTriangular>::InnerIterator it(m_matrix, col);
		int *outer = m_matrix._outerIndexPtr();
		double *data = m_matrix._valuePtr();
		// 1st element is the diagonal one
		/*if(it.value()<=0) {
			std::cout << ">>> " << it.value() << "  OMGWTF??????????";
			throw std::exception();
		}*/
		if(it.value()<0)
			std::cout << "uh-oh...";
		double isqrtDiagonal = 1/Eigen::ei_sqrt(it.value());
		// Multiply the whole column
		cblas_dscal(outer[col+1]-outer[col], isqrtDiagonal, data+outer[col], 1);
		// This is not unlike a sparse vector-vector multiplication
		while(++it) {
			Eigen::SparseMatrix<double, Eigen::LowerTriangular>::InnerIterator source(it);
			Eigen::SparseMatrix<double, Eigen::LowerTriangular>::InnerIterator target(m_matrix, it.row());
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
  const int size = m_matrix.rows();
  ei_assert(size==b.rows());
  m_matrix.solveTriangularInPlace(b);
  m_matrix.transpose().solveTriangularInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLT::halfSolveInPlace(Eigen::VectorXd &b)  const
{
  const int size = m_matrix.rows();
  ei_assert(size==b.rows());
  m_matrix.solveTriangularInPlace(b);

  return true;
}

/** Computes b = L^-T L^-1 b */
bool SparseIncompleteLLT::halfSolveInPlaceT(Eigen::VectorXd &b)  const
{
  const int size = m_matrix.rows();
  ei_assert(size==b.rows());
  m_matrix.transpose().solveTriangularInPlace(b);

  return true;
}

double SparseIncompleteLLT::getLINFinityNorm() const
{
	if(this->lInfNormCalc) {
		return this->lInfNorm;
	}

	Eigen::VectorXd w(Eigen::VectorXd::Ones(this->m_matrix.rows()));

	this->m_matrix.solveTriangularInPlace(w);

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
