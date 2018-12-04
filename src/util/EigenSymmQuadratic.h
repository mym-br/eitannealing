/*
  
  EigenSymmmQuadratic.hp
  quick-and-dirt Symmetric Quadratic Form evaluation for eigen3 

    By Thiago C. Martins
  
*/

#ifndef _EIGENSYMMQUADRATIC_H_
#define _EIGENSYMMQUADRATIC_H_

#include <Eigen/SparseCore>

// Lower triangular Symmetric view
template <typename _Scalar,typename t_vector > _Scalar EigenSymmQuadraticL(
  const typename Eigen::SparseMatrix<_Scalar, Eigen::ColMajor>::template SelfAdjointViewReturnType<Eigen::Lower>::Type &m,
  const Eigen::MatrixBase<t_vector>& x)
{
    typedef typename Eigen::SparseMatrix<_Scalar, Eigen::ColMajor>::Index Index;
    Index col, nextcol;
	_Scalar res = 0;
    col = m.matrix().outerIndexPtr()[0]; 
    for(Index j=0;j<m.matrix().outerSize();j++) {
	nextcol = m.matrix().outerIndexPtr()[j+1];
	if(col==nextcol) continue;
	_Scalar xi = x[m.matrix().innerIndexPtr()[col]];
	res += xi*xi*m.matrix().valuePtr()[col];
	col++;
	while(col!=nextcol) {
	  //res += 2*xi*x[m.matrix().innerIndexPtr()[col]]*m.matrix().valuePtr()[col];
		res += xi*x[m.matrix().innerIndexPtr()[col]] * m.matrix().valuePtr()[col];
		res += xi*x[m.matrix().innerIndexPtr()[col]] * m.matrix().valuePtr()[col];
	  col++;
	}
    }
    return res;
}

#endif	// _EIGENSYMMQUADRATIC_H_