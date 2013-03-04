#ifndef SPARSEINCOMPLETELQ_H
#define SPARSEINCOMPLETELQ_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class SparseIncompleteLQ
{
  public:
    typedef double Scalar;

    typedef Eigen::SparseMatrix<Scalar, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor > BaseMatrixType;
  // FIXME: This is actually a QR decomposition!
  protected:
    typedef Eigen::SparseMatrix<Scalar, Eigen::UpperTriangular | Eigen::ColMajor> LMatrixType;
  
    LMatrixType m_matrix;

  public:

    /** \returns the lower triangular matrix L */
    inline const LMatrixType& matrixL(void) const { return m_matrix; }

    bool solveInPlace(Eigen::VectorXd &b) const;
    bool solveInPlaceT(Eigen::VectorXd &b) const;
    
    void multInPlace(Eigen::VectorXd &b) const;

  
    SparseIncompleteLQ(const BaseMatrixType& matrix, unsigned int iq, unsigned int il);
};

#endif // SPARSEINCOMPLETELQ_H
