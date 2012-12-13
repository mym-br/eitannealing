#ifndef SPARSEINCOMPLETELQ_H
#define SPARSEINCOMPLETELQ_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class SparseIncompleteLQ
{
  public:
    typedef double Scalar;

    typedef Eigen::SparseMatrix<Scalar, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor > BaseMatrixType;
  
  protected:
    typedef Eigen::SparseMatrix<Scalar, Eigen::LowerTriangular | Eigen::ColMajor> LMatrixType;
  
    LMatrixType m_matrix;

  public:

    /** \returns the lower triangular matrix L */
    inline const LMatrixType& matrixL(void) const { return m_matrix; }

    bool solveInPlace(Eigen::VectorXd &b) const;
    bool halfSolveInPlace(Eigen::VectorXd &b) const;
    bool halfSolveInPlaceT(Eigen::VectorXd &b) const;

    void halfMultInPlace(Eigen::VectorXd &b) const;
    void multInPlace(Eigen::VectorXd &b) const;

  
    SparseIncompleteLQ(const BaseMatrixType& matrix, unsigned int iq, unsigned int il);
};

#endif // SPARSEINCOMPLETELQ_H
