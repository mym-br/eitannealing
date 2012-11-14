#ifndef SPARSEINCOMPLETELQ_H
#define SPARSEINCOMPLETELQ_H

class SparseIncompleteLQ
{
  protected:
    typedef double Scalar;

    typedef Eigen::SparseMatrix<Scalar, Eigen::LowerTriangular> LMatrixType;
    LMatrixType m_matrix;

  public:

    /** \returns the lower triangular matrix L */
    inline const LMatrixType& matrixL(void) const { return m_matrix; }

    bool solveInPlace(Eigen::VectorXd &b) const;
    bool halfSolveInPlace(Eigen::VectorXd &b) const;
    bool halfSolveInPlaceT(Eigen::VectorXd &b) const;

    void halfMultInPlace(Eigen::VectorXd &b) const;
    void multInPlace(Eigen::VectorXd &b) const;

  
    SparseIncompleteLQ(const Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >& matrix, unsigned int iq, unsigned int il);
};

#endif // SPARSEINCOMPLETELQ_H
