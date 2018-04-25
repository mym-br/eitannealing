#ifndef SOLVERCUDA_H_
#define SOLVERCUDA_H_

#include "settings.h"
#include "matrix-cpjds.h"

class PCGSolverCPJDS;
class MatrixCPJDSManager;
struct MatrixCPJDS;
namespace cgl { class Vector; };
namespace Eigen {
	//template<typename _Scalar, int _Flags = 0, typename _StorageIndex = int>  class SparseMatrix;
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows = _Rows, int _MaxCols = _Cols> class Matrix;
}

class CGCUDA_Solver {
	public:
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, Vector *x0, numType _LINFinityNorm);
		void do_iteration();
		//std::vector<numType> getX();
		Eigen::Matrix<numType, -1, 1, 0> getX();
		int getSize() { return this->size; }
		int getIteration();
		Vector *getCpjdsX();

		static numType createPreconditioner(MatrixCPJDS M, numType * pdata, numType * precond);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
		static cgl::Vector *createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int size, int n);
		
		numType getResidueSquaredNorm() const;
		numType getErrorl2Estimate() const;

	private:
		static void cblas_dscal(int n, numType alpha, numType *x, int inc);
		static numType m_preconditioner_eigen(MatrixCPJDS M, numType * pdata, numType * precond);
		std::vector<numType> transfX2Cpu();
		int size;
		PCGSolverCPJDS *solver;
		MatrixCPJDSManager *mgr;
		numType LINFinityNorm;
		
};

#endif // SOLVERCUDA_H_
