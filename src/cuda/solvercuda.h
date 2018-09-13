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
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm, double res, bool consolidatedKernels = false);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm, bool consolidatedKernels = false);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, Vector *x0, numType _LINFinityNorm, bool consolidatedKernels = false);
		void do_iteration();
		//std::vector<numType> getX();
		Eigen::Matrix<numType, -1, 1, 0> getX();
		int getSize() { return this->size; }
		int getIteration();
		Vector *getCpjdsX();

		static numType createPreconditioner(MatrixCPJDS &M, std::unique_ptr<numType[]> &pdata);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness);
		static cgl::Vector *createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int size, int n);
		
		double getResidueSquaredNorm() const;
		double getErrorl2Estimate() const;

		// Debug
		static Eigen::SparseMatrix<double, 0, int> getCpjdsStiffness(MatrixCPJDS &M, std::unique_ptr<numType[]> &pdata);
		static Eigen::Matrix<double, -1, 1, 0> getCpjdsCurrent(numType *vec, MatrixCPJDSManager &mgr, int size, int n);

	private:
		static void cblas_dscal(int n, numType alpha, numType *x, int inc);
		static numType m_preconditioner_eigen(MatrixCPJDS &M, std::unique_ptr<numType[]> &pdata, std::unique_ptr<numType[]> &precond);
		std::vector<numType> transfX2Cpu();
		int size;
		PCGSolverCPJDS *solver;
		MatrixCPJDSManager *mgr;
		numType LINFinityNorm;
		
};

#endif // SOLVERCUDA_H_
