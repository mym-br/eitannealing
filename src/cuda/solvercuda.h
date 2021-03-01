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
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> class Matrix;
}
enum CGSOLVERTYPE {DEFAULT, CONSOLIDATED, CONSOLIDATEDCG};

class CGCUDA_Solver {
	public:
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm, double res, CGSOLVERTYPE solverType = DEFAULT, bool init=true);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm, CGSOLVERTYPE = DEFAULT, bool init=true);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, Vector *x0, numType _LINFinityNorm, CGSOLVERTYPE = DEFAULT, bool init=true);
		void do_iteration();
		//std::vector<numType> getX();
		Eigen::Matrix<numType, -1, 1, 0,-1,1> getX();
		int getSize() { return this->size; }
		int getIteration();
		std::tuple<double, double, double> getIterationTimes();
		Vector *getCpjdsX();

		static numType createPreconditioner(MatrixCPJDS &M, std::unique_ptr<numType[]> &pdata);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<numType,0,int> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<numType, 0, int> *A, MatrixCPJDS *stiffness);
		static cgl::Vector *createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int size, int n);
		
		double getResidueSquaredNorm() const;
		double getErrorl2Estimate() const;

		// Debug
		static Eigen::SparseMatrix<numType, 0, int> getCpjdsStiffness(MatrixCPJDS &M, std::unique_ptr<numType[]> &pdata);
		static Eigen::Matrix<numType, -1, 1, 0, -1, 1> getCpjdsCurrent(numType *vec, MatrixCPJDSManager &mgr, int size, int n);
		void init(double res = -1);

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
