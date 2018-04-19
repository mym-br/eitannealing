#ifndef SOLVERCUDA_H_
#define SOLVERCUDA_H_

#include "settings.h"
#include "matrix-cpjds.h"

class PCGSolverCPJDS;
class MatrixCPJDSManager;
struct MatrixCPJDS;
namespace cgl { class Vector; };

class CGCUDA_Solver {
	public:
		CGCUDA_Solver(numType *A, numType *b, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients, numType *precond, int n);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec);
		void doIteration();
		std::vector<numType> getX();
		int getSize() { return this->size; }

		static bool m_preconditioner_eigen(MatrixCPJDS M, numType * pdata, numType * precond);
		static MatrixCPJDSManager *createManager(numType * A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients, int n);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
		static cgl::Vector *createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int size, int n);
	private:
		static void cblas_dscal(int n, numType alpha, numType *x, int inc);
		int size;
		PCGSolverCPJDS *solver;
		MatrixCPJDSManager *mgr;
};

#endif // SOLVERCUDA_H_
