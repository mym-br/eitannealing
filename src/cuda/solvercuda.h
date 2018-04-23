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
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm);
		void doIteration();
		std::vector<numType> getX();
		int getSize() { return this->size; }

		static numType createPreconditioner(MatrixCPJDS M, numType * pdata, numType * precond);
		static MatrixCPJDSManager *createManager(numType * A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients, int n);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
		static cgl::Vector *createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int size, int n);
		
		numType getResidueSquaredNorm() const;
		numType getErrorl2Estimate() const;

	private:
		static void cblas_dscal(int n, numType alpha, numType *x, int inc);
		static numType m_preconditioner_eigen(MatrixCPJDS M, numType * pdata, numType * precond);
		int size;
		PCGSolverCPJDS *solver;
		MatrixCPJDSManager *mgr;
		numType LINFinityNorm;
		
};

#endif // SOLVERCUDA_H_
