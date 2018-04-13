#ifndef SOLVERCUDA_H_
#define SOLVERCUDA_H_

#include "settings.h"
#include "matrix-cpjds.h"

class PCGSolverCPJDS;
class MatrixCPJDSManager;
namespace cgl { class Vector; };

class CGCUDA_Solver {
	public:
		CGCUDA_Solver(numType *A, numType *b, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients, numType *precond, int n);
		void doIteration();
		std::vector<numType> getX();
		int getSize() { return this->size; }
	private:
		void cblas_dscal(int n, numType alpha, numType *x, int inc);
		bool m_preconditioner_eigen(MatrixCPJDS M, numType * pdata, numType * precond);
		cgl::Vector *createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int n);
		int size;
		PCGSolverCPJDS *solver;
		MatrixCPJDSManager *mgr;
};

#endif // SOLVERCUDA_H_
