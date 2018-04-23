#include "solvercuda.h"
#include "vector.h"
#include "matrix-cpjds.h"
#include "solver-pcg.h"
#include "../nodecoefficients.h"

using namespace cgl;

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	solver = new PCGSolverCPJDS(mgr, stiffness, bVec);
	solver->init();
}

Vector *CGCUDA_Solver::createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int size, int n) {
	numType * vecArr = new numType[size];
	for (int i = 0; i < size; i++) {
		vecArr[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		vecArr[mgr.original2PaddedIdx[i]] = vec[i];
	}
	return new Vector(vecArr, size);
}

void CGCUDA_Solver::doIteration() {
	solver->doIteration();
	cudaDeviceSynchronize();
}

std::vector<numType> CGCUDA_Solver::getX() {
	Vector *x = solver->getX();
	
	return mgr->restore(x);
}

MatrixCPJDSManager *CGCUDA_Solver::createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients) {
	MatrixCPJDSManager *mgr = new MatrixCPJDSManager(A);
	mgr->buidMatrixCPJDS(stiffness, nodeCoef, nodesCount, numcoefficients);
	return mgr;
}

numType CGCUDA_Solver::getResidueSquaredNorm() const {
	return solver->getRmod();
}

numType CGCUDA_Solver::getErrorl2Estimate() const {
	numType r0norm = solver->getR0norm();
	return r0norm * r0norm * this->solver->getCurrentErr() * LINFinityNorm;
}

numType CGCUDA_Solver::createPreconditioner(MatrixCPJDS M, numType * pdata, numType * precond) {
	numType ans = m_preconditioner_eigen(M, M.cpuData.data, M.cpuData.precond); // FIXME: Use already implemented preconditioner
	cudaMemcpy(M.preconditionedData, M.cpuData.precond, (size_t)M.matrixData.elCount * sizeof(numType), cudaMemcpyHostToDevice);
	return ans;
}