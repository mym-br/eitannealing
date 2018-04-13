#include "solvercuda.h"
#include "vector.h"
#include "matrix-cpjds.h"
#include "solver-pcg.h"
#include "../nodecoefficients.h"

using namespace cgl;

CGCUDA_Solver::CGCUDA_Solver(numType *A, numType *b, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients, numType *precond, int n) {
	MatrixCPJDS *stiffness = new MatrixCPJDS;
	mgr = new MatrixCPJDSManager(A, n);
	mgr->buidMatrixCPJDS(stiffness, nodeCoef, nodesCount, numcoefficients);
	size = stiffness->matrixData.n;
	m_preconditioner_eigen(*stiffness, stiffness->cpuData.data, stiffness->cpuData.precond); // FIXME: Use already implemented preconditioner
	cudaMemcpy(stiffness->preconditionedData, stiffness->cpuData.precond, (size_t)stiffness->matrixData.elCount * sizeof(numType), cudaMemcpyHostToDevice);
	Vector *bVec = createCurrentVector(b, *mgr, n);
	solver = new PCGSolverCPJDS(mgr, stiffness, bVec);
	solver->init();
}

Vector *CGCUDA_Solver::createCurrentVector(numType *vec, MatrixCPJDSManager &mgr, int n) {
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