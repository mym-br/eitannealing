#include "solvercuda.h"
#include "vector.h"
#include "matrix-cpjds.h"
#include "solver-pcg.h"
#include "../nodecoefficients.h"

using namespace cgl;

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm, double res, bool consolidatedKernels) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	if (consolidatedKernels) solver = new PCGSolverCPJDS2(mgr, stiffness, bVec);
	else solver = new PCGSolverCPJDS(mgr, stiffness, bVec);
	solver->init(res);
}

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, numType _LINFinityNorm, bool consolidatedKernels) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	if(consolidatedKernels) solver = new PCGSolverCPJDS2(mgr, stiffness, bVec);
	else solver = new PCGSolverCPJDS(mgr, stiffness, bVec);
	solver->init();
}

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, Vector *x0, numType _LINFinityNorm, bool consolidatedKernels) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	if (consolidatedKernels) solver = new PCGSolverCPJDS2(mgr, stiffness, bVec);
	else solver = new PCGSolverCPJDS(mgr, stiffness, bVec);
	solver->init(x0);
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

void CGCUDA_Solver::do_iteration() {
	solver->doIteration();
	cudaDeviceSynchronize();
}

int CGCUDA_Solver::getIteration() {
	return solver->getIteration();
}

//std::vector<numType> CGCUDA_Solver::getX() {
//	Vector *x = solver->getX();
//	
//	return mgr->restore(x);
//}
//
//Eigen::VectorXd CGCUDA_Solver::getX() {
//	//	Vector *x = solver->getX();
//	//	
//	//	return mgr->restore(x);
//	//}
//}

std::vector<numType> CGCUDA_Solver::transfX2Cpu() {
	Vector *x = solver->getX();
		
	return mgr->restore(x);
}

MatrixCPJDSManager *CGCUDA_Solver::createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients) {
	MatrixCPJDSManager *mgr = new MatrixCPJDSManager(A);
	mgr->buidMatrixCPJDS(stiffness, nodeCoef, nodesCount, numcoefficients);
	return mgr;
}

MatrixCPJDSManager *CGCUDA_Solver::createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness) {
	MatrixCPJDSManager *mgr = new MatrixCPJDSManager(A);
	mgr->buidMatrixCPJDS(stiffness);
	return mgr;
}
double CGCUDA_Solver::getResidueSquaredNorm() const {
	return solver->getRmod();
}

double CGCUDA_Solver::getErrorl2Estimate() const {
	double r0norm = solver->getR0norm();
	return r0norm * r0norm * this->solver->getCurrentErr() * LINFinityNorm;
}

numType CGCUDA_Solver::createPreconditioner(MatrixCPJDS &M, std::unique_ptr<numType[]> &pdata) {
	numType ans = m_preconditioner_eigen(M, M.cpuData.data, M.cpuData.precond); // FIXME: Use already implemented preconditioner
	cudaMemcpy(M.preconditionedData.get(), M.cpuData.precond.get(), (size_t)M.matrixData.elCount * sizeof(numType), cudaMemcpyHostToDevice);
	return ans;
}

Vector *CGCUDA_Solver::getCpjdsX() {
	return solver->getX();
}