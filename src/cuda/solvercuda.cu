#include "solvercuda.h"
#include "vector.h"
#include "matrix-cpjds.h"
#include "solver-pcg.h"
#include "../nodecoefficients.h"

using namespace cgl;

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, double _LINFinityNorm, double res, CGSOLVERTYPE solverType, bool init) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	switch (solverType) {
	case DEFAULT: solver = new PCGSolverCPJDS(mgr, stiffness, bVec); break;
	case CONSOLIDATED: solver = new PCGSolverConsolidatedCPJDS(mgr, stiffness, bVec); break;
	#ifdef CGROUPS
	case CONSOLIDATEDCG: solver = new PCGSolverConsolidatedCPJDSCG(mgr, stiffness, bVec); break;
	#endif
	}
	if(init) solver->init(res);
}

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, double _LINFinityNorm, CGSOLVERTYPE solverType, bool init) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	switch (solverType) {
	case DEFAULT: solver = new PCGSolverCPJDS(mgr, stiffness, bVec); break;
	case CONSOLIDATED: solver = new PCGSolverConsolidatedCPJDS(mgr, stiffness, bVec); break;
	#ifdef CGROUPS
	case CONSOLIDATEDCG: solver = new PCGSolverConsolidatedCPJDSCG(mgr, stiffness, bVec); break;
	#endif
	}
	if(init) solver->init();
}

CGCUDA_Solver::CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, Vector *x0, double _LINFinityNorm, CGSOLVERTYPE solverType, bool init) : mgr(mgr), LINFinityNorm(_LINFinityNorm) {
	size = stiffness->matrixData.n;
	switch (solverType) {
	case DEFAULT: solver = new PCGSolverCPJDS(mgr, stiffness, bVec); break;
	case CONSOLIDATED: solver = new PCGSolverConsolidatedCPJDS(mgr, stiffness, bVec); break;
	#ifdef CGROUPS
	case CONSOLIDATEDCG: solver = new PCGSolverConsolidatedCPJDSCG(mgr, stiffness, bVec); break;
	#endif
	}
	if(init) solver->init(x0);
}

void CGCUDA_Solver::init(double res) {
	solver->init(res);
}

Vector *CGCUDA_Solver::createCurrentVector(double *vec, MatrixCPJDSManager &mgr, int size, int n) {
	double * vecArr = new double[size];
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

std::tuple<double, double, double> CGCUDA_Solver::getIterationTimes() {
	return solver->getAvgTimes();
}

std::vector<double> CGCUDA_Solver::transfX2Cpu() {
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

double CGCUDA_Solver::createPreconditioner(MatrixCPJDS &M, std::unique_ptr<double[]> &pdata) {
	double ans = m_preconditioner_eigen(M, M.cpuData.data, M.cpuData.precond); // FIXME: Use already implemented preconditioner
	cudaMemcpy(M.preconditionedData.get(), M.cpuData.precond.get(), (size_t)M.matrixData.elCount * sizeof(double), cudaMemcpyHostToDevice);
	return ans;
}

Vector *CGCUDA_Solver::getCpjdsX() {
	return solver->getX();
}