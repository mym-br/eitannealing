#include "directsolver.h"
#include "solver.h"
#include "problem.h"
#include "solution.h"

EitDirectSolver::EitDirectSolver(const char* meshfilename, const  char* currentfilename) {
	bool is2dProblem;
	input = problem::createNewProblem(meshfilename, &is2dProblem);
	input->initProblem(meshfilename, true);
	readings = std::make_unique<observations<double>>();
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	readings->initObs(&currentfilename, NULL, input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode(), input->getGroundNode() - input->getGenericElectrodesCount() + 1);
	m1 = NULL;
}

int EitDirectSolver::getCoeffCount() {
	return input->getNumCoefficients();
}

int EitDirectSolver::getCurrentPatternCount() {
	return readings->getCurrentsCount();
}

int EitDirectSolver::getElectrodesCount() {
	return input->getGenericElectrodesCount();
}
void EitDirectSolver::setconds(double* cond, int n) {
	if (n != input->getNumCoefficients()) { std::cout << "Wrong conductivities vector size " << n << " (should be " << input->getNumCoefficients() << ")" << std::endl; return; }
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[input->getNode2Coefficient(i)] = cond[i];
	if(m1) delete m1;
	input->assembleProblemMatrix(&v[0], &m1);
	input->postAssembleProblemMatrix(&m1);
	precond = std::shared_ptr< SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));
}

double* EitDirectSolver::solve(int patterno) {
	Eigen::VectorXd x, currents;
	currents = input->getCurrentVector(patterno, readings.get());
	CG_Solver solver(*m1, currents, *precond);
	for (int i = 0; i < 100; i++) solver.do_iteration();
	x = solver.getX();

	double *potentials = new double[input->getGenericElectrodesCount()];
	int firstElectrodeIdx = input->getGroundNode() - input->getGenericElectrodesCount() + 1;
	for (int i = 0; i < input->getGenericElectrodesCount(); i++) {
		if (firstElectrodeIdx + i == input->getGroundNode()) potentials[i] = 0;
		else potentials[i] = x[firstElectrodeIdx + i] * readings->getCurrentVal(patterno);
	}

	return potentials;
}