#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"

std::shared_ptr<problem> input;
observations<double>* readings;
matrix* m1;
std::shared_ptr<SparseIncompleteLLT> precond;

void init(const char* meshfilename, const  char* currentfilename)
{
	bool is2dProblem;
	input = problem::createNewProblem(meshfilename, &is2dProblem);
	input->initProblem(meshfilename, true);
	readings = new observations<double>;
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	readings->initObs(&currentfilename, NULL, input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode(), input->getGroundNode() - input->getGenericElectrodesCount() + 1);
}

void setconds(double* cond, int n) {
	if (n != input->getNumCoefficients()) { std::cout << "Wrong conductivities vector size " << n << " (should be " << input->getNumCoefficients() << ")" << std::endl; return;}
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[input->getNode2Coefficient(i)] = cond[i];
	input->assembleProblemMatrix(&v[0], &m1);
	input->postAssembleProblemMatrix(&m1);
	precond = std::shared_ptr< SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));
}

void solve(double* potentials, int n, int patterno) {
	if (n != input->getGenericElectrodesCount()) { std::cout << "Wrong potentials vector size " << n << " (should be " << input->getGenericElectrodesCount() << ")" << std::endl; return; }

	Eigen::VectorXd x, currents;
	currents = input->getCurrentVector(patterno, readings);
	CG_Solver solver(*m1, currents, *precond);
	for (int i = 0; i < 100; i++) solver.do_iteration();
	x = solver.getX();

	int firstElectrodeIdx = input->getGroundNode() - input->getGenericElectrodesCount() + 1;
	for (int i = 0; i < input->getGenericElectrodesCount(); i++) {
		if(firstElectrodeIdx + i == input->getGroundNode()) potentials[i] = 0;
		else potentials[i] = x[firstElectrodeIdx + i] * readings->getCurrentVal(patterno);
	}
}