#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "parameters/parametersparser.h"
#include "matrixview.h"
#include "observations.h"
#include <iostream>

int main(int argc, char *argv[])
{
	bool is2dProblem;
	std::shared_ptr<problem> input = problem::createNewProblem(argv[1], &is2dProblem);
	//input->setGroundNode(params.ground);
	input->initProblem(argv[1]);
	observations<double> *readings = new observations<double>;
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	const char* curfile = argv[2];
	int baseNode = argc > 3 ? std::atoi(argv[3]) : -1;
	bool clockwisePositive = argc > 4 ? std::atoi(argv[4]) : true;
	readings->initObs(&curfile, NULL, input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode(), baseNode, input->getGroundNode() - input->getGenericElectrodesCount() + 1, input->getGroundNode(), clockwisePositive);

	matrix *m1;
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[i] = 0.3815;
	input->assembleProblemMatrix(&v[0], &m1);
	input->postAssembleProblemMatrix(&m1);
	Eigen::VectorXd currents;
	Eigen::VectorXd x;
	SparseIncompleteLLT precond(*m1);

	std::vector<Eigen::VectorXd> solutions;
	for (int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
		currents = input->getCurrentVector(patterno, readings);
		CG_Solver solver(*m1, currents, precond);
		for (int i = 0; i < 100; i++) solver.do_iteration();
		x = solver.getX();

#ifdef ZEROELECSUM
		// Correct potentials
		double avg = 0;
		for (int i = input->getNodesCount() - input->getGenericElectrodesCount(); i < input->getNodesCount(); i++) avg += x[i];
		avg /= input->getGenericElectrodesCount();
		for (int i = 0; i < input->getNodesCount(); i++) x[i] -= avg;
#endif
		solutions.push_back(x);
		std::cout << "Finished solution " << patterno + 1 << " of " << readings->getCurrentsCount() << std::endl;
	}

	solution::savePotentials(solutions, "directsol.msh", input, readings);
}
