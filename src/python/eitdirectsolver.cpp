#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"

std::shared_ptr<problem> input;
observations<double>* readings;
matrix* m1;
std::shared_ptr <SparseIncompleteLLT> precond;

int init(const char* meshfilename, const  char* currentfilename)
{
	bool is2dProblem;
	input = problem::createNewProblem(meshfilename, &is2dProblem);
	//input->setGroundNode(params.ground);
	input->initProblem(meshfilename);
	readings = new observations<double>;
	readings->initObs(&currentfilename, NULL, input->getNodesCount(), input->getGenericElectrodesCount());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();


	
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[i] = 0.3815;
	input->assembleProblemMatrix(&v[0], &m1);
	input->postAssembleProblemMatrix(&m1);
	precond = std::shared_ptr< SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));

	return 2;
}

int solve() {
	Eigen::VectorXd currents;
	Eigen::VectorXd x;
	

	std::vector<Eigen::VectorXd> solutions;
	for (int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
		currents = input->getCurrentVector(patterno, readings);
		CG_Solver solver(*m1, currents, *precond);
		for (int i = 0; i < 100; i++) solver.do_iteration();
		x = solver.getX();

		solutions.push_back(x);
		std::cout << "Finished solution " << patterno + 1 << " of " << readings->getCurrentsCount() << std::endl;
	}

	solution::savePotentials(solutions, "directsol.msh", input, readings);

	return 3;
}