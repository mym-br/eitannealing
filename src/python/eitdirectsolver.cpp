#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"

int solve()
{
	bool is2dProblem;
	std::shared_ptr<problem> input = problem::createNewProblem("circular_A_2D.msh", &is2dProblem);
	//input->setGroundNode(params.ground);
	input->initProblem("circular_A_2D.msh");
	observations<double> *readings = new observations<double>;
	const char *curfile = "cuba_190ma_cp.txt";
	readings->initObs(&curfile, NULL, input->getNodesCount(), input->getGenericElectrodesCount());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();


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
		
		solutions.push_back(x);
		std::cout << "Finished solution " << patterno + 1 << " of " << readings->getCurrentsCount() << std::endl;
	}

	solution::savePotentials(solutions, "directsol.msh", input, readings);
	
	return 2;
}
