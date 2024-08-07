#include "pyeitsolver.h"
#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"
#include <memory>

std::shared_ptr<problem> input;
std::unique_ptr<observations<double>> readings;

std::map<std::string, int> init(const char* meshfilename, const  char* currentfilename) {
    // Load mesh file geometry into memory
	bool is2dProblem;
	input = problem::createNewProblem(meshfilename, &is2dProblem);
    input->setIgnoreouterring(true);
	input->initProblem(meshfilename);

    // Create FEM sparse matrix and auxiliary structures
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();

    // Read current pattern
    readings = std::make_unique<observations<double>>();
    //readings = std::make_unique<observations<std::complex<double>>(); //logica para saber qual é

	readings->initObs(&currentfilename, NULL, input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode());

    //TODO tirar argumentos que sobram >> ultimos 2

    // Return problem data
    return std::map<std::string, int> {
        {"coeffCount", input->getNumCoefficients()},
        {"currentsCount", readings->getCurrentsCount()},
        {"electrodesCount", input->getGenericElectrodesCount()},
        {"groundNode", input->getGroundNode()}
    };
}

std::pair<int, std::vector<double>> solveForwardProblem(std::vector<double> conds) {
    // Check conductivity vector size
    size_t n = conds.size();
	if (n != input->getNumCoefficients()) throw std::exception(("Wrong conductivities vector size " + std::to_string(n) + " (should be " + std::to_string(input->getNumCoefficients()) + ")").c_str());

    // Map node conductivities to coefficient indices
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[input->getNode2Coefficient(i)] = conds[i];

    // Create FEM conductivity matrix
	matrix* m1;
	input->assembleProblemMatrix(&v[0], &m1);
	input->postAssembleProblemMatrix(&m1);

    // Create preconditioner matrix
	std::shared_ptr<SparseIncompleteLLT> precond = std::shared_ptr<SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));

    // Solve forward problem to obtain potentials
    std::vector<double> potentials(input->getGenericElectrodesCount() * input->getGenericElectrodesCount());
    Eigen::VectorXd x, currents;
    int noIterations = 0;
    for(int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
        currents = input->getCurrentVector(patterno, readings.get());
        CG_Solver solver(*m1, currents, *precond);
        int i = 0; 
        for (; i < 1500 && solver.getResidueSquaredNorm() > 1e-19; i++) 
            solver.do_iteration();
        noIterations += i;
        //std::cout << "Pattern number: " << patterno << ". Total number of iterations: " << i << std::endl; 
        
        x = solver.getX();

        // Save results to appropriate index in the output vector
        int firstElectrodeIdx = input->getGroundNode() - input->getGenericElectrodesCount() + 1;
        for (int i = 0; i < input->getGenericElectrodesCount(); i++) {
            if (firstElectrodeIdx + i == input->getGroundNode()) potentials[patterno*input->getGenericElectrodesCount() + i] = 0;
            else potentials[patterno*input->getGenericElectrodesCount() +i] = x[firstElectrodeIdx + i] * readings->getCurrentVal(patterno);
        }
    }
    noIterations /= 32;

    // Return potentials
    delete m1;
    return std::pair<int, std::vector<double>>(noIterations, potentials);
}

std::vector<double> solveFullForwardProblem(std::vector<double> conds) {
    // Check conductivity vector size
    size_t n = conds.size();
	if (n != input->getNumCoefficients()) throw std::exception(("Wrong conductivities vector size " + std::to_string(n) + " (should be " + std::to_string(input->getNumCoefficients()) + ")").c_str());

    // Map node conductivities to coefficient indices
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[input->getNode2Coefficient(i)] = conds[i];

    // Create FEM conductivity matrix
	matrix* m1;
	input->assembleProblemMatrix(&v[0], &m1);
	input->postAssembleProblemMatrix(&m1);

    // Create preconditioner matrix
	std::shared_ptr<SparseIncompleteLLT> precond = std::shared_ptr<SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));

    // Solve forward problem to obtain potentials
    std::vector<double> fullPotentials(input->getGenericElectrodesCount() * input->getNumCoefficients());
    Eigen::VectorXd x, currents;
    for(int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
        currents = input->getCurrentVector(patterno, readings.get());
        CG_Solver solver(*m1, currents, *precond);
        for (int i = 0; i < 100; i++) solver.do_iteration();
        x = solver.getX();

        // Save results to appropriate index in the output vector
        for (int i = 0; i < input->getNumCoefficients(); i++) {
            if (i < input->getGroundNode()) fullPotentials[patterno*input->getNumCoefficients() +i] = x[i] * readings->getCurrentVal(patterno);  
            else if (i == input->getGroundNode()) fullPotentials[patterno*input->getNumCoefficients() +i] = 0;  
            else fullPotentials[patterno*input->getNumCoefficients() +i] = x[i-1] * readings->getCurrentVal(patterno);
        }
    }
    
    // Return potentials
    delete m1;
    return fullPotentials;
}