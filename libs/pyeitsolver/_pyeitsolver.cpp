#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

std::shared_ptr<problem> input;
std::unique_ptr<observations<double>> readings;

py::dict init(const char* meshfilename, const  char* currentfilename) {
    // Load mesh file geometry into memory
	bool is2dProblem;
	input = problem::createNewProblem(meshfilename, &is2dProblem);
	input->initProblem(meshfilename, true);

    // Create FEM sparse matrix and auxiliary structures
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();

    // Read current pattern
    readings = std::make_unique<observations<double>>();
	readings->initObs(&currentfilename, NULL, input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode(), input->getGroundNode() - input->getGenericElectrodesCount() + 1);

    // Return problem data
    return py::dict(
        "coeffCount"_a=input->getNumCoefficients(),
        "currentsCount"_a=readings->getCurrentsCount(),
        "electrodesCount"_a=input->getGenericElectrodesCount(),
        "groundNode"_a=input->getGroundNode()
        );
}

std::vector<double> solveForwardProblem(std::vector<double> conds) {
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
    for(int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
        currents = input->getCurrentVector(patterno, readings.get());
        CG_Solver solver(*m1, currents, *precond);
        for (int i = 0; i < 100; i++) solver.do_iteration();
        x = solver.getX();

        // Save results to appropriate index in the output vector
        int firstElectrodeIdx = input->getGroundNode() - input->getGenericElectrodesCount() + 1;
        for (int i = 0; i < input->getGenericElectrodesCount(); i++) {
            if (firstElectrodeIdx + i == input->getGroundNode()) potentials[patterno*input->getGenericElectrodesCount() + i] = 0;
            else potentials[patterno*input->getGenericElectrodesCount() +i] = x[firstElectrodeIdx + i] * readings->getCurrentVal(patterno);
        }
    }

    // Return potentials
    delete m1;
    return potentials;
}

PYBIND11_MODULE(_pyeitsolver, m) {
    m.doc() = "forward problem solver functions for eit problem"; // optional module docstring

    m.def("init", &init, "Initializes the solver with current and mesh data");
    m.def("solveForwardProblem", &solveForwardProblem, "Solves forward problem with provided conductivities");
}