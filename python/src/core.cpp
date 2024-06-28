#include <pybind11/pybind11.h>
#include <pybind11/eigen/matrix.h>
#include <pybind11/stl.h>
#include <Eigen/Core>
#include <memory>
#include <map>
#include <string>
#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"

#define STRINGIFY(x) #x
#define MAX_ITERATIONS 1500
#define RESIDUAL 1e-19

namespace pyeitsolver
{
    class EitSolver {
        public:
        EitSolver(std::string meshFilename, std::string currentFilename) {
            this->info = init(meshFilename, currentFilename);
        }

        // Function to create the indices vector
        Eigen::VectorXi createIndicesVector(std::shared_ptr<problem> input)
        {
            int numCoefficients = input->getNumCoefficients();
            Eigen::VectorXi indices(numCoefficients);
            std::vector<int> temp(numCoefficients);
            std::iota(temp.begin(), temp.end(), 0); // Fill with sequential integers starting from 0

            std::transform(temp.begin(), temp.end(), indices.data(), [&input](int i)
                        { return input->getNode2Coefficient(i); });

            return indices;
        }

        std::pair<int, Eigen::MatrixXd> solve_forward_problem(const Eigen::VectorXd &conds, bool meshPotentials = false)
        {
            // Check conductivity vector size
            auto n = conds.size();
            if (n != input->getNumCoefficients())
                throw std::exception(("Wrong conductivities vector size " + std::to_string(n) + " (should be " + std::to_string(input->getNumCoefficients()) + ")").c_str());

            // Map node conductivities to coefficient indices
            Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
            Eigen::VectorXi indices = createIndicesVector(input);
            v(indices) = conds;

            // Create FEM conductivity matrix
            matrix *m1;
            input->assembleProblemMatrix(&v[0], &m1);
            input->postAssembleProblemMatrix(&m1);

            // Create preconditioner matrix
            std::shared_ptr<SparseIncompleteLLT> precond = std::shared_ptr<SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));

            // Solve forward problem to obtain potentials
            int electrodeCount = input->getGenericElectrodesCount();
            Eigen::MatrixXd potentials(electrodeCount, meshPotentials ? n : electrodeCount);
            Eigen::VectorXd x, currents;
            int noIterations = 0;
            for (int patterno = 0; patterno < readings->getCurrentsCount(); patterno++)
            {
                currents = input->getCurrentVector(patterno, readings.get());
                CG_Solver solver(*m1, currents, *precond);
                int i = 0;
                for (; i < MAX_ITERATIONS && solver.getResidueSquaredNorm() > RESIDUAL; i++)
                    solver.do_iteration();
                noIterations += i;

                x = solver.getX();

                // Save results to appropriate index in the output vector
                int firstElectrodeIdx = input->getGroundNode() - electrodeCount + 1;
                potentials.row(patterno) = meshPotentials ? x : x.segment(firstElectrodeIdx, electrodeCount);
                potentials.row(patterno) *= readings->getCurrentVal(patterno);
            }
            noIterations /= 32;

            // Return potentials
            delete m1;
            return std::pair<int, Eigen::MatrixXd>(noIterations, potentials);
        }

        std::map<std::string, int> info;

        private:
        std::shared_ptr<problem> input;
        std::unique_ptr<observations<double>> readings;

        std::map<std::string, int> init(std::string meshFilename, std::string currentFilename)
        {
            // Load mesh file geometry into memory
            bool is2dProblem;
            input = problem::createNewProblem(meshFilename.c_str(), &is2dProblem);
            input->setIgnoreouterring(true);
            input->initProblem(meshFilename.c_str());

            // Create FEM sparse matrix and auxiliary structures
            input->buildNodeCoefficients();
            input->prepareSkeletonMatrix();
            input->createCoef2KMatrix();

            // Read current pattern
            readings = std::make_unique<observations<double>>();
            const char *currentFilenameCStar = currentFilename.c_str();
            readings->initObs(&currentFilenameCStar, NULL, input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode());

            // Return problem data
            return std::map<std::string, int>{
                {"coeff_count", input->getNumCoefficients()},
                {"currents_count", readings->getCurrentsCount()},
                {"electrodes_count", input->getGenericElectrodesCount()},
                {"ground_node", input->getGroundNode()}};
        }
    };
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        Pybind11 eitannealing solver plugin
        -----------------------

        .. currentmodule:: pyeitsolver
    )pbdoc";
    py::class_<pyeitsolver::EitSolver>(m, "EitSolver", "Loads mesh data, currents data and creates problem matrices")
        .def(py::init<std::string, std::string>())
        .def_readonly("info", &pyeitsolver::EitSolver::info)
        .def("solve_forward_problem", &pyeitsolver::EitSolver::solve_forward_problem, "Solves forward problem, returning potentials, for a given conductivity distribution", py::arg("conds"), py::kw_only(), py::arg("mesh_potentials") = false);

#ifdef VERSION_INFO
    m.attr("__version__") = STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}