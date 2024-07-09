#include <pybind11/pybind11.h>
#include <pybind11/eigen/matrix.h>
#include <pybind11/stl.h>
#include <filesystem>
#include <Eigen/Core>
#include <memory>
#include <map>
#include <string>
#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"

#define STRINGIFY(x) #x
#define DEFAULT_MAX_ITERATIONS 1500
#define DEFAULT_RESIDUAL 1e-19

namespace fs = std::filesystem;

namespace pyeitsolver
{
    class EitSolver {
        public:
        EitSolver(std::string meshFilename, std::string currentFilename) {
            init(meshFilename, currentFilename);
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

        std::pair<std::map<std::string, double>, Eigen::MatrixXd> solve_forward_problem(Eigen::VectorXd &conductivities, bool meshPotentials = false, int maxIterations = DEFAULT_MAX_ITERATIONS, double residual = DEFAULT_RESIDUAL)
        {
            // Check conductivity vector size
            auto n = conductivities.size();
            if (n != input->getNumCoefficients())
                throw std::runtime_error("Wrong conductivities vector size " + std::to_string(n) + " (should be " + std::to_string(input->getNumCoefficients()) + ")");

            // Map node conductivities to coefficient indices
            Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
            Eigen::VectorXi indices = createIndicesVector(input);
            v(indices) = conductivities;

            // Create FEM conductivity matrix
            matrix *m1;
            input->assembleProblemMatrix(&v[0], &m1);
            input->postAssembleProblemMatrix(&m1);

            // Create preconditioner matrix
            std::shared_ptr<SparseIncompleteLLT> precond = std::shared_ptr<SparseIncompleteLLT>(new SparseIncompleteLLT(*m1));

            // Solve forward problem to obtain potentials
            int electrodeCount = input->getGenericElectrodesCount();
            int groundNodeIdx = input->getGroundNode();

            Eigen::MatrixXd potentials(electrodeCount, meshPotentials ? n : electrodeCount);
            Eigen::VectorXd x, currents;

            double avgIterations = 0, avgResidual = 0;
            for (int patterno = 0; patterno < readings->getCurrentsCount(); patterno++)
            {
                currents = input->getCurrentVector(patterno, readings.get());
                CG_Solver solver(*m1, currents, *precond);
                int i = 0;
                for (; i < maxIterations && solver.getResidueSquaredNorm() > residual; i++)
                    solver.do_iteration();
                avgIterations += i;
                avgResidual += solver.getResidueSquaredNorm();

                x = solver.getX();

                // Create xExpanded with size n
                Eigen::VectorXd xExpanded(n);
                xExpanded << x.head(groundNodeIdx), 0, x.tail(n - 1 - groundNodeIdx);
                
                // Save results to appropriate index in the output vector
                potentials.row(patterno) = (meshPotentials ? xExpanded : xExpanded.segment(groundNodeIdx - electrodeCount + 1, electrodeCount)).transpose();
                potentials.row(patterno) *= readings->getCurrentVal(patterno);
            }
            avgIterations /= readings->getCurrentsCount();
            avgResidual /= readings->getCurrentsCount();

            std::map<std::string, double> executionInfo = {
                {"max_iterations", (double)maxIterations},
                {"residual", (double)residual},
                {"avg_iterations", avgIterations},
                {"avg_residual", avgResidual}};

            // Return potentials
            delete m1;
            return std::pair<std::map<std::string, double>, Eigen::MatrixXd>(executionInfo, potentials);
        }

        std::map<std::string, int> getProblemInfo() {
            return std::map<std::string, int>{
                {"nodes_count", input->getNodesCount()},
                {"currents_count", readings->getCurrentsCount()},
                {"electrodes_count", input->getGenericElectrodesCount()},
                {"ground_node", input->getGroundNode()},
                {"first_electrode_idx", input->getGroundNode() - input->getGenericElectrodesCount() + 1}};
        }

        int nodeCount;
        int electrodeCount;

        private:
        std::shared_ptr<problem> input;
        std::unique_ptr<observations<double>> readings;

        void init(std::string meshFilename, std::string currentFilename)
        {
            // Check if the mesh file exists
            if (!fs::exists(meshFilename)) {
                throw std::runtime_error("Mesh file does not exist: " + meshFilename);
            }

            // Check if the current file exists
            if (!fs::exists(currentFilename)) {
                throw std::runtime_error("Current file does not exist: " + currentFilename);
            }

            
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

            assert(input->getNodesCount() == input->getNumCoefficients() && "There should be a 1 to 1 correspondence between node and coefficient");

            this->nodeCount = input->getNodesCount();
            this->electrodeCount = input->getGenericElectrodesCount();
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
        .def_property_readonly("info", &pyeitsolver::EitSolver::getProblemInfo)
        .def("solve_forward_problem", &pyeitsolver::EitSolver::solve_forward_problem, "Solves forward problem, returning potentials, for a given conductivity distribution", py::arg("conds"), py::kw_only(), py::arg("mesh_potentials") = false, py::arg("max_iterations") = DEFAULT_MAX_ITERATIONS, py::arg("residual") = DEFAULT_RESIDUAL)
                .def("__repr__",
        [](const pyeitsolver::EitSolver &a) {
            return "<pyeitsolver.EitSolver of mesh with " + std::to_string(a.nodeCount) + " coefficients and " +  std::to_string(a.electrodeCount) + " electrodes>";
        });

#ifdef VERSION_INFO
    m.attr("__version__") = STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}