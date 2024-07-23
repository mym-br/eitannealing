#include <pybind11/pybind11.h>
#include <pybind11/eigen/matrix.h>
#include <pybind11/stl.h>
#include <filesystem>
#include <Eigen/Core>
#ifdef USE_PETSC
#include <petscksp.h>
#endif
#include <memory>
#include <tuple>
#include <map>
#include <string>

#include <mpi.h>

#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"

// Enable only one at a time.
//#define USE_MUMPS 1
//#define USE_MKL_PARDISO 1
//#define USE_CUDA_SOLVER 1

#define STRINGIFY(x) #x
#define CHKERRTHROW(ierr)                                     \
    if (PetscUnlikely(ierr))                                  \
    {                                                         \
        throw std::runtime_error(PetscErrorMessageStr(ierr)); \
    }

#define CURRENT_FREQUENCY 275000
#define PARASITE_CAPACITANCY 80E-12

namespace fs = std::filesystem;

namespace pyeitsolver
{
#ifdef USE_PETSC
    std::string PetscErrorMessageStr(PetscErrorCode ierr)
    {
        const char *msg;
        PetscErrorMessage(ierr, &msg, PETSC_NULLPTR);
        return std::string(msg);
    }
#endif

    class EitComplexSolver
    {
    public:
        EitComplexSolver(std::string meshFilename, std::string currentFilename)
        {
            init(meshFilename, currentFilename);
#ifdef USE_PETSC
            char a0[] = "";
#ifdef USE_MKL_PARDISO
            char pardiso_a1[] = "-pc_type";
            char pardiso_a2[] = "lu";
            char pardiso_a3[] = "-pc_factor_mat_solver_type";
            char pardiso_a4[] = "mkl_pardiso";
            char pardiso_a5[] = "-mat_mkl_pardiso_65"; // Suggested number of threads to use within MKL PARDISO
            char pardiso_a6[] = "4";
            char pardiso_a7[] = "-mat_mkl_pardiso_68"; // Message level information, use 1 to get detailed information on the solver options
            char pardiso_a8[] = "0";
#endif
#ifdef USE_CUDA_SOLVER
            char cuda_a1[] = "-vec_type";
            char cuda_a2[] = "cuda";
            char cuda_a3[] = "-mat_type";
            char cuda_a4[] = "aijcusparse";
#endif
            char petsc_a1[] = "";//"-info";
            char petsc_a2[] = "";//"-log_trace";
            char* argv[] = {
                &a0[0],
#ifdef USE_MKL_PARDISO
                &pardiso_a1[0],
                &pardiso_a2[0],
                &pardiso_a3[0],
                &pardiso_a4[0],
                &pardiso_a5[0],
                &pardiso_a6[0],
                &pardiso_a7[0],
                &pardiso_a8[0],
#endif
#ifdef USE_CUDA_SOLVER
                &cuda_a1[0],
                &cuda_a2[0],
                &cuda_a3[0],
                &cuda_a4[0],
#endif
                &petsc_a1[0],
                &petsc_a2[0],
                NULL
            };
            char** argv_ptr = &argv[0];
            int argc = (int) (sizeof(argv) / sizeof(argv[0])) - 1;
            PetscErrorCode ierr = PetscInitialize(&argc, &argv_ptr, PETSC_NULLPTR, PETSC_NULLPTR);
            //PetscErrorCode ierr = PetscInitialize(PETSC_NULLPTR, PETSC_NULLPTR, PETSC_NULLPTR, PETSC_NULLPTR);
            if (ierr)
            {
                throw std::runtime_error("Could not initialize PETSc");
            }
#endif
        }

        ~EitComplexSolver()
        {
#ifdef USE_PETSC
            PetscErrorCode ierr = PetscFinalize();
#endif
        }

        std::map<std::string, int> getProblemInfo()
        {
            return std::map<std::string, int>{
                {"nodes_count", input->getNodesCount()},
                {"currents_count", readings->getCurrentsCount()},
                {"electrodes_count", input->getGenericElectrodesCount()},
                {"ground_node", input->getGroundNode()},
                {"first_electrode_idx", input->getGroundNode() - input->getGenericElectrodesCount() + 1}};
        }

        std::tuple<Eigen::VectorXi, Eigen::VectorXi, Eigen::VectorXcd> getStiffnessCOO(Eigen::VectorXd &conductivities)
        {
            // Check conductivity vector size
            auto n = conductivities.size();
            if (n != input->getNumCoefficients())
                throw std::runtime_error("Wrong conductivities vector size " + std::to_string(n) + " (should be " + std::to_string(input->getNumCoefficients()) + ")");

            // Map node conductivities to coefficient indices
            Eigen::VectorXcd v = Eigen::VectorXcd::Zero(n);
            Eigen::VectorXi indices = createIndicesVector(input);
            v(indices) = conductivities;

            // Create FEM conductivity matrix
            matrixcomplex *m1;
            input->assembleProblemMatrix(&v[0], &m1);
            input->addMatrixCapacitances(&m1);
            // input->postAssembleProblemMatrix(&m1);

            Eigen::VectorXi rows(m1->nonZeros());
            Eigen::VectorXi cols(m1->nonZeros());
            Eigen::VectorXcd values(m1->nonZeros());

            int i = 0;
            for (int k = 0; k < m1->outerSize(); ++k)
                for (matrixcomplex::InnerIterator it(*m1, k); it; ++it)
                {
                    rows[i] = it.row();
                    cols[i] = it.col();
                    values[i] = it.value();
                    i++;
                }

            return std::tuple<Eigen::VectorXi, Eigen::VectorXi, Eigen::VectorXcd>(rows, cols, values);
        }

        Eigen::MatrixXd getCurrentsVectors()
        {
            int numCols = input->getCurrentVector(0, readings.get()).size();
            int numRows = readings->getCurrentsCount();
            Eigen::MatrixXd currentsMatrix(numRows, numCols);

            for (int pattern = 0; pattern < numRows; ++pattern)
            {
                currentsMatrix.row(pattern) = input->getCurrentVector(pattern, readings.get());
            }
            return currentsMatrix;
        }

#ifdef USE_PETSC
        Eigen::MatrixXcd solveForwardProblem(Eigen::VectorXcd &conductivities, bool meshPotentials = false)
        {
            Vec x, b; /* solution, right-hand-side */
            Mat A;    /* linear system matrix */
            KSP ksp;  /* linear solver context */
            PetscErrorCode ierr;

            int mpi_rank;
            int mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            if (mpi_err != MPI_SUCCESS) {
                throw std::runtime_error("Could not get the MPI rank (error: " + std::to_string(mpi_err) + ").\n");
            }
            int mpi_size;
            mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            if (mpi_err != MPI_SUCCESS) {
                throw std::runtime_error("Could not get the MPI size (error: " + std::to_string(mpi_err) + ").\n");
            }
            printf("MPI rank = %d, size = %d\n", mpi_rank, mpi_size);

            // Check conductivity vector size
            auto n = conductivities.size();
            if (n != input->getNumCoefficients())
                throw std::runtime_error("Wrong conductivities vector size " + std::to_string(n) + " (should be " + std::to_string(input->getNumCoefficients()) + ")");

            // Map node conductivities to coefficient indices
            Eigen::VectorXcd v = Eigen::VectorXcd::Zero(n);
            Eigen::VectorXi indices = createIndicesVector(input);
            v(indices) = conductivities;

            // Create FEM conductivity matrix
            matrixcomplex *m1;
            input->assembleProblemMatrix(&v[0], &m1);
            input->addMatrixCapacitances(&m1);

            // Create vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &x);
            CHKERRTHROW(ierr);
            ierr = VecSetSizes(x, PETSC_DECIDE, n);
            CHKERRTHROW(ierr);
            ierr = VecSetFromOptions(x);
            CHKERRTHROW(ierr);
            ierr = VecDuplicate(x, &b);
            CHKERRTHROW(ierr);

            // Create matrix
            ierr = MatCreate(PETSC_COMM_WORLD, &A);
            CHKERRTHROW(ierr);
            ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
            CHKERRTHROW(ierr);
            ierr = MatSetFromOptions(A);
            CHKERRTHROW(ierr);
            ierr = MatSetUp(A);
            CHKERRTHROW(ierr);

            // Assemble matrix A from Eigen::SparseMatrix
            for (int k = 0; k < m1->outerSize(); ++k)
            {
                for (matrixcomplex::InnerIterator it(*m1, k); it; ++it)
                {
                    PetscInt row = it.row();
                    PetscInt col = it.col();
                    PetscScalar value = it.value();
                    ierr = MatSetValues(A, 1, &row, 1, &col, &value, INSERT_VALUES);
                    CHKERRTHROW(ierr);
                    if (row == col)
                        continue;
                    ierr = MatSetValues(A, 1, &col, 1, &row, &value, INSERT_VALUES);
                    CHKERRTHROW(ierr);
                }
            }

            // Assemble the matrix
            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            CHKERRTHROW(ierr);
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            CHKERRTHROW(ierr);

            // Create linear solver context
            ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
            CHKERRTHROW(ierr);
            ierr = KSPSetOperators(ksp, A, A);
            CHKERRTHROW(ierr);

            PC pc; // (PETSc) preconditioner context
            ierr = KSPGetPC(ksp, &pc); // get preconditioner context
            CHKERRTHROW(ierr);

#ifdef USE_MUMPS
            ierr = KSPSetType(ksp, KSPPREONLY);
            CHKERRTHROW(ierr);

            // Cholesky factorization/decomposition requires a Hermitian, positive-definite matrix.
            // PETSc(PCSetType(pc, PCCHOLESKY)); // slow KSPSetUp()
            ierr = PCSetType(pc, PCLU);
            CHKERRTHROW(ierr);
            ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
            CHKERRTHROW(ierr);
#else
//            ierr = PCSetType(pc, PCJACOBI);
//            //ierr = PCSetType(pc, PCNONE);
//            CHKERRTHROW(ierr);
#endif

            // Set linear solver options
            ierr = KSPSetFromOptions(ksp);
            CHKERRTHROW(ierr);

            ierr = KSPSetUp(ksp);
            CHKERRTHROW(ierr);

            Eigen::MatrixXcd potentials(electrodeCount, n);
            for (int pattern = 0; pattern < readings->getCurrentsCount(); pattern++)
            {
                // Set right-hand side vector from Eigen::VectorXcd
                Eigen::VectorXd b_eigen = input->getCurrentVector(pattern, readings.get());
                for (PetscInt i = 0; i < n; ++i)
                {
                    std::complex<double> value = b_eigen[i];
                    ierr = VecSetValue(b, i, value, INSERT_VALUES);
                    CHKERRTHROW(ierr);
                }

                ierr = VecAssemblyBegin(b);
                CHKERRTHROW(ierr);
                ierr = VecAssemblyEnd(b);
                CHKERRTHROW(ierr);

                // Solve linear system
                ierr = KSPSolve(ksp, b, x);
                CHKERRTHROW(ierr);

                // // View the solution
                // ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);
                // CHKERRTHROW(ierr);

                // Convert PETSc vector x to Eigen::VectorXcd
                PetscInt size;
                ierr = VecGetSize(x, &size);
                CHKERRTHROW(ierr);

                const PetscScalar *x_array;
                ierr = VecGetArrayRead(x, &x_array);
                CHKERRTHROW(ierr);

                for (PetscInt i = 0; i < size; ++i)
                {
                    potentials.row(pattern)[i] = x_array[i];
                }
            }

            // Free work space
            ierr = VecDestroy(&x);
            CHKERRTHROW(ierr);
            ierr = VecDestroy(&b);
            CHKERRTHROW(ierr);
            ierr = MatDestroy(&A);
            CHKERRTHROW(ierr);
            ierr = KSPDestroy(&ksp);
            CHKERRTHROW(ierr);

            return potentials;
        }
#endif

        int nodeCount;
        int electrodeCount;
        double current;

    private:
        std::shared_ptr<problem> input;
        std::unique_ptr<observations<double>> readings;
        void init(std::string meshFilename, std::string currentFilename)
        {
            // Check if the mesh file exists
            if (!fs::exists(meshFilename))
            {
                throw std::runtime_error("Mesh file does not exist: " + meshFilename);
            }

            // Check if the current file exists
            if (!fs::exists(currentFilename))
            {
                throw std::runtime_error("Current file does not exist: " + currentFilename);
            }

            // Load mesh file geometry into memory
            bool is2dProblem;
            input = problem::createNewProblem(meshFilename.c_str(), &is2dProblem);
            input->setIgnoreouterring(true);
            input->setCapacitance(PARASITE_CAPACITANCY);
            input->setCurrentFreq(CURRENT_FREQUENCY);
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
            this->current = readings->getCurrentVal(0);
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
    };
}

namespace py = pybind11;

PYBIND11_MODULE(pyeitsolver_complex, m)
{
    m.doc() = R"pbdoc(
        Pybind11 eitannealing complex solver plugin
        -----------------------

        .. currentmodule:: pyeitsolver
    )pbdoc";
    py::class_<pyeitsolver::EitComplexSolver>(m, "EitComplexSolver", "Loads mesh data, currents data and creates problem matrices")
        .def(py::init<std::string, std::string>())
        .def_property_readonly("info", &pyeitsolver::EitComplexSolver::getProblemInfo)
        .def_readwrite("current", &pyeitsolver::EitComplexSolver::current)
        .def("getCOO_formatted_stiffness", &pyeitsolver::EitComplexSolver::getStiffnessCOO)
        .def("get_currents_vectors", &pyeitsolver::EitComplexSolver::getCurrentsVectors)
#ifdef USE_PETSC
        .def("solve_forward_problem", &pyeitsolver::EitComplexSolver::solveForwardProblem, "Solves forward problem, returning potentials, for a given conductivity distribution", py::arg("conds"), py::kw_only(), py::arg("mesh_potentials") = false)
#endif
        .def("__repr__",
             [](const pyeitsolver::EitComplexSolver &a)
             {
                 return "<pyeitsolver.EitComplexSolver of mesh with " + std::to_string(a.nodeCount) + " coefficients and " + std::to_string(a.electrodeCount) + " electrodes>";
             });

#ifdef VERSION_INFO
    m.attr("__version__") = STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
