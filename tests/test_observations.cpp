#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "parameters/parametersparser.h"
#include "matrixview.h"
#include "observations.h"
#include <iostream>
#include <gtest/gtest.h>

// Check if b and A matrix are compatible for CG solving
TEST(CgTwoDimTest, SizeAssertion)
{
    bool is2dProblem;
    const char *meshFile = "./data/circular_A_2D.msh";
    const char *currentsFile = "./data/cuba_190ma_cp.txt";
    const char *measurementsFile = "./data/tri02b.ampl_calibrado.txt";
    std::shared_ptr<problem> input = problem::createNewProblem(meshFile, &is2dProblem);
    input->initProblem(meshFile);
    observations<double> *readings = new observations<double>;
    readings->initObs(&currentsFile, measurementsFile, input->getNodesCount(), input->getGenericElectrodesCount());
    input->buildNodeCoefficients();
    input->prepareSkeletonMatrix();
    input->createCoef2KMatrix();

    matrix *m1;
    Eigen::VectorXd v = Eigen::VectorXd::Constant(input->getNumCoefficients(), 0.3815);
    input->assembleProblemMatrix(&v[0], &m1);
    input->postAssembleProblemMatrix(&m1);
    Eigen::VectorXd currents;
    Eigen::VectorXd x;
    SparseIncompleteLLT precond(*m1);

    currents = input->getCurrentVector(0, readings);
    EXPECT_EQ(currents.rows(), m1->rows());
}
