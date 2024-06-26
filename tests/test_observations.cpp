#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "parameters/parametersparser.h"
#include "matrixview.h"
#include "observations.h"
#include <iostream>
#include <gtest/gtest.h>

class CgTwoDimTest : public testing::Test
{
protected:
    CgTwoDimTest()
    {
        bool is2dProblem;
        const char *meshFile = "./data/circular_A_2D.msh";
        const char *currentsFile = "./data/cuba_190ma_cp.txt";
        const char *measurementsFile = "./data/tri02b.ampl_calibrado.txt";
        this->input = problem::createNewProblem(meshFile, &is2dProblem);
        this->input->initProblem(meshFile);
        this->readings = new observations<double>;
        this->readings->initObs(&currentsFile, measurementsFile, this->input->getNodesCount(), this->input->getGenericElectrodesCount());
        this->input->buildNodeCoefficients();
        this->input->prepareSkeletonMatrix();
        this->input->createCoef2KMatrix();
    }

    std::shared_ptr<problem> input;
    observations<double> *readings;
};

// Check if b and A matrix are compatible for CG solving
TEST_F(CgTwoDimTest, SizeAssertion)
{
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
