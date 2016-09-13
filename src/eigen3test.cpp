#include <iostream>
#include <vector>
#include "basematrix.h"
#include "solver.h"
#include "incomplete_cholesky.h"
#include "incomplete_LDLT.h"
#include "assembleMatrix.h"
#include "problemdescription.h"
typedef Eigen::Triplet<double> T;




int main(int argc, char *argv[])
{
    initProblem(argv[1]);
    buildNodeCoefficients();
    prepareSkeletonMatrix();
    createCoef2KMatrix();
    matrix *m1, *m2;
    Eigen::VectorXd v(numcoefficients);
    for(int i=0; i<v.rows(); i++)
      v[i] = 0.1 + 0.5*(i%5);
    assembleProblemMatrix_old(&v[0], &m1);
    
    SparseIncompleteLLT pre(*m1);
    
    Eigen::VectorXd b((m1->selfadjointView<Eigen::Lower>())*Eigen::VectorXd::Ones(m1->rows()));
    
    CG_Solver solver(*m1, b, pre);
    
    return 0;
}