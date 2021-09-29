#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/intcoef.h"
#include "../src/solver_lb.h"
#include "../src/incomplete_qr_builder.h"
#include "../src/sparse_incomplete_qr.h"

int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <4) {
    	 std::cerr << "need mesh filename, currents and tensions\n";
    	 return 1;
     }
     std::cout << "Parsing mesh file..." << std::flush;
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
     observations<double> *readingsScalar;
     input->initProblem(argv[1]);
     input->buildNodeCoefficients();
     input->prepareSkeletonMatrix();
     input->createCoef2KMatrix();
     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;

     std::cout << "Parsing observations..." << std::flush;
     readingsScalar = new observations<double>;
     readingsScalar->initObs((const char **)&argv[2], argv[3], input->getNodesCount(), input->getGenericElectrodesCount());
     std::cout << "Done\n" << std::flush;

     std::cout << "Preparing dummy solution..." << std::flush;
     double *sol = new double[input->getNumCoefficients()];
     for(int i=0;i<input->getNumCoefficients();i++) sol[i]=(mincond+maxcond)/2;
     std::cout << "Done\n" << std::flush;

     matrix *Aii, *Aic, *Acc;
     std::cout << "Assembling matrices..." << std::flush;
     assembleProblemMatrix_lb(sol, &Aii, &Aic, &Acc, *input);
     std::cout << "Done\n" << std::flush;

     std::unique_ptr<LB_Solver::Preconditioner> precond;
     std::cout << "Preparing preconditioner..." << std::flush;
     precond.reset(LB_Solver::makePreconditioner(*Aii, *Aic));
     std::cout << "Done\n" << std::flush;

     std::cout << "Preparing solver (eigenvalue estimator)..." << std::flush;
     std::unique_ptr<LB_Solver_EG_Estimate> solver_e;
     solver_e.reset(new LB_Solver_EG_Estimate(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, 600, (float)0.0000001
     ));
     std::cout << "Done\n" << std::flush;

     std::cout << "Preparing solver (new eigenvalue estimator)..." << std::flush;
     Eigen::VectorXd eigenvector;
     double eigenvalue;
     LB_Solver(Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, 600, (float)0.0000001, eigenvalue, eigenvector);
     std::cout << "Done\n" << std::flush;

     std::cout << "Eigenvalue (old): " << solver_e->getLeastEvEst() << " (new): " << eigenvalue << std::endl;
     std::cout << "norm of old-new eigenvector: " << (eigenvector-solver_e->getEvector()).norm() << std::endl;

     std::cout << "Testing emplace construction on a vector... ";
     {
        std::vector<LB_Solver> v;
        v.reserve(100);
        for(int i=0;i<100;i++) {
            v.emplace_back(Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, 600, (float)0.0000001, eigenvalue, eigenvector);
        }
        std::cout << "Created 100 instances on a vector... now deleting... ";
     }
     std::cout << "Done!\n";
     return 0;
 }
