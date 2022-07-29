#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"

#include "util/timestamp/timestamp.h"
#include "solver_lb.h"


int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <2) {
    	 std::cerr << "need mesh filename\n";
    	 return 1;
     }
     //Eigen::setNbThreads(16);
     std::cout << "Parsing mesh file..." << std::flush;
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
     input->initProblem(argv[1]);
     input->buildNodeCoefficients();
     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;
     std::cout << "Preparing dummy solution..." << std::flush;
     double *sol_R = new double[input->getNumCoefficients()];
     for(int i=0;i<input->getNumCoefficients();i++) sol_R[i]=1.0;


     double *sol_I = new double[input->getNumCoefficients()];
     for(int i=0;i<32;i++) sol_I[i]=0.0;
     for(int i=32;i<input->getNumCoefficients();i++) sol_I[i]=1.0;
     std::cout << "Done\n" << std::flush;
     matrix *Aii_R, *Aii_I, *Aic_R, *Aic_I, *Acc_R, *Acc_I;
     std::cout << "Assembling matrices..." << std::flush;
     assembleProblemMatrix_lb(sol_R, &Aii_R, &Aic_R, &Acc_R, *input);
     assembleProblemMatrix_lb(sol_I, &Aii_I, &Aic_I, &Acc_I, *input);
     std::cout << "Done\n" << std::flush;

     matrix ii(Aii_R->selfadjointView<Eigen::Lower>());
     matrix c(Aii_R->rows()+Aic_R->rows(), Aii_R->cols());
     
     c.reserve(ii.nonZeros() + Aic_R->nonZeros());
     for(int i = 0; i < Aii_R->cols(); i++) {
        c.startVec(i);
        for(Eigen::SparseMatrix<double>::InnerIterator iti(ii, i); iti; ++iti) {
            c.insertBack(iti.row(), i) = iti.value();
        }
        for(Eigen::SparseMatrix<double>::InnerIterator itc(*Aic_R, i); itc; ++itc) {
            c.insertBack(ii.rows()+itc.row(), i) = itc.value();
        }     
     }
     c.finalize();
     
     Eigen::VectorXd  x(Eigen::VectorXd::Random(c.cols())), y, y2, x2(Eigen::VectorXd::Random(c.rows()));
     
     long start, stop;
     for(int i = 0; i<100; i++)
            y.noalias() = c*x;
     start = get_usec_timestamp();
     for(int i = 0; i<800000; i++)
            y.noalias() = c*x;
     stop = get_usec_timestamp();
     std::cout << "full: "  <<  ((double)(stop - start))/800000 << std::endl;

     for(int i = 0; i<100; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x;
            y2.noalias() = (*Aic_R)*x;
     }
     start = get_usec_timestamp();
     for(int i = 0; i<800000; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x;
            y2.noalias() = (*Aic_R)*x;
     }
     stop = get_usec_timestamp();
     std::cout << "symmetric: "  <<  ((double)(stop - start))/800000 << std::endl;
     
     for(int i = 0; i<100; i++)
            y.noalias() = c.transpose()*x2;
     start = get_usec_timestamp();
     for(int i = 0; i<800000; i++)
            y.noalias() = c.transpose()*x2;
     stop = get_usec_timestamp();
     std::cout << "transposed full: "  <<  ((double)(stop - start))/800000 << std::endl;
     
     y2 = Eigen::VectorXd::Random(Aic_R->rows());
     for(int i = 0; i<100; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x + Aic_R->transpose()*y2;
            //y += Aic_R->transpose()*y2;
     }
     start = get_usec_timestamp();
     for(int i = 0; i<800000; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x;
            //y += Aic_R->transpose()*y2;
     }
     stop = get_usec_timestamp();
     std::cout << "transposed symmetric: "  <<  ((double)(stop - start))/800000 << std::endl;


    

     return 0;
 }
