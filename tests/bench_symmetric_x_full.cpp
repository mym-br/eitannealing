#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"

#include <sys/time.h>
#include <sys/resource.h>
#include "solver_lb.h"
unsigned long get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec*1e6 + t.tv_usec;
}


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
            c.insertBack(itc.row(), i) = itc.value();
        }     
     }
     c.finalize();
     
     Eigen::VectorXd  x(Eigen::VectorXd::Random(c.cols())), y, y2, x2(Eigen::VectorXd::Random(c.rows()));
     
     long start, stop;
     for(int i = 0; i<100; i++)
            y.noalias() = c*x;
     start = get_time();
     for(int i = 0; i<400000; i++)
            y.noalias() = c*x;
     stop = get_time();
     std::cout << "full: "  <<  ((double)(stop - start))/400000 << std::endl;

     for(int i = 0; i<100; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x;
            y2.noalias() = (*Aic_R)*x;
     }
     start = get_time();
     for(int i = 0; i<400000; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x;
            y2.noalias() = (*Aic_R)*x;
     }
     stop = get_time();
     std::cout << "symmetric: "  <<  ((double)(stop - start))/400000 << std::endl;
     
     for(int i = 0; i<100; i++)
            y.noalias() = c.transpose()*x2;
     start = get_time();
     for(int i = 0; i<400000; i++)
            y.noalias() = c.transpose()*x2;
     stop = get_time();
     std::cout << "transposed full: "  <<  ((double)(stop - start))/400000 << std::endl;
     
     y2 = Eigen::VectorXd::Random(Aic_R->rows());
     for(int i = 0; i<100; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x + Aic_R->transpose()*y2;
            //y += Aic_R->transpose()*y2;
     }
     start = get_time();
     for(int i = 0; i<400000; i++) {
            y.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x;
            //y += Aic_R->transpose()*y2;
     }
     stop = get_time();
     std::cout << "transposed symmetric: "  <<  ((double)(stop - start))/400000 << std::endl;


     /*

     Eigen::VectorXd V_R(32), V_I(32), J_R(32), J_I(32);

     V_R.fill(0.0);
     V_I.fill(0.0);
     J_R.fill(0.0);
     J_I.fill(0.0);
     V_R(4) = -30.738288239574686;
     V_R(0) = 30.733620796726793;
     V_R(5) = -12.495870755237416;
     V_R(31) = 12.49545757830731;
     V_R(1) = 8.627724263627385;
     V_R(3) = -8.624816166327843;

     V_I(4) = 27.881561780188846;
     V_I(0) = -27.87684498288155;
     V_I(5) = 12.393968674438954;
     V_I(31) = -12.3936017438714986;
     V_I(1) = -8.53400765435308;
     V_I(3) = 8.531044976628754;

     J_R(0) = 1.0;
     J_R(4) = -1.0;



     LB_Solver_Complex::Preconditioner *pre;
     pre = new SparseIncompleteQRComplex(8, 10, *Aii_R, *Aii_I, *Aic_R, *Aic_I);


     LB_Solver_Complex *solver;
     solver = new LB_Solver_Complex(Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_I, J_R, J_I, V_R, V_I, *pre, 0.002);

     for(int i=0; i<50;i++)
         solver->do_iteration();

     solver->do_iteration();
     */



     return 0;
 }
