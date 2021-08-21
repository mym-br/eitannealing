#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"
#include "../src/incomplete_qr_complex2.h"
#include "util/timestamp/timestamp.h"
#include "solver_lb.h"

int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <2) {
    	 std::cerr << "need mesh filename\n";
    	 return 1;
     }
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

     std::unique_ptr<LB_Solver_Complex::Preconditioner> pre =
        std::make_unique<LB_Solver_Complex::Preconditioner>(8, 10, *Aii_R, *Aii_I, *Aic_R, *Aic_I);

     std::unique_ptr<LB_Solver_Complex> solver =
        std::make_unique<LB_Solver_Complex>(Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_I, J_R, J_I, V_R, V_I, *pre, 0.002);

     std::cout << "Iterations:\n";
     for(int i=0; i<100;i++) {
         solver->do_iteration();
         std::cout << i << ":" << solver->getMinErrorl2Estimate() << ":" << solver->getMaxErrorl2Estimate() << std::endl;
     }

     std::cout << "Benchmarking old vs new...\n";
     matrixcomplex Aii = (Eigen::MatrixXcd(*Aii_R) + Complex(0,1)*Eigen::MatrixXcd(*Aii_I)).sparseView();
     matrixcomplex Aic = (Eigen::MatrixXcd(*Aic_R) + Complex(0,1)*Eigen::MatrixXcd(*Aic_I)).sparseView();
     std::unique_ptr<SparseIncompleteQRComplex2> pre2 = std::make_unique<SparseIncompleteQRComplex2>(8, 10, Aii, Aic);


     Eigen::VectorXd dummy_Re(Eigen::VectorXd::Random(Aii_R->cols()));
     Eigen::VectorXd dummy_Im(Eigen::VectorXd::Random(Aii_R->cols()));
     long start, stop;
     for(int i = 0; i<100; i++) {
         pre->solveInPlaceC(dummy_Re, dummy_Im);
         pre->solveInPlaceCT(dummy_Re, dummy_Im);
     }
     start = get_usec_timestamp();
     for(int i = 0; i<400000; i++) {
         pre->solveInPlaceC(dummy_Re, dummy_Im);
         pre->solveInPlaceCT(dummy_Re, dummy_Im);
     }
     stop = get_usec_timestamp();
     std::cout << "old: "  <<  ((double)(stop - start))/400000 << std::endl;

     Eigen::VectorXcd dummy_C(Eigen::VectorXcd::Random(Aii_R->cols()));
     for(int i = 0; i<100; i++) {
         pre2->solveInPlace(dummy_C);
         pre2->solveInPlaceT(dummy_C);
     }
     start = get_usec_timestamp();
     for(int i = 0; i<400000; i++) {
         pre2->solveInPlace(dummy_C);
         pre2->solveInPlaceT(dummy_C);
     }
     stop = get_usec_timestamp();
     std::cout << "new: "  <<  ((double)(stop - start))/400000 << std::endl;

     return 0;
 }
