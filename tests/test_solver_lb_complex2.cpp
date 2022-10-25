#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/solver_lb_complex2.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"


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
     matrixcomplex Aii = (Eigen::MatrixXcd(*Aii_R) + Complex(0,1)*Eigen::MatrixXcd(*Aii_I)).sparseView();
     matrixcomplex Aic = (Eigen::MatrixXcd(*Aic_R) + Complex(0,1)*Eigen::MatrixXcd(*Aic_I)).sparseView();
     matrixcomplex Acc = (Eigen::MatrixXcd(*Acc_R) + Complex(0,1)*Eigen::MatrixXcd(*Acc_I)).sparseView();
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


     eigen_complexdouble_engine::vector V = (eigen_complexdouble_engine::vector(V_R) + Complex(0,1)*eigen_complexdouble_engine::vector(V_I));
     eigen_complexdouble_engine::vector J = (eigen_complexdouble_engine::vector(J_R) + Complex(0,1)*eigen_complexdouble_engine::vector(J_I));

     std::unique_ptr<LB_Solver_Complex::Preconditioner> pre =
        std::make_unique<LB_Solver_Complex::Preconditioner>(8, 16, *Aii_R, *Aii_I, *Aic_R, *Aic_I);


     std::unique_ptr<eigen_complexdouble_engine::preconditioner> pre2(
        eigen_complexdouble_engine::make_new_preconditioner(Aii, Aic));

     Eigen::VectorXd r_R = -Aic_R->transpose()*V_R + Aic_I->transpose()*V_I;
     Eigen::VectorXd r_I = -Aic_R->transpose()*V_I - Aic_I->transpose()*V_R;



     Eigen::VectorXd rc_R = J_R - Acc_R->selfadjointView<Eigen::Lower>()*V_R + Acc_I->selfadjointView<Eigen::Lower>()*V_I;
     Eigen::VectorXd rc_I = J_I - Acc_R->selfadjointView<Eigen::Lower>()*V_I - Acc_I->selfadjointView<Eigen::Lower>()*V_R;

     double JhatNorm2 = r_R.squaredNorm()+r_I.squaredNorm()+rc_R.squaredNorm()+rc_I.squaredNorm();
     double delta = sqrt(JhatNorm2);

     Eigen::VectorXd p_R = r_R/delta;
     Eigen::VectorXd p_I = r_I/delta;
     Eigen::VectorXd pc_R = rc_R/delta;
     Eigen::VectorXd pc_I = rc_I/delta;


     Eigen::VectorXd s_R = Aii_R->selfadjointView<Eigen::Lower>()*p_R + Aii_I->selfadjointView<Eigen::Lower>()*p_I;
     Eigen::VectorXd s_I = Aii_R->selfadjointView<Eigen::Lower>()*p_I - Aii_I->selfadjointView<Eigen::Lower>()*p_R;
     s_R.noalias() += Aic_R->transpose()*pc_R;
     s_R.noalias() += Aic_I->transpose()*pc_I;
     s_I.noalias() += Aic_R->transpose()*pc_I;
     s_I.noalias() -= Aic_I->transpose()*pc_R;


     eigen_complexdouble_engine::vector r, rc;
     eigen_complexdouble_engine::j_minus_a_x_phi(r, rc, Aic, Acc, J, V);
     double JhatNorm2_2 = r.squaredNorm()+rc.squaredNorm();
     double delta_2 = sqrt(JhatNorm2_2);
     eigen_complexdouble_engine::vector p = r/delta_2;
     eigen_complexdouble_engine::vector pc = rc/delta_2;

     eigen_complexdouble_engine::vector s;

     eigen_complexdouble_engine::symmetric_conjtranspose_complex_mult_and_assign(Aii, p, s);
     s.noalias() += Aic.adjoint()*pc;

     //std::cout << s - s_R - s_I*Complex(0,1) << std::endl;
     pre->solveInPlaceCT(s_R, s_I);
     pre2->solveInPlaceT(s);


     //std::cout << s - s_R - s_I*Complex(0,1) << std::endl;

     //std::cout << Aic.transpose() - Aic.adjoint();
     //std::cout << "Ext:\n";
     //std::cout << s_R << "\n";
     //std::cout << "Int:\n";

     std::unique_ptr<LB_Solver_Complex> solver =
        std::make_unique<LB_Solver_Complex>(Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_I, J_R, J_I, V_R, V_I, *pre, 0.002);

     std::unique_ptr<LB_Solver_Complex2> solver2 =
        std::make_unique<LB_Solver_Complex2>(&Aii, &Aic, &Acc, J, V, *pre2, 0.002);

     Eigen::VectorXd X_R, X_I;
     solver->getX(X_R, X_I);
     //std::cout << X_R << "\n";

     //std::cout << solver2->getX() << "\n";

     std::cout << "Iterations:\n";
     for(int i=0; i<100;i++) {
         solver->do_iteration();
         solver2->do_iteration();
         std::cout << i << ":" << solver->getMinErrorl2Estimate() << ":" << solver->getMaxErrorl2Estimate() << std::endl;
         std::cout << i << ":" << solver2->getMinErrorl2Estimate() << ":" << solver2->getMaxErrorl2Estimate() << std::endl;
     }

     return 0;
 }
