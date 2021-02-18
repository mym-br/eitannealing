#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
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
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
     input->buildNodeCoefficients();

   
     double *sol_R = new double[input->getNumCoefficients()];
     for(int i=0;i<input->getNumCoefficients();i++) sol_R[i]=1.0;
     
     
     double *sol_I = new double[input->getNumCoefficients()];
     for(int i=0;i<32;i++) sol_I[i]=0.0;
     for(int i=32;i<input->getNumCoefficients();i++) sol_I[i]=1.0;
     
     matrix *Aii_R, *Aii_I, *Aic_R, *Aic_I, *Acc_R, *Acc_I;
     
     assembleProblemMatrix_lb(sol_R, &Aii_R, &Aic_R, &Acc_R, *input);
     assembleProblemMatrix_lb(sol_I, &Aii_I, &Aic_I, &Acc_I, *input);
             
  
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
     
     


     return 0;
 }
