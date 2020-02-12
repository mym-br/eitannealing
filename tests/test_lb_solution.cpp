#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"

#include <sys/time.h>
#include <sys/resource.h>

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
     std::cout << "Parsing mesh file..." << std::flush;
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
     input->initProblem(argv[1]);
     input->buildNodeCoefficients();
     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;
     std::cout << "Preparing dummy solution..." << std::flush;
     double *sol = new double[input->getNumCoefficients()];
     for(int i=0;i<input->getNumCoefficients();i++) sol[i]=1.0;


     matrix *Aii, *Aic, *Acc;
     std::cout << "Assembling matrices..." << std::flush;
     assembleProblemMatrix_lb(sol, &Aii, &Aic, &Acc, *input);
     std::cout << "Done\n" << std::flush;

     Eigen::VectorXd V(32), J(32);

     V.fill(0.0);
     V(4) = -30.738288239574686;
     V(0) = 30.733620796726793;
     V(5) = -12.495870755237416;
     V(31) = 12.49545757830731;
     V(1) = 8.627724263627385;
     V(3) = -8.624816166327843;

     J.fill(0.0);
     J(0) = 1.0;
     J(4) = -1.0;


     LB_Solver::Preconditioner *pre;
     pre = LB_Solver::makePreconditioner(*Aii, *Aic);

     LB_Solver *solver;
     solver = new LB_Solver(Aii, Aic, Acc, J, V, *pre, 0.002);

     std::cout << "Initial values: " << solver->getMinErrorl2Estimate() << " (min) " << solver->getMinErrorl2Estimate() << " (max)\n";
     for(int i=0; i<50;i++) {         
         solver->do_iteration();
         std::cout << i << ": " << solver->getMinErrorl2Estimate() << " (min) " << solver->getMinErrorl2Estimate() << " (max)\n";
     }
     solver->do_iteration();

     return 0;
 }
