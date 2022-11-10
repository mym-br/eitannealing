#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/intcoef.h"
#include "../src/solution_lb_complex2.h"
#include "../src/gradientnormregularisation.h"

int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <5) {
    	 std::cerr << "need mesh filename, currents (real), currents (imag) and tensions\n";
    	 return 1;
     }
     std::cout << "Parsing mesh file..." << std::flush;
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));

     input->initProblem(argv[1]);
     input->buildNodeCoefficients();
     input->prepareSkeletonMatrix();
     input->createCoef2KMatrix();

     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;

     std::cout << "Parsing observations..." << std::flush;
     std::unique_ptr<observations<std::complex<double> > > readings;
     readings = std::make_unique<observations<std::complex<double> > >();
     readings->initObs((const char **)&argv[2], argv[4], input->getNodesCount(), input->getGenericElectrodesCount());
     std::cout << "Done\n" << std::flush;


     std::cout << "Preparing dummy solution..." << std::flush;
     std::vector<double> sol_R(input->getNumCoefficients());
     for(int i=0;i<input->getNumCoefficients();i++) sol_R[i]=(mincond+maxcond)/2;
     std::vector<double> sol_I(input->getNumCoefficients());
     for(int i=0;i<32;i++) sol_I[i]=0.0;
     for(int i=32;i<input->getNumCoefficients();i++) sol_I[i]=0.05;

     std::vector<eigen_complexdouble_engine::scalar> sol;
     for(int i=0;i<input->getNumCoefficients();i++) {
        sol.emplace_back(sol_R[i], sol_I[i]);
     }
     std::cout << "Done.\n" << std::flush;

     init_genrand64(1);
     std::unique_ptr<solution_lb_complex2> current, next;
     shuffleData sdata;
     std::unique_ptr<shuffler> sh(new shuffler(input, readings.get()));
     std::shared_ptr<complexGradientNormRegularisation> reg(new complexGradientNormRegularisation(input));

     std::cout << "Evaluating initial solution..." << std::flush;
     current.reset(new solution_lb_complex2(input, *readings, reg, std::vector<std::complex<double> >(sol)));
     std::cout << "Initial values: " << current->getDMin() << " (min) " << current->getDMax() << " (max)\n";
     std::cout << "Saturating initial solution..." << std::flush;
     current->saturate();
     std::cout << "Done\n" << std::flush;
     std::cout << "Saturated values: " << current->getDMin() << " (min) " << current->getDMax() << " (max)\n";

     std::cout << "generating next solution..." << std::flush;
     next.reset(current->shuffle(&sdata, *sh));
     std::cout << "Done\n" << std::flush;
     std::cout << "Next solution Initial values: " << next->getDMin() << " (min) " << next->getDMax() << " (max)\n";
     std::cout << "improving next solution..." << std::flush;
     next->improve();
     std::cout << "Done\n" << std::flush;
     std::cout << "Next solution values: " << next->getDMin() << " (min) " << next->getDMax() << " (max)\n";
     std::cout << "improving next solution..." << std::flush;
     next->improve();
     std::cout << "Done\n" << std::flush;
     std::cout << "Next solution values: " << next->getDMin() << " (min) " << next->getDMax() << " (max)\n";
     std::cout << "improving next solution..." << std::flush;
     next->improve();
     std::cout << "Done\n" << std::flush;
     std::cout << "Next solution values: " << next->getDMin() << " (min) " << next->getDMax() << " (max)\n";

     return 0;
 }
