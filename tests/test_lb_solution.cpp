#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/intcoef.h"
#include "../src/solution_lb_real.h"
#include "../src/gradientnormregularisation.h"

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
     std::shared_ptr<gradientNormRegularisation> reg(new gradientNormRegularisation(input));
     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;
     
     std::cout << "Parsing observations..." << std::flush;
     readingsScalar = new observations<double>;
     readingsScalar->initObs((const char **)&argv[2], argv[3], input->getNodesCount(), input->getGenericElectrodesCount());
     std::cout << "Done\n" << std::flush;
     
     std::cout << "Preparing dummy solution..." << std::flush;
     std::vector<double> sol;
     for(int i=0;i<input->getNumCoefficients();i++) sol.push_back((mincond+maxcond)/2);
     std::cout << "Done\n" << std::flush;
     init_genrand64(1);
     std::unique_ptr<solution_lb> current, next;
     shuffleData sdata;
     std::unique_ptr<shuffler> sh;
     std::cout << "Evaluating initial solution..." << std::flush;     
     current.reset(new solution_lb(input, *readingsScalar, reg, std::vector<double>(sol)));
     std::cout << "Initial values: " << current->getDMin() << " (min) " << current->getDMax() << " (max)\n";
     std::cout << "Saturating initial solution..." << std::flush;
     current->saturate();
     std::cout << "Done\n" << std::flush;
     std::cout << "Saturated values: " << current->getDMin() << " (min) " << current->getDMax() << " (max)\n";
     
     std::cout << "generating next solution..." << std::flush;     
     sh.reset(new shuffler(input, readingsScalar));
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
