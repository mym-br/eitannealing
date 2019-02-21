#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problemdescription.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"
#include "../src/assembleMatrix.h"

int main(int argc, char *argv[])
 {
  //  struct timespec time;
  //  clock_gettime(CLOCK_REALTIME, &time);
  // init_genrand64(time.tv_nsec);
  initProblem(argv[1]);

	 buildNodeCoefficients();
	 prepareSkeletonMatrix();
   createCoef2KMatrix();
     double *sol = new double[numcoefficients];
     for(int i=0;i<numcoefficients;i++) sol[i]=1.0;


     return 0;
 }
