#include "solution_lb_impl.h"
#include "observations.h"
#include "gradientnormregularisation.h"
#include "solution_lb_real.h"

template class solution_lb_gen<LB_Solver, double, realobservations, gradientNormRegularisation, realMatrixBuilder, shuffleData, shuffler>;
