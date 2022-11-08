/*
 * solution_lb_real.h
 *
 *  Created on: Nov 07, 2022
 *      Author: thiago
 */

#ifndef SOLUTION_LB_REAL_H_
#define SOLUTION_LB_REAL_H_

#include "solution.h"
#include "solver_lb.h"
#include "solution_lb.h"
#include "gradientnormregularisation.h"

typedef observations<double> realobservations;
typedef LBMatrixBuilder<eigen_double_engine::scalar, eigen_double_engine::symMatrix, eigen_double_engine::matrix> realMatrixBuilder;
typedef solution_lb_gen<LB_Solver, double, realobservations, gradientNormRegularisation, realMatrixBuilder, shuffleData, shuffler> solution_lb;

#endif  // SOLUTION_LB_COMPLEX2_H_
