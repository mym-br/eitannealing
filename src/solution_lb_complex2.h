/*
 * solution_lb_complex2.h
 *
 *  Created on: Nov 07, 2022
 *      Author: thiago
 */

#ifndef SOLUTION_LB_COMPLEX2_H_
#define SOLUTION_LB_COMPLEX2_H_

#include "solutionbase.h"
#include "solver_lb.h"
#include "solution_lb.h"
#include "solver_lb_complex2.h"
#include "problem.h"

class complexGradientNormRegularisation
{
private:
    int electrodecoefficients;
    std::unique_ptr<matrix> regularizationMatrix;
    void buildMatrix();
    int coefficientMap(int node);
    std::set<int> lastElectrodeNodes;
	std::shared_ptr<problem> input;
public:
  	complexGradientNormRegularisation(std::shared_ptr<problem> _input);
    double getRegularisation(const std::complex<double> *sol) const;
};

typedef observations<std::complex<double> > complexobservations;
typedef LBMatrixBuilder<eigen_complexdouble_engine::scalar, eigen_complexdouble_engine::symMatrix, eigen_complexdouble_engine::matrix> complexMatrixBuilder;
typedef solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler> solution_lb_complex2;

#endif  // SOLUTION_LB_COMPLEX2_H_
