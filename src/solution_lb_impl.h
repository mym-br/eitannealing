/*
 * solution.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#include "solution_lb.h"
#include "random.h"
#include <iostream>
#include "util/standard_deviation.hpp"
#include "util/heap_siftdown.hpp"

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::improve()
{
	unsigned critical = solver_heap_ordered_by_errors[0];
	// Do another iteration on the critical solver
	simulations[critical].do_iteration();
	this->totalit++;
	// Recalcule expected distance and boundaries
	double d = simulations[critical].getErrorl2Estimate();
	double min_d = simulations[critical].getMinErrorl2Estimate();
	double max_d = simulations[critical].getMaxErrorl2Estimate();

	distance2[critical] = d*d;
	maxdist2[critical] = max_d*max_d;
	mindist2[critical] = min_d*min_d;
	err[critical] = max_d-min_d;
	err_x_dist[critical] = max_d*err[critical];
	totalDist = std::sqrt(distance2.sum())+regularisation;
	minTotalDist = std::sqrt(mindist2.sum())+regularisation;
	maxTotalDist = std::sqrt(maxdist2.sum())+regularisation;
	// reevaluate critical
	heap_sift_top_down(solver_heap_ordered_by_errors.begin(), solver_heap_ordered_by_errors.begin(), solver_heap_ordered_by_errors.end(),
		   [&errors = this->err_x_dist](unsigned i, unsigned j) { return errors[i]<errors[j]; });
	critErr = this->err_x_dist[solver_heap_ordered_by_errors[0]];
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
bool solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::
	compareWith(solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler>  &target, float kt, float prob)
{
	double delta, expdelta;
	// Ensure errors are within required margin
	while(true) {
		double min_delta = target.minTotalDist - this->maxTotalDist;
		double max_delta = target.maxTotalDist - this->minTotalDist;
		delta = target.totalDist - this->totalDist;
		expdelta = exp(-delta/kt);
		// Check boundary conditions:
		// Upper bound is negative, no doubt here
		if(max_delta < 0) break;
		// Estimate is negative, but upper bound is positive
		else if(delta <= 0) {
			if(exp(-max_delta/kt)>=(1-prob)) break; // Upper condition only
		}
		// Estimate and upper bounds are positive, lower bound is negative
		else if(min_delta <= 0) {
			if(expdelta >= 1 - prob) { // Lower condition
				if(expdelta <= prob) break;	// upper condition
				if(exp(-max_delta/kt)>=(expdelta-prob)) break;
			}
		}
		// Classic case, everything is positive
		else {
			if(exp(-min_delta/kt)<=prob+expdelta) { // lower condition
				if(expdelta <= prob) break;	// upper condition
				if(exp(-max_delta/kt)>=(expdelta-prob)) break;
			}
		}
		// Not there yet, improve boundaries
		// Select wich one to improve
		if(this->critErr > target.critErr)
			this->improve();
		else
			target.improve();
	}
	if(delta <= 0) {
		//std::cout << "+";
		return true;
	}
	//std::cout << "-" << expdelta << std::endl;


	// Now randomly accept based on the energy
	if(genreal()<expdelta) return true;
	return false;
}

// Construct from solution vector copy
template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler>::
	solution_lb_gen(std::shared_ptr<problem> p, const observations &o, std::shared_ptr<regularizer> r, std::vector<admittance> &&sol):
				p(p), matrixBuilder(new matBuilder(*p)), reg(r), o(o),
				sol(std::move(sol)),
				distance2(o.getNObs()),
				maxdist2(o.getNObs()),
				mindist2(o.getNObs()),
				err(o.getNObs()),
				err_x_dist(o.getNObs()),
				solver_heap_ordered_by_errors(o.getNObs())
{

    this->initMatrices();
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::
	solution_lb_gen(std::vector<admittance> &&sol, const solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler>  &base):
                p(base.p), matrixBuilder(base.matrixBuilder), reg(base.reg), o(base.o),
                sol(std::move(sol)),
                distance2(o.getNObs()),
                maxdist2(o.getNObs()),
                mindist2(o.getNObs()),
                err(o.getNObs()),
                err_x_dist(o.getNObs()),
				solver_heap_ordered_by_errors(o.getNObs())
{
        this->initMatrices();
        this->initSimulations(base);
        this->initErrors();
}


// New random solution
template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::
	solution_lb_gen(std::shared_ptr<problem> p, const observations &o, std::shared_ptr<regularizer> r):
                p(p), matrixBuilder(new matBuilder(*p)), reg(r), o(o),
		sol(solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::getNewRandomSolution(p->getNumCoefficients())),
		distance2(o.getNObs()),
		maxdist2(o.getNObs()),
		mindist2(o.getNObs()),
		err(o.getNObs()),
		err_x_dist(o.getNObs()),
		solver_heap_ordered_by_errors(o.getNObs())
{
      this->initMatrices();
      this->initSimulations();
	this->initErrors();
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::initMatrices()
{
      Aii = matrixBuilder->buildAiiMatrix(sol);
	  Aic = matrixBuilder->buildAicMatrix(sol);
	  Acc = matrixBuilder->buildAccMatrix(sol);
      precond.reset(solver::makePreconditioner(*Aii, *Aic));
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::initSimulations()
{
    // Prepare solvers
	int i;
	this->totalit = 0;

    // Reserve so we don't need to move solvers around
    simulations.reserve(o.getNObs());
    // 1st solution estimates also least eigenvalue
    simulations.emplace_back(Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[0].tail(p->getGenericElectrodesCount()),
                    typename solver::vector(o.getTensions()[0]), *precond,  100, (float)0.0001, this->eigenvalue_estimate, this->eigenvector_estimate);
    this->totalit += simulations[0].getIteration();
	for(i=1;i<o.getNObs();i++)
	{
        simulations.emplace_back(
            Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[i].tail(p->getGenericElectrodesCount()),
            typename solver::vector(o.getTensions()[i]), *precond, this->eigenvalue_estimate);
		simulations[i].do_iteration();
		this->totalit += simulations[i].getIteration();
	}
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::
	initSimulations(const solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler>  &base)
{
        // Prepare solvers
        int i;
        this->totalit = 0;

        // Reserve so we don't need to move solvers around
        simulations.reserve(o.getNObs());
        // 1st solution estimates also least eigenvalue

        simulations.emplace_back(Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[0].tail(p->getGenericElectrodesCount()),
                        typename solver::vector(o.getTensions()[0]), *precond,
                        base.simulations[0].getX(), base.getLeastEigenvectorEstimation(), 100, (float)0.0001, this->eigenvalue_estimate, this->eigenvector_estimate);
        this->totalit += simulations[0].getIteration();
        for(i=1;i<o.getNObs();i++)
        {
                simulations.emplace_back(
                        Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[i].tail(p->getGenericElectrodesCount()),
                        typename solver::vector(o.getTensions()[i]), *precond, this->eigenvalue_estimate,
					       base.simulations[i].getX());
                simulations[i].do_iteration();
                this->totalit += simulations[i].getIteration();
        }
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::initErrors()
{
	int i;
	// Calc regularisation value
  double electrodevariance = population_variance(this->sol.begin(), this->sol.begin()+this->p->getGenericElectrodesCount());
	this->regularisation = reg->getRegularisation(this->sol.data())*0.010+electrodevariance*20;
	// Retrieve distance estimates, errors and boundaries
	for(i=0;i<o.getNObs();i++) {
		double d = simulations[i].getErrorl2Estimate();
		double min_d = simulations[i].getMinErrorl2Estimate();
		double max_d = simulations[i].getMaxErrorl2Estimate();
		distance2[i] = d*d;
		maxdist2[i] = max_d*max_d;
		mindist2[i] = min_d*min_d;
		err[i] = max_d - min_d;
		err_x_dist[i] = max_d*err[i];
		solver_heap_ordered_by_errors[i] = i;
	}
	totalDist = std::sqrt(distance2.sum())+regularisation;
	minTotalDist = std::sqrt(mindist2.sum())+regularisation;
	maxTotalDist = std::sqrt(maxdist2.sum())+regularisation;
	// evaluate critical
	make_heap_down(solver_heap_ordered_by_errors.begin(), solver_heap_ordered_by_errors.end(),
				   [&errors = this->err_x_dist](unsigned i, unsigned j) { return errors[i]<errors[j]; });
	critErr = this->err_x_dist[solver_heap_ordered_by_errors[0]];
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler>  *solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::
	shuffle(shuffleData *data, const shuffler &sh) const
{
	solution_lb_gen *res;
	try {
		res = new solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> (getShuffledSolution(data, sh), *this);
	} catch(...) {

		for(int i=0;i<65;i++)
			std::cout << i << ":" << res->sol[i] << std::endl;
		exit(0);
	}
	return res;
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::saturate()
{
      ensureMinIt(p->getNodesCount()+30);
}

template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler>
void solution_lb_gen<solver, admittance, observations, regularizer, matBuilder, shuffleData, shuffler> ::ensureMinIt(unsigned int it)
{
      static Eigen::VectorXd aux(p->getGenericElectrodesCount()-1);
      for(unsigned int i = 0; i<o.getNObs();i++) {
            while(simulations[i].getIteration()<it) {
                simulations[i].do_iteration();
                this->totalit++;

				double d = simulations[i].getErrorl2Estimate();
				double min_d = simulations[i].getMinErrorl2Estimate();
				double max_d = simulations[i].getMaxErrorl2Estimate();

                distance2[i] = d*d;
                maxdist2[i] = max_d*max_d;
                mindist2[i] = min_d*min_d;
                err[i] = max_d-min_d;
                err_x_dist[i] = max_d*err[i];
                totalDist = std::sqrt(distance2.sum())+regularisation;
                minTotalDist = std::sqrt(mindist2.sum())+regularisation;
                maxTotalDist = std::sqrt(maxdist2.sum())+regularisation;

                // reevaluate critical
			    make_heap_down(solver_heap_ordered_by_errors.begin(), solver_heap_ordered_by_errors.end(),
				   [&errors = this->err_x_dist](unsigned i, unsigned j) { return errors[i]<errors[j]; });
				critErr = this->err_x_dist[solver_heap_ordered_by_errors[0]];
            }
      }
}

#ifdef __GENERATE_LB_BENCHMARKS
solution_lb::solution_lb(const float *sigma, char benchmarktag):
				sol(solution_lb::copySolution(sigma)),
				distance(o.getNObs()),
				maxdist(o.getNObs()),
				mindist(o.getNObs()),
				err(o.getNObs()),
				err_x_dist(o.getNObs())
{
}

solution_lb_benchmarked::solution_lb_benchmarked(const float *sigma, benchmark_entry *bench, int n):solution_lb(sigma), vector(bench), n(n)
{
	vector[0].timestamp = getTimestamp();
	vector[0].e_low = -1;
	vector[0].e_high = -1;
	assembleProblemMatrix_lb(sol, &Aii, &Aic, &Acc, 32);
	vector[1].timestamp = getTimestamp();
	vector[1].e_low = -2;
	vector[1].e_high = -2;
        precond.reset(LB_Solver::makePreconditioner(*Aii));
	vector[2].timestamp = getTimestamp();
	vector[2].e_low = -3;
	vector[2].e_high = -3;
	this->initSimulations();
	this->initErrors();
	vector[3].timestamp = getTimestamp();
	vector[3].e_low = getDMin();
	vector[3].e_high = getDMax();
	i = 4;
}

void solution_lb_benchmarked::performBenchmark()
{
    while(i<n) {
      this->improve();
      vector[i].timestamp = getTimestamp();
      vector[i].e_low = getDMin();
      vector[i].e_high = getDMax();
      i++;
    }
}


// Linux-only
#include <sys/time.h>

int solution_lb_benchmarked::getTimestamp()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec*1000000 + time.tv_usec;
}
#endif // __GENERATE_LB_BENCHMARKS
