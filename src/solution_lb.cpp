/*
 * solution.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#include "solution_lb.h"
#include "random.h"
#include "observations.h"
#include "problemdescription.h"
#include <iostream>
#include <boost/numeric/interval.hpp>

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

void solution_lb::improve()
{
	// Do another iteration on the critical solver
	simulations[critical]->do_iteration();
	this->totalit++;
	// Recalcule expected distance and boundaries
	
	distance[critical] = simulations[critical]->getErrorl2Estimate();
	maxdist[critical] = simulations[critical]->getMaxErrorl2Estimate();
	mindist[critical] = simulations[critical]->getMinErrorl2Estimate();
	err[critical] = maxdist[critical]-mindist[critical];
	err_x_dist[critical] = maxdist[critical]*err[critical];
	totalDist = distance.norm();
	minTotalDist = mindist.norm();
	maxTotalDist = maxdist.norm();
	// reevaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for(int i = 1; i<nobs;i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}

bool solution_lb::compareWith(solution_lb &target, float kt, float prob)
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


solution_lb::solution_lb(const float *sigma):
				sol(solution_lb::copySolution(sigma)),
				simulations(new LB_Solver *[nobs]),
				distance(nobs),
				maxdist(nobs),
				mindist(nobs),
				err(nobs),
				err_x_dist(nobs)
{
	assembleProblemMatrix_lb(sol, &Aii, &Aic, &Acc, 32);
    SparseIncompleteLLT precond(*Aii);
					
	this->initSimulations();
	this->initErrors();
}


// New random solution
solution_lb::solution_lb():
		sol(solution_lb::getNewRandomSolution()),
		simulations(new LB_Solver *[nobs]),
		distance(nobs),
		maxdist(nobs),
		mindist(nobs),
		err(nobs),
		err_x_dist(nobs)
{
	assembleProblemMatrix_lb(sol, &Aii, &Aic, &Acc, 32);
    SparseIncompleteLLT precond(*Aii);
	this->initSimulations();
	this->initErrors();
}

void solution_lb::initSimulations()
{
	// FIXME: Get a more acurate estimative of lowest eigenvalue!
	// Get a geometric average of coefficients
	double s = 1;
	for(int i=0; i< numcoefficients;i++)
		s*=sol[i];
	s = pow(s, 1/numcoefficients);
	double a = 0.01*s;
	// Prepare solvers
	int i;
	this->totalit = 0;
	for(i=0;i<nobs;i++)
	{
		simulations[i] = new LB_Solver(
			Aii, Aic, Acc, Eigen::VectorXd(currents[i].end(31)), 
			Eigen::VectorXd(tensions[i].end(31)), *precond, a);
		simulations[i]->do_iteration();
		this->totalit += simulations[i]->getIteration();
	}
}

void solution_lb::initErrors()
{
	int i;
	// Retrieve distance estimates, errors and boundaries
	for(i=0;i<nobs;i++) {
		// Compare with observation
		
		distance[i] = simulations[i]->getErrorl2Estimate();
		maxdist[i] = simulations[i]->getMaxErrorl2Estimate();
		mindist[i] = simulations[i]->getMinErrorl2Estimate();
		err[i] = maxdist[i]-mindist[i];
		err_x_dist[i] = maxdist[i]*err[i];
	}
	totalDist = distance.norm();
	minTotalDist = mindist.norm();
	maxTotalDist = maxdist.norm();
	// evaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for(i = 1; i<nobs;i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}


float *solution_lb::copySolution(const float *sol)
{
	float *res = new float[numcoefficients];

	for(int i=0;i<numcoefficients;i++)
		res[i] = sol[i];

	return res;
}

float *solution_lb::getNewRandomSolution()
{
	float *res = new float[numcoefficients];
	int i = 0;
/*
	res[i++] = 0.0795333;
	res[i++] = 0.154207;
	res[i++] = 0.10827;
	res[i++] = 0.107503;
	res[i++] = 0.120324;
	res[i++] = 0.115978;
	res[i++] = 0.112217;
	res[i++] = 0.109881;
	res[i++] = 0.103229;
	res[i++] = 0.0989397;
	res[i++] = 0.0964289;
	res[i++] = 0.0905536;
	res[i++] = 0.0856748;
	res[i++] = 0.0871923;
	res[i++] = 0.0870397;
	res[i++] = 0.0801445;
	res[i++] = 0.0874562;
	res[i++] = 0.0893944;
	res[i++] = 0.0892118;
	res[i++] = 0.0922962;
	res[i++] = 0.109619;
	res[i++] = 0.115378;
	res[i++] = 0.120788;
	res[i++] = 0.110558;
	res[i++] = 0.115594;
	res[i++] = 0.122183;
	res[i++] = 0.11994;
	res[i++] = 0.12327;
	res[i++] = 0.123105;
	res[i++] = 0.114889;
	res[i++] = 0.122205;
	res[i++] = 0.119641;
	res[i++] = 0.35731;*/

	for(i=0;i<numcoefficients;i++)
		res[i] = mincond+genreal()*(maxcond-mincond);

	return res;
}

float *solution_lb::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	float *res = solution_lb::copySolution(sol);
	// head or tails
	if(genint(2)) { // Normal
		int ncoef = genint(numcoefficients);	// Lower values fixed;

		if(sh.shuffleConsts[ncoef]==0) {
			res[ncoef] = mincond+genreal()*(maxcond-mincond);
		} else {
			float val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for(int i=0;i<sh.shuffleConsts[ncoef];i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoef];
				rnd -= 0.5;
				rnd *= (maxcond - mincond);
				val += rnd;
			} while((val < mincond) || (val > maxcond));
			res[ncoef] = val;
		}
		if(data) {
			data->swap = false;
			data->ncoef = ncoef;
		}
	} else { // swap
		int ncoef = genint(innerAdjacency.size());
		int node1, node2;

		node1 = node2coefficient[innerAdjacency[ncoef].first];
		node2 = node2coefficient[innerAdjacency[ncoef].second];
		
		// Order nodes
		if(res[node1]>res[node2]) {
			int aux = node1;
			node1 = node2;;
			node2 = aux;
		}
		float v1 = res[node1], v2 = res[node2];
		float a = max( min(v1-mincond, maxcond-v2), min(maxcond-v1, v2-mincond));

		float delta;
		do {
			if(sh.swapshuffleconsts[ncoef]==0) {
				delta = a*(genreal()*2 - 1);
			} else {
				double rnd = 0;
				for(int i=0;i<sh.swapshuffleconsts[ncoef];i++)
					rnd += genreal();
				rnd /= sh.swapshuffleconsts[ncoef];
				rnd -= 0.5;
				delta = a*rnd;
			}
			v1 = res[node1] - delta;
			v2 = res[node2] + delta;
		} while((v1 < mincond) || (v2 < mincond) || (v1 > maxcond) || (v2 > maxcond));
		res[node1] = v1;
		res[node2] = v2;
		if(data) {
			data->swap = true;
			data->ncoef = ncoef;
		}
	}
	return res;
}

solution_lb *solution_lb::shuffle(shuffleData *data, const shuffler &sh) const
{
	float *sigma = getShuffledSolution(data, sh);
	solution_lb *res;
	try {
		res = new solution_lb(sigma);
	} catch(...) {
		for(int i=0;i<65;i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
}


solution_lb::~solution_lb()
{
	delete[] sol;
	delete Aii;
	delete Aic;
	delete Acc;
	delete precond;
	for(int i=0;i<nobs;i++) {
		delete simulations[i];
	}
	delete[] simulations;
}

