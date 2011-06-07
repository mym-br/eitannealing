/*
 * solution.cpp
 *
 *  Created on: Sep 12, 2010
 *      Author: thiago
 */

#include "solution.h"
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

void solution::improve()
{
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXd aux(electrodes.size()-1);

	// Do another iteration on the critical solver
	simulations[critical]->do_iteration();
	this->totalit++;
	// Recalcule expected distance and boundaries
	
	aux = simulations[critical]->getX().end(electrodes.size()-1);
	aux -= tensions[critical];
					
	distance[critical] = aux.norm();
	err[critical] = sqrt(simulations[critical]->getErrorl2Estimate());
	maxdist[critical] = distance[critical] + err[critical];
	mindist[critical] = max(distance[critical] - err[critical],0);
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

bool solution::compareWith(solution &target, float kt, float prob)
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


solution::solution(const float *sigma):
				sol(solution::copySolution(sigma)),
				stiffness(solution::getNewStiffness(sol)),
				precond(new SparseIncompleteLLT(*stiffness)),
				simulations(new CG_Solver *[nobs]),
				distance(nobs),
				maxdist(nobs),
				mindist(nobs),
				err(nobs),
				err_x_dist(nobs)
{
	this->initSimulations();
	this->initErrors();
}


// New random solution
solution::solution():
		sol(solution::getNewRandomSolution()),
		stiffness(solution::getNewStiffness(sol)),
		precond(new SparseIncompleteLLT(*stiffness)),
		simulations(new CG_Solver *[nobs]),
		distance(nobs),
		maxdist(nobs),
		mindist(nobs),
		err(nobs),
		err_x_dist(nobs)
{
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solution::solution(float *sigma, const solution &base):
		sol(sigma),
		stiffness(solution::getNewStiffness(sol)),
		precond(new SparseIncompleteLLT(*stiffness)),
		simulations(new CG_Solver *[nobs]),
		distance(nobs),
		maxdist(nobs),
		mindist(nobs),
		err(nobs),
		err_x_dist(nobs)
{
	this->initSimulations(base);
	this->initErrors();
}

void solution::initSimulations(const solution &base)
{
	// Prepare solvers
	int i;
	this->totalit = 0;
	for(i=0;i<nobs;i++)
	{
		// Reuse previous solutions as initial values
		simulations[i] = new CG_Solver(*stiffness, currents[i], base.simulations[i]->getX(), *precond);
		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		simulations[i]->do_iteration();
		simulations[i]->do_iteration();
		simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while(ndecr<3) {
			simulations[i]->do_iteration();
			aux = simulations[i]->getErrorl2Estimate();
			if(aux>=err) ndecr = 0;
			else {
				ndecr++;
			}
			err = aux;
		}
		this->totalit += simulations[i]->getIteration();
	}
}

void solution::initSimulations()
{
	// Prepare solvers
	int i;
	this->totalit = 0;
	for(i=0;i<nobs;i++)
	{
		simulations[i] = new CG_Solver(*stiffness, currents[i], *precond);
		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		simulations[i]->do_iteration();
		simulations[i]->do_iteration();
		simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while(ndecr<3) {
			simulations[i]->do_iteration();
			aux = simulations[i]->getErrorl2Estimate();
			if(aux>=err) ndecr = 0;
			else {
				ndecr++;
			}
			err = aux;
		}
		this->totalit += simulations[i]->getIteration();
	}
}

void solution::initErrors()
{
	int i;
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXd aux(electrodes.size()-1);
	// Retrieve distance estimates, errors and boundaries
	for(i=0;i<nobs;i++) {
		// Compare with observation
		aux = simulations[i]->getX().end(aux.size());
		aux -= tensions[i];
		
		distance[i] = aux.norm();
		err[i] = sqrt(simulations[i]->getErrorl2Estimate());
		maxdist[i] = distance[i] + err[i];
		mindist[i] = max(distance[i] - err[i],0);
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


float *solution::copySolution(const float *sol)
{
	float *res = new float[numcoefficients];

	for(int i=0;i<numcoefficients;i++)
		res[i] = sol[i];

	return res;
}

float *solution::getNewRandomSolution()
{
	float *res = new float[numcoefficients];
/*
	res[0] = 0.118859;
	res[1] = 0.0911114;
	res[2] = 0.100251;
	res[3] = 0.106117;
	res[4] = 0.160247;
	res[5] = 0.105798;
	res[6] = 0.106278;
	res[7] = 0.104063;
	res[8] = 0.0882884;
	res[9] = 0.0978375;
	res[10] = 0.0876754;
	res[11] = 0.0793514;
	res[12] = 0.0787333;
	res[13] = 0.0845089;
	res[14] = 0.0753199;
	res[15] = 0.0786953;
	res[16] = 0.0791697;
	res[17] = 0.0832361;
	res[18] = 0.0828951;
	res[19] = 0.0927432;
	res[20] = 0.101078;
	res[21] = 0.100138;
	res[22] = 0.102309;
	res[23] = 0.112525;
	res[24] = 0.125682;
	res[25] = 0.0879498;
	res[26] = 0.121403;
	res[27] = 0.11835;
	res[28] = 0.120835;
	res[29] = 0.0954887;
	res[30] = 0.132771;
	res[31] = 0.118369;
	res[32] = 1/2.82314;
*/	
	for(int i=0;i<numcoefficients;i++)
		res[i] = mincond+genreal()*(maxcond-mincond);

	return res;
}

float *solution::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	float *res = solution::copySolution(sol);
	// head or tails
	//if(genint(2)) { // Normal
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
	/*} else { // swap
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
	}*/
	return res;
}

void shuffler::addShufflerFeedback(const shuffleData &data, bool pos)
{
	if(pos) { // positive feedback
		if(data.swap)
			this->swapshuffleconsts[data.ncoef] /= 2;
		else
			this->shuffleConsts[data.ncoef] /= 2;
	} else {
		if(data.swap)
			this->swapshuffleconsts[data.ncoef]++;
		else
			this->shuffleConsts[data.ncoef]++;
	}
}

solution *solution::shuffle(shuffleData *data, const shuffler &sh) const
{
	float *sigma = getShuffledSolution(data, sh);
	solution *res;
	try {
		res = new solution(sigma, *this);
	} catch(...) {
		for(int i=0;i<65;i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
}

solution::~solution()
{
	delete[] sol;
	delete stiffness;
	delete precond;
	for(int i=0;i<nobs;i++) {
		delete simulations[i];
	}
	delete[] simulations;
}

