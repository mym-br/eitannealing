/*
 * solution.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#include "solution_lb_complex.h"
#include "random.h"
#include "observations_complex.h"
//#include "problemdescription.h"
#include <iostream>
//#include <boost/numeric/interval.hpp>
#include "gradientnormregularisation.h"

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

float mincond_I = 0;
float maxcond_I = 0.1;

void solution_lb_complex::improve()
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
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
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

bool solution_lb_complex::compareWith(solution_lb_complex &target, float kt, float prob)
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


solution_lb_complex::solution_lb_complex(const double *s_R, const double *s_I):
				solRe(solution_lb_complex::copySolution(s_R)),
        solIm(solution_lb_complex::copySolution(s_I)),
				simulations(new LB_Solver_Complex *[nobs]),
				distance(nobs),
				maxdist(nobs),
				mindist(nobs),
				err(nobs),
				err_x_dist(nobs)
{
        assembleProblemMatrix_lb(solRe, &Aii_R, &Aic_R, &Acc_R, 32);
        assembleProblemMatrix_lb(solIm, &Aii_I, &Aic_I, &Acc_I, 32);
        precond.reset(LB_Solver_Complex::makePreconditioner(*Aii_R, *Aii_I, *Aic_R, *Aic_I));

	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solution_lb_complex::solution_lb_complex(double *s_R, double *s_I, const solution_lb_complex &base):
                solRe(s_R),
                solIm(s_I),
                simulations(new LB_Solver_Complex *[nobs]),
                distance(nobs),
                maxdist(nobs),
                mindist(nobs),
                err(nobs),
                err_x_dist(nobs)
{
      assembleProblemMatrix_lb(solRe, &Aii_R, &Aic_R, &Acc_R, 32);
      assembleProblemMatrix_lb(solIm, &Aii_I, &Aic_I, &Acc_I, 32);
      precond.reset(LB_Solver_Complex::makePreconditioner(*Aii_R, *Aii_I, *Aic_R, *Aic_I));

        this->initSimulations(base);
        this->initErrors();
}


// New random solution
solution_lb_complex::solution_lb_complex():
		solRe(solution_lb_complex::getNewRandomSolution_R()),
    solIm(solution_lb_complex::getNewRandomSolution_I()),
		simulations(new LB_Solver_Complex *[nobs]),
		distance(nobs),
		maxdist(nobs),
		mindist(nobs),
		err(nobs),
		err_x_dist(nobs)
{
  assembleProblemMatrix_lb(solRe, &Aii_R, &Aic_R, &Acc_R, 32);
  assembleProblemMatrix_lb(solIm, &Aii_I, &Aic_I, &Acc_I, 32);
  precond.reset(LB_Solver_Complex::makePreconditioner(*Aii_R, *Aii_I, *Aic_R, *Aic_I));
  this->initSimulations();
	this->initErrors();
}

void solution_lb_complex::initSimulations()
{

	// Prepare solvers
	int i;
	this->totalit = 0;
  Eigen::VectorXd J_R, J_I, V_R, V_I;
  J_R = currents[0].tail(31);
  J_I = currents_I[0].tail(31);
  J_R.conservativeResize(32);
  J_R(31) = -J_R.sum();
  J_I.conservativeResize(32);
  J_I(31) = -J_I.sum();
  V_R = tensions[0].tail(31);
  V_R.conservativeResize(32);
  V_R(31) = 0;
  V_I = tensions_I[0].tail(31);
  V_I.conservativeResize(32);
  V_I(31) = 0;

        // 1st solution estimates also least eigenvalue
        LB_Solver_EG_Complex_Estimate *solver = new LB_Solver_EG_Complex_Estimate(
                        Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_R, J_R, J_I, V_R, V_I, *precond,  100, 0.0001);
        double a = solver->getLeastEvEst();
        simulations[0] = solver;
	this->totalit += solver->getIteration();
	for(i=1;i<nobs;i++)
	{
    J_R = currents[i].tail(31);
    J_I = currents_I[i].tail(31);
    J_R.conservativeResize(32);
    J_R(31) = -J_R.sum();
    J_I.conservativeResize(32);
    J_I(31) = -J_I.sum();
    V_R = tensions[i].tail(31);
    V_R.conservativeResize(32);
    V_R(31) = 0;
    V_I = tensions_I[i].tail(31);
    V_I.conservativeResize(32);
    V_I(31) = 0;
                simulations[i] = new LB_Solver_Complex(
                        Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_R, J_R, J_I, V_R, V_I, *precond, a);
		simulations[i]->do_iteration();
		this->totalit += simulations[i]->getIteration();
	}
}

void solution_lb_complex::initSimulations(const solution_lb_complex &base)
{
        // Prepare solvers
        int i;
        this->totalit = 0;
        Eigen::VectorXd J_R, J_I, V_R, V_I;
        J_R = currents[0].tail(31);
        J_I = currents_I[0].tail(31);
        J_R.conservativeResize(32);
        J_R(31) = -J_R.sum();
        J_I.conservativeResize(32);
        J_I(31) = -J_I.sum();
        V_R = tensions[0].tail(31);
        V_R.conservativeResize(32);
        V_R(31) = 0;
        V_I = tensions[0].tail(31);
        V_I.conservativeResize(32);
        V_I(31) = 0;
        const LB_Solver_EG_Complex_Estimate *baseEVSolver = dynamic_cast<const LB_Solver_EG_Complex_Estimate *>(base.simulations[0]);
        // 1st solution estimates also least eigenvalue
        Eigen::VectorXd X_R, X_I;
        baseEVSolver->getX(X_R, X_I);
        LB_Solver_EG_Complex_Estimate *solver = new LB_Solver_EG_Complex_Estimate(
                        Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_R, J_R, J_I, V_R, V_I, *precond,
                        X_R, X_I, baseEVSolver->getEvector(), 100, 0.0001);
        double a = solver->getLeastEvEst();
        simulations[0] = solver;
	this->totalit += solver->getIteration();
        for(i=1;i<nobs;i++)
        {
                J_R = currents[i].tail(31);
                J_I = currents_I[i].tail(31);
                J_R.conservativeResize(32);
                J_R(31) = -J_R.sum();;
                J_I.conservativeResize(32);
                J_I(31) = -J_I.sum();
                V_R = tensions[0].tail(31);
                V_I = tensions_I[i].tail(31);
                V_R.conservativeResize(32);
                V_R(31) = 0;
                V_I.conservativeResize(32);
                V_I(31) = 0;
                base.simulations[i]->getX(X_R, X_I);
                simulations[i] = new LB_Solver_Complex(
                        Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_R, J_R, J_I, V_R, V_I, *precond, a,
					       X_R, X_I);

                simulations[i]->do_iteration();
                this->totalit += simulations[i]->getIteration();
        }
}


void solution_lb_complex::initErrors()
{
	int i;
	// Calc regularisation value
	this->regularisation = gradientNormRegularisation::getInstance()->getRegularisation(this->solRe)*0.0005+gradientNormRegularisation::getInstance()->getRegularisation(this->solIm)*0.0005;
	// Retrieve distance estimates, errors and boundaries
	for(i=0;i<nobs;i++) {
		// Compare with observation

		distance[i] = simulations[i]->getErrorl2Estimate();
		maxdist[i] = simulations[i]->getMaxErrorl2Estimate();
		mindist[i] = simulations[i]->getMinErrorl2Estimate();
		err[i] = maxdist[i]-mindist[i];
		err_x_dist[i] = maxdist[i]*err[i];
	}
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
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


double *solution_lb_complex::copySolution(const double *sol) const
{
	double *res = new double[numcoefficients];

	for(int i=0;i<numcoefficients;i++)
		res[i] = sol[i];

	return res;
}

double *solution_lb_complex::getNewRandomSolution_R()
{
	double *res = new double[numcoefficients];
	int i = 0;
	for(;i<numcoefficients;i++)
		res[i] = mincond+genreal()*(maxcond-mincond);

	return res;
}

double *solution_lb_complex::getNewRandomSolution_I()
{
	double *res = new double[numcoefficients];
	int i = 0;
	for(;i<numcoefficients;i++)
		res[i] = mincond_I+genreal()*(maxcond_I-mincond_I);
	return res;
}

double *solution_lb_complex::getShuffledSolution_I(shuffleData *data, const shuffler &sh) const
{
	double *res = solution_lb_complex::copySolution(solIm);
	// head or tails
	if(genint(2)) { // Normal
		int ncoef = genint(numcoefficients);
		float minc, maxc;
                  maxc = maxcond_I;
                  minc = mincond_I;

		if(sh.shuffleConsts[ncoef]==0) {
			res[ncoef] = minc+genreal()*(maxc-minc);
		} else {
			float val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for(int i=0;i<sh.shuffleConsts[ncoef];i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoef];
				rnd -= 0.5;
				rnd *= (maxc - minc);
				val += rnd;
			} while((val < minc) || (val > maxc));
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
		float a = max( min(v1-mincond_I, maxcond_I-v2), min(maxcond_I-v1, v2-mincond_I));

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
		} while((v1 < mincond_I) || (v2 < mincond_I) || (v1 > maxcond_I) || (v2 > maxcond_I));
		res[node1] = v1;
		res[node2] = v2;
		if(data) {
			data->swap = true;
			data->ncoef = ncoef;
		}
	}
	return res;
}


double *solution_lb_complex::getShuffledSolution_R(shuffleData *data, const shuffler &sh) const
{
	double *res = solution_lb_complex::copySolution(solRe);
	// head or tails
	if(genint(2)) { // Normal
		int ncoef = genint(numcoefficients);
		float minc, maxc;
                  maxc = maxcond;
                  minc = mincond;

		if(sh.shuffleConsts[ncoef]==0) {
			res[ncoef] = minc+genreal()*(maxc-minc);
		} else {
			float val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for(int i=0;i<sh.shuffleConsts[ncoef];i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoef];
				rnd -= 0.5;
				rnd *= (maxc - minc);
				val += rnd;
			} while((val < minc) || (val > maxc));
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

solution_lb_complex *solution_lb_complex::shuffle(shuffleData *data_R, const shuffler &sh_R, shuffleData *data_I, const shuffler &sh_I, bool *real) const
{
  double *s_R;
  double *s_I;
  if(genint(2)) { // real shuffle
    s_R = getShuffledSolution_R(data_R, sh_R);
    s_I = solution_lb_complex::copySolution(solIm);
    *real = true;
  } else {
    s_R = solution_lb_complex::copySolution(solRe);
    s_I = getShuffledSolution_R(data_I, sh_I);;
    *real = false;
  }
	solution_lb_complex *res;
	res = new solution_lb_complex(s_R, s_I, *this);
	return res;
}

void solution_lb_complex::saturate()
{
      ensureMinIt(nodes.size()+30);
}

void solution_lb_complex::ensureMinIt(unsigned int it)
{
      static Eigen::VectorXd aux(gelectrodes.size()-1);
      for(int i = 0; i<nobs;i++) {
            LB_Solver_Complex *sim = this->simulations[i];
            while(sim->getIteration()<it) {
                simulations[i]->do_iteration();
                this->totalit++;

                distance[i] = simulations[i]->getErrorl2Estimate();
                maxdist[i] = simulations[i]->getMaxErrorl2Estimate();
                mindist[i] = simulations[i]->getMinErrorl2Estimate();
                err[i] = maxdist[i]-mindist[i];
                err_x_dist[i] = maxdist[i]*err[i];
                totalDist = distance.norm()+regularisation;
                minTotalDist = mindist.norm()+regularisation;
                maxTotalDist = maxdist.norm()+regularisation;

               // reevaluate critical
                double max = err_x_dist[0];
                critical = 0;
                for(int j = 1; j<nobs;j++) {
                  if(max < err_x_dist[j]) {
                        max = err_x_dist[j];
                        critical = j;
                  }
                }
                critErr = err[critical];
            }
      }
}


solution_lb_complex::~solution_lb_complex()
{
	delete[] solRe;
  delete[] solIm;
	delete Aii_R;
  delete Aii_I;
	delete Aic_R;
  delete Aic_I;
	delete Acc_R;
  delete Acc_I;
	for(int i=0;i<nobs;i++) {
		delete simulations[i];
	}
	delete[] simulations;
}
