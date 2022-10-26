/*
 * solution.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#include "solution_lb.h"
#include "random.h"
#include "observations.h"
#include <iostream>
#include "gradientnormregularisation.h"
#include "util/standard_deviation.hpp"

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif
//#define USE_PREVIOUS_DATA
#ifdef USE_PREVIOUS_DATA

const float base[] = {
  // densermesh tri01a.ampl_calibrado
  //0.228468, 0.274277, 0.277748, 0.276384, 0.369586, 0.316081, 0.354511, 0.350053, 0.377487, 0.372943, 0.359078, 0.329245, 0.242468, 0.239446, 0.220835, 0.249638, 0.346346, 0.358612, 0.317659, 0.326589, 0.300691, 0.342779, 0.318762, 0.356848, 0.360735, 0.357942, 0.37484, 0.37047, 0.283102, 0.259762, 0.211597, 0.243358, 0.381399
  // Densermesh, tri01a.ampl_calibrado latter data
  //0.199859, 0.252996, 0.251331, 0.229678, 0.334548, 0.284054, 0.309068, 0.343423, 0.363469, 0.361536, 0.381482, 0.337849, 0.218437, 0.227098, 0.217082, 0.232915, 0.272974, 0.314924, 0.281226, 0.286308, 0.313496, 0.312229, 0.283623, 0.293245, 0.305662, 0.314062, 0.316468, 0.381498, 0.255632, 0.223154, 0.211192, 0.229608, 0.381494
  // Densermesh,  tri02b.ampl_calibrado.txt ,
  0.277581, 0.24364, 0.236519, 0.229084, 0.324995, 0.274411, 0.281469, 0.290241, 0.316553, 0.340117, 0.381493, 0.332735, 0.211533, 0.243014, 0.223994, 0.224349, 0.24471, 0.324983, 0.282342, 0.253281, 0.307704, 0.289102, 0.283188, 0.256502, 0.274983, 0.281164, 0.28214, 0.350598, 0.244814, 0.222791, 0.21463, 0.262809, 0.381499
  //0.204031 ,0.253713 ,0.258056 ,0.22902 ,0.334923 ,0.287182 ,0.309377 ,0.354451 ,0.370216 ,0.366888 ,0.378747 ,0.347058 ,0.222916 ,0.228123 ,0.218939 ,0.238365 ,0.278456 ,0.324395 ,0.288809 ,0.286103 ,0.324509 ,0.321599 ,0.284023 ,0.303139 ,0.312969 ,0.328236 ,0.320193 ,0.375504 ,0.256666 ,0.229927 ,0.213987 ,0.230564 ,0.381306

  // ueta_circuito_0XXX_AMPLITUDES_10mA.txt
 // 0.0962959, 0.202657, 0.136188, 0.205733, 0.171565, 0.165175, 0.163677, 0.151598, 0.151849, 0.136051, 0.137322, 0.123318, 0.114067, 0.123119, 0.118788, 0.112003, 0.121462, 0.126865, 0.121439, 0.130946, 0.158523, 0.153706, 0.167477, 0.151688, 0.176961, 0.165414, 0.170706, 0.183611, 0.176578, 0.147866, 0.14479, 0.143259, 0.381499

  //0.34648, 0.406421, 0.396935, 0.372705, 0.513649, 0.443898, 0.471985, 0.467473, 0.490282, 0.498275, 0.530917, 0.463241, 0.372647, 0.397438, 0.37826, 0.374614, 0.391259, 0.428486, 0.406962, 0.40457, 0.474347, 0.474314, 0.470977, 0.435977, 0.442748, 0.470965, 0.422019, 0.482234, 0.413204, 0.375263, 0.385004, 0.425495, 1.35
  //0.0607283, 0.0665641, 0.0672715, 0.0644222, 0.0831744, 0.0800821, 0.095222, 0.087998, 0.0795147, 0.07844, 0.0820809, 0.0738347, 0.0615422, 0.0666471, 0.0630172, 0.0626265, 0.0654864, 0.0707983, 0.067319, 0.0661263, 0.0802712, 0.0783579, 0.0762317, 0.0713599, 0.0815007, 0.0752741, 0.0709489, 0.0902495, 0.0724872, 0.062839, 0.0674892, 0.0733704
  //0.0607283, 0.0665641, 0.0672715, 0.0644222, 0.0831744, 0.0800821, 0.095222, 0.087998, 0.0795147, 0.07844, 0.0820809, 0.0738347, 0.0615422, 0.0666471, 0.0630172, 0.0626265, 0.0654864, 0.0707983, 0.067319, 0.0661263, 0.0802712, 0.0783579, 0.0762317, 0.0713599, 0.0815007, 0.0752741, 0.0709489, 0.0902495, 0.0724872, 0.062839, 0.0674892, 0.0733704, 0.261434
  //0.13807, 0.125363, 0.106238, 0.0994108, 0.144801, 0.12248, 0.133051, 0.128964, 0.13664, 0.138327, 0.15113, 0.123861, 0.0936534, 0.107083, 0.0969518, 0.0918484, 0.0992923, 0.118951, 0.10938, 0.102443, 0.136631, 0.124224, 0.132069, 0.11404, 0.119842, 0.11939, 0.111648, 0.131918, 0.107573, 0.0992579, 0.105823, 0.133198, 0.22967
/*
    0.09,//0.102741,
    0.24316,
    0.153582,
    0.242394,
    0.201954,
    0.1922,
    0.193071,
    0.176909,
    0.177522,
    0.156914,
    0.159333,
    0.140742,
    0.128991,
    0.141002,
    0.135405,
    0.126562,
    0.138666,
    0.145703,
    0.138299,
    0.150642,
    0.187453,
    0.180399,
    0.199719,
    0.176805,
    0.212709,
    0.196534,
    0.207548,
    0.221224,
    0.214724,
    0.170313,
    0.15991,
    0.16019,
    0.365//0.375*/
};

#endif /* USE_PREVIOUS_DATA */

void solution_lb::improve()
{
	// Do another iteration on the critical solver
	simulations[critical].do_iteration();
	this->totalit++;
	// Recalcule expected distance and boundaries

	distance[critical] = simulations[critical].getErrorl2Estimate();
	maxdist[critical] = simulations[critical].getMaxErrorl2Estimate();
	mindist[critical] = simulations[critical].getMinErrorl2Estimate();
	err[critical] = maxdist[critical]-mindist[critical];
	err_x_dist[critical] = maxdist[critical]*err[critical];
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
	// reevaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for(int i = 1; i<o.getNObs();i++) {
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

// Construct from solution vector copy
solution_lb::solution_lb(std::shared_ptr<problem> p, observations<double> &o, std::vector<double> &&sol):
                                p(p), o(o),
				sol(std::move(sol)),
				distance(o.getNObs()),
				maxdist(o.getNObs()),
				mindist(o.getNObs()),
				err(o.getNObs()),
				err_x_dist(o.getNObs())
{

    this->initMatrices();
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solution_lb::solution_lb(std::vector<double> &&sol, const solution_lb &base):
                p(base.p), o(base.o),
                sol(std::move(sol)),
                distance(o.getNObs()),
                maxdist(o.getNObs()),
                mindist(o.getNObs()),
                err(o.getNObs()),
                err_x_dist(o.getNObs())
{
        this->initMatrices();
        this->initSimulations(base);
        this->initErrors();
}


// New random solution
solution_lb::solution_lb(std::shared_ptr<problem> p, observations<double> &o):
                p(p), o(o),
		sol(solution_lb::getNewRandomSolution(p->getNumCoefficients())),
		distance(o.getNObs()),
		maxdist(o.getNObs()),
		mindist(o.getNObs()),
		err(o.getNObs()),
		err_x_dist(o.getNObs())
{
      this->initMatrices();
      this->initSimulations();
	this->initErrors();
}

void solution_lb::initMatrices()
{
      solver::matrix *_aii, *_aic, *_acc;
      assembleProblemMatrix_lb(sol.data(), &_aii, &_aic, &_acc, *p);
      Aii.reset(_aii);
      Aic.reset(_aic);
      Acc.reset(_acc);
      precond.reset(LB_Solver::makePreconditioner(*Aii, *Aic));
}

void solution_lb::initSimulations()
{
    // Prepare solvers
	int i;
	this->totalit = 0;

    // Reserve so we don't need to move solvers around
    simulations.reserve(o.getNObs());
    // 1st solution estimates also least eigenvalue
    simulations.emplace_back(Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[0].tail(p->getGenericElectrodesCount()),
                    Eigen::VectorXd(o.getTensions()[0]), *precond,  100, (float)0.0001, this->eigenvalue_estimate, this->eigenvector_estimate);
    this->totalit += simulations[0].getIteration();
	for(i=1;i<o.getNObs();i++)
	{
        simulations.emplace_back(
            Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[i].tail(p->getGenericElectrodesCount()),
            Eigen::VectorXd(o.getTensions()[i]), *precond, this->eigenvalue_estimate);
		simulations[i].do_iteration();
		this->totalit += simulations[i].getIteration();
	}
}

void solution_lb::initSimulations(const solution_lb &base)
{
        // Prepare solvers
        int i;
        this->totalit = 0;

        // Reserve so we don't need to move solvers around
        simulations.reserve(o.getNObs());
        // 1st solution estimates also least eigenvalue
        simulations.emplace_back(Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[0].tail(p->getGenericElectrodesCount()),
                        Eigen::VectorXd(o.getTensions()[0]), *precond,
                        base.simulations[0].getX(), base.getLeastEigenvectorEstimation(), 100, (float)0.0001, this->eigenvalue_estimate, this->eigenvector_estimate);
        this->totalit += simulations[0].getIteration();
        for(i=1;i<o.getNObs();i++)
        {
                simulations.emplace_back(
                        Aii.get(), Aic.get(), Acc.get(), o.getCurrents()[i].tail(p->getGenericElectrodesCount()),
                        Eigen::VectorXd(o.getTensions()[i]), *precond, this->eigenvalue_estimate,
					       base.simulations[i].getX());
                simulations[i].do_iteration();
                this->totalit += simulations[i].getIteration();
        }
}


void solution_lb::initErrors()
{
	int i;
	// Calc regularisation value
  double electrodevariance = population_variance(this->sol.begin(), this->sol.begin()+this->p->getGenericElectrodesCount());
	this->regularisation = gradientNormRegularisation::getInstance()->getRegularisation(this->sol.data())*0.010+electrodevariance*20;
	// Retrieve distance estimates, errors and boundaries
	for(i=0;i<o.getNObs();i++) {
		// Compare with observation

		distance[i] = simulations[i].getErrorl2Estimate();
		maxdist[i] = simulations[i].getMaxErrorl2Estimate();
		mindist[i] = simulations[i].getMinErrorl2Estimate();
		err[i] = maxdist[i]-mindist[i];
		err_x_dist[i] = maxdist[i]*err[i];
	}
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
	// evaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for(i = 1; i<o.getNObs();i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}

std::vector<double> solution_lb::getNewRandomSolution(int size)
{
	std::vector<double> res;
	res.reserve(size);
	int i = 0;
#ifdef USE_PREVIOUS_DATA
        for(;i<sizeof(base)/sizeof(*base);i++)
            res.push_back(base[i]);
#endif /* USE_PREVIOUS_DATA */
	for(;i<size;i++)
		res.push_back(mincond+genreal()*(maxcond-mincond));

	return res;
}

std::vector<double> solution_lb::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	std::vector<double> res(sol);
	// head or tails
	if(genint(2)) { // Normal
#ifndef USE_PREVIOUS_DATA
		int ncoef = genint(p->getNumCoefficients());

#else // USE_PREVIOUS_DATA
		int ncoef = genint(p->getNumCoefficients()-(sizeof(base)/sizeof(float)))+(sizeof(base)/sizeof(float));
#endif	// float minc, maxc;
		float minc, maxc;
                  maxc = (float)maxcond;
                  minc = (float)mincond;

		if(sh.shuffleConsts[ncoef]==0) {
			res[ncoef] = minc+genreal()*(maxc-minc);
		} else {
			double val;
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
		int ncoef = genint(p->getInnerAdjacencyCount());
		int node1, node2;

		std::pair<int, int> adj = p->getAdjacency(ncoef);
                node1 = p->getNode2Coefficient(adj.first);
		node2 = p->getNode2Coefficient(adj.second);

		// Order nodes
		if(res[node1]>res[node2]) {
			int aux = node1;
			node1 = node2;;
			node2 = aux;
		}
		double v1 = res[node1], v2 = res[node2];
		double a = max( min(v1-mincond, maxcond-v2), min(maxcond-v1, v2-mincond));

		double delta;
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
	solution_lb *res;
	try {
		res = new solution_lb(getShuffledSolution(data, sh), *this);
	} catch(...) {

		for(int i=0;i<65;i++)
			std::cout << i << ":" << res->sol[i] << std::endl;
		exit(0);
	}
	return res;
}

void solution_lb::saturate()
{
      ensureMinIt(p->getNodesCount()+30);
}

void solution_lb::ensureMinIt(unsigned int it)
{
      static Eigen::VectorXd aux(p->getGenericElectrodesCount()-1);
      for(unsigned int i = 0; i<o.getNObs();i++) {
            while(simulations[i].getIteration()<it) {
                simulations[i].do_iteration();
                this->totalit++;

                distance[i] = simulations[i].getErrorl2Estimate();
                maxdist[i] = simulations[i].getMaxErrorl2Estimate();
                mindist[i] = simulations[i].getMinErrorl2Estimate();
                err[i] = maxdist[i]-mindist[i];
                err_x_dist[i] = maxdist[i]*err[i];
                totalDist = distance.norm()+regularisation;
                minTotalDist = mindist.norm()+regularisation;
                maxTotalDist = maxdist.norm()+regularisation;

               // reevaluate critical
                double max = err_x_dist[0];
                critical = 0;
                for(int j = 1; j<o.getNObs();j++) {
                  if(max < err_x_dist[j]) {
                        max = err_x_dist[j];
                        critical = j;
                  }
                }
                critErr = err[critical];
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
