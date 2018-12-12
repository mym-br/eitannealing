/*
 * solution.cpp
 *
 *  Created on: Sep 12, 2010
 *      Author: thiago
 */

#include "solution.h"
#include "random.h"
//#include "observations.h"
//#include "problemdescription.h"
#include <iostream>
//#include <boost/numeric/interval.hpp>
#include "gradientnormregularisation.h"
#include "intcoef.h"
#include <numeric>

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
	static Eigen::VectorXd aux(input->getGenericElectrodesCount());

	// Do another iteration on the critical solver
	simulations[critical]->do_iteration();
	this->totalit++;
	// Recalcule expected distance and boundaries

	aux = simulations[critical]->getX().tail(input->getGenericElectrodesCount());
	#ifndef BLOCKGND
	// Rebase tension for zero sum
	zeroSumVector(aux);
	#endif
	aux -= readings->getTensions()[critical];

	distance[critical] = aux.norm();
	err[critical] = sqrt(simulations[critical]->getErrorl2Estimate());
	maxdist[critical] = distance[critical] + err[critical];
	mindist[critical] = max(distance[critical] - err[critical],0);
	err_x_dist[critical] = maxdist[critical]*err[critical];
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
	// reevaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for (int i = 1; i<readings->getNObs(); i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}

bool solution::compareWith(solutionbase &target, double kt, double prob)
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


bool solution::compareWithMinIt(solutionbase &target, double kt, int minit)
{
	double delta, expdelta;
	// Ensure errors are within required margin
	ensureMinIt(minit);
	target.ensureMinIt(minit);
	delta = target.totalDist - this->totalDist;
	expdelta = exp(-delta/kt);
	if(delta <= 0) {
	  return true;
	}
	if(genreal()<expdelta) return true;
	return false;
}

bool solution::compareWithMaxE2(solutionbase &target, double kt, double e2)
{
	double delta, expdelta;
	// Ensure errors are within required margin
	ensureMaxE2(e2);
	target.ensureMaxE2(e2);
	delta = target.totalDist - this->totalDist;
	expdelta = exp(-delta/kt);
	if(delta <= 0) {
	  return true;
	}
	if(genreal()<expdelta) return true;
	return false;
}


solution::solution(const double *sigma, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs) :
		sol(solution::copySolution(sigma, _input)),
		stiffness(solution::getNewStiffness(sol, _input)),
		precond(new SparseIncompleteLLT(*stiffness)),
		simulations(new CG_Solver *[_readings->getNObs()]),
		fixedCoeffs(_fixedCoeffs),
		solutionbase(sigma, _input, _readings),
		intcoef(new intCoef(*_input))
{
	this->initSimulations();
	this->initErrors();
}


// New random solution
solution::solution(std::shared_ptr<problem> _input, observations<double> *_readings, std::vector<double> &electrodesCoeffs) :
		sol(solution::getNewRandomSolution(_input, electrodesCoeffs)),
		stiffness(solution::getNewStiffness(sol,  _input)),
		precond(new SparseIncompleteLLT(*stiffness)),
		simulations(new CG_Solver *[_readings->getNObs()]),
		fixedCoeffs(electrodesCoeffs.size()),
		solutionbase(_input, _readings),
		intcoef(new intCoef(*_input))
{
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solution::solution(double *sigma, const solution &base, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs) :
		sol(sigma),
		stiffness(solution::getNewStiffness(sol, _input)),
		precond(new SparseIncompleteLLT(*stiffness)),
		simulations(new CG_Solver *[_readings->getNObs()]),
		fixedCoeffs(_fixedCoeffs),
		solutionbase(sigma, base, _input, _readings),
		intcoef(base.intcoef)
{
	this->initSimulations(base);
	this->initErrors();
}

void solution::initSimulations(const solution &base)
{
	// Prepare solvers
	int i;
	this->totalit = 0;
	for(i=0;i<readings->getNObs();i++)
	{
		// Reuse previous solutions as initial values
		simulations[i] = new CG_Solver(*stiffness, input->getCurrentVector(i, readings), base.simulations[i]->getX(), *precond);
		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while(ndecr<2) {
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
	for (i = 0; i<readings->getNObs(); i++)
	{
		simulations[i] = new CG_Solver(*stiffness, input->getCurrentVector(i, readings), *precond);
		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while(ndecr<2) {
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
	// Calc regularisation value
	this->regularisation = gradientNormRegularisation::getInstance()->getRegularisation(this->sol)*30;
	//for (int i = 0; i < input->getGenericElectrodesCount(); i++) { std::cout << this->sol[i] << " "; std::cout << std::endl; }
	if (std::abs(input->electrodevar) > 1e-6) {
		double sum = std::accumulate(this->sol, this->sol + input->getGenericElectrodesCount(), 0.0);
		double mean = sum / input->getGenericElectrodesCount();
		std::vector<double> diff(input->getGenericElectrodesCount());
		std::transform(this->sol, this->sol + input->getGenericElectrodesCount(), diff.begin(), [mean](double x) { return x - mean; });
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		double var = sq_sum / input->getGenericElectrodesCount();
		this->elecvariance = var*input->electrodevar;
		this->regularisation += this->elecvariance;
	}
	else
		this->elecvariance = 0;
	int i;
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXd aux(input->getGenericElectrodesCount());
	// Retrieve distance estimates, errors and boundaries
	for (i = 0; i<readings->getNObs(); i++) {
		// Compare with observation
		aux = simulations[i]->getX().tail(aux.size());
		#ifndef BLOCKGND
		// Rebase tension for zero sum
		zeroSumVector(aux);
		#endif
		aux -= readings->getTensions()[i];

		distance[i] = aux.norm();
		err[i] = sqrt(simulations[i]->getErrorl2Estimate());
		double distancei = distance[i];
		double erri = err[i];
		maxdist[i] = distance[i] + err[i];
		mindist[i] = max(distance[i] - err[i],0);
		err_x_dist[i] = maxdist[i]*err[i];


	}
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
	// evaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for (i = 1; i<readings->getNObs(); i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}


//double *solution::copySolution(const double *sol, std::shared_ptr<problem> input)
//{
//	double *res = new double[input->getNumCoefficients()];
//
//	for (int i = 0; i<input->getNumCoefficients(); i++)
//		res[i] = sol[i];
//
//	return res;
//}

double *solution::getNewRandomSolution(std::shared_ptr<problem> input, std::vector<double> &electrodesCoeffs)
{
	double *res = new double[input->getNumCoefficients()];
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
	for (i = 0; i < electrodesCoeffs.size(); i++)
		res[i] = electrodesCoeffs[i];
	for (; i<input->getNumCoefficients(); i++)
		res[i] = mincond+genreal()*(maxcond-mincond);

	return res;
}
//
//void solution::saveMesh(double *sol, const char *filename, const char *propertyname, std::shared_ptr<problem> input, int step) {
//	std::ofstream myfile;
//	myfile.open(filename);
//
//	std::ifstream inputfile(input->getMeshFilename());
//	for (int i = 0; inputfile.eof() != true; i++) {
//		std::string line;
//		std::getline(inputfile, line);
//		myfile << line << '\n';
//	}
//
//	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
//	myfile << "$NodeData\n1\n\"" << propertyname << "\"\n1\n0.0\n3\n" << step << "\n1\n" << input->getNodesCount() << "\n";
//	for (int j = 0; j < input->getNodesCount(); j++) {
//		myfile << (j + 1) << "\t" << sol[input->getNode2Coefficient(j)] << "\n";
//	}
//	myfile << "$EndNodeData\n";
//	myfile.flush();
//	myfile.close();
//}
//
//void solution::savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings) {
//	std::ofstream myfile;
//	myfile.open(filename);
//
//	std::ifstream inputfile(input->getMeshFilename());
//	for (int i = 0; inputfile.eof() != true; i++) {
//		std::string line;
//		std::getline(inputfile, line);
//		myfile << line << '\n';
//	}
//
//	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
//	for (int patterno = 0; patterno < sols.size(); patterno++) {
//		myfile << "$NodeData\n1\n\"Electric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
//		for (int j = 0; j < input->getNodesCount(); j++)
//			myfile << (j + 1) << "\t" << sols[patterno][j] * readings->getCurrentVal(patterno) << "\n";
//		myfile << "$EndNodeData\n";
//	}
//	myfile.flush();
//	myfile.close();
//}

double *solution::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	double *res = solutionbase<double>::copySolution(sol, input);
	// head or tails
	if(genint(2)) { // Normal
		int ncoef = fixedCoeffs + genint(input->getNumCoefficients() - fixedCoeffs);	// Lower values fixed;

		if(sh.shuffleConsts[ncoef]==0) {
			res[ncoef] = mincond+genreal()*(maxcond-mincond);
		} else {
			double val;
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
		int ncoef = genint(input->getInnerAdjacencyCount());
		int node1, node2;

		node1 = input->node2coefficient[input->innerAdjacency[ncoef].first];
		node2 = input->node2coefficient[input->innerAdjacency[ncoef].second];

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

//void shuffler::addShufflerFeedback(const shuffleData &data, bool pos)
//{
//	if(pos) { // positive feedback
//		if(data.swap)
//			this->swapshuffleconsts[data.ncoef] /= 2;
//		else
//			this->shuffleConsts[data.ncoef] /= 2;
//	} else {
//		if(data.swap)
//			this->swapshuffleconsts[data.ncoef]++;
//		else
//			this->shuffleConsts[data.ncoef]++;
//	}
//}

solution *solution::shuffle(shuffleData *data, const shuffler &sh) const
{
	double *sigma = getShuffledSolution(data, sh);
	solution *res;
	try {
		res = new solution(sigma, *this, input, readings, fixedCoeffs);
	} catch(...) {
		for(int i=0;i<65;i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
}

void solution::saturate()
{
      ensureMinIt(input->getNodesCount()+30);
}

void solution::ensureMinIt(unsigned int it)
{
	static Eigen::VectorXd aux(input->getGenericElectrodesCount());
	for (int i = 0; i<readings->getNObs(); i++) {
	    CG_Solver *sim = this->simulations[i];
	    while(sim->getIteration()<it) {
		simulations[i]->do_iteration();
		this->totalit++;
		// Recalcule expected distance and boundaries
		aux = simulations[i]->getX().tail(input->getGenericElectrodesCount());
		#ifndef BLOCKGND
		// Rebase tension for zero sum
		zeroSumVector(aux);
		#endif
		aux -= readings->getTensions()[i];
		distance[i] = aux.norm();
		err[i] = sqrt(simulations[i]->getErrorl2Estimate());
		maxdist[i] = distance[i] + err[i];
		mindist[i] = max(distance[i] - err[i],0);
		err_x_dist[i] = maxdist[i]*err[i];
		totalDist = distance.norm();
		minTotalDist = mindist.norm();
		maxTotalDist = maxdist.norm();
		// reevaluate critical
		double max = err_x_dist[0];
		critical = 0;
		for (int j = 1; j<readings->getNObs(); j++) {
		  if(max < err_x_dist[j]) {
			max = err_x_dist[j];
			critical = j;
		  }
		}
		critErr = err[critical];
	    }
      }
}

void solution::ensureMaxE2(double e2)
{
	static Eigen::VectorXd aux(input->getGenericElectrodesCount());
	for (int i = 0; i<readings->getNObs(); i++) {
	    CG_Solver *sim = this->simulations[i];
	    while(sim->getLastE2()>e2) {
		simulations[i]->do_iteration();
		this->totalit++;
		// Recalcule expected distance and boundaries
		aux = simulations[i]->getX().tail(input->getGenericElectrodesCount());
		#ifndef BLOCKGND
		// Rebase tension for zero sum
		zeroSumVector(aux);
		#endif
		aux -= readings->getTensions()[i];
		distance[i] = aux.norm();
		err[i] = sqrt(simulations[i]->getErrorl2Estimate());
		maxdist[i] = distance[i] + err[i];
		mindist[i] = max(distance[i] - err[i],0);
		err_x_dist[i] = maxdist[i]*err[i];
		totalDist = distance.norm();
		minTotalDist = mindist.norm();
		maxTotalDist = maxdist.norm();
		// reevaluate critical
		double max = err_x_dist[0];
		critical = 0;
		for (int j = 1; j<readings->getNObs(); j++) {
		  if(max < err_x_dist[j]) {
			max = err_x_dist[j];
			critical = j;
		  }
		}
		critErr = err[critical];
	    }
      }
}

solution::~solution()
{
	delete[] sol;
	delete stiffness;
	delete precond;
	for (int i = 0; i<readings->getNObs(); i++) {
		delete simulations[i];
	}
	delete[] simulations;
}

void solution::savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings) {
	std::ofstream myfile;
	myfile.open(filename);

	std::ifstream inputfile(input->getMeshFilename());
	for (int i = 0; inputfile.eof() != true; i++) {
		std::string line;
		std::getline(inputfile, line);
		myfile << line << '\n';
	}

	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
	for (int patterno = 0; patterno < sols.size(); patterno++) {
		myfile << "$NodeData\n1\n\"Electric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
		for (int j = 0; j < input->getNodesCount(); j++)
			myfile << (j + 1) << "\t" << sols[patterno][j] * readings->getCurrentVal(patterno) << "\n";
		myfile << "$EndNodeData\n";
	}
	myfile.flush();
	myfile.close();
}

void solution::savePotentials(std::vector<Eigen::VectorXf> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings) {
	std::ofstream myfile;
	myfile.open(filename);

	std::ifstream inputfile(input->getMeshFilename());
	for (int i = 0; inputfile.eof() != true; i++) {
		std::string line;
		std::getline(inputfile, line);
		myfile << line << '\n';
	}

	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
	for (int patterno = 0; patterno < sols.size(); patterno++) {
		myfile << "$NodeData\n1\n\"Electric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
		for (int j = 0; j < input->getNodesCount(); j++)
			myfile << (j + 1) << "\t" << sols[patterno][j] * readings->getCurrentVal(patterno) << "\n";
		myfile << "$EndNodeData\n";
	}
	myfile.flush();
	myfile.close();
}
