/*
 * solution.cpp
 *
 *  Created on: Sep 12, 2010
 *      Author: thiago
 */

#include "solutioncomplex.h"
#include "random.h"
//#include "observations.h"
//#include "problemdescription.h"
#include <iostream>
//#include <boost/numeric/interval.hpp>
#include "gradientnormregularisationcomplex.h"
#define _USE_MATH_DEFINES
#include <math.h>

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

void solutioncomplex::zeroSumVector(Eigen::VectorXcd &vec) {
	std::complex<double> avg = 0;
	for (int i = 0; i < vec.size(); i++) avg += vec[i];
	avg /= input->getGenericElectrodesCount();
	for (int i = 0; i < vec.size(); i++) vec[i] -= avg;
}

void solutioncomplex::improve()
{
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXcd aux(input->getGenericElectrodesCount());

	// Do another iteration on the critical solver
	double rnorm = simulations[critical]->getResidueSquaredNorm();
	if (simulations[critical]->getResidueSquaredNorm() == 0.0) return;
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
	for(int i = 1; i<readings->getNObs();i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}

bool solutioncomplex::compareWith(solutioncomplex &target, double kt, double prob)
{
	double delta, expdelta;
	//// Ensure errors are within required margin
	//while(true) {
	//	double min_delta = target.minTotalDist - this->maxTotalDist;
	//	double max_delta = target.maxTotalDist - this->minTotalDist;
	//	delta = target.totalDist - this->totalDist;
	//	expdelta = exp(-delta/kt);
	//	// Check boundary conditions:
	//	// Upper bound is negative, no doubt here
	//	if(max_delta < 0) break;
	//	// Estimate is negative, but upper bound is positive
	//	else if(delta <= 0) {
	//		if(exp(-max_delta/kt)>=(1-prob)) break; // Upper condition only
	//	}
	//	// Estimate and upper bounds are positive, lower bound is negative
	//	else if(min_delta <= 0) {
	//		if(expdelta >= 1 - prob) { // Lower condition
	//			if(expdelta <= prob) break;	// upper condition
	//			if(exp(-max_delta/kt)>=(expdelta-prob)) break;
	//		}
	//	}
	//	// Classic case, everything is positive
	//	else {
	//		if(exp(-min_delta/kt)<=prob+expdelta) { // lower condition
	//			if(expdelta <= prob) break;	// upper condition
	//			if(exp(-max_delta/kt)>=(expdelta-prob)) break;
	//		}
	//	}
	//	// Not there yet, improve boundaries
	//	// Select wich one to improve
	//	if(this->critErr > target.critErr)
	//		this->improve();
	//	else
	//		target.improve();
	//}

	for (int i = 0; i < 10 ; i++) {
		target.improve();
	}
	
	delta = target.totalDist - this->totalDist;
	expdelta = exp(-delta / kt);
	if(delta <= 0) {
		//std::cout << "+";
		return true;
	}
	//std::cout << "-" << expdelta << std::endl;


	// Now randomly accept based on the energy
	if(genreal()<expdelta) return true;
	return false;
}


bool solutioncomplex::compareWithMinIt(solutioncomplex &target, double kt, int minit)
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

bool solutioncomplex::compareWithMaxE2(solutioncomplex &target, double kt, double e2)
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


solutioncomplex::solutioncomplex(const std::complex<double> *sigma, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
				sol(solutioncomplex::copySolution(sigma, _input)),
				stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
				precond(new SparseIncompleteLLTComplex(*stiffness)),
				simulations(new CG_SolverComplex *[_readings->getNObs()]),
				distance(_readings->getNObs()),
				maxdist(_readings->getNObs()),
				mindist(_readings->getNObs()),
				err(_readings->getNObs()),
				err_x_dist(_readings->getNObs()),
				input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}


// New random solution
solutioncomplex::solutioncomplex(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
		sol(solutioncomplex::getNewRandomSolution(_input)),
		stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
		precond(new SparseIncompleteLLTComplex(*stiffness)),
		simulations(new CG_SolverComplex *[_readings->getNObs()]),
		distance(_readings->getNObs()),
		maxdist(_readings->getNObs()),
		mindist(_readings->getNObs()),
		err(_readings->getNObs()),
		err_x_dist(_readings->getNObs()),
		input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solutioncomplex::solutioncomplex(std::complex<double> *sigma, const solutioncomplex &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
		sol(sigma),
		stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
		precond(new SparseIncompleteLLTComplex(*stiffness)),
		simulations(new CG_SolverComplex *[_readings->getNObs()]),
		distance(_readings->getNObs()),
		maxdist(_readings->getNObs()),
		mindist(_readings->getNObs()),
		err(_readings->getNObs()),
		err_x_dist(_readings->getNObs()),
		input(_input), readings(_readings)
{
	this->initSimulations(base);
	this->initErrors();
}

void solutioncomplex::initSimulations(const solutioncomplex &base)
{
	// Prepare solvers
	int i;
	this->totalit = 0;
	for (i = 0; i<readings->getNObs(); i++)
	{
		// Reuse previous solutions as initial values
		simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), base.simulations[i]->getX(), *precond);
		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while (ndecr<2 && simulations[i]->getResidueSquaredNorm() != 0.0) {
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

void solutioncomplex::initSimulations()
{
	// Prepare solvers
	int i;
	this->totalit = 0;
	for (i = 0; i<readings->getNObs(); i++)
	{
		simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), *precond);
		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux = 1;
		int ndecr = 0;
		while(ndecr<2 && aux != 0.0) {
			simulations[i]->do_iteration();
			aux = simulations[i]->getErrorl2Estimate();
			if (aux >= err) ndecr = 0;
			else {
				ndecr++;
			}
			err = aux;
		}
		this->totalit += simulations[i]->getIteration();
	}
}

void solutioncomplex::initErrors()
{
	// Calc regularisation value
	this->regularisation = std::abs(gradientNormRegularisationComplex::getInstance()->getRegularisation(this->sol))*30;
	int i;
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXcd aux(input->getGenericElectrodesCount());
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
		maxdist[i] = distance[i] + err[i];
		mindist[i] = max(distance[i] - err[i],0);
		double distancei = distance[i];
		double erri = err[i];
		double maxdisti = maxdist[i];
		double mindisti = mindist[i];
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


std::complex<double> *solutioncomplex::copySolution(const std::complex<double> *sol, std::shared_ptr<problem> input)
{
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];

	for (int i = 0; i<input->getNumCoefficients(); i++)
		res[i] = sol[i];

	return res;
}

std::complex<double> *solutioncomplex::getNewRandomSolution(std::shared_ptr<problem> input)
{
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
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
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminperm = w*minperm, wmaxperm = w*maxperm;
	for (i = 0; i < input->getNumCoefficients(); i++)
		res[i] = std::complex<double>(mincond + genreal()*(maxcond - mincond), wminperm + genreal()*(wmaxperm - wminperm));

	return res;
}

void solutioncomplex::saveMesh(double *sol, const char *filename, std::shared_ptr<problem> input, int step) {
	std::ofstream myfile;
	myfile.open(filename);

	std::ifstream inputfile(input->getMeshFilename());
	for (int i = 0; inputfile.eof() != true; i++) {
		std::string line;
		std::getline(inputfile, line);
		myfile << line << '\n';
	}

	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
	myfile << "$NodeData\n1\n\"Conductivity\"\n1\n0.0\n3\n" << step << "\n1\n" << input->getNodesCount() << "\n";
	for (int j = 0; j < input->getNodesCount(); j++) {
		myfile << (j + 1) << "\t" << sol[input->getNode2Coefficient(j)] << "\n";
	}
	myfile << "$EndNodeData\n";
	myfile.flush();
	myfile.close();
}

void solutioncomplex::savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem> input, observations<std::complex<double>> *readings) {
	std::string refname(filename), imfname(filename), absfname(filename), angfname(filename);
	std::size_t dotfound = refname.find_last_of(".");
	refname.replace(dotfound, 1, "_re."); imfname.replace(dotfound, 1, "_im."); absfname.replace(dotfound, 1, "_abs."); angfname.replace(dotfound, 1, "_ang.");
	std::ofstream myfilereal(refname), myfileimag(imfname), myfileabs(absfname), myfileang(angfname);

	std::ifstream inputfile(input->getMeshFilename());
	for (int i = 0; inputfile.eof() != true; i++) {
		std::string line;
		std::getline(inputfile, line);
		myfilereal << line << '\n'; myfileimag << line << '\n';  myfileabs << line << '\n';  myfileang << line << '\n';
	}

	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
	for (int patterno = 0; patterno < sols.size(); patterno++) {
		myfilereal << "$NodeData\n1\n\"Real Eletric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
		myfileimag << "$NodeData\n1\n\"Imaginary Eletric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
		myfileabs << "$NodeData\n1\n\"Magnitude Eletric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
		myfileang << "$NodeData\n1\n\"Phase Angle Eletric Potential\"\n1\n0.0\n3\n" << patterno << "\n1\n" << input->getNodesCount() << "\n";
		for (int j = 0; j < input->getNodesCount(); j++) {
			myfilereal << (j + 1) << "\t" << sols[patterno][j].real() << "\n";
			myfileimag << (j + 1) << "\t" << sols[patterno][j].imag() << "\n";
			myfileabs << (j + 1) << "\t" << std::abs(sols[patterno][j]) << "\n";
			myfileang << (j + 1) << "\t" << std::arg(sols[patterno][j]) << "\n";
		}
		myfilereal << "$EndNodeData\n"; myfileimag << "$EndNodeData\n"; myfileabs << "$EndNodeData\n"; myfileang << "$EndNodeData\n";
	}
	myfilereal.flush(); myfileimag.flush(); myfileabs.flush(); myfileang.flush();
	myfilereal.close(); myfileimag.close(); myfileabs.close(); myfileang.close();
}

std::complex<double> *solutioncomplex::getShuffledSolution(shuffleData *data, const shufflercomplex &sh) const
{
	std::complex<double> *res = solutioncomplex::copySolution(sol, input);
	// Real or complex shuffle
	if (genint(2)) { // Real
		// head or tails
		if (genint(2)) { // Normal
			int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;

			if (sh.shuffleConsts[ncoef] == 0) {
				res[ncoef] = std::complex<double>(mincond + genreal()*(maxcond - mincond), res[ncoef].imag());
			}
			else {
				std::complex<double> val;
				do {
					val = res[ncoef];
					double rnd = 0;
					for (int i = 0; i < sh.shuffleConsts[ncoef]; i++)
						rnd += genreal();
					rnd /= sh.shuffleConsts[ncoef];
					rnd -= 0.5;
					rnd *= (maxcond - mincond);
					val += rnd;
				} while ((val.real() < mincond) || (val.real() > maxcond));
				res[ncoef] = val;
			}
			if (data) {
				data->swap = false;
				data->ncoef = ncoef;
			}
		}
		else { // swap
			int ncoef = genint(input->getInnerAdjacencyCount());
			int node1, node2;

			node1 = input->node2coefficient[input->innerAdjacency[ncoef].first];
			node2 = input->node2coefficient[input->innerAdjacency[ncoef].second];

			// Order nodes
			if (res[node1].real() > res[node2].real()) {
				int aux = node1;
				node1 = node2;;
				node2 = aux;
			}
			std::complex<double> v1 = res[node1], v2 = res[node2];
			double a = max(min(v1.real() - mincond, maxcond - v2.real()), min(maxcond - v1.real(), v2.real() - mincond));

			double delta;
			do {
				if (sh.swapshuffleconsts[ncoef] == 0) {
					delta = a*(genreal() * 2 - 1);
				}
				else {
					double rnd = 0;
					for (int i = 0; i < sh.swapshuffleconsts[ncoef]; i++)
						rnd += genreal();
					rnd /= sh.swapshuffleconsts[ncoef];
					rnd -= 0.5;
					delta = a*rnd;
				}
				v1 = res[node1] - delta;
				v2 = res[node2] + delta;
			} while ((v1.real() < mincond) || (v2.real() < mincond) || (v1.real() > maxcond) || (v2.real() > maxcond));
			res[node1] = v1;
			res[node2] = v2;
			if (data) {
				data->swap = true;
				data->ncoef = ncoef;
			}
		}
		return res;
	}

	// Imaginary
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminperm = w*minperm, wmaxperm = w*maxperm;
	if (genint(2)) { // Normal
		int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;

		if (sh.shuffleConsts[ncoef] == 0) {
			res[ncoef] = std::complex<double>(res[ncoef].real(), wminperm + genreal()*(wmaxperm - wminperm));
		}
		else {
			std::complex<double> val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for (int i = 0; i < sh.shuffleConsts[ncoef]; i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoef];
				rnd -= 0.5;
				rnd *= (wmaxperm - wminperm);
				val += std::complex<double>(0, rnd);
			} while ((val.imag() < wminperm) || (val.imag() > wmaxperm));
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
		if (res[node1].imag() > res[node2].imag()) {
			int aux = node1;
			node1 = node2;;
			node2 = aux;
		}
		std::complex<double> v1 = res[node1], v2 = res[node2];
		double a = max(min(v1.imag() - wminperm, wmaxperm - v2.imag()), min(wmaxperm - v1.imag(), v2.imag() - wminperm));

		double delta;
		do {
			if (sh.swapshuffleconsts[ncoef] == 0) {
				delta = a*(genreal() * 2 - 1);
			}
			else {
				double rnd = 0;
				for (int i = 0; i < sh.swapshuffleconsts[ncoef]; i++)
					rnd += genreal();
				rnd /= sh.swapshuffleconsts[ncoef];
				rnd -= 0.5;
				delta = a*rnd;
			}
			v1 = res[node1] - std::complex<double>(0, delta);
			v2 = res[node2] + std::complex<double>(0, delta);
		} while ((v1.imag() < wminperm) || (v2.imag() < wminperm) || (v1.imag() > wmaxperm) || (v2.imag() > wmaxperm));
		res[node1] = v1;
		res[node2] = v2;
		if(data) {
			data->swap = true;
			data->ncoef = ncoef;
		}
	}
	return res;
}

void shufflercomplex::addShufflerFeedback(const shuffleData &data, bool pos)
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

solutioncomplex *solutioncomplex::shuffle(shuffleData *data, const shufflercomplex &sh) const
{
	std::complex<double> *sigma = getShuffledSolution(data, sh);
	solutioncomplex *res;
	try {
		res = new solutioncomplex(sigma, *this, input, readings);
	} catch(...) {
		for(int i=0;i<65;i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
}

void solutioncomplex::saturate()
{
      ensureMinIt(input->getNodesCount()+30);
}

void solutioncomplex::ensureMinIt(unsigned int it)
{     
	static Eigen::VectorXcd aux(input->getGenericElectrodesCount());
	for (int i = 0; i<readings->getNObs(); i++) {
	    CG_SolverComplex *sim = this->simulations[i];
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
		//mindist[i] = max(distance[i] - err[i],0);
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

void solutioncomplex::ensureMaxE2(double e2)
{     
	static Eigen::VectorXcd aux(input->getGenericElectrodesCount());
	for (int i = 0; i<readings->getNObs(); i++) {
	    CG_SolverComplex *sim = this->simulations[i];
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
		//mindist[i] = max(distance[i] - err[i],0);
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

solutioncomplex::~solutioncomplex()
{
	delete[] sol;
	delete stiffness;
	delete stiffnessorig;
	delete precond;
	for (int i = 0; i<readings->getNObs(); i++) {
		delete simulations[i];
	}
	delete[] simulations;
}

