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

//void solutioncomplex::zeroSumVector(Eigen::VectorXcd &vec) {
//	std::complex<double> avg = 0;
//	for (int i = 0; i < vec.size(); i++) avg += vec[i];
//	avg /= input->getGenericElectrodesCount();
//	for (int i = 0; i < vec.size(); i++) vec[i] -= avg;
//}

bool solutioncomplex::compareWith(solutionbase &target, double kt, double prob)
{
	double delta, expdelta;
	
	delta = target.totalDist - this->totalDist;
	expdelta = exp(-delta / kt);
	if (delta <= 0) {
		//std::cout << "+";
		return true;
	}
	//std::cout << "-" << expdelta << std::endl;


	// Now randomly accept based on the energy
	if (genreal()<expdelta) return true;
	return false;
}

solutioncomplex::solutioncomplex(const std::complex<double> *sigma, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
solutionbase(sigma, _input, _readings),
sol(solutionbase::copySolution(sigma, _input)),
stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
simulations(new CG_SolverComplex *[_readings->getNObs()])
//input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}


// New random solution
solutioncomplex::solutioncomplex(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
solutionbase(_input, _readings),
sol(solutioncomplex::getNewRandomSolution(_input)),
stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
simulations(new CG_SolverComplex *[_readings->getNObs()])
//input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solutioncomplex::solutioncomplex(std::complex<double> *sigma, const solutioncomplex &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
solutionbase(sigma, base, _input, _readings),
sol(sigma),
stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
simulations(new CG_SolverComplex *[_readings->getNObs()])
//input(_input), readings(_readings)
{
	this->initSimulations(base);
	this->initErrors();
}

// New random solution
solutioncomplexcalibration::solutioncomplexcalibration(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
solutioncomplex(_input, _readings) {}

void solutioncomplex::initSimulations(const solutionbase<std::complex<double>> &base)
{
	// Prepare solvers
	for (int i = 0; i<readings->getNObs(); i++)
	{
		// Reuse previous solutions as initial values
		//simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), base.simulations[i]->getX(), *precond);
		simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), base.getSimulationX(i), *precond);
	}
}

void solutioncomplex::initSimulations()
{
	// Prepare solvers
	for (int i = 0; i<readings->getNObs(); i++)
	{
		simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), *precond);
		for (int k = 0; k < 5; k++) simulations[i]->do_iteration();
	}
}

void solutioncomplex::initErrors()
{
	// Calc regularisation value
	//Just some scrap space to avoid dynamic allocations
	//	WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXcd aux(input->getGenericElectrodesCount());
	// Retrieve distance estimates, errors and boundaries
	for (int i = 0; i<readings->getNObs(); i++) {
		// Compare with observation
		aux = simulations[i]->getX().tail(aux.size());
#ifndef BLOCKGND
		// Rebase tension for zero sum
		zeroSumVector(aux);
#endif
		aux -= readings->getTensions()[i];
		distance[i] = aux.norm();
	}
	totalDist = distance.norm();
}


//std::complex<double> *solutioncomplexcalibration::copySolution(const std::complex<double> *sol, std::shared_ptr<problem> input)
//{
//	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
//
//	for (int i = 0; i<input->getNumCoefficients(); i++)
//		res[i] = sol[i];
//
//	return res;
//}

std::complex<double> *solutioncomplex::getNewRandomSolution(std::shared_ptr<problem> input)
{
	if (input->getCalibrationMode() != 0) return solutioncomplexcalibration::getNewRandomSolution(input);
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
	int i = 0;
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminperm = w*minperm, wmaxperm = w*maxperm;
	for (i = 0; i < input->getNumCoefficients(); i++)
		res[i] = std::complex<double>(mincond + genreal()*(maxcond - mincond), wminperm + genreal()*(wmaxperm - wminperm));

	return res;
}


std::complex<double> *solutioncomplexcalibration::getNewRandomSolution(std::shared_ptr<problem> input)
{
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminpermint = w*minpermint, wmaxpermint = w*maxpermint;
	res[0] = std::complex<double>(mincondint + genreal()*(maxcondint - mincondint), wminpermint + genreal()*(wmaxpermint - wminpermint));
	double wminpermelec = w*minpermelec, wmaxpermelec = w*maxpermelec;
	for (int i = 1; i < input->getNumCoefficients(); i++)
		res[i] = std::complex<double>(mincondelec + genreal()*(maxcondelec - mincondelec), wminpermelec + genreal()*(wmaxpermelec - wminpermelec));

	return res;
}

std::complex<double> *solutioncomplex::getShuffledSolution(shuffleData *data, const shuffler &sh) const {
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
		if (data) {
			data->swap = true;
			data->ncoef = ncoef;
		}
	}
	return res;
}

std::complex<double> *solutioncomplexcalibration::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	std::complex<double> *res = solutioncomplexcalibration::copySolution(sol, input);
	int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;
	double minval, maxval;

	// Real or complex shuffle
	if (genint(2)) {
		minval = ncoef == 0 ? mincondint : mincondelec;
		maxval = ncoef == 0 ? maxcondint : maxcondelec;

		// Real
		double deltamax = (maxval - minval);
		int ncoefreal = 2 * ncoef;

		if (sh.shuffleConsts[ncoefreal] == 0) {
			res[ncoef] = std::complex<double>(minval + genreal()*deltamax, res[ncoef].imag());
		}
		else {
			std::complex<double> val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for (int i = 0; i < sh.shuffleConsts[ncoefreal]; i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoefreal];
				rnd -= 0.5;
				rnd *= deltamax / 4.0;
				val += rnd;
			} while ((val.real() < minval) || (val.real() > maxval));
			res[ncoef] = val;
		}
		if (data) data->ncoef = ncoefreal;
	}
	else {
		// Imaginary
		double w = 2 * M_PI * input->getCurrentFreq();
		minval = ncoef == 0 ? w*minpermint : w*minpermelec;
		maxval = ncoef == 0 ? w*maxpermint : w*maxpermelec;

		double deltamax = (maxval - minval);
		int ncoefimg = 2 * ncoef + 1;

		if (sh.shuffleConsts[ncoefimg] == 0) {
			res[ncoef] = std::complex<double>(res[ncoef].real(), minval + genreal()*deltamax);
		}
		else {
			std::complex<double> val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for (int i = 0; i < sh.shuffleConsts[ncoefimg]; i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoefimg];
				rnd -= 0.5;
				rnd *= deltamax / 4.0;
				val += std::complex<double>(0, rnd);
			} while ((val.imag() < minval) || (val.imag() > maxval));
			res[ncoef] = val;
		}
		if (data) data->ncoef = ncoefimg;
	}
	return res;
}

solutioncomplex *solutioncomplex::shuffle(shuffleData *data, const shuffler &sh) const
{
	std::complex<double> *sigma = getShuffledSolution(data, sh);
	solutioncomplex *res;
	try {
		res = new solutioncomplex(sigma, *this, input, readings);
	}
	catch (...) {
		for (int i = 0; i<65; i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
}

solutioncomplexcalibration *solutioncomplexcalibration::shuffle(shuffleData *data, const shuffler &sh) const
{
	std::complex<double> *sigma = getShuffledSolution(data, sh);
	solutioncomplexcalibration *res;
	try {
		res = new solutioncomplexcalibration(sigma, *this, input, readings);
	}
	catch (...) {
		for (int i = 0; i<65; i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
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
