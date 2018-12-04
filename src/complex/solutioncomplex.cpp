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

solutioncomplex::solutioncomplex(const std::complex<double> *sigma, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, int _fixedCoeffs) :
solutionbase(sigma, _input, _readings),
sol(solutionbase::copySolution(sigma, _input)),
stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
fixedCoeffs(_fixedCoeffs),
simulations(new CG_SolverComplex *[_readings->getNObs()])
//input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}


// New random solution
solutioncomplex::solutioncomplex(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, std::vector<std::complex<double>> &electrodesCoeffs) :
solutionbase(_input, _readings),
sol(solutioncomplex::getNewRandomSolution(_input, electrodesCoeffs)),
stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
fixedCoeffs(electrodesCoeffs.size()),
simulations(new CG_SolverComplex *[_readings->getNObs()])
//input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solutioncomplex::solutioncomplex(std::complex<double> *sigma, const solutioncomplex &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, int _fixedCoeffs) :
solutionbase(sigma, base, _input, _readings),
sol(sigma),
stiffness(solutioncomplex::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
fixedCoeffs(_fixedCoeffs),
simulations(new CG_SolverComplex *[_readings->getNObs()])
//input(_input), readings(_readings)
{
	this->initSimulations(base);
	this->initErrors();
}

// New random solution
solutioncomplexcalibration::solutioncomplexcalibration(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, std::vector<std::complex<double>> &electrodesCoeffs) :
solutioncomplex(_input, _readings, electrodesCoeffs) {}

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

std::complex<double> *solutioncomplex::getNewRandomSolution(std::shared_ptr<problem> input, std::vector<std::complex<double>> &electrodesCoeffs)
{
	if (input->getCalibrationMode() != 0) return solutioncomplexcalibration::getNewRandomSolution(input);
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
	int i = 0;
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminperm = w*minperm, wmaxperm = w*maxperm;
	for (i = 0; i < electrodesCoeffs.size(); i++)
		res[i] = std::complex<double>(electrodesCoeffs[i].real(), w*electrodesCoeffs[i].imag());
	std::complex<double> val  = std::complex<double>(mincond + genreal()*(maxcond - mincond), wminperm + genreal()*(wmaxperm - wminperm)); 
	for (; i < input->getNumCoefficients(); i++)
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
	std::complex<double> *res = solutionbase::copySolution(sol, input);

	// Real or complex shuffle
	if (genint(2)) { // Real
		// head or tails
		if (genint(2)) { // Normal
			int ncoef = fixedCoeffs + genint(input->getNumCoefficients() - fixedCoeffs);	// Lower values fixed;
			int ncoefidx = 2 * ncoef;
			res[ncoef] = std::complex<double>(solutionbase::calcNewShuffledValue(ncoefidx, res[ncoef].real(), data, sh, mincond, maxcond), res[ncoef].imag());
		}
		else { // swap
			int ncoef = genint(input->getInnerAdjacencyCount());	// Lower values fixed;
			int ncoefidx = 2 * ncoef;
			int node1, node2;
			node1 = input->node2coefficient[input->innerAdjacency[ncoef].first];
			node2 = input->node2coefficient[input->innerAdjacency[ncoef].second];

			std::pair<double, double> vals = solutionbase::calcNewSwappedValue(ncoefidx, node1, node2, res[node1].real(), res[node2].real(), data, sh, mincond, maxcond);
			res[node1] = std::complex<double>(vals.first, res[node1].imag());
			res[node2] = std::complex<double>(vals.second, res[node1].imag());
		}
		return res;
	}

	// Imaginary
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminperm = w*minperm, wmaxperm = w*maxperm;
	// head or tails
	if (genint(2)) { // Normal
		int ncoef = fixedCoeffs + genint(input->getNumCoefficients() - fixedCoeffs);	// Lower values fixed;
		int ncoefidx = 2 * ncoef;
		res[ncoef] = std::complex<double>(res[ncoef].real(), solutionbase::calcNewShuffledValue(ncoefidx, res[ncoef].imag(), data, sh, wminperm, wmaxperm));
	}
	else { // swap
		int ncoef = genint(input->getInnerAdjacencyCount());	// Lower values fixed;
		int ncoefidx = 2 * ncoef;
		int node1, node2;
		node1 = input->node2coefficient[input->innerAdjacency[ncoef].first];
		node2 = input->node2coefficient[input->innerAdjacency[ncoef].second];

		std::pair<double, double> vals = solutionbase::calcNewSwappedValue(ncoefidx, node1, node2, res[node1].imag(), res[node2].imag(), data, sh, wminperm, wmaxperm);
		res[node1] = std::complex<double>(res[node1].real(), vals.first);
		res[node2] = std::complex<double>(res[node1].real(), vals.second);
	}
	
	return res;
}

std::complex<double> *solutioncomplexcalibration::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	std::complex<double> *res = solutionbase::copySolution(sol, input);
	int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;
	int ncoefidx = 2 * ncoef;
	double minval, maxval;

	// Real or complex shuffle
	if (genint(2)) { // Real
		// head or tails
		minval = ncoef == 0 ? mincondint : mincondelec;
		maxval = ncoef == 0 ? maxcondint : maxcondelec;
		res[ncoef] = std::complex<double>(solutionbase::calcNewShuffledValue(ncoefidx, res[ncoef].real(), data, sh, minval, maxval), res[ncoef].imag());
		return res;
	}

	// Imaginary
	double w = 2 * M_PI * input->getCurrentFreq();
	minval = ncoef == 0 ? w*minpermint : w*minpermelec;
	maxval = ncoef == 0 ? w*maxpermint : w*maxpermelec;
	// head or tails
	res[ncoef] = std::complex<double>(res[ncoef].real(), solutionbase::calcNewShuffledValue(ncoefidx, res[ncoef].imag(), data, sh, minval, maxval));
	return res;

	//std::complex<double> *res = solutioncomplexcalibration::copySolution(sol, input);
	//int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;
	//double minval, maxval;

	//// Real or complex shuffle
	//if (genint(2)) {
	//	minval = ncoef == 0 ? mincondint : mincondelec;
	//	maxval = ncoef == 0 ? maxcondint : maxcondelec;

	//	// Real
	//	double deltamax = (maxval - minval);
	//	int ncoefreal = 2 * ncoef;

	//	if (sh.shuffleConsts[ncoefreal] == 0) {
	//		res[ncoef] = std::complex<double>(minval + genreal()*deltamax, res[ncoef].imag());
	//	}
	//	else {
	//		std::complex<double> val;
	//		do {
	//			val = res[ncoef];
	//			double rnd = 0;
	//			for (int i = 0; i < sh.shuffleConsts[ncoefreal]; i++)
	//				rnd += genreal();
	//			rnd /= sh.shuffleConsts[ncoefreal];
	//			rnd -= 0.5;
	//			rnd *= deltamax / 4.0;
	//			val += rnd;
	//		} while ((val.real() < minval) || (val.real() > maxval));
	//		res[ncoef] = val;
	//	}
	//	if (data) data->ncoef = ncoefreal;
	//}
	//else {
	//	// Imaginary
	//	double w = 2 * M_PI * input->getCurrentFreq();
	//	minval = ncoef == 0 ? w*minpermint : w*minpermelec;
	//	maxval = ncoef == 0 ? w*maxpermint : w*maxpermelec;

	//	double deltamax = (maxval - minval);
	//	int ncoefimg = 2 * ncoef + 1;

	//	if (sh.shuffleConsts[ncoefimg] == 0) {
	//		res[ncoef] = std::complex<double>(res[ncoef].real(), minval + genreal()*deltamax);
	//	}
	//	else {
	//		std::complex<double> val;
	//		do {
	//			val = res[ncoef];
	//			double rnd = 0;
	//			for (int i = 0; i < sh.shuffleConsts[ncoefimg]; i++)
	//				rnd += genreal();
	//			rnd /= sh.shuffleConsts[ncoefimg];
	//			rnd -= 0.5;
	//			rnd *= deltamax / 4.0;
	//			val += std::complex<double>(0, rnd);
	//		} while ((val.imag() < minval) || (val.imag() > maxval));
	//		res[ncoef] = val;
	//	}
	//	if (data) data->ncoef = ncoefimg;
	//}
	//return res;
}

solutioncomplex *solutioncomplex::shuffle(shuffleData *data, const shuffler &sh) const
{
	std::complex<double> *sigma = getShuffledSolution(data, sh);
	solutioncomplex *res;
	try {
		res = new solutioncomplex(sigma, *this, input, readings, fixedCoeffs);
	}
	catch (...) {
		for (int i = 0; i<65; i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
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
