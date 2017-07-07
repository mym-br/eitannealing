/*
* solution.cpp
*
*  Created on: Sep 12, 2010
*      Author: thiago
*/

#include "solutioncomplexcalibration.h"
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

void solutioncomplexcalibration::zeroSumVector(Eigen::VectorXcd &vec) {
	std::complex<double> avg = 0;
	for (int i = 0; i < vec.size(); i++) avg += vec[i];
	avg /= input->getGenericElectrodesCount();
	for (int i = 0; i < vec.size(); i++) vec[i] -= avg;
}

bool solutioncomplexcalibration::compareWith(solutioncomplexcalibration &target, double kt, double prob)
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

solutioncomplexcalibration::solutioncomplexcalibration(const std::complex<double> *sigma, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
sol(solutioncomplexcalibration::copySolution(sigma, _input)),
stiffness(solutioncomplexcalibration::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
simulations(new CG_SolverComplex *[_readings->getNObs()]),
distance(_readings->getNObs()),
input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}


// New random solution
solutioncomplexcalibration::solutioncomplexcalibration(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
sol(solutioncomplexcalibration::getNewRandomSolution(_input)),
stiffness(solutioncomplexcalibration::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
simulations(new CG_SolverComplex *[_readings->getNObs()]),
distance(_readings->getNObs()),
input(_input), readings(_readings)
{
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solutioncomplexcalibration::solutioncomplexcalibration(std::complex<double> *sigma, const solutioncomplexcalibration &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) :
sol(sigma),
stiffness(solutioncomplexcalibration::getNewStiffness(sol, &stiffnessorig, _input)),
precond(new SparseIncompleteLLTComplex(*stiffness)),
simulations(new CG_SolverComplex *[_readings->getNObs()]),
distance(_readings->getNObs()),
input(_input), readings(_readings)
{
	this->initSimulations(base);
	this->initErrors();
}


void solutioncomplexcalibration::initSimulations(const solutioncomplexcalibration &base)
{
	// Prepare solvers
	for (int i = 0; i<readings->getNObs(); i++)
	{
		// Reuse previous solutions as initial values
		simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), base.simulations[i]->getX(), *precond);
	}
}

void solutioncomplexcalibration::initSimulations()
{
	// Prepare solvers
	for (int i = 0; i<readings->getNObs(); i++)
	{
		simulations[i] = new CG_SolverComplex(*stiffness, input->getConjugatedCurrentVector(i, stiffnessorig, readings), *precond);
		for (int k = 0; k < 5; k++) simulations[i]->do_iteration();
	}
}

void solutioncomplexcalibration::initErrors()
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


std::complex<double> *solutioncomplexcalibration::copySolution(const std::complex<double> *sol, std::shared_ptr<problem> input)
{
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];

	for (int i = 0; i<input->getNumCoefficients(); i++)
		res[i] = sol[i];

	return res;
}


std::complex<double> *solutioncomplexcalibration::getNewRandomSolution(std::shared_ptr<problem> input)
{
	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
	double w = 2 * M_PI * input->getCurrentFreq();
	double wminperm = w*minperm, wmaxperm = w*maxperm;
	for (int i = 0; i < input->getNumCoefficients(); i++)
		res[i] = std::complex<double>(mincond + genreal()*(maxcond - mincond), wminperm + genreal()*(wmaxperm - wminperm));

	return res;
}

void solutioncomplexcalibration::saveMesh(double *sol, const char *filename, std::shared_ptr<problem> input, int step) {
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

void solutioncomplexcalibration::savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem> input, observations<std::complex<double>> *readings) {
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

std::complex<double> *solutioncomplexcalibration::getShuffledSolution(shuffleData *data, const shufflercomplexcalibration &sh) const
{
	std::complex<double> *res = solutioncomplexcalibration::copySolution(sol, input);
	// Real or complex shuffle
	if (genint(2)) {
		// Real
		int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;
		double deltamax = (maxcond - mincond);
		int ncoefreal = 2 * ncoef;

		if (sh.shuffleConsts[ncoefreal] == 0) {
			res[ncoef] = std::complex<double>(mincond + genreal()*deltamax, res[ncoef].imag());
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
			} while ((val.real() < mincond) || (val.real() > maxcond));
			res[ncoef] = val;
		}
		if (data) {
			data->ncoef = ncoefreal;
		}
	}
	else {
		// Imaginary
		double w = 2 * M_PI * input->getCurrentFreq();
		double wminperm = w*minperm, wmaxperm = w*maxperm;

		int ncoef = genint(input->getNumCoefficients());	// Lower values fixed;
		double deltamax = (wmaxperm - wminperm);
		int ncoefimg = 2 * ncoef + 1;

		if (sh.shuffleConsts[ncoefimg] == 0) {
			res[ncoef] = std::complex<double>(res[ncoef].real(), wminperm + genreal()*deltamax);
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
			} while ((val.imag() < wminperm) || (val.imag() > wmaxperm));
			res[ncoef] = val;
		}
		if (data) {
			data->ncoef = ncoefimg;
		}
	}
	return res;
}

void shufflercomplexcalibration::addShufflerFeedback(const shuffleData &data, bool pos)
{
	if (pos) { // positive feedback
		this->shuffleConsts[data.ncoef] /= 2;
	}
	else {
		this->shuffleConsts[data.ncoef]++;
	}
}

solutioncomplexcalibration *solutioncomplexcalibration::shuffle(shuffleData *data, const shufflercomplexcalibration &sh) const
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

solutioncomplexcalibration::~solutioncomplexcalibration()
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

