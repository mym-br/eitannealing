/*
* solutioncomplex.h
*
*  Created on: Sep 12, 2010
*      Author: aksato
*/

#ifndef SOLUTIONCOMPLEXCALIBRATION_H_
#define SOLUTIONCOMPLEXCALIBRATION_H_

class solution;

#include <vector>
#include "solvercomplex.h"
#include "problem.h"
#include "solutionbase.h"
#include <memory>

class solutioncomplex : public solutionbase<std::complex<double>> {
protected:
	std::complex<double> *sol;
private:
	matrixcomplex *stiffness, *stiffnessorig;
	SparseIncompleteLLTComplex *precond;


	CG_SolverComplex **simulations;
	void initSimulations();

	void initSimulations(const solutionbase<std::complex<double>> &base);
	void initErrors();

	int fixedCoeffs;
public:
	Eigen::VectorXcd getSimulationX(int i) const { return simulations[i]->getX(); }

	double *getShufledSolution();
	static std::complex<double> *getNewRandomSolution(std::shared_ptr<problem> input, const std::vector<std::complex<double>> &electrodesCoeffs);

	virtual std::complex<double> *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

	static matrixcomplex *getNewStiffness(std::complex<double> *sol, matrixcomplex **stiffnessorig, std::shared_ptr<problem> input) {
		matrixcomplex *aux = new matrixcomplex;
		input->assembleProblemMatrix(sol, stiffnessorig);
		input->addMatrixCapacitances(stiffnessorig);
		matrixcomplex Afull = (**stiffnessorig) + ((matrixcomplex)((matrixcomplex)((**stiffnessorig).selfadjointView<Eigen::Lower>())).triangularView<Eigen::StrictlyUpper>()).conjugate();
		matrixcomplex A_tfull = (**stiffnessorig).conjugate() + ((matrixcomplex)((**stiffnessorig).selfadjointView<Eigen::Lower>())).triangularView<Eigen::StrictlyUpper>();
		*aux = A_tfull * Afull;
		input->postAssembleProblemMatrix<Complex>(&aux);
		return aux;
	}


	// shuffle constructor
	solutioncomplex(std::complex<double> *sol, const solutioncomplex &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, int _fixedCoeffs);
	//std::shared_ptr<problem> input;
	//observations<std::complex<double>> *readings;

	solutioncomplex(const std::complex<double> *sol, std::shared_ptr<problem> input, observations<std::complex<double>> *_readings, int _fixedCoeffs);
	solutioncomplex(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, const std::vector<std::complex<double>> &electrodesCoeffs);	// New random solution
	bool compareWith(solutionbase &target, double kt, double prob);
	virtual solutioncomplex *shuffle(shuffleData *data, const shuffler &sh) const;

	static void savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem> input, observations<std::complex<double>> *readings);

	std::complex<double> *getSolution() { return this->sol; }

	~solutioncomplex();
};


class solutioncomplexcalibration : public solutioncomplex {
public:
	static std::complex<double> *getNewRandomSolution(std::shared_ptr<problem> input);
	std::complex<double> *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

	// shuffle constructor
	solutioncomplexcalibration(std::complex<double> *sol, const solutioncomplexcalibration &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) : solutioncomplex(sol, base, _input, _readings, 0) {};
	solutioncomplexcalibration(const std::complex<double> *sol, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings) : solutioncomplex(sol, _input, _readings, 0) {};
	solutioncomplexcalibration(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings, const std::vector<std::complex<double>> &electrodesCoeffs);	// New random solution
	solutioncomplexcalibration *shuffle(shuffleData *data, const shuffler &sh) const;

	static void saveMesh(double *sol, const char *filename, std::shared_ptr<problem> input, int step = 0);
	//static void savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem> input, observations<std::complex<double>> *readings);

	int getTotalIt() { return 32 * 5; }
	double getRegularisationValue() const { return 0.0; }
	~solutioncomplexcalibration() {}
};

const double mincondint = 0.005;
const double maxcondint = 1.000;
const double minpermint = 1.0e-14;
const double maxpermint = 1.0e-10;
const double mincondelec = 0.380;
//const double maxcondelec = 5000.000;
const double maxcondelec = 11000.000;
const double minpermelec = 1.0e-9;
const double maxpermelec = 1.09e-7;

#endif /* SOLUTIONCOMPLEXCALIBRATION_H_ */
