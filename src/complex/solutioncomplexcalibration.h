/*
* solutioncomplex.h
*
*  Created on: Sep 12, 2010
*      Author: aksato
*/

#ifndef SOLUTIONCOMPLEXCALIBRATION_H_
#define SOLUTIONCOMPLEXCALIBRATION_H_

class solution;

#include "solvercomplex.h"
#include "problem.h"
#include <memory>

struct shuffleData {
	int ncoef;
};

struct shufflercomplexcalibration {
	int * shuffleConsts;
	shufflercomplexcalibration(std::shared_ptr<problem> input) {
		shuffleConsts = new int[2 * input->getNumCoefficients()];

		for (int i = 0; i<2 * input->getNumCoefficients(); i++) shuffleConsts[i] = 0;
	}
	~shufflercomplexcalibration() {
		delete[] shuffleConsts;
	}
	void addShufflerFeedback(const shuffleData &data, bool pos);
};

class solutioncomplexcalibration {
private:

	std::complex<double> *sol;
	matrixcomplex *stiffness, *stiffnessorig;
	SparseIncompleteLLTComplex *precond;


	CG_SolverComplex **simulations;
	Eigen::VectorXd distance;
	double totalDist;

	void initSimulations();

	void initSimulations(const solutioncomplexcalibration &base);
	void initErrors();

public:

	double *getShufledSolution();
	static std::complex<double> *getNewRandomSolution(std::shared_ptr<problem> input);
	static std::complex<double> *copySolution(const std::complex<double> *sol, std::shared_ptr<problem> input);

	std::complex<double> *getShuffledSolution(shuffleData *data, const shufflercomplexcalibration &sh) const;

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
	solutioncomplexcalibration(std::complex<double> *sol, const solutioncomplexcalibration &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings);
	std::shared_ptr<problem> input;
	observations<std::complex<double>> *readings;
	void zeroSumVector(Eigen::VectorXcd &vec);

	solutioncomplexcalibration(const std::complex<double> *sol, std::shared_ptr<problem> input, observations<std::complex<double>> *_readings);
	solutioncomplexcalibration(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings);	// New random solution
	bool compareWith(solutioncomplexcalibration &target, double kt, double prob);
	solutioncomplexcalibration *shuffle(shuffleData *data, const shufflercomplexcalibration &sh) const;

	static void saveMesh(double *sol, const char *filename, std::shared_ptr<problem> input, int step = 0);
	static void savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem> input, observations<std::complex<double>> *readings);

	int getTotalIt() { return 32 * 5;  }

	double getDEstimate() const {
		return totalDist;
	}

	double getRegularisationValue() { return 0.0; }

	std::complex<double> *getSolution() {
		return this->sol;
	}

	~solutioncomplexcalibration();
};

const double mincondint = 0.005;
const double maxcondint = 0.380;
const double minpermint = 1.0e-9;
const double maxpermint = 7.09e-8;
const double mincondelec =  3.000;
const double maxcondelec = 10.000;
const double minpermelec = 1.0e-9;
const double maxpermelec = 7.09e-8;

#endif /* SOLUTIONCOMPLEXCALIBRATION_H_ */
