/*
 * solution.h
 *
 *  Created on: Sep 12, 2010
 *      Author: thiago
 */

#ifndef SOLUTIONCUDA_H_
#define SOLUTIONCUDA_H_

#include "../solution.h"
#include "settings.h"
#include "matrix-cpjds.h"
#include "solvercuda.h"

class solutionCuda : public solutionbase<double> {
//public:
//	solutionCuda(const double *sol, std::shared_ptr<problem> input, observations<double> *_readings, int _fixedCoeffs); //: solution(sol, input, _readings, _fixedCoeffs);
//	solutionCuda(std::shared_ptr<problem> _input, observations<double> *_readings, std::vector<double> &electrodesCoeffs); // New random solution //: solution(_input, _readings, electrodesCoeffs);	// New random solution
//	~solutionCuda() {};
//
//	CGCUDA_Solver **simulationsCjpds;
//private:
//	// shuffle constructor
//	solutionCuda(double *sol, const solutionCuda &base, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs);
//	void initSimulations();
//	void initSimulations(const solutionCuda &base);
//	void initErrors();
//
//	MatrixCPJDS *stiffnessCjpds;
//	MatrixCPJDSManager *mgrCjpds;
//	double lINFinityNorm;
//
//
//
//
//
//
private:
	double *sol;
	matrix *stiffness;
	SparseIncompleteLLT *precond;
	MatrixCPJDS *stiffnessCpjds;
	MatrixCPJDSManager *mgr;
	double lINFinityNorm;
	CGCUDA_Solver **simulations;
//
	void initSimulations();
	void initSimulations(const solutionCuda &base);
	void initErrors();
//
	int fixedCoeffs;
//
//	double *getShufledSolution();
	static double *getNewRandomSolution(std::shared_ptr<problem> input, std::vector<double> &electrodesCoeffs);
//	//static double *copySolution(const double *sol, std::shared_ptr<problem> input);
//
	double *getShuffledSolution(shuffleData *data, const shuffler &sh) const;
//
	static matrix *getNewStiffness(double *sol, std::shared_ptr<problem> input) {
		matrix *aux;
		input->assembleProblemMatrix(sol, &aux);
		input->postAssembleProblemMatrix(&aux);
		return aux;
	}
//
//	// shuffle constructor
	solutionCuda(double *sol, const solutionCuda &base, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs);
//
public:
	solutionCuda(const double *sol, std::shared_ptr<problem> input, observations<double> *_readings, int _fixedCoeffs);
	solutionCuda(std::shared_ptr<problem> _input, observations<double> *_readings, std::vector<double> &electrodesCoeffs);	// New random solution
	bool compareWith(solutionbase &target, double kt, double prob);
	solutionCuda *shuffle(shuffleData *data, const shuffler &sh) const;

	void improve();

	Eigen::VectorXd getSimulationX(int i) const { return (simulations[i]->getX()).cast<double>(); }
	double *getSolution() { return this->sol; }
//
//	static void savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings);
//
	~solutionCuda();
};

#endif /* SOLUTIONCUDA_H_ */
