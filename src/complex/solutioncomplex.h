/*
 * solutioncomplex.h
 *
 *  Created on: Sep 12, 2010
 *      Author: aksato
 */

#ifndef SOLUTIONCOMPLEX_H_
#define SOLUTIONCOMPLEX_H_

class solution;

#include "solvercomplex.h"
//#include "assembleMatrix.h"
//#include "problemdescription.h"
#include "problem.h"
#include <memory>

/*
 * std::autoptr<solution> current, *next;
 * while(kt > xxx) {
 * 	 while(totalit < xxx && acceptit < xxx) {
 *    next = solution->shuffle()
 *
 *    if(current.compareWith(*next, kt, 0.001)) {
 *    		current = next;
 *    		acceptit++;
 *    }
 *    totalit++;
 * 	}
 * 	kt *= 0.95
 * }
 *
 */

struct shuffleData {
	int ncoef;
	bool swap;
};

struct shufflercomplex {
	int * shuffleConsts;
	int * swapshuffleconsts;
	shufflercomplex(std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input) {
		shuffleConsts = new int[input->getNumCoefficients()];
		swapshuffleconsts = new int[input->getInnerAdjacencyCount()];

		for (int i = 0; i<input->getNumCoefficients(); i++) shuffleConsts[i] = 0;
		for (int i = 0; i<input->getInnerAdjacencyCount(); i++) swapshuffleconsts[i] = 0;
	}
	~shufflercomplex() {
		delete[] shuffleConsts;
		delete[] swapshuffleconsts;

	}
	void addShufflerFeedback(const shuffleData &data, bool pos);
};

class solutioncomplex {
	private:


			std::complex<double> *sol;
			matrixcomplex *stiffness, *stiffnessorig;
			SparseIncompleteLLTComplex *precond;


			CG_SolverComplex **simulations;
			Eigen::VectorXd distance;
			Eigen::VectorXd maxdist;
			Eigen::VectorXd mindist;
			Eigen::VectorXd err;
			Eigen::VectorXd err_x_dist;


			double totalDist;
			double minTotalDist;
			double maxTotalDist;
			int critical;
			double critErr;

			int totalit;

			void initSimulations();
			void initSimulations(const solutioncomplex &base);
			void initErrors();
public:
			double *getShufledSolution();
			static std::complex<double> *getNewRandomSolution(std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input);
			static std::complex<double> *copySolution(const std::complex<double> *sol, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input);

			std::complex<double> *getShuffledSolution(shuffleData *data, const shufflercomplex &sh) const;

			static matrixcomplex *getNewStiffness(std::complex<double> *sol, matrixcomplex **stiffnessorig, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input) {
				matrixcomplex *aux = new matrixcomplex;
				input->assembleProblemMatrix(sol, stiffnessorig);
				*aux = (**stiffnessorig).conjugate().selfadjointView<Eigen::Lower>() * (matrixcomplex)(**stiffnessorig).selfadjointView<Eigen::Lower>();
				input->postAssempleProblemMatrix(&aux);
				return aux;
			}

			// shuffle constructor
			solutioncomplex(std::complex<double> *sol, const solutioncomplex &base, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> _input);
			double regularisation;
			std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input;
			void zeroSumVector(Eigen::VectorXcd &vec);

	//public:

		solutioncomplex(const std::complex<double> *sol, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input);
		solutioncomplex(std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> _input);	// New random solution
		bool compareWith(solutioncomplex &target, double kt, double prob);
		bool compareWithMinIt(solutioncomplex &target, double kt, int minit);
		bool compareWithMaxE2(solutioncomplex &target, double kt, double e2);
		solutioncomplex *shuffle(shuffleData *data, const shufflercomplex &sh) const;

		static void saveMesh(double *sol, const char *filename, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input, int step = 0);
		static void savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input);
		static void savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> input);

		double getRegularisationValue() const {
		  return this->regularisation;
		  
		}

		void improve();

		double getDEstimate() const {
			return totalDist;
		}
		double getDMax() const {
			return maxTotalDist;
		}
		double getDMin() const {
			return minTotalDist;
		}

		std::complex<double> *getSolution() {
			return this->sol;
		}

		int getCritical() const {
			return this->critical;
		}

		double getCritErr() const {
			return this->critErr;
		}

		std::complex<double> getErrorAt(int sim) const {
			return this->distance[sim];
		}

		int getTotalIt() const {
			return this->totalit;
		}
		
		void saturate();
		
		void ensureMinIt(unsigned int it);
		
		void ensureMaxE2(double e2);

		~solutioncomplex();

};


#endif /* SOLUTIONCOMPLEX_H_ */
