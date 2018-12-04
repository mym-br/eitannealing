/*
 * solutioncomplexcomplete.h
 *
 *  Created on: Sep 12, 2010
 *      Author: aksato
 */

#ifndef solutioncomplexcompleteCOMPLETE_H_
#define solutioncomplexcompleteCOMPLETE_H_

class solution;

#include "solvercomplex.h"
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

struct shuffleDataComplete {
	int ncoef;
	bool swap;
};

struct shufflercomplexcomplete {
	int * shuffleConsts;
	int * swapshuffleconsts;
	shufflercomplexcomplete(std::shared_ptr<problem> input) {
		shuffleConsts = new int[input->getNumCoefficients()];
		swapshuffleconsts = new int[input->getInnerAdjacencyCount()];

		for (int i = 0; i<input->getNumCoefficients(); i++) shuffleConsts[i] = 0;
		for (int i = 0; i<input->getInnerAdjacencyCount(); i++) swapshuffleconsts[i] = 0;
	}
	~shufflercomplexcomplete() {
		delete[] shuffleConsts;
		delete[] swapshuffleconsts;

	}
	void addShufflerFeedback(const shuffleDataComplete &data, bool pos);
};

class solutioncomplexcomplete {
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
			void initSimulations(const solutioncomplexcomplete &base);
			void initErrors();
public:
			double *getShufledSolution();
			static std::complex<double> *getNewRandomSolution(std::shared_ptr<problem> input);
			static std::complex<double> *copySolution(const std::complex<double> *sol, std::shared_ptr<problem> input);

			std::complex<double> *getShuffledSolution(shuffleDataComplete *data, const shufflercomplexcomplete &sh) const;
			double calcNewShuffledValue(int ncoef, double curVal, shuffleDataComplete *data, const shufflercomplexcomplete &sh, double minval, double maxval) const;
			std::pair<double, double> calcNewSwappedValue(int ncoef, int &node1, int &node2, double v1, double v2, shuffleDataComplete *data, const shufflercomplexcomplete &sh, double minval, double maxval) const;

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
			solutioncomplexcomplete(std::complex<double> *sol, const solutioncomplexcomplete &base, std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings);
			double regularisation;
			std::shared_ptr<problem> input;
			observations<std::complex<double>> *readings;
			void zeroSumVector(Eigen::VectorXcd &vec);

	//public:

		solutioncomplexcomplete(const std::complex<double> *sol, std::shared_ptr<problem> input, observations<std::complex<double>> *_readings);
		solutioncomplexcomplete(std::shared_ptr<problem> _input, observations<std::complex<double>> *_readings);	// New random solution
		bool compareWith(solutioncomplexcomplete &target, double kt, double prob);
		bool compareWithMinIt(solutioncomplexcomplete &target, double kt, int minit);
		bool compareWithMaxE2(solutioncomplexcomplete &target, double kt, double e2);
		solutioncomplexcomplete *shuffle(shuffleDataComplete *data, const shufflercomplexcomplete &sh) const;

		static void saveMesh(double *sol, const char *filename, std::shared_ptr<problem> input, int step = 0);
		static void savePotentials(std::vector<Eigen::VectorXcd> &sols, const char *filename, std::shared_ptr<problem> input, observations<std::complex<double>> *readings);

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

		~solutioncomplexcomplete();

};


#endif /* solutioncomplexcompleteCOMPLETE_H_ */
