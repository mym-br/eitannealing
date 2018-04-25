/*
 * solution.h
 *
 *  Created on: Sep 12, 2010
 *      Author: thiago
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_

class solution;

#include "solver.h"
//#include "assembleMatrix.h"
//#include "problemdescription.h"
#include "problem.h"
#include "solutionbase.h"
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
//
//struct shuffleData {
//	int ncoef;
//	bool swap;
//};
//
//struct shuffler {
//	int * shuffleConsts;
//	int * swapshuffleconsts;
//	shuffler(std::shared_ptr<problem> input) {
//		shuffleConsts = new int[input->getNumCoefficients()];
//		swapshuffleconsts = new int[input->getInnerAdjacencyCount()];
//
//		for (int i = 0; i<input->getNumCoefficients(); i++) shuffleConsts[i] = 0;
//		for (int i = 0; i<input->getInnerAdjacencyCount(); i++) swapshuffleconsts[i] = 0;
//	}
//	~shuffler() {
//		delete[] shuffleConsts;
//		delete[] swapshuffleconsts;
//
//	}
//	void addShufflerFeedback(const shuffleData &data, bool pos);
//};
//
//class solution {
//	private:
//
//
//			double *sol;
//			matrix *stiffness;
//			SparseIncompleteLLT *precond;
//
//
//			CG_Solver **simulations;
//			Eigen::VectorXd distance;
//			Eigen::VectorXd maxdist;
//			Eigen::VectorXd mindist;
//			Eigen::VectorXd err;
//			Eigen::VectorXd err_x_dist;
//
//
//			double totalDist;
//			double minTotalDist;
//			double maxTotalDist;
//			int critical;
//			double critErr;
//
//			int totalit;
//
//
//
//
//			void initSimulations();
//			void initSimulations(const solution &base);
//			void initErrors();
//
//			double *getShufledSolution();
//			static double *getNewRandomSolution(std::shared_ptr<problem> input);
//			static double *copySolution(const double *sol, std::shared_ptr<problem> input);
//
//			double *getShuffledSolution(shuffleData *data, const shuffler &sh) const;
//
//			static matrix *getNewStiffness(double *sol, std::shared_ptr<problem> input) {
//				matrix *aux;
//				input->assembleProblemMatrix(sol, &aux);
//				input->postAssembleProblemMatrix(&aux);
//				return aux;
//			}
//
//			// shuffle constructor
//			solution(double *sol, const solution &base, std::shared_ptr<problem> _input, observations<double> *_readings);
//			double regularisation;
//			std::shared_ptr<problem> input;
//			observations<double> *readings;
//			void zeroSumVector(Eigen::VectorXd &vec);
//
//	public:
//
//		solution(const double *sol, std::shared_ptr<problem> input, observations<double> *_readings);
//		solution(std::shared_ptr<problem> _input, observations<double> *_readings);	// New random solution
//		bool compareWith(solution &target, double kt, double prob);
//		bool compareWithMinIt(solution &target, double kt, int minit);
//		bool compareWithMaxE2(solution &target, double kt, double e2);
//		solution *shuffle(shuffleData *data, const shuffler &sh) const;
//
//		static void saveMesh(double *sol, const char *filename, const char *propertyname, std::shared_ptr<problem> input, int step = 0);
//		static void savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings);
//
//		double getRegularisationValue() const {
//		  return this->regularisation;
//		  
//		}
//
//		void improve();
//
//		double getDEstimate() const {
//			return totalDist;
//		}
//		double getDMax() const {
//			return maxTotalDist;
//		}
//		double getDMin() const {
//			return minTotalDist;
//		}
//
//		double *getSolution() {
//			return this->sol;
//		}
//
//		int getCritical() const {
//			return this->critical;
//		}
//
//		double getCritErr() const {
//			return this->critErr;
//		}
//
//		double getErrorAt(int sim) const {
//			return this->distance[sim];
//		}
//
//		int getTotalIt() const {
//			return this->totalit;
//		}
//		
//		void saturate();
//		
//		void ensureMinIt(unsigned int it);
//		
//		void ensureMaxE2(double e2);
//
//		~solution();
//
//};



class solution : public solutionbase<double> {
protected:
	double *sol;
	matrix *stiffness, *stiffnessorig;
	SparseIncompleteLLT *precond;
	CG_Solver **simulations;

	void initSimulations();
	void initSimulations(const solution &base);
	void initErrors();

	int fixedCoeffs;

	double *getShufledSolution();
	static double *getNewRandomSolution(std::shared_ptr<problem> input, std::vector<double> &electrodesCoeffs);
	//static double *copySolution(const double *sol, std::shared_ptr<problem> input);

	double *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

	static matrix *getNewStiffness(double *sol, std::shared_ptr<problem> input) {
		matrix *aux;
		input->assembleProblemMatrix(sol, &aux);
		input->postAssembleProblemMatrix(&aux);
		return aux;
	}

	// shuffle constructor
	solution(double *sol, const solution &base, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs);

public:
	solution(const double *sol, std::shared_ptr<problem> input, observations<double> *_readings, int _fixedCoeffs);
	solution(std::shared_ptr<problem> _input, observations<double> *_readings, std::vector<double> &electrodesCoeffs);	// New random solution
	bool compareWith(solutionbase &target, double kt, double prob);
	bool compareWithMinIt(solutionbase &target, double kt, int minit);
	bool compareWithMaxE2(solutionbase &target, double kt, double e2);
	solution *shuffle(shuffleData *data, const shuffler &sh) const;

	void improve();
	void saturate();
	void ensureMinIt(unsigned int it);
	void ensureMaxE2(double e2);

	Eigen::VectorXd getSimulationX(int i) const { return simulations[i]->getX(); }
	double *getSolution() { return this->sol; }

	static void savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings);
	static void savePotentials(std::vector<Eigen::VectorXf> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings);

	~solution();
};

#endif /* SOLUTION_H_ */
