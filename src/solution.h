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
#include "assembleMatrix.h"
#include "problemdescription.h"


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

struct shuffler {
	int * shuffleConsts;
	int * swapshuffleconsts;
	shuffler() {
		shuffleConsts = new int[numcoefficients];
		swapshuffleconsts = new int[innerAdjacency.size()];

		for(int i=0;i<numcoefficients;i++) shuffleConsts[i] = 0;
		for(int i=0;i<innerAdjacency.size();i++) swapshuffleconsts[i] = 0;
	}
	~shuffler() {
		delete[] shuffleConsts;
		delete[] swapshuffleconsts;

	}
	void addShufflerFeedback(const shuffleData &data, bool pos);
};

class solution {
	private:


			double *sol;
			matrix *stiffness;
			SparseIncompleteLLT *precond;


			CG_Solver **simulations;
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
			void initSimulations(const solution &base);
			void initErrors();

			double *getShufledSolution();
			static double *getNewRandomSolution();
			static double *copySolution(const double *sol);

			double *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

			static matrix *getNewStiffness(double *sol) {
				matrix *aux;
				assembleProblemMatrix(sol, &aux);
				return aux;
			}

			// shuffle constructor
			solution(double *sol, const solution &base);
			double regularisation;


	public:

		solution(const double *sol);
		solution();	// New random solution
		bool compareWith(solution &target, double kt, double prob);
		bool compareWithMinIt(solution &target, double kt, int minit);
		bool compareWithMaxE2(solution &target, double kt, double e2);
		solution *shuffle(shuffleData *data, const shuffler &sh) const;
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

		double *getSolution() {
			return this->sol;
		}

		int getCritical() const {
			return this->critical;
		}

		double getCritErr() const {
			return this->critErr;
		}

		double getErrorAt(int sim) const {
			return this->distance[sim];
		}

		int getTotalIt() const {
			return this->totalit;
		}
		
		void saturate();
		
		void ensureMinIt(unsigned int it);
		
		void ensureMaxE2(double e2);

		~solution();

};


#endif /* SOLUTION_H_ */
