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
		shuffleConsts = new int[64];
		swapshuffleconsts = new int[2*7*8];

		for(int i=0;i<64;i++) shuffleConsts[i] = 0;
		for(int i=0;i<2*7*8;i++) swapshuffleconsts[i] = 0;
	}
	~shuffler() {
		delete[] shuffleConsts;
		delete[] swapshuffleconsts;

	}
	void addShufflerFeedback(const shuffleData &data, bool pos);
};

class solution {
	private:


			float *sol;
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

			float *getShufledSolution();
			static float *getNewRandomSolution();
			static float *copySolution(const float *sol);

			float *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

			static matrix *getNewStiffness(float *sol) {
				matrix *aux;
				assembleProblemMatrix(sol, &aux);
				return aux;
			}

			// shuffle constructor
			solution(float *sol, const solution &base);


	public:

		solution(const float *sol);
		solution();	// New random solution
		bool compareWith(solution &target, float kt, float prob);
		solution *shuffle(shuffleData *data, const shuffler &sh) const;

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

		float *getSolution() {
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

		~solution();

};


#endif /* SOLUTION_H_ */