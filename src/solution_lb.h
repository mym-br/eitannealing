/*
 * solution.h
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#ifndef SOLUTION_LB_H_
#define SOLUTION_LB_H_

class solution_lb;

#include "solution.h"
#include "solver_lb.h"
#include "problemdescription.h"
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


class solution_lb {
	private:


			float *sol;
			matrix *Aii, *Acc;
			matrix2 *Aic; 
			std::auto_ptr<LB_Solver::Preconditioner> precond;

			LB_Solver **simulations;
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
			double regularisation;
			int totalit;




			void initSimulations();
			void initSimulations(const solution_lb &base);
			void initErrors();

			float *getShufledSolution();
			static float *getNewRandomSolution();
			static float *copySolution(const float *sol);

			float *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

			
			// shuffle constructor
			solution_lb(float *sol, const solution_lb &base);


	public:

		solution_lb(const float *sol);
		solution_lb();	// New random solution
		bool compareWith(solution_lb &target, float kt, float prob);
		//bool compareWithMinIt(solution &target, float kt, int minit);
		//bool compareWithMaxE2(solution &target, float kt, double e2);
		solution_lb *shuffle(shuffleData *data, const shuffler &sh) const;

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
		
		double getRegularisationValue() const {
			    return this->regularisation;
			}
		
		void saturate();
		
		void ensureMinIt(unsigned int it);
		
		//void ensureMaxE2(double e2);

		~solution_lb();
};


#endif /* SOLUTION_H_ */
