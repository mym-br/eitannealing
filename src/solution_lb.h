/*
 * solution.h
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#ifndef SOLUTION_LB_H_
#define SOLUTION_LB_H_

#include "solution.h"
#include "solver_lb.h"
#include "problem.h"
#include "lbmatrixbuilder.h"
#include "observations.h"
#include <memory>
#include <vector>


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


template <class solver, class admittance, class observations, class regularizer, class matBuilder, class shuffleData, class shuffler> class solution_lb_gen {
	protected:

			std::vector<admittance> sol;
			std::unique_ptr<typename solver::matrix> Aii, Aic, Acc;
			std::unique_ptr<typename solver::Preconditioner> precond;
			std::shared_ptr<problem> p;
			std::shared_ptr<matBuilder> matrixBuilder;
			std::shared_ptr<regularizer> reg;

			const observations &o;

			std::vector<solver> simulations;
			Eigen::VectorXd distance2;
			Eigen::VectorXd maxdist2;
			Eigen::VectorXd mindist2;
			Eigen::VectorXd err;
			Eigen::VectorXd err_x_dist;

			double totalDist;
			double minTotalDist;
			double maxTotalDist;
			std::vector<unsigned > solver_heap_ordered_by_errors;
			double critErr;
			double regularisation;
			int totalit;

			// Least eigenvector and eigenvalue estimations
			double eigenvalue_estimate;
			typename solver::vector eigenvector_estimate;

			void initSimulations();
			void initSimulations(const solution_lb_gen &base);
			void initErrors();

			std::vector<admittance> getNewRandomSolution(int size);

			std::vector<admittance> getShuffledSolution(shuffleData *data, const shuffler &sh) const;

			// shuffle constructor
			solution_lb_gen(std::vector<admittance> &&sol, const solution_lb_gen &base);

			void initMatrices();


	public:

		solution_lb_gen(std::shared_ptr<problem> p, const observations &o, std::shared_ptr<regularizer> r, std::vector<admittance> &&sol);
		solution_lb_gen(std::shared_ptr<problem> p, const observations &o, std::shared_ptr<regularizer> r);	// New random solution
		bool compareWith(solution_lb_gen &target, float kt, float prob);
		//bool compareWithMinIt(solution &target, float kt, int minit);
		//bool compareWithMaxE2(solution &target, float kt, double e2);
		solution_lb_gen *shuffle(shuffleData *data, const shuffler &sh) const;

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

		const std::vector<admittance> &getSolution() {
			return this->sol;
		}

		int getCritical() const {
			return this->solver_heap_ordered_by_errors[0];
		}

		double getCritErr() const {
			return this->critErr;
		}

		double getErrorAt(int sim) const {
			return std::sqrt(this->distance2[sim]);
		}

		int getTotalIt() const {
			return this->totalit;
		}

		double getRegularisationValue() const {
			    return this->regularisation;
        }

        double getLeastEigenvalueEstimation() const {
            return this->eigenvalue_estimate;
        }

        const Eigen::VectorXd &getLeastEigenvectorEstimation() const {
            return this->eigenvector_estimate;
        }

		void saturate();

		void ensureMinIt(unsigned int it);

		//void ensureMaxE2(double e2);
};

#endif /* SOLUTION_H_ */
