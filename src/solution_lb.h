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


class solution_lb {
	protected:


			double *sol;
			matrix *Aii, *Acc;
			matrix2 *Aic; 
			std::unique_ptr<LB_Solver::Preconditioner> precond;
                        std::shared_ptr<problem> p;
                        // FIXME: o should be const!
                        observations<double> &o;

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

			double *getShufledSolution();
			static double *getNewRandomSolution(int size);
			static double *copySolution(const double *sol, unsigned int size);

			double *getShuffledSolution(shuffleData *data, const shuffler &sh) const;

			
			// shuffle constructor
			solution_lb(double *sol, const solution_lb &base);


	public:

		solution_lb(std::shared_ptr<problem> p, observations<double> &o, const double *sol);
		solution_lb(std::shared_ptr<problem> p, observations<double> &o);	// New random solution
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
		
		double getRegularisationValue() const {
			    return this->regularisation;
			}
		
		void saturate();
		
		void ensureMinIt(unsigned int it);
		
		//void ensureMaxE2(double e2);

		~solution_lb();
		
		
#ifdef __GENERATE_LB_BENCHMARKS
		solution_lb(const float *sol, char benchmarktag);
#endif	// __GENERATE_LB_BENCHMARKS
};


class solution_lb_benchmarked : public solution_lb
{
public: 
  struct benchmark_entry {
      unsigned int timestamp;
      double e_low;
      double e_high;
    } *vector;
    int i;
    int n;
    
    solution_lb_benchmarked(const float *sigma, benchmark_entry *bench, int n);
    
    void performBenchmark();
    
protected:
  static int getTimestamp();
};

#endif /* SOLUTION_H_ */
