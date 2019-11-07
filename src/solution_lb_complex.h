/*
 * solution.h
 *
 *  Created on: Feb 26, 2012
 *      Author: thiago
 */

#ifndef SOLUTION_LB_COMPLEX_H_
#define SOLUTION_LB_COMPLEX_H_

#include "solution.h"
#include "solver_lb.h"
#include "solver_lb_complex.h"
//#include "problemdescription.h"
#include <memory>
#include "observations.h"
#include "problem.h"

class solution_lb_complex {
	protected:


			double *solRe;
			double *solIm;
			matrix *Aii_R, *Aii_I, *Aic_R, *Aic_I, *Acc_R, *Acc_I;
			std::unique_ptr<LB_Solver_Complex::Preconditioner> precond;

			std::shared_ptr<problem> p;
			const observations<std::complex<double> > &o;

			LB_Solver_Complex **simulations;
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
			void initSimulations(const solution_lb_complex &base);
			void initErrors();

			void getShufledSolution(double **s_R, double **s_I);
			double *getNewRandomSolution_R();
			double *getNewRandomSolution_I();
			double *copySolution(const double *sol) const;

			double *getShuffledSolution_R(shuffleData *data, const shuffler &sh) const;
			double *getShuffledSolution_I(shuffleData *data, const shuffler &sh) const;


			// shuffle constructor
			solution_lb_complex(double *sol_R, double *sol_I, const solution_lb_complex &base);


	public:

		solution_lb_complex(std::shared_ptr<problem> p, const observations<std::complex<double> > &o, const double *sol_R, const double *sol_I);
		solution_lb_complex(std::shared_ptr<problem> p, const observations<std::complex<double> > &o);	// New random solution
		bool compareWith(solution_lb_complex &target, float kt, float prob);
		//bool compareWithMinIt(solution &target, float kt, int minit);
		//bool compareWithMaxE2(solution &target, float kt, double e2);
		solution_lb_complex *shuffle(shuffleData *data_R, const shuffler &sh_R, shuffleData *data_I, const shuffler &sh_I, bool *shufle_r) const;

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

		double *getSolution_R() {
			return this->solRe;
		}

		double *getSolution_I() {
			return this->solIm;
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

		~solution_lb_complex();
};

#endif /* SOLUTION_LB_COMPLEX_H_ */
