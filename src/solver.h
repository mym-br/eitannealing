/*
 * solver.h
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "basematrix.h"
#include "incomplete_cholesky.h"
#include "circularbuff.h"
//#include "nodecoefficients.h"

class CG_Solver {
	protected:
		int it;
		double totalItTime, totalTriangularTime, totalSpmvTime;
		int lastRefresh;
		bool refresh_at_next;
		double rmod, rmod_1;
		
		matrix::ConstSelfAdjointViewReturnType<Eigen::Lower>::Type A;
		const Eigen::VectorXd &b;
		Eigen::VectorXd x;		
		Eigen::VectorXd p;		
		Eigen::VectorXd r;
		Eigen::VectorXd q;

		// precond!
		const SparseIncompleteLLT	&precond;
		Eigen::VectorXd z;

		
		

		/* Lanczos
		Eigen::VectorXd v;
		Eigen::VectorXd v_1;
		Eigen::VectorXd vt;
		Eigen::VectorXd u;
		double lalpha;
		double leta;*/




		double beta, beta_1, gamma, gamma_1, gamma_2;
		double alpha, eta, eta_p1;
		double c, c_1, s, s_1;
		// Circular buffers
		circularbuff<double,8> w;
		circularbuff<double,8> wt;
		double rt1, r1, r2, r3;
		circularbuff<double,8> err;
		//double merr;
		double r0norm;
		double r0norm2;
		

		circularbuff<double,8> eta_;
		circularbuff<double,8> alpha_;
		


		void init(double res = -1);
						
		void saveVals(const char* fname, double val, bool app = true);

	public:
	
		void setrefresh() {
			this->refresh_at_next = true;			
		}
		
		int getIteration() const {
			return this->it;
		}
		
		std::tuple<double, double, double> getIterationTimes() {
			return{ totalItTime / (double)it, totalTriangularTime / (double)it, totalSpmvTime / (double)it };
		}

		double getResidueSquaredNorm() const {
			return this->rmod;
		}
		
		double getErrorl2Estimate() const {
			if(it>2) {
				return r0norm*r0norm*(this->err[it-1])*precond.getLINFinityNorm();
			} else return 0;
		}

		CG_Solver(matrix &A, Eigen::VectorXd &b, const SparseIncompleteLLT &precond, double res);
		CG_Solver(matrix &A, Eigen::VectorXd &b, const SparseIncompleteLLT &precond);
		CG_Solver(matrix &A, Eigen::VectorXd &b, const Eigen::VectorXd &x0, const SparseIncompleteLLT &precond);
		void do_iteration();
		
		const Eigen::VectorXd &getX() const {
			return this->x;
		}
		
		Eigen::VectorXd getY() const {
			Eigen::VectorXd y(x);
			this->precond.halfMultInPlace(y);
			return y;
		}

		const Eigen::VectorXd &getR() const {
			return this->r;
		}
		
		const SparseIncompleteLLT	&getPrecond() {
			return this->precond;
		}
		
		double getLastE2() const {
		  return rmod/r0norm2;
		}

		//matrix buildJacobiMatrx();

};


#endif /* SOLVER_H_ */
