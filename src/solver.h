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
//#include "nodecoefficients.h"

template<class base, size_t len> class circularbuff
{
	base _val[len];
	public:

		typedef base basetype;

		const base &operator[](int i) const {
			if(i>=0)
				return _val[i%len];
			return _val[(len - ((-i)%len))%len];
		}

		base &operator[](int i) {
			if(i>=0)
				return _val[i%len];
			return _val[(len - ((-i)%len))%len];
		}

		circularbuff(){};
};

class CG_Solver {
	protected:
		int it;
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
		


		void init();
						
	public:
	
		void setrefresh() {
			this->refresh_at_next = true;			
		}
		
		int getIteration() const {
			return this->it;
		}
		
		double getResidueSquaredNorm() const {
			return this->rmod;
		}
		
		double getErrorl2Estimate() const {
			if(it>2) {
				return r0norm*r0norm*(this->err[it-1])*precond.getLINFinityNorm();
			} else return 0;
		}

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
