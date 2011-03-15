/*
 * solver.h
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "incomplete_cholesky.h"
#include "partial_incomplete_cholesky.h"
#include "nodecoefficients.h"

// FIXME: IS Col-Major faster than Row-Major?
typedef Eigen::SparseMatrix<double, Eigen::ColMajor|Eigen::SelfAdjoint|Eigen::LowerTriangular> matrix;
void assembleProblemMatrix(float *cond, matrix **stiffnes);
void assembleProblemMatrix(float *cond, matrix **stiffnes, int numNodes, nodeCoefficients **nodeCoef);


class CG_Solver {
	protected:
		int it;
		int lastRefresh;
		bool refresh_at_next;
		double rmod, rmod_1;
		
		const matrix &A;
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
		double w[360];
		double wt[360];
		double rt1, r1, r2, r3;
		double err[360];
		//double merr;
		double r0norm;

		double eta_[360];
		double alpha_[360];
		


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

		//matrix buildJacobiMatrx();

};


#endif /* SOLVER_H_ */
