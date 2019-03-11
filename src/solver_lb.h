/*
 * solver.h
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#ifndef SOLVER_LB_H_
#define SOLVER_LB_H_

#include "solver.h"
#include <memory>
#include "problem.h"

// FIXME: IS Col-Major faster than Row-Major?

class LB_Solver {
	public:
		//typedef SparseIncompleteLQ Preconditioner;
		typedef SparseIncompleteLLT Preconditioner;
	protected:
                int it;

				matrix::ConstSelfAdjointViewReturnType<Eigen::Lower>::Type Aii;
                const matrix &Aic;
                const Preconditioner       &precond;

                double JhatNorm2;
                double ATJhatNorm2;

				double g;
                double pi, phi2, phi;
                double c, gamma_ip, delta;
				double si, psi_im;
                Eigen::VectorXd p, pc;
                Eigen::VectorXd r, rc;
                Eigen::VectorXd s, q, qaux;

				bool lowerSafe;
				double g_im;
				double pi_im;
				double a;
				double dt, at;
				double phi2t;
				double gr;
				double alpha, beta;

                Eigen::VectorXd x0, x, w;
                float fit;
                void init();

        public:


                int getIteration() const {
                  return it;
                }

                double getErrorl2Estimate() const ;
				double getMinErrorl2Estimate() const ;
				double getMaxErrorl2Estimate() const {
					return getErrorl2Estimate();
				}

                LB_Solver(matrix *Aii, matrix *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, double a);
                LB_Solver(matrix *Aii, matrix *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, double a, const Eigen::VectorXd &x0);
                void do_iteration();

                const Eigen::VectorXd getX() const {
			 Eigen::VectorXd xaux(this->x);
			 precond.solveInPlace(xaux);
                        return xaux+this->x0;
                }

                virtual ~LB_Solver() {};
		static Preconditioner *makePreconditioner(const matrix &A) {
		    return new Preconditioner(A);//,15,6);
		}

};

class LB_Solver_EG_Estimate : public LB_Solver
{
  double ev;
  Eigen::VectorXd evec;

  public:
    LB_Solver_EG_Estimate(matrix *Aii, matrix *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, int n, float e);
    LB_Solver_EG_Estimate(matrix *Aii, matrix *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, const Eigen::VectorXd &x0, const Eigen::VectorXd &egHint, int n, float e);
    double getLeastEvEst() const {
      return this->ev;
    }

    Eigen::VectorXd getEvector() const {
      return this->evec;
    }
};

void assembleProblemMatrix_lb(double *cond, matrix **Kii, matrix2 **Kic, matrix **Kcc, problem &p);

#endif  // SOLVER_LB_H_
