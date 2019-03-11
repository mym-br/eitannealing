/*
 * solver.h
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#ifndef SOLVER_LB_COMPLEX_H_
#define SOLVER_LB_COMPLEX_H_

#include <complex>
#include "incomplete_qr_complex.h"
#include <memory>

class LB_Solver_Complex {
	public:
		//typedef SparseIncompleteLQ Preconditioner;
		typedef SparseIncompleteQRComplex Preconditioner;
	protected:
        int it;

				matrix::ConstSelfAdjointViewReturnType<Eigen::Lower>::Type Aii_R, Aii_I;
        const matrix &Aic_R, &Aic_I;
        const Preconditioner       &precond;

        double JhatNorm2;
        double ATJhatNorm2;

				double g;
        double pi, phi2, phi;
        double c, gamma_ip, delta;
				double si, psi_im;
        Eigen::VectorXd p_R, p_I, pc_R, pc_I;
        Eigen::VectorXd r_R, r_I, rc_R, rc_I;
        Eigen::VectorXd s_R, s_I, q_R, q_I, qaux_R, qaux_I;

				bool lowerSafe;
				double g_im;
				double pi_im;
				double a;
				double dt, at;
				double phi2t;
				double gr;
				double alpha, beta;

        Eigen::VectorXd x0_R, x0_I, x_R, x_I, w_R, w_I;
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

        LB_Solver_Complex(matrix *Aii_R, matrix *Aii_I, matrix *Aic_R, matrix *Aic_I, matrix *Acc_R, matrix *Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, double a);
        LB_Solver_Complex(matrix *Aii_R, matrix *Aii_I, matrix *Aic_R, matrix *Aic_I, matrix *Acc_R, matrix *Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, double a, const Eigen::VectorXd &x0_R, const Eigen::VectorXd &x0_I);
        void do_iteration();

        void getX(Eigen::VectorXd &R, Eigen::VectorXd &I) const {
			 			R = this->x_R;
						I = this->x_I;
			 			precond.solveInPlaceC(R, I);
        }

        virtual ~LB_Solver_Complex() {};
				static Preconditioner *makePreconditioner(const matrix &Aii_R, const matrix &Aii_I, const matrix &Aic_R, const matrix &Aic_I) {
		    	return new Preconditioner(8, 16, Aii_R, Aii_I, Aic_R, Aic_I);
				}
};

class LB_Solver_EG_Complex_Estimate : public LB_Solver_Complex
{
  double ev;
  Eigen::VectorXd evec;

  public:
    LB_Solver_EG_Complex_Estimate(matrix *Aii_R, matrix *Aii_I, matrix *Aic_R, matrix *Aic_I, matrix *Acc_R, matrix *Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, int n, float e);
    LB_Solver_EG_Complex_Estimate(matrix *Aii_R, matrix *Aii_I, matrix *Aic_R, matrix *Aic_I, matrix *Acc_R, matrix *Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, const Eigen::VectorXd &x0_R, const Eigen::VectorXd &x0_I, const Eigen::VectorXd &egHint, int n, float e);
    double getLeastEvEst() const {
      return this->ev;
    }

    Eigen::VectorXd getEvector() const {
      return this->evec;
    }
};

void assembleProblemMatrix_lb(double *cond, matrix **Kii, matrix **Kic, matrix **Kcc,int numElect);

#endif  // SOLVER_LB_COMPLEX_H_
