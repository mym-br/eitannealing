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

// FIXME: IS Col-Major faster than Row-Major?

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> matrix2;

class LB_Solver {
        protected:
                int it;
                
                const matrix &Aii;
                const matrix2 &Aic;
                const SparseIncompleteLLT       &precond;
                
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

                LB_Solver(matrix *Aii, matrix2 *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const SparseIncompleteLLT &precond, double a);
                LB_Solver(matrix *Aii, matrix2 *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const SparseIncompleteLLT &precond, double a, const Eigen::VectorXd &x0);
                void do_iteration();
                
                const Eigen::VectorXd getX() const {
                        return this->x+this->x0;
                }
                
                virtual ~LB_Solver() {};
           
};

class LB_Solver_EG_Estimate : public LB_Solver
{
  typedef Eigen::SparseMatrix<double, Eigen::UpperTriangular> UMatrix;
  
  double ev;
  Eigen::VectorXd evec;
  
  public:
    LB_Solver_EG_Estimate(matrix *Aii, matrix2 *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const SparseIncompleteLLT &precond, int n, float e);
    LB_Solver_EG_Estimate(matrix *Aii, matrix2 *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const SparseIncompleteLLT &precond, const Eigen::VectorXd &x0, const Eigen::VectorXd &egHint, int n, float e); 
    double getLeastEvEst() const {
      return this->ev;
    }
    
    Eigen::VectorXd getEvector() const {
      return this->evec;
    }
};

void assembleProblemElectrodeIdentityMatrix(float *cond, matrix2 **Kic, int numElect);
void assembleProblemMatrix_lb(float *cond, matrix **Kii, matrix2 **Kic, matrix **Kcc,int numElect);

#endif  // SOLVER_LB_H_