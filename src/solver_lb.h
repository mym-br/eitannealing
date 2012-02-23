/*
 * solver.h
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#ifndef SOLVER_LB_H_
#define SOLVER_LB_H_
#endif

#include "solver.h"

// FIXME: IS Col-Major faster than Row-Major?

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> matrix2;


class LB_Solver {
        protected:
                int it;
          
                void init();
                                                
        public:
        
                
                int getIteration() const ;
                double getResidueSquaredNorm() const ;
                
                double getErrorl2Estimate() const ;

                LB_Solver(matrix &Aii, matrix *Aic, Eigen::VectorXd &b, const SparseIncompleteLLT &precond);
                void do_iteration();
           
           
           
};


void assembleProblemMatrix_lb(float *cond, matrix **Kii, matrix2 **Kic, matrix **Kcc,int numElect);