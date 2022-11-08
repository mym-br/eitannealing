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

template<class num_engine> class LB_Solver_A {
	public:
		//typedef SparseIncompleteLQ Preconditioner;
		typedef typename num_engine::preconditioner Preconditioner;
		typedef typename num_engine::vector vector;
		typedef typename num_engine::symMatrix symMatrix;
		typedef typename num_engine::matrix matrix;
		typedef typename num_engine::scalar scalar;
		typedef typename num_engine::real real;
	protected:
            unsigned int it;

            const symMatrix &Aii;
						//Eigen::SparseMatrix<double, Eigen::ColMajor>::ConstSelfAdjointViewReturnType<Eigen::Lower>::Type Aii;
            const matrix &Aic;
            const Preconditioner &precond;

            real JhatNorm2;
            real ATJhatNorm2;

			real g;
            real pi, phi2, phi;
            real c, gamma_ip, delta;
			real si, psi_im;
            vector p, pc;
            vector r, rc;
            vector s, q, qaux;

			bool lowerSafe;
			real g_im;
			real pi_im;
			real a;
			real dt, at;
			real phi2t;
			real gr;
			real alpha, beta;

            vector x0, x, w;
			real fit;
            void init(){
			    JhatNorm2 = r.squaredNorm()+rc.squaredNorm();
			    delta = sqrt(JhatNorm2);
			    p = r/delta; pc = rc/delta;
                num_engine::conjtranspose_product_ii_ic_vector(s, Aii, Aic, p, pc);
				precond.solveInPlaceT(s);
				ATJhatNorm2 = s.squaredNorm();
				gamma_ip = sqrt(ATJhatNorm2); // gamma of *NEXT* iteration is obtained here!
			    ATJhatNorm2*=JhatNorm2;

			    q = s/gamma_ip;              // uses gamma of *NEXT* iteration
				// *** Gauss
			    g=0;

			    // 1
			    w = q;
			    fit = delta;

				qaux = q;
			    precond.solveInPlace(qaux); // q  <- q*C^-1
			    num_engine::product_ii_ic_vector(r, rc, Aii, Aic, qaux);
				r -= gamma_ip*p; rc -= gamma_ip*pc;
			    delta = sqrt(r.squaredNorm()+rc.squaredNorm());
				p = r/delta; pc = rc/delta;
			    num_engine::conjtranspose_product_ii_ic_vector(s, Aii, Aic, p, pc);
				precond.solveInPlaceT(s);
				s -= delta*q;
			    // *** Gauss, as next value for gamma will be pertinent to iteration 2!
			    phi2 = gamma_ip*gamma_ip+delta*delta;
				phi = sqrt(phi2);
				c = -gamma_ip/phi;
				si = delta/phi;
				pi = 1/phi2;
				g+=pi;
				// This is due to gauss-radau
				alpha = gamma_ip*gamma_ip+delta*delta;
				gamma_ip = s.norm();
				q = s/gamma_ip;

				// Gauss-radau
				beta = gamma_ip*delta;
				at = a;
				dt = alpha - a;

			    x = - (c*fit/phi)*w;

			    it = 1;
			}

        public:


          unsigned int getIteration() const {
            return it;
          }

          real getErrorl2Estimate() const {
				    return sqrt(JhatNorm2 - ATJhatNorm2*g);
					}

				real getMinErrorl2Estimate() const {
				    if(!lowerSafe) return 0;
				        real v = JhatNorm2 - ATJhatNorm2*(gr);
						if(v<0) return 0;//this->getX().norm()*0.0005;
						return sqrt(v);//+this->getX().norm()*0.0005;
				}

				real getMaxErrorl2Estimate() const {
					return getErrorl2Estimate();
				}

	  /*
        // std::vector requires a copy constructor for resizes.
        //  This is just a non-working placeholder, call reserve() so resize is never used.
        /*LB_Solver_A ( const LB_Solver_A &s):
         Aii(s.Aii), Aic(s.Aic), precond(s.precond), lowerSafe(true), a(s.a), x0(vector::Zero(s.Aii.rows())){
            std::cerr << "LB_Solver_A copy constructor invoked. This should never happen (wrong container size?)\n";
            std::exit(-1);
        }

        LB_Solver_A& operator=(LB_Solver_A const&) = delete;*/
// FIXME: Why pointers on the constructors???
        // Specialized constructor for unknown initial guess.
	  //	Uses zero as an initial guess and skips the first subtraction on r
	  LB_Solver_A(symMatrix *_Aii, matrix *_Aic, symMatrix *_Acc, const vector &J, const vector &Phi, const Preconditioner &precond, real a):
				    Aii(*_Aii), Aic(*_Aic), precond(precond), lowerSafe(true), a(a), x0(vector::Zero(_Aii->rows()))
				{
				    it = 0;
				    // 0
				    num_engine::j_minus_a_x_phi(r, rc, Aic, *_Acc, J, Phi);
				    init();
				}

        LB_Solver_A(symMatrix *_Aii, matrix *_Aic, symMatrix *_Acc, const vector &J, const vector &Phi, const Preconditioner &precond, real a, const vector &x0):
				    Aii(*_Aii), Aic(*_Aic), precond(precond), lowerSafe(true), a(a), x0(x0)
				{
					it = 0;
					// 0
					num_engine::j_minus_a_x_phi(r, rc, Aic, *_Acc, J, Phi);
					Eigen::VectorXd xaux(x0);
					num_engine::subtract_a_x(r, rc, Aii, Aic, x0);
					init();
				}



        void do_iteration() {
				    qaux = q;
				    precond.solveInPlace(qaux); // q  <- q*C^-1
                    num_engine::product_ii_ic_vector(r, rc, Aii, Aic, qaux);
					r -= gamma_ip*p; rc -= gamma_ip*pc;
				    delta = sqrt(r.squaredNorm()+rc.squaredNorm());
					p = r/delta; pc = rc/delta;
				    num_engine::conjtranspose_product_ii_ic_vector(s, Aii, Aic, p, pc);
				    precond.solveInPlaceT(s);
					s -= delta*q;
                    // *** Gauss, as next value for gamma will be pertinent to next iteration!
                    psi_im = si*gamma_ip;

                    fit *= si;
                    w = q - (psi_im/phi)*w;
                    phi2 = c*c*gamma_ip*gamma_ip+delta*delta;
                    phi = sqrt(phi2);
                    c *= -gamma_ip/phi;
                    si = delta/phi;
                    pi_im = pi;
                    pi *= psi_im*psi_im/phi2;
                    g_im = g;
                    g+=pi;
                    // This is due to gauss-radau
                    alpha = gamma_ip*gamma_ip+delta*delta;

				    gamma_ip = s.norm();
					q = s/gamma_ip;

					// Gauss-radau
					at = a + beta*beta/dt;
					phi2t = at - psi_im*psi_im;
					dt = alpha - a - (beta*beta/dt);
					beta = gamma_ip*delta;

					gr = g_im + pi_im*(psi_im*psi_im/phi2t);

				     x -= (c*fit/phi)*w;
				     it++;
				}

                vector getX() const {
			 		vector xaux(this->x);
			 		precond.solveInPlace(xaux);
                        return xaux+this->x0;
                }

                //virtual ~LB_Solver_A() {};
				static Preconditioner *makePreconditioner(const symMatrix &A_ii, const matrix &A_ic) {
		    		return num_engine::make_new_preconditioner(A_ii, A_ic);
				}

    private:
        // This method only works if called *after* init and *before* any do_iteration!
        std::pair<scalar, vector> estimate_smallest_eigenvalue(int n, float e, vector &&temp_eigenvector) {
            Eigen::SparseMatrix<double, Eigen::ColMajor>  U(n,n);
            // Alpha and beta must be stored in order to recalculate dt
            std::vector<double> AlphaVector;
            std::vector<double> BetaVector;
            AlphaVector.push_back(this->alpha);
            BetaVector.push_back(this->beta);
            U.reserve(2*n-1);
            U.insert(0,0) = this->phi;
            for(int i=1;i<n;i++) {
                this->do_iteration();
                AlphaVector.push_back(this->alpha);
                BetaVector.push_back(this->beta);
                U.insert(i-1,i) = this->psi_im;
                U.insert(i,i) = this->phi;
            }
            U.makeCompressed();
            // Now calc eigenvalue
            scalar oev=0;
            scalar ev = 1.0;
            while(fabs(oev-ev)/ev > e) {
                U.triangularView<Eigen::Upper>().transpose().solveInPlace(temp_eigenvector);
                U.triangularView<Eigen::Upper>().solveInPlace(temp_eigenvector);
                oev = ev;
                ev = 1/temp_eigenvector.norm();
                temp_eigenvector *= ev;
            }
            this->a = ev;
            // We need now to recalculate Dt...
            this->dt = AlphaVector[0] - this->a;
            for(unsigned int i=1;i<this->it;i++) {
                this->dt = AlphaVector[i] - this->a - (BetaVector[i-1]*BetaVector[i-1])/this->dt;
	        }
            // Update Gauss-Radau values
            this->do_iteration();
            return std::make_pair(ev, std::move(temp_eigenvector));
        }
    public:
        // Constructor with eigevalue estimation
	  LB_Solver_A(symMatrix *Aii, matrix *Aic, symMatrix *Acc, const vector &J, const vector &Phi, const Preconditioner &precond,
                                int n, float e, scalar &res_smallest_eigenvalue, vector &res_eigenvector):
                                LB_Solver_A<num_engine>(Aii, Aic, Acc, J, Phi, precond, 0)
        {
            auto res = estimate_smallest_eigenvalue(n, e, Eigen::VectorXd::Constant(n,1/sqrt(n)));
            res_smallest_eigenvalue = res.first;
            res_eigenvector = std::move(res.second);
        }

        // Constructor with eigevalue estimation reusing previous eigenvector estimate
        LB_Solver_A(symMatrix *Aii, matrix *Aic, symMatrix *Acc, const vector &J, const vector &Phi, const Preconditioner &precond,
                    const vector &x0, const vector &egHint, unsigned int n, float e, scalar &res_smallest_eigenvalue, vector &res_eigenvector):
                    LB_Solver_A<num_engine>(Aii, Aic, Acc, J, Phi, precond, 0, x0)
        {
            auto res = estimate_smallest_eigenvalue(n, e, vector(egHint));
            res_smallest_eigenvalue = res.first;
            res_eigenvector = std::move(res.second);
        }
};

void assembleProblemMatrix_lb(double *cond, matrix **Kii, matrix **Kic, matrix **Kcc, problem &p);

struct eigen_double_engine {
	typedef double scalar;
	typedef double real;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> symMatrix;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> matrix;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector;
	typedef SparseIncompleteLLT preconditioner;

	inline static void product_ii_ic_vector(vector &dest_i, vector &dest_c, const symMatrix &ii, const symMatrix &ic, const vector &b) {
		dest_i.noalias() = ii.selfadjointView<Eigen::Lower>()*b;
		dest_c.noalias() = ic*b;
	}

	inline static void conjtranspose_product_ii_ic_vector(vector &dest, const symMatrix &ii, const symMatrix &ic, const vector &b_i, const vector &b_c) {
		dest.noalias() = ii.selfadjointView<Eigen::Lower>()*b_i;
		dest.noalias() += ic.transpose()*b_c;
	}

	inline static void j_minus_a_x_phi(vector &dest_i, vector &dest_c, const matrix &ic, const symMatrix &cc, const vector &j, const vector &phi) {
		dest_i.noalias() = -ic.transpose()*phi;
		dest_c.noalias() = j;
		dest_c.noalias() -= cc.selfadjointView<Eigen::Lower>()*phi;
	}

	inline static void subtract_a_x(vector &dest_i, vector &dest_c, const symMatrix &ii,  const symMatrix &ic, const vector &x) {
		dest_i.noalias() -= ii.selfadjointView<Eigen::Lower>()*x;
		dest_c.noalias() -= ic*x;
	}

	inline static preconditioner *make_new_preconditioner(const symMatrix &A_ii, const matrix &A_ic) {
		return new preconditioner(A_ii);
	}
};

typedef LB_Solver_A<eigen_double_engine> LB_Solver;
//typedef LB_Solver_EG_Estimate_A<eigen_double_engine> LB_Solver_EG_Estimate;


#endif  // SOLVER_LB_H_
