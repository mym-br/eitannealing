/*
 * solver.h
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#ifndef SOLVERCOMPLEX_H_
#define SOLVERCOMPLEX_H_

#include "basematrix.h"
#include "incomplete_choleskycomplex.h"
#include "solver.h"
//#include "nodecoefficients.h"

//template<class base, size_t len> class circularbuff
//{
//	base _val[len];
//	public:
//
//		typedef base basetype;
//
//		const base &operator[](int i) const {
//			if(i>=0)
//				return _val[i%len];
//			return _val[(len - ((-i)%len))%len];
//		}
//
//		base &operator[](int i) {
//			if(i>=0)
//				return _val[i%len];
//			return _val[(len - ((-i)%len))%len];
//		}
//
//		circularbuff(){};
//};

class CG_SolverComplex {
	protected:
		int it;
		int lastRefresh;
		bool refresh_at_next;
		std::complex<double> rmod, rmod_1;
		
		matrixcomplex::ConstSelfAdjointViewReturnType<Eigen::Lower>::Type A;
		const Eigen::VectorXcd &b;
		Eigen::VectorXcd x;		
		Eigen::VectorXcd p;		
		Eigen::VectorXcd r;
		Eigen::VectorXcd q;

		// precond!
		const SparseIncompleteLLTComplex &precond;
		Eigen::VectorXcd z;

		
		

		/* Lanczos
		Eigen::VectorXd v;
		Eigen::VectorXd v_1;
		Eigen::VectorXd vt;
		Eigen::VectorXd u;
		double lalpha;
		double leta;*/




		std::complex<double> beta, beta_1, gamma, gamma_1, gamma_2;
		std::complex<double> alpha, eta, eta_p1;
		std::complex<double> c, c_1, s, s_1;
		// Circular buffers
		circularbuff<std::complex<double>, 8> w;
		circularbuff<std::complex<double>, 8> wt;
		std::complex<double> rt1, r1, r2, r3;
		circularbuff<double, 8> err;
		//double merr;
		double r0norm;
		double r0norm2;
		

		circularbuff<double,8> eta_;
		circularbuff<double,8> alpha_;
		


		void init();

		//void saveVals(const char* fname, matrixcomplex &mat);
		//void saveVals(const char* fname, const Eigen::VectorXcd &vec);
		//void saveVals(const char* fname, std::complex<double> &val, bool app = true);
		//void saveVals(const char* fname, double val, bool app = true);

	public:
	
		void setrefresh() {
			this->refresh_at_next = true;			
		}
		
		int getIteration() const {
			return this->it;
		}
		
		double getResidueSquaredNorm() const {
			return this->r.norm()*this->r.norm();
		}
		
		double getErrorl2Estimate() const {
			if(it>2) {
				//double errit = this->err[it - 1];
				//double linfnorm = precond.getLINFinityNorm();
				return r0norm*r0norm*(this->err[it-1])*precond.getLINFinityNorm();
			} else return 0;
		}

		CG_SolverComplex(matrixcomplex &A, Eigen::VectorXcd &b, const SparseIncompleteLLTComplex &precond);
		CG_SolverComplex(matrixcomplex &A, Eigen::VectorXcd &b, const Eigen::VectorXcd &x0, const SparseIncompleteLLTComplex &precond);
		void do_iteration();
		
		const Eigen::VectorXcd &getX() const {
			return this->x;
		}

		Eigen::VectorXcd getY() const {
			Eigen::VectorXcd y(x);
			this->precond.halfMultInPlace(y);
			return y;
		}

		const Eigen::VectorXcd &getR() const {
			return this->r;
		}
		
		const SparseIncompleteLLTComplex &getPrecond() {
			return this->precond;
		}
		
		double getLastE2() const {
			return getResidueSquaredNorm() / r0norm2;
		}

		//matrix buildJacobiMatrx();

};


#endif /* SOLVERCOMPLEX_H_ */
