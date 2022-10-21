/*
 * solver.h
 *
 *  Created on: October 20, 2022
 *      Author: thiago
 */

#ifndef SOLVER_LB_COMPLEX2_H_
#define SOLVER_LB_COMPLEX2_H_

#include <complex>
#include "solver_lb.h"
#include "incomplete_qr_complex2.h"

#define LQ_COMPLEX_NQ (8)
#define LQ_COMPLEX_NR (16)

struct eigen_complexdouble_engine {
	typedef std::complex<double> scalar;
	typedef double real;
	typedef Eigen::SparseMatrix<scalar, Eigen::ColMajor> symMatrix;
	typedef Eigen::SparseMatrix<scalar, Eigen::ColMajor> matrix;
	typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector;
	typedef SparseIncompleteQRComplex2 preconditioner;

    struct __complex_wrapper_noconj : public std::complex<double> {
        using std::complex<double>::complex;
    };

    inline static void symmetric_complex_mult_and_assign(const symMatrix &m, const vector &x, vector &res)
    {
        const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *mm  = (const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *)&m;
        const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *xx  = (const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&x;
        Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *rr  = (Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&res;
        rr->noalias() = mm->selfadjointView<Eigen::Lower>()*(*xx);
    }

    inline static void symmetric_complex_mult_and_subtract(const symMatrix &m, const vector &x, vector &res)
    {
        const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *mm  = (const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *)&m;
        const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *xx  = (const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&x;
        Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *rr  = (Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&res;
        rr->noalias() -= mm->selfadjointView<Eigen::Lower>()*(*xx);
    }


	inline static void product_ii_ic_vector(vector &dest_i, vector &dest_c, const symMatrix &ii, const matrix &ic, const vector &b) {
		//dest_i.noalias() = ii.selfadjointView<Eigen::Lower>()*b;
        symmetric_complex_mult_and_assign(ii, b, dest_i);
		dest_c.noalias() = ic*b;
	}

	inline static void transposed_product_ii_ic_vector(vector &dest, const symMatrix &ii, const matrix &ic, const vector &b_i, const vector &b_c) {
		//dest.noalias() = ii.selfadjointView<Eigen::Lower>()*b_i;
        symmetric_complex_mult_and_assign(ii, b_i, dest);
		dest.noalias() += ic.transpose()*b_c;
	}

	inline static void j_minus_a_x_phi(vector &dest_i, vector &dest_c, const matrix &ic, const symMatrix &cc, const vector &j, const vector &phi) {
		dest_i.noalias() = -ic.transpose()*phi;
		dest_c.noalias() = j;
        //dest_c.noalias() -= cc.selfadjointView<Eigen::Lower>()*phi;
        symmetric_complex_mult_and_subtract(cc, phi, dest_c);
	}

	inline static void subtract_a_x(vector &dest_i, vector &dest_c, const symMatrix &ii,  const symMatrix &ic, const vector &x) {
		//dest_i.noalias() -= ii.selfadjointView<Eigen::Lower>()*x;
        symmetric_complex_mult_and_subtract(ii, x, dest_i);
		dest_c.noalias() -= ic*x;
	}

	inline static preconditioner *make_new_preconditioner(const symMatrix &A_ii, const matrix &A_ic) {
		return new preconditioner(LQ_COMPLEX_NQ, LQ_COMPLEX_NR, A_ii, A_ic);
	}
};

typedef LB_Solver_A<eigen_complexdouble_engine> LB_Solver_Complex2;

#endif  // SOLVER_LB_COMPLEX2_H_
