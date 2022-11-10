#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/intcoef.h"
#include "../src/solver_lb.h"
#include "../src/incomplete_qr_builder.h"
#include "../src/sparse_incomplete_qr.h"

#include "util/timestamp/timestamp.h"

struct eigen_double_qr_engine {
	typedef double scalar;
	typedef double real;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> symMatrix;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> matrix;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> real_vector;
	typedef SparseIncompleteQR<double> preconditioner;

    static SparseIncompleteQRBuilder<double> &getQRBuilder() {
        static SparseIncompleteQRBuilder<double> builder;
        return builder;
    }

	inline static void product_ii_ic_vector(vector &dest_i, vector &dest_c, const symMatrix &ii, const symMatrix &ic, const vector &b) {
		dest_i.noalias() = ii.selfadjointView<Eigen::Lower>()*b;
		dest_c.noalias() = ic*b;
	}

	inline static void transposed_product_ii_ic_vector(vector &dest, const symMatrix &ii, const symMatrix &ic, const vector &b_i, const vector &b_c) {
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
		return new preconditioner(8, 16, A_ii, A_ic, getQRBuilder());
	}
};


int main(int argc, char *argv[])
 {
     long start, stop;
     bool is2d;
     if(argc <2) {
    	 std::cerr << "need mesh filename\n";
    	 return 1;
     }
     std::cout << "Parsing mesh file..." << std::flush;
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
     observations<double> *readingsScalar;
     input->initProblem(argv[1]);
     input->buildNodeCoefficients();
     input->prepareSkeletonMatrix();
     input->createCoef2KMatrix();
     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;

     std::cout << "Preparing dummy solution..." << std::flush;
     std::vector<double> sol;
     for(int i=0;i<input->getNumCoefficients();i++) sol.push_back((mincond+maxcond)/2);
     std::cout << "Done\n" << std::flush;

     std::unique_ptr<matrix> Aii, Aic, Acc;
     std::cout << "Assembling matrices..." << std::flush;
     {
            matrix *_aii, *_aic, *_acc;
            assembleProblemMatrix_lb(sol.data(), &_aii, &_aic, &_acc, *input);
            Aii.reset(_aii);
            Aic.reset(_aic);
            Acc.reset(_acc);
     }
     std::cout << "Done\n" << std::flush;


     std::cout << "Building QR preconditioner...\n";
     std::unique_ptr<LB_Solver_A<eigen_double_qr_engine>::Preconditioner> precondqr;
     for(int i = 0; i<1000; i++) {
         precondqr.reset(LB_Solver_A<eigen_double_qr_engine>::makePreconditioner(*Aii, *Aic));
     }
     start = get_usec_timestamp();
     for(int i = 0; i<50000; i++) {
         precondqr.reset(LB_Solver_A<eigen_double_qr_engine>::makePreconditioner(*Aii, *Aic));
     }
     stop = get_usec_timestamp();
     std::cout << "QR preconditioner: "  <<  ((double)(stop - start))/50000 << std::endl;

     return 0;
 }
