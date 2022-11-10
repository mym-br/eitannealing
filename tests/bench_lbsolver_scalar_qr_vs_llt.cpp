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
		return new preconditioner(8, 16, A_ii, A_ic, getQRBuilder());
	}
};


int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <4) {
    	 std::cerr << "need mesh filename, currents and tensions\n";
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

     std::cout << "Parsing observations..." << std::flush;
     readingsScalar = new observations<double>;
     readingsScalar->initObs((const char **)&argv[2], argv[3], input->getNodesCount(), input->getGenericElectrodesCount());
     std::cout << "Done\n" << std::flush;

     std::cout << "Preparing dummy solution..." << std::flush;
     double *sol = new double[input->getNumCoefficients()];
     for(int i=0;i<input->getNumCoefficients();i++) sol[i]=(mincond+maxcond)/2;
     std::cout << "Done\n" << std::flush;

     matrix *Aii, *Aic, *Acc;
     std::cout << "Assembling matrices..." << std::flush;
     assembleProblemMatrix_lb(sol, &Aii, &Aic, &Acc, *input);
     std::cout << "Done\n" << std::flush;

     std::unique_ptr<LB_Solver::Preconditioner> precond;
     std::cout << "Preparing preconditioner..." << std::flush;
     precond.reset(LB_Solver::makePreconditioner(*Aii, *Aic));
     std::cout << "Done\n" << std::flush;

     std::cout << "Preparing solver (eigenvalue estimator)..." << std::flush;
     double a;
     LB_Solver::vector eigenvector;
     std::unique_ptr<LB_Solver> solver_e = std::make_unique<LB_Solver>(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, 800, (float)0.00000001,
        a, eigenvector
     );
     std::cout << "Done\n" << std::flush;

     std::cout << "Solved performed " << solver_e->getIteration() << " iterations\n";
     std::cout << "Least eigenvalue estimation: " << a << "\n";

     long start, stop;
     std::unique_ptr<LB_Solver_A<eigen_double_qr_engine>::Preconditioner> precondqr;

     std::cout << "Preparing solver..." << std::flush;
     std::unique_ptr<LB_Solver> solver = std::make_unique<LB_Solver>(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, a
     );
     std::cout << "Done\n" << std::flush;

     solver->do_iteration();
     for(unsigned i = 0; i<300; i++) {
         std::cout << solver->getIteration() << ": (min): " << solver->getMinErrorl2Estimate() << " (max): " << solver->getMaxErrorl2Estimate() << std::endl;
         solver->do_iteration();
     }

     std::cout << "Preparing preconditioner (QR)..." << std::flush;
     precondqr.reset(LB_Solver_A<eigen_double_qr_engine>::makePreconditioner(*Aii, *Aic));
     std::cout << "Done\n" << std::flush;

     std::cout << "Preparing solver (QR)..." << std::flush;
     std::unique_ptr<LB_Solver_A<eigen_double_qr_engine> > solverqr = std::make_unique<LB_Solver_A<eigen_double_qr_engine> >(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precondqr, a
     );
     std::cout << "Done\n" << std::flush;

     solverqr->do_iteration();
     for(unsigned i = 0; i<300; i++) {
         std::cout << solverqr->getIteration() << ": (min): " << solverqr->getMinErrorl2Estimate() << " (max): " << solverqr->getMaxErrorl2Estimate() << std::endl;
         solverqr->do_iteration();
     }

     std::cout << "Benchmarking iteration costs..." << std::flush;
     std::cout << "LLT preconditioner...\n";
     for(int i = 0; i<100; i++) {
         solver->do_iteration();
     }
     start = get_usec_timestamp();
     for(int i = 0; i<100000; i++) {
         solver->do_iteration();
     }
     stop = get_usec_timestamp();
     std::cout << "LLT preconditioner: "  <<  ((double)(stop - start))/100000 << std::endl;

     std::cout << "QR preconditioner...\n";
     for(int i = 0; i<100; i++) {
         solverqr->do_iteration();
     }
     start = get_usec_timestamp();
     for(int i = 0; i<100000; i++) {
         solverqr->do_iteration();
     }
     stop = get_usec_timestamp();
     std::cout << "QR preconditioner: "  <<  ((double)(stop - start))/100000 << std::endl;

     std::cout << "Benchmarking construction/destruction costs..." << std::flush;
     std::cout << "LLT preconditioner...\n";
     for(int i = 0; i<20; i++) {
         precond.reset(LB_Solver::makePreconditioner(*Aii, *Aic));
     }
     start = get_usec_timestamp();
     for(int i = 0; i<4000; i++) {
         precond.reset(LB_Solver::makePreconditioner(*Aii, *Aic));
     }
     stop = get_usec_timestamp();
     std::cout << "LLT preconditioner: "  <<  ((double)(stop - start))/4000 << std::endl;

     std::cout << "QR preconditioner...\n";
     for(int i = 0; i<20; i++) {
         precondqr.reset(LB_Solver_A<eigen_double_qr_engine>::makePreconditioner(*Aii, *Aic));
     }
     start = get_usec_timestamp();
     for(int i = 0; i<5000; i++) {
         precondqr.reset(LB_Solver_A<eigen_double_qr_engine>::makePreconditioner(*Aii, *Aic));
     }
     stop = get_usec_timestamp();
     std::cout << "QR preconditioner: "  <<  ((double)(stop - start))/5000 << std::endl;

     return 0;
 }
