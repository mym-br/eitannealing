#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/intcoef.h"
#include "../src/solver_lb.h"
#include "../src/incomplete_qr_builder.h"

#include <sys/time.h>
#include <sys/resource.h>

unsigned long get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec*1e6 + t.tv_usec;
}

class SparseIncompleteQR
{
  protected:
    vectorx idiagonal;
    matrix rmatrix;

    struct MatricesStorageAdaptor {
        const matrix &ii;
        const matrix &ic;
        unsigned long square_size;
        // FIXME: This map can be reused!
        std::vector<std::vector<std::pair<unsigned long, double > > > upperMap;
        MatricesStorageAdaptor(const matrix &Aii_low, const matrix &Aic):
             ii(Aii_low), ic(Aic), square_size(Aii_low.cols()), upperMap(Aii_low.rows())  {}
        void iterateOverColumn(unsigned long j, std::function<void(unsigned long, double)> &&f) {
            // FIXME: This assumes that the columns will be iterated in order, so the current column has already
            //  an appropriate upper map
            // First iterate over upper elements
            for(auto [j, x] : upperMap[j]) f(j, x);
            // Now, lower elements and fill map for next columns
            for(typename matrix::InnerIterator it(ii, j); it; ++it) {
                if(it.index()>j) {
                  upperMap[it.index()].push_back(std::make_pair(j, it.value()));
                }
                f(it.index(), it.value());
            }
            // Finally, ic elements
            for(typename matrix::InnerIterator it(ic, j); it; ++it) {
                f(it.index() + square_size, it.value());
            }
        }
        unsigned long rows() const { return square_size+ic.rows(); }
        unsigned long cols() const { return square_size; }
    };

  public:

    SparseIncompleteQR(unsigned long nr, unsigned long nq, const matrix &Aii_low, const matrix &Aic):
        idiagonal(Aii_low.rows()), rmatrix(Aii_low.rows(), Aii_low.rows()) {

        SparseIncompleteQRBuilder<double> builder;

        this->rmatrix.reserve(Aii_low.rows()*nr);

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, double x) {
            this->idiagonal(j) = x;
         },
         [this](unsigned long i, unsigned long j, double x) {
             this->rmatrix.insert(i,j) = x;
         });
        rmatrix.makeCompressed();
        idiagonal = idiagonal.cwiseInverse();
        for(int i = 0; i<rmatrix.outerSize(); i++)
            rmatrix.col(i) *= idiagonal(i);
    }

    void solveInPlace(vectorx &b) const {
      rmatrix.triangularView<Eigen::UnitUpper>().solveInPlace(b);
      b = b.cwiseProduct(idiagonal);
    }

     // conjugated transpose
    void solveInPlaceT(vectorx &b) const {
        b = b.cwiseProduct(idiagonal);
        rmatrix.triangularView<Eigen::UnitUpper>().transpose().solveInPlace(b);
    }
};

struct eigen_double_qr_engine {
	typedef double scalar;
	typedef double real;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> symMatrix;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> matrix;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector;
	typedef SparseIncompleteQR preconditioner;

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
		return new preconditioner(16, 8, A_ii, A_ic);
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
     std::unique_ptr<LB_Solver_EG_Estimate> solver_e;
     solver_e.reset(new LB_Solver_EG_Estimate(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, 600, (float)0.00001
     ));
     std::cout << "Done\n" << std::flush;
     
     double a = solver_e->getLeastEvEst();
     std::cout << "Solved performed " << solver_e->getIteration() << " iterations\n";
     std::cout << "Least eigenvalue estimation: " << a << "\n";
     
     std::cout << "Preparing solver..." << std::flush;
     std::unique_ptr<LB_Solver> solver;
     solver.reset(new LB_Solver(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precond, a
     ));
     std::cout << "Done\n" << std::flush;
     
     solver->do_iteration();
     for(unsigned i = 0; i<100; i++) {
         std::cout << solver->getIteration() << ": (min): " << solver->getMinErrorl2Estimate() << " (max): " << solver->getMaxErrorl2Estimate() << std::endl;
         solver->do_iteration();
     }
     
     std::unique_ptr<LB_Solver_A<eigen_double_qr_engine>::Preconditioner> precondqr;
     std::cout << "Preparing preconditioner (QR)..." << std::flush;
     precondqr.reset(LB_Solver_A<eigen_double_qr_engine>::makePreconditioner(*Aii, *Aic));
     std::cout << "Done\n" << std::flush;
     
     std::cout << "Preparing solver (QR)..." << std::flush;
     std::unique_ptr<LB_Solver_A<eigen_double_qr_engine> > solverqr;
     solverqr.reset(new LB_Solver_A<eigen_double_qr_engine>(
        Aii, Aic, Acc, Eigen::VectorXd(readingsScalar->getCurrents()[0].tail(input->getGenericElectrodesCount())), Eigen::VectorXd(readingsScalar->getTensions()[0]), *precondqr, a
     ));
     std::cout << "Done\n" << std::flush;
     
     solverqr->do_iteration();
     for(unsigned i = 0; i<100; i++) {
         std::cout << solverqr->getIteration() << ": (min): " << solverqr->getMinErrorl2Estimate() << " (max): " << solverqr->getMaxErrorl2Estimate() << std::endl;
         solverqr->do_iteration();
     }
    
     return 0;
 }
