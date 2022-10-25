#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"

#include "util/timestamp/timestamp.h"
#include "solver_lb.h"
#include <complex>

struct __complex_wrapper_alwaysconj;

struct __complex_wrapper_noconj : public std::complex<double> {
    using std::complex<double>::complex;

    inline __complex_wrapper_noconj operator *(__complex_wrapper_alwaysconj other) const;
};

struct __complex_wrapper_alwaysconj : public std::complex<double> {
    using std::complex<double>::complex;

    __complex_wrapper_noconj operator *(__complex_wrapper_noconj other) const {
        return __complex_wrapper_noconj(std::conj(*this)*other);
    }
};

inline __complex_wrapper_noconj __complex_wrapper_noconj::operator *(__complex_wrapper_alwaysconj other) const {
    return other*(*this);
}


template<>
struct Eigen::ScalarBinaryOpTraits<__complex_wrapper_alwaysconj, __complex_wrapper_noconj, Eigen::internal::scalar_product_op<__complex_wrapper_alwaysconj, __complex_wrapper_noconj> > {
    typedef __complex_wrapper_noconj ReturnType;
};

inline void symmetric_complex_mult_and_assign(const matrixcomplex &m, const vectorxcomplex &x, vectorxcomplex &res)
{
    const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *mm  = (const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *)&m;
    const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *xx  = (const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&x;
    Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *rr  = (Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&res;
    rr->noalias() = mm->selfadjointView<Eigen::Lower>()*(*xx);
}

inline void symmetric_conjtranspose_complex_mult_and_assign(const matrixcomplex &m, const vectorxcomplex &x, vectorxcomplex &res)
{
    const Eigen::SparseMatrix<__complex_wrapper_alwaysconj, Eigen::ColMajor> *mm  = (const Eigen::SparseMatrix<__complex_wrapper_alwaysconj, Eigen::ColMajor> *)&m;
    const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *xx  = (const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&x;
    Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *rr  = (Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&res;
    rr->noalias() = mm->selfadjointView<Eigen::Lower>()*(*xx);
}


int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <2) {
    	 std::cerr << "need mesh filename\n";
    	 return 1;
     }
     std::cout << "Parsing mesh file..." << std::flush;
     std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
     input->initProblem(argv[1]);
     input->buildNodeCoefficients();
     std::cout << "Done\n" << std::flush;
     std::cout << "Solution with " << input->getNumCoefficients() << " coefficients\n" << std::flush;
     std::cout << "Preparing dummy solution..." << std::flush;
     double *sol_R = new double[input->getNumCoefficients()];
     for(int i=0;i<input->getNumCoefficients();i++) sol_R[i]=1.0;


     double *sol_I = new double[input->getNumCoefficients()];
     for(int i=0;i<32;i++) sol_I[i]=0.0;
     for(int i=32;i<input->getNumCoefficients();i++) sol_I[i]=1.0;
     std::cout << "Done\n" << std::flush;
     matrix *Aii_R, *Aii_I, *Aic_R, *Aic_I, *Acc_R, *Acc_I;
     std::cout << "Assembling matrices..." << std::flush;
     assembleProblemMatrix_lb(sol_R, &Aii_R, &Aic_R, &Acc_R, *input);
     assembleProblemMatrix_lb(sol_I, &Aii_I, &Aic_I, &Acc_I, *input);
     std::cout << "Done\n" << std::flush;
     
     matrixcomplex Aii = (Eigen::MatrixXcd(*Aii_R) + Complex(0,1)*Eigen::MatrixXcd(*Aii_I)).sparseView();
     matrixcomplex Aic = (Eigen::MatrixXcd(*Aic_R) + Complex(0,1)*Eigen::MatrixXcd(*Aic_I)).sparseView();
     
     

     matrix ii(Aii_R->selfadjointView<Eigen::Lower>());
     matrix c(Aii_R->rows()+Aic_R->rows(), Aii_R->cols());
     
     c.reserve(ii.nonZeros() + Aic_R->nonZeros());
     for(int i = 0; i < Aii_R->cols(); i++) {
        c.startVec(i);
        for(Eigen::SparseMatrix<double>::InnerIterator iti(ii, i); iti; ++iti) {
            c.insertBack(iti.row(), i) = iti.value();
        }
        for(Eigen::SparseMatrix<double>::InnerIterator itc(*Aic_R, i); itc; ++itc) {
            c.insertBack(ii.rows()+itc.row(), i) = itc.value();
        }     
     }
     c.finalize();
     
     Eigen::VectorXcd  x(Eigen::VectorXd::Random(c.cols())*Complex(1,0)+Eigen::VectorXd::Random(c.cols())*Complex(0,1)), y, y2, x2(Eigen::VectorXcd::Random(c.rows()));
     
     symmetric_complex_mult_and_assign(Aii, x, y);
     
     std::cout << "Residual (real): " << (y.real() - (Aii_R->selfadjointView<Eigen::Lower>()*x.real() - Aii_I->selfadjointView<Eigen::Lower>()*x.imag())).squaredNorm() << std::endl;
     std::cout << "Residual (imag): " << (y.imag() - (Aii_R->selfadjointView<Eigen::Lower>()*x.imag() + Aii_I->selfadjointView<Eigen::Lower>()*x.real())).squaredNorm() << std::endl;
     
     Eigen::VectorXcd  y_r, y_i;
     long start, stop;
     
     
     
     for(int i = 0; i<100; i++)
        symmetric_complex_mult_and_assign(Aii, x, y);
     start = get_usec_timestamp();
     for(int i = 0; i<40000; i++)
            symmetric_complex_mult_and_assign(Aii, x, y);
     stop = get_usec_timestamp();
     std::cout << "Eigen complex symmetric: "  <<  ((double)(stop - start))/40000 << std::endl;
     
     
     for(int i = 0; i<100; i++) {
        y_r.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x.real() - Aii_I->selfadjointView<Eigen::Lower>()*x.imag();
        y_i.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x.imag() + Aii_I->selfadjointView<Eigen::Lower>()*x.real();
     }
     start = get_usec_timestamp();
     for(int i = 0; i<4000; i++)  {
        y_r.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x.real() - Aii_I->selfadjointView<Eigen::Lower>()*x.imag();
        y_i.noalias() = Aii_R->selfadjointView<Eigen::Lower>()*x.imag() + Aii_I->selfadjointView<Eigen::Lower>()*x.real();
     }
     stop = get_usec_timestamp();
     std::cout << "Split components symmetric: "  <<  ((double)(stop - start))/4000 << std::endl;

     std::cout << "Conjugated product:\n";

     symmetric_conjtranspose_complex_mult_and_assign(Aii, x, y);

     std::cout << "Residual (real): " << (y.real() - (Aii_R->selfadjointView<Eigen::Lower>()*x.real() + Aii_I->selfadjointView<Eigen::Lower>()*x.imag())).squaredNorm() << std::endl;
     std::cout << "Residual (imag): " << (y.imag() - (Aii_R->selfadjointView<Eigen::Lower>()*x.imag() - Aii_I->selfadjointView<Eigen::Lower>()*x.real())).squaredNorm() << std::endl;

     for(int i = 0; i<100; i++)
        symmetric_conjtranspose_complex_mult_and_assign(Aii, x, y);
     start = get_usec_timestamp();
     for(int i = 0; i<40000; i++)
            symmetric_conjtranspose_complex_mult_and_assign(Aii, x, y);
     stop = get_usec_timestamp();
     std::cout << "Eigen complex conjugated symmetric: "  <<  ((double)(stop - start))/40000 << std::endl;

     
     
     
    
     
     
     return 0;
 }
