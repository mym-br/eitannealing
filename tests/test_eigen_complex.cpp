#include <memory>
#include <ctime>
#include <iostream>
#include "../src/problem.h"
#include "../src/nodecoefficients.h"
#include "../src/solver_lb_complex.h"
#include "../src/solver_lb.h"
#include "../src/observations_complex.h"
#include "../src/intcoef.h"

#include <sys/time.h>
#include <sys/resource.h>
#include "solver_lb.h"
#include <complex>

unsigned long get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec*1e6 + t.tv_usec;
}



struct __complex_wrapper_noconj : public std::complex<double> {
    using std::complex<double>::complex;
};
inline const __complex_wrapper_noconj & conj(const __complex_wrapper_noconj &x) { return x; }

inline void symmetric_complex_mult_and_assign(const matrixcomplex &m, const vectorxcomplex &x, vectorxcomplex &res)
{
    const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *mm  = (const Eigen::SparseMatrix<__complex_wrapper_noconj, Eigen::ColMajor> *)&m;
    const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *xx  = (const Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&x;
    Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *rr  = (Eigen::Matrix<__complex_wrapper_noconj, Eigen::Dynamic, 1> *)&res;
    rr->noalias() = (*mm)*(*xx);
}


int main(int argc, char *argv[])
 {
     bool is2d;
     if(argc <2) {
    	 std::cerr << "need mesh filename\n";
    	 return 1;
     }
     //Eigen::setNbThreads(16);
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
     
     Eigen::VectorXd  x(Eigen::VectorXd::Random(c.cols())), y, y2, x2(Eigen::VectorXd::Random(c.rows()));
    
     return 0;
 }
