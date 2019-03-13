#include <iostream>
#include "../src/incomplete_qr_builder.h"
#include "../src/util/fill_with_smallest.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include <complex>

int test1(void) {
  double simpleKh[] = {
    -3.0, -3.0, 0.0, 0.0,
    0.0, 0.0, -3.0, -3.0,
    6.0, -1.0, -2.0, 0.0,
    -1.0, 6.0, 0.0, -2.0,
    -2.0, 0.0, 6.0, -1.0,
    0.0, -2.0, -1.0, 6.0
  };

  double testR[] = {
    7.0, 0.42857142857142855,-3.4285714285714284,0.0,
    0.0,6.923271995742359,0.2888807490489974,-3.466568988587969,
    0.,0.,6.002671099874577,0.6629701628693009,
    0.,0.,0.,5.875706681229786
  };

  Eigen::SparseMatrix<double, Eigen::ColMajor> a =
    Eigen::Map<Eigen::MatrixXd>(simpleKh,4,6).sparseView(0.0, 0.00001).transpose();
  SparseIncompleteQRBuilder<double> builder;
  Eigen::SparseMatrix<double, Eigen::ColMajor> RMatrix(4, 4);
  RMatrix.reserve(make_sizestype_adaptor([](unsigned long i){return (i+1)>4?4:(i+1);}));

  builder.buildRMatrix(a, 3, 3, [&RMatrix](unsigned long j, double x) {
    RMatrix.insert(j,j) = x;
  }, [&RMatrix](unsigned long i, unsigned long j, double x) {
    RMatrix.insert(i,j) = x;
  });

  RMatrix.makeCompressed();
  if(
    (RMatrix - Eigen::SparseMatrix<double, Eigen::ColMajor>(Eigen::Map<Eigen::MatrixXd>(testR,4,4).sparseView(0.0, 0.001).transpose())).squaredNorm()
      > 1e-7)
      return 1;
  return 0;
}

int test2(void) {
  std::complex<double> simpleKh[] = {
    {-3.0, -1.0}, {-3.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
    {0.0, 0.0}, {0.0, 0.0}, {-3.0, 0.0}, {-3.0, 0.0},
    {6.0, 1.0}, {-1.0, 0.0}, {-2.0, 0.0}, {0.0, 0.0},
    {-1.0, 0.0}, {6.0, 0.0}, {0.0, 0.0}, {-2.0, 0.0},
    {-2.0, 0.0}, {0.0, 0.0}, {6.0, 0.0}, {-1.0, 0.0},
    {0.0, 0.0}, {-2.0, 0.0}, {-1.0, 0.0}, {6.0, 0.0}
  };

  std::complex<double> testR[] = {
    {7.211102550927978, 0.0}, {-0.41602514716892186,-0.27735009811261467},{-3.328201177351375,0.2773500981126146},{0.5547001962252291, 0.0},
    {0.0, 0.0}, {7.053367989832942, 0.0}, {0.3817052642352578, 0.14722917334788516}, {-3.369912189962704, -0.021811729384871872},
    {0.0, 0.0}, {0.0, 0.0}, {6.169058052718128, 0.0}, {-0.23215226567515956,-0.027068766200180533},
    {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {6.160058221029499, 0.0}
  };
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> a =
    Eigen::Map<Eigen::MatrixXcd>(simpleKh,4,6).sparseView(0.0, 0.00001).transpose();

  SparseIncompleteQRBuilder<std::complex<double> > builder;

  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> RMatrix(4, 4);
  RMatrix.reserve(make_sizestype_adaptor([](unsigned long i){return (i+1)>4?4:(i+1);}));

  builder.buildRMatrix(a, 4, 5, [&RMatrix](unsigned long j, double x) {
    RMatrix.insert(j,j) = x;
  }, [&RMatrix](unsigned long i, unsigned long j, std::complex<double> x) {
    RMatrix.insert(i,j) = x;
  });

  RMatrix.makeCompressed();

  if(
    (RMatrix - Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor>(Eigen::Map<Eigen::MatrixXcd>(testR,4,4).sparseView(0.0, 0.001).transpose())).squaredNorm()
      > 1e-7)
      return 1;
  return 0;
}

int main(int argc, char **argv) {
  int r;
  r = test1();
  if(r) return r;
  r = test2();
  if(r) return r;
  return 0;
}
