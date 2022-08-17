#include <iostream>
#include "../src/incomplete_qr_builder.h"
#include "../src/util/fill_with_smallest.hpp"
#include "../incomplete_qr_complex2.h"
#include "../incomplete_qr_complex.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include <complex>


std::complex<double> simpleKh3[] = {
  {6.0, 0.0}, {-1.0, 0.0}, {-2.0, 0.0}, {0.0, 0.0},
  {-1.0, 0.0}, {6.0, 0.0}, {0.0, 0.0}, {-2.0, 0.0},
  {-2.0, 0.0}, {0.0, 0.0}, {6.0, 0.0}, {-1.0, 0.0},
  {0.0, 0.0}, {-2.0, 0.0}, {-1.0, 0.0}, {6.0, 1.0},
  {-3.0, 0.0}, {-3.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
  {0.0, 0.0}, {0.0, 0.0}, {-3.0, 0.0}, {-3.0, -1.0},
};


std::complex<double> simpleKh2[] = {
  {-3.0, -1.0}, {-3.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
  {0.0, 0.0}, {0.0, 0.0}, {-3.0, 0.0}, {-3.0, 0.0},
  {6.0, 1.0}, {-1.0, 0.0}, {-2.0, 0.0}, {0.0, 0.0},
  {-1.0, 0.0}, {6.0, 0.0}, {0.0, 0.0}, {-2.0, 0.0},
  {-2.0, 0.0}, {0.0, 0.0}, {6.0, 0.0}, {-1.0, 0.0},
  {0.0, 0.0}, {-2.0, 0.0}, {-1.0, 0.0}, {6.0, 0.0}
};

std::complex<double> testR2[] = {
  {7.211102550927978, 0.0}, {-0.41602514716892186,-0.27735009811261467},{-3.328201177351375,0.2773500981126146},{0.5547001962252291, 0.0},
  {0.0, 0.0}, {7.053367989832942, 0.0}, {0.3817052642352578, 0.14722917334788516}, {-3.369912189962704, -0.021811729384871872},
  {0.0, 0.0}, {0.0, 0.0}, {6.169058052718128, 0.0}, {-0.23215226567515956,-0.027068766200180533},
  {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {6.160058221029499, 0.0}
};

double simpleKh1[] = {
  -3.0, -3.0, 0.0, 0.0,
  0.0, 0.0, -3.0, -3.0,
  6.0, -1.0, -2.0, 0.0,
  -1.0, 6.0, 0.0, -2.0,
  -2.0, 0.0, 6.0, -1.0,
  0.0, -2.0, -1.0, 6.0
};

double testR1[] = {
  7.0, 0.42857142857142855,-3.4285714285714284,0.0,
  0.0,6.923271995742359,0.2888807490489974,-3.466568988587969,
  0.,0.,6.002671099874577,0.6629701628693009,
  0.,0.,0.,5.875706681229786
};


int test1(void) {

  Eigen::SparseMatrix<double, Eigen::ColMajor> a =
    Eigen::Map<Eigen::MatrixXd>(simpleKh1,4,6).sparseView(0.0, 0.00001).transpose();
  SparseIncompleteQRBuilder<double> builder;
  Eigen::SparseMatrix<double, Eigen::ColMajor> RMatrix(4, 4);
  RMatrix.reserve(make_sizestype_adaptor([](unsigned long i){return (i+1)>4?4:(i+1);}));

  builder.buildRMatrix(a, 3, 3, [&RMatrix](unsigned long j, double x) {
    RMatrix.insert(j,j) = 1/x;
  }, [&RMatrix](unsigned long i, unsigned long j, double x) {
    RMatrix.insert(i,j) = x;
  });

  RMatrix.makeCompressed();
  if(
    (RMatrix - Eigen::SparseMatrix<double, Eigen::ColMajor>(Eigen::Map<Eigen::MatrixXd>(testR1,4,4).sparseView(0.0, 0.001).transpose())).squaredNorm()
      > 1e-7)
      return 1;
  return 0;
}

int test2(void) {

  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> a =
    Eigen::Map<Eigen::MatrixXcd>(simpleKh2,4,6).sparseView(0.0, 0.00001).transpose();

  SparseIncompleteQRBuilder<std::complex<double> > builder;

  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> RMatrix(4, 4);
  RMatrix.reserve(make_sizestype_adaptor([](unsigned long i){return (i+1)>4?4:(i+1);}));

  builder.buildRMatrix(a, 4, 5, [&RMatrix](unsigned long j, double x) {
    RMatrix.insert(j,j) = 1/x;
  }, [&RMatrix](unsigned long i, unsigned long j, std::complex<double> x) {
    RMatrix.insert(i,j) = x;
  });

  RMatrix.makeCompressed();

  if(
    (RMatrix - Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor>(Eigen::Map<Eigen::MatrixXcd>(testR2,4,4).sparseView(0.0, 0.001).transpose())).squaredNorm()
      > 1e-7)
      return 1;
  return 0;
}


int test3(void) {
    matrixcomplex Kh(Eigen::Map<Eigen::MatrixXcd>(simpleKh3,4,6).sparseView(0.0, 0.00001).transpose());
    matrixcomplex Kii(Kh.block(0,0,4,4).triangularView<Eigen::Lower>());
    matrixcomplex Kic(Kh.block(4,0,2,4));
    matrix Kii_R(Kii.real()), Kii_I(Kii.imag()), Kic_R(Kic.real()), Kic_I(Kic.imag());

    std::cout << "Khat:" << Kh << std::endl;

    SparseIncompleteQRComplex2 precond2(4, 4, Kii, Kic);
    vectorxcomplex b(vectorxcomplex::Ones(4));
    precond2.solveInPlace(b);
    std::cout << "QRComplex2:\n" << b.real() << std::endl << b.imag() << std::endl;

    vectorx b_R(vectorx::Ones(4)), b_I(vectorx::Zero(4));
    SparseIncompleteQRComplex precond(4, 4, Kii_R, Kii_I, Kic_R, Kic_I);
    precond.solveInPlaceC(b_R, b_I);
    std::cout << "QRComplex:\n" << b_R << std::endl << b_I << std::endl;

    SparseIncompleteQRBuilder<std::complex<double> > builder;

    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> a(Eigen::Map<Eigen::MatrixXcd>(simpleKh3,4,6).sparseView(0.0, 0.00001).transpose());

    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> RMatrix(4, 4);
    RMatrix.reserve(make_sizestype_adaptor([](unsigned long i){return (i+1)>4?4:(i+1);}));
    builder.buildRMatrix(a, 4, 4, [&RMatrix](unsigned long j, double x) {
      RMatrix.insert(j,j) = 1/x;
    }, [&RMatrix](unsigned long i, unsigned long j, std::complex<double> x) {
      RMatrix.insert(i,j) = x;
    });
    RMatrix.makeCompressed();

    b = vectorx::Ones(4);

    RMatrix.triangularView<Eigen::Upper>().solveInPlace(b);

    std::cout << "Builder:\n" << b_R << std::endl << b_I << std::endl;

    b = vectorx::Ones(4);
    precond2.solveInPlaceT(b);
    std::cout << "QRComplex2 Transposed:\n" << b.real() << std::endl << b.imag() << std::endl;

    b_R = vectorx::Ones(4);
    b_I = vectorx::Zero(4);
    precond.solveInPlaceCT(b_R, b_I);

    std::cout << "QRComplex Transposed:\n" << b_R << std::endl << b_I << std::endl;

    b = vectorx::Ones(4);
    RMatrix.triangularView<Eigen::Upper>().transpose().solveInPlace(b);

    std::cout << "Builder Transposed:\n" << b_R << std::endl << b_I << std::endl;



    return 0;
}

int main(int argc, char **argv) {
  int r;
  std::cout << "Test 1... ";
  r = test1();
  if(r) {
    std::cout << "Failed!\n";
    return r;
  }
  std::cout << "Succeeded.\nTest 2... ";
  r = test2();
  if(r) {
    std::cout << "Failed!\n";
    return r;
  }
  std::cout << "Succeeded.\nTest 3... ";
  r = test3();
  if(r) {
    std::cout << "Failed!\n";
    return r;
  }
  std::cout << "Succeeded.\n";
  return 0;
}
