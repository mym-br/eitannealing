#include <iostream>
#include "../src/incomplete_qr_builder.h"
#include "../src/incomplete_qr_complex.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include <complex>


int test1(void) {
  std::vector<std::pair< std::pair<int, int>, std::complex<double> > > AiiV = {
    {{0,0}, {6.0, 1.0}},
    {{1,0}, {-1.0, 0.0}},
    {{2,0}, {-2.0, 0.0}},
    {{1,1}, {6.0, 0.0}},
    {{3,1}, {-2.0, 0.0}},
    {{2,2}, {6.0, 0.0}},
    {{3,2}, {-1.0, 0.0}},
    {{3,3}, {6.0, 0.0}}
  };

  std::vector<std::pair< std::pair<int, int>, std::complex<double> > > AicV = {
    {{0,0}, {-3.0, -1.0}},
    {{0,1}, {-3.0, 0.0}},
    {{1,2}, {-3.0, 0.0}},
    {{1,3}, {-3.0, 0.0}}
  };

  Eigen::SparseMatrix<double, Eigen::ColMajor> AiiR(4, 4), AiiI(4, 4), AicR(2, 4), AicI(2, 4);

  for(auto [ij, x] : AiiV) {
    AiiR.insert(ij.first, ij.second) = x.real();
    AiiI.insert(ij.first, ij.second) = x.imag();
  }

  for(auto [ij, x] : AicV) {
    AicR.insert(ij.first, ij.second) = x.real();
    AicI.insert(ij.first, ij.second) = x.imag();
  }

  SparseIncompleteQRComplex<double> qr(4, 5, AiiR, AiiI, AicR, AicI);

  Eigen::VectorXd x_R(4), x_I(4);
  Eigen::VectorXd xi_R(4), xi_I(4);

  x_R << 1, 1, 1, 1;
  x_I << 1, 1, 1, 1;

  qr.solveInPlaceC(x_R, x_I);

  xi_R << 0.21432671692081437, 0.2132957768126756, 0.16749598599918508, 0.16233612802329572;
  xi_I << 0.21786640837878715, 0.2072004687291166, 0.1689205920436586, 0.16233612802329572;

  if(((xi_R-x_R).squaredNorm()+(xi_I-x_I).squaredNorm())>1e-7)
    return 1;

  return 0;
}

int main(int argc, char **argv) {
  int r;
  r = test1();
  if(r) return r;
  return 0;
}
