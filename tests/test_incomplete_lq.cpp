#include <iostream>
#include "../src/incomplete_qr_builder.h"
#include "../src/util/fill_with_smallest.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

int main(int argc, char **argv) {
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
    Eigen::Map<Eigen::MatrixXd>(simpleKh,4,6).sparseView(0.0, 0.001).transpose();
  SparseIncompleteQRBuilder<double> builder;

  Eigen::SparseMatrix<double, Eigen::ColMajor> R = builder.buildRMatrix(a, 3, 3);
  if(
    (R - Eigen::SparseMatrix<double, Eigen::ColMajor>(Eigen::Map<Eigen::MatrixXd>(testR,4,4).sparseView(0.0, 0.001).transpose())).squaredNorm()
      < 1e-7)
      return 0;
  else return 1;
}
