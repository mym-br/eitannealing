/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_QR_COMPLEX_HPP_
#define INCOMPLETE_QR_COMPLEX_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <unordered_map>

#include "basematrix.h"
#include "incomplete_qr_builder.h"

template<class scalar> class SparseIncompleteQRComplex
{
  protected:
    std::vector<double> idiagonal;
    std::vector<std::unordered_map<unsigned long, std::complex<double> > > rows;
    std::vector<std::unordered_map<unsigned long, std::complex<double> > > cols;

    struct MatricesStorageAdaptor {
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &iiR;
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &iiI;
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &icR;
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &icI;
        // FIXME: This map can be reused!
        std::vector<std::vector<std::pair<unsigned long, std::complex<double> > > > upperMap;

        MatricesStorageAdaptor(
          const matrix &iiUpper_R, const matrix &iiUpper_I,
          const matrix &Aic_R, const matrix &Aic_I):
          iiR(iiUpper_R), iiI(iiUpper_I), icR(Aic_R), icI(Aic_I), upperMap(iiUpper_R.rows())
           {}
        void iterateOverColumn(unsigned long j, std::function<void(unsigned long, std::complex<double>)> &&f) {
            // FIXME: This assumes that the columns will be iterated in order, so the current column has already
            //  an appropriate upper map
            // First iterate over upper elements
            for(auto [j, x] : upperMap[j]) f(j, x);
            // Now, lower elements and fill map for next columns
            for(typename Eigen::SparseMatrix<scalar>::InnerIterator itR(iiR, j), itI(iiI, j); itR; ++itR, ++itI) {
                std::complex<double> val(itR.value(), itI.value());
                if(itR.index()>j) {
                  upperMap[itR.index()].push_back(std::pair<unsigned long, std::complex<double> >(j, val));
                }
                f(itR.index(), val);
            }
            // Finally, ic elements
            for(typename Eigen::SparseMatrix<scalar>::InnerIterator itR(icR, j), itI(icI, j); itR; ++itR, ++itI) {
                std::complex<double> val(itR.value(), itI.value());
                f(itR.index(), val);
            }
        }
        unsigned long rows() const { return iiR.rows(); }
        unsigned long cols() const { return iiR.cols(); }
    };

  public:

    SparseIncompleteQRComplex(
      const matrix &iiUpper_R, const matrix &iiUpper_I,
      const matrix &Aic_R, const matrix &Aic_I
    ): idiagonal(iiUpper_R.rows()), rows(iiUpper_R.rows()), cols(iiUpper_R.rows()) {

      SparseIncompleteQRBuilder<std::complex<double> > builder;

      builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(iiUpper_R, iiUpper_I, Aic_R, Aic_I), 8, 16, [this](unsigned long j, double x) {
        this->idiagonal[j] = 1/x;
      }, [this](unsigned long i, unsigned long j, std::complex<double> x) {
        this->rows[i][j] = x;
        this->cols[j][i] = x;
      });
    }

    void solveInPlaceC(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const {
      unsigned long i = idiagonal.size();
      for(; i>=0; i--) {
        for(auto [j, x] : rows[i]) {
          bR[i] -= bR[j]*x.real() - bI[j]*x.imag();
          bI[i] -= bR[j]*x.imag() + bI[j]*x.real();
        }
        bR[i] *= idiagonal[i];
        bI[i] *= idiagonal[i];
      }
    };

     // conjugated transpose
    void solveInPlaceCT(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const {
      unsigned long n = idiagonal.size();
      for(unsigned int i = 0; i<n; i++) {
        for(auto [j, x] : cols[i]) {
          bR[i] -= bR[j]*x.real() + bI[j]*x.imag();
          bI[i] -= bI[j]*x.real() - bR[j]*x.imag();
        }
        bR[i] *= idiagonal[i];
        bI[i] *= idiagonal[i];
      }
    };
};



#endif /* INCOMPLETE_QR_COMPLEX_HPP_ */
