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
#include <iostream>

#include "basematrix.h"
#include "incomplete_qr_builder.h"

class SparseIncompleteQRComplex
{
  protected:
    std::vector<double> idiagonal;
    std::vector<std::vector<std::pair<unsigned long, std::complex<double> > > > rows;
    std::vector<std::vector<std::pair<unsigned long, std::complex<double> > > > cols;


    struct MatricesStorageAdaptor {
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &iiR;
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &iiI;
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &icR;
        const Eigen::SparseMatrix<double, Eigen::ColMajor> &icI;
        unsigned long square_size;
        // FIXME: This map can be reused!
        std::vector<std::vector<std::pair<unsigned long, std::complex<double> > > > upperMap;

        MatricesStorageAdaptor(
          const matrix &iiUpper_R, const matrix &iiUpper_I,
          const matrix &Aic_R, const matrix &Aic_I):
          iiR(iiUpper_R), iiI(iiUpper_I), icR(Aic_R), icI(Aic_I), square_size((unsigned long)iiUpper_R.cols()), upperMap(iiUpper_R.rows())
           {}
        void iterateOverColumn(unsigned long j, std::function<void(unsigned long, std::complex<double>)> &&f) {
            // FIXME: This assumes that the columns will be iterated in order, so the current column has already
            //  an appropriate upper map
            // First iterate over upper elements
            for(auto [i, x] : upperMap[j]) f(i, x);
            // Now, lower elements and fill map for next columns
            for(typename Eigen::SparseMatrix<double>::InnerIterator itR(iiR, j), itI(iiI, j); itR; ++itR, ++itI) {
                std::complex<double> val(itR.value(), itI.value());
                if((unsigned)itR.index()>j) {
                  upperMap[itR.index()].push_back(std::pair<unsigned long, std::complex<double> >(j, val));
                }
                f(itR.index(), val);
            }
            // Finally, ic elements
            for(typename Eigen::SparseMatrix<double>::InnerIterator itR(icR, j), itI(icI, j); itR; ++itR, ++itI) {
                std::complex<double> val(itR.value(), itI.value());
                f(itR.index() + square_size, val);
            }
        }
        unsigned long rows() const { return square_size + (unsigned long)icR.rows(); }
        unsigned long cols() const { return (unsigned long)square_size; }
    };

  public:

    SparseIncompleteQRComplex(unsigned long nr, unsigned long nq,
      const matrix &iiUpper_R, const matrix &iiUpper_I,
      const matrix &Aic_R, const matrix &Aic_I
    ): idiagonal(iiUpper_R.rows()), rows(iiUpper_R.rows()), cols(iiUpper_R.rows()) {

      SparseIncompleteQRBuilder<std::complex<double> > builder;

      builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(iiUpper_R, iiUpper_I, Aic_R, Aic_I), nr, nq, [this](unsigned long j, double invx) {
        this->idiagonal[j] = invx;
      }, [this](unsigned long i, unsigned long j, std::complex<double> x) {
        this->rows[i].push_back(std::make_pair(j, x));
        this->cols[j].push_back(std::make_pair(i, x));
      });
    }

    void solveInPlaceC(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const {
      int i = (int)(idiagonal.size() - 1);
      for(; i>=0; i--) {
        double vr = bR[i];
        double vi = bI[i];
        for(auto [j, x] : rows[i]) {
          vr -= bR[j]*x.real() - bI[j]*x.imag();
          vi -= bR[j]*x.imag() + bI[j]*x.real();
        }
        vr *= idiagonal[i];
        vi *= idiagonal[i];
        bR[i] = vr;
        bI[i] = vi;
      }
    }

     // conjugated transpose
    void solveInPlaceCT(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const {
	  size_t n = idiagonal.size();
      for(unsigned int i = 0; i<n; i++) {
        double vr = bR[i];
        double vi = bI[i];
        for(auto [j, x] : cols[i]) {
          vr -= bR[j]*x.real() + bI[j]*x.imag();
          vi -= bI[j]*x.real() - bR[j]*x.imag();
        }
        vr *= idiagonal[i];
        vi *= idiagonal[i];
        bR[i] = vr;
        bI[i] = vi;
      }
    };
};



#endif /* INCOMPLETE_QR_COMPLEX_HPP_ */
