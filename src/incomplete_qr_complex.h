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
    std::vector<std::unordered_map<unsigned long, double> > rowsR;
    std::vector<std::unordered_map<unsigned long, double> > rowsI;
    std::vector<std::unordered_map<unsigned long, double> > colsR;
    std::vector<std::unordered_map<unsigned long, double> > colsI;

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
    ): idiagonal(iiUpper_R.rows()), rowsR(iiUpper_R.rows()),
      rowsI(iiUpper_R.rows()), colsR(iiUpper_R.rows()), colsI(iiUpper_R.rows()) {
      idiagonal.reserve(iiUpper_R.size());
      SparseIncompleteQRBuilder<std::complex<double> > builder;


      builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(iiUpper_R, iiUpper_I, Aic_R, Aic_I), 8, 16, [this](unsigned long j, double x) {
        this->idiagonal[j] = x;
      }, [this](unsigned long i, unsigned long j, std::complex<double> x) {
        this->rowsR[i][j] = x.real();
        this->rowsI[i][j] = x.imag();
        this->colsR[j][i] = x.real();
        this->colsI[j][i] = x.imag();
      });
    }

    bool solveInPlaceC(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const;
    bool solveInPlaceCT(Eigen::VectorXd &bR, Eigen::VectorXd &bI) const;
};



#endif /* INCOMPLETE_QR_COMPLEX_HPP_ */
