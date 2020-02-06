/*
 * incomplete_cholesky.hpp
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_QR_COMPLEX2_HPP_
#define INCOMPLETE_QR_COMPLEX2_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <unordered_map>
#include <iostream>

#include "basematrix.h"
#include "incomplete_qr_builder.h"

class SparseIncompleteQRComplex2
{
  protected:
    vectorx idiagonal;
    matrixcomplex rmatrix;

    struct MatricesStorageAdaptor {
        const matrixcomplex &ii;
        const matrixcomplex &ic;
        unsigned long square_size;
        // FIXME: This map can be reused!
        std::vector<std::vector<std::pair<unsigned long, std::complex<double> > > > upperMap;
        MatricesStorageAdaptor(const matrixcomplex &Aii_low, const matrixcomplex &Aic):
             ii(Aii_low), ic(Aic), square_size(Aii_low.cols()), upperMap(Aii_low.rows())  {}
        void iterateOverColumn(unsigned long j, std::function<void(unsigned long, std::complex<double>)> &&f) {
            // FIXME: This assumes that the columns will be iterated in order, so the current column has already
            //  an appropriate upper map
            // First iterate over upper elements
            for(auto [j, x] : upperMap[j]) f(j, x);
            // Now, lower elements and fill map for next columns
            for(typename matrixcomplex::InnerIterator it(ii, j); it; ++it) {
                if(it.index()>j) {
                  upperMap[it.index()].push_back(std::make_pair(j, it.value()));
                }
                f(it.index(), it.value());
            }
            // Finally, ic elements
            for(typename matrixcomplex::InnerIterator it(ic, j); it; ++it) {
                f(it.index() + square_size, it.value());
            }
        }
        unsigned long rows() const { return square_size+ic.rows(); }
        unsigned long cols() const { return square_size; }
    };

  public:

    SparseIncompleteQRComplex2(unsigned long nr, unsigned long nq, const matrixcomplex &Aii_low, const matrixcomplex &Aic):
        idiagonal(Aii_low.rows()), rmatrix(Aii_low.rows(), Aii_low.rows()) {

        SparseIncompleteQRBuilder<std::complex<double> > builder;

        this->rmatrix.reserve(Aii_low.rows()*nr);

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, double x) {
            this->idiagonal(j) = x;
         },
         [this](unsigned long i, unsigned long j, std::complex<double> x) {
             this->rmatrix.insert(i,j) = x;
         });
        rmatrix.makeCompressed();
        idiagonal = idiagonal.cwiseInverse();
        for(int i = 0; i<rmatrix.outerSize(); i++)
            rmatrix.col(i) *= idiagonal(i);
    }

    void solveInPlaceC(vectorxcomplex &b) const {
      rmatrix.triangularView<Eigen::UnitUpper>().solveInPlace(b);
      b = b.cwiseProduct(idiagonal);
    }

     // conjugated transpose
    void solveInPlaceCT(vectorxcomplex &b) const {
        b = b.cwiseProduct(idiagonal);
        rmatrix.triangularView<Eigen::UnitUpper>().transpose().solveInPlace(b);
    }
};



#endif /* INCOMPLETE_QR_COMPLEX2_HPP_ */
