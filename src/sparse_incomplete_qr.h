#ifndef _SPARSE_INCOMPLETE_QR_H_
#define _SPARSE_INCOMPLETE_QR_H_

#include "incomplete_qr_builder.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>


template <class scalar> class SparseIncompleteQR
{
  protected:
    typedef Eigen::SparseMatrix<scalar, Eigen::ColMajor> basematrix;
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> basevector;
    typedef typename Eigen::NumTraits<scalar>::Real real;
    // idiagonal is always real, even when scalar is complex!
    Eigen::Matrix<real, Eigen::Dynamic, 1> idiagonal;
    basematrix rmatrix;

    struct MatricesStorageAdaptor {
        
        const basematrix &ii;
        const basematrix &ic;
        unsigned long square_size;
        // FIXME: This map can be reused!
        std::vector<std::vector<std::pair<unsigned long, scalar > > > upperMap;
        MatricesStorageAdaptor(const basematrix &Aii_low, const basematrix &Aic):
             ii(Aii_low), ic(Aic), square_size((unsigned long)Aii_low.cols()), upperMap(Aii_low.rows())  {}
        void iterateOverColumn(unsigned long j, std::function<void(unsigned long, scalar)> &&f) {
            // FIXME: This assumes that the columns will be iterated in order, so the current column has already
            //  an appropriate upper map
            // First iterate over upper elements
            for(auto [j, x] : upperMap[j]) f(j, x);
            // Now, lower elements and fill map for next columns
            for(typename basematrix::InnerIterator it(ii, j); it; ++it) {
                if(it.index()>j) {
                  upperMap[it.index()].push_back(std::make_pair(j, it.value()));
                }
                f(it.index(), it.value());
            }
            // Finally, ic elements
            for(typename basematrix::InnerIterator it(ic, j); it; ++it) {
                f(it.index() + square_size, it.value());
            }
        }
        unsigned long rows() const { return square_size+(unsigned long)ic.rows(); }
        unsigned long cols() const { return square_size; }
    };



  public:

    SparseIncompleteQR(unsigned long nr, unsigned long nq, const basematrix &Aii_low, const basematrix &Aic):
        idiagonal(Aii_low.rows()), rmatrix(Aii_low.rows(), Aii_low.rows()) {

        SparseIncompleteQRBuilder<scalar> builder;

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, real invx) {
            this->idiagonal(j) = invx;
         },
         [this](unsigned long i, unsigned long j, scalar x) {
             this->rmatrix.insert(i,j) = x;
         });
         rmatrix.makeCompressed();
         idiagonal = idiagonal.cwiseInverse();
         for(int j = 1; j<rmatrix.outerSize(); j++) {
             rmatrix.col(j) *= idiagonal(j);
         }
    }

    SparseIncompleteQR(unsigned long nr, unsigned long nq, const basematrix &Aii_low, const basematrix &Aic,  SparseIncompleteQRBuilder<double> &builder):
        idiagonal(Aii_low.rows()), rmatrix(Aii_low.rows(), Aii_low.rows()) {

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, real invx) {
            this->idiagonal(j) = invx;
         },
         [this](unsigned long i, unsigned long j, scalar x) {
             this->rmatrix.insert(i,j) = x;
         });
         rmatrix.makeCompressed();
         for(int j = 1; j<rmatrix.outerSize(); j++) {
             rmatrix.col(j) *= idiagonal(j);
         }
    }

    void solveInPlace(basevector &b) const {
        //FIXME: On gcc without fast-math complex x complex product is very slow (Something to do with Inf checks)
        rmatrix.template triangularView<Eigen::UnitUpper>().solveInPlace(b);
        b = b.cwiseProduct(idiagonal);
    }

     // conjugated transpose
    void solveInPlaceT(basevector &b) const {
        b = b.cwiseProduct(idiagonal);
        rmatrix.template triangularView<Eigen::UnitUpper>().transpose().solveInPlace(b);
    }
};

#endif
