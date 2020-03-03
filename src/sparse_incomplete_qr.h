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
    basevector idiagonal;
    basematrix rmatrix;

    struct MatricesStorageAdaptor {
        
        const basematrix &ii;
        const basematrix &ic;
        unsigned long square_size;
        // FIXME: This map can be reused!
        std::vector<std::vector<std::pair<unsigned long, scalar > > > upperMap;
        MatricesStorageAdaptor(const basematrix &Aii_low, const basematrix &Aic):
             ii(Aii_low), ic(Aic), square_size(Aii_low.cols()), upperMap(Aii_low.rows())  {}
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
        unsigned long rows() const { return square_size+ic.rows(); }
        unsigned long cols() const { return square_size; }
    };

  public:

    SparseIncompleteQR(unsigned long nr, unsigned long nq, const basematrix &Aii_low, const basematrix &Aic):
        idiagonal(Aii_low.rows()), rmatrix(Aii_low.rows(), Aii_low.rows()) {

        SparseIncompleteQRBuilder<scalar> builder;

        this->rmatrix.reserve(Aii_low.rows()*nr);

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, scalar x) {
            this->idiagonal(j) = x;
         },
         [this](unsigned long i, unsigned long j, scalar x) {
             this->rmatrix.insert(i,j) = x;
         });
        rmatrix.makeCompressed();
        idiagonal = idiagonal.cwiseInverse();
        for(int i = 0; i<rmatrix.outerSize(); i++)
            rmatrix.col(i) *= idiagonal(i);
    }

    void solveInPlace(basevector &b) const {
      this->rmatrix.template triangularView<Eigen::UnitUpper>().solveInPlace(b);
      b = b.cwiseProduct(this->idiagonal);
    }

     // conjugated transpose
    void solveInPlaceT(basevector &b) const {
        b = b.cwiseProduct(idiagonal);
        rmatrix.template triangularView<Eigen::UnitUpper>().transpose().solveInPlace(b);
    }
};

#endif