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
    // FIXME: NO! idiagonal is always real, even when scalar is complex!
    Eigen::Matrix<double, Eigen::Dynamic, 1> idiagonal;
    basematrix rmatrix;
    std::vector<std::vector<std::pair<unsigned, scalar > > > cols;
    std::vector<std::vector<std::pair<unsigned, scalar > > > rows;

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
        idiagonal(Aii_low.rows()), rmatrix(Aii_low.rows(), Aii_low.rows()), cols(Aii_low.cols()), rows(Aii_low.cols()) {

        SparseIncompleteQRBuilder<scalar> builder;

        this->rmatrix.reserve(Aii_low.rows()*nr);

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, double x) {
            this->idiagonal(j) = x;
         },
         [this](unsigned long i, unsigned long j, scalar x) {
             this->rmatrix.insert(i,j) = x;
             this->cols[j].push_back(std::pair<long, scalar>(i, x));
             this->rows[i].push_back(std::pair<long, scalar>(j, x));
         });
        rmatrix.makeCompressed();
        idiagonal = idiagonal.cwiseInverse();
        for(int j = 1; j<rmatrix.outerSize(); j++) {
            rmatrix.col(j) *= idiagonal(j);
            for(auto &x : cols[j])
                x.second *= idiagonal(j);
        }
    }

    void solveInPlace(basevector &b) const {
        //FIXME: For some reason, triangularView::solveInPlace is slower than this?
        for(auto j = cols.size()-1; j >0; j--) { // 1st column is empty
            scalar coef = b[j];
            for(auto [i, x] : cols[j]) {
                b[i] -= coef*x;
            }
        }
        //this->rmatrix.template triangularView<Eigen::UnitUpper>().solveInPlace(b);
        b = b.cwiseProduct(this->idiagonal);
        /*for(auto i = rows.size(); i > 0; i--) {
            scalar tot = b[i-1];
            for(auto [j, x] : rows[i-1]) {
                tot -= b[j]*x;
            }
            b[i-1] = tot*idiagonal[i-1];
        }
        int i = (int)(idiagonal.size() - 1);
        for(; i>=0; i--) {
            scalar v = b[i];
            for(auto [j, x] : rows[i]) {
                v -= b[j]*x;
            }
            v *= idiagonal[i];
            b[i] = v;
        }*/
    }

     // conjugated transpose
    void solveInPlaceT(basevector &b) const {
        b = b.cwiseProduct(idiagonal);
        rmatrix.template triangularView<Eigen::UnitUpper>().transpose().solveInPlace(b);
    }
};

#endif
