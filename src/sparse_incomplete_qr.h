#ifndef _SPARSE_INCOMPLETE_QR_H_
#define _SPARSE_INCOMPLETE_QR_H_

#include "incomplete_qr_builder.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>

// Conjugate on complex types
template <class base, int iscomplex> struct transpconjimpl;

template <class base> struct transpconjimpl<base, 1> {
    static base get_res(base &&x) {
        return std::conj(x);
    }
};

template <class base> struct transpconjimpl<base, 0> {
    static base get_res(base &&x) {
        return x;
    }
};

template <class base> base transpconj(base &&x) {
    return transpconjimpl<base, Eigen::NumTraits<base>::IsComplex>::get_res(x);
}

template <class scalar> class SparseIncompleteQR
{
  protected:
    typedef Eigen::SparseMatrix<scalar, Eigen::ColMajor> basematrix;
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> basevector;
    typedef typename Eigen::NumTraits<scalar>::Real real;
    // idiagonal is always real, even when scalar is complex!
    Eigen::Matrix<real, Eigen::Dynamic, 1> idiagonal;
    std::vector<std::vector<std::pair<unsigned, scalar > > > rows;
    std::vector<std::vector<std::pair<unsigned, scalar > > > trows;

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
        idiagonal(Aii_low.rows()), rows(Aii_low.cols()), trows(Aii_low.cols()) {

        SparseIncompleteQRBuilder<scalar> builder;

        builder.buildRMatrixFromColStorage(MatricesStorageAdaptor(Aii_low, Aic), nr, nq,
         [this](unsigned long j, real x) {
            this->idiagonal(j) = x;
         },
         [this](unsigned long i, unsigned long j, scalar x) {
             this->trows[j].push_back(std::pair<long, scalar>(i, transpconj(x)));
             this->rows[i].push_back(std::pair<long, scalar>(j, x));
         });
        idiagonal = idiagonal.cwiseInverse();
    }

    void solveInPlace(basevector &b) const {
        //FIXME: For some reason, triangularView::solveInPlace is slower than this?
        for(int i = rows.size()-1; i >= 0; i--) {
            scalar tot = b[i];
            for(auto [j, x] : rows[i]) {
                tot -= b[j]*x;
            }
            b[i] = tot*idiagonal[i];
        }
    }

     // conjugated transpose
    void solveInPlaceT(basevector &b) const {
        for(int i = 0; i<trows.size()-1; i++) {
            scalar tot = b[i];
            for(auto [j, x] : trows[i]) {
                tot -= b[j]*x;
            }
            b[i] = tot*idiagonal[i];
        }
    }
};

// Specialization for complex. For some reason, code is much faster
// by splitting the acumulator into real and imaginary components.
// I'm not really sure why.
template <> void SparseIncompleteQR<std::complex<double> >::solveInPlace(basevector &b) const {
    int i = (int)(idiagonal.size() - 1);
    for(; i>=0; i--) {
        double vr = b[i].real();
        double vi = b[i].imag();
        for(auto [j, x] : rows[i]) {
            vr -= b[j].real()*x.real() - b[j].imag()*x.imag();
            vi -= b[j].real()*x.imag() + b[j].imag()*x.real();
        }
        vr *= idiagonal[i];
        vi *= idiagonal[i];
        b[i] = std::complex(vr, vi);
    }
}

template <> void SparseIncompleteQR<std::complex<double> >::solveInPlaceT(basevector &b) const {
    for(int i = 0; i<trows.size()-1; i++) {
        double vr = b[i].real();
        double vi = b[i].imag();
        for(auto [j, x] : trows[i]) {
            vr -= b[j].real()*x.real() - b[j].imag()*x.imag();
            vi -= b[j].real()*x.imag() + b[j].imag()*x.real();
        }
        vr *= idiagonal[i];
        vi *= idiagonal[i];
        b[i] = std::complex(vr, vi);
    }
}

#endif
