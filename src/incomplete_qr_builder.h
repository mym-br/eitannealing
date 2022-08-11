/*
 * incomplete_qr_builder.h
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_QR_BUILDER_H_
#define INCOMPLETE_QR_BUILDER_H_

#include <utility>
#include <vector>
#include <algorithm>
#include <complex>
#include "util/fill_with_smallest.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "util/eigen_sizestype_adaptor.hpp"

// ABSL Flat hash seems to be faster than std::unordered_map
#ifdef USE_ABSL_FLAT_HASH
#include <absl/container/flat_hash_map.h>
#define INCOMPLETEQRBUILDER_MAP absl::flat_hash_map
#else
#include <unordered_map>
#define INCOMPLETEQRBUILDER_MAP std::unordered_map
#endif

template<class scalar> class SparseIncompleteQRBuilder
{
    protected:
        typedef typename Eigen::NumTraits<scalar>::Real real;
        typedef std::pair<unsigned long, scalar> i_c; // Pair index coefficient

        // Those persistent objects should provide storage for
        // temporary results, so that new calls to buildLMatrix()
        // won't incur into memory allocation
        std::vector<std::vector<i_c> > qrows;
        std::vector<std::vector<i_c> > qcols;
        INCOMPLETEQRBUILDER_MAP<unsigned long, scalar> buildingR;
        INCOMPLETEQRBUILDER_MAP<unsigned long, scalar> buildingQ;
        std::vector<i_c> selectedR;
        std::vector<i_c> selectedQ;

        static real sqnorm(const scalar &x);
        static scalar inner(const scalar &x, const scalar &y);
        static inline bool cmp_larger_abs_coef(const i_c &a, const i_c &b);

    public:
        SparseIncompleteQRBuilder(){};
        // This works for a generic "columnMajorStorage" concept.
        // This concept must implement the methods:
        //  iterateOverColumn(unsigned long i, std::function<void(unsigned long, scalar)> &f) const
        //    applies f to each nonzero element of a's i-th column, passing its row number and value.
        // Unsigned long rows() and unsigned long cols()
        //  Attention: diagonalInsertFunction receives the *Reciprocal* diagonal value
        template <class columnMajorStorage, class diagonalInsertFunction, class upperElementsInsertFunction> void
        buildRMatrixFromColStorage(columnMajorStorage &&a, unsigned long nr, unsigned long nq, diagonalInsertFunction &&insert_diagonal, upperElementsInsertFunction &&insert_upper)
        {
                unsigned long m = a.rows();
                unsigned long n = a.cols();

                // Initialize temp q storage
                qcols.resize(n);
                for(auto &x : this->qcols) {
                    x.clear();
                    x.reserve(nq);
                }
                qrows.resize(m);
                for(auto &x : this->qrows) {
                    x.clear();
                    x.reserve(3*nq);  // That's actually somehow unpredictable, but 3*nq should be enought
                }
                selectedR.reserve(nr-1);
                selectedQ.reserve(nq);
                for(unsigned long j = 0; j<n ; j++) {
                    // Calc the current L vector, product of i-th row of a
                    //  and the previous q matrix
                    buildingR.clear();
                    buildingQ.clear();  // As long as we'll be iterating over A[:j] column, initialize next Q vector
                    // The following should have ~O(nq^2) complexity
                    a.iterateOverColumn(j,[this](unsigned long i, scalar v){
                        this->buildingQ[i] = v;
                        for(auto [qj, qv] : qrows[i])
                            this->buildingR[qj] += inner(v, qv);
                    });
                    // Get nr-1 *largest* elements
                    //  -1 accounts for the diagonal
                    selectedR.clear();
                    fillWithNSmallest(selectedR, buildingR, nr - 1, cmp_larger_abs_coef);
                    // Sort it according to index
                    std::sort(selectedR.begin(), selectedR.end(), [](const i_c &a, const i_c &b){return a.first<b.first;});
                    // Now fill R matrix column and finalize Q calculation
                    for(auto [ri, rv] : selectedR) {
                        insert_upper(ri, j, rv);
                        for(auto [qj, qv] : qcols[ri])
                            buildingQ[qj] -= rv*qv;
                    }
                    // Get ql largest elements, same procedure as for L above
                    selectedQ.clear();
                    fillWithNSmallest(selectedQ, buildingQ, nq, cmp_larger_abs_coef);
                    // Renormalize
                    real qnorm2 = 0;
                    for(auto [i, v] : selectedQ) qnorm2 += sqnorm(v);
                    real qnorm = std::sqrt(qnorm2);
                    real inorm = 1/qnorm;
                    // Final element of R is the norm
                    insert_diagonal(j, inorm);
                    // Now update q storage
                    for(auto [i, v] : selectedQ) {
                        scalar nv = v*inorm;
                        qrows[i].emplace_back(std::move(std::pair(j, nv)));
                        qcols[j].emplace_back(std::move(std::pair(i, nv)));
                    }
                }
        };
        struct columnMajorStorageAdaptor {
            const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m;
            columnMajorStorageAdaptor(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m):m(m){}
            void iterateOverColumn(unsigned long j, std::function<void(unsigned long, scalar)> &&f) const {
                for(typename Eigen::SparseMatrix<scalar>::InnerIterator it(m, j); it; ++it)
                    f(it.index(), it.value());
            }
            unsigned long rows() const { return m.rows(); }
            unsigned long cols() const { return m.cols(); }
        };
        template <class diagonalInsertFunction, class upperElementsInsertFunction> void
        buildRMatrix(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &a, unsigned long nr, unsigned long nq, diagonalInsertFunction &&insert_diagonal, upperElementsInsertFunction &&insert_upper) {
            buildRMatrixFromColStorage(columnMajorStorageAdaptor(a), nr, nq, insert_diagonal, insert_upper);
        }

        typedef decltype(sqnorm) norm_type;
};

template<> inline double SparseIncompleteQRBuilder<double>::sqnorm(const double &x) {
  return x*x;
}

template<> inline double SparseIncompleteQRBuilder<double>::inner(const double &x, const double &y) {
  return x*y;
}

template<> inline bool SparseIncompleteQRBuilder<double>::cmp_larger_abs_coef(
  const SparseIncompleteQRBuilder<double>::i_c &a,
  const SparseIncompleteQRBuilder<double>::i_c &b) {
  return std::abs(a.second) > std::abs(b.second);
}

template<> inline double SparseIncompleteQRBuilder<typename std::complex<double> >::sqnorm(const std::complex<double> &x) {
  return std::norm(x);
}

template<> inline std::complex<double> SparseIncompleteQRBuilder<std::complex<double> >::inner(const std::complex<double> &x, const std::complex<double> &y) {
  return x*std::conj(y);
}
// FIXME: Can we cache results of l2 squared norm?
template<> inline bool SparseIncompleteQRBuilder<std::complex<double> >::cmp_larger_abs_coef(
  const SparseIncompleteQRBuilder<std::complex<double> >::i_c &a,
  const SparseIncompleteQRBuilder<std::complex<double> >::i_c &b) {
  return std::norm(a.second) > std::norm(b.second);
}


#endif  // INCOMPLETE_LQ_BUILDER_H_
