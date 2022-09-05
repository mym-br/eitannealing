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
template<class i, class v> using SparseIncompleteQRBuilder_map = absl::flat_hash_map<i, v>;
#else
#include <unordered_map>
template<class i, class v> using SparseIncompleteQRBuilder_map = std::unordered_map<i, v>;
#endif


/* Scalar-specific operations
 * Handles special cases for complex numbers,
 * such as product by complex conjugate and cache of squared norm
 */
template<class scalar> struct scalar_ops;

template<> struct scalar_ops<double> {
    typedef std::pair<unsigned long, double> i_c;
    typedef typename Eigen::NumTraits<double>::Real real;

    static inline double inner(const double &a, const double &b) {
      return a*b;
    }

    class coefficient_vector {
          std::vector<i_c> vector;
    public:
        void reserve(unsigned n) {
            vector.reserve(n);
        }

        inline void fill_with_selected_coeficients(const SparseIncompleteQRBuilder_map<unsigned long, double> &building, unsigned n) {
            vector.clear();
            fillWithNSmallest(vector, building, n, [](const i_c &a, const i_c &b){return std::abs(a.second) > std::abs(b.second);});
        }

        inline real inv_total_coefficients_norm() {
            real qnorm2 = 0;
            for(auto [i, v] : vector) qnorm2 += v*v;
            real qnorm = std::sqrt(qnorm2);
            return 1/qnorm;
        }

        inline void sort_by_index() {
            std::sort(vector.begin(), vector.end(), [](const i_c &a, const i_c &b){return a.first<b.first;});
        }
        template <class func> inline void iterate_over_coefficients(func f) {
            for(auto [ri, rv] : vector) f(ri, rv);
        }
    };
};

template<> struct scalar_ops<std::complex<double> > {
    typedef std::pair<unsigned long, std::complex<double> > i_c;
    typedef typename Eigen::NumTraits<double>::Real real;

    static inline std::complex<double> inner(const std::complex<double> &a, const std::complex<double> &b) {
        return a*std::conj(b);
    }

    class coefficient_vector {
          struct coefficient_with_cached_norm {
              unsigned long i;
              std::complex<double> v;
              real n2;  // cached squared norm

              // Construct from index/value pair
              coefficient_with_cached_norm(const std::pair<unsigned long, std::complex<double> > &x):
                i(x.first), v(x.second), n2(std::norm(x.second)) {}
          };

          std::vector<coefficient_with_cached_norm> vector;
    public:
        void reserve(unsigned n) {
            vector.reserve(n);
        }

        inline void fill_with_selected_coeficients(const SparseIncompleteQRBuilder_map<unsigned long, std::complex<double>> &building, unsigned n) {
            vector.clear();
            fillWithNSmallest(vector, building, n, [](const coefficient_with_cached_norm &a, const coefficient_with_cached_norm &b) {
                return a.n2 > b.n2;
            });
        }

        inline real inv_total_coefficients_norm() {
            real qnorm2 = 0;
            for(auto [i, v, n2] : vector) qnorm2 += n2;
            real qnorm = std::sqrt(qnorm2);
            return 1/qnorm;
        }

        inline void sort_by_index() {
            std::sort(vector.begin(), vector.end(), [](const coefficient_with_cached_norm &a, const coefficient_with_cached_norm &b){return a.i<b.i;});
        }
        template <class func> inline void iterate_over_coefficients(func f) {
            for(auto [ri, rv, n2] : vector) f(ri, rv);
        }
    };
};


template<class scalar> class SparseIncompleteQRBuilder
{
    protected:

        typedef scalar_ops<scalar> ops;
        typedef typename ops::real real;
        typedef typename ops::i_c i_c; // Pair index coefficient
        typedef typename ops::coefficient_vector selected_coefficients_vector;
        typedef SparseIncompleteQRBuilder_map<unsigned long, scalar> building_coefficients_map;

        // Those persistent objects should provide storage for
        // temporary results, so that new calls to buildLMatrix()
        // won't incur into memory allocation
        std::vector<std::vector<i_c> > qrows;
        std::vector<std::vector<i_c> > qcols;

        building_coefficients_map buildingR;
        building_coefficients_map buildingQ;

        selected_coefficients_vector selectedR;
        selected_coefficients_vector selectedQ;

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
                            this->buildingR[qj] += ops::inner(v, qv);
                    });
                    // Get nr-1 *largest* elements
                    //  -1 accounts for the diagonal
                    selectedR.fill_with_selected_coeficients(buildingR, nr-1);
                    // Sort it according to index
                    selectedR.sort_by_index();
                    // Now fill R matrix column and finalize Q calculation
                    selectedR.iterate_over_coefficients([&](unsigned long ri, scalar rv) {
                         insert_upper(ri, j, rv);
                         for(auto [qj, qv] : qcols[ri])
                              buildingQ[qj] -= rv*qv;
                    });
                    // Get ql largest elements, same procedure as for L above
                    selectedQ.fill_with_selected_coeficients(buildingQ, nq);
                    // Renormalize
                    real inorm = selectedQ.inv_total_coefficients_norm();
                    // Final element of R is the norm
                    insert_diagonal(j, inorm);
                    // Now update q storage
                    selectedQ.iterate_over_coefficients([&](unsigned long i, scalar v) {
                        scalar nv = v*inorm;
                        qrows[i].emplace_back(j, nv);
                        qcols[j].emplace_back(i, nv);
                    });
                }
        };
        struct columnMajorStorageAdaptor {
            const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m;
            columnMajorStorageAdaptor(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m):m(m){}
            template <class func> void iterateOverColumn(unsigned long j, func &&f) const {
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
};


#endif  // INCOMPLETE_LQ_BUILDER_H_
