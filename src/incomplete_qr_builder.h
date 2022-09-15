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
    typedef double scalar;

    static inline double inner(const double &a, const double &b) {
      return a*b;
    }

    class coefficient {
    protected:
        unsigned long i;
        double val;
    public:
        coefficient(const std::pair<unsigned long, scalar> &p):
            i(p.first), val(p.second) {}

        unsigned long index() const { return i; }
        double value() const { return val; }
        real sqnorm() const { return val*val; }
        bool norm_greater_than(const coefficient &other) const {
            return std::abs(val) > std::abs(other.val);
        }
    };
};

template<> struct scalar_ops<std::complex<double> > {
    typedef std::pair<unsigned long, std::complex<double> > i_c;
    typedef typename Eigen::NumTraits<double>::Real real;
    typedef std::complex<double> scalar;

    static inline std::complex<double> inner(const std::complex<double> &a, const std::complex<double> &b) {
        return a*std::conj(b);
    }

    class coefficient {
    protected:
        unsigned long i;
        std::complex<double> val;
        real n2;    // Cached squared norm
    public:
        coefficient(const std::pair<unsigned long, scalar> &p):
            i(p.first), val(p.second), n2(std::norm(p.second)) {}

        unsigned long index() const { return i; }
        std::complex<double> value() const { return val; }
        real sqnorm() const { return n2; }
        bool norm_greater_than(const coefficient &other) const {
            return n2 > other.n2;
        }
    };
};


template<class scalar> class SparseIncompleteQRBuilder
{
    protected:

        typedef scalar_ops<scalar> ops;
        typedef typename ops::real real;
        typedef typename ops::i_c i_c; // Pair index coefficient
        typedef std::vector<typename ops::coefficient> selected_coefficients_vector;

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
                    selectedR.clear();
                    fillWithNSmallest(selectedR, buildingR, nr-1, [](const auto &a, const auto &b){return a.norm_greater_than(b);});
                    // Sort it according to index
                    std::sort(selectedR.begin(), selectedR.end(), [](const auto &a, const auto &b){return a.index()<b.index();});
                    // Now fill R matrix column and finalize Q calculation
                    for(auto c : selectedR) {
                        insert_upper(c.index(), j, c.value());
                        for(auto [qj, qv] : qcols[c.index()])
                            buildingQ[qj] -= c.value()*qv;
                    }
                    // Get ql largest elements, same procedure as for L above
                    selectedQ.clear();
                    fillWithNSmallest(selectedQ, buildingQ, nq, [](const auto &a, const auto &b){return a.norm_greater_than(b);});
                    // Renormalize
                    real qnorm2 = 0;
                    for(auto c : selectedQ) qnorm2 += c.sqnorm();
                    real inorm = 1/std::sqrt(qnorm2);
                    // Final element of R is the norm
                    insert_diagonal(j, inorm);
                    // Now update q storage
                    for(auto c : selectedQ) {
                        scalar nv = c.value()*inorm;
                        qrows[c.index()].emplace_back(j, nv);
                        qcols[j].emplace_back(c.index(), nv);
                    }
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
