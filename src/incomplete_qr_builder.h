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
#include <unordered_map>
#include <algorithm>
#include "util/fill_with_smallest.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

template<class scalar> class SparseIncompleteQRBuilder
{
    protected:

        typedef std::pair<unsigned long, scalar> i_c; // Pair index coefficient

        // Those persistent objects should provide storage for
        // temporary results, so that new calls to buildLMatrix()
        // won't incur into memory allocation
        std::vector<std::vector<i_c> > qrows;
        std::vector<std::vector<i_c> > qcols;
        std::unordered_map<unsigned long, scalar> buildingR;
        std::unordered_map<unsigned long, scalar> buildingQ;
        std::vector<i_c> selectedR;
        std::vector<i_c> selectedQ;

    public:
        SparseIncompleteQRBuilder(){};

        Eigen::SparseMatrix<scalar, Eigen::ColMajor>
        buldRMatrix(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &a, unsigned long nr, unsigned long nq)
        {
                unsigned long m = a.rows();
                unsigned long n = a.cols();
                Eigen::SparseMatrix<scalar, Eigen::ColMajor> RMatrix(n, n);
                struct _calc_q_cols_size { 
                    unsigned long mr; _calc_q_cols_size(unsigned long n):mr(n){}; 
                    typedef unsigned long value_type;
                    value_type operator[](unsigned long i) const { return (i+1)>mr?mr:(i+1); }
                };
                RMatrix.reserve(_calc_q_cols_size(nr));
                // Initialize temp q storage
                qcols.resize(n);
                for(auto &x : this->qcols) {
                    x.clear();
                    x.reserve(nq); 
                }
                qrows.resize(m);
                for(auto &x : this->qrows) {
                    x.clear();
                    x.reserve(2*nq);  // That's actually somehow unpredictable, but 2*nq should be enought
                }
                for(unsigned long j = 0; j<n ; j++) {
                    // Calc the current L vector, product of i-th row of a
                    //  and the previous q matrix
                    buildingR.clear();           
                    buildingQ.clear();  // As long as we'll be iterating over A[:j] column, initialize next Q vector
                    // The following should have ~O(nq^2) complexity
                    for(typename Eigen::SparseMatrix<scalar>::InnerIterator it(a,j); it; ++it) {
                        buildingQ[it.index()] = it.value(); // Used on the 2nd step
                        for(auto [k, v] : qrows[it.index()]) {
                            buildingR[k] += v * it.value();
                        }
                    }
                    auto cmp_larger_abs_coef = [](const i_c &a, i_c const &b) {return std::abs(a.second) > std::abs(b.second);};
                    // Get nr *largest* elements, notice the reversed comparator above
                    fillWithNSmallest(selectedR, buildingR, nr, cmp_larger_abs_coef);
                    // Sort it according to index
                    std::sort(selectedR.begin(), selectedR.end(), [](const i_c &a, const i_c &b){return a.first<b.first;});
                    // Now fill R matrix column and finalize Q calculation
                    for(auto [ri, rv] : selectedR) {
                        RMatrix.insert(ri, j) = rv;
                        for(auto [qj, qv] : qcols[ri])
                            buildingQ[qj] -= rv*qv;
                    }
                    // Get ql largest elements, same procedure as for L above
                    fillWithNSmallest(selectedQ, buildingQ, nq, cmp_larger_abs_coef);
                    // Now update q storage
                    for(auto [i, v] : selectedQ) {
                        qrows[i].push_back(std::pair(j, v));
                        qcols[j].push_back(std::pair(i, v));
                    }
                }
                // Finish
                RMatrix.makeCompressed();
                return RMatrix;
        };


};

#endif  // INCOMPLETE_LQ_BUILDER_H_
