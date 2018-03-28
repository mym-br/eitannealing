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
#include "util/heap_siftdown.hpp"
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
        std::unordered_map<unsigned long, scalar> buildingL;
        std::unordered_map<unsigned long, scalar> buildingQ;
        std::vector<i_c> selectedL;
        std::vector<i_c> selectedQ;

    public:
        SparseIncompleteQRBuilder(){};

        Eigen::SparseMatrix<scalar, Eigen::ColMajor>
        buldRMatrix(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &a, unsigned long nl, unsigned long nq)
        {
                unsigned long m = a.rows();
                unsigned long n = a.cols();
                Eigen::SparseMatrix<scalar, Eigen::ColMajor> RMatrix(n, n);
                RMatrix.reserve(Eigen::VectorXi::Constant(n, nl));
                for(auto &x : this->qrows) {
                    x.clear();
                    x.reserve(2*nq);  // That's actually somehow unpredictable, but 2*nq should be enought
                }
                for(auto &x : this->qcols) {
                    x.clear();
                    x.reserve(nq);  // That's actually somehow unpredictable, but 2*nq should be enought
                }

                for(unsigned long j = 0; j<n ; j++) {
                    // Calc the current L vector, product of i-th row of a
                    //  and the previous q matrix
                    buildingL.clear();
                        
                    for(unsigned long i = 0; i < j; i++) {
                        for(typename Eigen::SparseMatrix<scalar>::InnerIterator it(a,i); it; ++it) {
                            for(auto [k, v] : qrows[it.index()]) {
                                buildingL[k] += v * it.value();
                            }
                            // Get nl elements
                            selectedL.clear();
                            auto ll = buildingL.begin();
                            unsigned long count = 0;
                            while(count < nl && ll != buildingL.end()) {
                                selectedL.push_back(std::pair(ll->first, ll->second));
                                count++;
                                ll++;
                            }
                            if(ll != buildingL.end()) { // There are remaining elements, get the nl largest
                                auto cmp = [](const i_c &a, i_c const &b) { return a.second > b.second; };
                                std::make_heap(selectedL.begin(), selectedL.end(), cmp); // min heap
                                while(ll != buildingL.end()) {
                                    if(ll->second > selectedL[0].second) {
                                        // Replace smallest element and fix the heap
                                        selectedL[0].first = ll->first;
                                        selectedL[0].second = ll->second;
                                        heap_sift_top_down(selectedL.begin(), selectedL.end(), cmp);
                                    }
                                    ll++;
                                }
                        }
                    }
                }

        };


};

#endif  // INCOMPLETE_LQ_BUILDER_H_
