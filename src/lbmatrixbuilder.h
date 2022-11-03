/*
 * lbmatrixbuilder.h
 *
 *  Created on: Nov 1, 2022
 *      Author: Thiago C. Martins
 */

#ifndef LBMATRIXBUILDER_H_
#define LBMATRIXBUILDER_H_

#include <vector>
#include <cstring>
#include <utility>
#include "problem.h"

template<class scalar, class symmetricMatrix, class matrix> class LBMatrixBuilder {
    protected:

        symmetricMatrix aii_skel;
        matrix aic_skel;
        symmetricMatrix acc_skel;

        typedef decltype(std::declval<problem>().getNodeCoefficients()[0]->coefficient) coefficient;

        Eigen::SparseMatrix<coefficient> aii_coef2k;
        Eigen::SparseMatrix<coefficient> aic_coef2k;
        Eigen::SparseMatrix<coefficient> acc_coef2k;

    public:

        LBMatrixBuilder(const problem &p):
            aii_skel(p.getNodesCount()-p.getGenericElectrodesCount(), p.getNodesCount()-p.getGenericElectrodesCount()),
            aic_skel(p.getGenericElectrodesCount(), p.getNodesCount()-p.getGenericElectrodesCount()),
            acc_skel(p.getGenericElectrodesCount(), p.getGenericElectrodesCount())
        {
            int n = p.getNodesCount();
            int numCoefficients = p.getNumCoefficients();
            int iiLimit = p.getNodesCount()-p.getGenericElectrodesCount();
            int i;

            // Phase 1: Initialize skeleton matrices

            // First iiLimit columns: coefficients of aii and aic
            for (i=0; i<iiLimit; ++i) {
                const nodeCoefficients *aux = p.getNodeCoefficients()[i];
                while(aux && aux->node < i) aux = aux->next;  // skip upper triangular
                while(aux && aux->node < iiLimit) { // Insert on aii
                     int row = aux->node;
                     while(aux && aux->node==row)  aux = aux->next; // Combine all coefficients with same row number
                     aii_skel.insert(row, i) = 1.0; // Dummy coefficient
                }
                while(aux) { // Insert on aic
                     int row = aux->node;
                     while(aux && aux->node==row)  aux = aux->next; // Combine all coefficients with same row number
                     aic_skel.insert(row - iiLimit, i) = 1.0; // Dummy coefficient
                }
            }
            // Last columns: coefficients of acc
            for (; i<p.getNodesCount(); ++i) {
                const nodeCoefficients *aux = p.getNodeCoefficients()[i];
                while(aux && aux->node < i) aux = aux->next;  // skip upper triangular
                int row = aux->node;
                while(aux && aux->node==row)  aux = aux->next; // Combine all coefficients with same row number
                acc_skel.insert(row-iiLimit,i-iiLimit) = 1.0;
            }
            aii_skel.makeCompressed();
            aic_skel.makeCompressed();
            acc_skel.makeCompressed();

            // Phase 2: Initialize coefficient2matrix
            scalar *iibase = aii_skel.valuePtr();
            scalar *icbase = aic_skel.valuePtr();
            scalar *ccbase = acc_skel.valuePtr();

            std::vector<Eigen::Triplet<coefficient>> iiTripletList;
            std::vector<Eigen::Triplet<coefficient>> icTripletList;
            std::vector<Eigen::Triplet<coefficient>> ccTripletList;

            // Aii and Aic
            for (i = 0; i<iiLimit; ++i) {
                for (const nodeCoefficients *aux = p.getNodeCoefficients()[i]; aux != NULL; aux = aux->next) {
                    if (aux->node < i) continue; // skip upper triangular
                    if(aux->node<iiLimit) { // aii
                        typename symmetricMatrix::StorageIndex iiIndex = (typename symmetricMatrix::StorageIndex)(&aii_skel.coeffRef(aux->node, i) - iibase);
                        iiTripletList.push_back(Eigen::Triplet<coefficient>(iiIndex, aux->condIndex, aux->coefficient));
                    } else {    // aic
                        typename matrix::StorageIndex icIndex = (typename matrix::StorageIndex)(&aic_skel.coeffRef(aux->node - iiLimit, i) - icbase);
                        icTripletList.push_back(Eigen::Triplet<coefficient>(icIndex, aux->condIndex, aux->coefficient));
                    }
                }
            }
            // Acc
            for (; i<n; ++i) {
                for (const nodeCoefficients *aux = p.getNodeCoefficients()[i]; aux != NULL; aux = aux->next) {
                    if (aux->node < i) continue; // skip upper triangular
                    typename symmetricMatrix::StorageIndex ccIndex = (typename symmetricMatrix::StorageIndex)(&acc_skel.coeffRef(aux->node - iiLimit, i - iiLimit) - ccbase);
                    ccTripletList.push_back(Eigen::Triplet<coefficient>(ccIndex, aux->condIndex, aux->coefficient));
                }
            }

            aii_coef2k.resize(aii_skel.nonZeros(), numCoefficients);
            aic_coef2k.resize(aic_skel.nonZeros(), numCoefficients);
            acc_coef2k.resize(acc_skel.nonZeros(), numCoefficients);

            aii_coef2k.setFromTriplets(iiTripletList.begin(), iiTripletList.end());
            aic_coef2k.setFromTriplets(icTripletList.begin(), icTripletList.end());
            acc_coef2k.setFromTriplets(ccTripletList.begin(), ccTripletList.end());

            aii_coef2k.makeCompressed();
            aic_coef2k.makeCompressed();
            acc_coef2k.makeCompressed();
        }

        std::unique_ptr<symmetricMatrix> buildAiiMatrix(const std::vector<scalar> &c) {
            std::unique_ptr<symmetricMatrix> ret = std::make_unique<symmetricMatrix>(aii_skel.rows(), aii_skel.cols());
            ret->reserve(aii_skel.nonZeros());

            std::memcpy(ret->outerIndexPtr(), aii_skel.outerIndexPtr(), (aii_skel.outerSize() + 1)*sizeof(typename symmetricMatrix::StorageIndex));
            std::memcpy(ret->innerIndexPtr(), aii_skel.innerIndexPtr(), aii_skel.nonZeros()*sizeof(typename symmetricMatrix::StorageIndex));

            Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1> >(ret->valuePtr(), aii_skel.nonZeros()).noalias() =
                aii_coef2k*Eigen::Map<const Eigen::Matrix<scalar, Eigen::Dynamic, 1> >(c.data(), c.size());
            ret->resizeNonZeros(aii_skel.nonZeros());

            return ret;
        }

        std::unique_ptr<symmetricMatrix> buildAccMatrix(const std::vector<scalar> &c) {
            std::unique_ptr<symmetricMatrix> ret = std::make_unique<symmetricMatrix>(acc_skel.rows(), acc_skel.cols());
            ret->reserve(acc_skel.nonZeros());

            std::memcpy(ret->outerIndexPtr(), acc_skel.outerIndexPtr(), (acc_skel.outerSize() + 1)*sizeof(typename symmetricMatrix::StorageIndex));
            std::memcpy(ret->innerIndexPtr(), acc_skel.innerIndexPtr(), acc_skel.nonZeros()*sizeof(typename symmetricMatrix::StorageIndex));

            Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1> >(ret->valuePtr(), acc_skel.nonZeros()).noalias() =
                acc_coef2k*Eigen::Map<const Eigen::Matrix<scalar, Eigen::Dynamic, 1> >(c.data(), c.size());
            ret->resizeNonZeros(acc_skel.nonZeros());

            return ret;
        }

        std::unique_ptr<matrix> buildAicMatrix(const std::vector<scalar> &c) {
            std::unique_ptr<matrix> ret = std::make_unique<matrix>(aic_skel.rows(), aic_skel.cols());
            ret->reserve(aic_skel.nonZeros());


            std::memcpy(ret->outerIndexPtr(), aic_skel.outerIndexPtr(), (aic_skel.outerSize() + 1)*sizeof(typename symmetricMatrix::StorageIndex));
            std::memcpy(ret->innerIndexPtr(), aic_skel.innerIndexPtr(), aic_skel.nonZeros()*sizeof(typename symmetricMatrix::StorageIndex));

            Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1> >(ret->valuePtr(), aic_skel.nonZeros()).noalias() =
                aic_coef2k*Eigen::Map<const Eigen::Matrix<scalar, Eigen::Dynamic, 1> >(c.data(), c.size());

            return ret;
        }
};

#endif  // LBMATRIXBUILDER_H_
