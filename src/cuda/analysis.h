#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "settings.h"
#include <Eigen/Core>
#include "../basematrix.h"
#include <memory>

/** data must represent a square complete(filled) matrix - default reorders by lower triangular */
std::unique_ptr<int[]> colorSort(Eigen::SparseMatrix<double> * data, int &colorCount, std::unique_ptr<int[]> &newIdx);
/** data must represent a square complete(filled) matrix */
std::unique_ptr<int[]> colorSort(Eigen::SparseMatrix<double> * data, int &colorCount, std::unique_ptr<int[]> &newIdx, bool reorderByLowerTriangular);

/** data must represent a square complete(filled) matrix */
void swap(Eigen::SparseMatrix<double> * data, int a, int b, std::unique_ptr<int[]> &newIdx);

/** data must represent a square complete(filled) matrix */
std::unique_ptr<Eigen::SparseMatrix<double>> fillPadding(Eigen::SparseMatrix<double> * data, int colorCount, std::unique_ptr<int[]> &colorOff, std::shared_ptr<int[]> colorsOffPadded, int &sizePadding);

/** both reorderIdx and unorderIdx must be changed to take padding in effect */
void fixIndicesMapPadding(int size, std::unique_ptr<int[]> &reorderIdx, std::unique_ptr<int[]> &unorderIdx, int sizePadded, std::unique_ptr<int[]> &reorderIdxFixed, std::unique_ptr<int[]> &unorderIdxFixed, int colorCount, std::unique_ptr<int[]> &colorOff, std::shared_ptr<int[]> colorsOffPadded);

#endif /* ANALYSIS_H */
