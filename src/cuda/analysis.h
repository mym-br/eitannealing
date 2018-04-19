#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "settings.h"
#include <Eigen/Core>
#include "../basematrix.h"

/** data must represent a square complete(filled) matrix - default reorders by lower triangular */
int * colorSort(int n, numType * data, int * colorCount, int * newIdx);
int * colorSort(Eigen::SparseMatrix<float, Eigen::ColMajor> * data, int * colorCount, int * newIdx);
/** data must represent a square complete(filled) matrix */
int * colorSort(int n, numType * data, int * colorCount, int * newIdx, bool reorderByLowerTriangular);
int * colorSort(matrix * data, int * colorCount, int * newIdx, bool reorderByLowerTriangular);

/** data must represent a square complete(filled) matrix */
void swap(numType * data, int size, int a, int b, int * newIdx);
void swap(matrix * data, int a, int b, int * newIdx);

/** data must represent a square complete(filled) matrix */
numType * fillPadding(int size, numType * data, int colorCount, int * colorOff, int * sizePadding);
numType * fillPadding(matrix * data, int colorCount, int * colorOff, int * sizePadding);

/** data must represent an array */
numType * fillPaddingV(int size, numType * data, int sizePadded, numType * dataPadded, int colorCount, int * colorOff);

/** both reorderIdx and unorderIdx must be changed to take padding in effect */
void fixIndicesMapPadding(int size, int * reorderIdx, int * unorderIdx, int sizePadded, int * reorderIdxFixed, int * unorderIdxFixed, int colorCount, int * colorOff);

#endif /* ANALYSIS_H */