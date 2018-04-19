#include "matrix-cpjds.h"

#include "settings.h"

#include "utils.h"
#include <vector>
#include <deque>

#include <fstream>
#include <iostream>

#include "color.h"
#include "analysis.h"

#include <Eigen/Core>
#include "../basematrix.h"

/*
 * This constructor does require to have been pre-processed. Thus, it is still
 * not color-sorted nor does it have padding. Also, it is not filled with padding
 * and n is not necessarily a multiple of WARP_SIZE.
 */
MatrixCPJDSManager::MatrixCPJDSManager(numType * pdata, int n) {
	int colorCount = 0;
	int * reorderIdx = new int[n]; // map for original index -> color-sorted index
	// color-sort stiffness matrix
	LOG("--Color sorting--");
	int * colorsOff = colorSort(n, pdata, &colorCount, reorderIdx, false);
	if (colorsOff == NULL) {
		LOG("Erro no color-sort.");
		delete reorderIdx;
		return;
	}

	int nPadded = n;
	// padding-fill for warps
	//std::cout << "--Padding-filling--" << std::endl;
	this->data = fillPadding(n, pdata, colorCount, colorsOff, &nPadded);
	this->n = nPadded;
	this->nOrig = n;

	// corrigir offsets no vetor de cores
	int * colorsOffPadded = new int[colorCount + 1];
	int newSize = 0;
	for (int i = 0; i < colorCount; i++) {
		colorsOffPadded[i] = newSize;
		newSize += (int)ceil((double)(colorsOff[i + 1] - colorsOff[i]) / WARP_SIZE) * WARP_SIZE;
	}
	colorsOffPadded[colorCount] = newSize;

	this->colorCount = colorCount;
	this->colors = colorsOffPadded;

	// map for color-sorted index -> original index
	int * unorderIdx = new int[n];
	// fill padding reorder indices array
	int * reorderIdxPadded = new int[nPadded];
	int * unorderIdxPadded = new int[n];
	fixIndicesMapPadding(n, reorderIdx, unorderIdx, nPadded, reorderIdxPadded, unorderIdxPadded, colorCount, colorsOff);

	delete reorderIdx;
	delete unorderIdx;

	this->padded2OriginalIdx = reorderIdxPadded;
	this->original2PaddedIdx = unorderIdxPadded;
	//LOGV2(padded2OriginalIdx, nPadded, "Padded 2 Original", LOGCPU);
	//LOGV2(original2PaddedIdx, n, "Original 2 Padded", LOGCPU);


	this->auxv = new numType[this->n];
	this->auxi = new int[this->n];

	blocks = (int)ceil((double)n / BLOCKSIZE);
}

/*
* data: full symmetric matrix (n-by-n), with zeroes - data is row-major!!!
* n: matrix size - must be multiple of WARP_SIZE (32)
* colors: array of each colors offset (size is colorCount + 1, last position being equal to n)
* colorCount: number of colors
*/
MatrixCPJDSManager::MatrixCPJDSManager(numType * data, int n, int * colors, int colorCount, int nOrig) {
	this->data = data;
	this->n = n;
	this->colors = colors;
	this->colorCount = colorCount;
	this->nOrig = nOrig;
};

MatrixCPJDSManager::~MatrixCPJDSManager() {
	delete data;
	delete colors;
	delete padded2OriginalIdx;
	delete original2PaddedIdx;
	delete auxv;
	delete auxi;
}

/* this method allows the conversion os (row, col) coordinates to the index in the data and indices arrays */
int MatrixCPJDSManager::coordinates2Index(int row, int col) {
	std::map<int, int> columnsMap = coord2IndexMap[row];
	if (columnsMap.find(col) == columnsMap.end()) {
		return -1;
	}
	else {
		return columnsMap[col];
	}
}

/* method for rearranging a regular vector to follow the CPJDS transformations */
Vector * MatrixCPJDSManager::transform(std::vector<numType> v, bool removeGround) {
	//FIXME: transform precisa usar mapeamento!

	//Vector * transformed = new Vector(this->n);
	numType aux = 0;
	if (removeGround) {
		aux = GROUNDNODE < 0 || GROUNDNODE >= v.size() ? 0 : v[GROUNDNODE];
	}
	for (int i = 0; i < this->n; i++) {
		auxv[i] = 0;
	}

	// removing ground node and changing to numtype array
	numType * vecOrig = new numType[v.size() - 1];
	for (int i = 0; i < this->n; i++) {
		vecOrig[i] = 0;
	}
	for (int i = 0; i < v.size(); i++) {
		int pos = i;
		if (pos == GROUNDNODE) { // skipping ground node
			continue;
		}
		if (pos > GROUNDNODE) { // as ground node was skipped, fix row coordinates
			pos--;
		}

		vecOrig[pos] = v[i] - aux;
	}

	// mapping to padded format
	for (int i = 0; i < v.size() - 1; i++) {
		auxv[original2PaddedIdx[i]] = vecOrig[i];
	}
	delete vecOrig;

	return new Vector(auxv, this->n);
}

/* method for creating an "electrode" mask vector */
Vector * MatrixCPJDSManager::mask() {
	numType * mask_h = new numType[n];
	for (int i = 0; i < n; i++) {
		mask_h[i] = 0;
	}
	// electrodes are last 32 (31, if one of them is used as ground) elements in the matrix
	// electrodes are first 32 only on the solution coefficientes vector
	int size = nOrig; // must consider the size of the original matrix, where electrodes are the last elements
	int firstElectrode = size - 31; // FIXME: considering 31 electrodes for now, as the last one is used as ground
	for (int i = firstElectrode; i < nOrig; i++) {
		mask_h[original2PaddedIdx[i]] = 1;
	}
	LOGV2(mask_h, n, "mask vector", LOGCPU);
	Vector * mask = new Vector(mask_h, n);
	delete mask_h;
	return mask;
}

MatrixCPJDSManager::MatrixCPJDSManager(Eigen::SparseMatrix<double> *pdata) {
	int n = pdata->cols();
	int colorCount = 0;
	int * reorderIdx = new int[n]; // map for original index -> color-sorted index
								   // color-sort stiffness matrix
	LOG("--Color sorting--");
	int * colorsOff = colorSort(pdata, &colorCount, reorderIdx, false);
	if (colorsOff == NULL) {
		LOG("Erro no color-sort.");
		delete reorderIdx;
		return;
	}

	int nPadded = n;
	// padding-fill for warps
	//std::cout << "--Padding-filling--" << std::endl;
	this->data = fillPadding(pdata, colorCount, colorsOff, &nPadded);
	this->n = nPadded;
	this->nOrig = n;

	// corrigir offsets no vetor de cores
	int * colorsOffPadded = new int[colorCount + 1];
	int newSize = 0;
	for (int i = 0; i < colorCount; i++) {
		colorsOffPadded[i] = newSize;
		newSize += (int)ceil((double)(colorsOff[i + 1] - colorsOff[i]) / WARP_SIZE) * WARP_SIZE;
	}
	colorsOffPadded[colorCount] = newSize;

	this->colorCount = colorCount;
	this->colors = colorsOffPadded;

	// map for color-sorted index -> original index
	int * unorderIdx = new int[n];
	// fill padding reorder indices array
	int * reorderIdxPadded = new int[nPadded];
	int * unorderIdxPadded = new int[n];
	fixIndicesMapPadding(n, reorderIdx, unorderIdx, nPadded, reorderIdxPadded, unorderIdxPadded, colorCount, colorsOff);

	delete reorderIdx;
	delete unorderIdx;

	this->padded2OriginalIdx = reorderIdxPadded;
	this->original2PaddedIdx = unorderIdxPadded;
	//LOGV2(padded2OriginalIdx, nPadded, "Padded 2 Original", LOGCPU);
	//LOGV2(original2PaddedIdx, n, "Original 2 Padded", LOGCPU);


	this->auxv = new numType[this->n];
	this->auxi = new int[this->n];

	blocks = (int)ceil((double)n / BLOCKSIZE);
}