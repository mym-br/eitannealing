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

	return new Vector(auxv.get(), this->n);
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

MatrixCPJDSManager::MatrixCPJDSManager(Eigen::SparseMatrix<numType> *pdata) {
	int n = pdata->cols();
	int colorCount = 0;
	std::unique_ptr<int[]> reorderIdx(new int[n]); // map for original index -> color-sorted index
								   // color-sort stiffness matrix
	LOG("--Color sorting--");
	std::unique_ptr<int[]> colorsOff = colorSort(pdata, colorCount, reorderIdx, false);
	if (colorsOff == NULL) {
		LOG("Erro no color-sort.");
		return;
	}

	int nPadded = n;
	// padding-fill for warps
	//std::cout << "--Padding-filling--" << std::endl;
	std::vector<int> colorsOffPaddedVec; this->data2 = fillPadding2(pdata, colorCount, colorsOff, colorsOffPaddedVec, nPadded);
	this->n = nPadded;
	this->nOrig = n;

	// corrigir offsets no vetor de cores
	std::shared_ptr<int[]> colorsOffPadded(new int[colorCount + 1]);
	int newSize = 0;
	for (int i = 0; i < colorCount; i++) {
		colorsOffPadded[i] = newSize;
		newSize += (int)ceil((double)(colorsOff[i + 1] - colorsOff[i]) / WARP_SIZE) * WARP_SIZE;
	}
	colorsOffPadded[colorCount] = newSize;
	this->colorCount = colorCount;
	this->colors = colorsOffPadded;

	// map for color-sorted index -> original index
	std::unique_ptr<int[]> unorderIdx(new int[n]);
	// fill padding reorder indices array
	std::unique_ptr<int[]> reorderIdxPadded(new int[nPadded]);
	std::unique_ptr<int[]> unorderIdxPadded(new int[n]);
	fixIndicesMapPadding(n, reorderIdx, unorderIdx, nPadded, reorderIdxPadded, unorderIdxPadded, colorCount, colorsOff);

	this->padded2OriginalIdx = std::move(reorderIdxPadded);
	this->original2PaddedIdx = std::move(unorderIdxPadded);
	//LOGV2(padded2OriginalIdx, nPadded, "Padded 2 Original", LOGCPU);
	//LOGV2(original2PaddedIdx, n, "Original 2 Padded", LOGCPU);


	this->auxv = std::unique_ptr<numType[]>(new numType[this->n]);
	this->auxi = std::unique_ptr<int[]>(new int[this->n]);

	blocks = (int)ceil((double)this->n / BLOCKSIZE);
}

void MatrixCPJDSManager::leftShiftMatrix(std::vector<std::deque<int>> &rowsL, std::vector<std::deque<int>> &rowsU) {
	for (int k = 0; k < data2->outerSize(); ++k)
		for (Eigen::SparseMatrix<numType>::InnerIterator it(*data2, k); it; ++it)
		{
			if (it.row() == it.col()) continue;
			rowsL[it.row()].push_back(it.col()); // adding non zero Aij to each row array
			rowsU[it.col()].push_back(it.row()); // adding non zero Aij to each row array
		}
}

void MatrixCPJDSManager::createDataAndIndicesVectors(numType *mdata, int *indices, int *colOffset, int *colorOffsetCount, std::vector<std::deque<int>> &rowsL, std::vector<std::deque<int>> &rowsU, std::vector<std::deque<int>> &padding) {
	int pos = 0;
	// each block is built in color order
	for (int k = 0; k < colorCount; k++) {
		bool hasElements = true;
		int col = 0;
		// insert diagonal in first column
		colOffset[colorOffsetCount[k] + col] = pos;
		for (int i = colors[k]; i < colors[k + 1]; i++) {
			indices[pos] = i;
			mdata[pos] = data2->coeff(i,i);
			pos++;
		}
		col++;
		while (true) { // iterate over all rows in color block for as long as it has elements
			hasElements = false;
			int startPos = pos;
			for (int i = colors[k]; i < colors[k + 1]; i++) {
				// fill data array in COL-MAJOR format
				if (rowsL[i].size() > 0) {// lower triangular
					int idx = rowsL[i].front(); // first element
					rowsL[i].pop_front(); // remove element

					indices[pos] = idx;
					mdata[pos] = data2->coeff(i, idx);

					pos++;
					hasElements = true;
				}
				else if (rowsU[i].size() > 0) {// upper triangular
					int idx = rowsU[i].front(); // first element
					rowsU[i].pop_front(); // remove element

					indices[pos] = idx;
					mdata[pos] = data2->coeff(idx, i);

					pos++;
					hasElements = true;
				}
				else if (padding[i].size() > 0) {// padding zeroes
					indices[pos] = -1;
					mdata[pos] = 0;
					padding[i].pop_front();

					pos++;
					hasElements = true;
				}
			}
			if (hasElements) {
				colOffset[colorOffsetCount[k] + col] = startPos;
				//std::cout << " colOfsset idx: " << colorOffsetCount[k] + col << std::endl;
				col++;
			}
			else {
				break;
			}
		}
	}
}

void MatrixCPJDSManager::dependencies_analysis2(int n, DependeciesMap * dependenciesMap) {
	// TODO: implement dependencie analysis on the eigen sparse matrix
}