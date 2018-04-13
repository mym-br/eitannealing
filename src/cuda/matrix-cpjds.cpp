#include "matrix-cpjds.h"

#include "settings.h"

#include "utils.h"
#include <vector>
#include <deque>

#include <fstream>
#include <iostream>

#include "color.h"
#include "analysis.h"
#include "model.h"
#include "stiffness.h"

#ifndef USECUDA

struct Dependencies {
	std::vector<int> dependencies;
	int dependenciesSize;
	std::pair <int, int> lower; // upper can be calculated from lower
};

typedef std::vector<std::vector<Dependencies>> DependeciesMap;

void dependencies_analysis(numType * data, int n, DependeciesMap * dependenciesMap) {
	for (int i = 0; i < n; i++) { // for each row...
		for (int j = i + 1; j < n; j++) { // ... for every (upper triangle) column in the row
		//for (int j = 0; j < i; j++) { // ... for every (lower triangle) column in the row
			if (MOD(data[i * n + j]) > EPS) { // ... if element (i, j) is non-zero, search dependencies
				Dependencies d;
				for (int k = 0; k < i; k++) { // check all rows before current row
				//for (int k = 0; k < j; k++) { // check all columns before current column
					// check if (k, i) and (k, j) are both non-zero
					if (MOD(data[k * n + i]) > EPS && MOD(data[k * n + j]) > EPS) {
						d.dependencies.push_back(k); // push dependency (only row in lower triangular is needed)
					}
				}
				d.lower = std::pair<int, int>(i, j);
				d.dependenciesSize = d.dependencies.size();
				(*dependenciesMap)[i].push_back(d);
			}
		}
	}
}

/*
* data: full symmetric matrix (n-by-n), with zeroes - data is row-major!!!
* n: matrix size - must be multiple of WARP_SIZE (32)
* colors: array of each colors offset (size is colorCount + 1, last position being equal to n)
* colorCount: number of colors
* M: resulting matrix in colored pJDS format
*/
int MatrixCPJDSManager::buidMatrixCPJDS(MatrixCPJDS * M) {
	// checking size (must be multiple of WARP_SIZE)
	if (n % WARP_SIZE != 0) {
		LOG("Invalid matrix size for pJDS format!");
		return -1;
	}
	// checking colors offsets (must be multiple of WARP_SIZE)
	for (int i = 0; i <= colorCount; i++) {
		if (colors[i] % WARP_SIZE != 0) {
			LOG("Invalid color offset for colored pJDS format!");
			return -1;
		}
	}

	// left-shift each rows' non-zeroes
	std::vector<std::deque<int>> rowsL(n), rowsU(n), padding(n);
	for (int i = 0; i < n; i++) {
		// saving diagonal Aii's index is not needed (Aii is expected to be positive)
		// lower triangular: add rest of the line with column < row
		for (int j = 0; j < i; j++) {
			if (MOD(data[i * n + j]) > EPS) { // DATA IS ROW-MAJOR!!!
				rowsL[i].push_back(j); // adding Aij to each row array, if non-zero
			}
		}
		// upper triangular: add rest of the line with column > row
		for (int j = i + 1; j < n; j++) {
			if (MOD(data[i * n + j]) > EPS) { // DATA IS ROW-MAJOR!!!
				rowsU[i].push_back(j); // adding Aij to each row array, if non-zero
			}
		}
	}

	// row's length (equivalent to number of non-zeroes)
	int * rowLength = new int[n];
	// row's size, including padding zeroes
	int * rowSize = new int[n];

	int nnz = 0,
		total = 0;
	int blocks = n / WARP_SIZE;
	for (int k = 0; k < blocks; k++) {
		int blockSize = 0;
		for (int i = 0; i < WARP_SIZE; i++) {
			int row = k * WARP_SIZE + i;
			// row length = lower triangular row's nnz + main diagonal + upper triangular row's nnz
			rowLength[row] = rowsL[row].size() + 1 + rowsU[row].size();
			nnz += rowLength[row]; // total number of non-zeroes

			if (rowLength[row] > blockSize) { // lower triangular warp-block size
				blockSize = rowLength[row];
			}
		}
		for (int i = 0; i < WARP_SIZE; i++) {
			int row = k * WARP_SIZE + i;
			rowSize[row] = blockSize;
			// fill zeroes so that all rows in a warp have the same size
			for (int j = rowLength[row]; j < blockSize; j++) {
				padding[row].push_back(-1);
			}
			// add to total element count
			total += blockSize;
		}
	}

	// preparing columns' offsets (COL-MAJOR!!)
	int * colorOffsetCount = new int[colorCount + 1];
	int maxSize = 0,
		offsetSize = 0;
	// find largest row size in each color (should be first row!)
	for (int k = 0; k < colorCount; k++) {
		colorOffsetCount[k] = offsetSize;
		maxSize = rowsL[colors[k]].size() + 1 + rowsU[colors[k]].size();
		// step is warp size because all rows in a warp block have same size
		for (int i = colors[k]; i < colors[k + 1]; i += WARP_SIZE) {
			int size = rowsL[i].size() + 1 + rowsU[i].size();
			if (size > maxSize) {
				std::cout << "something is not right... " << size << " (" << i << ") > " << maxSize << " (" << k << ")" << std::endl;
				maxSize = size;
			}
		}
		// lower triangular
		offsetSize += maxSize;
	}
	colorOffsetCount[colorCount] = offsetSize;

	// offset arrays (each color has a different size, computed above)
	int * colOffset = new int[offsetSize];

	// data and indices arrays
	numType * mdata = new numType[total];
	int * indices = new int[total];

	// values
	//std::cout << "values:" << std::endl;
	//for (int i = 0; i < n; i++) {
	//	std::cout << data[i * n + i] << "\t";
	//	for (int j = 0; j < rowsL[i].size(); j++) {
	//		std::cout << data[i * n + rowsL[i][j]] << "\t";
	//	}
	//	for (int j = 0; j < rowsU[i].size(); j++) {
	//		std::cout << data[i * n + rowsU[i][j]] << "\t";
	//	}
	//	for (int j = 0; j < padding[i].size(); j++) {
	//		std::cout << "0\t";
	//	}
	//	std::cout << std::endl;
	//}
	//// indices
	//std::cout << "\nindices:" << std::endl;
	//for (int i = 0; i < n; i++) {
	//	std::cout << i << "\t";
	//	for (int j = 0; j < rowsL[i].size(); j++) {
	//		std::cout << rowsL[i][j] << "\t";
	//	}
	//	for (int j = 0; j < rowsU[i].size(); j++) {
	//		std::cout << rowsU[i][j] << "\t";
	//	}
	//	for (int j = 0; j < padding[i].size(); j++) {
	//		std::cout << "x\t";
	//	}
	//	std::cout << std::endl;
	//}

	int pos = 0;
	// each block is built in color order
	for (int k = 0; k < colorCount; k++) {
		bool hasElements = true;
		int col = 0;
		// insert diagonal in first column
		colOffset[colorOffsetCount[k] + col] = pos;
		for (int i = colors[k]; i < colors[k + 1]; i++) {
			indices[pos] = i;
			mdata[pos] = data[i * n + i];
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
					mdata[pos] = data[i * n + idx];

					pos++;
					hasElements = true;
				}
				else if (rowsU[i].size() > 0) {// upper triangular
					int idx = rowsU[i].front(); // first element
					rowsU[i].pop_front(); // remove element

					indices[pos] = idx;
					mdata[pos] = data[i * n + idx];

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

	/* computing (x,y)=>IDX map */
	std::vector<std::map<int, int>> rowCol2IdxMap(n);
	for (int k = 0; k < colorCount; k++) { // color block
		for (int i = colors[k]; i < colors[k + 1]; i++) { // iterating rows in such color
			for (int j = 0; j < rowLength[i]; j++) {
				int idx = colOffset[colorOffsetCount[k] + j] + i - colors[k];
				rowCol2IdxMap[i].insert(std::pair<int, int>(indices[idx], idx));
			}
		}
	}
	this->coord2IndexMap = rowCol2IdxMap;

	// ---- DEBUG ----
	//for (int i = 0; i < n; i++) { // iterating rows in such color
	//	for (std::map<int, int>::iterator it = rowCol2IdxMap[i].begin(); it != rowCol2IdxMap[i].end(); ++it) {
	//		std::cout << it->first << ">" << it->second << '\t';
	//	}
	//	std::cout << std::endl;
	//}
	// ---- DEBUG END ----

	/* computing dependencies */
	DependeciesMap dependencies(n);
	dependencies_analysis(data, n, &dependencies);

	// fix dependencies' indices
	//std::cout << "\nFixing dependencies..." << std::endl;
	DependeciesMap dependenciesCPJDS(n);
	for (int i = 0; i < n; i++) {
		std::vector<Dependencies> dv = dependencies[i];

		for (int j = 0; j < dv.size(); j++) {
			Dependencies d = dv[j];
			// this code was removed, as all non-zeroes must be listed, regardless of having dependencies or not
			//if (d.dependencies.empty()) {
			//	std::cout << "skipping dependency..." << std::endl;
			//	continue;
			//}
			Dependencies dCPJDS;
			int coli = d.lower.first;
			int colj = d.lower.second;
			//dCPJDS.lower.first = coli % n;
			//dCPJDS.lower.second = colj % n;
			for (int k = 0; k < d.dependencies.size(); k++) {
				int row = d.dependencies[k];
				// convert to data array index
				int idxi = coordinates2Index(row, coli);
				int idxj = coordinates2Index(row, colj);
				dCPJDS.dependencies.push_back(idxi);
				dCPJDS.dependencies.push_back(idxj);
			}
			dCPJDS.lower = std::pair<int, int>(coli, colj);
			dCPJDS.dependenciesSize = d.dependenciesSize;
			dependenciesCPJDS[i].push_back(dCPJDS);
		}
	}

	//// ---- DEBUG ----
	//std::cout << "Dependencies" << std::endl;
	//for (int i = 0; i < n; i++) { // iterating rows
	//	std::vector<Dependencies> dvCPJDS = dependenciesCPJDS[i];
	//	for (int j = 0; j < dvCPJDS.size(); j++) {
	//		Dependencies d = dvCPJDS[j];
	//		std::cout << "row " << i << ": (" << d.lower.first << ", " << d.lower.second << ") - " << d.dependenciesSize << std::endl;
	//		for (int k = 0; k < d.dependencies.size(); k += 2) {
	//			std::cout << mdata[d.dependencies[k]] << ", " << mdata[d.dependencies[k + 1]] << std::endl;
	//		}
	//		std::cout << std::endl;
	//	}
	//}
	//std::cout << std::endl;
	//// ---- DEBUG END ----

	// make all dependecies' arrays for the same row have the same size
	// (it is supposed that dependencies's array size is smaller than warp size - however, warp will not be filled)
	//std::vector<std::vector<std::deque<int>>> rowsDeps(n);
	for (int i = 0; i < n; i++) { // iterating rows
		std::vector<Dependencies> * dvCPJDS = &(dependenciesCPJDS[i]);

		int maxSize = 0;
		for (int j = 0; j < (*dvCPJDS).size(); j++) { // for each non-zero element with dependency..
			int thisDepSize = (*dvCPJDS)[j].dependencies.size();
			if (thisDepSize > maxSize) { // ...check size, to find largest...
				maxSize = thisDepSize;
			}
		}
		for (int j = 0; j < (*dvCPJDS).size(); j++) { // ... and pad the others
			int thisDepSize = (*dvCPJDS)[j].dependencies.size();
			for (int k = thisDepSize; k < maxSize; k++) {
				(*dvCPJDS)[j].dependencies.push_back(-1);
			}
		}
	}

	//// ---- DEBUG ----
	//std::cout << "Dependencies (padded)" << std::endl;
	//for (int i = 0; i < n; i++) { // iterating rows
	//	std::vector<Dependencies> dvCPJDS = dependenciesCPJDS[i];
	//	for (int j = 0; j < dvCPJDS.size(); j++) {
	//		Dependencies d = dvCPJDS[j];
	//		std::cout << "row " << i << ": (" << d.lower.first << ", " << d.lower.second << ") - " << d.dependenciesSize << std::endl;
	//		for (int k = 0; k < d.dependencies.size(); k += 2) {
	//			if (mdata[d.dependencies[k]] > -1) {
	//				std::cout << mdata[d.dependencies[k]];
	//			}
	//			else {
	//				std::cout << "x";
	//			}
	//			std::cout << ", ";
	//			if (mdata[d.dependencies[k + 1]] > -1) {
	//				std::cout << mdata[d.dependencies[k + 1]];
	//			}
	//			else {
	//				std::cout << "x";
	//			}
	//			std::cout << std::endl;
	//		}
	//		std::cout << std::endl;
	//	}
	//}
	//std::cout << std::endl;
	//// ---- DEBUG END ----

	// create dependencies's arrays in CPJDS format
	std::vector<int> depsArrDiag;
	std::vector<int> depsArrRow;
	std::vector<int> depsIdxArr(n);
	std::vector<int> depsLowerArr;
	std::vector<int> depsUpperArr;
	std::vector<int> depsNNZArr(n + 1);
	std::vector<int> depsSizeArr(n);
	for (int i = 0; i < n; i++) { // iterating rows
		std::vector<Dependencies> dvCPJDS = dependenciesCPJDS[i];

		// initial index
		depsIdxArr[i] = depsArrRow.size(); // size of depsArrRow = size of depsArrDiag

		int depsSize = 0;
		if (dvCPJDS.size() > 0) {
			depsSize = dvCPJDS[0].dependencies.size() / 2; // all dvCPJDS have same size (padded in previous step)
			// also, dependenciesCPJDS.dependencies array is doubled as it contains both row and diagonal elements
		}
		// size
		depsSizeArr[i] = depsSize;
		
		// adding elements to "global" array, column-major
		for (int j = 0; j < depsSize; j++) {
			for (int k = 0; k < dvCPJDS.size(); k++) {
				depsArrDiag.push_back(dvCPJDS[k].dependencies[2 * j]);
				depsArrRow.push_back(dvCPJDS[k].dependencies[2 * j + 1]);
			}
		}

		depsNNZArr[i] = depsLowerArr.size(); // size of depsLowerArr = size of depsUpperArr
		// element's indices (lower and upper)
		for (int j = 0; j < dvCPJDS.size(); j++) {
			depsLowerArr.push_back(coord2IndexMap[dvCPJDS[j].lower.first][dvCPJDS[j].lower.second]);
			depsUpperArr.push_back(coord2IndexMap[dvCPJDS[j].lower.second][dvCPJDS[j].lower.first]);
		}
	}
	depsNNZArr[n] = depsLowerArr.size(); // size of depsLowerArr = size of depsUpperArr

	//// ---- DEBUG ----
	//std::cout << "Dependencies (padded, single array)" << std::endl;
	//for (int i = 0; i < n; i++) { // iterating rows
	//	int depsSize = depsSizeArr[i];
	//	if (depsSize > 0) {
	//		std::cout << "row " << i << std::endl;
	//		// adding elements to "global" array, column-major
	//		for (int j = 0; j < depsSize; j++) {
	//			if (depsArrDiag[depsIdxArr[i] + j] > -1) {
	//				std::cout << mdata[depsArrDiag[depsIdxArr[i] + j]];
	//			}
	//			else {
	//				std::cout << "x";
	//			}
	//			std::cout << "\t";
	//		}
	//		std::cout << std::endl;
	//
	//		for (int j = 0; j < depsSize; j++) {
	//			if (depsArrRow[depsIdxArr[i] + j] > -1) {
	//				std::cout << mdata[depsArrRow[depsIdxArr[i] + j]];
	//			}
	//			else {
	//				std::cout << "x";
	//			}
	//			std::cout << "\t";
	//		}
	//		std::cout << std::endl;
	//
	//		// element's indices (lower and upper)
	//		for (int j = 0; j < dependenciesCPJDS[i].size(); j++) {
	//			depsLowerArr[depsNNZArr[i] + j];
	//			depsUpperArr[depsNNZArr[i] + j];
	//		}
	//	}
	//}
	//std::cout << std::endl;
	//std::cout << "Dependencies (dependencies)" << std::endl;
	//for (int i = 0; i < n; i++) { // iterating rows
	//	if (dependenciesCPJDS[i].size() > 0) {
	//		std::cout << "row " << i << std::endl;
	//		for (int j = 0; j < dependenciesCPJDS[i].size(); j++) {
	//			std::cout << depsLowerArr[depsNNZArr[i] + j] << "(" << mdata[depsLowerArr[depsNNZArr[i] + j]] << ")" << "\t";
	//		}
	//		std::cout << std::endl;
	//		for (int j = 0; j < dependenciesCPJDS[i].size(); j++) {
	//			std::cout << depsUpperArr[depsNNZArr[i] + j] << "(" << mdata[depsUpperArr[depsNNZArr[i] + j]] << ")" << "\t";
	//		}
	//		std::cout << std::endl;
	//	}
	//}
	//std::cout << std::endl;
	//// ---- DEBUG END ----

	// create arrays in specified formats (MatrixDependencies e ElementDependencies)
	// element and diagonal rows dependency's index in data array
	int * dependencyRowDataIndex = new int[depsArrRow.size()];
	int * dependencyDiagDataIndex = new int[depsArrDiag.size()];
	// dependency's array initial index (index in dependencyRowDataIndex and dependencyDiagDataIndex)
	int * dependencyArrayInitialIndex = new int[depsIdxArr.size()];
	// lower and upper triangulars non-zero elements (with dependencies) index in data array
	int * dependencyLowerDataIndex = new int[depsLowerArr.size()];
	int * dependencyUpperDataIndex = new int[depsUpperArr.size()];
	// data array index for lower and upper triangular's elements (index in dependencyLowerDataIndex and dependencyUpperDataIndex)
	int * nnzElementDataIndex = new int[depsNNZArr.size()];
	// dependency's count
	int * dependenciesSize = new int[depsSizeArr.size()];

	// copy arrays
	for (int i = 0; i < depsArrRow.size(); i++) dependencyRowDataIndex[i] = depsArrRow[i];
	for (int i = 0; i < depsArrDiag.size(); i++) dependencyDiagDataIndex[i] = depsArrDiag[i];
	for (int i = 0; i < depsIdxArr.size(); i++) dependencyArrayInitialIndex[i] = depsIdxArr[i];
	for (int i = 0; i < depsLowerArr.size(); i++) dependencyLowerDataIndex[i] = depsLowerArr[i];
	for (int i = 0; i < depsUpperArr.size(); i++) dependencyUpperDataIndex[i] = depsUpperArr[i];
	for (int i = 0; i < depsNNZArr.size(); i++) nnzElementDataIndex[i] = depsNNZArr[i];
	for (int i = 0; i < depsSizeArr.size(); i++) dependenciesSize[i] = depsSizeArr[i];

	int depsSize = depsArrRow.size(), // = depsArrDiag.size()
		idxSize = depsLowerArr.size(); // = depsUpperArr.size()

	//// sanity check
	//std::cout << "dependency lower and upper indices sanity check...." << std::endl;
	//for (int i = 0; i < depsLowerArr.size(); i++) {
	//	for (int j = i + 1; j < depsLowerArr.size(); j++) {
	//		if (depsLowerArr[i] == depsLowerArr[j]) {
	//			std::cout << "repeated lower indices..." << std::endl;
	//		}
	//	}
	//}
	//for (int i = 0; i < depsUpperArr.size(); i++) {
	//	for (int j = i + 1; j < depsUpperArr.size(); j++) {
	//		if (depsUpperArr[i] == depsUpperArr[j]) {
	//			std::cout << "repeated upper indices..." << std::endl;
	//		}
	//	}
	//}
	//for (int i = 0; i < depsNNZArr.size(); i++) {
	//	for (int j = i + 1; j < depsNNZArr.size(); j++) {
	//		if (depsNNZArr[i] == depsNNZArr[j]) {
	//			std::cout << "repeated depsNNZArr indices..." << std::endl;
	//		}
	//	}
	//}
	//std::cout << "dependency lower and upper indices sanity check....done!" << std::endl;

	/* set matrix data */
	(*M).matrixData.data = mdata;
	(*M).matrixData.indices = indices;
	(*M).matrixData.rowLength = rowLength;
	(*M).matrixData.rowSize = rowSize;
	(*M).matrixData.colOffset = colOffset;
	(*M).matrixData.n = n;
	(*M).matrixData.nnz = nnz;
	(*M).matrixData.elCount = total;
	(*M).matrixData.offsetCount = offsetSize;
	/* set matrix colors */
	(*M).matrixColors.colors = colors;
	(*M).matrixColors.colorCount = colorCount;
	(*M).matrixColors.colorsColOffsetSize = colorOffsetCount;
	/* set matrix dependencies*/
	(*M).matrixDependencies.dependencyRowDataIndex = dependencyRowDataIndex;
	(*M).matrixDependencies.dependencyDiagDataIndex = dependencyDiagDataIndex;
	(*M).matrixDependencies.dependencyLowerDataIndex = dependencyLowerDataIndex;
	(*M).matrixDependencies.dependencyUpperDataIndex = dependencyUpperDataIndex;
	(*M).matrixDependencies.dependencyArrayInitialIndex = dependencyArrayInitialIndex;
	(*M).matrixDependencies.dependenciesSize = dependenciesSize;
	(*M).matrixDependencies.nnzElementDataIndex = nnzElementDataIndex;
	(*M).matrixDependencies.depsSize = depsSize;
	(*M).matrixDependencies.idxSize = idxSize;
	//(*M).dependenciesSupportData = dependenciesSupportData; // -- removed

	(*M).preconditionedData = new numType[total];
	for (int i = 0; i < total; i++) {
		(*M).preconditionedData[i] = 0;
	}

	return 0;
}

void MatrixCPJDSManager::deleteMatrixCPJDS(MatrixCPJDS M) {
	delete M.matrixData.data;
	delete M.matrixData.indices;
	delete M.matrixData.rowLength;
	delete M.matrixData.rowSize;
	delete M.matrixData.colOffset;
	//delete M.matrixColors.colors; // deleted in destructor
	delete M.matrixColors.colorsColOffsetSize;

	delete M.matrixDependencies.dependencyRowDataIndex;
	delete M.matrixDependencies.dependencyDiagDataIndex;
	delete M.matrixDependencies.dependencyLowerDataIndex;
	delete M.matrixDependencies.dependencyUpperDataIndex;
	delete M.matrixDependencies.dependencyArrayInitialIndex;
	delete M.matrixDependencies.dependenciesSize;
	delete M.matrixDependencies.nnzElementDataIndex;
	
	delete M.preconditionedData;
};

/* sets an element value according to row and column indexes */
void MatrixCPJDSManager::set(MatrixCPJDS M, int row, int col, numType val) {
	int rowP = original2PaddedIdx[row];
	int colP = original2PaddedIdx[col];
	int idx = coordinates2Index(rowP, colP);
	M.matrixData.data[idx] = val;
}
/* increments an element value according to row and column indexes */
void MatrixCPJDSManager::increment(MatrixCPJDS M, int row, int col, numType val) {
	int rowP = original2PaddedIdx[row];
	int colP = original2PaddedIdx[col];
	int idx = coordinates2Index(rowP, colP);
	M.matrixData.data[idx] += val;
}
/* method for restoring CPJDS-transformed vector to its original size and indexes */
std::vector<numType> MatrixCPJDSManager::restore(Vector * v, int groundNode){
	// groundNode = -1: no ground node
	std::vector<numType> restored(nOrig);

	for (int i = 0; i < nOrig; i++) {
		restored[i] = v->getData()[original2PaddedIdx[i]];
	}

	return restored;
}

#endif

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
		newSize += ceil((double)(colorsOff[i + 1] - colorsOff[i]) / WARP_SIZE) * WARP_SIZE;
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

	blocks = ceil((double)n / BLOCKSIZE);
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