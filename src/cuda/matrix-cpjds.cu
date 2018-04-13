#include "matrix-cpjds.h"

#include "settings.h"
#ifdef USECUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "utils.h"
#include <vector>
#include <deque>

#include <fstream>
#include <iostream>

#ifdef DEBUG
#include <sstream>
#endif

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

/* CPJDS matrix element setter */
__global__ void cm_set_cpjds(int idx, numType * data, numType val);

/* CPJDS matrix element incrementer */
__global__ void cm_increment_cpjds(int idx, numType * data, numType val);

/* CPJDS matrix batch incrementer */
__global__ void cm_increment_cpjds(int size, numType * data, numType * vals, int * indices);

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * bData);

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds2(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * bData);

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_inner_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * yData, numType * partial);

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult(numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorColOffset, numType * xData, numType * yData, numType * partial);

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_2(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * yData, numType * partial);

/* CPJDS lower triangular solver */
__global__ void cmv_solve_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * bData);

/* CPJDS upper triangular solver */
__global__ void cmv_solve_t_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * bData);

/* CPJDS lower and upper triangular solver */
__global__ void cmv_solve(numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * bData);

/* CPJDS lower and upper triangular solver */
__global__ void cmv_solve_inner(numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * bData, numType * interm, numType * partial);

/*
* data: full symmetric matrix (n-by-n), with zeroes - data is row-major!!!
* n: matrix size - must be multiple of WARP_SIZE (32)
* colors: array of each colors offset (size is colorCount + 1, last position being equal to n)
* colorCount: number of colors
* M: resulting matrix in colored pJDS format
*/
int MatrixCPJDSManager::buidMatrixCPJDS(MatrixCPJDS * M, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients) {
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

	std::ofstream myfile;
	myfile.open("E:\\Temp\\adata.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			numType val = data[i * n + j];
			if (MOD(val) < EPS) {
				val = 0;
			}
			myfile << val << "\t";
		}
		myfile << "\n";
	}
	myfile.flush();
	myfile.close();

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
	// preconditioner data
	int * pdata = new int[total];

	// copy arrays
	for (int i = 0; i < depsArrRow.size(); i++) dependencyRowDataIndex[i] = depsArrRow[i];
	for (int i = 0; i < depsArrDiag.size(); i++) dependencyDiagDataIndex[i] = depsArrDiag[i];
	for (int i = 0; i < depsIdxArr.size(); i++) dependencyArrayInitialIndex[i] = depsIdxArr[i];
	for (int i = 0; i < depsLowerArr.size(); i++) dependencyLowerDataIndex[i] = depsLowerArr[i];
	for (int i = 0; i < depsUpperArr.size(); i++) dependencyUpperDataIndex[i] = depsUpperArr[i];
	for (int i = 0; i < depsNNZArr.size(); i++) nnzElementDataIndex[i] = depsNNZArr[i];
	for (int i = 0; i < depsSizeArr.size(); i++) dependenciesSize[i] = depsSizeArr[i];
	for (int i = 0; i < total; i++) pdata[i] = 0;

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

	/* matrix data */
	numType * c_mdata;
	int * c_indices;
	int * c_rowLength;
	int * c_rowSize;
	int * c_colOffset;
	/* matrix colors */
	int * c_colors;
	int * c_colorsColOffsetSize;
	/* matrix dependencies*/
	int * c_dependencyRowDataIndex;
	int * c_dependencyDiagDataIndex;
	int * c_dependencyLowerDataIndex;
	int * c_dependencyUpperDataIndex;
	int * c_dependencyArrayInitialIndex;
	int * c_dependenciesSize;
	int * c_nnzElementDataIndex;
	/* preconditioner data */
	numType * c_pdata;

	int depsSize = depsArrRow.size(), // = depsArrDiag.size()
		idxSize = depsLowerArr.size(); // = depsUpperArr.size()

	// CUDA device memory allocation
	/* matrix data */
	cudaMalloc((void**)& c_mdata,						total * sizeof (numType));
	cudaMalloc((void**)& c_indices,						total * sizeof (int));
	cudaMalloc((void**)& c_rowLength,					n * sizeof (int));
	cudaMalloc((void**)& c_rowSize,						n * sizeof (int));
	cudaMalloc((void**)& c_colOffset,					offsetSize * sizeof (int));
	/* matrix colors */
	cudaMalloc((void**)& c_colors,						(colorCount + 1) * sizeof (int));
	cudaMalloc((void**)& c_colorsColOffsetSize,			colorCount * sizeof(int));
	/* matrix dependencies*/
	cudaMalloc((void**)& c_dependencyRowDataIndex,		depsSize * sizeof (int));
	cudaMalloc((void**)& c_dependencyDiagDataIndex,		depsSize * sizeof (int));
	cudaMalloc((void**)& c_dependencyLowerDataIndex,	idxSize * sizeof (int));
	cudaMalloc((void**)& c_dependencyUpperDataIndex,	idxSize * sizeof (int));
	cudaMalloc((void**)& c_dependencyArrayInitialIndex,	n * sizeof (int));
	cudaMalloc((void**)& c_dependenciesSize,			n * sizeof (int));
	cudaMalloc((void**)& c_nnzElementDataIndex,			(n + 1) * sizeof (int));
	/* preconditioner data */
	cudaMalloc((void**)& c_pdata,						total * sizeof (numType));

	// CUDA device memory transfer
	/* matrix data */
	cudaMemcpy(c_mdata,							mdata,							(size_t)total * sizeof (numType),	cudaMemcpyHostToDevice);
	cudaMemcpy(c_indices,						indices,						(size_t)total * sizeof (int),		cudaMemcpyHostToDevice);
	cudaMemcpy(c_rowLength,						rowLength,						(size_t)n * sizeof (int),			cudaMemcpyHostToDevice);
	cudaMemcpy(c_rowSize,						rowSize,						(size_t)n * sizeof (int),			cudaMemcpyHostToDevice);
	cudaMemcpy(c_colOffset,						colOffset,						(size_t)offsetSize * sizeof (int),	cudaMemcpyHostToDevice);
	/* matrix colors */
	cudaMemcpy(c_colors,						colors,							(size_t)(colorCount + 1) * sizeof (int),cudaMemcpyHostToDevice);
	cudaMemcpy(c_colorsColOffsetSize,			colorOffsetCount,				(size_t)colorCount * sizeof (int),	cudaMemcpyHostToDevice);
	/* matrix dependencies*/
	cudaMemcpy(c_dependencyRowDataIndex,		dependencyRowDataIndex,			(size_t)depsSize * sizeof (int),	cudaMemcpyHostToDevice);
	cudaMemcpy(c_dependencyDiagDataIndex,		dependencyDiagDataIndex,		(size_t)depsSize * sizeof (int),	cudaMemcpyHostToDevice);
	cudaMemcpy(c_dependencyLowerDataIndex,		dependencyLowerDataIndex,		(size_t)idxSize * sizeof (int),		cudaMemcpyHostToDevice);
	cudaMemcpy(c_dependencyUpperDataIndex,		dependencyUpperDataIndex,		(size_t)idxSize * sizeof (int),		cudaMemcpyHostToDevice);
	cudaMemcpy(c_dependencyArrayInitialIndex,	dependencyArrayInitialIndex,	(size_t)n * sizeof (int),			cudaMemcpyHostToDevice);
	cudaMemcpy(c_dependenciesSize,				dependenciesSize,				(size_t)n * sizeof (int),			cudaMemcpyHostToDevice);
	cudaMemcpy(c_nnzElementDataIndex,			nnzElementDataIndex,			(size_t)(n + 1) * sizeof (int),		cudaMemcpyHostToDevice);
	/* preconditioner data */
	cudaMemcpy(c_pdata,							pdata,							(size_t)total * sizeof (numType),	cudaMemcpyHostToDevice);

	/* set matrix data */
	(*M).matrixData.data = c_mdata;
	(*M).matrixData.indices = c_indices;
	(*M).matrixData.rowLength = c_rowLength;
	(*M).matrixData.rowSize = c_rowSize;
	(*M).matrixData.colOffset = c_colOffset;
	(*M).matrixData.n = n;
	(*M).matrixData.nnz = nnz;
	(*M).matrixData.elCount = total;
	(*M).matrixData.offsetCount = offsetSize;
	/* set matrix colors (CPU MEMORY!) */
	(*M).matrixColors.colors = colors;
	(*M).matrixColors.colorCount = colorCount;
	(*M).matrixColors.colorsColOffsetSize = colorOffsetCount;
	(*M).matrixColors.colors_d = c_colors;
	(*M).matrixColors.colorsColOffsetSize_d = c_colorsColOffsetSize;
	/* set matrix dependencies*/
	(*M).matrixDependencies.dependencyRowDataIndex = c_dependencyRowDataIndex;
	(*M).matrixDependencies.dependencyDiagDataIndex = c_dependencyDiagDataIndex;
	(*M).matrixDependencies.dependencyLowerDataIndex = c_dependencyLowerDataIndex;
	(*M).matrixDependencies.dependencyUpperDataIndex = c_dependencyUpperDataIndex;
	(*M).matrixDependencies.dependencyArrayInitialIndex = c_dependencyArrayInitialIndex;
	(*M).matrixDependencies.dependenciesSize = c_dependenciesSize;
	(*M).matrixDependencies.nnzElementDataIndex = c_nnzElementDataIndex;
	(*M).matrixDependencies.depsSize = depsSize;
	(*M).matrixDependencies.idxSize = idxSize;
	//(*M).dependenciesSupportData = dependenciesSupportData; // -- removed
	/* set preconditioner data */
	(*M).preconditionedData = c_pdata;
	// OBS: is color-data needed in GPU memory?

	// setting CPU aux data
	(*M).cpuData.data = mdata;
	(*M).cpuData.indices = indices;
	(*M).cpuData.precond = new numType[total]; for (int i = 0; i < total; i++) (*M).cpuData.precond[i] = 0;
	
	// free CPU memory used for initializing GPU data
	/* matrix data */
	delete rowLength;
	delete rowSize;
	delete colOffset;
	/* matrix dependencies*/
	delete dependencyRowDataIndex;
	delete dependencyDiagDataIndex;
	delete dependencyArrayInitialIndex;
	delete dependencyLowerDataIndex;
	delete dependencyUpperDataIndex;
	delete nnzElementDataIndex;
	delete dependenciesSize;
	/* preconditioner data */
	delete pdata;

	std::vector<int> coefDepCount(numcoefficients);
	//for (int i = 0; i < numcoefficients; ++i) { // computing coefficients
	//	coefDepCount[i] = 0;
	//	CondutanceCoefficientsVector ccv = condutanceCoefficients[i];
	//	for (int j = 0; j < ccv.size(); j++) {
	//		int row = ccv[j].row;
	//		int col = ccv[j].col;
	//		if (row == GROUNDNODE || col == GROUNDNODE) { // skipping ground node
	//			continue;
	//		}
	//		coefDepCount[i]++;
	//	}
	//}
	for (int i = 0; i<nodesCount; ++i) {
		for (nodeCoefficients *aux = nodeCoef[i]; aux != NULL; aux = aux->next) { 
			coefDepCount[aux->condIndex]++;
		}
	}
	// finding maximum dependency's count
	maxDepCount = -1;
	for (int i = 0; i < numcoefficients; i++) {
		if (coefDepCount[i] > maxDepCount) {
			maxDepCount = coefDepCount[i];
		}
	}
	cudaMalloc((void**)& auxv_d, maxDepCount * sizeof(numType));
	cudaMalloc((void**)& auxi_d, maxDepCount * sizeof(int));

	std::vector<int> csr2cpjds_map;
	std::vector<int> csr2cpjds_map_upper;
	std::vector<int> csr2cpjds_idx;
	std::vector<int> csr2cpjds_row;
	int elCount = 0;
	for (int col = 0; col < n; col++) { // column
		for (int row = col; row < n; row++) { // row
			int dataIdx = coordinates2Index(row, col);
			if (dataIdx >= 0) {
				csr2cpjds_map.push_back(dataIdx);
				csr2cpjds_idx.push_back(col); // nao esta sendo usado
				csr2cpjds_row.push_back(row);
				elCount++;
			}

			// upper triangular
			int dataIdxUpper = coordinates2Index(col, row);
			if (dataIdxUpper >= 0) {
				csr2cpjds_map_upper.push_back(dataIdxUpper);
			}
		}
	}

	MatrixCPJDS2CSR csrMap;
	csrMap.n = n;
	csrMap.nnz = elCount;
	csrMap.csr2cpjds = csr2cpjds_map;
	csrMap.csr2cpjds_upper = csr2cpjds_map_upper;
	csrMap.indices = csr2cpjds_idx;
	csrMap.row = csr2cpjds_row;
	(*M).csrMap = csrMap;

	return 0;
}

void MatrixCPJDSManager::deleteMatrixCPJDS(MatrixCPJDS M) {
	/* matrix data */
	cudaFree(M.matrixData.data);
	cudaFree(M.matrixData.indices);
	cudaFree(M.matrixData.rowLength);
	cudaFree(M.matrixData.rowSize);
	cudaFree(M.matrixData.colOffset);
	/* matrix colors (CPU MEMORY!) */
	//delete M.matrixColors.colors; // deleted in destructor
	delete M.matrixColors.colorsColOffsetSize;
	cudaFree(M.matrixColors.colors_d);
	cudaFree(M.matrixColors.colorsColOffsetSize_d);
	/* matrix dependencies*/
	cudaFree(M.matrixDependencies.dependencyRowDataIndex);
	cudaFree(M.matrixDependencies.dependencyDiagDataIndex);
	cudaFree(M.matrixDependencies.dependencyArrayInitialIndex);
	cudaFree(M.matrixDependencies.dependencyLowerDataIndex);
	cudaFree(M.matrixDependencies.dependencyUpperDataIndex);
	cudaFree(M.matrixDependencies.nnzElementDataIndex);
	cudaFree(M.matrixDependencies.dependenciesSize);
	/* preconditioner data */
	cudaFree(M.preconditionedData);
	/* aux data */
	delete M.cpuData.data;
	delete M.cpuData.precond;
	delete M.cpuData.indices;
};

/* sets an element value according to row and column indexes */
void MatrixCPJDSManager::set(MatrixCPJDS M, int row, int col, numType val) {
	int rowP = original2PaddedIdx[row];
	int colP = original2PaddedIdx[col];
	int idx = coordinates2Index(rowP, colP);
	cm_set_cpjds<<<1,1>>>(idx, M.matrixData.data, val);
}

/* increments an element value according to row and column indexes */
void MatrixCPJDSManager::increment(MatrixCPJDS M, int row, int col, numType val) {
	int rowP = original2PaddedIdx[row];
	int colP = original2PaddedIdx[col];
	int idx = coordinates2Index(rowP, colP);
	cm_increment_cpjds<<<1,1>>>(idx, M.matrixData.data, val);
}

/* increments an array of elements value according to elements' indexes */
void MatrixCPJDSManager::pushIncrements(MatrixCPJDS M, int size, numType * vals, int * indices) {
	cudaMemcpy(auxv_d, vals, (size_t)size * sizeof(numType), cudaMemcpyHostToDevice);
	cudaMemcpy(auxi_d, indices, (size_t)size * sizeof(int), cudaMemcpyHostToDevice);
	int blocksi = 1,
		blocksizes = size;
#ifdef DEBUG
	if (size > BLOCKSIZE) {
		blocksi = ceil((float)size / BLOCKSIZE);
		blocksizes = BLOCKSIZE;
	}
#endif
	cm_increment_cpjds <<<blocksi, blocksizes >>>(size, M.matrixData.data, auxv_d, auxi_d);
}

/* increments an array of elements value according to elements' indexes */
void MatrixCPJDSManager::pushIncrements(MatrixCPJDS M, int size, numType * vals, int * indices, cudaStream_t stream) {
	cudaMemcpyAsync(auxv_d, vals, (size_t)size * sizeof(numType), cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(auxi_d, indices, (size_t)size * sizeof(int), cudaMemcpyHostToDevice, stream);
	int blocksi = 1,
		blocksizes = size;
#ifdef DEBUG
	if (size > BLOCKSIZE) {
		blocksi = ceil((float)size / BLOCKSIZE);
		blocksizes = BLOCKSIZE;
	}
#endif
	cm_increment_cpjds <<<blocksi, blocksizes, 0, stream >>>(size, M.matrixData.data, auxv_d, auxi_d);
}

/* method for restoring CPJDS-transformed vector to its original size and indexes */
std::vector<numType> MatrixCPJDSManager::restore(Vector * v){
	// groundNode = -1: no ground node
	std::vector<numType> restored(nOrig);

	numType * data_h = new numType[this->n];
	cudaMemcpy(data_h, v->getData(), (size_t)this->n * sizeof(numType), cudaMemcpyDeviceToHost);

	for (int i = 0; i < nOrig; i++) {
		restored[i] = data_h[original2PaddedIdx[i]];
	}

	return restored;
}

/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
void MatrixCPJDSManager::mult(MatrixCPJDS M, Vector * x, Vector * b) {
#ifdef DEBUG
	if (x == NULL || b == NULL) {
		printf("mv_mult error: vector(s) null!\n");
		return;
	}
	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
		printf("mv_mult error: vector(s) of different size!\n");
		return;
	}
#endif
	//numType * mData = M.matrixData.data;
	//int * mIndices = M.matrixData.indices;
	//int * mRowLength = M.matrixData.rowLength;
	//int * mRowSize = M.matrixData.rowSize;
	//int * mColOffset = M.matrixData.colOffset;

	//int colorCount = M.matrixColors.colorCount;
	//int * colors = M.matrixColors.colors;
	//int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	//numType * xData = x->getData();
	//numType * bData = b->getData();

	//for (int k = 0; k < colorCount; k++) {

	//	int colorStart = colors[k];
	//	int colorEnd = colors[k + 1];
	//	int colorColOffset = mColorsColOffsetSize[k];

	//	int size = colorEnd - colorStart;

	//	// Launch CUDA matrix-vector multiplication kernel
	//	cmv_mult_cpjds <<<blocks, BLOCKSIZE>>>(size, mData, mIndices, mRowLength, mRowSize,
	//											mColOffset, colorStart, colorColOffset, xData, bData);
	cudaError_t cudaStatus;

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "mv_mult kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
	//}

	int size = M.matrixData.n;
	numType * mData = M.matrixData.data;//M.preconditionedData;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors_d = M.matrixColors.colors_d;
	int * mColorsColOffsetSize_d = M.matrixColors.colorsColOffsetSize_d;

	numType * xData = x->getData();
	numType * byData = b->getData();

	//// Launch CUDA matrix-vector multiplication kernel
	//cmv_mult <<<blocks, BLOCKSIZE, 0, stream >>>(mData, mIndices, mRowLength, mRowSize,
	//	mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, yData, pData);

	// Launch CUDA matrix-vector multiplication kernel
	cmv_mult_cpjds2 << <blocks, BLOCKSIZE>> >(size, mData, mIndices, mRowLength, mRowSize,
		mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, byData);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cout << "mv_mult kernel failed: " << cudaGetErrorString(cudaStatus) << std::endl;
	}
}

/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
void MatrixCPJDSManager::mult(MatrixCPJDS M, Vector * x, Vector * b, cudaStream_t * streams) {
#ifdef DEBUG
	if (x == NULL || b == NULL) {
		printf("mv_mult error: vector(s) null!\n");
		return;
	}
	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
		printf("mv_mult error: vector(s) of different size!\n");
		return;
	}
#endif
	numType * mData = M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors = M.matrixColors.colors;
	int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	numType * xData = x->getData();
	numType * bData = b->getData();

	for (int k = 0; k < colorCount; k++) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = mColorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA matrix-vector multiplication kernel
		cmv_mult_cpjds <<<blocks, BLOCKSIZE, 0, streams[k]>>>(size, mData, mIndices, mRowLength, mRowSize,
															mColOffset, colorStart, colorColOffset, xData, bData);

#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			printf("%s", cudaGetErrorString(cudaStatus));
			/*std::ostringstream msg;
			msg << "mv_mult kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
			msg.flush(); LOG(&msg.str()[0]);*/
		}
#endif
	}
}

/* Matrix-vector multiplication (stored) followed by an inner produtct: k = xt * M * x
* (totalization among blocks is NOT performed)
* where y = M * x, M: matrix; x, y: vectors, k: number (not provided) */
void MatrixCPJDSManager::multInner(MatrixCPJDS M, Vector * x, Vector * y, cudaStream_t * streams) {
//#ifdef DEBUG
//	if (x == NULL || b == NULL) {
//		printf("mv_mult error: vector(s) null!\n");
//		return;
//	}
//	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
//		printf("mv_mult error: vector(s) of different size!\n");
//		return;
//	}
//#endif
	numType * mData = M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors = M.matrixColors.colors;
	int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	numType * xData = x->getData();
	numType * yData = y->getData();
	numType * pData = y->getPartial();

	for (int k = 0; k < colorCount; k++) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = mColorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA matrix-vector multiplication kernel
		cmv_mult_inner_cpjds <<<blocks, BLOCKSIZE, 0, streams[k] >>>(size, mData, mIndices, mRowLength, mRowSize,
			mColOffset, colorStart, colorColOffset, xData, yData, pData);

#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			printf("%s", cudaGetErrorString(cudaStatus));
			/*std::ostringstream msg;
			msg << "mv_mult kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
			msg.flush(); LOG(&msg.str()[0]);*/
		}
#endif
	}
}

/* Matrix-vector multiplication (stored) followed by an inner produtct: k = xt * M * x
* (totalization among blocks is NOT performed)
* where y = M * x, M: matrix; x, y: vectors, k: number (not provided) */
void MatrixCPJDSManager::multInner2(MatrixCPJDS M, Vector * x, Vector * y, cudaStream_t stream) {
//#ifdef DEBUG
//	if (x == NULL || b == NULL) {
//		printf("mv_mult error: vector(s) null!\n");
//		return;
//	}
//	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
//		printf("mv_mult error: vector(s) of different size!\n");
//		return;
//	}
//#endif
	int size = M.matrixData.n;
	numType * mData = M.preconditionedData;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors_d = M.matrixColors.colors_d;
	int * mColorsColOffsetSize_d = M.matrixColors.colorsColOffsetSize_d;

	numType * xData = x->getData();
	numType * yData = y->getData();
	numType * pData = x->getPartial();

	//// Launch CUDA matrix-vector multiplication kernel
	//cmv_mult <<<blocks, BLOCKSIZE, 0, stream >>>(mData, mIndices, mRowLength, mRowSize,
	//	mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, yData, pData);

	// Launch CUDA matrix-vector multiplication kernel
	cmv_mult_2 <<<blocks, BLOCKSIZE, 0, stream >>>(size, mData, mIndices, mRowLength, mRowSize,
		mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, yData, pData);

//#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			printf("%s", cudaGetErrorString(cudaStatus));
			/*std::ostringstream msg;
			msg << "mv_mult kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
			msg.flush(); LOG(&msg.str()[0]);*/
		}
//#endif
}

/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
// L must be lower triangular
void MatrixCPJDSManager::solve(MatrixCPJDS M, Vector * b, Vector * x) {
#ifdef DEBUG
	if (x == NULL || b == NULL) {
		printf("mv_solve error: vector(s) null!\n");
		return;
	}
	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
		printf("mv_solve error: vector(s) of different size!\n");
		return;
	}
#endif

	//int size = M.n;
	numType * mData = M.preconditionedData;//; M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors = M.matrixColors.colors;
	int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	numType * xData = x->getData();
	numType * bData = b->getData();

	for (int k = 0; k < colorCount; k++) {
		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = mColorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_cpjds << <blocks, BLOCKSIZE >> >(size, mData, mIndices, mRowLength, mRowSize, 
													mColOffset, colorStart, colorColOffset, xData, bData);

#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			std::ostringstream msg;
			msg << "mv_solve kernel failed: " << cudaGetErrorString(cudaStatus) << ", it:  " << k << "\n";
			msg.flush(); LOG(&msg.str()[0]);
		}
#endif
	}
}

/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
// L must be lower triangular
void MatrixCPJDSManager::solve(MatrixCPJDS M, Vector * b, Vector * x, cudaStream_t stream) {
#ifdef DEBUG
	if (x == NULL || b == NULL) {
		printf("mv_solve error: vector(s) null!\n");
		return;
	}
	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
		printf("mv_solve error: vector(s) of different size!\n");
		return;
	}
#endif

	//int size = M.n;
	numType * mData = M.preconditionedData; // M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors = M.matrixColors.colors;
	int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	numType * xData = x->getData();
	numType * bData = b->getData();

	for (int k = 0; k < colorCount; k++) {
		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = mColorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_cpjds <<<blocks, BLOCKSIZE, 0, stream>>>(size, mData, mIndices, mRowLength, mRowSize, 
													mColOffset, colorStart, colorColOffset, xData, bData);

#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			std::ostringstream msg;
			msg << "mv_solve kernel failed: " << cudaGetErrorString(cudaStatus) << ", it:  " << k << "\n";
			msg.flush(); LOG(&msg.str()[0]);
		}
#endif
	}
}

/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
// U must be an upper triangular matrix
void MatrixCPJDSManager::solve_t(MatrixCPJDS M, Vector * b, Vector * x) {
#ifdef DEBUG
	if (x == NULL || b == NULL) {
		printf("mv_solve error: vector(s) null!\n");
		return;
	}
	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
		printf("mv_solve error: vector(s) of different size!\n");
		return;
	}
#endif

	//int size = M.n;
	numType * mData = M.preconditionedData; //M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors = M.matrixColors.colors;
	int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	numType * xData = x->getData();
	numType * bData = b->getData();

	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = mColorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_t_cpjds << <blocks, BLOCKSIZE >> >(size, mData, mIndices, mRowLength, mRowSize, 
													mColOffset, colorStart, colorColOffset, xData, bData);

#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			std::ostringstream msg;
			msg << "mv_solve_t kernel failed: " << cudaGetErrorString(cudaStatus) << ", it:  " << k << "\n";
			msg.flush(); LOG(&msg.str()[0]);
		}
#endif
	}
}

/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
// U must be an upper triangular matrix
void MatrixCPJDSManager::solve_t(MatrixCPJDS M, Vector * b, Vector * x, cudaStream_t stream) {
#ifdef DEBUG
	if (x == NULL || b == NULL) {
		printf("mv_solve error: vector(s) null!\n");
		return;
	}
	if (M.matrixData.n != x->getSize() || M.matrixData.n != b->getSize()) {
		printf("mv_solve error: vector(s) of different size!\n");
		return;
	}
#endif

	//int size = M.n;
	numType * mData = M.preconditionedData; //M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors = M.matrixColors.colors;
	int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	numType * xData = x->getData();
	numType * bData = b->getData();

	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = mColorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_t_cpjds <<<blocks, BLOCKSIZE, 0, stream >>>(size, mData, mIndices, mRowLength, mRowSize, 
													mColOffset, colorStart, colorColOffset, xData, bData);

#ifdef DEBUG
		// Check for any errors launching the kernel
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			std::ostringstream msg;
			msg << "mv_solve_t kernel failed: " << cudaGetErrorString(cudaStatus) << ", it:  " << k << "\n";
			msg.flush(); LOG(&msg.str()[0]);
		}
#endif
	}
}

/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
void MatrixCPJDSManager::solve_complete(MatrixCPJDS M, Vector * b, Vector * x, cudaStream_t stream) {
	numType * mData = M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors_d = M.matrixColors.colors_d;
	int * mColorsColOffsetSize_d = M.matrixColors.colorsColOffsetSize_d;

	numType * xData = x->getData();
	numType * bData = b->getData();

	// Launch CUDA (lower) triangular solver kernel
	cmv_solve <<<blocks, BLOCKSIZE, 0, stream >>>(mData, mIndices, mRowLength, mRowSize,
		mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, bData);
}

/* multiple operations
* full solver (lower + upper triangulars) M * x = b => x = inv(M) * b
* inner product b.x (x's partials are filled) */
void MatrixCPJDSManager::solve_and_inner(MatrixCPJDS M, Vector * b, Vector * x, Vector * u, cudaStream_t stream) {
	numType * mData = M.matrixData.data;
	int * mIndices = M.matrixData.indices;
	int * mRowLength = M.matrixData.rowLength;
	int * mRowSize = M.matrixData.rowSize;
	int * mColOffset = M.matrixData.colOffset;

	int colorCount = M.matrixColors.colorCount;
	int * colors_d = M.matrixColors.colors_d;
	int * mColorsColOffsetSize_d = M.matrixColors.colorsColOffsetSize_d;

	numType * xData = x->getData();
	numType * xParcial = x->getPartial();
	numType * bData = b->getData();
	numType * uData = u->getData();

	// Launch CUDA (lower) triangular solver kernel
	cmv_solve_inner <<<blocks, BLOCKSIZE, 0, stream >>>(mData, mIndices, mRowLength, mRowSize,
		mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, bData, uData, xParcial);
}

void MatrixCPJDSManager::saveToFile(char * filename, MatrixCPJDS M, numType * data, bool isCPU) {

	std::ofstream file(filename);
	if (!isCPU) {
		// copy from CUDA device memory to host memory
		numType * data_h = new numType[M.matrixData.elCount];
		cudaMemcpy(data_h, data, (size_t)M.matrixData.elCount * sizeof(numType), cudaMemcpyDeviceToHost);
		data = data_h;
	}

	for (int row = 0; row < M.matrixData.n; row++) {
		for (int col = 0; col < M.matrixData.n; col++) {
			int idx = coordinates2Index(row, col);
			numType val = 0;
			if (idx > -1) {
				int dcol = M.cpuData.indices[idx];
				val = M.cpuData.data[idx];
			}
			file << val << " ";
		}
		file << "\n";
	}
	file.close();

	if (!isCPU) {
		delete data;
	}
}

/*****************************************************************************************/
/*****************************************************************************************/
/***********************************    kernel    ****************************************/
/*****************************************************************************************/
/*****************************************************************************************/

/* CPJDS matrix element setter */
__global__ void cm_set_cpjds(int idx, numType * data, numType val) {
	// row index
	data[idx] = val;
}

/* CPJDS matrix element incrementer */
__global__ void cm_increment_cpjds(int idx, numType * data, numType val) {
	// row index
	data[idx] += val;
}

/* CPJDS matrix batch incrementer */
__global__ void cm_increment_cpjds(int size, numType * data, numType * vals, int * indices) {
	// row index
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < size) {
		data[indices[idx]] = vals[idx];
	}
}

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * bData) {
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int rowSize = aRowSize[row];

		numType sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * xData[idx]; // NOT coalesced!
			//sum += mData[offset] * xData[mIndices[offset]];

			__syncthreads(); // synchronization so all threads load from memory
		}
		bData[row] = sum; // coalesced
	}

	__syncthreads();
}

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds2(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * bData) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	if (row < dim) {
		int colorIdx = 0;
		for (; colorIdx < colorCount; colorIdx++) {
			if (row < colors[colorIdx]) {
				break;
			}
		}
		colorIdx--;
		int colorStart = colors[colorIdx];
		int colorColOffset = colorsColOffset[colorIdx];

		// row size (length + padding zeros)
		int rowSize = aRowLength[row];

		numType sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * xData[idx]; // NOT coalesced!
			//if (row == 5) {
			//	printf("%d\t%d\t%g\n", row, j, xData[idx]);
			//}
			//sum += rowData;

			__syncthreads(); // synchronization so all threads load from memory
		}
		bData[row] = sum; // coalesced
	}

	__syncthreads();
}

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_inner_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * yData, numType * partial) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int rowSize = aRowSize[row];

		numType sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			sum += rowData * xData[idx]; // NOT coalesced!
			//sum += mData[offset] * xData[mIndices[offset]];

			__syncthreads(); // synchronization so all threads load from memory
		}
		yData[tidx] = sum; // coalesced
		cache[cacheIndex] = sum * xData[row];
	}
	else {
		cache[cacheIndex] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (cacheIndex < half) {
			cache[cacheIndex] += cache[cacheIndex + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (cacheIndex == 0) {
		partial[blockIdx.x] = cache[0];
	}

	__syncthreads();

	// isto precisa ser enviado para outro ponto, para que as somas parciais sejam totalizadas
	//if (row == 0) {
	//	numType sum = 0;
	//	for (int i = 0; i < blocks; i++) {
	//		sum += r[i];
	//	}
	//	val[0] = sum;
	//}
}

// backup do kernel cmv_mult, modificado abaixo
/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult(numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * yData, numType * partial) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	cache[cacheIndex] = 0.0;
	for (int k = 0; k < colorCount; k++) {
		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];

		if (row >= colorStart && row < colorEnd) {
			// row size (length + padding zeros)
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();
			for (int j = 0; j < rowSize; j++) {
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

				numType rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				sum += rowData * xData[idx]; // NOT coalesced!

				__syncthreads(); // synchronization so all threads load from memory
			}
			yData[row] = sum; // coalesced
			cache[cacheIndex] = sum * xData[row];
		}

		__syncthreads();

		int half = BLOCKSIZE >> 1;
		while (half != 0){
			if (cacheIndex < half) {
				cache[cacheIndex] += cache[cacheIndex + half];
			}

			__syncthreads();

			half >>= 1;
		}

		if (cacheIndex == 0) {
			partial[blockIdx.x] = cache[0];
		}

		__syncthreads();
	}
}

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_2(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * yData, numType * partial) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	cache[cacheIndex] = 0.0;

	if (row < dim) {
		int colorIdx = 0;
		for (; colorIdx < colorCount; colorIdx++) {
			if (row < colors[colorIdx]) {
				break;
			}
		}
		int colorStart = colors[colorIdx];
		int colorColOffset = colorsColOffset[colorIdx];

		// row size (length + padding zeros)
		int rowSize = aRowSize[row];

		numType sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
		//	sum += rowData * xData[idx]; // NOT coalesced!
			//if (row == 0 || row == 1045 || row == 1151) {
			//	printf("%d\t%d\t%g\n", row, j, xData[idx]);
			//}
			sum += rowData;

			__syncthreads(); // synchronization so all threads load from memory
		}
		yData[row] = sum; // coalesced
		cache[cacheIndex] = sum * xData[row];
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (cacheIndex < half) {
			cache[cacheIndex] += cache[cacheIndex + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (cacheIndex == 0) {
		partial[blockIdx.x] = cache[0];
	}

	__syncthreads();
}

/* CPJDS lower triangular solver */
// does not work per se, in parallel
__global__ void cmv_solve_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * bData) {
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int rowSize = aRowLength[row];

		numType sum = 0;
		__syncthreads();

		for (int j = 1; j < rowSize; j++) { // first element is main diagonal
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			if (idx < row) { // main diagonal can be skiped
				sum += rowData * xData[idx];
			}
			__syncthreads();
		}
		xData[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + tidx];
	}

	__syncthreads();
}

/* CPJDS upper triangular solver */
// does not work per se, in parallel
__global__ void cmv_solve_t_cpjds(int dim, numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, numType * xData, numType * bData) {
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int rowSize = aRowLength[row];

		numType sum = 0;
		__syncthreads();

		for (int j = 1; j < rowSize; j++) { // first idx is main diagonal
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

			numType rowData = aData[offset]; // coalesced
			int idx = aIndices[offset]; // coalesced
			if (idx > row && idx > -1) { // main diagonal can be skiped
				sum += rowData * xData[idx];
			}
			__syncthreads();
		}
		xData[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + tidx];
	}

	__syncthreads();
}

/* CPJDS lower and upper triangular solver */
__global__ void cmv_solve(numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * bData) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	for (int k = 0; k < colorCount; k++) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];
		int inColorIdx = (row - colorStart);

		if (row >= colorStart && row < colorEnd) {
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				numType rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
				__syncthreads();
			}
			xData[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + inColorIdx];
		}
	}

	__syncthreads();

	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];
		int inColorIdx = (row - colorStart);

		if (row >= colorStart && row < colorEnd) {
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				numType rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
				__syncthreads();
			}
			xData[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + inColorIdx];
		}
	}

	__syncthreads();

}

/* CPJDS lower and upper triangular solver followed by inner product */
__global__ void cmv_solve_inner(numType * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, numType * xData, numType * bData, numType * interm, numType * partial) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	cache[cacheIndex] = 0.0;
	interm[cacheIndex] = 0.0;
	for (int k = 0; k < colorCount; k++) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];
		int inColorIdx = (row - colorStart);

		if (row >= colorStart && row < colorEnd) {
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				numType rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
				__syncthreads();
			}
			interm[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + inColorIdx];
		}
	}

	__syncthreads();

	for (int k = colorCount - 1; k >= 0; k--) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];
		int inColorIdx = (row - colorStart);

		if (row >= colorStart && row < colorEnd) {
			int rowSize = aRowSize[row];

			numType sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				numType rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
				__syncthreads();
			}
			//numType bVal = bData[row];
			numType bVal = interm[row]; // result from previous (lower triangular) solver
			numType xVal = (bVal - sum) / aData[aColOffset[colorColOffset] + inColorIdx];
			xData[row] = xVal;
			cache[cacheIndex] = xVal * bVal;
		}
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (cacheIndex < half) {
			cache[cacheIndex] += cache[cacheIndex + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (cacheIndex == 0) {
		partial[blockIdx.x] = cache[0];
	}

	__syncthreads();

}

#endif