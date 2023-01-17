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

void DeleterCudaIntPtr::operator()(int* ptr) { cudaFree(ptr); };
void DeleterCudaIntPtr::operator()(double* ptr) { cudaFree(ptr); };

//void dependencies_analysis(double * data, int n, DependeciesMap * dependenciesMap) {
void dependencies_analysis(std::unique_ptr<double[]> &data, int n, DependeciesMap * dependenciesMap) {
	//for (int i = 0; i < n; i++) { // for each row...
	//	for (int j = i + 1; j < n; j++) { // ... for every (upper triangle) column in the row
	//		//for (int j = 0; j < i; j++) { // ... for every (lower triangle) column in the row
	//		if (MOD(data[i * n + j]) > EPS) { // ... if element (i, j) is non-zero, search dependencies
	//			Dependencies d;
	//			for (int k = 0; k < i; k++) { // check all rows before current row
	//				//for (int k = 0; k < j; k++) { // check all columns before current column
	//				// check if (k, i) and (k, j) are both non-zero
	//				if (MOD(data[k * n + i]) > EPS && MOD(data[k * n + j]) > EPS) {
	//					d.dependencies.push_back(k); // push dependency (only row in lower triangular is needed)
	//				}
	//			}
	//			d.lower = std::pair<int, int>(i, j);
	//			d.dependenciesSize = d.dependencies.size();
	//			(*dependenciesMap)[i].push_back(d);
	//		}
	//	}
	//}
}

/* CPJDS matrix element setter */
__global__ void cm_set_cpjds(int idx, double * data, double val);

/* CPJDS matrix element incrementer */
__global__ void cm_increment_cpjds(int idx, double * data, double val);

/* CPJDS matrix batch incrementer */
__global__ void cm_increment_cpjds(int size, double * data, double * vals, int * indices);

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * bData);

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds2(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * bData);

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_inner_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * yData, double * partial);

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult(double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorColOffset, double * xData, double * yData, double * partial);

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_2(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * yData, double * partial);

/* CPJDS lower triangular solver */
__global__ void cmv_solve_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * bData);

/* CPJDS upper triangular solver */
__global__ void cmv_solve_t_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * bData);

/* CPJDS lower and upper triangular solver */
__global__ void cmv_solve(double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * bData);

/* CPJDS lower and upper triangular solver */
__global__ void cmv_solve_inner(double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * bData, double * interm, double * partial);

/*
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
	leftShiftMatrix(rowsL, rowsU);

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
	std::unique_ptr<int[]> colorOffsetCount(new int[colorCount + 1]);
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
	std::unique_ptr<int[]> colOffset(new int[offsetSize]);

	// data and indices arrays
	std::unique_ptr<double[]> mdata(new double[total]);
	std::unique_ptr<int[]> indices(new int[total]);
	createDataAndIndicesVectors(mdata.get(), indices.get(), colOffset.get(), colorOffsetCount.get(), rowsL, rowsU, padding);

	/* computing (x,y)=>IDX map */
	//std::vector<std::unique_ptr<std::map<int, int>>> rowCol2IdxMap(n);
	this->coord2IndexMap.resize(n);
	for (int i = 0; i < n; i++) this->coord2IndexMap[i] = std::unique_ptr<std::map<int, int>>(new std::map<int, int>());
	for (int k = 0; k < colorCount; k++) { // color block
		for (int i = colors[k]; i < colors[k + 1]; i++) { // iterating rows in such color
			for (int j = 0; j < rowLength[i]; j++) {
				int idx = colOffset[colorOffsetCount[k] + j] + i - colors[k];
				this->coord2IndexMap[i]->insert(std::pair<int, int>(indices[idx], idx));
			}
		}
	}

	///* computing dependencies */
	//DependeciesMap dependencies(n);
	//dependencies_analysis2(n, &dependencies);

	//// fix dependencies' indices
	////std::cout << "\nFixing dependencies..." << std::endl;
	//DependeciesMap dependenciesCPJDS(n);
	//for (int i = 0; i < n; i++) {
	//	std::vector<Dependencies> dv = dependencies[i];

	//	for (int j = 0; j < dv.size(); j++) {
	//		Dependencies d = dv[j];
	//		// this code was removed, as all non-zeroes must be listed, regardless of having dependencies or not
	//		//if (d.dependencies.empty()) {
	//		//	std::cout << "skipping dependency..." << std::endl;
	//		//	continue;
	//		//}
	//		Dependencies dCPJDS;
	//		int coli = d.lower.first;
	//		int colj = d.lower.second;
	//		//dCPJDS.lower.first = coli % n;
	//		//dCPJDS.lower.second = colj % n;
	//		for (int k = 0; k < d.dependencies.size(); k++) {
	//			int row = d.dependencies[k];
	//			// convert to data array index
	//			int idxi = coordinates2Index(row, coli);
	//			int idxj = coordinates2Index(row, colj);
	//			dCPJDS.dependencies.push_back(idxi);
	//			dCPJDS.dependencies.push_back(idxj);
	//		}
	//		dCPJDS.lower = std::pair<int, int>(coli, colj);
	//		dCPJDS.dependenciesSize = d.dependenciesSize;
	//		dependenciesCPJDS[i].push_back(dCPJDS);
	//	}
	//}

	//// make all dependecies' arrays for the same row have the same size
	//// (it is supposed that dependencies's array size is smaller than warp size - however, warp will not be filled)
	////std::vector<std::vector<std::deque<int>>> rowsDeps(n);
	//for (int i = 0; i < n; i++) { // iterating rows
	//	std::vector<Dependencies> * dvCPJDS = &(dependenciesCPJDS[i]);

	//	int maxSize = 0;
	//	for (int j = 0; j < (*dvCPJDS).size(); j++) { // for each non-zero element with dependency..
	//		int thisDepSize = (*dvCPJDS)[j].dependencies.size();
	//		if (thisDepSize > maxSize) { // ...check size, to find largest...
	//			maxSize = thisDepSize;
	//		}
	//	}
	//	for (int j = 0; j < (*dvCPJDS).size(); j++) { // ... and pad the others
	//		int thisDepSize = (*dvCPJDS)[j].dependencies.size();
	//		for (int k = thisDepSize; k < maxSize; k++) {
	//			(*dvCPJDS)[j].dependencies.push_back(-1);
	//		}
	//	}
	//}

	//// create dependencies's arrays in CPJDS format
	//std::vector<int> depsArrDiag;
	//std::vector<int> depsArrRow;
	//std::vector<int> depsIdxArr(n);
	//std::vector<int> depsLowerArr;
	//std::vector<int> depsUpperArr;
	//std::vector<int> depsNNZArr(n + 1);
	//std::vector<int> depsSizeArr(n);
	//for (int i = 0; i < n; i++) { // iterating rows
	//	std::vector<Dependencies> dvCPJDS = dependenciesCPJDS[i];

	//	// initial index
	//	depsIdxArr[i] = depsArrRow.size(); // size of depsArrRow = size of depsArrDiag

	//	int depsSize = 0;
	//	if (dvCPJDS.size() > 0) {
	//		depsSize = dvCPJDS[0].dependencies.size() / 2; // all dvCPJDS have same size (padded in previous step)
	//													   // also, dependenciesCPJDS.dependencies array is doubled as it contains both row and diagonal elements
	//	}
	//	// size
	//	depsSizeArr[i] = depsSize;

	//	// adding elements to "global" array, column-major
	//	for (int j = 0; j < depsSize; j++) {
	//		for (int k = 0; k < dvCPJDS.size(); k++) {
	//			depsArrDiag.push_back(dvCPJDS[k].dependencies[2 * j]);
	//			depsArrRow.push_back(dvCPJDS[k].dependencies[2 * j + 1]);
	//		}
	//	}

	//	depsNNZArr[i] = depsLowerArr.size(); // size of depsLowerArr = size of depsUpperArr
	//										 // element's indices (lower and upper)
	//	for (int j = 0; j < dvCPJDS.size(); j++) {
	//		depsLowerArr.push_back((*coord2IndexMap[dvCPJDS[j].lower.first])[dvCPJDS[j].lower.second]);
	//		depsUpperArr.push_back((*coord2IndexMap[dvCPJDS[j].lower.second])[dvCPJDS[j].lower.first]);
	//	}
	//}
	//depsNNZArr[n] = depsLowerArr.size(); // size of depsLowerArr = size of depsUpperArr

	//int * dependencyRowDataIndex = new int[depsArrRow.size()];
	//int * dependencyDiagDataIndex = new int[depsArrDiag.size()];
	//// dependency's array initial index (index in dependencyRowDataIndex and dependencyDiagDataIndex)
	//int * dependencyArrayInitialIndex = new int[depsIdxArr.size()];
	//// lower and upper triangulars non-zero elements (with dependencies) index in data array
	//int * dependencyLowerDataIndex = new int[depsLowerArr.size()];
	//int * dependencyUpperDataIndex = new int[depsUpperArr.size()];
	//// data array index for lower and upper triangular's elements (index in dependencyLowerDataIndex and dependencyUpperDataIndex)
	//int * nnzElementDataIndex = new int[depsNNZArr.size()];
	//// dependency's count
	//int * dependenciesSize = new int[depsSizeArr.size()];
	// preconditioner data
	double * pdata = new double[total];

	//// copy arrays
	//for (int i = 0; i < depsArrRow.size(); i++) dependencyRowDataIndex[i] = depsArrRow[i];
	//for (int i = 0; i < depsArrDiag.size(); i++) dependencyDiagDataIndex[i] = depsArrDiag[i];
	//for (int i = 0; i < depsIdxArr.size(); i++) dependencyArrayInitialIndex[i] = depsIdxArr[i];
	//for (int i = 0; i < depsLowerArr.size(); i++) dependencyLowerDataIndex[i] = depsLowerArr[i];
	//for (int i = 0; i < depsUpperArr.size(); i++) dependencyUpperDataIndex[i] = depsUpperArr[i];
	//for (int i = 0; i < depsNNZArr.size(); i++) nnzElementDataIndex[i] = depsNNZArr[i];
	//for (int i = 0; i < depsSizeArr.size(); i++) dependenciesSize[i] = depsSizeArr[i];
	for (int i = 0; i < total; i++) pdata[i] = 0;

	//int depsSize = depsArrRow.size(), // = depsArrDiag.size()
	//	idxSize = depsLowerArr.size(); // = depsUpperArr.size()

	/* matrix data */
	std::unique_ptr<double[], DeleterCudaIntPtr> c_mdata(new double[total]);
	std::unique_ptr<int[], DeleterCudaIntPtr> c_indices(new int[total]);
	std::unique_ptr<int[], DeleterCudaIntPtr> c_rowLength(new int[n]);
	std::unique_ptr<int[], DeleterCudaIntPtr> c_rowSize(new int[n]);
	std::unique_ptr<int[], DeleterCudaIntPtr> c_colOffset(new int[offsetSize]);
	/* matrix colors */
	std::unique_ptr<int[], DeleterCudaIntPtr> c_colors(new int[colorCount + 1]);
	std::unique_ptr<int[], DeleterCudaIntPtr> c_colorsColOffsetSize(new int[colorCount]);
	///* matrix dependencies*/
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_dependencyRowDataIndex(new int[depsSize]);
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_dependencyDiagDataIndex(new int[depsSize]);
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_dependencyLowerDataIndex(new int[idxSize]);
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_dependencyUpperDataIndex(new int[idxSize]);
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_dependencyArrayInitialIndex(new int[n]);
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_dependenciesSize(new int[n]);
	//std::unique_ptr<int[], DeleterCudaIntPtr> c_nnzElementDataIndex(new int[n + 1]);
	/* preconditioner data */
	std::unique_ptr<double[], DeleterCudaIntPtr> c_pdata(new double[total]);

	// CUDA device memory allocation
	/* matrix data */
	cudaMalloc((void**)& c_mdata, total * sizeof(double));
	cudaMalloc((void**)& c_indices, total * sizeof(int));
	cudaMalloc((void**)& c_rowLength, n * sizeof(int));
	cudaMalloc((void**)& c_rowSize, n * sizeof(int));
	cudaMalloc((void**)& c_colOffset, offsetSize * sizeof(int));
	/* matrix colors */
	cudaMalloc((void**)& c_colors, (colorCount + 1) * sizeof(int));
	cudaMalloc((void**)& c_colorsColOffsetSize, colorCount * sizeof(int));
	///* matrix dependencies*/
	//cudaMalloc((void**)& c_dependencyRowDataIndex, depsSize * sizeof(int));
	//cudaMalloc((void**)& c_dependencyDiagDataIndex, depsSize * sizeof(int));
	//cudaMalloc((void**)& c_dependencyLowerDataIndex, idxSize * sizeof(int));
	//cudaMalloc((void**)& c_dependencyUpperDataIndex, idxSize * sizeof(int));
	//cudaMalloc((void**)& c_dependencyArrayInitialIndex, n * sizeof(int));
	//cudaMalloc((void**)& c_dependenciesSize, n * sizeof(int));
	//cudaMalloc((void**)& c_nnzElementDataIndex, (n + 1) * sizeof(int));
	/* preconditioner data */
	cudaMalloc((void**)& c_pdata, total * sizeof(double));

	// CUDA device memory transfer
	/* matrix data */
	cudaMemcpy(c_mdata.get(), mdata.get(), (size_t)total * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(c_indices.get(), indices.get(), (size_t)total * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(c_rowLength.get(), rowLength, (size_t)n * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(c_rowSize.get(), rowSize, (size_t)n * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(c_colOffset.get(), colOffset.get(), (size_t)offsetSize * sizeof(int), cudaMemcpyHostToDevice);
	/* matrix colors */
	cudaMemcpy(c_colors.get(), colors.get(), (size_t)(colorCount + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(c_colorsColOffsetSize.get(), colorOffsetCount.get(), (size_t)colorCount * sizeof(int), cudaMemcpyHostToDevice);
	///* matrix dependencies*/
	//cudaMemcpy(c_dependencyRowDataIndex.get(), dependencyRowDataIndex, (size_t)depsSize * sizeof(int), cudaMemcpyHostToDevice);
	//cudaMemcpy(c_dependencyDiagDataIndex.get(), dependencyDiagDataIndex, (size_t)depsSize * sizeof(int), cudaMemcpyHostToDevice);
	//cudaMemcpy(c_dependencyLowerDataIndex.get(), dependencyLowerDataIndex, (size_t)idxSize * sizeof(int), cudaMemcpyHostToDevice);
	//cudaMemcpy(c_dependencyUpperDataIndex.get(), dependencyUpperDataIndex, (size_t)idxSize * sizeof(int), cudaMemcpyHostToDevice);
	//cudaMemcpy(c_dependencyArrayInitialIndex.get(), dependencyArrayInitialIndex, (size_t)n * sizeof(int), cudaMemcpyHostToDevice);
	//cudaMemcpy(c_dependenciesSize.get(), dependenciesSize, (size_t)n * sizeof(int), cudaMemcpyHostToDevice);
	//cudaMemcpy(c_nnzElementDataIndex.get(), nnzElementDataIndex, (size_t)(n + 1) * sizeof(int), cudaMemcpyHostToDevice);
	/* preconditioner data */
	cudaMemcpy(c_pdata.get(), pdata, (size_t)total * sizeof(double), cudaMemcpyHostToDevice);

	/* set matrix data */
	(*M).matrixData.data = std::move(c_mdata);
	(*M).matrixData.indices = std::move(c_indices);
	(*M).matrixData.rowLength = std::move(c_rowLength);
	(*M).matrixData.rowSize = std::move(c_rowSize);
	(*M).matrixData.colOffset = std::move(c_colOffset);
	(*M).matrixData.n = n;
	(*M).matrixData.nnz = nnz;
	(*M).matrixData.elCount = total;
	(*M).matrixData.offsetCount = offsetSize;
	/* set matrix colors (CPU MEMORY!) */
	(*M).matrixColors.colors = colors;
	(*M).matrixColors.colorCount = colorCount;
	(*M).matrixColors.colorsColOffsetSize = std::move(colorOffsetCount);
	(*M).matrixColors.colors_d = std::move(c_colors);
	(*M).matrixColors.colorsColOffsetSize_d = std::move(c_colorsColOffsetSize);
	///* set matrix dependencies*/
	//(*M).matrixDependencies.dependencyRowDataIndex = std::move(c_dependencyRowDataIndex);
	//(*M).matrixDependencies.dependencyDiagDataIndex = std::move(c_dependencyDiagDataIndex);
	//(*M).matrixDependencies.dependencyLowerDataIndex = std::move(c_dependencyLowerDataIndex);
	//(*M).matrixDependencies.dependencyUpperDataIndex = std::move(c_dependencyUpperDataIndex);
	//(*M).matrixDependencies.dependencyArrayInitialIndex = std::move(c_dependencyArrayInitialIndex);
	//(*M).matrixDependencies.dependenciesSize = std::move(c_dependenciesSize);
	//(*M).matrixDependencies.nnzElementDataIndex = std::move(c_nnzElementDataIndex);
	//(*M).matrixDependencies.depsSize = depsSize;
	//(*M).matrixDependencies.idxSize = idxSize;
	/* set preconditioner data */
	(*M).preconditionedData = std::move(c_pdata);
	// OBS: is color-data needed in GPU memory?

	// setting CPU aux data
	(*M).cpuData.data = std::move(mdata);
	(*M).cpuData.indices = std::move(indices);//indices;
	(*M).cpuData.precond = std::unique_ptr<double[]>(new double[total]); for (int i = 0; i < total; i++) (*M).cpuData.precond[i] = 0;

	// free CPU memory used for initializing GPU data
	/* matrix data */
	delete rowLength;
	delete rowSize;
	///* matrix dependencies*/
	//delete dependencyRowDataIndex;
	//delete dependencyDiagDataIndex;
	//delete dependencyArrayInitialIndex;
	//delete dependencyLowerDataIndex;
	//delete dependencyUpperDataIndex;
	//delete nnzElementDataIndex;
	//delete dependenciesSize;
	/* preconditioner data */
	delete pdata;

	//std::vector<int> csr2cpjds_map;
	//std::vector<int> csr2cpjds_map_upper;
	//std::vector<int> csr2cpjds_idx;
	//std::vector<int> csr2cpjds_row;
	//int elCount = 0;
	//for (int col = 0; col < n; col++) { // column
	//	for (int row = col; row < n; row++) { // row
	//		int dataIdx = coordinates2Index(row, col);
	//		if (dataIdx >= 0) {
	//			csr2cpjds_map.push_back(dataIdx);
	//			csr2cpjds_idx.push_back(col); // nao esta sendo usado
	//			csr2cpjds_row.push_back(row);
	//			elCount++;
	//		}

	//		// upper triangular
	//		int dataIdxUpper = coordinates2Index(col, row);
	//		if (dataIdxUpper >= 0) {
	//			csr2cpjds_map_upper.push_back(dataIdxUpper);
	//		}
	//	}
	//}

	MatrixCPJDS2CSR csrMap;
	createCsr2CpjdsMap(csrMap);
	//csrMap.n = n;
	//csrMap.nnz = elCount;
	//csrMap.csr2cpjds = csr2cpjds_map;
	//csrMap.csr2cpjds_upper = csr2cpjds_map_upper;
	//csrMap.indices = csr2cpjds_idx;
	//csrMap.row = csr2cpjds_row;
	(*M).csrMap = csrMap;

	return 0;
}

/*
* data: full symmetric matrix (n-by-n), with zeroes - data is row-major!!!
* n: matrix size - must be multiple of WARP_SIZE (32)
* colors: array of each colors offset (size is colorCount + 1, last position being equal to n)
* colorCount: number of colors
* M: resulting matrix in colored pJDS format
*/
int MatrixCPJDSManager::buidMatrixCPJDS(MatrixCPJDS * M, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients) {
	int ans;
	if ((ans = buidMatrixCPJDS(M)) != 0) return ans;

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
	cudaMalloc((void**)& auxv_d, maxDepCount * sizeof(double));
	cudaMalloc((void**)& auxi_d, maxDepCount * sizeof(int));
	return 0;
}

/* sets an element value according to row and column indexes */
void MatrixCPJDSManager::set(MatrixCPJDS &M, int row, int col, double val) {
	int rowP = original2PaddedIdx[row];
	int colP = original2PaddedIdx[col];
	int idx = coordinates2Index(rowP, colP);
	cm_set_cpjds<<<1,1>>>(idx, M.matrixData.data.get(), val);
}

/* increments an element value according to row and column indexes */
void MatrixCPJDSManager::increment(MatrixCPJDS &M, int row, int col, double val) {
	int rowP = original2PaddedIdx[row];
	int colP = original2PaddedIdx[col];
	int idx = coordinates2Index(rowP, colP);
	cm_increment_cpjds<<<1,1>>>(idx, M.matrixData.data.get(), val);
}

/* increments an array of elements value according to elements' indexes */
void MatrixCPJDSManager::pushIncrements(MatrixCPJDS &M, int size, double * vals, int * indices) {
	cudaMemcpy(auxv_d.get(), vals, (size_t)size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(auxi_d.get(), indices, (size_t)size * sizeof(int), cudaMemcpyHostToDevice);
	int blocksi = 1,
		blocksizes = size;
#ifdef DEBUG
	if (size > BLOCKSIZE) {
		blocksi = ceil((float)size / BLOCKSIZE);
		blocksizes = BLOCKSIZE;
	}
#endif
	cm_increment_cpjds <<<blocksi, blocksizes >>>(size, M.matrixData.data.get(), auxv_d.get(), auxi_d.get());
}

/* increments an array of elements value according to elements' indexes */
void MatrixCPJDSManager::pushIncrements(MatrixCPJDS &M, int size, double * vals, int * indices, cudaStream_t stream) {
	cudaMemcpyAsync(auxv_d.get(), vals, (size_t)size * sizeof(double), cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(auxi_d.get(), indices, (size_t)size * sizeof(int), cudaMemcpyHostToDevice, stream);
	int blocksi = 1,
		blocksizes = size;
#ifdef DEBUG
	if (size > BLOCKSIZE) {
		blocksi = ceil((float)size / BLOCKSIZE);
		blocksizes = BLOCKSIZE;
	}
#endif
	cm_increment_cpjds <<<blocksi, blocksizes, 0, stream >>>(size, M.matrixData.data.get(), auxv_d.get(), auxi_d.get());
}

/* method for restoring CPJDS-transformed vector to its original size and indexes */
std::vector<double> MatrixCPJDSManager::restore(Vector * v){
	// groundNode = -1: no ground node
	std::vector<double> restored(nOrig);

	double * data_h = new double[this->n];
	cudaMemcpy(data_h, v->getData(), (size_t)this->n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < nOrig; i++) {
		restored[i] = data_h[original2PaddedIdx[i]];
	}

	return restored;
}

/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
void MatrixCPJDSManager::mult(MatrixCPJDS &M, Vector * x, Vector * b) {
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
	//double * mData = M.matrixData.data;
	//int * mIndices = M.matrixData.indices;
	//int * mRowLength = M.matrixData.rowLength;
	//int * mRowSize = M.matrixData.rowSize;
	//int * mColOffset = M.matrixData.colOffset;

	//int colorCount = M.matrixColors.colorCount;
	//int * colors = M.matrixColors.colors;
	//int * mColorsColOffsetSize = M.matrixColors.colorsColOffsetSize;

	//double * xData = x->getData();
	//double * bData = b->getData();

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

	//// Launch CUDA matrix-vector multiplication kernel
	//cmv_mult <<<blocks, BLOCKSIZE, 0, stream >>>(mData, mIndices, mRowLength, mRowSize,
	//	mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, yData, pData);

	// Launch CUDA matrix-vector multiplication kernel
	cmv_mult_cpjds2 << <blocks, BLOCKSIZE>> >(M.matrixData.n, M.matrixData.data.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get() , M.matrixData.rowSize.get(),
		M.matrixData.colOffset.get(), M.matrixColors.colorCount, M.matrixColors.colors_d.get(), M.matrixColors.colorsColOffsetSize_d.get(), x->getData(), b->getData());

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cout << "mv_mult kernel failed: " << cudaGetErrorString(cudaStatus) << std::endl;
	}
}

/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
void MatrixCPJDSManager::mult(MatrixCPJDS &M, Vector * x, Vector * b, cudaStream_t * streams) {
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
	for (int k = 0; k < colorCount; k++) {

		int colorStart = M.matrixColors.colors[k];
		int colorEnd = M.matrixColors.colors[k + 1];
		int colorColOffset = M.matrixColors.colorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA matrix-vector multiplication kernel
		cmv_mult_cpjds <<<blocks, BLOCKSIZE, 0, streams[k]>>>(size, M.matrixData.data.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
			M.matrixData.colOffset.get(), colorStart, colorColOffset, x->getData(), b->getData());

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
void MatrixCPJDSManager::multInner(MatrixCPJDS &M, Vector * x, Vector * y, cudaStream_t * streams) {
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
	for (int k = 0; k < M.matrixColors.colorCount; k++) {

		int colorStart = M.matrixColors.colors[k];
		int colorEnd = M.matrixColors.colors[k + 1];
		int colorColOffset = M.matrixColors.colorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA matrix-vector multiplication kernel
		cmv_mult_inner_cpjds <<<blocks, BLOCKSIZE, 0, streams[k] >>>(size, M.matrixData.data.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
			M.matrixData.colOffset.get(), colorStart, colorColOffset, x->getData(), y->getData(), y->getPartial());

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
void MatrixCPJDSManager::multInner2(MatrixCPJDS &M, Vector * x, Vector * y, cudaStream_t stream) {
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
	//double * mData = M.preconditionedData;
	//// Launch CUDA matrix-vector multiplication kernel
	//cmv_mult <<<blocks, BLOCKSIZE, 0, stream >>>(mData, mIndices, mRowLength, mRowSize,
	//	mColOffset, colorCount, colors_d, mColorsColOffsetSize_d, xData, yData, pData);

	// Launch CUDA matrix-vector multiplication kernel
	cmv_mult_2 <<<blocks, BLOCKSIZE, 0, stream >>>(M.matrixData.n, M.preconditionedData.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
		M.matrixData.colOffset.get(), M.matrixColors.colorCount, M.matrixColors.colors_d.get(), M.matrixColors.colorsColOffsetSize_d.get(), x->getData(), y->getData(), x->getPartial());

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
void MatrixCPJDSManager::solve(MatrixCPJDS &M, Vector * b, Vector * x) {
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
	//double * mData = M.preconditionedData;//; M.matrixData.data;
	for (int k = 0; k < M.matrixColors.colorCount; k++) {
		int colorStart = M.matrixColors.colors[k];
		int colorEnd = M.matrixColors.colors[k + 1];
		int colorColOffset = M.matrixColors.colorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_cpjds << <blocks, BLOCKSIZE >> >(size, M.preconditionedData.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
			M.matrixData.colOffset.get(), colorStart, colorColOffset, x->getData(), b->getData());

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
void MatrixCPJDSManager::solve(MatrixCPJDS &M, Vector * b, Vector * x, cudaStream_t stream) {
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
	//double * mData = M.preconditionedData; // M.matrixData.data;
	for (int k = 0; k < M.matrixColors.colorCount; k++) {
		int colorStart = M.matrixColors.colors[k];
		int colorEnd = M.matrixColors.colors[k + 1];
		int colorColOffset = M.matrixColors.colorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_cpjds <<<blocks, BLOCKSIZE, 0, stream>>>(size, M.preconditionedData.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
			M.matrixData.colOffset.get(), colorStart, colorColOffset, x->getData(), b->getData());

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
void MatrixCPJDSManager::solve_t(MatrixCPJDS &M, Vector * b, Vector * x) {
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
	//double * mData = M.preconditionedData; //M.matrixData.data;
	for (int k = M.matrixColors.colorCount - 1; k >= 0; k--) {

		int colorStart = M.matrixColors.colors[k];
		int colorEnd = M.matrixColors.colors[k + 1];
		int colorColOffset = M.matrixColors.colorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_t_cpjds << <blocks, BLOCKSIZE >> >(size, M.preconditionedData.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
			M.matrixData.colOffset.get(), colorStart, colorColOffset, x->getData(), b->getData());

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
void MatrixCPJDSManager::solve_t(MatrixCPJDS &M, Vector * b, Vector * x, cudaStream_t stream) {
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
	//double * mData = M.preconditionedData; //M.matrixData.data;
	for (int k = M.matrixColors.colorCount - 1; k >= 0; k--) {

		int colorStart = M.matrixColors.colors[k];
		int colorEnd = M.matrixColors.colors[k + 1];
		int colorColOffset = M.matrixColors.colorsColOffsetSize[k];

		int size = colorEnd - colorStart;

		// Launch CUDA (lower) triangular solver kernel
		cmv_solve_t_cpjds <<<blocks, BLOCKSIZE, 0, stream >>>(size, M.preconditionedData.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
			M.matrixData.colOffset.get(), colorStart, colorColOffset, x->getData(), b->getData());

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
void MatrixCPJDSManager::solve_complete(MatrixCPJDS &M, Vector * b, Vector * x, cudaStream_t stream) {
	// Launch CUDA (lower) triangular solver kernel
	cmv_solve <<<blocks, BLOCKSIZE, 0, stream >>>(M.matrixData.data.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
		M.matrixData.colOffset.get(), M.matrixColors.colorCount, M.matrixColors.colors_d.get(), M.matrixColors.colorsColOffsetSize_d.get(), x->getData(), b->getData());
}

/* multiple operations
* full solver (lower + upper triangulars) M * x = b => x = inv(M) * b
* inner product b.x (x's partials are filled) */
void MatrixCPJDSManager::solve_and_inner(MatrixCPJDS &M, Vector * b, Vector * x, Vector * u, cudaStream_t stream) {
	// Launch CUDA (lower) triangular solver kernel
	cmv_solve_inner <<<blocks, BLOCKSIZE, 0, stream >>>(M.matrixData.data.get(), M.matrixData.indices.get(), M.matrixData.rowLength.get(), M.matrixData.rowSize.get(),
		M.matrixData.colOffset.get(), M.matrixColors.colorCount, M.matrixColors.colors_d.get(), M.matrixColors.colorsColOffsetSize_d.get(), x->getData(), b->getData(), u->getData(), x->getPartial());
}

void MatrixCPJDSManager::saveToFile(char * filename, MatrixCPJDS &M, double * data, bool isCPU) {

	std::ofstream file(filename);
	if (!isCPU) {
		// copy from CUDA device memory to host memory
		double * data_h = new double[M.matrixData.elCount];
		cudaMemcpy(data_h, data, (size_t)M.matrixData.elCount * sizeof(double), cudaMemcpyDeviceToHost);
		data = data_h;
	}

	for (int row = 0; row < M.matrixData.n; row++) {
		for (int col = 0; col < M.matrixData.n; col++) {
			int idx = coordinates2Index(row, col);
			double val = 0;
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
__global__ void cm_set_cpjds(int idx, double * data, double val) {
	// row index
	data[idx] = val;
}

/* CPJDS matrix element incrementer */
__global__ void cm_increment_cpjds(int idx, double * data, double val) {
	// row index
	data[idx] += val;
}

/* CPJDS matrix batch incrementer */
__global__ void cm_increment_cpjds(int size, double * data, double * vals, int * indices) {
	// row index
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < size) {
		data[indices[idx]] = vals[idx];
	}
}

/* CPJDS matrix-vector multiplication */
__global__ void cmv_mult_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * bData) {
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int rowSize = aRowSize[row];

		double sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

			double rowData = aData[offset]; // coalesced
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
__global__ void cmv_mult_cpjds2(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * bData) {
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
		int maxRowSize = aRowLength[0];
		for(int j = 1; j < dim; j++)
			if(aRowLength[j] > maxRowSize)
				maxRowSize = aRowLength[j];
		int rowSize = aRowLength[row];

		double sum = 0;
		__syncthreads();
		for (int j = 0; j < maxRowSize; j++) {
			if(j < rowSize) {
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

				double rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				sum += rowData * xData[idx]; // NOT coalesced!
				//if (row == 5) {
				//	printf("%d\t%d\t%g\n", row, j, xData[idx]);
				//}
				//sum += rowData;
			}

			__syncthreads(); // synchronization so all threads load from memory
		}
		bData[row] = sum; // coalesced
	}

	__syncthreads();
}

/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult_inner_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * yData, double * partial) {
	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
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

		double sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

			double rowData = aData[offset]; // coalesced
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
	//	double sum = 0;
	//	for (int i = 0; i < blocks; i++) {
	//		sum += r[i];
	//	}
	//	val[0] = sum;
	//}
}

// backup do kernel cmv_mult, modificado abaixo
/* CPJDS matrix-vector multiplication, followed by vector inner product */
__global__ void cmv_mult(double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * yData, double * partial) {
	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
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

			double sum = 0;
			__syncthreads();
			for (int j = 0; j < rowSize; j++) {
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

				double rowData = aData[offset]; // coalesced
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
__global__ void cmv_mult_2(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * yData, double * partial) {
	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
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

		double sum = 0;
		__syncthreads();
		for (int j = 0; j < rowSize; j++) {
			// colorColOffset already includes colorOffset (thus, color's first row)
			int offset = aColOffset[colorColOffset + j] + (row - colorStart); // coalesced?

			double rowData = aData[offset]; // coalesced
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
__global__ void cmv_solve_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * bData) {
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int maxRowSize = aRowLength[0];
		for(int j = 1; j < dim; j++)
			if(aRowLength[j] > maxRowSize)
				maxRowSize = aRowLength[j];
		int rowSize = aRowLength[row];

		double sum = 0;
		__syncthreads();
		for (int j = 1; j < maxRowSize; j++) {  // first element is main diagonal
			if(j < rowSize) {
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

				double rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
			}
			__syncthreads();
		}
		xData[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + tidx];
	}

	__syncthreads();
}

/* CPJDS upper triangular solver */
// does not work per se, in parallel
__global__ void cmv_solve_t_cpjds(int dim, double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorOffset, int colorColOffset, double * xData, double * bData) {
	// thread index
	int tidx = blockDim.x * blockIdx.x + threadIdx.x;

	if (tidx < dim) {
		// row index
		int row = tidx + colorOffset;
		// row length
		//int rowLength = aRowLength[row];
		// row size (length + padding zeros)
		int maxRowSize = aRowLength[0];
		for(int j = 1; j < dim; j++)
			if(aRowLength[j] > maxRowSize)
				maxRowSize = aRowLength[j];
		int rowSize = aRowLength[row];

		double sum = 0;
		__syncthreads();
		for (int j = 1; j < maxRowSize; j++) {  // first element is main diagonal
			if(j < rowSize) {
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + tidx; // coalesced?

				double rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx > row && idx > -1) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
			}
			__syncthreads();
		}
		xData[row] = (bData[row] - sum) / aData[aColOffset[colorColOffset] + tidx];
	}

	__syncthreads();
}

/* CPJDS lower and upper triangular solver */
__global__ void cmv_solve(double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * bData) {
	// thread index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	for (int k = 0; k < colorCount; k++) {

		int colorStart = colors[k];
		int colorEnd = colors[k + 1];
		int colorColOffset = colorsColOffset[k];
		int inColorIdx = (row - colorStart);

		if (row >= colorStart && row < colorEnd) {
			int rowSize = aRowSize[row];

			double sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				double rowData = aData[offset]; // coalesced
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

			double sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				double rowData = aData[offset]; // coalesced
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
__global__ void cmv_solve_inner(double * aData, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset,
	int colorCount, int * colors, int * colorsColOffset, double * xData, double * bData, double * interm, double * partial) {
	// shared memory for reduction
	__shared__ double cache[BLOCKSIZE];
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

			double sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				double rowData = aData[offset]; // coalesced
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

			double sum = 0;
			__syncthreads();

			for (int j = 1; j < rowSize; j++) { // first element is main diagonal
				// colorColOffset already includes colorOffset (thus, color's first row)
				int offset = aColOffset[colorColOffset + j] + inColorIdx; // coalesced?

				double rowData = aData[offset]; // coalesced
				int idx = aIndices[offset]; // coalesced
				if (idx < row) { // main diagonal can be skiped
					sum += rowData * xData[idx];
				}
				__syncthreads();
			}
			//double bVal = bData[row];
			double bVal = interm[row]; // result from previous (lower triangular) solver
			double xVal = (bVal - sum) / aData[aColOffset[colorColOffset] + inColorIdx];
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