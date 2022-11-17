#ifndef CPJDS_H
#define CPJDS_H

#include "settings.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include <deque>
#include <map>
#include<memory>

#include <device_launch_parameters.h>

#include "vector.h"
#include "../nodecoefficients.h"

using namespace cgl;

struct DeleterCudaIntPtr { void operator()(int* ptr); void operator()(numType* ptr); };

struct MatrixData {
	int n; // n is multiple of WARP_SIZE (32)
	// TODO: convert ALL to smart cuda ptr
	std::unique_ptr<numType[], DeleterCudaIntPtr> data;
	std::unique_ptr<int[], DeleterCudaIntPtr> indices;
	std::unique_ptr<int[], DeleterCudaIntPtr> rowLength; // row length (non-zeros)
	std::unique_ptr<int[], DeleterCudaIntPtr> rowSize; // row max (includes trailing padding-zeros) - all rows in the warp have same size
	std::unique_ptr<int[], DeleterCudaIntPtr> colOffset; // row start (index of the first element of the row)
	int offsetCount; // merely for debug, size of colOffset array
	int nnz; // number of non-zero elements
	int elCount; // total elements (including padding-zeroes)
};

struct MatrixColors {
	std::shared_ptr<int[]> colors; // color groups/offsets
	std::unique_ptr<int[]> colorsColOffsetSize; // each color's colOffset block position in the colOffset array
	int colorCount;

	std::unique_ptr<int[], DeleterCudaIntPtr> colors_d; // color groups/offsets
	std::unique_ptr<int[], DeleterCudaIntPtr> colorsColOffsetSize_d; // each color's colOffset block position in the colOffset array
};

struct MatrixDependencies {
	// TODO: convert ALL to smart cuda ptr
	std::unique_ptr<int[], DeleterCudaIntPtr> dependencyRowDataIndex; // element row dependency's index in data array
	std::unique_ptr<int[], DeleterCudaIntPtr> dependencyDiagDataIndex; // diagonal row dependency's index in data array
	std::unique_ptr<int[], DeleterCudaIntPtr> dependencyLowerDataIndex; // lower triangular non-zero elements (with dependencies) index in data array
	std::unique_ptr<int[], DeleterCudaIntPtr> dependencyUpperDataIndex; // upper triangular non-zero elements (with dependencies) index in data array

	std::unique_ptr<int[], DeleterCudaIntPtr> dependencyArrayInitialIndex; // dependency's array initial index (index in dependencyRowDataIndex and dependencyDiagDataIndex)
	std::unique_ptr<int[], DeleterCudaIntPtr> dependenciesSize; // dependency's count
	std::unique_ptr<int[], DeleterCudaIntPtr> nnzElementDataIndex; // data array index for lower and upper triangular's elements (index in dependencyLowerDataIndex and dependencyUpperDataIndex)

	int depsSize; // size of dependencyRowDataIndex and dependencyDiagDataIndex arrays (just for debugging)
	int idxSize; // size of dependencyLowerDataIndex and dependencyUpperDataIndex arrays (just for debugging)
};

struct MatrixCPJDS2CSR {
	int n;
	int nnz;// this is actually element count in both triangular matrices, not
	std::vector<int> csr2cpjds;
	std::vector<int> csr2cpjds_upper;
	std::vector<int> indices;
	std::vector<int> row;
};

struct CPUData {
	std::unique_ptr<numType[]> data;
	std::unique_ptr<numType[]> precond;
	std::unique_ptr<int[]> indices;
};

struct MatrixCPJDS {
	MatrixData matrixData;
	MatrixColors matrixColors;
	MatrixDependencies matrixDependencies;
	std::unique_ptr<numType[], DeleterCudaIntPtr> preconditionedData; // TODO: convert to smart cuda ptr

	/* aux map for preconditioner computation */
	CPUData cpuData;
	MatrixCPJDS2CSR csrMap;
};

struct Dependencies {
	std::vector<int> dependencies;
	int dependenciesSize;
	std::pair <int, int> lower; // upper can be calculated from lower
};

typedef std::vector<std::vector<Dependencies>> DependeciesMap;

// column-major!
class MatrixCPJDSManager {
private:
	std::unique_ptr<Eigen::SparseMatrix<double>> data2;

	int n;
	std::shared_ptr<int[]> colors;
	int colorCount;

	/* map (row,col)=>data_array_idx */
	std::vector<std::unique_ptr<std::map<int, int>>> coord2IndexMap;
	/* original matrix size, before padding */
	int nOrig;

	/* aux vector for transforming, incrementing, etc */
	std::unique_ptr<numType[]> auxv;
	/* aux indices vector for transforming, incrementing, etc */
	std::unique_ptr<int[]> auxi;
	/* maximum coefficient dependencies' count */
	int maxDepCount;
	/* aux device array for incrementing */
	std::unique_ptr<numType[], DeleterCudaIntPtr> auxv_d; // TODO: convert to smart cuda ptr
	/* aux device array for incrementing */
	std::unique_ptr<int[], DeleterCudaIntPtr> auxi_d; // TODO: convert to smart cuda ptr

	/* cached number of blocks so that it does not need to be computed on demand */
	int blocks;

	void leftShiftMatrix(std::vector<std::deque<int>> &rowsL, std::vector<std::deque<int>> &rowsU);
	void createDataAndIndicesVectors(numType *mdata, int *indices, int *colOffset, int *colorOffsetCount, std::vector<std::deque<int>> &rowsL, std::vector<std::deque<int>> &rowsU, std::vector<std::deque<int>> &padding);
	void dependencies_analysis2(int n, DependeciesMap * dependenciesMap);
	void createCsr2CpjdsMap(MatrixCPJDS2CSR &csrMap);
public:
	/* map (original index)=>(color-sorted-padded index) [size N] */
	std::unique_ptr<int[]> original2PaddedIdx;
	/* map (color-sorted-padded index)=>(original index) [size N padded] */
	std::unique_ptr<int[]> padded2OriginalIdx;
	// data must be processed and color-sorted
	MatrixCPJDSManager(Eigen::SparseMatrix<double> *data);

	/* provided M matrix is filled with a complete CPJDS matrix */
	int buidMatrixCPJDS(MatrixCPJDS * M, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
	int buidMatrixCPJDS(MatrixCPJDS * M);

	/* this method allows the conversion os (row, col) coordinates to the index in the data and indices arrays */
	int coordinates2Index(int row, int col);

	/* sets an element value according to row and column indexes */
	void set(MatrixCPJDS &M, int row, int col, numType val);

	/* increments an element value according to row and column indexes */
	void increment(MatrixCPJDS &M, int row, int col, numType val);

	/* increments an array of elements value according to elements' indexes */
	void pushIncrements(MatrixCPJDS &M, int size, numType * vals, int * indices);
	/* increments an array of elements value according to elements' indexes */
	void pushIncrements(MatrixCPJDS &M, int size, numType * vals, int * indices, cudaStream_t stream);

	/* method for rearranging a regular vector to follow the CPJDS transformations */
	Vector * transform(std::vector<numType> v, bool removeGround);

	/* method for restoring CPJDS-transformed vector to its original size and indexes */
	std::vector<numType> restore(Vector * v);

	/* method for creating an "electrode" mask vector */
	Vector * mask();

	/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
	void mult(MatrixCPJDS &M, Vector * x, Vector * b);
	/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
	void mult(MatrixCPJDS &M, Vector * x, Vector * b, cudaStream_t * streams);

	/* Matrix-vector multiplication (stored) followed by an inner produtct: k = xt * M * x
	 * (totalization among blocks is NOT performed!)
	 * where y = M * x, M: matrix; x, y: vectors, k: number */
	void multInner(MatrixCPJDS &M, Vector * x, Vector * y, cudaStream_t * streams);
	/* Matrix-vector multiplication (stored) followed by an inner produtct: k = xt * M * x
	* (totalization among blocks is NOT performed!)
	* where y = M * x, M: matrix; x, y: vectors, k: number
	* (uses dynamic parallelism instead of multiple streams) */
	void multInner2(MatrixCPJDS &M, Vector * x, Vector * y, cudaStream_t stream);

	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// L must be lower triangular
	void solve(MatrixCPJDS &M, Vector * b, Vector * x);
	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// L must be lower triangular
	void solve(MatrixCPJDS &M, Vector * b, Vector * x, cudaStream_t stream);

	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// U must be an upper triangular matrix
	void solve_t(MatrixCPJDS &M, Vector * b, Vector * x);
	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// U must be an upper triangular matrix
	void solve_t(MatrixCPJDS &M, Vector * b, Vector * x, cudaStream_t stream);

	/* full solver (lower + upper triangulars) M * x = b => x = inv(M) * b */
	void solve_complete(MatrixCPJDS &M, Vector * b, Vector * x, cudaStream_t stream);

	/* multiple operations
	 * full solver (lower + upper triangulars) M * x = b => x = inv(M) * b
	 * inner product b.x (x's partials are filled) */
	void solve_and_inner(MatrixCPJDS &M, Vector * b, Vector * x, Vector * u, cudaStream_t stream);

	void saveToFile(char * filename, MatrixCPJDS &M, numType * data, bool isCPU);
};

#endif /* CPJDS_H */
