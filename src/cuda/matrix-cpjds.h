#ifndef CPJDS_H
#define CPJDS_H

#include "settings.h"

#include <vector>
#include <map>

#include <device_launch_parameters.h>

#include "vector.h"
#include "../nodecoefficients.h"

namespace Eigen {
	template<typename _Scalar, int _Flags = 0, typename _StorageIndex = int>  class SparseMatrix;
}
using namespace cgl;

struct MatrixData {
	int n; // n is multiple of WARP_SIZE (32)
	numType * data;
	int * indices;
	int * rowLength; // row length (non-zeros)
	int * rowSize; // row max (includes trailing padding-zeros) - all rows in the warp have same size
	int * colOffset; // row start (index of the first element of the row)
	int offsetCount; // merely for debug, size of colOffset array
	int nnz; // number of non-zero elements
	int elCount; // total elements (including padding-zeroes)
};

struct MatrixColors {
	int * colors; // color groups/offsets
	int * colorsColOffsetSize; // each color's colOffset block position in the colOffset array
	int colorCount;

	int * colors_d; // color groups/offsets
	int * colorsColOffsetSize_d; // each color's colOffset block position in the colOffset array
};

struct MatrixDependencies {
	int * dependencyRowDataIndex; // element row dependency's index in data array
	int * dependencyDiagDataIndex; // diagonal row dependency's index in data array
	int * dependencyLowerDataIndex; // lower triangular non-zero elements (with dependencies) index in data array
	int * dependencyUpperDataIndex; // upper triangular non-zero elements (with dependencies) index in data array

	int * dependencyArrayInitialIndex; // dependency's array initial index (index in dependencyRowDataIndex and dependencyDiagDataIndex)
	int * dependenciesSize; // dependency's count
	int * nnzElementDataIndex; // data array index for lower and upper triangular's elements (index in dependencyLowerDataIndex and dependencyUpperDataIndex)

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
	numType * data;
	numType * precond;
	int * indices;
};

struct MatrixCPJDS {
	MatrixData matrixData;
	MatrixColors matrixColors;
	MatrixDependencies matrixDependencies;
	numType * preconditionedData;
	
	/* aux map for preconditioner computation */
	CPUData cpuData;
	MatrixCPJDS2CSR csrMap;
};

// column-major!
class MatrixCPJDSManager {
private:
	numType * data;
	int n;
	int * colors;
	int colorCount;

	/* map (row,col)=>data_array_idx */
	std::vector<std::map<int, int>> coord2IndexMap;
	/* original matrix size, before padding */
	int nOrig;

	/* aux vector for transforming, incrementing, etc */
	numType * auxv;
	/* aux indices vector for transforming, incrementing, etc */
	int * auxi;
	/* maximum coefficient dependencies' count */
	int maxDepCount;
	/* aux device array for incrementing */
	numType * auxv_d;
	/* aux device array for incrementing */
	int * auxi_d;

	/* cached number of blocks so that it does not need to be computed on demand */
	int blocks;

public:
	/* map (original index)=>(color-sorted-padded index) [size N] */
	int * original2PaddedIdx;
	/* map (color-sorted-padded index)=>(original index) [size N padded] */
	int * padded2OriginalIdx;
	// data must be processed and color-sorted
	MatrixCPJDSManager(numType * data, int n);
	MatrixCPJDSManager(Eigen::SparseMatrix<double> *data);
	// data has been pre-processed and color-sorted, n already includes padding
	MatrixCPJDSManager(numType * data, int n, int * colors, int colorCount, int nOrig);
	~MatrixCPJDSManager();

	/* provided M matrix is filled with a complete CPJDS matrix */
	int buidMatrixCPJDS(MatrixCPJDS * M, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
	/* provided L and U matrices are filled with corresponding triangular CPJDS matrices */
	int buidMatrixCPJDS(MatrixCPJDS * L, MatrixCPJDS * U);

	/* clear memory */
	void deleteMatrixCPJDS(MatrixCPJDS M);

	/* this method allows the conversion os (row, col) coordinates to the index in the data and indices arrays */
	int coordinates2Index(int row, int col);

	/* sets an element value according to row and column indexes */
	void set(MatrixCPJDS M, int row, int col, numType val);
	
	/* increments an element value according to row and column indexes */
	void increment(MatrixCPJDS M, int row, int col, numType val);
	
	/* increments an array of elements value according to elements' indexes */
	void pushIncrements(MatrixCPJDS M, int size, numType * vals, int * indices);
	/* increments an array of elements value according to elements' indexes */
	void pushIncrements(MatrixCPJDS M, int size, numType * vals, int * indices, cudaStream_t stream);

	/* method for rearranging a regular vector to follow the CPJDS transformations */
	Vector * transform(std::vector<numType> v, bool removeGround);
	
	/* method for restoring CPJDS-transformed vector to its original size and indexes */
	std::vector<numType> restore(Vector * v);
	
	/* method for creating an "electrode" mask vector */
	Vector * mask();

	/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
	void mult(MatrixCPJDS M, Vector * x, Vector * b);
	/* Matrix-vector multiplication: b = A * x, A: matrix; b, x: vectors */
	void mult(MatrixCPJDS M, Vector * x, Vector * b, cudaStream_t * streams);

	/* Matrix-vector multiplication (stored) followed by an inner produtct: k = xt * M * x
	 * (totalization among blocks is NOT performed!)
	 * where y = M * x, M: matrix; x, y: vectors, k: number */
	void multInner(MatrixCPJDS M, Vector * x, Vector * y, cudaStream_t * streams);
	/* Matrix-vector multiplication (stored) followed by an inner produtct: k = xt * M * x
	* (totalization among blocks is NOT performed!)
	* where y = M * x, M: matrix; x, y: vectors, k: number
	* (uses dynamic parallelism instead of multiple streams) */
	void multInner2(MatrixCPJDS M, Vector * x, Vector * y, cudaStream_t stream);

	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// L must be lower triangular
	void solve(MatrixCPJDS M, Vector * b, Vector * x);
	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// L must be lower triangular
	void solve(MatrixCPJDS M, Vector * b, Vector * x, cudaStream_t stream);

	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// U must be an upper triangular matrix
	void solve_t(MatrixCPJDS M, Vector * b, Vector * x);
	/* Linear system solver: A * x = b => x = inv(A) * b; A: matrix; x, b: vectors */
	// U must be an upper triangular matrix
	void solve_t(MatrixCPJDS M, Vector * b, Vector * x, cudaStream_t stream);

	/* full solver (lower + upper triangulars) M * x = b => x = inv(M) * b */
	void solve_complete(MatrixCPJDS M, Vector * b, Vector * x, cudaStream_t stream);

	/* multiple operations
	 * full solver (lower + upper triangulars) M * x = b => x = inv(M) * b
	 * inner product b.x (x's partials are filled) */
	void solve_and_inner(MatrixCPJDS M, Vector * b, Vector * x, Vector * u, cudaStream_t stream);

	void saveToFile(char * filename, MatrixCPJDS M, numType * data, bool isCPU);
};

#endif /* CPJDS_H */