#include "solvercuda.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

void CGCUDA_Solver::cblas_dscal(int n, numType alpha, numType *x, int inc) {
	for (int i = 0; i < n; i++) {
		*x *= alpha;
		x += inc;
	}
}

typedef Eigen::SparseMatrix<numType, Eigen::ColMajor> matrix;
typedef Eigen::SparseMatrix<numType, Eigen::ColMajor> matrixLower;
bool CGCUDA_Solver::m_preconditioner_eigen(MatrixCPJDS M, numType * pdata, numType * precond) {

	int size = M.csrMap.n;
	int nnz = M.csrMap.nnz;
	MatrixCPJDS2CSR csrMap = M.csrMap;

	//matrix stiff(size, size);
	//stiff.startFill(nnz); // estimate of the number of nonzeros (optional)
	std::vector<Eigen::Triplet<numType>> tripletList;
	for (int i = 0; i < nnz; ++i) {
		int row = M.csrMap.row[i];
		int dataIdx = csrMap.csr2cpjds[i];

		numType val = pdata[dataIdx];
		int col = M.csrMap.indices[i];

		//stiff.fill(row, col) = val;
		tripletList.push_back(Eigen::Triplet<numType>(row, col, val));
	}
	//stiff.endFill();
	matrix stiff(size, size);
	stiff.setFromTriplets(tripletList.begin(), tripletList.end());
	stiff.makeCompressed();

	matrixLower chol_mat(stiff);
	numType *data = chol_mat.valuePtr();
	for (int col = 0; col < size; col++) {
		matrixLower::InnerIterator it(chol_mat, col);
		int *outer = chol_mat.outerIndexPtr();

		if (it.value() < 0) {
			assert(false);
			return false;
		}
		numType isqrtDiagonal = 1 / std::sqrt(it.value());
		// Multiply the whole column
		cblas_dscal(outer[col + 1] - outer[col], isqrtDiagonal, data + outer[col], 1);
		// This is not unlike a sparse vector-vector multiplication
		while (++it) {
			matrixLower::InnerIterator source(it);
			matrixLower::InnerIterator target(chol_mat, it.row());
			while (target && source) {
				// Sweep and subtract on coincident rows
				//	This should be relatively quick, as both target and source have very few
				//		non-zero entries
				if (target.row() == source.row()) {
					target.valueRef() -= source.value()*it.value();
					++target; ++source;
				}
				while (target && target.row() < source.row()) {
					++target;
				}
				while (source && source.row() < target.row()) {
					++source;
				}
			}
		}
	}

	// fill precond
	for (int i = 0; i < nnz; ++i) {
		int dataIdx = csrMap.csr2cpjds[i];
		int dataIdxUper = csrMap.csr2cpjds_upper[i];
		precond[dataIdx] = data[i];
		precond[dataIdxUper] = data[i];
	}

	return true;
}