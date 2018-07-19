#include "solvercuda.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "../incomplete_cholesky.h"
#include "../basematrix.h"

vectorx CGCUDA_Solver::getCpjdsCurrent(numType *vec, MatrixCPJDSManager &mgr, int size, int n) {
	vectorx vecArr(size);
	for (int i = 0; i < size; i++) {
		vecArr[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		vecArr[mgr.original2PaddedIdx[i]] = vec[i];
	}
	return vecArr;
}

matrix CGCUDA_Solver::getCpjdsStiffness(MatrixCPJDS M, numType * pdata) {
	int size = M.csrMap.n;
	int nnz = M.csrMap.nnz;
	MatrixCPJDS2CSR csrMap = M.csrMap;

	std::vector<Eigen::Triplet<double>> tripletList;
	for (int i = 0; i < nnz; ++i) {
		int row = M.csrMap.row[i];
		int dataIdx = csrMap.csr2cpjds[i];

		double val = pdata[dataIdx];
		int col = M.csrMap.indices[i];

		tripletList.push_back(Eigen::Triplet<double>(row, col, val));
	}
	matrix stiff(size, size);
	stiff.setFromTriplets(tripletList.begin(), tripletList.end());
	stiff.makeCompressed();
	return stiff;
}

void CGCUDA_Solver::cblas_dscal(int n, numType alpha, numType *x, int inc) {
	for (int i = 0; i < n; i++) {
		*x *= alpha;
		x += inc;
	}
}

numType CGCUDA_Solver::m_preconditioner_eigen(MatrixCPJDS &M, std::shared_ptr<numType> pdata, std::shared_ptr<numType> precond) {

	int size = M.csrMap.n;
	int nnz = M.csrMap.nnz;
	MatrixCPJDS2CSR csrMap = M.csrMap;

	std::vector<Eigen::Triplet<numType>> tripletList;
	for (int i = 0; i < nnz; ++i) {
		int row = M.csrMap.row[i];
		int dataIdx = csrMap.csr2cpjds[i];

		numType val = pdata.get()[dataIdx];
		int col = M.csrMap.indices[i];

		tripletList.push_back(Eigen::Triplet<numType>(row, col, val));
	}
	Eigen::SparseMatrix<double> stiff(size, size);
	stiff.setFromTriplets(tripletList.begin(), tripletList.end());
	stiff.makeCompressed();

	SparseIncompleteLLT chol_mat(stiff);
	const Scalar *data = chol_mat.matrixL().valuePtr();

	// fill precond
	for (int i = 0; i < nnz; ++i) {
		int dataIdx = csrMap.csr2cpjds[i];
		int dataIdxUper = csrMap.csr2cpjds_upper[i];
		precond.get()[dataIdx] = data[i];
		precond.get()[dataIdxUper] = data[i];
	}

	return chol_mat.getLINFinityNorm();
}


Eigen::VectorXf CGCUDA_Solver::getX() {
	std::vector<numType> xCpu = this->transfX2Cpu();
	if(std::isnan(xCpu[0])) throw std::exception();
	float *data = &(*xCpu.begin());
	Eigen::Map<Eigen::VectorXf> ans(data, xCpu.size());
	return ans;
}
