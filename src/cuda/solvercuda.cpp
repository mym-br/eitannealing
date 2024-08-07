#include "solvercuda.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "../incomplete_cholesky.h"
#include "../basematrix.h"

vectorx CGCUDA_Solver::getCpjdsCurrent(double *vec, MatrixCPJDSManager &mgr, int size, int n) {
	vectorx vecArr(size);
	for (int i = 0; i < size; i++) {
		vecArr[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		vecArr[mgr.original2PaddedIdx[i]] = vec[i];
	}
	return vecArr;
}

matrix CGCUDA_Solver::getCpjdsStiffness(MatrixCPJDS &M, std::unique_ptr<double[]> &pdata) {
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

void CGCUDA_Solver::cblas_dscal(int n, double alpha, double *x, int inc) {
	for (int i = 0; i < n; i++) {
		*x *= alpha;
		x += inc;
	}
}

double CGCUDA_Solver::m_preconditioner_eigen(MatrixCPJDS &M, std::unique_ptr<double[]> &pdata, std::unique_ptr<double[]> &precond) {

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
	Eigen::SparseMatrix<double> stiff(size, size);
	stiff.setFromTriplets(tripletList.begin(), tripletList.end());
	stiff.makeCompressed();

	SparseIncompleteLLT chol_mat(stiff);
	const Scalar *data = chol_mat.matrixL().valuePtr();

	// fill precond
	for (int i = 0; i < nnz; ++i) {
		int dataIdx = csrMap.csr2cpjds[i];
		int dataIdxUper = csrMap.csr2cpjds_upper[i];
		precond[dataIdx] = data[i];
		precond[dataIdxUper] = data[i];
	}

	return chol_mat.getLINFinityNorm();
}


Eigen::Matrix<double, -1, 1, 0> CGCUDA_Solver::getX() {
	std::vector<double> xCpu = this->transfX2Cpu();
	if(std::isnan(xCpu[0])) throw std::exception();
	double *data = &(*xCpu.begin());
	Eigen::Map<Eigen::Matrix<double, -1, 1, 0>> ans(data, xCpu.size());
	return ans;
}
