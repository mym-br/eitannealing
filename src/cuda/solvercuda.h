#ifndef SOLVERCUDA_H_
#define SOLVERCUDA_H_

#include "settings.h"
#include "matrix-cpjds.h"

class PCGSolverCPJDS;
class MatrixCPJDSManager;
struct MatrixCPJDS;
namespace cgl { class Vector; };
enum CGSOLVERTYPE {DEFAULT, CONSOLIDATED, CONSOLIDATEDCG};

class CGCUDA_Solver {
	public:
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, double _LINFinityNorm, double res, CGSOLVERTYPE solverType = DEFAULT, bool init=true);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, double _LINFinityNorm, CGSOLVERTYPE = DEFAULT, bool init=true);
		CGCUDA_Solver(MatrixCPJDS *stiffness, MatrixCPJDSManager *mgr, Vector *bVec, Vector *x0, double _LINFinityNorm, CGSOLVERTYPE = DEFAULT, bool init=true);
		void do_iteration();
		Eigen::Matrix<double, -1, 1, 0> getX();
		int getSize() { return this->size; }
		int getIteration();
		std::tuple<double, double, double> getIterationTimes();
		Vector *getCpjdsX();

		static double createPreconditioner(MatrixCPJDS &M, std::unique_ptr<double[]> &pdata);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness, nodeCoefficients **nodeCoef, int nodesCount, int numcoefficients);
		static MatrixCPJDSManager *createManager(Eigen::SparseMatrix<double> *A, MatrixCPJDS *stiffness);
		static cgl::Vector *createCurrentVector(double *vec, MatrixCPJDSManager &mgr, int size, int n);
		
		double getResidueSquaredNorm() const;
		double getErrorl2Estimate() const;

		// Debug
		static Eigen::SparseMatrix<double, 0, int> getCpjdsStiffness(MatrixCPJDS &M, std::unique_ptr<double[]> &pdata);
		static Eigen::Matrix<double, -1, 1, 0> getCpjdsCurrent(double *vec, MatrixCPJDSManager &mgr, int size, int n);
		void init(double res = -1);

	private:
		static void cblas_dscal(int n, double alpha, double *x, int inc);
		static double m_preconditioner_eigen(MatrixCPJDS &M, std::unique_ptr<double[]> &pdata, std::unique_ptr<double[]> &precond);
		std::vector<double> transfX2Cpu();
		int size;
		PCGSolverCPJDS *solver;
		MatrixCPJDSManager *mgr;
		double LINFinityNorm;

};

#endif // SOLVERCUDA_H_
