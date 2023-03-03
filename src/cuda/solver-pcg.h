#ifndef SOLVERPCG_H
#define SOLVERPCG_H

#include <device_launch_parameters.h>

#include "matrix-cpjds.h"
#include "vector.h"
#include "number.h"
#include "../circularbuff.h"

using namespace cgl;

class PCGSolverCPJDS {
protected:
	int it;
	Number * rmod, * rmod_prev, * rmod_aux;
	Number *gamma;

	MatrixCPJDSManager * mgr;
	MatrixCPJDS *A;
	Vector * b;
	Vector * x;

	Vector * r;
	Vector * z;
	Vector * p;
	Vector * q;
	Vector * u;
	Vector * partial, *partial2;

	/* array of cuda streams for parallel queueing */

	int size, blocks;
	double totalItTime, totalTriangularTime, totalSpmvTime;

public:
	cudaStream_t stream;

	int getIteration() const {
		return this->it;
	}
	std::tuple<double, double, double> getAvgTimes() {
		return{ totalItTime/(double)it, totalTriangularTime/(double)it, totalSpmvTime/(double)it };
	}

	PCGSolverCPJDS(MatrixCPJDSManager * mgr, MatrixCPJDS *M, Vector * b);
	~PCGSolverCPJDS() {
		delete x;
		delete r;
		delete z;
		delete p;
		delete q;
		delete u;
		delete partial;

		delete rmod, rmod_prev, rmod_aux;
		delete gamma;

		streamDestroy();
	}

	void init(double res = -1);
	virtual void init(Vector *x0, double res = -1);
	virtual  void doIteration(int iteration = -1);

	virtual Vector * getX() {
		return this->x;
	}

	void streamInit();
	void streamDestroy();

	virtual double getRmod() { 
		#ifndef CALCULATE_ERRORS
		return rmod->transf2CPU();
		#else
		return rmod2;
		#endif
	}
	double getR0norm() { return r0norm; }
	double getCurrentErr() {
		if (it>2) return this->err[it - 1];
		return 0;
	}
protected:
	void checkedCudaEventRecord(cudaEvent_t &event);
	double rmod2, rmod2_1, gamma2, gamma2_1, beta, r0norm2, r0norm, alpha, eta, eta_p1, rt1, r1, c, s, r2, c_1, s_1, r3;
	// Circular buffers
	circularbuff<double, 8> w;
	circularbuff<double, 8> wt;
	circularbuff<double, 8> err;
};

class PCGSolverConsolidatedCPJDS : public PCGSolverCPJDS {
public:
	PCGSolverConsolidatedCPJDS(MatrixCPJDSManager * mgr, MatrixCPJDS *M, Vector * b) : PCGSolverCPJDS(mgr, M, b) { this->x_1 = new Vector(M->matrixData.n); }

	void init(Vector *x0, double res = -1);
	void doIteration(int iteration = -1);
	void doIteration0(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	void doIteration1(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	void doIteration2(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	void doIteration3(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	double getRmod() { 
		#ifndef CALCULATE_ERRORS
		return rmod->transf2CPU(); 
		#else
		return rmod2;
		#endif
	}
	Vector * getX() {
		return this->x_1;
	}
private:
	Vector * x_1;
};

#ifdef CGROUPS
class PCGSolverConsolidatedCPJDSCG : public PCGSolverCPJDS {
public:
	PCGSolverConsolidatedCPJDSCG(MatrixCPJDSManager * mgr, MatrixCPJDS *M, Vector * b) : PCGSolverCPJDS(mgr, M, b) { this->x_1 = new Vector(M->matrixData.n); }

	void init(Vector *x0, double res = -1);
	void doIteration(int iteration = -1);
	void doIteration0(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	void doIteration1(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	void doIteration2(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	void doIteration3(double * aData, double * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, double * zData, double * rData, double * xData, double * pData, double * qData, double * partialData, double * partialData2);
	double getRmod() {
#ifndef CALCULATE_ERRORS
		return rmod->transf2CPU();
#else
		return rmod2;
#endif
	}
	Vector * getX() {
		return this->x_1;
	}
private:
	Vector * x_1;
};
#endif

#endif /* SOLVERPCG_H */