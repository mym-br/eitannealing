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
	Vector * partial;

	/* array of cuda streams for parallel queueing */

	int size, blocks;

public:
	cudaStream_t stream;

	int getIteration() const {
		return this->it;
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

	void init();
	virtual void init(Vector *x0);
	virtual  void doIteration(int iteration = -1);
	//void doIteration0(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
	//void doIteration1(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
	//void doIteration2(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
	//void doIteration3(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);

	Vector * getX() {
		return this->x;
	}

	void streamInit();
	void streamDestroy();

	double getRmod() { return rmod_prev->transf2CPU(); }
	double getR0norm() { return r0norm; }
	double getCurrentErr() {
		if (it>2) return this->err[it - 1];
		return 0;
	}
private:
	double rmod2, rmod2_1, gamma2, gamma2_1, beta, r0norm2, r0norm, alpha, eta, eta_p1, rt1, r1, c, s, r2, c_1, s_1, r3;
	// Circular buffers
	circularbuff<double, 8> w;
	circularbuff<double, 8> wt;
	circularbuff<double, 8> err;
};

class PCGSolverCPJDS2 : public PCGSolverCPJDS {
public:
	PCGSolverCPJDS2(MatrixCPJDSManager * mgr, MatrixCPJDS *M, Vector * b) : PCGSolverCPJDS(mgr, M, b) {}

	void init(Vector *x0);
	void doIteration(int iteration = -1);
	void doIteration0(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
	void doIteration1(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
	void doIteration2(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
	void doIteration3(numType * aData, numType * precond, int * aIndices, int * aRowLength, int * aRowSize, int * aColOffset, int colorCount, int * colors, int * colorsColOffset, numType * zData, numType * rData, numType * xData, numType * pData, numType * qData, numType * partialData);
};
#endif /* SOLVERPCG_H */