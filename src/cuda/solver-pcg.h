#ifndef SOLVERPCG_H
#define SOLVERPCG_H

#include <device_launch_parameters.h>

#include "matrix-cpjds.h"
#include "vector.h"
#include "number.h"

using namespace cgl;

class PCGSolverCPJDS {
protected:
	int it;
	//Number * rmod, * rmod_prev, * rmod_aux;

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

	//Number * alpha, *beta, *gamma;

	/* array of cuda streams for parallel queueing */

	int size, blocks;

public:

	Number * rmod, *rmod_prev, *rmod_aux;
	Number * alpha, *beta, *gamma;

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
		delete alpha, beta, gamma;

		streamDestroy();
	}

	void init();
	void doIteration(int iteration = -1);

	Vector * getX() {
		return this->x;
	}

	void streamInit();
	void streamDestroy();
};

#endif /* SOLVERPCG_H */