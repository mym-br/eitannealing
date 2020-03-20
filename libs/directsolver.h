#ifndef DIRECTSOLVER_H
#define DIRECTSOLVER_H
#include <memory>
#include "observations.h"

class problem;
class SparseIncompleteLLT;

class EitDirectSolver {
public:
	EitDirectSolver(const char* meshfilename, const  char* currentfilename);
	int getCoeffCount();
	int getCurrentPatternCount();
	int getElectrodesCount();
	void setconds(double* cond, int n);
	double *solve(int patterno);

private:
	std::shared_ptr<problem> input;
	std::unique_ptr<observations<double>> readings;
	matrix *m1;
	std::shared_ptr<SparseIncompleteLLT> precond;
};

#endif // DIRECTSOLVER_H
