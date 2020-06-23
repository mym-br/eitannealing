#ifndef DIRECTSOLVER_H
#define DIRECTSOLVER_H
#include <memory>
#include "observations.h"

class problem;
class SparseIncompleteLLT;

class EitDirectSolver {
public:
	EitDirectSolver(const char* meshfilename, const  char* currentfilename);
	EitDirectSolver(const EitDirectSolver &other);
	void loadMesh(const char* meshfilename);
	int getCoeffCount();
	int getCurrentPatternCount();
	int getElectrodesCount();
	int getGroundNode();
	void setconds(double* cond, int n);
	double* solve(int patterno);

private:
	void initialize(const char* meshfilename, const  char* currentfilename);

	std::shared_ptr<problem> input;
	std::unique_ptr<observations<double>> readings;
	matrix* m1;
	std::shared_ptr<SparseIncompleteLLT> precond;
	const char* m_meshfilename, * m_currentfilename;
};

#endif // DIRECTSOLVER_H