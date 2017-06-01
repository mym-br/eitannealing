/*
* problem.h
*
*  Created on: Feb 23, 2017
*      Author: aksato
*/

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include <vector>
#include <map>
#include <fstream>
#include <functional>
#include <memory>
#include <Eigen/Core>
#include "basematrix.h"
class solution;

struct nodeCoefficients {
	int node;
	int condIndex;
	double coefficient;
	nodeCoefficients *next;

	nodeCoefficients(int node, int index, double coefficient) :
		node(node), condIndex(index), coefficient(coefficient), next(NULL) {}

	nodeCoefficients(int node, int index, double coefficient, nodeCoefficients *next) :
		node(node), condIndex(index), coefficient(coefficient), next(next) {}
};

class problem {
	friend class gradientNormRegularisation;
	friend class gradientNormRegularisation_old;
	friend class gradientNormRegularisationComplex;
	friend class solution;
	friend class solutioncomplex;

	private:
	int groundNode;

	protected:
	// Node data
	float electrodeh;
	std::map<int, int> node2coefficient;
	int numcoefficients;
	int nodeCount;
	// Observations
	int nobs;
	Eigen::VectorXd *tensions;
	Eigen::VectorXd *rebased_tensions;
	Eigen::VectorXd *currents;
	Eigen::VectorXd currentVals;
	nodeCoefficients **nodeCoef;
	std::vector<std::pair<int, int> > innerAdjacency;
	matrix *skeleton;
	matrix *coef2KMatrix;
	const char *filename;

	Eigen::VectorXcd *tensionscomplex;
	Eigen::VectorXcd *rebased_tensionscomplex;
	Eigen::VectorXcd *currentscomplex;
	Eigen::VectorXcd currentValscomplex;

public:
	static std::shared_ptr<problem> createNewProblem(const char *meshfilename, bool &is2D);
	virtual void initProblem(const char *meshfilename) = 0;
	void initObs(const char *filecurrents, const char* filename);
	void initObsComplex(const char *filecurrents, const char* filename);
	virtual void buildNodeCoefficients() = 0;
	virtual int getGenericElectrodesCount() = 0;
	virtual int getNodesCount() { return nodeCount; }
	virtual int getInnerAdjacencyCount() = 0;
	int getNumCoefficients() { return numcoefficients; }
	nodeCoefficients **getNodeCoefficients() { return nodeCoef; }
	int getNode2Coefficient(int id) { return node2coefficient[id]; }
	problem(const char *meshfilename) : filename(meshfilename), groundNode(-1), nobs(-1), 
		tensions(nullptr), rebased_tensions(nullptr), currents(nullptr), 
		tensionscomplex(nullptr), rebased_tensionscomplex(nullptr), currentscomplex(nullptr), 
		skeleton(nullptr), coef2KMatrix(nullptr), nodeCoef(nullptr) {};
	virtual ~problem();
	int getNObs() { return nobs; }
	Eigen::VectorXd *getTensions() { return tensions; }
	Eigen::VectorXd *getCurrents() { return currents; }
	Eigen::VectorXcd *getTensionsComplex() { return tensionscomplex; }
	Eigen::VectorXcd getCurrentsComplex(matrixcomplex *m, int i) {
		Eigen::VectorXcd vec = currentscomplex[i];
		vec = (*m).conjugate().selfadjointView<Eigen::Lower>() * currentscomplex[i];
		#ifndef BLOCKGND
		vec[getGroundNode()] = 0.0;
		#endif
		return vec;
	}
	std::complex<double> getCurrentValComplex(int i) { return currentValscomplex[i]; }
	int getCurrentsComplexCount() { return (int)currentValscomplex.size(); }
	void prepareSkeletonMatrix();
	void createCoef2KMatrix();
	void assembleProblemMatrix(double *cond, matrix **stiffnes);
	void postAssempleProblemMatrix(matrix **stiffnes);
	void assembleProblemMatrix(std::complex<double> *cond, matrixcomplex **stiffnes);
	void postAssempleProblemMatrix(matrixcomplex **stiffnes);
	const char* getMeshFilename() { return filename; }
	double getCurrentVal(int i) { return currentVals[i]; }
	int getCurrentsCount() { return (int)currentVals.size(); }
	void setGroundNode(int nodeid);
	int getGroundNode() { return this->groundNode; }
};

const double mincond = 0.005;
const double maxcond = 0.3815;

#endif // PROBLEM_H_