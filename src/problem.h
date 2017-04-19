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
	friend class solution;

	protected:
	// Node data
	float electrodeh;
	std::map<int, int> node2coefficient;
	int numcoefficients;
	int groundNode;
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
	
public:
	static std::shared_ptr<problem> createNewProblem(const char *meshfilename, bool &is2D);
	virtual void initProblem(const char *meshfilename) = 0;
	void initObs(const char *filecurrents, const char* filename);
	virtual void buildNodeCoefficients() = 0;
	virtual int getGenericElectrodesCount() = 0;
	virtual int getNodesCount() = 0;
	virtual int getInnerAdjacencyCount() = 0;
	int getNumCoefficients() { return numcoefficients; }
	nodeCoefficients **getNodeCoefficients() { return nodeCoef; }
	int getNode2Coefficient(int id) { return node2coefficient[id]; }
	problem(const char *meshfilename) : filename(meshfilename), groundNode(-1) {};
	virtual ~problem(){};
	int getNObs() { return nobs; }
	Eigen::VectorXd *getTensions() { return tensions; }
	Eigen::VectorXd *getCurrents() { return currents; }
	void prepareSkeletonMatrix();
	void createCoef2KMatrix();
	void assembleProblemMatrix(double *cond, matrix **stiffnes);
	const char* getMeshFilename() { return filename; }
	double getCurrentVal(int i) { return currentVals[i]; }
	int getCurrentsCount() { return (int)currentVals.size(); }
	void setGroundNode(int nodeid) { this->groundNode = nodeid; }
	int getGroundNode() { return this->groundNode; }
};

const double mincond = 0.005;
const double maxcond = 0.3815;

#endif // PROBLEM_H_