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
	nodeCoefficients **nodeCoef;
	std::vector<std::pair<int, int> > innerAdjacency;
	
public:
	virtual void initProblem(char *meshfilename) = 0;
	void initObs(char *filecurrents, char* filename);
	virtual void buildNodeCoefficients() = 0;
	virtual int getGenericElectrodesCount() = 0;
	virtual int getNodesCount() = 0;
	virtual int getInnerAdjacencyCount() = 0;
	int getNumCoefficients() { return numcoefficients; }
	virtual void prepareSkeletonMatrix() = 0;
	virtual void createCoef2KMatrix() = 0;
	virtual void assembleProblemMatrix(double *cond, matrix **stiffnes) = 0;
	nodeCoefficients **getNodeCoefficients() { return nodeCoef; }
	problem() {};
	virtual ~problem(){};
	int getNObs() { return nobs; }
	Eigen::VectorXd *getTensions() { return tensions; }
	Eigen::VectorXd *getCurrents() { return currents; }
};

const double mincond = 0.005;
const double maxcond = 0.3815;

#endif // PROBLEM_H_