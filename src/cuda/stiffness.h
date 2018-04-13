#ifndef STIFFNESS_H
#define STIFFNESS_H

#include <vector>
#include "settings.h"

struct NodeCoefficients {
	int node; // node index
	int condutanceIndex; // condutance (solution) vector index
	numType coefficient; // condutance coefficient
	NodeCoefficients *next;

	NodeCoefficients(int node, int solIndex, numType coefficient) :
		node(node), condutanceIndex(solIndex), coefficient(coefficient), next(NULL) {}

	NodeCoefficients(int node, int solIndex, numType coefficient, NodeCoefficients *next) :
		node(node), condutanceIndex(solIndex), coefficient(coefficient), next(next) {}
};

struct CondutanceCoefficients {
	int row; // node's row
	int col; // node's column
	numType coefficient; // this node's (row, column) condutance coefficient

	CondutanceCoefficients(int row, int col, numType coefficient) :
		row(row), col(col), coefficient(coefficient) {}
};

//typedef std::vector<NodeCoefficients> NodeCoefficientsVector;

//extern NodeCoefficientsVector * nodeCoefficients;

extern NodeCoefficients **nodeCoef;

typedef std::vector<CondutanceCoefficients> CondutanceCoefficientsVector;

extern CondutanceCoefficientsVector * condutanceCoefficients; // list of nodes that are affected by a given condutance (solution)

void buildNodeCoefficients();

#endif /* STIFFNESS_H */