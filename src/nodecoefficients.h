#ifndef NODECOEFFICIENTS_H_
#define NODECOEFFICIENTS_H_

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

#endif // NODECOEFFICIENTS_H_