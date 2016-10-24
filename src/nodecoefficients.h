/*
 * nodecoefficients.h
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#ifndef NODECOEFFICIENTS_H_
#define NODECOEFFICIENTS_H_

#include <stdlib.h>
#include "problemdescription.h"

struct nodeCoefficients {
	int node;
	int condIndex;
	double coefficient;
	nodeCoefficients *next;

	nodeCoefficients(int node, int index, double coefficient):
		node(node), condIndex(index), coefficient(coefficient), next(NULL) {}

	nodeCoefficients(int node, int index, double coefficient, nodeCoefficients *next):
			node(node), condIndex(index), coefficient(coefficient), next(next) {}
};

struct elementCoefficients {
  double aa, bb, cc, ab, ac, bc; 
  inline elementCoefficients operator *=(const double x) {
      aa *= x; bb *= x; cc *= x;
      ab *= x; ac *= x; bc *= x;
      return *this;
  }
};

extern nodeCoefficients **nodeCoef;

void buildNodeCoefficients();
void calcGenericElectrodeCoefficients(int electrode, const std::vector<genericEletrode> &electrodes);
void insertNewElementCoefficient(nodeCoefficients **target, int node, int element, double coefficient, std::map<int, int> &coefficientMap);

#endif /* NODECOEFFICIENTS_H_ */
