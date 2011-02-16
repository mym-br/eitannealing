/*
 * problemdescription.h
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#ifndef PROBLEMDESCRIPTION_H_
#define PROBLEMDESCRIPTION_H_


struct node {
	double x, y;
};

struct triangularElement {
	int n1, n2, n3;
	int condIndex;
};

struct triangularEletrode {
	int baseNode;
	int n1, n2, n3;
};

extern node *nodes;
extern int numNodes;
extern triangularElement *elements;
extern int numElements;
extern triangularEletrode *electrodes;
extern int numElectrodes;

void initProblem();

#endif /* PROBLEMDESCRIPTION_H_ */
