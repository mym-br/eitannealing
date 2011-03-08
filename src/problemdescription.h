/*
 * problemdescription.h
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#ifndef PROBLEMDESCRIPTION_H_
#define PROBLEMDESCRIPTION_H_

#include <vector>

struct node {
	double x, y;
};

struct triangularElement {
	int n1, n2, n3;
	int condIndex;
};

structt riangularEletrode {
	int baseNode;
	int n1, n2, n3;
};

extern std::vector<node> nodes;
extern std::vector<triangularElement> elements;
extern std::vector<triangularEletrode> electrodes;

void initProblem(char *meshfilename);

#endif /* PROBLEMDESCRIPTION_H_ */
