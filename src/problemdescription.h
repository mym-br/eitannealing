/*
 * problemdescription.h
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#ifndef PROBLEMDESCRIPTION_H_
#define PROBLEMDESCRIPTION_H_

#include <vector>
#include <map>

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
	//int getBaseNode() const {return this->baseNode;}
};

extern float electrodeh;
extern float totalheight;
extern std::vector<node> nodes;
extern std::vector<triangularElement> elements;
extern std::vector<triangularEletrode> electrodes;

extern std::map<int, int> node2coefficient;
extern int numcoefficients;
extern int groundNode;


void initProblem(char *meshfilename);

#endif /* PROBLEMDESCRIPTION_H_ */
