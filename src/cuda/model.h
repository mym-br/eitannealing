/*
* model.h
*
*  Created on: Jun 4, 2015
*      Author: renato
*/

#ifndef MODEL_H_
#define MODEL_H_

#include <vector>
#include <map>

#include "settings.h"

struct node {
	double x, y;
};

struct triangularElement {
	int n1, n2, n3;
	int condIndex;
};

struct triangularEletrode {
	int baseNode;
	std::vector<std::pair<int, int>> nodesPairs;
};

extern numType electrodeh;
extern numType totalheight;

extern std::vector<node> nodes;
extern std::vector<triangularElement> elements;
extern std::vector<triangularEletrode> electrodes;

extern std::vector<std::pair<int, int> > innerAdjacency;

extern std::map<int, int> node2coefficient;
extern int numcoefficients;

const numType mincond = 0.005;
const numType maxcond = 0.3815;

void initProblem(char *meshfilename);

#endif MODEL_H_