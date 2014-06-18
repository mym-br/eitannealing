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

struct genericEletrode {
	int baseNode;
	std::vector<std::pair<int, int> > nodesPairs;
};

extern float electrodeh;
extern float totalheight;
extern std::vector<node> nodes;
extern std::vector<triangularElement> elements;
//extern std::vector<triangularEletrode> electrodes;
extern std::vector<genericEletrode> gelectrodes;
extern std::vector<std::pair<int, int> > innerAdjacency;

// For displaying only
extern std::vector<std::pair<int, int> > perimeter;

extern std::map<int, int> node2coefficient;
extern int numcoefficients;
extern int groundNode;

const double mincond = 0.005;
const double maxcond = 0.3815;
void initProblem(char *meshfilename);

#endif /* PROBLEMDESCRIPTION_H_ */
