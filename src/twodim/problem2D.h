/*
* problem2D.h
*
*  Created on: Feb 23, 2017
*      Author: aksato
*/

#ifndef PROBLEM2D_H_
#define PROBLEM2D_H_

#include <vector>
#include <map>
#include <set>
#include <fstream>
#include "../problem.h"
class solution;
class viewport;

class problem2D : public problem {
	friend class viewport;
	friend class viewportcomplex;

	struct node {
		double x, y;
	};

	struct triangularElement {
		int a, b, c;
		int condIndex;
	};

	struct genericEletrode {
		int baseNode;
		std::vector<std::pair<int, int> > nodesPairs;
	};

	struct elementCoefficients {
		double aa, bb, cc, ab, ac, bc;
		inline elementCoefficients operator *=(const double x) {
			aa *= x; bb *= x; cc *= x;
			ab *= x; ac *= x; bc *= x;
			return *this;
		}
	};

private:
	void fillNodes();
	void fillElementsGenericElectrode();
	void addToGenericElectrode(int n1, int n2, int n3);
	void preparePerimeter();

	void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient);
	elementCoefficients calcElementCoefficients(const triangularElement &e);
	void insertNewElementCoefficient(nodeCoefficients **target, int node, const triangularElement &e, double coefficient, const std::map<int, int> &coefficientMap);
	void calcAndInsertGenericElectrodeCoefficients(const genericEletrode &e, const std::vector<node> &nodes, double electrodeh, double totalheight,
		const std::map<int, int> &coefficientMap,
		const std::function<void(int, int, int, double)> &insert);

	float totalheight;
	std::ifstream file;
	std::vector<node> nodes;
	std::vector<triangularElement> elements;
	std::vector<genericEletrode> gelectrodes;
	std::vector<std::pair<int, int> > perimeter;

public:
	void initProblem(const char *meshfilename);
	void buildNodeCoefficients();
	int getGenericElectrodesCount() { return (int)gelectrodes.size(); }
	int getInnerAdjacencyCount() { return (int)innerAdjacency.size(); }
	problem2D(const char *meshfilename) : problem(meshfilename) {};
	~problem2D(){};
};

#endif // PROBLEM2D_H_