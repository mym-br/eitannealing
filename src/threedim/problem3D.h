/*
* problem2D.h
*
*  Created on: Feb 23, 2017
*      Author: aksato
*/

#ifndef PROBLEM3D_H_
#define PROBLEM3D_H_

#include <vector>
#include <map>
#include <fstream>
#include "../problem.h"
class solution;
class viewport;

class problem3D : public problem {
	friend class viewport;

	struct node {
		double x, y, z;
	};

	struct tetrahedralElement {
		int a, b, c, d;
		int condIndex;
	};

	struct triangularElement {
		int a, b, c;
	};

	struct genericEletrode {
		int baseNode;
		std::vector<triangularElement > nodesTriangles;
	};

	struct elementCoefficients {
		double aa, bb, cc, dd, ab, ac, ad, bc, bd, cd;
		inline elementCoefficients operator *=(const double x) {
			aa *= x; bb *= x; cc *= x; dd *= x;
			ab *= x; ac *= x; bc *= x; bd*= x;
			ad *= x; cd *= x; 
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
	matrix *skeleton;
	matrix *coef2KMatrix;

public:
	void initProblem(char *meshfilename);
	void buildNodeCoefficients();
	void prepareSkeletonMatrix();
	void createCoef2KMatrix();
	void assembleProblemMatrix(double *cond, matrix **stiffnes);
	int getGenericElectrodesCount() { return (int)gelectrodes.size(); }
	int getNodesCount() { return (int)nodes.size(); }
	int getInnerAdjacencyCount() { return (int)innerAdjacency.size(); }
	~problem3D(){};
};

#endif // PROBLEM3D_H_