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
#include <set>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <cmath> 
#include "../problem.h"
class solution;
class viewport;

class problem3D : public problem {
	friend class viewport;

	struct node {
		node() {}
		node(const node &n) : x(n.x), y(n.y), z(n.z) {}
		double x, y, z;
	};

	struct tetrahedralElement {
		tetrahedralElement() {}
		tetrahedralElement(int _a, int _b, int _c, int _d) : a(_a), b(_b), c(_c), d(_d) {}
		int a, b, c, d;
		int condIndex;
	};

	struct triangularElement {
		triangularElement() {}
		triangularElement(int _a, int _b, int _c) : a(_a), b(_b), c(_c) {}
		int a, b, c;
	};

	struct genericEletrode {
		int baseNode = -1;
		std::vector<triangularElement> nodesTriangles;
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

	struct gmeshElement {
		int elm_number;
		int elm_type;
		int number_of_tags;
		std::vector<int> tags;
		std::vector<int> node_number_list;
	};

private:
	void fillNodes();
	gmeshElement getNextElement();
	void fillElementsGenericElectrode(bool ignoreouterring);
	void addToGenericElectrode(triangularElement base, int eletrodeTag, std::set<int> &baseNodes);
	void preparePerimeter();

	void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient);
	elementCoefficients calcElementCoefficients(const tetrahedralElement &e);
	void insertNewElementCoefficient(nodeCoefficients **target, int node, const tetrahedralElement &e, double coefficient, const std::map<int, int> &coefficientMap);
	void calcAndInsertGenericElectrodeCoefficients(const genericEletrode &e, const std::vector<node> &nodes, double electrodeh,
		const std::map<int, int> &coefficientMap,
		const std::function<void(int, int, int, double)> &insert);

	std::ifstream file;
	std::vector<node> nodes;
	std::vector<tetrahedralElement> elements;
	std::map<int, genericEletrode> gelectrodes;
	std::vector<std::pair<int, int> > perimeter;

public:
	void initProblem(const char *meshfilename, bool ignoreouterring = false);
	void setCalibrationMode(bool individualcoeffs = false);
	void buildNodeCoefficients();
	int getGenericElectrodesCount() { return (int)gelectrodes.size(); }
	int getInnerAdjacencyCount() { return (int)innerAdjacency.size(); }
	int getGenericElectrodesCoeffCount() { return (int)gelectrodes.size(); };
	problem3D(const char *meshfilename) : problem(meshfilename) {};
	~problem3D(){};
};

#endif // PROBLEM3D_H_