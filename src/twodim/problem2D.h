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

template <typename _Scalar, typename t_vector, typename t_matrix >
class problem2D : public problem <_Scalar, t_vector, t_matrix> {
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
	void fillNodes() {
		//file.seekg(0)
		std::string aux;
		do {
			std::getline(file, aux);
		} while (aux.compare("$NOD") != 0);

		int numNodes;
		file >> numNodes;

		nodes.reserve(numNodes);

		int i;
		for (i = 0; i<numNodes; i++) {
			problem2D::node n;

			file.ignore(256, ' ');
			file >> n.x;
			file >> n.y;
			file.ignore(256, ' ');
			nodes.push_back(n);
		}

		nodeCount = nodes.size();
	}
	void fillElementsGenericElectrode() {
		std::set<int> outerRingNodes;
		std::set<int> innerNodes;
		std::string aux;
		do {
			std::getline(file, aux);
		} while (aux.compare("$ELM") != 0);

		triangularElement temp;
		int numElements;	// Notice this includes also electrode elements
		file >> numElements;

		int i;

		elements.reserve(numElements);

		for (i = 0; i<numElements; i++) {
			int id;
			int na, nb, nc;
			file.ignore(256, ' ');
			file.ignore(256, ' ');
			file >> id;
			file.ignore(256, ' ');
			file.ignore(256, ' ');
			file.ignore(256, ' ');
			file >> na >> nb >> nc;
			// 1 based
			na--;
			nb--;
			nc--;
			switch (id) {
			case 1001:	// external ring
				//case 2001:
				//case 3001:
				innerNodes.erase(na);
				innerNodes.erase(nb);
				innerNodes.erase(nc);
				outerRingNodes.insert(na);
				outerRingNodes.insert(nb);
				outerRingNodes.insert(nc);
				temp.a = na;
				temp.b = nb;
				temp.c = nc;
				elements.push_back(temp);
				break;

			case 2001:
			case 3001:	// internal elements
			case 4001:
				if (!outerRingNodes.count(na))
					innerNodes.insert(na);
				if (!outerRingNodes.count(nb))
					innerNodes.insert(nb);
				if (!outerRingNodes.count(nc))
					innerNodes.insert(nc);
				temp.a = na;
				temp.b = nb;
				temp.c = nc;
				elements.push_back(temp);
				break;

			case 10000:	// electrode
				// For this to work, electrodes bust be the last
				//	entity declared
				// FIXME: Add consistency tests (i.e.: exactly two nodes
				//	should be present in outterRingNodes)
				int baseNode = -1;
				int e1, e2;
				if (outerRingNodes.count(na) == 0) {
					baseNode = na;
					e1 = nb;
					e2 = nc;
				}
				else if (outerRingNodes.count(nb) == 0) {
					baseNode = nb;
					e1 = na;
					e2 = nc;
				}
				else if (outerRingNodes.count(nc) == 0) {
					baseNode = nc;
					e1 = na;
					e2 = nb;
				}
				addToGenericElectrode(baseNode, e1, e2);
				break;
			}
		}

		// Prepare node <-> condindex map
		int condIndex = 0;
		// Electrode coefficients
		for (auto e : gelectrodes) {
			node2coefficient[e.baseNode] = condIndex++;
		}

		// Outter ring coefficient
		for (auto e : outerRingNodes) {
			node2coefficient[e] = condIndex;
		}

		condIndex++;

		// Inner coefficients
		for (auto i : innerNodes) {
			node2coefficient[i] = condIndex++;
		}

		numcoefficients = condIndex;
		if (getGroundNode() == -1) setGroundNode(gelectrodes.back().baseNode);

		// Prepare inner nodes adjacency map
		//	Adjacency is established between nodes that are NOT in the outter ring
		//	set AND share at least one element
		typedef std::set<std::pair<int, int> > adjacencySet;
		adjacencySet auxAdjacency;
		// Functor that adds an ordered pair to innerAdjacency
		//  Notice the pair is ordered first, so the pair (2,1) is translated
		//	to (1,2)
		auto insertAdjNodePair = [&auxAdjacency](int n1, int n2) {
			if (n1<n2) {
				auxAdjacency.insert(std::pair<int, int>(n1, n2));
			}
			else {
				auxAdjacency.insert(std::pair<int, int>(n2, n1));
			}
		};
		// For each element, add its node pairs in order
		for (auto e : elements) {
			bool naok = (outerRingNodes.find(e.a) == outerRingNodes.end());
			bool nbok = (outerRingNodes.find(e.b) == outerRingNodes.end());
			bool ncok = (outerRingNodes.find(e.c) == outerRingNodes.end());

			if (naok && nbok) insertAdjNodePair(e.a, e.b);
			if (naok && ncok) insertAdjNodePair(e.a, e.c);
			if (nbok && ncok) insertAdjNodePair(e.b, e.c);
		}

		innerAdjacency.resize(auxAdjacency.size());
		std::copy(auxAdjacency.begin(), auxAdjacency.end(), innerAdjacency.begin());
	}

	void addToGenericElectrode(int n1, int n2, int n3) {
		// Boost Lambda for locally-defined functor
		// STL has no support for pointer to member :(

		// find the electrode node
		std::vector<genericEletrode>::iterator electrode =
			std::find_if(gelectrodes.begin(), gelectrodes.end(),
			[n1](genericEletrode e){return e.baseNode == n1; });

		if (electrode == gelectrodes.end()) {	// New electrode
			genericEletrode aux;
			aux.baseNode = n1;
			std::pair<int, int>nodes(n2, n3);
			aux.nodesPairs.push_back(nodes);
			gelectrodes.push_back(aux);
		}
		else {
			// add new nodes pair to electrode nodes pair list
			std::pair<int, int>nodes(n2, n3);
			electrode->nodesPairs.push_back(nodes);
		}
	}

	void preparePerimeter() {

		// Search for elements that are fully inserted in the outter ring
		//  Since we no longer have access to the InnerNodes/OutterNodes sets,
		//	we must attempt to find it from their condition index
		//	Sadly, this has QUADRATIC complexity in the element count, which is quite bad

		// Yes, there's no clean way to find what is the condition index for external
		//	nodes either...
		int index = (int)gelectrodes.size();	// This should be the index number for the external nodes
		// FIXME: Is there's an actual risk of repeating edges here?

		auto CheckInternal = [index, this](int n) {
			if (node2coefficient[n]>index) return true;
			// Check if there's ANY adjacent node with coefficient > index
			//	1. The element contains the node
			//	2. The element contains ANY node with coefficient > index
			for (auto const &e : elements) {
				if (e.a == n || e.b == n || e.c == n) {
					if (node2coefficient[e.a]>index) return true;
					if (node2coefficient[e.b]>index) return true;
					if (node2coefficient[e.c]>index) return true;
				}
			}
			return false;
		};

		for (auto const &elem : elements)
		{
			// We want eleemnts that match the following criteria:
			// Either the 1st and 2nd nodes are external and the 3rd is internal,
			//  or the 1st and 3rd are external and the 2nd is internal
			//  or the 1st is internal and both 2nd and 3rd are external

			if (!CheckInternal(elem.a)) { // 1st is external...
				if (!CheckInternal(elem.b)) { // 2nd is external...
					if (CheckInternal(elem.c)) { // ... and 3rd is internal, our pair is a and b
						perimeter.push_back(std::make_pair(elem.a, elem.b));
					}
				}
				else
				if (!CheckInternal(elem.c)) { // 2nd is internal, 3rd is external, pair is a and c
					perimeter.push_back(std::make_pair(elem.a, elem.c));
				}
			}
			else
			if (!CheckInternal(elem.b)) { // 1st is interal, 2nd is external, check 3rd
				if (!CheckInternal(elem.c)) { // 3rd is external, pair is b and c
					perimeter.push_back(std::make_pair(elem.b, elem.c));
				}
			}
		}
	}

	void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient) {
		while (*target && (*target)->node < node) target = &(*target)->next;
		while (*target && (*target)->node == node)  {// special case, possible merge
			if ((*target)->condIndex == index) { // Merge
				(*target)->coefficient += coefficient;
				return;
			}
			target = &(*target)->next;
		}
		// Insert node
		*target = new nodeCoefficients(node, index, coefficient, *target);
	}

	elementCoefficients calcElementCoefficients(const triangularElement &e) {
		elementCoefficients c; // Hopefully copy-elision happens here

		// Vectors of each node are the opposed edge rotated by 90 degrees
		// Notice in the end it doesn't matter if the nodes are clockwise or not!
		Eigen::Vector2d
			va(nodes[e.b].y - nodes[e.c].y,
			nodes[e.c].x - nodes[e.b].x),
			vb(nodes[e.c].y - nodes[e.a].y,
			nodes[e.a].x - nodes[e.c].x),
			vc(nodes[e.a].y - nodes[e.b].y,
			nodes[e.b].x - nodes[e.a].x);

		double areaFactor = 2 * fabs(va.x()*vb.y() - va.y()*vb.x());

		c.ab = va.dot(vb) / areaFactor;
		c.ac = va.dot(vc) / areaFactor;
		c.bc = vb.dot(vc) / areaFactor;
		c.aa = -c.ab - c.ac;
		c.bb = -c.ab - c.bc;
		c.cc = -c.ac - c.bc;

		return c;
	}

	void insertNewElementCoefficient(nodeCoefficients **target, int node, const triangularElement &e, double coefficient, const std::map<int, int> &coefficientMap) {
		// Coefficient is spread among the 3 nodes
		coefficient /= 3;
		std::map<int, int>::const_iterator i;
		if ((i = coefficientMap.find(e.a)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
		if ((i = coefficientMap.find(e.b)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
		if ((i = coefficientMap.find(e.c)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
	}

	void calcAndInsertGenericElectrodeCoefficients(const genericEletrode &e, const std::vector<node> &nodes, double electrodeh, double totalheight,
		const std::map<int, int> &coefficientMap,
		const std::function<void(int, int, int, double)> &insert) {
		for (auto p : e.nodesPairs) {
			// Vij: nj -> ni
			Eigen::Vector2d	vij(nodes[p.second].x - nodes[p.first].x,
				nodes[p.second].y - nodes[p.first].y);
			// A = h / (2*l); B = l / (2*h) => A = 1 / (4*B) = 0.25/B
			double B = vij.norm() / (2 * electrodeh);
			double A = 0.25 / B;

			double kiie = A + B;	// (h+l)/(2hl)
			double kjje = kiie; // (h+l)/(2hl)
			double kije = -A;	// -h/(2hl)
			double kbie = -B;	// -l/(2hl)
			double kbje = -B;	// -l/(2hl)
			double kbbe = 2 * B;	// 2*l/(2hl) = l/(hl) = l/h

			// multiplicar pela altura (totalheight) e somar aos acumuladores de coeficiente
			std::map<int, int>::const_iterator ii = coefficientMap.find(e.baseNode);
			if (ii != coefficientMap.end()) {
				int index = ii->second;
				// Add to the base node...
				insert(e.baseNode, e.baseNode, index, kbbe * totalheight);
				insert(e.baseNode, p.first, index, kbie * totalheight);
				insert(e.baseNode, p.second, index, kbje * totalheight);

				// the i-th node...
				insert(p.first, p.first, index, kiie * totalheight);
				insert(p.first, p.second, index, kije * totalheight);	// FIXME: Necessary?
				insert(p.first, e.baseNode, index, kbie * totalheight);

				// ... and the j-th node
				insert(p.second, p.second, index, kjje * totalheight);
				insert(p.second, p.first, index, kije * totalheight);	// FIXME: Necessary?
				insert(p.second, e.baseNode, index, kbje * totalheight);
			}
		}
	}

	float totalheight;
	std::ifstream file;
	std::vector<node> nodes;
	std::vector<triangularElement> elements;
	std::vector<genericEletrode> gelectrodes;
	std::vector<std::pair<int, int> > perimeter;

public:
	void initProblem(const char *meshfilename) {
		file.open(meshfilename);

		fillNodes();
		fillElementsGenericElectrode();
		preparePerimeter();
		this->electrodeh = 0.0004f;
		this->totalheight = 0.020f;

		file.close();
	}

	void buildNodeCoefficients() {
		// Init coefficients;
		nodeCoef = new nodeCoefficients *[nodes.size()];
		for (int i = 0; i<nodes.size(); i++) nodeCoef[i] = NULL;
		// Build electrodes
		auto insertElectrodeCoefficient = [this](int node_a, int node_b, int condIndex, double cab) {
			insertNewCoefficient(&nodeCoef[node_a], node_b, condIndex, cab);
		};
		for (const genericEletrode &e : gelectrodes) calcAndInsertGenericElectrodeCoefficients(
			e, nodes, electrodeh, totalheight, node2coefficient, insertElectrodeCoefficient);

			// Now prepare the coefficients due to the elements
		for (const triangularElement e : elements) {
			elementCoefficients c = calcElementCoefficients(e);
			c *= totalheight;
			// Node 1
			insertNewElementCoefficient(&nodeCoef[e.a], e.a, e, c.aa, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.a], e.b, e, c.ab, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.a], e.c, e, c.ac, node2coefficient);
			// Node 2
			insertNewElementCoefficient(&nodeCoef[e.b], e.b, e, c.bb, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.b], e.a, e, c.ab, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.b], e.c, e, c.bc, node2coefficient);
			// Node 3
			insertNewElementCoefficient(&nodeCoef[e.c], e.c, e, c.cc, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.c], e.a, e, c.ac, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.c], e.b, e, c.bc, node2coefficient);
		}
	}

	int getGenericElectrodesCount() { return (int)gelectrodes.size(); }
	int getInnerAdjacencyCount() { return (int)innerAdjacency.size(); }
	problem2D(const char *meshfilename) : problem(meshfilename) {};
	~problem2D(){};
};

#endif // PROBLEM2D_H_