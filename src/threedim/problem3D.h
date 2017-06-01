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

template <typename _Scalar, typename t_vector, typename t_matrix >
class problem3D : public problem <_Scalar, t_vector, t_matrix> {
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
	void fillNodes() {
		//file.seekg(0)
		std::string aux;
		do {
			std::getline(file, aux);
		} while (aux.compare("$Nodes") != 0);

		int numNodes;
		file >> numNodes;

		nodes.reserve(numNodes);

		int i;
		for (i = 0; i<numNodes; i++) {
			problem3D::node n;

			file.ignore(256, ' ');
			file >> n.x;
			file >> n.y;
			file >> n.z;
			nodes.push_back(n);
		}

		nodeCount = nodes.size();
	}
	gmeshElement getNextElement() {
		gmeshElement element;
		int id, eltype, tagCount, temp;
		file >> id >> eltype >> tagCount;
		element.elm_number = id;
		element.elm_type = eltype;
		element.number_of_tags = tagCount;
		for (int i = 0; i < tagCount; i++) {
			file >> temp;
			element.tags.push_back(temp);
		}
		int nodeCount;
		switch (eltype) {
		case 15: nodeCount = 1; break; // 1-node point 
		case  2: nodeCount = 3; break; // 3-node triangle
		case  4: nodeCount = 4; break; // 4-node tetrahedron
		}
		for (int i = 0; i < nodeCount; i++) {
			file >> temp;
			element.node_number_list.push_back(--temp);
		}
		return element;
	}

	void fillElementsGenericElectrode() {
		std::set<int> outerRingNodes;
		std::set<int> innerNodes;
		std::set<int> baseElectrodeNodes;
		std::string aux;
		problem3D::gmeshElement currElement;
		std::vector<problem3D::gmeshElement> gmeshGeElectrodes;

		do {
			std::getline(file, aux);
		} while (aux.compare("$Elements") != 0);

		int numElements;	// Notice this includes also electrode elements
		file >> numElements;

		int i;

		elements.reserve(numElements);

		for (i = 0; i<numElements; i++) {
			currElement = getNextElement();
			int firstTag = *currElement.tags.begin();
			if (firstTag == 1001) {
				// external ring
				for (std::vector<int>::iterator it = currElement.node_number_list.begin(); it != currElement.node_number_list.end(); it++) {
					innerNodes.erase(*it); outerRingNodes.insert(*it);
				}
				elements.push_back(tetrahedralElement(currElement.node_number_list[0], currElement.node_number_list[1], currElement.node_number_list[2], currElement.node_number_list[3]));
			}
			if (firstTag == 2001) {
				// internal elements
				for (std::vector<int>::iterator it = currElement.node_number_list.begin(); it != currElement.node_number_list.end(); it++)
				if (!outerRingNodes.count(*it))
					innerNodes.insert(*it);
				elements.push_back(tetrahedralElement(currElement.node_number_list[0], currElement.node_number_list[1], currElement.node_number_list[2], currElement.node_number_list[3]));
			}
			if (firstTag > 3000 && firstTag < 4000) {
				// Add to a list to finish processing afterwards
				gmeshGeElectrodes.push_back(currElement);
			}
			if (firstTag == 9999) {
				// Electrode base node
				baseElectrodeNodes.insert(currElement.node_number_list[0]);
			}
		}

		// Process electrodes
		for (auto ge : gmeshGeElectrodes)
			addToGenericElectrode(triangularElement(ge.node_number_list[0], ge.node_number_list[1], ge.node_number_list[2]), *ge.tags.begin(), baseElectrodeNodes);

		// Prepare node <-> condindex map
		int condIndex = 0;
		// Electrode coefficients
		for (auto e : gelectrodes) {
			node2coefficient[e.second.baseNode] = condIndex++;
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
		if (getGroundNode() == -1) setGroundNode((*(--gelectrodes.end())).second.baseNode);

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
			bool ndok = (outerRingNodes.find(e.d) == outerRingNodes.end());

			if (naok && nbok) insertAdjNodePair(e.a, e.b);
			if (naok && ncok) insertAdjNodePair(e.a, e.c);
			if (naok && ndok) insertAdjNodePair(e.a, e.d);
			if (nbok && ncok) insertAdjNodePair(e.b, e.c);
			if (nbok && ndok) insertAdjNodePair(e.b, e.d);
			if (ncok && ndok) insertAdjNodePair(e.c, e.d);
		}

		innerAdjacency.resize(auxAdjacency.size());
		std::copy(auxAdjacency.begin(), auxAdjacency.end(), innerAdjacency.begin());
	}

	void addToGenericElectrode(triangularElement base, int eletrodeTag, std::set<int> &baseNodes) {
		// Attempt to get existing electrode
		if (gelectrodes.find(eletrodeTag) == gelectrodes.end()) {
			// Create new electrode
			gelectrodes[eletrodeTag] = genericEletrode();
		}
		// Check if any triangle node is in the base node list
		std::set<int>::iterator it;
		if (gelectrodes[eletrodeTag].baseNode == -1 && (
			(it = baseNodes.find(base.a)) != baseNodes.end() || (it = baseNodes.find(base.b)) != baseNodes.end() || (it = baseNodes.find(base.c)) != baseNodes.end())) {
			problem3D::node n(nodes[*it]);
			nodes.push_back(n);
			gelectrodes[eletrodeTag].baseNode = (int)nodes.size() - 1;
			baseNodes.erase(it); baseNodes.insert((int)nodes.size() - 1);
		}
		// Add base triangle part to electrode
		gelectrodes[eletrodeTag].nodesTriangles.push_back(base);
		// Search for base node!
	}

	void preparePerimeter();

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

	elementCoefficients calcElementCoefficients(const tetrahedralElement &e) {
		elementCoefficients c; // Hopefully copy-elision happens here

		// Edges vectors
		Eigen::Vector3d
			vjl(nodes[e.d].x - nodes[e.b].x, nodes[e.d].y - nodes[e.b].y, nodes[e.d].z - nodes[e.b].z),
			vjk(nodes[e.c].x - nodes[e.b].x, nodes[e.c].y - nodes[e.b].y, nodes[e.c].z - nodes[e.b].z),
			vik(nodes[e.c].x - nodes[e.a].x, nodes[e.c].y - nodes[e.a].y, nodes[e.c].z - nodes[e.a].z),
			vil(nodes[e.d].x - nodes[e.a].x, nodes[e.d].y - nodes[e.a].y, nodes[e.d].z - nodes[e.a].z),
			vij(nodes[e.b].x - nodes[e.a].x, nodes[e.b].y - nodes[e.a].y, nodes[e.b].z - nodes[e.a].z);

		// Vectors of each node are the opposed edge rotated to coincide with normal of opposed face,
		// order of vertexes does not impact the result
		Eigen::Vector3d va(vjl.cross(vjk)), vb(vik.cross(vil)), vc(vil.cross(vij)), vd(vij.cross(vik));

		double areaFactor = 6 * fabs((vij.cross(vik)).dot(vil));
		c.ab = va.dot(vb) / areaFactor;
		c.ac = va.dot(vc) / areaFactor;
		c.ad = va.dot(vd) / areaFactor;
		c.bc = vb.dot(vc) / areaFactor;
		c.bd = vb.dot(vd) / areaFactor;
		c.cd = vc.dot(vd) / areaFactor;
		c.aa = -c.ab - c.ac - c.ad;
		c.bb = -c.ab - c.bc - c.bd;
		c.cc = -c.ac - c.bc - c.cd;
		c.dd = -c.ad - c.bd - c.cd;

		return c;
	}

	void insertNewElementCoefficient(nodeCoefficients **target, int node, const tetrahedralElement &e, double coefficient, const std::map<int, int> &coefficientMap) {
		// Coefficient is spread among the 4 nodes
		coefficient /= 4;
		std::map<int, int>::const_iterator i;
		if ((i = coefficientMap.find(e.a)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
		if ((i = coefficientMap.find(e.b)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
		if ((i = coefficientMap.find(e.c)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
		if ((i = coefficientMap.find(e.d)) != coefficientMap.end())
			insertNewCoefficient(target, node, i->second, coefficient);
	}

	void calcAndInsertGenericElectrodeCoefficients(const genericEletrode &e, const std::vector<node> &nodes, double electrodeh,
		const std::map<int, int> &coefficientMap,
		const std::function<void(int, int, int, double)> &insert) {
		for (auto p : e.nodesTriangles) {
			// Get triangle base parameters
			double hz = electrodeh;
			node vi(nodes[p.a]), vj(nodes[p.b]), vk(nodes[p.c]);

			Eigen::Vector3d	vik(vk.x - vi.x, vk.y - vi.y, vk.z - vi.z), vij(vj.x - vi.x, vj.y - vi.y, vj.z - vi.z), vjk(vk.x - vj.x, vk.y - vj.y, vk.z - vj.z);
			Eigen::Vector3d projLjLk = (vik.dot(vij) / vij.dot(vij)) * vij;
			Eigen::Vector3d projLiLk = (vjk.dot(vij) / vij.dot(vij)) * vij;
			double xji = vij.norm();
			double hy = (vik - projLjLk).norm();
			double xki = std::copysign(projLjLk.norm(), projLjLk.dot(vij));
			double xkj = std::copysign(projLiLk.norm(), projLiLk.dot(vij));
			double areaFactor = 6 * hy*hz*xji;

			double A, B, C, D, E, F;
			A = hy*hz;
			B = hz*xkj;
			C = hy*xji;
			D = hz*xki;
			E = hz*xji;
			F = hy*xji;

			double kbbe = (B*B + 3 * C * C - 2 * B * (D - E) + (D - E) * (D - E)) / areaFactor;
			double kbie = -(C*C) / areaFactor;
			double kbje = kbie;
			double kbke = kbie;
			double kiie = (A*A + B*B + C*C) / areaFactor;
			double kije = -(A*A + B*D) / areaFactor;
			double kike = (B*E) / areaFactor;
			double kjje = (A*A + C*C + D*D + E*E) / areaFactor;
			double kjke = -(E*(D + E)) / areaFactor;
			double kkke = (C*C + 2 * E*E) / areaFactor;

			// somar aos acumuladores de coeficiente
			std::map<int, int>::const_iterator ii = coefficientMap.find(e.baseNode);
			if (ii != coefficientMap.end()) {
				int index = ii->second;
				// Add to the base node...
				insert(e.baseNode, e.baseNode, index, kbbe);
				insert(e.baseNode, p.a, index, kbie);
				insert(e.baseNode, p.b, index, kbje);
				insert(e.baseNode, p.c, index, kbke);

				// the i-th node...
				insert(p.a, e.baseNode, index, kbie);
				insert(p.a, p.a, index, kiie);
				insert(p.a, p.b, index, kije);
				insert(p.a, p.c, index, kike);

				// the j-th node...
				insert(p.b, e.baseNode, index, kbje);
				insert(p.b, p.a, index, kije);
				insert(p.b, p.b, index, kjje);
				insert(p.b, p.c, index, kjke);

				// ... and the k-th node
				insert(p.c, e.baseNode, index, kbke);
				insert(p.c, p.a, index, kike);
				insert(p.c, p.b, index, kjke);
				insert(p.c, p.c, index, kkke);
			}
		}
	}

	std::ifstream file;
	std::vector<node> nodes;
	std::vector<tetrahedralElement> elements;
	std::map<int, genericEletrode> gelectrodes;
	std::vector<std::pair<int, int> > perimeter;

public:
	void initProblem(const char *meshfilename) {
		file.open(meshfilename);

		fillNodes();
		fillElementsGenericElectrode();
		this->electrodeh = 0.0004f;

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
		for (const std::pair<int, genericEletrode> &e : gelectrodes) calcAndInsertGenericElectrodeCoefficients(
			e.second, nodes, electrodeh, node2coefficient, insertElectrodeCoefficient);

			// Now prepare the coefficients due to the elements
		for (const tetrahedralElement e : elements) {
			elementCoefficients c = calcElementCoefficients(e);

			// Node 1
			insertNewElementCoefficient(&nodeCoef[e.a], e.a, e, c.aa, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.a], e.b, e, c.ab, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.a], e.c, e, c.ac, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.a], e.d, e, c.ad, node2coefficient);
			// Node 2
			insertNewElementCoefficient(&nodeCoef[e.b], e.b, e, c.bb, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.b], e.a, e, c.ab, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.b], e.c, e, c.bc, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.b], e.d, e, c.bd, node2coefficient);
			// Node 3
			insertNewElementCoefficient(&nodeCoef[e.c], e.c, e, c.cc, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.c], e.a, e, c.ac, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.c], e.b, e, c.bc, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.c], e.d, e, c.cd, node2coefficient);
			// Node 4
			insertNewElementCoefficient(&nodeCoef[e.d], e.d, e, c.dd, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.d], e.a, e, c.ad, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.d], e.b, e, c.bd, node2coefficient);
			insertNewElementCoefficient(&nodeCoef[e.d], e.c, e, c.cd, node2coefficient);
		}
	}
	int getGenericElectrodesCount() { return (int)gelectrodes.size(); }
	int getInnerAdjacencyCount() { return (int)innerAdjacency.size(); }
	problem3D(const char *meshfilename) : problem(meshfilename) {};
	~problem3D(){};
};

#endif // PROBLEM3D_H_