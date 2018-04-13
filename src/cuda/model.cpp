#include "model.h"

#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include <iostream>

std::ifstream file;

std::vector<node> nodes;
std::vector<triangularElement> elements;
std::vector<triangularEletrode> electrodes;

std::vector<std::pair<int, int> > innerAdjacency;

std::map<int, int> node2coefficient;
int numcoefficients;

numType electrodeh;
numType totalheight;

void addToElectrode(int n1, int n2, int n3) {
	// find the electrode node
	std::vector<triangularEletrode>::iterator electrode =
		std::find_if(electrodes.begin(), electrodes.end(), [n1](triangularEletrode e) {
			return e.baseNode == n1;
		});

	if (electrode == electrodes.end()) {	// New electrode
		triangularEletrode aux;
		aux.baseNode = n1;
		std::pair<int, int>nodes(n2, n3);
		aux.nodesPairs.push_back(nodes);
		electrodes.push_back(aux);
	}
	else {
		// add new nodes pair to electrode nodes pair list
		std::pair<int, int>nodes(n2, n3);
		electrode->nodesPairs.push_back(nodes);
	}
}

void fillNodes() {
	//file.seekg(0)
	std::string aux;
	do {
		std::getline(file, aux);
	} while (aux.compare("$NOD") != 0);

	int numNodes;
	file >> numNodes;

	nodes.reserve(numNodes);

	for (int i = 0; i < numNodes; i++) {
		node n;

		file.ignore(256, ' ');
		file >> n.x;
		file >> n.y;
		file.ignore(256, ' ');
		nodes.push_back(n);
	}
}

void fillElementsElectrode() {
	std::set<int> outerRingNodes;
	std::set<int> innerNodes;
	std::string aux;
	do {
		std::getline(file, aux);
	} while (aux.compare("$ELM") != 0);

	triangularElement temp;
	int numElements;	// Notice this includes also electrode elements
	file >> numElements;

	elements.reserve(numElements);

	for (int i = 0; i < numElements; i++) {
		int id;
		int n1, n2, n3;
		file.ignore(256, ' ');
		file.ignore(256, ' ');
		file >> id;
		file.ignore(256, ' ');
		file.ignore(256, ' ');
		file.ignore(256, ' ');
		file >> n1 >> n2 >> n3;
		// 1 based
		n1--;
		n2--;
		n3--;
		switch (id) {
		case 1001:	// external ring
			//case 2001:
			//case 3001:
			innerNodes.erase(n1);
			innerNodes.erase(n2);
			innerNodes.erase(n3);
			outerRingNodes.insert(n1);
			outerRingNodes.insert(n2);
			outerRingNodes.insert(n3);
			temp.n1 = n1;
			temp.n2 = n2;
			temp.n3 = n3;
			elements.push_back(temp);
			break;

		case 2001:
		case 3001:	// internal elements
			if (!outerRingNodes.count(n1)) {
				innerNodes.insert(n1);
			}
			if (!outerRingNodes.count(n2)) {
				innerNodes.insert(n2);
			}
			if (!outerRingNodes.count(n3)) {
				innerNodes.insert(n3);
			}
			temp.n1 = n1;
			temp.n2 = n2;
			temp.n3 = n3;
			elements.push_back(temp);
			break;

		case 10000:	// electrode
			// For this to work, electrodes bust be the last
			//	entity declared
			// FIXME: Add consistency tests (i.e.: exactly two nodes
			//	should be present in outterRingNodes)
			int baseNode = -1;
			int e1, e2;
			if (outerRingNodes.count(n1) == 0) {
				baseNode = n1;
				e1 = n2;
				e2 = n3;
			}
			else if (outerRingNodes.count(n2) == 0) {
				baseNode = n2;
				e1 = n1;
				e2 = n3;
			}
			else if (outerRingNodes.count(n3) == 0) {
				baseNode = n3;
				e1 = n1;
				e2 = n2;
			}
			addToElectrode(baseNode, e1, e2);
			break;
		}
	}

	// Prepare node <-> condindex map
	int condIndex = 0;
	// Electrode coefficients
	for (auto e : electrodes) {
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
		bool n1ok = (outerRingNodes.find(e.n1) == outerRingNodes.end());
		bool n2ok = (outerRingNodes.find(e.n2) == outerRingNodes.end());
		bool n3ok = (outerRingNodes.find(e.n3) == outerRingNodes.end());

		if (n1ok && n2ok) insertAdjNodePair(e.n1, e.n2);
		if (n1ok && n3ok) insertAdjNodePair(e.n1, e.n3);
		if (n2ok && n3ok) insertAdjNodePair(e.n2, e.n3);
	}

	innerAdjacency.resize(auxAdjacency.size());
	std::copy(auxAdjacency.begin(), auxAdjacency.end(), innerAdjacency.begin());
}

void initProblem(char *filename) {
	file.open(filename);

	fillNodes();
	fillElementsElectrode();
	electrodeh = 0.0004;
	totalheight = 0.020;

	file.close();
}