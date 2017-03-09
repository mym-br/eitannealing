#include "problem3D.h"
#include <set>

void problem3D::initProblem(char *meshfilename) {
	file.open(meshfilename);

	fillNodes();
	fillElementsGenericElectrode();
	this->electrodeh = 0.0004f;

	file.close();
}

void problem3D::fillNodes() {
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
}

problem3D::gmeshElement problem3D::getNextElement() {
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

void problem3D::addToGenericElectrode(triangularElement base, int eletrodeTag, std::set<int> &baseNodes) {
	// Attempt to get existing electrode
	if (gelectrodes.find(eletrodeTag) == gelectrodes.end()) {
		// Create new electrode
		gelectrodes[eletrodeTag] = genericEletrode();
	}
	// Check if any triangle node is in the base node list
	std::set<int>::iterator it;
	if (gelectrodes[eletrodeTag].baseNode == -1 && (
		(it = baseNodes.find(base.a)) != baseNodes.end() || (it = baseNodes.find(base.b)) != baseNodes.end() || (it = baseNodes.find(base.c)) != baseNodes.end()) ) {
		problem3D::node n(nodes[*it]);
		nodes.push_back(n);
		gelectrodes[eletrodeTag].baseNode = (int)nodes.size() - 1;
		baseNodes.erase(it); baseNodes.insert((int)nodes.size() - 1);
	}
	// Add base triangle part to electrode
	gelectrodes[eletrodeTag].nodesTriangles.push_back(base);
	// Search for base node!
}

void problem3D::fillElementsGenericElectrode() {
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
	groundNode = (*(--gelectrodes.end())).second.baseNode;

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