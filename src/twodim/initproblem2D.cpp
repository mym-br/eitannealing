#include "problem2D.h"

void problem2D::initProblem(const char *meshfilename) {
	file.open(meshfilename);

	// Check mesh file version
	int meshVer = 1;
	std::string aux;
	std::getline(file, aux);
	if(aux.compare("$MeshFormat") == 0) {
		std::string numVer; file >> numVer;
		if (numVer[0] == '2') meshVer = 2;
		else { std::cerr << "Incompatible mesh file version: " << numVer << std::endl; return; }
		int binary; file >> binary;
		if(binary != 0) { std::cerr << "Only compatible with ASCII mesh file" << std::endl; return; }
	}
	else file.seekg(0, std::ios::beg);

	fillNodes(meshVer);
	fillElementsGenericElectrode(meshVer);
	preparePerimeter();
	//this->electrodeh = 0.023f;
	//this->totalheight = 0.016f;
	this->electrodeh = 0.0004f;
	this->totalheight = 0.020f;

	file.close();
}

void problem2D::fillNodes(int meshVer) {
	//file.seekg(0)
	std::string nodeStartSection(meshVer == 1 ? "$NOD" : "$Nodes");
	std::string aux;
	do {
		std::getline(file, aux);
	} while (aux.compare(nodeStartSection) != 0);

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

void problem2D::fillElementsGenericElectrode(int meshVer) {
	std::set<int> outerRingNodes;
	std::set<int> innerNodes;
	std::string elmStartSection(meshVer == 1 ? "$ELM" : "$Elements");
	std::string aux;
	do {
		std::getline(file, aux);
	} while (aux.compare(elmStartSection) != 0);

	triangularElement temp;
	int numElements;	// Notice this includes also electrode elements
	file >> numElements;

	int i;

	elements.reserve(numElements);

	int lastBaseNode = 0;
	for (i = 0; i<numElements; i++) {
		int id, numTags;
		int na, nb, nc;
		file.ignore(256, ' ');
		file.ignore(256, ' ');
		if (meshVer == 2) file >> numTags;
		else numTags = 3;
		file >> id;
		for(int i = 0; i < numTags; i++) file.ignore(256, ' ');
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

		case 10001: // multilayered electrode middle layer
			if (!outerRingNodes.count(na))
				gelectrodesNonBaseNodes.insert(na);
			if (!outerRingNodes.count(nb))
				gelectrodesNonBaseNodes.insert(nb);
			if (!outerRingNodes.count(nc))
				gelectrodesNonBaseNodes.insert(nc);
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
			if (outerRingNodes.count(na) == 0 && gelectrodesNonBaseNodes.count(na) == 0) {
				baseNode = na;
				e1 = nb;
				e2 = nc;
			}
			else if (outerRingNodes.count(nb) == 0 && gelectrodesNonBaseNodes.count(nb) == 0) {
				baseNode = nb;
				e1 = na;
				e2 = nc;
			}
			else if (outerRingNodes.count(nc) == 0 && gelectrodesNonBaseNodes.count(nc) == 0) {
				baseNode = nc;
				e1 = na;
				e2 = nb;
			}
			if (baseNode > lastBaseNode) lastBaseNode = baseNode;
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
	for (auto e : gelectrodesNonBaseNodes) {
		node2coefficient[e] = condIndex++;
	}

	// Outter ring coefficient
	for (auto e : outerRingNodes) {
		node2coefficient[e] = condIndex;
		if (this->ignoreouterring) condIndex++;
	}

	if (!outerRingNodes.empty() && !this->ignoreouterring) condIndex++;

	// Inner coefficients
	for (auto i : innerNodes) {
		node2coefficient[i] = condIndex++;
	}

	numcoefficients = condIndex;
	if (getGroundNode() == -1) setGroundNode(lastBaseNode);

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

void problem2D::setCalibrationMode(bool individualcoeffs) {
	// Set all coeffs to 0
	for (int i = 0; i < getNodesCount(); i++) node2coefficient[i] = 0;
	// Electrode coefficients are 1
	int k = 1;
	std::map<int, int> baseNodeCoeffs;
	for (auto e : gelectrodes) {
		if (individualcoeffs) {
			if (baseNodeCoeffs.find(e.baseNode) == baseNodeCoeffs.end())
				baseNodeCoeffs.insert(std::pair<int, int>(e.baseNode, k++));
			node2coefficient[e.baseNode] = baseNodeCoeffs[e.baseNode];
		}
		else node2coefficient[e.baseNode] = 1;
	}
	// Update coefficients count
	numcoefficients = individualcoeffs ? baseNodeCoeffs.size() + 1 : 2;
	// Set calibration mode
	calibrationMode = individualcoeffs ? 2 : 1;
}

void problem2D::addToGenericElectrode(int n1, int n2, int n3) {
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

void problem2D::preparePerimeter() {

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
