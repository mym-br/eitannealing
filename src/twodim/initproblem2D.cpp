#include "problem2D.h"
#include <set>

void problem2D::initProblem(const char *meshfilename) {
	file.open(meshfilename);

	fillNodes();
	fillElementsGenericElectrode();
	preparePerimeter();
	this->electrodeh = 0.0004f;
	this->totalheight = 0.027f;

	file.close();
}

void problem2D::fillNodes() {
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
}

void problem2D::fillElementsGenericElectrode() {
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
	groundNode = gelectrodes.back().baseNode;

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
		if (nbok && naok) insertAdjNodePair(e.a, e.b);
		if (nbok && ncok) insertAdjNodePair(e.b, e.c);
	}

	innerAdjacency.resize(auxAdjacency.size());
	std::copy(auxAdjacency.begin(), auxAdjacency.end(), innerAdjacency.begin());
}

void problem2D::addToGenericElectrode(int n1, int n2, int n3)
{
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

///
// Prepare the perimeter for display porposes
//	This is a hack that presumes the outter elements
//	are composed of JUST A SINGLE LAYER.
void problem2D::preparePerimeter()
{

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