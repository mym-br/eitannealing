/*
* initproblem.cpp
*
* 	Read from file
*
*  Created on: Feb 17, 2011
*      Author: thiago
*/

#include "problemdescription.h"

#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include <iostream>
#include <boost/concept_check.hpp>

//#include "gradientnormregularisation.h"


std::ifstream file;

// Get the index of a node in the grid (no electrodes then) from
//	its x, y coordinates
inline int getElementIndex(int x, int y) { return 31+17*x+y; }

std::vector<node> nodes;
std::vector<triangularElement> elements;
std::vector<genericEletrode> gelectrodes;
std::vector<std::pair<int, int> > innerAdjacency;
std::vector<std::pair<int, int> > perimeter;
std::map<int, int> node2coefficient;
int numcoefficients;
float electrodeh;
float totalheight;

int groundNode;

void addToGenericElectrode(int n1, int n2, int n3)
{
	// Boost Lambda for locally-defined functor
	// STL has no support for pointer to member :(

	// find the electrode node
	std::vector<genericEletrode>::iterator electrode = 
		std::find_if(gelectrodes.begin(), gelectrodes.end(),
		[n1](genericEletrode e){return e.baseNode==n1;});

	if(electrode==gelectrodes.end()) {	// New electrode
		genericEletrode aux;
		aux.baseNode = n1;
		std::pair<int, int>nodes(n2, n3);
		aux.nodesPairs.push_back(nodes);
		gelectrodes.push_back(aux);
	} else {
		// add new nodes pair to electrode nodes pair list
		std::pair<int, int>nodes(n2, n3);
		electrode->nodesPairs.push_back(nodes);
	}
}

void fillNodes() {
	//file.seekg(0)
	std::string aux;
	do {
		std::getline(file,aux);
	} while(aux.compare("$NOD")!=0);

	int numNodes;
	file >> numNodes;

	nodes.reserve(numNodes);

	int i;
	for(i=0;i<numNodes;i++) {
		node n;

		file.ignore(256, ' ');
		file >> n.x;
		file >> n.y;
		file.ignore(256, ' ');
		nodes.push_back(n);
	}
}

void fillElementsGenericElectrode() {
	std::set<int> outerRingNodes;
	std::set<int> innerNodes;
	std::string aux;
	do {
		std::getline(file,aux);
	} while(aux.compare("$ELM")!=0);

	triangularElement temp;
	int numElements;	// Notice this includes also electrode elements
	file >> numElements;

	int i;

	elements.reserve(numElements);

	for(i=0;i<numElements;i++) {
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
		switch(id) {
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
		case 4001:
			if(!outerRingNodes.count(n1))
				innerNodes.insert(n1);
			if(!outerRingNodes.count(n2))
				innerNodes.insert(n2);
			if(!outerRingNodes.count(n3))
				innerNodes.insert(n3);
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
			if(outerRingNodes.count(n1)==0) {
				baseNode = n1;
				e1 = n2;
				e2 = n3;
			} else if(outerRingNodes.count(n2)==0) {
				baseNode = n2;
				e1 = n1;
				e2 = n3;
			} 
			else if(outerRingNodes.count(n3)==0) {
				baseNode = n3;
				e1 = n1;
				e2 = n2;
			} 			 			 
			addToGenericElectrode(baseNode, e1, e2);
			break;
		}
	}

	// Prepare node <-> condindex map
	int condIndex = 0;
	// Electrode coefficients
	for(auto e:gelectrodes) {
		node2coefficient[e.baseNode] = condIndex++;
	}

	// Outter ring coefficient
	for(auto e:outerRingNodes) {
		node2coefficient[e] = condIndex;
	}

	condIndex++;

	// Inner coefficients
	for(auto i:innerNodes) {
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
		if(n1<n2) {
			auxAdjacency.insert(std::pair<int,int>(n1, n2));   
		} else {
			auxAdjacency.insert(std::pair<int,int>(n2, n1));
		}
	};
	// For each element, add its node pairs in order
	for(auto e: elements) {
		bool n1ok = (outerRingNodes.find(e.n1)==outerRingNodes.end());
		bool n2ok = (outerRingNodes.find(e.n2)==outerRingNodes.end());
		bool n3ok = (outerRingNodes.find(e.n3)==outerRingNodes.end());

		if(n1ok && n2ok) insertAdjNodePair(e.n1, e.n2);
		if(n1ok && n3ok) insertAdjNodePair(e.n1, e.n3);
		if(n2ok && n3ok) insertAdjNodePair(e.n2, e.n3);
	}

	innerAdjacency.resize(auxAdjacency.size());
	std::copy(auxAdjacency.begin(), auxAdjacency.end(), innerAdjacency.begin());
}



///
// Prepare the perimeter for display porposes
//	This is a hack that presumes the outter elements
//	are composed of JUST A SINGLE LAYER.
void preparePerimeter()
{
    
    // Search for elements that are fully inserted in the outter ring
    //  Since we no longer have access to the InnerNodes/OutterNodes sets,
    //	we must attempt to find it from their condition index
    //	Sadly, this has QUADRATIC complexity in the element count, which is quite bad
        
    // Yes, there's no clean way to find what is the condition index for external
    //	nodes either...
    int index = gelectrodes.size();	// This should be the index number for the external nodes
    // FIXME: Is there's an actual risk of repeating edges here?
    
    auto CheckInternal = [index] (int n) {
      if(node2coefficient[n]>index) return true;
      // Check if there's ANY adjacent node with coefficient > index
      //	1. The element contains the node
      //	2. The element contains ANY node with coefficient > index
      for(auto const &e : elements) {
      	  if(e.n1==n ||e.n2==n ||  e.n3==n) {
	    if(node2coefficient[e.n1]>index) return true;
	    if(node2coefficient[e.n2]>index) return true;
	    if(node2coefficient[e.n3]>index) return true;
	  }
      }
      return false;
    };
    
    for(auto const &elem : elements)
    {
	// We want eleemnts that match the following criteria:
	// Either the 1st and 2nd nodes are external and the 3rd is internal,
	//  or the 1st and 3rd are external and the 2nd is internal
	//  or the 1st is internal and both 2nd and 3rd are external
      
        if(!CheckInternal(elem.n1)) { // 1st is external...
	  if(!CheckInternal(elem.n2)) { // 2nd is external...
	     if(CheckInternal(elem.n3)) { // ... and 3rd is internal, our pair is n1 and n2
		perimeter.push_back(std::make_pair(elem.n1,elem.n2));
	     }
	  } else 
	      if(!CheckInternal(elem.n3)) { // 2nd is internal, 3rd is external, pair is n1 and n3
		perimeter.push_back(std::make_pair(elem.n1,elem.n3));
	      }
	} else
	  if(!CheckInternal(elem.n2)) { // 1st is interal, 2nd is external, check 3rd
	      if(!CheckInternal(elem.n3)) { // 3rd is external, pair is n2 and n3
		perimeter.push_back(std::make_pair(elem.n2,elem.n3));
	      }
	  }
    }
}

void initProblem(char *filename)
{
	file.open(filename);

	fillNodes();
	//fillElements();//	fillElectrodes();
	fillElementsGenericElectrode();
	preparePerimeter();
	electrodeh = 0.0004;
	totalheight = 0.020;

	file.close();
}


