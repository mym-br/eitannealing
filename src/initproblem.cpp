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
#include <boost/lambda/lambda.hpp>


std::ifstream file;

// Get the index of a node in the grid (no electrodes then) from
//	its x, y coordinates
inline int getElementIndex(int x, int y) { return 31+17*x+y; }

std::vector<node> nodes;
std::vector<triangularElement> elements;
std::vector<triangularEletrode> electrodes;
std::map<int, int> node2coefficient;


void addToElectrode(int n1, int n2, int n3)
{
	using namespace boost::lambda;
	// Boost Lambda for locally-defined functor
	// STL has no support for pointer to member :(
		
	// find the electrode node
	std::vector<triangularEletrode>::iterator electrode = 
	    std::find_if(electrodes.begin(), electrodes.end(),
		      ((&_1 ->* &triangularEletrode::baseNode)==n1));
	
	if(electrode==electrodes.end()) {	// New electrode
	      triangularEletrode aux;
	      aux.baseNode = n1;
	      aux.n1 = n2;
	      aux.n2 = n3;
	      electrodes.push_back(aux);
	} else {
	      // Find the missing node
	      if(electrode->n1 == n3 ||
		 electrode->n2 == n3) 
		    electrode->n3 = n2;
	      else
		    electrode->n3 = n3;
		
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

void fillElements() {
	std::set<int> outerRingNodes;
	std::set<int> innerNodes;
	std::string aux;
	do {
	  std::getline(file,aux);
	} while(aux.compare("$ELM")!=0);
  
	triangularElement temp;
	int numElements;	// Notice this includes also electrode elements
	file >> numElements;
	
	elements.reserve(numElements);
	int i;
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
	      
	      case 2001:	// internal elements
			if(!innerNodes.count(n1))
				outerRingNodes.insert(n1);
			if(!innerNodes.count(n2))
				outerRingNodes.insert(n2);
			if(!innerNodes.count(n3))
				outerRingNodes.insert(n3);
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
			 addToElectrode(baseNode, e1, e2);

			 break;
	  }
	}
	// Prepare node <-> condindex map
	using namespace boost::lambda;	// Lambda black magic
	int condIndex = 0;
	// Electrode coefficients
	std::for_each(electrodes.begin(), electrodes.end(),
		var(node2coefficient)[(&_1 ->* &triangularEletrode::baseNode)]
			=var(condIndex)++);

	// Outter ring coefficient
	std::for_each(electrodes.begin(), electrodes.end(),
		var(node2coefficient)[(&_1 ->* &triangularEletrode::baseNode)]
			=condIndex);
	condIndex++;
	// Outter ring coefficient
	std::for_each(electrodes.begin(), electrodes.end(),
		var(node2coefficient)[(&_1 ->* &triangularEletrode::baseNode)]
			=var(condIndex)++);
}

void initProblem(char *filename)
{
	file.open(filename);
	
  
	fillNodes();
	fillElements();//	fillElectrodes();
}
