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
#include <boost/lambda/lambda.hpp>


std::ifstream file;

// Get the index of a node in the grid (no electrodes then) from
//	its x, y coordinates
inline int getElementIndex(int x, int y) { return 31+17*x+y; }

std::vector<node> nodes;
std::vector<triangularElement> elements;
std::vector<triangularEletrode> electrodes;



void addToElectrode(int n1, int n2, int n3)
{
	
	int tmp;
	// find lesser node (suposedly the base one
	if(n3 < n2) {
	  tmp = n2; n2 = n3; n3 = tmp;
	}
	if(n2 < n1) {
	  tmp = n1; n1 = n2; n2 = tmp;
	}
	using namespace boost::lambda;
	// Boost Lambda for locally-defined functor
	// STL has no support for pointer to member :(
		
	//int e = &(*electrodes.begin()) ->* &triangularEletrode::baseNode;
	
	//bool res = ((&(*_1) ->* &triangularEletrode::baseNode)==n1)(electrodes.begin());
	
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

	std::string aux;
	do {
	  std::getline(file,aux);
	} while(aux.compare("$ELM")!=0);
  
	
	int numElements;	// Notice this includes also electrode elements
	file >> numElements;
	
	elements.reserve(numElements);
	int i;
	for(i=0;i<numElements;i++) {
	  node n;
	  int id;
	  int n1, n2, n3;
	  file.ignore(256, ' ');
	  file.ignore(256, ' ');
	  file >> id;	  
	  file.ignore(256, ' ');
	  file.ignore(256, ' ');
	  file >> n1 >> n2 >> n3;
	  
	  switch(id) {
	      case 10000:	// electrode
		 addToElectrode(n1, n2, n3); 
		break;
	      case 2001:	// internal elements
		  
		break;
	      case 1001:	// external ring
		break;
		
		
	      
	      
		
	  }
	
	}
}
/*
triangularEletrode *electrodes;
int numElectrodes;

void fillElectrodes()
{
	numElectrodes = 32;
	electrodes = new triangularEletrode[numElectrodes];
	int i;
	// Left side (1st is the ground node
	electrodes[0].baseNode = numNodes-1;
	electrodes[0].n1 = getElementIndex(0, 0);
	electrodes[0].n2 = getElementIndex(0, 1);
	electrodes[0].n3 = getElementIndex(0, 2);
	for(i=1;i<8;i++) {
		electrodes[i].baseNode = i-1;
		electrodes[i].n1 = getElementIndex(0, 2*i);
		electrodes[i].n2 = getElementIndex(0, 2*i+1);
		electrodes[i].n3 = getElementIndex(0, 2*i+2);
	}
	// Bottom
	for(i=0;i<8;i++) {
		electrodes[i+8].baseNode = i+7;
		electrodes[i+8].n1 = getElementIndex(2*i, 16);
		electrodes[i+8].n2 = getElementIndex(2*i+1, 16);
		electrodes[i+8].n3 = getElementIndex(2*i+2, 16);
	}
	// Right
	for(i=0;i<8;i++) {
		electrodes[i+16].baseNode = i+15;
		electrodes[i+16].n1 = getElementIndex(16, 16-2*i);
		electrodes[i+16].n2 = getElementIndex(16, 15-2*i);
		electrodes[i+16].n3 = getElementIndex(16, 14-2*i);
	}
	// Top
	for(i=0;i<8;i++) {
		electrodes[i+24].baseNode = i+23;
		electrodes[i+24].n1 = getElementIndex(16-2*i, 0);
		electrodes[i+24].n2 = getElementIndex(15-2*i, 0);
		electrodes[i+24].n3 = getElementIndex(14-2*i, 0);
	}
}
*/
void initProblem(char *filename)
{
	file.open(filename);
	
  
	fillNodes();
	fillElements();//	fillElectrodes();
}
