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

std::ifstream file;

// Get the index of a node in the grid (no electrodes then) from
//	its x, y coordinates
inline int getElementIndex(int x, int y) { return 31+17*x+y; }

std::vector<node> nodes;
std::vector<triangularElement> elements;
std::vector<triangularEletrode> electrodes;


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
	  float x, y;
	  
	  file.ignore(256, ' ');
	  file >> x;
	  file >> y;
	}
}
/*
triangularElement *elements;
int numElements;
void fillElements() {

	numElements = 512;
	elements = new triangularElement[512];

	int elementCount = 0;
	int i, j;
	for(i=0;i<16;i++) { // Column
		for(j=0;j<16;j++) { // Row
			if((i+j)%2) { // odd elements
				elements[elementCount].n1 = getElementIndex(i,j);
				elements[elementCount].n2 = getElementIndex(i,j+1);
				elements[elementCount].n3 = getElementIndex(i+1,j+1);
				elements[elementCount++].condIndex = 1+8*(i/2) + j/2;
				elements[elementCount].n1 = getElementIndex(i,j);
				elements[elementCount].n2 = getElementIndex(i+1,j);
				elements[elementCount].n3 = getElementIndex(i+1,j+1);
				elements[elementCount++].condIndex = 1+8*(i/2) + j/2;
			} else { // Even elements
				elements[elementCount].n1 = getElementIndex(i,j);
				elements[elementCount].n2 = getElementIndex(i,j+1);
				elements[elementCount].n3 = getElementIndex(i+1,j);
				elements[elementCount++].condIndex = 1+8*(i/2) + j/2;
				elements[elementCount].n1 = getElementIndex(i+1,j);
				elements[elementCount].n2 = getElementIndex(i,j+1);
				elements[elementCount].n3 = getElementIndex(i+1,j+1);
				elements[elementCount++].condIndex = 1+8*(i/2) + j/2;
			}
		}
	}
}

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
//	fillElements();
//	fillElectrodes();
}
