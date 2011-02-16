/*
 * initproblem.cpp
 *
 * 	Prepare the problem instance
 *
 * 		A 8x8 grid, each cell composed of 8 triangular elements
 * 			The grid is surrounded by 32 "electodes"
 *
 *        *   *
 *       /|\ /|\
 *		+-+-+-+-+...
 *     /|/|\|/|\|
 *    *-+-+-+-+-+
 *     \|\|/|\|/|
 *      +-+-+-+-+
 *     /|/|\|/|\|
 *    *-+-+-+-+-+
 *     \|\|/|\|/|
 *      +-+-+-+-+
 *		.        .
 *		.         .
 *      .          .
 *
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#include "problemdescription.h"

// Get the index of a node in the grid (no electrodes then) from
//	its x, y coordinates
inline int getElementIndex(int x, int y) { return 31+17*x+y; }

node *nodes;
int numNodes;

void fillNodes() {
	numNodes = 32 + (8*2+1)*(8*2+1);
	nodes = new node[numNodes];
	int i,j;
	// Electrodes 1st:
	// Left side (notice the top left electrode is the ground node
	//	and as such is added only at the end
	for(i=1;i<8;i++) { nodes[i-1].x = -1; nodes[i-1].y = 15 - 2*i; }
	// Bottom
	for(i=0;i<8;i++) { nodes[i+7].x = 1+2*i; nodes[i+7].y = -1; }
	// Right
	for(i=0;i<8;i++) { nodes[i+15].x = 17; nodes[i+15].y = 1+2*i; }
	// Top
	for(i=0;i<8;i++) { nodes[i+23].x = 15 - 2*i; nodes[i+23].y = 17; }
	// Now main nodes
	for(i=0;i<17;i++) { // Column
		for(j=0;j<17;j++) { // Row
			int index = getElementIndex(i,j);
			nodes[index].x = i;
			nodes[index].y = 16 - j;
		}
	}
	// and finally, the top left electrode (Ground node)
	nodes[numNodes-1].x = -1.0;
	nodes[numNodes-1].y = 15.0;
}

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

void initProblem()
{
	fillNodes();
	fillElements();
	fillElectrodes();
}
