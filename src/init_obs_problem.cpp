/*
 * init_obs_problem.cpp
 *
 *  Created on: Oct 10, 2010
 *      Author: thiago
 *
 * 	Prepare the problem instance for the observations
 *
 * 		A 8x8 grid, each cell composed of 32 triangular elements
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
#include "nodecoefficients.h"
#include "solver.h"

namespace obs {

// Get the index of a node in the grid (no electrodes then) from
//	its x, y coordinates
static inline int getElementIndex(int x, int y) { return 31+32*5 + 33*x+y; }

node *nodes;
int numNodes;

static void fillNodes() {
	numNodes = 32 + // electrode entry
			   32*5 + // electrode contact
			(8*4+1)*(8*4+1);
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

	// Electrode contact nodes
	// Left
	for(i=0;i<8;i++) {
		 for(j=0;j<5;j++) {
			 nodes[31+5*i+j].x = -0.5;
			 nodes[31+5*i+j].y = 16 - 2*i - 0.5*j;
		 }
	}

	// Bottom
	for(i=0;i<8;i++) {
		for(j=0;j<5;j++) {
			nodes[31+40+5*i+j].x = 2*i+0.5*j;
			nodes[31+40+5*i+j].y = -0.5;
		}
	}

	// Right
	for(i=0;i<8;i++) {
		for(j=0;j<5;j++) {
			nodes[31+80+5*i+j].x = 16.5;
			nodes[31+80+5*i+j].y = 2*i+0.5*j;
		}
	}
	// Top
	for(i=0;i<8;i++) {
		for(j=0;j<5;j++) {
			nodes[31+120+5*i+j].x = 16-2*i-0.5*j;
			nodes[31+120+5*i+j].y = 16.5;
		}
	}


	// Now main nodes
	for(i=0;i<33;i++) { // Column
		for(j=0;j<33;j++) { // Row
			int index = getElementIndex(i,j);
			nodes[index].x = (float)i/2;
			nodes[index].y = 16 - (float)j/2;
		}
	}
	// and finally, the top left electrode (Ground node)
	nodes[numNodes-1].x = -1.0;
	nodes[numNodes-1].y = 15.0;
}

static triangularElement *elements;
static int numElements;
static void fillElements() {

	numElements = 4*512+8*32;
	elements = new triangularElement[numElements];

	int elementCount = 0;
	int i, j;
	// Main body
	for(i=0;i<32;i++) { // Column
		for(j=0;j<32;j++) { // Row
			if((i+j)%2) { // odd elements
				elements[elementCount].n1 = getElementIndex(i,j);
				elements[elementCount].n2 = getElementIndex(i+1,j);
				elements[elementCount].n3 = getElementIndex(i+1,j+1);
				elements[elementCount++].condIndex = 1+8*(i/4) + j/4;
				elements[elementCount].n1 = getElementIndex(i,j);
				elements[elementCount].n2 = getElementIndex(i,j+1);
				elements[elementCount].n3 = getElementIndex(i+1,j+1);
				elements[elementCount++].condIndex = 1+8*(i/4) + j/4;
			} else { // Even elements
				elements[elementCount].n1 = getElementIndex(i,j);
				elements[elementCount].n2 = getElementIndex(i,j+1);
				elements[elementCount].n3 = getElementIndex(i+1,j);
				elements[elementCount++].condIndex = 1+8*(i/4) + j/4;
				elements[elementCount].n1 = getElementIndex(i+1,j);
				elements[elementCount].n2 = getElementIndex(i,j+1);
				elements[elementCount].n3 = getElementIndex(i+1,j+1);
				elements[elementCount++].condIndex = 1+8*(i/4) + j/4;
			}
		}
	}
	// Electrode body
	// Left column
	for(i=0;i<8;i++) {
		for(j=0;j<4;j+=2) {
			elements[elementCount].n1 = 31+5*i+j;
			elements[elementCount].n2 = getElementIndex(0, j+4*i);
			elements[elementCount].n3 = getElementIndex(0, j+1+4*i);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+5*i+j;
			elements[elementCount].n2 = 31+5*i+j+1;
			elements[elementCount].n3 = getElementIndex(0, j+1+4*i);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+5*i+j+1;
			elements[elementCount].n2 = 31+5*i+j+2;
			elements[elementCount].n3 = getElementIndex(0, j+1+4*i);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+5*i+j+2;
			elements[elementCount].n2 = getElementIndex(0, j+1+4*i);
			elements[elementCount].n3 = getElementIndex(0, j+2+4*i);
			elements[elementCount++].condIndex = 0;
		}
	}
	// bottom row
	for(i=0;i<8;i++) {
		for(j=0;j<4;j+=2) {
			elements[elementCount].n1 = 31+40+5*i+j;
			elements[elementCount].n2 = getElementIndex(j+4*i, 32);
			elements[elementCount].n3 = getElementIndex(j+1+4*i,32);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+40+5*i+j;
			elements[elementCount].n2 = 31+40+5*i+j+1;
			elements[elementCount].n3 = getElementIndex(j+1+4*i, 32);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+40+5*i+j+1;
			elements[elementCount].n2 = 31+40+5*i+j+2;
			elements[elementCount].n3 = getElementIndex(j+1+4*i, 32);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+40+5*i+j+2;
			elements[elementCount].n2 = getElementIndex(j+1+4*i, 32);
			elements[elementCount].n3 = getElementIndex(j+2+4*i, 32);
			elements[elementCount++].condIndex = 0;
		}
	}
	// right column
	for(i=0;i<8;i++) {
		for(j=0;j<4;j+=2) {
			elements[elementCount].n1 = 31+80+5*i+j;
			elements[elementCount].n2 = getElementIndex(32, 32-j-4*i);
			elements[elementCount].n3 = getElementIndex(32, 31-j-4*i);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+80+5*i+j;
			elements[elementCount].n2 = 31+80+5*i+j+1;
			elements[elementCount].n3 = getElementIndex(32, 31-j-4*i);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+80+5*i+j+1;
			elements[elementCount].n2 = 31+80+5*i+j+2;
			elements[elementCount].n3 = getElementIndex(32, 31-j-4*i);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+80+5*i+j+2;
			elements[elementCount].n2 = getElementIndex(32, 31-j-4*i);
			elements[elementCount].n3 = getElementIndex(32, 30-j-4*i);
			elements[elementCount++].condIndex = 0;
		}
	}
	// top row
	for(i=0;i<8;i++) {
		for(j=0;j<4;j+=2) {
			elements[elementCount].n1 = 31+120+5*i+j;
			elements[elementCount].n2 = getElementIndex(32-j-4*i, 0);
			elements[elementCount].n3 = getElementIndex(31-j-4*i, 0);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+120+5*i+j;
			elements[elementCount].n2 = 31+120+5*i+j+1;
			elements[elementCount].n3 = getElementIndex(31-j-4*i, 0);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+120+5*i+j+1;
			elements[elementCount].n2 = 31+120+5*i+j+2;
			elements[elementCount].n3 = getElementIndex(31-j-4*i, 0);
			elements[elementCount++].condIndex = 0;
			elements[elementCount].n1 = 31+120+5*i+j+2;
			elements[elementCount].n2 = getElementIndex(31-j-4*i, 0);
			elements[elementCount].n3 = getElementIndex(30-j-4*i, 0);
			elements[elementCount++].condIndex = 0;
		}
	}
}

static nodeCoefficients **nodeCoef;

int numElectrodes = 32;

void buildNodeCoefficients();

void initObsProblem()
{
	fillNodes();
	fillElements();
	buildNodeCoefficients();
	// partial Cleanup
	delete[] elements;
	delete[] nodes;

	//fillElectrodes();
}

matrix *buildObsProblemMatrix(float *coeff)
{
	matrix *res;
	assembleProblemMatrix(coeff, &res, numNodes, nodeCoef);
	return res;
}

static void calcElementCoefficients(int element,
		double &c11, double &c22, double &c33,
		double &c12, double &c13, double &c23)
{
	// Vectors of each node are the opposed edge rotated by 90 degrees
	// Notice in the end it doesn't matter if the nodes are clockwise or not!
	Eigen::Vector2d
		v1(	nodes[elements[element].n2].y - nodes[elements[element].n3].y,
			nodes[elements[element].n3].x - nodes[elements[element].n2].x),
		v2(	nodes[elements[element].n3].y - nodes[elements[element].n1].y,
			nodes[elements[element].n1].x - nodes[elements[element].n3].x),
		v3(	nodes[elements[element].n1].y - nodes[elements[element].n2].y,
			nodes[elements[element].n2].x - nodes[elements[element].n1].x);

	double areaFactor = 2*fabs(v1.x()*v2.y()-v1.y()*v2.x());

	c11 = v1.squaredNorm()/areaFactor;
	c22 = v2.squaredNorm()/areaFactor;
	c33 = v3.squaredNorm()/areaFactor;
	c12 = v1.dot(v2)/areaFactor;
	c13 = v1.dot(v3)/areaFactor;
	c23 = v2.dot(v3)/areaFactor;
}

#define _COEFFICIENT_EPSLON_ (0.00001)

static void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient)
{
	if(fabs(coefficient)<_COEFFICIENT_EPSLON_) return;
	while(*target && (*target)->node < node) target = &(*target)->next;
	while(*target && (*target)->node == node)  {// special case, possible merge
		if((*target)->condIndex == index) { // Merge
			(*target)->coefficient += coefficient;
			return;
		}
		target = &(*target)->next;
	}
	// Insert node
	*target = new nodeCoefficients(node, index, coefficient, *target);
}

void buildNodeCoefficients()
{
	int i;
	// Init coefficients;
	nodeCoef = new nodeCoefficients *[numNodes-1];
	for(int i = 0; i<numNodes-1; i++) nodeCoef[i] = NULL;
	// Build electrodes
	// First for the ground node
	// 1
	insertNewCoefficient(&nodeCoef[31], 31, 0, 1.0);
	insertNewCoefficient(&nodeCoef[31], 32, 0, -0.5);
	// 2
	insertNewCoefficient(&nodeCoef[32], 32, 0, 2.0);
	insertNewCoefficient(&nodeCoef[32], 31, 0, -0.5);
	insertNewCoefficient(&nodeCoef[32], 33, 0, -0.5);
	// 3
	insertNewCoefficient(&nodeCoef[33], 33, 0, 2.0);
	insertNewCoefficient(&nodeCoef[33], 32, 0, -0.5);
	insertNewCoefficient(&nodeCoef[33], 34, 0, -0.5);
	// 4
	insertNewCoefficient(&nodeCoef[34], 34, 0, 2.0);
	insertNewCoefficient(&nodeCoef[34], 33, 0, -0.5);
	insertNewCoefficient(&nodeCoef[34], 35, 0, -0.5);
	// 5
	insertNewCoefficient(&nodeCoef[35], 35, 0, 1.0);
	insertNewCoefficient(&nodeCoef[35], 34, 0, -0.5);

	for(i = 0; i < numElectrodes-1; i++) {
		// Add to the base node...
		insertNewCoefficient(&nodeCoef[i], i, 0, 4);
		insertNewCoefficient(&nodeCoef[i], 5*i+36, 0, -0.5);
		insertNewCoefficient(&nodeCoef[i], 5*i+37, 0, -1.0);
		insertNewCoefficient(&nodeCoef[i], 5*i+38, 0, -1.0);
		insertNewCoefficient(&nodeCoef[i], 5*i+39, 0, -1.0);
		insertNewCoefficient(&nodeCoef[i], 5*i+40, 0, -0.5);
		// 1
		insertNewCoefficient(&nodeCoef[5*i+36], 5*i+36, 0, 1.0);
		insertNewCoefficient(&nodeCoef[5*i+36], i,      0, -0.5);
		insertNewCoefficient(&nodeCoef[5*i+36], 5*i+37, 0, -0.5);
		// 2
		insertNewCoefficient(&nodeCoef[5*i+37], 5*i+37, 0, 2.0);
		insertNewCoefficient(&nodeCoef[5*i+37], i,      0, -1.0);
		insertNewCoefficient(&nodeCoef[5*i+37], 5*i+36, 0, -0.5);
		insertNewCoefficient(&nodeCoef[5*i+37], 5*i+38, 0, -0.5);
		// 3
		insertNewCoefficient(&nodeCoef[5*i+38], 5*i+38, 0, 2.0);
		insertNewCoefficient(&nodeCoef[5*i+38], i,      0, -1.0);
		insertNewCoefficient(&nodeCoef[5*i+38], 5*i+37, 0, -0.5);
		insertNewCoefficient(&nodeCoef[5*i+38], 5*i+39, 0, -0.5);
		// 4
		insertNewCoefficient(&nodeCoef[5*i+39], 5*i+39, 0, 2.0);
		insertNewCoefficient(&nodeCoef[5*i+39], i,      0, -1.5);
		insertNewCoefficient(&nodeCoef[5*i+39], 5*i+38, 0, -0.5);
		insertNewCoefficient(&nodeCoef[5*i+39], 5*i+40, 0, -0.5);
		// 5
		insertNewCoefficient(&nodeCoef[5*i+40], 5*i+40, 0, 1.0);
		insertNewCoefficient(&nodeCoef[5*i+40], i,      0, -0.5);
		insertNewCoefficient(&nodeCoef[5*i+40], 5*i+39, 0, -0.5);
	}
	// Now prepare the coefficients due to the elements
	double c11, c22, c33, c12, c23, c13;
	for(i = 0; i < numElements; i++) {
		calcElementCoefficients(i, c11, c22, c33, c12, c13, c23);
		// Node 1
		insertNewCoefficient(&nodeCoef[elements[i].n1],
				elements[i].n1, elements[i].condIndex, c11);
		insertNewCoefficient(&nodeCoef[elements[i].n1],
				elements[i].n2, elements[i].condIndex, c12);
		insertNewCoefficient(&nodeCoef[elements[i].n1],
				elements[i].n3, elements[i].condIndex, c13);
		// Node 2
		insertNewCoefficient(&nodeCoef[elements[i].n2],
				elements[i].n2, elements[i].condIndex, c22);
		insertNewCoefficient(&nodeCoef[elements[i].n2],
				elements[i].n1, elements[i].condIndex, c12);
		insertNewCoefficient(&nodeCoef[elements[i].n2],
				elements[i].n3, elements[i].condIndex, c23);
		// Node 3
		insertNewCoefficient(&nodeCoef[elements[i].n3],
				elements[i].n3, elements[i].condIndex, c33);
		insertNewCoefficient(&nodeCoef[elements[i].n3],
				elements[i].n1, elements[i].condIndex, c13);
		insertNewCoefficient(&nodeCoef[elements[i].n3],
				elements[i].n2, elements[i].condIndex, c23);
	}
}

}

