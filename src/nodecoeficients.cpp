/*
 * nodecoeficients.cpp
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#include "nodecoefficients.h"
#include "problemdescription.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

nodeCoefficients **nodeCoef;

#define _COEFFICIENT_EPSLON_ (0.00001)

void calcElementCoefficients(int element,
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

void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient)
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
	nodeCoef = new nodeCoefficients *[nodes.size()-1];
	for(int i = 0; i<nodes.size()-1; i++) nodeCoef[i] = NULL;
	// Build electrodes
	for(i = 0; i < electrodes.size(); i++) {
		if(i==0) { // Ground node, coefficients of the base node are zeroed
		    // the 1st node...
		    insertNewCoefficient(&nodeCoef[electrodes[i].n1],
						    electrodes[i].n1, 0.0, 1.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n1],
						    electrodes[i].n2, 0.0,-0.5);
		    // the 2nd...
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].n2, 0.0, 2.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].n1, 0.0,-0.5);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].n3, 0.0,-0.5);
		    // ... and the 3rd
		    insertNewCoefficient(&nodeCoef[electrodes[i].n3],
						    electrodes[i].n3, 0.0, 1.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n3],
						    electrodes[i].n2, 0.0,-0.5);
		} else {
		    // Add to the base node...
		    insertNewCoefficient(&nodeCoef[electrodes[i].baseNode],
				    electrodes[i].baseNode, 0.0, 2.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].baseNode],
						    electrodes[i].n1, 0.0,-0.5);
		    insertNewCoefficient(&nodeCoef[electrodes[i].baseNode],
						    electrodes[i].n2, 0.0,-1.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].baseNode],
						    electrodes[i].n3, 0.0,-0.5);
		    // And to the 1st node...
		    insertNewCoefficient(&nodeCoef[electrodes[i].n1],
						    electrodes[i].n1, 0.0, 1.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n1],
						    electrodes[i].baseNode, 0.0,-0.5);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n1],
						    electrodes[i].n2, 0.0,-0.5);
		    // the 2nd...
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].n2, 0.0, 2.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].baseNode, 0.0,-1.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].n1, 0.0,-0.5);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n2],
						    electrodes[i].n3, 0.0,-0.5);
		    // ... and the 3rd
		    insertNewCoefficient(&nodeCoef[electrodes[i].n3],
						    electrodes[i].n3, 0.0, 1.0);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n3],
						    electrodes[i].baseNode, 0.0,-0.5);
		    insertNewCoefficient(&nodeCoef[electrodes[i].n3],
						    electrodes[i].n2, 0.0,-0.5);
		}
	}
	// Now prepare the coefficients due to the elements
	double c11, c22, c33, c12, c23, c13;
	for(i = 0; i < elements.size(); i++) {
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
