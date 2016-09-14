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

// removido para usar novos eletrodos
//void calcELectrodeCoefficients(int electrode,
//							   double &cbb, double &c11, double &c22, double &c33,
//							   double &cb1, double &cb2, double &cb3,
//							   double &c12, double &c23)
//{
//	// FIXME: Approximate model, presumes 
//	//	the outter edges are approximately paralel
//	Eigen::Vector2d	// V3: n2 -> n1 v3: n2 -> n3
//		v1(	nodes[electrodes[electrode].n1].x - nodes[electrodes[electrode].n2].x,
//		nodes[electrodes[electrode].n1].y - nodes[electrodes[electrode].n2].y),
//		v3(	nodes[electrodes[electrode].n3].x - nodes[electrodes[electrode].n2].x,
//		nodes[electrodes[electrode].n3].y - nodes[electrodes[electrode].n2].y);
//
//	cb1 = -v1.norm()/(2*electrodeh);
//	c12 = .25/cb1;
//	c11 = -(cb1 + c12);
//
//	cb3 = -v3.norm()/(2*electrodeh);
//	c23 = .25/cb3;
//	c33 = -(cb3 + c23);
//
//	cb2 = cb3 + cb1;
//	c22 = -(cb2 + c12 + c23);
//
//	cbb = -(cb1 + cb2 + cb3);  
//}

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

	//c11 = v1.squaredNorm()/areaFactor;
	//c22 = v2.squaredNorm()/areaFactor;
	//c33 = v3.squaredNorm()/areaFactor;
	c12 = v1.dot(v2)/areaFactor;
	c13 = v1.dot(v3)/areaFactor;
	c23 = v2.dot(v3)/areaFactor;
	c11 = -c12 - c13;
	c22 = -c12 - c23;
	c33 = -c13 - c23;	
}

void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient)
{
	//if(fabs(coefficient)<_COEFFICIENT_EPSLON_) return;
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

void insertNewElementCoefficient(nodeCoefficients **target, int node, int element, double coefficient)
{
	// Coefficient is spread among the 3 nodes
	coefficient /= 3;
	insertNewCoefficient(target, node, node2coefficient[elements[element].n1], coefficient);
	insertNewCoefficient(target, node, node2coefficient[elements[element].n2], coefficient);
	insertNewCoefficient(target, node, node2coefficient[elements[element].n3], coefficient);
}

void calcGenericElectrodeCoefficients(int electrode)
{
	// FIXME: Approximate model, presumes 
	//	the outter edges are approximately paralel
	genericEletrode thisElectrode = gelectrodes[electrode];
	for (std::vector<std::pair<int, int>>::iterator it = thisElectrode.nodesPairs.begin() ; it != thisElectrode.nodesPairs.end(); ++it) {
		// Vij: nj -> ni
		Eigen::Vector2d	vij(nodes[it->second].x - nodes[it->first].x,
							nodes[it->second].y - nodes[it->first].y);
		// A = h / (2*l); B = l / (2*h) => A = 1 / (4*B) = 0.25/B
		double B = vij.norm()/(2*electrodeh);
		double A = 0.25/B;

		double kiie = A+B;	// (h�+l�)/(2hl)
		double kjje = kiie; // (h�+l�)/(2hl)
		double kije = -A;	// -h�/(2hl)
		double kbie = -B;	// -l�/(2hl)
		double kbje = -B;	// -l�/(2hl)
		double kbbe = 2*B;	// 2*l�/(2hl) = l�/(hl) = l/h

		// multiplicar pela altura (totalheight) e somar aos acumuladores de coeficiente
		int index = node2coefficient[thisElectrode.baseNode]; 
		// Add to the base node...
		insertNewCoefficient(&nodeCoef[thisElectrode.baseNode], thisElectrode.baseNode, index, kbbe * totalheight);
		insertNewCoefficient(&nodeCoef[thisElectrode.baseNode], it->first, index,kbie * totalheight);
		insertNewCoefficient(&nodeCoef[thisElectrode.baseNode], it->second, index, kbje * totalheight);

		// the i-th node...
		insertNewCoefficient(&nodeCoef[it->first], it->first, index, kiie * totalheight);
		insertNewCoefficient(&nodeCoef[it->first], it->second, index, kije * totalheight);
		insertNewCoefficient(&nodeCoef[it->first], thisElectrode.baseNode, index, kbie * totalheight);

		// ... and the j-th node
		insertNewCoefficient(&nodeCoef[it->second], it->second, index, kjje * totalheight);
		insertNewCoefficient(&nodeCoef[it->second], it->first, index, kije * totalheight);
		insertNewCoefficient(&nodeCoef[it->second], thisElectrode.baseNode, index, kbje * totalheight);
	} 
}

void buildNodeCoefficients()
{
	int i;
	// Init coefficients;
	nodeCoef = new nodeCoefficients *[nodes.size()];
	for(int i = 0; i<nodes.size(); i++) nodeCoef[i] = NULL;
	// Build electrodes
	for(i = 0; i < gelectrodes.size(); i++) {
		calcGenericElectrodeCoefficients(i);
	}
	// Now prepare the coefficients due to the elements
	double c11, c22, c33, c12, c23, c13;
	for(i = 0; i < elements.size(); i++) {
		calcElementCoefficients(i, c11, c22, c33, c12, c13, c23);
		c11 *= totalheight; c22 *= totalheight; c33 *= totalheight;
		c12 *= totalheight; c23 *= totalheight; c13 *= totalheight;
		// Node 1
		insertNewElementCoefficient(&nodeCoef[elements[i].n1],
			elements[i].n1, i, c11);
		insertNewElementCoefficient(&nodeCoef[elements[i].n1],
			elements[i].n2, i, c12);
		insertNewElementCoefficient(&nodeCoef[elements[i].n1],
			elements[i].n3, i, c13);
		// Node 2
		insertNewElementCoefficient(&nodeCoef[elements[i].n2],
			elements[i].n2, i, c22);
		insertNewElementCoefficient(&nodeCoef[elements[i].n2],
			elements[i].n1, i, c12);
		insertNewElementCoefficient(&nodeCoef[elements[i].n2],
			elements[i].n3, i, c23);
		// Node 3
		insertNewElementCoefficient(&nodeCoef[elements[i].n3],
			elements[i].n3, i, c33);
		insertNewElementCoefficient(&nodeCoef[elements[i].n3],
			elements[i].n1, i, c13);
		insertNewElementCoefficient(&nodeCoef[elements[i].n3],
			elements[i].n2, i, c23);
	}
}
