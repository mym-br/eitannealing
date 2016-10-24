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

elementCoefficients calcElementCoefficients(const element &e)
{
	elementCoefficients c; // Hopefully copy-elision happens here
  
	// Vectors of each node are the opposed edge rotated by 90 degrees
	// Notice in the end it doesn't matter if the nodes are clockwise or not!
	Eigen::Vector2d
		va(nodes[e.b].y - nodes[e.c].y,
		   nodes[e.c].x - nodes[e.b].x),
		vb(nodes[e.c].y - nodes[e.a].y,
		   nodes[e.a].x - nodes[e.c].x),
		vc(nodes[e.a].y - nodes[e.b].y,
		   nodes[e.b].x - nodes[e.a].x);

	double areaFactor = 2*fabs(va.x()*vb.y()-va.y()*vb.x());

	c.ab = va.dot(vb)/areaFactor;
	c.ac = va.dot(vc)/areaFactor;
	c.bc = vb.dot(vc)/areaFactor;
	c.aa = -c.ab - c.ac;
	c.bb = -c.ab - c.bc;
	c.cc = -c.ac - c.bc;	

	return c;
}

void insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient)
{
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

void insertNewElementCoefficient(nodeCoefficients **target, int node, const triangularElement &e, double coefficient, const std::map<int, int> &coefficientMap)
{
	// Coefficient is spread among the 3 nodes
	coefficient /= 3;
	std::map<int, int>::const_iterator i;
	if((i=coefficientMap.find(e.a))!=coefficientMap.end())
	  insertNewCoefficient(target, node, i->second, coefficient);
	if((i=coefficientMap.find(e.b))!=coefficientMap.end())
	  insertNewCoefficient(target, node, i->second, coefficient);
	if((i=coefficientMap.find(e.c))!=coefficientMap.end())
	  insertNewCoefficient(target, node, i->second, coefficient);
}

void calcGenericElectrodeCoefficients(int electrode, const std::vector<genericEletrode> &electrodes)
{
	// FIXME: Approximate model, presumes 
	//	the outter edges are approximately paralel
	genericEletrode thisElectrode = electrodes[electrode];
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
		calcGenericElectrodeCoefficients(i, gelectrodes);
	}
	// Now prepare the coefficients due to the elements
	for(const element e: elements) {		
		elementCoefficients c = calcElementCoefficients(e);
		c *= totalheight;
		// Node 1
		insertNewElementCoefficient(&nodeCoef[e.a], e.a, e, c.aa, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.a], e.b, e, c.ab, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.a], e.c, e, c.ac, node2coefficient);
		// Node 2
		insertNewElementCoefficient(&nodeCoef[e.b], e.b, e, c.bb, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.b], e.a, e, c.ab, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.b], e.c, e, c.bc, node2coefficient);
		// Node 3
		insertNewElementCoefficient(&nodeCoef[e.c], e.c, e, c.cc, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.c], e.a, e, c.ac, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.c], e.b, e, c.bc, node2coefficient);
	}
}
