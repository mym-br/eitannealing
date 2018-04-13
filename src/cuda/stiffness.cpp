#include "stiffness.h"

#include <cmath>

#include "model.h"

//NodeCoefficientsVector * nodeCoefficients;
NodeCoefficients **nodeCoef;

CondutanceCoefficientsVector * condutanceCoefficients;

struct Vector2D {
	numType x;
	numType y;
	Vector2D() {
		x = 0; y = 0;
	}
	Vector2D(numType px, numType py) {
		x = px; y = py;
	}
};

// old method ----------------------------------------------------------------------------
void insertNewCoefficient(NodeCoefficients **target, int node, int index, numType coefficient) {
	if (fabs(coefficient) < EPS) return;
	while (*target && (*target)->node < node) target = &(*target)->next;
	while (*target && (*target)->node == node)  {// special case, possible merge
		if ((*target)->condutanceIndex == index) { // Merge
			(*target)->coefficient += coefficient;
			return;
		}
		target = &(*target)->next;
	}
	// Insert node
	*target = new NodeCoefficients(node, index, coefficient, *target);
}
// also using old method --------------------------------------------------------END

void insertNewCoefficient(int nodeBase, int nodeDependency, int solIndex, numType coefficient) {
	if (MOD(coefficient) < EPS) { // zero
		return;
	}

	// condutance map
	CondutanceCoefficientsVector * ccv = &(condutanceCoefficients[solIndex]);
	for (int i = 0; i < ccv->size(); i++) {
		if ((*ccv)[i].row == nodeBase && (*ccv)[i].col == nodeDependency) { // merge
			(*ccv)[i].coefficient += coefficient;
			return;
		}
	}
	// if it got here, not found, add!
	CondutanceCoefficients temp_c(nodeBase, nodeDependency, coefficient);
	condutanceCoefficients[solIndex].push_back(temp_c);
}

void insertNewElementCoefficient(int nodeBase, int nodeDependency, int element, double coefficient) {
	// Coefficient is spread among the 3 nodes
	coefficient /= 3;
	insertNewCoefficient(nodeBase, nodeDependency, node2coefficient[elements[element].n1], coefficient);
	insertNewCoefficient(nodeBase, nodeDependency, node2coefficient[elements[element].n2], coefficient);
	insertNewCoefficient(nodeBase, nodeDependency, node2coefficient[elements[element].n3], coefficient);
}

// also using old method -----------------------------------------------------------
void insertNewElementCoefficient(NodeCoefficients **target, int node, int element, double coefficient) {
	// Coefficient is spread among the 3 nodes
	coefficient /= 3;
	insertNewCoefficient(target, node, node2coefficient[elements[element].n1], coefficient);
	insertNewCoefficient(target, node, node2coefficient[elements[element].n2], coefficient);
	insertNewCoefficient(target, node, node2coefficient[elements[element].n3], coefficient);
}
// also using old method --------------------------------------------------------END

void calcElementCoefficients(int element,
	double &c11, double &c22, double &c33,
	double &c12, double &c13, double &c23) {
	// Vectors of each node are the opposed edge rotated by 90 degrees
	// Notice in the end it doesn't matter if the nodes are clockwise or not!
	Vector2D
		v1(nodes[elements[element].n2].y - nodes[elements[element].n3].y,
		   nodes[elements[element].n3].x - nodes[elements[element].n2].x),
		v2(nodes[elements[element].n3].y - nodes[elements[element].n1].y,
		   nodes[elements[element].n1].x - nodes[elements[element].n3].x),
		v3(nodes[elements[element].n1].y - nodes[elements[element].n2].y,
		   nodes[elements[element].n2].x - nodes[elements[element].n1].x);

	double areaFactor = 2 * abs(v1.x * v2.y - v1.y * v2.x);

	c11 = (v1.x * v1.x + v1.y * v1.y) / areaFactor;
	c22 = (v2.x * v2.x + v2.y * v2.y) / areaFactor;
	c33 = (v3.x * v3.x + v3.y * v3.y) / areaFactor;

	c12 = (v1.x * v2.x + v1.y * v2.y) / areaFactor;
	c13 = (v3.x * v1.x + v3.y * v1.y) / areaFactor;
	c23 = (v2.x * v3.x + v2.y * v3.y) / areaFactor;
}

void calcElectrodeCoefficients(int electrode) {
	// FIXME: Approximate model, presumes 
	//	the outter edges are approximately paralel
	triangularEletrode thisElectrode = electrodes[electrode];

	std::vector<std::pair<int, int>>::iterator it;
	for (it = thisElectrode.nodesPairs.begin(); it != thisElectrode.nodesPairs.end(); it++) {
		// Vij: nj -> ni
		Vector2D vij(nodes[it->second].x - nodes[it->first].x,
					 nodes[it->second].y - nodes[it->first].y);
		// A = h / (2*l); B = l / (2*h) => A = 1 / (4*B) = 0.25/B
		double B = sqrt(vij.x * vij.x + vij.y * vij.y) / (2 * electrodeh);
		double A = 0.25 / B;

		double kiie = A + B;	// (h²+l²)/(2hl)
		double kjje = kiie;		// (h²+l²)/(2hl)
		double kije = -A;		// -h²/(2hl)
		double kbie = -B;		// -l²/(2hl)
		double kbje = -B;		// -l²/(2hl)
		double kbbe = 2 * B;	// 2*l²/(2hl) = l²/(hl) = l/h

		// multiplicar pela altura (totalheight) e somar aos acumuladores de coeficiente
		int index = node2coefficient[thisElectrode.baseNode];
		// Add to the base node...
		insertNewCoefficient(thisElectrode.baseNode, thisElectrode.baseNode,	index, kbbe * totalheight);
		insertNewCoefficient(thisElectrode.baseNode, it->first,					index, kbie * totalheight);
		insertNewCoefficient(thisElectrode.baseNode, it->second,				index, kbje * totalheight);

		// the i-th node...
		insertNewCoefficient(it->first, it->first,				index, kiie * totalheight);
		insertNewCoefficient(it->first, it->second,				index, kije * totalheight);
		insertNewCoefficient(it->first, thisElectrode.baseNode,	index, kbie * totalheight);

		// ... and the j-th node
		insertNewCoefficient(it->second, it->second,				index, kjje * totalheight);
		insertNewCoefficient(it->second, it->first,					index, kije * totalheight);
		insertNewCoefficient(it->second, thisElectrode.baseNode,	index, kbje * totalheight);

		// also using old method --------------------------------------------------------
		// Add to the base node...
		insertNewCoefficient(&nodeCoef[thisElectrode.baseNode], thisElectrode.baseNode, index, kbbe * totalheight);
		insertNewCoefficient(&nodeCoef[thisElectrode.baseNode], it->first, index, kbie * totalheight);
		insertNewCoefficient(&nodeCoef[thisElectrode.baseNode], it->second, index, kbje * totalheight);

		// the i-th node...
		insertNewCoefficient(&nodeCoef[it->first], it->first, index, kiie * totalheight);
		insertNewCoefficient(&nodeCoef[it->first], it->second, index, kije * totalheight);
		insertNewCoefficient(&nodeCoef[it->first], thisElectrode.baseNode, index, kbie * totalheight);

		// ... and the j-th node
		insertNewCoefficient(&nodeCoef[it->second], it->second, index, kjje * totalheight);
		insertNewCoefficient(&nodeCoef[it->second], it->first, index, kije * totalheight);
		insertNewCoefficient(&nodeCoef[it->second], thisElectrode.baseNode, index, kbje * totalheight);
		// also using old method --------------------------------------------------------END
	}
}

void buildNodeCoefficients() {
	// Init coefficients;
	nodeCoef = new NodeCoefficients *[nodes.size()];
	for (int i = 0; i<nodes.size(); i++) nodeCoef[i] = NULL;
	//nodeCoefficients = new NodeCoefficientsVector[nodes.size()];
	condutanceCoefficients = new CondutanceCoefficientsVector[numcoefficients];

	// Build electrodes
	for (int i = 0; i < electrodes.size(); i++) {
		calcElectrodeCoefficients(i);
	}
	// Now prepare the coefficients due to the elements
	double c11, c22, c33, c12, c23, c13;
	for (int i = 0; i < elements.size(); i++) {
		calcElementCoefficients(i, c11, c22, c33, c12, c13, c23);
		c11 *= totalheight; c22 *= totalheight; c33 *= totalheight;
		c12 *= totalheight; c23 *= totalheight; c13 *= totalheight;
		// Node 1
		insertNewElementCoefficient(elements[i].n1, elements[i].n1, i, c11);
		insertNewElementCoefficient(elements[i].n1, elements[i].n2, i, c12);
		insertNewElementCoefficient(elements[i].n1, elements[i].n3, i, c13);
		// Node 2
		insertNewElementCoefficient(elements[i].n2, elements[i].n2, i, c22);
		insertNewElementCoefficient(elements[i].n2, elements[i].n1, i, c12);
		insertNewElementCoefficient(elements[i].n2, elements[i].n3, i, c23);
		// Node 3
		insertNewElementCoefficient(elements[i].n3, elements[i].n3, i, c33);
		insertNewElementCoefficient(elements[i].n3, elements[i].n1, i, c13);
		insertNewElementCoefficient(elements[i].n3, elements[i].n2, i, c23);

		// also using old method -------------------------------------------------------
		// Node 1
		insertNewElementCoefficient(&nodeCoef[elements[i].n1], elements[i].n1, i, c11);
		insertNewElementCoefficient(&nodeCoef[elements[i].n1], elements[i].n2, i, c12);
		insertNewElementCoefficient(&nodeCoef[elements[i].n1], elements[i].n3, i, c13);
		// Node 2
		insertNewElementCoefficient(&nodeCoef[elements[i].n2], elements[i].n2, i, c22);
		insertNewElementCoefficient(&nodeCoef[elements[i].n2], elements[i].n1, i, c12);
		insertNewElementCoefficient(&nodeCoef[elements[i].n2], elements[i].n3, i, c23);
		// Node 3
		insertNewElementCoefficient(&nodeCoef[elements[i].n3], elements[i].n3, i, c33);
		insertNewElementCoefficient(&nodeCoef[elements[i].n3], elements[i].n1, i, c13);
		insertNewElementCoefficient(&nodeCoef[elements[i].n3], elements[i].n2, i, c23);
		// also using old method --------------------------------------------------------END
	}
}
