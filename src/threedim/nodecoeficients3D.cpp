#include "problem3D.h"
#include <Eigen/Dense>
#include <iostream>

void problem3D::buildNodeCoefficients()
{
	// Init coefficients;
	nodeCoef = new nodeCoefficients *[nodes.size()];
	for (int i = 0; i<nodes.size(); i++) nodeCoef[i] = NULL;
	//// Build electrodes
	//auto insertElectrodeCoefficient = [this](int node_a, int node_b, int condIndex, double cab) {
	//	insertNewCoefficient(&nodeCoef[node_a], node_b, condIndex, cab);
	//};
	//for (const genericEletrode &e : gelectrodes) calcAndInsertGenericElectrodeCoefficients(
	//	e, nodes, electrodeh, totalheight, node2coefficient, insertElectrodeCoefficient);

		// Now prepare the coefficients due to the elements
	for (const tetrahedralElement e : elements) {
		elementCoefficients c = calcElementCoefficients(e);
		c *= totalheight;
		// Node 1
		insertNewElementCoefficient(&nodeCoef[e.a], e.a, e, c.aa, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.a], e.b, e, c.ab, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.a], e.c, e, c.ac, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.d], e.d, e, c.ad, node2coefficient);
		// Node 2
		insertNewElementCoefficient(&nodeCoef[e.b], e.b, e, c.bb, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.b], e.a, e, c.ab, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.b], e.c, e, c.bc, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.b], e.d, e, c.bd, node2coefficient);
		// Node 3
		insertNewElementCoefficient(&nodeCoef[e.c], e.c, e, c.cc, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.c], e.a, e, c.ac, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.c], e.b, e, c.bc, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.c], e.d, e, c.cd, node2coefficient);
	}
}

problem3D::elementCoefficients problem3D::calcElementCoefficients(const tetrahedralElement &e)
{
	elementCoefficients c; // Hopefully copy-elision happens here

	// Edges vectors
	Eigen::Vector3d
		vjl(nodes[e.d].x - nodes[e.b].x, nodes[e.d].y - nodes[e.b].y, nodes[e.d].z - nodes[e.b].z),
		vjk(nodes[e.c].x - nodes[e.b].x, nodes[e.c].y - nodes[e.b].y, nodes[e.c].z - nodes[e.b].z),
		vik(nodes[e.c].x - nodes[e.a].x, nodes[e.c].y - nodes[e.a].y, nodes[e.c].z - nodes[e.a].z),
		vil(nodes[e.d].x - nodes[e.a].x, nodes[e.d].y - nodes[e.a].y, nodes[e.d].z - nodes[e.a].z),
		vij(nodes[e.b].x - nodes[e.a].x, nodes[e.b].y - nodes[e.a].y, nodes[e.b].z - nodes[e.a].z);

	// Vectors of each node are the opposed edge rotated to coincide with normal of opposed face
	// Check if order of vertexes impacts the result, it appears not
	Eigen::Vector3d va(vjl.cross(vjk)), vb(vik.cross(vil)), vc(vil.cross(vij)), vd(vij.cross(vik));

	double areaFactor = 6 * fabs((vij.cross(vik)).dot(vil));
	c.ab = va.dot(vb) / areaFactor;
	c.ac = va.dot(vc) / areaFactor;
	c.ad = va.dot(vd) / areaFactor;
	c.bc = vb.dot(vc) / areaFactor;
	c.bd = vb.dot(vd) / areaFactor;
	c.cd = vc.dot(vd) / areaFactor;
	c.aa = -c.ab - c.ac - c.ad;
	c.bb = -c.ab - c.bc - c.bd;
	c.cc = -c.ac - c.bc - c.cd;
	c.dd = -c.ad - c.bd - c.cd;

	std::cout << c.aa << " " << c.ab << " " << c.ac << " " << c.ad << std::endl;
	std::cout << c.ab << " " << c.bb << " " << c.bc << " " << c.bd << std::endl;
	std::cout << c.ac << " " << c.bc << " " << c.cc << " " << c.cd << std::endl;
	std::cout << c.ad << " " << c.bd << " " << c.cd << " " << c.dd << std::endl;

	return c;
}

void problem3D::insertNewElementCoefficient(nodeCoefficients **target, int node, const tetrahedralElement &e, double coefficient, const std::map<int, int> &coefficientMap)
{
	// Coefficient is spread among the 4 nodes
	coefficient /= 4;
	std::map<int, int>::const_iterator i;
	if ((i = coefficientMap.find(e.a)) != coefficientMap.end())
		insertNewCoefficient(target, node, i->second, coefficient);
	if ((i = coefficientMap.find(e.b)) != coefficientMap.end())
		insertNewCoefficient(target, node, i->second, coefficient);
	if ((i = coefficientMap.find(e.c)) != coefficientMap.end())
		insertNewCoefficient(target, node, i->second, coefficient);
	if ((i = coefficientMap.find(e.d)) != coefficientMap.end())
		insertNewCoefficient(target, node, i->second, coefficient);
}

void problem3D::insertNewCoefficient(nodeCoefficients **target, int node, int index, double coefficient)
{
	while (*target && (*target)->node < node) target = &(*target)->next;
	while (*target && (*target)->node == node)  {// special case, possible merge
		if ((*target)->condIndex == index) { // Merge
			(*target)->coefficient += coefficient;
			return;
		}
		target = &(*target)->next;
	}
	// Insert node
	*target = new nodeCoefficients(node, index, coefficient, *target);
}