#include "problem3D.h"
#include <Eigen/Dense>
#include <iostream>

void problem3D::buildNodeCoefficients()
{
	// Init coefficients;
	nodeCoef = new nodeCoefficients *[nodes.size()];
	for (int i = 0; i<nodes.size(); i++) nodeCoef[i] = NULL;
	// Build electrodes
	auto insertElectrodeCoefficient = [this](int node_a, int node_b, int condIndex, double cab) {
		insertNewCoefficient(&nodeCoef[node_a], node_b, condIndex, cab);
	};
	for (const std::pair<int,genericEletrode> &e : gelectrodes) calcAndInsertGenericElectrodeCoefficients(
		e.second, nodes, electrodeh, totalheight, node2coefficient, insertElectrodeCoefficient);

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
		// Node 4
		insertNewElementCoefficient(&nodeCoef[e.d], e.d, e, c.dd, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.d], e.a, e, c.ad, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.d], e.b, e, c.bd, node2coefficient);
		insertNewElementCoefficient(&nodeCoef[e.d], e.c, e, c.cd, node2coefficient);
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

	// Vectors of each node are the opposed edge rotated to coincide with normal of opposed face,
	// order of vertexes does not impact the result
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

	//std::cout << c.aa << " " << c.ab << " " << c.ac << " " << c.ad << std::endl;
	//std::cout << c.ab << " " << c.bb << " " << c.bc << " " << c.bd << std::endl;
	//std::cout << c.ac << " " << c.bc << " " << c.cc << " " << c.cd << std::endl;
	//std::cout << c.ad << " " << c.bd << " " << c.cd << " " << c.dd << std::endl;

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

void problem3D::calcAndInsertGenericElectrodeCoefficients(const genericEletrode &e, const std::vector<node> &nodes, double electrodeh, double totalheight,
	const std::map<int, int> &coefficientMap,
	const std::function<void(int, int, int, double)> &insert)
{
	for (auto p : e.nodesTriangles) {
		// Get triangle base parameters
		double hz = electrodeh;
		node vi(nodes[p.a]), vj(nodes[p.b]), vk(nodes[p.c]);
		Eigen::Vector3d	vji(vj.x - vi.x, vj.y - vi.y, vj.z - vi.z), vkj(vk.x - vj.x, vk.y - vj.y, vk.z - vj.z), vki(vk.x - vi.x, vk.y - vi.y, vk.z - vi.z);
		double lk(vji.norm()), li(vkj.norm()), lj(vki.norm());
		double S = (li + lj + lk) / 2;
		double xji = lk;
		double hy = (2 * sqrt(S*(S - li)*(S - lj)*(S - lk))) / lk;
		double xki = sqrt(lj*lj - hy*hy);
		double xkj = sqrt(li*li - hy*hy);
		xki = xkj > xji ? -xki : xki;
		xkj = xki > xji ? xkj : -xkj;
		double areaFactor = 6 * hy*hz*xji;
		
		//std::cout << "vi = (" << vi.x << ", " << vi.y << ", " << vi.z << ")." << std::endl;
		//std::cout << "vj = (" << vj.x << ", " << vj.y << ", " << vj.z << ")." << std::endl;
		//std::cout << "vk = (" << vk.x << ", " << vk.y << ", " << vk.z << ")." << std::endl;
		//std::cout << "areaFactor = " << areaFactor << std::endl;

		double A, B, C, D, E, F;
		A = hy*hz;
		B = hz*xkj;
		C = hy*xji;
		D = hz*xki;
		E = hz*xji;
		F = hy*xji;

		//std::cout << "A = " << A << std::endl << "B = " << B << std::endl << "C = " << C << std::endl << "D = " << D << std::endl << "E = " << E << std::endl << "F = " << F << std::endl;

		double kbbe = (B*B + 3 * C*C - 2 * B*B*(D - E) + (D - E)*(D - E)) / areaFactor;
		double kbie = -(C*C) / areaFactor;
		double kbje = kbie;
		double kiie = (A*A + B*B + C*C) / areaFactor;
		double kije = -(A*A + B*D) / areaFactor;
		double kike = kije;
		double kjje = (A*A + C*C + D*D + E*E) / areaFactor;
		double kjke = -(E*(D + E)) / areaFactor;
		double kkke = (C*C + 2*E*E) / areaFactor;

		//std::cout << "kbbe = " << kbbe << std::endl << "kbie = " << kbie << std::endl << "kbje = " << kbje << std::endl << "kiie = " << kiie << std::endl << "kije = " << kije << std::endl << "kike = " << kike << std::endl << "kjje = " << kjje << std::endl << "kjke = " << kjke << std::endl << "kkke = " << kkke << std::endl;

		// somar aos acumuladores de coeficiente
		std::map<int, int>::const_iterator ii = coefficientMap.find(e.baseNode);
		if (ii != coefficientMap.end()) {
			int index = ii->second;
			// Add to the base node...
			insert(e.baseNode, e.baseNode, index, kbbe);
			insert(e.baseNode, p.a, index, kbie);
			insert(e.baseNode, p.b, index, kbje);
			// insert(e.baseNode, p.b, index, kbke); = 0

			// the i-th node...
			insert(p.a, e.baseNode, index, kbie);
			insert(p.a, p.a, index, kiie);
			insert(p.a, p.b, index, kije);
			insert(p.a, p.c, index, kike);

			// the j-th node...
			insert(p.b, e.baseNode, index, kbje);
			insert(p.b, p.a, index, kije);
			insert(p.b, p.b, index, kjje);
			insert(p.b, p.c, index, kjke);

			// ... and the k-th node
			//insert(p.c, e.baseNode, index, kbke); = 0
			insert(p.c, p.a, index, kike);
			insert(p.c, p.b, index, kjke);
			insert(p.c, p.c, index, kkke);
		}
	}
}