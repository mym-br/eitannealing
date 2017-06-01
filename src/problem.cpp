#include "problem.h"
#include "twodim/problem2D.h"
#include "threedim/problem3D.h"
#include <memory>
#include <iostream>
#include <fstream>

#define TETRAHEDRON_TYPE 4

std::shared_ptr<problem> problem::createNewProblem(const char *meshfilename, bool &is2D) {
	std::ifstream file;
	file.open(meshfilename);

	std::string aux;
	do {
		std::getline(file, aux);
	} while (aux.compare("$ELM") != 0 && aux.compare("$Elements") != 0);
	// Get number of elements
	int elementCout;
	file >> elementCout;
	// Check if there are any thetahedral elements
	int id, eltype;
	for (int i = 0; i < elementCout; i++) {
		file >> id >> eltype;
		if (eltype == TETRAHEDRON_TYPE)  {
			is2D = false;
			return std::shared_ptr<problem>(new problem3D(meshfilename));
		}
		std::getline(file, aux);
	}
	file.close();

	is2D = true;
	return std::shared_ptr<problem>(new problem2D(meshfilename));
}

void problem::setGroundNode(int nodeid) {
	this->groundNode = nodeid;
	#ifndef BLOCKGND
	// Zeroing current vector for ground node. FIXME: Ground not can only be set once!
	for (int i = 0; i < nobs; i++) currents[i][this->groundNode] = 0;
	#endif
}

problem::~problem(){
	delete[] tensions;
	delete[] rebased_tensions;
	delete[] currents;
	delete[] tensionscomplex;
	delete[] rebased_tensionscomplex;
	delete[] currentscomplex;
	delete coef2KMatrix;
	delete skeleton;
	for (int i = 0; i < getNodesCount(); i++) {
		nodeCoefficients *node = nodeCoef[i];
		while (node != NULL)
		{
			nodeCoefficients* tmpNode = node->next;
			delete node;
			node = tmpNode;
		}
	}
	delete[] nodeCoef;
};