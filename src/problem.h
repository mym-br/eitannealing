/*
* problem.h
*
*  Created on: Feb 23, 2017
*      Author: aksato
*/

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include <vector>
#include <map>
#include <fstream>
#include <functional>
#include <memory>
#include <Eigen/Core>
#include "basematrix.h"
class solution;

#define TETRAHEDRON_TYPE 4

struct nodeCoefficients {
	int node;
	int condIndex;
	double coefficient;
	nodeCoefficients *next;

	nodeCoefficients(int node, int index, double coefficient) :
		node(node), condIndex(index), coefficient(coefficient), next(NULL) {}

	nodeCoefficients(int node, int index, double coefficient, nodeCoefficients *next) :
		node(node), condIndex(index), coefficient(coefficient), next(next) {}
};

template <typename _Scalar, typename t_vector, typename t_matrix >
class problem {
	friend class gradientNormRegularisation;
	friend class gradientNormRegularisation_old;
	friend class gradientNormRegularisationComplex;
	friend class solution;
	friend class solutioncomplex;

	private:
	int groundNode;

	protected:
	// Node data
	float electrodeh;
	std::map<int, int> node2coefficient;
	int numcoefficients;
	int nodeCount;
	// Observations
	int nobs;
	t_vector *tensions;
	t_vector *currents;
	t_vector currentVals;
	nodeCoefficients **nodeCoef;
	std::vector<std::pair<int, int> > innerAdjacency;
	matrix *skeleton;
	matrix *coef2KMatrix;
	const char *filename;
	double currentFreq;

public:

	//static std::shared_ptr<problem> createNewProblem(const char *meshfilename, bool &is2D) {
		//std::ifstream file;
		//file.open(meshfilename);

		//std::string aux;
		//do {
		//	std::getline(file, aux);
		//} while (aux.compare("$ELM") != 0 && aux.compare("$Elements") != 0);
		//// Get number of elements
		//int elementCout;
		//file >> elementCout;
		//// Check if there are any thetahedral elements
		//int id, eltype;
		//for (int i = 0; i < elementCout; i++) {
		//	file >> id >> eltype;
		//	if (eltype == TETRAHEDRON_TYPE)  {
		//		is2D = false;
		//		return std::shared_ptr<problem>(new problem3D(meshfilename));
		//	}
		//	std::getline(file, aux);
		//}
		//file.close();

		//is2D = true;
		//return std::shared_ptr<problem>(new problem2D(meshfilename));
	//}
	static bool isProblem2D(const char *meshfilename) {
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
				return false;
			}
			std::getline(file, aux);
		}
		file.close();

		return true;
	}
	
	// Virtual functions
	virtual void initProblem(const char *meshfilename) = 0;
	virtual void buildNodeCoefficients() = 0;
	virtual int getGenericElectrodesCount() = 0;
	virtual int getNodesCount() { return nodeCount; }
	virtual int getInnerAdjacencyCount() = 0;

	// Getters and setters
	int getNumCoefficients() { return numcoefficients; }
	nodeCoefficients **getNodeCoefficients() { return nodeCoef; }
	int getNode2Coefficient(int id) { return node2coefficient[id]; }
	int getNObs() { return nobs; }
	t_vector *getTensions() { return tensions; }
	t_vector *getCurrents() { return currents; }
	const char* getMeshFilename() { return filename; }
	_Scalar getCurrentVal(int i) { return currentVals[i]; }
	int getCurrentsCount() { return (int)currentVals.size(); }
	void setGroundNode(int nodeid) { this->groundNode = nodeid; }
	int getGroundNode() { return this->groundNode; }
	void setCurrentFreq(double _currentFreq) { this->currentFreq = _currentFreq; }
	double getCurrentFreq() { return this->currentFreq; }

	// Contructor and destructors
	problem(const char *meshfilename) : filename(meshfilename), groundNode(-1), nobs(-1),
		tensions(nullptr), currents(nullptr),
		skeleton(nullptr), coef2KMatrix(nullptr), nodeCoef(nullptr) {};
	virtual ~problem() {
		delete[] tensions;
		delete[] currents;
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
	}

	void initObs(const char *filecurrents, const char* filename) {
		std::ifstream file;
		std::ifstream filec;

		filec.open(filecurrents);
		file.open(filename);

		int n = getGenericElectrodesCount() - 1;
		int valuesCount = std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());
		nobs = valuesCount / getGenericElectrodesCount();

		file.clear();
		file.seekg(0, std::ios::beg);
		tensions = new t_vector[nobs];
		currents = new t_vector[nobs];
		currentVals = t_vector(nobs);
		t_vector current(getNodesCount());
		current.fill(0);
		int baseIndex = (int)current.size() - n - 1;
		for (int i = 0; i<nobs; i++) {
			double c;
			int entry, exit;
			filec >> entry;
			filec >> exit;
			filec >> c;
			entry--; exit--;	// zero-based
			currentVals[i] = c;
			currents[i] = current;
			currents[i][baseIndex + entry] = 1;
			currents[i][baseIndex + exit] = -1;
			#ifndef BLOCKGND
			currents[i][getGroundNode()] = 0;
			#endif

			// read tensions from file
			tensions[i].resize(getGenericElectrodesCount());
			_Scalar val, avg;
			avg = 0;
			for (int j = 0; j<getGenericElectrodesCount(); j++) {
				file >> val;
				tensions[i][j] = val / c;  // Values are normalized by current
				avg += tensions[i][j];
			}
			// rebase tensions, apply offset for zero sum on electrodes
			avg /= (double)getGenericElectrodesCount();
			for (int j = 0; j < getGenericElectrodesCount(); j++) {
				tensions[i][j] -= avg;
			}
		}

		filec.close();
		file.close();
	}

	void prepareSkeletonMatrix() {
		std::vector<Eigen::Triplet<Scalar>> tripletList;
		for (int i = 0; i<getNodesCount(); ++i) {
			nodeCoefficients *aux = nodeCoef[i];
			while (aux) { // Col-major storage
				while (aux->node < i) aux = aux->next; // skip upper triangular
				int row = aux->node;
				while (aux && aux->node == row) aux = aux->next;
				// 1.0 value as placeholder
				tripletList.push_back(Eigen::Triplet<Scalar>(row, i, 1.0));
			}
		}
		skeleton = new matrix(getNodesCount(), getNodesCount());
		skeleton->setFromTriplets(tripletList.begin(), tripletList.end());
		skeleton->makeCompressed();
	}

	void createCoef2KMatrix() {
		Scalar *base = skeleton->valuePtr();// coeffRef(0, 0);
		std::vector<Eigen::Triplet<Scalar> > tripletList;
		bool groundNodeSet = false;
		for (int i = 0; i<getNodesCount(); ++i) {
			for (nodeCoefficients *aux = nodeCoef[i]; aux != NULL; aux = aux->next) {
				if (aux->node < i) continue; // skip upper triangular
				// Find index 
				matrix::StorageIndex row = (matrix::StorageIndex)(&skeleton->coeffRef(aux->node, i) - base);
				tripletList.push_back(Eigen::Triplet<Scalar>(row, aux->condIndex, aux->coefficient));
			}
		}
		coef2KMatrix = new matrix(skeleton->nonZeros(), numcoefficients);
		coef2KMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
		coef2KMatrix->makeCompressed();
	}

	void assembleProblemMatrix(_Scalar *cond, t_matrix **stiffnes) {
		// FIXME: Is there any way to get rid of the initial memset of zeros
		//	on the outer index?
		t_matrix *m = new t_matrix(skeleton->rows(), skeleton->cols());
		m->reserve(skeleton->nonZeros());
		// Now byte-copy outer and inner vectors from skeleton
		memcpy(m->outerIndexPtr(), skeleton->outerIndexPtr(), (m->rows() + 1)*sizeof(t_matrix::StorageIndex));
		memcpy(m->innerIndexPtr(), skeleton->innerIndexPtr(), m->nonZeros()*sizeof(t_matrix::StorageIndex));

		// Final coefficient vector is the Sparse x Dense product of coef2KMatrix times coefficients
		Eigen::Map<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(m->valuePtr(), m->nonZeros()).noalias() = (*coef2KMatrix)*Eigen::Map<t_vector>(cond, numcoefficients);
		m->resizeNonZeros(skeleton->nonZeros());
		*stiffnes = m;
	}

	void postAssempleProblemMatrix(t_matrix **stiffnes) {
		#ifdef BLOCKGND
		for (int i = getNodesCount() - nobs; i < getNodesCount(); i++)
		for (int j = i; j < getNodesCount(); j++) {
			std::complex<double> *val = &m->coeffRef(j, i);
			*val = *val + 1 / 32.0;
		}
		#else
		for (int i = 0; i < getGroundNode(); i++) *(&(*stiffnes)->coeffRef(getGroundNode(), i)) = 0.0;
		*(&(*stiffnes)->coeffRef(getGroundNode(), getGroundNode())) = 1.0;
		#endif
	}

	vectorxcomplex getConjugatedCurrentVector(int i, matrixcomplex *stiffnes) {
		vectorxcomplex current = getCurrents()[i];
		return stiffnes->conjugate().selfadjointView<Eigen::Lower>() * current;
	}
};

const double mincond = 0.005;
const double maxcond = 0.3815;
const double minperm = 0.000000000070922044418976;//0.00000005;
const double maxperm = 0.0000000070922044418976;// 0.05;

#endif // PROBLEM_H_