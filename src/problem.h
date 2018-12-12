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
#define _USE_MATH_DEFINES
#include <math.h>
#include "observations.h"
#include "nodecoefficients.h"
class solution;

//#define TETRAHEDRON_TYPE 4

class problem {
        // FIXME: Rework this interface?
	friend class gradientNormRegularisation;
	friend class gradientNormRegularisation_old;
	friend class gradientNormRegularisationComplex;
	friend class solution;
	friend class solutioncomplex;
	friend class solutioncomplexcomplete;
	friend class solutionCuda;

	private:
	int groundNode;

	protected:
	// Node data
	float electrodeh;
	std::map<int, int> node2coefficient;
	int numcoefficients;
	int nodeCount;
	nodeCoefficients **nodeCoef;
	std::vector<std::pair<int, int> > innerAdjacency;
	matrix *skeleton;
	matrix *coef2KMatrix;
	const char *filename;
	double currentFreq;
	double capacitance;
	bool isCapacitive;
	int calibrationMode;

public:

	static std::shared_ptr<problem> createNewProblem(const char *meshfilename, bool *is2D);

	// Virtual functions
	virtual void initProblem(const char *meshfilename) = 0;
	virtual void setCalibrationMode(bool individualcoeffs = false) = 0;
	virtual void buildNodeCoefficients() = 0;
	virtual int getGenericElectrodesCount() = 0;
	virtual int getGenericElectrodesCoeffCount() = 0;
	virtual int getNodesCount() { return nodeCount; }
	virtual int getInnerAdjacencyCount() = 0;

	// Getters and setters
	int getNumCoefficients() { return numcoefficients; }
	nodeCoefficients **getNodeCoefficients() { return nodeCoef; }
	int getNode2Coefficient(int id) { return node2coefficient[id]; }
	const char* getMeshFilename() { return filename; }
	void setGroundNode(int nodeid) { this->groundNode = nodeid; }
	int getGroundNode() { return this->groundNode; }
	void setCurrentFreq(double _currentFreq) { this->currentFreq = _currentFreq; }
	double getCurrentFreq() { return this->currentFreq; }
	void setCapacitance(double _capacitance) { this->capacitance = _capacitance; isCapacitive = true; }
	int getCalibrationMode() { return this->calibrationMode; }
	std::pair<int, int> getAdjacency(int index) { return this->innerAdjacency[index]; }

	// Contructor and destructors
	problem(const char *meshfilename) : filename(meshfilename), groundNode(-1),
		skeleton(nullptr), coef2KMatrix(nullptr), nodeCoef(nullptr),
		capacitance(0.0), isCapacitive(false), calibrationMode(0) {};
	virtual ~problem();

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

	template <typename  _Scalar> void assembleProblemMatrix(_Scalar *cond, Eigen::SparseMatrix<_Scalar, Eigen::ColMajor> **stiffnes) {
		// FIXME: Is there any way to get rid of the initial memset of zeros
		//	on the outer index?
		Eigen::SparseMatrix<_Scalar, Eigen::ColMajor> *m = new Eigen::SparseMatrix<_Scalar, Eigen::ColMajor>(skeleton->rows(), skeleton->cols());
		m->reserve(skeleton->nonZeros());
		// Now byte-copy outer and inner vectors from skeleton
		memcpy(m->outerIndexPtr(), skeleton->outerIndexPtr(), (m->rows() + 1)*sizeof(typename Eigen::SparseMatrix<_Scalar, Eigen::ColMajor>::StorageIndex));
		memcpy(m->innerIndexPtr(), skeleton->innerIndexPtr(), m->nonZeros()*sizeof(typename Eigen::SparseMatrix<_Scalar, Eigen::ColMajor>::StorageIndex));

		// Final coefficient vector is the Sparse x Dense product of coef2KMatrix times coefficients
		Eigen::Map<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(m->valuePtr(), m->nonZeros()).noalias() = (*coef2KMatrix)*Eigen::Map<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> >(cond, numcoefficients);
		m->resizeNonZeros(skeleton->nonZeros());
		*stiffnes = m;
	}

	template <typename  _Scalar> void postAssembleProblemMatrix(Eigen::SparseMatrix<_Scalar, Eigen::ColMajor> **stiffnes) {
		if (!isCapacitive) {
			#ifdef BLOCKGND
			for (int i = getNodesCount() - nobs; i < getNodesCount(); i++)
			for (int j = i; j < getNodesCount(); j++) {
				std::complex<double> *val = &m->coeffRef(j, i);
				*val = *val + 1 / 32.0;
			}
			#else
			for (int i = 0; i < getGroundNode(); i++) *(&(*stiffnes)->coeffRef(getGroundNode(), i)) = 0.0;
			for (int i = getGroundNode() + 1; i < getNodesCount(); i++) *(&(*stiffnes)->coeffRef(i, getGroundNode())) = 0.0;
			*(&(*stiffnes)->coeffRef(getGroundNode(), getGroundNode())) = 1.0;
			#endif
		}
	}

	void addMatrixCapacitances(matrixcomplex **stiffnes) {
		for (int i = getNodesCount() - getGenericElectrodesCount(); i < getNodesCount(); i++) {
			Complex *val = &(*stiffnes)->coeffRef(i, i);
			Complex jwc = std::complex<double>(0, 2 * M_PI * getCurrentFreq() * capacitance);
			*val = *val + jwc;
		}
	}

	vectorx getCurrentVector(int i, observations<double> *obs) {
		vectorx current = obs->getCurrents()[i];
		#ifndef BLOCKGND
		current[groundNode] = 0;
		#endif
		return current;
	}

	vectorxcomplex getConjugatedCurrentVector(int i, matrixcomplex *stiffnes, observations<std::complex<double>> *obs) {
		vectorxcomplex current = obs->getCurrents()[i];
		matrixcomplex Atransconj = stiffnes->conjugate() + ((matrixcomplex)(stiffnes->selfadjointView<Eigen::Lower>())).triangularView<Eigen::StrictlyUpper>();
		current = Atransconj * current;
		#ifndef BLOCKGND
		if (!isCapacitive) current[groundNode] = 0;
		#endif
		return current;
	}

	double electrodevar, regularizationFactor;
};

const double mincond = 0.001;
const double maxcond = 0.3815;
//const double maxcond = 0.3;
const double minperm = 0.000000000070922044418976;//0.00000005;
const double maxperm = 0.0000000070922044418976;// 0.05;

#endif // PROBLEM_H_
