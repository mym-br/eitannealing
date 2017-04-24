#include "problem.h"

void problem::prepareSkeletonMatrix()
{
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

void problem::createCoef2KMatrix()
{
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


void problem::assembleProblemMatrix(double *cond, matrix **stiffnes)
{
	// FIXME: Is there any way to get rid of the initial memset of zeros
	//	on the outer index?
	matrix *m = new matrix(skeleton->rows(), skeleton->cols());
	m->reserve(skeleton->nonZeros());
	// Now byte-copy outer and inner vectors from skeleton
	memcpy(m->outerIndexPtr(), skeleton->outerIndexPtr(), (m->rows() + 1)*sizeof(matrix::StorageIndex));
	memcpy(m->innerIndexPtr(), skeleton->innerIndexPtr(), m->nonZeros()*sizeof(matrix::StorageIndex));

	// Final coefficient vector is the Sparse x Dense product of coef2KMatrix times coefficients
	Eigen::Map<vectorx>(m->valuePtr(), m->nonZeros()).noalias() = (*coef2KMatrix)*Eigen::Map<Eigen::VectorXd>(cond, numcoefficients);
	m->resizeNonZeros(skeleton->nonZeros());
	*stiffnes = m;

	for (int i = getNodesCount() - nobs; i < getNodesCount(); i++)
	for (int j = getNodesCount() - nobs; j < getNodesCount(); j++) {
		double *val = &(*stiffnes)->coeffRef(i, j);
		*val = *val + 1/32.0;
	}
}