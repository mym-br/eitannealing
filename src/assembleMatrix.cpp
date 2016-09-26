#include "assembleMatrix.h"
#include "basematrix.h"
#include "problemdescription.h"
#include "nodecoefficients.h"


matrix *skeleton;
matrix *coef2KMatrix;

void prepareSkeletonMatrix()
{
    std::vector<Eigen::Triplet<Scalar>> tripletList;    
    for (int i=0; i<nodes.size()-1; ++i) {
	nodeCoefficients *aux = nodeCoef[i];
	while(aux) { // Col-major storage
	    while(aux->node < i) aux = aux->next; // skip upper triangular
	    int row = aux->node;
	    if(row==groundNode) {
	      aux = aux->next;
	      continue;   // Skip ground node
	    }
	    while(aux && aux->node==row) aux = aux->next;
	    // 1.0 value as placeholder
	    tripletList.push_back(Eigen::Triplet<Scalar>(row, i, 1.0));
	}
    }
    skeleton = new matrix(nodes.size()-1, nodes.size()-1);
    skeleton->setFromTriplets(tripletList.begin(), tripletList.end());
    skeleton->makeCompressed();
}

void createCoef2KMatrix()
{
    Scalar *base = &skeleton->coeffRef(0,0);
    std::vector<Eigen::Triplet<Scalar> > tripletList;  
    for (int i=0; i<nodes.size()-1; ++i) {
        for(nodeCoefficients *aux = nodeCoef[i]; aux!=NULL; aux = aux->next) {
	    if(aux->node < i) continue; // skip upper triangular
	    if(aux->node==groundNode) continue;   // Skip ground node
	    // Find index 
	    matrix::StorageIndex row = (matrix::StorageIndex)(&skeleton->coeffRef(aux->node, i) - base);	    
	    tripletList.push_back(Eigen::Triplet<Scalar>(row, aux->condIndex, aux->coefficient));
	}
    }
    coef2KMatrix = new matrix(skeleton->nonZeros(), numcoefficients);
    coef2KMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
    coef2KMatrix->makeCompressed();
}


void assembleProblemMatrix(double *cond, matrix **stiffnes)
{
	// FIXME: Is there any way to get rid of the initial memset of zeros
        //	on the outer index?
	matrix *m = new matrix(skeleton->rows(), skeleton->cols());
	m->reserve(skeleton->nonZeros());
	// Now byte-copy outer and inner vectors from skeleton
	memcpy(m->outerIndexPtr(), skeleton->outerIndexPtr(), (m->rows()+1)*sizeof(matrix::StorageIndex));
	memcpy(m->innerIndexPtr(), skeleton->innerIndexPtr(), m->nonZeros()*sizeof(matrix::StorageIndex));
	
	// Final coefficient vector is the Sparse x Dense product of coef2KMatrix times coefficients
	Eigen::Map<vectorx>(m->valuePtr(), m->nonZeros()).noalias() = (*coef2KMatrix)*Eigen::Map<Eigen::VectorXd>(cond, numcoefficients);
	m->resizeNonZeros(skeleton->nonZeros());
	*stiffnes = m;
}


// Old algorithm, just as a baseline for tests
void assembleProblemMatrix_old(double *cond, matrix **stiffnes)
{
	matrix *out = new matrix(nodes.size()-1, nodes.size()-1);
	double val;
	std::vector<Eigen::Triplet<Scalar>> tripletList; 
	for (int i=0; i<nodes.size()-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux) { // Col-major storage
			while(aux->node < i) aux = aux->next; // skip upper triangular
			int row = aux->node;
                        if(row==groundNode) {
                          aux = aux->next;
                          continue;   // Skip ground node
                        }
                        val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			tripletList.push_back(Eigen::Triplet<Scalar>(row, i, val));
		}
	}

	out->setFromTriplets(tripletList.begin(), tripletList.end());
	out->makeCompressed();
	*stiffnes = out;
}

// Build full symetric matrix
void assembleProblemMatrix_old_sim(double *cond, matrix **stiffnes)
{
	matrix *out = new matrix(nodes.size()-1, nodes.size()-1);
	double val;
	std::vector<Eigen::Triplet<Scalar>> tripletList; 
	for (int i=0; i<nodes.size()-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux) { // Col-major storage
			int row = aux->node;
                        if(row==groundNode) {
                          aux = aux->next;
                          continue;   // Skip ground node
                        }
                        val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			tripletList.push_back(Eigen::Triplet<Scalar>(row, i, val));
		}
	}

	out->setFromTriplets(tripletList.begin(), tripletList.end());
	out->makeCompressed();
	*stiffnes = out;
}