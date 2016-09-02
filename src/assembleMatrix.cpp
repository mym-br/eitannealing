#include "basematrix.h"


// FIXME: It should be slightly faster to have an "skeleton"
//		matrix (with ones for instance), copy it and iterate over its
//		coefficients.

void assembleProblemMatrix(float *cond, matrix **stiffnes)
{

	/*
	// Jacobi preconditioning
	*sqDiagonal = new Eigen::VectorXd(numNodes-1);
	Eigen::VectorXd &sqrtdiagonal = **sqDiagonal;
	int i;
	// Prepare diagonal square roots
	for(i=0;i<numNodes-1;i++) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux && aux->node < i) aux = aux->next;
		double val = 0;
		while(aux && aux->node==i) {
			val += aux->coefficient*cond[aux->condIndex];
			aux = aux->next;
		}
		sqrtdiagonal[i] = sqrt(val);
	}


	matrix *out = new matrix(numNodes-1, numNodes-1);
	out->startFill(3*(numNodes-1)); // estimate of the number of nonzeros (optional)
	for (i=0; i<numNodes-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux && aux->node <= i) aux = aux->next; // skip upper triangular
		// Diagonal is 1
		out->fill(i,i) = 1;
		while(aux) { // Col-major storage
			int row = aux->node;
			double val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			out->fill(row,i) = val/(sqrtdiagonal[row]*sqrtdiagonal[i]);
		}
	}
	out->endFill();

	*stiffnes = out;*/


	matrix *out = new matrix(nodes.size()-1, nodes.size()-1);
	double val;
	out->startFill(3*(nodes.size()-1)); // estimate of the number of nonzeros (optional)
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
			out->fill(row,i) = val;
		}
	}
	out->endFill();

	*stiffnes = out;
}


void assembleProblemMatrix(float *cond, matrix **stiffnes, int numNodes, nodeCoefficients **nodeCoef)
{

	/*
	// Jacobi preconditioning
	*sqDiagonal = new Eigen::VectorXd(numNodes-1);
	Eigen::VectorXd &sqrtdiagonal = **sqDiagonal;
	int i;
	// Prepare diagonal square roots
	for(i=0;i<numNodes-1;i++) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux && aux->node < i) aux = aux->next;
		double val = 0;
		while(aux && aux->node==i) {
			val += aux->coefficient*cond[aux->condIndex];
			aux = aux->next;
		}
		sqrtdiagonal[i] = sqrt(val);
	}


	matrix *out = new matrix(numNodes-1, numNodes-1);
	out->startFill(3*(numNodes-1)); // estimate of the number of nonzeros (optional)
	for (i=0; i<numNodes-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux && aux->node <= i) aux = aux->next; // skip upper triangular
		// Diagonal is 1
		out->fill(i,i) = 1;
		while(aux) { // Col-major storage
			int row = aux->node;
			double val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			out->fill(row,i) = val/(sqrtdiagonal[row]*sqrtdiagonal[i]);
		}
	}
	out->endFill();

	*stiffnes = out;*/


	matrix *out = new matrix(numNodes-1, numNodes-1);
	double val;
	out->startFill(3*(numNodes-1)); // estimate of the number of nonzeros (optional)
	for (int i=0; i<numNodes-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux) { // Col-major storage
			while(aux->node < i) aux = aux->next; // skip upper triangular
			int row = aux->node;
			val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			out->fill(row,i) = val;
		}
	}
	out->endFill();

	*stiffnes = out;
}
