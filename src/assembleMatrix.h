#ifndef ASSEMBLEMATRIX_H_
#define ASSEMBLEMATRIX_H_


#include "basematrix.h"

void prepareSkeletonMatrix();
void createCoef2KMatrix();
void assembleProblemMatrix(double *cond, matrix **stiffnes);
void assembleProblemMatrix_old(double *cond, matrix **stiffnes);


#endif	// ASSEMBLEMATRIX_H_