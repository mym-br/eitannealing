#include "solver_lb.h"
#include "problem.h"
#include <iostream>

// FIXME: Use new implementation
void assembleProblemMatrix_lb(double *cond, matrix **Kii, matrix **Kic, matrix **Kcc, problem &p)
{
      int iiLimit = p.getNodesCount()-p.getGenericElectrodesCount();

        matrix *out = new matrix(iiLimit, iiLimit);
        matrix *outBottom = new matrix(p.getGenericElectrodesCount(), iiLimit);
        double val;
        out->reserve(7*iiLimit); // estimate of the number of nonzeros (optional)
        int i;
        for (i=0; i<iiLimit; ++i) {
                nodeCoefficients *aux = p.getNodeCoefficients()[i];
                while(aux) { // Col-major storage
                        while(aux->node < i) aux = aux->next; // skip upper triangular
                        int row = aux->node;
                        val = 0.0;
                        while(aux && aux->node==row) {
                                val += aux->coefficient*cond[aux->condIndex];
                                aux = aux->next;
                        }
                        if(row>=iiLimit) {
                          outBottom->insert(row - iiLimit, i) = val;
                        }
                        else
                          out->insert(row,i) = val;
                }
        }
        out->makeCompressed();
        outBottom->makeCompressed();
        *Kii = out;
        *Kic = outBottom;
        // Now Kcc
        matrix *out2 = new matrix(p.getGenericElectrodesCount(), p.getGenericElectrodesCount());
        out2->reserve((p.getGenericElectrodesCount())*4);
        for (; i<p.getNodesCount(); ++i) {
                nodeCoefficients *aux =  p.getNodeCoefficients()[i];
                while(aux) { // Col-major storage in Kcc
                        while(aux->node < i) aux = aux->next; // skip upper triangular
                        int row = aux->node;
                        val = 0.0;
                        while(aux && aux->node==row) {
                                val += aux->coefficient*cond[aux->condIndex];
                                aux = aux->next;
                        }
                        out2->insert(row-iiLimit,i-iiLimit) = val;
                }
        }
        out2->makeCompressed();
        *Kcc = out2;
}
