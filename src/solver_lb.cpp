#include "solver_lb.h"
#include "nodecoefficients.h"
#include "problemdescription.h"

void assembleProblemMatrix_lb(float *cond, matrix **Kii, matrix2 **Kic, matrix **Kcc, int numElect)
{
      int iiLimit = nodes.size()-numElect;

        matrix *out = new matrix(iiLimit, iiLimit);
        double val;
        out->startFill(); // estimate of the number of nonzeros (optional)
        int i;
        for (i=0; i<iiLimit; ++i) {
                nodeCoefficients *aux = nodeCoef[i];
                while(aux) { // Col-major storage
                        while(aux->node < i) aux = aux->next; // skip upper triangular
                        int row = aux->node;
                        if(row>=iiLimit) break; // Cut lower segment
                        val = 0.0;
                        while(aux && aux->node==row) {
                                val += aux->coefficient*cond[aux->condIndex];
                                aux = aux->next;
                        }
                        out->fill(row,i) = val;
                }
        }
        out->endFill();
        *Kii = out;
        // Now Kii and Kic
        matrix *out2 = new matrix(numElect-1, numElect-1);
        // Row major! Built as the transpose
        matrix2 *outLeft = new matrix2(numElect-1, iiLimit);
        out2->startFill();
        outLeft->startFill();
        for (; i<nodes.size()-1; ++i) {
                nodeCoefficients *aux = nodeCoef[i];
                while(aux) { // Col-major storage in Kcc, row major in Kic
                        while(aux->node > iiLimit && aux->node < i) aux = aux->next; // skip small upper triangular section at the bottom left
                        int row = aux->node;
                        val = 0.0;
                        while(aux && aux->node==row) {
                                val += aux->coefficient*cond[aux->condIndex];
                                aux = aux->next;
                        }
                        if(row<iiLimit)
                          outLeft->fill(i-iiLimit,row) = val; // As noted previously, outLeft is filled sideways row-wise
                        else 
                          out2->fill(row-iiLimit,i-iiLimit) = val;
                }
        }
        out->endFill();
        outLeft->endFill();
        *Kcc = out2;
        *Kic = outLeft;        
}
