#include "solver_lb.h"
#include "nodecoefficients.h"
#include "problemdescription.h"

LB_Solver::LB_Solver(matrix *_Aii, matrix2 *_Aic, matrix *_Acc, Eigen::VectorXd &J, Eigen::VectorXd &Phi, const SparseIncompleteLLT &precond):
    Aii(*_Aii), Aic(*_Aic), precond(precond)
{
    
    // 0
    r = -_Aic->transpose()*Phi;
    rc = J.end(_Acc->rows()) - *_Acc*Phi;
    JhatNorm2 = r.squaredNorm()+rc.squaredNorm();
    delta = sqrt(JhatNorm2);
    
    p = r/delta; pc = rc/delta;
    
    // S = (ACi)T*p
    s = Aii.transpose()*p.lazy(); s += Aic.transpose()*pc.lazy();
    precond.solveInPlace(s);
    ATJhatNorm2 = s.squaredNorm();
    gamma_ip = sqrt(ATJhatNorm2); // gamma of *NEXT* iteration is obtained here!
    
    ATJhatNorm2*=JhatNorm2;
    
    q = s/gamma_ip;              // uses gamma of *NEXT* iteration
    std::cout << "Gamma:" << gamma_ip << std::endl;
    // *** Gauss
    g=0;
    
    // 1
    precond.solveInPlace(q); // q  <- q*C^-1
    r = Aii*q.lazy(); rc = Aic*q.lazy();
    delta = sqrt(r.squaredNorm()+rc.squaredNorm());
    p = r/delta; pc = rc/delta;
    s = Aii.transpose()*p.lazy(); s += Aic.transpose()*pc.lazy();
    precond.solveInPlace(s);
        // *** Gauss, as next value for gamma will be pertinent to iteration 2!
        phi2 = gamma_ip*gamma_ip+delta*delta;
    gamma_ip = s.norm();
    q = s/gamma_ip;

    // *** Gauss
    pi = 1/phi2;
    g+=pi;
    
    
    
    it = 1;        
}

void LB_Solver::do_iteration()
{
    precond.solveInPlace(q); // q  <- q*C^-1
    r = Aii*q.lazy(); rc = Aic*q.lazy();
    delta = sqrt(r.squaredNorm()+rc.squaredNorm());
    p = r/delta; pc = rc/delta;
    s = Aii.transpose()*p.lazy(); s += Aic.transpose()*pc.lazy();
    precond.solveInPlace(s);
        // *** Gauss, as next value for gamma will be pertinent to *NEXT*iteration!
        phi2 = gamma_ip*gamma_ip+delta*delta;    
    gamma_ip = s.norm();
    q = s/gamma_ip;
        
    // *** Gauss
    pi /= phi2;
    g+=pi;
    
    it++;
}

double LB_Solver::getErrorl2Estimate() const
{
    return sqrt(JhatNorm2 - ATJhatNorm2*g);
}




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
