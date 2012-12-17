#include "solver_lb.h"
#include "nodecoefficients.h"
#include "problemdescription.h"

LB_Solver::LB_Solver(matrix *_Aii, matrix2 *_Aic, matrix *_Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, double a):
    Aii(*_Aii), Aic(*_Aic), precond(precond), a(a), lowerSafe(true), x0(Eigen::VectorXd::Zero(_Aii->rows()))
{
    it = 0;
    // 0
    r = -_Aic->transpose()*Phi;
    rc = J.end(_Acc->rows()) - *_Acc*Phi;      
    init();  
}

LB_Solver::LB_Solver(matrix *_Aii, matrix2 *_Aic, matrix *_Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, double a, const Eigen::VectorXd &x0):
    Aii(*_Aii), Aic(*_Aic), precond(precond), a(a), lowerSafe(true), x0(x0)
{
    it = 0;
    // 0
    r = -_Aic->transpose()*Phi;
    rc = J.end(_Acc->rows()) - *_Acc*Phi;      
    //std::cout << "Before:" << sqrt(r.squaredNorm()+rc.squaredNorm()); 
    Eigen::VectorXd xaux(x0);
    //precond.solveInPlace(xaux);
    r -=  Aii*xaux;
    rc -= Aic*xaux;
    //std::cout << "After:" << sqrt(r.squaredNorm()+rc.squaredNorm()) << std::endl;
    
    init();  
}
    
void LB_Solver::init()
{
    JhatNorm2 = r.squaredNorm()+rc.squaredNorm();
    delta = sqrt(JhatNorm2);
    p = r/delta; pc = rc/delta;
    
    // S = (ACi)T*p
    s = Aii.transpose()*p.lazy(); s += Aic.transpose()*pc.lazy();
    precond.solveInPlaceT(s);
	ATJhatNorm2 = s.squaredNorm();
    gamma_ip = sqrt(ATJhatNorm2); // gamma of *NEXT* iteration is obtained here!
    	ATJhatNorm2*=JhatNorm2;
	   
    q = s/gamma_ip;              // uses gamma of *NEXT* iteration
	// *** Gauss
    g=0;
    
    // 1
    w = q;
    fit = delta;
    
	qaux = q;
    precond.solveInPlace(qaux); // q  <- q*C^-1
    r = Aii*qaux.lazy(); rc = Aic*qaux.lazy();
	r -= gamma_ip*p; rc -= gamma_ip*pc;
    delta = sqrt(r.squaredNorm()+rc.squaredNorm());
	p = r/delta; pc = rc/delta;
	s = Aii.transpose()*p.lazy(); s += Aic.transpose()*pc.lazy();
	precond.solveInPlaceT(s);
	s -= delta*q;
	    // *** Gauss, as next value for gamma will be pertinent to iteration 2!
        phi2 = gamma_ip*gamma_ip+delta*delta;
		phi = sqrt(phi2);
		c = -gamma_ip/phi;
		si = delta/phi;
		pi = 1/phi2;
		g+=pi;    
		// This is due to gauss-radau
		alpha = gamma_ip*gamma_ip+delta*delta;
	gamma_ip = s.norm();
	q = s/gamma_ip;

	// Gauss-radau
	beta = gamma_ip*delta;
	at = a;
	dt = alpha - a;
    
    x = - (c*fit/phi)*w;
        
    it = 1;        
}

void LB_Solver::do_iteration()
{
    qaux = q;
    precond.solveInPlace(qaux); // q  <- q*C^-1
    r = Aii*qaux.lazy(); rc = Aic*qaux.lazy();
	r -= gamma_ip*p; rc -= gamma_ip*pc;
    delta = sqrt(r.squaredNorm()+rc.squaredNorm());
	p = r/delta; pc = rc/delta;
    s = Aii.transpose()*p.lazy(); s += Aic.transpose()*pc.lazy();
    precond.solveInPlaceT(s);
	s -= delta*q;
        // *** Gauss, as next value for gamma will be pertinent to next iteration!
		psi_im = si*gamma_ip;
                
                fit *= si;
                w = q - (psi_im/phi)*w;
                phi2 = c*c*gamma_ip*gamma_ip+delta*delta;
		phi = sqrt(phi2);
		c *= -gamma_ip/phi;
		si = delta/phi;
		pi_im = pi;
		pi *= psi_im*psi_im/phi2;
		g_im = g;
		g+=pi;    
		// This is due to gauss-radau
		alpha = gamma_ip*gamma_ip+delta*delta;
		
    gamma_ip = s.norm();
	q = s/gamma_ip;

	// Gauss-radau
	at = a + beta*beta/dt;
	phi2t = at - psi_im*psi_im;
	dt = alpha - a - (beta*beta/dt);
	beta = gamma_ip*delta;
		
	gr = g_im + pi_im*(psi_im*psi_im/phi2t);
    
     x -= (c*fit/phi)*w;
     //std::cout << "x_1[" << it+1 << "]:" << x[0] << std::endl;
     //std::cout << "fit[" << it+1 << "]:" << fit << std::endl; 
     it++;
}

double LB_Solver::getErrorl2Estimate() const
{
    return sqrt(JhatNorm2 - ATJhatNorm2*g);//+this->getX().norm()*0.0005;
}


double LB_Solver::getMinErrorl2Estimate() const
{
    if(!lowerSafe) return 0;
        double v = JhatNorm2 - ATJhatNorm2*(gr);
	if(v<0) return 0;//this->getX().norm()*0.0005;
	return sqrt(v);//+this->getX().norm()*0.0005;
}

LB_Solver_EG_Estimate::LB_Solver_EG_Estimate(matrix* Aii, matrix2* Aic, matrix* Acc, const Eigen::VectorXd& J, const Eigen::VectorXd& Phi, const Preconditioner& precond, int n, float e): 
        LB_Solver(Aii, Aic, Acc, J, Phi, precond, 0)
{
      UMatrix U(n,n);
      // Alpha and beta must be stored in order to recalculate dt
      std::vector<double> AlphaVector;    
      std::vector<double> BetaVector;
      AlphaVector.push_back(alpha);
      BetaVector.push_back(beta);
      U.startFill(2*n-1);
      U.fill(0,0) = phi;
      for(int i=1;i<n;i++) {
        this->do_iteration();
        AlphaVector.push_back(alpha);
        BetaVector.push_back(beta);      
        U.fill(i-1,i) = psi_im;
        U.fill(i,i) = phi;
      }
      U.endFill();
      // Now calc eigenvalue     
      double oev=0;
      ev = 1.0;
      evec = Eigen::VectorXd::Constant(n,1/sqrt(n));
      while(fabs(oev-ev)/ev > e) {
        U.transpose().solveTriangularInPlace(evec);
        U.solveTriangularInPlace(evec);
        oev = ev;
        ev = 1/evec.norm(); 
        evec *= ev;
      }
      this->a = ev;
      // We need now to recalculate Dt...
      //std::cout << "Alpha[1]:"<<AlphaVector[0]<<std::endl; 
      dt = AlphaVector[0] - a;
      //std::cout << "dt[1]:"<<dt<<std::endl; 
      for(int i=1;i<this->it;i++) {
        dt = AlphaVector[i] - a - (BetaVector[i-1]*BetaVector[i-1])/dt;
        //std::cout << "dt[" << i+1 << "]:" << dt << std::endl;
      }
      // Update Gauss-Radau values 
      this->do_iteration();      
}


LB_Solver_EG_Estimate::LB_Solver_EG_Estimate(matrix *Aii, matrix2 *Aic, matrix *Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, const Eigen::VectorXd &x0, const Eigen::VectorXd &egHint, int n, float e):
 LB_Solver(Aii, Aic, Acc, J, Phi, precond, 0, x0)
{
    UMatrix U(n,n);
      // Alpha and beta must be stored in order to recalculate dt
      std::vector<double> AlphaVector;    
      std::vector<double> BetaVector;
      AlphaVector.push_back(alpha);
      BetaVector.push_back(beta);
      U.startFill(2*n-1);
      U.fill(0,0) = phi;
      for(int i=1;i<n;i++) {
        this->do_iteration();
        AlphaVector.push_back(alpha);
        BetaVector.push_back(beta);      
        U.fill(i-1,i) = psi_im;
        U.fill(i,i) = phi;
      }
      U.endFill();
      // Now calc eigenvalue     
      double oev=0;
      ev = 1.0;
      evec = egHint;
      while(fabs(oev-ev)/ev > e) {
        U.transpose().solveTriangularInPlace(evec);
        U.solveTriangularInPlace(evec);
        oev = ev;
        ev = 1/evec.norm(); 
        evec *= ev;
      }
      this->a = ev;
      // We need now to recalculate Dt...
      //std::cout << "Alpha[1]:"<<AlphaVector[0]<<std::endl; 
      dt = AlphaVector[0] - a;
      //std::cout << "dt[1]:"<<dt<<std::endl; 
      for(int i=1;i<this->it;i++) {
        dt = AlphaVector[i] - a - (BetaVector[i-1]*BetaVector[i-1])/dt;
        //std::cout << "dt[" << i+1 << "]:" << dt << std::endl;
      }
      // Update Gauss-Radau values 
      this->do_iteration();
    
}


void assembleProblemElectrodeIdentityMatrix(float *cond, matrix2 **Kic, int numElect)
{
        int iiLimit = nodes.size()-numElect;      
        matrix2 *out = new matrix2(numElect-1, nodes.size()-1);
        double val;
        out->startFill(numElect-1); // estimate of the number of nonzeros (optional)
        for (int i=iiLimit; i<nodes.size()-1; ++i) {
              out->fill(i-iiLimit, i) = totalheight*mincond/2;
        }
        out->endFill();

        *Kic = out;
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
