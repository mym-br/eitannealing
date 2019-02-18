#include "solver_lb_complex.h"
#include "nodecoefficients.h"
#include "problemdescription.h"

LB_Solver_Complex::LB_Solver_Complex(matrix *_Aii, matrix2 *_Aic, matrix *_Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, double a):
    Aii(*_Aii), Aic(*_Aic), precond(precond), a(a), lowerSafe(true), x0(Eigen::VectorXd::Zero(_Aii->rows()))
{
    it = 0;
    // 0
    r = -_Aic->transpose()*Phi;
    rc = J.tail(_Acc->rows()) - _Acc->selfadjointView<Eigen::Lower>()*Phi;
    init();
}

LB_Solver_Complex::LB_Solver_Complex(matrix *_Aii, matrix2 *_Aic, matrix *_Acc, const Eigen::VectorXd &J, const Eigen::VectorXd &Phi, const Preconditioner &precond, double a, const Eigen::VectorXd &x0):
    Aii(*_Aii), Aic(*_Aic), precond(precond), a(a), lowerSafe(true), x0(x0)
{
    it = 0;
    // 0
    r = -_Aic->transpose()*Phi;
    rc = J.tail(_Acc->rows()) - _Acc->selfadjointView<Eigen::Lower>()*Phi;
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
    s.noalias() = Aii*p;	// Aii^T = Aii
    s.noalias() += Aic.transpose()*pc;
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
        r.noalias() = Aii*qaux;
        rc.noalias() = Aic*qaux;
	r -= gamma_ip*p; rc -= gamma_ip*pc;
        delta = sqrt(r.squaredNorm()+rc.squaredNorm());
	p = r/delta; pc = rc/delta;
	s.noalias() = Aii*p; 	// Aii^T = Aii
        s.noalias() += Aic.transpose()*pc;
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

void LB_Solver_Complex::do_iteration()
{
    qaux_R = q_R; qaux_I = q_I;
    precond.solveInPlaceC(qaux_R, qaux_I);
    r_R.noalias() = Aii_R*qaux_R - Aii_I*qaux_I;
		r_I.noalias() = Aii_R*qaux_I + Aii_I*qaux_R;
		rc_R.noalias() = Aic_R*qaux_R - Aic_I*qaux_I;
		rc_I.noalias() = Aic_R*qaux_I + Aic_I*qaux_R;

		r_R -= gamma_ip*p_R; r_i -= gamma_ip*p_I;
		rc_R -= gamma_ip*pc_R; rc_I -= gamma_ip*pc_I;
    delta = sqrt(r_R.squaredNorm()+r_I.squaredNorm()+rc_R.squaredNorm()+rc_I.squaredNorm());
		p_R = r_R/delta; p_I = r_I/delta; pc_R = rc_R/delta; pc_I = rc_I/delta;
    s_R.noalias() = Aii_R*p_R + Aii_I*p_I; 	// Aii^* = Aii
		s_I.noalias() = Aii_R*p_I - Aii_I*p_R; 	// Aii^* = Aii
		s_R.noalias() += Aic_R.transpose()*pc_R;
		s_R.noalias() += Aic_I.transpose()*pc_I;
		s_I.noalias() += Aic_R.transpose()*pc_I;
		s_R.noalias() -= Aic_I.transpose()*pc_R;

		precond.solveInPlaceCT(s_R, s_I);
		s_R -= delta*q_R;
		s_I -= delta*q_i;
        // *** Gauss, as next value for gamma will be pertinent to next iteration!
		psi_im = si*gamma_ip;

    fit *= si;
    w_R = q_R - (psi_im/phi)*w_R;
		w_I = q_I - (psi_im/phi)*w_I;
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

    gamma_ip = sqrt(s_R.squaredNorm()+s_I.squaredNorm());
		q_R = s_R/gamma_ip;
		q_I = s_I/gamma_ip;


	// Gauss-radau
	at = a + beta*beta/dt;
	phi2t = at - psi_im*psi_im;
	dt = alpha - a - (beta*beta/dt);
	beta = gamma_ip*delta;

	gr = g_im + pi_im*(psi_im*psi_im/phi2t);

     x_R -= (c*fit/phi)*w_R;
		 x_I -= (c*fit/phi)*w_I;
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
      Eigen::SparseMatrix<double, Eigen::ColMajor>  U(n,n);
      // Alpha and beta must be stored in order to recalculate dt
      std::vector<double> AlphaVector;
      std::vector<double> BetaVector;
      AlphaVector.push_back(alpha);
      BetaVector.push_back(beta);
      U.reserve(2*n-1);
      U.insert(0,0) = phi;
      for(int i=1;i<n;i++) {
        this->do_iteration();
        AlphaVector.push_back(alpha);
        BetaVector.push_back(beta);
        U.insert(i-1,i) = psi_im;
        U.insert(i,i) = phi;
      }
      U.makeCompressed();
      // Now calc eigenvalue
      double oev=0;
      ev = 1.0;
      evec = Eigen::VectorXd::Constant(n,1/sqrt(n));
      while(fabs(oev-ev)/ev > e) {
        U.triangularView<Eigen::Upper>().transpose().solveInPlace(evec);
        U.triangularView<Eigen::Upper>().solveInPlace(evec);
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
      Eigen::SparseMatrix<double, Eigen::ColMajor>  U(n,n);
      // Alpha and beta must be stored in order to recalculate dt
      std::vector<double> AlphaVector;
      std::vector<double> BetaVector;
      AlphaVector.push_back(alpha);
      BetaVector.push_back(beta);
      U.reserve(2*n-1);
      U.insert(0,0) = phi;
      for(int i=1;i<n;i++) {
        this->do_iteration();
        AlphaVector.push_back(alpha);
        BetaVector.push_back(beta);
        U.insert(i-1,i) = psi_im;
        U.insert(i,i) = phi;
      }
      U.makeCompressed();
      // Now calc eigenvalue
      double oev=0;
      ev = 1.0;
      evec = egHint;
      while(fabs(oev-ev)/ev > e) {
        U.triangularView<Eigen::Upper>().transpose().solveInPlace(evec);
        U.triangularView<Eigen::Upper>().solveInPlace(evec);
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

// FIXME: Use new implementation
void assembleProblemMatrix_lb(double *cond, matrix **Kii, matrix2 **Kic, matrix **Kcc, int numElect)
{
      int iiLimit = nodes.size()-numElect;

        matrix *out = new matrix(iiLimit, iiLimit);
        double val;
        out->reserve(7*iiLimit); // estimate of the number of nonzeros (optional)
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
                        out->insert(row,i) = val;
                }
        }
        out->makeCompressed();
        *Kii = out;
        // Now Kcc and Kic
        matrix *out2 = new matrix(numElect-1, numElect-1);
        // Row major! Built as the transpose
        matrix2 *outLeft = new matrix2(numElect-1, iiLimit);
        out2->reserve((numElect-1)*4);
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
                          outLeft->insert(i-iiLimit,row) = val; // As noted previously, outLeft is filled sideways row-wise
                        else
                          out2->insert(row-iiLimit,i-iiLimit) = val;
                }
        }
        out->makeCompressed();
        outLeft->makeCompressed();
        *Kcc = out2;
        *Kic = outLeft;
}
