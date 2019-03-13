#include "solver_lb_complex.h"
#include "nodecoefficients.h"
#include <iostream>

LB_Solver_Complex::LB_Solver_Complex(matrix *_Aii_R, matrix *_Aii_I, matrix *_Aic_R, matrix *_Aic_I, matrix *_Acc_R, matrix *_Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, double a):
    Aii_R(*_Aii_R), Aii_I(*_Aii_I), Aic_R(*_Aic_R), Aic_I(*_Aic_I), precond(precond), lowerSafe(true), a(a), x0_R(Eigen::VectorXd::Zero(_Aii_R->rows())), x0_I(Eigen::VectorXd::Zero(_Aii_R->rows()))
{
    it = 0;
    // 0
    r_R = -Aic_R.transpose()*Phi_R + Aic_I.transpose()*Phi_I;
    r_I = -Aic_R.transpose()*Phi_I - Aic_I.transpose()*Phi_R;
    rc_R = J_R - _Acc_R->selfadjointView<Eigen::Lower>()*Phi_R + _Acc_I->selfadjointView<Eigen::Lower>()*Phi_I;
    rc_I = J_I - _Acc_R->selfadjointView<Eigen::Lower>()*Phi_I - _Acc_I->selfadjointView<Eigen::Lower>()*Phi_R;

    init();
}

LB_Solver_Complex::LB_Solver_Complex(matrix *_Aii_R, matrix *_Aii_I, matrix *_Aic_R, matrix *_Aic_I, matrix *_Acc_R, matrix *_Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, double a, const Eigen::VectorXd &_x0_R, const Eigen::VectorXd &_x0_I):
    Aii_R(*_Aii_R), Aii_I(*_Aii_I), Aic_R(*_Aic_R), Aic_I(*_Aic_I), precond(precond), lowerSafe(true), a(a), x0_R(_x0_R), x0_I(_x0_I)
{
    it = 0;
    // 0
    r_R = -Aic_R.transpose()*Phi_R + Aic_I.transpose()*Phi_I;
    r_I = -Aic_R.transpose()*Phi_I - Aic_I.transpose()*Phi_R;
    rc_R = J_R - _Acc_R->selfadjointView<Eigen::Lower>()*Phi_R + _Acc_I->selfadjointView<Eigen::Lower>()*Phi_I;
    rc_I = J_I - _Acc_R->selfadjointView<Eigen::Lower>()*Phi_I - _Acc_I->selfadjointView<Eigen::Lower>()*Phi_R;

    r_R -=  Aii_R*x0_R - Aii_I*x0_I;
    r_I -=  Aii_R*x0_I + Aii_I*x0_R;
    rc_R -= Aic_R*x0_R - Aic_I*x0_I;
    rc_I -= Aic_R*x0_I + Aic_I*x0_R;

    init();
}

void LB_Solver_Complex::init()
{
    JhatNorm2 = r_R.squaredNorm()+r_I.squaredNorm()+rc_R.squaredNorm()+rc_I.squaredNorm();
    delta = sqrt(JhatNorm2);
    p_R = r_R/delta; p_I = r_I/delta;
    pc_R = rc_R/delta; pc_I = rc_I/delta;
    // S = (ACi)T*p
    s_R.noalias() = Aii_R*p_R;
    s_R.noalias() += Aii_I*p_I;
    s_I.noalias() = Aii_R*p_I - Aii_I*p_R;
    s_R.noalias() += Aic_R.transpose()*pc_R;
    s_R.noalias() += Aic_I.transpose()*pc_I;
    s_I.noalias() += Aic_R.transpose()*pc_I;
    s_I.noalias() -= Aic_I.transpose()*pc_R;

    precond.solveInPlaceCT(s_R, s_I);

    std::cout << s_R << "\n\n" << s_I << "\n";

    ATJhatNorm2 = s_R.squaredNorm();
    ATJhatNorm2 += s_I.squaredNorm();
    gamma_ip = sqrt(ATJhatNorm2); // gamma of *NEXT* iteration is obtained here!
    ATJhatNorm2*=JhatNorm2;

    q_R = s_R/gamma_ip; q_I = s_I/gamma_ip;               // uses gamma of *NEXT* iteration
	// *** Gauss
    g=0;

        // 1
        w_R = q_R; w_I = q_I;
        fit = delta;

	      qaux_R = q_R; qaux_I = q_I;
        precond.solveInPlaceC(qaux_R, qaux_I); // q  <- q*C^-1

        r_R.noalias() = Aii_R*qaux_R - Aii_I*qaux_I;
    		r_I.noalias() = Aii_R*qaux_I + Aii_I*qaux_R;
    		rc_R.noalias() = Aic_R*qaux_R - Aic_I*qaux_I;
    		rc_I.noalias() = Aic_R*qaux_I + Aic_I*qaux_R;

        r_R -= gamma_ip*p_R; r_I -= gamma_ip*p_I;
        rc_R -= gamma_ip*pc_R; rc_I -= gamma_ip*pc_I;
        delta = sqrt(r_R.squaredNorm()+r_I.squaredNorm()+rc_R.squaredNorm()+rc_I.squaredNorm());
        p_R = r_R/delta; p_I = r_I/delta; pc_R = rc_R/delta; pc_I = rc_I/delta;

        s_R.noalias() = Aii_R*p_R + Aii_I*p_I; 	// Aii^* = Aii
        s_I.noalias() = Aii_R*p_I - Aii_I*p_R; 	// Aii^* = Aii
        s_R.noalias() += Aic_R.transpose()*pc_R;
        s_R.noalias() += Aic_I.transpose()*pc_I;
        s_I.noalias() += Aic_R.transpose()*pc_I;
        s_I.noalias() -= Aic_I.transpose()*pc_R;

	      precond.solveInPlaceCT(s_R, s_I);
        s_R -= delta*q_R;
        s_I -= delta*q_I;
	       // *** Gauss, as next value for gamma will be pertinent to iteration 2!
        phi2 = gamma_ip*gamma_ip+delta*delta;
		    phi = sqrt(phi2);
		    c = -gamma_ip/phi;
		    si = delta/phi;
		    pi = 1/phi2;
		    g+=pi;
		    // This is due to gauss-radau
		    alpha = gamma_ip*gamma_ip+delta*delta;
	      gamma_ip = sqrt(s_R.squaredNorm()+s_I.squaredNorm());
	      q_R = s_R/gamma_ip;
        q_I = s_I/gamma_ip;

	      // Gauss-radau
	      beta = gamma_ip*delta;
	      at = a;
	      dt = alpha - a;

        x_R = - (c*fit/phi)*w_R;
        x_I = - (c*fit/phi)*w_I;

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

		r_R -= gamma_ip*p_R; r_I -= gamma_ip*p_I;
		rc_R -= gamma_ip*pc_R; rc_I -= gamma_ip*pc_I;
    delta = sqrt(r_R.squaredNorm()+r_I.squaredNorm()+rc_R.squaredNorm()+rc_I.squaredNorm());
		p_R = r_R/delta; p_I = r_I/delta; pc_R = rc_R/delta; pc_I = rc_I/delta;
    s_R.noalias() = Aii_R*p_R + Aii_I*p_I; 	// Aii^* = Aii
		s_I.noalias() = Aii_R*p_I - Aii_I*p_R; 	// Aii^* = Aii
		s_R.noalias() += Aic_R.transpose()*pc_R;
		s_R.noalias() += Aic_I.transpose()*pc_I;
		s_I.noalias() += Aic_R.transpose()*pc_I;
		s_I.noalias() -= Aic_I.transpose()*pc_R;

		precond.solveInPlaceCT(s_R, s_I);
		s_R -= delta*q_R;
		s_I -= delta*q_I;
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

double LB_Solver_Complex::getErrorl2Estimate() const
{
    return sqrt(JhatNorm2 - ATJhatNorm2*g);//+this->getX().norm()*0.0005;
}


double LB_Solver_Complex::getMinErrorl2Estimate() const
{
    if(!lowerSafe) return 0;
        double v = JhatNorm2 - ATJhatNorm2*(gr);
	if(v<0) return 0;//this->getX().norm()*0.0005;
	return sqrt(v);//+this->getX().norm()*0.0005;
}

LB_Solver_EG_Complex_Estimate::LB_Solver_EG_Complex_Estimate(matrix *Aii_R, matrix *Aii_I, matrix *Aic_R, matrix *Aic_I, matrix *Acc_R, matrix *Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, int n, float e):
        LB_Solver_Complex(Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_I, J_R, J_I, Phi_R, Phi_I, precond, 0)
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


LB_Solver_EG_Complex_Estimate::LB_Solver_EG_Complex_Estimate(matrix *Aii_R, matrix *Aii_I, matrix *Aic_R, matrix *Aic_I, matrix *Acc_R, matrix *Acc_I, const Eigen::VectorXd &J_R, const Eigen::VectorXd &J_I, const Eigen::VectorXd &Phi_R, const Eigen::VectorXd &Phi_I, const Preconditioner &precond, const Eigen::VectorXd &x0_R, const Eigen::VectorXd &x0_I, const Eigen::VectorXd &egHint, int n, float e):
 LB_Solver_Complex(Aii_R, Aii_I, Aic_R, Aic_I, Acc_R, Acc_I, J_R, J_I, Phi_R, Phi_I, precond, 0, x0_R, x0_I)
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

/*// FIXME: Use new implementation
void assembleProblemMatrix_lb(double *cond, matrix **Kii, matrix **Kic, matrix **Kcc, int numElect)
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
        matrix *outLeft = new matrix(numElect-1, iiLimit);
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
}*/
