/*
 * solver.cpp
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */


#include "solvercomplex.h"
//#include "nodecoefficients.h"
#include <fstream>

CG_SolverComplex::CG_SolverComplex(matrixcomplex &_A, const Eigen::VectorXcd &b, const SparseIncompleteLLTComplex &pre) :
	A(_A),
	b(b),
	x(Eigen::VectorXd::Zero(_A.rows())),
	precond(pre)
{
	this->init();
}

CG_SolverComplex::CG_SolverComplex(matrixcomplex &_A,  matrixcomplex &A_H, const Eigen::VectorXcd &b, const SparseIncompleteLLTComplex &pre) :
	A(_A),
	Afull(_A),
	b(b),
	x(Eigen::VectorXd::Zero(_A.rows())),
	precond(pre)
	{
	this->init(A_H);
}

CG_SolverComplex::CG_SolverComplex(matrixcomplex &A_, const Eigen::VectorXcd &b, const Eigen::VectorXcd &x0, const SparseIncompleteLLTComplex &pre) :
	A(A_), b(b), x(x0),
	// Precond!
	precond(pre)  {
	this->init();
}

// Setup and calculate the 1st iteraction
void CG_SolverComplex::init()
{
	r.resize(A.rows()); // matrixcomplex Aprint = A.matrix(); saveVals("Aprint.txt", Aprint); saveVals("bprint.txt", b);
	beta = 0;
	r = b - A*x; //saveVals("rnorm.txt", r.norm(), false);
	//p = r;

	z = r;
	precond.solveInPlace(z);
	rmod_1 = rmod = z.dot(r);
	p = z;

	//rmod_1 = rmod = r.squaredNorm();
	q = A*p;
	gamma = rmod / q.dot(p);
	r0norm = r.norm(); //r0norm = sqrt(rmod);
	r0norm2 = r0norm*r0norm; //r0norm2 = rmod;
	alpha = std::complex<double>(1.0, 0.0) / gamma;

	// ########## LANCZOS
	/*v = r/r0norm;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v;
	alpha_[0] = lalpha;
	std::cout << "Alpha by lanczos:" << lalpha << " by CG:" << alpha << std::endl;*/


	// Now do the first 3 iterations
	it = 0;
	x += gamma*p;
	r = r - gamma*q; //saveVals("rnorm.txt", r.norm());
	z = r;
	precond.solveInPlace(z);
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;

	eta_p1 = sqrt(beta)/gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = alpha;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	c = rt1/r1;
	s = eta_p1/r1;
	wt[0] = std::complex<double>(1.0, 0.0) / rt1;
	w[0] = std::complex<double>(1.0, 0.0) / r1;

	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A*p;
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	it++;
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[1] = lalpha;
	eta_[1] = leta;*/
	
	if (getResidueSquaredNorm() < 1e-19) return;
	x += gamma*p;	// 2...
	r = r - gamma*q; //saveVals("rnorm.txt", r.norm());
	z = r;
	precond.solveInPlace(z);
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta)/gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c*alpha - s*eta;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	r2 = c*eta + s*alpha;	// r_2,2 = c_1*eta_2
	c_1 = c;
	c = rt1/r1;
	s_1 = s;
	s = eta_p1/r1;
	w[1] = wt[1] = -r2*w[0];
	wt[1] /= rt1;
	w[1] /= r1;

	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A*p;
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	it++;
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[2] = lalpha;
	eta_[2] = leta;*/

	if (getResidueSquaredNorm() < 1e-19) return;
	x += gamma*p;	// ...and 3!
	r = r - gamma*q; //saveVals("rnorm.txt", r.norm());
	z = r;
	precond.solveInPlace(z);
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);// saveVals("rmod.txt", rmod);
	beta = rmod/rmod_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta)/gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c*alpha - s*c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	r2 = c_1*c*eta + s*alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1*eta;
	c_1 = c;
	c = rt1/r1;
	s_1 = s;
	s = eta_p1/r1;
	w[2] = wt[2] = -(r3*w[0]+r2*w[1]);
	wt[2] /= rt1;
	w[2] /= r1;

	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A*p;
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	it++;
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[3] = lalpha;
	eta_[3] = leta;*/

	err[0] = (std::conj(wt[0]) * wt[0]).real();
	err[1] = (std::conj(w[0]) * w[0] + std::conj(wt[1]) * wt[1]).real();
	err[2] = (std::conj(w[0]) * w[0] + std::conj(w[1] * w[1]) + std::conj(wt[2] * wt[2])).real();

	//err[0] = (wt[0] * wt[0]).real();
	//err[1] = (w[0] * w[0] + wt[1] * wt[1]).real();
	//err[2] = (w[0]*w[0]+w[1]*w[1]+wt[2]*wt[2]).real();
}

// Setup and calculate the 1st iteraction
void CG_SolverComplex::init(matrixcomplex &A_H)
{
	r.resize(Afull.rows()); //saveVals("Aprint.txt", Afull); saveVals("Atprint.txt", A_H); saveVals("bprint.txt", b);
	beta = 0;
	r = b - A_H*(Afull*x); //saveVals("rnorm2.txt", r.norm(), false);
	//p = r;

	z = r;
	precond.solveInPlace2(z); //saveVals("rprint.txt", r); saveVals("zprint.txt", z);
	rmod_1 = rmod = z.dot(r);
	p = z;

	//rmod_1 = rmod = r.squaredNorm();
	q = A_H*(Afull*p);
	gamma = rmod / q.dot(p);
	r0norm = r.norm(); //r0norm = sqrt(rmod);
	r0norm2 = r0norm*r0norm; //r0norm2 = rmod;
	alpha = std::complex<double>(1.0, 0.0) / gamma;

	// ########## LANCZOS
	/*v = r/r0norm;
	vt = Afull*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v;
	alpha_[0] = lalpha;
	std::cout << "Alpha by lanczos:" << lalpha << " by CG:" << alpha << std::endl;*/


	// Now do the first 3 iterations

	x += gamma*p;
	r = r - gamma*q; //saveVals("rnorm2.txt", r.norm());
	z = r;
	precond.solveInPlace2(z);
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod / rmod_1;

	eta_p1 = sqrt(beta) / gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = alpha;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	c = rt1 / r1;
	s = eta_p1 / r1;
	wt[0] = std::complex<double>(1.0, 0.0) / rt1;
	w[0] = std::complex<double>(1.0, 0.0) / r1;

	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A_H*(Afull*p);
	gamma_1 = gamma;
	gamma = rmod / q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1

	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = Afull*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[1] = lalpha;
	eta_[1] = leta;*/

	x += gamma*p;	// 2...
	r = r - gamma*q; //saveVals("rnorm2.txt", r.norm());
	z = r;
	precond.solveInPlace2(z);
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod / rmod_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c*alpha - s*eta;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	r2 = c*eta + s*alpha;	// r_2,2 = c_1*eta_2
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[1] = wt[1] = -r2*w[0];
	wt[1] /= rt1;
	w[1] /= r1;

	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A_H*(Afull*p);
	gamma_1 = gamma;
	gamma = rmod / q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = Afull*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[2] = lalpha;
	eta_[2] = leta;*/

	x += gamma*p;	// ...and 3!
	r = r - gamma*q; //saveVals("rnorm2.txt", r.norm());
	z = r;
	precond.solveInPlace2(z);
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);// saveVals("rmod.txt", rmod);
	beta = rmod / rmod_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c*alpha - s*c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	r2 = c_1*c*eta + s*alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1*eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[2] = wt[2] = -(r3*w[0] + r2*w[1]);
	wt[2] /= rt1;
	w[2] /= r1;

	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A_H*(Afull*p);
	gamma_1 = gamma;
	gamma = rmod / q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	it = 3;
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = Afull*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[3] = lalpha;
	eta_[3] = leta;*/

	err[0] = (std::conj(wt[0]) * wt[0]).real();
	err[1] = (std::conj(w[0]) * w[0] + std::conj(wt[1]) * wt[1]).real();
	err[2] = (std::conj(w[0]) * w[0] + std::conj(w[1] * w[1]) + std::conj(wt[2] * wt[2])).real();

	//err[0] = (wt[0] * wt[0]).real();
	//err[1] = (w[0] * w[0] + wt[1] * wt[1]).real();
	//err[2] = (w[0]*w[0]+w[1]*w[1]+wt[2]*wt[2]).real();
}

void CG_SolverComplex::do_iteration(matrixcomplex &A_H) {
	it++;

	x += gamma*p;
	if ((it % 500) == 0) {
		r = b - A_H*(Afull*x); //saveVals("rnorm2.txt", r.norm());
	}
	else {
		r = r - gamma*q; //saveVals("rnorm2.txt", r.norm());
	}
	z = r;
	precond.solveInPlace2(z);

	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod / rmod_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c*alpha - s*c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	r2 = c_1*c*eta + s*alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1*eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;


	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A_H*(Afull*p);
	gamma_1 = gamma;
	gamma = rmod / q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1

	/*// Brenziski:
	rA = Afull*(b - Afull*x);
	c0 = rA.squaredNorm();
	c1 = rA.dot(r);
	c2 = rA.squaredNorm();*/

	//std::cout << it << "=>" << r3 << " : " << r2 << " : " << r1 << "[" << rt1 << "]\n";

	// FIXME: GET RID OF THOSE STUPID BUFFERS!!!!!!!!!!!!!!!!!!!!!
	if (it<360) {
		w[it - 1] = wt[it - 1] = -(r3*w[it - 3] + r2*w[it - 2]);
		wt[it - 1] /= rt1;
		w[it - 1] /= r1;
	}
	else {
		w[0] = it;
	}

	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = Afull*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[it] = lalpha;
	eta_[it] = leta;*/

	//err[it - 1] = (w[it - 2] * w[it - 2] + wt[it - 1] * wt[it - 1] - wt[it - 2] * wt[it - 2]).real();
	err[it - 1] = (std::conj(w[it - 2]) * w[it - 2] + std::conj(wt[it - 1]) * wt[it - 1] - std::conj(wt[it - 2]) * wt[it - 2]).real();


	//std::cout << it << ":"  << x.squaredNorm() << std::endl;
}

void CG_SolverComplex::do_iteration() {

	it++;

	x += gamma*p;
	if((it % 500)==0) {
		r = b - A*x; //saveVals("rnorm.txt", r.norm());
	}
	else {
		r = r - gamma*q; //saveVals("rnorm.txt", r.norm());
	}
	z = r;
	precond.solveInPlace(z);

	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta)/gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c*alpha - s*c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	r2 = c_1*c*eta + s*alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1*eta;
	c_1 = c;
	c = rt1/r1;
	s_1 = s;
	s = eta_p1/r1;


	//p = r + beta*p;
	p = z + beta*p;
	q.noalias() = A*p;
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = std::complex<double>(1.0, 0.0) / gamma + beta / gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1

	/*// Brenziski:
	rA = A*(b - A*x);
	c0 = rA.squaredNorm();
	c1 = rA.dot(r);
	c2 = rA.squaredNorm();*/

	//std::cout << it << "=>" << r3 << " : " << r2 << " : " << r1 << "[" << rt1 << "]\n";

	// FIXME: GET RID OF THOSE STUPID BUFFERS!!!!!!!!!!!!!!!!!!!!!
	if(it<360) {
		w[it-1] = wt[it-1] = -(r3*w[it-3]+r2*w[it-2]);
		wt[it-1] /= rt1;
		w[it-1] /= r1;
	} else {
		w[0] = it;
	}

	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[it] = lalpha;
	eta_[it] = leta;*/

	//err[it - 1] = (w[it - 2] * w[it - 2] + wt[it - 1] * wt[it - 1] - wt[it - 2] * wt[it - 2]).real();
	err[it - 1] = (std::conj(w[it - 2]) * w[it - 2] + std::conj(wt[it - 1]) * wt[it - 1] - std::conj(wt[it - 2]) * wt[it - 2]).real();


	//std::cout << it << ":"  << x.squaredNorm() << std::endl;
}
//
//void CG_SolverComplex::saveVals(const char* fname, matrixcomplex &mat) {
//	std::ofstream myfile;
//	myfile.open(fname, std::ios::binary);
//	for (int i = 0; i < mat.rows(); i++) {
//		for (int j = 0; j < mat.cols(); j++) {
//			double valre = mat.coeff(j, i).real();
//			double valim = mat.coeff(j, i).imag();
//			myfile.write((char*)&valre, sizeof(double));
//			myfile.write((char*)&valim, sizeof(double));
//		}
//	}
//	myfile.close();
//}
//
//void CG_SolverComplex::saveVals(const char* fname, const Eigen::VectorXcd &vec) {
//	std::ofstream myfile;
//	myfile.open(fname, std::ios::binary); myfile;
//	for (int i = 0; i < vec.size(); i++) {
//		double valre = vec.coeff(i).real();
//		double valim = vec.coeff(i).imag();
//		myfile.write((char*)&valre, sizeof(double));
//		myfile.write((char*)&valim, sizeof(double));
//	}
//	myfile.close();
//}
//
//void CG_SolverComplex::saveVals(const char* fname, std::complex<double> &val, bool app) {
//	std::ofstream myfile;
//	if (app) myfile.open(fname, std::ios::binary | std::ofstream::app);
//	else myfile.open(fname, std::ios::binary);
//	double valre = val.real();
//	double valim = val.imag();
//	myfile.write((char*)&valre, sizeof(double));
//	myfile.write((char*)&valim, sizeof(double));
//	myfile.close();
//}
//
//void CG_SolverComplex::saveVals(const char* fname, double val, bool app) {
//	std::ofstream myfile;
//	if (app) myfile.open(fname, std::ios::binary | std::ofstream::app);
//	else myfile.open(fname, std::ios::binary);
//	myfile.write((char*)&val, sizeof(double));
//	myfile.close();
//}
