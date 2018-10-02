/*
 * solver.cpp
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */


#include "solver.h"
//#include "nodecoefficients.h"
#include <fstream>
#include <chrono>

CG_Solver::CG_Solver(matrix &_A, Eigen::VectorXd &b, const SparseIncompleteLLT &pre, double res) :
	A(_A),
	b(b),
	x(Eigen::VectorXd::Zero(_A.rows())),
	precond(pre),
	totalItTime(0.0), totalTriangularTime(0.0), totalSpmvTime(0.0)
{
	this->init(res);
}

CG_Solver::CG_Solver(matrix &_A, Eigen::VectorXd &b,  const SparseIncompleteLLT &pre):
	A(_A),
	b(b), 
	x(Eigen::VectorXd::Zero(_A.rows())), 
	precond(pre),
	totalItTime(0.0), totalTriangularTime(0.0), totalSpmvTime(0.0)
{
	this->init();
}

CG_Solver::CG_Solver(matrix &A_, Eigen::VectorXd &b, const Eigen::VectorXd &x0, const SparseIncompleteLLT &pre):
	A(A_), b(b), x(x0),
	// Precond!
	precond(pre),
	totalItTime(0.0), totalTriangularTime(0.0), totalSpmvTime(0.0) 
{
	this->init();
}

// Setup and calculate the 1st iteraction	
void CG_Solver::init(double res)
{
	r.resize(A.rows());
	beta = 0;
	r = b - A*x;
	//p = r;

	z = r;
	precond.solveInPlace(z);
	rmod_1 = rmod = z.dot(r);
	p = z;

	//rmod_1 = rmod = r.squaredNorm();
	q = A*p;
	gamma = rmod/q.dot(p);
	r0norm2 = rmod;
	r0norm = sqrt(rmod);
	alpha = 1/gamma;
	
	if (rmod < res) { it = 0; return; }
	// ########## LANCZOS
	/*v = r/r0norm;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v;
	alpha_[0] = lalpha;
	std::cout << "Alpha by lanczos:" << lalpha << " by CG:" << alpha << std::endl;*/


	// Now do the first 3 iterations
#ifdef CGTIMING
	auto t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	x += gamma*p;	// 1...
	if (rmod < res) { it = 1; return; }
	r = r - gamma*q;
	z = r;
#ifdef CGTIMING
	auto tri_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	precond.solveInPlace(z);
#ifdef CGTIMING
	auto tri_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;
#ifdef CALCULATE_ERRORS
	eta_p1 = sqrt(beta)/gamma; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = alpha;
	r1 = sqrt(rt1*rt1 + eta_p1*eta_p1);
	c = rt1/r1;
	s = eta_p1/r1;
	wt[0] = 1/rt1;
	w[0] = 1/r1;
#endif // CALCULATE_ERRORS
	//p = r + beta*p;
	p = z + beta*p;
#ifdef CGTIMING
	auto spmv_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	q.noalias() = A*p;
#ifdef CGTIMING
	auto spmv_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
#ifdef CGTIMING
	auto t2 = std::chrono::high_resolution_clock::now();
	totalItTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	totalTriangularTime += std::chrono::duration_cast<std::chrono::microseconds>(tri_t2 - tri_t1).count();
	totalSpmvTime += std::chrono::duration_cast<std::chrono::microseconds>(spmv_t2 - spmv_t1).count();
#endif // CGTIMING
	if (rmod < res) { it = 1; return; }
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[1] = lalpha;
	eta_[1] = leta;*/
#ifdef CGTIMING
	t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	x += gamma*p;	// 2...
	if (rmod < res) { it = 2; return; }
	r = r - gamma*q;
	z = r;
#ifdef CGTIMING
	tri_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	precond.solveInPlace(z);
#ifdef CGTIMING
	tri_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;
#ifdef CALCULATE_ERRORS
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
#endif // CALCULATE_ERRORS

	//p = r + beta*p;
	p = z + beta*p;
#ifdef CGTIMING
	spmv_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	q.noalias() = A*p;
#ifdef CGTIMING
	spmv_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
#ifdef CGTIMING
	t2 = std::chrono::high_resolution_clock::now();
	totalItTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	totalTriangularTime += std::chrono::duration_cast<std::chrono::microseconds>(tri_t2 - tri_t1).count();
	totalSpmvTime += std::chrono::duration_cast<std::chrono::microseconds>(spmv_t2 - spmv_t1).count();
#endif // CGTIMING
	if (rmod < res) { it = 2; return; }

	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[2] = lalpha;
	eta_[2] = leta;*/

#ifdef CGTIMING
	t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	x += gamma*p;	// ...and 3!
	r = r - gamma*q;
	z = r;
#ifdef CGTIMING
	tri_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	precond.solveInPlace(z);
#ifdef CGTIMING
	tri_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;
#ifdef CALCULATE_ERRORS
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
#endif // CALCULATE_ERRORS
	//p = r + beta*p;
	p = z + beta*p;
#ifdef CGTIMING
	spmv_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	q.noalias() = A*p;
#ifdef CGTIMING
	spmv_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	it = 3;
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[3] = lalpha;
	eta_[3] = leta;*/
#ifdef CALCULATE_ERRORS
	err[0] = wt[0]*wt[0];
	err[1] = w[0]*w[0]+wt[1]*wt[1];
	err[2] = w[0]*w[0]+w[1]*w[1]+wt[2]*wt[2];
#endif // CALCULATE_ERRORS
#ifdef CGTIMING
	t2 = std::chrono::high_resolution_clock::now();
	totalItTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	totalTriangularTime += std::chrono::duration_cast<std::chrono::microseconds>(tri_t2 - tri_t1).count();
	totalSpmvTime += std::chrono::duration_cast<std::chrono::microseconds>(spmv_t2 - spmv_t1).count();
#endif // CGTIMING
}
		
void CG_Solver::do_iteration() {
	auto t1 = std::chrono::high_resolution_clock::now();
	it++;

	x += gamma*p;
	if((it % 500)==0) {
		r = b - A*x;
	}
	else {
		r = r - gamma*q;
	}
	z = r;
#ifdef CGTIMING
	auto tri_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	precond.solveInPlace(z);
#ifdef CGTIMING
	auto tri_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	rmod_1 = rmod;
	//rmod = r.squaredNorm();
	rmod = z.dot(r);
	beta = rmod/rmod_1;
#ifdef CALCULATE_ERRORS
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
#endif // CALCULATE_ERRORS

	//p = r + beta*p;
	p = z + beta*p;
#ifdef CGTIMING
	auto spmv_t1 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	q.noalias() = A*p;
#ifdef CGTIMING
	auto spmv_t2 = std::chrono::high_resolution_clock::now();
#endif // CGTIMING
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
#ifdef CGTIMING
	auto t2 = std::chrono::high_resolution_clock::now();
	totalItTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	totalTriangularTime += std::chrono::duration_cast<std::chrono::microseconds>(tri_t2 - tri_t1).count();
	totalSpmvTime += std::chrono::duration_cast<std::chrono::microseconds>(spmv_t2 - spmv_t1).count();
#endif // CGTIMING
	/*// Brenziski:
	rA = A*(b - A*x);
	c0 = rA.squaredNorm();
	c1 = rA.dot(r);
	c2 = rA.squaredNorm();*/

	//std::cout << it << "=>" << r3 << " : " << r2 << " : " << r1 << "[" << rt1 << "]\n";

#ifdef CALCULATE_ERRORS
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

	err[it-1] = w[it-2]*w[it-2]+wt[it-1]*wt[it-1] - wt[it-2]*wt[it-2];
#endif // CALCULATE_ERRORS
	//std::cout << it << ":"  << x.squaredNorm() << std::endl;
}

void CG_Solver::saveVals(const char* fname, double val, bool app) {
	std::ofstream myfile;
	if (app) myfile.open(fname, std::ios::binary | std::ofstream::app);
	else myfile.open(fname, std::ios::binary);
	myfile.write((char*)&val, sizeof(double));
	myfile.close();
}