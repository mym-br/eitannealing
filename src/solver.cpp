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

void checkedRegisterTimestamp(std::chrono::high_resolution_clock::time_point &t) {
#ifdef CGTIMING 
	t = std::chrono::high_resolution_clock::now();
#endif
}

CG_Solver::CG_Solver(matrix &_A, const Eigen::VectorXd &b, const SparseIncompleteLLT &pre, double res, bool init) :
	A(_A),
	b(b),
	x(Eigen::VectorXd::Zero(_A.rows())),
	precond(pre),
	totalItTime(0.0), totalTriangularTime(0.0), totalSpmvTime(0.0)
{
	if(init) this->init(res);
}

CG_Solver::CG_Solver(matrix &_A, const Eigen::VectorXd &b,  const SparseIncompleteLLT &pre, bool init) :
	A(_A),
	b(b), 
	x(Eigen::VectorXd::Zero(_A.rows())), 
	precond(pre),
	totalItTime(0.0), totalTriangularTime(0.0), totalSpmvTime(0.0)
{
	if(init) this->init();
}

CG_Solver::CG_Solver(matrix &_A, const Eigen::VectorXd &b, const Eigen::VectorXd &x0, const SparseIncompleteLLT &pre, bool init) :
	A(_A), b(b), x(x0),
	// Precond!
	precond(pre),
	totalItTime(0.0), totalTriangularTime(0.0), totalSpmvTime(0.0) 
{
	if(init) this->init();
}

// Setup and calculate the 1st iteraction	
void CG_Solver::init(double res)
{
	std::chrono::high_resolution_clock::time_point t1, t2, tri_t1, tri_t2, spmv_t1, spmv_t2;

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
	checkedRegisterTimestamp(t1);
	x += gamma*p;	// 1...
	if (rmod < res) { it = 1; return; }
	r = r - gamma*q;
	z = r;
	checkedRegisterTimestamp(tri_t1);
	precond.solveInPlace(z);
	checkedRegisterTimestamp(tri_t2);
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
	checkedRegisterTimestamp(spmv_t1);
	q.noalias() = A*p;
	checkedRegisterTimestamp(spmv_t2);
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	checkedRegisterTimestamp(t2);

#ifdef CGTIMING
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
	checkedRegisterTimestamp(t1);
	x += gamma*p;	// 2...
	if (rmod < res) { it = 2; return; }
	r = r - gamma*q;
	z = r;
	checkedRegisterTimestamp(tri_t1);
	precond.solveInPlace(z);
	checkedRegisterTimestamp(tri_t2);
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
	checkedRegisterTimestamp(spmv_t1);
	q.noalias() = A*p;
	checkedRegisterTimestamp(spmv_t2);
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	checkedRegisterTimestamp(t2);

#ifdef CGTIMING
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

	checkedRegisterTimestamp(t1);
	x += gamma*p;	// ...and 3!
	r = r - gamma*q;
	z = r;
	checkedRegisterTimestamp(tri_t1);
	precond.solveInPlace(z);
	checkedRegisterTimestamp(tri_t2);
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
	checkedRegisterTimestamp(spmv_t1);
	q.noalias() = A*p;
	checkedRegisterTimestamp(spmv_t2);
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	it = 3;
	checkedRegisterTimestamp(t2);
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
	totalItTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	totalTriangularTime += std::chrono::duration_cast<std::chrono::microseconds>(tri_t2 - tri_t1).count();
	totalSpmvTime += std::chrono::duration_cast<std::chrono::microseconds>(spmv_t2 - spmv_t1).count();
#endif // CGTIMING
}

void CG_Solver::do_iteration() {
	std::chrono::high_resolution_clock::time_point t1, t2, tri_t1, tri_t2, spmv_t1, spmv_t2;

	checkedRegisterTimestamp(t1);
	it++;

	x += gamma*p;
	if((it % 500)==0) {
		r = b - A*x;
	}
	else {
		r = r - gamma*q;
	}
	z = r;
	checkedRegisterTimestamp(tri_t1);
	precond.solveInPlace(z);
	checkedRegisterTimestamp(tri_t2);
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
	checkedRegisterTimestamp(spmv_t1);
	q.noalias() = A*p;
	checkedRegisterTimestamp(spmv_t2);
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	checkedRegisterTimestamp(t2);

#ifdef CGTIMING
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
    w[it-1] = wt[it-1] = -(r3*w[it-3]+r2*w[it-2]);
    wt[it-1] /= rt1;
    w[it-1] /= r1;

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
