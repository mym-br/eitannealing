/*
 * solver.cpp
 *
 *  Created on: Jun 27, 2010
 *      Author: thiago
 */

#include "solver.h"
#include "nodecoefficients.h"
#include "problemdescription.h"



CG_Solver::CG_Solver(matrix &A, Eigen::VectorXd &b,  const SparseIncompleteLLT &pre):
	A(A), b(b), x(Eigen::VectorXd::Zero(A.rows())),
	// Precond!
	precond(pre)	{
	this->init();
}

CG_Solver::CG_Solver(matrix &A, Eigen::VectorXd &b, const Eigen::VectorXd &x0, const SparseIncompleteLLT &pre):
	A(A), b(b), x(x0),
	// Precond!
	precond(pre)  {
	this->init();
}

// Setup and calculate the 1st iteraction	
void CG_Solver::init()
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
	r0norm = sqrt(rmod);
	alpha = 1/gamma;

	// ########## LANCZOS
	/*v = r/r0norm;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v;
	alpha_[0] = lalpha;
	std::cout << "Alpha by lanczos:" << lalpha << " by CG:" << alpha << std::endl;*/


	// Now do the first 3 iterations

	x += gamma*p;	// 1...
	r = r - gamma*q;
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
	wt[0] = 1/rt1;
	w[0] = 1/r1;

	//p = r + beta*p;
	p = z + beta*p;
	q = A*p.lazy();
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1

	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[1] = lalpha;
	eta_[1] = leta;*/

	x += gamma*p;	// 2...
	r = r - gamma*q;
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
	q = A*p.lazy();
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1
	/* ########## LANCZOS
	leta = vt.norm();
	v_1 = v;
	v = vt/leta;
	vt = A*v;
	lalpha = vt.dot(v);
	vt -= lalpha*v + leta*v_1;
	alpha_[2] = lalpha;
	eta_[2] = leta;*/

	x += gamma*p;	// ...and 3!
	r = r - gamma*q;
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
	w[2] = wt[2] = -(r3*w[0]+r2*w[1]);
	wt[2] /= rt1;
	w[2] /= r1;

	//p = r + beta*p;
	p = z + beta*p;
	q = A*p.lazy();
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

	err[0] = wt[0]*wt[0];
	err[1] = w[0]*w[0]+wt[1]*wt[1];
	err[2] = w[0]*w[0]+w[1]*w[1]+wt[2]*wt[2];
}
		
void CG_Solver::do_iteration() {

	it++;

	x += gamma*p;
	if((it % 500)==0) {
		r = b - A*x;
	}
	else {
		r = r - gamma*q;
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
	q = A*p.lazy();
	gamma_1 = gamma;
	gamma = rmod/q.dot(p);
	alpha = 1/gamma + beta/gamma_1;			// alpha_k+1 = 1/gamma_k + beta_k/gamma_k-1

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

	if(it<360)
		err[it-1] = w[it-2]*w[it-2]+wt[it-1]*wt[it-1] - wt[it-2]*wt[it-2];
	//err += wt[it-1]*wt[it-1] - wt[it-2]*wt[it-2] + w[it-2]*w[it-2];

	//std::cout << it << ":"  << x.squaredNorm() << std::endl;
}


/*matrix CG_Solver::buildJacobiMatrx()
{
	matrix result(it, it);
	result.startFill(2*it-1);
	int i;
	// Fill column-wise
	result.fill(0,0) = alpha_[0];
	result.fill(1,0) = eta_[1];
	for(i=1;i<it-1;i++) {
		result.fill(i-1,i) = eta_[i];
		result.fill(i,i) = alpha_[i];
		result.fill(i+1,i) = eta_[i+1];
	}
	result.fill(it-2,it-1) = eta_[it-1];
	result.fill(it-1,it-1) = alpha_[it-1];
	result.endFill();
	return result;
}*/

// FIXME: It should be slightly faster to have an "skeleton"
//		matrix (with ones for instance), copy it and iterate over its
//		coefficients.

void assembleProblemMatrix(float *cond, matrix **stiffnes)
{

	/*
	// Jacobi preconditioning
	*sqDiagonal = new Eigen::VectorXd(numNodes-1);
	Eigen::VectorXd &sqrtdiagonal = **sqDiagonal;
	int i;
	// Prepare diagonal square roots
	for(i=0;i<numNodes-1;i++) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux && aux->node < i) aux = aux->next;
		double val = 0;
		while(aux && aux->node==i) {
			val += aux->coefficient*cond[aux->condIndex];
			aux = aux->next;
		}
		sqrtdiagonal[i] = sqrt(val);
	}


	matrix *out = new matrix(numNodes-1, numNodes-1);
	out->startFill(3*(numNodes-1)); // estimate of the number of nonzeros (optional)
	for (i=0; i<numNodes-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux && aux->node <= i) aux = aux->next; // skip upper triangular
		// Diagonal is 1
		out->fill(i,i) = 1;
		while(aux) { // Col-major storage
			int row = aux->node;
			double val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			out->fill(row,i) = val/(sqrtdiagonal[row]*sqrtdiagonal[i]);
		}
	}
	out->endFill();

	*stiffnes = out;*/


	matrix *out = new matrix(nodes.size()-1, nodes.size()-1);
	double val;
	out->startFill(3*(nodes.size()-1)); // estimate of the number of nonzeros (optional)
	for (int i=0; i<nodes.size()-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux) { // Col-major storage
			while(aux->node < i) aux = aux->next; // skip upper triangular
			int row = aux->node;
			val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			out->fill(row,i) = val;
		}
	}
	out->endFill();

	*stiffnes = out;
}



void assembleProblemMatrix(float *cond, matrix **stiffnes, int numNodes, nodeCoefficients **nodeCoef)
{
	matrix *out = new matrix(numNodes-1, numNodes-1);
	double val;
	out->startFill(3*(numNodes-1)); // estimate of the number of nonzeros (optional)
	for (int i=0; i<numNodes-1; ++i) {
		nodeCoefficients *aux = nodeCoef[i];
		while(aux) { // Col-major storage
			while(aux->node < i) aux = aux->next; // skip upper triangular
			int row = aux->node;
			val = 0.0;
			while(aux && aux->node==row) {
				val += aux->coefficient*cond[aux->condIndex];
				aux = aux->next;
			}
			out->fill(row,i) = val;
		}
	}
	out->endFill();

	*stiffnes = out;
}
