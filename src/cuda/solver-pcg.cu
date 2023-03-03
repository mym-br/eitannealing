#include "solver-pcg.h"

#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <iomanip>

#include "settings.h"
#include "utils.h"
//#include "preconditioner.h"

void PCGSolverCPJDS::checkedCudaEventRecord(cudaEvent_t &event) {
#ifdef CGTIMING 
	cudaEventRecord(event);
#endif
}

using namespace cgl;

PCGSolverCPJDS::PCGSolverCPJDS(MatrixCPJDSManager * mgr, MatrixCPJDS *M, Vector * b) {
	this->mgr = mgr;
	this->A = M;
	this->b = b;
	this->x = new Vector(M->matrixData.n);

	r = new Vector(M->matrixData.n);
	z = new Vector(M->matrixData.n);
	p = new Vector(M->matrixData.n);
	q = new Vector(M->matrixData.n);
	u = new Vector(M->matrixData.n);
	partial = new Vector(M->matrixData.n);

	rmod = new Number(0);
	rmod_prev = new Number(1);
	rmod_aux = new Number(1);
	gamma = new Number(0);

	size = A->matrixData.n;
	blocks = ceil((double)size / BLOCKSIZE);
	partial2 = new Vector(this->blocks);

	totalItTime = totalTriangularTime = totalSpmvTime = 0;
	streamInit();
}

void PCGSolverCPJDS::init(double res) {
	x->reset(stream);
	this->init(x, res);
}

// Setup and calculate the 1st iteration	
void PCGSolverCPJDS::init(Vector *x0, double res) {
	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
#endif // CGTIMING

	double *data_h = new double[1];
	//m_preconditioner(A, PCGSolverCPJDS::streams[PCGSolverCPJDS::mainStream]);
	x->reset();
	r->reset();
	z->reset();
	p->reset();
	q->reset();
	u->reset();
	partial->reset();

	double * zData = z->getData();
	double * rData = r->getData();
	double * xData = x->getData();
	double * bData = b->getData();
	double * partialData = partial->getData();

	it = 0;

	// Initial x value
	x0->copyTo(x);

	//r = b - A * x;
	mgr->mult(*A, x, u); // u = A * x
	b->subtr(u, r); // r = b - u = b - A * x
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	z->copyTo(p); //p = z;
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	mgr->mult(*A, p, q); // q = A * p;
	p->inner(q, gamma); // gamma = q.dot(p);
	#ifdef CALCULATE_ERRORS
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	rmod2_1 = rmod2 = *data_h;
	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;
	if (rmod2 < res) { it = 0; return; }

	// Error calculations
	r0norm2 = rmod2;
	r0norm = sqrt(r0norm2);
	beta = 0;
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2;
	#endif

	// Now do the first 3 iterations
	checkedCudaEventRecord(startTotal);
	x->scalarAdd(rmod, gamma, p, NULL); // 1...
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = alpha;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	c = rt1 / r1;
	s = eta_p1 / r1;
	wt[0] = 1 / rt1;
	w[0] = 1 / r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;
	if (rmod2 < res) { it = 1; return; }

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING

	checkedCudaEventRecord(startTotal);
	x->scalarAdd(rmod, gamma, p, NULL); // 2...
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = c * alpha - s * eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c * eta + s * alpha;	// r_2,2 = c_1*eta_2
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[1] = wt[1] = -r2 * w[0];
	wt[1] /= rt1;
	w[1] /= r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;
	if (rmod2 < res) { it = 2; return; }

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING

	checkedCudaEventRecord(startTotal);
	x->scalarAdd(rmod, gamma, p, NULL); // 3...
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2;
	rt1 = c * alpha - s * c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c_1 * c*eta + s * alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1 * eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;
	w[2] = wt[2] = -(r3*w[0] + r2 * w[1]);
	wt[2] /= rt1;
	w[2] /= r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif //CGTIMING

#ifdef CALCULATE_ERRORS
	err[0] = wt[0] * wt[0];
	err[1] = w[0] * w[0] + wt[1] * wt[1];
	err[2] = w[0] * w[0] + w[1] * w[1] + wt[2] * wt[2];
#endif
	it = 3;

	delete data_h;
}

void PCGSolverCPJDS::doIteration(int iteration) {
	cudaEvent_t startTotal, stopTotal, startTri, stopTri, startSpmv, stopSpmv;
#ifdef CGTIMING
	cudaEventCreate(&startTotal); cudaEventCreate(&stopTotal);
	cudaEventCreate(&startTri); cudaEventCreate(&stopTri);
	cudaEventCreate(&startSpmv); cudaEventCreate(&stopSpmv);
#endif // CGTIMING

	checkedCudaEventRecord(startTotal);
	double *data_h = new double[1];
	it++;
	x->scalarAdd(rmod, gamma, p, NULL);
	r->scalarSubtr(rmod, gamma, q, NULL); //r -= alpha * q; alpha = rmod / gamma
	// M * z = r -> solve for z, having M = L * Lt
	// L * Lt * z = r
	checkedCudaEventRecord(startTri);
	mgr->solve(*A, r, u); // solving L * u = r
	mgr->solve_t(*A, u, z); // now solve Lt * z = u
	checkedCudaEventRecord(stopTri);
	std::swap(rmod_prev, rmod); //rmod_1 = rmod;
	z->inner(r, rmod); // rmod = z.dot(r);
	p->scalar(rmod, rmod_prev, u); // p = z + beta * p; beta = rmod / rmod_prev
	u->sum(z, p);
	r->inner(z, rmod); // rmod = z.dot(r); isto nao deveria ser recalculado
	checkedCudaEventRecord(startSpmv);
	mgr->mult(*A, p, q); // q = A * p;
	checkedCudaEventRecord(stopSpmv);
	p->inner(q, gamma); // gamma = q.dot(p);

#ifdef CALCULATE_ERRORS
	rmod2_1 = rmod2;
	cudaMemcpy(data_h, rmod->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	rmod2 = *data_h;

	// Error calculations
	beta = rmod2 / rmod2_1;
	eta = eta_p1;
	eta_p1 = sqrt(beta) / gamma2; 	// eta_k+1 = sqrt(beta_k)/gamma_k-1
	rt1 = c * alpha - s * c_1*eta;
	r1 = sqrt(rt1*rt1 + eta_p1 * eta_p1);
	r2 = c_1 * c*eta + s * alpha; // r_2,k = c_k-2*c_k-1*eta_k + s_k-1*alpha_k
	r3 = s_1 * eta;
	c_1 = c;
	c = rt1 / r1;
	s_1 = s;
	s = eta_p1 / r1;

	gamma2_1 = gamma2;
	cudaMemcpy(data_h, gamma->getData(), (size_t)1 * sizeof(double), cudaMemcpyDeviceToHost);
	gamma2 = *data_h;

	// FIXME: GET RID OF THOSE STUPID BUFFERS!!!!!!!!!!!!!!!!!!!!!
	if (it<360) {
		w[it - 1] = wt[it - 1] = -(r3*w[it - 3] + r2 * w[it - 2]);
		wt[it - 1] /= rt1;
		w[it - 1] /= r1;
	}
	else w[0] = it;
	err[it - 1] = w[it - 2] * w[it - 2] + wt[it - 1] * wt[it - 1] - wt[it - 2] * wt[it - 2];

	// Update values for next iteration
	gamma2 = rmod2 / gamma2;
	alpha = 1 / gamma2 + beta / gamma2_1;
#endif
	delete data_h;
	checkedCudaEventRecord(stopTotal);
#ifdef CGTIMING
	cudaEventSynchronize(stopTotal); cudaEventSynchronize(stopTri); cudaEventSynchronize(stopSpmv);
	float msTotal, msTri, msSpmv;  msTotal = msTri = msSpmv = 0;
	cudaEventElapsedTime(&msTotal, startTotal, stopTotal); cudaEventElapsedTime(&msTri, startTri, stopTri); cudaEventElapsedTime(&msSpmv, startSpmv, stopSpmv);
	totalItTime += (float)(1e3 * msTotal);
	totalTriangularTime += (float)(1e3 * msTri);
	totalSpmvTime += (float)(1e3 *  msSpmv);
#endif // CGTIMING
}

void PCGSolverCPJDS::streamInit() {
	cudaStreamCreate(&this->stream);
}

void PCGSolverCPJDS::streamDestroy() {
	cudaStreamDestroy(this->stream);
}