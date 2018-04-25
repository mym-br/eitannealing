#include "solutioncuda.h"
#include "gradientnormregularisation.h"
#include <numeric>
#include <algorithm>

solutionCuda::solutionCuda(const double *sigma, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs) :
	sol(solutionCuda::copySolution(sigma, _input)),
	stiffness(solutionCuda::getNewStiffness(sol, _input)),
	precond(new SparseIncompleteLLT(*stiffness)),
	simulations(new CGCUDA_Solver *[_readings->getNObs()]),
	fixedCoeffs(_fixedCoeffs),
	solutionbase(sigma, _input, _readings)
{
	stiffnessCpjds = new MatrixCPJDS;
	mgr = CGCUDA_Solver::createManager(stiffness, stiffnessCpjds, input->getNodeCoefficients(), input->getNodesCount(), input->getNumCoefficients());
	lINFinityNorm = CGCUDA_Solver::createPreconditioner(*stiffnessCpjds, stiffnessCpjds->cpuData.data, stiffnessCpjds->cpuData.precond);
	this->initSimulations();
	this->initErrors();
}


// New random solution
solutionCuda::solutionCuda(std::shared_ptr<problem> _input, observations<double> *_readings, std::vector<double> &electrodesCoeffs) :
	sol(solutionCuda::getNewRandomSolution(_input, electrodesCoeffs)),
	stiffness(solutionCuda::getNewStiffness(sol, _input)),
	precond(new SparseIncompleteLLT(*stiffness)),
	simulations(new CGCUDA_Solver *[_readings->getNObs()]),
	fixedCoeffs(electrodesCoeffs.size()),
	solutionbase(_input, _readings)
{
	stiffnessCpjds = new MatrixCPJDS;
	mgr = CGCUDA_Solver::createManager(stiffness, stiffnessCpjds, input->getNodeCoefficients(), input->getNodesCount(), input->getNumCoefficients());
	lINFinityNorm = CGCUDA_Solver::createPreconditioner(*stiffnessCpjds, stiffnessCpjds->cpuData.data, stiffnessCpjds->cpuData.precond);
	this->initSimulations();
	this->initErrors();
}

// New randomly modified solution
solutionCuda::solutionCuda(double *sigma, const solutionCuda &base, std::shared_ptr<problem> _input, observations<double> *_readings, int _fixedCoeffs) :
	sol(sigma),
	stiffness(solutionCuda::getNewStiffness(sol, _input)),
	precond(new SparseIncompleteLLT(*stiffness)),
	simulations(new CGCUDA_Solver *[_readings->getNObs()]),
	fixedCoeffs(_fixedCoeffs),
	solutionbase(sigma, base, _input, _readings)
{
	stiffnessCpjds = new MatrixCPJDS;
	mgr = CGCUDA_Solver::createManager(stiffness, stiffnessCpjds, input->getNodeCoefficients(), input->getNodesCount(), input->getNumCoefficients());
	lINFinityNorm = CGCUDA_Solver::createPreconditioner(*stiffnessCpjds, stiffnessCpjds->cpuData.data, stiffnessCpjds->cpuData.precond);
	this->initSimulations(base);
	this->initErrors();
}

void solutionCuda::initSimulations() {
	// Prepare solvers
	int i;
	this->totalit = 0;
	for (i = 0; i<readings->getNObs(); i++)
	{
		numType *currentsData = new numType[input->getNodesCount()];
		for (int j = 0; j < input->getNodesCount(); j++) currentsData[j] = input->getCurrentVector(i, readings)[j];
		Vector *bVec = CGCUDA_Solver::createCurrentVector(currentsData, *mgr, stiffnessCpjds->matrixData.n, input->getNodesCount());
		simulations[i] = new CGCUDA_Solver(stiffnessCpjds, mgr, bVec, lINFinityNorm);

		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while (ndecr<2) {
			simulations[i]->do_iteration();
			aux = simulations[i]->getErrorl2Estimate();
			if (aux >= err) ndecr = 0;
			else {
				ndecr++;
			}
			err = aux;
		}
		this->totalit += simulations[i]->getIteration();
	}
}

void solutionCuda::initSimulations(const solutionCuda &base) {
	// Prepare solvers
	int i;
	this->totalit = 0;
	for (i = 0; i<readings->getNObs(); i++)
	{
		numType *currentsData = new numType[input->getNodesCount()];
		for (int j = 0; j < input->getNodesCount(); j++) currentsData[j] = input->getCurrentVector(i, readings)[j];
		Vector *bVec = CGCUDA_Solver::createCurrentVector(currentsData, *mgr, stiffnessCpjds->matrixData.n, input->getNodesCount());
		// Reuse previous solutions as initial values
		//simulations[i] = new CGCUDA_Solver(stiffnessCpjds, mgr, bVec, base.simulations[i]->getCpjdsX(), lINFinityNorm);
		simulations[i] = new CGCUDA_Solver(stiffnessCpjds, mgr, bVec, lINFinityNorm);

		// Run three iterations, then wait for 3 consecutive decreasing error estimates
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		//simulations[i]->do_iteration();
		double err = simulations[i]->getErrorl2Estimate();
		double aux;
		int ndecr = 0;
		while (ndecr<2) {
			simulations[i]->do_iteration();
			aux = simulations[i]->getErrorl2Estimate();
			if (aux >= err) ndecr = 0;
			else {
				ndecr++;
			}
			err = aux;
		}
		this->totalit += simulations[i]->getIteration();
	}
}

void solutionCuda::initErrors() {
	// Calc regularisation value
	this->regularisation = gradientNormRegularisation::getInstance()->getRegularisation(this->sol)*input->regularizationFactor;
	// Calc electrode contact condutivity variance
	//for (int i = 0; i < input->getGenericElectrodesCount(); i++) { std::cout << this->sol[i] << " "; std::cout << std::endl; }
	if (std::abs(input->electrodevar) > 1e-6) {
		double sum = std::accumulate(this->sol, this->sol + input->getGenericElectrodesCount(), 0.0);
		double mean = sum / input->getGenericElectrodesCount();
		std::vector<double> diff(input->getGenericElectrodesCount());
		std::transform(this->sol, this->sol + input->getGenericElectrodesCount(), diff.begin(), [mean](double x) { return x - mean; });
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		double var = sq_sum / input->getGenericElectrodesCount();
		this->elecvariance = var*input->electrodevar;
		this->regularisation += this->elecvariance;
	}
	else 
		this->elecvariance = 0;
	int i;
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXd aux(input->getGenericElectrodesCount());
	// Retrieve distance estimates, errors and boundaries
	for (i = 0; i<readings->getNObs(); i++) {
		// Compare with observation
		Eigen::VectorXd xcudavec = (simulations[i]->getX()).cast<double>();
		aux = xcudavec.tail(aux.size());
		#ifndef BLOCKGND
		// Rebase tension for zero sum
		zeroSumVector(aux);
		#endif
		aux -= readings->getTensions()[i];

		distance[i] = aux.norm();
		err[i] = sqrt(simulations[i]->getErrorl2Estimate());
		double distancei = distance[i];
		double erri = err[i];
		maxdist[i] = distance[i] + err[i];
		mindist[i] = std::max(distance[i] - err[i],0.0);
		err_x_dist[i] = maxdist[i]*err[i];
		

	}
	totalDist = distance.norm()+regularisation;
	minTotalDist = mindist.norm()+regularisation;
	maxTotalDist = maxdist.norm()+regularisation;
	// evaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for (i = 1; i<readings->getNObs(); i++) {
		if(max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}

void solutionCuda::improve()
{
	// Just some scrap space to avoid dynamic allocations
	//		WARNING: Obviously thread-unsafe!!!!
	static Eigen::VectorXd aux(input->getGenericElectrodesCount());

	// Do another iteration on the critical solver
	simulations[critical]->do_iteration();
	this->totalit++;
	// Recalcule expected distance and boundaries
	Eigen::VectorXd xcudavec = (simulations[critical]->getX()).cast<double>();
	aux = (simulations[critical]->getX()).cast<double>().tail(input->getGenericElectrodesCount());
	#ifndef BLOCKGND
	// Rebase tension for zero sum
	zeroSumVector(aux);
	#endif
	aux -= readings->getTensions()[critical];

	distance[critical] = aux.norm();
	err[critical] = sqrt(simulations[critical]->getErrorl2Estimate());
	maxdist[critical] = distance[critical] + err[critical];
	mindist[critical] = std::max(distance[critical] - err[critical], 0.0);
	err_x_dist[critical] = maxdist[critical] * err[critical];
	totalDist = distance.norm() + regularisation;
	minTotalDist = mindist.norm() + regularisation;
	maxTotalDist = maxdist.norm() + regularisation;
	// reevaluate critical
	double max = err_x_dist[0];
	critical = 0;
	for (int i = 1; i<readings->getNObs(); i++) {
		if (max < err_x_dist[i]) {
			max = err_x_dist[i];
			critical = i;
		}
	}
	critErr = err[critical];
}

bool solutionCuda::compareWith(solutionbase &target, double kt, double prob)
{
	double delta, expdelta;
	// Ensure errors are within required margin
	while (true) {
		double min_delta = target.minTotalDist - this->maxTotalDist;
		double max_delta = target.maxTotalDist - this->minTotalDist;
		delta = target.totalDist - this->totalDist;
		expdelta = exp(-delta / kt);
		// Check boundary conditions:
		// Upper bound is negative, no doubt here
		if (max_delta < 0) break;
		// Estimate is negative, but upper bound is positive
		else if (delta <= 0) {
			if (exp(-max_delta / kt) >= (1 - prob)) break; // Upper condition only
		}
		// Estimate and upper bounds are positive, lower bound is negative
		else if (min_delta <= 0) {
			if (expdelta >= 1 - prob) { // Lower condition
				if (expdelta <= prob) break;	// upper condition
				if (exp(-max_delta / kt) >= (expdelta - prob)) break;
			}
		}
		// Classic case, everything is positive
		else {
			if (exp(-min_delta / kt) <= prob + expdelta) { // lower condition
				if (expdelta <= prob) break;	// upper condition
				if (exp(-max_delta / kt) >= (expdelta - prob)) break;
			}
		}
		// Not there yet, improve boundaries
		// Select wich one to improve
		if (this->critErr > target.critErr)
			this->improve();
		else
			target.improve();
	}
	if (delta <= 0) {
		//std::cout << "+";
		return true;
	}
	//std::cout << "-" << expdelta << std::endl;


	// Now randomly accept based on the energy
	if (genreal()<expdelta) return true;
	return false;
}

double *solutionCuda::getNewRandomSolution(std::shared_ptr<problem> input, std::vector<double> &electrodesCoeffs)
{
	double *res = new double[input->getNumCoefficients()];
	int i = 0;
	for (i = 0; i < electrodesCoeffs.size(); i++)
		res[i] = electrodesCoeffs[i];
	for (; i<input->getNumCoefficients(); i++)
		res[i] = mincond + genreal()*(maxcond - mincond);

	return res;
}

solutionCuda *solutionCuda::shuffle(shuffleData *data, const shuffler &sh) const
{
	double *sigma = getShuffledSolution(data, sh);
	solutionCuda *res;
	try {
		res = new solutionCuda(sigma, *this, input, readings, fixedCoeffs);
	}
	catch (...) {
		for (int i = 0; i<65; i++)
			std::cout << i << ":" << sigma[i] << std::endl;
		exit(0);
	}
	return res;
}

double *solutionCuda::getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	double *res = solutionbase<double>::copySolution(sol, input);
	// head or tails
	if (genint(2)) { // Normal
		int ncoef = fixedCoeffs + genint(input->getNumCoefficients() - fixedCoeffs);	// Lower values fixed;

		if (sh.shuffleConsts[ncoef] == 0) {
			res[ncoef] = mincond + genreal()*(maxcond - mincond);
		}
		else {
			double val;
			do {
				val = res[ncoef];
				double rnd = 0;
				for (int i = 0; i<sh.shuffleConsts[ncoef]; i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoef];
				rnd -= 0.5;
				rnd *= (maxcond - mincond);
				val += rnd;
			} while ((val < mincond) || (val > maxcond));
			res[ncoef] = val;
		}
		if (data) {
			data->swap = false;
			data->ncoef = ncoef;
		}
	}
	else { // swap
		int ncoef = genint(input->getInnerAdjacencyCount());
		int node1, node2;

		node1 = input->node2coefficient[input->innerAdjacency[ncoef].first];
		node2 = input->node2coefficient[input->innerAdjacency[ncoef].second];

		// Order nodes
		if (res[node1]>res[node2]) {
			int aux = node1;
			node1 = node2;;
			node2 = aux;
		}
		double v1 = res[node1], v2 = res[node2];
		double a = std::max(std::min(v1 - mincond, maxcond - v2), std::min(maxcond - v1, v2 - mincond));

		double delta;
		do {
			if (sh.swapshuffleconsts[ncoef] == 0) {
				delta = a * (genreal() * 2 - 1);
			}
			else {
				double rnd = 0;
				for (int i = 0; i<sh.swapshuffleconsts[ncoef]; i++)
					rnd += genreal();
				rnd /= sh.swapshuffleconsts[ncoef];
				rnd -= 0.5;
				delta = a * rnd;
			}
			v1 = res[node1] - delta;
			v2 = res[node2] + delta;
		} while ((v1 < mincond) || (v2 < mincond) || (v1 > maxcond) || (v2 > maxcond));
		res[node1] = v1;
		res[node2] = v2;
		if (data) {
			data->swap = true;
			data->ncoef = ncoef;
		}
	}
	return res;
}

solutionCuda::~solutionCuda()
{
	delete[] sol;
	delete stiffness;
	delete precond;
	for (int i = 0; i<readings->getNObs(); i++) {
		delete simulations[i];
	}
	delete[] simulations;
	delete stiffnessCpjds;
	delete mgr;
}