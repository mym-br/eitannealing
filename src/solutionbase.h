#ifndef SOLUTIONBASE_H_
#define SOLUTIONBASE_H_

#include "solver.h"
#include "problem.h"
#include "random.h"
#include <memory>
#include <algorithm>

struct shuffleData {
	int ncoef;
	bool swap;
};

struct shuffler {
	int * shuffleConsts;
	int * swapshuffleconsts;
	shuffler(std::shared_ptr<problem> input, observations<double> *readings) {
		shuffleConsts = new int[input->getNumCoefficients()];
		swapshuffleconsts = new int[input->getInnerAdjacencyCount()];

		for (int i = 0; i<input->getNumCoefficients(); i++) shuffleConsts[i] = 0;
		for (int i = 0; i<input->getInnerAdjacencyCount(); i++) swapshuffleconsts[i] = 0;
	}
	shuffler(std::shared_ptr<problem> input, observations<std::complex<double>> *readings) {
		shuffleConsts = new int[2 * input->getNumCoefficients()];
		swapshuffleconsts = new int[2 * input->getInnerAdjacencyCount()];

		for (int i = 0; i<2 * input->getNumCoefficients(); i++) shuffleConsts[i] = 0;
		for (int i = 0; i<2 * input->getInnerAdjacencyCount(); i++) swapshuffleconsts[i] = 0;
	}
	~shuffler() {
		delete[] shuffleConsts;
		delete[] swapshuffleconsts;

	}
	void addShufflerFeedback(const shuffleData &data, bool pos);
};

template <class T>
class solutionbase {
	friend class solution;
	friend class solutioncomplex;
protected:
	//T* sol;
	//Eigen::SparseMatrix<T, Eigen::ColMajor> *stiffness, *stiffnessorig;

	Eigen::VectorXd distance;
	Eigen::VectorXd maxdist;
	Eigen::VectorXd mindist;
	Eigen::VectorXd err;
	Eigen::VectorXd err_x_dist;

	double totalDist;
	double minTotalDist;
	double maxTotalDist;
	int critical;
	double critErr;
	int totalit;

	static double calcNewShuffledValue(int ncoef, double curVal, shuffleData *data, const shuffler &sh, double minval, double maxval) {
		double val;
		if (sh.shuffleConsts[ncoef] == 0) {
			val = minval + genreal()*(maxval - minval);
		}
		else {
			do {
				val = curVal;
				double rnd = 0;
				for (int i = 0; i < sh.shuffleConsts[ncoef]; i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[ncoef];
				rnd -= 0.5;
				rnd *= (maxval - minval);
				val += rnd;
			} while ((val < minval) || (val > maxval));
		}
		if (data) {
			data->swap = false;
			data->ncoef = ncoef;
		}
		return val;
	}

	static std::pair<double, double> calcNewSwappedValue(int ncoef, int &node1, int &node2, double v1, double v2, shuffleData *data, const shuffler &sh, double minval, double maxval) {
		// Order nodes
		if (v1 > v2) {
			int aux = node1;
			node1 = node2;
			node2 = aux;
		}
		double a = std::max(std::min(v1 - minval, maxval - v2), std::min(maxval - v1, v2 - minval));

		std::pair<double, double> newVals(v1, v2);
		double delta;
		do {
			if (sh.swapshuffleconsts[ncoef] == 0) {
				delta = a*(genreal() * 2 - 1);
			}
			else {
				double rnd = 0;
				for (int i = 0; i < sh.swapshuffleconsts[ncoef]; i++)
					rnd += genreal();
				rnd /= sh.swapshuffleconsts[ncoef];
				rnd -= 0.5;
				delta = a*rnd;
			}
			newVals.first = v1 - delta;
			newVals.second = v2 + delta;
		} while ((v1 < minval) || (v2 < minval) || (v1 > maxval) || (v2 > maxval));
		if (data) {
			data->swap = true;
			data->ncoef = ncoef;
		}

		return newVals;
	}

	static T *copySolution(const T *sol, std::shared_ptr<problem> input)
	{
		T *res = new T[input->getNumCoefficients()];

		for (int i = 0; i<input->getNumCoefficients(); i++)
			res[i] = sol[i];

		return res;
	}
	//static matrix *getNewStiffness(double *sol, std::shared_ptr<problem> input);

	// shuffle constructor
	double regularisation;
	std::shared_ptr<problem> input;
	observations<T> *readings;
	//void zeroSumVector(Eigen::VectorXd &vec);
	
	void zeroSumVector(Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
		T avg = 0;
		for (int i = 0; i < vec.size(); i++) avg += vec[i];
		avg /= input->getGenericElectrodesCount();
		for (int i = 0; i < vec.size(); i++) vec[i] -= avg;
	}

	solutionbase(const T *sigma, std::shared_ptr<problem> _input, observations<T> *_readings) :
		//sol(solutionbase::copySolution(sigma, _input)),
		input(_input), readings(_readings),
		distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	solutionbase(std::shared_ptr<problem> _input, observations<T> *_readings) :
		//sol(solutionbase::getNewRandomSolution(_input, _readings)),
		input(_input), readings(_readings),
		distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	solutionbase(T *sigma, const solutionbase &base, std::shared_ptr<problem> _input, observations<T> *_readings) :
		//sol(sigma),
		input(_input), readings(_readings),
		distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};

	//solutionbase(const double *sigma, std::shared_ptr<problem> _input, observations<double> *_readings) :
	//	//sol(solutionbase::copySolution(sigma, _input)),
	//	distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	//solutionbase(std::shared_ptr<problem> _input, observations<double> *_readings) :
	//	//sol(solutionbase::getNewRandomSolution(_input, _readings)), 
	//	distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	//solutionbase(double *sigma, const solutionbase &base, observations<double> *_readings) :
	//	//sol(sigma),
	//	distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	virtual void improve() {}
	virtual void ensureMinIt(unsigned int it) {}
	virtual void ensureMaxE2(double e2) {}

public:
	virtual Eigen::Matrix<T, Eigen::Dynamic, 1> getSimulationX(int i) const = 0;

	virtual bool compareWith(solutionbase &target, double kt, double prob) = 0;
	virtual solutionbase *shuffle(shuffleData *data, const shuffler &sh) const = 0;
	//virtual void improve() = 0;

	static void saveMesh(double *sol, const char *filename, const char *propertyname, std::shared_ptr<problem> input, int step = 0) {
		std::ofstream myfile;
		myfile.open(filename);

		std::ifstream inputfile(input->getMeshFilename());
		for (int i = 0; inputfile.eof() != true; i++) {
			std::string line;
			std::getline(inputfile, line);
			myfile << line << '\n';
		}

		//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
		myfile << "$NodeData\n1\n\"" << propertyname << "\"\n1\n0.0\n3\n" << step << "\n1\n" << input->getNodesCount() << "\n";
		for (int j = 0; j < input->getNodesCount(); j++) {
			myfile << (j + 1) << "\t" << sol[input->getNode2Coefficient(j)] << "\n";
		}
		myfile << "$EndNodeData\n";
		myfile.flush();
		myfile.close();
	}

	//static void savePotentials(std::vector<Eigen::VectorXd> &sols, const char *filename, std::shared_ptr<problem> input, observations<double> *readings);

	virtual double getRegularisationValue() const { return this->regularisation; }
	double getDEstimate() const { return totalDist; }
	double getDMax() const { return maxTotalDist; }
	double getDMin() const { return minTotalDist; }
	int getCritical() const { return this->critical; }
	double getCritErr() const { return this->critErr; }
	double getErrorAt(int sim) const { return this->distance[sim]; }
	int getTotalIt() const { return this->totalit; }
	virtual  T* getSolution() = 0;// { return this->sol; }
	virtual ~solutionbase() {};

};

#endif