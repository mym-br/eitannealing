#ifndef SOLUTIONBASE_H_
#define SOLUTIONBASE_H_

#include "solver.h"
#include "problem.h"
#include "random.h"
#include <memory>

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
protected:
	//T* sol;
	//Eigen::SparseMatrix<T, Eigen::ColMajor> *stiffness, *stiffnessorig;

	Eigen::VectorXd distance;
	Eigen::VectorXd maxdist;
	Eigen::VectorXd mindist;
	Eigen::VectorXd err;
	Eigen::VectorXd err_x_dist;

	double minTotalDist;
	double maxTotalDist;
	int critical;
	double critErr;

	int totalit;

	//virtual void initSimulations() = 0;
	//virtual void initSimulations(const solutionbase &base) = 0;
	//virtual void initErrors() = 0;
	//virtual T *getShuffledSolution(shuffleData *data, const shuffler &sh) const = 0;

	//static double *getNewRandomSolution(std::shared_ptr<problem> input, observations<double> *readings)
	//{
	//	double *res = new double[input->getNumCoefficients()];
	//	for (int i = 0; i<input->getNumCoefficients(); i++) res[i] = mincond + genreal()*(maxcond - mincond);
	//	return res;
	//}

	//static std::complex<double> *getNewRandomSolution(std::shared_ptr<problem> input, observations<std::complex<double>> *readings)
	//{
	//	std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];
	//	double w = 2 * M_PI * input->getCurrentFreq();
	//	double wminperm = w*minperm, wmaxperm = w*maxperm;
	//	for (int i = 0; i < input->getNumCoefficients(); i++)
	//		res[i] = std::complex<double>(mincond + genreal()*(maxcond - mincond), wminperm + genreal()*(wmaxperm - wminperm));

	//	return res;
	//}

	static T *copySolution(const std::complex<double> *sol, std::shared_ptr<problem> input)
	{
		std::complex<double> *res = new std::complex<double>[input->getNumCoefficients()];

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

	solutionbase(const std::complex<double> *sigma, std::shared_ptr<problem> _input, observations<T> *_readings) :
		//sol(solutionbase::copySolution(sigma, _input)),
		input(_input), readings(_readings),
		distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	solutionbase(std::shared_ptr<problem> _input, observations<T> *_readings) :
		//sol(solutionbase::getNewRandomSolution(_input, _readings)),
		input(_input), readings(_readings),
		distance(_readings->getNObs()), maxdist(_readings->getNObs()), mindist(_readings->getNObs()), err(_readings->getNObs()), err_x_dist(_readings->getNObs()) {};
	solutionbase(std::complex<double> *sigma, const solutionbase &base, std::shared_ptr<problem> _input, observations<T> *_readings) :
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
	

public:
	double totalDist;
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