#include "problem.h"

void problem::initObs(const char *filecurrents, const char* filename)
{
	std::ifstream file;
	std::ifstream filec;

	filec.open(filecurrents);
	file.open(filename);

	int n = getGenericElectrodesCount() - 1;
	int valuesCount = std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());
	nobs = valuesCount / getGenericElectrodesCount();

	file.clear();
	file.seekg(0, std::ios::beg);
	tensions = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	currentVals = Eigen::VectorXd(nobs);
	//Eigen::VectorXd current(getNodesCount() - 1);
	Eigen::VectorXd current(getNodesCount());
	current.fill(0);
	int baseIndex = (int)current.size() - n - 1;
	for (int i = 0; i<nobs; i++) {
		double c;
		int entry, exit;
		filec >> entry;
		filec >> exit;
		filec >> c;
		entry--; exit--;	// zero-based
		currentVals[i] = c;
		currents[i] = current;
		currents[i][baseIndex + entry] = 1;
		currents[i][baseIndex + exit] = -1;

		// read tensions from file
		tensions[i].resize(getGenericElectrodesCount() - 1);
		double val;
		for (int j = 0; j<getGenericElectrodesCount() - 1; j++) {
			file >> val;
			tensions[i][j] = val / c;  // Values are normalized by current
		}
		// rebase tensions, as our last electrode is always the ground
		file >> val;
		for (int j = 0; j<getGenericElectrodesCount() - 1; j++) {
			tensions[i][j] -= val / c;
		}
	}

	filec.close();
	file.close();
}