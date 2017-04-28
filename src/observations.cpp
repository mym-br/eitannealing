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
		#ifndef BLOCKGND
		// Zero ground node current
		if (groundNode > 0 && groundNode < getNodesCount()) currents[i][groundNode] = 0;
		#endif

		// read tensions from file
		tensions[i].resize(getGenericElectrodesCount());
		double val, avg;
		avg = 0;
		for (int j = 0; j<getGenericElectrodesCount(); j++) {
			file >> val;
			tensions[i][j] = val / c;  // Values are normalized by current
			avg += tensions[i][j];
		}
		// rebase tensions, apply offset for zero sum on electrodes
		avg /= (double)getGenericElectrodesCount();
		for (int j = 0; j < getGenericElectrodesCount(); j++) {
			tensions[i][j] -= avg;
		}
	}

	filec.close();
	file.close();
}