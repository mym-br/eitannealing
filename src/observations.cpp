#include "problem.h"

void problem::initObs(const char *filecurrents, const char* filename)
{
	std::ifstream file;
	std::ifstream filec;

	filec.open(filecurrents);
	file.open(filename);

	// FIXME: Read from file!
	nobs = 32;
	int n = getGenericElectrodesCount() - 1;

	tensions = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	Eigen::VectorXd current(getNodesCount() - 1);
	current.fill(0);
	int baseIndex = (int)current.size() - n;
	for (int i = 0; i<nobs; i++) {
		double c;
		int entry, exit;
		filec >> entry;
		filec >> exit;
		filec >> c;
		entry--; exit--;	// zero-based
		currents[i] = current;
		if (entry != n)
			currents[i][baseIndex + entry] = 1;
		if (exit != n)
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