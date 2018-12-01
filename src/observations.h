/*
* observations.h
*
*  Created on: Jun 14, 2017
*      Author: aksato
*/
#ifndef OBSERVATIONS_H_
#define OBSERVATIONS_H_

#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include "problem.h"
#include "basematrix.h"

template <typename _Scalar>
class observations {
private:
	int nobs;
	std::string mesh_file;
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *tensions;
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *currents;
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> currentVals;

public:
	observations() : nobs(-1), tensions(nullptr), currents(nullptr) {}
	~observations() {
		delete[] tensions;
		delete[] currents;
	};
	int getNObs() { return nobs; }
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *getTensions() { return tensions; }
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *getCurrents() { return currents; }
	const char* getMeshFilename() { return mesh_file.c_str(); }
	_Scalar getCurrentVal(int i) { return currentVals[i]; }
	int getCurrentsCount() { return (int)currentVals.size(); }

	void initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount) {};
};

template<>
inline void observations<double>::initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount) {
	std::ifstream file;
	std::ifstream filec;
  mesh_file = filename;
	filec.open(*filecurrents);
	file.open(mesh_file);

	int n = electrodesCount; // TODO: read number of measurements from file
	int valuesCount = std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());
	nobs = valuesCount / n;

	file.clear();
	file.seekg(0, std::ios::beg);
	tensions = new vectorx[nobs];
	currents = new vectorx[nobs];
	currentVals = vectorx(nobs);
	vectorx current(nodesCount);
	current.fill(0);
	int baseIndex = (int)current.size() - n;
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
		tensions[i].resize(n);
		Scalar val, avg;
		avg = 0;
		for (int j = 0; j<n; j++) {
			file >> val;
			tensions[i][j] = val / c;  // Values are normalized by current
			avg += tensions[i][j];
		}
		// rebase tensions, apply offset for zero sum on electrodes
		avg /= (double)n;
		for (int j = 0; j < n; j++) {
			tensions[i][j] -= avg;
		}
	}

	filec.close();
	file.close();
};

template<>
inline void observations<std::complex<double>>::initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount) {
	std::ifstream file;
	std::ifstream filecin, filecout;

	filecin.open(filecurrents[0]);
	filecout.open(filecurrents[1]);
	file.open(filename);

	int n = electrodesCount;
	int valuesCount = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	nobs = valuesCount / n;

	file.clear();
	file.seekg(0, std::ios::beg);
	tensions = new vectorxcomplex[nobs];
	currents = new vectorxcomplex[nobs];
	currentVals = vectorxcomplex(nobs);
	vectorxcomplex current(nodesCount);
	current.fill(0);
	int baseIndex = (int)current.size() - n;
	for (int i = 0; i<nobs; i++) {
		char ch;
		double valreal, valimag;
		currents[i] = current;
		filecin >> ch >> valreal >> ch >> valimag >> ch;
		currents[i][baseIndex + i] = std::complex<double>(valreal, valimag);
		filecout >> ch >> valreal >> ch >> valimag >> ch;
		currents[i][i == nobs - 1 ? baseIndex : baseIndex + i + 1] = std::complex<double>(valreal, valimag);

		// read tensions from file
		tensions[i].resize(n);
		std::complex<double> avg;
		avg = 0;
		for (int j = 0; j<n; j++) {
			file >> ch >> valreal >> ch >> valimag >> ch;
			tensions[i][j] = std::complex<double>(valreal, valimag);
		}
	}

	filecin.close();
	filecout.close();
	file.close();
};

#endif // OBSERVATIONS_H_
