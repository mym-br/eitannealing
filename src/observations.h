/*
* observations.h
*
*  Created on: Jun 14, 2017
*      Author: aksato
*/
#ifndef OBSERVATIONS_H_
#define OBSERVATIONS_H_

#include <iostream>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
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
	int getNObs() const { return nobs; }
	const Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *getTensions() const { return tensions; }
	const Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *getCurrents() const { return currents; }
	const char* getMeshFilename() const { return mesh_file.c_str(); }
	_Scalar getCurrentVal(int i) const { return currentVals[i]; }
	int getCurrentsCount() const { return (int)currentVals.size(); }

	void initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount, int groundNode = -1) {};
};

template<>
inline void observations<double>::initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount, int groundNode) {
	std::ifstream file;
	std::ifstream filec;
        //mesh_file = filename;
	filec.open(*filecurrents);
	file.open(filename);

	int n = electrodesCount; // TODO: read number of measurements from file
	long valuesCount = (long)std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());
	nobs = valuesCount / electrodesCount;

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
	#ifdef ZEROELECSUM
		// Zero ground node current
		if (groundNode > 0 && groundNode < nodesCount) currents[i][groundNode] = 0;
	#endif
		// read tensions from file
		tensions[i].resize(n);
		Scalar val, avg;
		avg = 0;
		for (int j = 0; j<n; j++) {
			file >> val;
			tensions[i][j] = val / c;  // Values are normalized by current
			avg += tensions[i][j];
		}
#if defined(ZEROELECSUM) || defined(BLOCKGND)
		// rebase tensions, apply offset for zero sum on electrodes
		avg /= (double)n;
		for (int j = 0; j < n; j++) {
			tensions[i][j] -= avg;
		}
#else
		// rebase tensions, as our last electrode is always the ground
		val = tensions[i][n-1];
		for (int j = 0; j < n-1; j++) {
			tensions[i][j] -= val;
		}
		tensions[i][n-1] = 0.0f;
#endif
	}

	filec.close();
	file.close();
};

template<>
inline void observations<std::complex<double>>::initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount, int groundNode) {
	std::ifstream file;
	std::ifstream filecin, filecout;

	filecin.open(filecurrents[0]);
	filecout.open(filecurrents[1]);
	file.open(filename);

	int n = electrodesCount;
	long valuesCount = (long)std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
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
