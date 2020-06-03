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
#include <fstream>
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

	void initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount, int groundNode = -1, int baseIndex = -1, int minIndex = -1, int maxIndex = -1, bool posDirection = true) {};
};

template<>
inline void observations<double>::initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount, int groundNode, int baseIndex, int minIndex, int maxIndex, bool posDirection) {
	std::ifstream file;
	std::ifstream filec;
	bool readPotentials = (filename != NULL);
	filec.open(*filecurrents);
	if(readPotentials) file.open(filename);

#if defined(ZEROELECSUM) || defined(BLOCKGND)
	int n = electrodesCount; // TODO: read number of measurements from file
#else
	int n = electrodesCount - 1; // TODO: read number of measurements from file
#endif
	//int valuesCount = std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());
	int valuesCount = std::distance(std::istream_iterator<double>(filec), std::istream_iterator<double>());
	nobs = valuesCount / 3;

	filec.clear();
	filec.seekg(0, std::ios::beg);
	tensions = new vectorx[nobs];
	currents = new vectorx[nobs];
	currentVals = vectorx(nobs);
#if defined(ZEROELECSUM) || defined(BLOCKGND)
	vectorx current(nodesCount);
#else
	vectorx current(nodesCount-1);
#endif
	current.fill(0);
	if(baseIndex < 0) baseIndex = (int)current.size() - n;
	for (int i = 0; i<nobs; i++) {
		double c;
		int entry, exit;
		filec >> entry;
		filec >> exit;
		filec >> c;
		entry--; exit--;	// zero-based
		int entryIndex;
		int exitIndex;
		if (minIndex > 0 && maxIndex > 0) {
			entryIndex = baseIndex + (posDirection ? entry : -entry);
			exitIndex = baseIndex + (posDirection ? exit : -exit);
			if (entryIndex > maxIndex) entryIndex = minIndex + (entryIndex - maxIndex); if (entryIndex < minIndex) entryIndex = maxIndex - (minIndex - entryIndex) + 1;
			if (exitIndex > maxIndex) exitIndex = minIndex + (exitIndex - maxIndex); if (exitIndex < minIndex) exitIndex = maxIndex - (minIndex - exitIndex) + 1;
		}
		else {
			entryIndex = baseIndex + entry;
			exitIndex = baseIndex + exit;
		}
		if (groundNode < 0) groundNode = n;
		currentVals[i] = c;
		currents[i] = current;
#if defined(ZEROELECSUM) || defined(BLOCKGND)
		currents[i][entryIndex] = 1;
		currents[i][exitIndex] = -1;
	#ifndef BLOCKGND
		// Zero ground node current
		if (groundNode > 0 && groundNode < nodesCount) currents[i][groundNode] = 0;
	#endif
#else
		if (entryIndex != groundNode) currents[i][entryIndex] = 1;
		if (exitIndex != groundNode) currents[i][exitIndex] = -1;
#endif
		// read tensions from file
		if (!readPotentials) continue;
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
		file >> val;
		for (int j = 0; j < n; j++) {
			tensions[i][j] -= val / c;
		}
#endif
	}

	filec.close();
	if (readPotentials) file.close();
};

template<>
inline void observations<std::complex<double>>::initObs(const char **filecurrents, const char* filename, int nodesCount, int electrodesCount, int groundNode, int baseIndex, int minIndex, int maxIndex, bool posDirection) {
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
	if (baseIndex < 0) baseIndex = (int)current.size() - n;
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
