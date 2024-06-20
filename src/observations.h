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
class observations
{
private:
	int nobs;
	std::string mesh_file;
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *tensions;
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *currents;
	Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> currentVals;

public:
	observations() : nobs(-1), tensions(nullptr), currents(nullptr) {}
	~observations()
	{
		delete[] tensions;
		delete[] currents;
	};
	int getNObs() const { return nobs; }
	const Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *getTensions() const { return tensions; }
	const Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> *getCurrents() const { return currents; }
	const char *getMeshFilename() const { return mesh_file.c_str(); }
	_Scalar getCurrentVal(int i) const { return currentVals[i]; }
	int getCurrentsCount() const { return (int)currentVals.size(); }

	void initObs(const char **filecurrents, const char *filename, int nodesCount, int electrodesCount, int groundNode = -1, int baseIndex = -1) {};
};

template <>
inline void observations<double>::initObs(const char **filecurrents, const char *filename, int nodesCount, int electrodesCount, int groundNode, int baseIndex)
{
	std::ifstream file;
	std::ifstream filec;
	bool readPotentials = (filename != NULL);
	filec.open(*filecurrents);
	if (readPotentials)
		file.open(filename);

	int n = electrodesCount - 1; // TODO: read number of measurements from file
	long valuesCount = std::distance(std::istream_iterator<double>(filec), std::istream_iterator<double>());
	nobs = valuesCount / 3;

	filec.clear();
	filec.seekg(0, std::ios::beg);
	tensions = new vectorx[nobs];
	currents = new vectorx[nobs];
	currentVals = vectorx(nobs);
	vectorx current(nodesCount - 1);
	current.fill(0);
	if (baseIndex < 0)
		baseIndex = (int)current.size() - n;
	for (int i = 0; i < nobs; i++)
	{
		double c;
		int entry, exit;
		filec >> entry;
		filec >> exit;
		filec >> c;
		entry--;
		exit--; // zero-based
		currentVals[i] = c;
		currents[i] = current;
		if (entry != n)
			currents[i][baseIndex + entry] = 1;
		if (exit != n)
			currents[i][baseIndex + exit] = -1;

		// read tensions from file
		if (!readPotentials)
			continue;
		tensions[i].resize(n);
		Scalar val;
		for (int j = 0; j < n; j++)
		{
			file >> val;
			tensions[i][j] = val / c; // Values are normalized by current
		}
		// rebase tensions, as our last electrode is always the ground
		file >> val;
		for (int j = 0; j < n; j++)
		{
			tensions[i][j] -= val / c;
		}
	}

	filec.close();
	file.close();
};

// FIXME: This actually ignores the 2nd current file, and currents can only be real
template <>
inline void observations<std::complex<double>>::initObs(const char **filecurrents, const char *filename, int nodesCount, int electrodesCount, int groundNode, int baseIndex)
{
	std::ifstream file;
	std::ifstream filec;
	// mesh_file = filename;
	filec.open(*filecurrents);
	file.open(filename);

	int n = electrodesCount; // TODO: read number of measurements from file
	long valuesCount = (long)std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());
	nobs = valuesCount / (2 * electrodesCount);

	file.clear();
	file.seekg(0, std::ios::beg);
	tensions = new vectorxcomplex[nobs];
	currents = new vectorxcomplex[nobs];
	vectorx current(nodesCount);
	current.fill(0);
	if (baseIndex < 0)
		baseIndex = (int)current.size() - n;
	for (int i = 0; i < nobs; i++)
	{
		double c;
		int entry, exit;
		filec >> entry;
		filec >> exit;
		filec >> c;
		entry--;
		exit--; // zero-based
		currents[i] = current;
		currents[i][baseIndex + entry] = 1;
		currents[i][baseIndex + exit] = -1;
		// read tensions from file
		tensions[i].resize(n);
		Scalar val_re, val_im;
		for (int j = 0; j < n; j++)
		{
			file >> val_re; // 2 coefficients per tension
			file >> val_im;
			tensions[i][j] = std::complex<double>(val_re, val_im) / c; // Values are normalized by current
		}
		// rebase tensions, as our last electrode is always the ground
		std::complex<double> val = tensions[i][n - 1];
		for (int j = 0; j < n - 1; j++)
		{
			tensions[i][j] -= val;
		}
		tensions[i][n - 1] = 0.0f;
	}

	filec.close();
	file.close();
};

#endif // OBSERVATIONS_H_
