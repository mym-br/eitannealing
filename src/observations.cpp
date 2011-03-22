/*
 * observations.cpp
 *
 *  Created on: Sep 10, 2010
 *      Author: thiago
 */

#include "observations.h"
#include "solver.h"
#include "graphics.h"
#include <fstream>
#include <string>
#include "observations.h"
#include "problemdescription.h"

int nobs;

Eigen::VectorXd *tensions;
Eigen::VectorXd *currents;

viewport *exact;


float *solution;

void initObs(char *filename)
{
	std::ifstream file;

	file.open(filename);

	// FIXME: Read from file!
	nobs = 32;

	tensions = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	Eigen::VectorXd current(electrodes.size()-1);
	current.fill(0);
	for(int i=0;i<nobs;i++) {
		if(i!=current.size())
			current[i] = 0.001;
		if(i+4!=current.size())
			current[(i+4)%current.size()] = -0.001;

		// read tensions from file
		tensions[i].resize(electrodes.size());
		for(unsigned int j=0;j<electrodes.size();j++) {
			double val;
			file >> val;
			tensions[i][j] = val;
		}
	}
}


