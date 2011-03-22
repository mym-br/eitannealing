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
	int n = electrodes.size()-1;
	

	tensions = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	Eigen::VectorXd current(nodes.size()-1);
	current.fill(0);
	int baseIndex = current.size()-n;
	for(int i=0;i<nobs;i++) {
		currents[i] = current;
		if(i!=n)
			currents[i][baseIndex+i] = 0.001;
		if(i+4<n)
			currents[i][baseIndex+i+4] = -0.001;
		if(i+4>n)
			currents[i][baseIndex+i+3-n] = -0.001;
			

		// read tensions from file
		tensions[i].resize(electrodes.size());
		for(unsigned int j=0;j<electrodes.size();j++) {
			double val;
			file >> val;
			tensions[i][j] = val;
		}
	}
}


