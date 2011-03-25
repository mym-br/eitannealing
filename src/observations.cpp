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

void initObs(char *filecurrents, char* filename)
{
	std::ifstream file;
	std::ifstream filec;

	filec.open(filecurrents);
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
		double c;
		filec >> c;
		currents[i] = current;
		if(i!=n)
			currents[i][baseIndex+i] = 1;
		if(i+4<n)
			currents[i][baseIndex+i+4] = -1;
		if(i+4>n)
			currents[i][baseIndex+i+3-n] = -1;
			

		// read tensions from file
		tensions[i].resize(electrodes.size()-1);
		double val;
		for(unsigned int j=0;j<electrodes.size()-1;j++) {			
			file >> val;
			tensions[i][j] = val/c;  // Values are normalized by current
		}
		// rebase tensions, as our last electrode is always the ground
		file >> val;
		for(unsigned int j=0;j<electrodes.size()-1;j++) {			
			tensions[i][j] -= val/c;
		}
	}
}


