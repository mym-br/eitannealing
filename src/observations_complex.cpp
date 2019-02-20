/*
* observations.cpp
*
*  Created on: Sep 10, 2010
*      Author: thiago
*/

#include "observations_complex.h"
#include "solver.h"
//#include "graphics.h"
#include <fstream>
#include <string>
#include "observations.h"
#include "problemdescription.h"

Eigen::VectorXd *tensions_I;
Eigen::VectorXd *currents_I;

//viewport *exact;

void initObsComplex(char *filecurrents, char* filetensions)
{
	std::ifstream filec;
	std::ifstream filet;

	filec.open(filecurrents);
	filet.open(filetensions);

	// FIXME: Read from file!
	nobs = 32;
	int n = gelectrodes.size()-1;

	tensions = new Eigen::VectorXd[nobs];
	tensions_I = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	currents_I = new Eigen::VectorXd[nobs];
	Eigen::VectorXd current(nodes.size()-1);
	current.fill(0);
	Eigen::VectorXd current_I(current);
	int baseIndex = current.size()-n;
	for(int i=0;i<nobs;i++) {
		double c;
		int entry, exit;
		filec >> entry;
		filec >> exit;
		filec >> c;
		entry--; exit--;	// zero-based
		currents[i] = current;
		currents_I[i] = current_I;
		if(entry!=n)
			currents[i][baseIndex+entry] = 1;
		if(exit!=n)
			currents[i][baseIndex+exit] = -1;

		// read tensions from file
		tensions[i].resize(gelectrodes.size()-1);
		double val_R, val_I;
		for(unsigned int j=0;j<gelectrodes.size()-1;j++) {
			filet >> val_R;
			filet >> val_I;
			tensions[i][j] = val_R/c;  // Values are normalized by current
			tensions_I[i][j] = val_I/c;  // Values are normalized by current
		}
		// rebase tensions, as our last electrode is always the ground
		filet >> val_R;
		filet >> val_I;
		for(unsigned int j=0;j<gelectrodes.size()-1;j++) {
			tensions[i][j] -= val_R/c;
			tensions_I[i][j] -= val_I/c;
		}
	}

	filec.close();
	filet.close();
}
