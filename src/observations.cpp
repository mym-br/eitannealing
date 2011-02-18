/*
 * observations.cpp
 *
 *  Created on: Sep 10, 2010
 *      Author: thiago
 */

#include "observations.h"
#include "solver.h"
#include "graphics.h"

#include "problemdescription.h"

int nobs;

Eigen::VectorXd *tensions;
Eigen::VectorXd *currents;

viewport *exact;


float *solution;

// Just simulate the system to get an "observation"
void initObs()
{
	
	float *solution = new float[65];
	int i;
	solution[0] = 1.0;
	for(i=1;i<65;i++) {
	//	int x = (i-1)%8;
	//	int y = (i-1)/8;
	//	solution[i] = 1+ ((x/2)+(y/2))%2;
		solution[i] = 2.0;
	}
	
	solution[28] = 1.0;
	solution[29] = 1.0;
	solution[36] = 1.0;
	solution[37] = 1.0;


	// Prepare a problem matrix
	matrix *stiffness;
	assembleProblemMatrix(solution, &stiffness);
	
	SparseIncompleteLLT pre(*stiffness);

	double exitcurrent = -(double)1/(numElectrodes-1);

	nobs = numElectrodes - 1;
	tensions = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	Eigen::VectorXd current(numNodes-1);
	current.fill(0);
	current.start(nobs).fill(exitcurrent);
	for(i=0;i<nobs;i++) {
		current[i] = 1;
		// Prepare result
		currents[i] = Eigen::VectorXd::Zero(numNodes-1);
		currents[i].start(nobs).fill(exitcurrent);
		currents[i][i]=1;

		// solve the system
		CG_Solver solver(*stiffness, current, pre);
		for(int j=0;j<numNodes+20;j++)
			solver.do_iteration();
		tensions[i] = solver.getX().start(numElectrodes-1);
		// Reset the current
		current[i] = exitcurrent;
	}

	exact =new viewport(400, 400, "Solution");
	exact->setCurrentSolution(solution);
	exact->show();
}


