/*
 * observations.cpp
 *
 *  Created on: Sep 10, 2010
 *      Author: thiago
 */

#include "observations.h"
#include "solver.h"
#include "graphics.h"

#include "observations.h"
#include "init_obs_problem.h"
#include "problemdescription.h"

int nobs;

Eigen::VectorXd *tensions;
Eigen::VectorXd *currents;

viewport *exact;


float *solution;

// Just simulate the system to get an "observation"
void initObs()
{
	//obs::initObsProblem();

	float *solution = new float[65];
	int i;
	solution[0] = 1.0;
	for(i=1;i<65;i++) {
	//	int x = (i-1)%8;
	//	int y = (i-1)/8;
	//	solution[i] = 1+ ((x/2)+(y/2))%2;
		solution[i] = 1.0;
	}
	solution[19] = 2.0;
	solution[20] = 2.0;
	solution[21] = 2.0;
	solution[22] = 2.0;

	solution[27] = 2.0;
	solution[28] = 2.0;
	solution[29] = 2.0;
	solution[30] = 2.0;
	solution[35] = 2.0;
	solution[36] = 2.0;
	solution[37] = 2.0;
	solution[38] = 2.0;
	solution[43] = 2.0;
	solution[44] = 2.0;
	solution[45] = 2.0;
	solution[46] = 2.0;



	// Prepare a problem matrix
	matrix *stiffness = obs::buildObsProblemMatrix(solution);
	SparseIncompleteLLT pre(*stiffness);

	double exitcurrent = -(double)1/(obs::numElectrodes-1);

	nobs = obs::numElectrodes - 1;
	tensions = new Eigen::VectorXd[nobs];
	currents = new Eigen::VectorXd[nobs];
	Eigen::VectorXd current(obs::numNodes-1);
	current.fill(0);
	current.start(nobs).fill(exitcurrent);
	for(i=0;i<nobs;i++) {
		current[i] = 1;
		// Prepare result
		currents[i] = Eigen::VectorXd::Zero(nodes.size()-1);
		currents[i].start(nobs).fill(exitcurrent);
		currents[i][i]=1;

		// solve the system
		CG_Solver solver(*stiffness, current, pre);
		for(int j=0;j<obs::numNodes+20;j++)
			solver.do_iteration();
		tensions[i] = solver.getX().start(obs::numElectrodes-1);
		// Reset the current
		current[i] = exitcurrent;
	}

	exact =new viewport(400, 400, "Solution");
	exact->setCurrentSolution(solution);
	exact->show();
}


