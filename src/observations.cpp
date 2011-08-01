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
	float va[65];
	int i;
	solution[0] = 1.0;
	for(i=1;i<65;i++) {
	//	int x = (i-1)%8;
	//	int y = (i-1)/8;
	//	solution[i] = 1+ ((x/2)+(y/2))%2;
		solution[i] = 2.0;
		va[i] = 0;
	}
	
	solution[1] = 1.0;
	solution[2] = 1.0;
	solution[9] = 1.0;
	solution[10] = 1.0;
	
	solution[5] = 1.0;
	solution[6] = 1.0;
	solution[13] = 1.0;
	solution[14] = 1.0;
	
	solution[19] = 1.0;
	solution[20] = 1.0;
	solution[27] = 1.0;
	solution[28] = 1.0;
	
	solution[23] = 1.0;
	solution[24] = 1.0;
	solution[31] = 1.0;
	solution[32] = 1.0;
	
	solution[33] = 1.0;
	solution[34] = 1.0;
	solution[41] = 1.0;
	solution[42] = 1.0;
	
	solution[37] = 1.0;
	solution[38] = 1.0;
	solution[45] = 1.0;
	solution[46] = 1.0;
	
	solution[51] = 1.0;
	solution[52] = 1.0;
	solution[59] = 1.0;
	solution[60] = 1.0;
	
	solution[55] = 1.0;
	solution[56] = 1.0;
	solution[63] = 1.0;
	solution[64] = 1.0;

	
	
	


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
	exact->setCurrentSolution(solution,va);
	exact->show();
}


