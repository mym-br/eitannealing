/*
 * init_obs_problem.h
 *
 *  Created on: Oct 10, 2010
 *      Author: thiago
 */

#ifndef INIT_OBS_PROBLEM_H_
#define INIT_OBS_PROBLEM_H_

namespace obs {
	void initObsProblem();
	matrix *buildObsProblemMatrix(float *coeff);

	extern int numNodes;
	extern int numElements;
	extern int numElectrodes;
}

#endif /* INIT_OBS_PROBLEM_H_ */
