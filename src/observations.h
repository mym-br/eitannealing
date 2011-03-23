/*
 * observations.h
 *
 *  Created on: Sep 10, 2010
 *      Author: thiago
 */

#ifndef OBSERVATIONS_H_
#define OBSERVATIONS_H_

#include <Eigen/Core>

extern int nobs;

extern Eigen::VectorXd *tensions;
extern Eigen::VectorXd *currents;

void initObs(char *filecurrents, char* filename);

#endif /* OBSERVATIONS_H_ */
