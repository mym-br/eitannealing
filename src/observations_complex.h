/*
 * observations.h
 *
 *  Created on: Feb 20, 2019
 *      Author: thiago
 */

#ifndef OBSERVATIONS_COMPLEX_H_
#define OBSERVATIONS_COMPLEX_H_

#include <Eigen/Core>
#include "observations.h"

extern Eigen::VectorXd *tensions_I;
extern Eigen::VectorXd *currents_I;

void initObsComplex(char *filecurrents, char* filetensions_R, char* filetensions_I);

#endif /* OBSERVATIONS_COMPLEX_H_ */
