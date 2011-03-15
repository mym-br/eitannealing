/*
 * random.h
 *
 *  Created on: Sep 14, 2010
 *      Author: thiago
 */

#ifndef RANDOM_H_
#define RANDOM_H_

//#ifdef _cpluscplus
extern "C" {
//#endif

#include "mt19937-64/mt64.h"

//#ifdef _cpluscplus
}
//#endif

#define initrandom() (init_genrand64(0))

#define genreal() (genrand64_real1())
#define genint(n) (genrand64_int64()%(n))

#endif /* RANDOM_H_ */
