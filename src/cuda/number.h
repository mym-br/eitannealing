#ifndef NUMBER_H
#define NUMBER_H

#include "settings.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

class Number {
protected:
	double * data;

public:
	Number(double n);
	~Number();

	double * getData() {
		return data;
	}
	void copy(Number * src);
	void copy(Number * src, cudaStream_t stream);

	double transf2CPU();
};

#endif /* NUMBER_H */
