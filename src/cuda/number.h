#ifndef NUMBER_H
#define NUMBER_H

#include "settings.h"

#include <cuda_runtime.h> 
#include <device_launch_parameters.h>

class Number {
protected:
	numType * data;

public:
	Number(numType n);
	~Number();

	numType * getData() {
		return data;
	}
	void copy(Number * src);
	void copy(Number * src, cudaStream_t stream);

	numType transf2CPU();
};

#endif /* NUMBER_H */