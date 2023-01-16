#include "number.h"

#ifdef USECUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/*****************************************************************************************/
/*****************************************************************************************/
/**********************************    functions    **************************************/
/*****************************************************************************************/
/*****************************************************************************************/

Number::Number(double val) {
	// CUDA device memory allocation
	cudaMalloc((void**)& this->data, (size_t) sizeof (double));
	// CUDA memory copy
	cudaMemcpy(this->data, &val, (size_t) sizeof (double), cudaMemcpyHostToDevice);
}
Number::~Number() {
	cudaFree(this->data);
}
void Number::copy(Number * src) {
	// CUDA memory copy
	cudaMemcpy(this->data, src->data, (size_t) sizeof (double), cudaMemcpyHostToDevice);
}
void Number::copy(Number * src, cudaStream_t stream) {
	// CUDA memory copy
	cudaMemcpyAsync(this->data, src->data, (size_t) sizeof (double), cudaMemcpyHostToDevice, stream);
}

double Number::transf2CPU() {
	double val;
	cudaMemcpy(&val, this->data, (size_t)sizeof(double), cudaMemcpyDeviceToHost);
	return val;
}

#endif