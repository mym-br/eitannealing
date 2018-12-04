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

Number::Number(numType val) {
	// CUDA device memory allocation
	cudaMalloc((void**)& this->data, sizeof (numType));
	// CUDA memory copy
	cudaMemcpy(this->data, &val, (size_t) sizeof (numType), cudaMemcpyHostToDevice);
}
Number::~Number() {
	cudaFree(data);
}
void Number::copy(Number * src) {
	// CUDA memory copy
	cudaMemcpy(this->data, src->data, (size_t) sizeof (numType), cudaMemcpyHostToDevice);
}
void Number::copy(Number * src, cudaStream_t stream) {
	// CUDA memory copy
	cudaMemcpyAsync(this->data, src->data, (size_t) sizeof (numType), cudaMemcpyHostToDevice, stream);
}

numType Number::transf2CPU() {
	numType *val = new numType[1];
	cudaMemcpy(val, data, (size_t)sizeof(numType), cudaMemcpyDeviceToHost);
	return val[0];
}

#endif