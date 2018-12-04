#include "vector.h"

#include "settings.h"
#ifdef USECUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>

#include "utils.h"

#include "vector.h"
//#include "matrix-op.h"

#ifdef DEBUG
#include <sstream>
#endif

/*****************************************************************************************/
/*****************************************************************************************/
/***********************************    kernel    ****************************************/
/*****************************************************************************************/
/*****************************************************************************************/

__global__ void cv_init(int dim, numType * v) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	if (row < dim) {
		v[row] = 0;
	}

	__syncthreads();
}

__global__ void cv_sum(int dim, numType * v, numType * u, numType * r) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	if (row < dim) {
		r[row] = v[row] + u[row];
	}

	__syncthreads();
}

__global__ void cv_subtraction(int dim, numType * v, numType * u, numType * r) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	if (row < dim) {
		r[row] = v[row] - u[row];
	}

	__syncthreads();
}

__global__ void cv_scalar(int dim, numType * v, numType * a, numType * r) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	const numType val = a[0];

	if (row < dim) {
		r[row] = v[row] * val;
	}

	__syncthreads();
}

__global__ void cv_scalar(int dim, numType * v, numType * a, numType * b, numType * r) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	numType val;
	val = a[0] / b[0];

	if (row < dim) {
		r[row] = v[row] * val;
	}

	__syncthreads();
}

__global__ void cv_scalar_add(int dim, numType * r, numType * a, numType * b, numType * v) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	numType val;
	val = a[0] / b[0];

	if (row < dim) {
		r[row] += v[row] * val;
	}

	__syncthreads();
}

__global__ void cv_scalar_add_totalize(int dim, numType * r, numType * a, numType * b, numType * v,
	int blocks, numType * partials) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	// totalization
	numType aVal;
	if (tidx == 0) {
		numType sum = 0;
		for (int i = 0; i < blocks; i++) {
			sum += partials[i];
		}
		aVal = sum;
		if (row == 0) {
			a[0] = sum;
		}
	}
	__syncthreads();

	const numType val = aVal / b[0];

	if (row < dim) {
		r[row] += v[row] * val;
	}

	__syncthreads();
}

__global__ void cv_scalar_subtr(int dim, numType * r, numType * a, numType * b, numType * v) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	numType val;

	val = a[0] / b[0];
	if (row < dim) {
		r[row] -= v[row] * val;
	}

	__syncthreads();
}

__global__ void cv_scalar_add_subtr_totalize(int dim, numType * v, numType * k, numType * m,
	numType * r, numType * s, numType * u, int blocks, numType * partials) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int tidx = threadIdx.x;

	numType mVal;
	// partial totalization
	if (tidx == 0) {

		numType sum = 0;
		for (int i = 0; i < blocks; i++) {
			sum += partials[i];
		}
		mVal = sum;
		if (row == 0) {
			m[0] = sum;
		}
	}

	__syncthreads();

	// r += (k/m) * v
	const numType val = k[0] / m[0];
	if (row < dim) {
		// r += (k/m) * v
		r[row] += v[row] * val;
		// s -= (k/m) * u
		s[row] -= u[row] * val;
	}

	__syncthreads();
}

__global__ void cv_inner(int dim, numType * v, numType * u, numType * r, int blocks, numType * val) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	if (row < dim) {
		cache[cacheIndex] = v[row] * u[row];
	}
	else {
		cache[cacheIndex] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (cacheIndex < half) {
			cache[cacheIndex] += cache[cacheIndex + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (cacheIndex == 0) {
		r[blockIdx.x] = cache[0];
	}

	__syncthreads();

	if (row == 0) {
		numType sum = 0;
		for (int i = 0; i < blocks; i++) {
			sum += r[i];
		}
		val[0] = sum;
	}
}

// probably not used
__global__ void cv_inner(int dim, numType * v, numType * u, numType * r) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	if (row < dim) {
		cache[cacheIndex] = v[row] * u[row];
	}
	else {
		cache[cacheIndex] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (cacheIndex < half) {
			cache[cacheIndex] += cache[cacheIndex + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (cacheIndex == 0) {
		r[blockIdx.x] = cache[0];
	}

	__syncthreads();
}

__global__ void cv_mask(int dim, numType * v, numType * m) {
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;

	if (row < dim) {
		v[row] *= m[row];
	}

	__syncthreads();
}

// probably not used
__global__ void cv_reduction_sum(int dim, numType * v, numType * r) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// row index
	int row = threadIdx.x;

	if (row < dim) {
		cache[row] = v[row];
	}
	else {
		cache[row] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (row < half) {
			cache[row] += cache[row + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (row == 0) {
		r[0] = cache[0];
	}

	__syncthreads();
}

/*****************************************************************************************/
/*****************************************************************************************/
/***********************************    kernel    ****************************************/
/*****************************************************************************************/
/*****************************************************************************************/

// can be replaced by inner product with itself
__global__ void cv_norm(int dim, numType * v, numType * norm) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// row index
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int cacheIndex = threadIdx.x;

	if (row < dim) {
		cache[cacheIndex] = v[row] * v[row];
	}
	else {
		cache[cacheIndex] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (cacheIndex < half) {
			cache[cacheIndex] += cache[cacheIndex + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (cacheIndex == 0) {
		norm[blockIdx.x] = cache[0];
	}

	__syncthreads();
}

// can be replaced by inner product with itself
__global__ void cv_reduction(int dim, numType * v) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// row index
	int row = threadIdx.x;

	if (row < dim) {
		cache[row] = v[row];
	}
	else {
		cache[row] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (row < half) {
			cache[row] += cache[row + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (row == 0) {
		v[0] = cache[0];
	}

	__syncthreads();
}

__global__ void cv_reduction2(int dim, numType * v, numType * n) {
	// shared memory for reduction
	__shared__ numType cache[BLOCKSIZE];
	// row index
	int row = threadIdx.x;

	if (row < dim) {
		cache[row] = v[row];
	}
	else {
		cache[row] = 0.0;
	}

	__syncthreads();

	int half = BLOCKSIZE >> 1;
	while (half != 0){
		if (row < half) {
			cache[row] += cache[row + half];
		}

		__syncthreads();

		half >>= 1;
	}

	if (row == 0) {
		n[0] = cache[0];
	}

	__syncthreads();
}

__global__ void cv_add(int dim, numType * v, int idx, numType val) {
	if (idx < dim) {
		v[idx] += val;
	}
}

__global__ void cv_set(int dim, numType * v, int idx, numType val) {
	if (idx < dim) {
		v[idx] = val;
	}
}

__global__ void cv_set2(int dim, numType * v, numType * target, int * indices) {
	// row index
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < dim) {
		target[indices[idx]] = v[idx];
	}
}

/*****************************************************************************************/
/*****************************************************************************************/
/**********************************    functions    **************************************/
/*****************************************************************************************/
/*****************************************************************************************/

using namespace cgl;
Vector::Vector(int n) {
	size = n;
	blocks = ceil((double)size / BLOCKSIZE);

	// CUDA device memory allocation
	cudaMalloc((void**)& this->data, size * sizeof (numType));

	// CUDA device memory allocation
	cudaMalloc((void**)& this->partial_d, blocks * sizeof(numType));

	// initialization
	cv_init<<<blocks, BLOCKSIZE >>>(n, this->data);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_init kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

Vector::Vector(numType * v, int n) {
	size = n;
	blocks = ceil((double)size / BLOCKSIZE);

	// CUDA device memory allocation
	cudaMalloc((void**)& partial_d, blocks * sizeof(numType));

	// CUDA device memory allocation
	cudaMalloc((void**)& this->data, size * sizeof (numType));
	// CUDA memory copy
	cudaMemcpy(this->data, v, (size_t)size * sizeof (numType), cudaMemcpyHostToDevice);
}

Vector::~Vector() {
	cudaFree(data);
	cudaFree(partial_d);
}

/* makes all elements 0 */
void Vector::reset() {
	cv_init <<<blocks, BLOCKSIZE >>>(size, this->data);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_init kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* makes all elements 0 */
void Vector::reset(cudaStream_t stream) {
	cv_init <<<blocks, BLOCKSIZE, 0, stream>>>(size, this->data);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_init kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
}
#endif
}

/* sets a single element */
void Vector::set(int idx, numType val) {
#ifdef DEBUG
	if (idx < 0 || idx >= size) {
		LOG("set: Invalid vector index!");
		printf("set: Invalid vector index! %d", idx);
		return;
	}
#endif

	//cv_set <<<1, 1 >>>(size, data, idx, val);
	cudaMemcpy(&this->data[idx], &val, (size_t) sizeof(numType), cudaMemcpyHostToDevice);

#ifdef DEBUG
	// TODO: should check for kernel errors
#endif
}

/* sets a single element */
void Vector::set(int idx, numType val, cudaStream_t stream) {
#ifdef DEBUG
	if (idx < 0 || idx >= size) {
		LOG("set: Invalid vector index!");
		printf("set: Invalid vector index! %d", idx);
		return;
	}
#endif

	cudaMemcpyAsync(&(this->data[idx]), &(val), (size_t)size * sizeof(numType), cudaMemcpyHostToDevice, stream);

#ifdef DEBUG
	// TODO: should check for kernel errors
#endif
}

/* sets elements in batch */
void Vector::set(int size, numType * src, int * indices) {
	int sblocks = 1;
//#ifdef DEBUG
//	if (size < 0) {
//		LOG("set2: Invalid vector index!");
//		printf("set2: Invalid vector index! %d", idx);
//		return;
//	} else if (size > BLOCKSIZE) {
//		LOG("set2: Invalid vector index!");
//		printf("set2: Invalid vector index! %d", idx);
//		// recompute size?
//		sblocks = bceil((double)size / BLOCKSIZE);
//		size = BLOCKSIZE;
//	}
//#endif

	cv_set2 <<<sblocks, size >>>(this->size, data, src, indices);

#ifdef DEBUG
	// TODO: should check for kernel errors
#endif
}

/* sets elements in batch */
void Vector::set(int size, numType * src, int * indices, cudaStream_t stream) {
	int sblocks = 1;
//#ifdef DEBUG
//	if (size < 0) {
//		LOG("set2: Invalid vector index!");
//		printf("set2: Invalid vector index! %d", idx);
//		return;
//	} else if (size > BLOCKSIZE) {
//		LOG("set2: Invalid vector index!");
//		printf("set2: Invalid vector index! %d", idx);
//		// recompute size?
//		sblocks = bceil((double)size / BLOCKSIZE);
//		size = BLOCKSIZE;
//	}
//#endif

	cv_set2 <<<sblocks, size, 0, stream>>>(this->size, data, src, indices);

#ifdef DEBUG
	// TODO: should check for kernel errors
#endif
}

/* copies this vector to another */
void Vector::copyTo(Vector * target) {
#ifdef DEBUG
	// instead of erroring, it just warns
	if (size != this->size) {
		LOG("copy: Invalid source array size!");
		printf("copy: Invalid source array size! %d is different from %d", size, this->size);

		cudaFree(target->data);
		// TODO: should check for different sizes
		target->size = this->size;
		// CUDA device memory allocation
		cudaMalloc((void**)& target->data, size * sizeof(numType));
	}
#endif
	// CUDA memory copy
	cudaMemcpy(target->data, this->data, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice);
}

/* copies this vector to another */
void Vector::copyTo(Vector * target, cudaStream_t stream) {
#ifdef DEBUG
	// instead of erroring, it just warns
	if (size != this->size) {
		LOG("copy: Invalid source array size!");
		printf("copy: Invalid source array size! %d is different from %d", size, this->size);

		cudaFree(target->data);
		// TODO: should check for different sizes
		target->size = this->size;
		// CUDA device memory allocation
		cudaMalloc((void**)& target->data, size * sizeof(numType));
	}
#endif
	// CUDA memory copy
	cudaMemcpyAsync(target->data, this->data, (size_t)size * sizeof(numType), cudaMemcpyDeviceToDevice, stream);
}

/* copy values from another array to this vector */
void Vector::copy(int size, numType * src) {
	// TODO: should check for different sizes
#ifdef DEBUG
	if (size != this->size) {
		LOG("copy2: Invalid source array size!");
		printf("copy2: Invalid source array size! %d is different from %d", size, this->size);
	}
	return;
#endif
	// CUDA memory copy
	cudaMemcpy(this->data, src, (size_t)size * sizeof(numType), cudaMemcpyHostToDevice);
}

/* copy values from another array to this vector */
void Vector::copy(int size, numType * src, cudaStream_t stream) {
	// TODO: should check for different sizes
#ifdef DEBUG
	if (size != this->size) {
		LOG("copy2: Invalid source array size!");
		printf("copy2: Invalid source array size! %d is different from %d", size, this->size);
	}
	return;
#endif
	// CUDA memory copy
	cudaMemcpyAsync(this->data, src, (size_t)size * sizeof(numType), cudaMemcpyHostToDevice, stream);
}

/* computes vector norm */
numType Vector::norm() {
	// TODO: should we be transfering and allocing all the time?
	numType inner = 0;

	// Launch CUDA vector scalar multiplication kernel
	cv_norm << <blocks, BLOCKSIZE >> >(size, data, partial_d);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_norm kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif

	//LOGV(partial_d, blocks, "partial sum");
	cv_reduction <<<1, BLOCKSIZE >>>(blocks, partial_d);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_reduction kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif

	//LOGV(partial_d, blocks, "partial sum");
	cudaMemcpy(&inner, partial_d, (size_t) sizeof(numType), cudaMemcpyDeviceToHost);

#ifdef DEBUG
	// Check for any errors
#endif

	return sqrt(inner);
}

/* computes vector norm */
numType Vector::norm(cudaStream_t stream) {
	// TODO: should we be transfering and allocing all the time?
	numType inner = 0;

	// Launch CUDA vector scalar multiplication kernel
	cv_norm <<<blocks, BLOCKSIZE, 0, stream >>>(size, data, partial_d);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_norm kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif

	//LOGV(partial_d, blocks, "partial sum");
	cv_reduction <<<1, BLOCKSIZE, 0, stream>>>(blocks, partial_d);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_reduction kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif

	//LOGV(partial_d, blocks, "partial sum");
	cudaMemcpyAsync(&inner, partial_d, (size_t) sizeof(numType), cudaMemcpyDeviceToHost, stream);

#ifdef DEBUG
	// Check for any errors
#endif

	return sqrt(inner);
}

/* Sum of two vectors: r = a + b */
void Vector::sum(Vector * b, Vector * r) {
#ifdef DEBUG
	if (b == NULL || r == NULL) {
		LOG("v_add error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != b->getSize() || a->getSize() != r->getSize()) {
	//	LOG("v_add error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector sum kernel
	cv_sum << <blocks, BLOCKSIZE >> >(size, data, b->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_sum kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Sum of two vectors: r = a + b */
void Vector::sum(Vector * b, Vector * r, cudaStream_t stream) {
#ifdef DEBUG
	if (b == NULL || r == NULL) {
		LOG("v_add error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != b->getSize() || a->getSize() != r->getSize()) {
	//	LOG("v_add error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector sum kernel
	cv_sum <<<blocks, BLOCKSIZE, 0, stream >>>(size, data, b->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_sum kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Subtraction of two vectors: r = a - b */
void Vector::subtr(Vector * b, Vector * r) {
#ifdef DEBUG
	if (b == NULL || r == NULL) {
		LOG("v_subtr error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != b->getSize() || a->getSize() != r->getSize()) {
	//	LOG("v_subtr error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector subtraction kernel
	cv_subtraction << <blocks, BLOCKSIZE >> >(size, data, b->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_subtraction kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Subtraction of two vectors: r = a - b */
void Vector::subtr(Vector * b, Vector * r, cudaStream_t stream) {
#ifdef DEBUG
	if (b == NULL || r == NULL) {
		LOG("v_subtr error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != b->getSize() || a->getSize() != r->getSize()) {
	//	LOG("v_subtr error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector subtraction kernel
	cv_subtraction <<<blocks, BLOCKSIZE, 0, stream >>>(size, data, b->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_subtraction kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling a vector: r = k * a; k: scalar */
void Vector::scalar(Number * k, Vector * r) {
#ifdef DEBUG
	if (r == NULL) {
		LOG("v_scalar error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != r->getSize()) {
	//	LOG("v_scalar error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector scalar multiplication kernel
	cv_scalar << <blocks, BLOCKSIZE >> >(size, data, k->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling a vector: r = k * a; k: scalar */
void Vector::scalar(Number * k, Vector * r, cudaStream_t stream) {
#ifdef DEBUG
	if (r == NULL) {
		LOG("cv_scalar error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != r->getSize()) {
	//	LOG("cv_scalar error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector scalar multiplication kernel
	cv_scalar <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, k->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling a vector: r = (k/m) * a; k, m: scalar */
void Vector::scalar(Number * k, Number * m, Vector * r) {
#ifdef DEBUG
	if (r == NULL) {
		LOG("v_scalar error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != r->getSize()) {
	//	LOG("v_scalar error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector scalar multiplication kernel
	cv_scalar << <blocks, BLOCKSIZE >> >(size, data, k->getData(), m->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling a vector: r = (k/m) * a; k, m: scalar */
void Vector::scalar(Number * k, Number * m, Vector * r, cudaStream_t stream) {
#ifdef DEBUG
	if (r == NULL) {
		LOG("cv_scalar error: vector(s) null!\n");
		return;
	}
	if (size != r->getSize()) {
		LOG("cv_scalar error: vector(s) of different size!\n");
		return;
	}
#endif

	// Launch CUDA vector scalar multiplication kernel
	cv_scalar <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, k->getData(), m->getData(), r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling and adding a vector: r += (k/m) * a; k, m: scalar */
void Vector::scalarAdd(Number * k, Number * m, Vector * a, cudaStream_t stream) {
#ifdef DEBUG
	if (a == NULL) {
		LOG("cv_scalar_add error: vector(s) null!\n");
		return;
	}
	if (a->getSize() != size) {
		LOG("cv_scalar_add error: vector(s) of different size!\n");
		return;
	}
#endif

	// Launch CUDA vector scalar multiplication and sum kernel
	cv_scalar_add <<<blocks, BLOCKSIZE>>>(size, data, k->getData(), m->getData(), a->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar_add kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling a vector: r = (k/m) * a; k, m: scalar, m is totalized (and saved) from partials before-hand */
void Vector::scalarAddTotalization(Number * k, Number * m, Vector * r, cudaStream_t stream) {
#ifdef DEBUG
	//if (a == NULL) {
	//	LOG("cv_scalar_add_totalize error: vector(s) null!\n");
	//	return;
	//}
	//if (a->getSize() != size) {
	//	LOG("cv_scalar_add_totalize error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// Launch CUDA vector scalar multiplication and sum kernel
	cv_scalar_add_totalize <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, k->getData(), m->getData(), 
		r->getData(), blocks, partial_d);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar_add_totalize kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Scaling and subtracting a vector: r -= (k/m) * a; k, m: scalar */
void Vector::scalarSubtr(Number * k, Number * m, Vector * a, cudaStream_t stream) {
#ifdef DEBUG
	if (a == NULL) {
		LOG("cv_scalar_add error: vector(s) null!\n");
		return;
	}
	if (a->getSize() != size) {
		LOG("cv_scalar_add error: vector(s) of different size!\n");
		return;
	}
#endif

	// Launch CUDA vector scalar multiplication and subtraction kernel
	cv_scalar_subtr <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, k->getData(), m->getData(), a->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar_add kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Multiple vector operations:
* Scales and adds a vector: r += (k/m) * v;
*    where k, m: scalar, m is totalized (and saved) from partials before-hand (belongs to v)
* Scales and subtracts a vector: s -= (k/m) * u
/    where k, m: scalar, m was previously totalized */
void Vector::scalarAddSubtrTotalization(Number * k, Number * m, Vector * r, Vector * s, Vector * u, cudaStream_t stream) {
#ifdef DEBUG
	//if (a == NULL) {
	//	LOG("cv_scalar_add_totalize error: vector(s) null!\n");
	//	return;
	//}
	//if (a->getSize() != size) {
	//	LOG("cv_scalar_add_totalize error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	// partial totalization
	// r += (k/m) * v
	// s -= (k/m) * u
	cv_scalar_add_subtr_totalize <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, k->getData(), m->getData(), 
		r->getData(), s->getData(), u->getData(), blocks, partial_d);

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "cv_scalar_add_totalize kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Inner product two vectors: r = a . b (dot product) */
void Vector::inner(Vector * b, Number * r) {
#ifdef DEBUG
	if (b == NULL || r == NULL) {
		LOG("v_inner error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != b->getSize()) {
	//	LOG("v_inner error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	cv_inner << <blocks, BLOCKSIZE >> >(size, data, b->getData(), partial_d, blocks, r->getData());
	cv_reduction2 <<<1, BLOCKSIZE >>>(blocks, partial_d, r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "v_inner kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Inner product two vectors: r = a . b (dot product) */
void Vector::inner(Vector * b, Number * r, cudaStream_t stream) {
#ifdef DEBUG
	if (b == NULL || r == NULL) {
		LOG("v_inner error: vector(s) null!\n");
		return;
	}
	//if (a->getSize() != b->getSize()) {
	//	LOG("v_inner error: vector(s) of different size!\n");
	//	return;
	//}
#endif

	cv_inner <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, b->getData(), partial_d, blocks, r->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "v_inner kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Apply mask to vector: a = a .* mask (in place!) */
void Vector::mask(Vector * mask) {
//#ifdef DEBUG
//	if (b == NULL) {
//		LOG("v_mask error: vector(s) null!\n");
//		return;
//	}
//	if (a->getSize() != b->getSize()) {
//		LOG("v_inner error: vector(s) of different size!\n");
//		return;
//	}
//#endif

	cv_mask << <blocks, BLOCKSIZE >> >(size, data, mask->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "v_mask kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}

/* Apply mask to vector: a = a .* mask (in place!) */
void Vector::mask(Vector * mask, cudaStream_t stream) {
//#ifdef DEBUG
//	if (b == NULL) {
//		LOG("v_mask error: vector(s) null!\n");
//		return;
//	}
//	if (a->getSize() != b->getSize()) {
//		LOG("v_mask error: vector(s) of different size!\n");
//		return;
//	}
//#endif

	cv_mask <<<blocks, BLOCKSIZE, 0, stream>>>(size, data, mask->getData());

#ifdef DEBUG
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::ostringstream msg;
		msg << "v_mask kernel failed: " << cudaGetErrorString(cudaStatus) << "\n";
		msg.flush(); LOG(&msg.str()[0]);
	}
#endif
}


#endif