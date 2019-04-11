/*
* vector.h
*
*  Created on: Jun 4, 2015
*      Author: renato
*/

#ifndef VECTOR_H
#define VECTOR_H

#include <vector>

#include <device_launch_parameters.h>

#include "settings.h"
#include "number.h"

namespace cgl {

	class Vector {
	private:
		/* data array */
		numType * data;
		/* vector size */
		int size;

		/* cached number of blocks so that it does not need to be computed on demand */
		int blocks;
		/* pre-alloced memory for norm computation */
		numType * partial_d;

	public:
		Vector(int n);
		Vector(numType * v, int n);
		~Vector();

		/* returns vector size */
		int getSize() {
			return size;
		}

		/* returns vector's inner data array (memory location is not defined) */
		numType * getData() {
			return data;
		}

		/* returns vector's inner data array (memory location is not defined) */
		numType * getPartial() {
			return partial_d;
		}

		/* makes all elements 0 */
		void reset();
		/* makes all elements 0 */
		void reset(cudaStream_t stream);

		/* sets a single element */
		void set(int idx, numType val);
		/* sets a single element */
		void set(int idx, numType val, cudaStream_t stream);

		/* sets elements in batch */
		void set(int size, numType * src, int * indices);
		/* sets elements in batch */
		void set(int size, numType * src, int * indices, cudaStream_t stream);

		/* copies this vector to another */
		void copyTo(Vector * target);
		/* copies this vector to another */
		void copyTo(Vector * target, cudaStream_t stream);

		/* copy values from another array to this vector */
		void copy(int size, numType * src);
		/* copy values from another array to this vector */
		void copy(int size, numType * src, cudaStream_t stream);

		/* computes vector norm */
		numType norm();
		/* computes vector norm */
		numType norm(cudaStream_t stream);

		//void swap(int a, int b);

		/* Sum of two vectors: r = a + b */
		void sum(Vector * b, Vector * r);
		/* Sum of two vectors: r = a + b */
		void sum(Vector * b, Vector * r, cudaStream_t stream);

		/* Subtraction of two vectors: r = a - b */
		void subtr(Vector * b, Vector * r);
		/* Subtraction of two vectors: r = a - b */
		void subtr(Vector * b, Vector * r, cudaStream_t stream);

		/* Scaling a vector: r = k * a; k: scalar */
		void scalar(Number * k, Vector * r);
		/* Scaling a vector: r = (k/m) * a; k, m: scalar */
		void scalar(Number * k, Number * m, Vector * r);
		/* Scaling a vector: r = k * a; k: scalar */
		void scalar(Number * k, Vector * r, cudaStream_t stream);
		/* Scaling a vector: r = (k/m) * a; k, m: scalar */
		void scalar(Number * k, Number * m, Vector * r, cudaStream_t stream);
		/* Scaling a vector: r = (k/m) * a; k, m: scalar */
		void scalarAdd(Number * k, Number * m, Vector * r, cudaStream_t stream);
		/* Scaling a vector: r = (k/m) * a; k, m: scalar, m is totalized (and saved) from partials before-hand */
		void scalarAddTotalization(Number * k, Number * m, Vector * r, cudaStream_t stream);
		/* Scaling a vector: r = (k/m) * a; k, m: scalar */
		void scalarSubtr(Number * k, Number * m, Vector * r, cudaStream_t stream);
		/* Multiple vector operations:
		 * Scales and adds a vector: r += (k/m) * v;
		 *    where k, m: scalar, m is totalized (and saved) from partials before-hand (belongs to v)
		 * Scales and subtracts a vector: s -= (k/m) * u
		 /    where k, m: scalar, m was previously totalized */
		void scalarAddSubtrTotalization(Number * k, Number * m, Vector * r, Vector * s, Vector * u, cudaStream_t stream);

		/* Inner product two vectors: r = a . b (dot product) */
		void inner(Vector * b, Number * r);
		/* Inner product two vectors: r = a . b (dot product) */
		void inner(Vector * b, Number * r, cudaStream_t stream);

		/* Apply mask to vector: a = a .* mask (in place!) */
		void mask(Vector * mask);
		/* Apply mask to vector: a = a .* mask (in place!) */
		void mask(Vector * mask, cudaStream_t stream);
	};

	typedef std::vector<std::vector<numType>> VectorArray;
}

#endif /* VECTOR_H */