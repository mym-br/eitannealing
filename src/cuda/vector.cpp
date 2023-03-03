#include "vector.h"

#include "settings.h"
#ifndef USECUDA

#include <math.h>

#include "utils.h"

using namespace cgl;

Vector::Vector(int n) {
	data = new double[n];
	for (int i = 0; i < n; i++) {
		data[i] = 0;
	}
	size = n;
}
Vector::Vector(double * v, int n) {
	this->data = new double[n];
	for (int i = 0; i < n; i++) {
		this->data[i] = v[i];
	}
	size = n;
}
Vector::~Vector() {
	delete data;
}
double * Vector::getData() {
	return data;
}
void Vector::reset() {
	for (int i = 0; i < size; i++) {
		data[i] = 0;
	}
}
void Vector::add(int idx, double val) {
	if (idx < 0 || idx >= size) {
		LOG("Invalid vector index!");
		return;
	}
	data[idx] += val;
}
void Vector::set(int idx, double val) {
	if (idx < 0 || idx >= size) {
		LOG("Invalid vector index!");
		return;
	}
	data[idx] = val;
}
void Vector::set(int size, double * src, int * indices) {
	if (idx < 0 || idx >= size) {
		LOG("Invalid vector index!");
		return;
	}
	for (int i = 0; i < size; i++) {
		data[indices[i]] = src[i];
	}
}
void Vector::copyTo(Vector * target) {
	for (int i = 0; i < this->size; i++) {
		target->data[i] = this->data[i];
	}
	target->size = this->size;
	// TODO: should check for different sizes
}
void Vector::copy(int size, double * src) {
	for (int i = 0; i < this->size; i++) {
		this->data[i] = src[i];
	}
	this->size = size;
	// TODO: should check for different sizes
}
double Vector::norm() {
	double sum = 0;
	for (int i = 0; i < size; i++) {
		sum += data[i] * data[i];
	}
	return sqrt(sum);
}
void Vector::swap(int a, int b) {
	double aux = data[a];
	data[a] = data[b];
	data[b] = aux;
}

#endif
