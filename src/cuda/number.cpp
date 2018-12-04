#include "number.h"

#ifndef USECUDA

Number::Number(numType val) {
	data = new numType[1];
	data[0] = val;
}
Number::~Number() {
	delete data;
}
void Number::copy(Number * src) {
	this->data[0] = src->data[0];
}

#endif
