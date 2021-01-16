#ifndef CIRCULARBUFF_H_
#define CIRCULARBUFF_H_

#include <cstddef>

template<class base, size_t len> class circularbuff
{
	base _val[len];
public:

	typedef base basetype;

	const base &operator[](int i) const {
		if (i >= 0)
			return _val[i%len];
		return _val[(len - ((-i) % len)) % len];
	}

	base &operator[](int i) {
		if (i >= 0)
			return _val[i%len];
		return _val[(len - ((-i) % len)) % len];
	}

	circularbuff() {};
};

#endif // CIRCULARBUFF_H_
