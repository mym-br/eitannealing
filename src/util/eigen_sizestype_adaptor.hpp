/*
 * braket_functor_adaptor.hpp
 *
 *  Created on: Aug 05, 2018
 *      Author: thiago
 */


#ifndef _EIGEN_SIZESTYPE_ADAPTOR_HPP
#define _EIGEN_SIZESTYPE_ADAPTOR_HPP


template<class fn> class sizestype_adaptor {
protected:
    fn &_f;
public:
    sizestype_adaptor(fn &f): _f(f) {}
    typedef decltype(_f(0)) value_type;
    value_type operator[](unsigned long i) const { return _f(i); }
};

template<class fn> sizestype_adaptor<fn> make_sizestype_adaptor(fn &&f) {
    return sizestype_adaptor<fn>(f);
}

#endif  // _EIGEN_SIZESTYPE_ADAPTOR_HPP
