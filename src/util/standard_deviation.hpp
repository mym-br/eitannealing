/*
 * standard_deviation.hpp
 *
 *  Created on: Dec 19, 2018
 *      Author: thiago
 */
#ifndef _STANDARD_DEVIATION_HPP
#define _STANDARD_DEVIATION_HPP

#include <numeric>
#include <iterator>

template<class InputIterator> typename std::iterator_traits<InputIterator>::value_type population_variance(InputIterator start, InputIterator end)
{
    struct moments {
        typename std::iterator_traits<InputIterator>::value_type s;
        typename std::iterator_traits<InputIterator>::value_type ss;
        moments operator+(typename std::iterator_traits<InputIterator>::value_type x) {
            return { this->s + x, this->ss + x*x };
        }
    };
    typename std::iterator_traits<InputIterator>::difference_type n = end - start;
    moments m = std::accumulate(start, end, moments{0, 0});
    return (n*m.ss - m.s*m.s)/(n*n);
}

#endif  // _STANDARD_DEVIATION_HPP

