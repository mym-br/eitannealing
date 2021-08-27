/*
 * heap_sifdown.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: thiago
 */
#ifndef _HEAP_SIFTDOWN_HPP
#define _HEAP_SIFTDOWN_HPP

#include <iterator>

// Just a simple sifdown algorithm for incomplete_lq_builder
// This probably works only with std::vector iterators
template<class iterator, class compare> void heap_sift_top_down(iterator s, const iterator &end, const compare &cmp)
{
    iterator c1(s);
    typename std::iterator_traits<iterator>::difference_type a = std::distance(s, ++c1); 
    iterator c2 = c1;
    c2++;
    if(c1 >= end) return;
    while(end - c2 > 0) {
        if(cmp(*c2,*c1) && cmp(*s, *c1)) {
            std::iter_swap(s, c1);
            s = c1;
        } else if(cmp(*s, *c2)) {
            std::iter_swap(s, c2);
            s = c2;
        } else return;
        a += a;
        c1 = s; std::advance(c1,a);
        c2 = c1; c2++;
    }
    if((c1 < end) && cmp(*s, *c1)) std::iter_swap(s, c1);
}

template<class iterator, class compare> void make_heap_down(const iterator &start, const iterator &end, const compare &cmp)
{
    iterator last_parent(start);
    typename std::iterator_traits<iterator>::difference_type length = std::distance(start, end);
    std::advance(last_parent, length/2 - 1);
    while(last_parent >= start) {
        heap_sift_top_down(last_parent, end, cmp);
        std::advance(last_parent, -1);
    }
}

#endif  // _HEAP_SIFTDOWN_HPP
