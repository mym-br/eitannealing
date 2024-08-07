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
template<class iterator, class compare> void heap_sift_top_down(iterator s, const iterator &start, const iterator &end, compare &&cmp)
{
    typename std::iterator_traits<iterator>::difference_type a = std::distance(start, s);
    iterator c1(s);
    std::advance(c1, a+1);
    iterator c2 = c1;
    c2++;
    if(c1 >= end) return;
    while(end > c2) {
        iterator t = c1 + cmp(*c1, *c2);
        if(cmp(*s,*t)) {
            std::iter_swap(s, t);
            s = t;
        } else return;
        c1 = s;
        a = std::distance(start, s);
        std::advance(c1, a+1);
        c2 = c1;
        c2++;
    }
    if((c1 < end) && cmp(*s, *c1)) std::iter_swap(s, c1);
}

template<class iterator, class compare> void make_heap_down(const iterator &start, const iterator &end, compare &&cmp)
{
    iterator last_parent(start);
    typename std::iterator_traits<iterator>::difference_type length = std::distance(start, end);
    std::advance(last_parent, length/2 - 1);
    while(last_parent >= start) {
        heap_sift_top_down(last_parent, start, end, std::forward<compare>(cmp));
        std::advance(last_parent, -1);
    }
}

#endif  // _HEAP_SIFTDOWN_HPP
