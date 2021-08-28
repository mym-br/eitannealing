/*
 * fill_with_smallest.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: thiago
 */


#ifndef _FILL_WITH_SMALLEST_HPP
#define _FILL_WITH_SMALLEST_HPP

#include <algorithm>
#include "heap_siftdown.hpp"

// Fill dest with n smallest elements from source.
//  cmp(a,b) must return a<b
// container1 should be a storage class heap-friendly (is there anything else than vectors?)
template <class container1, class container2, class comparator> 
void fillWithNSmallest(container1 &dest, const container2 &orig, unsigned long n, const comparator &cmp) {
    dest.clear();
    auto oo = orig.begin();
    while(dest.size() < n && oo != orig.end()) {
        dest.push_back(*oo);
        oo++;
    }
    if(oo != orig.end()) { // There are remaining elements, get the n largest
        make_heap_down(dest.begin(), dest.end(), cmp); // max heap
        while(oo != orig.end()) {
            if(cmp(*oo, dest.front())) {
                // Replace smallest element and fix the heap
                dest.front() = *oo;
                heap_sift_top_down(dest.begin(), dest.begin(), dest.end(), cmp);
            }
            oo++;
        }
    }
    
}

#endif  // FILL_WITH_SMALLEST
