#include <vector>
#include <memory>
#include <ctime>
#include <iostream>
#include "../src/util/fill_with_smallest.hpp"
#include "util/timestamp/timestamp.h"
#include <random>

// Just a simple sifdown algorithm for incomplete_lq_builder
// This probably works only with std::vector iterators
template<class iterator, class compare> void heap_sift_top_down2(iterator s, const iterator &start, const iterator &end, const compare &cmp)
{
    typename std::iterator_traits<iterator>::difference_type a = std::distance(start, s);
    iterator c1(s);
    std::advance(c1,a+1);
    iterator c2 = c1;
    c2++;
    if(c1 >= end) return;
    while(end > c2) {
        iterator t = cmp(*c2,*c1)?c1:c2;
        std::cout << "hsfd p: " << *s << " c1: " << *c1 << " c2: " << *c2 << " t: " << *t;
        if(cmp(*s,*t)) {
            std::iter_swap(s, t);
            s = t;
            std::cout << "!\n";
        } else {
            std::cout << "\n";
            return;
        }
        a = std::distance(start, s);
        c1 = s;
        std::advance(c1,a+1);
        c2 = c1;
        c2++;
    }
    if(c1 < end) {
        std::cout << "hsfd p: " << *s << " c: " << *c1;
        if(cmp(*s, *c1)) {
            std::iter_swap(s, c1);
            std::cout << "!\n";
        } else {
            std::cout << "\n";
        }
    }
}


template<class iterator, class compare> void make_heap_down2(const iterator &start, const iterator &end, const compare &cmp)
{
    std::cout<< "make_heap_down2\nVector:\n";
    for(iterator i = start; i!= end; i++)
        std::cout << *i << std::endl;

    iterator last_parent(start);
    typename std::iterator_traits<iterator>::difference_type length = std::distance(start, end);
    std::cout << "Length: " << length << std::endl;
    for(std::advance(last_parent, length/2 - 1); last_parent >= start; std::advance(last_parent, -1)) {
        std::cout << "Last parent value: " << *last_parent << std::endl;
        std::cout<< "Vector before heap_sift_top_down:\n";
        for(iterator i = start; i!= end; i++)
            std::cout << *i << std::endl;
        heap_sift_top_down2(last_parent, start, end, cmp);
        std::cout<< "Vector after heap_sift_top_down:\n";
        for(iterator i = start; i!= end; i++)
            std::cout << *i << std::endl;
    }
}


int main()
{
    std::mt19937 gen32;


    std::vector<int> origin;
    std::vector<int> dest;

    std::vector<int> temp;
    std::vector<int> prev;
    gen32.seed(get_usec_timestamp());

    for(int j = 0; j < 100; j++) {
        temp.clear();
        prev.clear();
        for(int i = 0; i < 5; i++) {
            int x = gen32()%20;
            temp.push_back(x);
            prev.push_back(x);
        }

        make_heap_down(temp.begin(), temp.end(), [](int a, int b) {
            return a>b;
        });

        if(!std::is_heap(temp.begin(), temp.end(), [](int a, int b) {
            return a>b;
        })) {
            std::cout << "make_heap_down error: Result is not a heap!\n";
            std::cout << "Vector before make_heap_down:\n";
            for(auto x : prev)
                std::cout << x << std::endl;
            std::cout << "Vector after make_heap_down:\n";
            for(auto x : temp)
                std::cout << x << std::endl;
            std::cout << "Running make_heap_down2\n";
            make_heap_down2(prev.begin(), prev.end(), [](int a, int b) {
                return a>b;
            });

            return -1;
        }

        for(int i=0;i<100;i++) {
            temp[0]=gen32();
            std::vector aux(temp);
            heap_sift_top_down(temp.begin(), temp.begin(), temp.end(), [](int a, int b) {
                return a>b;
            });
            if(!std::is_heap(temp.begin(), temp.end(), [](int a, int b) {
                return a>b;
            })) {
                std::cout << "Error!\n";
                std::cout << "Vector before sift_down:\n";
                for(auto x : aux) {
                    std::cout << i << std::endl;
                }
                std::cout << "Vector after sift_down:\n";
                for(auto x : temp) {
                    std::cout << i << std::endl;
                }
                return -1;
            }
        }
    }

    for(int i = 0; i < 100; i++) {
        origin.push_back(gen32());
    }

    fillWithNSmallest(dest, origin, 10, [](int a, int b) {
        return a>b;
    });
    std::cout << "dest vector:\n";
    for(auto x : dest)
        std::cout << x << std::endl;

    if(std::is_heap(dest.begin(), dest.end(), [](int a, int b) {
        return a>b;
    })) {
        std::cout << "ok!\n";
    } else {
        std::cout << "not a heap!\n";
        return -1;
    }


    return 0;
}
