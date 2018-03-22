/*
 * incomplete_lq_builder.h
 *
 *  Created on: Jul 27, 2010
 *      Author: thiago
 */

#ifndef INCOMPLETE_LQ_BUILDER_H_
#define INCOMPLETE_LQ_BUILDER_H_

#include <utility>
#include <vector>
#include <queue>
#include <map>

template<class scalar> class SparseIncompleteLQBuilder
{
  protected:


    typedef std::pair<unsigned long, scalar> i_c; // Pair index coefficient

    // Those persistent objects should provide storage for
    // temporary results, so that new calls to buildLMatrix()
    // won't incur into memory allocation
    std::vector<std::vector<i_c> > rows;
    std::map<unsigned long, scalar> buildingL;
    static const auto compare_coefficients = [](const i_c &a, const i_c &b){
       return a.second > b.second;
    };
    std::priority_queue<i_c,
      decltype(compare_coefficients)> selected_l;




  public:

    SparseIncompleteLQBuilder(){};


};

#endif  // INCOMPLETE_LQ_BUILDER_H_
