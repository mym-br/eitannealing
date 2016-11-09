#ifndef _INTCOEF_H_
#define _INTCOEF_H_

#include <memory>

#include "basematrix.h"
#include "problemdescription.h"
#include "nodecoefficients.h"

class intCoef {
  protected:
    vectorx x;
    int numcoefficients;
  public:
    intCoef(
      const std::vector<node> &nodes, 	       
      const std::vector<triangularElement> &elements,
      const std::map<int, int> &node2coefficient,
      int numcoefficients,
      float totalheight); 
    
    double getInt(double *coefficients) const;
};

#endif	// _INTCOEF_H_