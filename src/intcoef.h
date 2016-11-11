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
    
    intCoef(
      const std::vector<node> &nodes, 	       
      const std::vector<triangularElement> &elements,
      const std::map<int, int> &node2coefficient,
      int numcoefficients); 
  public:
   
    
    double getInt(double *coefficients) const;
    
    static void initInstance();
    static intCoef *getInstance();
};

#endif	// _INTCOEF_H_