#ifndef _INTCOEF_H_
#define _INTCOEF_H_

#include <memory>

#include "basematrix.h"
#include "problem.h"
class intCoef {
  protected:
    vectorx x;
    int numcoefficients;
    
  public:
    intCoef(const problem &problem);    
    double getInt(double *coefficients) const;
};

#endif	// _INTCOEF_H_
