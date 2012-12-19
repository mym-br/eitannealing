#ifndef GRADIENTNORMREGULARISATION_H
#define GRADIENTNORMREGULARISATION_H
#include <memory>
#include "solver.h"
#include "problemdescription.h"

class gradientNormRegularisation
{

private:
    gradientNormRegularisation();
    static std::unique_ptr<gradientNormRegularisation> instance;
    
    
    std::unique_ptr<matrix> regularizationMatrix;
    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> n2cmatrix;
    std::unique_ptr<n2cmatrix> adaptMatrix;
    static n2cmatrix *buildCoefficient2NodeMatrix(triangularEletrode *);
    triangularEletrode lastElectrode;
public:
    double getRegularisation(const float *sol) const;
  
   static void initInstance();
   static gradientNormRegularisation *getInstance() {
   
     return instance.get();
  }
};

#endif // GRADIENTNORMREGULARISATION_H
