#ifndef GRADIENTNORMREGULARISATION_H
#define GRADIENTNORMREGULARISATION_H
#include <memory>
#include "solver.h"
//#include "problemdescription.h"
#include "problem.h"
#include <set>
class gradientNormRegularisation
{

private:
    static std::unique_ptr<gradientNormRegularisation> instance;
    int electrodecoefficients;
    std::unique_ptr<matrix> regularizationMatrix;
    void buildMatrix();
    int coefficientMap(int node);
    std::set<int> lastElectrodeNodes;    
	std::shared_ptr<problem> input;
public:
  	gradientNormRegularisation(std::shared_ptr<problem> _input);
    double getRegularisation(const double *sol) const;
  
	static void initInstance(std::shared_ptr<problem> _input);
   static gradientNormRegularisation *getInstance() {
   
     return instance.get();
  }
};

//class gradientNormRegularisation_old
//{
//
//private:
//    gradientNormRegularisation_old();
//    static std::unique_ptr<gradientNormRegularisation_old> instance;
//    
//    
//    std::unique_ptr<matrix> regularizationMatrix;
//    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> n2cmatrix;
//    std::unique_ptr<n2cmatrix> adaptMatrix;
//    static n2cmatrix *buildCoefficient2NodeMatrix(genericEletrode **);
//    std::set<int> lastElectrodeNodes;    
//public:
//    double getRegularisation(const double *sol) const;
//  
//   static void initInstance();
//   static gradientNormRegularisation_old *getInstance() {
//   
//     return instance.get();
//  }
//};

#endif // GRADIENTNORMREGULARISATION_H
