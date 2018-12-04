#ifndef GRADIENTNORMREGULARISATIONCOMPLEX_H
#define GRADIENTNORMREGULARISATIONCOMPLEX_H
#include <memory>
#include "solver.h"
//#include "problemdescription.h"
#include "problem.h"
#include <set>
class gradientNormRegularisationComplex
{

private:
	gradientNormRegularisationComplex(std::shared_ptr<problem> _input, bool _calibrationmode);
	static std::unique_ptr<gradientNormRegularisationComplex> instance;
    int electrodecoefficients;
    std::unique_ptr<matrixcomplex> regularizationMatrix;
    void buildMatrix();
    int coefficientMap(int node);
    std::set<int> lastElectrodeNodes;    
	std::shared_ptr<problem> input;
	bool calibrationmode;
public:
    std::complex<double> getRegularisation(const std::complex<double> *sol) const;
  
	static void initInstance(std::shared_ptr<problem> _input);
	static void initCalibrationInstance(std::shared_ptr<problem> _input);
	static gradientNormRegularisationComplex *getInstance() {
   
     return instance.get();
  }
};

#endif // GRADIENTNORMREGULARISATIONCOMPLEX_H
