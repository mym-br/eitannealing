#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/control_structures.hpp>
#include <boost/function.hpp>
#include <iostream>

#include "gradientnormregularisation.h"
#include "problemdescription.h"
#include "solver.h"

std::auto_ptr<gradientNormRegularisation> gradientNormRegularisation::instance;

gradientNormRegularisation::n2cmatrix *gradientNormRegularisation::buildCoefficient2NodeMatrix(triangularEletrode *last)
{
    using namespace boost::lambda;
    n2cmatrix *m = new n2cmatrix(nodes.size()-1,numcoefficients);
    
    m->startFill();
    for(int i = 0; i<nodes.size()-1; i++)  {
      int coefficient = node2coefficient[i];
      if(coefficient <32) {// electrode node
	 std::vector<triangularEletrode>::iterator ee = 
	    std::find_if(electrodes.begin(), electrodes.end(),
			 ((&_1 ->* &triangularEletrode::baseNode)==i));
	if(ee!= electrodes.end())
	  coefficient = node2coefficient[ee->n2];
	else continue;
      }
      m->fill(i,coefficient)=1;
    }
    m->endFill();
    
    // now get last electrode
    *last = *std::find_if(electrodes.begin(), electrodes.end(),
		 ((&_1 ->* &triangularEletrode::baseNode)==nodes.size()-1));
    return m;
}

gradientNormRegularisation::gradientNormRegularisation()
{
    std::vector<float> sol(numcoefficients, 1.0);
    matrix *base;
    assembleProblemMatrix(&sol[0], &base);
    this->regularizationMatrix.reset(base);
    this->adaptMatrix.reset(gradientNormRegularisation::buildCoefficient2NodeMatrix(&this->lastElectrode));
}

void gradientNormRegularisation::initInstance()
{
    instance.reset(new gradientNormRegularisation());
}

double gradientNormRegularisation::getRegularisation(const float *sol) const
{
    Eigen::VectorXd s(Eigen::VectorXf::Map(sol, numcoefficients).cast<double>());
    Eigen::VectorXd c(*adaptMatrix*s);
    Eigen::VectorXd d(*regularizationMatrix*c);
    d[lastElectrode.n1] = 0;
    d[lastElectrode.n2] = 0;
    d[lastElectrode.n3] = 0;
    return c.dot(d);
}