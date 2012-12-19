#include "gradientnormregularisation.h"
#include "problemdescription.h"
#include "solver.h"
#include <algorithm>
#include <iostream>
#include <boost/concept_check.hpp>

std::unique_ptr<gradientNormRegularisation> gradientNormRegularisation::instance;

gradientNormRegularisation::n2cmatrix *gradientNormRegularisation::buildCoefficient2NodeMatrix(triangularEletrode *last)
{
    n2cmatrix *m = new n2cmatrix(nodes.size()-1,numcoefficients);
    
    m->startFill();
    for(int i = 0; i<nodes.size()-1; i++)  {
      int coefficient = node2coefficient[i];
      if(coefficient <32) {// electrode node
	auto ee = std::find_if(electrodes.begin(), electrodes.end(), [i](const triangularEletrode &e) {
	  return e.baseNode == i;  
	});
	if(ee!= electrodes.end())
	  coefficient = node2coefficient[ee->n2];
	else continue;
      }
      m->fill(i,coefficient)=1;
    }
    m->endFill();
    
    // now get last electrode
    *last = *std::find_if(electrodes.begin(), electrodes.end(), [](const triangularEletrode &e) {
	  return e.baseNode == nodes.size()-1;  
	});
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