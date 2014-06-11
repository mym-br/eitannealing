#include "gradientnormregularisation.h"
#include "problemdescription.h"
#include "solver.h"
#include <algorithm>
#include <iostream>
#include <boost/concept_check.hpp>

std::unique_ptr<gradientNormRegularisation> gradientNormRegularisation::instance;

gradientNormRegularisation::n2cmatrix *gradientNormRegularisation::buildCoefficient2NodeMatrix(genericEletrode **last)
{
    n2cmatrix *m = new n2cmatrix(nodes.size()-1,numcoefficients);
    
    m->startFill();
    for(int i = 0; i<nodes.size()-1; i++)  {
      int coefficient = node2coefficient[i];
      if(coefficient <32) {// electrode node
	auto ee = std::find_if(gelectrodes.begin(), gelectrodes.end(), [i](const genericEletrode &e) {
	  return e.baseNode == i;  
	});
	if(ee!= gelectrodes.end())
	  coefficient = node2coefficient[ee->nodesPairs.begin()->first];
	else continue;
      }
      m->fill(i,coefficient)=1;
    }
    m->endFill();
    
    // now get last electrode
    *last = &(*std::find_if(gelectrodes.begin(), gelectrodes.end(), [](const genericEletrode &e) {
	  return e.baseNode == nodes.size()-1;  
	}));
    return m;
}

gradientNormRegularisation::gradientNormRegularisation()
{
    std::vector<float> sol(numcoefficients, 1.0);
    matrix *base;
    assembleProblemMatrix(&sol[0], &base);
    this->regularizationMatrix.reset(base);
    genericEletrode *last;
    this->adaptMatrix.reset(gradientNormRegularisation::buildCoefficient2NodeMatrix(&last));
    for(std::pair<int, int> & nodes : last->nodesPairs) {
      this->lastElectrodeNodes.insert(nodes.first);
      this->lastElectrodeNodes.insert(nodes.second);
    }
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
    for(int node : lastElectrodeNodes) {
      d[node] = 0;
    }
    return c.dot(d);
}