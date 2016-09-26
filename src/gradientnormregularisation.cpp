#include "gradientnormregularisation.h"
#include "problemdescription.h"
#include "solver.h"
#include <algorithm>
#include <iostream>
#include "assembleMatrix.h"

#include "util/EigenSymmQuadratic.h"

std::unique_ptr<gradientNormRegularisation_old> gradientNormRegularisation_old::instance;

gradientNormRegularisation_old::n2cmatrix *gradientNormRegularisation_old::buildCoefficient2NodeMatrix(genericEletrode **last)
{
    std::vector<Eigen::Triplet<Scalar>> tripletList; 
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
      tripletList.push_back(Eigen::Triplet<Scalar>(i,coefficient,1));
    }
    n2cmatrix *m = new n2cmatrix(nodes.size()-1,numcoefficients);
    m->setFromTriplets(tripletList.begin(), tripletList.end());
    m->makeCompressed();
    
    // now get last electrode
    *last = &(*std::find_if(gelectrodes.begin(), gelectrodes.end(), [](const genericEletrode &e) {
	  return e.baseNode == nodes.size()-1;  
	}));
    return m;
}

gradientNormRegularisation_old::gradientNormRegularisation_old()
{
    std::vector<double> sol(numcoefficients, 1.0);
    matrix *base;
    assembleProblemMatrix(&sol[0], &base);
    this->regularizationMatrix.reset(base);
    genericEletrode *last;
    this->adaptMatrix.reset(gradientNormRegularisation_old::buildCoefficient2NodeMatrix(&last));
    for(std::pair<int, int> & nodes : last->nodesPairs) {
      this->lastElectrodeNodes.insert(nodes.first);
      this->lastElectrodeNodes.insert(nodes.second);
    }
}

void gradientNormRegularisation_old::initInstance()
{
    instance.reset(new gradientNormRegularisation_old());
}

double gradientNormRegularisation_old::getRegularisation(const double *sol) const
{
    Eigen::VectorXd s(Eigen::VectorXd::Map(sol, numcoefficients));
    Eigen::VectorXd c(*adaptMatrix*s);
    Eigen::VectorXd d(regularizationMatrix->selfadjointView<Eigen::Lower>()*c);
    for(int node : lastElectrodeNodes) {
      d[node] = 0;
    }
    return c.dot(d);
}


std::unique_ptr<gradientNormRegularisation> gradientNormRegularisation::instance;

int gradientNormRegularisation::coefficientMap(int node)
{
    int c = node2coefficient[node];
    return c>32?c:32;
}


void gradientNormRegularisation::buildMatrix()
{
    matrix *out = new matrix(numcoefficients, numcoefficients);
    std::vector<Eigen::Triplet<Scalar>> tripletList; 
    for (int i=0; i<nodes.size()-1; ++i) {
	int ci = coefficientMap(i);

	for(nodeCoefficients *aux = nodeCoef[i]; aux; aux = aux->next) {
	    int cj = coefficientMap(aux->node);
	    if(ci>cj) continue; // skip upper triangular
	    tripletList.push_back(Eigen::Triplet<Scalar>(cj, ci, aux->coefficient));
	}
    }
    out->setFromTriplets(tripletList.begin(), tripletList.end());
    out->makeCompressed();
    regularizationMatrix.reset(out);
}

gradientNormRegularisation::gradientNormRegularisation()
{
    this->buildMatrix();
}

void gradientNormRegularisation::initInstance()
{
    instance.reset(new gradientNormRegularisation());
}

double gradientNormRegularisation::getRegularisation(const double *sol) const
{
    return EigenSymmQuadraticL<Scalar>(regularizationMatrix->selfadjointView<Eigen::Lower>(),
				       Eigen::VectorXd::Map(sol, numcoefficients));
}
