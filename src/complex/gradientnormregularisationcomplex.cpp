#include "gradientnormregularisationcomplex.h"
//#include "problemdescription.h"
#include "solver.h"
#include <algorithm>
#include <iostream>
//#include "assembleMatrix.h"

#include "util/EigenSymmQuadratic.h"

std::unique_ptr<gradientNormRegularisationComplex> gradientNormRegularisationComplex::instance;

int gradientNormRegularisationComplex::coefficientMap(int node)
{
	int c = input->node2coefficient[node] - electrodecoefficients;
    return c>0?c:0;
}


void gradientNormRegularisationComplex::buildMatrix()
{
    
	matrixcomplex *out = new matrixcomplex(input->numcoefficients - electrodecoefficients, input->numcoefficients - electrodecoefficients);
    std::vector<Eigen::Triplet<Scalar>> tripletList; 
	for (int i = 0; i<input->getNodesCount() - 1; ++i) {
	int ci = coefficientMap(i);

	for (nodeCoefficients *aux = input->getNodeCoefficients()[i]; aux; aux = aux->next) {
	    int cj = coefficientMap(aux->node);
	    if(ci>cj) continue; // skip upper triangular
	    tripletList.push_back(Eigen::Triplet<Scalar>(cj, ci, aux->coefficient));
	}
    }
    out->setFromTriplets(tripletList.begin(), tripletList.end());
    out->makeCompressed();
    regularizationMatrix.reset(out);
}

gradientNormRegularisationComplex::gradientNormRegularisationComplex(std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> _input) : input(_input)
{
	electrodecoefficients = input->getGenericElectrodesCount();
    this->buildMatrix();
}

void gradientNormRegularisationComplex::initInstance(std::shared_ptr<problem<Complex, Eigen::VectorXcd, matrixcomplex>> _input)
{
	instance.reset(new gradientNormRegularisationComplex(_input));
}

std::complex<double> gradientNormRegularisationComplex::getRegularisation(const std::complex<double> *sol) const
{
	return EigenSymmQuadraticL<Complex>(regularizationMatrix->selfadjointView<Eigen::Lower>(),
		Eigen::VectorXcd::Map(sol + electrodecoefficients, input->numcoefficients - electrodecoefficients));
}
