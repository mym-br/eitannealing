#include "observations.h"
#include "gradientnormregularisation.h"
#include "solution_lb_complex2.h"
#include "util/EigenSymmQuadratic.h"
#include "util/standard_deviation.hpp"

#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

// Population variance specialization for complex values
template<> double population_variance<std::vector<std::complex<double> >::iterator>(std::vector<std::complex<double> >::iterator start, std::vector<std::complex<double> >::iterator end)
{
    struct moments {
        double re_s;
        double re_ss;
        double im_s;
        double im_ss;

        moments operator+(std::complex<double> x) {
            return { this->re_s + x.real(), this->re_ss + x.real()*x.real(), this->im_s + x.imag(), this->im_ss + x.imag()*x.imag()};
        }
    };
    unsigned n = end - start;
    moments m = std::accumulate(start, end, moments{0, 0, 0, 0});
    return (n*(m.re_ss + m.im_ss) - m.re_s*m.re_s - m.im_s*m.im_s)/(n*n);
}


static float mincond_I = (float)minperm;
static float maxcond_I = (float)maxperm;

// Shuffler specialization for real admitance
template <>
std::vector<std::complex<double> > solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::
	getShuffledSolution(shuffleData *data, const shuffler &sh) const
{
	std::vector<std::complex<double> > res(sol);
	bool shufflereal = genint(2);
	float minc, maxc;
	if(shufflereal) {
		maxc = (float)maxcond;
		minc = (float)mincond;
	} else {
		maxc = maxcond_I;
		minc = mincond_I;
	}
	// head or tails
	if(genint(2)) { // Normal
		int ncoef = genint(p->getNumCoefficients());
		int pcoef = ncoef + (shufflereal?0:p->getNumCoefficients());

		if(sh.shuffleConsts[pcoef]==0) {
			if(shufflereal)
				res[ncoef].real(minc+genreal()*(maxc-minc));
			else
				res[ncoef].imag(minc+genreal()*(maxc-minc));
		} else {
			double val;
			do {
				if(shufflereal)
					val = res[ncoef].real();
				else
					val = res[ncoef].imag();
				double rnd = 0;
				for(int i=0;i<sh.shuffleConsts[pcoef];i++)
					rnd += genreal();
				rnd /= sh.shuffleConsts[pcoef];
				rnd -= 0.5;
				rnd *= (maxc - minc);
				val += rnd;
			} while((val < minc) || (val > maxc));
			if(shufflereal)
				res[ncoef].real(val);
			else
				res[ncoef].imag(val);
		}
		if(data) {
			data->swap = false;
			data->ncoef = pcoef;
		}
	} else { // swap
		int ncoef = genint(p->getInnerAdjacencyCount());
		int pcoef = ncoef + (shufflereal?0:p->getInnerAdjacencyCount());
		int node1, node2;

		std::pair<int, int> adj = p->getAdjacency(ncoef);
                node1 = p->getNode2Coefficient(adj.first);
		node2 = p->getNode2Coefficient(adj.second);

		// Order nodes
		if(shufflereal)
			if(res[node1].real()>res[node2].real()) {
				int aux = node1;
				node1 = node2;;
				node2 = aux;
			}
		else
			if(res[node1].imag()>res[node2].imag()) {
				int aux = node1;
				node1 = node2;;
				node2 = aux;
			}
		double v1p, v2p, v1, v2;
		if(shufflereal) {
			v1p = v1 = res[node1].real();
			v2p = v2 = res[node2].real();
		} else {
			v1p = v1 = res[node1].imag();
			v2p = v2 = res[node2].imag();
		}
		double a = max( min(v1-minc, maxc-v2), min(maxc-v1, v2-minc));

		double delta;
		do {
			if(sh.swapshuffleconsts[pcoef]==0) {
				delta = a*(genreal()*2 - 1);
			} else {
				double rnd = 0;
				for(int i=0;i<sh.swapshuffleconsts[pcoef];i++)
					rnd += genreal();
				rnd /= sh.swapshuffleconsts[pcoef];
				rnd -= 0.5;
				delta = a*rnd;
			}
			v1 = v1p - delta;
			v2 = v2p + delta;
		} while((v1 < minc) || (v2 < minc) || (v1 > maxc) || (v2 > maxc));
		if(shufflereal) {
			res[node1].real(v1);
			res[node2].real(v2);
		} else {
			res[node1].imag(v1);
			res[node2].imag(v2);
		}
		if(data) {
			data->swap = true;
			data->ncoef = pcoef;
		}
	}
	return res;
}

template <>
std::vector<std::complex<double> > solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::
    getNewRandomSolution(int size)
{
	std::vector<std::complex<double> > res;
	res.reserve(size);
	int i = 0;
	for(;i<size;i++) {
		double re, im;
		re = mincond+genreal()*(maxcond-mincond);
		im = mincond_I+genreal()*(maxcond_I-mincond_I);
		res.emplace_back(re, im);
	}

	return res;
}

// Complex regularization

int complexGradientNormRegularisation::coefficientMap(int node)
{
	int c = input->getNode2Coefficient(node) - electrodecoefficients;
    return c>0?c:0;
}


void complexGradientNormRegularisation::buildMatrix()
{

	matrix *out = new matrix(input->getNumCoefficients() - electrodecoefficients, input->getNumCoefficients() - electrodecoefficients);
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

complexGradientNormRegularisation::complexGradientNormRegularisation(std::shared_ptr<problem> _input) : input(_input)
{
	electrodecoefficients = input->getGenericElectrodesCoeffCount();
    this->buildMatrix();
}

double complexGradientNormRegularisation::getRegularisation(const std::complex<double> *sol) const
{
    auto v = Eigen::VectorXcd::Map(sol + electrodecoefficients, input->getNumCoefficients() - electrodecoefficients);

	double re = EigenSymmQuadraticL<Scalar>(regularizationMatrix->selfadjointView<Eigen::Lower>(), v.real());
	double im = EigenSymmQuadraticL<Scalar>(regularizationMatrix->selfadjointView<Eigen::Lower>(), v.imag());

	return re + im;
}

// Declare template specialization
template class solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>;
// FIXME: Intel Compiler somehow needs explicit declarations of those methods. GCC doesn't. Which one is right?
template solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::solution_lb_gen(std::shared_ptr<problem>, complexobservations const&, std::shared_ptr<complexGradientNormRegularisation>, std::vector< std::complex<double> >&&);
template solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::solution_lb_gen(std::shared_ptr<problem>, complexobservations const&, std::shared_ptr<complexGradientNormRegularisation>);
template void solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::saturate();
template solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler> *solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::shuffle(shuffleData *data, const shuffler &sh) const;
template void solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::improve();
template bool solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler>::compareWith(solution_lb_gen<LB_Solver_Complex2, std::complex<double>, complexobservations, complexGradientNormRegularisation, complexMatrixBuilder, shuffleData, shuffler> &, float, float);

#include "solution_lb_impl.h"

