#include <iostream>
#include "../src/incomplete_qr_builder.h"
#include "../src/problem.h"
#include "../src/util/fill_with_smallest.hpp"
#include <vector>

class upperTriangularElementsMap {
    std::vector<std::vector<std::pair<Eigen::Index, Eigen::Index>>> UTCols;
    // FIXME: assume compressed!!!!!
public:
    template<class scalar> upperTriangularElementsMap(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m) {
        UTCols.resize(m.cols());
        for(Eigen::Index j = 0; j < m.outerSize(); j++) 
            for(Eigen::Index idx = m.outerIndexPtr()[j]+1; idx < m.outerIndexPtr()[j+1]; idx++)
                UTCols[m.innerIndexPtr()[idx]].push_back(std::make_pair(j, idx));
    }
    
    void iterateOverUpperElements(unsigned int j, std::function<void(Eigen::Index,Eigen::Index)> &&f) const {
        for(auto [i, offset] : UTCols[j]) f(i, offset);
    }
};

template<class scalar> class symmetric2ColStorage {
    const upperTriangularElementsMap &utmap;
    const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m;
public:
    symmetric2ColStorage(const Eigen::SparseMatrix<scalar, Eigen::ColMajor> &m, const upperTriangularElementsMap &map): utmap(map), m(m) {}

    void iterateOverColumn(unsigned long j, std::function<void(unsigned long, scalar)> &&f) const {
        // 1st iterate over upper elements
        utmap.iterateOverUpperElements(j, [f, this](Eigen::Index i, Eigen::Index offset) {        
            f(i, m.valuePtr()[offset]);
        });
        (typename SparseIncompleteQRBuilder<scalar>::columnMajorStorageAdaptor(m)).iterateOverColumn(j,std::move(f));
    }
    unsigned long rows() const { return m.rows(); }
    unsigned long cols() const { return m.cols(); }
};


int main(int argc, char **argv) {
  if(argc!=2) {
    std::cerr << "Parameters must be mesh file name\n";
    return 1;
  }
  std::shared_ptr<problem> prob = problem::createNewProblem(argv[1], true);
  prob->initProblem(argv[1]);
  prob->buildNodeCoefficients();
  prob->prepareSkeletonMatrix();
  prob->createCoef2KMatrix();
  std::vector<double> cond(prob->getNumCoefficients(), 1.0);
  Eigen::SparseMatrix<double, Eigen::ColMajor> *a;
  prob->assembleProblemMatrix(&cond[0], &a);
  SparseIncompleteQRBuilder<double> builder;

  Eigen::SparseMatrix<double> R = builder.buildRMatrixFromColStorage(symmetric2ColStorage<double>(*a, upperTriangularElementsMap(*a)), 6, 12);
  return 0;
}
