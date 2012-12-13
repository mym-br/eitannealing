#include "sparseincompletelq.h"
#include <Eigen/Sparse>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

using namespace boost::lambda;	// FIXME: Replace boost lambda with C++11 lambda
SparseIncompleteLQ::SparseIncompleteLQ(
  const Eigen::SparseMatrix<double, Eigen::SelfAdjoint | Eigen::LowerTriangular | Eigen::ColMajor >& A, unsigned int iq, unsigned int il)
{
    typedef std::pair<int, Scalar> element;
    typedef std::vector<std::vector<element> > VSV;
    // Q Matrix (As a vector of sparse vectors)
    VSV QColumns;
    // Q rows
    VSV QRows;
    // building Q vector
    std::map<int, double> buildingQ;
    std::vector<double> buildingSQ;
    // building L vector
    std::map<int, double> buildingL;
    std::vector<element> buildingSL;
    
    
    // first step: extract first rows
    
    for(int col=0;col<A.cols();col++) {
      buildingQ.clear();
      buildingL.clear();
      for(BaseMatrixType::InnerIterator it(A, col);it;++it) {
	// insert everything in current building column
	buildingQ[it.row()] = it.value();
	// And do a scalar product with previous columns... put result in building L
	for(std::vector<element>::iterator qit = QRows[it.row()].begin(); qit!=QRows[it.row()].end(); qit++) {
	      buildingL[qit->first] += it.value()*qit->second;
	}
      }
      // Extract l1 largest L elements: first add all elements to a vector
      buildingSL.clear();
      std::copy(buildingL.begin(), buildingL.end(), buildingSL.begin());
      // Get l1 largest elements
      if(il<buildingSL.size()) {
	//std::partial_sort(buildingSL.begin(),buildingSL.begin()+il, buildingSL.end(), 
	//   (bind<double>(fabs,&_1 ->* &element::second)>bind<double>(fabs,&_2 ->* &element::second)));
	std::partial_sort(buildingSL.begin(),buildingSL.begin()+il, buildingSL.end(), 
	   [](const element &e1, const element &e2){return fabs(e1.second)>fabs(e2.second);});
	buildingSL.resize(il);
      }

      // And push them back to building matrix
	
      // Now subtract projections with previous vectors
      for(VSV::iterator i=QColumns.begin(); i!=QColumns.end(); i++) {
	
      }
      // Extract largest elements
      // And normalize
      
      
      
    }
    
    // FIXME: should we SORT L vector?
    
    
}