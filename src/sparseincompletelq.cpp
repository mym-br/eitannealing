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
    std::vector<element> buildingSQ;
    // building L vector
    std::map<int, double> buildingL;
    std::vector<element> buildingSL;
    
    this->m_matrix.startFill(il*A.cols());
    
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
      // Extract il largest L elements: first add all elements to a vector
      buildingSL.resize(il <= buildingL.size()?il:buildingL.size());
      std::partial_sort_copy(buildingL.begin(), buildingL.end(),buildingSL.begin(), buildingSL.end(),
	[](const element &e1, const element &e2){return fabs(e1.second)>fabs(e2.second);}
      );
      // Now re-sort them according to index FIXME: Is this step actually necessary?
      std::sort(buildingSL.begin(),buildingSL.end(),
	[](const element &e1, const element &e2){return e1.first<e2.first;});
      // And push them back to building matrix, while subtracting projections from previous vectors
      for(auto e:buildingSL) {
	// At most iq complexity
	for(auto q:QColumns[e.first]) {
	    buildingQ[q.first] -= e.second*q.second;	  
	}
	this->m_matrix.fill(e.first,col) = e.second;
      }
      // Extract largest elements from Q
      buildingSQ.resize(std::min(iq,(unsigned int)buildingQ.size()));
      std::partial_sort_copy(buildingQ.begin(), buildingQ.end(),buildingSQ.begin(), buildingSQ.end(),
	[](const element &e1, const element &e2){return fabs(e1.second)>fabs(e2.second);}
      );
      // Normalize it and update rows
      double normq2 = 0;
      for(auto q:buildingSQ){
	normq2 += q.second * q.second;	
      }
      double inormq = 1/sqrt(normq2);
      for(auto q:buildingSQ){
	q.second*=inormq;	
	// Update rows here
	QRows[q.first].push_back(element(col, q.second));
      }
      // And add it to built columns
      QColumns.push_back(buildingSQ);
      
    }
    //this->m_matrix.endFill();
    
    
}