#include <iostream>
#include <vector>
#include "basematrix.h"
#include "solver.h"
#include "incomplete_cholesky.h"
#include "incomplete_LDLT.h"


typedef Eigen::Triplet<double> T;


const matrix buildInPlace() {
  matrix::StorageIndex *outer;
  matrix::StorageIndex *inner;
  Scalar *val;
  
  
  val = new Scalar[5];
  val[0] = 2; val[1] = -1; val[2] = 2; val[3] = 2; val[4] = 2;
  
  outer = new matrix::StorageIndex[4];
  outer[0] = 1; outer[1] = 2; outer[2] = 3; outer[3] = 4;
  
  inner = new matrix::StorageIndex[5];
  inner[0] = 0; inner[1] = 3; inner[2] = 1; inner[3] = 2; inner[4] = 3;
  
  return std::move(Eigen::Map<matrix>(4, 4, 5, outer, inner, val));    
}

int main(int argc, char *argv[])
{
    /*
  
    std::vector<T> tripletList;
    tripletList.push_back(T(0,0,2));
    tripletList.push_back(T(1,1,2));
    tripletList.push_back(T(2,2,2));
    tripletList.push_back(T(3,3,2));
    tripletList.push_back(T(2,0,-1));
    
    
    matrix m(4,4);
    m.setFromTriplets(tripletList.begin(), tripletList.end());
    
    SparseIncompleteLLT precond(m);
        
    
    Eigen::VectorXd w(Eigen::VectorXd::Ones(4));
    
    precond.solveInPlace(w);
    
    std::cout << w << std::endl << std::endl;
    
    SparseIncompleteLDLT precond2(m);
    
    w = Eigen::VectorXd::Ones(4);
    
    precond2.solveInPlace(w);
    
    std::cout << w << std::endl;
    */
    
    matrix m = buildInPlace();
    matrix n(m);
    
    std::cout << n.selfadjointView<Eigen::Lower>() << std::endl;
    return 0;
}