#include <iostream>
#include <vector>
#include "basematrix.h"
#include "solver.h"
#include "incomplete_cholesky.h"


typedef Eigen::Triplet<double> T;


int main(int argc, char *argv[])
{
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
    
    precond.halfSolveInPlaceT(w);
    
    std::cout << w << std::endl;
    
    return 0;
}