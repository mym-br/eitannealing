# distutils: language=c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair

cdef extern from "pyeitsolver.h":
    map[string, int] init(const char* meshfilename, const  char* currentfilename) except +
    pair[int, vector[double]] solveForwardProblem(vector[double] conds) except +
    vector[double] solveFullForwardProblem(vector[double] conds) except +