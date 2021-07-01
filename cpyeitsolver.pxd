# distutils: language=c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map

cdef extern from "pyeitsolver.h":
    map[string, int] init(const char* meshfilename, const  char* currentfilename) except +
    vector[double] solveForwardProblem(vector[double] conds) except +
    vector[double] solveFullForwardProblem(vector[double] conds) except +