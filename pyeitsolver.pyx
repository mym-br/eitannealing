# distutils: language=c++
 
cimport cpyeitsolver
import pathlib
from libcpp.string cimport string
from libcpp.map cimport map
from cython.operator cimport dereference, postincrement

cpdef init(meshFilename, currentsFilename):
    if(pathlib.Path(meshFilename).exists() == False):
        raise Exception('Input mesh file not found. The value was: {}'.format(meshFilename))
    if(pathlib.Path(currentsFilename).exists() == False):
        raise Exception('Input currents file not found. The value was: {}'.format(currentsFilename))
    modelInfo = cpyeitsolver.init(meshFilename.encode(), currentsFilename.encode())
    
    # Convert map to dict with string keys
    modelInfoStr = {}
    cdef map[string,int].iterator it = modelInfo.begin()
    while(it != modelInfo.end()):
        modelInfoStr[dereference(it).first.decode()] = dereference(it).second
        postincrement(it)
    return modelInfoStr

cpdef forwardSolve(conductivities):
    return cpyeitsolver.solveForwardProblem(conductivities)

cpdef fullForwardSolve(conductivities):
    return cpyeitsolver.solveFullForwardProblem(conductivities)
