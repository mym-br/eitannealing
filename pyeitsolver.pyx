# distutils: language=c++
 
cimport cpyeitsolver
import pathlib

cpdef init(meshFilename, currentsFilename):
    if(pathlib.Path(meshFilename).exists() == False):
        raise Exception('Input mesh file not found. The value was: {}'.format(meshFilename))
    if(pathlib.Path(currentsFilename).exists() == False):
        raise Exception('Input currents file not found. The value was: {}'.format(currentsFilename))
    return cpyeitsolver.init(meshFilename, currentsFilename)

cpdef forwardSolve(conductivities):
    return cpyeitsolver.solveForwardProblem(conductivities)

cpdef fullForwardSolve(conductivities):
    return cpyeitsolver.solveFullForwardProblem(conductivities)
