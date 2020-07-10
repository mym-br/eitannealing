from . import _pyeitsolver
import pathlib

def init(meshFilename, currentsFilename):
    if(pathlib.Path(meshFilename).exists() == False):
        raise Exception('Input mesh file not found. The value was: {}'.format(meshFilename))
    if(pathlib.Path(currentsFilename).exists() == False):
        raise Exception('Input currents file not found. The value was: {}'.format(currentsFilename))

    return _pyeitsolver.init(meshFilename, currentsFilename)

def forwardSolve(conds):
    return _pyeitsolver.solveForwardProblem(conds)