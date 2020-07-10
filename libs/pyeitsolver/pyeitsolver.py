from . import _pyeitsolver

def init(meshFilename, currentsFilename):
    return _pyeitsolver.init(meshFilename, currentsFilename)

def forwardSolve(conds):
    return _pyeitsolver.solveForwardProblem(conds)