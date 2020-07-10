from . import _pyeitsolver
import pathlib

def forwardSolve(conductivities, meshFilename, currentsFilename):
    if(pathlib.Path(meshFilename).exists() == False):
        raise Exception('Input mesh file not found. The value was: {}'.format(meshFilename))
    if(pathlib.Path(currentsFilename).exists() == False):
        raise Exception('Input currents file not found. The value was: {}'.format(currentsFilename))

    problemData = _pyeitsolver.init(meshFilename, currentsFilename)
    for cond in conductivities:
        if(len(cond) != problemData['coeffCount']):
            raise Exception('Incorrect number of conductivities, should be {}. The value was: {}'.format(problemData['coeffCount'], len(cond)))
        yield _pyeitsolver.solveForwardProblem(cond)