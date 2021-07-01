from setuptools import Extension, setup

USE_CYTHON = ...  # command line option, try-import, ...

ext = '.pyx' if USE_CYTHON else '.cpp'

sourcefiles = [
    "pyeitsolver" + ext,
    "src/pyeitsolver.cpp",
    "src/gradientnormregularisation.cpp",
    "src/intcoef.cpp",
    "src/incomplete_cholesky.cpp",
    "src/incomplete_ldlt.cpp",
    "src/problem.cpp",
    "src/solution.cpp",
    "src/solutionbase.cpp",
    "src/solver.cpp",
    "src/mt19937-64/mt19937-64.c",
    "src/threedim/initproblem3D.cpp",
    "src/solver_lb.cpp",
    "src/solution_lb.cpp",
    "src/threedim/nodecoeficients3D.cpp",
    "src/twodim/initproblem2D.cpp",
    "src/twodim/nodecoeficients2D.cpp",
    "src/complex/incomplete_choleskycomplex.cpp",
    "src/complex/solvercomplex.cpp",
    "src/complex/solutioncomplex.cpp",
    "src/complex/gradientnormregularisationcomplex.cpp",
]

extensions = [
    Extension("pyeitsolver", sourcefiles, include_dirs=['src', 'eigen-3.3.5'])
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(ext_modules=extensions)