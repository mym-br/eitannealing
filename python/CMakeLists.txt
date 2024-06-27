# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 3.3.5 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR} ../src)

set(eit_core_SRCS
  ../src/gradientnormregularisation.cpp
  ../src/intcoef.cpp
  ../src/incomplete_cholesky.cpp
  ../src/incomplete_ldlt.cpp
  ../src/problem.cpp
  ../src/solution.cpp
  ../src/solutionbase.cpp
  ../src/solver.cpp
  ../src/mt19937-64/mt19937-64.c
  ../src/threedim/initproblem3D.cpp
  ../src/solver_lb.cpp
  ../src/solution_lb_real.cpp
  ../src/threedim/nodecoeficients3D.cpp
  ../src/twodim/initproblem2D.cpp
  ../src/twodim/nodecoeficients2D.cpp
)

add_library(eit_core OBJECT ${eit_core_SRCS})
set_property(TARGET eit_core PROPERTY CXX_STANDARD 17)

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(_core MODULE src/core.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)

# This is passing in the version as a define just as an example
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION pyeitsolver)
