# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES C CXX)

set(PYBIND11_NOPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 3.3.5 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR} ${PYTHON_INCLUDE_DIR} ../src)

set(eit_python_core_SRCS
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

add_library(eit_python_core OBJECT ${eit_python_core_SRCS})
set_property(TARGET eit_python_core PROPERTY CXX_STANDARD 17)

# Create the Python extension module manually
add_library(_core MODULE src/core.cpp $<TARGET_OBJECTS:eit_python_core>)
set_target_properties(_core PROPERTIES PREFIX "" SUFFIX ".pyd" CXX_STANDARD 17)
target_link_libraries(_core PRIVATE pybind11::module ${PYTHON_LIBRARY})

# This is passing in the version as a define just as an example
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION pyeitsolver)

# Create the Python extension module manually
add_library(_complex MODULE src/complex.cpp $<TARGET_OBJECTS:eit_python_core>)
set_target_properties(_complex PROPERTIES PREFIX "" SUFFIX ".pyd" CXX_STANDARD 17)
target_link_libraries(_complex PRIVATE pybind11::module ${PYTHON_LIBRARY})

# This is passing in the version as a define just as an example
target_compile_definitions(_complex PRIVATE VERSION_INFO=${PROJECT_VERSION} BLOCKGND=1)

# The install directory is the output (wheel) directory
install(TARGETS _complex DESTINATION pyeitsolver)
