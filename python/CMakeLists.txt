
# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  pyeitsolver
  VERSION 0.1
  LANGUAGES C CXX)

# Option to select between PYBIND11_NOPYTHON and the new version
option(USE_PYBIND11_NOPYTHON "Use PYBIND11_NOPYTHON for building" OFF)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

if (USE_PYBIND11_NOPYTHON)
  set(PYBIND11_NOPYTHON ON)
else()
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
endif()

set(PETSC_DIR "" CACHE STRING "PETSc dir.")
set(PETSC_ARCH "" CACHE STRING "PETSc arch.")

find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(PETSc)

# No validation required as PETSc does this internally
message(STATUS "Found PETSc ${PETSc_VERSION}")
if (PETSC_FOUND)
  list(APPEND PETSC_DEFINITIONS -DUSE_PETSC=1)
endif()

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR} ${PYTHON_INCLUDE_DIR} ${PETSC_INCLUDES} ../src)
add_definitions(${PETSC_DEFINITIONS})

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

if (USE_PYBIND11_NOPYTHON)
  # Create the Python extension core module manually
  add_library(_core MODULE src/core.cpp $<TARGET_OBJECTS:eit_python_core>)
  set_target_properties(_core PROPERTIES PREFIX "" SUFFIX ".pyd" CXX_STANDARD 17)
  target_link_libraries(_core PRIVATE pybind11::module ${PYTHON_LIBRARY})

  # Create the Python extension complex module manually
  add_library(_complex MODULE src/complex.cpp $<TARGET_OBJECTS:eit_python_core>)
  set_target_properties(_complex PROPERTIES PREFIX "" SUFFIX ".pyd" CXX_STANDARD 17)
  target_link_libraries(_complex PRIVATE pybind11::module ${PYTHON_LIBRARY} ${PETSC_LIBRARIES})
else()
  # Add core library using FindPython's tooling
  python_add_library(_core MODULE src/core.cpp $<TARGET_OBJECTS:eit_python_core> WITH_SOABI)
  target_link_libraries(_core PRIVATE pybind11::headers)
  
  # Add complex library using FindPython's tooling
  python_add_library(_complex MODULE src/complex.cpp $<TARGET_OBJECTS:eit_python_core> WITH_SOABI)
  target_link_libraries(_complex PRIVATE pybind11::headers ${PETSC_LIBRARIES})
endif()

set(CMAKE_CXX_FLAGS "-fPIC")

# This is passing in the version as a define just as an example
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION pyeitsolver)

# This is passing in the version as a define just as an example
target_compile_definitions(_complex PRIVATE VERSION_INFO=${PROJECT_VERSION} BLOCKGND=1)

target_link_options(_complex PUBLIC -fPIC)

# The install directory is the output (wheel) directory
install(TARGETS _complex DESTINATION pyeitsolver)
