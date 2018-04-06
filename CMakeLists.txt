project(eitannealingtest)
cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)

if(NOT DEFINED GMSH_INCLUDE_DIR)
if(NOT DEFINED DO_NOT_USE_NATIVE_GMSH)
    find_path(GMSH_INCLUDE_DIR NAMES onelab.h PATH_SUFFIXES gmsh)
    if (GMSH_INCLUDE_DIR)
        message(STATUS "LibGmsh includes found at ${GMSH_INCLUDE_DIR}.")
    else()
        message(STATUS "LibGmsh includes not found, falling back to redist.")
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/gmsh) 
    endif()
else()
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/gmsh) 
endif()
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GMSH_INCLUDE_DIR} src src/complex)

set(eitannealingtest_SRCS
  src/gradientnormregularisation.cpp
  src/graphics.cpp
  src/incomplete_cholesky.cpp
  src/incomplete_ldlt.cpp
  tests/test_incomplete_lq.cpp
  src/problem.cpp
  src/solution.cpp
  src/solutionbase.cpp
  src/solver.cpp
  src/mt19937-64/mt19937-64.c
  src/gmsh/gmshgraphics.cpp
  src/parameters/parametersparser.cpp
  src/threedim/initproblem3D.cpp
  src/solver_lb.cpp
  src/solution_lb.cpp
  src/threedim/nodecoeficients3D.cpp
  src/twodim/initproblem2D.cpp
  src/twodim/nodecoeficients2D.cpp
  src/complex/incomplete_choleskycomplex.cpp
  src/complex/solvercomplex.cpp
  src/complex/solutioncomplex.cpp
  src/complex/gradientnormregularisationcomplex.cpp
)

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
if(WIN32)
    set(ADDITIONAL_LIBRARIES wsock32)
endif()
set(CMAKE_CXX_FLAGS ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

add_executable(eitannealingtest ${eitannealingtest_SRCS})
set_property(TARGET eitannealingtest PROPERTY CXX_STANDARD 17)
target_link_libraries(eitannealingtest ${libs} ${ADDITIONAL_LIBRARIES})
