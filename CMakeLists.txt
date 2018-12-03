project(eitannealingtest)
cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)
find_package(GMSH REQUIRED)

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GMSH_INCLUDE_DIRS} src src/complex)

set(eit_core_SRCS
  src/gradientnormregularisation.cpp
  src/intcoef.cpp
  src/incomplete_cholesky.cpp
  src/incomplete_ldlt.cpp
  src/problem.cpp
  src/solution.cpp
  src/solutionbase.cpp
  src/solver.cpp
  src/mt19937-64/mt19937-64.c
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

set(eit_gui_SRCS
  src/graphics.cpp
  src/gmsh/gmshgraphics.cpp
)

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GMSH_LIBRARIES})
set(CMAKE_CXX_FLAGS ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

add_library(eit_core OBJECT ${eit_core_SRCS})
set_property(TARGET eit_core PROPERTY CXX_STANDARD 17)
set_property(TARGET eit_core PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(eitannealingtest tests/test_incomplete_lq.cpp ${eit_gui_SRCS})
set_property(TARGET eitannealingtest PROPERTY CXX_STANDARD 17)
target_link_libraries(eitannealingtest $<TARGET_OBJECTS:eit_core> ${libs} ${ADDITIONAL_LIBRARIES})
