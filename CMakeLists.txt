 project(eitannealingtests)
cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} src src/complex)

set(eitannealingtest_SRCS 
  src/gradientnormregularisation.cpp 
  src/graphics.cpp
  src/incomplete_cholesky.cpp
  src/incomplete_ldlt.cpp
  src/maincomplex.cpp
  src/problem.cpp
  src/solution.cpp
  src/solver.cpp
  src/mt19937-64/mt19937-64.c
  src/gmsh/gmshgraphics.cpp
  src/parameters/parametersparser.cpp
   src/threedim/initproblem3D.cpp
  src/threedim/nodecoeficients3D.cpp
  src/twodim/initproblem2D.cpp
  src/twodim/nodecoeficients2D.cpp
  src/complex/incomplete_choleskycomplex.cpp
  src/complex/solvercomplex.cpp
  src/complex/solutioncomplex.cpp
  src/complex/gradientnormregularisationcomplex.cpp
  src/complex/graphicscomplex.cpp
)

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
if(WIN32)
    set(ADDITIONAL_LIBRARIES wsock32)
endif()
set(CMAKE_CXX_FLAGS ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

add_executable(eitannealingtest ${eitannealingtest_SRCS})
set_property(TARGET eitannealingtest PROPERTY CXX_STANDARD 11)
target_link_libraries(eitannealingtest ${libs} ${ADDITIONAL_LIBRARIES})
