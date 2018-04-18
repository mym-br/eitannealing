cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(eitannealingtests LANGUAGES C CXX CUDA)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)
find_package(CUDA QUIET REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} src src/complex)

set(eitannealingtest_SRCS 
  src/gradientnormregularisation.cpp 
  src/graphics.cpp
  src/incomplete_cholesky.cpp
  src/incomplete_ldlt.cpp
  src/eigen3test.cpp
  src/problem.cpp
  src/solution.cpp
  src/solutionbase.cpp
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
  src/cuda/matrix-cpjds.cpp
  src/cuda/matrix-cpjds.cu
  src/cuda/number.cpp
  src/cuda/number.cu
  src/cuda/vector.cpp
  src/cuda/vector.cu
  src/cuda/utils.cpp
  src/cuda/color.cpp
  src/cuda/analysis.cpp
  src/cuda/solver-pcg.cu
  src/cuda/HighResClock.cpp
  src/cuda/solvercuda.cu
  src/cuda/solvercuda.cpp
  src/cuda/solvercublas.cu
  src/cuda/solvercublas.cpp
)

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
if(WIN32)
    set(ADDITIONAL_LIBRARIES wsock32)
endif()
set(CMAKE_CXX_FLAGS ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

add_executable(eitannealingtest ${eitannealingtest_SRCS})
target_compile_features(eitannealingtest PUBLIC cxx_std_11)
set_property(TARGET eitannealingtest PROPERTY CXX_STANDARD 11)
set_property(TARGET eitannealingtest PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(eitannealingtest ${libs} ${ADDITIONAL_LIBRARIES})
