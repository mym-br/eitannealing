cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(directsolver LANGUAGES C CXX CUDA)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)
find_package(CUDA QUIET REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${EIGEN3_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} src ../src ../src/complex)

set(directsolver_SRCS 
  src/main.cpp
  src/conversions.cpp
  src/mm/mmio.c
  ../src/incomplete_cholesky.cpp
  ../src/solver.cpp
  ../src/cuda/solvercuda.cu
  ../src/cuda/solvercuda.cpp
  ../src/cuda/matrix-cpjds.cpp
  ../src/cuda/matrix-cpjds.cu
  ../src/cuda/number.cpp
  ../src/cuda/number.cu
  ../src/cuda/vector.cpp
  ../src/cuda/vector.cu
  ../src/cuda/utils.cpp
  ../src/cuda/color.cpp
  ../src/cuda/analysis.cpp
  ../src/cuda/solver-pcg.cu
  ../src/cuda/solvercuda.cu
  ../src/cuda/solvercuda.cpp
  ../src/cuda/solvercublas.cu
  ../src/cuda/solvercublas.cpp
  ../src/problem.cpp
  ../src/threedim/initproblem3D.cpp
  ../src/threedim/nodecoeficients3D.cpp
  ../src/twodim/initproblem2D.cpp
  ../src/twodim/nodecoeficients2D.cpp

)

set(libs ${CMAKE_THREAD_LIBS_INIT} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
if(WIN32)
    set(ADDITIONAL_LIBRARIES wsock32)
endif()

add_executable(directsolver ${directsolver_SRCS})
set_property(TARGET directsolver PROPERTY CXX_STANDARD 17)
set_property(TARGET directsolver PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET directsolver PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(directsolver ${libs} ${ADDITIONAL_LIBRARIES})
