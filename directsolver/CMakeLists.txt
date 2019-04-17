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
  #src/cudacg.cu
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
  ../src/cuda/solver-cpcg.cu
  ../src/cuda/solver-cpcgcg.cu
  ../src/cuda/solvercuda.cu
  ../src/cuda/solvercuda.cpp
  #../src/cuda/solvercublas.cu
  #../src/cuda/solvercublas.cpp
  ../src/problem.cpp
  ../src/threedim/initproblem3D.cpp
  ../src/threedim/nodecoeficients3D.cpp
  ../src/twodim/initproblem2D.cpp
  ../src/twodim/nodecoeficients2D.cpp

)

if(WIN32)
    set(ADDITIONAL_LIBRARIES wsock32)
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
endif()

if(CUSPARSE)
	set(libs ${CMAKE_THREAD_LIBS_INIT} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
	add_executable(directsolver ${directsolver_SRCS} ../src/cuda/solvercublas.cu ../src/cuda/solvercublas.cpp src/cudacg.cu)
	target_compile_definitions(directsolver PRIVATE CUBLASCUSPARSE)
else()
	set(libs ${CMAKE_THREAD_LIBS_INIT})
	add_executable(directsolver ${directsolver_SRCS})
endif()
set_property(TARGET directsolver PROPERTY CXX_STANDARD 17)
set_property(TARGET directsolver PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET directsolver PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(directsolver ${libs} ${ADDITIONAL_LIBRARIES})
if(KERNELTIMING)
	target_compile_definitions(directsolver PRIVATE CGTIMING)
endif()
if(COOPERATIVEGROUPS)
	target_compile_definitions(directsolver PRIVATE CGROUPS)
	set(CMAKE_CUDA_FLAGS "-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=compute_70 -gencode arch=compute_75,code=compute_75 -rdc=true")
endif()