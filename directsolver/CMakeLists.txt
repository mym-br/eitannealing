cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(cufppcgsolver LANGUAGES C CXX CUDA)
find_package(Eigen3 REQUIRED)
find_package(CUDAToolkit REQUIRED)
find_package(MKL CONFIG REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} src/common ../src ../src/complex)

set(common_SRCS
  src/common/mm/mmio.c
)

set(cufppcgsolver_SRCS 
  src/cufppcgsolver/main.cpp
  src/cufppcgsolver/conversions.cpp
  ../src/incomplete_cholesky.cpp
  ../src/solver.cpp
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
  ../src/problem.cpp
  ../src/threedim/initproblem3D.cpp
  ../src/threedim/nodecoeficients3D.cpp
  ../src/twodim/initproblem2D.cpp
  ../src/twodim/nodecoeficients2D.cpp
)

set(cusolver_SRCS 
  src/cusolver/main.cpp
  src/cusolver/mmio_wrapper.cpp
)

set(pardiso_SRCS 
  src/pardiso/main.cpp
)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
endif()

if(CUSPARSE)
	set(libs CUDA::cudart CUDA::cublas CUDA::cusparse)
	add_executable(cufppcgsolver ${cufppcgsolver_SRCS} ${common_SRCS} ../src/cuda/solvercublas.cu ../src/cuda/solvercublas.cpp src/cufppcgsolver/cudacg.cu)
	target_compile_definitions(cufppcgsolver PRIVATE CUBLASCUSPARSE)
else()
	add_executable(cufppcgsolver ${cufppcgsolver_SRCS})
endif()
set_property(TARGET cufppcgsolver PROPERTY CXX_STANDARD 17)
set_property(TARGET cufppcgsolver PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET cufppcgsolver PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(cufppcgsolver ${libs})
if(KERNELTIMING)
	target_compile_definitions(cufppcgsolver PRIVATE CGTIMING)
endif()
if(COOPERATIVEGROUPS)
	target_compile_definitions(cufppcgsolver PRIVATE CGROUPS)
  set_target_properties(cufppcgsolver PROPERTIES CUDA_ARCHITECTURES "60;61;70;75")
	set(CMAKE_CUDA_FLAGS "-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=compute_70 -gencode arch=compute_75,code=compute_75")
else()
  set_target_properties(cufppcgsolver PROPERTIES CUDA_ARCHITECTURES OFF)
endif()

if(CUSOLVER)
  add_executable(cusolver ${cusolver_SRCS} ${common_SRCS})
  target_link_libraries(cusolver PRIVATE CUDA::cudart CUDA::cusolver)
endif()
if(PARDISO)
  add_executable(pardiso ${pardiso_SRCS} ${common_SRCS})
  set_property(TARGET pardiso PROPERTY CXX_STANDARD 17)
  set_property(TARGET pardiso PROPERTY CXX_STANDARD_REQUIRED ON)
  target_include_directories(pardiso PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_options(pardiso PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_link_libraries(pardiso PUBLIC $<LINK_ONLY:MKL::MKL>)
endif()