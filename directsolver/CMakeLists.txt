cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(cufppcgsolver LANGUAGES C CXX CUDA)
find_package(Eigen3 REQUIRED)
find_package(CUDAToolkit)
find_package(MKL CONFIG)

include_directories(${EIGEN3_INCLUDE_DIR} ${CUDAToolkit_INCLUDE_DIRS} src/common ../src ../src/complex)

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

set(cufppcgsolver_cuda_SRCS 
  ../src/cuda/solvercublas.cu
  ../src/cuda/solvercublas.cpp
  src/cufppcgsolver/cudacg.cu
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


if(CUDAToolkit_FOUND)
  set(cufppcgsolver_SRCS ${cufppcgsolver_SRCS} ${cufppcgsolver_cuda_SRCS})
endif()

add_executable(cufppcgsolver ${cufppcgsolver_SRCS} ${common_SRCS})
add_executable(cufppcgsolver_cgtiming ${cufppcgsolver_SRCS} ${common_SRCS})

set_target_properties(
  cufppcgsolver
  cufppcgsolver_cgtiming
    PROPERTIES 
    CXX_STANDARD 17
      CXX_STANDARD_REQUIRED ON
      CUDA_SEPARABLE_COMPILATION ON
)
if(CUDAToolkit_FOUND)
  set(libs CUDA::cudart CUDA::cublas CUDA::cusparse)
  target_compile_definitions(cufppcgsolver PRIVATE CUBLASCUSPARSE)
  target_compile_definitions(cufppcgsolver_cgtiming PRIVATE CUBLASCUSPARSE)
  if(NOCOOPERATIVEGROUPS)
      set_target_properties(cufppcgsolver cufppcgsolver_cgtiming PROPERTIES CUDA_ARCHITECTURES OFF)
    else()
      target_compile_definitions(cufppcgsolver PRIVATE CGROUPS)
      target_compile_definitions(cufppcgsolver_cgtiming PRIVATE CGROUPS)
      set_target_properties(cufppcgsolver cufppcgsolver_cgtiming PROPERTIES CUDA_ARCHITECTURES "60;61;70;75")
  endif()
endif()
target_compile_definitions(cufppcgsolver_cgtiming PRIVATE CGTIMING)

target_link_libraries(cufppcgsolver ${libs})
target_link_libraries(cufppcgsolver_cgtiming ${libs})

if(CUDAToolkit_FOUND)
  add_executable(cusolver ${cusolver_SRCS} ${common_SRCS})
  target_link_libraries(cusolver PRIVATE CUDA::cudart CUDA::cusolver)
endif()
if(MKL_FOUND)
  add_executable(pardiso ${pardiso_SRCS} ${common_SRCS})
  set_property(TARGET pardiso PROPERTY CXX_STANDARD 17)
  set_property(TARGET pardiso PROPERTY CXX_STANDARD_REQUIRED ON)
  target_include_directories(pardiso PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_options(pardiso PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_link_libraries(pardiso PUBLIC $<LINK_ONLY:MKL::MKL>)
endif()