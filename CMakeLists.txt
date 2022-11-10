cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(eitannealingtests)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 3.3.5 REQUIRED)  # We only need this version for CUDA
find_package(Threads REQUIRED)
find_package(absl COMPONENTS container)
find_package(CUDA 10)
if(CUDA_FOUND)
  if(${CUDA_VERSION_MAJOR} GREATER 10) # FIXME: Cuda support is rotting fast
    message(WARNING "Cuda version greater than 10 is not supported. Cuda Target will be skipped.")
    unset(CUDA_FOUND)
  endif()
endif()
find_package(GMSH REQUIRED)

include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GMSH_INCLUDE_DIRS} src src/complex)

if(absl_FOUND)
  message("Using absl flat_hash_map")
  add_compile_definitions(USE_ABSL_FLAT_HASH)
  set(corelibs absl::flat_hash_map)
endif()

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
  src/threedim/initproblem3D.cpp
  src/solver_lb.cpp
  src/solution_lb_real.cpp
  src/threedim/nodecoeficients3D.cpp
  src/twodim/initproblem2D.cpp
  src/twodim/nodecoeficients2D.cpp
  src/complex/incomplete_choleskycomplex.cpp
  src/complex/solvercomplex.cpp
  src/complex/solutioncomplex.cpp
  src/complex/gradientnormregularisationcomplex.cpp
)

set(eit_new_complex_engine_SRCS
  src/solver_lb_complex.cpp
  src/solution_lb_complex.cpp
  src/solution_lb_complex2.cpp
)

set(eit_cuda_SRCS
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
  src/cuda/solvercuda.cu
  src/cuda/solvercuda.cpp
  src/cuda/solvercublas.cu
  src/cuda/solvercublas.cpp
  src/cuda/solutioncuda.cpp
)

set(eit_gui_SRCS
  src/graphics.cpp
  src/gmsh/gmshgraphics.cpp
)

IF (WIN32)
  set(eit_time_SRCS src/util/timestamp/win/timestamp.cpp)
ELSE()
  set(eit_time_SRCS src/util/timestamp/ux/timestamp.cpp)
ENDIF()

set(corelibs ${corelibs})
set(libs ${corelibs} ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GMSH_LIBRARIES})

add_library(eit_core OBJECT ${eit_core_SRCS})
set_property(TARGET eit_core PROPERTY CXX_STANDARD 17)

add_library(eit_new_complex_engine OBJECT ${eit_new_complex_engine_SRCS})
set_property(TARGET eit_new_complex_engine PROPERTY CXX_STANDARD 17)

add_library(eit_gui OBJECT ${eit_gui_SRCS})
target_compile_options(eit_gui PRIVATE ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})
set_property(TARGET eit_gui PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET eit_gui PROPERTY AUTOMOC ON)

add_library(eit_time OBJECT ${eit_time_SRCS})
set_property(TARGET eit_core PROPERTY CXX_STANDARD 17)

add_executable(eitmain src/main.cpp src/parameters/parametersparser.cpp $<TARGET_OBJECTS:eit_gui> $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitmain  ${libs} ${ADDITIONAL_LIBRARIES} )
set_property(TARGET eitmain PROPERTY CXX_STANDARD 17)
set_property(TARGET eitmain PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(eitwosteps src/maintwosteps.cpp src/parameters/parametersparser.cpp $<TARGET_OBJECTS:eit_gui> $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitwosteps  ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitwosteps PROPERTY CXX_STANDARD 17)
set_property(TARGET eitwosteps PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(eitlb src/mainlb.cpp src/parameters/parametersparser.cpp $<TARGET_OBJECTS:eit_gui> $<TARGET_OBJECTS:eit_core> $<TARGET_OBJECTS:eit_new_complex_engine>)
target_link_libraries(eitlb ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitlb PROPERTY CXX_STANDARD 17)
set_property(TARGET eitlb PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(eitdirect src/maindirect.cpp $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitdirect  ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitdirect PROPERTY CXX_STANDARD 17)
set_property(TARGET eitdirect PROPERTY CXX_STANDARD_REQUIRED ON)

if(${CUDA_FOUND})
    cuda_add_executable(eitcuda src/maincuda.cpp src/parameters/parametersparser.cpp ${eit_cuda_SRCS} $<TARGET_OBJECTS:eit_gui> $<TARGET_OBJECTS:eit_core>)
    target_link_libraries(eitcuda ${libs} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${ADDITIONAL_LIBRARIES})
    set_property(TARGET eitmain PROPERTY CXX_STANDARD 17)
    set_property(TARGET eitcuda PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

add_subdirectory(tests)
