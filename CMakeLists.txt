cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(eitannealingtests)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 3.3.5 REQUIRED)  # We only need this version for CUDA
find_package(Threads REQUIRED)
find_package(CUDA 10 REQUIRED)
find_package(GMSH REQUIRED)
find_package(SWIG REQUIRED)
find_package(Python3 COMPONENTS NumPy)

set(CMAKE_AUTOMOC ON)
include(${SWIG_USE_FILE})
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GMSH_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIR} src src/complex)

set(eit_core_SRCS
  src/gradientnormregularisation.cpp
  src/intcoef.cpp
  src/incomplete_cholesky.cpp
  src/incomplete_ldlt.cpp
  #src/eigen3test.cpp
  src/problem.cpp
  src/solution.cpp
  src/solutionbase.cpp
  src/solver.cpp
  src/mt19937-64/mt19937-64.c
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
  src/mm/mmio.c
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

set(CMAKE_SWIG_FLAGS "")

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GMSH_LIBRARIES})

add_library(eit_core OBJECT ${eit_core_SRCS})
set_property(TARGET eit_core PROPERTY CXX_STANDARD 17)

add_executable(eitmain src/main.cpp src/parameters/parametersparser.cpp ${eit_gui_SRCS} $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitmain  ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitmain PROPERTY CXX_STANDARD 17)
set_property(TARGET eitmain PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(eitwosteps src/maintwosteps.cpp src/parameters/parametersparser.cpp ${eit_gui_SRCS} $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitwosteps  ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitwosteps PROPERTY CXX_STANDARD 17)
set_property(TARGET eitwosteps PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(eitlb src/mainlb.cpp src/parameters/parametersparser.cpp ${eit_gui_SRCS} $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitlb  ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitlb PROPERTY CXX_STANDARD 17)
set_property(TARGET eitlb PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(eitdirect src/maindirect.cpp $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitdirect  ${libs} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitdirect PROPERTY CXX_STANDARD 17)
set_property(TARGET eitdirect PROPERTY CXX_STANDARD_REQUIRED ON)

cuda_add_executable(eitcuda src/maincuda.cpp src/parameters/parametersparser.cpp ${eit_cuda_SRCS} ${eit_gui_SRCS} $<TARGET_OBJECTS:eit_core>)
target_link_libraries(eitcuda ${libs} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${ADDITIONAL_LIBRARIES})
set_property(TARGET eitmain PROPERTY CXX_STANDARD 17)
set_property(TARGET eitcuda PROPERTY CUDA_SEPARABLE_COMPILATION ON)

IF(CMAKE_BUILD_TYPE MATCHES RELEASE)
    set_source_files_properties(src/python/eitdirectsolver.i PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(src/python/eitdirectsolver.i PROPERTIES SWIG_FLAGS "-includeall")
    SWIG_ADD_LIBRARY(eitdirectsolver LANGUAGE python SOURCES src/python/eitdirectsolver.i src/python/eitdirectsolver.cpp $<TARGET_OBJECTS:eit_core>)
    SWIG_LINK_LIBRARIES(eitdirectsolver ${Python3_LIBRARIES} ${libs} ${ADDITIONAL_LIBRARIES})
ENDIF(CMAKE_BUILD_TYPE MATCHES RELEASE)

add_subdirectory(tests)