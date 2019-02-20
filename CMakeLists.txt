project(eitannealingtests)
cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_AUTOMOC ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

set(eit_core_SRCS
  src/gradientnormregularisation.cpp
  src/initproblem.cpp
  src/nodecoeficients.cpp
  src/solver.cpp
  src/incomplete_cholesky.cpp
  src/incomplete_ldlt.cpp
  src/incomplete_qr_builder.cpp
  src/solution.cpp
  src/mt19937-64/mt19937-64.c
  src/observations.cpp
  src/solver_lb.cpp
  src/solution_lb.cpp
  src/intcoef.cpp
  src/assembleMatrix.cpp
)

set(eit_gui_SRCS
  src/graphics.cpp
  src/gmsh/gmshgraphics.cpp
)

set(eitannealing_gui_SRCS
  src/graphics.cpp
)

add_subdirectory(tests)

add_library(eit_core OBJECT ${eit_core_SRCS})
set_property(TARGET eit_core PROPERTY CXX_STANDARD 17)

add_library(eit_gui OBJECT ${eitannealing_gui_SRCS})
target_compile_options(eit_gui PRIVATE ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})
set_property(TARGET eit_gui PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(eitannealingtest src/main.cpp $<TARGET_OBJECTS:eit_gui> $<TARGET_OBJECTS:eit_core>)
set_property(TARGET eitannealingtest PROPERTY CXX_STANDARD 17)
target_compile_definitions(eitannealingtest PRIVATE ${Qt5Widgets_DEFINITIONS})
target_link_libraries(eitannealingtest ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${ADDITIONAL_LIBRARIES})
