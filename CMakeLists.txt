 project(eitannealingtests)
cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Qt5Widgets REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

set(eitannealingtest_SRCS 
#  src/gradientnormregularisation.cpp 
#  src/sparseincompletelq.cpp
#  src/graphics.cpp
#  src/initproblem.cpp
#  src/main.cpp
#  src/nodecoeficients.cpp
#  src/solver.cpp
  src/incomplete_cholesky.cpp
#  src/solution.cpp
#  src/mt19937-64/mt19937-64.c
#  src/observations.cpp
#  src/solver_lb.cpp
#  src/solution_lb.cpp
  src/eigen3test.cpp
)

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set(CMAKE_CXX_FLAGS ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

add_executable(eitannealingtest ${eitannealingtest_SRCS})
set_property(TARGET eitannealingtest PROPERTY CXX_STANDARD 11)
target_link_libraries(eitannealingtest ${libs})
