 project(eitannealingtests)
cmake_minimum_required(VERSION 2.8.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Qt5Widgets REQUIRED)
#find_package(Boost REQUIRED COMPONENTS thread)
find_package(Eigen2 REQUIRED)
find_package(CXX11)

MESSAGE("Qt include:" ${Qt5Widgets_INCLUDE_DIRS})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
include_directories(${Qt5Widgets_INCLUDE_DIRS} ${EIGEN2_INCLUDES})

set(eitannealingtest_SRCS src/gradientnormregularisation.cpp src/sparseincompletelq.cpp
  src/graphics.cpp
  src/initproblem.cpp
  src/main.cpp
  src/nodecoeficients.cpp
  src/solver.cpp
  src/incomplete_cholesky.cpp
  src/solution.cpp
  src/mt19937-64/mt19937-64.c
  src/observations.cpp
  src/solver_lb.cpp
  src/solution_lb.cpp
)

add_definitions(${Qt5Widgets_DEFINITIONS})
set(libs ${Qt5Widgets_LIBRARIES})
set(CMAKE_CXX_FLAGS ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})


IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if(${MSVC}) #Under Visual Studio it is better to specify the path and let pragmas guide the linker
	if(NOT MSVC_VERSION LESS 1600) 
	   add_definitions(-D_HAS_CPP0X=0)
    endif(NOT MSVC_VERSION LESS 1600)
else(${MSVC})
    set(libs ${libs})
endif(${MSVC})

add_executable(eitannealingtest ${eitannealingtest_SRCS})
target_link_libraries(eitannealingtest ${libs})
