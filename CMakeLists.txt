 project(eitannealingtests)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Qt4 REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(Eigen2 REQUIRED)
find_package(CXX11)

IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "-std=c++11")
endif()

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR} ${Boost_INCLUDE_DIR} ${EIGEN2_INCLUDES})

set(eitannealingtest_SRCS src/sparseincompletelq.cpp
  src/graphics.cpp
  src/initproblem.cpp
  src/main.cpp
  src/nodecoeficients.cpp
  src/solver.cpp
  src/incomplete_cholesky.cpp
  src/init_obs_problem.cpp
  src/solution.cpp
  src/mt19937-64/mt19937-64.c
  src/observations.cpp
  src/solver_lb.cpp
  src/solution_lb.cpp
)
set(libs ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})

qt4_automoc(${eitannealingtest_SRCS})

if(${MSVC}) #Under Visual Studio it is better to specify the path and let pragmas guide the linker
	if(NOT MSVC_VERSION LESS 1600) 
	   add_definitions(-D_HAS_CPP0X=0)
    endif(NOT MSVC_VERSION LESS 1600)
	link_directories(${Boost_LIBRARY_DIRS})
else(${MSVC})
    set(libs ${libs} ${Boost_LIBRARIES})
endif(${MSVC})

add_executable(eitannealingtest ${eitannealingtest_SRCS})
target_link_libraries(eitannealingtest ${libs})
