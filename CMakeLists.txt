project(eitannealingtests)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) 
find_package(Qt4 REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(Eigen2 REQUIRED)

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR} ${EIGEN2_INCLUDES})

set(eitannealingtest_SRCS
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
)
qt4_automoc(${eitannealingtest_SRCS})
add_executable(eitannealingtest ${eitannealingtest_SRCS})
target_link_libraries(eitannealingtest ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${Boost_LIBRARIES} cblas)
