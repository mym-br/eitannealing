find_package(SWIG)
find_package(Python3 COMPONENTS NumPy)

include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
include(${SWIG_USE_FILE})

set(CMAKE_SWIG_FLAGS "")

add_library(eitannealinglib STATIC directsolver.cpp $<TARGET_OBJECTS:eit_core>)
set_target_properties(eitannealinglib PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1 DEBUG_POSTFIX _d)
install(FILES directsolver.h ../src/observations.h ../src/basematrix.h DESTINATION EitAnnealing/include)
install(TARGETS eitannealinglib EXPORT eitannealinglib-targets DESTINATION EitAnnealing/lib)
install(EXPORT eitannealinglib-targets DESTINATION EitAnnealing/lib)
install(FILES EitAnnealing-config.cmake DESTINATION EitAnnealing)

# FIXME: compile python debug library
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
else()
    if(CMAKE_BUILD_TYPE EQUAL "RELEASE" AND SWIG_FOUND AND Python3_FOUND AND Python3_NumPy_FOUND)
        set_source_files_properties(python/eitdirectsolver.i PROPERTIES CPLUSPLUS ON)
        set_source_files_properties(python/eitdirectsolver.i PROPERTIES SWIG_FLAGS "-includeall")
        SWIG_ADD_LIBRARY(eitdirectsolver LANGUAGE python SOURCES python/eitdirectsolver.i python/eitdirectsolver.cpp $<TARGET_OBJECTS:eit_core>)
        SWIG_LINK_LIBRARIES(eitdirectsolver ${Python3_LIBRARIES} ${libs} ${ADDITIONAL_LIBRARIES})
    endif()
endif()