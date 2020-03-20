find_package(SWIG)
find_package(Python3 COMPONENTS NumPy)

include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
include(${SWIG_USE_FILE})

set(CMAKE_SWIG_FLAGS "")

if(SWIG_FOUND AND Python3_FOUND AND Python3_NumPy_FOUND)
    set_source_files_properties(python/eitdirectsolver.i PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(python/eitdirectsolver.i PROPERTIES SWIG_FLAGS "-includeall")
    SWIG_ADD_LIBRARY(eitdirectsolver LANGUAGE python SOURCES python/eitdirectsolver.i python/eitdirectsolver.cpp $<TARGET_OBJECTS:eit_core>)
    SWIG_LINK_LIBRARIES(eitdirectsolver ${Python3_LIBRARIES} ${libs} ${ADDITIONAL_LIBRARIES})
endif()
