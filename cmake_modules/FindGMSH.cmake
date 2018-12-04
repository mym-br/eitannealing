option(DO_NOT_USE_REDIST_GMSH "Avoid using redistributed gmsh includes")
option(DO_NOT_USE_NATIVE_GMSH "Avoid using local gmsh includes")
if(NOT GMSH_INCLUDE_DIRS)
    if(NOT DO_NOT_USE_NATIVE_GMSH)
        message(STATUS "Looking for gmsh includes")
        find_path(GMSH_INCLUDE_DIRS NAMES onelab.h PATH_SUFFIXES gmsh)
        if (EXISTS ${GMSH_INCLUDE_DIRS})
            message(STATUS "Looking for gmsh includes - found at ${GMSH_INCLUDE_DIRS}")
        else()
            if(NOT DO_NOT_USE_REDIST_GMSH)
                message(STATUS "Looking for gmsh includes - not found, falling back to redist")
                find_path(GMSH_INCLUDE_DIRS NAMES onelab.h PATHS ${CMAKE_SOURCE_DIR}/thirdparty PATH_SUFFIXES gmsh)
                if (EXISTS ${GMSH_INCLUDE_DIRS})
                    message(STATUS "Looking for gmsh includes - using redist at " ${GMSH_INCLUDE_DIRS})
                else()
                    message(STATUS "Looking for gmsh includes - redist missing")
                    message(STATUS "Looking for gmsh includes - not found")
                endif()
            else()
                message(STATUS "Looking for gmsh includes - not found")
            endif()
        endif()
    else()
        if(NOT DO_NOT_USE_REDIST_GMSH)
            find_path(GMSH_INCLUDE_DIRS NAMES onelab.h PATHS ${CMAKE_SOURCE_DIR}/thirdparty PATH_SUFFIXES gmsh)
            if (EXISTS ${GMSH_INCLUDE_DIRS})
                message(STATUS "Looking for gmsh includes - using redist at ${GMSH_INCLUDE_DIRS}")
            else()
                message(STATUS "Looking for gmsh includes - redist missing")
                message(STATUS "Looking for gmsh includes - not found")
            endif()
        endif()
    endif()
endif()
find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH_INCLUDE_DIRS)
# win32 requires wsock32
if(WIN32)
    set(GMSH_LIBRARIES wsock32)
endif()