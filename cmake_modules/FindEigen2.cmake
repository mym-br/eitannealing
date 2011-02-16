# - Try to find Eigen2
# Input variables
#
# Once done this will define
#
#  EIGEN2_FOUND - system has libqalculate
#  EIGEN2_INCLUDES - eigen include path
#

if(EIGEN2_INCLUDES)

  # in cache already
  set(EIGEN2_FOUND)

else(EIGEN2_INCLUDES)
  if(NOT WIN32)
    include(UsePkgConfig)

    exec_program(${PKGCONFIG_EXECUTABLE} ARGS eigen2 --exists RETURN_VALUE _return_VALUE OUTPUT_VARIABLE _pkgconfigDevNull)
    
    if(_return_VALUE STREQUAL "0")
      exec_program(${PKGCONFIG_EXECUTABLE} ARGS eigen2 --cflags-only-I OUTPUT_VARIABLE EIGEN2_INCLUDES)    
      string(REPLACE "-I" "" EIGEN2_INCLUDES ${EIGEN2_INCLUDES})
      set(EIGEN2_FOUND TRUE)
    endif(_return_VALUE STREQUAL "0")

  else(NOT WIN32)
    # XXX: currently no libqalculate on windows
    set(EIGEN2_FOUND FALSE)

  endif(NOT WIN32)

  include(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen2 DEFAULT_MSG EIGEN2_INCLUDES )

  mark_as_advanced(EIGEN2_INCLUDES)

endif(EIGEN2_INCLUDES)

