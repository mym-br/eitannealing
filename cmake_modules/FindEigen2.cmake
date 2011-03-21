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
    # Build a list of possible locations
    GET_FILENAME_COMPONENT(W32PROGRAMFILESDIR  "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion;ProgramFilesDir]" ABSOLUTE CACHE)

    GET_FILENAME_COMPONENT(W32COMMONFILESDIR  "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion;CommonFilesDir]" ABSOLUTE CACHE)

    find_path(EIGEN2_found_INCLUDES "Eigen\\Eigen" HINTS W32PROGRAMFILESDIR W32COMMONFILESDIR $ENV{ProgramData} PATH_SUFFIXES "Eigen2" "Eigen")
    if(${EIGEN2_found_INCLUDES} STREQUAL "EIGEN2_found_INCLUDES-NOTFOUND")
      set(EIGEN2_FOUND FALSE)
    else(${EIGEN2_found_INCLUDES} STREQUAL "EIGEN2_found_INCLUDES-NOTFOUND")
      set(EIGEN2_INCLUDES ${EIGEN2_found_INCLUDES})
      set(EIGEN2_FOUND TRUE)
    endif(${EIGEN2_found_INCLUDES} STREQUAL "EIGEN2_found_INCLUDES-NOTFOUND")

  endif(NOT WIN32)

  include(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen2 DEFAULT_MSG EIGEN2_INCLUDES )

  mark_as_advanced(EIGEN2_INCLUDES)

endif(EIGEN2_INCLUDES)

