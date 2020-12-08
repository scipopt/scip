# ${MPFR_INCLUDE_DIRS} contains the paths to MPFR.h 
# ${MPFR_LIBRARIES} contains libMPFR.

# Check whether environment variable MPFR_DIR was set.
if(NOT MPFR_DIR)
  set(ENV_MPFR_DIR $ENV{MPFR_DIR})
  if(ENV_MPFR_DIR)
    set(MPFR_DIR $ENV{MPFR_DIR} CACHE PATH "Path to MPFR directory")
  endif()
endif()

find_path(MPFR_INCLUDE_DIRS
    NAMES mpfr.h
    HINTS ${MPFR_DIR}
    PATH_SUFFIXES include)

if(STATIC_GMP)
#  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    find_library(MPFR_LIBRARY
        NAMES libmpfr.a MPFR
        HINTS ${MPFR_DIR}
        NO_DEFAULT_PATH
        PATH_SUFFIXES lib)
else()
    find_library(MPFR_LIBRARY
        NAMES libmpfr.so MPFR
        HINTS ${MPFR_DIR}
        PATH_SUFFIXES lib)
endif()

SET(MPFR_LIBRARIES ${MPFR_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDE_DIRS MPFR_LIBRARIES)
