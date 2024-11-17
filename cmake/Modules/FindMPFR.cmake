# ${MPFR_INCLUDE_DIRS} contains the paths to mpfr.h (and mpfrxx.h) if MPFR is found.
# ${MPFR_LIBRARIES} contains libmpfr and libmpfrxx if MPFR is found.

# Check whether environment variable MPFR_DIR was set.
if(NOT MPFR_DIR)
  set(ENV_MPFR_DIR $ENV{MPFR_DIR})
  if(ENV_MPFR_DIR)
    set(MPFR_DIR $ENV{MPFR_DIR} CACHE PATH "Path to mpfr directory")
  endif()
endif()

find_path(MPFR_INCLUDE_DIRS
    NAMES mpfr.h
    HINTS ${MPFR_DIR}
    PATH_SUFFIXES include)

if(STATIC_MPFR)
    find_library(MPFR_LIBRARY
        NAMES libmpfr.a mpfr
        HINTS ${MPFR_DIR}
        PATH_SUFFIXES lib)
else()
    find_library(MPFR_LIBRARY
        NAMES libmpfr.so mpfr
        HINTS ${MPFR_DIR}
        PATH_SUFFIXES lib)
endif()

SET(MPFR_LIBRARIES ${MPFR_LIBRARY} ${MPFRXX_LIBRARY})

# look for mpir library and include files when mpfr could not be found
if(NOT MPFR_LIBRARIES)
    find_path(MPFR_INCLUDE_DIRS
       NAMES mpir.h
       HINTS ${MPFR_DIR}
       PATH_SUFFIXES include)

   find_library(MPFR_LIBRARY
      NAMES mpir
      HINTS ${MPFR_DIR}
      PATH_SUFFIXES lib)

   SET(MPFR_LIBRARIES ${MPFR_LIBRARY})
endif()

file(GLOB MPFR_HEADERS "${MPFR_INCLUDE_DIRS}/mpfr.h" "${MPFR_INCLUDE_DIRS}/mpfr-*.h")
foreach (mpfr_header_filename ${MPFR_HEADERS})
    file(READ "${mpfr_header_filename}" _mpfr_version_header)
    string(REGEX MATCH
            "define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" _mpfr_major_version_match
            "${_mpfr_version_header}")
    if (_mpfr_major_version_match)
        set(MPFR_MAJOR_VERSION "${CMAKE_MATCH_1}")
        string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)"
                _mpfr_minor_version_match "${_mpfr_version_header}")
        set(MPFR_MINOR_VERSION "${CMAKE_MATCH_1}")
        string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
                _mpfr_patchlevel_version_match "${_mpfr_version_header}")
        set(MPFR_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")
        set(MPFR_VERSION
                ${MPFR_MAJOR_VERSION}.${MPFR_MINOR_VERSION}.${MPFR_PATCHLEVEL_VERSION})
    endif ()
endforeach ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDE_DIRS MPFR_LIBRARIES)
