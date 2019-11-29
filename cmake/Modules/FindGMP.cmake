# ${GMP_INCLUDE_DIRS} contains the paths to gmp.h (and gmpxx.h) if GMP is found.
# ${GMP_LIBRARIES} contains libgmp and libgmpxx if GMP is found.

# Check whether environment variable GMP_DIR was set.
if(NOT GMP_DIR)
  set(ENV_GMP_DIR $ENV{GMP_DIR})
  if(ENV_GMP_DIR)
    set(GMP_DIR $ENV{GMP_DIR} CACHE PATH "Path to gmp directory")
  endif()
endif()

find_path(GMP_INCLUDE_DIRS
    NAMES gmp.h gmpxx.h
    HINTS ${GMP_DIR}
    PATH_SUFFIXES include)

if(STATIC_GMP)
    find_library(GMP_LIBRARY
        NAMES libgmp.a gmp
        HINTS ${GMP_DIR}
        PATH_SUFFIXES lib)

    find_library(GMPXX_LIBRARY
        NAMES libgmpxx.a gmpxx
        HINTS ${GMP_DIR}
        PATH_SUFFIXES lib)
else()
    find_library(GMP_LIBRARY
        NAMES gmp
        HINTS ${GMP_DIR}
        PATH_SUFFIXES lib)

    find_library(GMPXX_LIBRARY
        NAMES gmpxx
        HINTS ${GMP_DIR}
        PATH_SUFFIXES lib)
endif()

SET(GMP_LIBRARIES ${GMP_LIBRARY} ${GMPXX_LIBRARY})

# look for mpir library and include files when gmp could not be found
if(NOT GMP_LIBRARIES)
    find_path(GMP_INCLUDE_DIRS
       NAMES mpir.h
       HINTS ${GMP_DIR}
       PATH_SUFFIXES include)

   find_library(GMP_LIBRARY
      NAMES mpir
      HINTS ${GMP_DIR}
      PATH_SUFFIXES lib)

   SET(GMP_LIBRARIES ${GMP_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIRS GMP_LIBRARIES)
