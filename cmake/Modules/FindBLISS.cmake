find_path(BLISS_INCLUDE_DIR
    NAMES bliss/graph.hh
    HINTS ${BLISS_DIR}
    PATH_SUFFIXES include)

find_library(BLISS_LIBRARY
    NAMES bliss
    HINTS ${BLISS_DIR}
    PATH_SUFFIXES lib)

set(BLISS_LIBRARIES ${BLISS_LIBRARY})
set(BLISS_INCLUDE_DIRS ${BLISS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLISS DEFAULT_MSG BLISS_INCLUDE_DIRS BLISS_LIBRARIES)
