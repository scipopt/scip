find_path(BLISS_INCLUDE_DIRS
    NAMES graph.hh
    HINTS ${BLISS_DIR}
    PATH_SUFFIXES hh)

find_library(BLISS_LIBRARY
    NAMES bliss
    HINTS ${BLISS_DIR})

set(BLISS_LIBRARIES ${BLISS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLISS DEFAULT_MSG BLISS_INCLUDE_DIRS BLISS_LIBRARIES)
