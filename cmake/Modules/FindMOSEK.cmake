find_path(MOSEK_INCLUDE_DIRS
    NAMES mosek.h
    HINTS ${MOSEK_DIR} $ENV{MOSEK_DIR}
    PATH_SUFFIXES h)

find_library(MOSEK_LIBRARY
    NAMES mosek64
    HINTS ${MOSEK_DIR} $ENV{MOSEK_DIR}
    PATH_SUFFIXES bin)

set(MOSEK_LIBRARIES ${MOSEK_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK DEFAULT_MSG MOSEK_INCLUDE_DIRS MOSEK_LIBRARIES)
