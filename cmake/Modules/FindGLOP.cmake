find_path(GLOP_INCLUDE_DIRS
    NAMES ortools/glop/lp_solver.h
    HINTS ${GLOP_DIR} $ENV{GLOP_DIR}
    PATH_SUFFIXES include)

set(GLOP_INCLUDE_DIRS
    ${GLOP_INCLUDE_DIRS}
    ${GLOP_INCLUDE_DIRS}/ortools/gen
    ${GLOP_INCLUDE_DIRS}/dependencies/install/include
    )

# todo: enable recursive search
find_library(GLOP_LIBRARY
    NAMES ortools
    HINTS ${GLOP_DIR} $ENV{GLOP_DIR}
    PATH_SUFFIXES lib)

set(GLOP_LIBRARIES ${GLOP_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLOP DEFAULT_MSG GLOP_INCLUDE_DIRS GLOP_LIBRARIES)
