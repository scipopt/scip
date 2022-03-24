find_path(GLOP_INCLUDE_DIRS
    NAMES ortools/glop/lp_solver.h
    HINTS ${GLOP_DIR} $ENV{GLOP_DIR}
    PATH_SUFFIXES include)

set(GLOP_INCLUDE_DIRS
    ${GLOP_INCLUDE_DIRS}
    ${GLOP_INCLUDE_DIRS}/ortools/gen
    ${GLOP_INCLUDE_DIRS}/dependencies/install/include
    ${GLOP_INCLUDE_DIRS}/include
    )

# todo: enable recursive search
find_library(GLOP_LIBRARY
    NAMES ortools
    HINTS ${GLOP_DIR} $ENV{GLOP_DIR}
    PATH_SUFFIXES lib)

find_library(GLOG_LIBRARY
    NAMES glog
    HINTS ${GLOP_DIR}/dependencies/install/lib $ENV{GLOP_DIR}/dependencies/install/lib
    HINTS ${GLOP_DIR} $ENV{GLOP_DIR}
    PATH_SUFFIXES lib)

if(GLOG_LIBRARY)
    set(GLOP_LIBRARIES ${GLOP_LIBRARY} ${GLOG_LIBRARY})
else()
    set(GLOP_LIBRARIES ${GLOP_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLOP DEFAULT_MSG GLOP_INCLUDE_DIRS GLOP_LIBRARIES)
