find_path(PDLP_INCLUDE_DIRS
    NAMES ortools/glop/lp_solver.h
    HINTS ${PDLP_DIR} $ENV{PDLP_DIR}
    PATH_SUFFIXES include)

set(PDLP_INCLUDE_DIRS
    ${PDLP_INCLUDE_DIRS}
    ${PDLP_INCLUDE_DIRS}/ortools/gen
    ${PDLP_INCLUDE_DIRS}/dependencies/install/include
    ${PDLP_INCLUDE_DIRS}/include
    ${PDLP_INCLUDE_DIRS}/eigen3
    )

# todo: enable recursive search
find_library(PDLP_LIBRARY
    NAMES ortools
    HINTS ${PDLP_DIR} $ENV{PDLP_DIR}
    PATH_SUFFIXES lib)

find_library(GLOG_LIBRARY
    NAMES glog
    HINTS ${PDLP_DIR}/dependencies/install/lib $ENV{PDLP_DIR}/dependencies/install/lib
    HINTS ${PDLP_DIR} $ENV{PDLP_DIR}
    PATH_SUFFIXES lib)

if(GLOG_LIBRARY)
    set(PDLP_LIBRARIES ${PDLP_LIBRARY} ${GLOG_LIBRARY})
else()
    set(PDLP_LIBRARIES ${PDLP_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PDLP DEFAULT_MSG PDLP_INCLUDE_DIRS PDLP_LIBRARIES)

