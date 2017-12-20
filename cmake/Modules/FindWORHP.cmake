find_path(WORHP_INCLUDE_DIRS
    NAMES worhp/worhp.h
    HINTS ${WORHP_DIR} $ENV{WORHP_DIR}
    PATH_SUFFIXES include)

# todo: enable recursive search
find_library(WORHP_LIBRARY
    NAMES worhp
    HINTS ${WORHP_DIR} $ENV{WORHP_DIR}
    PATH_SUFFIXES lib)

set(WORHP_LIBRARIES ${WORHP_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(WORHP DEFAULT_MSG WORHP_INCLUDE_DIRS WORHP_LIBRARIES)
