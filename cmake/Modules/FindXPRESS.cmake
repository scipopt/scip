find_path(XPRESS_INCLUDE_DIRS
    NAMES xprs.h
    HINTS ${XPRESS_DIR} $ENV{XPRESS_DIR}
    PATH_SUFFIXES include)

# todo: enable recursive search
find_library(XPRESS_LIBRARY
    NAMES xprs
    HINTS ${XPRESS_DIR} $ENV{XPRESS_DIR}
    PATH_SUFFIXES lib)

# todo properly check when pthread is necessary
set(XPRESS_LIBRARIES ${XPRESS_LIBRARY} pthread)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XPRESS DEFAULT_MSG XPRESS_INCLUDE_DIRS XPRESS_LIBRARIES)
