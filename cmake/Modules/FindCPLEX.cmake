find_path(CPLEX_INCLUDE_DIRS
    NAMES cplex.h
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES include/ilcplex include)

# todo: enable recursive search
find_library(CPLEX_LIBRARY
    NAMES cplex
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES lib/x86-64_linux/static_pic lib)

# todo properly check when pthread is necessary
set(CPLEX_LIBRARIES ${CPLEX_LIBRARY} pthread)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_INCLUDE_DIRS CPLEX_LIBRARIES)
