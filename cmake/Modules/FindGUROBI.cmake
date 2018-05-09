find_path(GUROBI_INCLUDE_DIRS
    NAMES gurobi_c.h
    HINTS ${GUROBI_DIR} $ENV{GUROBI_DIR}
    PATH_SUFFIXES linux64/include include)

# todo: enable recursive search
find_library(GUROBI_LIBRARY
    NAMES gurobi gurobi70 gurobi75 gurobi80
    HINTS ${GUROBI_DIR} $ENV{GUROBI_DIR}
    PATH_SUFFIXES linux64/lib lib)

set(GUROBI_LIBRARIES ${GUROBI_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_INCLUDE_DIRS GUROBI_LIBRARIES)
