SET(SOPLEX_ROOT_DIR ${PROJECT_SOURCE_DIR}/../../../soplex-3.0.0 CACHE PATH "SoPlex root directory")

FIND_LIBRARY(SOPLEX_LIB soplex
  HINTS ${SOPLEX_ROOT_DIR}/lib ${PROJECT_SOURCE_DIR}/../../../soplex/lib
)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(soplex DEFAULT_MSG SOPLEX_LIB)
