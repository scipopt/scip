SET(SOPLEX_ROOT_DIR ${PROJECT_SOURCE_DIR}/../../../soplex-2.2.1 CACHE PATH "SoPlex root directory")

FIND_LIBRARY(SOPLEX_LIB soplex
  HINTS ${SOPLEX_ROOT_DIR}/lib ${PROJECT_SOURCE_DIR}/../../../soplex/lib
)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(soplex DEFAULT_MSG SOPLEX_LIB)
