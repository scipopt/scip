SET(SCIP_ROOT_DIR ${PROJECT_SOURCE_DIR}/../.. CACHE PATH "SCIP root directory")

find_path(SCIP_INC scip/scip.h
  HINTS ${SCIP_ROOT_DIR}/src
)

find_library(OBJSCIP_LIB objscip
  HINTS ${SCIP_ROOT_DIR}/lib/static
)

find_library(SCIP_LIB scip
  HINTS ${SCIP_ROOT_DIR}/lib/static
)

find_library(NLPI_LIB NAMES nlpi.cppad nlpi
  HINTS ${SCIP_ROOT_DIR}/lib/static
)

find_library(LPS_LIB NAMES lpispx2
  HINTS ${SCIP_ROOT_DIR}/lib/static
)

find_library(TPINONE_LIB NAMES tpinone
        HINTS ${SCIP_ROOT_DIR}/lib/static
        )

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(scip REQUIRED SCIP_INC OBJSCIP_LIB SCIP_LIB NLPI_LIB LPS_LIB TPINONE_LIB)
