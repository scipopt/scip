cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Coloring)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(coloring
    src/branch_coloring.c
    src/branch_strongcoloring.c
    src/coloringplugins.c
    src/cons_storeGraph.c
    src/heur_init.c
    src/main.c
    src/pricer_coloring.c
    src/probdata_coloring.c
    src/reader_col.c
    src/reader_csol.c
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(coloring ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(coloring PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET applications )
    add_dependencies( applications coloring )
endif()

add_subdirectory(check)

enable_testing()
