cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(CycleClustering)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(BEFORE PUBLIC ${SCIP_INCLUDE_DIRS})

add_executable(cycleclustering
    src/main.c
    src/cycplugins.c
    src/sepa_partition.c
    src/sepa_edge.c
    src/sepa_subtour.c
    src/heur_cycgreedy.c
    src/heur_fuzzyround.c
    src/heur_cyckerlin.c
    src/heur_redsize.c
    src/branch_multinode.c
    src/event_newsol.c
    src/reader_cyc.c
    src/probdata_cyc.c
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(cycleclustering ${SCIP_LIBRARIES} ${LIBM})
target_compile_options(cycleclustering PRIVATE ${SCIP_COMPILE_FLAGS})

if( TARGET applications )
    add_dependencies( applications cycleclustering)
endif()

add_subdirectory(check)

enable_testing()
