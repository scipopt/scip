cmake_minimum_required(VERSION 3.0)

project(CycleClustering)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

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

if( TARGET applications )
    add_dependencies( applications cycleclustering)
endif()

add_subdirectory(check)

enable_testing()
