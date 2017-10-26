cmake_minimum_required(VERSION 3.0)

project(CycleClustering)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(cycleclustering
    src/main.c
    src/spaplugins.c
    src/sepa_partition.c
    src/sepa_edge.c
    src/sepa_subtour.c
    src/heur_spagreedy.c
    src/heur_fuzzyround.c
    src/heur_spakerlin.c
    src/branch_multinode.c
    src/reader_spa.c
    src/probdata_spa.c
)

target_link_libraries(cycleclustering ${SCIP_LIBRARIES} m)

if( TARGET applications )
    add_dependencies( applications ccscip )
endif()

add_subdirectory(check)

enable_testing()
