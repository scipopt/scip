cmake_minimum_required(VERSION 3.0)

project(ccscip)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../lib)

set (SOURCE_FILES
        src/main.c
        src/spaplugins.c
        src/sepa_partition.c
        src/sepa_edge.c
        src/sepa_subtour.c
        src/heur_spagreedy.c
        src/heur_fuzzyround.c
        src/heur_spakerlin.c
        src/branch_multinode.c
        src/cons_triangle.c
        src/reader_spa.c
        src/probdata_spa.c )

add_executable(ccscip
        ${SOURCE_FILES}
        )

        target_link_libraries(ccscip libscip.so)
