project(Coloring)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

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

target_link_libraries(coloring ${SCIP_LIBRARIES} m)
                                                 
add_dependencies(applications coloring)

add_subdirectory(check)

enable_testing()
