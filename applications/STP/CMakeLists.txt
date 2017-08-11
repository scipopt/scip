project(scipstp)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scipstp 
    src/cmain.c
    src/branch_stp.c
    src/cons_stp.c
    src/dialog_stp.c
    src/event_bestsol.c
    src/grphbase.c
    src/grphload.c
    src/grphmcut.c
    src/grphpath.c
    src/grphsave.c
    src/heur_ascendprune.c
    src/heur_local.c
    src/heur_prune.c
    src/heur_rec.c
    src/heur_slackprune.c
    src/heur_tm.c
    src/misc_stp.c
    src/pricer_stp.c
    src/probdata_stp.c
    src/prop_stp.c
    src/reader_stp.c
    src/reduce_alt.c
    src/reduce_bnd.c
    src/reduce.c
    src/reduce_simple.c
    src/validate.c
)

target_link_libraries(scipstp ${SCIP_LIBRARIES} m)
                                                 
add_dependencies(applications scipstp)

add_subdirectory(check)

enable_testing()
