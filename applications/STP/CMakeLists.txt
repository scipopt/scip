cmake_minimum_required(VERSION 3.0)

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
    src/reader_gr.c
    src/reduce_alt.c
    src/reduce_bnd.c
    src/reduce.c
    src/reduce_simple.c
    src/validate.c
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(scipstp ${SCIP_LIBRARIES} ${LIBM})

set_property(TARGET scipstp PROPERTY C_STANDARD 99)

if( TARGET applications )
    add_dependencies( applications scipstp )
endif()

add_subdirectory(check)

enable_testing()
