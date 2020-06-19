cmake_minimum_required(VERSION 3.3)

project(scipstp)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scipstp
    src/cmain.c
    src/branch_stp.c
    src/completegraph.c
    src/cons_stp.c
    src/dialog_stp.c
    src/dualascent.c
    src/enumeration.c
    src/event_bestsol.c
    src/extreduce_base.c
    src/extreduce_core.c
    src/extreduce_data.c
    src/extreduce_dbg.c
    src/extreduce_extmst.c
    src/extreduce_redcosts.c
    src/extreduce_util.c
    src/graph_history.c
    src/graph_base.c
    src/graph_load.c
    src/graph_mcut.c
    src/graph_path.c
    src/graph_sdpath.c
    src/graph_pcbase.c
    src/graph_save.c
    src/graph_sol.c
    src/graph_trans.c
    src/graph_util.c
    src/graph_vnoi.c
    src/heur_ascendprune.c
    src/heur_local.c
    src/heur_prune.c
    src/heur_rec.c
    src/heur_slackprune.c
    src/heur_tm.c
    src/misc_stp.c
    src/mst.c
    src/pricer_stp.c
    src/probdata_stp.c
    src/prop_stp.c
    src/reader_stp.c
    src/reader_gr.c
    src/reduce_alt.c
    src/reduce_base.c
    src/reduce_bnd.c
    src/reduce_da.c
    src/reduce_ext.c
    src/reduce_pcsimple.c
    src/reduce_sd.c
    src/reduce_sdcomp.c
    src/reduce_sdutil.c
    src/reduce_simple.c
    src/reduce_util.c
    src/shortestpath.c
    src/stptest_base.c
    src/stptest_da.c
    src/stptest_extreduce.c
    src/stptest_extutils.c
    src/stptest_graph.c
    src/stptest_graphpath.c
    src/stptest_heurlocal.c
    src/stptest_heurtm.c
    src/stptest_misc.c
    src/stptest_pcreduce.c
    src/stptest_reducesd.c
    src/stptest_reduceutils.c
    src/solpool.c
    src/solstp.c
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
