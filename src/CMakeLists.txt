include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(scipsources
    scip/branch_allfullstrong.c
    scip/branch_cloud.c
    scip/branch_distribution.c
    scip/branch_fullstrong.c
    scip/branch_inference.c
    scip/branch_leastinf.c
    scip/branch_mostinf.c
    scip/branch_multaggr.c
    scip/branch_nodereopt.c
    scip/branch_pscost.c
    scip/branch_random.c
    scip/branch_relpscost.c
    scip/cons_abspower.c
    scip/compr_largestrepr.c
    scip/compr_weakcompr.c
    scip/concsolver_scip.c
    scip/cons_and.c
    scip/cons_bivariate.c
    scip/cons_bounddisjunction.c
    scip/cons_cardinality.c
    scip/cons_components.c
    scip/cons_conjunction.c
    scip/cons_countsols.c
    scip/cons_cumulative.c
    scip/cons_disjunction.c
    scip/cons_indicator.c
    scip/cons_integral.c
    scip/cons_knapsack.c
    scip/cons_linear.c
    scip/cons_linking.c
    scip/cons_logicor.c
    scip/cons_nonlinear.c
    scip/cons_or.c
    scip/cons_orbitope.c
    scip/cons_pseudoboolean.c
    scip/cons_quadratic.c
    scip/cons_setppc.c
    scip/cons_soc.c
    scip/cons_sos1.c
    scip/cons_sos2.c
    scip/cons_superindicator.c
    scip/cons_varbound.c
    scip/cons_xor.c
    scip/dialog_default.c
    scip/event_globalbnd.c
    scip/event_softtimelimit.c
    scip/event_solvingphase.c
    scip/disp_default.c
    scip/heur_actconsdiving.c
    scip/heur_bound.c
    scip/heur_clique.c
    scip/heur_coefdiving.c
    scip/heur_completesol.c
    scip/heur_crossover.c
    scip/heur_dins.c
    scip/heur_distributiondiving.c
    scip/heur_dualval.c
    scip/heur_feaspump.c
    scip/heur_fixandinfer.c
    scip/heur_fracdiving.c
    scip/heur_gins.c
    scip/heur_guideddiving.c
    scip/heur_indicator.c
    scip/heur_intdiving.c
    scip/heur_intshifting.c
    scip/heur_linesearchdiving.c
    scip/heur_localbranching.c
    scip/heur_lpface.c
    scip/heur_locks.c
    scip/heur_mutation.c
    scip/heur_multistart.c
    scip/heur_mpec.c
    scip/heur_nlpdiving.c
    scip/heur_objpscostdiving.c
    scip/heur_octane.c
    scip/heur_ofins.c
    scip/heur_oneopt.c
    scip/heur_proximity.c
    scip/heur_pscostdiving.c
    scip/heur_reoptsols.c
    scip/heur_repair.c
    scip/heur_randrounding.c
    scip/heur_rens.c
    scip/heur_rins.c
    scip/heur_rootsoldiving.c
    scip/heur_rounding.c
    scip/heur_shiftandpropagate.c
    scip/heur_shifting.c
    scip/heur_simplerounding.c
    scip/heur_subnlp.c
    scip/heur_sync.c
    scip/heur_trivial.c
    scip/heur_trivialnegation.c
    scip/heur_trysol.c
    scip/heur_twoopt.c
    scip/heur_undercover.c
    scip/heur_vbounds.c
    scip/heur_veclendiving.c
    scip/heur_zeroobj.c
    scip/heur_zirounding.c
    scip/message_default.c
    scip/nodesel_bfs.c
    scip/nodesel_breadthfirst.c
    scip/nodesel_dfs.c
    scip/nodesel_estimate.c
    scip/nodesel_hybridestim.c
    scip/nodesel_restartdfs.c
    scip/nodesel_uct.c
    scip/presol_boundshift.c
    scip/presol_convertinttobin.c
    scip/presol_domcol.c
    scip/presol_dualagg.c
    scip/presol_dualcomp.c
    scip/presol_dualinfer.c
    scip/presol_gateextraction.c
    scip/presol_implfree.c
    scip/presol_implics.c
    scip/presol_inttobinary.c
    scip/presol_qpkktref.c
    scip/presol_redvub.c
    scip/presol_trivial.c
    scip/presol_tworowbnd.c
    scip/presol_stuffing.c
    scip/prop_dualfix.c
    scip/prop_genvbounds.c
    scip/prop_nlobbt.c
    scip/prop_obbt.c
    scip/prop_probing.c
    scip/prop_pseudoobj.c
    scip/prop_redcost.c
    scip/prop_rootredcost.c
    scip/prop_vbounds.c
    scip/prop_sync.c
    scip/reader_bnd.c
    scip/reader_ccg.c
    scip/reader_cip.c
    scip/reader_cnf.c
    scip/reader_diff.c
    scip/reader_fix.c
    scip/reader_fzn.c
    scip/reader_gms.c
    scip/reader_lp.c
    scip/reader_mps.c
    scip/reader_mst.c
    scip/reader_opb.c
    scip/reader_osil.c
    scip/reader_pip.c
    scip/reader_pbm.c
    scip/reader_ppm.c
    scip/reader_rlp.c
    scip/reader_sol.c
    scip/reader_wbo.c
    scip/reader_zpl.c
    scip/sepa_cgmip.c
    scip/sepa_clique.c
    scip/sepa_closecuts.c
    scip/sepa_cmir.c
    scip/sepa_convexproj.c
    scip/sepa_disjunctive.c
    scip/sepa_eccuts.c
    scip/sepa_flowcover.c
    scip/sepa_gauge.c
    scip/sepa_gomory.c
    scip/sepa_impliedbounds.c
    scip/sepa_intobj.c
    scip/sepa_mcf.c
    scip/sepa_oddcycle.c
    scip/sepa_rapidlearning.c
    scip/sepa_strongcg.c
    scip/sepa_zerohalf.c
    scip/boundstore.c
    scip/branch.c
    scip/clock.c
    scip/concsolver.c
    scip/concurrent.c
    scip/conflict.c
    scip/conflictstore.c
    scip/cons.c
    scip/cutpool.c
    scip/cuts.c
    scip/debug.c
    scip/dialog.c
    scip/disp.c
    scip/event.c
    scip/fileio.c
    scip/heur.c
    scip/heuristics.c
    scip/compr.c
    scip/history.c
    scip/implics.c
    scip/interrupt.c
    scip/intervalarith.c
    scip/lp.c
    scip/matrix.c
    scip/mem.c
    scip/misc.c
    scip/nlp.c
    scip/nodesel.c
    scip/paramset.c
    scip/presol.c
    scip/presolve.c
    scip/pricestore.c
    scip/pricer.c
    scip/primal.c
    scip/prob.c
    scip/prop.c
    scip/reader.c
    scip/relax.c
    scip/reopt.c
    scip/retcode.c
    scip/scip.c
    scip/scipbuildflags.c
    scip/scipdefplugins.c
    scip/scipgithash.c
    scip/scipshell.c
    scip/sepa.c
    scip/sepastore.c
    scip/set.c
    scip/sol.c
    scip/solve.c
    scip/stat.c
    scip/syncstore.c
    scip/tree.c
    scip/var.c
    scip/visual.c
    tclique/tclique_branch.c
    tclique/tclique_coloring.c
    tclique/tclique_graph.c
    dijkstra/dijkstra.c
    xml/xmlparse.c
)

set(objscipsources
    objscip/objbranchrule.cpp
    objscip/objconshdlr.cpp
    objscip/objdialog.cpp
    objscip/objdisp.cpp
    objscip/objeventhdlr.cpp
    objscip/objheur.cpp
    objscip/objmessagehdlr.cpp
    objscip/objnodesel.cpp
    objscip/objpresol.cpp
    objscip/objpricer.cpp
    objscip/objprobdata.cpp
    objscip/objprop.cpp
    objscip/objreader.cpp
    objscip/objrelax.cpp
    objscip/objsepa.cpp
    objscip/objvardata.cpp
)

set(nlpisources
    nlpi/nlpi.c
    nlpi/nlpi_all.c
    nlpi/nlpioracle.c
    nlpi/expr.c
    nlpi/intervalarithext.cpp
    blockmemshell/memory.c
    scip/rbtree.c
    scip/misc.c
    scip/intervalarith.c
    scip/interrupt.c
    scip/message.c
    ${nlpi}
    ${exprinterpret}
)

set(lpisources ${lpi} scip/bitencode.c blockmemshell/memory.c scip/rbtree.c scip/message.c)

set(blockmemshellheader
    blockmemshell/memory.h
)

set(lpiheaders
    lpi/lpi.h
    lpi/type_lpi.h
)

set(dijkstraheaders
    dijkstra/dijkstra.h
)

set(nlpiheaders
    nlpi/exprinterpret.h
    nlpi/intervalarithext.h
    nlpi/nlpi.h
    nlpi/nlpi_ipopt.h
    nlpi/nlpi_worhp.h
    nlpi/nlpioracle.h
    nlpi/pub_expr.h
    nlpi/struct_expr.h
    nlpi/struct_nlpi.h
    nlpi/type_expr.h
    nlpi/type_exprinterpret.h
    nlpi/type_nlpi.h
)

set(objscipheaders
    objscip/objbranchrule.h
    objscip/objcloneable.h
    objscip/objconshdlr.h
    objscip/objdialog.h
    objscip/objdisp.h
    objscip/objeventhdlr.h
    objscip/objheur.h
    objscip/objmessagehdlr.h
    objscip/objnodesel.h
    objscip/objpresol.h
    objscip/objpricer.h
    objscip/objprobcloneable.h
    objscip/objprobdata.h
    objscip/objprop.h
    objscip/objreader.h
    objscip/objrelax.h
    objscip/objscipdefplugins.h
    objscip/objscip.h
    objscip/objsepa.h
    objscip/objvardata.h
    objscip/type_objcloneable.h
    objscip/type_objprobcloneable.h
)

set(scipheaders
    scip/bitencode.h
    scip/boundstore.h
    scip/branch_allfullstrong.h
    scip/branch_cloud.h
    scip/branch_distribution.h
    scip/branch_fullstrong.h
    scip/branch.h
    scip/branch_inference.h
    scip/branch_leastinf.h
    scip/branch_mostinf.h
    scip/branch_multaggr.h
    scip/branch_nodereopt.h
    scip/branch_pscost.h
    scip/branch_random.h
    scip/branch_relpscost.h
    scip/clock.h
    scip/compr.h
    scip/compr_largestrepr.h
    scip/compr_weakcompr.h
    scip/concsolver.h
    scip/concsolver_scip.h
    scip/concurrent.h
    scip/conflict.h
    scip/conflictstore.h
    scip/cons_abspower.h
    scip/cons_and.h
    scip/cons_bivariate.h
    scip/cons_bounddisjunction.h
    scip/cons_cardinality.h
    scip/cons_components.h
    scip/cons_conjunction.h
    scip/cons_countsols.h
    scip/cons_cumulative.h
    scip/cons_disjunction.h
    scip/cons.h
    scip/cons_indicator.h
    scip/cons_integral.h
    scip/cons_knapsack.h
    scip/cons_linear.h
    scip/cons_linking.h
    scip/cons_logicor.h
    scip/cons_nonlinear.h
    scip/cons_orbitope.h
    scip/cons_or.h
    scip/cons_pseudoboolean.h
    scip/cons_quadratic.h
    scip/cons_setppc.h
    scip/cons_soc.h
    scip/cons_sos1.h
    scip/cons_sos2.h
    scip/cons_superindicator.h
    scip/cons_varbound.h
    scip/cons_xor.h
    scip/cutpool.h
    scip/cuts.h
    scip/debug.h
    scip/def.h
    scip/dialog_default.h
    scip/dialog.h
    scip/disp_default.h
    scip/disp.h
    scip/event_globalbnd.h
    scip/event.h
    scip/event_softtimelimit.h
    scip/event_solvingphase.h
    scip/heur_actconsdiving.h
    scip/heur_bound.h
    scip/heur_clique.h
    scip/heur_coefdiving.h
    scip/heur_completesol.h
    scip/heur_crossover.h
    scip/heur_dins.h
    scip/heur_distributiondiving.h
    scip/heur_dualval.h
    scip/heur_feaspump.h
    scip/heur_fixandinfer.h
    scip/heur_fracdiving.h
    scip/heur_gins.h
    scip/heur_guideddiving.h
    scip/heur.h
    scip/heur_indicator.h
    scip/heur_intdiving.h
    scip/heur_intshifting.h
    scip/heuristics.h
    scip/heur_linesearchdiving.h
    scip/heur_localbranching.h
    scip/heur_locks.h
    scip/heur_lpface.h
    scip/heur_multistart.h
    scip/heur_mutation.h
    scip/heur_mpec.h
    scip/heur_nlpdiving.h
    scip/heur_objpscostdiving.h
    scip/heur_octane.h
    scip/heur_ofins.h
    scip/heur_oneopt.h
    scip/heur_proximity.h
    scip/heur_pscostdiving.h
    scip/heur_randrounding.h
    scip/heur_rens.h
    scip/heur_reoptsols.h
    scip/heur_repair.h
    scip/heur_rins.h
    scip/heur_rootsoldiving.h
    scip/heur_rounding.h
    scip/heur_shiftandpropagate.h
    scip/heur_shifting.h
    scip/heur_simplerounding.h
    scip/heur_subnlp.h
    scip/heur_sync.h
    scip/heur_trivial.h
    scip/heur_trivialnegation.h
    scip/heur_trysol.h
    scip/heur_twoopt.h
    scip/heur_undercover.h
    scip/heur_vbounds.h
    scip/heur_veclendiving.h
    scip/heur_zeroobj.h
    scip/heur_zirounding.h
    scip/history.h
    scip/implics.h
    scip/interrupt.h
    scip/intervalarith.h
    scip/lp.h
    scip/mem.h
    scip/message_default.h
    scip/message.h
    scip/misc.h
    scip/nlp.h
    scip/nodesel_bfs.h
    scip/nodesel_breadthfirst.h
    scip/nodesel_dfs.h
    scip/nodesel_estimate.h
    scip/nodesel.h
    scip/nodesel_hybridestim.h
    scip/nodesel_restartdfs.h
    scip/nodesel_uct.h
    scip/paramset.h
    scip/presol_boundshift.h
    scip/presol_convertinttobin.h
    scip/presol_domcol.h
    scip/presol_dualagg.h
    scip/presol_dualcomp.h
    scip/presol_dualinfer.h
    scip/presol_gateextraction.h
    scip/presol.h
    scip/presol_implfree.h
    scip/presol_implics.h
    scip/presol_inttobinary.h
    scip/presol_qpkktref.h
    scip/presol_redvub.h
    scip/presol_stuffing.h
    scip/presol_trivial.h
    scip/presol_tworowbnd.h
    scip/presolve.h
    scip/pricer.h
    scip/pricestore.h
    scip/primal.h
    scip/prob.h
    scip/prop_dualfix.h
    scip/prop_genvbounds.h
    scip/prop.h
    scip/prop_nlobbt.h
    scip/prop_obbt.h
    scip/prop_probing.h
    scip/prop_pseudoobj.h
    scip/prop_redcost.h
    scip/prop_rootredcost.h
    scip/prop_sync.h
    scip/prop_vbounds.h
    scip/pub_branch.h
    scip/pub_compr.h
    scip/pub_conflict.h
    scip/pub_cons.h
    scip/pub_cutpool.h
    scip/pub_dialog.h
    scip/pub_disp.h
    scip/pub_event.h
    scip/pub_fileio.h
    scip/pub_heur.h
    scip/pub_history.h
    scip/pub_implics.h
    scip/pub_lp.h
    scip/pub_matrix.h
    scip/pub_message.h
    scip/pub_misc.h
    scip/pub_misc_select.h
    scip/pub_misc_sort.h
    scip/pub_nlp.h
    scip/pub_nodesel.h
    scip/pub_paramset.h
    scip/pub_presol.h
    scip/pub_pricer.h
    scip/pub_prop.h
    scip/pub_reader.h
    scip/pub_relax.h
    scip/pub_reopt.h
    scip/pub_sepa.h
    scip/pub_sol.h
    scip/pub_tree.h
    scip/pub_var.h
    scip/rbtree.h
    scip/reader_bnd.h
    scip/reader_ccg.h
    scip/reader_cip.h
    scip/reader_cnf.h
    scip/reader_diff.h
    scip/reader_fix.h
    scip/reader_fzn.h
    scip/reader_gms.h
    scip/reader.h
    scip/reader_lp.h
    scip/reader_mps.h
    scip/reader_mst.h
    scip/reader_opb.h
    scip/reader_osil.h
    scip/reader_pbm.h
    scip/reader_pip.h
    scip/reader_ppm.h
    scip/reader_rlp.h
    scip/reader_sol.h
    scip/reader_wbo.h
    scip/reader_zpl.h
    scip/relax.h
    scip/reopt.h
    scip/retcode.h
    scip/scipbuildflags.h
    scip/scipdefplugins.h
    scip/scipgithash.h
    scip/scip.h
    scip/scipshell.h
    scip/sepa_cgmip.h
    scip/sepa_clique.h
    scip/sepa_closecuts.h
    scip/sepa_cmir.h
    scip/sepa_convexproj.h
    scip/sepa_disjunctive.h
    scip/sepa_eccuts.h
    scip/sepa_flowcover.h
    scip/sepa_gauge.h
    scip/sepa_gomory.h
    scip/sepa.h
    scip/sepa_impliedbounds.h
    scip/sepa_intobj.h
    scip/sepa_mcf.h
    scip/sepa_oddcycle.h
    scip/sepa_rapidlearning.h
    scip/sepastore.h
    scip/sepa_strongcg.h
    scip/sepa_zerohalf.h
    scip/set.h
    scip/sol.h
    scip/solve.h
    scip/stat.h
    scip/struct_branch.h
    scip/struct_clock.h
    scip/struct_compr.h
    scip/struct_concsolver.h
    scip/struct_concurrent.h
    scip/struct_conflict.h
    scip/struct_conflictstore.h
    scip/struct_cons.h
    scip/struct_cutpool.h
    scip/struct_dialog.h
    scip/struct_disp.h
    scip/struct_event.h
    scip/struct_heur.h
    scip/struct_history.h
    scip/struct_implics.h
    scip/struct_lp.h
    scip/struct_matrix.h
    scip/struct_mem.h
    scip/struct_message.h
    scip/struct_misc.h
    scip/struct_nlp.h
    scip/struct_nodesel.h
    scip/struct_paramset.h
    scip/struct_presol.h
    scip/struct_pricer.h
    scip/struct_pricestore.h
    scip/struct_primal.h
    scip/struct_prob.h
    scip/struct_prop.h
    scip/struct_reader.h
    scip/struct_relax.h
    scip/struct_reopt.h
    scip/struct_scip.h
    scip/struct_sepa.h
    scip/struct_sepastore.h
    scip/struct_set.h
    scip/struct_sol.h
    scip/struct_stat.h
    scip/struct_syncstore.h
    scip/struct_tree.h
    scip/struct_var.h
    scip/struct_visual.h
    scip/syncstore.h
    scip/tree.h
    scip/type_branch.h
    scip/type_clock.h
    scip/type_compr.h
    scip/type_concsolver.h
    scip/type_concurrent.h
    scip/type_conflict.h
    scip/type_conflictstore.h
    scip/type_cons.h
    scip/type_cutpool.h
    scip/type_dialog.h
    scip/type_disp.h
    scip/type_event.h
    scip/type_heur.h
    scip/type_history.h
    scip/type_implics.h
    scip/type_interrupt.h
    scip/type_lp.h
    scip/type_matrix.h
    scip/type_mem.h
    scip/type_message.h
    scip/type_misc.h
    scip/type_nlp.h
    scip/type_nodesel.h
    scip/type_paramset.h
    scip/type_presol.h
    scip/type_pricer.h
    scip/type_pricestore.h
    scip/type_primal.h
    scip/type_prob.h
    scip/type_prop.h
    scip/type_reader.h
    scip/type_relax.h
    scip/type_reopt.h
    scip/type_result.h
    scip/type_retcode.h
    scip/type_scip.h
    scip/type_sepa.h
    scip/type_sepastore.h
    scip/type_set.h
    scip/type_sol.h
    scip/type_stat.h
    scip/type_syncstore.h
    scip/type_timing.h
    scip/type_tree.h
    scip/type_var.h
    scip/type_visual.h
    scip/var.h
    scip/visual.h
)

set(tcliqueheaders
    tclique/tclique_coloring.h
    tclique/tclique_def.h
    tclique/tclique.h
)

set(tinycthreadheader
    tinycthread/tinycthread.h
   )

set(tpiheaders
    tpi/def_openmp.h
    tpi/tpi.h
    tpi/tpi_none.h
    tpi/tpi_openmp.h
    tpi/tpi_tnycthrd.h
    tpi/type_tpi.h
    tpi/type_tpi_none.h
    tpi/type_tpi_openmp.h
    tpi/type_tpi_tnycthrd.h
)

set(xmlheaders
    xml/xmldef.h
    xml/xml.h
)


add_library(scip OBJECT ${scipsources})
add_dependencies(scip update_githash update_buildflags)
add_library(objscip STATIC ${objscipsources})
add_library(nlpi STATIC ${nlpisources})
add_library(lpi STATIC ${lpisources})
add_library(tpi STATIC ${tpisources})
target_link_libraries(lpi ${GMP_LIBRARIES} ${LPS_LIBRARIES})
target_link_libraries(nlpi ${NLPI_LIBRARIES})

if(SHARED)
    #in the shared case the scipopt needs to compile the sources relocatable again
    add_library(scipopt ${scipsources} ${objscipsources} ${nlpisources} ${lpisources} ${tpisources})
    add_dependencies(scipopt update_githash update_buildflags)
else()
    #in the static case scipopt can reuse the non relocatable objects
    add_library(scipopt $<TARGET_OBJECTS:scip>)
    target_link_libraries(scipopt lpi nlpi tpi objscip)
endif()

target_link_libraries(scipopt ${ZLIB_LIBRARIES} ${Readline_LIBRARY} ${GMP_LIBRARIES} ${ZIMPL_PIC_LIBRARIES} ${LPS_PIC_LIBRARIES} ${THREAD_LIBRARIES} ${NLPI_LIBRARIES})

add_executable(scip-bin main.c $<TARGET_OBJECTS:scip>)
target_link_libraries(scip-bin objscip nlpi lpi tpi ${ZLIB_LIBRARIES} ${LPS_LIBRARIES} ${Readline_LIBRARY} ${GMP_LIBRARIES} ${ZIMPL_LIBRARIES} ${THREAD_LIBRARIES})

# we cannot have the same name for library and binary, so we redefine it here
set_target_properties(scip-bin PROPERTIES OUTPUT_NAME scip)

# install the header files of scip
install(FILES ${lpiheaders} DESTINATION include/lpi)
install(FILES ${dijkstraheaders} DESTINATION include/dijkstra)
install(FILES ${nlpiheaders} DESTINATION include/nlpi)
install(FILES ${objscipheaders} DESTINATION include/objscip)
install(FILES ${scipheaders} DESTINATION include/scip)
install(FILES ${tcliqueheaders} DESTINATION include/tclique)
install(FILES ${tinycthreadheader} DESTINATION include/tinycthread)
install(FILES ${tpiheaders} DESTINATION include/tpi)
install(FILES ${xmlheaders} DESTINATION include/xml)
install(FILES ${blockmemshellheader} DESTINATION include/blockmemshell)

# install the binary and the library to appropriate lcoations and add them to an export group
install(TARGETS scip-bin scipopt EXPORT scip-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

# Add all targets to the build-tree export set
export(TARGETS scip-bin scipopt
  FILE "${PROJECT_BINARY_DIR}/scip-targets.cmake")

#configure the config file for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/scip-config.cmake.in
  "${PROJECT_BINARY_DIR}/scip-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${PROJECT_SOURCE_DIR}/scip-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/scip-config.cmake" @ONLY)

# install the targets of the scip export group and the config file so that other projects
# can link easily against scip
install(EXPORT scip-targets FILE scip-targets.cmake DESTINATION lib/cmake/scip)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/scip-config.cmake" DESTINATION lib/cmake/scip)
