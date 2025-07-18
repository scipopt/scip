#
# interface function for setting common library properties
#
include(GNUInstallDirs)

function(setLibProperties targetname outputname)
    set_target_properties(${targetname} PROPERTIES
        OUTPUT_NAME ${outputname})
endfunction(setLibProperties)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT WIN32)
    set(SCIP_COMP_OPTIONS ${SCIP_COMP_OPTIONS} "$<$<COMPILE_LANGUAGE:C>:-D_XOPEN_SOURCE=600>")
endif()

#configure the scipbuildflags.c source file for the build tree
configure_file(${PROJECT_SOURCE_DIR}/src/scipbuildflags.c.in
  "${PROJECT_BINARY_DIR}/scipbuildflag.c" @ONLY)

set(scipsources
    ${PROJECT_BINARY_DIR}/scipbuildflag.c
    blockmemshell/memory.c
    scip/message.c
    scip/bitencode.c
    scip/bandit.c
    scip/bandit_epsgreedy.c
    scip/bandit_exp3.c
    scip/bandit_exp3ix.c
    scip/bandit_ucb.c
    scip/branch_allfullstrong.c
    scip/branch_cloud.c
    scip/branch_distribution.c
    scip/branch_fullstrong.c
    scip/branch_gomory.c
    scip/branch_inference.c
    scip/branch_leastinf.c
    scip/branch_lookahead.c
    scip/branch_mostinf.c
    scip/branch_multaggr.c
    scip/branch_nodereopt.c
    scip/branch_pscost.c
    scip/branch_random.c
    scip/branch_relpscost.c
    scip/branch_vanillafullstrong.c
    scip/certificate.cpp
    scip/compr_largestrepr.c
    scip/compr_weakcompr.c
    scip/concsolver_scip.c
    scip/cons_abspower.c
    scip/cons_and.c
    scip/cons_benders.c
    scip/cons_benderslp.c
    scip/cons_bounddisjunction.c
    scip/cons_cardinality.c
    scip/cons_components.c
    scip/cons_conjunction.c
    scip/cons_countsols.c
    scip/cons_cumulative.c
    scip/cons_disjunction.c
    scip/cons_exactlinear.c
    scip/cons_exactsol.c
    scip/cons_fixedvar.c
    scip/cons_indicator.c
    scip/cons_integral.c
    scip/cons_knapsack.c
    scip/cons_linear.c
    scip/cons_linking.c
    scip/cons_logicor.c
    scip/cons_nonlinear.c
    scip/cons_or.c
    scip/cons_orbisack.c
    scip/cons_orbitope.c
    scip/cons_orbitope_full.c
    scip/cons_orbitope_pp.c
    scip/cons_pseudoboolean.c
    scip/cons_quadratic.c
    scip/cons_setppc.c
    scip/cons_soc.c
    scip/cons_sos1.c
    scip/cons_sos2.c
    scip/cons_superindicator.c
    scip/cons_symresack.c
    scip/cons_varbound.c
    scip/cons_xor.c
    scip/cutsel_ensemble.c
    scip/cutsel_hybrid.c
    scip/cutsel_dynamic.c
    scip/dcmp.c
    scip/dialog_default.c
    scip/event_globalbnd.c
    scip/event_estim.c
    scip/event_shadowtree.c
    scip/event_softtimelimit.c
    scip/event_solvingphase.c
    scip/expr_abs.c
    scip/expr_entropy.c
    scip/expr_erf.c
    scip/expr_exp.c
    scip/expr_log.c
    scip/expr_pow.c
    scip/expr_product.c
    scip/expr_sum.c
    scip/expr_trig.c
    scip/expr_value.c
    scip/expr_var.c
    scip/expr_varidx.c
    ${exprinterpret}
    scip/disp_default.c
    scip/heur_actconsdiving.c
    scip/heur_adaptivediving.c
    scip/heur_bound.c
    scip/heur_clique.c
    scip/heur_coefdiving.c
    scip/heur_completesol.c
    scip/heur_conflictdiving.c
    scip/heur_crossover.c
    scip/heur_dins.c
    scip/heur_distributiondiving.c
    scip/heur_dks.c
    scip/heur_dps.c
    scip/heur_dualval.c
    scip/heur_farkasdiving.c
    scip/heur_feaspump.c
    scip/heur_fixandinfer.c
    scip/heur_fracdiving.c
    scip/heur_gins.c
    scip/heur_guideddiving.c
    scip/heur_indicator.c
    scip/heur_indicatordiving.c
    scip/heur_intdiving.c
    scip/heur_intshifting.c
    scip/heur_linesearchdiving.c
    scip/heur_localbranching.c
    scip/heur_lpface.c
    scip/heur_locks.c
    scip/heur_alns.c
    scip/heur_mutation.c
    scip/heur_multistart.c
    scip/heur_mpec.c
    scip/heur_nlpdiving.c
    scip/heur_objpscostdiving.c
    scip/heur_octane.c
    scip/heur_ofins.c
    scip/heur_oneopt.c
    scip/heur_padm.c
    scip/heur_proximity.c
    scip/heur_pscostdiving.c
    scip/heur_reoptsols.c
    scip/heur_repair.c
    scip/heur_randrounding.c
    scip/heur_rens.c
    scip/heur_rins.c
    scip/heur_rootsoldiving.c
    scip/heur_rounding.c
    scip/heur_scheduler.c
    scip/heur_shiftandpropagate.c
    scip/heur_shifting.c
    scip/heur_simplerounding.c
    scip/heur_subnlp.c
    scip/heur_sync.c
    scip/heur_trivial.c
    scip/heur_trivialnegation.c
    scip/heur_trustregion.c
    scip/heur_trysol.c
    scip/heur_twoopt.c
    scip/heur_undercover.c
    scip/heur_vbounds.c
    scip/heur_veclendiving.c
    scip/heur_zeroobj.c
    scip/heur_zirounding.c
    scip/hypergraph.c
    scip/iisfinder_greedy.c
    scip/message_default.c
    scip/network.c
    scip/nlhdlr_bilinear.c
    scip/nlhdlr_convex.c
    scip/nlhdlr_default.c
    scip/nlhdlr_perspective.c
    scip/nlhdlr_quadratic.c
    scip/nlhdlr_quotient.c
    scip/nlhdlr_signomial.c
    scip/nlhdlr_soc.c
    scip/nlpi_all.c
    ${nlpi}
    scip/nodesel_bfs.c
    scip/nodesel_breadthfirst.c
    scip/nodesel_dfs.c
    scip/nodesel_estimate.c
    scip/nodesel_hybridestim.c
    scip/nodesel_restartdfs.c
    scip/nodesel_uct.c
    scip/presol_milp.cpp
    scip/presol_boundshift.c
    scip/presol_convertinttobin.c
    scip/presol_domcol.c
    scip/presol_dualagg.c
    scip/presol_dualcomp.c
    scip/presol_dualinfer.c
    scip/presol_gateextraction.c
    scip/presol_implics.c
    scip/presol_implint.c
    scip/presol_inttobinary.c
    scip/presol_qpkktref.c
    scip/presol_redvub.c
    scip/presol_sparsify.c
    scip/presol_dualsparsify.c
    scip/presol_stuffing.c
    scip/presol_trivial.c
    scip/presol_tworowbnd.c
    scip/prop_dualfix.c
    scip/prop_genvbounds.c
    scip/prop_nlobbt.c
    scip/prop_obbt.c
    scip/prop_probing.c
    scip/prop_pseudoobj.c
    scip/prop_redcost.c
    scip/prop_rootredcost.c
    scip/prop_symmetry.c
    scip/prop_vbounds.c
    scip/prop_symmetry.c
    scip/prop_sync.c
    scip/reader_bnd.c
    scip/reader_ccg.c
    scip/reader_cip.c
    scip/reader_cnf.c
    scip/reader_cor.c
    scip/reader_dec.c
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
    scip/reader_sto.c
    scip/reader_smps.c
    scip/reader_tim.c
    scip/reader_wbo.c
    scip/reader_zpl.c
    scip/sepa_cgmip.c
    scip/sepa_clique.c
    scip/sepa_closecuts.c
    scip/sepa_aggregation.c
    scip/sepa_convexproj.c
    scip/sepa_disjunctive.c
    scip/sepa_eccuts.c
    scip/sepa_flower.c
    scip/sepa_gauge.c
    scip/sepa_gomory.c
    scip/sepa_impliedbounds.c
    scip/sepa_interminor.c
    scip/sepa_intobj.c
    scip/sepa_lagromory.c
    scip/sepa_mcf.c
    scip/sepa_minor.c
    scip/sepa_mixing.c
    scip/sepa_oddcycle.c
    scip/sepa_rapidlearning.c
    scip/sepa_rlt.c
    scip/sepa_zerohalf.c
    scip/treemodel.c
    scip/benders_default.c
    scip/benderscut_feas.c
    scip/benderscut_feasalt.c
    scip/benderscut_int.c
    scip/benderscut_nogood.c
    scip/benderscut_opt.c
    scip/table_default.c
    scip/benders.c
    scip/benderscut.c
    scip/bendersdefcuts.c
    scip/relax_benders.c
    scip/boundstore.c
    scip/branch.c
    scip/clock.c
    scip/concsolver.c
    scip/concurrent.c
    scip/conflict_general.c
    scip/conflict_graphanalysis.c
    scip/conflict_dualproofanalysis.c
    scip/conflict_resolution.c
    scip/conflictstore.c
    scip/cons.c
    scip/cutpool.c
    scip/cuts.c
    scip/cutsel.c
    scip/debug.c
    scip/dialog.c
    scip/disp.c
    scip/event.c
    scip/expr.c
    scip/exprcurv.c
    scip/expriter.c
    scip/fileio.c
    scip/heur.c
    scip/heuristics.c
    scip/compr.c
    scip/history.c
    scip/iisfinder.c
    scip/implics.c
    scip/interrupt.c
    scip/intervalarith.c
    scip/lapack_calls.c
    scip/lp.c
    scip/lpexact.c
    scip/lpexact_bounding.c
    scip/matrix.c
    scip/mem.c
    scip/misc.c
    scip/misc_linear.c
    scip/misc_rowprep.c
    scip/nlhdlr.c
    scip/nlp.c
    scip/nlpi.c
    scip/nlpioracle.c
    scip/nodesel.c
    scip/paramset.c
    scip/presol.c
    scip/presolve.c
    scip/pricestore.c
    scip/pricer.c
    scip/primal.c
    scip/prob.c
    scip/prop.c
    scip/rational.cpp
    scip/reader.c
    scip/relax.c
    scip/reopt.c
    scip/retcode.c
    scip/scip_bandit.c
    scip/scip_benders.c
    scip/scip_branch.c
    scip/scip_certificate.c
    scip/scip_compr.c
    scip/scip_concurrent.c
    scip/scip_conflict.c
    scip/scip_cons.c
    scip/scip_copy.c
    scip/scip_cut.c
    scip/scip_cutsel.c
    scip/scip_datastructures.c
    scip/scip_datatree.c
    scip/scip_debug.c
    scip/scip_dcmp.c
    scip/scip_dialog.c
    scip/scip_disp.c
    scip/scip_event.c
    scip/scip_exact.c
    scip/scip_expr.c
    scip/scip_general.c
    scip/scip_heur.c
    scip/scip_iisfinder.c
    scip/scip_lp.c
    scip/scip_lpexact.c
    scip/scip_mem.c
    scip/scip_message.c
    scip/scip_nlp.c
    scip/scip_nlpi.c
    scip/scip_nodesel.c
    scip/scip_numerics.c
    scip/scip_param.c
    scip/scip_presol.c
    scip/scip_pricer.c
    scip/scip_prob.c
    scip/scip_probing.c
    scip/scip_prop.c
    scip/scip_randnumgen.c
    scip/scip_reader.c
    scip/scip_relax.c
    scip/scip_reopt.c
    scip/scip_sepa.c
    scip/scip_sol.c
    scip/scip_solve.c
    scip/scip_solvingstats.c
    scip/scip_table.c
    scip/scip_timing.c
    scip/scip_tree.c
    scip/scip_validation.c
    scip/scip_var.c
    scip/datatree.c
    scip/scipdefplugins.c
    scip/scipcoreplugins.c
    scip/scipgithash.c
    scip/scipshell.c
    scip/sepa.c
    scip/sepastore.c
    scip/sepastoreexact.c
    scip/set.c
    scip/sol.c
    scip/solve.c
    scip/stat.c
    scip/symmetry.c
    scip/symmetry_graph.c
    scip/symmetry_orbitopal.c
    scip/symmetry_orbital.c
    scip/symmetry_lexred.c
    scip/syncstore.c
    scip/table.c
    scip/tree.c
    scip/var.c
    scip/visual.c
    scip/rbtree.c
    tclique/tclique_branch.c
    tclique/tclique_coloring.c
    tclique/tclique_graph.c
    dijkstra/dijkstra.c
    xml/xmlparse.c
    ${amplsrc}
)

set(objscipsources
    objscip/objbenders.cpp
    objscip/objbenderscut.cpp
    objscip/objbranchrule.cpp
    objscip/objconshdlr.cpp
    objscip/objcutsel.cpp
    objscip/objdialog.cpp
    objscip/objdisp.cpp
    objscip/objeventhdlr.cpp
    objscip/objexprhdlr.cpp
    objscip/objheur.cpp
    objscip/objiisfinder.cpp
    objscip/objmessagehdlr.cpp
    objscip/objnodesel.cpp
    objscip/objpresol.cpp
    objscip/objpricer.cpp
    objscip/objprobdata.cpp
    objscip/objprop.cpp
    objscip/objreader.cpp
    objscip/objrelax.cpp
    objscip/objsepa.cpp
    objscip/objtable.cpp
    objscip/objvardata.cpp
)

set(lpisources
   ${lpi}
   )

set(symsources
   ${sym}
   )

set(lpiexactsources
   ${lpiexact}
   )

if(SCIP_WITH_EXACTSOLVE)
set(rectlusources
   rectlu/rectlu_factor.c
   rectlu/rectlu_num.c
   )

set(rectluheaders
    rectlu/rectlu.h
    rectlu/rectlu_num.h
    rectlu/rectlu_factor.h)
endif()

set(blockmemshellheader
    blockmemshell/memory.h
)

set(lpiheaders
    lpi/lpi.h
    lpi/type_lpi.h
)

set(lpiexactheaders
    lpiexact/lpiexact.h
    lpiexact/type_lpiexact.h
)

set(dijkstraheaders
    dijkstra/dijkstra.h
)

set(objscipheaders
    objscip/objbenders.h
    objscip/objbenderscut.h
    objscip/objbranchrule.h
    objscip/objcloneable.h
    objscip/objconshdlr.h
    objscip/objcutsel.h
    objscip/objdialog.h
    objscip/objdisp.h
    objscip/objeventhdlr.h
    objscip/objexprhdlr.h
    objscip/objheur.h
    objscip/objiisfinder.h
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
    objscip/objtable.h
    objscip/objvardata.h
    objscip/type_objcloneable.h
    objscip/type_objprobcloneable.h
)

set(scipheaders
    scip/bandit.h
    scip/bandit_epsgreedy.h
    scip/bandit_exp3.h
    scip/bandit_exp3ix.h
    scip/bandit_ucb.h
    scip/benders.h
    scip/benders_default.h
    scip/benderscut.h
    scip/benderscut_feas.h
    scip/benderscut_feasalt.h
    scip/benderscut_int.h
    scip/benderscut_nogood.h
    scip/benderscut_opt.h
    scip/bendersdefcuts.h
    scip/relax_benders.h
    scip/bitencode.h
    scip/boundstore.h
    scip/branch_allfullstrong.h
    scip/branch_cloud.h
    scip/branch_distribution.h
    scip/branch_fullstrong.h
    scip/branch_gomory.h
    scip/branch.h
    scip/branch_inference.h
    scip/branch_leastinf.h
    scip/branch_lookahead.h
    scip/branch_mostinf.h
    scip/branch_multaggr.h
    scip/branch_nodereopt.h
    scip/branch_pscost.h
    scip/branch_random.h
    scip/branch_relpscost.h
    scip/branch_vanillafullstrong.h
    scip/certificate.h
    scip/clock.h
    scip/compr.h
    scip/compr_largestrepr.h
    scip/compr_weakcompr.h
    scip/concsolver.h
    scip/concsolver_scip.h
    scip/concurrent.h
    scip/conflict.h
    scip/conflict_graphanalysis.h
    scip/conflict_dualproofanalysis.h
    scip/conflict_resolution.h
    scip/conflict_general.h
    scip/conflictstore.h
    scip/cons_abspower.h
    scip/cons_and.h
    scip/cons_benders.h
    scip/cons_benderslp.h
    scip/cons_bounddisjunction.h
    scip/cons_cardinality.h
    scip/cons_components.h
    scip/cons_conjunction.h
    scip/cons_countsols.h
    scip/cons_cumulative.h
    scip/cons_disjunction.h
    scip/cons_exactlinear.h
    scip/cons_exactsol.h
    scip/cons_fixedvar.h
    scip/cons.h
    scip/cons_indicator.h
    scip/cons_integral.h
    scip/cons_knapsack.h
    scip/cons_linear.h
    scip/cons_linking.h
    scip/cons_logicor.h
    scip/cons_nonlinear.h
    scip/cons_orbisack.h
    scip/cons_orbitope.h
    scip/cons_orbitope_full.h
    scip/cons_orbitope_pp.h
    scip/cons_or.h
    scip/cons_pseudoboolean.h
    scip/cons_quadratic.h
    scip/cons_setppc.h
    scip/cons_soc.h
    scip/cons_sos1.h
    scip/cons_sos2.h
    scip/cons_superindicator.h
    scip/cons_symresack.h
    scip/cons_varbound.h
    scip/cons_xor.h
    scip/cutpool.h
    scip/cuts.h
    scip/cutsel.h
    scip/cutsel_ensemble.h
    scip/cutsel_hybrid.h
    scip/cutsel_dynamic.h
    scip/dbldblarith.h
    scip/debug.h
    scip/dcmp.h
    scip/def.h
    scip/dialog_default.h
    scip/dialog.h
    scip/disp_default.h
    scip/disp.h
    scip/event_globalbnd.h
    scip/event.h
    scip/event_estim.h
    scip/event_shadowtree.h
    scip/event_softtimelimit.h
    scip/event_solvingphase.h
    scip/expr.h
    scip/expr_abs.h
    scip/expr_entropy.h
    scip/expr_erf.h
    scip/expr_exp.h
    scip/expr_log.h
    scip/expr_pow.h
    scip/expr_product.h
    scip/expr_sum.h
    scip/expr_trig.h
    scip/expr_value.h
    scip/expr_var.h
    scip/expr_varidx.h
    scip/exprinterpret.h
    scip/heur_actconsdiving.h
    scip/heur_adaptivediving.h
    scip/heur_bound.h
    scip/heur_clique.h
    scip/heur_coefdiving.h
    scip/heur_completesol.h
    scip/heur_conflictdiving.h
    scip/heur_crossover.h
    scip/heur_dins.h
    scip/heur_distributiondiving.h
    scip/heur_dks.h
    scip/heur_dps.h
    scip/heur_dualval.h
    scip/heur_farkasdiving.h
    scip/heur_feaspump.h
    scip/heur_fixandinfer.h
    scip/heur_fracdiving.h
    scip/heur_gins.h
    scip/heur_guideddiving.h
    scip/heur.h
    scip/heur_indicator.h
    scip/heur_indicatordiving.h
    scip/heur_intdiving.h
    scip/heur_intshifting.h
    scip/heuristics.h
    scip/heur_linesearchdiving.h
    scip/heur_localbranching.h
    scip/heur_locks.h
    scip/heur_alns.h
    scip/heur_lpface.h
    scip/heur_multistart.h
    scip/heur_mutation.h
    scip/heur_mpec.h
    scip/heur_nlpdiving.h
    scip/heur_objpscostdiving.h
    scip/heur_octane.h
    scip/heur_ofins.h
    scip/heur_oneopt.h
    scip/heur_padm.h
    scip/heur_proximity.h
    scip/heur_pscostdiving.h
    scip/heur_randrounding.h
    scip/heur_rens.h
    scip/heur_reoptsols.h
    scip/heur_repair.h
    scip/heur_rins.h
    scip/heur_rootsoldiving.h
    scip/heur_rounding.h
    scip/heur_scheduler.h
    scip/heur_shiftandpropagate.h
    scip/heur_shifting.h
    scip/heur_simplerounding.h
    scip/heur_subnlp.h
    scip/heur_sync.h
    scip/heur_trivial.h
    scip/heur_trivialnegation.h
    scip/heur_trustregion.h
    scip/heur_trysol.h
    scip/heur_twoopt.h
    scip/heur_undercover.h
    scip/heur_vbounds.h
    scip/heur_veclendiving.h
    scip/heur_zeroobj.h
    scip/heur_zirounding.h
    scip/history.h
    scip/hypergraph.h
    scip/iisfinder.h
    scip/iisfinder_greedy.h
    scip/implics.h
    scip/interrupt.h
    scip/intervalarith.h
    scip/lapack_calls.h
    scip/lp.h
    scip/lpexact.h
    scip/lpexact_bounding.h
    scip/mem.h
    scip/message_default.h
    scip/message.h
    scip/multiprecision.hpp
    scip/misc.h
    scip/nlhdlr_bilinear.h
    scip/nlhdlr_convex.h
    scip/nlhdlr_default.h
    scip/nlhdlr_perspective.h
    scip/nlhdlr_quadratic.h
    scip/nlhdlr_quotient.h
    scip/nlhdlr_signomial.h
    scip/nlhdlr_soc.h
    scip/nlhdlr.h
    scip/nlp.h
    scip/nlpi.h
    scip/nlpioracle.h
    scip/nlpi_all.h
    scip/nlpi_filtersqp.h
    scip/nlpi_ipopt.h
    scip/nlpi_worhp.h
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
    scip/presol_milp.h
    scip/presol_convertinttobin.h
    scip/presol_domcol.h
    scip/presol_dualagg.h
    scip/presol_dualcomp.h
    scip/presol_dualinfer.h
    scip/presol_gateextraction.h
    scip/presol.h
    scip/presol_implics.h
    scip/presol_implint.h
    scip/presol_inttobinary.h
    scip/presol_qpkktref.h
    scip/presol_redvub.h
    scip/presol_sparsify.h
    scip/presol_dualsparsify.h
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
    scip/prop_symmetry.h
    scip/prop_sync.h
    scip/prop_vbounds.h
    scip/pub_branch.h
    scip/pub_bandit.h
    scip/pub_bandit_epsgreedy.h
    scip/pub_bandit_exp3.h
    scip/pub_bandit_exp3ix.h
    scip/pub_bandit_ucb.h
    scip/pub_benders.h
    scip/pub_benderscut.h
    scip/pub_compr.h
    scip/pub_conflict.h
    scip/pub_cons.h
    scip/pub_cutpool.h
    scip/pub_cutsel.h
    scip/pub_datatree.h
    scip/pub_dcmp.h
    scip/pub_dialog.h
    scip/pub_disp.h
    scip/pub_event.h
    scip/pub_expr.h
    scip/pub_fileio.h
    scip/pub_heur.h
    scip/pub_history.h
    scip/pub_iisfinder.h
    scip/pub_implics.h
    scip/pub_lp.h
    scip/pub_lpexact.h
    scip/pub_matrix.h
    scip/pub_message.h
    scip/pub_misc.h
    scip/pub_misc_linear.h
    scip/pub_misc_rowprep.h
    scip/pub_misc_select.h
    scip/pub_misc_sort.h
    scip/pub_network.h
    scip/pub_nlhdlr.h
    scip/pub_nlp.h
    scip/pub_nlpi.h
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
    scip/pub_table.h
    scip/pub_tree.h
    scip/pub_var.h
    scip/rational.h
    scip/rbtree.h
    scip/reader_bnd.h
    scip/reader_ccg.h
    scip/reader_cip.h
    scip/reader_cnf.h
    scip/reader_cor.h
    scip/reader_dec.h
    scip/reader_diff.h
    scip/reader_fix.h
    scip/reader_fzn.h
    scip/reader_gms.h
    scip/reader.h
    scip/reader_lp.h
    scip/reader_mps.h
    scip/reader_mst.h
    scip/reader_nl.h
    scip/reader_opb.h
    scip/reader_osil.h
    scip/reader_pbm.h
    scip/reader_pip.h
    scip/reader_ppm.h
    scip/reader_rlp.h
    scip/reader_sol.h
    scip/reader_smps.h
    scip/reader_sto.h
    scip/reader_tim.h
    scip/reader_wbo.h
    scip/reader_zpl.h
    scip/relax.h
    scip/reopt.h
    scip/retcode.h
    scip/scipbuildflags.h
    scip/scipcoreplugins.h
    scip/scipdefplugins.h
    scip/scipgithash.h
    scip/scip.h
    scip/scip_bandit.h
    scip/scip_benders.h
    scip/scip_branch.h
    scip/scip_certificate.h
    scip/scip_compr.h
    scip/scip_concurrent.h
    scip/scip_conflict.h
    scip/scip_cons.h
    scip/scip_copy.h
    scip/scip_cut.h
    scip/scip_cutsel.h
    scip/scip_datastructures.h
    scip/scip_datatree.h
    scip/scip_debug.h
    scip/scip_dcmp.h
    scip/scip_dialog.h
    scip/scip_disp.h
    scip/scip_event.h
    scip/scip_exact.h
    scip/scip_expr.h
    scip/scip_general.h
    scip/scip_heur.h
    scip/scip_iisfinder.h
    scip/scip_lp.h
    scip/scip_lpexact.h
    scip/scip_mem.h
    scip/scip_message.h
    scip/scip_nlp.h
    scip/scip_nlpi.h
    scip/scip_nodesel.h
    scip/scip_numerics.h
    scip/scip_param.h
    scip/scip_presol.h
    scip/scip_pricer.h
    scip/scip_prob.h
    scip/scip_probing.h
    scip/scip_prop.h
    scip/scip_randnumgen.h
    scip/scip_reader.h
    scip/scip_relax.h
    scip/scip_reopt.h
    scip/scip_sepa.h
    scip/scip_sol.h
    scip/scip_solve.h
    scip/scip_solvingstats.h
    scip/scip_table.h
    scip/scip_timing.h
    scip/scip_tree.h
    scip/scip_validation.h
    scip/scip_var.h
    scip/scipshell.h
    scip/sepa_cgmip.h
    scip/sepa_clique.h
    scip/sepa_closecuts.h
    scip/sepa_aggregation.h
    scip/sepa_convexproj.h
    scip/sepa_disjunctive.h
    scip/sepa_eccuts.h
    scip/sepa_flower.h
    scip/sepa_gauge.h
    scip/sepa_gomory.h
    scip/sepa.h
    scip/sepa_impliedbounds.h
    scip/sepa_interminor.h
    scip/sepa_intobj.h
    scip/sepa_lagromory.h
    scip/sepa_mcf.h
    scip/sepa_minor.h
    scip/sepa_mixing.h
    scip/sepa_oddcycle.h
    scip/sepa_rapidlearning.h
    scip/sepa_rlt.h
    scip/sepastore.h
    scip/sepastoreexact.h
    scip/sepa_zerohalf.h
    scip/set.h
    scip/sol.h
    scip/solve.h
    scip/stat.h
    scip/struct_bandit.h
    scip/struct_benders.h
    scip/struct_benderscut.h
    scip/struct_branch.h
    scip/struct_clock.h
    scip/struct_compr.h
    scip/struct_concsolver.h
    scip/struct_concurrent.h
    scip/struct_conflict.h
    scip/struct_conflictstore.h
    scip/struct_cons.h
    scip/struct_cutpool.h
    scip/struct_cuts.h
    scip/struct_cutsel.h
    scip/struct_datatree.h
    scip/struct_dcmp.h
    scip/struct_dialog.h
    scip/struct_disp.h
    scip/struct_event.h
    scip/struct_expr.h
    scip/struct_heur.h
    scip/struct_history.h
    scip/struct_hypergraph.h
    scip/struct_iisfinder.h
    scip/struct_implics.h
    scip/struct_lp.h
    scip/struct_lpexact.h
    scip/struct_matrix.h
    scip/struct_mem.h
    scip/struct_message.h
    scip/struct_misc.h
    scip/struct_nlhdlr.h
    scip/struct_nlp.h
    scip/struct_nlpi.h
    scip/struct_nodesel.h
    scip/struct_paramset.h
    scip/struct_presol.h
    scip/struct_pricer.h
    scip/struct_pricestore.h
    scip/struct_primal.h
    scip/struct_prob.h
    scip/struct_prop.h
    scip/struct_rational.h
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
    scip/struct_table.h
    scip/struct_tree.h
    scip/struct_var.h
    scip/struct_visual.h
    scip/datatree.h
    scip/symmetry.h
    scip/symmetry_graph.h
    scip/symmetry_orbitopal.h
    scip/symmetry_orbital.h
    scip/symmetry_lexred.h
    scip/syncstore.h
    scip/table_default.h
    scip/table.h
    scip/tree.h
    scip/treemodel.h
    scip/type_bandit.h
    scip/type_benders.h
    scip/type_benderscut.h
    scip/type_branch.h
    scip/type_certificate.h
    scip/type_clock.h
    scip/type_compr.h
    scip/type_concsolver.h
    scip/type_concurrent.h
    scip/type_conflict.h
    scip/type_conflictstore.h
    scip/type_cons.h
    scip/type_cutpool.h
    scip/type_cuts.h
    scip/type_cutsel.h
    scip/type_datatree.h
    scip/type_dcmp.h
    scip/type_dialog.h
    scip/type_disp.h
    scip/type_event.h
    scip/type_expr.h
    scip/type_exprinterpret.h
    scip/type_heur.h
    scip/type_history.h
    scip/type_hypergraph.h
    scip/type_iisfinder.h
    scip/type_implics.h
    scip/type_interrupt.h
    scip/type_lp.h
    scip/type_lpexact.h
    scip/type_matrix.h
    scip/type_mem.h
    scip/type_message.h
    scip/type_misc.h
    scip/type_nlhdlr.h
    scip/type_nlp.h
    scip/type_nlpi.h
    scip/type_nodesel.h
    scip/type_paramset.h
    scip/type_presol.h
    scip/type_pricer.h
    scip/type_pricestore.h
    scip/type_primal.h
    scip/type_prob.h
    scip/type_prop.h
    scip/type_rational.h
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
    scip/type_table.h
    scip/type_timing.h
    scip/type_tree.h
    scip/type_var.h
    scip/type_visual.h
    scip/var.h
    scip/visual.h
)

set(symheaders
  symmetry/build_dejavu_graph.h
  symmetry/compute_symmetry.h
  symmetry/struct_symmetry.h
  symmetry/type_symmetry.h
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
    tpi/type_tpi.h
)

set(xmlheaders
    xml/xmldef.h
    xml/xml.h
)

# all source files should be compiled with a c++ compiler
if(CXXONLY)
    set_source_files_properties(main.c ${scipsources} ${objscipsources} ${lpisources} ${lpiexactsources} ${rectlusources} ${tpisources} ${symsources} PROPERTIES LANGUAGE CXX)

    # for the clang compiler this suppresses the warnings about treating 'c' input as 'c++' when CXXONLY is enabled
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
       set(SCIP_COMP_OPTIONS ${SCIP_COMP_OPTIONS} -x c++)
    endif()
endif()

add_library(libobjscip EXCLUDE_FROM_ALL ${objscipsources})
setLibProperties(libobjscip "objscip")

add_library(liblpi EXCLUDE_FROM_ALL
    ${lpisources}
    scip/bitencode.c
    blockmemshell/memory.c
    scip/message.c)
setLibProperties(liblpi "lpi")

add_library(liblpiexact EXCLUDE_FROM_ALL
    ${lpiexactsources})
setLibProperties(liblpiexact "lpiexact")

add_library(librectlu EXCLUDE_FROM_ALL
    ${rectlusources}
    scip/bitencode.c)
setLibProperties(librectlu "rectlu")

add_library(libscip ${scipsources} ${objscipsources} ${lpisources} ${lpiexactsources} ${rectlusources} ${tpisources} ${symsources})
if(MSVC)
   # msvc otherwise is not smart enough and tries to link the scip.exe binary as a library
   setLibProperties(libscip "libscip")
else()
   setLibProperties(libscip "scip")
endif()

target_include_directories(libscip BEFORE PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(libscip
    PRIVATE
    ${SYM_LIBRARIES}
    $<$<BOOL:${SCIP_WITH_ZLIB}>:ZLIB::ZLIB>
    ${Readline_LIBRARY}
    ${THREAD_LIBRARIES}
    ${NLPI_LIBRARIES}
    ${PAPILO_IMPORTED_TARGETS}
)
target_link_libraries(libscip
    PUBLIC
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

set(SOPLEX_NEEDED 0)
set(ZIMPL_NEEDED 0)
if(SHARED)
    target_link_libraries(libscip PRIVATE ${ZIMPL_PIC_LIBRARIES} ${LPS_PIC_LIBRARIES})
    add_executable(scip main.c ${scipsources} ${objscipsources} ${lpisources} ${lpiexactsources} ${rectlusources} ${tpisources} ${symsources})
else()
    if(SOPLEX_FOUND)
        set(SOPLEX_NEEDED 1)
    endif()
    if(ZIMPL_FOUND)
        set(ZIMPL_NEEDED 1)
    endif()
    target_link_libraries(libscip PRIVATE ${ZIMPL_LIBRARIES} ${LPS_LIBRARIES})
    add_executable(scip main.c)
    target_link_libraries(scip PUBLIC libscip)
endif()

add_executable(SCIP::SCIP ALIAS scip)
add_library(SCIP::libscip ALIAS libscip)

target_compile_options(scip PRIVATE ${SCIP_COMP_OPTIONS})
target_compile_options(libscip PRIVATE ${SCIP_COMP_OPTIONS})

include(GenerateExportHeader)
generate_export_header(libscip BASE_NAME scip EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/scip/scip_export.h)
target_compile_definitions(scip PRIVATE SCIP_STATIC_DEFINE)

if(SANITIZE)
    find_package(Sanitizers)
    add_sanitizers(libscip)
    add_sanitizers(scip)
    get_target_property(CONF_SANITIZE_FLAGS libscip SANITIZE_FLAGS)
endif(SANITIZE)

target_link_libraries(scip
    PUBLIC
    ${SYM_LIBRARIES}
    $<$<BOOL:${SCIP_WITH_ZLIB}>:ZLIB::ZLIB>
    ${Readline_LIBRARY}
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${ZIMPL_LIBRARIES}
    ${THREAD_LIBRARIES}
    ${NLPI_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${PAPILO_IMPORTED_TARGETS}
    ${LPS_LIBRARIES}
)

add_dependencies(libscip scip_update_githash)
add_dependencies(scip scip_update_githash)

set_target_properties(libscip PROPERTIES
    VERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}
    SOVERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)

target_include_directories(scip BEFORE PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

# set the install rpath to the installed destination
set_target_properties(scip PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH TRUE)

# install the header files of scip
install(FILES ${lpiheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lpi)
install(FILES ${lpiexactheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lpiexact)
install(FILES ${rectluheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rectlu)
install(FILES ${dijkstraheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dijkstra)
install(FILES ${objscipheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/objscip)
install(FILES ${scipheaders} ${PROJECT_BINARY_DIR}/scip/config.h ${PROJECT_BINARY_DIR}/scip/scip_export.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/scip)
install(FILES ${tcliqueheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tclique)
install(FILES ${tinycthreadheader} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tinycthread)
install(FILES ${tpiheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tpi)
install(FILES ${xmlheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/xml)
install(FILES ${symheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/symmetry)
install(FILES ${blockmemshellheader} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/blockmemshell)

# install the binary and the library to appropriate locations and add them to an export group
install(TARGETS scip libscip EXPORT scip-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Add all targets to the build-tree export set
export(TARGETS scip libscip
  FILE "${PROJECT_BINARY_DIR}/scip-targets.cmake")

# make soplex and zimpl dir absolute for the config file
if(SOPLEX_NEEDED)
    get_filename_component(CONF_SOPLEX_DIR ${SOPLEX_DIR} REALPATH BASE_DIR ${CMAKE_SOURCE_DIR})
endif()

if(ZIMPL_NEEDED)
    get_filename_component(CONF_ZIMPL_DIR ${ZIMPL_DIR} REALPATH BASE_DIR ${CMAKE_SOURCE_DIR})
endif()

# Configure scip-config.cmake for the build-tree:
# We add the CMake module path from the sources.
set(EXTRA_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/scip-config.cmake.in
  "${PROJECT_BINARY_DIR}/scip-config.cmake" @ONLY)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
${PROJECT_BINARY_DIR}/scip-config-version.cmake
    VERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
    )

# Configure the scip-config.cmake for the installation.
# We add the current directory of the installed scip-config.cmake as a CMake module path.
set(EXTRA_CMAKE_MODULE_PATH "\${CMAKE_CURRENT_LIST_DIR}")
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
if(SOPLEX_NEEDED)
    set(CONF_SOPLEX_DIR "\${CMAKE_CURRENT_LIST_DIR}/../soplex")
endif()
if(ZIMPL_NEEDED)
    set(CONF_ZIMPL_DIR "\${CMAKE_CURRENT_LIST_DIR}/../zimpl")
endif()
configure_file(${PROJECT_SOURCE_DIR}/scip-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/scip-config.cmake" @ONLY)

# install the targets of the scip export group and the config file so that other projects
# can link easily against scip
install(EXPORT scip-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scip)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/scip-config.cmake"
              ${PROJECT_BINARY_DIR}/scip-config-version.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scip)

# Install FindGMP.cmake module.
if(SCIP_WITH_GMP)
   install(FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindGMP.cmake
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scip)
endif()

# Install FindMPFR.cmake module.
if(SCIP_WITH_MPFR)
   install(FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindMPFR.cmake
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scip)
endif()

