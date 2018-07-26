#! /bin/bash

#
# brief descriptions in an array called 'descriptions'
#
declare -A descriptions

descriptions[scip_bandit]="public methods for bandit algorithms"
descriptions[scip_benders]="public methods for Benders decomposition"
descriptions[scip_branch]="public methods for branching rule plugins and branching"
descriptions[scip_compr]="public methods for compression plugins"
descriptions[scip_concurrent]="public methods for concurrent solving mode"
descriptions[scip_conflict]="public methods for conflict handler plugins and conflict analysis"
descriptions[scip_cons]="public methods for constraint handler plugins and constraints"
descriptions[scip_copy]="public methods for problem copies"
descriptions[scip_cut]="public methods for cuts and aggregation rows"
descriptions[scip_datastructures]="public methods for data structures"
descriptions[scip_debug]="public methods for debugging"
descriptions[scip_dialog]="public methods for dialog handler plugins"
descriptions[scip_disp]="public methods for display handler plugins"
descriptions[scip_event]="public methods for event handler plugins and event handlers"
descriptions[scip_expr]="public methods for expression handlers"
descriptions[scip_general]="general public methods"
descriptions[scip_heur]="public methods for primal heuristic plugins and divesets"
descriptions[scip_lp]="public methods for the LP relaxation, rows and columns"
descriptions[scip_mem]="public methods for memory management"
descriptions[scip_message]="public methods for message handling"
descriptions[scip_nlp]="public methods for nonlinear relaxations"
descriptions[scip_nodesel]="public methods for node selector plugins"
descriptions[scip_nonlinear]="public methods for nonlinear functions"
descriptions[scip_numerics]="public methods for numerical tolerances"
descriptions[scip_param]="public methods for SCIP parameter handling"
descriptions[scip_presol]="public methods for presolving plugins"
descriptions[scip_pricer]="public methods for variable pricer plugins"
descriptions[scip_prob]="public methods for global and local (sub)problems"
descriptions[scip_probing]="public methods for the probing mode"
descriptions[scip_prop]="public methods for propagator plugins"
descriptions[scip_randnumgen]="public methods for random numbers"
descriptions[scip_reader]="public methods for reader plugins"
descriptions[scip_relax]="public methods for relaxator plugins"
descriptions[scip_reopt]="public methods for reoptimization"
descriptions[scip_sepa]="public methods for separator plugins"
descriptions[scip_sol]="public methods for solutions"
descriptions[scip_solve]="public solving methods"
descriptions[scip_solvingstats]="public methods for querying solving statistics"
descriptions[scip_table]="public methods for statistics table plugins"
descriptions[scip_timing]="public methods for timing"
descriptions[scip_tree]="public methods for the branch-and-bound tree"
descriptions[scip_validation]="public methods for validation"
descriptions[scip_var]="public methods for SCIP variables"

echo "substituting doxygen brief descriptions..."
for i in newfiles/scip/*
do
    filename_no_extension=$i
    for extension in c h
    do
        filename_no_extension=$(basename $filename_no_extension .${extension})
    done
    description=${descriptions[${filename_no_extension}]}
    sed -i "s/ \\* @brief.*/ * @brief  $description/g" $i
    echo $i
done
echo "finished substituting doxygen brief descriptions"