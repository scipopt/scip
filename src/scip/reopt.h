/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reopt.h
 * @brief  data structures and methods for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_REOPT_H__
#define __SCIP_REOPT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_primal.h"
#include "scip/type_retcode.h"
#include "scip/type_reopt.h"

#include "scip/struct_reopt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates reopt data */
extern
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** frees reopt data */
extern
SCIP_RETCODE SCIPreoptFree(
   SCIP_REOPT**          reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,              /**< original primal */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** returns the number constraints added by reoptimization plug-in */
extern
int SCIPreoptGetNAddedConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** returns the number of stored nodes providing an improving feasible LP solution in the current run */
extern
int SCIPreoptGetNFeasNodes(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   );

/** returns the number of stored nodes providing an improving feasible LP solution over all runs */
extern
int SCIPreoptGetNTotalfeasNodes(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   );

/* returns the number of stored nodes that exceeded the cutoff bound in the current run */
extern
int SCIPreoptGetNPrunedNodes(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   );

/* returns the number of stored nodes that exceeded the cutoff bound over all runs */
extern
int SCIPreoptGetNTotalprunedNodes(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   );

/** rerturns the number of reoptimized nodes that were cutoff in the same iteration in the current run */
extern
int SCIPreoptGetNRediedNodes(
   SCIP_REOPT*           reopt                    /*< reoptimization data structure */
   );

/** rerturns the number of reoptimized nodes that were cutoff in the same iteration over all runs */
extern
int SCIPreoptGetNTotalrediedNodes(
   SCIP_REOPT*           reopt                    /*< reoptimization data structure */
   );

/** returns the number of stored nodes with an infeasible LP in the current run */
extern
int SCIPreoptGetNInfNodes(
   SCIP_REOPT*           reopt                    /*< reoptimization data structure */
   );

/** returns the number of stored nodes with an infeasible LP over all runs */
extern
int SCIPreoptGetNTotalinfNodes(
   SCIP_REOPT*           reopt                    /*< reoptimization data structure */
   );

/** returns the number of found infeasible subtrees */
extern
int SCIPreoptGetNInfsubtrees(
   SCIP_REOPT*           reopt                    /*< reoptimization data structure */
   );

/** add a solution to sols */
extern
SCIP_RETCODE SCIPreoptAddSol(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,              /**< original primal */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_SOL*             sol,                     /**< solution to add */
   SCIP_Bool             bestsol,                 /**< is the current solution an optimal solution? */
   SCIP_Bool*            added,                   /**< pointer to store the information if the soltion was added */
   SCIP_VAR**            vars,                    /**< variable array */
   int                   nvars,                   /**< number of variables */
   int                   run                      /**< number of the current run (1,2,...) */
   );

/** add optimal solution */
extern
SCIP_RETCODE SCIPreoptAddOptSol(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SOL*             sol,                     /**< solution to add */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal               /**< original primal */
   );

/** add a run */
extern
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_VAR**            vars,                    /**< trnasformed variables */
   int                   nvars,                   /**< number of transformed variables */
   int                   size                     /**< number of expected solutions */
   );

/** get the number of checked during the reoptimization process */
extern
int SCIPreoptGetNCheckedsols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** update the number of checked during the reoptimization process */
extern
void SCIPreoptSetNCheckedsols(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   ncheckedsols             /**< number of updated solutions */
   );

/** get the number of checked during the reoptimization process */
extern
int SCIPreoptGetNImprovingsols(
   SCIP_REOPT*           reopt                    /**< reoptimization data */
   );

/** update the number of checked during the reoptimization process */
extern
void SCIPreoptSetNImprovingsols(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   nimprovingsols           /**< number of improving solutions */
   );

/** returns number of solution */
extern
int SCIPreoptGetNSolsRun(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   run                      /**< number of the run (1,2,..) */
   );

/** returns number of all solutions */
extern
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** return the stored solutions of a given run */
extern
SCIP_RETCODE SCIPreoptGetSolsRun(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run,                     /**< number of the run (1,2,...) */
   SCIP_SOL**            sols,                    /**< array of solutions to fill */
   int                   allocmem,                /**< length of the array */
   int*                  nsols                    /**< pointer to store the number of added solutions */
   );

/** returns the number of saved solutions overall runs */
extern
int SCIPreoptNSavedSols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** returns the number of reused sols over all runs */
extern
int SCIPreoptNUsedSols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** save objective function */
extern
SCIP_RETCODE SCIPreoptSaveNewObj(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_VAR**            vars,                    /**< array od variables */
   int                   nvars                    /**< number of variables */
   );

/** check if the current and the previous objective are similar enough
 * returns TRUE if we want to restart, otherwise FALSE */
extern
SCIP_RETCODE SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** returns the similarity to the previous objective function */
extern
SCIP_Real SCIPreoptGetSimToPrevious(
      SCIP_REOPT*        reopt                   /**< reoptimization data */
   );

/** returns the similarity to the first objective functions */
extern
SCIP_Real SCIPreoptGetSimToFirst(
      SCIP_REOPT*        reopt                   /**< reoptimization data */
   );

/** return the similarity between two of objective functions of two given runs */
extern
SCIP_Real SCIPreoptGetSim(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run1,                    /**< number of the first run */
   int                   run2                     /**< number of the second run */
   );

/** returns the best solution of the last run */
extern
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** returns the coefficent of variable with index @param idx in run @param run */
extern
SCIP_Real SCIPreoptGetObjCoef(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run,                     /**< number of the run */
   int                   idx                      /**< index of variable */
   );

/** checks the changes of the objective coefficients */
extern
void SCIPreoptGetVarCoefChg(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   varidx,                  /**< index of the variable */
   SCIP_Bool*            negated,                 /**< variable switches the impact, i.e., coefficient changed the sign */
   SCIP_Bool*            entering,                /**< variable enters the objective, i.e., coefficient gets different to 0 */
   SCIP_Bool*            leaving                  /**< variable leaves the objective, i.e., coefficient gets 0 */
   );

/** print the optimal solutions of all previous runs */
extern
SCIP_RETCODE SCIPreoptPrintOptSols(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt                    /**< reoptimization data */
   );

/** return all optimal solutions of the previous runs
 *  depending on the current stage the method copies the solutions into
 *  the origprimal or primal space. That means, all solutions need to be
 *  freed before starting a new iteration!!!
 */
extern
SCIP_RETCODE SCIPreoptGetOptSols(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SOL**            sols,                    /**< array to store solution pointers */
   int*                  nsols                    /**< pointer to store number of solutions */
   );

/** reset marks of stored solutions to not updated */
extern
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** returns the number of stored nodes */
extern
int SCIPreoptGetNNodes(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** returns the number of child nodes of @param node */
int SCIPreoptGetNChildren(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** Save information if infeasible nodes */
extern
SCIP_RETCODE SCIPreoptAddInfNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** check the reason for cut off a node and if necessary store the node */
extern
SCIP_RETCODE SCIPreoptCheckCutoff(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memery */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_EVENTTYPE        eventtype,               /**< eventtype */
   SCIP_LPSOLSTAT        lpsolstat,               /**< solution status of the LP */
   SCIP_Bool             isrootnode,              /**< the node is the root */
   SCIP_Bool             isfocusnode,             /**< the node is the current focus node */
   SCIP_Real             lowerbound,              /**< lower bound of the node */
   int                   effectiverootdepth       /**< effective root depth */
   );

/** store bound changes based on dual information */
extern
SCIP_RETCODE SCIPreoptAddDualBndchg(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_VAR*             var,                     /**< variables */
   SCIP_Real             newval,                  /**< new bound */
   SCIP_Real             oldval                   /**< old bound */
   );

/** returns the number of bound changes based on dual information */
extern
int SCIPreoptGetNDualBndchs(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** returns the number of leave nodes of the subtree induced
 *  by @param node (of the whole tree of node == NULL)
 */
extern
int SCIPreoptGetNLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** returns the child nodes of @param node that need to be
 *  reoptimized next or NULL if @param node is a leaf
 */
extern
SCIP_RETCODE SCIPreoptGetChildIDs(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int*                  childs,                  /**< array to store the child ids */
   int                   mem,                     /**< allocated memory */
   int*                  nchilds                  /**< pointer to store the number of child nodes */
   );

/** returns all leave if of the subtree induced by @param node */
extern
SCIP_RETCODE SCIPreoptGetLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int*                  leaves,                  /**< array to the the ids */
   int                   mem,                     /**< allocated memory */
   int*                  nleaves                  /**< pointer to store the number of leav node */
   );

/** returns the time needed to store the nodes */
extern
SCIP_Real SCIPreoptGetSavingtime(
   SCIP_REOPT*           reopt                    /**< reoptimization data */
   );

/** store a global constraint that should be added at the beginning of the next iteration */
extern
SCIP_RETCODE SCIPreoptAddGlbCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   LOGICORDATA*          consdata,                /**< logicor data */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** add the stored constraints globally to the problem */
extern
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** check if the LP of the given node should be solved or not */
extern
SCIP_RETCODE SCIPreoptGetSolveLP(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_NODE*            node,                    /**< node of the current search tree */
   SCIP_Bool*            solvelp                  /**< pointer to store whether the LP of the node should be solved */
   );

/** reoptimize the node stored at ID id */
extern
SCIP_RETCODE SCIPreoptApply(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PROB*            transprob,               /**< transformed problem */
   SCIP_PROB*            origprob,                /**< original problem */
   SCIP_TREE*            tree,                    /**< branching tree */
   SCIP_LP*              lp,                      /**< current LP */
   SCIP_BRANCHCAND*      branchcand,              /**< branching candidated */
   SCIP_EVENTQUEUE*      eventqueue,              /**< event queue */
   SCIP_NODE*            node_fix,                /**< node curresponding to the fixed part */
   SCIP_NODE*            node_cons,               /**< node corresponding to the pruned part */
   int                   id,                      /**< id of a stored node which should be reoptimized */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** reoptimize the node stored at ID id in the fashion of interdiction branching */
SCIP_RETCODE SCIPreoptApplyInterdiction(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PROB*            transprob,               /**< transformed problem */
   SCIP_PROB*            origprob,                /**< original problem */
   SCIP_TREE*            tree,                    /**< branching tree */
   SCIP_LP*              lp,                      /**< current LP */
   SCIP_BRANCHCAND*      branchcand,              /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,              /**< event queue */
   SCIP_NODE**           nodes,                   /**< array to store created nodes */
   int                   nnodes,                  /**< size of the array */
   int                   id,                      /**< id of a stored node which should be reoptimized */
   int*                  permutation,             /**< permutation of the variable order (within the constraint) */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** delete a node stored in the reopttree */
extern
SCIP_RETCODE SCIPreopttreeDeleteNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id,                      /**< id of the node */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** replace the node stored at ID id by its child nodes */
extern
SCIP_RETCODE SCIPreoptShrinkNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   id                       /**< id of the node */
   );

/** reset the stored information abound bound changes based on dual information */
extern
SCIP_RETCODE SCIPreoptResetDualcons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node,                    /**< node of the search tree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** split the root node and move all children to one of the two resulting nodes */
extern
SCIP_RETCODE SCIPreoptSplitRoot(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** use some kind of interdiction branching to handle dual information with the root node */
extern
SCIP_RETCODE SCIPreoptInterdictRoot(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_VAR**            vars,                    /**< ordered array of variables */
   SCIP_Real*            vals,                    /**< values of variables */
   int                   nvars,                   /**< number of variables */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** create the constraint which splits the node stored at ID id on the basis of
 *  the stored dual information.
 */
extern
void SCIPreoptCreateSplitCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id,                      /**< id of the node */
   LOGICORDATA*          consdata                 /**< data of the split constraint */
   );

/** returns if a node need to be split because some bound changes
 * were based on dual information */
extern
SCIP_Bool SCIPreoptSplitNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node                     /**< node of the search tree */
   );

/** calculates a local similarity of a given node and returns if the subproblem
 * should be solved from scratch */
extern
SCIP_RETCODE SCIPreoptCheckLocalRestart(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_VAR**            origvars,                /**< original variable */
   int                   norigvars,               /**< number of original variables */
   SCIP_Bool*            localrestart             /**< pointer to store whether the induced subtree should solved from
                                                    *  scratch */
);

/** reset the complete tree an set a new search frontier */
extern
SCIP_RETCODE SCIPreoptApplyCompression(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_REOPTNODE**      representatives,         /**< array of representatives */
   int                   nrepresentatives,        /**< number of representatives */
   SCIP_Bool*            success                  /**< pointer to store if the method was successful */
   );

/** add the node @param node to the reopttree */
extern
SCIP_RETCODE SCIPreoptAddNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_REOPTTYPE        reopttype,               /**< reopttype of the node */
   SCIP_LPSOLSTAT        lpsolstat,               /**< solution status of the LP */
   SCIP_Bool             saveafterduals,          /**< store decisions directly after the first dual decision */
   SCIP_Bool             isrootnode,              /**< the node is the root */
   SCIP_Bool             isfocusnode,             /**< the node is the current focus node */
   SCIP_Real             lowerbound               /**< lower bound of the node */
   );

/** add all unprocessed nodes to the reoptimization tree */
extern
SCIP_RETCODE SCIPreoptSaveOpenNodes(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE**           leaves,                  /**< array of open leave nodes */
   int                   nleaves,                 /**< number of open leave nodes */
   SCIP_NODE**           childs,                  /**< array of open children nodes */
   int                   nchilds,                 /**< number of open leave nodes */
   SCIP_NODE**           siblings,                /**< array of open sibling nodes */
   int                   nsiblings                /**< number of open leave nodes */
   );

/** returns the number of restarts */
extern
int SCIPreoptGetNRestarts(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/** return the last lower bound of the stored node */
extern
SCIP_Real SCIPreoptGetNodeLb(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id                       /**< id of the node */
   );

/*
 * methods for reoptnode
 */

/** initialize an empty node */
extern
void SCIPreoptnodeCreateClear(
   SCIP_REOPTNODE*       reoptnode                /**< node of the reopttree */
   );

/** clear the node */
extern
SCIP_RETCODE SCIPreoptnodeClear(
   SCIP_REOPTNODE**      reoptnode,               /**< pointer of reoptnode */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** add a variable to a given reoptnode */
extern
SCIP_RETCODE SCIPreoptnodeAddVar(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   SCIP_VAR*             var,                     /**< variable to add */
   SCIP_Real             val,                     /**< value of the variable */
   SCIP_BOUNDTYPE        boundtype,               /**< boundtype of the variable */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** add a constraint to a given reoptnode */
extern
SCIP_RETCODE SCIPreoptnodeAddCons(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   SCIP_VAR**            consvars,                /**< variables which are part of the constraint */
   SCIP_Real*            consvals,                /**< values of the variables */
   int                   nvars,                   /**< number of variables */
   REOPT_CONSTYPE        constype,                /**< type of the constraint */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/** returns the number of bound changes stored in the reoptnode */
extern
int SCIPreoptnodeGetNVars(
   SCIP_REOPTNODE*       reoptnode                /**< node of the roepttree */
   );

/** set the parent id */
extern
void SCIPreoptnodeSetParentID(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   int                   parentid                 /**< id of the parent node */
   );


/** return the branching path stored at ID id */
extern
void SCIPreoptnodeGetPath(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id,                      /**< id of the stored node */
   SCIP_VAR**            vars,                    /**< array for variables */
   SCIP_Real*            vals,                    /**< array for values */
   SCIP_BOUNDTYPE*       boundtypes,              /**< array for bound types */
   int                   mem,                     /**< allocated memory */
   int*                  nvars,                   /**< pointer to store the number of variables */
   int*                  nafterdualvars           /**< pointer to store the number of variables stored after
                                                    *  the first dual decision */
   );

/** returns all added constraints at ID id */
extern
void SCIPreoptnodeGetConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id,                      /**< id of a stored node */
   SCIP_VAR***           vars,                    /**< 2-dim array of variables */
   SCIP_Real**           vals,                    /**< 2-dim array of values */
   int                   mem,                     /**< allocated memory for constraints */
   int*                  nconss,                  /**< pointer to store the number of constraints */
   int*                  nvars                    /**< pointer to store the number of variables */
   );

/** returns the number of bound changes at the node stored at ID id */
extern
int SCIPreoptnodeGetNConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id                       /**< id of the node */
   );

/** returns the number of stored bound changes based on dual information in the reopttree at ID id */
int SCIPreoptnodeGetNDualBoundChgs(
   SCIP_REOPT*           reopt,                  /**< reoptimization data structure */
   int                   id                      /**< id of a stored node */
   );

/** returns the reopttype of a node stored at ID id */
extern
SCIP_REOPTTYPE SCIPreoptnodeGetType(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id                       /**< id of the node */
   );

#ifdef NDEBUG

/** In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPreoptnodeGetNVars(reoptnode)           ((reoptnode)->nvars + (reoptnode)->nafterdualvars)
#define SCIPreoptnodeGetNConss(reopt, id)          ((((reopt)->reopttree)->reoptnodes[id])->nconss)
#define SCIPreoptnodeGetType(reopt, id)            ((((reopt)->reopttree)->reoptnodes[id])->reopttype)

#define SCIPreoptGetNodeLb(reopt, id)              ((((reopt)->reopttree)->reoptnodes[id])->lowerbound)
#define SCIPreoptGetNRestarts(reopt)               ((reopt)->nrestarts)
#define SCIPreoptGetSimToFirst(reopt)              ((reopt)->simtofirstobj)
#define SCIPreoptGetSimToPrevious(reopt)           ((reopt)->simtolastobj)
#define SCIPreoptGetNImprovingsols(reopt)          ((reopt)->nimprovingsols)
#define SCIPreoptGetNCheckedsols(reopt)            ((reopt)->ncheckedsols)
#define SCIPreoptGetNInfNodes(reopt)               (((reopt)->reopttree)->ninfnodes)
#define SCIPreoptGetNTotalinfNodes(reopt)          (((reopt)->reopttree)->ntotalinfnodes)
#define SCIPreoptGetNRediedNodes(reopt)            (((reopt)->reopttree)->nrediednodes)
#define SCIPreoptGetNTotalrediedNodes(reopt)       (((reopt)->reopttree)->ntotalrediednodes)
#define SCIPreoptGetNPrunedNodes(reopt)            (((reopt)->reopttree)->nprunednodes)
#define SCIPreoptGetNTotalprunedNodes(reopt)       (((reopt)->reopttree)->ntotalprunednodes)
#define SCIPreoptGetNFeasNodes(reopt)              (((reopt)->reopttree)->nfeasnodes)
#define SCIPreoptGetNTotalfeasNodes(reopt)         (((reopt)->reopttree)->ntotalfeasnodes)
#endif

#ifdef __cplusplus
}
#endif

#endif
