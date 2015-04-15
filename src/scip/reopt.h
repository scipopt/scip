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
   SCIP_PRIMAL*          origprimal,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* returns the number constraints added by reoptimization plug-in */
extern
int SCIPreoptGetNAddedConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node
   );

/* returns the number of stored nodes providing an improving feasible LP solution */
extern
int SCIPreoptGetNFeasNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   );

/* returns the number of stored nodes that exceeded the cutoff bound */
extern
int SCIPreoptGetNPrunedNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   );

/* rerturns the number of reoptimized nodes that were cutoff in the same iteration */
extern
int SCIPreoptGetNRediedNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   );

/* returns the number of stored nodes with an infeasible LP */
extern
int SCIPreoptGetNInfNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   );

/** add a solution to sols */
extern
SCIP_RETCODE SCIPreoptAddSol(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_SOL*             sol,                     /**< solution to add */
   SCIP_Bool             bestsol,                 /**< is the current solution an optimal solution? */
   SCIP_Bool*            added,                   /**< pointer to store the information if the soltion was added */
   SCIP_VAR**            vars,
   int                   nvars,
   int                   run                      /**< number of the current run (1,2,...) */
   );

/* add optimal solution */
extern
SCIP_RETCODE SCIPreoptAddOptSol(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SOL*             sol,                     /**< solution to add */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,
   SCIP_PRIMAL*          origprimal
   );

/* add a run */
extern
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_VAR**            origvars,
   int                   norigvars,
   int                   size                     /**< number of expected solutions */
   );

/* get the number of checked during the reoptimization process */
extern
int SCIPreoptGetNCheckedsols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/* update the number of checked during the reoptimization process */
extern
void SCIPreoptSetNCheckedsols(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   ncheckedsols
   );

/* get the number of checked during the reoptimization process */
extern
int SCIPreoptGetNImprovingsols(
   SCIP_REOPT*           reopt
   );

/* update the number of checked during the reoptimization process */
extern
void SCIPreoptSetNImprovingsols(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   nimprovingsols
   );

/* returns number of solution */
extern
int SCIPreoptGetNSolsRun(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   run
   );

/* returns number of all solutions */
extern
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/* return the stored solutions of a given run */
extern
SCIP_RETCODE SCIPreoptGetSolsRun(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   run,
   SCIP_SOL**            sols,
   int                   allocmem,
   int*                  nsols
   );

/* returns the number of saved solutions overall runs */
extern
int SCIPreoptNSavedSols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/* returns the number of reused sols over all runs */
extern
int SCIPreoptNUsedSols(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/* save objective function */
extern
SCIP_RETCODE SCIPreoptSaveNewObj(
      SCIP_REOPT*           reopt,                   /**< reoptimization data */
      SCIP_SET*             set,                     /**< global SCIP settings */
      BMS_BLKMEM*           blkmem,                  /**< block memory */
      SCIP_VAR**            vars,
      int                   nvars
   );

/* check if the current and the previous objective are similar enough
 * returns TRUE if we want to restart, otherwise FALSE */
extern
SCIP_RETCODE SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* returns an array of indices with similar objective functions
 * to obj_idx */
extern
int* SCIPreoptGetSimilarityIdx(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   obj_id,
   int*                  sim_ids,
   int*                  nids
   );

/*
 * returns the similarity to the previous objective function
 */
extern
SCIP_Real SCIPreoptGetSimToPrevious(
      SCIP_REOPT*        reopt                   /**< reoptimization data */
   );

/*
 * returns the similarity to the first objective function
 */
extern
SCIP_Real SCIPreoptGetSimToFirst(
      SCIP_REOPT*        reopt                   /**< reoptimization data */
   );

/*
 * return the similarity between two of objective functions of two given runs
 */
extern
SCIP_Real SCIPreoptGetSim(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   run1,
   int                   run2
   );

/*
 * returns the best solution of the last run
 */
extern
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/*
 * returns the coefficent of variable with index @param idx in run @param run
 */
extern
SCIP_Real SCIPreoptGetObjCoef(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   run,
   int                   idx
   );

/* checks the changes of the objective coefficients */
extern
void SCIPreoptGetVarCoefChg(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   varidx,
   SCIP_Bool*            negated,
   SCIP_Bool*            entering,
   SCIP_Bool*            leaving
   );

/*
 * print the optimal solutions of all previous runs
 */
extern
SCIP_RETCODE SCIPreoptPrintOptSols(
   SCIP*                 scip,
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/*
 * return all optimal solutions of the previous runs
 * depending on the current stage the method copies the solutions into
 * the origprimal or primal space. That means, all solutions need to be
 * freed before starting a new iteration!!!
 */
extern
SCIP_RETCODE SCIPreoptGetOptSols(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SOL**            sols,
   int*                  nsols
   );

/* reset marks of stored solutions to not updated */
extern
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/* returns the number of stored nodes */
extern
int SCIPreoptGetNNodes(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node
   );

/* returns the number of child nodes of @param node */
int SCIPreoptGetNChildren(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node
   );

/*
 *  Save information if infeasible nodes
 */
extern
SCIP_RETCODE SCIPreoptAddInfNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node
   );

/**
 * check the reason for cut off a node and if necessary store the node
 */
extern
SCIP_RETCODE SCIPreoptCheckCutoff(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   SCIP_NODE*            node,
   SCIP_EVENT*           event,
   SCIP_LPSOLSTAT        lpsolstat,
   SCIP_Bool             isrootnode,
   SCIP_Bool             isfocusnode,
   SCIP_Real             lowerbound,
   int                   effectiverootdepth
   );

/** store bound changes based on dual information */
extern
SCIP_RETCODE SCIPreoptAddDualBndchg(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   SCIP_NODE*            node,
   SCIP_VAR*             var,
   SCIP_Real             newval,
   SCIP_Real             oldval
   );

/* returns the number of bound changes based on dual information */
extern
int SCIPreoptGetNDualBndchs(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node
   );

/* returns the number of leave nodes of the subtree induced
 * by @param node (of the whole tree of node == NULL) */
extern
int SCIPreoptGetNLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node
   );

/* returns the child nodes of @param node that need to be
 * reoptimized next or NULL if @param node is a leaf */
extern
SCIP_RETCODE SCIPreoptGetChildIDs(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   SCIP_NODE*            node,
   int*                  childs,
   int                   mem,
   int*                  nchilds
   );

/* returns all leave if of the subtree induced by @param node */
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

/* returns the time needed to store the nodes */
extern
SCIP_Real SCIPreoptGetSavingtime(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/* store a global constraint that should be added at the beginning of the next iteration */
extern
SCIP_RETCODE SCIPreoptAddGlbCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   LOGICORDATA*          consdata,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* add the stored constraints globally to the problem */
extern
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

extern
SCIP_RETCODE SCIPreoptAddGlbSolCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SOL*             sol,
   SCIP_VAR**            vars,
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   int                   nvars
   );

extern
SCIP_RETCODE SCIPreoptGetSolveLP(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node,
   SCIP_Bool*            solvelp
   );

/* returns the reopttype of a node stored at ID id */
extern
SCIP_REOPTTYPE SCIPreoptnodeGetType(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id
   );

/* reoptimize the node stored at ID id */
extern
SCIP_RETCODE SCIPreoptApply(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP*                 scip,
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,
   SCIP_PROB*            transprob,
   SCIP_PROB*            origprob,
   SCIP_TREE*            tree,
   SCIP_LP*              lp,
   SCIP_BRANCHCAND*      branchcand,
   SCIP_EVENTQUEUE*      eventqueue,
   SCIP_NODE*            node_fix,
   SCIP_NODE*            node_cons,
   int                   id,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* delete a node stored in the reopttree */
extern
SCIP_RETCODE SCIPreopttreeDeleteNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* replace the node stored at ID id by its child nodes */
extern
SCIP_RETCODE SCIPreoptShrinkNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   int                   id
   );

/* return the branching path stored at ID id */
extern
void SCIPreoptnodeGetPath(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       boundtypes,
   int                   mem,
   int*                  nvars,
   int*                  nafterdualvars
   );

/* returns all added constraints at ID id */
extern
void SCIPreoptnodeGetConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id,
   SCIP_VAR***           vars,
   SCIP_Real**           vals,
   int                   mem,
   int*                  nconss,
   int*                  nvars
   );

/* returns the number of bound changes at the node
 * stored at ID id */
extern
int SCIPreoptnodeGetNConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id
   );

/* returns the number of bound changes based on primal information including bound
 * changes directly after the first bound change based on dual information at the node
 * stored at ID id */
extern
int SCIPreoptnodeGetIdNVars(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id
   );

/* reset the stored information abound bound changes based on dual information */
extern
void SCIPreoptResetDualcons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* split the root node and move all children to one of the two resulting nodes */
extern
SCIP_RETCODE SCIPreoptSplitRoot(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

extern
void SCIPreoptCreateSplitCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id,
   LOGICORDATA*          consdata
   );

/* returns if a node need to be split because some bound changes
 * were based on dual information */
extern
SCIP_Bool SCIPreoptSplitNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_NODE*            node
   );

/* calculates a local similarity of a given node and returns if the subproblem
 * should be solved from scratch */
extern
SCIP_RETCODE SCIPreoptCheckLocalRestart(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   SCIP_NODE*            node,
   SCIP_VAR**            origvars,
   int                   norigvars,
   SCIP_Bool*            localrestart
);

/* reset the complete tree an set a new search frontier */
SCIP_RETCODE SCIPreoptApplyCompression(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   SCIP_REOPTNODE**      representatives,
   int                   nrepresentatives,
   SCIP_Bool*            success
   );

/* add the node @param node to the reopttree */
extern
SCIP_RETCODE SCIPreoptAddNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                   /**< block memory */
   SCIP_NODE*            node,
   SCIP_REOPTTYPE        reopttype,
   SCIP_LPSOLSTAT        lpsolstat,
   SCIP_Bool             saveafterduals,
   SCIP_Bool             isrootnode,
   SCIP_Bool             isfocusnode,
   SCIP_Real             lowerbound
   );

/* add all unprocessed nodes to the reoptimization tree */
extern
SCIP_RETCODE SCIPreoptSaveOpenNodes(
   SCIP_REOPT*           reopt,                   /* reoptimization data structure */
   SCIP_SET*             set,                     /* global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /* block memory */
   SCIP_NODE**           leaves,                  /* array of open leave nodes */
   int                   nleaves,                 /* number of open leave nodes */
   SCIP_NODE**           childs,                  /* array of open children nodes */
   int                   nchilds,                 /* number of open leave nodes */
   SCIP_NODE**           siblings,                /* array of open sibling nodes */
   int                   nsiblings                /* number of open leave nodes */
   );

/* returns the number of restarts */
extern
int SCIPreoptGetNRestarts(
   SCIP_REOPT*           reopt                   /**< reoptimization data */
   );

/*
 * return the last lower bound of the stored node
 */
extern
SCIP_Real SCIPreoptGetNodeLb(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   id                       /**< id of the node */
   );

/*
 * initialize an empty node
 */
extern
void SCIPreoptnodeCreateClear(
   SCIP_REOPTNODE*       reoptnode                /**< node of the reopttree */
   );

/*
 * clear the node
 */
extern
SCIP_RETCODE SCIPreoptnodeClear(
   SCIP_REOPTNODE**      reoptnode,               /**< pointer of reoptnode */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* add a variable to a given reoptnode */
extern
SCIP_RETCODE SCIPreoptnodeAddVar(
   SCIP_REOPTNODE*       reoptnode,
   SCIP_VAR*             var,
   SCIP_Real             val,
   SCIP_BOUNDTYPE        bound,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* add a constraint to a given reoptnode */
extern
SCIP_RETCODE SCIPreoptnodeAddCons(
   SCIP_REOPTNODE*       reoptnode,
   SCIP_VAR**            consvars,
   SCIP_Real*            consvals,
   int                   nvars,
   REOPT_CONSTYPE        constype,
   BMS_BLKMEM*           blkmem                   /**< block memory */
   );

/* returns the number of bound changes stored in the reopttree at ID id*/
extern
int SCIPreoptnodeGetNVars(
   SCIP_REOPTNODE*       reoptnode
   );

/* set the parent id */
extern
void SCIPreoptnodeSetParentID(
   SCIP_REOPTNODE*       reoptnode,
   int                   parentid
   );

#ifdef __cplusplus
}
#endif

#endif
