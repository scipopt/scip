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

/**@file   struct_reopt.h
 * @brief  data structures for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_REOPT_H__
#define __SCIP_STRUCT_REOPT_H__


#include "scip/def.h"
#include "scip/type_reopt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** nodes of SCIP_SolTree */
struct SCIP_SolNode
{
   SCIP_SOL*             sol;                /**< the stored solution */
   SCIP_SOLNODE*         father;             /**< pointer to the parent node */
   SCIP_SOLNODE*         rchild;             /**< pointer to the right child node (0-branch) */
   SCIP_SOLNODE*         lchild;             /**< pointer to the left child node (1-branch) */
   SCIP_Bool             updated;            /**< flag if the solution is already updated
                                              *   w.r.t. the new objective function */
};

/** tree for solution */
struct SCIP_SolTree
{
   SCIP_SOLNODE*         root;               /**< root node of the solution tree */
   SCIP_SOLNODE***       sols;               /**< array of arrays of solutions of the reoptimization runs */
   int*                  solssize;           /**< size of sols[x] arrays */
   int*                  nsols;              /**< number of solutions stored in sols[x] array */
};

/** data for constraints to split nodes during reoptimization */
struct LogicOrData
{
   SCIP_VAR**            vars;               /**< array of variables */
   SCIP_Real*            vals;               /**< array of variable bounds */
   REOPT_CONSTYPE        constype;           /**< type of the constraint */
   int                   varssize;           /**< available sitze in the arrays */
   int                   nvars;              /**< number of enties in the arrays */
};

/** nodes of SCIP_ReoptTree */
struct SCIP_ReoptNode
{
   /* data to store the branching path */
   SCIP_VAR**            vars;                    /**< variables along the branching path up to the next stored node */
   SCIP_Real*            varbounds;               /**< bounds along the branching path up to the next stored node */
   SCIP_BOUNDTYPE*       varboundtypes;           /**< boundtypes along the branching path up to the next stored node */
   int                   nvars;                   /**< number of branching decisions up to the next stored node */
   int                   varssize;                /**< size of allocated memory */

   SCIP_VAR**            afterdualvars;           /**< variables along the branching path after the first decision based on dual information */
   SCIP_Real*            afterdualvarbounds;      /**< bounds along the branching path after the first decision based on dual information */
   SCIP_BOUNDTYPE*       afterdualvarboundtypes;  /**< boundtypes along the branching path after the first dual information */
   int                   nafterdualvars;          /**< number of branching decisions after the first dual information */
   int                   afterdualvarssize;       /**< size of allocated memory */

   SCIP_Bool             dualfixing;              /**< flag whether bound changes based on dual decisions exist */
   SCIP_Real             lowerbound;              /**< the last lowerbound of this node in the previous iteration */

   /* other information: child nodes, added constraints, ... */
   int*                  childids;                /**< array of child nodes that need to be reoptimized */
   int                   nchilds;                 /**< number of child nodes */
   int                   allocchildmem;           /**< allocated memory for child nodes */
   LOGICORDATA**         conss;                   /**< array of constraints added to the node, i.e., logic-or constraints */
   int                   nconss;                  /**< number of added constraints */
   int                   consssize;               /**< allocated memory for constraints */
   LOGICORDATA*          dualconscur;             /**< dual constraint that need to be added the current round */
   LOGICORDATA*          dualconsnex;             /**< dual constraint that need to be added the next round */

   unsigned int          parentID:29;             /**< id of the stored parent node */
   unsigned int          reopttype:3;             /**< reason for storing the node */
};

/* tree to store the current search tree */
struct SCIP_ReoptTree
{
   SCIP_Bool             initialized;             /**< is the data structure initialized? */

   SCIP_REOPTNODE**      reoptnodes;              /**< array of SCIP_REOPTNODE */
   SCIP_QUEUE*           openids;                 /**< queue of open positions in the reoptnodes array */
   int                   reoptnodessize;          /**< size of allocated memory for the reoptnodes array and the openid queue */
   int                   nreoptnodes;             /**< number of saved nodes */

   /* statistics */
   int                   nfeasnodes;              /**< number of feasible nodes in the current run */
   int                   ntotalfeasnodes;         /**< number of feasible nodes over all runs */
   int                   ninfnodes;               /**< number of (LP-)infeasible nodes in the current run */
   int                   ntotalinfnodes;          /**< number of (LP-)infeasible nodes over all runs */
   int                   nprunednodes;            /**< number of pruned nodes in the current run */
   int                   ntotalprunednodes;       /**< number of pruned nodes over all runs */
   int                   nrediednodes;            /**< number of pruned reoptimized nodes in the current run */
   int                   ntotalrediednodes;       /**< number of pruned reoptimized nodes over all runs */
   int                   ninfsubtrees;            /**< number of found infeasible subtrees */
};

/** reoptimization data and solution storage */
struct SCIP_Reopt
{
   int                   run;                     /**< number of the current reoptimization run */
   int                   runsize;                 /**< allocated memory for runs */
   int                   firstobj;                /**< first non empty objective function */
   int                   nobjvars;                /**< number of variables */

   SCIP_Real**           objs;                    /**< list of objective coefficients */
   SCIP_SOL**            prevbestsols;            /**< list of best solutions of all previous rounds */
   SCIP_SOLTREE*         soltree;                 /**< tree to handle all saved solutions */

   SCIP_Real             simtolastobj;            /**< similarity to the last objective function */
   SCIP_Real             simtofirstobj;           /**< similarity to the first objective function */

   SCIP_REOPTTREE*       reopttree;               /**< data structure to store the current reoptimization search tree */

   SCIP_Longint          lastbranched;            /**< number of the last branched node */

   int                   noptsolsbyreoptsol;      /**< number of successive optimal solutions found by heur_reoptsols */

   SCIP_Longint          lastseennode;            /**< node number of the last caught event */

   /* data structure to track decisions based on dual information */
   SCIP_Longint          currentnode;             /**< number of the current node */
   LOGICORDATA*          dualcons;                /**< constraint describing bound changes based on dual information */

   /* global constraints */
   LOGICORDATA**         glbconss;                /**< global constraints that need to be added at the beginning of the next iteration */
   int                   nglbconss;               /**< number of stored global constraints */
   int                   allocmemglbconss;        /**< allocated memory for global constraints */

   /* statistics */
   int                   nlocalrestarts;          /**< number of local restarts in the current iteration */
   int                   ncheckedsols;            /**< number of updated solutions by reoptsols */
   int                   nimprovingsols;          /**< number of improving solutions found by reoptsols */
   int                   nrestarts;               /**< number of restarts */

   /* clocks */
   SCIP_CLOCK*           savingtime;              /**< time needed to store the nodes */
};

#ifdef __cplusplus
}
#endif

#endif
