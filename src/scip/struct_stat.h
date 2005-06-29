/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_stat.h,v 1.36 2005/06/29 11:08:07 bzfpfend Exp $"

/**@file   struct_stat.h
 * @brief  datastructures for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_STAT_H__
#define __STRUCT_STAT_H__


#include "scip/def.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_vbc.h"
#include "scip/type_history.h"


/** problem and runtime specific statistics */
struct Stat
{
   Longint          nlpiterations;      /**< total number of LP iterations */
   Longint          nprimallpiterations;/**< number of iterations in primal simplex */
   Longint          nduallpiterations;  /**< number of iterations in dual simplex */
   Longint          nbarrierlpiterations;/**< number of iterations in barrier algorithm */
   Longint          nprimalresolvelpiterations;  /**< number of primal LP iterations with advanced start basis */
   Longint          ndualresolvelpiterations;    /**< number of dual LP iterations with advanced start basis */
   Longint          nnodelpiterations;  /**< number of iterations for totally solving node relaxations */
   Longint          ninitlpiterations;  /**< number of iterations for solving nodes' initial relaxations */
   Longint          ndivinglpiterations;/**< number of iterations in diving and probing */
   Longint          nsblpiterations;    /**< number of simplex iterations used in strong branching */
   Longint          nrootsblpiterations;/**< number of simplex iterations used in strong branching at the root node */
   Longint          nconflictlpiterations;/**< number of simplex iterations used in conflict analysis */
   Longint          nredcoststrcalls;   /**< number of times, reduced cost strengthening was called */
   Longint          nredcoststrfound;   /**< number of reduced cost strengthenings found */
   Longint          nnodes;             /**< number of nodes processed in current run (including focus node) */
   Longint          ntotalnodes;        /**< total number of nodes processed in all runs (including focus node) */
   Longint          ncreatednodes;      /**< total number of nodes created */
   Longint          ncreatednodesrun;   /**< number of nodes created in current run */
   Longint          nactivatednodes;    /**< number of times, a node got activated in current run */
   Longint          ndeactivatednodes;  /**< number of times, a node got deactivated in current run */
   Longint          nbacktracks;        /**< number of times, the new node was chosen from the leaves queue */
   Longint          ndelayedcutoffs;    /**< number of times, the selected node was from a cut off subtree */
   Longint          nreprops;           /**< number of times, a solved node is repropagated again */
   Longint          nlpsolsfound;       /**< number of CIP-feasible LP solutions found so far */
   Longint          npssolsfound;       /**< number of CIP-feasible pseudo solutions found so far */
   Longint          lastdispnode;       /**< last node for which an information line was displayed */
   Longint          lastdivenode;       /**< last node where LP diving was applied */
   Longint          domchgcount;        /**< internal counter, where all domain changes are counted */
   Longint          nrootboundchgs;     /**< total number of bound changes generated in the root node */
   Longint          nrepropboundchgs;   /**< total number of bound changes generated in repropagating nodes */
   Longint          nboundchgs;         /**< total number of bound changes generated in the tree */
   Longint          nholechgs;          /**< total number of hole changes generated in the tree */
   Real             rootlowerbound;     /**< lower bound of root node */
   CLOCK*           solvingtime;        /**< total time used for solving (including presolving) the current problem */
   CLOCK*           presolvingtime;     /**< total time used for presolving the current problem */
   CLOCK*           primallptime;       /**< primal LP solution time */
   CLOCK*           duallptime;         /**< dual LP solution time */
   CLOCK*           barrierlptime;      /**< barrier LP solution time */
   CLOCK*           divinglptime;       /**< diving and probing LP solution time */
   CLOCK*           strongbranchtime;   /**< strong branching time */
   CLOCK*           conflictlptime;     /**< conflict analysis LP solution time */
   CLOCK*           lpsoltime;          /**< time needed for storing feasible LP solutions */
   CLOCK*           pseudosoltime;      /**< time needed for storing feasible pseudo solutions */
   CLOCK*           redcoststrtime;     /**< time needed for reduced cost strengthening */
   CLOCK*           nodeactivationtime; /**< time needed for path switching and activating nodes */
   HISTORY*         glbhistory;         /**< global history information over all variables */
   HISTORY*         glbhistorycrun;     /**< global history information over all variables for current run */
   VAR*             lastbranchvar;      /**< last variable, that was branched on */
   VBC*             vbc;                /**< VBC Tool information */
   STATUS           status;             /**< SCIP solving status */
   BRANCHDIR        lastbranchdir;      /**< direction of the last branching */
   int              nruns;              /**< number of branch and bound runs on current problem, including current run */
   int              nrootboundchgsrun;  /**< total number of bound changes generated in the root node of current run */
   int              nvaridx;            /**< number of used variable indices */
   int              ncolidx;            /**< number of used column indices */
   int              nrowidx;            /**< number of used row indices */
   int              marked_nvaridx;     /**< number of used variable indices before solving started */
   int              marked_ncolidx;     /**< number of used column indices before solving started */
   int              marked_nrowidx;     /**< number of used row indices before solving started */
   int              lpcount;            /**< internal counter, where all simplex calls are counted */
   int              nlps;               /**< total number of LPs solved with at least 1 iteration */
   int              nprimallps;         /**< number of primal LPs solved */
   int              nduallps;           /**< number of dual LPs solved */
   int              nbarrierlps;        /**< number of barrier LPs solved */
   int              nprimalresolvelps;  /**< number of primal LPs solved with advanced start basis and at least 1 iteration */
   int              ndualresolvelps;    /**< number of dual LPs solved with advanced start basis and at least 1 iteration */
   int              nnodelps;           /**< number of LPs solved for node relaxations */
   int              ninitlps;           /**< number of LPs solved for nodes' initial relaxations */
   int              ndivinglps;         /**< number of LPs solved during diving and probing */
   int              nstrongbranchs;     /**< number of strong branching calls */
   int              nrootstrongbranchs; /**< number of strong branching calls at the root node */
   int              nconflictlps;       /**< number of LPs solved during conflict analysis */
   int              npricerounds;       /**< number of pricing rounds performed in current node */
   int              nseparounds;        /**< number of separation rounds performed in current node */
   int              ndisplines;         /**< number of displayed information lines */
   int              maxdepth;           /**< maximal depth of all processed nodes in current run */
   int              maxtotaldepth;      /**< maximal depth of all processed nodes over all runs */
   int              plungedepth;        /**< current plunging depth (successive times, a child was selected as next node) */
   int              nactiveconss;       /**< total number of currently active constraints */
   int              nenabledconss;      /**< total number of currently enabled constraints */
   int              nimplications;      /**< total number of implications stored in the implication graph */
   int              npresolrounds;      /**< number of presolving rounds in current run */
   int              npresolfixedvars;   /**< number of presolving fixings in current run */
   int              npresolaggrvars;    /**< number of presolving aggregations in current run */
   int              npresolchgvartypes; /**< number of presolving variable type changes in current run */
   int              npresolchgbds;      /**< number of presolving bound changes in current run */
   int              npresoladdholes;    /**< number of presolving hole additions in current run */
   int              npresoldelconss;    /**< number of presolving constraint deletions in current run */
   int              npresolupgdconss;   /**< number of presolving constraint upgrades in current run */
   int              npresolchgcoefs;    /**< number of presolving coefficient changes in current run */
   int              npresolchgsides;    /**< number of presolving side changes in current run */
   Bool             memsavemode;        /**< should algorithms be switched to memory saving mode? */
};


#endif
