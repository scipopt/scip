/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_stat.h,v 1.45 2006/03/16 14:43:07 bzfpfend Exp $"

/**@file   struct_stat.h
 * @brief  datastructures for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_STAT_H__
#define __SCIP_STRUCT_STAT_H__


#include "scip/def.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_vbc.h"
#include "scip/type_history.h"


/** problem and runtime specific statistics */
struct SCIP_Stat
{
   SCIP_Longint          nlpiterations;      /**< total number of LP iterations */
   SCIP_Longint          nprimallpiterations;/**< number of iterations in primal simplex */
   SCIP_Longint          nduallpiterations;  /**< number of iterations in dual simplex */
   SCIP_Longint          nbarrierlpiterations;/**< number of iterations in barrier algorithm */
   SCIP_Longint          nprimalresolvelpiterations;  /**< number of primal LP iterations with advanced start basis */
   SCIP_Longint          ndualresolvelpiterations;    /**< number of dual LP iterations with advanced start basis */
   SCIP_Longint          nnodelpiterations;  /**< number of iterations for totally solving node relaxations */
   SCIP_Longint          ninitlpiterations;  /**< number of iterations for solving nodes' initial relaxations */
   SCIP_Longint          ndivinglpiterations;/**< number of iterations in diving and probing */
   SCIP_Longint          nsblpiterations;    /**< number of simplex iterations used in strong branching */
   SCIP_Longint          nrootsblpiterations;/**< number of simplex iterations used in strong branching at the root node */
   SCIP_Longint          nconflictlpiterations;/**< number of simplex iterations used in conflict analysis */
   SCIP_Longint          nnodes;             /**< number of nodes processed in current run (including focus node) */
   SCIP_Longint          ntotalnodes;        /**< total number of nodes processed in all runs (including focus node) */
   SCIP_Longint          ncreatednodes;      /**< total number of nodes created */
   SCIP_Longint          ncreatednodesrun;   /**< number of nodes created in current run */
   SCIP_Longint          nactivatednodes;    /**< number of times, a node got activated in current run */
   SCIP_Longint          ndeactivatednodes;  /**< number of times, a node got deactivated in current run */
   SCIP_Longint          nbacktracks;        /**< number of times, the new node was chosen from the leaves queue */
   SCIP_Longint          ndelayedcutoffs;    /**< number of times, the selected node was from a cut off subtree */
   SCIP_Longint          nreprops;           /**< number of times, a solved node is repropagated again */
   SCIP_Longint          nrepropboundchgs;   /**< number of bound changes generated in repropagating nodes */
   SCIP_Longint          nrepropcutoffs;     /**< number of times, a repropagated node was cut off */
   SCIP_Longint          nlpsolsfound;       /**< number of CIP-feasible LP solutions found so far */
   SCIP_Longint          npssolsfound;       /**< number of CIP-feasible pseudo solutions found so far */
   SCIP_Longint          lastdispnode;       /**< last node for which an information line was displayed */
   SCIP_Longint          lastdivenode;       /**< last node where LP diving was applied */
   SCIP_Longint          lastconflictnode;   /**< last node where conflict analysis was applied */
   SCIP_Longint          bestsolnode;        /**< node number where the last incumbent solution was found */
   SCIP_Longint          domchgcount;        /**< internal counter, where all domain changes are counted */
   SCIP_Longint          nboundchgs;         /**< total number of bound changes generated in the tree */
   SCIP_Longint          nholechgs;          /**< total number of hole changes generated in the tree */
   SCIP_Real             rootlowerbound;     /**< lower bound of root node */
   SCIP_Real             conflictscoreweight;/**< current weight to use for updating conflict scores in history */
   SCIP_CLOCK*           solvingtime;        /**< total time used for solving (including presolving) the current problem */
   SCIP_CLOCK*           presolvingtime;     /**< total time used for presolving the current problem */
   SCIP_CLOCK*           primallptime;       /**< primal LP solution time */
   SCIP_CLOCK*           duallptime;         /**< dual LP solution time */
   SCIP_CLOCK*           barrierlptime;      /**< barrier LP solution time */
   SCIP_CLOCK*           divinglptime;       /**< diving and probing LP solution time */
   SCIP_CLOCK*           strongbranchtime;   /**< strong branching time */
   SCIP_CLOCK*           conflictlptime;     /**< conflict analysis LP solution time */
   SCIP_CLOCK*           lpsoltime;          /**< time needed for storing feasible LP solutions */
   SCIP_CLOCK*           pseudosoltime;      /**< time needed for storing feasible pseudo solutions */
   SCIP_CLOCK*           nodeactivationtime; /**< time needed for path switching and activating nodes */
   SCIP_HISTORY*         glbhistory;         /**< global history information over all variables */
   SCIP_HISTORY*         glbhistorycrun;     /**< global history information over all variables for current run */
   SCIP_VAR*             lastbranchvar;      /**< last variable, that was branched on */
   SCIP_VBC*             vbc;                /**< VBC Tool information */
   SCIP_STATUS           status;             /**< SCIP solving status */
   SCIP_BRANCHDIR        lastbranchdir;      /**< direction of the last branching */
   int                   nruns;              /**< number of branch and bound runs on current problem, including current run */
   int                   nrootboundchgs;     /**< total number of bound changes generated in the root node */
   int                   nrootboundchgsrun;  /**< total number of bound changes generated in the root node of current run */
   int                   nrootintfixings;    /**< total number of global fixings of integer variables */
   int                   nrootintfixingsrun; /**< total number of global fixings of integer variables of current run */
   int                   nvaridx;            /**< number of used variable indices */
   int                   ncolidx;            /**< number of used column indices */
   int                   nrowidx;            /**< number of used row indices */
   int                   marked_nvaridx;     /**< number of used variable indices before solving started */
   int                   marked_ncolidx;     /**< number of used column indices before solving started */
   int                   marked_nrowidx;     /**< number of used row indices before solving started */
   int                   lpcount;            /**< internal counter, where all simplex calls are counted */
   int                   nlps;               /**< total number of LPs solved with at least 1 iteration */
   int                   nprimallps;         /**< number of primal LPs solved */
   int                   nduallps;           /**< number of dual LPs solved */
   int                   nbarrierlps;        /**< number of barrier LPs solved */
   int                   nprimalresolvelps;  /**< number of primal LPs solved with advanced start basis and at least 1 iteration */
   int                   ndualresolvelps;    /**< number of dual LPs solved with advanced start basis and at least 1 iteration */
   int                   nnodelps;           /**< number of LPs solved for node relaxations */
   int                   ninitlps;           /**< number of LPs solved for nodes' initial relaxations */
   int                   ndivinglps;         /**< number of LPs solved during diving and probing */
   int                   nstrongbranchs;     /**< number of strong branching calls */
   int                   nrootstrongbranchs; /**< number of strong branching calls at the root node */
   int                   nconflictlps;       /**< number of LPs solved during conflict analysis */
   int                   npricerounds;       /**< number of pricing rounds performed in current node */
   int                   nseparounds;        /**< number of separation rounds performed in current node */
   int                   ndisplines;         /**< number of displayed information lines */
   int                   maxdepth;           /**< maximal depth of all processed nodes in current run */
   int                   maxtotaldepth;      /**< maximal depth of all processed nodes over all runs */
   int                   plungedepth;        /**< current plunging depth (successive times, a child was selected as next node) */
   int                   nactiveconss;       /**< total number of currently active constraints */
   int                   nenabledconss;      /**< total number of currently enabled constraints */
   int                   nimplications;      /**< total number of implications stored in the implication graph */
   int                   npresolrounds;      /**< number of presolving rounds in current run */
   int                   npresolfixedvars;   /**< number of presolving fixings in current run */
   int                   npresolaggrvars;    /**< number of presolving aggregations in current run */
   int                   npresolchgvartypes; /**< number of presolving variable type changes in current run */
   int                   npresolchgbds;      /**< number of presolving bound changes in current run */
   int                   npresoladdholes;    /**< number of presolving hole additions in current run */
   int                   npresoldelconss;    /**< number of presolving constraint deletions in current run */
   int                   npresolupgdconss;   /**< number of presolving constraint upgrades in current run */
   int                   npresolchgcoefs;    /**< number of presolving coefficient changes in current run */
   int                   npresolchgsides;    /**< number of presolving side changes in current run */
   int                   solindex;           /**< consecutively numbered solution index */
   SCIP_Bool             memsavemode;        /**< should algorithms be switched to memory saving mode? */
   SCIP_Bool             userinterrupt;      /**< has the user asked to interrupt the solving process? */
};


#endif
