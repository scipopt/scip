/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_stat.h,v 1.20 2004/09/07 18:22:21 bzfpfend Exp $"

/**@file   struct_stat.h
 * @brief  datastructures for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_STAT_H__
#define __STRUCT_STAT_H__


#include "def.h"
#include "type_stat.h"
#include "type_clock.h"
#include "type_vbc.h"
#include "type_history.h"


/** problem and runtime specific statistics */
struct Stat
{
   Longint          nlpiterations;      /**< number of simplex iterations (primal + dual) */
   Longint          nprimallpiterations;/**< number of iterations in primal simplex */
   Longint          nduallpiterations;  /**< number of iterations in dual simplex */
   Longint          nnodelpiterations;  /**< number of iterations for solving node relaxations */
   Longint          ndivinglpiterations;/**< number of iterations in diving */
   Longint          nsblpiterations;    /**< number of simplex iterations used in strong branching */
   Longint          nconflictlpiterations;/**< number of simplex iterations used in conflict analysis */
   Longint          nredcoststrcalls;   /**< number of times, reduced cost strengthening was called */
   Longint          nredcoststrfound;   /**< number of reduced cost strengthenings found */
   Longint          nnodes;             /**< number of nodes processed in current run (including active node) */
   Longint          ntotalnodes;        /**< total number of nodes processed in all runs (including active node) */
   Longint          ncreatednodes;      /**< number of nodes created */
   Longint          nbacktracks;        /**< number of times, the new node was chosen from the leaves queue */
   Longint          nlpsolsfound;       /**< number of CIP-feasible LP solutions found so far */
   Longint          npssolsfound;       /**< number of CIP-feasible pseudo solutions found so far */
   Longint          lastdispnode;       /**< last node for which an information line was displayed */
   Longint          lastdivenode;       /**< last node where LP diving was applied */
   Longint          domchgcount;        /**< internal counter, where all domain changes are counted */
   Longint          nrootboundchgs;     /**< total number of bound changes generated in the root node */
   Longint          nboundchgs;         /**< total number of bound changes generated in the tree */
   Longint          nholechgs;          /**< total number of hole changes generated in the tree */
   Real             rootlowerbound;     /**< lower bound of root node */
   CLOCK*           solvingtime;        /**< total time used for solving (including presolving) the current problem */
   CLOCK*           presolvingtime;     /**< total time used for presolving the current problem */
   CLOCK*           primallptime;       /**< primal LP solution time */
   CLOCK*           duallptime;         /**< dual LP solution time */
   CLOCK*           divinglptime;       /**< diving LP solution time (primal + dual) */
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
   BRANCHDIR        lastbranchdir;      /**< direction of the last branching */
   int              nruns;              /**< number of branch and bound runs on current problem, including current run */
   int              nrootboundchgsrun;  /**< total number of bound changes generated in the root node of current run */
   int              nvaridx;            /**< number of used variable indices */
   int              ncolidx;            /**< number of used column indices */
   int              nrowidx;            /**< number of used row indices */
   int              marked_nvaridx;     /**< number of used variable indices before solving started */
   int              marked_ncolidx;     /**< number of used column indices before solving started */
   int              marked_nrowidx;     /**< number of used row indices before solving started */
   int              lpcount;            /**< internal counter, where all SCIPlpSolve() calls are counted */
   int              nlps;               /**< number of LPs solved (primal + dual) with at least 1 iteration */
   int              nprimallps;         /**< number of primal LPs solved */
   int              nduallps;           /**< number of dual LPs solved */
   int              nnodelps;           /**< number of LPs solved for node relaxations */
   int              ndivinglps;         /**< number of LPs solved during diving */
   int              nstrongbranchs;     /**< number of strong branching calls */
   int              nconflictlps;       /**< number of LPs solved during conflict analysis */
   int              npricerounds;       /**< number of pricing rounds performed in current node */
   int              nseparounds;        /**< number of separation rounds performed in current node */
   int              ndisplines;         /**< number of displayed information lines */
   int              maxdepth;           /**< maximal depth of all processed nodes in current run */
   int              maxtotaldepth;      /**< maximal depth of all processed nodes over all runs */
   int              plungedepth;        /**< current plunging depth (successive times, a child was selected as next node) */
   int              nactiveconss;       /**< total number of currently active constraints */
   int              nenabledconss;      /**< total number of currently enabled constraints */
   Bool             memsavemode;        /**< should algorithms be switched to memory saving mode? */
};


#endif
