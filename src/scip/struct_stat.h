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

/**@file   struct_stat.h
 * @brief  datastructures for problem statistics
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_STAT_H__
#define __SCIP_STRUCT_STAT_H__


#include "scip/def.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_vbc.h"
#include "scip/type_history.h"

#ifdef __cplusplus
extern "C" {
#endif

/** problem and runtime specific statistics */
struct SCIP_Stat
{
   SCIP_Longint          nlpiterations;      /**< total number of LP iterations */
   SCIP_Longint          nrootlpiterations;  /**< total number of LP iterations in root node */
   SCIP_Longint          nprimallpiterations;/**< number of iterations in primal simplex */
   SCIP_Longint          nduallpiterations;  /**< number of iterations in dual simplex */
   SCIP_Longint          nlexduallpiterations;/**< number of iterations in lexicographic dual simplex */
   SCIP_Longint          nbarrierlpiterations;/**< number of iterations in barrier algorithm */
   SCIP_Longint          nprimalresolvelpiterations;  /**< number of primal LP iterations with advanced start basis */
   SCIP_Longint          ndualresolvelpiterations;    /**< number of dual LP iterations with advanced start basis */
   SCIP_Longint          nlexdualresolvelpiterations; /**< number of lexicographic dual LP iterations with advanced start basis */
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
   SCIP_Longint          nprobboundchgs;     /**< total number of bound changes generated in the tree during probing */
   SCIP_Longint          nprobholechgs;      /**< total number of hole changes generated in the tree  during probing */
   SCIP_Longint          nnodesbeforefirst;  /**< number of nodes before first primal solution */   
   SCIP_Real             rootlowerbound;     /**< lower bound of root node */
   SCIP_Real             vsidsweight;        /**< current weight to use for updating VSIDS in history */
   SCIP_Real             firstprimalbound;   /**< objective value of first primal solution */
   SCIP_Real             firstprimaltime;    /**< time (in seconds) needed for first primal solution */
   SCIP_Real             primalzeroittime;   /**< time used in primal simplex calls without iterations */
   SCIP_Real             dualzeroittime;     /**< time used in dual simplex calls without iterations */
   SCIP_Real             barrierzeroittime;  /**< time used in barrier calls without iterations */
   SCIP_Real             maxcopytime;        /**< maxmimal time needed for copying a problem */
   SCIP_Real             mincopytime;        /**< minimal time needed for copying a problem */
   SCIP_CLOCK*           solvingtime;        /**< total time used for solving (including presolving) the current problem */
   SCIP_CLOCK*           presolvingtime;     /**< total time used for presolving the current problem */
   SCIP_CLOCK*           primallptime;       /**< primal LP solution time */
   SCIP_CLOCK*           duallptime;         /**< dual LP solution time */
   SCIP_CLOCK*           lexduallptime;      /**< lexicographic dual LP solution time */
   SCIP_CLOCK*           barrierlptime;      /**< barrier LP solution time */
   SCIP_CLOCK*           divinglptime;       /**< diving and probing LP solution time */
   SCIP_CLOCK*           strongbranchtime;   /**< strong branching time */
   SCIP_CLOCK*           conflictlptime;     /**< conflict analysis LP solution time */
   SCIP_CLOCK*           lpsoltime;          /**< time needed for storing feasible LP solutions */
   SCIP_CLOCK*           pseudosoltime;      /**< time needed for storing feasible pseudo solutions */
   SCIP_CLOCK*           nodeactivationtime; /**< time needed for path switching and activating nodes */
   SCIP_CLOCK*           nlpsoltime;         /**< time needed for solving NLPs */
   SCIP_CLOCK*           copyclock;          /**< time needed for copying problems */
   SCIP_HISTORY*         glbhistory;         /**< global history information over all variables */
   SCIP_HISTORY*         glbhistorycrun;     /**< global history information over all variables for current run */
   SCIP_VAR*             lastbranchvar;      /**< last variable, that was branched on */
   SCIP_VBC*             vbc;                /**< VBC Tool information */
   SCIP_HEUR*            firstprimalheur;    /**< heuristic which found the first primal solution */     
   SCIP_STATUS           status;             /**< SCIP solving status */
   SCIP_BRANCHDIR        lastbranchdir;      /**< direction of the last branching */
   SCIP_Longint          lpcount;            /**< internal counter, where all lp calls are counted; this includes the restored lps after diving and probing */
   SCIP_Longint          nlps;               /**< total number of LPs solved with at least 1 iteration */
   SCIP_Longint          nrootlps;           /**< number of LPs solved at the root node with at least 1 iteration */
   SCIP_Longint          nprimallps;         /**< number of primal LPs solved with at least 1 iteration */
   SCIP_Longint          nprimalzeroitlps;   /**< number of primal LPs with 0 iterations */
   SCIP_Longint          nduallps;           /**< number of dual LPs solved with at least 1 iteration */
   SCIP_Longint          ndualzeroitlps;     /**< number of dual LPs with 0 iterations */
   SCIP_Longint          nlexduallps;        /**< number of lexicographic dual LPs solved */
   SCIP_Longint          nbarrierlps;        /**< number of barrier LPs solved with at least 1 iteration */
   SCIP_Longint          nbarrierzeroitlps;  /**< number of barrier LPs with 1 iteration */
   SCIP_Longint          nprimalresolvelps;  /**< number of primal LPs solved with advanced start basis and at least 1 iteration */
   SCIP_Longint          ndualresolvelps;    /**< number of dual LPs solved with advanced start basis and at least 1 iteration */
   SCIP_Longint          nlexdualresolvelps; /**< number of lexicographic dual LPs solved with advanced start basis and at least 1 iteration */
   SCIP_Longint          nnodelps;           /**< number of LPs solved for node relaxations */
   SCIP_Longint          ninitlps;           /**< number of LPs solved for nodes' initial relaxations */
   SCIP_Longint          ndivinglps;         /**< number of LPs solved during diving and probing */
   SCIP_Longint          nstrongbranchs;     /**< number of strong branching calls */
   SCIP_Longint          nrootstrongbranchs; /**< number of strong branching calls at the root node */
   SCIP_Longint          nconflictlps;       /**< number of LPs solved during conflict analysis */
   SCIP_Longint          nnlps;              /**< number of NLPs solved */
   int                   subscipdepth;       /**< depth of current scip instance (increased by each copy call) */
   int                   nruns;              /**< number of branch and bound runs on current problem, including current run */
   int                   nconfrestarts;      /**< number of restarts performed due to conflict analysis */
   int                   nrootboundchgs;     /**< total number of bound changes generated in the root node */
   int                   nrootboundchgsrun;  /**< total number of bound changes generated in the root node of current run */
   int                   nrootintfixings;    /**< total number of global fixings of integer variables */
   int                   nrootintfixingsrun; /**< total number of global fixings of integer variables of current run */
   int                   prevrunnvars;       /**< number of variables in the previous run */
   int                   nvaridx;            /**< number of used variable indices */
   int                   ncolidx;            /**< number of used column indices */
   int                   nrowidx;            /**< number of used row indices */
   int                   marked_nvaridx;     /**< number of used variable indices before solving started */
   int                   marked_ncolidx;     /**< number of used column indices before solving started */
   int                   marked_nrowidx;     /**< number of used row indices before solving started */
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
   int                   npresoladdconss;    /**< number of presolving constraint additions in current run */
   int                   npresolupgdconss;   /**< number of presolving constraint upgrades in current run */
   int                   npresolchgcoefs;    /**< number of presolving coefficient changes in current run */
   int                   npresolchgsides;    /**< number of presolving side changes in current run */
   int                   lastnpresolfixedvars;/**< number of presolving fixings in current run */
   int                   lastnpresolaggrvars;/**< number of presolving aggregations in current run */
   int                   lastnpresolchgvartypes;/**< number of presolving variable type changes in current run */
   int                   lastnpresolchgbds;  /**< number of presolving bound changes in current run */
   int                   lastnpresoladdholes;/**< number of presolving hole additions in current run */
   int                   lastnpresoldelconss;/**< number of presolving constraint deletions in current run */
   int                   lastnpresoladdconss;/**< number of presolving constraint additions in current run */
   int                   lastnpresolupgdconss;/**< number of presolving constraint upgrades in current run */
   int                   lastnpresolchgcoefs;/**< number of presolving coefficient changes in current run */
   int                   lastnpresolchgsides;/**< number of presolving side changes in current run */
   int                   solindex;           /**< consecutively numbered solution index */
   int                   nrunsbeforefirst;   /**< number of runs until first primal solution */
   int                   firstprimaldepth;   /**< depth in which first primal solution was found */
   int                   ncopies;            /**< counter how often SCIPcopy() was performed */
   SCIP_Bool             memsavemode;        /**< should algorithms be switched to memory saving mode? */
   SCIP_Bool             userinterrupt;      /**< has the user asked to interrupt the solving process? */
   SCIP_Bool             userrestart;        /**< has the user asked to restart the solving process? */
   SCIP_Bool             inrestart;          /**< are we currently restarting the system? */
   SCIP_Bool             collectvarhistory;  /**< should variable history statistics be collected */
};

#ifdef __cplusplus
}
#endif

#endif
