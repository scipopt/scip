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
#pragma ident "@(#) $Id: set.c,v 1.89 2004/02/04 17:27:42 bzfpfend Exp $"

/**@file   set.c
 * @brief  methods for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "def.h"
#include "set.h"
#include "stat.h"
#include "misc.h"
#include "event.h"
#include "lp.h"
#include "paramset.h"
#include "scip.h"
#include "branch.h"
#include "conflict.h"
#include "cons.h"
#include "disp.h"
#include "heur.h"
#include "nodesel.h"
#include "presol.h"
#include "pricer.h"
#include "reader.h"
#include "sepa.h"



/*
 * Default settings
 */

/* Message Output */

#define SCIP_DEFAULT_VERBLEVEL    SCIP_VERBLEVEL_NORMAL


/* CTRL-C interrupt */

#define SCIP_DEFAULT_CATCHCTRLC       TRUE /**< should the CTRL-C interrupt be catched by SCIP? */


/* Dynamic Memory */

#define SCIP_DEFAULT_MEMLIMIT        1e+20 /**< maximal memory usage in MB */
#define SCIP_DEFAULT_MEMSAVEFAC        0.8 /**< fraction of maximal mem usage when switching to memory saving mode */
#define SCIP_DEFAULT_MEMGROWFAC        1.2 /**< memory growing factor for dynamically allocated arrays */
#define SCIP_DEFAULT_MEMGROWINIT         4 /**< initial size of dynamically allocated arrays */
#define SCIP_DEFAULT_TREEGROWFAC       2.0 /**< memory growing factor for tree array */
#define SCIP_DEFAULT_TREEGROWINIT    65536 /**< initial size of tree array */
#define SCIP_DEFAULT_PATHGROWFAC       2.0 /**< memory growing factor for path array */
#define SCIP_DEFAULT_PATHGROWINIT      256 /**< initial size of path array */


/* Branching */
#define SCIP_DEFAULT_BRANCHSCOREFAC  0.167 /**< branching score factor to weigh downward and upward gain prediction */


/* Presolving */
#define SCIP_DEFAULT_MAXPRESOLROUNDS    -1 /**< maximal number of presolving rounds (-1: unlimited) */
#define SCIP_DEFAULT_PRESOLABORTFAC  1e-04 /**< abort presolve, if l.t. frac of the problem was changed in last round */


/* LP Solving */

#define SCIP_DEFAULT_CHECKLPFEAS      TRUE /**< should LP solutions be checked, resolving LP when numerical troubles occur? */
#define SCIP_DEFAULT_FASTMIP          TRUE /**< should FASTMIP setting of LP solver be used? */
#define SCIP_DEFAULT_SCALING          TRUE /**< should scaling of LP solver be used? */
#define SCIP_DEFAULT_LPSOLVEFREQ         1 /**< frequency for solving LP at the nodes; -1: never; 0: only root LP */
#define SCIP_DEFAULT_LPSOLVEDEPTH       -1 /**< maximal depth for solving LPs (-1: no depth limit) */
#define SCIP_DEFAULT_REDCOSTFREQ         5 /**< frequency for applying reduced cost fixing (-1: never; 0: only root LP) */
#define SCIP_DEFAULT_COLAGELIMIT        10 /**< maximum age a dynamic column can reach before it is deleted from the LP */
#define SCIP_DEFAULT_ROWAGELIMIT        10 /**< maximum age a dynamic row can reach before it is deleted from the LP */


/* Pricing */

#define SCIP_DEFAULT_MAXPRICEVARS      100 /**< maximal number of variables priced in per pricing round */
#define SCIP_DEFAULT_MAXPRICEVARSROOT 2000 /**< maximal number of priced variables at the root node */
#define SCIP_DEFAULT_ABORTPRICEVARSFAC 2.0 /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
#define SCIP_DEFAULT_CLEANUPCOLS     FALSE /**< should new non-basic columns be removed after LP solving? */


/* Cut Separation */

#define SCIP_DEFAULT_MAXSEPACUTS       100 /**< maximal number of cuts separated per separation round */
#define SCIP_DEFAULT_MAXSEPACUTSROOT  2000 /**< maximal separated cuts at the root node */
#define SCIP_DEFAULT_CUTAGELIMIT       100 /**< maximum age a cut can reach before it is deleted from global pool, or -1 */
#define SCIP_DEFAULT_CLEANUPROWS      TRUE /**< should new basic rows be removed after LP solving? */


/* Constraint Settings */

#define SCIP_DEFAULT_CONSAGELIMIT      200 /**< maximum age an unnecessary constr. can reach before it is deleted, or -1 */
#define SCIP_DEFAULT_CONSOBSOLETEAGE   100 /**< age of a constraint after which it is marked obsolete */


/* History Settings */

#define SCIP_DEFAULT_HISTORYEPS      1e-01 /**< default minimal variable distance value to use for history updates */
#define SCIP_DEFAULT_HISTORYDELTA    1e-04 /**< default minimal objective distance value to use for history updates */


/* Conflict Analysis */
#define SCIP_DEFAULT_USEPROPCONFLICT  TRUE /**< should propagation conflict analysis be used? */
#define SCIP_DEFAULT_USELPCONFLICT   FALSE /**< should infeasible LP conflict analysis be used? */
#define SCIP_DEFAULT_MAXCONFVARSFAC   0.02 /**< maximal fraction of binary variables involved in a conflict clause */
#define SCIP_DEFAULT_MINMAXCONFVARS     20 /**< minimal absolute maximum of variables involved in a conflict clause */


/* Primal Solutions */

#define SCIP_DEFAULT_GAPLIMIT          0.0 /**< solving stops, if the gap is below the given value */
#define SCIP_DEFAULT_SOLLIMIT           -1 /**< solving stops, if the given number of solutions were found (-1: no limit) */
#define SCIP_DEFAULT_MAXSOL            100 /**< maximal number of solutions to store in the solution storage */


/* Tree */

#define SCIP_DEFAULT_NODELIMIT        -1LL /**< maximal number of nodes to process (-1: no limit) */


/* Display */

#define SCIP_DEFAULT_DISPWIDTH         139 /**< maximal number of characters in a node information line */
#define SCIP_DEFAULT_DISPFREQ          100 /**< frequency for displaying node information lines */
#define SCIP_DEFAULT_DISPHEADERFREQ     15 /**< frequency for displaying header lines (every n'th node information line) */


/* Timing */

#define SCIP_DEFAULT_CLOCKTYPE  SCIP_CLOCKTYPE_CPU  /**< default clock type for timing */
#define SCIP_DEFAULT_CLOCKSENABLED    TRUE /**< is timing enabled? */
#define SCIP_DEFAULT_TIMELIMIT       1e+20 /**< maximal time in seconds to run */



/** calculate memory size for dynamically allocated arrays */
static
int calcGrowSize(
   int              initsize,           /**< initial size of array */
   Real             growfac,            /**< growing factor of array */
   int              num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(initsize >= 0);
   assert(growfac >= 1.0);

   if( growfac == 1.0 )
      size = MAX(initsize, num);
   else
   {
      /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
      size = initsize;
      while( size < num )
         size = (int)(growfac * size + 1);
   }

   return size;
}


/** information method for a parameter change of feastol */
static
DECL_PARAMCHGD(paramChgdFeastol)
{  /*lint --e{715}*/
   Real newfeastol;

   newfeastol = SCIPparamGetReal(param);
   
   /* change the feastol through the SCIP call in order to mark the LP unsolved */
   CHECK_OKAY( SCIPsetFeastol(scip, newfeastol) );

   return SCIP_OKAY;
}

/** information method for a parameter change of dualfeastol */
static
DECL_PARAMCHGD(paramChgdDualfeastol)
{  /*lint --e{715}*/
   Real newdualfeastol;

   newdualfeastol = SCIPparamGetReal(param);
   
   /* change the dualfeastol through the SCIP call in order to mark the LP unsolved */
   CHECK_OKAY( SCIPsetDualfeastol(scip, newdualfeastol) );

   return SCIP_OKAY;
}

/** creates global SCIP settings */
RETCODE SCIPsetCreate(
   SET**            set,                /**< pointer to SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   Real machineeps = SCIPcalcMachineEpsilon();

   assert(set != NULL);
   assert(scip != NULL);
   assert(machineeps < SCIP_MAXEPSILON);

   ALLOC_OKAY( allocMemory(set) );

   (*set)->scip = scip;

   CHECK_OKAY( SCIPparamsetCreate(&(*set)->paramset) );
   CHECK_OKAY( SCIPbufferCreate(&(*set)->buffer) );

   (*set)->readers = NULL;
   (*set)->nreaders = 0;
   (*set)->readerssize = 0;
   (*set)->pricers = NULL;
   (*set)->npricers = 0;
   (*set)->nactivepricers = 0;
   (*set)->pricerssize = 0;
   (*set)->pricerssorted = FALSE;
   (*set)->conshdlrs = NULL;
   (*set)->nconshdlrs = 0;
   (*set)->conshdlrssize = 0;
   (*set)->conflicthdlrs = NULL;
   (*set)->nconflicthdlrs = 0;
   (*set)->conflicthdlrssize = 0;
   (*set)->conflicthdlrssorted = FALSE;
   (*set)->presols = NULL;
   (*set)->npresols = 0;
   (*set)->presolssize = 0;
   (*set)->presolssorted = FALSE;
   (*set)->sepas = NULL;
   (*set)->nsepas = 0;
   (*set)->sepassize = 0;
   (*set)->sepassorted = FALSE;
   (*set)->heurs = NULL;
   (*set)->nheurs = 0;
   (*set)->heurssize = 0;
   (*set)->heurssorted = FALSE;
   (*set)->eventhdlrs = NULL;
   (*set)->neventhdlrs = 0;
   (*set)->eventhdlrssize = 0;
   (*set)->nodesels = NULL;
   (*set)->nnodesels = 0;
   (*set)->nodeselssize = 0;
   (*set)->actnodesel = NULL;
   (*set)->branchrules = NULL;
   (*set)->nbranchrules = 0;
   (*set)->branchrulessize = 0;
   (*set)->branchrulessorted = FALSE;
   (*set)->disps = NULL;
   (*set)->ndisps = 0;
   (*set)->dispssize = 0;

   /* SCIP parameters */
   assert(sizeof(int) == sizeof(VERBLEVEL));
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "display/verblevel",
                  "verbosity level of output",
                  (int*)&(*set)->verblevel, (int)SCIP_DEFAULT_VERBLEVEL, (int)SCIP_VERBLEVEL_NONE, (int)SCIP_VERBLEVEL_FULL,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "display/dispwidth",
                  "maximal number of characters in a node information line",
                  &(*set)->dispwidth, SCIP_DEFAULT_DISPWIDTH, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "display/dispfreq",
                  "frequency for displaying node information lines",
                  &(*set)->dispfreq, SCIP_DEFAULT_DISPFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "display/dispheaderfreq",
                  "frequency for displaying header lines (every n'th node information line)",
                  &(*set)->dispheaderfreq, SCIP_DEFAULT_DISPHEADERFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/infinity",
                  "values larger than this are considered infinity",
                  &(*set)->infinity, SCIP_DEFAULT_INFINITY, 1e+10, SCIP_INVALID/10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/epsilon",
                  "absolute values smaller than this are considered zero",
                  &(*set)->epsilon, SCIP_DEFAULT_EPSILON, machineeps, SCIP_MAXEPSILON,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/sumepsilon",
                  "absolute values of sums smaller than this are considered zero",
                  &(*set)->sumepsilon, SCIP_DEFAULT_SUMEPSILON, machineeps*1e+03, SCIP_MAXEPSILON,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/feastol",
                  "LP feasibility tolerance for constraints",
                  &(*set)->feastol, SCIP_DEFAULT_FEASTOL, machineeps*1e+03, SCIP_MAXEPSILON,
                  paramChgdFeastol, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/dualfeastol",
                  "LP feasibility tolerance for reduced costs",
                  &(*set)->dualfeastol, SCIP_DEFAULT_DUALFEASTOL, machineeps*1e+03, SCIP_MAXEPSILON,
                  paramChgdDualfeastol, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/boundstreps",
                  "minimal improve for strengthening bounds",
                  &(*set)->boundstreps, SCIP_DEFAULT_BOUNDSTREPS, machineeps*1e+03, SCIP_INVALID/10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/cutvioleps",
                  "epsilon for deciding if a cut is violated",
                  &(*set)->cutvioleps, SCIP_DEFAULT_CUTVIOLEPS, machineeps*1e+03, SCIP_INVALID/10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/cutviolepsroot",
                  "epsilon for deciding if a cut is violated in the root node",
                  &(*set)->cutviolepsroot, 0.05*SCIP_DEFAULT_CUTVIOLEPS, machineeps*1e+03, SCIP_INVALID/10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/historyeps",
                  "minimal variable distance value to use for branching history updates",
                  &(*set)->historyeps, SCIP_DEFAULT_HISTORYEPS, machineeps*1e+03, 1.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "numerics/historydelta",
                  "minimal objective distance value to use for branching history updates",
                  &(*set)->historydelta, SCIP_DEFAULT_HISTORYDELTA, 0.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "memory/memsavefac",
                  "fraction of maximal memory usage resulting in switch to memory saving mode",
                  &(*set)->memsavefac, SCIP_DEFAULT_MEMSAVEFAC, 0.0, 1.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "memory/memgrowfac",
                  "memory growing factor for dynamically allocated arrays",
                  &(*set)->memgrowfac, SCIP_DEFAULT_MEMGROWFAC, 1.0, 10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "memory/memgrowinit",
                  "initial size of dynamically allocated arrays",
                  &(*set)->memgrowinit, SCIP_DEFAULT_MEMGROWINIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "memory/treegrowfac",
                  "memory growing factor for tree array",
                  &(*set)->treegrowfac, SCIP_DEFAULT_TREEGROWFAC, 1.0, 10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "memory/treegrowinit",
                  "initial size of tree array",
                  &(*set)->treegrowinit, SCIP_DEFAULT_TREEGROWINIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "memory/pathgrowfac",
                  "memory growing factor for path array",
                  &(*set)->pathgrowfac, SCIP_DEFAULT_PATHGROWFAC, 1.0, 10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "memory/pathgrowinit",
                  "initial size of path array",
                  &(*set)->pathgrowinit, SCIP_DEFAULT_PATHGROWINIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "branching/branchscorefac",
                  "branching score factor to weigh downward and upward gain prediction",
                  &(*set)->branchscorefac, SCIP_DEFAULT_BRANCHSCOREFAC, 0.0, 1.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "presolving/maxpresolrounds",
                  "maximal number of presolving rounds (-1: unlimited)",
                  &(*set)->maxpresolrounds, SCIP_DEFAULT_MAXPRESOLROUNDS, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "presolving/presolabortfac",
                  "abort presolve, if less than this fraction of the problem was changed in last presolve round",
                  &(*set)->presolabortfac, SCIP_DEFAULT_PRESOLABORTFAC, 0.0, 1.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "pricing/maxpricevars",
                  "maximal number of variables priced in per pricing round",
                  &(*set)->maxpricevars, SCIP_DEFAULT_MAXPRICEVARS, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "pricing/maxpricevarsroot",
                  "maximal number of priced variables at the root node",
                  &(*set)->maxpricevarsroot, SCIP_DEFAULT_MAXPRICEVARSROOT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "pricing/abortpricevarsfac",
                  "pricing is aborted, if fac * maxpricevars pricing candidates were found",
                  &(*set)->abortpricevarsfac, SCIP_DEFAULT_ABORTPRICEVARSFAC, 1.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "separating/maxsepacuts",
                  "maximal number of cuts separated per separation round",
                  &(*set)->maxsepacuts, SCIP_DEFAULT_MAXSEPACUTS, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "separating/maxsepacutsroot",
                  "maximal number of separated cuts at the root node",
                  &(*set)->maxsepacutsroot, SCIP_DEFAULT_MAXSEPACUTSROOT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "separating/cutagelimit",
                  "maximum age a cut can reach before it is deleted from the global cut pool, or -1 to keep all cuts",
                  &(*set)->cutagelimit, SCIP_DEFAULT_CUTAGELIMIT, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "constraints/consagelimit",
                  "maximum age an unnecessary constraint can reach before it is deleted, or -1 to keep all constraints",
                  &(*set)->consagelimit, SCIP_DEFAULT_CONSAGELIMIT, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "constraints/consobsoleteage",
                  "age of a constraint after which it is marked obsolete, or -1 to not mark constraints obsolete",
                  &(*set)->consobsoleteage, SCIP_DEFAULT_CONSOBSOLETEAGE, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "conflict/usepropconflict",
                  "should propagation conflict analysis be used?",
                  &(*set)->usepropconflict, SCIP_DEFAULT_USEPROPCONFLICT,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "conflict/uselpconflict",
                  "should infeasible LP conflict analysis be used?",
                  &(*set)->uselpconflict, SCIP_DEFAULT_USELPCONFLICT,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "conflict/maxconfvarsfac",
                  "maximal fraction of binary variables involved in a conflict clause",
                  &(*set)->maxconfvarsfac, SCIP_DEFAULT_MAXCONFVARSFAC, 0.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "conflict/minmaxconfvars",
                  "minimal absolute maximum of variables involved in a conflict clause",
                  &(*set)->minmaxconfvars, SCIP_DEFAULT_MINMAXCONFVARS, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddLongintParam(*set, memhdr,
                  "limits/nodelimit",
                  "maximal number of nodes to process (-1: no limit)",
                  &(*set)->nodelimit, SCIP_DEFAULT_NODELIMIT, -1LL, LONGINT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "limits/timelimit",
                  "maximal time in seconds to run",
                  &(*set)->timelimit, SCIP_DEFAULT_TIMELIMIT, 0.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "limits/memlimit",
                  "maximal memory usage in MB; reported memory usage is lower than real memory usage!",
                  &(*set)->memlimit, SCIP_DEFAULT_MEMLIMIT, 0.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "limits/gaplimit",
                  "solving stops, if the gap = |(primalbound - dualbound)/dualbound| is below the given value",
                  &(*set)->gaplimit, SCIP_DEFAULT_GAPLIMIT, 0.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "limits/sollimit",
                  "solving stops, if the given number of solutions were found (-1: no limit)",
                  &(*set)->sollimit, SCIP_DEFAULT_SOLLIMIT, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "limits/maxsol",
                  "maximal number of solutions to store in the solution storage",
                  &(*set)->maxsol, SCIP_DEFAULT_MAXSOL, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "lp/lpsolvefreq",
                  "frequency for solving LP at the nodes (-1: never; 0: only root LP)",
                  &(*set)->lpsolvefreq, SCIP_DEFAULT_LPSOLVEFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "lp/lpsolvedepth",
                  "maximal depth for solving LP at the nodes (-1: no depth limit)",
                  &(*set)->lpsolvedepth, SCIP_DEFAULT_LPSOLVEDEPTH, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "lp/redcostfreq",
                  "frequency for applying reduced cost fixing (-1: never; 0: only root LP)",
                  &(*set)->redcostfreq, SCIP_DEFAULT_REDCOSTFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "lp/checklpfeas",
                  "should LP solutions be checked, resolving LP when numerical troubles occur?",
                  &(*set)->checklpfeas, SCIP_DEFAULT_CHECKLPFEAS,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "lp/fastmip",
                  "should FASTMIP setting of LP solver be used?",
                  &(*set)->fastmip, SCIP_DEFAULT_FASTMIP,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "lp/scaling",
                  "should scaling of LP solver be used?",
                  &(*set)->scaling, SCIP_DEFAULT_SCALING,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "lp/cleanupcols",
                  "should new non-basic columns be removed after LP solving?",
                  &(*set)->cleanupcols, SCIP_DEFAULT_CLEANUPCOLS,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "lp/cleanuprows",
                  "should new basic rows be removed after LP solving?",
                  &(*set)->cleanuprows, SCIP_DEFAULT_CLEANUPROWS,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "lp/colagelimit",
                  "maximum age a dynamic column can reach before it is deleted from the LP",
                  &(*set)->colagelimit, SCIP_DEFAULT_COLAGELIMIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "lp/rowagelimit",
                  "maximum age a dynamic row can reach before it is deleted from the LP",
                  &(*set)->rowagelimit, SCIP_DEFAULT_ROWAGELIMIT, 0, INT_MAX,
                  NULL, NULL) );
   assert(sizeof(int) == sizeof(CLOCKTYPE));
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "timing/clocktype",
                  "default clock type (1: CPU user seconds, 2: wall clock time)",
                  (int*)&(*set)->clocktype, (int)SCIP_DEFAULT_CLOCKTYPE, 1, 2,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "timing/clocksenabled",
                  "is timing enabled?",
                  &(*set)->clocksenabled, SCIP_DEFAULT_CLOCKSENABLED,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "misc/catchctrlc",
                  "should the CTRL-C interrupt be catched by SCIP?",
                  &(*set)->catchctrlc, SCIP_DEFAULT_CATCHCTRLC,
                  NULL, NULL) );

   return SCIP_OKAY;
}

/** frees global SCIP settings */
RETCODE SCIPsetFree(
   SET**            set,                /**< pointer to SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   int i;

   assert(set != NULL);

   /* free parameter set */
   SCIPparamsetFree(&(*set)->paramset, memhdr);

   /* free memory buffers */
   SCIPbufferFree(&(*set)->buffer);

   /* free file readers */
   for( i = 0; i < (*set)->nreaders; ++i )
   {
      CHECK_OKAY( SCIPreaderFree(&(*set)->readers[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->readers);

   /* free variable pricers */
   for( i = 0; i < (*set)->npricers; ++i )
   {
      CHECK_OKAY( SCIPpricerFree(&(*set)->pricers[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->pricers);

   /* free constraint handlers */
   for( i = 0; i < (*set)->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrFree(&(*set)->conshdlrs[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->conshdlrs);

   /* free conflict handlers */
   for( i = 0; i < (*set)->nconflicthdlrs; ++i )
   {
      CHECK_OKAY( SCIPconflicthdlrFree(&(*set)->conflicthdlrs[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->conflicthdlrs);

   /* free presolvers */
   for( i = 0; i < (*set)->npresols; ++i )
   {
      CHECK_OKAY( SCIPpresolFree(&(*set)->presols[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->presols);

   /* free separators */
   for( i = 0; i < (*set)->nsepas; ++i )
   {
      CHECK_OKAY( SCIPsepaFree(&(*set)->sepas[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->sepas);

   /* free primal heuristics */
   for( i = 0; i < (*set)->nheurs; ++i )
   {
      CHECK_OKAY( SCIPheurFree(&(*set)->heurs[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->heurs);

   /* free event handlers */
   for( i = 0; i < (*set)->neventhdlrs; ++i )
   {
      CHECK_OKAY( SCIPeventhdlrFree(&(*set)->eventhdlrs[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->eventhdlrs);

   /* free node selectors */
   for( i = 0; i < (*set)->nnodesels; ++i )
   {
      CHECK_OKAY( SCIPnodeselFree(&(*set)->nodesels[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->nodesels);

   /* free branching methods */
   for( i = 0; i < (*set)->nbranchrules; ++i )
   {
      CHECK_OKAY( SCIPbranchruleFree(&(*set)->branchrules[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->branchrules);

   /* free display columns */
   for( i = 0; i < (*set)->ndisps; ++i )
   {
      CHECK_OKAY( SCIPdispFree(&(*set)->disps[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->disps);

   freeMemory(set);

   return SCIP_OKAY;
}

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPsetAddBoolParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetAddBool(set->paramset, memhdr, name, desc, valueptr, defaultvalue, paramchgd, paramdata) );
   
   return SCIP_OKAY;
}

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPsetAddIntParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue,       /**< default value of the parameter */
   int              minvalue,           /**< minimum value for parameter */
   int              maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetAddInt(set->paramset, memhdr, name, desc, valueptr, defaultvalue, minvalue, maxvalue, 
                  paramchgd, paramdata) );
   
   return SCIP_OKAY;
}

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPsetAddLongintParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue,       /**< default value of the parameter */
   Longint          minvalue,           /**< minimum value for parameter */
   Longint          maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetAddLongint(set->paramset, memhdr, name, desc, valueptr, defaultvalue, minvalue, maxvalue, 
                  paramchgd, paramdata) );
   
   return SCIP_OKAY;
}

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPsetAddRealParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue,       /**< default value of the parameter */
   Real             minvalue,           /**< minimum value for parameter */
   Real             maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetAddReal(set->paramset, memhdr, name, desc, valueptr, defaultvalue, minvalue, maxvalue, 
                  paramchgd, paramdata) );
   
   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPsetAddCharParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue,       /**< default value of the parameter */
   const char*      allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetAddChar(set->paramset, memhdr, name, desc, valueptr, defaultvalue, allowedvalues, 
                  paramchgd, paramdata) );
   
   return SCIP_OKAY;
}

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPsetAddStringParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetAddString(set->paramset, memhdr, name, desc, valueptr, defaultvalue, paramchgd, paramdata) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing Bool parameter */
RETCODE SCIPsetGetBoolParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Bool*            value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetGetBool(set->paramset, name, value) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing Int parameter */
RETCODE SCIPsetGetIntParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   int*             value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetGetInt(set->paramset, name, value) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing Longint parameter */
RETCODE SCIPsetGetLongintParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Longint*         value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetGetLongint(set->paramset, name, value) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing Real parameter */
RETCODE SCIPsetGetRealParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Real*            value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetGetReal(set->paramset, name, value) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing Char parameter */
RETCODE SCIPsetGetCharParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   char*            value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetGetChar(set->paramset, name, value) );
   
   return SCIP_OKAY;
}

/** gets the value of an existing String parameter */
RETCODE SCIPsetGetStringParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   char**           value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetGetString(set->paramset, name, value) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing Bool parameter */
RETCODE SCIPsetSetBoolParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetSetBool(set->paramset, set->scip, name, value) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing Int parameter */
RETCODE SCIPsetSetIntParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetSetInt(set->paramset, set->scip, name, value) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing Longint parameter */
RETCODE SCIPsetSetLongintParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetSetLongint(set->paramset, set->scip, name, value) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing Real parameter */
RETCODE SCIPsetSetRealParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetSetReal(set->paramset, set->scip, name, value) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing Char parameter */
RETCODE SCIPsetSetCharParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetSetChar(set->paramset, set->scip, name, value) );
   
   return SCIP_OKAY;
}

/** changes the value of an existing String parameter */
RETCODE SCIPsetSetStringParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetSetString(set->paramset, set->scip, name, value) );
   
   return SCIP_OKAY;
}

/** reads parameters from a file */
RETCODE SCIPsetReadParams(
   SET*             set,                /**< global SCIP settings */
   const char*      filename            /**< file name */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetRead(set->paramset, set->scip, filename) );

   return SCIP_OKAY;
}

/** writes all parameters in the parameter set to a file */
RETCODE SCIPsetWriteParams(
   SET*             set,                /**< global SCIP settings */
   const char*      filename,           /**< file name, or NULL for stdout */
   Bool             comments,           /**< should parameter descriptions be written as comments? */
   Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetWrite(set->paramset, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** returns the array of all available SCIP parameters */
PARAM** SCIPsetGetParams(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return SCIPparamsetGetParams(set->paramset);
}

/** returns the total number of all available SCIP parameters */
int SCIPsetGetNParams(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return SCIPparamsetGetNParams(set->paramset);
}

/** inserts file reader in file reader list */
RETCODE SCIPsetIncludeReader(
   SET*             set,                /**< global SCIP settings */
   READER*          reader              /**< file reader */
   )
{
   assert(set != NULL);
   assert(reader != NULL);

   if( set->nreaders >= set->readerssize )
   {
      set->readerssize = SCIPsetCalcMemGrowSize(set, set->nreaders+1);
      ALLOC_OKAY( reallocMemoryArray(&set->readers, set->readerssize) );
   }
   assert(set->nreaders < set->readerssize);
   
   set->readers[set->nreaders] = reader;
   set->nreaders++;

   return SCIP_OKAY;
}   

/** returns the file reader of the given name, or NULL if not existing */
READER* SCIPsetFindReader(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of file reader */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nreaders; ++i )
   {
      if( strcmp(SCIPreaderGetName(set->readers[i]), name) == 0 )
         return set->readers[i];
   }

   return NULL;
}

/** inserts variable pricer in variable pricer list */
RETCODE SCIPsetIncludePricer(
   SET*             set,                /**< global SCIP settings */
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(set != NULL);
   assert(pricer != NULL);

   if( set->npricers >= set->pricerssize )
   {
      set->pricerssize = SCIPsetCalcMemGrowSize(set, set->npricers+1);
      ALLOC_OKAY( reallocMemoryArray(&set->pricers, set->pricerssize) );
   }
   assert(set->npricers < set->pricerssize);
   
   set->pricers[set->npricers] = pricer;
   set->npricers++;
   set->pricerssorted = FALSE;

   return SCIP_OKAY;
}   

/** returns the variable pricer of the given name, or NULL if not existing */
PRICER* SCIPsetFindPricer(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of variable pricer */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->npricers; ++i )
   {
      if( strcmp(SCIPpricerGetName(set->pricers[i]), name) == 0 )
         return set->pricers[i];
   }

   return NULL;
}

/** sorts pricers by priorities */
void SCIPsetSortPricers(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->pricerssorted )
   {
      SCIPbsortPtr((void**)set->pricers, set->npricers, SCIPpricerComp);
      set->pricerssorted = TRUE;
   }
}

/** inserts constraint handler in constraint handler list */
RETCODE SCIPsetIncludeConshdlr(
   SET*             set,                /**< global SCIP settings */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   int checkpriority;
   int i;

   assert(set != NULL);
   assert(conshdlr != NULL);
   assert(!SCIPconshdlrIsInitialized(conshdlr));

   if( set->nconshdlrs >= set->conshdlrssize )
   {
      set->conshdlrssize = SCIPsetCalcMemGrowSize(set, set->nconshdlrs+1);
      ALLOC_OKAY( reallocMemoryArray(&set->conshdlrs, set->conshdlrssize) );
   }
   assert(set->nconshdlrs < set->conshdlrssize);

   checkpriority = SCIPconshdlrGetCheckPriority(conshdlr);
   for( i = set->nconshdlrs; i > 0 && SCIPconshdlrGetCheckPriority(set->conshdlrs[i-1]) < checkpriority; --i )
   {
      set->conshdlrs[i] = set->conshdlrs[i-1];
   }

   set->conshdlrs[i] = conshdlr;
   set->nconshdlrs++;

   return SCIP_OKAY;
}   

/** returns the constraint handler of the given name, or NULL if not existing */
CONSHDLR* SCIPsetFindConshdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of constraint handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nconshdlrs; ++i )
   {
      if( strcmp(SCIPconshdlrGetName(set->conshdlrs[i]), name) == 0 )
         return set->conshdlrs[i];
   }

   return NULL;
}

/** inserts conflict handler in conflict handler list */
RETCODE SCIPsetIncludeConflicthdlr(
   SET*             set,                /**< global SCIP settings */
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(set != NULL);
   assert(conflicthdlr != NULL);
   assert(!SCIPconflicthdlrIsInitialized(conflicthdlr));

   if( set->nconflicthdlrs >= set->conflicthdlrssize )
   {
      set->conflicthdlrssize = SCIPsetCalcMemGrowSize(set, set->nconflicthdlrs+1);
      ALLOC_OKAY( reallocMemoryArray(&set->conflicthdlrs, set->conflicthdlrssize) );
   }
   assert(set->nconflicthdlrs < set->conflicthdlrssize);

   set->conflicthdlrs[set->nconflicthdlrs] = conflicthdlr;
   set->nconflicthdlrs++;
   set->conflicthdlrssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the conflict handler of the given name, or NULL if not existing */
CONFLICTHDLR* SCIPsetFindConflicthdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of conflict handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      if( strcmp(SCIPconflicthdlrGetName(set->conflicthdlrs[i]), name) == 0 )
         return set->conflicthdlrs[i];
   }

   return NULL;
}

/** sorts conflict handlers by priorities */
void SCIPsetSortConflicthdlrs(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->conflicthdlrssorted )
   {
      SCIPbsortPtr((void**)set->conflicthdlrs, set->nconflicthdlrs, SCIPconflicthdlrComp);
      set->conflicthdlrssorted = TRUE;
   }
}

/** inserts presolver in presolver list */
RETCODE SCIPsetIncludePresol(
   SET*             set,                /**< global SCIP settings */
   PRESOL*          presol              /**< presolver */
   )
{
   assert(set != NULL);
   assert(presol != NULL);

   if( set->npresols >= set->presolssize )
   {
      set->presolssize = SCIPsetCalcMemGrowSize(set, set->npresols+1);
      ALLOC_OKAY( reallocMemoryArray(&set->presols, set->presolssize) );
   }
   assert(set->npresols < set->presolssize);
   
   set->presols[set->npresols] = presol;
   set->npresols++;
   set->presolssorted = FALSE;

   return SCIP_OKAY;
}   

/** returns the presolver of the given name, or NULL if not existing */
PRESOL* SCIPsetFindPresol(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of presolver */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->npresols; ++i )
   {
      if( strcmp(SCIPpresolGetName(set->presols[i]), name) == 0 )
         return set->presols[i];
   }

   return NULL;
}

/** sorts presolvers by priorities */
void SCIPsetSortPresols(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->presolssorted )
   {
      SCIPbsortPtr((void**)set->presols, set->npresols, SCIPpresolComp);
      set->presolssorted = TRUE;
   }
}

/** inserts separator in separator list */
RETCODE SCIPsetIncludeSepa(
   SET*             set,                /**< global SCIP settings */
   SEPA*            sepa                /**< separator */
   )
{
   assert(set != NULL);
   assert(sepa != NULL);
   assert(!SCIPsepaIsInitialized(sepa));

   if( set->nsepas >= set->sepassize )
   {
      set->sepassize = SCIPsetCalcMemGrowSize(set, set->nsepas+1);
      ALLOC_OKAY( reallocMemoryArray(&set->sepas, set->sepassize) );
   }
   assert(set->nsepas < set->sepassize);
   
   set->sepas[set->nsepas] = sepa;
   set->nsepas++;
   set->sepassorted = FALSE;

   return SCIP_OKAY;
}   

/** returns the separator of the given name, or NULL if not existing */
SEPA* SCIPsetFindSepa(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of separator */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nsepas; ++i )
   {
      if( strcmp(SCIPsepaGetName(set->sepas[i]), name) == 0 )
         return set->sepas[i];
   }

   return NULL;
}

/** sorts separators by priorities */
void SCIPsetSortSepas(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->sepassorted )
   {
      SCIPbsortPtr((void**)set->sepas, set->nsepas, SCIPsepaComp);
      set->sepassorted = TRUE;
   }
}

/** inserts primal heuristic in primal heuristic list */
RETCODE SCIPsetIncludeHeur(
   SET*             set,                /**< global SCIP settings */
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(set != NULL);
   assert(heur != NULL);
   assert(!SCIPheurIsInitialized(heur));

   if( set->nheurs >= set->heurssize )
   {
      set->heurssize = SCIPsetCalcMemGrowSize(set, set->nheurs+1);
      ALLOC_OKAY( reallocMemoryArray(&set->heurs, set->heurssize) );
   }
   assert(set->nheurs < set->heurssize);
   
   set->heurs[set->nheurs] = heur;
   set->nheurs++;
   set->heurssorted = FALSE;

   return SCIP_OKAY;
}   

/** returns the primal heuristic of the given name, or NULL if not existing */
HEUR* SCIPsetFindHeur(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of primal heuristic */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nheurs; ++i )
   {
      if( strcmp(SCIPheurGetName(set->heurs[i]), name) == 0 )
         return set->heurs[i];
   }

   return NULL;
}

/** sorts heuristics by priorities */
void SCIPsetSortHeurs(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->heurssorted )
   {
      SCIPbsortPtr((void**)set->heurs, set->nheurs, SCIPheurComp);
      set->heurssorted = TRUE;
   }
}

/** inserts event handler in event handler list */
RETCODE SCIPsetIncludeEventhdlr(
   SET*             set,                /**< global SCIP settings */
   EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(set != NULL);
   assert(eventhdlr != NULL);
   assert(!SCIPeventhdlrIsInitialized(eventhdlr));

   if( set->neventhdlrs >= set->eventhdlrssize )
   {
      set->eventhdlrssize = SCIPsetCalcMemGrowSize(set, set->neventhdlrs+1);
      ALLOC_OKAY( reallocMemoryArray(&set->eventhdlrs, set->eventhdlrssize) );
   }
   assert(set->neventhdlrs < set->eventhdlrssize);
   
   set->eventhdlrs[set->neventhdlrs] = eventhdlr;
   set->neventhdlrs++;

   return SCIP_OKAY;
}   

/** returns the event handler of the given name, or NULL if not existing */
EVENTHDLR* SCIPsetFindEventhdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->neventhdlrs; ++i )
   {
      if( strcmp(SCIPeventhdlrGetName(set->eventhdlrs[i]), name) == 0 )
         return set->eventhdlrs[i];
   }

   return NULL;
}

/** inserts node selector in node selector list */
RETCODE SCIPsetIncludeNodesel(
   SET*             set,                /**< global SCIP settings */
   NODESEL*         nodesel             /**< node selector */
   )
{
   int i;

   assert(set != NULL);
   assert(nodesel != NULL);
   assert(!SCIPnodeselIsInitialized(nodesel));

   if( set->nnodesels >= set->nodeselssize )
   {
      set->nodeselssize = SCIPsetCalcMemGrowSize(set, set->nnodesels+1);
      ALLOC_OKAY( reallocMemoryArray(&set->nodesels, set->nodeselssize) );
   }
   assert(set->nnodesels < set->nodeselssize);
   
   for( i = set->nnodesels; i > 0 && SCIPnodeselGetStdPriority(nodesel) > SCIPnodeselGetStdPriority(set->nodesels[i-1]);
        --i )
   {
      set->nodesels[i] = set->nodesels[i-1];
   }
   set->nodesels[i] = nodesel;
   set->nnodesels++;

   return SCIP_OKAY;
}   

/** returns the node selector of the given name, or NULL if not existing */
NODESEL* SCIPsetFindNodesel(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nnodesels; ++i )
   {
      if( strcmp(SCIPnodeselGetName(set->nodesels[i]), name) == 0 )
         return set->nodesels[i];
   }

   return NULL;
}

/** returns node selector with highest priority in the current mode */
NODESEL* SCIPsetGetActNodesel(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(set != NULL);
   assert(stat != NULL);

   /* check, if old node selector is still valid */
   if( set->actnodesel == NULL && set->nnodesels > 0 )
   {
      int i;

      set->actnodesel = set->nodesels[0];

      /* search highest priority node selector */
      if( stat->memsavemode )
      {
         for( i = 1; i < set->nnodesels; ++i )
         {
            if( SCIPnodeselGetMemsavePriority(set->nodesels[i]) > SCIPnodeselGetMemsavePriority(set->actnodesel) )
               set->actnodesel = set->nodesels[i];
         }
      }
      else
      {
         for( i = 1; i < set->nnodesels; ++i )
         {
            if( SCIPnodeselGetStdPriority(set->nodesels[i]) > SCIPnodeselGetStdPriority(set->actnodesel) )
               set->actnodesel = set->nodesels[i];
         }
      }
   }
   
   return set->actnodesel;
}

/** inserts branching rule in branching rule list */
RETCODE SCIPsetIncludeBranchrule(
   SET*             set,                /**< global SCIP settings */
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(set != NULL);
   assert(branchrule != NULL);
   assert(!SCIPbranchruleIsInitialized(branchrule));

   if( set->nbranchrules >= set->branchrulessize )
   {
      set->branchrulessize = SCIPsetCalcMemGrowSize(set, set->nbranchrules+1);
      ALLOC_OKAY( reallocMemoryArray(&set->branchrules, set->branchrulessize) );
   }
   assert(set->nbranchrules < set->branchrulessize);

   set->branchrules[set->nbranchrules] = branchrule;
   set->nbranchrules++;
   set->branchrulessorted = FALSE;

   return SCIP_OKAY;
}   

/** returns the branching rule of the given name, or NULL if not existing */
BRANCHRULE* SCIPsetFindBranchrule(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nbranchrules; ++i )
   {
      if( strcmp(SCIPbranchruleGetName(set->branchrules[i]), name) == 0 )
         return set->branchrules[i];
   }

   return NULL;
}

/** sorts branching rules by priorities */
void SCIPsetSortBranchrules(
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->branchrulessorted )
   {
      SCIPbsortPtr((void**)set->branchrules, set->nbranchrules, SCIPbranchruleComp);
      set->branchrulessorted = TRUE;
   }
}

/** inserts display column in display column list */
RETCODE SCIPsetIncludeDisp(
   SET*             set,                /**< global SCIP settings */
   DISP*            disp                /**< display column */
   )
{
   int i;

   assert(set != NULL);
   assert(disp != NULL);
   assert(!SCIPdispIsInitialized(disp));

   if( set->ndisps >= set->dispssize )
   {
      set->dispssize = SCIPsetCalcMemGrowSize(set, set->ndisps+1);
      ALLOC_OKAY( reallocMemoryArray(&set->disps, set->dispssize) );
   }
   assert(set->ndisps < set->dispssize);

   for( i = set->ndisps; i > 0 && SCIPdispGetPosition(disp) < SCIPdispGetPosition(set->disps[i-1]); --i )
   {
      set->disps[i] = set->disps[i-1];
   }
   set->disps[i] = disp;
   set->ndisps++;

   return SCIP_OKAY;
}   

/** returns the display column of the given name, or NULL if not existing */
DISP* SCIPsetFindDisp(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->ndisps; ++i )
   {
      if( strcmp(SCIPdispGetName(set->disps[i]), name) == 0 )
         return set->disps[i];
   }

   return NULL;
}

/** initializes all user callback functions */
RETCODE SCIPsetInitCallbacks(
   SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* active variable pricers */
   SCIPsetSortPricers(set);
   for( i = 0; i < set->nactivepricers; ++i )
   {
      CHECK_OKAY( SCIPpricerInit(set->pricers[i], set->scip) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrInit(set->conshdlrs[i], set->scip) );
   }

   /* conflict handlers */
   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      CHECK_OKAY( SCIPconflicthdlrInit(set->conflicthdlrs[i], set->scip) );
   }

   /* presolvers */
   for( i = 0; i < set->npresols; ++i )
   {
      CHECK_OKAY( SCIPpresolInit(set->presols[i], set->scip) );
   }

   /* separators */
   for( i = 0; i < set->nsepas; ++i )
   {
      CHECK_OKAY( SCIPsepaInit(set->sepas[i], set->scip) );
   }

   /* primal heuristics */
   for( i = 0; i < set->nheurs; ++i )
   {
      CHECK_OKAY( SCIPheurInit(set->heurs[i], set->scip) );
   }

   /* event handlers */
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      CHECK_OKAY( SCIPeventhdlrInit(set->eventhdlrs[i], set->scip) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      CHECK_OKAY( SCIPnodeselInit(set->nodesels[i], set->scip) );
   }

   /* branching rules */
   for( i = 0; i < set->nbranchrules; ++i )
   {
      CHECK_OKAY( SCIPbranchruleInit(set->branchrules[i], set->scip) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      CHECK_OKAY( SCIPdispInit(set->disps[i], set->scip) );
   }
   CHECK_OKAY( SCIPdispAutoActivate(set) );

   return SCIP_OKAY;
}

/** calls exit methods of all user callback functions */
RETCODE SCIPsetExitCallbacks(
   SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* active variable pricers */
   SCIPsetSortPricers(set);
   for( i = 0; i < set->nactivepricers; ++i )
   {
      CHECK_OKAY( SCIPpricerExit(set->pricers[i], set->scip) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrExit(set->conshdlrs[i], set->scip) );
   }

   /* conflict handlers */
   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      CHECK_OKAY( SCIPconflicthdlrExit(set->conflicthdlrs[i], set->scip) );
   }

   /* presolvers */
   for( i = 0; i < set->npresols; ++i )
   {
      CHECK_OKAY( SCIPpresolExit(set->presols[i], set->scip) );
   }

   /* separators */
   for( i = 0; i < set->nsepas; ++i )
   {
      CHECK_OKAY( SCIPsepaExit(set->sepas[i], set->scip) );
   }

   /* primal heuristics */
   for( i = 0; i < set->nheurs; ++i )
   {
      CHECK_OKAY( SCIPheurExit(set->heurs[i], set->scip) );
   }

   /* event handlers */
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      CHECK_OKAY( SCIPeventhdlrExit(set->eventhdlrs[i], set->scip) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      CHECK_OKAY( SCIPnodeselExit(set->nodesels[i], set->scip) );
   }

   /* branchruleing rules */
   for( i = 0; i < set->nbranchrules; ++i )
   {
      CHECK_OKAY( SCIPbranchruleExit(set->branchrules[i], set->scip) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      CHECK_OKAY( SCIPdispExit(set->disps[i], set->scip) );
   }

   return SCIP_OKAY;
}

/** calculate memory size for dynamically allocated arrays */
int SCIPsetCalcMemGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->memgrowinit, set->memgrowfac, num);
}

/** calculate memory size for tree array */
int SCIPsetCalcTreeGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->treegrowinit, set->treegrowfac, num);
}

/** calculate memory size for path array */
int SCIPsetCalcPathGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->pathgrowinit, set->pathgrowfac, num);
}

/** sets verbosity level for message output */
RETCODE SCIPsetSetVerbLevel(
   SET*             set,                /**< global SCIP settings */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   )
{
   assert(set != NULL);

   if( verblevel > SCIP_VERBLEVEL_FULL )
   {
      errorMessage("Invalid verbosity level <%d>, maximum is <%d>\n", verblevel, SCIP_VERBLEVEL_FULL);
      return SCIP_INVALIDCALL;
   }
   
   set->verblevel = verblevel;

   return SCIP_OKAY;
}

/** sets LP feasibility tolerance */
RETCODE SCIPsetSetFeastol(
   SET*             set,                /**< global SCIP settings */
   Real             feastol             /**< new feasibility tolerance */
   )
{
   assert(set != NULL);

   set->feastol = feastol;

   return SCIP_OKAY;
}

/** sets LP feasibility tolerance for reduced costs */
RETCODE SCIPsetSetDualfeastol(
   SET*             set,                /**< global SCIP settings */
   Real             dualfeastol         /**< new reduced costs feasibility tolerance */
   )
{
   assert(set != NULL);

   set->dualfeastol = dualfeastol;

   return SCIP_OKAY;
}

/** returns the maximal number of variables priced into the LP per round */
int SCIPsetGetMaxpricevars(
   const SET*       set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   )
{
   assert(set != NULL);

   if( root )
      return set->maxpricevarsroot;
   else
      return set->maxpricevars;
}

/** returns the maximal number of cuts separated per round */
int SCIPsetGetMaxsepacuts(
   const SET*       set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   )
{
   assert(set != NULL);

   if( root )
      return set->maxsepacutsroot;
   else
      return set->maxsepacuts;
}

   

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
Real SCIPsetRelDiff(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real absval1;
   Real absval2;
   Real quot;

   assert(set != NULL);

   absval1 = ABS(val1);
   absval2 = ABS(val2);
   quot = MAX(absval1, absval2);
   quot = MAX(quot, 1.0);
   
   return (val1-val2)/quot;
}

/** checks, if values are in range of epsilon */
Bool SCIPsetIsEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSEQ(val1, val2, set->epsilon);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
Bool SCIPsetIsLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSLT(val1, val2, set->epsilon);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
Bool SCIPsetIsLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSLE(val1, val2, set->epsilon);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
Bool SCIPsetIsGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSGT(val1, val2, set->epsilon);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
Bool SCIPsetIsGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSGE(val1, val2, set->epsilon);
}

/** checks, if value is in range epsilon of 0.0 */
Bool SCIPsetIsZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->epsilon);
}

/** checks, if value is greater than epsilon */
Bool SCIPsetIsPositive(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSP(val, set->epsilon);
}

/** checks, if value is lower than -epsilon */
Bool SCIPsetIsNegative(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSN(val, set->epsilon);
}

/** checks, if values are in range of sumepsilon */
Bool SCIPsetIsSumEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSEQ(val1, val2, set->sumepsilon);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
Bool SCIPsetIsSumLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSLT(val1, val2, set->sumepsilon);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
Bool SCIPsetIsSumLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSLE(val1, val2, set->sumepsilon);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
Bool SCIPsetIsSumGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSGT(val1, val2, set->sumepsilon);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
Bool SCIPsetIsSumGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSGE(val1, val2, set->sumepsilon);
}

/** checks, if value is in range sumepsilon of 0.0 */
Bool SCIPsetIsSumZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->sumepsilon);
}

/** checks, if value is greater than sumepsilon */
Bool SCIPsetIsSumPositive(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSP(val, set->sumepsilon);
}

/** checks, if value is lower than -sumepsilon */
Bool SCIPsetIsSumNegative(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSN(val, set->sumepsilon);
}

/** checks, if relative difference of values is in range of feastol */
Bool SCIPsetIsFeasEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSZ(diff, set->feastol);
}

/** checks, if relative difference of val1 and val2 is lower than feastol */
Bool SCIPsetIsFeasLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSN(diff, set->feastol);
}

/** checks, if relative difference of val1 and val2 is not greater than feastol */
Bool SCIPsetIsFeasLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return !EPSP(diff, set->feastol);
}

/** checks, if relative difference of val1 and val2 is greater than feastol */
Bool SCIPsetIsFeasGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSP(diff, set->feastol);
}

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
Bool SCIPsetIsFeasGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return !EPSN(diff, set->feastol);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
Bool SCIPsetIsFeasZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->feastol);
}

/** checks, if value is greater than feasibility tolerance */
Bool SCIPsetIsFeasPositive(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSP(val, set->feastol);
}

/** checks, if value is lower than -feasibility tolerance */
Bool SCIPsetIsFeasNegative(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSN(val, set->feastol);
}

/** checks, if the first given lower bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
Bool SCIPsetIsLbBetter(
   const SET*       set,                /**< global SCIP settings */
   Real             lb1,                /**< first lower bound to compare */
   Real             lb2                 /**< second lower bound to compare */
   )
{
   assert(set != NULL);

   return EPSGT(lb1, lb2, set->boundstreps);
}

/** checks, if the first given upper bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
Bool SCIPsetIsUbBetter(
   const SET*       set,                /**< global SCIP settings */
   Real             ub1,                /**< first upper bound to compare */
   Real             ub2                 /**< second upper bound to compare */
   )
{
   assert(set != NULL);

   return EPSLT(ub1, ub2, set->boundstreps);
}

/** checks, if the cut's activity is more then cutvioleps larger than the given right hand side;
 *  both, the activity and the rhs, should be normed
 */
Bool SCIPsetIsCutViolated(
   const SET*       set,                /**< global SCIP settings */
   Bool             root,               /**< should the root's cutvioleps be used? */
   Real             cutactivity,        /**< activity of the cut */
   Real             cutrhs              /**< right hand side value of the cut */
   )
{
   assert(set != NULL);

   if( root )
      return EPSGT(cutactivity, cutrhs, set->cutviolepsroot);
   else
      return EPSGT(cutactivity, cutrhs, set->cutvioleps);
}

/** checks, if relative difference of values is in range of epsilon */
Bool SCIPsetIsRelEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSZ(diff, set->epsilon);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
Bool SCIPsetIsRelLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSN(diff, set->epsilon);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
Bool SCIPsetIsRelLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return !EPSP(diff, set->epsilon);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
Bool SCIPsetIsRelGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSP(diff, set->epsilon);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
Bool SCIPsetIsRelGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return !EPSN(diff, set->epsilon);
}

/** checks, if relative difference of values is in range of sumepsilon */
Bool SCIPsetIsSumRelEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSZ(diff, set->sumepsilon);
}

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
Bool SCIPsetIsSumRelLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSN(diff, set->sumepsilon);
}

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
Bool SCIPsetIsSumRelLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return !EPSP(diff, set->sumepsilon);
}

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
Bool SCIPsetIsSumRelGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return EPSP(diff, set->sumepsilon);
}

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
Bool SCIPsetIsSumRelGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   Real diff;

   assert(set != NULL);

   diff = SCIPsetRelDiff(set, val1, val2);

   return !EPSN(diff, set->sumepsilon);
}

/** checks, if value is (positive) infinite */
Bool SCIPsetIsInfinity(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   )
{
   assert(set != NULL);

   return (val >= set->infinity);
}

/** checks, if value is non-negative within the LP feasibility bounds */
Bool SCIPsetIsFeasible(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return (val >= -set->feastol);
}

/** checks, if value is integral within the LP feasibility bounds */
Bool SCIPsetIsIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSISINT(val, set->feastol);
}

/** checks, if given fractional part is smaller than feastol */
Bool SCIPsetIsFracIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);
   assert(SCIPsetIsGE(set, val, -set->feastol));
   assert(SCIPsetIsLE(set, val, 1.0+set->feastol));

   return (val <= set->feastol);
}

/** rounds value + feasibility tolerance down to the next integer */
Real SCIPsetFloor(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSFLOOR(val, set->feastol);
}

/** rounds value - feasibility tolerance up to the next integer */
Real SCIPsetCeil(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return EPSCEIL(val, set->feastol);
}

/** returns fractional part of value, i.e. x - floor(x) */
Real SCIPsetFrac(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to return fractional part for */
   )
{
   assert(set != NULL);

   return EPSFRAC(val, set->feastol);
}


#endif

