/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   set.c
 * @brief  global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "memory.h"
#include "tree.h"
#include "set.h"


/*
 * Default settings
 */

/* Message Output */

#define SCIP_DEFAULT_VERBLEVEL    SCIP_VERBLEVEL_NORMAL


/* Dynamic Memory */

#define SCIP_DEFAULT_MEMGROWFAC        1.2
#define SCIP_DEFAULT_MEMGROWINIT         4
#define SCIP_DEFAULT_BUFGROWFAC        2.0
#define SCIP_DEFAULT_BUFGROWINIT     65536
#define SCIP_DEFAULT_TREEGROWFAC       2.0
#define SCIP_DEFAULT_TREEGROWINIT    65536
#define SCIP_DEFAULT_PATHGROWFAC       2.0
#define SCIP_DEFAULT_PATHGROWINIT      256


/* Branching */
#define SCIP_DEFAULT_BRANCHSCOREFAC  0.167 /**< branching score factor to weigh downward and upward gain prediction */


/* Presolving */
#define SCIP_DEFAULT_MAXPRESOLROUNDS    -1 /**< maximal number of presolving rounds (-1: unlimited) */
#define SCIP_DEFAULT_PRESOLABORTFAC  1e-04 /**< abort presolve, if l.t. frac of the problem was changed in last round */


/* LP Solving */

#define SCIP_DEFAULT_LPSOLVEFREQ         4 /**< frequency for solving LP at the nodes; -1: never; 0: only root LP */
#define SCIP_DEFAULT_LPSOLVEDEPTH       -1 /**< maximal depth for solving LPs (-1: no depth limit) */
#define SCIP_DEFAULT_COLAGELIMIT         8 /**< maximum age a column can reach before it is deleted from the LP */
#define SCIP_DEFAULT_ROWAGELIMIT         8 /**< maximum age a row can reach before it is deleted from the LP */


/* Pricing */

#define SCIP_DEFAULT_MAXPRICEVARS      128 /**< maximal number of variables priced in per pricing round */
#define SCIP_DEFAULT_MAXPRICEVARSROOT 2048 /**< maximal number of priced variables at the root node */
#define SCIP_DEFAULT_ABORTPRICEVARSFAC 2.0 /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
#define SCIP_DEFAULT_CLEANUPCOLS     FALSE /**< should new non-basic columns be removed after LP solving? */


/* Cut Separation */

#define SCIP_DEFAULT_MAXSEPACUTS       128 /**< maximal number of cuts separated per separation round */
#define SCIP_DEFAULT_MAXSEPACUTSROOT  2048 /**< maximal separated cuts at the root node */
#define SCIP_DEFAULT_CUTAGELIMIT       128 /**< maximum age a cut can reach before it is deleted from global cut pool */
#define SCIP_DEFAULT_CLEANUPROWS      TRUE /**< should new basic rows be removed after LP solving? */


/* Constraint Settings */

#define SCIP_DEFAULT_CONSAGELIMIT      128 /**< maximum age an unnecessary constraint can reach before it is deleted */


/* Primal Solutions */

#define SCIP_DEFAULT_MAXSOL            256 /**< maximal number of solutions to store in the solution storage */


/* Tree */

#define SCIP_DEFAULT_NODELIMIT          -1 /**< maximal number of nodes to process (-1: no limit) */


/* Display */

#define SCIP_DEFAULT_DISPWIDTH         138 /**< maximal number of characters in a node information line */
#define SCIP_DEFAULT_DISPFREQ         1000 /**< frequency for displaying node information lines */
#define SCIP_DEFAULT_DISPHEADERFREQ     15 /**< frequency for displaying header lines (every n'th node information line) */




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
{
   Real newfeastol;

   newfeastol = SCIPparamGetReal(param);
   
   /* change the feastol through the SCIP call in order to mark the LP unsolved */
   CHECK_OKAY( SCIPsetFeastol(scip, newfeastol) );

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
   (*set)->pricerssize = 0;
   (*set)->conshdlrs = NULL;
   (*set)->nconshdlrs = 0;
   (*set)->conshdlrssize = 0;
   (*set)->presols = NULL;
   (*set)->npresols = 0;
   (*set)->presolssize = 0;
   (*set)->sepas = NULL;
   (*set)->nsepas = 0;
   (*set)->sepassize = 0;
   (*set)->heurs = NULL;
   (*set)->nheurs = 0;
   (*set)->heurssize = 0;
   (*set)->eventhdlrs = NULL;
   (*set)->neventhdlrs = 0;
   (*set)->eventhdlrssize = 0;
   (*set)->nodesels = NULL;
   (*set)->nnodesels = 0;
   (*set)->nodeselssize = 0;
   (*set)->nodesel = NULL;
   (*set)->branchrules = NULL;
   (*set)->nbranchrules = 0;
   (*set)->branchrulessize = 0;
   (*set)->disps = NULL;
   (*set)->ndisps = 0;
   (*set)->dispssize = 0;

   /* SCIP parameters */
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_VerbLevel", "verbosity level of output",
                  &(*set)->verblevel, SCIP_DEFAULT_VERBLEVEL, SCIP_VERBLEVEL_NONE, SCIP_VERBLEVEL_FULL,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_Infinity", "values larger than this are considered infinity",
                  &(*set)->infinity, SCIP_DEFAULT_INFINITY, 1e+10, SCIP_INVALID/10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_Epsilon", "absolute values smaller than this are considered zero",
                  &(*set)->epsilon, SCIP_DEFAULT_EPSILON, machineeps, SCIP_MAXEPSILON,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_SumEpsilon", "absolute values of sums smaller than this are considered zero",
                  &(*set)->sumepsilon, SCIP_DEFAULT_SUMEPSILON, machineeps*1e+03, SCIP_MAXEPSILON,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_Feastol", "LP feasibility tolerance",
                  &(*set)->feastol, SCIP_DEFAULT_FEASTOL, machineeps*1e+03, SCIP_MAXEPSILON,
                  paramChgdFeastol, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_CutViolEps", "epsilon for deciding if a cut is violated",
                  &(*set)->cutvioleps, SCIP_DEFAULT_CUTVIOLEPS, machineeps*1e+03, SCIP_MAXEPSILON,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_MemGrowFac", "memory growing factor for dynamically allocated arrays",
                  &(*set)->memgrowfac, SCIP_DEFAULT_MEMGROWFAC, 1.0, 10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_MemGrowInit", "initial size of dynamically allocated arrays",
                  &(*set)->memgrowinit, SCIP_DEFAULT_MEMGROWINIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_TreeGrowFac", "memory growing factor for tree array",
                  &(*set)->treegrowfac, SCIP_DEFAULT_TREEGROWFAC, 1.0, 10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_TreeGrowInit", "initial size of tree array",
                  &(*set)->treegrowinit, SCIP_DEFAULT_TREEGROWINIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "global_PathGrowFac", "memory growing factor for path array",
                  &(*set)->pathgrowfac, SCIP_DEFAULT_PATHGROWFAC, 1.0, 10.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_PathGrowInit", "initial size of path array",
                  &(*set)->pathgrowinit, SCIP_DEFAULT_PATHGROWINIT, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "global_BranchScoreFac", "branching score factor to weigh downward and upward gain prediction",
                  &(*set)->branchscorefac, SCIP_DEFAULT_BRANCHSCOREFAC, 0.0, 1.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "global_DispWidth", "maximal number of characters in a node information line",
                  &(*set)->dispwidth, SCIP_DEFAULT_DISPWIDTH, 0, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "global_DispFreq", "frequency for displaying node information lines",
                  &(*set)->dispfreq, SCIP_DEFAULT_DISPFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "global_DispHeaderFreq", "frequency for displaying header lines (every n'th node information line)",
                  &(*set)->dispheaderfreq, SCIP_DEFAULT_DISPHEADERFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "global_MaxPresolRounds", "maximal number of presolving rounds (-1: unlimited)",
                  &(*set)->maxpresolrounds, SCIP_DEFAULT_MAXPRESOLROUNDS, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr, 
                  "global_PresolAbortFac",
                  "abort presolve, if less than this fraction of the problem was changed in last presolve round",
                  &(*set)->presolabortfac, SCIP_DEFAULT_PRESOLABORTFAC, 0.0, 1.0,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr, 
                  "global_MaxPriceVars", "maximal number of variables priced in per pricing round",
                  &(*set)->maxpricevars, SCIP_DEFAULT_MAXPRICEVARS, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_MaxPriceVarsRoot", "maximal number of priced variables at the root node",
                  &(*set)->maxpricevarsroot, SCIP_DEFAULT_MAXPRICEVARSROOT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddRealParam(*set, memhdr,
                  "global_AbortPriceVarsFac", "pricing is aborted, if fac * maxpricevars pricing candidates were found",
                  &(*set)->abortpricevarsfac, SCIP_DEFAULT_ABORTPRICEVARSFAC, 1.0, REAL_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_MaxSepaCuts", "maximal number of cuts separated per separation round",
                  &(*set)->maxsepacuts, SCIP_DEFAULT_MAXSEPACUTS, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_MaxSepaCutsRoot", "maximal number of separated cuts at the root node",
                  &(*set)->maxsepacutsroot, SCIP_DEFAULT_MAXSEPACUTSROOT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_ColAgeLimit", "maximum age a column can reach before it is deleted from the LP",
                  &(*set)->colagelimit, SCIP_DEFAULT_COLAGELIMIT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_RowAgeLimit", "maximum age a row can reach before it is deleted from the LP",
                  &(*set)->rowagelimit, SCIP_DEFAULT_ROWAGELIMIT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_CutAgeLimit", "maximum age a cut can reach before it is deleted from the global cut pool",
                  &(*set)->cutagelimit, SCIP_DEFAULT_CUTAGELIMIT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_ConsAgeLimit", "maximum age an unnecessary constraint can reach before it is deleted",
                  &(*set)->consagelimit, SCIP_DEFAULT_CONSAGELIMIT, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_MaxSol", "maximal number of solutions to store in the solution storage",
                  &(*set)->maxsol, SCIP_DEFAULT_MAXSOL, 1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddLongintParam(*set, memhdr,
                  "global_NodeLimit", "maximal number of nodes to process (-1: no limit)",
                  &(*set)->nodelimit, SCIP_DEFAULT_NODELIMIT, -1, LONGINT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_LpSolveFreq", "frequency for solving LP at the nodes (-1: never; 0: only root LP)",
                  &(*set)->lpsolvefreq, SCIP_DEFAULT_LPSOLVEFREQ, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddIntParam(*set, memhdr,
                  "global_LpSolveDepth", "maximal depth for solving LP at the nodes (-1: no depth limit)",
                  &(*set)->lpsolvedepth, SCIP_DEFAULT_LPSOLVEDEPTH, -1, INT_MAX,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "global_CleanupCols", "should new non-basic columns be removed after LP solving?",
                  &(*set)->cleanupcols, SCIP_DEFAULT_CLEANUPCOLS,
                  NULL, NULL) );
   CHECK_OKAY( SCIPsetAddBoolParam(*set, memhdr,
                  "global_CleanupRows", "should new basic rows be removed after LP solving?",
                  &(*set)->cleanuprows, SCIP_DEFAULT_CLEANUPROWS,
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
      errorMessage("pricers not implemented yet");
      abort();
#if 0
      CHECK_OKAY( SCIPpricerFree(&(*set)->pricers[i], (*set)->scip) );
#endif
   }
   freeMemoryArrayNull(&(*set)->pricers);

   /* free constraint handlers */
   for( i = 0; i < (*set)->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrFree(&(*set)->conshdlrs[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->conshdlrs);

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
   Bool             comments            /**< should parameter descriptions be written as comments? */
   )
{
   assert(set != NULL);

   CHECK_OKAY( SCIPparamsetWrite(set->paramset, filename, comments) );

   return SCIP_OKAY;
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
      errorMessage("pricers not yet implemented");
      abort();
#if 0
      if( strcmp(SCIPpricerGetName(set->pricers[i]), name) == 0 )
         return set->pricers[i];
#endif
   }

   return NULL;
}

/** inserts constraint handler in constraint handler list */
RETCODE SCIPsetIncludeConsHdlr(
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
CONSHDLR* SCIPsetFindConsHdlr(
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

/** inserts event handler in event handler list */
RETCODE SCIPsetIncludeEventHdlr(
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
EVENTHDLR* SCIPsetFindEventHdlr(
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
   assert(set != NULL);
   assert(nodesel != NULL);
   assert(!SCIPnodeselIsInitialized(nodesel));

   if( set->nnodesels >= set->nodeselssize )
   {
      set->nodeselssize = SCIPsetCalcMemGrowSize(set, set->nnodesels+1);
      ALLOC_OKAY( reallocMemoryArray(&set->nodesels, set->nodeselssize) );
   }
   assert(set->nnodesels < set->nodeselssize);
   
   set->nodesels[set->nnodesels] = nodesel;
   set->nnodesels++;

   if( set->nodesel == NULL )
      set->nodesel = nodesel;

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

/** inserts branching rule in branching rule list */
RETCODE SCIPsetIncludeBranchrule(
   SET*             set,                /**< global SCIP settings */
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   int i;

   assert(set != NULL);
   assert(branchrule != NULL);
   assert(!SCIPbranchruleIsInitialized(branchrule));

   if( set->nbranchrules >= set->branchrulessize )
   {
      set->branchrulessize = SCIPsetCalcMemGrowSize(set, set->nbranchrules+1);
      ALLOC_OKAY( reallocMemoryArray(&set->branchrules, set->branchrulessize) );
   }
   assert(set->nbranchrules < set->branchrulessize);

   for( i = set->nbranchrules; i > 0
           && SCIPbranchruleGetPriority(branchrule) > SCIPbranchruleGetPriority(set->branchrules[i-1]); --i )
   {
      set->branchrules[i] = set->branchrules[i-1];
   }
   set->branchrules[i] = branchrule;
   set->nbranchrules++;

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
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* variable pricers */
   for( i = 0; i < set->npricers; ++i )
   {
      errorMessage("pricers not implemented yet");
      abort();
#if 0
      CHECK_OKAY( SCIPpricerInit(set->pricers[i], set->scip) );
#endif
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrInit(set->conshdlrs[i], set->scip) );
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
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* variable pricers */
   for( i = 0; i < set->npricers; ++i )
   {
      errorMessage("pricers not implemented yet");
      abort();
#if 0
      CHECK_OKAY( SCIPpricerExit(set->pricers[i], set->scip) );
#endif
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrExit(set->conshdlrs[i], set->scip) );
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
      char s[MAXSTRLEN];
      sprintf(s, "Invalid verbosity level <%d>, maximum is <%d>", verblevel, SCIP_VERBLEVEL_FULL);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }
   
   set->verblevel = verblevel;

   return SCIP_OKAY;
}

/** sets LP feasibility tolerance */
RETCODE SCIPsetSetFeastol(
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data (or NULL) */
   Real             feastol             /**< new feasibility tolerance */
   )
{
   assert(set != NULL);

   set->feastol = feastol;
   if( lp != NULL )
   {
      CHECK_OKAY( SCIPlpSetFeastol(lp, feastol) );
   }

   return SCIP_OKAY;
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

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
Bool SCIPsetRealToRational(
   const SET*       set,                /**< global SCIP settings */
   Real             val,                /**< real value to convert into rational number */
   Longint          maxdnom,            /**< maximal denominator allowed */
   Longint*         nominator,          /**< pointer to store the nominator of the rational number */
   Longint*         denominator         /**< pointer to store the denominator of the rational number */
   )
{
   Real a;
   Real b;
   Real g0;
   Real g1;
   Real gx;
   Real h0;
   Real h1;
   Real hx;

   assert(nominator != NULL);
   assert(denominator != NULL);

   b = val;
   a = SCIPsetFloor(set, b);
   g0 = a;
   g1 = 1.0;
   h0 = 1.0;
   h1 = 0.0;

   while( !SCIPsetIsFeasZero(set, val - g0/h0) )
   {
      assert(SCIPsetIsGT(set, b, a));
      assert(h0 >= 0.0);
      assert(h1 >= 0.0);

      b = 1.0 / (b - a);
      a = SCIPsetFloor(set, b);

      assert(a >= 0.0);

      gx = g0;
      hx = h0;

      g0 = a * g0 + g1;
      h0 = a * h0 + h1;

      g1 = gx;
      h1 = hx;
      
      if( h0 > maxdnom )
         return FALSE;
   }

   if( ABS(g0) > (LONGINT_MAX << 4) || h0 > (LONGINT_MAX << 4) )
      return FALSE;

   assert(h0 >= 0.0);

   *nominator = (Longint)g0;
   *denominator = (Longint)h0;

   return TRUE;
}

   

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

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

/** checks, if the cut's activity is more then cutvioleps larger than the given right hand side;
 *  both, the activity and the rhs, should be normed
 */
Bool SCIPsetIsCutViolated(
   const SET*       set,                /**< global SCIP settings */
   Real             cutactivity,        /**< activity of the cut */
   Real             cutrhs              /**< right hand side value of the cut */
   )
{
   assert(set != NULL);

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
   assert(val >= -set->feastol);
   assert(val < 1.0);

   return (val <= set->feastol);
}


#endif

