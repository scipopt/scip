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


/* LP solving */

#define SCIP_DEFAULT_LPSOLVEFREQ         3 /**< frequency for solving LP at the nodes */
#define SCIP_DEFAULT_LPSOLVEDEPTH  INT_MAX /**< maximal depth for solving LPs */


/* Pricing */

#define SCIP_DEFAULT_USEPRICING       TRUE /**< activate pricing of variables */
#define SCIP_DEFAULT_MAXPRICEVARS       32 /**< maximal number of variables priced in per pricing round */
#define SCIP_DEFAULT_MAXPRICEVARSROOT 1024 /**< maximal number of priced variables at the root node */
#define SCIP_DEFAULT_ABORTPRICEVARSFAC 2.0 /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */


/* Cut Separation */

#define SCIP_DEFAULT_MAXSEPACUTS       128 /**< maximal number of cuts separated per separation round */
#define SCIP_DEFAULT_MAXSEPACUTSROOT  4092 /**< maximal separated cuts at the root node */
#define SCIP_DEFAULT_AGELIMIT          128 /**< maximum age a cut can reach before it is deleted from global cut pool */


/* Primal Solutions */

#define SCIP_DEFAULT_MAXSOL            256 /**< maximal number of solutions to store in the solution storage */


/* Tree */

#define SCIP_DEFAULT_NODELIMIT /*LONGINT_MAX*/ 10000 /**< maximal number of nodes to create */


/* Display */

#define SCIP_DEFAULT_DISPWIDTH         140 /**< maximal number of characters in a node information line */
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

   /* calculate the size with this loop, such that the resulting numbers are allways the same (-> block memory) */
   size = initsize;
   while( size < num )
      size = growfac * size + 1;

   return size;
}



/** creates global SCIP settings */
RETCODE SCIPsetCreate(
   SET**            set,                /**< pointer to SCIP settings */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(set != NULL);
   assert(scip != NULL);

   ALLOC_OKAY( allocMemory(set) );

   (*set)->scip = scip;
   (*set)->verblevel = SCIP_DEFAULT_VERBLEVEL;
   (*set)->epsilon = SCIP_DEFAULT_EPSILON;
   (*set)->sumepsilon = SCIP_DEFAULT_SUMEPSILON;
   (*set)->infinity = SCIP_DEFAULT_INFINITY;
   (*set)->feastol = SCIP_DEFAULT_FEASTOL;
   (*set)->memGrowFac = SCIP_DEFAULT_MEMGROWFAC;
   (*set)->memGrowInit = SCIP_DEFAULT_MEMGROWINIT;
   (*set)->treeGrowFac = SCIP_DEFAULT_TREEGROWFAC;
   (*set)->treeGrowInit = SCIP_DEFAULT_TREEGROWINIT;
   (*set)->pathGrowFac = SCIP_DEFAULT_PATHGROWFAC;
   (*set)->pathGrowInit = SCIP_DEFAULT_PATHGROWINIT;

   CHECK_OKAY( SCIPbufferCreate(&(*set)->buffer) );

   (*set)->readers = NULL;
   (*set)->nreaders = 0;
   (*set)->readerssize = 0;
   (*set)->conshdlrs = NULL;
   (*set)->nconshdlrs = 0;
   (*set)->conshdlrssize = 0;
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
   (*set)->dispwidth = SCIP_DEFAULT_DISPWIDTH;
   (*set)->dispfreq = SCIP_DEFAULT_DISPFREQ;
   (*set)->dispheaderfreq = SCIP_DEFAULT_DISPHEADERFREQ;
   (*set)->maxpricevars = SCIP_DEFAULT_MAXPRICEVARS;
   (*set)->maxpricevarsroot = SCIP_DEFAULT_MAXPRICEVARSROOT;
   (*set)->abortpricevarsfac = SCIP_DEFAULT_ABORTPRICEVARSFAC;
   (*set)->maxsepacuts = SCIP_DEFAULT_MAXSEPACUTS;
   (*set)->maxsepacutsroot = SCIP_DEFAULT_MAXSEPACUTSROOT;
   (*set)->agelimit = SCIP_DEFAULT_AGELIMIT;
   (*set)->maxsol = SCIP_DEFAULT_MAXSOL;
   (*set)->nodelimit = SCIP_DEFAULT_NODELIMIT;
   (*set)->lpsolvefreq = SCIP_DEFAULT_LPSOLVEFREQ;
   (*set)->lpsolvedepth = SCIP_DEFAULT_LPSOLVEDEPTH;
   (*set)->usepricing = SCIP_DEFAULT_USEPRICING;

   return SCIP_OKAY;
}

/** frees global SCIP settings */
RETCODE SCIPsetFree(
   SET**            set                 /**< pointer to SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* free memory buffers */
   SCIPbufferFree(&(*set)->buffer);

   /* free file readers */
   for( i = 0; i < (*set)->nreaders; ++i )
   {
      CHECK_OKAY( SCIPreaderFree(&(*set)->readers[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->readers);

   /* free constraint handlers */
   for( i = 0; i < (*set)->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrFree(&(*set)->conshdlrs[i], (*set)->scip) );
   }
   freeMemoryArrayNull(&(*set)->conshdlrs);

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

/** finds the file reader of the given name */
RETCODE SCIPsetFindReader(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of file reader */
   READER**         reader              /**< pointer for storing the file reader (returns NULL, if not found) */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);
   assert(reader != NULL);

   *reader = NULL;
   for( i = 0; i < set->nreaders; ++i )
   {
      if( strcmp(SCIPreaderGetName(set->readers[i]), name) == 0 )
      {
         *reader = set->readers[i];
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** inserts constraint handler in constraint handler list */
RETCODE SCIPsetIncludeConsHdlr(
   SET*             set,                /**< global SCIP settings */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   int chckpriority;
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

   chckpriority = SCIPconshdlrGetChckPriority(conshdlr);
   for( i = set->nconshdlrs; i > 0 && SCIPconshdlrGetChckPriority(set->conshdlrs[i-1]) < chckpriority; --i )
   {
      set->conshdlrs[i] = set->conshdlrs[i-1];
   }

   set->conshdlrs[i] = conshdlr;
   set->nconshdlrs++;

   return SCIP_OKAY;
}   

/** finds the constraint handler of the given name */
RETCODE SCIPsetFindConsHdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   *conshdlr = NULL;
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      if( strcmp(SCIPconshdlrGetName(set->conshdlrs[i]), name) == 0 )
      {
         *conshdlr = set->conshdlrs[i];
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
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

/** finds the primal heuristic of the given name */
RETCODE SCIPsetFindHeur(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of primal heuristic */
   HEUR**           heur                /**< pointer for storing the primal heuristic (returns NULL, if not found) */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);
   assert(heur != NULL);

   *heur = NULL;
   for( i = 0; i < set->nheurs; ++i )
   {
      if( strcmp(SCIPheurGetName(set->heurs[i]), name) == 0 )
      {
         *heur = set->heurs[i];
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
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

/** finds the event handler of the given name */
RETCODE SCIPsetFindEventHdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of event handler */
   EVENTHDLR**      eventhdlr           /**< pointer for storing the event handler (returns NULL, if not found) */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);
   assert(eventhdlr != NULL);

   *eventhdlr = NULL;
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      if( strcmp(SCIPeventhdlrGetName(set->eventhdlrs[i]), name) == 0 )
      {
         *eventhdlr = set->eventhdlrs[i];
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
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

/** initializes all user callback functions */
RETCODE SCIPsetInitCallbacks(
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrInit(set->conshdlrs[i], set->scip) );
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

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrExit(set->conshdlrs[i], set->scip) );
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
   return calcGrowSize(set->memGrowInit, set->memGrowFac, num);
}

/** calculate memory size for tree array */
int SCIPsetCalcTreeGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->treeGrowInit, set->treeGrowFac, num);
}

/** calculate memory size for path array */
int SCIPsetCalcPathGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->pathGrowInit, set->pathGrowFac, num);
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
      char s[255];
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

/** checks, if values are in range of feasibility tolerance */
Bool SCIPsetIsFeasEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSEQ(val1, val2, set->feastol);
}

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
Bool SCIPsetIsFeasLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSLT(val1, val2, set->feastol);
}

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
Bool SCIPsetIsFeasLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSLE(val1, val2, set->feastol);
}

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
Bool SCIPsetIsFeasGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSGT(val1, val2, set->feastol);
}

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
Bool SCIPsetIsFeasGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   return EPSGE(val1, val2, set->feastol);
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

/** rounds value down to the next integer */
Real SCIPsetFloor(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return floor(val + set->feastol);
}

/** rounds value up to the next integer */
Real SCIPsetCeil(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return ceil(val - set->feastol);
}

/** returns fractional part of value, i.e. x - floor(x) */
Real SCIPsetFrac(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to return fractional part for */
   )
{
   assert(set != NULL);

   return val - SCIPsetFloor(set, val);
}

/** checks, if value is integral within the LP feasibility bounds */
Bool SCIPsetIsIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(set != NULL);

   return (SCIPsetCeil(set, val) - val <= set->feastol);
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

/** checks, if the given integer bounds correspond to a fixed interval */
Bool SCIPsetIsFixed(
   const SET*       set,                /**< global SCIP settings */
   Real             lb,                 /**< lower integer bound */
   Real             ub                  /**< upper integer bound */
   )
{
   assert(set != NULL);
   assert(SCIPsetIsIntegral(set, lb));
   assert(SCIPsetIsIntegral(set, ub));

   return SCIPsetIsEQ(set, lb, ub);
}


#endif

