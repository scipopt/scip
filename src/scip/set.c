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


/* LP Solving */

#define SCIP_DEFAULT_LPSOLVEFREQ         4 /**< frequency for solving LP at the nodes */
#define SCIP_DEFAULT_LPSOLVEDEPTH  INT_MAX /**< maximal depth for solving LPs */
#define SCIP_DEFAULT_COLAGELIMIT         8 /**< maximum age a column can reach before it is deleted from the LP */
#define SCIP_DEFAULT_ROWAGELIMIT         8 /**< maximum age a row can reach before it is deleted from the LP */


/* Pricing */

#define SCIP_DEFAULT_USEPRICING       TRUE /**< activate pricing of variables */
#define SCIP_DEFAULT_MAXPRICEVARS      128 /**< maximal number of variables priced in per pricing round */
#define SCIP_DEFAULT_MAXPRICEVARSROOT 2048 /**< maximal number of priced variables at the root node */
#define SCIP_DEFAULT_ABORTPRICEVARSFAC 2.0 /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
#define SCIP_DEFAULT_CLEANUPCOLS     FALSE /**< should new non-basic columns be removed after LP solving? */


/* Cut Separation */

#define SCIP_DEFAULT_MAXSEPACUTS       256 /**< maximal number of cuts separated per separation round */
#define SCIP_DEFAULT_MAXSEPACUTSROOT  4092 /**< maximal separated cuts at the root node */
#define SCIP_DEFAULT_CUTAGELIMIT       128 /**< maximum age a cut can reach before it is deleted from global cut pool */
#define SCIP_DEFAULT_CLEANUPROWS      TRUE /**< should new basic rows be removed after LP solving? */


/* Constraint Settings */

#define SCIP_DEFAULT_CONSAGELIMIT      128 /**< maximum age an unnecessary constraint can reach before it is deleted */


/* Primal Solutions */

#define SCIP_DEFAULT_MAXSOL            256 /**< maximal number of solutions to store in the solution storage */


/* Tree */

#define SCIP_DEFAULT_NODELIMIT LONGINT_MAX /**< maximal number of nodes to create */


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

   /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
   size = initsize;
   while( size < num )
      size = (int)(growfac * size + 1);

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
   (*set)->infinity = SCIP_DEFAULT_INFINITY;
   (*set)->epsilon = SCIP_DEFAULT_EPSILON;
   (*set)->sumepsilon = SCIP_DEFAULT_SUMEPSILON;
   (*set)->feastol = SCIP_DEFAULT_FEASTOL;
   (*set)->cutvioleps = SCIP_DEFAULT_CUTVIOLEPS;
   (*set)->memgrowfac = SCIP_DEFAULT_MEMGROWFAC;
   (*set)->memgrowinit = SCIP_DEFAULT_MEMGROWINIT;
   (*set)->treegrowfac = SCIP_DEFAULT_TREEGROWFAC;
   (*set)->treegrowinit = SCIP_DEFAULT_TREEGROWINIT;
   (*set)->pathgrowfac = SCIP_DEFAULT_PATHGROWFAC;
   (*set)->pathgrowinit = SCIP_DEFAULT_PATHGROWINIT;
   (*set)->branchscorefac = SCIP_DEFAULT_BRANCHSCOREFAC;

   CHECK_OKAY( SCIPbufferCreate(&(*set)->buffer) );

   (*set)->readers = NULL;
   (*set)->nreaders = 0;
   (*set)->readerssize = 0;
   (*set)->conshdlrs = NULL;
   (*set)->nconshdlrs = 0;
   (*set)->conshdlrssize = 0;
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
   (*set)->dispwidth = SCIP_DEFAULT_DISPWIDTH;
   (*set)->dispfreq = SCIP_DEFAULT_DISPFREQ;
   (*set)->dispheaderfreq = SCIP_DEFAULT_DISPHEADERFREQ;
   (*set)->maxpricevars = SCIP_DEFAULT_MAXPRICEVARS;
   (*set)->maxpricevarsroot = SCIP_DEFAULT_MAXPRICEVARSROOT;
   (*set)->abortpricevarsfac = SCIP_DEFAULT_ABORTPRICEVARSFAC;
   (*set)->maxsepacuts = SCIP_DEFAULT_MAXSEPACUTS;
   (*set)->maxsepacutsroot = SCIP_DEFAULT_MAXSEPACUTSROOT;
   (*set)->colagelimit = SCIP_DEFAULT_COLAGELIMIT;
   (*set)->rowagelimit = SCIP_DEFAULT_ROWAGELIMIT;
   (*set)->cutagelimit = SCIP_DEFAULT_CUTAGELIMIT;
   (*set)->consagelimit = SCIP_DEFAULT_CONSAGELIMIT;
   (*set)->maxsol = SCIP_DEFAULT_MAXSOL;
   (*set)->nodelimit = SCIP_DEFAULT_NODELIMIT;
   (*set)->lpsolvefreq = SCIP_DEFAULT_LPSOLVEFREQ;
   (*set)->lpsolvedepth = SCIP_DEFAULT_LPSOLVEDEPTH;
   (*set)->usepricing = SCIP_DEFAULT_USEPRICING;
   (*set)->cleanupcols = SCIP_DEFAULT_CLEANUPCOLS;
   (*set)->cleanuprows = SCIP_DEFAULT_CLEANUPROWS;

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

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrInit(set->conshdlrs[i], set->scip) );
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

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrExit(set->conshdlrs[i], set->scip) );
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

/** calculates the greatest common divisor of the two given values */
Longint SCIPsetGreComDiv(
   const SET*       set,                /**< global SCIP settings */
   Longint          val1,               /**< first value of greatest common devisor calculation */
   Longint          val2                /**< second value of greatest common devisor calculation */
   )
{
   Longint t;
   Longint gcd;

   assert(val1 >= 0);
   assert(val2 >= 0);

   /* extract all prime factors 2 */
   gcd = 1;
   while( !(val1 & 1) && !(val2 & 1) )
   {
      val1 /= 2;
      val2 /= 2;
      gcd *= 2;
   }

   t = val1 & 1 ? -val2 : val1;
   do
   {
      while( !(t & 1) )
	 t /= 2;

      if( t > 0 )
	 val1 = t;
      else
	 val2 = -t;

      t = val1 - val2;
   }
   while( t != 0 );
   gcd *= val1;

   return gcd;
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

