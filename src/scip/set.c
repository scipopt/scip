/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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


static
int calcGrowSize(                       /**< calculate memory size for dynamically allocated arrays */
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



RETCODE SCIPsetCreate(                  /**< creates global SCIP settings */
   SET**            set,                /**< pointer to SCIP settings */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(set != NULL);
   assert(scip != NULL);

   ALLOC_OKAY( allocMemory(*set) );

   (*set)->scip = scip;
   (*set)->verblevel = SCIP_DEFAULT_VERBLEVEL;
   (*set)->epsilon = SCIP_DEFAULT_EPSILON;
   (*set)->infinity = SCIP_DEFAULT_INFINITY;
   (*set)->feastol = SCIP_DEFAULT_FEASTOL;
   (*set)->memGrowFac = SCIP_DEFAULT_MEMGROWFAC;
   (*set)->memGrowInit = SCIP_DEFAULT_MEMGROWINIT;
   (*set)->treeGrowFac = SCIP_DEFAULT_TREEGROWFAC;
   (*set)->treeGrowInit = SCIP_DEFAULT_TREEGROWINIT;
   (*set)->pathGrowFac = SCIP_DEFAULT_PATHGROWFAC;
   (*set)->pathGrowInit = SCIP_DEFAULT_PATHGROWINIT;
   (*set)->readers = NULL;
   (*set)->nreaders = 0;
   (*set)->readerssize = 0;
   (*set)->conshdlrs = NULL;
   (*set)->nconshdlrs = 0;
   (*set)->conshdlrssize = 0;
   (*set)->nodesels = NULL;
   (*set)->nnodesels = 0;
   (*set)->nodeselssize = 0;
   (*set)->nodesel = NULL;
   (*set)->disps = NULL;
   (*set)->ndisps = 0;
   (*set)->dispssize = 0;
   (*set)->dispwidth = SCIP_DEFAULT_DISPWIDTH;
   (*set)->dispfreq = SCIP_DEFAULT_DISPFREQ;
   (*set)->dispheaderfreq = SCIP_DEFAULT_DISPHEADERFREQ;
   (*set)->maxpricevars = SCIP_DEFAULT_MAXPRICEVARS;
   (*set)->maxsepacuts = SCIP_DEFAULT_MAXSEPACUTS;
   (*set)->maxsol = SCIP_DEFAULT_MAXSOL;
   (*set)->nodelimit = SCIP_DEFAULT_NODELIMIT;
   (*set)->usepricing = SCIP_DEFAULT_USEPRICING;

   return SCIP_OKAY;
}

RETCODE SCIPsetFree(                    /**< frees global SCIP settings */
   SET**            set                 /**< pointer to SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* free file readers */
   for( i = 0; i < (*set)->nreaders; ++i )
   {
      CHECK_OKAY( SCIPreaderFree(&(*set)->readers[i]) );
   }
   freeMemoryArray((*set)->readers);

   /* free constraint handlers */
   for( i = 0; i < (*set)->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrFree(&(*set)->conshdlrs[i]) );
   }
   freeMemoryArray((*set)->conshdlrs);

   /* free node selectors */
   for( i = 0; i < (*set)->nnodesels; ++i )
   {
      CHECK_OKAY( SCIPnodeselFree(&(*set)->nodesels[i]) );
   }
   freeMemoryArray((*set)->nodesels);

   /* free display columns */
   for( i = 0; i < (*set)->ndisps; ++i )
   {
      CHECK_OKAY( SCIPdispFree(&(*set)->disps[i]) );
   }
   freeMemoryArray((*set)->disps);

   freeMemory(*set);

   return SCIP_OKAY;
}

RETCODE SCIPsetIncludeReader(           /**< inserts file reader in file reader list */
   SET*             set,                /**< global SCIP settings */
   READER*          reader              /**< file reader */
   )
{
   assert(set != NULL);
   assert(reader != NULL);
   assert(!SCIPreaderIsInitialized(reader));

   if( set->nreaders >= set->readerssize )
   {
      set->readerssize = SCIPsetCalcMemGrowSize(set, set->nreaders+1);
      ALLOC_OKAY( reallocMemoryArray(set->readers, set->readerssize) );
   }
   assert(set->nreaders < set->readerssize);
   
   set->readers[set->nreaders] = reader;
   set->nreaders++;

   return SCIP_OKAY;
}   

RETCODE SCIPsetFindReader(              /**< finds the file reader of the given name */
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

RETCODE SCIPsetIncludeConsHdlr(         /**< inserts constraint handler in constraint handler list */
   SET*             set,                /**< global SCIP settings */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(set != NULL);
   assert(conshdlr != NULL);
   assert(!SCIPconshdlrIsInitialized(conshdlr));

   if( set->nconshdlrs >= set->conshdlrssize )
   {
      set->conshdlrssize = SCIPsetCalcMemGrowSize(set, set->nconshdlrs+1);
      ALLOC_OKAY( reallocMemoryArray(set->conshdlrs, set->conshdlrssize) );
   }
   assert(set->nconshdlrs < set->conshdlrssize);
   
   set->conshdlrs[set->nconshdlrs] = conshdlr;
   set->nconshdlrs++;

   return SCIP_OKAY;
}   

RETCODE SCIPsetFindConsHdlr(            /**< finds the constraint handler of the given name */
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

RETCODE SCIPsetIncludeNodesel(          /**< inserts node selector in node selector list */
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
      ALLOC_OKAY( reallocMemoryArray(set->nodesels, set->nodeselssize) );
   }
   assert(set->nnodesels < set->nodeselssize);
   
   set->nodesels[set->nnodesels] = nodesel;
   set->nnodesels++;

   if( set->nodesel == NULL )
      set->nodesel = nodesel;

   return SCIP_OKAY;
}   

RETCODE SCIPsetIncludeDisp(             /**< inserts display column in display column list */
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
      ALLOC_OKAY( reallocMemoryArray(set->disps, set->dispssize) );
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

RETCODE SCIPsetInitCallbacks(           /**< initializes all user callback functions */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* file readers */
   for( i = 0; i < set->nreaders; ++i )
   {
      CHECK_OKAY( SCIPreaderInit(set->readers[i], set->scip) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrInit(set->conshdlrs[i], set->scip) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      CHECK_OKAY( SCIPnodeselInit(set->nodesels[i], set->scip) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      CHECK_OKAY( SCIPdispInit(set->disps[i], set->scip) );
   }
   CHECK_OKAY( SCIPdispAutoActivate(set) );

   return SCIP_OKAY;
}

RETCODE SCIPsetExitCallbacks(           /**< calls exit methods of all user callback functions */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   assert(set != NULL);

   /* file readers */
   for( i = 0; i < set->nreaders; ++i )
   {
      CHECK_OKAY( SCIPreaderExit(set->readers[i], set->scip) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrExit(set->conshdlrs[i], set->scip) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      CHECK_OKAY( SCIPnodeselExit(set->nodesels[i], set->scip) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      CHECK_OKAY( SCIPdispExit(set->disps[i], set->scip) );
   }

   return SCIP_OKAY;
}

int SCIPsetCalcMemGrowSize(             /**< calculate memory size for dynamically allocated arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->memGrowInit, set->memGrowFac, num);
}

int SCIPsetCalcTreeGrowSize(            /**< calculate memory size for tree array */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->treeGrowInit, set->treeGrowFac, num);
}

int SCIPsetCalcPathGrowSize(            /**< calculate memory size for path array */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->pathGrowInit, set->pathGrowFac, num);
}

RETCODE SCIPsetSetVerbLevel(            /**< sets verbosity level for message output */
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

RETCODE SCIPsetSetFeastol(              /**< sets LP feasibility tolerance */
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

Bool SCIPsetIsEQ(                       /**< checks, if values are in range of epsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( ABS(val1-val2) < set->epsilon );
}

Bool SCIPsetIsL(                        /**< checks, if val1 is (more than epsilon) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 < val2 - set->epsilon );
}

Bool SCIPsetIsLE(                       /**< checks, if val1 is not (more than epsilon) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 <= val2 + set->epsilon );
}

Bool SCIPsetIsG(                        /**< checks, if val1 is (more than epsilon) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 > val2 + set->epsilon );
}

Bool SCIPsetIsGE(                       /**< checks, if val1 is not (more than epsilon) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 >= val2 - set->epsilon );
}

Bool SCIPsetIsInfinity(                 /**< checks, if value is (positive) infinite */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   )
{
   return( val >= set->infinity );
}

Bool SCIPsetIsZero(                     /**< checks, if value is in range epsilon of 0.0 */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( ABS(val) <= set->epsilon );
}

Bool SCIPsetIsPos(                      /**< checks, if value is greater than epsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( val > set->epsilon );
}

Bool SCIPsetIsNeg(                      /**< checks, if value is lower than -epsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( val < -set->epsilon );
}

Bool SCIPsetIsFeasible(                 /**< checks, if value is non-negative within the LP feasibility bounds */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( val >= -set->feastol );
}

Real SCIPsetFloor(                      /**< rounds value down to the next integer */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return floor(val + set->feastol);
}

Real SCIPsetCeil(                       /**< rounds value up to the next integer */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return ceil(val - set->feastol);
}

Bool SCIPsetIsIntegral(                 /**< checks, if value is integral within the LP feasibility bounds */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
#if 0
   return (SCIPsetCeil(set, val) - SCIPsetFloor(set, val) < 0.5);
#else
   return (SCIPsetCeil(set, val) - val <= set->feastol);
#endif
}
