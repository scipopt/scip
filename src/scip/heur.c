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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur.c,v 1.39 2005/01/18 09:26:45 bzfpfend Exp $"

/**@file   heur.c
 * @brief  methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "clock.h"
#include "paramset.h"
#include "primal.h"
#include "scip.h"
#include "heur.h"

#include "struct_heur.h"



/** compares two heuristics w. r. to their delay positions and their priority */
DECL_SORTPTRCOMP(SCIPheurComp)
{  /*lint --e{715}*/
   HEUR* heur1 = (HEUR*)elem1;
   HEUR* heur2 = (HEUR*)elem2;

   assert(heur1 != NULL);
   assert(heur2 != NULL);

   if( heur1->delaypos == heur2->delaypos )
   {
      assert(heur1->delaypos == -1);
      assert(heur2->delaypos == -1);
      return heur2->priority - heur1->priority; /* prefer higher priorities */
   }
   else if( heur1->delaypos == -1 )
      return +1;                                /* prefer delayed heuristics */
   else if( heur2->delaypos == -1 )
      return -1;                                /* prefer delayed heuristics */
   else if( heur1->ncalls * heur1->freq != heur2->ncalls * heur2->freq )
      return heur1->ncalls * heur1->freq - heur2->ncalls * heur2->freq;
   else
      return heur1->delaypos - heur2->delaypos; /* prefer lower delay positions */
}

/** method to call, when the priority of a heuristic was changed */
static
DECL_PARAMCHGD(paramChgdHeurPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetHeurPriority() to mark the heurs unsorted */
   CHECK_OKAY( SCIPsetHeurPriority(scip, (HEUR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a primal heuristic */
RETCODE SCIPheurCreate(
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   char             dispchar,           /**< display character of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   int              freqofs,            /**< frequency offset for calling primal heuristic */
   int              maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   Bool             pseudonodes,        /**< call heuristic at nodes where only a pseudo solution exist? */
   Bool             duringplunging,     /**< call heuristic during plunging? */
   DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(heur != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(freqofs >= 0);
   assert(heurexec != NULL);

   ALLOC_OKAY( allocMemory(heur) );
   ALLOC_OKAY( duplicateMemoryArray(&(*heur)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*heur)->desc, desc, strlen(desc)+1) );
   (*heur)->dispchar = dispchar;
   (*heur)->priority = priority;
   (*heur)->freq = freq;
   (*heur)->freqofs = freqofs;
   (*heur)->maxdepth = maxdepth;
   (*heur)->delaypos = -1;
   (*heur)->pseudonodes = pseudonodes;
   (*heur)->duringplunging = duringplunging;
   (*heur)->heurfree = heurfree;
   (*heur)->heurinit = heurinit;
   (*heur)->heurexit = heurexit;
   (*heur)->heurexec = heurexec;
   (*heur)->heurdata = heurdata;
   CHECK_OKAY( SCIPclockCreate(&(*heur)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*heur)->ncalls = 0;
   (*heur)->nsolsfound = 0;
   (*heur)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "heuristics/%s/priority", name);
   sprintf(paramdesc, "priority of heuristic <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*heur)->priority, priority, INT_MIN, INT_MAX, 
                  paramChgdHeurPriority, (PARAMDATA*)(*heur)) ); /*lint !e740*/
   sprintf(paramname, "heuristics/%s/freq", name);
   sprintf(paramdesc, "frequency for calling primal heuristic <%s> (-1: never, 0: only in root node)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*heur)->freq, freq, -1, INT_MAX, NULL, NULL) );
   sprintf(paramname, "heuristics/%s/freqofs", name);
   sprintf(paramdesc, "frequency offset for calling primal heuristic <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*heur)->freqofs, freqofs, 0, INT_MAX, NULL, NULL) );
   sprintf(paramname, "heuristics/%s/maxdepth", name);
   sprintf(paramdesc, "maximal depth level to call primal heuristic <%s> (-1: no limit)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*heur)->maxdepth, maxdepth, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of primal heuristic */
RETCODE SCIPheurFree(
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(heur != NULL);
   assert(*heur != NULL);
   assert(!(*heur)->initialized);
   assert(scip != NULL);

   /* call destructor of primal heuristic */
   if( (*heur)->heurfree != NULL )
   {
      CHECK_OKAY( (*heur)->heurfree(scip, *heur) );
   }

   SCIPclockFree(&(*heur)->clock);
   freeMemoryArray(&(*heur)->name);
   freeMemoryArray(&(*heur)->desc);
   freeMemory(heur);

   return SCIP_OKAY;
}

/** initializes primal heuristic */
RETCODE SCIPheurInit(
   HEUR*            heur,               /**< primal heuristic */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(heur != NULL);
   assert(scip != NULL);

   if( heur->initialized )
   {
      errorMessage("primal heuristic <%s> already initialized\n", heur->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(heur->clock);

   heur->ncalls = 0;
   heur->nsolsfound = 0;

   if( heur->heurinit != NULL )
   {
      CHECK_OKAY( heur->heurinit(scip, heur) );
   }
   heur->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of primal heuristic */
RETCODE SCIPheurExit(
   HEUR*            heur,               /**< primal heuristic */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(heur != NULL);
   assert(scip != NULL);

   if( !heur->initialized )
   {
      errorMessage("primal heuristic <%s> not initialized\n", heur->name);
      return SCIP_INVALIDCALL;
   }

   if( heur->heurexit != NULL )
   {
      CHECK_OKAY( heur->heurexit(scip, heur) );
   }
   heur->initialized = FALSE;

   return SCIP_OKAY;
}

/** initializes solution process data of primal heuristic */
RETCODE SCIPheurInitsol(
   HEUR*            heur,               /**< primal heuristic */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(heur != NULL);
   assert(set != NULL);

   if( heur->delaypos != -1 )
   {
      heur->delaypos = -1;
      set->heurssorted = FALSE;
   }

   return SCIP_OKAY;
}

/** calls execution method of primal heuristic */
RETCODE SCIPheurExec(
   HEUR*            heur,               /**< primal heuristic */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   int              depth,              /**< depth of current node */
   int              lpforkdepth,        /**< depth of the last node with solved LP */
   Bool             currentnodehaslp,   /**< is LP being processed in the current node? */
   Bool             plunging,           /**< is the next node to be processed a child or sibling? */
   int*             ndelayedheurs,      /**< pointer to count the number of delayed heuristics */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   Bool execute;

   assert(heur != NULL);
   assert(heur->heurexec != NULL);
   assert(heur->freq >= -1);
   assert(heur->freqofs >= 0);
   assert(heur->maxdepth >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(primal != NULL);
   assert(depth >= 0);
   assert(ndelayedheurs != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( heur->pseudonodes )
   {
      /* heuristic may be executed on every node: check, if the current depth matches the execution frequency and offset */
      execute = (heur->freq > 0 && (depth - heur->freqofs) % heur->freq == 0);
   }
   else
   {
      /* heuristic may only be executed on LP nodes: check, if a node matching the execution frequency lies between the
       * current node and the last LP node of the path
       */
      execute = (heur->freq > 0
         && ((depth + heur->freq - heur->freqofs) / heur->freq
            != (lpforkdepth + heur->freq - heur->freqofs) / heur->freq));
   }

   /* if frequency is zero, execute heuristic at root node */
   execute = execute || (depth == 0 && heur->freq == 0);

   /* compare current depth against heuristic's maximal depth level */
   execute = execute && (heur->maxdepth == -1 || depth <= heur->maxdepth);
   
   /* if the heuristic was delayed, execute it anywas */
   execute = execute || (heur->delaypos >= 0);
   
   /* execute LP heuristics only at LP nodes */
   execute = execute && (heur->pseudonodes || currentnodehaslp);

   if( execute )
   {
      if( plunging && !heur->duringplunging && depth > 0 )
      {
         /* the heuristic should be delayed until plunging is finished */
         *result = SCIP_DELAYED;
      }
      else
      {
         Longint oldnsolsfound;
         
         debugMessage("executing primal heuristic <%s> in depth %d (delaypos: %d)\n", heur->name, depth, heur->delaypos);

         oldnsolsfound = primal->nsolsfound;

         /* start timing */
         SCIPclockStart(heur->clock, set);

         /* call external method */
         CHECK_OKAY( heur->heurexec(set->scip, heur, result) );

         /* stop timing */
         SCIPclockStop(heur->clock, set);

         /* evaluate result */
         if( *result != SCIP_FOUNDSOL
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            errorMessage("execution method of primal heuristic <%s> returned invalid result <%d>\n", 
               heur->name, *result);
            return SCIP_INVALIDRESULT;
         }
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            heur->ncalls++;
         heur->nsolsfound += primal->nsolsfound - oldnsolsfound;

         /* update delay position of heuristic */
         if( *result != SCIP_DELAYED && heur->delaypos != -1 )
         {
            heur->delaypos = -1;
            set->heurssorted = FALSE;
         }
      }
   }
   assert(*result == SCIP_DIDNOTRUN || *result == SCIP_DELAYED || heur->delaypos == -1);

   /* check if the heuristic was (still) delayed */
   if( *result == SCIP_DELAYED || heur->delaypos >= 0 )
   {
      debugMessage("delaying execution of primal heuristic <%s> in depth %d (delaypos: %d)\n", 
         heur->name, depth, *ndelayedheurs);

      /* mark the heuristic delayed */
      if( heur->delaypos != *ndelayedheurs )
      {
         heur->delaypos = *ndelayedheurs;
         set->heurssorted = FALSE;
      }
      (*ndelayedheurs)++;
   }

   return SCIP_OKAY;
}

/** gets user data of primal heuristic */
HEURDATA* SCIPheurGetData(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->heurdata;
}

/** sets user data of primal heuristic; user has to free old data in advance! */
void SCIPheurSetData(
   HEUR*            heur,               /**< primal heuristic */
   HEURDATA*        heurdata            /**< new primal heuristic user data */
   )
{
   assert(heur != NULL);

   heur->heurdata = heurdata;
}

/** gets name of primal heuristic */
const char* SCIPheurGetName(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->name;
}

/** gets description of primal heuristic */
const char* SCIPheurGetDesc(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->desc;
}

/** gets display character of primal heuristic */
char SCIPheurGetDispchar(
   HEUR*            heur                /**< primal heuristic */
   )
{
   if( heur == NULL )
      return '*';
   else
      return heur->dispchar;
}

/** gets priority of primal heuristic */
int SCIPheurGetPriority(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->priority;
}

/** sets priority of primal heuristic */
void SCIPheurSetPriority(
   HEUR*            heur,               /**< primal heuristic */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the primal heuristic */
   )
{
   assert(heur != NULL);
   assert(set != NULL);
   
   heur->priority = priority;
   set->heurssorted = FALSE;
}

/** gets frequency of primal heuristic */
int SCIPheurGetFreq(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->freq;
}

/** gets frequency offset of primal heuristic */
int SCIPheurGetFreqofs(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->freqofs;
}

/** gets maximal depth level for calling primal heuristic (returns -1, if no depth limit exists) */
int SCIPheurGetMaxdepth(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->maxdepth;
}

/** gets the number of times, the heuristic was called and tried to find a solution */
Longint SCIPheurGetNCalls(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->ncalls;
}

/** gets the number of primal feasible solutions found by this heuristic */
Longint SCIPheurGetNSolsFound(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->nsolsfound;
}

/** is primal heuristic initialized? */
Bool SCIPheurIsInitialized(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->initialized;
}

/** gets time in seconds used in this heuristic */
Real SCIPheurGetTime(
   HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return SCIPclockGetTime(heur->clock);
}

