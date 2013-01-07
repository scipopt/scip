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

/**@file   heur.c
 * @brief  methods for primal heuristics
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/primal.h"
#include "scip/scip.h"
#include "scip/heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_heur.h"



/** compares two heuristics w. r. to their delay positions and their priority */
SCIP_DECL_SORTPTRCOMP(SCIPheurComp)
{  /*lint --e{715}*/
   SCIP_HEUR* heur1 = (SCIP_HEUR*)elem1;
   SCIP_HEUR* heur2 = (SCIP_HEUR*)elem2;

   assert(heur1 != NULL);
   assert(heur2 != NULL);

   if( heur1->delaypos == heur2->delaypos )
      return heur2->priority - heur1->priority; /* prefer higher priorities */
   else if( heur1->delaypos == -1 )
      return +1;                                /* prefer delayed heuristics */
   else if( heur2->delaypos == -1 )
      return -1;                                /* prefer delayed heuristics */
   else if( heur1->ncalls * heur1->freq > heur2->ncalls * heur2->freq )
      return +1;
   else if( heur1->ncalls * heur1->freq < heur2->ncalls * heur2->freq )
      return -1;
   else
      return heur1->delaypos - heur2->delaypos; /* prefer lower delay positions */
}


/** comparison method for sorting heuristics w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPheurCompName)
{
   return strcmp(SCIPheurGetName((SCIP_HEUR*)elem1), SCIPheurGetName((SCIP_HEUR*)elem2));
}

/** method to call, when the priority of a heuristic was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdHeurPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetHeurPriority() to mark the heuristics unsorted */
   SCIP_CALL( SCIPsetHeurPriority(scip, (SCIP_HEUR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given primal heuristic to a new scip */
SCIP_RETCODE SCIPheurCopyInclude(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(heur != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( heur->heurcopy != NULL )
   {
      SCIPdebugMessage("including heur %s in subscip %p\n", SCIPheurGetName(heur), (void*)set->scip);
      SCIP_CALL( heur->heurcopy(set->scip, heur) );
   }

   return SCIP_OKAY;
}

/** creates a primal heuristic */
SCIP_RETCODE SCIPheurCreate(
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   unsigned int          timingmask,         /**< positions in the node solving loop where heuristic should be executed */
   SCIP_Bool             usessubscip,        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_DECL_HEURCOPY    ((*heurcopy)),      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol)),   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol)),   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(heur != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(freqofs >= 0);
   assert(heurexec != NULL);

   SCIP_ALLOC( BMSallocMemory(heur) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*heur)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*heur)->desc, desc, strlen(desc)+1) );
   (*heur)->dispchar = dispchar;
   (*heur)->priority = priority;
   (*heur)->freq = freq;
   (*heur)->freqofs = freqofs;
   (*heur)->maxdepth = maxdepth;
   (*heur)->delaypos = -1;
   (*heur)->timingmask = timingmask;
   (*heur)->usessubscip = usessubscip;
   (*heur)->heurcopy = heurcopy;
   (*heur)->heurfree = heurfree;
   (*heur)->heurinit = heurinit;
   (*heur)->heurexit = heurexit;
   (*heur)->heurinitsol = heurinitsol;
   (*heur)->heurexitsol = heurexitsol;
   (*heur)->heurexec = heurexec;
   (*heur)->heurdata = heurdata;
   SCIP_CALL( SCIPclockCreate(&(*heur)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*heur)->heurclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*heur)->ncalls = 0;
   (*heur)->nsolsfound = 0;
   (*heur)->nbestsolsfound = 0;
   (*heur)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of heuristic <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*heur)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdHeurPriority, (SCIP_PARAMDATA*)(*heur)) ); /*lint !e740*/
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/freq", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "frequency for calling primal heuristic <%s> (-1: never, 0: only at depth freqofs)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*heur)->freq, FALSE, freq, -1, INT_MAX, NULL, NULL) );
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/freqofs", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "frequency offset for calling primal heuristic <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*heur)->freqofs, FALSE, freqofs, 0, INT_MAX, NULL, NULL) );
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxdepth", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "maximal depth level to call primal heuristic <%s> (-1: no limit)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*heur)->maxdepth, TRUE, maxdepth, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of primal heuristic */
SCIP_RETCODE SCIPheurFree(
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(heur != NULL);
   assert(*heur != NULL);
   assert(!(*heur)->initialized);
   assert(set != NULL);

   /* call destructor of primal heuristic */
   if( (*heur)->heurfree != NULL )
   {
      SCIP_CALL( (*heur)->heurfree(set->scip, *heur) );
   }

   SCIPclockFree(&(*heur)->heurclock);
   SCIPclockFree(&(*heur)->setuptime);
   BMSfreeMemoryArray(&(*heur)->name);
   BMSfreeMemoryArray(&(*heur)->desc);
   BMSfreeMemory(heur);

   return SCIP_OKAY;
}

/** initializes primal heuristic */
SCIP_RETCODE SCIPheurInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(heur != NULL);
   assert(set != NULL);

   if( heur->initialized )
   {
      SCIPerrorMessage("primal heuristic <%s> already initialized\n", heur->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(heur->setuptime);
      SCIPclockReset(heur->heurclock);
      
      heur->delaypos = -1;
      heur->ncalls = 0;
      heur->nsolsfound = 0;
      heur->nbestsolsfound = 0;
   }

   if( heur->heurinit != NULL )
   {
      /* start timing */
      SCIPclockStart(heur->setuptime, set);

      SCIP_CALL( heur->heurinit(set->scip, heur) );

      /* stop timing */
      SCIPclockStop(heur->setuptime, set);
   }
   heur->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of primal heuristic */
SCIP_RETCODE SCIPheurExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(heur != NULL);
   assert(set != NULL);

   if( !heur->initialized )
   {
      SCIPerrorMessage("primal heuristic <%s> not initialized\n", heur->name);
      return SCIP_INVALIDCALL;
   }

   if( heur->heurexit != NULL )
   {
      /* start timing */
      SCIPclockStart(heur->setuptime, set);

      SCIP_CALL( heur->heurexit(set->scip, heur) );

      /* stop timing */
      SCIPclockStop(heur->setuptime, set);
   }
   heur->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs primal heuristic that the branch and bound process is being started */
SCIP_RETCODE SCIPheurInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(heur != NULL);
   assert(set != NULL);

   if( heur->delaypos != -1 )
   {
      heur->delaypos = -1;
      set->heurssorted = FALSE;
   }

   /* call solving process initialization method of primal heuristic */
   if( heur->heurinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(heur->setuptime, set);

      SCIP_CALL( heur->heurinitsol(set->scip, heur) );

      /* stop timing */
      SCIPclockStop(heur->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs primal heuristic that the branch and bound process data is being freed */
SCIP_RETCODE SCIPheurExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(heur != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of primal heuristic */
   if( heur->heurexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(heur->setuptime, set);

      SCIP_CALL( heur->heurexitsol(set->scip, heur) );

      /* stop timing */
      SCIPclockStop(heur->setuptime, set);
   }

   return SCIP_OKAY;
}

/** should the heuristic be executed at the given depth, frequency, timing, ... */
SCIP_Bool SCIPheurShouldBeExecuted(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   depth,              /**< depth of current node */
   int                   lpstateforkdepth,   /**< depth of the last node with solved LP */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   SCIP_Bool*            delayed             /**< pointer to store whether the heuristic should be delayed */
   )
{
   SCIP_Bool execute;

   if( ((heur->timingmask & SCIP_HEURTIMING_BEFOREPRESOL) && heurtiming == SCIP_HEURTIMING_BEFOREPRESOL)
       || ((heur->timingmask & SCIP_HEURTIMING_DURINGPRESOLLOOP) && heurtiming == SCIP_HEURTIMING_DURINGPRESOLLOOP) )
   {
      /* heuristic may be executed before/during presolving. Do so, if it was not disabled by setting the frequency to -1 */
      execute = heur->freq >= 0; 
   } 
   else if( (heur->timingmask & SCIP_HEURTIMING_AFTERPSEUDONODE) == 0
      && (heurtiming == SCIP_HEURTIMING_AFTERLPNODE || heurtiming == SCIP_HEURTIMING_AFTERLPPLUNGE) )
   {
      /* heuristic was skipped on intermediate pseudo nodes: check, if a node matching the execution frequency lies
       * between the current node and the last LP node of the path
       */
      execute = (heur->freq > 0 && depth >= heur->freqofs 
         && ((depth + heur->freq - heur->freqofs) / heur->freq
            != (lpstateforkdepth + heur->freq - heur->freqofs) / heur->freq));
   }
   else
   {
      /* heuristic may be executed on every node: check, if the current depth matches the execution frequency and offset */
      execute = (heur->freq > 0 && depth >= heur->freqofs && (depth - heur->freqofs) % heur->freq == 0);
   }

   /* if frequency is zero, execute heuristic only at the depth level of the frequency offset */
   execute = execute || (depth == heur->freqofs && heur->freq == 0);

   /* compare current depth against heuristic's maximal depth level */
   execute = execute && (heur->maxdepth == -1 || depth <= heur->maxdepth);

   /* if the heuristic was delayed, execute it anyway */
   execute = execute || (heur->delaypos >= 0);

   /* if the heuristic should be called after plunging but not during plunging, delay it if we are in plunging */
   if( execute
      && ((heurtiming == SCIP_HEURTIMING_AFTERLPNODE
            && (heur->timingmask & SCIP_HEURTIMING_AFTERLPNODE) == 0
            && (heur->timingmask & SCIP_HEURTIMING_AFTERLPPLUNGE) > 0)
         || (heurtiming == SCIP_HEURTIMING_AFTERPSEUDONODE
            && (heur->timingmask & SCIP_HEURTIMING_AFTERPSEUDONODE) == 0
            && (heur->timingmask & SCIP_HEURTIMING_AFTERPSEUDOPLUNGE) > 0)) )
   {
      /* the heuristic should be delayed until plunging is finished */
      execute = FALSE;
      *delayed = TRUE;
   }

   /* execute heuristic only if its timing mask fits the current point in the node solving process */
   execute = execute && (heur->timingmask & heurtiming) > 0;

   return execute;
}

/** calls execution method of primal heuristic */
SCIP_RETCODE SCIPheurExec(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   int                   depth,              /**< depth of current node */
   int                   lpstateforkdepth,   /**< depth of the last node with solved LP */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   int*                  ndelayedheurs,      /**< pointer to count the number of delayed heuristics */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_Bool execute;
   SCIP_Bool delayed;

   assert(heur != NULL);
   assert(heur->heurexec != NULL);
   assert(heur->freq >= -1);
   assert(heur->freqofs >= 0);
   assert(heur->maxdepth >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(primal != NULL);
   assert(depth >= 0 || heurtiming == SCIP_HEURTIMING_BEFOREPRESOL || heurtiming == SCIP_HEURTIMING_DURINGPRESOLLOOP);
   assert(ndelayedheurs != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   delayed = FALSE;
   execute = SCIPheurShouldBeExecuted(heur, depth, lpstateforkdepth, heurtiming, &delayed);

   if( delayed )
   {
      assert(!execute);
      *result = SCIP_DELAYED;
   }

   if( execute )
   {
      SCIP_Longint oldnsolsfound;
      SCIP_Longint oldnbestsolsfound;
         
      SCIPdebugMessage("executing primal heuristic <%s> in depth %d (delaypos: %d)\n", heur->name, depth, heur->delaypos);

      oldnsolsfound = primal->nsolsfound;
      oldnbestsolsfound = primal->nbestsolsfound;

      /* start timing */
      SCIPclockStart(heur->heurclock, set);

      /* call external method */
      SCIP_CALL( heur->heurexec(set->scip, heur, heurtiming, result) );

      /* stop timing */
      SCIPclockStop(heur->heurclock, set);

      /* evaluate result */
      if( *result != SCIP_FOUNDSOL
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN
         && *result != SCIP_DELAYED )
      {
         SCIPerrorMessage("execution method of primal heuristic <%s> returned invalid result <%d>\n", 
            heur->name, *result);
         return SCIP_INVALIDRESULT;
      }
      if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
         heur->ncalls++;
      heur->nsolsfound += primal->nsolsfound - oldnsolsfound;
      heur->nbestsolsfound += primal->nbestsolsfound - oldnbestsolsfound;

      /* update delay position of heuristic */
      if( *result != SCIP_DELAYED && heur->delaypos != -1 )
      {
         heur->delaypos = -1;
         set->heurssorted = FALSE;
      }
   }
   assert(*result == SCIP_DIDNOTRUN || *result == SCIP_DELAYED || heur->delaypos == -1);

   /* check if the heuristic was (still) delayed */
   if( *result == SCIP_DELAYED || heur->delaypos >= 0 )
   {
      SCIPdebugMessage("delaying execution of primal heuristic <%s> in depth %d (delaypos: %d), heur was%s delayed before, had delaypos %d\n",
         heur->name, depth, *ndelayedheurs, heur->delaypos >= 0 ? "" : " not", heur->delaypos);

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
SCIP_HEURDATA* SCIPheurGetData(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->heurdata;
}

/** sets user data of primal heuristic; user has to free old data in advance! */
void SCIPheurSetData(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< new primal heuristic user data */
   )
{
   assert(heur != NULL);

   heur->heurdata = heurdata;
}

/* new callback setter methods */

/** sets copy callback of primal heuristic */
void SCIPheurSetCopy(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURCOPY    ((*heurcopy))       /**< copy callback of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(heur != NULL);

   heur->heurcopy = heurcopy;
}

/** sets destructor callback of primal heuristic */
void SCIPheurSetFree(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURFREE    ((*heurfree))       /**< destructor of primal heuristic */
   )
{
   assert(heur != NULL);

   heur->heurfree = heurfree;
}

/** sets initialization callback of primal heuristic */
void SCIPheurSetInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit))       /**< initialize primal heuristic */
   )
{
   assert(heur != NULL);

   heur->heurinit = heurinit;
}

/** sets deinitialization callback of primal heuristic */
void SCIPheurSetExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit))       /**< deinitialize primal heuristic */
   )
{
   assert(heur != NULL);

   heur->heurexit = heurexit;
}

/** sets solving process initialization callback of primal heuristic */
void SCIPheurSetInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol))    /**< solving process initialization callback of primal heuristic */
   )
{
   assert(heur != NULL);

   heur->heurinitsol = heurinitsol;
}

/** sets solving process deinitialization callback of primal heuristic */
void SCIPheurSetExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol))    /**< solving process deinitialization callback of primal heuristic */
   )
{
   assert(heur != NULL);

   heur->heurexitsol = heurexitsol;
}

/** gets name of primal heuristic */
const char* SCIPheurGetName(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->name;
}

/** gets description of primal heuristic */
const char* SCIPheurGetDesc(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->desc;
}

/** gets display character of primal heuristic */
char SCIPheurGetDispchar(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->dispchar;
}

/** returns the timing mask of the heuristic */
SCIP_HEURTIMING SCIPheurGetTimingmask(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->timingmask;
}

/** sets new timing mask for heuristic */
void SCIPheurSetTimingmask(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_HEURTIMING       timingmask          /**< new timing mask of heuristic */
   )
{
   assert(heur != NULL);

   heur->timingmask = timingmask;
}

/** does the heuristic use a secondary SCIP instance? */
SCIP_HEURTIMING SCIPheurUsesSubscip(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->usessubscip;
}

/** gets priority of primal heuristic */
int SCIPheurGetPriority(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->priority;
}

/** sets priority of primal heuristic */
void SCIPheurSetPriority(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the primal heuristic */
   )
{
   assert(heur != NULL);
   assert(set != NULL);
   
   heur->priority = priority;
   set->heurssorted = FALSE;
}

/** gets frequency of primal heuristic */
int SCIPheurGetFreq(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->freq;
}

/** sets frequency of primal heuristic */
void SCIPheurSetFreq(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   freq                /**< new frequency of heuristic */
   )
{
   assert(heur != NULL);
   
   heur->freq = freq;
}

/** gets frequency offset of primal heuristic */
int SCIPheurGetFreqofs(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->freqofs;
}

/** gets maximal depth level for calling primal heuristic (returns -1, if no depth limit exists) */
int SCIPheurGetMaxdepth(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->maxdepth;
}

/** gets the number of times, the heuristic was called and tried to find a solution */
SCIP_Longint SCIPheurGetNCalls(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->ncalls;
}

/** gets the number of primal feasible solutions found by this heuristic */
SCIP_Longint SCIPheurGetNSolsFound(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->nsolsfound;
}

/** gets the number of new best primal feasible solutions found by this heuristic */
SCIP_Longint SCIPheurGetNBestSolsFound(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->nbestsolsfound;
}

/** is primal heuristic initialized? */
SCIP_Bool SCIPheurIsInitialized(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->initialized;
}

/** gets time in seconds used in this heuristic for setting up for next stages */
SCIP_Real SCIPheurGetSetupTime(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return SCIPclockGetTime(heur->setuptime);
}

/** gets time in seconds used in this heuristic */
SCIP_Real SCIPheurGetTime(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return SCIPclockGetTime(heur->heurclock);
}
