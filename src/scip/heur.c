/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
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

/* resets diving settings counters */
void SCIPdivesetReset(
   SCIP_DIVESET*         diveset             /**< diveset to be reset */
   )
{
   assert(diveset != NULL);

   diveset->nlpiterations = 0L;
   diveset->totaldepth = 0L;
   diveset->totalsoldepth = 0L;
   diveset->totalnnodes = 0L;
   diveset->totalnbacktracks = 0L;
   diveset->minsoldepth = INT_MAX;
   diveset->maxsoldepth = -1;
   diveset->mindepth = INT_MAX;
   diveset->maxdepth = -1;
   diveset->nlps = 0;
   diveset->nsolsfound = 0;
   diveset->nbestsolsfound = 0;
   diveset->ncalls = 0;
   diveset->nsolcalls = 0;
}

/** update diveset statistics and global diveset statistics */
void SCIPdivesetUpdateStats(
   SCIP_DIVESET*         diveset,            /**< diveset to be reset */
   SCIP_STAT*            stat,               /**< global SCIP statistics */
   int                   depth,              /**< the depth reached this time */
   int                   nprobingnodes,      /**< the number of probing nodes explored this time */
   int                   nbacktracks,        /**< the number of backtracks during probing this time */
   SCIP_Longint          nsolsfound,         /**< number of new solutions found this time */
   SCIP_Longint          nbestsolsfound,     /**< number of new best solutions found this time */
   SCIP_Bool             leavesol            /**< has the diving heuristic reached a feasible leaf */
   )
{
   assert(diveset != NULL);

   diveset->totaldepth += depth;
   diveset->mindepth = MIN(diveset->mindepth, depth);
   diveset->maxdepth = MAX(diveset->maxdepth, depth);
   diveset->totalnnodes += nprobingnodes;
   diveset->totalnbacktracks += nbacktracks;
   diveset->ncalls++;

   /* update solution statistics only if a solution was found */
   if( leavesol )
   {
      diveset->totalsoldepth += depth;
      diveset->minsoldepth = MIN(diveset->minsoldepth, depth);
      diveset->maxsoldepth = MAX(diveset->maxsoldepth, depth);
      diveset->nsolcalls++;
   }

   diveset->nsolsfound += nsolsfound;
   diveset->nbestsolsfound += nbestsolsfound;

   stat->totaldivesetdepth += depth;
   stat->ndivesetcalls++;
}

/** append diveset to heuristic array of divesets */
static
SCIP_RETCODE heurAddDiveset(
   SCIP_HEUR*            heur,               /**< the heuristic to which this dive setting belongs */
   SCIP_DIVESET*         diveset             /**< pointer to the freshly created diveset */
   )
{
   assert(heur != NULL);
   assert(diveset != NULL);
   assert(diveset->heur == NULL);

   diveset->heur = heur;

   if( heur->divesets == NULL )
   {
      assert(heur->ndivesets == 0);
      SCIP_ALLOC( BMSallocMemoryArray(&heur->divesets, 1) );
   }
   else
   {
      assert(heur->ndivesets > 0);
      SCIP_ALLOC( BMSreallocMemoryArray(&heur->divesets, heur->ndivesets + 1) ); /*lint !e776 I expect no overflow here */
   }

   /* append diveset to the end of the array */
   heur->divesets[heur->ndivesets] = diveset;
   heur->ndivesets++;

   return SCIP_OKAY;
}

/** create a set of diving heuristic settings */
SCIP_RETCODE SCIPdivesetCreate(
   SCIP_DIVESET**        diveset,            /**< pointer to the freshly created diveset */
   SCIP_HEUR*            heur,               /**< the heuristic to which this dive setting belongs */
   const char*           name,               /**< name for the diveset, or NULL if the name of the heuristic should be used */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_Real             minreldepth,        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth,        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot,      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   SCIP_Real             maxdiveubquot,      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot,     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol, /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol,/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             lpresolvedomchgquot,/**< percentage of immediate domain changes during probing to trigger LP resolve */
   int                   lpsolvefreq,        /**< LP solve frequency for (0: only if enough domain reductions are found by propagation)*/
   int                   maxlpiterofs,       /**< additional number of allowed LP iterations */
   SCIP_Bool             backtrack,          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Bool             onlylpbranchcands,  /**< should only LP branching candidates be considered instead of the slower but
                                              *   more general constraint handler diving variable selection? */
   SCIP_DIVETYPE         divetypemask,       /**< bit mask that represents the supported dive types by this dive set */
   SCIP_DECL_DIVESETGETSCORE((*divesetgetscore))  /**< method for candidate score and rounding direction */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   const char* divesetname;

   assert(diveset != NULL);
   assert(set != NULL);
   assert(divesetgetscore != NULL);
   assert(heur != NULL);

   SCIP_ALLOC( BMSallocMemory(diveset) );

   /* for convenience, the name gets inferred from the heuristic to which the diveset is added if no name is provided */
   divesetname = (name == NULL ? SCIPheurGetName(heur) : name);
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*diveset)->name, divesetname, strlen(divesetname)+1) );
   (*diveset)->heur = NULL;

   /* copy callbacks */
   (*diveset)->divesetgetscore = divesetgetscore;

   SCIP_CALL( heurAddDiveset(heur, *diveset) );
   (*diveset)->sol = NULL;
   (*diveset)->divetypemask = divetypemask;

   /* add collection of diving heuristic specific parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/minreldepth", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem,
         paramname, "minimal relative depth to start diving",
         &(*diveset)->minreldepth, TRUE, minreldepth, 0.0, 1.0, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxreldepth", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "maximal relative depth to start diving",
         &(*diveset)->maxreldepth, TRUE, maxreldepth, 0.0, 1.0, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxlpiterquot", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem,
         paramname,
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &(*diveset)->maxlpiterquot, FALSE, maxlpiterquot, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxlpiterofs", (*diveset)->name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem,
         paramname,
         "additional number of allowed LP iterations",
         &(*diveset)->maxlpiterofs, FALSE, maxlpiterofs, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxdiveubquot", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem,
         paramname,
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &(*diveset)->maxdiveubquot, TRUE, maxdiveubquot, 0.0, 1.0, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxdiveavgquot", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem,
         paramname,
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &(*diveset)->maxdiveavgquot, TRUE, maxdiveavgquot, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxdiveubquotnosol", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem,
         paramname,
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &(*diveset)->maxdiveubquotnosol, TRUE, maxdiveubquotnosol, 0.0, 1.0, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/maxdiveavgquotnosol", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem,
         paramname,
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &(*diveset)->maxdiveavgquotnosol, TRUE, maxdiveavgquotnosol, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/backtrack", (*diveset)->name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem,
         paramname,
         "use one level of backtracking if infeasibility is encountered?",
         &(*diveset)->backtrack, FALSE, backtrack, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/lpresolvedomchgquot", (*diveset)->name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "percentage of immediate domain changes during probing to trigger LP resolve",
         &(*diveset)->lpresolvedomchgquot, FALSE, lpresolvedomchgquot,  0.0, SCIP_REAL_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/lpsolvefreq", (*diveset)->name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem,
         paramname,
         "LP solve frequency for diving heuristics (0: only after enough domain changes have been found)",
         &(*diveset)->lpsolvefreq, FALSE, lpsolvefreq, 0, INT_MAX,
         NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "heuristics/%s/onlylpbranchcands", (*diveset)->name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem,
            paramname,
            "should only LP branching candidates be considered instead of the slower but "
            "more general constraint handler diving variable selection?",
            &(*diveset)->onlylpbranchcands, FALSE, onlylpbranchcands, NULL, NULL) );

   SCIPdivesetReset(*diveset);

   return SCIP_OKAY;
}

/** get the heuristic to which this diving setting belongs */
SCIP_HEUR* SCIPdivesetGetHeur(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->heur;
}

/** get the working solution of this dive set */
SCIP_SOL* SCIPdivesetGetWorkSolution(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->sol;
}

/** set the working solution for this dive set */
void SCIPdivesetSetWorkSolution(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_SOL*             sol                 /**< new working solution for this dive set, or NULL */
   )
{
   assert(diveset != NULL);

   diveset->sol = sol;
}

/** get the name of the dive set */
const char* SCIPdivesetGetName(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->name;
}

/** get the minimum relative depth of the diving settings */
SCIP_Real SCIPdivesetGetMinRelDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->minreldepth;
}

/** get the maximum relative depth of the diving settings */
SCIP_Real SCIPdivesetGetMaxRelDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxreldepth;
}

/** get the number of successful runs of the diving settings */
SCIP_Longint SCIPdivesetGetSolSuccess(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return 10 * diveset->nbestsolsfound + diveset->nsolsfound;
}

/** get the number of calls to this dive set */
int SCIPdivesetGetNCalls(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->ncalls;
}

/** get the number of calls successfully terminated at a feasible leaf node */
int SCIPdivesetGetNSolutionCalls(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->nsolcalls;
}

/** get the minimum depth reached by this dive set */
int SCIPdivesetGetMinDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->mindepth;
}

/** get the maximum depth reached by this dive set */
int SCIPdivesetGetMaxDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->maxdepth;
}

/** get the average depth this dive set reached during execution */
SCIP_Real SCIPdivesetGetAvgDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return (diveset->ncalls == 0 ? 0.0 : diveset->totaldepth / (SCIP_Real)diveset->ncalls);
}

/** get the minimum depth at which this dive set found a solution */
int SCIPdivesetGetMinSolutionDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->minsoldepth;
}

/** get the maximum depth at which this dive set found a solution */
int SCIPdivesetGetMaxSolutionDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->maxsoldepth;
}

/** get the average depth at which this dive set found a solution */
SCIP_Real SCIPdivesetGetAvgSolutionDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return (diveset->nsolcalls == 0 ? 0.0 : diveset->totalsoldepth / (SCIP_Real)diveset->nsolcalls);
}

/** get the total number of LP iterations used by this dive set */
SCIP_Longint SCIPdivesetGetNLPIterations(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->nlpiterations;
}

/** get the total number of probing nodes used by this dive set */
SCIP_Longint SCIPdivesetGetNProbingNodes(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->totalnnodes;
}

/** get the total number of backtracks performed by this dive set */
SCIP_Longint SCIPdivesetGetNBacktracks(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->totalnbacktracks;
}

/** get the maximum LP iterations quotient of the diving settings */
SCIP_Real SCIPdivesetGetMaxLPIterQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxlpiterquot;
}

/** get the maximum LP iterations offset of the diving settings */
int SCIPdivesetGetMaxLPIterOffset(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxlpiterofs;
}

/** get the maximum upper bound quotient parameter of the diving settings if no solution is available */
SCIP_Real SCIPdivesetGetUbQuotNoSol(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxdiveubquotnosol;
}

/** get the average quotient parameter of the diving settings if no solution is available */
SCIP_Real SCIPdivesetGetAvgQuotNoSol(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxdiveavgquotnosol;
}
/** get the maximum upper bound quotient parameter of the diving settings if an incumbent solution exists */
SCIP_Real SCIPdivesetGetUbQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxdiveubquot;
}

/** get the average upper bound quotient parameter of the diving settings if an incumbent solution exists */
SCIP_Real SCIPdivesetGetAvgQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->maxdiveavgquot;
}

/** should backtracking be applied? */
SCIP_Bool SCIPdivesetUseBacktrack(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   return diveset->backtrack;
}

/** returns the LP solve frequency for diving LPs (0: dynamically based on number of intermediate domain reductions) */
int SCIPdivesetGetLPSolveFreq(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->lpsolvefreq;
}

/** returns the domain reduction quotient for triggering an immediate resolve of the diving LP (0.0: always resolve)*/
SCIP_Real SCIPdivesetGetLPResolveDomChgQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->lpresolvedomchgquot;
}

/** should only LP branching candidates be considered instead of the slower but
 *  more general constraint handler diving variable selection?
 */
SCIP_Bool SCIPdivesetUseOnlyLPBranchcands(
   SCIP_DIVESET*         diveset             /**< diving settings */
   )
{
   assert(diveset != NULL);

   return diveset->onlylpbranchcands;
}

/** returns TRUE if dive set supports diving of the specified type */
SCIP_Bool SCIPdivesetSupportsType(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVETYPE         divetype            /**< bit mask that represents the supported dive types by this dive set */
   )
{
   assert(diveset != NULL);

   return (divetype & diveset->divetypemask);
}

/** update diveset LP statistics, should be called after every LP solved by this diving heuristic */
void SCIPdivesetUpdateLPStats(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_STAT*            stat,               /**< global SCIP statistics */
   SCIP_Longint          niterstoadd         /**< additional number of LP iterations to be added */
   )
{
   diveset->nlpiterations += niterstoadd;
   stat->ndivesetlpiterations += niterstoadd;
   diveset->nlps++;
   stat->ndivesetlps++;
}

/** frees memory of a diveset */
static
void divesetFree(
   SCIP_DIVESET**        diveset             /**< general diving settings */
   )
{
   assert(*diveset != NULL);
   assert((*diveset)->name != NULL);

   BMSfreeMemoryArray(&(*diveset)->name);
   BMSfreeMemory(diveset);
}

/** get the candidate score and preferred rounding direction for a candidate variable */
SCIP_RETCODE SCIPdivesetGetScore(
   SCIP_DIVESET*         diveset,            /**< general diving settings */
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_DIVETYPE         divetype,           /**< the type of diving that should be applied */
   SCIP_VAR*             divecand,           /**< the candidate for which the branching direction is requested */
   SCIP_Real             divecandsol,        /**< LP solution value of the candidate */
   SCIP_Real             divecandfrac,       /**< fractionality of the candidate */
   SCIP_Real*            candscore,          /**< pointer to store the candidate score */
   SCIP_Bool*            roundup             /**< pointer to store whether preferred direction for diving is upwards */
   )
{
   assert(diveset->divesetgetscore != NULL);
   assert(candscore != NULL);
   assert(roundup != NULL);
   assert(divecand != NULL);
   assert(divetype & diveset->divetypemask);

   SCIP_CALL( diveset->divesetgetscore(set->scip, diveset, divetype, divecand, divecandsol, divecandfrac, candscore, roundup) );

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
   (*heur)->divesets = NULL;
   (*heur)->ndivesets = 0;

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
   int d;
   assert(heur != NULL);
   assert(*heur != NULL);
   assert(!(*heur)->initialized);
   assert(set != NULL);
   assert((*heur)->divesets != NULL || (*heur)->ndivesets == 0);

   /* call destructor of primal heuristic */
   if( (*heur)->heurfree != NULL )
   {
      SCIP_CALL( (*heur)->heurfree(set->scip, *heur) );
   }

   for( d = 0; d < (*heur)->ndivesets; ++d )
   {
      assert((*heur)->divesets[d] != NULL);
      divesetFree(&((*heur)->divesets[d]));
   }
   BMSfreeMemoryArrayNull(&(*heur)->divesets);
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
   int d;
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

   /* reset dive sets */
   for( d = 0; d < heur->ndivesets; ++d )
   {
      assert(heur->divesets[d] != NULL);
      SCIPdivesetReset(heur->divesets[d]);
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
   SCIP_Bool             nodeinfeasible,     /**< was the current node already detected to be infeasible? */
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
      SCIP_CALL( heur->heurexec(set->scip, heur, heurtiming, nodeinfeasible, result) );

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
SCIP_Bool SCIPheurUsesSubscip(
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

/** enables or disables all clocks of \p heur, depending on the value of the flag */
void SCIPheurEnableOrDisableClocks(
   SCIP_HEUR*            heur,               /**< the heuristic for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the heuristic be enabled? */
   )
{
   assert(heur != NULL);

   SCIPclockEnableOrDisable(heur->setuptime, enable);
   SCIPclockEnableOrDisable(heur->heurclock, enable);
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

/** returns array of divesets of this primal heuristic, or NULL if it has no divesets */
SCIP_DIVESET** SCIPheurGetDivesets(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->divesets;
}

/** returns the number of divesets of this primal heuristic */
int SCIPheurGetNDivesets(
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(heur != NULL);

   return heur->ndivesets;
}
