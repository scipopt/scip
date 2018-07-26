/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_rpa.c
 * @brief  constraint handler for recursive circle packing
 * @author Benjamin Mueller
 *
 * This constraint handler is used to store information about which (not verified) rectangular patterns have been locked
 * and which circular patterns have not been tried to be verified yet.
 *
 * @todo Is it enough the lock the unverified circular pattern variables only in the positive direction?
 * @todo remove all unnecessary callbacks and defines
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_rpa.h"
#include "probdata_rpa.h"
#include "pattern.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "rpa"
#define CONSHDLR_DESC          "ringpacking constraint handler"
#define CONSHDLR_ENFOPRIORITY         1 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY       -1 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */

/* new best solution event handler properties */
#define EVENTHDLR_NAME         "bestsol"
#define EVENTHDLR_DESC         "best solution event handler"

/* default values of parameters */
#define DEFAULT_VERIFICATION_NLPTILIM        10.0      /**< time limit for each verification NLP */
#define DEFAULT_VERIFICATION_NLPNODELIM      10000L    /**< node limit for each verification NLP */
#define DEFAULT_VERIFICATION_HEURTILIM       10.0      /**< time limit for heuristic verification */
#define DEFAULT_VERIFICATION_HEURITERLIM     1000      /**< iteration limit for each heuristic verification */
#define DEFAULT_VERIFICATION_TOTALTILIM      3600.0    /**< total time limit for all verification problems during the solving process */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;         /**< event handler */

   SCIP_Bool*            locked;            /**< array to remember which (not verified) patterns have been locked */
   int                   lockedsize;        /**< size of locked array */
   SCIP_Bool*            tried;             /**< array to mark circular patterns that have been tried to be verified */\

   /* parameter for verification */
   SCIP_Real             timeleft;          /**< time left for solving verification problem during the solving process */
   SCIP_Real             nlptilim;          /**< hard time limit for verification nlp */
   SCIP_Real             heurtilim;         /**< hard time limit for verification heuristic*/
   SCIP_Longint          nlpnodelim;        /**< hard node limit for verification nlp */
   int                   heuriterlim;       /**< hard iteration limit for verification heuristic*/
};

/*
 * Local methods
 */

/** auxiliary function to decide whether a proposed solution is feasible; a solution is called feasible if and only if
 *  z*_C = 0 holds for all circular patterns that are either not packable, i.e., SCIP_PACKABLE_NO or SCIP_PACKABLE_UNKNOWN
 */
static
SCIP_Bool isSolFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution (NULL for LP solution) */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_PATTERN** cpatterns;
   SCIP_VAR** cvars;
   int ncpatterns;
   int p;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get information about circular patterns and their corresponding variables */
   SCIPprobdataGetCInfos(probdata, &cpatterns, &cvars, &ncpatterns);
   assert(ncpatterns > 0);

   for( p = 0; p < ncpatterns; ++p )
   {
      assert(cpatterns[p] != NULL);
      assert(cvars[p] != NULL);

      /* check only circular patterns which might not be packable */
      if( SCIPpatternGetPackableStatus(cpatterns[p]) != SCIP_PACKABLE_YES )
      {
         SCIP_Real solval = SCIPgetSolVal(scip, sol, cvars[p]);

         if( !SCIPisFeasZero(scip, solval) )
         {
            SCIPdebugMsg(scip, "solution might be infeasible because of circular pattern %d = (%g,%d)\n", p,
               SCIPgetSolVal(scip, sol, cvars[p]), SCIPpatternGetPackableStatus(cpatterns[p]));
            return FALSE;
         }
      }
   }

   return TRUE;
}

/** tries to verify a circular pattern; it first tries to call heuristic(s) and afterwards uses a verification NLP */
static
SCIP_RETCODE verifyCircularPattern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern             /**< circular pattern */
   )
{
   SCIP_Real timelimit;
   SCIP_Real nlptimelimit;
   SCIP_Real heurtimelimit;

   assert(probdata != NULL);
   assert(pattern != NULL);
   assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR);
   assert(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* get the total time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* verify heuristically */
   heurtimelimit = MIN(timelimit - SCIPgetSolvingTime(scip), conshdlrdata->heurtilim); /*lint !e666*/

   SCIPdebugMsg(scip, "call verification heuristic (%g,%d)\n", heurtimelimit, conshdlrdata->heuriterlim);
   conshdlrdata->timeleft += SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, heurtimelimit, conshdlrdata->heuriterlim) );
   conshdlrdata->timeleft -= SCIPgetSolvingTime(scip);

   /* use verification NLP if pattern is still not verified */
   if( SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN )
   {
      nlptimelimit = MIN3(conshdlrdata->timeleft, timelimit - SCIPgetSolvingTime(scip), conshdlrdata->nlptilim); /*lint !e666*/

      SCIPdebugMsg(scip, "call verification NLP (%g,%lld)\n", nlptimelimit, conshdlrdata->nlpnodelim);
      conshdlrdata->timeleft += SCIPgetSolvingTime(scip);
      SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, nlptimelimit, conshdlrdata->nlpnodelim) );
      conshdlrdata->timeleft -= SCIPgetSolvingTime(scip);
   }

   SCIPdebugMsg(scip, "packable status? %d\n", SCIPpatternGetPackableStatus(pattern));
   SCIPcheckPattern(scip, probdata, pattern);

   return SCIP_OKAY;
}

/** auxiliary function for enforcing ringpacking constraint; the function checks whether
 *
 *  1. the solution is feasible; if yes -> skip
 *  2. tries to verify an unverified circular pattern C with z*_c > 0
 *     2a. case packable or unknown: go to 2.
 *     2b. case not packable: fix z_C to 0 -> skip
 *  3. fix all unverified circular patterns to 0
 *
 *  Note that after step 3. the dual bound is not valid anymore.
 */
static
SCIP_RETCODE enforceSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< solution (NULL for LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;
   SCIP_PATTERN** cpatterns;
   SCIP_VAR** cvars;
   int ncpatterns;
   int p;

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "enforce solution:\n");
      SCIP_CALL( SCIPprintSol(scip, sol, NULL, TRUE) );
#endif

   *result = SCIP_FEASIBLE;

   /* (1.) check whether the solution is already feasible */
   if( isSolFeasible(scip, sol) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get circular pattern information */
   SCIPprobdataGetCInfos(probdata, &cpatterns, &cvars, &ncpatterns);
   assert(cpatterns != NULL);
   assert(cvars != NULL);
   assert(ncpatterns > 0);

   /* (2.) try to verify a pattern */
   for( p = 0; p < ncpatterns; ++p )
   {
      SCIP_Real solval;
      SCIP_Bool infeasible;
      SCIP_Bool success;

      assert(cpatterns[p] != NULL);
      assert(cvars[p] != NULL);

      solval = SCIPgetSolVal(scip, sol, cvars[p]);

      /* skip unused circular patterns */
      if( SCIPisFeasZero(scip, solval) )
         continue;

      /* try to verify an unknown circular pattern */
      if( SCIPpatternGetPackableStatus(cpatterns[p]) == SCIP_PACKABLE_UNKNOWN && !conshdlrdata->tried[p] )
      {
         SCIP_CALL( verifyCircularPattern(scip, conshdlrdata, probdata, cpatterns[p]) );
         conshdlrdata->tried[p] = TRUE;
      }

      /* (2a.) fix corresponding variable to zero if pattern is not packable */
      if( SCIPpatternGetPackableStatus(cpatterns[p]) == SCIP_PACKABLE_NO )
      {
         SCIP_CALL( SCIPfixVar(scip, cvars[p], 0.0, &infeasible, &success) );
         SCIPdebugMsg(scip, "fix unpackable pattern %d\n", p);

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if( success )
         {
            /* stop since we cutoff the LP solution */
            *result = SCIP_REDUCEDDOM;
            return SCIP_OKAY;
         }
      }
   }

   SCIPdebugMsg(scip, "fix all tested but still unverified circular patterns\n");

   /* (3.) fix an unverified patterns that is used */
   for( p = 0; p < ncpatterns; ++p )
   {
      if( SCIPpatternGetPackableStatus(cpatterns[p]) == SCIP_PACKABLE_UNKNOWN )
      {
         SCIP_Bool infeasible;
         SCIP_Bool success;

         SCIP_CALL( SCIPfixVar(scip, cvars[p], 0.0, &infeasible, &success) );
         SCIPdebugMsg(scip, "fix unknown pattern %d in [%g,%g] (success=%u)\n", p, SCIPvarGetLbLocal(cvars[p]),
            SCIPvarGetUbLocal(cvars[p]), success);

         /* dual bound is not valid anymore */
         SCIPprobdataInvalidateDualbound(scip, probdata);

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if( success )
            *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}

/** get shading of a given pattern type */
static
int getShadingVal(
   int                   type,              /**< pattern type */
   int                   ntypes             /**< total number of patterns */
   )
{
   assert(type >= 0);
   assert(type < ntypes);

   return (int)MIN(100, MAX(ntypes, 100) / (type+1));
}

/*
 * Callback methods of event handler
 */

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{  /*lint --e{715}*/
   SCIP_PATTERN** patterns;
   SCIP_VAR** vars;
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   char* filename;
   int npatterns;
   int p;

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   /* check whether new solution is indeed feasible */
#ifndef NDEBUG
   {
      /* check circular patterns */
      SCIPprobdataGetCInfos(probdata, &patterns, &vars, &npatterns);
      assert(npatterns > 0);

      for( p = 0; p < npatterns; ++p )
      {
         SCIP_Real val;

         assert(patterns[p] != NULL);
         assert(vars[p] != NULL);

         val = SCIPgetSolVal(scip, sol, vars[p]);

         /* pattern is either not used or packable */
         assert(SCIPisFeasZero(scip, val) || SCIPpatternGetPackableStatus(patterns[p]) == SCIP_PACKABLE_YES);
         SCIPcheckPattern(scip, probdata, patterns[p]);
      }

      /* check rectangular patterns */
      SCIPprobdataGetRInfos(probdata, &patterns, &vars, &npatterns);
      assert(npatterns > 0);

      for( p = 0; p < npatterns; ++p )
         SCIPcheckPattern(scip, probdata, patterns[p]);
   }
#endif

   /* write best solution information into a tex file */
   SCIP_CALL( SCIPgetStringParam(scip, "ringpacking/texfilename", &filename) );

   if( strcmp(filename, "") != 0 )
   {
      FILE* file;
      SCIP_Real* rexts;
      SCIP_Real* rints;
      int* demands;
      SCIP_Real width;
      SCIP_Real height;
      int ntypes;

      rexts = SCIPprobdataGetRexts(probdata);
      rints = SCIPprobdataGetRints(probdata);
      demands = SCIPprobdataGetDemands(probdata);
      width = SCIPprobdataGetWidth(probdata);
      height = SCIPprobdataGetHeight(probdata);
      ntypes = SCIPprobdataGetNTypes(probdata);

      SCIPdebugMsg(scip, "write best solution into %s\n", filename);

      file = fopen(filename, "w");
      assert(file != NULL);

      /* latex header */
      SCIPinfoMessage(scip, file, "\\documentclass[preview]{standalone}\n");
      SCIPinfoMessage(scip, file, "\\usepackage{tikz}\n");
      SCIPinfoMessage(scip, file, "\\usepackage{xstring}\n\n");
      SCIPinfoMessage(scip, file, "\\begin{document}\n\n");

      /* circular patterns */
      SCIPinfoMessage(scip, file, "\\section*{circular patterns}\n\n");
      SCIPprobdataGetCInfos(probdata, &patterns, &vars, &npatterns);
      for( p = 0; p < npatterns; ++p )
      {
         if( SCIPpatternGetPackableStatus(patterns[p]) == SCIP_PACKABLE_YES )
         {
            int type = SCIPpatternGetCircleType(patterns[p]);
            int i;

            SCIPinfoMessage(scip, file, "\\StrSubstitute{%s}{_}{-}[\\pname]\n", SCIPvarGetName(vars[p]));
            SCIPinfoMessage(scip, file, "\\subsection*{\\texttt{\\pname} = %g demand=%d}\n",
               SCIPgetSolVal(scip, sol, vars[p]), demands[type]);
            SCIPinfoMessage(scip, file, "\\resizebox{0.3\\textwidth}{!}{\n");
            SCIPinfoMessage(scip, file, "\\begin{tikzpicture}\n");
            SCIPinfoMessage(scip, file, "\\draw[draw=none,fill=black!%d!white] (%g,%g) circle (%g);\n",
               getShadingVal(type, ntypes), 0.0, 0.0, rexts[type]);
            SCIPinfoMessage(scip, file, "\\draw[draw=none,fill=white] (%g,%g) circle (%g);\n", 0.0, 0.0, rints[type]);

            for( i = 0; i < SCIPpatternGetNElemens(patterns[p]); ++i )
            {
               int elemtype = SCIPpatternGetElementType(patterns[p], i);
               SCIP_Real x = SCIPpatternGetElementPosX(patterns[p], i);
               SCIP_Real y = SCIPpatternGetElementPosY(patterns[p], i);
               SCIP_Real _rext = rexts[elemtype];
               SCIP_Real _rint = rints[elemtype];

               SCIPinfoMessage(scip, file, "\\draw[draw=none,fill=black!%d!white] (%g,%g) circle (%g);\n",
                  getShadingVal(elemtype, ntypes), x, y, _rext);
               SCIPinfoMessage(scip, file, "\\draw[draw=none,fill=white] (%g,%g) circle (%g);\n", x, y, _rint);
            }

            SCIPinfoMessage(scip, file, "\\end{tikzpicture}\n");
            SCIPinfoMessage(scip, file, "}\n\n");
         }
      }

      /* rectangular patterns */
      SCIPinfoMessage(scip, file, "\\section*{rectangular patterns}\n\n");
      SCIPprobdataGetRInfos(probdata, &patterns, &vars, &npatterns);
      for( p = 0; p < npatterns; ++p )
      {
         int i;

         assert(SCIPpatternGetPackableStatus(patterns[p]) == SCIP_PACKABLE_YES);

         SCIPinfoMessage(scip, file, "\\StrSubstitute{%s}{_}{-}[\\pname]\n", SCIPvarGetName(vars[p]));
         SCIPinfoMessage(scip, file, "\\subsection*{\\texttt{\\pname} = %g}\n", SCIPgetSolVal(scip, sol, vars[p]));
         SCIPinfoMessage(scip, file, "\\resizebox{0.3\\textwidth}{!}{\n");
         SCIPinfoMessage(scip, file, "\\begin{tikzpicture}\n");

         for( i = 0; i < SCIPpatternGetNElemens(patterns[p]); ++i )
         {
            int elemtype = SCIPpatternGetElementType(patterns[p], i);
            SCIP_Real x = SCIPpatternGetElementPosX(patterns[p], i);
            SCIP_Real y = SCIPpatternGetElementPosY(patterns[p], i);
            SCIP_Real _rext = rexts[elemtype];

            SCIPinfoMessage(scip, file, "\\draw[draw=none,fill=black!%d!white] (%g,%g) circle (%g);\n",
               getShadingVal(elemtype, ntypes), x, y, _rext);
            /* SCIPinfoMessage(scip, file, "\\draw[draw=none,fill=white] (%g,%g) circle (%g);\n", x, y, _rint); */
         }

         SCIPinfoMessage(scip, file, "\\draw[] (%g,%g) -- (%g,%g) -- (%g,%g) -- (%g,%g) -- cycle;\n",
            0.0, 0.0, 0.0, height, width, height, width, 0.0);

         SCIPinfoMessage(scip, file, "\\end{tikzpicture}\n");
         SCIPinfoMessage(scip, file, "}\n\n");
      }

      SCIPinfoMessage(scip, file, "\\end{document}\n");

      fclose(file);
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeRpa)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->locked != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->locked, conshdlrdata->lockedsize);
      conshdlrdata->lockedsize = 0;
   }

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolRpa)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;
   int ncpatterns;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->tried == NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIPprobdataGetCInfos(probdata, NULL, NULL, &ncpatterns);
   assert(ncpatterns > 0);

   /* allocate memory to remember which patterns have been tested to be packable */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->tried, ncpatterns) );
   BMSclearMemoryArray(conshdlrdata->tried, ncpatterns);

   /* catch events for new solutions */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, conshdlrdata->eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolRpa)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;
   int ncpatterns;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->tried != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIPprobdataGetCInfos(probdata, NULL, NULL, &ncpatterns);
   assert(ncpatterns > 0);

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->tried, ncpatterns);

   /* free events for new solutions */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, conshdlrdata->eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpRpa)
{  /*lint --e{715}*/
   SCIP_CALL( enforceSol(scip, conshdlr, NULL, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxRpa)
{  /*lint --e{715}*/
   SCIP_CALL( enforceSol(scip, conshdlr, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsRpa)
{  /*lint --e{715}*/
   SCIP_CALL( enforceSol(scip, conshdlr, NULL, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckRpa)
{  /*lint --e{715}*/
   *result = isSolFeasible(scip, sol) ? SCIP_FEASIBLE : SCIP_INFEASIBLE;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockRpa)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;
   SCIP_PATTERN** cpatterns;
   SCIP_VAR** cvars;
   SCIP_Bool first;
   int ncpatterns;
   int p;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get circular patterns and corresponding variables */
   SCIPprobdataGetCInfos(probdata, &cpatterns, &cvars, &ncpatterns);
   assert(cpatterns != NULL);
   assert(cvars != NULL);
   assert(ncpatterns > 0);

   /* remember whether we have locked the variables for the first time */
   if( conshdlrdata->locked == NULL )
   {
      first = TRUE;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->locked, ncpatterns) );
      BMSclearMemoryArray(conshdlrdata->locked, ncpatterns);
      conshdlrdata->lockedsize = ncpatterns;
   }
   else
      first = FALSE;

   /*
    * A pattern might change its status during a later verification step and we only need to lock patterns with a
    * SCIP_PACKABLE_UNKNOWN status. For this reason, we keep track of patterns that have been locked. The CONSLOCK
    * callback should only be called twice because the constraint handler does not have constraints.
    */
   for( p = 0; p < ncpatterns; ++p )
   {
      assert(cpatterns[p] != NULL);
      assert(cvars[p] != NULL);

      if( first && SCIPpatternGetPackableStatus(cpatterns[p]) == SCIP_PACKABLE_UNKNOWN )
      {
         assert(!conshdlrdata->locked[p]);
         assert(nlocksneg + nlockspos > 0);

         SCIP_CALL( SCIPaddVarLocksType(scip, cvars[p], SCIP_LOCKTYPE_MODEL, nlocksneg + nlockspos,
               nlocksneg + nlockspos) );
         conshdlrdata->locked[p] = TRUE;
         SCIPdebugMsg(scip, "lock %s\n", SCIPvarGetName(cvars[p]));
      }
      else if( !first && conshdlrdata->locked[p] )
      {
         assert(nlocksneg + nlockspos < 0);

         SCIP_CALL( SCIPaddVarLocksType(scip, cvars[p], SCIP_LOCKTYPE_MODEL, nlocksneg + nlockspos,
               nlocksneg + nlockspos) );
         conshdlrdata->locked[p] = FALSE;
         SCIPdebugMsg(scip, "unlock %s\n", SCIPvarGetName(cvars[p]));
      }
   }

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */


/** creates the handler for ringpacking */
SCIP_RETCODE SCIPincludeConshdlrRpa(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpRpa, consEnfopsRpa, consCheckRpa, consLockRpa,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolRpa) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeRpa) );;
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolRpa) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxRpa) );

   /* add event handler for new solutios */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         processNewSolutionEvent, NULL) );

   /* add verification parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
      "ringpacking/verification/nlptilim",
      "time limit for verification NLP",
      &conshdlrdata->nlptilim, FALSE, DEFAULT_VERIFICATION_NLPTILIM, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
      "ringpacking/verification/nlpnodelim",
      "node limit for verification NLP",
      &conshdlrdata->nlpnodelim, FALSE, DEFAULT_VERIFICATION_NLPNODELIM, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
      "ringpacking/verification/heurtilim",
      "time limit for heuristic verification",
      &conshdlrdata->heurtilim, FALSE, DEFAULT_VERIFICATION_HEURTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "ringpacking/verification/heuriterlim",
      "iteration limit for heuristic verification",
      &conshdlrdata->heuriterlim, FALSE, DEFAULT_VERIFICATION_HEURITERLIM, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
      "ringpacking/verification/totaltilim",
      "total time limit for all verification problems during the solving process",
      &conshdlrdata->timeleft, FALSE, DEFAULT_VERIFICATION_TOTALTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
