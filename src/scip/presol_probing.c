/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_probing.c
 * @ingroup PRESOLVERS
 * @brief  probing presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_probing.h"


#define PRESOL_NAME            "probing"
#define PRESOL_DESC            "probing presolver on binary variables"
#define PRESOL_PRIORITY         -100000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               TRUE /**< should presolver be delayed, if other presolvers found reductions? */

#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */




/*
 * Default parameter settings
 */

#define DEFAULT_MAXRUNS              1  /**< maximal number of runs, probing participates in (-1: no limit) */
#define DEFAULT_PROPROUNDS          -1  /**< maximal number of propagation rounds in probing subproblems */
#define DEFAULT_MAXFIXINGS          50  /**< maximal number of fixings found, until probing is interrupted
                                         *   (0: don't interrupt) */
#define DEFAULT_MAXUSELESS        2000  /**< maximal number of successive probings without fixings,
                                         *   until probing is aborted (0: don't abort) */
#define DEFAULT_MAXTOTALUSELESS    100  /**< maximal number of successive probings without fixings, bound changes,
                                         *   and implications, until probing is aborted (0: don't abort) */
#define DEFAULT_MAXSUMUSELESS        0  /**< maximal number of probings without fixings, until probing is aborted
                                         *   (0: don't abort) */




/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   SCIP_VAR**            sortedvars;         /**< problem variables sorted by number of rounding locks */
   int                   nsortedvars;        /**< number of problem variables */
   int                   nsortedbinvars;     /**< number of binary problem variables */
   int                   maxruns;            /**< maximal number of runs, probing participates in (-1: no limit) */
   int                   proprounds;         /**< maximal number of propagation rounds in probing subproblems */
   int                   maxfixings;         /**< maximal number of fixings found, until probing is interrupted
                                              *   (0: don't interrupt) */
   int                   maxuseless;         /**< maximal number of successive probings without fixings,
                                              *   until probing is aborted (0: don't abort) */
   int                   maxtotaluseless;    /**< maximal number of successive probings without fixings, bound changes,
                                              *   and implications, until probing is aborted (0: don't abort) */
   int                   maxsumuseless;      /**< maximal number of probings without fixings, until probing is aborted
                                              *   (0: don't abort) */
   int                   startidx;           /**< starting variable index of next call */
   int                   lastsortstartidx;   /**< last starting variable index where the variables have been sorted */
   int                   nfixings;           /**< total number of fixings found in probing */
   int                   naggregations;      /**< total number of aggregations found in probing */
   int                   nimplications;      /**< total number of implications found in probing */
   int                   nbdchgs;            /**< total number of bound changes found in probing */
   int                   nuseless;           /**< current number of successive useless probings */
   int                   ntotaluseless;      /**< current number of successive totally useless probings */
   int                   nsumuseless;        /**< current number of useless probings */
   SCIP_Bool             called;             /**< was probing applied at least once? */
};




/*
 * Local methods
 */

/** initializes the presolver data */
static
void initPresoldata(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   presoldata->sortedvars = NULL;
   presoldata->nsortedvars = 0;
   presoldata->nsortedbinvars = 0;
   presoldata->startidx = 0;
   presoldata->lastsortstartidx = -1;
   presoldata->nfixings = 0;
   presoldata->naggregations = 0;
   presoldata->nimplications = 0;
   presoldata->nbdchgs = 0;
   presoldata->nuseless = 0;
   presoldata->ntotaluseless = 0;
   presoldata->nsumuseless = 0;
}

/** frees the sorted vars array */
static
SCIP_RETCODE freeSortedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   if( presoldata->sortedvars != NULL )
   {
      int i;

      /* release variables */
      for( i = 0; i < presoldata->nsortedvars; ++i )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &presoldata->sortedvars[i]) );
      }
      SCIPfreeMemoryArray(scip, &presoldata->sortedvars);
      presoldata->nsortedvars = 0;
      presoldata->nsortedbinvars = 0;
   }

   return SCIP_OKAY;
}

/** sorts the binary variables starting with the given index by rounding locks and implications */
static
SCIP_RETCODE sortVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   int                   firstidx            /**< first index that should be subject to sorting */
   )
{
   SCIP_VAR** vars;
   int* scores;
   int nvars;
   int i;

   assert(presoldata != NULL);
   assert(presoldata->sortedvars != NULL);

   vars = &presoldata->sortedvars[firstidx];
   nvars = presoldata->nsortedbinvars - firstidx;
   if( nvars <= 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("resorting probing variables %d to %d\n", firstidx, presoldata->nsortedbinvars-1);

   /* sort the variables by number of rounding locks and implications */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      var = vars[i];

      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      if( SCIPvarIsActive(var) )
         scores[i] = SCIPvarGetNLocksDown(var) + SCIPvarGetNLocksUp(var)
            + SCIPvarGetNImpls(var, FALSE) + SCIPvarGetNImpls(var, TRUE)
            + 5*(SCIPvarGetNCliques(var, FALSE) + SCIPvarGetNCliques(var, TRUE));
      else
         scores[i] = -1;
   }

   SCIPsortDownIntPtr(scores, (void**) vars, nvars);

   SCIPfreeBufferArray(scip, &scores);

   return SCIP_OKAY;
}

/** applies and evaluates probing of a single variable in the given direction */
static
SCIP_RETCODE applyProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   probingpos,         /**< variable number to apply probing on */
   SCIP_Bool             probingdir,         /**< value to fix probing variable to */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   assert(presoldata != NULL);
   assert(impllbs != NULL);
   assert(implubs != NULL);
   assert(proplbs != NULL);
   assert(propubs != NULL);
   assert(cutoff != NULL);
   assert(0 <= probingpos && probingpos < nvars);
   assert(SCIPvarGetType(vars[probingpos]) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarGetLbLocal(vars[probingpos]) < 0.5);
   assert(SCIPvarGetUbLocal(vars[probingpos]) > 0.5);

   SCIPdebugMessage("applying probing on variable <%s> == %u (nlocks=%d/%d, impls=%d/%d, clqs=%d/%d)\n",
      SCIPvarGetName(vars[probingpos]), probingdir,
      SCIPvarGetNLocksDown(vars[probingpos]), SCIPvarGetNLocksUp(vars[probingpos]),
      SCIPvarGetNImpls(vars[probingpos], FALSE), SCIPvarGetNImpls(vars[probingpos], TRUE),
      SCIPvarGetNCliques(vars[probingpos], FALSE), SCIPvarGetNCliques(vars[probingpos], TRUE));

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* fix variable */
   if( probingdir == FALSE )
   {
      SCIP_CALL( SCIPchgVarUbProbing(scip, vars[probingpos], 0.0) );
   }
   else
   {
      SCIP_CALL( SCIPchgVarLbProbing(scip, vars[probingpos], 1.0) );
   }

   /* apply propagation of implication graph and clique table */
   SCIP_CALL( SCIPpropagateProbingImplications(scip, cutoff) );
   if( !(*cutoff) )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         impllbs[i] = SCIPvarGetLbLocal(vars[i]);
         implubs[i] = SCIPvarGetUbLocal(vars[i]);
      }
   }

   /* apply propagation */
   if( !(*cutoff) )
   {
      SCIP_CALL( SCIPpropagateProbing(scip, presoldata->proprounds, cutoff, NULL) );
   }

   /* evaluate propagation */
   if( !(*cutoff) )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         proplbs[i] = SCIPvarGetLbLocal(vars[i]);
         propubs[i] = SCIPvarGetUbLocal(vars[i]);
#if 0
#ifdef SCIP_DEBUG
         if( SCIPisGT(scip, proplbs[i], SCIPvarGetLbGlobal(vars[i])) )
         {
            SCIPdebugMessage(" -> <%s>[%g,%g] >= %g\n", SCIPvarGetName(vars[i]), 
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), proplbs[i]);
         }
         if( SCIPisLT(scip, propubs[i], SCIPvarGetUbGlobal(vars[i])) )
         {
            SCIPdebugMessage(" -> <%s>[%g,%g] <= %g\n", SCIPvarGetName(vars[i]), 
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), propubs[i]);
         }
#endif
#endif
      }
   }

   /* exit probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}




/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyProbing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolProbing(scip) );
 
   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeProbing)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   assert(presoldata->sortedvars == NULL);
   assert(presoldata->nsortedvars == 0);
   assert(presoldata->nsortedbinvars == 0);

   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitProbing)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   initPresoldata(presoldata);

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitProbing)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIP_CALL( freeSortedvars(scip, presoldata) );
   assert(presoldata->sortedvars == NULL);
   assert(presoldata->nsortedvars == 0);
   assert(presoldata->nsortedbinvars == 0);

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreProbing)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   presoldata->called = FALSE;

   return SCIP_OKAY;
}


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreProbing)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* delete the vars array, if the maximal number of runs are exceeded */
   if( presoldata->maxruns >= 0 && SCIPgetNRuns(scip) >= presoldata->maxruns )
   {
      SCIP_CALL( freeSortedvars(scip, presoldata) );
      assert(presoldata->sortedvars == NULL);
      assert(presoldata->nsortedvars == 0);
      assert(presoldata->nsortedbinvars == 0);
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecProbing)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_VAR** vars;
   SCIP_Real* zeroimpllbs;
   SCIP_Real* zeroimplubs;
   SCIP_Real* zeroproplbs;
   SCIP_Real* zeropropubs;
   SCIP_Real* oneimpllbs;
   SCIP_Real* oneimplubs;
   SCIP_Real* oneproplbs;
   SCIP_Real* onepropubs;
   int nvars;
   int nbinvars;
   int maxfixings;
   int maxuseless;
   int maxtotaluseless;
   int maxsumuseless;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldnimplications;
   int localnfixedvars;
   int localnaggrvars;
   int localnchgbds;
   int localnimplications;
   int i;
   SCIP_Bool delay;
   SCIP_Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* check, if probing should be applied in the current run */
   if( presoldata->maxruns >= 0 && SCIPgetNRuns(scip) > presoldata->maxruns )
      return SCIP_OKAY;

   /* if no domains changed since the last call, we don't need to probe */
   if( presoldata->called && nnewfixedvars == 0 && nnewaggrvars == 0 && nnewchgbds == 0 && nnewholes == 0 )
      return SCIP_OKAY;

   /**@todo currently we only perform probing on variables which are of binary type; there are, however, more
    *       implicit binary variables which can be found with the method SCIPvarIsBinary(); for these variable we could
    *       also perform probing */

   SCIPdebugMessage("executing probing (used %.1f sec)\n", SCIPpresolGetTime(presol));

   *result = SCIP_DIDNOTFIND;

   /* allow some additional probings */
   presoldata->nuseless -= presoldata->nuseless/10;
   presoldata->ntotaluseless -= presoldata->ntotaluseless/10;

   /* get variable data */
   if( presoldata->sortedvars == NULL )
   {
      SCIP_VAR** probvars;

      assert(presoldata->startidx == 0);

      SCIP_CALL( SCIPgetVarsData(scip, &probvars, &nvars, &nbinvars, NULL, NULL, NULL) );
      if( nbinvars == 0 )
         return SCIP_OKAY;

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &presoldata->sortedvars, probvars, nvars) );
      presoldata->nsortedvars = nvars;
      presoldata->nsortedbinvars = nbinvars;

      /* capture variables to make sure, the variables are not deleted */
      for( i = 0; i < presoldata->nsortedvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, presoldata->sortedvars[i]) );
      }
   }
   else
      nbinvars = SCIPgetNBinVars(scip);

   /* if we probed all binary variables in previous runs, start again with the first one */
   if( !presoldata->called && presoldata->startidx >= nbinvars )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "   (%.1fs) probing cycle finished: starting next cycle\n", SCIPgetSolvingTime(scip));
      presoldata->startidx = 0;
      presoldata->lastsortstartidx = -1;
      presoldata->nuseless = 0;
      presoldata->ntotaluseless = 0;
   }
   presoldata->called = TRUE;

   /* sort the binary variables by number of rounding locks, if at least 100 variables were probed since last sort */
   if( presoldata->lastsortstartidx < 0 || presoldata->startidx - presoldata->lastsortstartidx >= 100 )
   {
      SCIP_CALL( sortVariables(scip, presoldata, presoldata->startidx) );
      presoldata->lastsortstartidx = presoldata->startidx;
   }

   vars = presoldata->sortedvars;
   nvars = presoldata->nsortedvars;
   nbinvars = presoldata->nsortedbinvars;
   assert(vars != NULL);
   assert(nbinvars > 0);

   maxfixings = (presoldata->maxfixings > 0 ? presoldata->maxfixings : INT_MAX);
   maxuseless = (presoldata->maxuseless > 0 ? presoldata->maxuseless : INT_MAX);
   maxtotaluseless = (presoldata->maxtotaluseless > 0 ? presoldata->maxtotaluseless : INT_MAX);
   maxsumuseless = (presoldata->maxsumuseless > 0 ? presoldata->maxsumuseless : INT_MAX);
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldnimplications = presoldata->nimplications;

   /* get temporary memory for storing probing results */
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimpllbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimplubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroproplbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeropropubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimpllbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimplubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneproplbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &onepropubs, nvars) );

   /* for each binary variable, probe fixing the variable to zero and one */
   delay = FALSE;
   cutoff = FALSE;
   for( i = presoldata->startidx; i < nbinvars && !cutoff; ++i )
   {
      SCIP_Bool localcutoff;
      SCIP_Bool probingzero;
      SCIP_Bool probingone;

      /* check whether probing should be aborted */
      if( SCIPisStopped(scip) || presoldata->nuseless >= maxuseless
         || presoldata->ntotaluseless >= maxtotaluseless || presoldata->nsumuseless >= maxsumuseless )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n", 
            SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(SCIP_Real)(i+1)/(SCIP_Real)nbinvars,
            presoldata->nfixings, presoldata->naggregations, presoldata->nimplications, presoldata->nbdchgs);

         if( SCIPisStopped(scip) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: solving stopped\n", SCIPgetSolvingTime(scip));
         }
         else if( presoldata->nuseless >= maxuseless )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: %d/%d successive useless probings\n", SCIPgetSolvingTime(scip), 
               presoldata->nuseless, maxuseless);
         }
         else if( presoldata->ntotaluseless >= maxtotaluseless )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: %d/%d successive totally useless probings\n", SCIPgetSolvingTime(scip), 
               presoldata->ntotaluseless, maxtotaluseless);
         }
         else if( presoldata->nsumuseless >= maxsumuseless )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: %d/%d useless probings in total\n", SCIPgetSolvingTime(scip), 
               presoldata->nsumuseless, maxsumuseless);
         }
         break;
      }

      /* check if we already fixed enough variables for this round */
      if( *nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars >= maxfixings )
      {
         delay = TRUE;
         break;
      }

      /* display probing status */
      if( (i+1) % 100 == 0 )
      {
         SCIP_VERBLEVEL verblevel;

         verblevel = ((i+1) % 1000 == 0 ? SCIP_VERBLEVEL_HIGH : SCIP_VERBLEVEL_FULL);
         SCIPverbMessage(scip, verblevel, NULL,
            "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n", 
            SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(SCIP_Real)(i+1)/(SCIP_Real)nbinvars,
            presoldata->nfixings, presoldata->naggregations, presoldata->nimplications, presoldata->nbdchgs);
      }

      /* ignore variables, that were fixed, aggregated, or deleted in prior probings */
      if( !SCIPvarIsActive(vars[i]) || SCIPvarIsDeleted(vars[i])
         || SCIPvarGetLbGlobal(vars[i]) > 0.5 || SCIPvarGetUbGlobal(vars[i]) < 0.5 )
         continue;

      if( presoldata->nuseless > 0 )
         presoldata->nsumuseless++;
      else
         presoldata->nsumuseless = MAX(presoldata->nsumuseless-1, 0);
      presoldata->nuseless++;
      presoldata->ntotaluseless++;

      /* determine whether zero probing should happen */
      probingzero = TRUE;
      if ( SCIPvarGetNLocksDown(vars[i]) == 0 )
         probingzero = FALSE;

      if ( probingzero )
      {
         /* apply probing for fixing the variable to zero */
         SCIP_CALL( applyProbing(scip, presoldata, vars, nvars, i, FALSE, zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs,
               &localcutoff) );
         
         if( localcutoff )
         {
            SCIP_Bool fixed;
            
            /* the variable can be fixed to TRUE */
            SCIPdebugMessage("fixing probing variable <%s> to 1.0, nlocks=(%d/%d)\n",
               SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
            SCIP_CALL( SCIPfixVar(scip, vars[i], 1.0, &cutoff, &fixed) );
            assert(fixed);
            (*nfixedvars)++;
            presoldata->nfixings++;
            presoldata->nuseless = 0;
            presoldata->ntotaluseless = 0;
            continue; /* don't try upwards direction, because the variable is already fixed */
         }

         /* ignore variables, that were fixed, aggregated, or deleted in prior probings
          * (propagators in zero-probe might have found global fixings but did not trigger the localcutoff)
          */
         if( !SCIPvarIsActive(vars[i]) || SCIPvarIsDeleted(vars[i])
            || SCIPvarGetLbGlobal(vars[i]) > 0.5 || SCIPvarGetUbGlobal(vars[i]) < 0.5 )
            continue;
      }

      /* determine whether one probing should happen */
      probingone = TRUE;
      if ( SCIPvarGetNLocksUp(vars[i]) == 0 )
         probingone = FALSE;

      if ( probingone )
      {
         /* apply probing for fixing the variable to one */
         SCIP_CALL( applyProbing(scip, presoldata, vars, nvars, i, TRUE, oneimpllbs, oneimplubs, oneproplbs, onepropubs,
               &localcutoff) );
         
         if( localcutoff )
         {
            SCIP_Bool fixed;
            
            /* the variable can be fixed to FALSE */
            SCIPdebugMessage("fixing probing variable <%s> to 0.0, nlocks=(%d/%d)\n", 
               SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
            SCIP_CALL( SCIPfixVar(scip, vars[i], 0.0, &cutoff, &fixed) );
            assert(fixed);
            (*nfixedvars)++;
            presoldata->nfixings++;
            presoldata->nuseless = 0;
            presoldata->ntotaluseless = 0;
            continue; /* don't analyze probing deductions, because the variable is already fixed */
         }
      }

      /* not have to check deductions if only one probing direction has been checked */
      if ( ! probingzero || ! probingone )
         continue;

      /* analyze probing deductions */
      localnfixedvars    = 0;
      localnaggrvars     = 0;
      localnimplications = 0;
      localnchgbds       = 0;
      SCIP_CALL( SCIPanalyzeDeductionsProbing(scip, vars[i], 0.0, 1.0,
         nvars, vars, zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs, oneimpllbs, oneimplubs, oneproplbs, onepropubs,
         &localnfixedvars, &localnaggrvars, &localnimplications, &localnchgbds, &cutoff) );

      *nfixedvars += localnfixedvars;
      *naggrvars  += localnaggrvars;
      *nchgbds    += localnchgbds;
      presoldata->nfixings      += localnfixedvars;
      presoldata->naggregations += localnaggrvars;
      presoldata->nbdchgs       += localnchgbds;
      presoldata->nimplications += localnimplications;

      if( localnfixedvars > 0 || localnaggrvars > 0 )
      {
         presoldata->nuseless = 0;
         presoldata->ntotaluseless = 0;
      }
      if( localnimplications > 0 || localnchgbds > 0 )
         presoldata->ntotaluseless = 0;
   }

   presoldata->startidx = i;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &onepropubs);
   SCIPfreeBufferArray(scip, &oneproplbs);
   SCIPfreeBufferArray(scip, &oneimplubs);
   SCIPfreeBufferArray(scip, &oneimpllbs);
   SCIPfreeBufferArray(scip, &zeropropubs);
   SCIPfreeBufferArray(scip, &zeroproplbs);
   SCIPfreeBufferArray(scip, &zeroimplubs);
   SCIPfreeBufferArray(scip, &zeroimpllbs);

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
      *result = SCIP_DELAYED;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds
      || presoldata->nimplications > oldnimplications )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the probing presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create probing presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
   initPresoldata(presoldata);

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolCopyProbing,
         presolFreeProbing, presolInitProbing, presolExitProbing, 
         presolInitpreProbing, presolExitpreProbing, presolExecProbing,
         presoldata) );

   /* add probing presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/probing/maxruns", 
         "maximal number of runs, probing participates in (-1: no limit)",
         &presoldata->maxruns, FALSE, DEFAULT_MAXRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/probing/proprounds", 
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &presoldata->proprounds, TRUE, DEFAULT_PROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/probing/maxfixings",
         "maximal number of fixings found, until probing is interrupted (0: don't iterrupt)",
         &presoldata->maxfixings, TRUE, DEFAULT_MAXFIXINGS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/probing/maxuseless",
         "maximal number of successive probings without fixings, until probing is aborted (0: don't abort)",
         &presoldata->maxuseless, TRUE, DEFAULT_MAXUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/probing/maxtotaluseless",
         "maximal number of successive probings without fixings, bound changes, and implications, until probing is aborted (0: don't abort)",
         &presoldata->maxtotaluseless, TRUE, DEFAULT_MAXTOTALUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/probing/maxsumuseless",
         "maximal number of probings without fixings, until probing is aborted (0: don't abort)",
         &presoldata->maxsumuseless, TRUE, DEFAULT_MAXSUMUSELESS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** analyses boundchanges resulting from probing on a variable and performs deduced fixations, aggregations, and domain tightenings
 * Given a variable probingvar with domain [l,u] and bound tightening results from reducing the domain
 * once to [l,leftub] and once to [rightlb,u], the method computes and applies resulting variable fixations, aggregations,
 * implications, and bound changes. Variable probingvar does not need to be binary.
 * The whole domain of probingvar need to be covered by the left and right branches, i.e.,
 * we assume leftub >= rightlb for continuous variables or floor(leftub) >= ceil(rightlb)-1 for discrete variables.
 * Bounds after applying implications and cliques do not need to be provided, but if they are omitted and probingvar is a binary variable,
 * then already existing implications may be added.
 */
SCIP_RETCODE SCIPanalyzeDeductionsProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             probingvar,         /**< the probing variable */
   SCIP_Real             leftub,             /**< upper bound of probing variable in left branch */
   SCIP_Real             rightlb,            /**< lower bound of probing variable in right branch */
   int                   nvars,              /**< number of variables which bound changes should be analyzed */
   SCIP_VAR**            vars,               /**< variables which bound changes should be analyzed */
   SCIP_Real*            leftimpllbs,        /**< lower bounds after applying implications and cliques in left branch, or NULL */
   SCIP_Real*            leftimplubs,        /**< upper bounds after applying implications and cliques in left branch, or NULL */
   SCIP_Real*            leftproplbs,        /**< lower bounds after applying domain propagation in left branch */
   SCIP_Real*            leftpropubs,        /**< upper bounds after applying domain propagation in left branch */
   SCIP_Real*            rightimpllbs,       /**< lower bounds after applying implications and cliques in right branch, or NULL */
   SCIP_Real*            rightimplubs,       /**< upper bounds after applying implications and cliques in right branch, or NULL */
   SCIP_Real*            rightproplbs,       /**< lower bounds after applying domain propagation in right branch */
   SCIP_Real*            rightpropubs,       /**< upper bounds after applying domain propagation in right branch */
   int*                  nfixedvars,         /**< pointer to counter which is increased by the number of deduced variable fixations */
   int*                  naggrvars,          /**< pointer to counter which is increased by the number of deduced variable aggregations */
   int*                  nimplications,      /**< pointer to counter which is increased by the number of deduced implications */
   int*                  nchgbds,            /**< pointer to counter which is increased by the number of deduced bound tightenings */
   SCIP_Bool*            cutoff              /**< buffer to store whether a cutoff is detected */
   )
{
   SCIP_Bool fixedleft;
   SCIP_Bool fixedright;
   int j;

   assert(scip != NULL);
   assert(probingvar != NULL);
   assert(SCIPisGE(scip, leftub,  SCIPvarGetLbGlobal(probingvar))); /* left  branch should not be empty by default */
   assert(SCIPisLE(scip, rightlb, SCIPvarGetUbGlobal(probingvar))); /* right branch should not be empty by default */
   assert(vars != NULL || nvars == 0);
   assert(leftproplbs != NULL);
   assert(leftpropubs != NULL);
   assert(rightproplbs != NULL);
   assert(rightpropubs != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nimplications != NULL);
   assert(nchgbds != NULL);
   assert(cutoff != NULL);

   /* @todo the asserts below could be relaxed by taking domain holes into account */
   if( SCIPvarGetType(probingvar) != SCIP_VARTYPE_CONTINUOUS )
   {
      /* adjust bounds to actually used ones */
      leftub  = SCIPfloor(scip, leftub);
      rightlb = SCIPceil(scip, rightlb);
      /* assert dichotomy in case of discrete var: leftub >= rightlb - 1.0 */
      assert(SCIPisGE(scip, leftub, rightlb - 1.0));
   }
   else
   {
      /* assert dichotomy in case of continuous var: leftub >= rightlb */
      assert(SCIPisGE(scip, leftub, rightlb));
   }

   /* check if probing variable was fixed in the branches */
   fixedleft  = SCIPisEQ(scip, SCIPvarGetLbGlobal(probingvar), leftub);
   fixedright = SCIPisEQ(scip, SCIPvarGetUbGlobal(probingvar), rightlb);

   *cutoff = FALSE;

   for( j = 0; j < nvars && !*cutoff; ++j )
   {
      SCIP_Real newlb;
      SCIP_Real newub;

      /* @todo why not also look at probingvar? */
      if( vars[j] == probingvar )
         continue;

      /* new bounds of the variable is the union of the propagated bounds of the left and right case */
      newlb = MIN(leftproplbs[j], rightproplbs[j]);
      newub = MAX(leftpropubs[j], rightpropubs[j]);

      /* check for fixed variables */
      if( SCIPisFeasEQ(scip, newlb, newub) )
      {
         SCIP_Real fixval;
         SCIP_Bool fixed;

         /* in both probings, variable j is deduced to the same value: fix variable to this value */
         fixval = SCIPselectSimpleValue(newlb - SCIPepsilon(scip), newub + SCIPepsilon(scip), MAXDNOM);
         SCIP_CALL( SCIPfixVar(scip, vars[j], fixval, cutoff, &fixed) );
         if( fixed )
         {
            SCIPdebugMessage("fixed variable <%s> to %g due to probing on <%s> with nlocks=(%d/%d)\n",
               SCIPvarGetName(vars[j]), fixval,
               SCIPvarGetName(probingvar), SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
            (*nfixedvars)++;
         }
         continue;
      }
      else
      {
         /* check for bound tightenings */
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Bool tightened;

         oldlb = SCIPvarGetLbGlobal(vars[j]);
         oldub = SCIPvarGetUbGlobal(vars[j]);
         if( SCIPisLbBetter(scip, newlb, oldlb, oldub) )
         {
            /* in both probings, variable j is deduced to be at least newlb: tighten lower bound */
            SCIP_CALL( SCIPtightenVarLb(scip, vars[j], newlb, FALSE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage("tightened lower bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), oldlb, oldub, newlb,
                  SCIPvarGetName(probingvar), SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
               (*nchgbds)++;
            }
         }
         if( SCIPisUbBetter(scip, newub, oldlb, oldub) && !*cutoff )
         {
            /* in both probings, variable j is deduced to be at most newub: tighten upper bound */
            SCIP_CALL( SCIPtightenVarUb(scip, vars[j], newub, FALSE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage("tightened upper bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), oldlb, oldub, newub,
                  SCIPvarGetName(probingvar), SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
               (*nchgbds)++;
            }
         }
         if( *cutoff )
            break;
      }

      /* check for aggregations and implications */
      if( fixedleft  && SCIPisEQ(scip, leftproplbs[j],  leftpropubs[j]) &&
          fixedright && SCIPisEQ(scip, rightproplbs[j], rightpropubs[j]) )
      {
         /* vars[j] is fixed whenever probingvar is fixed, i.e.,
          *   vars[j] = leftproplbs[j] + (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub) * (probingvar - leftub)
          * -> both variables can be aggregated:
          *    (rightlb - leftub) * (vars[j] - leftproplbs[j]) = (rightproplbs[j] - leftproplbs[j]) * (probingvar - leftub)
          * -> (rightlb - leftub) * vars[j] - (rightproplbs[j] - leftproplbs[j]) * probingvar = leftproplbs[j] * rightlb - rightproplbs[j] * leftub
          *
          * check for case where both variables are binary: leftub = 1, rightlb = 0
          * case leftproplbs[j] = 0, rightproplbs[j] = 1, i.e., vars[j] and probingvar are fixed to same value
          *    -> aggregation is 1 * vars[j] - 1 * probingvar = 0 * 1 - 1 * 0 = 0 -> correct
          * case leftproplbs[j] = 1, rightproblbs[j] = 0, i.e., vars[j] and probingvar are fixed to oppositve values
          *    -> aggregation is 1 * vars[j] + 1 * probingvar = 1 * 1 - 0 * 0 = 0 -> correct
          */
         SCIP_Bool aggregated;
         SCIP_Bool redundant;

         SCIP_CALL( SCIPaggregateVars(scip, vars[j], probingvar,
            rightlb - leftub, -(rightproplbs[j] - leftproplbs[j]), leftproplbs[j] * rightlb - rightproplbs[j] * leftub,
            cutoff, &redundant, &aggregated) );

         if( aggregated )
         {
            SCIPdebugMessage("aggregated variables %g<%s> - %g<%s> == %g, nlocks=(%d/%d)\n",
               rightlb - leftub, SCIPvarGetName(vars[j]),
               rightproplbs[j] - leftproplbs[j], SCIPvarGetName(probingvar),
               leftproplbs[j] * rightlb - rightproplbs[j] * leftub,
               SCIPvarGetNLocksDown(vars[j]), SCIPvarGetNLocksUp(probingvar));
            (*naggrvars)++;
         }
      }
      else if( SCIPvarGetType(probingvar) == SCIP_VARTYPE_BINARY )
      {
         /* implications can be added only for binary variables */
         int nboundchanges;

         /* since probing var is binary variable, probing should have fixed variable in both branches,
          * which is to 0.0 in the left branch and to 1.0 in the right branch */
         assert(fixedleft);
         assert(fixedright);
         assert(SCIPisZero(scip, leftub));
         assert(SCIPisEQ(scip, rightlb, 1.0));

         if( SCIPisEQ(scip, newlb, leftpropubs[j]) && (leftimplubs == NULL || leftimplubs[j] > leftpropubs[j]) )
         {
            /* vars[j] is fixed to lower bound whenever probingvar is fixed to 0.0
             * and implication is not already known
             * -> insert implication: probingvar == 0  =>  vars[j] <= leftpropubs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), leftpropubs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, leftpropubs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         else if( SCIPisEQ(scip, newub, leftproplbs[j]) && (leftimpllbs == NULL || leftimpllbs[j] < leftproplbs[j]) )
         {
            /* vars[j] is fixed to upper bound whenever probingvar is fixed to 0.0
             * and implication is not already known
             * -> insert implication: probingvar == 0  =>  vars[j] >= leftproplbs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), leftproplbs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, leftproplbs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         /* we can do an else here, since the case where vars[j] is fixed for both fixings of probingvar had been handled as aggregation */
         else if( SCIPisEQ(scip, newlb, rightpropubs[j]) && (rightimplubs == NULL || rightimplubs[j] > rightpropubs[j]) )
         {
            /* vars[j] is fixed to lower bound whenever probingvar is fixed to 1.0
             * and implication is not already known
             * -> insert implication: probingvar == 1  =>  vars[j] <= rightpropubs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), rightpropubs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, rightpropubs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         else if( SCIPisEQ(scip, newub, rightproplbs[j]) && (rightimpllbs == NULL || rightimpllbs[j] < rightproplbs[j]) )
         {
            /* vars[j] is fixed to upper bound whenever probingvar is fixed to 1.0
             * and implication is not already known
             * -> insert implication: probingvar == 1  =>  vars[j] >= leftproplbs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), rightproplbs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, rightproplbs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         else if( SCIPvarGetType(vars[j]) != SCIP_VARTYPE_BINARY )
         {
            /* check for implications for lower or upper bounds (only store implications with bounds tightened at least by 0.5)
             * in case of binary variables, this should have been handled in the previous cases, since every boundchange also fixes the variable
             */
            if( leftpropubs[j] < newub - 0.5 && (leftimplubs == NULL || leftpropubs[j] < leftimplubs[j]) )
            {
               /* insert implication: probingvar == 0  =>  vars[j] <= leftpropubs[j] */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] <= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, leftpropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, leftpropubs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
            if( leftproplbs[j] > newlb + 0.5 && (leftimpllbs == NULL || leftproplbs[j] > leftimpllbs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 0  =>  vars[j] >= leftproplbs[j] */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] >= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, leftproplbs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, leftproplbs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
            if( rightpropubs[j] < newub - 0.5 && (rightpropubs == NULL || rightpropubs[j] < rightimplubs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 1  =>  vars[j] <= rightpropubs[j] */
               /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] <= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, rightpropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, rightpropubs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
            if( rightproplbs[j] > newlb + 0.5 && (rightproplbs == NULL || rightproplbs[j] > rightimpllbs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 1  =>  vars[j] >= rightproplbs[j] */
               /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] >= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, rightproplbs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, rightproplbs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
         }
      }
   }

   return SCIP_OKAY;
}
