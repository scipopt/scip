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
#pragma ident "@(#) $Id: presol_probing.c,v 1.63 2011/02/24 16:55:38 bzfpfets Exp $"

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
      int j;

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

      /* not not have to check deductions if only one probing direction has been checked */
      if ( ! probingzero || ! probingone )
         continue;

      /* analyze probing deductions */
      for( j = 0; j < nvars && !cutoff; ++j )
      {
         SCIP_Real newlb;
         SCIP_Real newub;
         SCIP_Bool fixed;

         if( j == i )
            continue;

         /* new bounds of the variable is the union of the propagated bounds of the zero and one case */
         newlb = MIN(zeroproplbs[j], oneproplbs[j]);
         newub = MAX(zeropropubs[j], onepropubs[j]);

         /* check for fixed variables */
         if( SCIPisFeasEQ(scip, newlb, newub) )
         {
            SCIP_Real fixval;

            /* in both probings, variable j is deduced to 0: fix variable to 0 */
            fixval = SCIPselectSimpleValue(newlb - SCIPepsilon(scip), newub + SCIPepsilon(scip), MAXDNOM);
            SCIP_CALL( SCIPfixVar(scip, vars[j], fixval, &cutoff, &fixed) );
            if( fixed )
            {
               SCIPdebugMessage("fixed variable <%s> to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), fixval, 
                  SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               (*nfixedvars)++;
               presoldata->nfixings++;
               presoldata->nuseless = 0;
               presoldata->ntotaluseless = 0;
            }
            continue;
         }

         if( j < nbinvars )
         {
            SCIP_Bool redundant;
            SCIP_Bool aggregated;
            int nboundchanges;

            assert(SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY);

            /* check for aggregations and implications on binary variables */
            if( zeropropubs[j] < 0.5 && oneproplbs[j] > 0.5 )
            {
               /* variable j is always deduced to the same value as probing variable i:
                * both variables can be aggregated with x_i - x_j == 0
                */
               SCIP_CALL( SCIPaggregateVars(scip, vars[i], vars[j], 1.0, -1.0, 0.0, &cutoff, &redundant, &aggregated) );
               if( aggregated )
               {
                  SCIPdebugMessage("aggregated variables <%s> == <%s>, nlocks=(%d/%d)\n",
                     SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), 
                     SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
                  (*naggrvars)++;
                  presoldata->naggregations++;
                  presoldata->nuseless = 0;
                  presoldata->ntotaluseless = 0;
               }
            }
            else if( zeroproplbs[j] > 0.5 && onepropubs[j] < 0.5 )
            {
               /* variable j is always deduced to the opposite value of probing variable i:
                * both variables can be aggregated with x_i + x_j == 1
                */
               SCIP_CALL( SCIPaggregateVars(scip, vars[i], vars[j], 1.0, 1.0, 1.0, &cutoff, &redundant, &aggregated) );
               if( aggregated )
               {
                  SCIPdebugMessage("aggregated variables <%s> == 1 - <%s>, nlocks=(%d/%d)\n",
                     SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]),
                     SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
                  (*naggrvars)++;
                  presoldata->naggregations++;
                  presoldata->nuseless = 0;
                  presoldata->ntotaluseless = 0;
               }
            }
            else if( zeropropubs[j] < 0.5 && zeroimplubs[j] > 0.5 ) /* ignore already existing implications */
            {
               /* insert implication: x_i == 0  =>  x_j == 0 */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s> == 0\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));*/
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, 0.0,
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
            else if( zeroproplbs[j] > 0.5 && zeroimpllbs[j] < 0.5 ) /* ignore already existing implications */
            {
               /* insert implication: x_i == 0  =>  x_j == 1 */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s> == 1\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));*/
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, 1.0,
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
            else if( onepropubs[j] < 0.5 && oneimplubs[j] > 0.5 ) /* ignore already existing implications */
            {
               /* insert implication: x_i == 1  =>  x_j == 0 */
               /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s> == 0\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));*/
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, 0.0,
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
            else if( oneproplbs[j] > 0.5 && oneimpllbs[j] < 0.5 ) /* ignore already existing implications */
            {
               /* insert implication: x_i == 1  =>  x_j == 1 */
               /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s> == 1\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));*/
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, 1.0, 
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
         }
         else
         {
            SCIP_Real oldlb;
            SCIP_Real oldub;
            SCIP_Bool tightened;
            int nboundchanges;

            assert(SCIPvarGetType(vars[j]) != SCIP_VARTYPE_BINARY);

            /* check for bound tightenings */
            oldlb = SCIPvarGetLbGlobal(vars[j]);
            oldub = SCIPvarGetUbGlobal(vars[j]);
            if( SCIPisLbBetter(scip, newlb, oldlb, oldub) )
            {
               /* in both probings, variable j is deduced to be at least newlb: tighten lower bound */
               SCIP_CALL( SCIPtightenVarLb(scip, vars[j], newlb, FALSE, &cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMessage("tightened lower bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                     SCIPvarGetName(vars[j]), oldlb, oldub, newlb,
                     SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
                  (*nchgbds)++;
                  presoldata->nbdchgs++;
                  /*presoldata->nuseless = 0;*/
                  presoldata->ntotaluseless = 0;
               }
            }
            if( SCIPisUbBetter(scip, newub, oldlb, oldub) && !cutoff )
            {
               /* in both probings, variable j is deduced to be at most newub: tighten upper bound */
               SCIP_CALL( SCIPtightenVarUb(scip, vars[j], newub, FALSE, &cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMessage("tightened upper bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                     SCIPvarGetName(vars[j]), oldlb, oldub, newub,
                     SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
                  (*nchgbds)++;
                  presoldata->nbdchgs++;
                  /*presoldata->nuseless = 0;*/
                  presoldata->ntotaluseless = 0;
               }
            }

            /* check for implications (only store implications with bounds tightened at least by 0.5) */
            if( zeropropubs[j] < newub - 0.5 && zeropropubs[j] < zeroimplubs[j] && !cutoff )
            {
               /* insert implication: x_i == 0  =>  x_j <= zeropropubs[j] */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] <= %g\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, zeropropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, zeropropubs[j],
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
            if( zeroproplbs[j] > newlb + 0.5 && zeroproplbs[j] > zeroimpllbs[j] && !cutoff )
            {
               /* insert implication: x_i == 0  =>  x_j >= zeroproplbs[j] */
               SCIPdebugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] >= %g\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, zeroproplbs[j]);
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, zeroproplbs[j],
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
            if( onepropubs[j] < newub - 0.5 && onepropubs[j] < oneimplubs[j] && !cutoff )
            {
               /* insert implication: x_i == 1  =>  x_j <= onepropubs[j] */
               SCIPdebugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] <= %g\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, onepropubs[j]);
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, onepropubs[j],
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
            if( oneproplbs[j] > newlb + 0.5 && oneproplbs[j] > oneimpllbs[j] && !cutoff )
            {
               /* insert implication: x_i == 1  =>  x_j >= oneproplbs[j] */
               SCIPdebugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] >= %g\n", 
                 SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, oneproplbs[j]);
               SCIP_CALL( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, oneproplbs[j],
                     &cutoff, &nboundchanges) );
               presoldata->nimplications++;
               (*nchgbds) += nboundchanges;
               presoldata->nbdchgs += nboundchanges;
               presoldata->ntotaluseless = 0;
            }
         }
      }
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
