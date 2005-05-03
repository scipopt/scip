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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_probing.c,v 1.15 2005/05/03 14:48:03 bzfpfend Exp $"

/**@file   presol_probing.c
 * @brief  probing presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

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
#define DEFAULT_MAXFIXINGS          10  /**< maximal number of fixings found, until probing is aborted (0: don't abort) */
#define DEFAULT_MAXUSELESS        2000  /**< maximal number of successive useless probings, until probing is aborted
                                         *   (0: don't abort) */




/*
 * Data structures
 */

/** presolver data */
struct PresolData
{
   VAR**            sortedvars;         /**< problem variables sorted by number of rounding locks */
   int              nsortedvars;        /**< number of problem variables */
   int              nsortedbinvars;     /**< number of binary problem variables */
   int              maxruns;            /**< maximal number of runs, probing participates in (-1: no limit) */
   int              proprounds;         /**< maximal number of propagation rounds in probing subproblems */
   int              maxfixings;         /**< maximal number of fixings found, until probing is aborted (0: don't abort) */
   int              maxuseless;         /**< maximal number of successive useless probings, until probing is aborted
                                         *   (0: don't abort) */
   int              startidx;           /**< starting variable index of next call */
   int              nfixings;           /**< total number of fixings found in probing */
   int              naggregations;      /**< total number of aggregations found in probing */
   int              nimplications;      /**< total number of implications found in probing */
   int              nbdchgs;            /**< total number of bound changes found in probing */
   Bool             called;             /**< was probing applied at least once? */
};




/*
 * Local methods
 */

/** applies and evaluates probing of a single variable in the given direction */
static
RETCODE applyProbing(
   SCIP*            scip,               /**< SCIP data structure */
   PRESOLDATA*      presoldata,         /**< presolver data */
   VAR**            vars,               /**< problem variables */
   int              nvars,              /**< number of problem variables */
   int              probingpos,         /**< variable number to apply probing on */
   Bool             probingdir,         /**< value to fix probing variable to */
   Real*            lbs,                /**< array to store propagated lower bounds */
   Real*            ubs,                /**< array to store propagated upper bounds */
   Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   assert(presoldata != NULL);
   assert(lbs != NULL);
   assert(ubs != NULL);
   assert(cutoff != NULL);
   assert(0 <= probingpos && probingpos < nvars);
   assert(SCIPvarGetType(vars[probingpos]) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarGetLbLocal(vars[probingpos]) < 0.5);
   assert(SCIPvarGetUbLocal(vars[probingpos]) > 0.5);

   debugMessage("applying probing on variable <%s> == %d\n", SCIPvarGetName(vars[probingpos]), probingdir);

   /* start probing mode */
   CHECK_OKAY( SCIPstartProbing(scip) );

   /* fix variable */
   if( probingdir == FALSE )
   {
      CHECK_OKAY( SCIPchgVarUbProbing(scip, vars[probingpos], 0.0) );
   }
   else
   {
      CHECK_OKAY( SCIPchgVarLbProbing(scip, vars[probingpos], 1.0) );
   }

   /* apply propagation */
   CHECK_OKAY( SCIPpropagateProbing(scip, presoldata->proprounds, cutoff) );
   
   /* evaluate propagation */
   if( !(*cutoff) )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         lbs[i] = SCIPvarGetLbLocal(vars[i]);
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
#ifdef DEBUG
         if( SCIPisGT(scip, lbs[i], SCIPvarGetLbGlobal(vars[i])) )
         {
            debugMessage(" -> <%s>[%g,%g] >= %g\n", SCIPvarGetName(vars[i]), 
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), lbs[i]);
         }
         if( SCIPisLT(scip, ubs[i], SCIPvarGetUbGlobal(vars[i])) )
         {
            debugMessage(" -> <%s>[%g,%g] <= %g\n", SCIPvarGetName(vars[i]), 
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), ubs[i]);
         }
#endif
      }
   }

   /* exit probing mode */
   CHECK_OKAY( SCIPendProbing(scip) );

   return SCIP_OKAY;
}




/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
DECL_PRESOLFREE(presolFreeProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
DECL_PRESOLINIT(presolInitProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   presoldata->sortedvars = NULL;
   presoldata->nsortedvars = 0;
   presoldata->nsortedbinvars = 0;
   presoldata->startidx = 0;
   presoldata->nfixings = 0;
   presoldata->naggregations = 0;
   presoldata->nimplications = 0;
   presoldata->nbdchgs = 0;

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
DECL_PRESOLEXIT(presolExitProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   if( presoldata->sortedvars != NULL )
   {
      SCIPfreeMemoryArray(scip, &presoldata->sortedvars);
      presoldata->nsortedvars = 0;
      presoldata->nsortedbinvars = 0;
   }
   assert(presoldata->sortedvars == NULL);
   assert(presoldata->nsortedvars == 0);
   assert(presoldata->nsortedbinvars == 0);

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
DECL_PRESOLINITPRE(presolInitpreProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   presoldata->called = FALSE;

   return SCIP_OKAY;
}


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DECL_PRESOLEXITPRE(presolExitpreProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* delete the vars array, if the maximal number of runs are exceeded */
   if( presoldata->maxruns >= 0 && SCIPgetNRuns(scip) >= presoldata->maxruns )
   {
      if( presoldata->sortedvars != NULL )
      {
         SCIPfreeMemoryArray(scip, &presoldata->sortedvars);
         presoldata->nsortedvars = 0;
         presoldata->nsortedbinvars = 0;
      }
      assert(presoldata->sortedvars == NULL);
      assert(presoldata->nsortedvars == 0);
      assert(presoldata->nsortedbinvars == 0);
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
DECL_PRESOLEXEC(presolExecProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;
   VAR** vars;
   Real* zerolbs;
   Real* zeroubs;
   Real* onelbs;
   Real* oneubs;
   int nvars;
   int nbinvars;
   int maxfixings;
   int maxuseless;
   int nuseless;
   int oldnfixedvars;
   int oldnaggrvars;
   int i;
   Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* check, if probing should be applied in the current run */
   if( presoldata->maxruns >= 0 && SCIPgetNRuns(scip) > presoldata->maxruns )
      return SCIP_OKAY;

   /* if no domains, sides, or coefficients changed since the last call, we don't need to probe */
   if( presoldata->called && nnewfixedvars == 0 && nnewaggrvars == 0 && nnewchgbds == 0 && nnewholes == 0
      && nnewchgcoefs == 0 && nnewchgsides == 0 )
      return SCIP_OKAY;

   /* get variable data */
   if( presoldata->sortedvars == NULL )
   {
      VAR** probvars;
      int* varsnlocks;

      CHECK_OKAY( SCIPgetVarsData(scip, &probvars, &nvars, &nbinvars, NULL, NULL, NULL) );
      if( nbinvars == 0 )
         return SCIP_OKAY;

      CHECK_OKAY( SCIPallocMemoryArray(scip, &presoldata->sortedvars, nvars) );
      presoldata->nsortedvars = nvars;
      presoldata->nsortedbinvars = nbinvars;
      vars = presoldata->sortedvars;

      /* sort the variables by number of rounding locks */
      CHECK_OKAY( SCIPallocBufferArray(scip, &varsnlocks, nbinvars) );
      for( i = 0; i < nbinvars; ++i )
      {
         int j;
         int nlocks;
         
         nlocks = SCIPvarGetNLocksDown(probvars[i]) + SCIPvarGetNLocksUp(probvars[i]);
         for( j = i; j > 0 && nlocks > varsnlocks[j-1]; --j )
         {
            vars[j] = vars[j-1];
            varsnlocks[j] = varsnlocks[j-1];
         }
         vars[j] = probvars[i];
         varsnlocks[j] = nlocks;
      }
      for( i = nbinvars; i < nvars; ++i )
         vars[i] = probvars[i];
      SCIPfreeBufferArray(scip, &varsnlocks);
   }
   else
   {
      vars = presoldata->sortedvars;
      nvars = presoldata->nsortedvars;
      nbinvars = presoldata->nsortedbinvars;
   }
   assert(vars != NULL);
   assert(nbinvars > 0);

   *result = SCIP_DIDNOTFIND;

   maxfixings = (presoldata->maxfixings > 0 ? presoldata->maxfixings : INT_MAX);
   maxuseless = (presoldata->maxuseless > 0 ? presoldata->maxuseless : INT_MAX);
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;

   /* get temporary memory for storing probing results */
   CHECK_OKAY( SCIPallocBufferArray(scip, &zerolbs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &zeroubs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &onelbs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &oneubs, nvars) );

   /* for each binary variable, probe fixing the variable to zero and one */
   cutoff = FALSE;
   nuseless = 0;
   for( i = presoldata->startidx; i < nbinvars && !cutoff && !SCIPpressedCtrlC(scip)
           && *nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars < maxfixings
           && nuseless < maxuseless; ++i )
   {
      Bool localcutoff;
      int j;

      /* display probing status */
      if( (i+1) % 1000 == 0 )
      {
         SCIPmessage(scip, SCIP_VERBLEVEL_FULL, "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n", 
            SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(Real)(i+1)/(Real)nbinvars,
            presoldata->nfixings, presoldata->naggregations, presoldata->nimplications, presoldata->nbdchgs);
      }

      /* ignore variables, that were fixed or aggregated in prior probings */
      if( !SCIPvarIsActive(vars[i]) )
         continue;

      nuseless++;

      /* apply probing for fixing the variable to zero */
      CHECK_OKAY( applyProbing(scip, presoldata, vars, nvars, i, FALSE, zerolbs, zeroubs, &localcutoff) );
      if( localcutoff )
      {
         Bool fixed;
         
         /* the variable can be fixed to TRUE */
         debugMessage("fixing probing variable <%s> to 1.0, nlocks=(%d/%d)\n",
            SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
         CHECK_OKAY( SCIPfixVar(scip, vars[i], 1.0, &cutoff, &fixed) );
         assert(fixed);
         (*nfixedvars)++;
         presoldata->nfixings++;
         nuseless = 0;
         if( cutoff )
            break;
         continue; /* don't try upwards direction, because the variable is already fixed */
      }

      /* apply probing for fixing the variable to one */
      CHECK_OKAY( applyProbing(scip, presoldata, vars, nvars, i, TRUE, onelbs, oneubs, &localcutoff) );
      if( localcutoff )
      {
         Bool fixed;
         
         /* the variable can be fixed to FALSE */
         debugMessage("fixing probing variable <%s> to 0.0, nlocks=(%d/%d)\n", 
            SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
         CHECK_OKAY( SCIPfixVar(scip, vars[i], 0.0, &cutoff, &fixed) );
         assert(fixed);
         (*nfixedvars)++;
         presoldata->nfixings++;
         nuseless = 0;
         if( cutoff )
            break;
         continue; /* don't analyze probing deductions, because the variable is already fixed */
      }

      /* analyze probing deductions */
      for( j = 0; j < nvars && !cutoff; ++j )
      {
         Real newlb;
         Real newub;
         Bool fixed;

         if( j == i )
            continue;

         /* new bounds of the variable is the union of the propagated bounds of the zero and one case */
         newlb = MIN(zerolbs[j], onelbs[j]);
         newub = MAX(zeroubs[j], oneubs[j]);

         /* check for fixed variables */
         if( SCIPisFeasEQ(scip, newlb, newub) )
         {
            Real fixval;

            /* in both probings, variable j is deduced to 0: fix variable to 0 */
            fixval = SCIPselectSimpleValue(newlb, newub, MAXDNOM);
            debugMessage("fixing variable <%s> to %g due to probing on <%s> with nlocks=(%d/%d)\n",
               SCIPvarGetName(vars[j]), fixval, 
               SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
            CHECK_OKAY( SCIPfixVar(scip, vars[j], fixval, &cutoff, &fixed) );
            if( fixed )
            {
               (*nfixedvars)++;
               presoldata->nfixings++;
               nuseless = 0;
            }
            continue;
         }

         if( j < nbinvars )
         {
            Bool redundant;
            Bool aggregated;

            assert(SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY);

            /* check for aggregations and implications on binary variables */
            if( zeroubs[j] < 0.5 && onelbs[j] > 0.5 )
            {
               /* variable j is always deduced to the same value as probing variable i:
                * both variables can be aggregated with x_i - x_j == 0
                */
               debugMessage("aggregating variables <%s> == <%s>, nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), 
                  SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               CHECK_OKAY( SCIPaggregateVars(scip, vars[i], vars[j], 1.0, -1.0, 0.0, &cutoff, &redundant, &aggregated) );
               if( aggregated )
               {
                  (*naggrvars)++;
                  presoldata->naggregations++;
                  nuseless = 0;
               }
            }
            else if( zerolbs[j] > 0.5 && oneubs[j] < 0.5 )
            {
               /* variable j is always deduced to the opposite value of probing variable i:
                * both variables can be aggregated with x_i + x_j == 1
                */
               debugMessage("aggregating variables <%s> == 1 - <%s>, nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]),
                  SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               CHECK_OKAY( SCIPaggregateVars(scip, vars[i], vars[j], 1.0, 1.0, 1.0, &cutoff, &redundant, &aggregated) );
               if( aggregated )
               {
                  (*naggrvars)++;
                  presoldata->naggregations++;
                  nuseless = 0;
               }
            }
            else if( zeroubs[j] < 0.5 )
            {
               /* insert implication: x_i == 0  =>  x_j == 0 */
               debugMessage("found implication <%s> == 0  =>  <%s> == 0\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff) );
               presoldata->nimplications++;
            }
            else if( zerolbs[j] > 0.5 )
            {
               /* insert implication: x_i == 0  =>  x_j == 1 */
               debugMessage("found implication <%s> == 0  =>  <%s> == 1\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff) );
               presoldata->nimplications++;
            }
            else if( oneubs[j] < 0.5 )
            {
               /* insert implication: x_i == 1  =>  x_j == 0 */
               debugMessage("found implication <%s> == 1  =>  <%s> == 0\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff) );
               presoldata->nimplications++;
            }
            else if( onelbs[j] > 0.5 )
            {
               /* insert implication: x_i == 1  =>  x_j == 1 */
               debugMessage("found implication <%s> == 1  =>  <%s> == 1\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff) );
               presoldata->nimplications++;
            }
         }
         else
         {
            Real oldlb;
            Real oldub;
            Bool tightened;

            assert(SCIPvarGetType(vars[j]) != SCIP_VARTYPE_BINARY);

            /* check for bound tightenings */
            oldlb = SCIPvarGetLbGlobal(vars[j]);
            oldub = SCIPvarGetUbGlobal(vars[j]);
            if( SCIPisLbBetter(scip, newlb, oldlb) )
            {
               /* in both probings, variable j is deduced to be at least newlb: tighten lower bound */
               debugMessage("tighten lower bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), oldlb, oldub, newlb,
                  SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               CHECK_OKAY( SCIPtightenVarLb(scip, vars[j], newlb, &cutoff, &tightened) );
               if( tightened )
               {
                  (*nchgbds)++;
                  presoldata->nbdchgs++;
                  nuseless = 0;
               }
            }
            if( SCIPisUbBetter(scip, newub, oldub) )
            {
               /* in both probings, variable j is deduced to be at most newub: tighten upper bound */
               debugMessage("tighten upper bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), oldlb, oldub, newub,
                  SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               CHECK_OKAY( SCIPtightenVarUb(scip, vars[j], newub, &cutoff, &tightened) );
               if( tightened )
               {
                  (*nchgbds)++;
                  presoldata->nbdchgs++;
                  nuseless = 0;
               }
            }

            /* check for implications (only store implications with bounds tightened at least by 0.5) */
            if( zeroubs[j] < newub - 0.5 )
            {
               /* insert implication: x_i == 0  =>  x_j <= zeroubs[j] */
               debugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] <= %g\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, zeroubs[j]);
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, zeroubs[j],
                     &cutoff) );
               presoldata->nimplications++;
            }
            if( zerolbs[j] > newlb + 0.5 )
            {
               /* insert implication: x_i == 0  =>  x_j >= zerolbs[j] */
               debugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] >= %g\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, zerolbs[j]);
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, zerolbs[j],
                     &cutoff) );
               presoldata->nimplications++;
            }
            if( oneubs[j] < newub - 0.5 )
            {
               /* insert implication: x_i == 1  =>  x_j <= oneubs[j] */
               debugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] <= %g\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, oneubs[j]);
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, oneubs[j],
                     &cutoff) );
               presoldata->nimplications++;
            }
            if( onelbs[j] > newlb + 0.5 )
            {
               /* insert implication: x_i == 1  =>  x_j >= onelbs[j] */
               debugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] >= %g\n", 
                  SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]), newlb, newub, onelbs[j]);
               CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, onelbs[j],
                     &cutoff) );
               presoldata->nimplications++;
            }
         }
      }
   }
   presoldata->startidx = i;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &oneubs);
   SCIPfreeBufferArray(scip, &onelbs);
   SCIPfreeBufferArray(scip, &zeroubs);
   SCIPfreeBufferArray(scip, &zerolbs);

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the probing presolver and includes it in SCIP */
RETCODE SCIPincludePresolProbing(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PRESOLDATA* presoldata;

   /* create probing presolver data */
   CHECK_OKAY( SCIPallocMemory(scip, &presoldata) );

   /* include presolver */
   CHECK_OKAY( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolFreeProbing, presolInitProbing, presolExitProbing, 
         presolInitpreProbing, presolExitpreProbing, presolExecProbing,
         presoldata) );

   /* add probing presolver parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
         "presolving/probing/maxruns", 
         "maximal number of runs, probing participates in (-1: no limit)",
         &presoldata->maxruns, DEFAULT_MAXRUNS, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "presolving/probing/proprounds", 
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &presoldata->proprounds, DEFAULT_PROPROUNDS, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "presolving/probing/maxfixings",
         "maximal number of fixings found, until probing is aborted (0: don't abort)",
         &presoldata->maxfixings, DEFAULT_MAXFIXINGS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "presolving/probing/maxuseless",
         "maximal number of successive useless probings, until probing is aborted (0: don't abort)",
         &presoldata->maxuseless, DEFAULT_MAXUSELESS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
