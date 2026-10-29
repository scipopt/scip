/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_probing.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  probing propagator
 * @author Tobias Achterberg
 * @author Matthias Miltenberger
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/prop_probing.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_prop.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_prop.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"


#define PROP_NAME               "probing"
#define PROP_DESC               "probing propagator on binary variables"
#define PROP_TIMING             SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           -100000 /**< propagation priority */
#define PROP_FREQ                    -1 /**< propagation frequency */
#define PROP_DELAY                 TRUE /**< should propagation method be delayed, if other propagators found
                                         *   reductions? */
#define PROP_PRESOL_PRIORITY    -100000 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */
#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */


/* @todo check for restricting the maximal number of implications that can be added by probing */

/* sorting of probing variables, two different variants are implemeneted */
/* #define VARIANT_B */


/*
 * Default parameter settings
 */

#define DEFAULT_MAXRUNS              1  /**< maximal number of runs, probing participates in (-1: no limit) */
#define DEFAULT_PROPROUNDS          -1  /**< maximal number of propagation rounds in probing subproblems */
#define DEFAULT_MAXFIXINGS          25  /**< maximal number of fixings found, until probing is interrupted
                                         *   (0: don't interrupt) */
#define DEFAULT_MAXUSELESS        1000  /**< maximal number of successive probings without fixings,
                                         *   until probing is aborted (0: don't abort) */
#define DEFAULT_MAXTOTALUSELESS     50  /**< maximal number of successive probings without fixings, bound changes,
                                         *   and implications, until probing is aborted (0: don't abort) */
#define DEFAULT_MAXSUMUSELESS        0  /**< maximal number of probings without fixings, until probing is aborted
                                         *   (0: don't abort) */
#define DEFAULT_MAXDEPTH            -1  /**< maximal depth until propagation is executed(-1: no limit) */
#define DEFAULT_RANDSEED            59  /**< random initial seed */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_VAR**            sortedvars;         /**< problem variables sorted by number of rounding locks, used in presolving */
   int*                  nprobed;            /**< array of numbers how often we already probed on each variables */
   int                   noldtotalvars;      /**< number of total variables in problem */
   int                   nsortedvars;        /**< number of problem variables, used in presolving */
   int                   nsortedbinvars;     /**< number of binary problem variables, used in presolving */
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
   int                   startidx;           /**< starting variable index of next call, used in presolving */
   int                   lastsortstartidx;   /**< last starting variable index where the variables have been sorted, used in presolving */
   int                   nfixings;           /**< total number of fixings found in probing */
   int                   naggregations;      /**< total number of aggregations found in probing */
   int                   nimplications;      /**< total number of implications found in probing */
   int                   nbdchgs;            /**< total number of bound changes found in probing */
   int                   nuseless;           /**< current number of successive useless probings */
   int                   ntotaluseless;      /**< current number of successive totally useless probings */
   int                   nsumuseless;        /**< current number of useless probings */
   int                   maxdepth;           /**< maximal depth until propagation is executed */
   SCIP_Longint          lastnode;           /**< last node where probing was applied, or -1 for presolving, and -2 for not applied yet */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
};


/** sparse bookkeeping of probing in presolving: only the probed variables whose bounds a probing changed */
struct ProbingSparse
{
   int*                  posof;              /**< by variable index: position in the probed variables, or -1 */
   int*                  depidx;             /**< sorted indices of active variables, with ... */
   int*                  deppos;             /**< ... the position of a probed variable whose bounds derive from them */
   int*                  eqpos;              /**< positions of the unfixed probed variables with equal bounds */
   int*                  implmark[2];        /**< by position, per direction: stamp when implications changed it */
   int*                  propmark[2];        /**< by position, per direction: stamp when propagation changed it */
   int*                  proppos[2];         /**< per direction: positions that the last probing changed */
   int*                  setmark;            /**< by position: stamp of the last analysis that included it */
   int*                  setpos;             /**< positions of the variables of the current analysis */
   int                   implstamp[2];       /**< per direction: stamp after implications and cliques */
   int                   propstamp[2];       /**< per direction: stamp after propagation */
   int                   nproppos[2];        /**< per direction: number of positions in proppos */
   int                   nposof;             /**< length of posof */
   int                   ndeps;              /**< number of entries in depidx and deppos */
   int                   depssize;           /**< size of depidx and deppos */
   int                   neqpos;             /**< number of entries in eqpos */
   int                   stamp;              /**< last stamp given out */
   SCIP_Bool             dirty;              /**< do eqpos and the derived bounds need a new scan? */
};
typedef struct ProbingSparse PROBINGSPARSE;


/*
 * Local methods
 */
/** initializes the propagator data */
static
SCIP_RETCODE initPropdata(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   propdata->sortedvars = NULL;
   propdata->nprobed = NULL;
   propdata->noldtotalvars = 0;
   propdata->nsortedvars = 0;
   propdata->nsortedbinvars = 0;
   propdata->startidx = 0;
   propdata->lastsortstartidx = -1;
   propdata->nfixings = 0;
   propdata->naggregations = 0;
   propdata->nimplications = 0;
   propdata->nbdchgs = 0;
   propdata->nuseless = 0;
   propdata->ntotaluseless = 0;
   propdata->nsumuseless = 0;
   propdata->lastnode = -2;
   propdata->randnumgen = NULL;

   return SCIP_OKAY;
}

/** frees the sorted vars array */
static
SCIP_RETCODE freeSortedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   if( propdata->sortedvars != NULL )
   {
      int i;

      /* release variables */
      for( i = 0; i < propdata->nsortedvars; ++i )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &propdata->sortedvars[i]) );
      }
      SCIPfreeMemoryArray(scip, &propdata->sortedvars);
      propdata->nsortedvars = 0;
      propdata->nsortedbinvars = 0;
   }

   SCIPfreeMemoryArrayNull(scip, &propdata->nprobed);
   propdata->noldtotalvars = 0;

   return SCIP_OKAY;
}

/** sorts the binary variables starting with the given index by rounding locks and implications */
static
SCIP_RETCODE sortVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR**            vars,               /**< problem variables to be sorted */
   int                   nvars,              /**< number of problem variables to be sorted */
   int                   firstidx            /**< first index that should be subject to sorting */
   )
{
   SCIP_VAR** sortedvars;
   int nsortedvars;
   SCIP_Real* scores;
   int i;
   int minnprobings;
   SCIP_Real maxscore;
   int nlocksdown;
   int nlocksup;
   int nimplzero;
   int nimplone;
   int nclqzero;
   int nclqone;

   assert(propdata != NULL);
   assert(propdata->nprobed != NULL);

   assert(vars != NULL || nvars == 0);

   nsortedvars = nvars - firstidx;
   if( nsortedvars <= 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   sortedvars = &(vars[firstidx]);

   SCIPdebugMsg(scip, "resorting probing variables %d to %d\n", firstidx, nvars-1);

   /* sort the variables by number of rounding locks and implications */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, nsortedvars) );

   maxscore = -1.0;
   minnprobings = INT_MAX;

   /* determine maximal possible score and minimal number of probings over all variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real tmp;

      var = vars[i];

      assert(SCIPvarIsBinary(var));
      assert(propdata->noldtotalvars > SCIPvarGetIndex(var));
      assert(propdata->nprobed[SCIPvarGetIndex(var)] >= 0);

      if( SCIPvarIsActive(var) )
      {
         nlocksdown = SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL);
         nlocksup = SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL);
         nimplzero = SCIPvarGetNImpls(var, FALSE);
         nimplone = SCIPvarGetNImpls(var, TRUE);
         nclqzero = SCIPvarGetNCliques(var, FALSE);
         nclqone = SCIPvarGetNCliques(var, TRUE);

#ifndef VARIANT_B
         tmp = -MAX(nlocksdown, nlocksup)
            + 10.0 * MIN(nimplzero, nimplone)
            + 100.0 * MIN(nclqzero, nclqone);
#else
         tmp = - ABS(nlocksdown - nlocksup)
            + MIN(nlocksdown, nlocksup)
            + 500.0 * nimplzero + 50.0 * nimplone
            + 50000.0 * nclqzero + 5000.0 * nclqone;
#endif

         if( tmp > maxscore )
            maxscore = tmp;
         if( propdata->nprobed[SCIPvarGetIndex(var)] < minnprobings )
            minnprobings = propdata->nprobed[SCIPvarGetIndex(var)];
      }
   }

   /* correct number of probings on each variable by minimal number of probings */
   if( minnprobings > 0 )
   {
      for( i = 0; i < nvars; ++i )
      {
         SCIP_VAR* var;

         var = vars[i];

         if( SCIPvarIsActive(var) )
            propdata->nprobed[SCIPvarGetIndex(var)] -= minnprobings;
      }
   }

   for( i = 0; i < nsortedvars; ++i )
   {
      SCIP_VAR* var;
      var = sortedvars[i];

      assert(SCIPvarIsBinary(var));

      /* prefer variables that we did not already probe on */
      if( SCIPvarIsActive(var) )
      {
         SCIP_Real randomoffset;
         nlocksdown = SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL);
         nlocksup = SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL);
         nimplzero = SCIPvarGetNImpls(var, FALSE);
         nimplone = SCIPvarGetNImpls(var, TRUE);
         nclqzero = SCIPvarGetNCliques(var, FALSE);
         nclqone = SCIPvarGetNCliques(var, TRUE);

         assert(propdata->noldtotalvars > SCIPvarGetIndex(var));
         assert(propdata->nprobed[SCIPvarGetIndex(var)] >= 0);

         /* use a random offset to break possible ties arbitrarily */
         randomoffset = SCIPrandomGetReal(propdata->randnumgen, 0.0, 0.5);

#ifndef VARIANT_B
         scores[i] = -maxscore * propdata->nprobed[SCIPvarGetIndex(var)]
            - MAX(nlocksdown, nlocksup)
            + 10.0 * MIN(nimplzero, nimplone)
            + 100.0 * MIN(nclqzero, nclqone)  /*lint !e790*/
            - randomoffset; /* to break ties randomly */
#else
         scores[i] = -maxscore * propdata->nprobed[SCIPvarGetIndex(var)]
         - ABS(nlocksdown - nlocksup)
         + MIN(nlocksdown, nlocksup)
         + 500.0 * nimplzero + 50.0 * nimplone  /*lint !e790*/
         + 50000.0 * nclqzero + 5000.0 * nclqone  /*lint !e790*/
         - randomoffset; /* to break ties randomly */
#endif
      }
      else
         scores[i] = -SCIPinfinity(scip);
   }

   SCIPsortDownRealPtr(scores, (void**) sortedvars, nsortedvars);

   SCIPfreeBufferArray(scip, &scores);

   return SCIP_OKAY;
}

/** creates the sparse bookkeeping for probing on the given variables */
static
SCIP_RETCODE sparseCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping */
   SCIP_VAR**            vars,               /**< probed variables */
   int                   nvars               /**< number of probed variables */
   )
{
   int d;
   int p;

   assert(sparse != NULL);
   assert(vars != NULL || nvars == 0);

   sparse->nposof = SCIPgetNTotalVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &sparse->posof, sparse->nposof) );
   for( p = 0; p < sparse->nposof; ++p )
      sparse->posof[p] = -1;
   for( p = 0; p < nvars; ++p )
   {
      assert(SCIPvarGetIndex(vars[p]) < sparse->nposof);
      sparse->posof[SCIPvarGetIndex(vars[p])] = p;
   }

   for( d = 0; d < 2; ++d )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &sparse->implmark[d], nvars) ); /*lint !e866*/
      SCIP_CALL( SCIPallocClearBufferArray(scip, &sparse->propmark[d], nvars) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &sparse->proppos[d], nvars) );
      sparse->implstamp[d] = 0;
      sparse->propstamp[d] = 0;
      sparse->nproppos[d] = 0;
   }
   SCIP_CALL( SCIPallocClearBufferArray(scip, &sparse->setmark, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sparse->setpos, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sparse->eqpos, nvars) );

   sparse->depssize = MAX(nvars, 1);
   SCIP_CALL( SCIPallocBufferArray(scip, &sparse->depidx, sparse->depssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sparse->deppos, sparse->depssize) );

   sparse->ndeps = 0;
   sparse->neqpos = 0;
   sparse->stamp = 0;
   sparse->dirty = TRUE;

   return SCIP_OKAY;
}

/** frees the sparse bookkeeping */
static
void sparseFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBINGSPARSE*        sparse              /**< sparse bookkeeping */
   )
{
   int d;

   assert(sparse != NULL);

   SCIPfreeBufferArray(scip, &sparse->deppos);
   SCIPfreeBufferArray(scip, &sparse->depidx);
   SCIPfreeBufferArray(scip, &sparse->eqpos);
   SCIPfreeBufferArray(scip, &sparse->setpos);
   SCIPfreeBufferArray(scip, &sparse->setmark);
   for( d = 1; d >= 0; --d )
   {
      SCIPfreeBufferArray(scip, &sparse->proppos[d]);
      SCIPfreeBufferArray(scip, &sparse->propmark[d]);
      SCIPfreeBufferArray(scip, &sparse->implmark[d]);
   }
   SCIPfreeBufferArray(scip, &sparse->posof);
}

/** records that the bounds of the probed variable at pos derive from the active variables that var represents */
static
SCIP_RETCODE sparseAddDeps(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping */
   SCIP_VAR*             var,                /**< inactive variable */
   int                   pos                 /**< position of the probed variable */
   )
{
   SCIP_VAR** aggrvars;
   int naggrvars;
   int a;

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      aggrvars = SCIPvarGetMultaggrVars(var);
      naggrvars = SCIPvarGetMultaggrNVars(var);
   }
   else
   {
      aggrvars = &var;
      naggrvars = 1;
   }

   for( a = 0; a < naggrvars; ++a )
   {
      SCIP_VAR* probvar;

      probvar = SCIPvarGetProbvar(aggrvars[a]);
      if( probvar == NULL )
         continue;

      switch( SCIPvarGetStatus(probvar) )
      {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         if( sparse->ndeps == sparse->depssize )
         {
            sparse->depssize = SCIPcalcMemGrowSize(scip, sparse->ndeps + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &sparse->depidx, sparse->depssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &sparse->deppos, sparse->depssize) );
         }
         sparse->depidx[sparse->ndeps] = SCIPvarGetIndex(probvar);
         sparse->deppos[sparse->ndeps] = pos;
         ++sparse->ndeps;
         break;
      case SCIP_VARSTATUS_MULTAGGR:
         if( probvar != var )
         {
            SCIP_CALL( sparseAddDeps(scip, sparse, probvar, pos) );
         }
         break;
      case SCIP_VARSTATUS_FIXED:
         /* the bounds that derive from a fixed variable do not change */
         break;
      case SCIP_VARSTATUS_ORIGINAL:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_NEGATED:
      default:
         SCIPerrorMessage("invalid status of an active representative\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** scans the probed variables for unfixed ones with equal bounds (SCIPanalyzeDeductionsProbing() fixes them) and
 *  inactive ones (their bounds derive from active variables)
 */
static
SCIP_RETCODE sparseScan(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping */
   SCIP_VAR**            vars,               /**< probed variables */
   int                   nvars               /**< number of probed variables */
   )
{
   int p;

   assert(sparse != NULL);

   sparse->neqpos = 0;
   sparse->ndeps = 0;

   for( p = 0; p < nvars; ++p )
   {
      SCIP_VARSTATUS status;

      status = SCIPvarGetStatus(vars[p]);
      if( status == SCIP_VARSTATUS_FIXED )
         continue;

      if( SCIPisEQ(scip, SCIPvarGetLbLocal(vars[p]), SCIPvarGetUbLocal(vars[p])) )
         sparse->eqpos[sparse->neqpos++] = p;

      if( status != SCIP_VARSTATUS_LOOSE && status != SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( sparseAddDeps(scip, sparse, vars[p], p) );
      }
   }

   SCIPsortIntInt(sparse->depidx, sparse->deppos, sparse->ndeps);
   sparse->dirty = FALSE;

   return SCIP_OKAY;
}

/** marks and stores the positions of the probed variables that the bound changes of the probing node changed */
static
void sparseMark(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping */
   int*                  mark,               /**< by position: marks */
   int                   stamp,              /**< stamp to mark with */
   int*                  pos,                /**< array to store the positions of the marked variables */
   int*                  npos                /**< pointer to store the number of marked variables */
   )
{
   SCIP_DOMCHG* domchg;
   int nboundchgs;
   int b;

   assert(sparse != NULL);
   assert(!sparse->dirty);

   domchg = SCIPnodeGetDomchg(SCIPgetCurrentNode(scip));
   nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
   *npos = 0;

   for( b = 0; b < nboundchgs; ++b )
   {
      int idx;
      int first;
      int d;

      idx = SCIPvarGetIndex(SCIPboundchgGetVar(SCIPdomchgGetBoundchg(domchg, b)));
      if( idx >= sparse->nposof )
         continue;

      d = sparse->posof[idx];
      if( d >= 0 && mark[d] != stamp )
      {
         mark[d] = stamp;
         pos[(*npos)++] = d;
      }

      /* the inactive probed variables whose bounds derive from this variable */
      if( sparse->ndeps > 0 && SCIPsortedvecFindInt(sparse->depidx, idx, sparse->ndeps, &first) )
      {
         while( first > 0 && sparse->depidx[first - 1] == idx )
            --first;

         for( ; first < sparse->ndeps && sparse->depidx[first] == idx; ++first )
         {
            d = sparse->deppos[first];
            if( mark[d] != stamp )
            {
               mark[d] = stamp;
               pos[(*npos)++] = d;
            }
         }
      }
   }
}

/** SCIPanalyzeDeductionsProbing() for probing a binary variable to 0 and 1, only on the variables it can deduce
 *  something for: those that probing changed in either direction and the unfixed ones with equal bounds; in the order
 *  of the probed variables and with the bounds from before probing where a direction did not change a variable, it
 *  makes the same deductions as on all variables
 */
static
SCIP_RETCODE analyzeDeductionsSparse(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping of both directions */
   SCIP_VAR**            vars,               /**< probed variables */
   int                   nvars,              /**< number of probed variables */
   int                   probingpos,         /**< position of the probing variable */
   SCIP_Real*            zeroimpllbs,        /**< lower bounds after implications and cliques when probing to 0 */
   SCIP_Real*            zeroimplubs,        /**< upper bounds after implications and cliques when probing to 0 */
   SCIP_Real*            zeroproplbs,        /**< lower bounds after propagation when probing to 0 */
   SCIP_Real*            zeropropubs,        /**< upper bounds after propagation when probing to 0 */
   SCIP_Real*            oneimpllbs,         /**< lower bounds after implications and cliques when probing to 1 */
   SCIP_Real*            oneimplubs,         /**< upper bounds after implications and cliques when probing to 1 */
   SCIP_Real*            oneproplbs,         /**< lower bounds after propagation when probing to 1 */
   SCIP_Real*            onepropubs,         /**< upper bounds after propagation when probing to 1 */
   int*                  nfixedvars,         /**< pointer to counter of deduced variable fixations */
   int*                  naggrvars,          /**< pointer to counter of deduced variable aggregations */
   int*                  nimplications,      /**< pointer to counter of deduced implications */
   int*                  nchgbds,            /**< pointer to counter of deduced bound tightenings */
   SCIP_Bool*            cutoff              /**< buffer to store whether a cutoff is detected */
   )
{
   /* the 4 bound arrays of a probing direction */
   enum ProbingBounds
   {
      IMPL_LB = 0,                           /**< lower bounds after implications and cliques */
      IMPL_UB = 1,                           /**< upper bounds after implications and cliques */
      PROP_LB = 2,                           /**< lower bounds after propagation */
      PROP_UB = 3                            /**< upper bounds after propagation */
   };
   SCIP_VAR** setvars;
   SCIP_Real* src[2][4];
   SCIP_Real* dst[2][4];
   int* setpos;
   int nset;
   int stamp;
   int a;
   int d;
   int k;

   assert(sparse != NULL);
   assert(!sparse->dirty);

   /* the variables that probing changed in either direction and the unfixed ones with equal bounds */
   setpos = sparse->setpos;
   stamp = ++sparse->stamp;
   nset = 0;
   for( d = 0; d < 3; ++d )
   {
      int* pos = (d < 2 ? sparse->proppos[d] : sparse->eqpos);
      int npos = (d < 2 ? sparse->nproppos[d] : sparse->neqpos);

      for( k = 0; k < npos; ++k )
      {
         if( sparse->setmark[pos[k]] != stamp )
         {
            sparse->setmark[pos[k]] = stamp;
            setpos[nset++] = pos[k];
         }
      }
   }

   /* in the order of the probed variables: for a large set, scanning the marks of all variables is cheaper than sorting
    * (a probing can change thousands of variables, and sorting them cost about as much as the dense comparison)
    */
   if( nset > nvars / 16 )
   {
      int p;

      nset = 0;
      for( p = 0; p < nvars; ++p )
      {
         if( sparse->setmark[p] == stamp )
            setpos[nset++] = p;
      }
   }
   else
      SCIPsortInt(setpos, nset);

   /* per direction (0: probing to 0, 1: to 1) */
   src[0][IMPL_LB] = zeroimpllbs;
   src[0][IMPL_UB] = zeroimplubs;
   src[0][PROP_LB] = zeroproplbs;
   src[0][PROP_UB] = zeropropubs;
   src[1][IMPL_LB] = oneimpllbs;
   src[1][IMPL_UB] = oneimplubs;
   src[1][PROP_LB] = oneproplbs;
   src[1][PROP_UB] = onepropubs;
   SCIP_CALL( SCIPallocBufferArray(scip, &setvars, MAX(nset, 1)) );
   for( d = 0; d < 2; ++d )
   {
      for( a = 0; a < 4; ++a )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &dst[d][a], MAX(nset, 1)) ); /*lint !e866*/
      }
   }

   /* their bounds in the order of the probed variables: as probing stored them for a direction that changed the
    * variable, the bounds from before probing otherwise
    */
   for( k = 0; k < nset; ++k )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      int p;

      p = setpos[k];
      setvars[k] = vars[p];
      lb = SCIPvarGetLbLocal(vars[p]);
      ub = SCIPvarGetUbLocal(vars[p]);

      for( d = 0; d < 2; ++d )
      {
         SCIP_Bool implchanged;
         SCIP_Bool propchanged;

         implchanged = (sparse->implmark[d][p] == sparse->implstamp[d]);
         propchanged = (sparse->propmark[d][p] == sparse->propstamp[d]);
         dst[d][IMPL_LB][k] = (implchanged ? src[d][IMPL_LB][p] : lb);
         dst[d][IMPL_UB][k] = (implchanged ? src[d][IMPL_UB][p] : ub);
         dst[d][PROP_LB][k] = (propchanged ? src[d][PROP_LB][p] : lb);
         dst[d][PROP_UB][k] = (propchanged ? src[d][PROP_UB][p] : ub);
      }
   }

#ifndef NDEBUG
   /* debug mode: the set has the bounds of the dense arrays, and on all other variables the analysis does nothing */
   {
      int p;

      for( k = 0, p = 0; p < nvars; ++p )
      {
         if( k < nset && setpos[k] == p )
         {
            for( d = 0; d < 2; ++d )
            {
               for( a = 0; a < 4; ++a )
                  assert(dst[d][a][k] == src[d][a][p]); /*lint !e777*/
            }
            ++k;
         }
         else
         {
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbLocal(vars[p]);
            ub = SCIPvarGetUbLocal(vars[p]);
            for( d = 0; d < 2; ++d )
            {
               assert(src[d][IMPL_LB][p] == lb && src[d][IMPL_UB][p] == ub); /*lint !e777*/
               assert(src[d][PROP_LB][p] == lb && src[d][PROP_UB][p] == ub); /*lint !e777*/
            }
            assert(SCIPvarGetStatus(vars[p]) == SCIP_VARSTATUS_FIXED || !SCIPisEQ(scip, lb, ub));
         }
      }
   }
#endif

   SCIP_CALL( SCIPanalyzeDeductionsProbing(scip, vars[probingpos], 0.0, 1.0, nset, setvars,
         dst[0][IMPL_LB], dst[0][IMPL_UB], dst[0][PROP_LB], dst[0][PROP_UB],
         dst[1][IMPL_LB], dst[1][IMPL_UB], dst[1][PROP_LB], dst[1][PROP_UB],
         nfixedvars, naggrvars, nimplications, nchgbds, cutoff) );

   for( d = 1; d >= 0; --d )
   {
      for( a = 3; a >= 0; --a )
      {
         SCIPfreeBufferArray(scip, &dst[d][a]);
      }
   }
   SCIPfreeBufferArray(scip, &setvars);

   return SCIP_OKAY;
}

/* forward declaration: SCIPapplyProbingVar() with sparse bookkeeping */
static
SCIP_RETCODE applyProbingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   probingpos,         /**< variable number to apply probing on */
   SCIP_BOUNDTYPE        boundtype,          /**< which bound should be changed */
   SCIP_Real             bound,              /**< which bound should be set */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping, or NULL to store the bounds of all variables */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   );

/** the main probing loop */
static
SCIP_RETCODE applyProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   nbinvars,           /**< number of binary variables */
   int*                  startidx,           /**< pointer to store starting variable index of next call */
   int*                  nfixedvars,         /**< pointer to store number of fixed variables */
   int*                  naggrvars,          /**< pointer to store number of aggregated variables */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   int                   oldnfixedvars,      /**< number of previously fixed variables */
   int                   oldnaggrvars,       /**< number of previously aggregated variables */
   SCIP_Bool*            delay,              /**< pointer to store whether propagator should be delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether cutoff occurred */
   )
{
   SCIP_Real* zeroimpllbs;
   SCIP_Real* zeroimplubs;
   SCIP_Real* zeroproplbs;
   SCIP_Real* zeropropubs;
   SCIP_Real* oneimpllbs;
   SCIP_Real* oneimplubs;
   SCIP_Real* oneproplbs;
   SCIP_Real* onepropubs;
   PROBINGSPARSE sparsedata;
   PROBINGSPARSE* sparse;
   int localnfixedvars;
   int localnaggrvars;
   int localnchgbds;
   int localnimplications;
   int maxfixings;
   int maxuseless;
   int maxtotaluseless;
   int maxsumuseless;
   int i;
   int oldstartidx;
   SCIP_Bool aborted;
   SCIP_Bool looped;

   assert(vars != NULL);
   assert(nbinvars > 0);

   maxfixings = (propdata->maxfixings > 0 ? propdata->maxfixings : INT_MAX);
   maxuseless = (propdata->maxuseless > 0 ? propdata->maxuseless : INT_MAX);
   maxtotaluseless = (propdata->maxtotaluseless > 0 ? propdata->maxtotaluseless : INT_MAX);
   maxsumuseless = (propdata->maxsumuseless > 0 ? propdata->maxsumuseless : INT_MAX);
   aborted = FALSE;
   looped = FALSE;
   oldstartidx = *startidx;
   i = *startidx;

   /* get temporary memory for storing probing results */
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimpllbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimplubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroproplbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeropropubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimpllbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimplubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneproplbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &onepropubs, nvars) );

   /* in presolving, record per probing only the variables whose bounds it changed */
   sparse = NULL;
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      sparse = &sparsedata;
      SCIP_CALL( sparseCreate(scip, sparse, vars, nvars) );
   }

   /* for each binary variable, probe fixing the variable to zero and one */
   *delay = FALSE;
   *cutoff = FALSE;
   do
   {
      for( ; i < nbinvars && !(*cutoff); ++i )
      {
         SCIP_Bool localcutoff;
         SCIP_Bool probingzero;
         SCIP_Bool probingone;

         /* check whether probing should be aborted */
         if( propdata->nuseless >= maxuseless || propdata->ntotaluseless >= maxtotaluseless || propdata->nsumuseless >= maxsumuseless || SCIPisStopped(scip) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n",
               SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(SCIP_Real)(i+1)/(SCIP_Real)nbinvars,
               propdata->nfixings, propdata->naggregations, propdata->nimplications, propdata->nbdchgs);

            aborted = TRUE;

            if( propdata->nuseless >= maxuseless )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "   (%.1fs) probing aborted: %d/%d successive useless probings\n", SCIPgetSolvingTime(scip),
                  propdata->nuseless, maxuseless);
            }
            else if( propdata->ntotaluseless >= maxtotaluseless )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "   (%.1fs) probing aborted: %d/%d successive totally useless probings\n", SCIPgetSolvingTime(scip),
                  propdata->ntotaluseless, maxtotaluseless);
            }
            else if( propdata->nsumuseless >= maxsumuseless )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "   (%.1fs) probing aborted: %d/%d useless probings in total\n", SCIPgetSolvingTime(scip),
                  propdata->nsumuseless, maxsumuseless);
            }
            else
            {
               assert(SCIPisStopped(scip));
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "   (%.1fs) probing aborted: solving stopped\n", SCIPgetSolvingTime(scip));
            }
            break;
         }

         /* check if we already fixed enough variables for this round, or probed on all variables */
         if( *nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars >= maxfixings || (looped && oldstartidx == i) )
         {
            if( *nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars > 0 )
               *delay = TRUE;
            else
               aborted = TRUE;
            break;
         }

         /* display probing status */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && (i+1) % 100 == 0 )
         {
            SCIP_VERBLEVEL verblevel;

            verblevel = ((i+1) % 1000 == 0 ? SCIP_VERBLEVEL_HIGH : SCIP_VERBLEVEL_FULL);
            SCIPverbMessage(scip, verblevel, NULL,
               "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n",
               SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(SCIP_Real)(i+1)/(SCIP_Real)nbinvars,
               propdata->nfixings, propdata->naggregations, propdata->nimplications, propdata->nbdchgs);
         }

         /* ignore variables, that were fixed, aggregated, or deleted in prior probings */
         if( !SCIPvarIsActive(vars[i]) || SCIPvarIsDeleted(vars[i])
            || SCIPvarGetLbLocal(vars[i]) > 0.5 || SCIPvarGetUbLocal(vars[i]) < 0.5 )
            continue;

         if( propdata->nuseless > 0 )
            propdata->nsumuseless++;
         else
            propdata->nsumuseless = MAX(propdata->nsumuseless-1, 0);
         propdata->nuseless++;
         propdata->ntotaluseless++;

         /* bounds or statuses changed since the last scan: rescan for the sparse bookkeeping */
         if( sparse != NULL && sparse->dirty )
         {
            SCIP_CALL( sparseScan(scip, sparse, vars, nvars) );
         }

         /* determine whether one probing should happen */
         probingone = TRUE;
         if( SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL) == 0 )
            probingone = FALSE;

         if( probingone )
         {
            /* apply probing for fixing the variable to one */
            SCIP_CALL( applyProbingVar(scip, vars, nvars, i, SCIP_BOUNDTYPE_LOWER, 1.0, propdata->proprounds,
                  oneimpllbs, oneimplubs, oneproplbs, onepropubs, sparse, &localcutoff) );

            if( localcutoff )
            {
               SCIP_Bool fixed;

               if( sparse != NULL )
                  sparse->dirty = TRUE;

               if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
               {
                  /* the variable can be fixed to FALSE */
                  SCIP_CALL( SCIPfixVar(scip, vars[i], 0.0, cutoff, &fixed) );
                  assert(fixed);
               }
               else
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, vars[i], 0.0, TRUE, cutoff, &fixed) );
               }

               if( fixed )
               {
                  SCIPdebugMsg(scip, "fixed probing variable <%s> to 0.0, nlocks=(%d/%d)\n",
                     SCIPvarGetName(vars[i]), SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL),
                     SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL));
                  (*nfixedvars)++;
                  propdata->nfixings++;
                  propdata->nuseless = 0;
                  propdata->ntotaluseless = 0;
               }
               else if( *cutoff )
               {
                  SCIPdebugMsg(scip, "tightening upper bound of probing variable <%s> to 0.0 led to a cutoff\n",
                     SCIPvarGetName(vars[i]));
               }
               continue; /* don't try downwards direction, because the variable is already fixed */
            }

            /* ignore variables, that were fixed, aggregated, or deleted in prior probings
             * (propagators in one-probe might have found global fixings but did not trigger the localcutoff)
             */
            if( !SCIPvarIsActive(vars[i]) || SCIPvarIsDeleted(vars[i])
               || SCIPvarGetLbLocal(vars[i]) > 0.5 || SCIPvarGetUbLocal(vars[i]) < 0.5 )
            {
               if( sparse != NULL )
                  sparse->dirty = TRUE;
               continue;
            }
         }

         /* determine whether zero probing should happen */
         probingzero = TRUE;
         if( SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL) == 0 )
            probingzero = FALSE;

         if( probingzero )
         {
            /* apply probing for fixing the variable to zero */
            SCIP_CALL( applyProbingVar(scip, vars, nvars, i, SCIP_BOUNDTYPE_UPPER, 0.0, propdata->proprounds,
                  zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs, sparse, &localcutoff) );

            if( localcutoff )
            {
               SCIP_Bool fixed;

               if( sparse != NULL )
                  sparse->dirty = TRUE;

               if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
               {
                  /* the variable can be fixed to TRUE */
                  SCIP_CALL( SCIPfixVar(scip, vars[i], 1.0, cutoff, &fixed) );
               }
               else
               {
                  SCIP_CALL( SCIPtightenVarLb(scip, vars[i], 1.0, TRUE, cutoff, &fixed) );
               }

               if( fixed )
               {
                  SCIPdebugMsg(scip, "fixed probing variable <%s> to 1.0, nlocks=(%d/%d)\n",
                     SCIPvarGetName(vars[i]), SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL),
                     SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL));
                  (*nfixedvars)++;
                  propdata->nfixings++;
                  propdata->nuseless = 0;
                  propdata->ntotaluseless = 0;
               }
               else if( *cutoff )
               {
                  SCIPdebugMsg(scip, "tightening lower bound of probing variable <%s> to 1.0 led to a cutoff\n",
                     SCIPvarGetName(vars[i]));
               }
               continue; /* don't analyze probing deductions, because the variable is already fixed */
            }
         }

         /* not have to check deductions if only one probing direction has been checked */
         if( !probingzero || !probingone )
            continue;

         assert(propdata->noldtotalvars > SCIPvarGetIndex(vars[i]));

         /* count number of probings on each variable */
         propdata->nprobed[SCIPvarGetIndex(vars[i])] += 1;

         /* analyze probing deductions */
         localnfixedvars    = 0;
         localnaggrvars     = 0;
         localnimplications = 0;
         localnchgbds       = 0;
         if( sparse != NULL )
         {
            SCIP_CALL( analyzeDeductionsSparse(scip, sparse, vars, nvars, i, zeroimpllbs, zeroimplubs, zeroproplbs,
                  zeropropubs, oneimpllbs, oneimplubs, oneproplbs, onepropubs, &localnfixedvars, &localnaggrvars,
                  &localnimplications, &localnchgbds, cutoff) );

            /* deductions changed bounds or statuses */
            if( localnfixedvars > 0 || localnaggrvars > 0 || localnchgbds > 0 || *cutoff )
               sparse->dirty = TRUE;
         }
         else
         {
            SCIP_CALL( SCIPanalyzeDeductionsProbing(scip, vars[i], 0.0, 1.0,
                  nvars, vars, zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs, oneimpllbs, oneimplubs, oneproplbs, onepropubs,
                  &localnfixedvars, &localnaggrvars, &localnimplications, &localnchgbds, cutoff) );
         }

         *nfixedvars += localnfixedvars;
         *naggrvars  += localnaggrvars;
         *nchgbds    += localnchgbds;
         propdata->nfixings      += localnfixedvars;
         propdata->naggregations += localnaggrvars;
         propdata->nbdchgs       += localnchgbds;
         propdata->nimplications += localnimplications;

         if( localnfixedvars > 0 || localnaggrvars > 0 )
         {
            SCIPdebugMsg(scip, "probing on <%s> led to %d fixed and %d aggregated variables\n", SCIPvarGetName(vars[i]),
               localnfixedvars, localnaggrvars);
            propdata->nuseless = 0;
            propdata->ntotaluseless = 0;
         }
         if( localnimplications > 0 || localnchgbds > 0 )
            propdata->ntotaluseless = 0;
      }

      looped = TRUE;

      /* check if we reached the end of all binary variables but did not stop, so we start from the beginning */
      if( i == nbinvars && !(*cutoff) && !(*delay) && !aborted )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) probing cycle finished: starting next cycle\n", SCIPgetSolvingTime(scip));
         i = 0;

         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
         {
            SCIP_VAR** newvars;
            int nnewvars;
            int nnewbinvars;
            int nnewintvars;
            int nnewimplvars;
            int lastidx;
            int v;

            assert(vars == propdata->sortedvars);
            assert(nbinvars == propdata->nsortedbinvars);

            /* release old variables and free memory */
            for( v = propdata->nsortedvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &propdata->sortedvars[v]) );
            }
            SCIPfreeMemoryArray(scip, &propdata->sortedvars);
            propdata->nsortedvars = 0;
            propdata->nsortedbinvars = 0;

            /* get new variables */
            nnewvars = SCIPgetNVars(scip);
            newvars = SCIPgetVars(scip);
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &(propdata->sortedvars), newvars, nnewvars) ); /*lint !e666*/
            propdata->nsortedvars = nnewvars;

            nnewbinvars = SCIPgetNBinVars(scip);
            nnewintvars = SCIPgetNIntVars(scip);
            nnewimplvars = SCIPgetNImplVars(scip);

            /* determine implicit binary variables */
            lastidx = nnewbinvars + nnewintvars + nnewimplvars;
            for( v = nnewbinvars; v < lastidx; ++v )
            {
               if( SCIPvarIsBinary(propdata->sortedvars[v]) )
               {
                  SCIPswapPointers((void**) &(propdata->sortedvars[nnewbinvars]), (void**) &(propdata->sortedvars[v]));
                  ++nnewbinvars;
               }
            }
            propdata->nsortedbinvars = nnewbinvars;

            nbinvars = nnewbinvars;
            vars = propdata->sortedvars;
            nvars = propdata->nsortedvars;

            /* the probed variables change: the sparse bookkeeping is created anew after sorting them */
            if( sparse != NULL )
            {
               sparseFree(scip, sparse);
               sparse = NULL;
            }

            SCIP_CALL( SCIPreallocBufferArray(scip, &zeroimpllbs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &zeroimplubs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &zeroproplbs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &zeropropubs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &oneimpllbs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &oneimplubs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &oneproplbs, nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &onepropubs, nvars) );

            /* correct oldstartidx which is used for early termination */
            if( oldstartidx >= nbinvars )
               oldstartidx = nbinvars - 1;

            /* capture variables to make sure, the variables are not deleted */
            for( v = propdata->nsortedvars - 1; v >= 0; --v )
            {
               SCIP_CALL( SCIPcaptureVar(scip, propdata->sortedvars[v]) );
            }

            if( nnewbinvars == 0 )
            {
               *startidx = 0;
               propdata->lastsortstartidx = -1;
               propdata->nuseless = 0;
               propdata->ntotaluseless = 0;

               goto TERMINATE;
            }

            /* resorting here might lead to probing a second time on the same variable */
            SCIP_CALL( sortVariables(scip, propdata, propdata->sortedvars, propdata->nsortedbinvars, 0) );
            propdata->lastsortstartidx = 0;

            sparse = &sparsedata;
            SCIP_CALL( sparseCreate(scip, sparse, vars, nvars) );
         }
      }
   }
   while( i == 0 && !(*cutoff) && !(*delay) && !aborted );

   *startidx = i;

 TERMINATE:
   /* free temporary memory */
   if( sparse != NULL )
      sparseFree(scip, sparse);
   SCIPfreeBufferArray(scip, &onepropubs);
   SCIPfreeBufferArray(scip, &oneproplbs);
   SCIPfreeBufferArray(scip, &oneimplubs);
   SCIPfreeBufferArray(scip, &oneimpllbs);
   SCIPfreeBufferArray(scip, &zeropropubs);
   SCIPfreeBufferArray(scip, &zeroproplbs);
   SCIPfreeBufferArray(scip, &zeroimplubs);
   SCIPfreeBufferArray(scip, &zeroimpllbs);

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyProbing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   /* call inclusion method for propagator */
   SCIP_CALL( SCIPincludePropProbing(scip) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->sortedvars == NULL);
   assert(propdata->nsortedvars == 0);
   assert(propdata->nsortedbinvars == 0);

   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( initPropdata(propdata) );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &propdata->randnumgen,
      DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( freeSortedvars(scip, propdata) );
   assert(propdata->sortedvars == NULL);
   assert(propdata->nsortedvars == 0);
   assert(propdata->nsortedbinvars == 0);

   /* free random number generator */
   SCIPfreeRandom(scip, &propdata->randnumgen);

   return SCIP_OKAY;
}

/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->lastnode = -2;

   return SCIP_OKAY;
}


/** presolving deinitialization method of propagator (called after presolving has been finished) */
static
SCIP_DECL_PROPEXITPRE(propExitpreProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* delete the vars array, if the maximal number of runs are exceeded */
   if( propdata->maxruns >= 0 && SCIPgetNRuns(scip) >= propdata->maxruns )
   {
      SCIP_CALL( freeSortedvars(scip, propdata) );
      assert(propdata->sortedvars == NULL);
      assert(propdata->nsortedvars == 0);
      assert(propdata->nsortedbinvars == 0);
   }

   return SCIP_OKAY;
}


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolProbing)
{
  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* reset all propdata elements for stopping propagation earlier */
   propdata->nuseless = 0;
   propdata->ntotaluseless = 0;
   propdata->nsumuseless = 0;

   return SCIP_OKAY;
}


/** presolve method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int nvars;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldnimplications;
   int ntotalvars;
   SCIP_Bool delay;
   SCIP_Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);
   nimplvars = SCIPgetNImplVars(scip);

   /* if we have no binary variable anymore, we stop probing */
   if( nbinvars + nintvars + nimplvars == 0 )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* check, if probing should be applied in the current run */
   if( propdata->maxruns >= 0 && SCIPgetNRuns(scip) > propdata->maxruns )
      return SCIP_OKAY;

   /* if no domains changed since the last call, we don't need to probe */
   if( propdata->lastnode == -1 && nnewfixedvars == 0 && nnewaggrvars == 0 && nnewchgbds == 0 && nnewholes == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "executing probing (used %.1f sec)\n", SCIPpropGetTime(prop));

   *result = SCIP_DIDNOTFIND;

   /* allow some additional probings */
   propdata->nuseless -= propdata->nuseless/10;
   propdata->ntotaluseless -= propdata->ntotaluseless/10;

   /* get variable data */
   if( propdata->sortedvars == NULL )
   {
      SCIP_VAR** vars = SCIPgetVars(scip);
      int lastidx;
      int v;

      assert(propdata->startidx == 0);

      nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(propdata->sortedvars), vars, nvars) );
      propdata->nsortedvars = nvars;

      /* determine implicit binary variables */
      lastidx = nbinvars + nintvars + nimplvars;
      for( v = nbinvars; v < lastidx; ++v )
      {
         if( SCIPvarIsBinary(propdata->sortedvars[v]) )
         {
            SCIPswapPointers((void**) &(propdata->sortedvars[nbinvars]), (void**) &(propdata->sortedvars[v]));
            ++nbinvars;
         }
      }
      propdata->nsortedbinvars = nbinvars;

      /* capture variables to make sure, the variables are not deleted */
      for( v = propdata->nsortedvars - 1; v >= 0 ; --v )
      {
         SCIP_CALL( SCIPcaptureVar(scip, propdata->sortedvars[v]) );
      }
   }

   if( propdata->nsortedbinvars == 0 )
      return SCIP_OKAY;

   /* number of total variables is not decreasing, and we can identify every variable by their index, so allocate
    * enough space
    */
   ntotalvars = SCIPgetNTotalVars(scip);
   if( propdata->noldtotalvars < ntotalvars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->nprobed, ntotalvars) );
      BMSclearMemoryArray(&(propdata->nprobed[propdata->noldtotalvars]), ntotalvars - propdata->noldtotalvars); /*lint !e866*/
      propdata->noldtotalvars = ntotalvars;
   }

   propdata->lastnode = -1;

   /* sort the binary variables by number of rounding locks, if at least 100 variables were probed since last sort */
   if( propdata->lastsortstartidx < 0 || propdata->startidx - propdata->lastsortstartidx >= 100 )
   {
      SCIP_CALL( sortVariables(scip, propdata, propdata->sortedvars, propdata->nsortedbinvars, propdata->startidx) );
      propdata->lastsortstartidx = propdata->startidx;
   }

   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldnimplications = propdata->nimplications;

   /* start probing on variables */
   SCIP_CALL( applyProbing(scip, propdata, propdata->sortedvars, propdata->nsortedvars, propdata->nsortedbinvars,
         &(propdata->startidx), nfixedvars, naggrvars, nchgbds, oldnfixedvars, oldnaggrvars, &delay, &cutoff) );

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else
   {
      if( delay )
      {
         /* probing was interrupted because it reached the maximal fixings parameter, so we want to rerun it at the next call */
         propdata->lastnode = -2;
      }

      if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds
         || propdata->nimplications > oldnimplications )
         *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR** binvars;
   int nvars;
   int nbinvars;
   int i;
   int nfixedvars;
   int naggrvars;
   int nchgbds;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   SCIPdebug( int oldnimplications; )
   int startidx;
   int ntotalvars;
   SCIP_Bool delay;
   SCIP_Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* avoid recursive infinity loop */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* only call propagation on branching candidates, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* if already called stop */
   if( propdata->lastnode == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
      return SCIP_OKAY;

   /* if maximal depth for propagation is reached, stop */
   if( propdata->maxdepth >= 0 && propdata->maxdepth < SCIPgetDepth(scip) )
      return SCIP_OKAY;

   propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* get (number of) fractional variables that should be integral */
   /* todo check if integrating fractional implicit integer variables is beneficial for probing */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &vars, NULL, NULL, &nvars, NULL, NULL) );
   nbinvars = 0;

   /* alloc array for fractional binary variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nvars) );

   /* copy binary variables to array binvars */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      var = vars[i];

      assert(var != NULL);
      if( SCIPvarIsBinary(var) )
      {
         assert(SCIPvarGetLbLocal(var) < 0.5);
         assert(SCIPvarGetUbLocal(var) > 0.5);

         binvars[nbinvars] = var;
         ++nbinvars;
      }
   }
   SCIPdebugMsg(scip, "problem <%s> node %" SCIP_LONGINT_FORMAT " probing propagation found %d of %d possible probing candidates\n", SCIPgetProbName(scip), SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), nbinvars, nvars);

   if( nbinvars == 0 )
   {
      *result = SCIP_DIDNOTFIND;
      goto TERMINATE;
   }

   /* number of total variables is not decreasing, and we can identify every variable by their index, so allocate
    * enough space
    */
   ntotalvars = SCIPgetNTotalVars(scip);
   if( propdata->noldtotalvars < ntotalvars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->nprobed, ntotalvars) );
      BMSclearMemoryArray(&(propdata->nprobed[propdata->noldtotalvars]), ntotalvars - propdata->noldtotalvars); /*lint !e866*/
      propdata->noldtotalvars = ntotalvars;
   }

   /* sort binary variables */
   SCIP_CALL( sortVariables(scip, propdata, binvars, nbinvars, 0) );

   oldnfixedvars = 0;
   oldnaggrvars = 0;
   oldnchgbds = 0;
   nfixedvars = 0;
   naggrvars = 0;
   nchgbds = 0;
   startidx = 0;
   SCIPdebug( oldnimplications = propdata->nimplications; )

   /* start probing on found variables */
   SCIP_CALL( applyProbing(scip, propdata, binvars, nbinvars, nbinvars, &startidx, &nfixedvars, &naggrvars, &nchgbds, oldnfixedvars, oldnaggrvars, &delay, &cutoff) );
   SCIPdebug( SCIPdebugMsg(scip, "probing propagation found %d fixings, %d aggregation, %d nchgbds, and %d implications\n",
      nfixedvars, naggrvars, nchgbds, (propdata->nimplications) - oldnimplications); )

   if( delay )
   {
      /* probing was interrupted because it reached the maximal fixings parameter, so we want to rerun it at the next call */
      propdata->lastnode = -2;
   }

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > oldnfixedvars || naggrvars > oldnaggrvars || nchgbds > oldnchgbds )
      *result = SCIP_REDUCEDDOM;

 TERMINATE:
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY; /*lint !e438*/
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropProbing)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the probing propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create probing propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   SCIP_CALL( initPropdata(propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecProbing, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyProbing) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeProbing) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitProbing) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitProbing) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolProbing) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreProbing) );
   SCIP_CALL( SCIPsetPropExitpre(scip, prop, propExitpreProbing) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolProbing, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS,
         PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropProbing) );

   /* add probing propagator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxruns",
         "maximal number of runs, probing participates in (-1: no limit)",
         &propdata->maxruns, FALSE, DEFAULT_MAXRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/proprounds",
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &propdata->proprounds, TRUE, DEFAULT_PROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxfixings",
         "maximal number of fixings found, until probing is interrupted (0: don't iterrupt)",
         &propdata->maxfixings, TRUE, DEFAULT_MAXFIXINGS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxuseless",
         "maximal number of successive probings without fixings, until probing is aborted (0: don't abort)",
         &propdata->maxuseless, TRUE, DEFAULT_MAXUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxtotaluseless",
         "maximal number of successive probings without fixings, bound changes, and implications, until probing is aborted (0: don't abort)",
         &propdata->maxtotaluseless, TRUE, DEFAULT_MAXTOTALUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxsumuseless",
         "maximal number of probings without fixings, until probing is aborted (0: don't abort)",
         &propdata->maxsumuseless, TRUE, DEFAULT_MAXSUMUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxdepth",
         "maximal depth until propagation is executed(-1: no limit)",
         &propdata->maxdepth, TRUE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}


/** applies and evaluates probing of a single variable in the given direction and bound; with sparse bookkeeping, it
 *  stores the bounds only of the variables whose bounds probing changed, at their positions in the arrays, and records
 *  these positions for the direction; all other variables keep the bounds from before probing
 */
static
SCIP_RETCODE applyProbingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   probingpos,         /**< variable number to apply probing on */
   SCIP_BOUNDTYPE        boundtype,          /**< which bound should be changed */
   SCIP_Real             bound,              /**< which bound should be set */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   PROBINGSPARSE*        sparse,             /**< sparse bookkeeping, or NULL to store the bounds of all variables */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   SCIP_Bool dense;
   int direction;
   int i;

   assert(impllbs != NULL);
   assert(implubs != NULL);
   assert(proplbs != NULL);
   assert(propubs != NULL);
   assert(cutoff != NULL);
   assert(0 <= probingpos && probingpos < nvars);
   assert(SCIPisGE(scip, bound, SCIPvarGetLbLocal(vars[probingpos])));
   assert(SCIPisLE(scip, bound, SCIPvarGetUbLocal(vars[probingpos])));
   assert(sparse == NULL || SCIPvarIsBinary(vars[probingpos]));

   SCIPdebugMsg(scip, "applying probing on variable <%s> %s %g (nlocks=%d/%d, impls=%d/%d, clqs=%d/%d)\n",
      SCIPvarGetName(vars[probingpos]), boundtype == SCIP_BOUNDTYPE_UPPER ? "<=" : ">=", bound,
      SCIPvarGetNLocksDownType(vars[probingpos], SCIP_LOCKTYPE_MODEL),
      SCIPvarGetNLocksUpType(vars[probingpos], SCIP_LOCKTYPE_MODEL),
      SCIPvarGetNImpls(vars[probingpos], FALSE), SCIPvarGetNImpls(vars[probingpos], TRUE),
      SCIPvarGetNCliques(vars[probingpos], FALSE), SCIPvarGetNCliques(vars[probingpos], TRUE));

   /* debug mode also stores all bounds: analyzeDeductionsSparse() checks the sparse ones against them */
#ifdef NDEBUG
   dense = (sparse == NULL);
#else
   dense = TRUE;
#endif

   /* the sparse bookkeeping probes a binary variable: to 1 by its lower bound, to 0 by its upper bound */
   direction = (boundtype == SCIP_BOUNDTYPE_LOWER ? 1 : 0);
   if( sparse != NULL )
      sparse->nproppos[direction] = 0;

   /* in debug mode we assert above that this trivial infeasibility does not occur (for performance reasons), but in
    * optimized mode we return safely
    */
   if( SCIPisLT(scip, bound, SCIPvarGetLbLocal(vars[probingpos]))
         || SCIPisGT(scip, bound, SCIPvarGetUbLocal(vars[probingpos])) )
   {
      SCIPdebugMsg(scip, " -> trivial infeasibility detected\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* fix variable */
   if( boundtype == SCIP_BOUNDTYPE_UPPER )
   {
      SCIP_CALL( SCIPchgVarUbProbing(scip, vars[probingpos], bound) );
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      SCIP_CALL( SCIPchgVarLbProbing(scip, vars[probingpos], bound) );
   }

   /* apply propagation of implication graph and clique table */
   SCIP_CALL( SCIPpropagateProbingImplications(scip, cutoff) );

   if( !(*cutoff) )
   {
      if( sparse != NULL )
      {
         int* pos;
         int npos;

         /* proppos serves as scratch here: the propagation below runs on the same probing node, so it changes a superset
          * of these variables and overwrites the list
          */
         sparse->implstamp[direction] = ++sparse->stamp;
         pos = sparse->proppos[direction];
         sparseMark(scip, sparse, sparse->implmark[direction], sparse->implstamp[direction], pos, &npos);
         for( i = 0; i < npos; ++i )
         {
            impllbs[pos[i]] = SCIPvarGetLbLocal(vars[pos[i]]);
            implubs[pos[i]] = SCIPvarGetUbLocal(vars[pos[i]]);
         }
      }

      if( dense ) /*lint !e774*/
      {
         for( i = 0; i < nvars; ++i )
         {
            impllbs[i] = SCIPvarGetLbLocal(vars[i]);
            implubs[i] = SCIPvarGetUbLocal(vars[i]);
         }
      }

      /* apply propagation */
      SCIP_CALL( SCIPpropagateProbing(scip, maxproprounds, cutoff, NULL) );
   }
   else
   {
      SCIPdebugMsg(scip, "propagating probing implications after <%s> to %g led to a cutoff\n",
         SCIPvarGetName(vars[probingpos]), bound);
   }

   /* evaluate propagation */
   if( !(*cutoff) )
   {
      if( sparse != NULL )
      {
         int* pos;

         sparse->propstamp[direction] = ++sparse->stamp;
         pos = sparse->proppos[direction];
         sparseMark(scip, sparse, sparse->propmark[direction], sparse->propstamp[direction], pos,
            &sparse->nproppos[direction]);
         for( i = 0; i < sparse->nproppos[direction]; ++i )
         {
            proplbs[pos[i]] = SCIPvarGetLbLocal(vars[pos[i]]);
            propubs[pos[i]] = SCIPvarGetUbLocal(vars[pos[i]]);
         }
      }

      if( dense ) /*lint !e774*/
      {
         for( i = 0; i < nvars; ++i )
         {
            proplbs[i] = SCIPvarGetLbLocal(vars[i]);
            propubs[i] = SCIPvarGetUbLocal(vars[i]);
         }
      }
   }

   /* exit probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}

/** applies and evaluates probing of a single variable in the given direction and bound */
SCIP_RETCODE SCIPapplyProbingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   probingpos,         /**< variable number to apply probing on */
   SCIP_BOUNDTYPE        boundtype,          /**< which bound should be changed */
   SCIP_Real             bound,              /**< which bound should be set */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   SCIP_CALL( applyProbingVar(scip, vars, nvars, probingpos, boundtype, bound, maxproprounds, impllbs, implubs, proplbs,
         propubs, NULL, cutoff) );

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
   SCIP_Bool probingvarisbinary;
   SCIP_Bool probingvarisinteger;
   int j;

   assert(scip != NULL);
   assert(probingvar != NULL);
   assert(SCIPisGE(scip, leftub,  SCIPvarGetLbLocal(probingvar))); /* left  branch should not be empty by default */
   assert(SCIPisLE(scip, rightlb, SCIPvarGetUbLocal(probingvar))); /* right branch should not be empty by default */
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
   if( SCIPvarIsIntegral(probingvar) )
   {
      /* adjust bounds to actually used ones */
      leftub  = SCIPfloor(scip, leftub);
      rightlb = SCIPceil(scip, rightlb);

      probingvarisinteger = TRUE;
      probingvarisbinary = SCIPvarIsBinary(probingvar);
   }
   else
   {
      /* assert dichotomy in case of continuous var: leftub >= rightlb */
      assert(SCIPisGE(scip, leftub, rightlb));
      probingvarisbinary = FALSE;
      probingvarisinteger = FALSE;
   }

   /* check if probing variable was fixed in the branches */
   fixedleft  = SCIPisEQ(scip, SCIPvarGetLbLocal(probingvar), leftub);
   fixedright = SCIPisEQ(scip, SCIPvarGetUbLocal(probingvar), rightlb);

   *cutoff = FALSE;

   for( j = 0; j < nvars && !*cutoff; ++j )
   {
      SCIP_VAR* var;
      SCIP_Bool varisinteger;
      SCIP_Real newlb;
      SCIP_Real newub;

      assert(vars != NULL); /* for flexelint */

      var = vars[j];
      assert(var != NULL);

      /* @todo: add holes, and even add holes if x was the probing variable and it followed a better bound on x itself */
      /* @todo: check if we probed on an integer variable, that this maybe led to aggregation on two other variables, i.e
       *        probing on x <= 1 and x >= 2 led to y = 1, z = 1 and y = 0, z = 0 resp., which means y = Z
       */

      /* if probing variable is binary, then there is nothing we could deduce here (variable should be fixed in both branches)
       * if it is not binary, we want to see if we found bound tightenings, even though it seems quite unlikely */
      if( var == probingvar && probingvarisbinary )
         continue;

      /* new bounds of the variable is the union of the propagated bounds of the left and right case */
      newlb = MIN(leftproplbs[j], rightproplbs[j]);
      newub = MAX(leftpropubs[j], rightpropubs[j]);
      varisinteger = SCIPvarIsIntegral(var);

      /* check for fixed variables */
      if( SCIPisEQ(scip, newlb, newub) )
      {
         SCIP_Real fixval;
         SCIP_Bool fixed;

         if( !varisinteger )
         {
            /* in both probings, variable j is deduced to the same value: fix variable to this value */
            fixval = SCIPselectSimpleValue(newlb - 0.9 * SCIPepsilon(scip), newub + 0.9 * SCIPepsilon(scip), MAXDNOM);
         }
         else
         {
            fixval = newlb;
         }

         if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
         {
            SCIP_CALL( SCIPfixVar(scip, var, fixval, cutoff, &fixed) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, fixval, TRUE, cutoff, &fixed) );
            if( !*cutoff )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarUb(scip, var, fixval, TRUE, cutoff, &tightened) );
               fixed &= tightened;
            }
         }

         if( fixed )
         {
            SCIPdebugMsg(scip, "fixed variable <%s> to %g due to probing on <%s> with nlocks=(%d/%d)\n",
               SCIPvarGetName(var), fixval,
               SCIPvarGetName(probingvar), SCIPvarGetNLocksDownType(probingvar, SCIP_LOCKTYPE_MODEL),
               SCIPvarGetNLocksUpType(probingvar, SCIP_LOCKTYPE_MODEL));
            (*nfixedvars)++;
         }
         else if( *cutoff )
         {
            SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible fixing of variable <%s> to %g\n",
               SCIPvarGetName(probingvar), SCIPvarGetName(var), fixval);
         }

         continue;
      }
      else
      {
         /* check for bound tightenings */
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Bool tightened;
         SCIP_Bool tightenlb;
         SCIP_Bool tightenub;
         SCIP_Bool force;

         oldlb = SCIPvarGetLbLocal(var);
         oldub = SCIPvarGetUbLocal(var);

         if( varisinteger )
         {
            force = TRUE;
            tightenlb = (newlb > oldlb + 0.5);
            tightenub = (newub < oldub - 0.5);
         }
         else
         {
            force = TRUE;
            tightenlb = SCIPisLbBetter(scip, newlb, oldlb, oldub);
            tightenub = SCIPisUbBetter(scip, newub, oldlb, oldub);
         }

         if( tightenlb )
         {
            /* in both probings, variable j is deduced to be at least newlb: tighten lower bound */
            SCIP_CALL( SCIPtightenVarLb(scip, var, newlb, force, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMsg(scip, "tightened lower bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(var), oldlb, oldub, newlb,
                  SCIPvarGetName(probingvar), SCIPvarGetNLocksDownType(probingvar, SCIP_LOCKTYPE_MODEL),
                  SCIPvarGetNLocksUpType(probingvar, SCIP_LOCKTYPE_MODEL));
               (*nchgbds)++;
            }
            else if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible new lower bound of variable <%s> to %g\n",
                  SCIPvarGetName(probingvar), SCIPvarGetName(var), newlb);
            }
         }

         if( tightenub && !*cutoff )
         {
            /* in both probings, variable j is deduced to be at most newub: tighten upper bound */
            SCIP_CALL( SCIPtightenVarUb(scip, var, newub, force, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMsg(scip, "tightened upper bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(var), oldlb, oldub, newub,
                  SCIPvarGetName(probingvar), SCIPvarGetNLocksDownType(probingvar, SCIP_LOCKTYPE_MODEL),
                  SCIPvarGetNLocksUpType(probingvar, SCIP_LOCKTYPE_MODEL));
               (*nchgbds)++;
            }
            else if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible new lower bound of variable <%s> to %g\n",
                  SCIPvarGetName(probingvar), SCIPvarGetName(var), newub);
            }
         }
         if( *cutoff )
            break;
      }

      /* below we add aggregations and implications between probingvar and var,
       * we don't want this if both variables are the same
       */
      if( var == probingvar )
         continue;

      /* check for aggregations and implications */
      if( fixedleft && fixedright &&
          SCIPisEQ(scip, leftproplbs[j], leftpropubs[j]) && SCIPisEQ(scip, rightproplbs[j], rightpropubs[j]) )
      {
         /* var is fixed whenever probingvar is fixed, i.e.,
          *   var = leftproplbs[j] + (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub) * (probingvar - leftub)
          * -> both variables can be aggregated:
          *    (rightlb - leftub) * (var - leftproplbs[j]) = (rightproplbs[j] - leftproplbs[j]) * (probingvar - leftub)
          * -> (rightlb - leftub) * var - (rightproplbs[j] - leftproplbs[j]) * probingvar = leftproplbs[j] * rightlb - rightproplbs[j] * leftub
          *
          * check for case where both variables are binary: leftub = 1, rightlb = 0
          * case leftproplbs[j] = 0, rightproplbs[j] = 1, i.e., var and probingvar are fixed to same value
          *    -> aggregation is 1 * var - 1 * probingvar = 0 * 1 - 1 * 0 = 0 -> correct
          * case leftproplbs[j] = 1, rightproblbs[j] = 0, i.e., var and probingvar are fixed to opposite values
          *    -> aggregation is 1 * var + 1 * probingvar = 1 * 1 - 0 * 0 = 0 -> correct
          */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
         {
            SCIP_Bool aggregated;
            SCIP_Bool redundant;

            SCIP_CALL( SCIPaggregateVars(scip, var, probingvar,
               rightlb - leftub, -(rightproplbs[j] - leftproplbs[j]), leftproplbs[j] * rightlb - rightproplbs[j] * leftub,
               cutoff, &redundant, &aggregated) );

            if( aggregated )
            {
               SCIPdebugMsg(scip, "aggregated variables %g<%s> - %g<%s> == %g, nlocks=(%d/%d)\n",
                  rightlb - leftub, SCIPvarGetName(var),
                  rightproplbs[j] - leftproplbs[j], SCIPvarGetName(probingvar),
                  leftproplbs[j] * rightlb - rightproplbs[j] * leftub,
                  SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL),
                  SCIPvarGetNLocksUpType(probingvar, SCIP_LOCKTYPE_MODEL));
               (*naggrvars)++;
            }
            if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible aggregation: %g<%s> - %g<%s> == %g\n",
                  SCIPvarGetName(probingvar), rightlb - leftub, SCIPvarGetName(var),
                  rightproplbs[j] - leftproplbs[j], SCIPvarGetName(probingvar),
                  leftproplbs[j] * rightlb - rightproplbs[j] * leftub);
            }
         }
         else if( probingvarisinteger && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
         {
            /* if we are not in presolving, then we cannot do aggregations
             * but we can use variable bounds to code the same equality
             * var == ((leftproplbs[j] * rightlb - rightproplbs[j] * leftub) + (rightproplbs[j] - leftproplbs[j]) * probingvar) / (rightlb - leftub)
             */
            int nboundchanges;

            assert(!SCIPisEQ(scip, leftub, rightlb));

            SCIP_CALL( SCIPaddVarVlb(scip, var, probingvar, (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub), (leftproplbs[j] * rightlb - rightproplbs[j] * leftub) / (rightlb - leftub), cutoff, &nboundchanges) );
            (*nchgbds) += nboundchanges;

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible vlb: %g<%s> - %g<%s> == %g\n",
                  SCIPvarGetName(probingvar), rightlb - leftub, SCIPvarGetName(var),
                  rightproplbs[j] - leftproplbs[j], SCIPvarGetName(probingvar),
                  leftproplbs[j] * rightlb - rightproplbs[j] * leftub);
            }
            else
            {
               SCIP_CALL( SCIPaddVarVub(scip, var, probingvar, (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub), (leftproplbs[j] * rightlb - rightproplbs[j] * leftub) / (rightlb - leftub), cutoff, &nboundchanges) );
               (*nchgbds) += nboundchanges;

               if( *cutoff )
               {
                  SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible vub: %g<%s> - %g<%s> == %g\n",
                     SCIPvarGetName(probingvar), rightlb - leftub, SCIPvarGetName(var),
                     rightproplbs[j] - leftproplbs[j], SCIPvarGetName(probingvar),
                     leftproplbs[j] * rightlb - rightproplbs[j] * leftub);
               }
            }
            (*nimplications)++;
         }
         /* if probingvar is continuous and we are in solving stage, then we do nothing, but it's unlikely that we get
          * here (fixedleft && fixedright) with a continuous variable
          */
      }
      /* @todo: check if we can add variable lowerbounds/upperbounds on integer variables */
      /* can only add implications on binary variables which are globally valid */
      else if( probingvarisbinary && (SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0) )
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
            /* var is fixed to lower bound whenever probingvar is fixed to 0.0
             * and implication is not already known
             * -> insert implication: probingvar == 0  =>  var <= leftpropubs[j]
             */
            /*SCIPdebugMsg(scip, "found implication <%s> == 0  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(var), leftpropubs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, var, SCIP_BOUNDTYPE_UPPER, leftpropubs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 0  =>  <%s> == %g\n",
                  SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), leftpropubs[j]);
            }
         }
         else if( SCIPisEQ(scip, newub, leftproplbs[j]) && (leftimpllbs == NULL || leftimpllbs[j] < leftproplbs[j]) )
         {
            /* var is fixed to upper bound whenever probingvar is fixed to 0.0
             * and implication is not already known
             * -> insert implication: probingvar == 0  =>  var >= leftproplbs[j]
             */
            /*SCIPdebugMsg(scip, "found implication <%s> == 0  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(var), leftproplbs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, var, SCIP_BOUNDTYPE_LOWER, leftproplbs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 0  =>  <%s> == %g\n",
                  SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), leftproplbs[j]);
            }
         }
         /* we can do an else here, since the case where var is fixed for both fixings of probingvar had been handled as aggregation */
         else if( SCIPisEQ(scip, newlb, rightpropubs[j]) && (rightimplubs == NULL || rightimplubs[j] > rightpropubs[j]) )
         {
            /* var is fixed to lower bound whenever probingvar is fixed to 1.0
             * and implication is not already known
             * -> insert implication: probingvar == 1  =>  var <= rightpropubs[j]
             */
            /*SCIPdebugMsg(scip, "found implication <%s> == 1  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(var), rightpropubs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, var, SCIP_BOUNDTYPE_UPPER, rightpropubs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 1  =>  <%s> == %g\n",
                  SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), rightpropubs[j]);
            }
         }
         else if( SCIPisEQ(scip, newub, rightproplbs[j]) && (rightimpllbs == NULL || rightimpllbs[j] < rightproplbs[j]) )
         {
            /* var is fixed to upper bound whenever probingvar is fixed to 1.0
             * and implication is not already known
             * -> insert implication: probingvar == 1  =>  var >= leftproplbs[j]
             */
            /*SCIPdebugMsg(scip, "found implication <%s> == 1  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(var), rightproplbs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, var, SCIP_BOUNDTYPE_LOWER, rightproplbs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 1  =>  <%s> == %g\n",
                  SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), rightproplbs[j]);
            }
         }
         else if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPvarIsImpliedIntegral(var) )
         {
            /* check for implications for lower or upper bounds (only store implications with bounds tightened at least by 0.5)
             * in case of binary variables, this should have been handled in the previous cases, since every boundchange also fixes the variable
             */
            if( leftpropubs[j] < newub - 0.5 && (leftimplubs == NULL || leftpropubs[j] < leftimplubs[j]) )
            {
               /* insert implication: probingvar == 0  =>  var <= leftpropubs[j] */
               /*SCIPdebugMsg(scip, "found implication <%s> == 0  =>  <%s>[%g,%g] <= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(var), newlb, newub, leftpropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, var, SCIP_BOUNDTYPE_UPPER, leftpropubs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;

               if( *cutoff )
               {
                  SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 0  =>  <%s> <= %g\n",
                     SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), leftpropubs[j]);
               }
            }
            if( leftproplbs[j] > newlb + 0.5 && (leftimpllbs == NULL || leftproplbs[j] > leftimpllbs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 0  =>  var >= leftproplbs[j] */
               /*SCIPdebugMsg(scip, "found implication <%s> == 0  =>  <%s>[%g,%g] >= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(var), newlb, newub, leftproplbs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, var, SCIP_BOUNDTYPE_LOWER, leftproplbs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;

               if( *cutoff )
               {
                  SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 0  =>  <%s> >= %g\n",
                     SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), leftproplbs[j]);
               }
            }
            if( rightpropubs[j] < newub - 0.5 && (rightimplubs == NULL || rightpropubs[j] < rightimplubs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 1  =>  var <= rightpropubs[j] */
               /*SCIPdebugMsg(scip, "found implication <%s> == 1  =>  <%s>[%g,%g] <= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(var), newlb, newub, rightpropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, var, SCIP_BOUNDTYPE_UPPER, rightpropubs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;

               if( *cutoff )
               {
                  SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 1  =>  <%s> <= %g\n",
                     SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), rightpropubs[j]);
               }
            }
            if( rightproplbs[j] > newlb + 0.5 && (rightimpllbs == NULL || rightproplbs[j] > rightimpllbs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 1  =>  var >= rightproplbs[j] */
               /*SCIPdebugMsg(scip, "found implication <%s> == 1  =>  <%s>[%g,%g] >= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(var), newlb, newub, rightproplbs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, var, SCIP_BOUNDTYPE_LOWER, rightproplbs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;

               if( *cutoff )
               {
                  SCIPdebugMsg(scip, "analyzing probing deduction of <%s> led to an infeasible implication <%s> == 1  =>  <%s> <= %g\n",
                     SCIPvarGetName(probingvar), SCIPvarGetName(probingvar), SCIPvarGetName(var), rightproplbs[j]);
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}
