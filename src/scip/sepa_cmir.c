/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_cmir.c,v 1.58 2006/08/31 08:27:28 bzfpfend Exp $"

/**@file   sepa_cmir.c
 * @brief  complemented mixed integer rounding cuts separator (Marchand's version)
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_cmir.h"


#define SEPA_NAME              "cmir"
#define SEPA_DESC              "complemented mixed integer rounding cuts separator (Marchand's version)"
#define SEPA_PRIORITY             -3000
#define SEPA_FREQ                     0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             3 /**< maximal number of cmir separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        10 /**< maximal number of cmir separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXTRIES            100 /**< maximal number of rows to start aggregation with per separation round
                                         *   (-1: unlimited) */
#define DEFAULT_MAXTRIESROOT         -1 /**< maximal number of rows to start aggregation with per round in the root node
                                         *   (-1: unlimited) */
#define DEFAULT_MAXFAILS             20 /**< maximal number of consecutive unsuccesful aggregation tries (-1: unlimited) */
#define DEFAULT_MAXFAILSROOT        100 /**< maximal number of consecutive unsuccesful aggregation tries in the root node
                                         *   (-1: unlimited) */
#define DEFAULT_MAXAGGRS              3 /**< maximal number of aggregations for each row per separation round */
#define DEFAULT_MAXAGGRSROOT          6 /**< maximal number of aggreagtions for each row per round in the root node */
#define DEFAULT_MAXSEPACUTS         100 /**< maximal number of cmir cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of cmir cuts separated per separation round in root node */
#define DEFAULT_MAXSLACK            0.0 /**< maximal slack of rows to be used in aggregation */
#define DEFAULT_MAXSLACKROOT        0.1 /**< maximal slack of rows to be used in aggregation in the root node */
#define DEFAULT_DENSITYSCORE      1e-04 /**< weight of row density in the aggregation scoring of the rows */
#define DEFAULT_SLACKSCORE        1e-03 /**< weight of slack in the aggregation scoring of the rows */
#define DEFAULT_MAXAGGDENSITY      0.20 /**< maximal density of aggregated row */
#define DEFAULT_MAXROWDENSITY      0.05 /**< maximal density of row to be used in aggregation */
#define DEFAULT_MAXROWFAC          1e+4 /**< maximal row aggregation factor */
#define DEFAULT_MAXTESTDELTA         -1 /**< maximal number of different deltas to try (-1: unlimited) */
#define DEFAULT_MAXTESTDELTAROOT     -1 /**< maximal number of different deltas to try in the root node (-1: unlimited) */
#define DEFAULT_MAXCONTS             10 /**< maximal number of active continuous variables in aggregated row */
#define DEFAULT_MAXCONTSROOT         10 /**< maximal number of active continuous variables in aggregated row in the root */
#define DEFAULT_AGGRTOL             0.1 /**< aggregation heuristic: tolerance for bounddistances used to select real 
                                         *   variable in current aggregated constraint to be eliminated */
#define DEFAULT_TRYNEGSCALING      TRUE /**< should negative values also be tested in scaling? */
#define DEFAULT_FIXINTEGRALRHS     TRUE /**< should an additional variable be complemented if f0 = 0? */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */

#define BOUNDSWITCH                 0.5
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE
#define MINFRAC                    0.05
#define MAXFRAC                    1.00
#define MAKECONTINTEGRAL          FALSE




/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_Real             maxslack;           /**< maximal slack of rows to be used in aggregation */
   SCIP_Real             maxslackroot;       /**< maximal slack of rows to be used in aggregation in the root node */
   SCIP_Real             densityscore;       /**< weight of row density in the aggregation scoring of the rows */
   SCIP_Real             slackscore;         /**< weight of slack in the aggregation scoring of the rows */
   SCIP_Real             maxaggdensity;      /**< maximal density of aggregated row */
   SCIP_Real             maxrowdensity;      /**< maximal density of row to be used in aggregation */
   SCIP_Real             maxrowfac;          /**< maximal row aggregation factor */
   SCIP_Real             aggrtol;            /**< tolerance for bounddistance used in aggregation heuristic */
   int                   maxrounds;          /**< maximal number of cmir separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of cmir separation rounds in the root node (-1: unlimited) */
   int                   maxtries;           /**< maximal number of rows to start aggregation with per separation round
                                              *   (-1: unlimited) */
   int                   maxtriesroot;       /**< maximal number of rows to start aggregation with per round in the root node
                                              *   (-1: unlimited) */
   int                   maxfails;           /**< maximal number of consecutive unsuccesful aggregation tries
                                              *   (-1: unlimited) */
   int                   maxfailsroot;       /**< maximal number of consecutive unsuccesful aggregation tries in the root
                                              *   node (-1: unlimited) */
   int                   maxaggrs;           /**< maximal number of aggregations for each row per separation round */
   int                   maxaggrsroot;       /**< maximal number of aggregations for each row per round in the root node */
   int                   maxsepacuts;        /**< maximal number of cmir cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cmir cuts separated per separation round in root node */
   int                   maxtestdelta;	     /**< maximal number of different deltas to try (-1: unlimited) */
   int                   maxtestdeltaroot;   /**< maximal number of different deltas to try in the root node (-1: unlimited) */
   int                   maxconts;	     /**< maximal number of active continuous variables in aggregated row */
   int                   maxcontsroot;       /**< maximal number of active continuous variables in aggregated row in the root */
   SCIP_Bool             trynegscaling;      /**< should negative values also be tested in scaling? */
   SCIP_Bool             fixintegralrhs;     /**< should an additional variable be complemented if f0 = 0? */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};




/*
 * Local methods
 */

/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   SCIP_VAR**            cutvars,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real val;
   SCIP_Real absval;
   SCIP_Real cutsqrnorm;
   SCIP_Real act;
   SCIP_Real norm;
   int len;
   int v;

   assert(nvars == 0 || cutcoefs != NULL);
   assert(nvars == 0 || varsolvals != NULL);
   assert(cutvars != NULL);
   assert(cutvals != NULL);
   assert(cutlen != NULL);
   assert(cutact != NULL);
   assert(cutnorm != NULL);

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch( normtype )
   {
   case 'e':
      cutsqrnorm = 0.0;
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            absval = REALABS(val);
            norm = MAX(norm, absval);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 's':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 'd':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm = 1.0;
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}

/** adds given cut to LP if violated */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            varsolvals,         /**< solution values of active variables */
   SCIP_Real*            cutcoefs,           /**< coefficients of active variables in cut */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Bool             cutislocal,         /**< is the cut only locally valid? */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   int*                  ncuts               /**< pointer to count the number of added cuts */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_Real cutact;
   SCIP_Real cutnorm;
   int nvars;
   int cutlen;
   SCIP_Bool success;
   
   assert(scip != NULL);
   assert(sepadata != NULL);      
   assert(varsolvals != NULL);
   assert(cutcoefs != NULL);
   assert(ncuts != NULL);

   /* get active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars == 0 || vars != NULL);

   /* get temporary memory for storing the cut as sparse row */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, nvars) );
   
   /* store the cut as sparse row, calculate activity and norm of cut */
   SCIP_CALL( storeCutInArrays(scip, nvars, vars, cutcoefs, varsolvals, normtype,
         cutvars, cutvals, &cutlen, &cutact, &cutnorm) );

   if( SCIPisPositive(scip, cutnorm) && SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm) )
   {
      SCIP_ROW* cut;
      char cutname[SCIP_MAXSTRLEN];
      
      /* create the cut */
      sprintf(cutname, "cmir%d_%d", SCIPgetNLPs(scip), *ncuts);
      SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), cutrhs, 
            cutislocal, FALSE, sepadata->dynamiccuts) );
      SCIP_CALL( SCIPaddVarsToRow(scip, cut, cutlen, cutvars, cutvals) );

      SCIPdebugMessage(" -> found potential c-mir cut <%s>: activity=%f, rhs=%f, norm=%f, eff=%f\n",
         cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut));
      SCIPdebug(SCIPprintRow(scip, cut, NULL));
      
      /* try to scale the cut to integral values, but only if the scaling is small; otherwise keep the fractional cut */
      SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
            30, 100.0, MAKECONTINTEGRAL, &success) );
      if( success && !SCIPisCutEfficacious(scip, sol, cut) )
      {
         SCIPdebugMessage(" -> c-mir cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
            cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut));
         SCIPdebug(SCIPprintRow(scip, cut, NULL));
         success = FALSE;
      }
      else
         success = TRUE; /* also use cut if scaling failed */

      /* if scaling was successful, add the cut */
      if( success ) /*lint !e774*/ /* Boolean within 'if' always evaluates to True */
      {
         SCIPdebugMessage(" -> found c-mir cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%g)\n",
            cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut),
            SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
            SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
         SCIPdebug(SCIPprintRow(scip, cut, NULL));
         SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE) );
         if( !cutislocal )
         {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         }
         (*ncuts)++;
      }
      
      /* release the row */
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }
   
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cutvals);
   SCIPfreeBufferArray(scip, &cutvars);

   return SCIP_OKAY;   
}

/** adds delta to active continuous variables counter */
static
void updateNActiveConts(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_Real*            varsolvals,         /**< LP solution value of all variables in LP */
   SCIP_Real*            bestcontlbs,        /**< best lower (variable or standard) bounds of continuous variables */
   SCIP_Real*            bestcontubs,        /**< best upper (variable or standard) bounds of continuous variables */
   int                   nintvars,           /**< number of integer variables */
   SCIP_VAR*             var,                /**< continuous variable */
   int                   delta,              /**< delta value of counters */
   int*                  nactiveconts        /**< pointer to count number of active continuous variabls */
   )
{
   SCIP_Real primsol;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(nactiveconts != NULL);

   primsol = varsolvals[SCIPvarGetProbindex(var)];
   lb = bestcontlbs[SCIPvarGetProbindex(var) - nintvars];
   ub = bestcontubs[SCIPvarGetProbindex(var) - nintvars];

   if( SCIPisLT(scip, lb, primsol) && SCIPisLT(scip, primsol, ub) )
      (*nactiveconts) += delta;
}

/** decreases the score of a row in order to not aggregate it again too soon */
static
void decreaseRowScore(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_Real*            rowlhsscores,       /**< aggregation scores for left hand sides of row */
   SCIP_Real*            rowrhsscores,       /**< aggregation scores for right hand sides of row */
   int                   rowidx              /**< index of row to decrease score for */
   )
{
   assert(rowlhsscores != NULL);
   assert(rowrhsscores != NULL);
   assert(rowlhsscores[rowidx] >= 0.0);
   assert(rowrhsscores[rowidx] >= 0.0);

   rowlhsscores[rowidx] *= 0.9;
   rowrhsscores[rowidx] *= 0.9;
}

/** calculates efficacy of the given cut */
static
SCIP_Real calcEfficacy(
   int                   nvars,              /**< number of variables in the problem */
   SCIP_Real*            cutcoefs,           /**< dense vector of cut coefficients */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Real             cutact              /**< activity of cut */
   )
{
   SCIP_Real sqrnorm;
   int i;
   
   assert(cutcoefs != NULL);

   if( cutact <= cutrhs )
      return -1.0;

   sqrnorm = 0.0;
   for( i = 0; i < nvars; ++i )
      sqrnorm += SQR(cutcoefs[i]);
   sqrnorm = MAX(sqrnorm, 1e-12);

   return (cutact - cutrhs)/SQRT(sqrnorm);
}

/** calculates the c-MIR cut for the given rowweights and delta value, and updates testeddeltas, bestdelta, and
 *  bestefficacy
 */
static
SCIP_RETCODE tryDelta(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   int                   nvars,              /**< number of problem variables */
   SCIP_Real*            rowweights,         /**< weight of rows in aggregated row */ 
   SCIP_Real*            cutcoefs,           /**< array to store the cut coefficients */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Real*            testeddeltas,       /**< array with already tested deltas */
   int*                  ntesteddeltas,      /**< pointer to the number of elements in testeddeltas */
   SCIP_Real             delta,              /**< delta value to scale mixed knapsack equation with */
   SCIP_Real*            bestdelta,          /**< pointer to the currently best delta value */
   SCIP_Real*            bestefficacy        /**< pointer to the currently best efficacy */
   )
{
   SCIP_Bool tested;
   int i;

   assert(sepadata != NULL);
   assert(testeddeltas != NULL);
   assert(ntesteddeltas != NULL);
   assert(bestdelta != NULL);
   assert(bestefficacy != NULL);

   /* do not use too small deltas */
   if( SCIPisFeasZero(scip, delta) )
      return SCIP_OKAY;

   /* check, if delta with mult was already tested */
   tested = FALSE;
   for( i = 0; i < *ntesteddeltas && !tested; i++ )
      tested = SCIPisEQ(scip, testeddeltas[i], delta);
   if( !tested )
   {
      SCIP_Real cutact;
      SCIP_Real cutrhs;
      SCIP_Bool success;
      SCIP_Bool cutislocal;

      testeddeltas[*ntesteddeltas] = delta;
      (*ntesteddeltas)++;

      /* create a MIR cut out of the weighted LP rows */
      SCIP_CALL( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, sepadata->fixintegralrhs, sepadata->maxrowfac,
            MINFRAC, MAXFRAC, rowweights, delta, mksetcoefs, cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
      assert(ALLOWLOCAL || !cutislocal);
      SCIPdebugMessage("delta = %g  -> success: %d\n", delta, success);
               
      /* check if delta generates cut which is more violated */
      if( success )
      {
         SCIP_Real efficacy;

         efficacy = calcEfficacy(nvars, cutcoefs, cutrhs, cutact);
         SCIPdebugMessage("act = %g  rhs = %g  eff = %g, old besteff = %g, old bestdelta=%g\n", 
            cutact, cutrhs, efficacy, *bestefficacy, *bestdelta);
         if( efficacy > *bestefficacy )
         {
            *bestdelta = delta;
            *bestefficacy = efficacy;
         }
      }
   }

   return SCIP_OKAY;
}

/** Performs the cut generation heuristic of the c-MIR separation algorithm, i.e., tries to generate a c-MIR cut which is
 *  valid for the mixed knapsack set corresponding to the current aggregated constraint.
 */
static
SCIP_RETCODE cutGenerationHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_COL**            cols,               /**< columns of the LP */
   int                   ncols,              /**< number of columns in the LP */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   nintvars,           /**< number of integer variables */
   SCIP_Real*            varsolvals,         /**< LP solution value of all variables in LP */
   int                   maxtestdelta,       /**< maximal number of different deltas to try (-1: unlimited) */
   SCIP_Real*            rowweights,         /**< weight of rows in aggregated row */ 
   SCIP_Real*            aggrcoefs,          /**< coefficients of all variables in aggregated row */
   int*                  aggrintnonzposs,    /**< positions of non-zero coefficients for integer variables in aggrcoefs */
   int*                  aggrnonzidxs,       /**< position of column in aggrintnonzposs array, or -1 */
   int                   naggrintnonzs,      /**< number of non-zero coefficients for integer variables in aggrcoefs */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   int*                  ncuts               /**< pointer to count the number of generated cuts */
)
{
   SCIP_Real* cutcoefs;        
   SCIP_Real* mksetcoefs;
   SCIP_Real* testeddeltas;
   SCIP_Real bestdelta;
   SCIP_Real bestefficacy; 
   SCIP_Real maxabsmksetcoef;
   int ntesteddeltas;
   int vi;

   if( maxtestdelta == -1 )
      maxtestdelta = INT_MAX;

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &mksetcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &testeddeltas, 2*(nintvars+1)) );

   /* As in Marchand's version. Use the absolute value of the coefficients of the integer variables (lying 
    * stricly between its bounds) in the constructed mixed knapsack set, i.e., 
    * N* = { |alpha'_j| : j in N, alpha'_j != 0 and l_j < x*_j < u_j }
    */ 
         
   /* search delta for generating a cut with maximum efficacy: 
    * delta = coefficient of integer variable in constructed mixed knapsack set which lies between its bounds
    */ 
   ntesteddeltas = 0;
   bestdelta = 0.0;
   bestefficacy = 0.0;
   maxabsmksetcoef = 0.0;

   /* try delta = 1 and get the coefficients of all variables in the constructed mixed knapsack set */
   SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, mksetcoefs, testeddeltas, &ntesteddeltas, 1.0,
         &bestdelta, &bestefficacy) );
   if( sepadata->trynegscaling )
   {
      SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, NULL, testeddeltas, &ntesteddeltas, -1.0,
            &bestdelta, &bestefficacy) );
   }

   /* find mult in { +1, -1 } and delta in the corresponding set N* leading to the most violated c-MIR cut */
   for( vi = 0; vi < nintvars; vi++ )
   {
      SCIP_VAR* var;
      SCIP_Real primsol;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real absmksetcoef;
            
      var = vars[vi];
      assert(vi == SCIPvarGetProbindex(var));
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
      assert(SCIPvarIsActive(var));
      assert(SCIPvarIsIntegral(var));

      /* update maximum coefficient of integer variables in constructed mixed knapsack set for 
       *   mult = +1 and delta = 1 and 
       *   mult = -1 and delta = 1 
       */
      absmksetcoef = REALABS(mksetcoefs[vi]);
      maxabsmksetcoef = MAX(maxabsmksetcoef, absmksetcoef);

      if( ntesteddeltas >= maxtestdelta )
         continue; /* remaining loop is only for maxabsmksetcoef calculation */

      /* ignore variables with current solution value on its bounds */
      primsol = varsolvals[vi];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      if( SCIPisEQ(scip, primsol, lb) || SCIPisEQ(scip, primsol, ub) )
         continue;

      /* try to divide aggregated row by absmksetcoef */
      if( !SCIPisFeasZero(scip, absmksetcoef) )
      {
         SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, NULL, testeddeltas, &ntesteddeltas,
               1.0/absmksetcoef, &bestdelta, &bestefficacy) );
         if( sepadata->trynegscaling )
         {
            SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, NULL, testeddeltas, &ntesteddeltas,
                  -1.0/absmksetcoef, &bestdelta, &bestefficacy) );
         }
      }
   }

   /* additionally try delta = maxabscoef+1 */
   if( !SCIPisFeasZero(scip, maxabsmksetcoef) )
   {
      SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, NULL, testeddeltas, &ntesteddeltas,
            1.0/(maxabsmksetcoef+1.0), &bestdelta, &bestefficacy) );
      if( sepadata->trynegscaling )
      {
         SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, NULL, testeddeltas, &ntesteddeltas,
               -1.0/(maxabsmksetcoef+1.0), &bestdelta, &bestefficacy) );
      }
   }

   /* delta found */
   if( SCIPisEfficacious(scip, bestefficacy) )
   {
      SCIP_Real delta;
      SCIP_Real cutrhs;           
      SCIP_Real cutact;
      SCIP_Bool success;
      SCIP_Bool cutislocal;
      int i;

      assert(!SCIPisFeasZero(scip, bestdelta));

      /* Try to improve efficacy by multiplying delta with 2, 4 and 8 */
      for( i = 0, delta = 2.0 * bestdelta; i < 3; i++, delta *= 2.0 )
      {
         SCIP_CALL( tryDelta(scip, sepadata, nvars, rowweights, cutcoefs, NULL, testeddeltas, &ntesteddeltas,
               delta, &bestdelta, &bestefficacy) );
      }

      /* generate cut with bestdelta and best boundswitch value */
      SCIP_CALL( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, sepadata->fixintegralrhs, sepadata->maxrowfac,
            MINFRAC, MAXFRAC, rowweights, bestdelta, NULL, cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
      assert(ALLOWLOCAL || !cutislocal);
      assert(success); 
#ifndef NDEBUG
      {
         SCIP_Real recalcedefficacy;
         recalcedefficacy = calcEfficacy(nvars, cutcoefs, cutrhs, cutact);
         assert(SCIPisEQ(scip, bestefficacy, recalcedefficacy));
      }
#endif
      SCIP_CALL( addCut(scip, sepadata, sol, varsolvals, cutcoefs, cutrhs, cutislocal, normtype, ncuts) );
   }

   /* free datastructures */
   SCIPfreeBufferArray(scip, &testeddeltas);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &mksetcoefs);

   return SCIP_OKAY; 
}

/** aggregates different single mixed integer constraints by taking linear combinations of the rows of the LP  */
static
SCIP_RETCODE aggregation(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            varsolvals,         /**< LP solution value of all variables in LP */
   SCIP_Real*            bestcontlbs,        /**< best lower (variable or standard) bounds of continuous variables */
   SCIP_Real*            bestcontubs,        /**< best upper (variable or standard) bounds of continuous variables */
   SCIP_Real*            rowlhsscores,       /**< aggregation scores for left hand sides of row */
   SCIP_Real*            rowrhsscores,       /**< aggregation scores for right hand sides of row */
   int                   startrow,           /**< index of row to start aggregation */ 
   int                   maxaggrs,           /**< maximal number of aggregations */
   SCIP_Real             maxslack,           /**< maximal slack of rows to be used in aggregation */
   int                   maxtestdelta,       /**< maximal number of different deltas to try (-1: unlimited) */
   int                   maxconts,           /**< maximal number of active continuous variables in aggregated row */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   SCIP_Bool*            wastried,           /**< pointer to store whether the given startrow was actually tried */
   int*                  ncuts               /**< pointer to count the number of generated cuts */
   )
{
   SCIP_COL** startnonzcols;   
   SCIP_COL** cols;
   SCIP_VAR** vars;
   SCIP_ROW** rows;
   SCIP_COL* bestcol;          
   SCIP_Real* startnonzcoefs;      
   SCIP_Real* aggrcoefs;       
   SCIP_Real* rowweights;       
   int* aggrnonzidxs;
   int* aggrintnonzposs;
   int* aggrcontnonzposs;
   SCIP_Real maxweight;
   SCIP_Real minweight;
   SCIP_Real startrowact;      
   SCIP_Bool hasfractional;
   int naggrintnonzs;
   int naggrcontnonzs;
   int maxaggrnonzs;
   int nstartnonzcols;         
   int naggrs;
   int nactiveconts;
   int nvars;
   int nintvars;
   int ncontvars;
   int ncols;
   int nrows;
   int c;
   int r;

   assert(scip != NULL);
   assert(sepadata != NULL);      
   assert(varsolvals != NULL);
   assert(rowlhsscores != NULL);
   assert(rowrhsscores != NULL);
   assert(wastried != NULL);
   assert(ncuts != NULL);

   *wastried = FALSE;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );
   nintvars = nvars - ncontvars;
   assert((nvars == 0 && nintvars == 0 && ncontvars == 0) || vars != NULL);
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   assert(ncols == 0 || cols != NULL);
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert(nrows == 0 || rows != NULL);
   assert(0 <= startrow && startrow < nrows);

   SCIPdebugMessage("start c-MIR aggregation with row <%s> (%d/%d)\n", SCIProwGetName(rows[startrow]), startrow, nrows);

   /* calculate maximal number of non-zeros in aggregated row */
   maxaggrnonzs = (int)(sepadata->maxaggdensity * ncols);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrnonzidxs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrintnonzposs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrcontnonzposs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowweights, nrows) );

   /* initialize weights of rows in aggregation */
   BMSclearMemoryArray(rowweights, nrows);
   startrowact = SCIPgetRowSolActivity(scip, rows[startrow], sol);
   if( startrowact <= 0.5 * SCIProwGetLhs(rows[startrow]) + 0.5 * SCIProwGetRhs(rows[startrow]) )
      rowweights[startrow] = -1.0;
   else 
      rowweights[startrow] = 1.0;
   maxweight = 1.0;
   minweight = 1.0;
   
   /* decrease score of startrow in order to not aggregate it again too soon */
   decreaseRowScore(scip, rowlhsscores, rowrhsscores, startrow);
   
   /* get nonzero columns and coefficients of startrow */
   startnonzcols =  SCIProwGetCols(rows[startrow]);
   nstartnonzcols = SCIProwGetNLPNonz(rows[startrow]);
   startnonzcoefs = SCIProwGetVals(rows[startrow]);
   
   /* for all columns of startrow store coefficient as coefficient in aggregated row */ 
   BMSclearMemoryArray(aggrcoefs, ncols);
   naggrintnonzs = 0;
   naggrcontnonzs = 0;
   nactiveconts = 0;
   hasfractional = FALSE;
   for( c = 0; c < nstartnonzcols; c++ )
   {
      SCIP_VAR* var;
      int pos;

      var = SCIPcolGetVar(startnonzcols[c]);
      pos = SCIPcolGetLPPos(startnonzcols[c]);
      assert(pos >= 0); 
      assert(!SCIPisZero(scip, startnonzcoefs[c]));
      aggrcoefs[pos] = rowweights[startrow] * startnonzcoefs[c];
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         updateNActiveConts(scip, varsolvals, bestcontlbs, bestcontubs, nintvars, var, +1, &nactiveconts);
         aggrnonzidxs[pos] = naggrcontnonzs;
         aggrcontnonzposs[naggrcontnonzs] = pos;
         naggrcontnonzs++;
      }
      else
      {
         aggrnonzidxs[pos] = naggrintnonzs;
         aggrintnonzposs[naggrintnonzs] = pos;
         naggrintnonzs++;
         if( !hasfractional )
         {
            SCIP_Real primsol;

            primsol = varsolvals[SCIPvarGetProbindex(var)];
            hasfractional = !SCIPisIntegral(scip, primsol);
         }
      }
   }

   /* don't try aggregation if there is no integer variable with fractional value */
   if( !hasfractional )
      maxaggrs = -1;

   /* try to generate cut from the current aggregated row 
    * add cut if found, otherwise add another row to aggregated row 
    * in order to get rid of a continuous variable
    */
   naggrs = 0;
   while( nactiveconts <= maxconts && naggrs <= maxaggrs && naggrcontnonzs + naggrintnonzs <= maxaggrnonzs )
   {
      SCIP_ROW* bestrow;          
      SCIP_COL** bestrownonzcols;     /* columns with nonzero coefficients in best row to add */
      SCIP_Real* bestrownonzcoefs;    /* nonzero coefficients of columns in best row to add */
      int nbestrownonzcols;           /* number of columns with nonzero coefficients in best row to add */
      SCIP_Real bestbounddist;
      SCIP_Real bestscore;
      SCIP_Real aggrfac;
      SCIP_Real absaggrfac;
      int nzi;
      int oldncuts;

      *wastried = TRUE;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("aggregation of startrow %d and %d additional rows with %d integer and %d continuous variables (%d active):\n",
         startrow, naggrs, naggrintnonzs, naggrcontnonzs, nactiveconts);
      for( c = 0; c < ncols; ++c )
      {
         if( aggrcoefs[c] != 0.0 )
            SCIPdebugPrintf(" %+g<%s>(%g)", aggrcoefs[c], SCIPvarGetName(SCIPcolGetVar(cols[c])),
               varsolvals[SCIPvarGetProbindex(SCIPcolGetVar(cols[c]))]);
      }
      SCIPdebugPrintf("\n");
#endif

      /* Step 1: 
       * try to generate a MIR cut out of the current aggregation 
       */
      oldncuts = *ncuts;
      SCIP_CALL( cutGenerationHeuristic(scip, sepadata, sol, cols, ncols, vars, nvars, nintvars, varsolvals, 
            maxtestdelta, rowweights, aggrcoefs, aggrintnonzposs, aggrnonzidxs, naggrintnonzs, normtype, ncuts) );
      
      /* if the cut was successfully added, abort the aggregation of further rows */
      if( *ncuts > oldncuts )
      {
         SCIPdebugMessage(" -> abort aggregation: cut found\n");
         break;
      }

      /* Step 2: 
       * aggregate an additional row in order to remove a continuous variable 
       */

      /* abort, if we reached the maximal number of aggregations */
      if( naggrs == maxaggrs )
      {
         SCIPdebugMessage(" -> abort aggregation: %s\n", nactiveconts == 0 ? "no more active continuous variables"
            : "maximal number of aggregations reached");
         break;
      }

      SCIPdebugMessage(" -> search column to eliminate\n");
      
      /* search for "best" continuous variable in aggregated row:
       * - solution value is strictly between lower and upper bound
       * - it exists a not yet aggregated row with nonzero coefficient in this column
       * out of these variables:
       * - prefer variables with larger distance of current solution value to its bounds
       * - of those with large bound distance, prefer variables that can be eliminated with a row of high score
       */
      bestcol = NULL;
      bestbounddist = -1.0;
      bestscore = 0.0;
      bestrow = NULL;
      aggrfac = 0.0;
      for( nzi = 0; nzi < naggrcontnonzs; ++nzi )
      {
         SCIP_COL* col;
         SCIP_VAR* var;
         SCIP_Real primsol;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real distlower;
         SCIP_Real distupper;
         SCIP_Real bounddist;

         c = aggrcontnonzposs[nzi];
         assert(0 <= c && c < ncols);
         assert(aggrnonzidxs[c] == nzi);
         assert(!SCIPisZero(scip, aggrcoefs[c]));

         col = cols[c];
         var = SCIPcolGetVar(col);
         assert(!SCIPvarIsIntegral(var));

         /* get minimum distance of LP solution value of variable to its bounds */
         primsol = varsolvals[SCIPvarGetProbindex(var)];
         lb = bestcontlbs[SCIPvarGetProbindex(var) - nintvars];
         ub = bestcontubs[SCIPvarGetProbindex(var) - nintvars];
         assert(SCIPisGE(scip, lb, SCIPvarGetLbGlobal(var)));
         assert(SCIPisLE(scip, ub, SCIPvarGetUbGlobal(var)));
         distlower = primsol - lb;
         distupper = ub - primsol;
         bounddist = MIN(distlower, distupper);

         /* check, if variable is candidate to be the new best variable */
         if( bounddist >= bestbounddist - sepadata->aggrtol )
         {
            SCIP_ROW** nonzrows;
            SCIP_Real* nonzcoefs;
            int nnonzrows;
            
            SCIPdebugMessage("     -> col <%s>[%g,%g]: sol=%g, dist=%g\n", 
               SCIPvarGetName(var), lb, ub, primsol, bounddist);

            /* look for "best" row to add (minimal slack), but don't add rows again,
             * that are already involved in aggregation
             */
            nnonzrows = SCIPcolGetNLPNonz(col);
            nonzrows = SCIPcolGetRows(col);
            nonzcoefs = SCIPcolGetVals(col);
            
            for( r = 0; r < nnonzrows; r++ )
            {
               SCIP_Real score;
               SCIP_Real factor;
               SCIP_Real absfactor;
               SCIP_Real activity;
               SCIP_Real lhs;
               SCIP_Real rhs;
               int lppos;
               
               lppos = SCIProwGetLPPos(nonzrows[r]);
               assert(0 <= lppos && lppos < nrows);
               
               SCIPdebugMessage("        -> r=%d row <%s>: weight=%g, pos=%d, alpha_j=%g, a^r_j=%g, factor=%g, %g <= %g <= %g\n",
                  r, SCIProwGetName(nonzrows[r]), rowweights[lppos], lppos, aggrcoefs[c], nonzcoefs[r], 
                  - aggrcoefs[c] / nonzcoefs[r], SCIProwGetLhs(nonzrows[r]), 
                  SCIPgetRowSolActivity(scip, nonzrows[r], sol), SCIProwGetRhs(nonzrows[r]));

               /* take only unmodifiable LP rows, that are not yet aggregated */
               if( rowweights[lppos] != 0.0 || SCIProwIsModifiable(nonzrows[r]) )
                  continue;
               
               /* don't aggregate rows that would lead to a too extreme aggregation factor */
               factor = - aggrcoefs[c] / nonzcoefs[r]; 
               absfactor = REALABS(factor);
               if( !SCIPisPositive(scip, absfactor)
                  || absfactor > sepadata->maxrowfac * minweight
                  || maxweight > sepadata->maxrowfac * absfactor )
                  continue;
               
               /* check, if the row's slack multiplied with the aggregation factor is too large */
               activity = SCIPgetRowSolActivity(scip, nonzrows[r], sol);
               lhs = SCIProwGetLhs(nonzrows[r]);
               rhs = SCIProwGetRhs(nonzrows[r]);
               if( (factor < 0.0 && SCIPisGT(scip, factor * (lhs - activity), maxslack))
                  || (factor > 0.0 && SCIPisGT(scip, factor * (rhs - activity), maxslack)) )
                  continue;
               
               /* for selected real variable y_k, select constraint r with best score SCORE_r with r in P\Q, 
                * where P\Q is the set of constraints not yet involved in the aggregation zet
                */
               assert(!SCIPisInfinity(scip, -SCIProwGetLhs(nonzrows[r])) || rowlhsscores[lppos] == 0.0);
               assert(!SCIPisInfinity(scip, SCIProwGetRhs(nonzrows[r])) || rowrhsscores[lppos] == 0.0);
               score = (factor < 0.0 ? rowlhsscores[lppos] : rowrhsscores[lppos]);
               if( score > 0.0 && (bounddist > bestbounddist + sepadata->aggrtol || score > bestscore) )
               {
                  bestbounddist = bounddist;
                  bestscore = score; 
                  bestcol = col;
                  bestrow = nonzrows[r];
                  aggrfac = factor;
                  SCIPdebugMessage("     -> col <%s>: %g * row <%s>, bounddist=%g, score=%g\n",
                     SCIPvarGetName(SCIPcolGetVar(bestcol)), aggrfac, SCIProwGetName(bestrow), bestbounddist, score);
               }
            }
         }
      }
      assert((bestcol == NULL) == (bestrow == NULL));

      /* abort, if no row can be added to remove an additional active continuous variable */
      if( bestcol == NULL )
      {
         SCIPdebugMessage(" -> abort aggregation: no removable column found\n");
         break;
      }

      /* Step 3: add row to aggregation */
      SCIPdebugMessage(" -> adding %+g<%s> to eliminate variable <%s> (aggregation %d)\n", 
         aggrfac, SCIProwGetName(bestrow), SCIPvarGetName(SCIPcolGetVar(bestcol)), naggrs+1);
      assert(rowweights[SCIProwGetLPPos(bestrow)] == 0.0);
      assert(!SCIPisZero(scip, aggrfac));

      /* change row's aggregation weight */
      rowweights[SCIProwGetLPPos(bestrow)] = aggrfac;
      absaggrfac = REALABS(aggrfac);
      maxweight = MAX(maxweight, absaggrfac);
      minweight = MIN(minweight, absaggrfac);

      /* decrease score of aggregation row in order to not aggregate it again too soon */
      decreaseRowScore(scip, rowlhsscores, rowrhsscores, SCIProwGetLPPos(bestrow));

      /* change coefficients of aggregation and update the number of continuous variables */
      bestrownonzcols = SCIProwGetCols(bestrow);
      bestrownonzcoefs = SCIProwGetVals(bestrow);
      nbestrownonzcols = SCIProwGetNLPNonz(bestrow);
      for( c = 0; c < nbestrownonzcols; c++ )
      {
         SCIP_VAR* var;
         int pos;
         SCIP_Bool iscont;
         SCIP_Bool waszero;
         SCIP_Bool iszero;

         var = SCIPcolGetVar(bestrownonzcols[c]);
         pos = SCIPcolGetLPPos(bestrownonzcols[c]);
         assert(pos >= 0);
         assert(!SCIPisZero(scip, bestrownonzcoefs[c]));

         iscont = (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
         waszero = (aggrcoefs[pos] == 0.0);
         aggrcoefs[pos] += bestrownonzcoefs[c] * aggrfac;
         iszero = SCIPisZero(scip, aggrcoefs[pos]);

         if( iszero )
         {
            aggrcoefs[pos] = 0.0;
            if( !waszero )
            {
               /* coefficient switched from non-zero to zero */
               if( iscont )
               {
                  nzi = aggrnonzidxs[pos];
                  assert(0 <= nzi && nzi < naggrcontnonzs);
                  assert(aggrcontnonzposs[nzi] == pos);
                  aggrcontnonzposs[nzi] = aggrcontnonzposs[naggrcontnonzs-1];
                  aggrnonzidxs[aggrcontnonzposs[nzi]] = nzi;
                  naggrcontnonzs--;
                  updateNActiveConts(scip, varsolvals, bestcontlbs, bestcontubs, nintvars, var, -1, &nactiveconts);
               }
               else
               {
                  nzi = aggrnonzidxs[pos];
                  assert(0 <= nzi && nzi < naggrintnonzs);
                  assert(aggrintnonzposs[nzi] == pos);
                  aggrintnonzposs[nzi] = aggrintnonzposs[naggrintnonzs-1];
                  aggrnonzidxs[aggrintnonzposs[nzi]] = nzi;
                  naggrintnonzs--;
               }
            }
         }
         else if( waszero )
         {
            /* coefficient switched from zero to non-zero */
            if( iscont )
            {
               assert(naggrcontnonzs < ncols);
               aggrnonzidxs[pos] = naggrcontnonzs;
               aggrcontnonzposs[naggrcontnonzs] = pos;
               naggrcontnonzs++;
               updateNActiveConts(scip, varsolvals, bestcontlbs, bestcontubs, nintvars, var, +1, &nactiveconts);
            }
            else
            {
               assert(naggrintnonzs < ncols);
               aggrnonzidxs[pos] = naggrintnonzs;
               aggrintnonzposs[naggrintnonzs] = pos;
               naggrintnonzs++;
            }
         }
      }
      naggrs++;

      SCIPdebugMessage(" -> %d continuous variables left (%d/%d active), %d/%d nonzeroes, %d/%d aggregations\n", 
         naggrcontnonzs, nactiveconts, maxconts, naggrcontnonzs + naggrintnonzs, maxaggrnonzs, naggrs, maxaggrs);
   }
#ifdef SCIP_DEBUG
   if( nactiveconts > maxconts )
      SCIPdebugMessage(" -> abort aggregation: %d/%d active continuous variables\n", nactiveconts, maxconts);
   if( naggrcontnonzs + naggrintnonzs > maxaggrnonzs )
      SCIPdebugMessage(" -> abort aggregation: %d/%d nonzeros\n", naggrcontnonzs + naggrintnonzs, maxaggrnonzs);
#endif

   /* free datastructures */
   SCIPfreeBufferArray(scip, &rowweights);
   SCIPfreeBufferArray(scip, &aggrcontnonzposs);
   SCIPfreeBufferArray(scip, &aggrintnonzposs);
   SCIPfreeBufferArray(scip, &aggrnonzidxs);
   SCIPfreeBufferArray(scip, &aggrcoefs);

   return SCIP_OKAY; 
}

/** searches and adds c-MIR cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPA*            sepa,               /**< the c-MIR separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** vars;
   SCIP_Real* varsolvals;
   SCIP_Real* bestcontlbs;
   SCIP_Real* bestcontubs;
   SCIP_ROW** rows;     
   SCIP_Real* rowlhsscores;
   SCIP_Real* rowrhsscores;
   SCIP_Real* rowscores;
   int* roworder;
   SCIP_Real maxslack;
   SCIP_Real objnorm;
   int nvars;
   int nintvars;
   int ncontvars;
   int nrows;
   int ntries;
   int nfails;
   int depth;
   int ncalls;
   int maxtries;
   int maxfails;
   int maxaggrs;
   int maxsepacuts;
   int maxtestdelta;
   int maxconts;
   int ncuts;
   int r;
   int v;
   char normtype;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the cmir cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* get all rows and number of columns */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) ); 
   assert(nrows == 0 || rows != NULL);

   /* nothing to do, if LP is empty */
   if( nrows == 0 )
      return SCIP_OKAY;

   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   ncontvars = SCIPgetNContVars(scip);
   nintvars = nvars-ncontvars;
   assert(nvars == 0 || vars != NULL);

   /* nothing to do, if problem has no variables */
   if( nvars == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* get data structure */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowlhsscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowrhsscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &roworder, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcontlbs, ncontvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcontubs, ncontvars) );
  
   /* get the solution values for all active variables */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, varsolvals) );

   /* calculate the tightest bounds w.r.t. current solution for the continuous variables */
   for( v = nintvars; v < nvars; ++v )
   {
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Real bestvlb;
      SCIP_Real bestvub;
      int bestvlbidx;
      int bestvubidx;

      if( ALLOWLOCAL )
      {
         bestlb = SCIPvarGetLbLocal(vars[v]);
         bestub = SCIPvarGetUbLocal(vars[v]);
      }
      else
      {
         bestlb = SCIPvarGetLbGlobal(vars[v]);
         bestub = SCIPvarGetUbGlobal(vars[v]);
      }
      SCIPvarGetClosestVlb(vars[v], &bestvlb, &bestvlbidx);
      SCIPvarGetClosestVub(vars[v], &bestvub, &bestvubidx);
      if( bestvlbidx >= 0 )
         bestlb = MAX(bestlb, bestvlb);
      if( bestvubidx >= 0 )
         bestub = MIN(bestub, bestvub);

      bestcontlbs[v-nintvars] = bestlb;
      bestcontubs[v-nintvars] = bestub;
   }

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
   {
      maxtries = sepadata->maxtriesroot;
      maxfails = sepadata->maxfailsroot;
      maxaggrs = sepadata->maxaggrsroot;
      maxsepacuts = sepadata->maxsepacutsroot;
      maxslack = sepadata->maxslackroot;
      maxtestdelta = sepadata->maxtestdeltaroot;
      maxconts = sepadata->maxcontsroot;
   }   
   else
   {
      maxtries = sepadata->maxtries;
      maxfails = sepadata->maxfails;
      maxaggrs = sepadata->maxaggrs;
      maxsepacuts = sepadata->maxsepacuts;
      maxslack = sepadata->maxslack;
      maxtestdelta = sepadata->maxtestdelta;
      maxconts = sepadata->maxconts;
   }

   /* calculate aggregation scores for both sides of all rows, and sort rows by nonincreasing maximal score */
   objnorm = SCIPgetObjNorm(scip);
   objnorm = MAX(objnorm, 1.0);
   for( r = 0; r < nrows; r++ )
   {
      int nnonz;
      int i;

      assert(SCIProwGetLPPos(rows[r]) == r);

      nnonz = SCIProwGetNNonz(rows[r]);
      if( nnonz == 0 )
      {
         /* ignore empty rows */
         rowlhsscores[r] = 0.0;
         rowrhsscores[r] = 0.0;
      }
      else
      {
         SCIP_Real activity;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real dualsol;
         SCIP_Real dualscore;
         SCIP_Real rowdensity;
         SCIP_Real rownorm;
         SCIP_Real slack;

         dualsol = (sol == NULL ? SCIProwGetDualsol(rows[r]) : 1.0);
         activity = SCIPgetRowSolActivity(scip, rows[r], sol);
         lhs = SCIProwGetLhs(rows[r]);
         rhs = SCIProwGetRhs(rows[r]);
         rownorm = SCIProwGetNorm(rows[r]);
         rownorm = MAX(rownorm, 0.1);
         rowdensity = (SCIP_Real)nnonz/(SCIP_Real)nvars;
         assert(SCIPisPositive(scip, rownorm));

         slack = (activity - lhs)/rownorm;
         dualscore = MAX(dualsol/objnorm, 0.0001);
         if( !SCIPisInfinity(scip, -lhs) && SCIPisLE(scip, slack, maxslack)
            && (ALLOWLOCAL || !SCIProwIsLocal(rows[r])) /*lint !e506 !e774*/
            && rowdensity <= sepadata->maxrowdensity
            && rowdensity <= sepadata->maxaggdensity )  /*lint !e774*/
         {
            rowlhsscores[r] = dualscore + sepadata->densityscore * (1.0-rowdensity)
               + sepadata->slackscore * MAX(1.0 - slack, 0.0);
            assert(rowlhsscores[r] > 0.0);
         }
         else
            rowlhsscores[r] = 0.0;

         slack = (rhs - activity)/rownorm;
         dualscore = MAX(-dualsol/objnorm, 0.0001);
         if( !SCIPisInfinity(scip, rhs) && SCIPisLE(scip, slack, maxslack)
            && (ALLOWLOCAL || !SCIProwIsLocal(rows[r])) /*lint !e506 !e774*/
            && rowdensity <= sepadata->maxrowdensity
            && rowdensity <= sepadata->maxaggdensity )  /*lint !e774*/
         {
            rowrhsscores[r] = dualscore + sepadata->densityscore * (1.0-rowdensity)
               + sepadata->slackscore * MAX(1.0 - slack, 0.0);
            assert(rowrhsscores[r] > 0.0);
         }
         else
            rowrhsscores[r] = 0.0;
      }
      rowscores[r] = MAX(rowlhsscores[r], rowrhsscores[r]);
      for( i = r; i > 0 && rowscores[r] > rowscores[roworder[i-1]]; --i )
         roworder[i] = roworder[i-1];
      assert(0 <= i && i <= r);
      roworder[i] = r;
   }

   /* start aggregation heuristic for each row in the LP */
   ncuts = 0;
   if( maxtries < 0 )
      maxtries = INT_MAX;
   if( maxfails < 0 )
      maxfails = INT_MAX;
   else if( depth == 0 && 2*SCIPgetNSepaRounds(scip) < maxfails )
      maxfails += maxfails - 2*SCIPgetNSepaRounds(scip); /* allow up to double as many fails in early separounds of root node */
   ntries = 0;
   nfails = 0;
   for( r = 0; r < nrows && ntries < maxtries && ncuts < maxsepacuts && rowscores[roworder[r]] > 0.0; r++ )
   {
      SCIP_Bool wastried;
      int oldncuts;

      oldncuts = ncuts;
      SCIP_CALL( aggregation(scip, sepadata, sol, varsolvals, bestcontlbs, bestcontubs, rowlhsscores, rowrhsscores, roworder[r], 
            maxaggrs, maxslack, maxtestdelta, maxconts, normtype, &wastried, &ncuts) );
      if( !wastried )
         continue;
      ntries++;
      if( ncuts == oldncuts )
      {
         nfails++;
         if( nfails >= maxfails )
            break;
      }
      else
         nfails = 0;
   }

   /* free data structure */
   SCIPfreeBufferArray(scip, &bestcontubs);
   SCIPfreeBufferArray(scip, &bestcontlbs);
   SCIPfreeBufferArray(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &roworder);
   SCIPfreeBufferArray(scip, &rowscores);
   SCIPfreeBufferArray(scip, &rowrhsscores);
   SCIPfreeBufferArray(scip, &rowlhsscores);

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;
    
   return SCIP_OKAY;
}




/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeCmir)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;

}


/** initialization method of separator (called when problem solving starts) */
#define sepaInitCmir NULL


/** deinitialization method of separator (called when problem solving exits) */
#define sepaExitCmir NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolCmir NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolCmir NULL


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpCmir)
{  /*lint --e{715}*/
   SCIP_CALL( separateCuts(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolCmir)
{  /*lint --e{715}*/
   SCIP_CALL( separateCuts(scip, sepa, sol, result) );

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the cmir separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaCmir(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create cmir separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_DELAY,
         sepaFreeCmir, sepaInitCmir, sepaExitCmir, 
         sepaInitsolCmir, sepaExitsolCmir, 
         sepaExeclpCmir, sepaExecsolCmir,
         sepadata) );

   /* add cmir separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxrounds",
         "maximal number of cmir separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxroundsroot",
         "maximal number of cmir separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxtries",
         "maximal number of rows to start aggregation with per separation round (-1: unlimited)",
         &sepadata->maxtries, DEFAULT_MAXTRIES, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxtriesroot",
         "maximal number of rows to start aggregation with per separation round in the root node (-1: unlimited)",
         &sepadata->maxtriesroot, DEFAULT_MAXTRIESROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxfails",
         "maximal number of consecutive unsuccesful aggregation tries (-1: unlimited)",
         &sepadata->maxfails, DEFAULT_MAXFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxfailsroot",
         "maximal number of consecutive unsuccesful aggregation tries in the root node (-1: unlimited)",
         &sepadata->maxfailsroot, DEFAULT_MAXFAILSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxaggrs",
         "maximal number of aggregations for each row per separation round",
         &sepadata->maxaggrs, DEFAULT_MAXAGGRS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxaggrsroot",
         "maximal number of aggregations for each row per separation round in the root node",
         &sepadata->maxaggrsroot, DEFAULT_MAXAGGRSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxsepacuts",
         "maximal number of cmir cuts separated per separation round",
         &sepadata->maxsepacuts, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxsepacutsroot",
         "maximal number of cmir cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/maxslack",
         "maximal slack of rows to be used in aggregation",
         &sepadata->maxslack, DEFAULT_MAXSLACK, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/maxslackroot",
         "maximal slack of rows to be used in aggregation in the root node",
         &sepadata->maxslackroot, DEFAULT_MAXSLACKROOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/densityscore",
         "weight of row density in the aggregation scoring of the rows",
         &sepadata->densityscore, DEFAULT_DENSITYSCORE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/slackscore",
         "weight of slack in the aggregation scoring of the rows",
         &sepadata->slackscore, DEFAULT_SLACKSCORE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/maxaggdensity",
         "maximal density of aggregated row",
         &sepadata->maxaggdensity, DEFAULT_MAXAGGDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/maxrowdensity",
         "maximal density of row to be used in aggregation",
         &sepadata->maxrowdensity, DEFAULT_MAXROWDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/maxrowfac",
         "maximal row aggregation factor",
         &sepadata->maxrowfac, DEFAULT_MAXROWFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxtestdelta",
         "maximal number of different deltas to try (-1: unlimited)",
         &sepadata->maxtestdelta, DEFAULT_MAXTESTDELTA, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxtestdeltaroot",
         "maximal number of different deltas to try in the root node (-1: unlimited)",
         &sepadata->maxtestdeltaroot, DEFAULT_MAXTESTDELTAROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxconts",
         "maximal number of active continuous variables in aggregated row", 
         &sepadata->maxconts, DEFAULT_MAXCONTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cmir/maxcontsroot",
         "maximal number of active continuous variables in aggregated row in the root node", 
         &sepadata->maxcontsroot, DEFAULT_MAXCONTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cmir/aggrtol",
         "tolerance for bounddistances used to select continuous variable in current aggregated constraint to be eliminated",
         &sepadata->aggrtol, DEFAULT_AGGRTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cmir/trynegscaling",
         "should negative values also be tested in scaling?",
         &sepadata->trynegscaling, DEFAULT_TRYNEGSCALING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cmir/fixintegralrhs",
         "should an additional variable be complemented if f0 = 0?",
         &sepadata->fixintegralrhs, DEFAULT_FIXINTEGRALRHS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cmir/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
