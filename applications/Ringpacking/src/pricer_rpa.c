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

/**@file   pricer_rpa.c
 * @brief  Ringpacking variable pricer
 * @author Benjamin Mueller
 *
 * This file implements the variable pricer which checks if variables negative reduced cost exist. See
 * @ref RINGPACKING_PRICER for more details.
 *
 * @page RINGPACKING_PRICER Pricing new variables
 *
 * The task of the pricer is to search for new variables with negative reduced costs. For this, the following non-linear
 * program is solved:
 *
 * \f[
 *  \begin{equation}
 *    \min_{P \in \mathcal{RP}} \left\{1 - \sum_{t \in \mathcal{T}} \lambda_t P_t\right\},
 *  \end{equation}
 * \f]
 *
 * where \f$\lambda\f$ is given by the dual solution of the restricted master problem. See the
 * \ref RINGPACKING_PROBLEM "problem description" for more details.
 *
 * This problem is very hard, but can be modeled as a weighted circle packing problem for a single rectangle. Therefore,
 * we first use a simple greedy heuristic to solve the problem. If the heuristic fails, the MINLP is solved with
 * conventional methods on a new \SCIP instance and a given time limit. If the problem can be solved and the optimal
 * value is non-negative, the LP relaxation has been solved to optimality and what remains is ensuring integrality of
 * the solution by the normal \SCIP framework. If, on the other hand, the best solution found by both methods is negative,
 * we have found an improving pattern, whose corresponding variable needs to be added to the restricted master problem.
 * It is possible (and not unlikely) that neither method succeeds in finding a pattern with negative solution value. In
 * that case, we also exit the pricing loop, just as if we had found an optimal solution, and proceed with enforcing
 * integrality resulting in a feasible primal solution to the whole problem.
 *
 * In case that the pricing problem cannot be solved to optimality, we cannot directly deduce a lower bound for the
 * master problem. However, the following theorem by Farley shows how a valid dual bound can be computed from the
 * LP solution and the pricing solution, see
 * <a href="https://doi.org/10.1287/opre.38.5.922">A Note on Bounding a Class of Linear Programming Problems, Including
 * Cutting Stock Problems</a> for more details.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "scip/scipdefplugins.h"
#include "scip/scip.h"
#include "pricer_rpa.h"
#include "probdata_rpa.h"
#include "pattern.h"

/**@name Pricer properties
 *
 * @{
 */

#define PRICER_NAME            "ringpacking"
#define PRICER_DESC            "pricer for ringpacking"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE          /* only call pricer if all problem variables have non-negative reduced costs */

/* default values of pricing parameters */
#define DEFAULT_PRICING_NLPTILIM             120.0     /**< time limit for each pricing NLP */
#define DEFAULT_PRICING_NLPNODELIM           SCIP_LONGINT_MAX /**< node limit for each pricing NLP */
#define DEFAULT_PRICING_HEURTILIM            60.0      /**< time limit for each heuristic pricing */
#define DEFAULT_PRICING_HEURITERLIM          1000      /**< iteration limit for each heuristic pricing */
#define DEFAULT_PRICING_TOTALTILIM           3600.0    /**< total time limit for all pricing NLPs and heuristic calls */

/**@} */

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif

/*
 * Data structures
 */

/** @brief Variable pricer data used in the \ref pricer_ringpacking.c "pricer" */
struct SCIP_PricerData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Real             timeleft;           /**< time left for solving pricing problems (with NLP or heuristic) */

   /* parameters */
   SCIP_Real             nlptilim;           /**< time limit for pricing NLP */
   SCIP_Real             heurtilim;          /**< time limit for pricing heuristic */
   SCIP_Longint          nlpnodelim;         /**< node limit for pricing NLP */
   int                   heuriterlim;        /**< iteration limit for pricing heuristic */
};


/**@name Local methods
 *
 * @{
 */

/** returns an upper bound on the density for n equal circles in a square (holds also for rectangles); this is a result
 * of Groemer's formula (see, 'Ueber die Einlagerung von Kreisen in einen konvexen Bereich')
 */
static
SCIP_Real getDensityUb(int n)
{
   assert(n > 0);
   if( n == 1 )
      return M_PI / 4.0;
   return (n * M_PI) / SQR(2.0 - SQRT(3.0) + SQRT(7.0 - M_PI + SQRT(3.0)*(2.0*n - 6.0 + M_PI)) );/*lint !e666*/
}

/** helper function to count how many circles of a given type are needed */
static
int getNCircles(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             rext,               /**< external radius */
   int                   demand,             /**< demand */
   SCIP_Real             width,              /**< width of the rectangle */
   SCIP_Real             height,             /**< height of the rectangle */
   SCIP_Real             lambda              /**< objective coefficient of each circle of the given type */
   )
{
   SCIP_Real volrect;
   SCIP_Real volcircle;
   int ncircles;

   assert(!SCIPisFeasNegative(scip, lambda));

   /* objective coefficient of this circle type is zero -> not needed */
   if( !SCIPisFeasPositive(scip, lambda) )
      return 0;

   volrect =  width * height;
   volcircle = M_PI * SQR(rext);
   assert(volcircle != 0.0);

   ncircles = (int)SCIPfeasFloor(scip, (getDensityUb(demand) * volrect) / volcircle);

   /* special cases where too big circles have a minimum distance to each other (in x direction) */
   if( SCIPisGT(scip, rext, height / 4.0) )
   {
      SCIP_Real c = SQRT(4.0 * rext * height - SQR(height));
      ncircles = (int)MIN(ncircles, 1 + (int)SCIPfloor(scip, (width - 2.0*rext) / c)); /*lint !e666*/
   }
   if( SCIPisGT(scip, rext, height / 6.0) && SCIPisLE(scip, rext, height / 4.0) )
   {
      SCIP_Real c = MIN(SQRT(3.0*rext*rext + rext * height - height * height / 4.0),
         SQRT(8.0 * rext * height - height * height - 12.0 * rext * rext)); /*lint !e666*/
      SCIP_Real k = SCIPfloor(scip, height / (2.0 * rext)) + 1;
      SCIP_Real l = (width - 2.0 * rext) / c;
      ncircles = (int)MIN(ncircles, k + l*(k-1) - 1);
   }
   assert(ncircles > 0);

   return MIN(ncircles, demand);
}

/** adds a variable to the restricted master problem */
static
SCIP_RETCODE addVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  types,              /**< types of elements */
   SCIP_Real*            xs,                 /**< x coordinate of each element */
   SCIP_Real*            ys,                 /**< y coordinate of each element */
   SCIP_Bool*            ispacked,           /**< checks whether an element has been packed (might be NULL) */
   int                   nelems              /**< total number of elements */
   )
{
   SCIP_CONS** conss;
   SCIP_PATTERN* pattern;
   SCIP_VAR* var;
   char name[SCIP_MAXSTRLEN];
   char strtmp[SCIP_MAXSTRLEN];
   int i;

   assert(probdata != NULL);
   assert(types != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(nelems > 0);

   conss = SCIPprobdataGetPatternConss(probdata);
   assert(conss != NULL);

   /* create rectangular pattern */
   SCIP_CALL( SCIPpatternCreateRectangular(scip, &pattern) );

   /* create variable name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "r");
   for( i = 0; i < nelems; ++i )
   {
      if( ispacked == NULL || ispacked[i] )
      {
         (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", types[i]);
         (void) strcat(name, strtmp);
         SCIP_CALL( SCIPpatternAddElement(pattern, types[i], xs[i], ys[i]) );
      }
   }

   /* mark pattern to be packable */
   SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

   /* create and add variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPvarSetInitial(var, TRUE) );
   SCIP_CALL( SCIPvarSetRemovable(var, TRUE) );
   SCIPvarMarkDeletable(var);
   SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
   SCIPdebugMsg(scip, "added variable %s\n", name);

   /* add coefficients */
   for( i = 0; i < nelems; ++i )
   {
      if( ispacked == NULL || ispacked[i] )
      {
         assert(types[i] >= 0 && types[i] < SCIPprobdataGetNTypes(probdata));
         SCIP_CALL( SCIPaddCoefLinear(scip, conss[types[i]], var, 1.0) );
      }
   }

   /* add pattern and variable to the problem data */
   SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, var) );

   /* release memory */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );
   SCIPpatternRelease(scip, &pattern);

   return SCIP_OKAY;
}

/* extracts and adds a variable with (hopefully) negative reduced costs */
static
SCIP_RETCODE extractVariablesMINLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_VAR**            x,                  /**< x variables of sub-SCIP */
   SCIP_VAR**            y,                  /**< y variables of sub-SCIP */
   SCIP_VAR**            w,                  /**< w variables of sub-SCIP */
   int*                  types,              /**< type corresponding to each variable */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            success             /**< pointer to store if we could add at least one variable with negative reduced costs */
   )
{
   SCIP_SOL* sol;
   SCIP_Real* xs;
   SCIP_Real* ys;
   int* selectedtypes;
   int nselected;
   int i;

   assert(success != NULL);

   if( SCIPgetNSols(subscip) == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   sol = SCIPgetBestSol(subscip);
   assert(sol != NULL);
   SCIPdebugMsg(scip, "found column with reduced cost = %f\n", 1.0 + SCIPgetSolOrigObj(subscip, sol));

   /* reduced cost is non-negative */
   if( SCIPisFeasGE(subscip, 1.0 + SCIPgetSolOrigObj(subscip, sol), 0.0) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
      *success = TRUE;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedtypes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, nvars) );

   /* scan which circles have been selected */
   nselected = 0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real solval = SCIPgetSolVal(subscip, sol, w[i]);

      if( solval >= 0.5 )
      {
         selectedtypes[nselected] = types[i];
         xs[nselected] = SCIPgetSolVal(subscip, sol, x[i]);
         ys[nselected] = SCIPgetSolVal(subscip, sol, y[i]);
         ++nselected;
      }
   }
   assert(nselected > 0); /* otherwise the reduced cost can not be negative */

   /* add variable to main SCIP */
   SCIP_CALL( addVariable(scip, probdata, selectedtypes, xs, ys, NULL, nselected) );

   /* free memory */
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);
   SCIPfreeBufferArray(scip, &selectedtypes);

   return SCIP_OKAY;
}

/** array to compute the score of each element */
static
void computeScores(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            rexts,              /**< external radii for each type */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   int*                  elements,           /**< type of each element */
   int                   nelements,          /**< total number of elements */
   SCIP_Real*            lambdas,            /**< dual multipliers for each type */
   SCIP_Real*            scores,             /**< array to store the score of each element */
   int                   iter,               /**< iteration round */
   int                   iterlim             /**< total iteration limit */
   )
{
   int i;

   assert(iter < iterlim);

   for( i = 0; i < nelements; ++i )
   {
      int elemtype = elements[i];
      SCIP_Real rext = rexts[elemtype];

      /* use items with largest multipliers first */
      if( iter == 0 )
         scores[i] = lambdas[elemtype];

      /* use largest elements first */
      else if( iter == 1 )
         scores[i] = rext;

      /* use smallest elements first */
      else if( iter == 2 )
         scores[i] = -rext;

      /* use [0,1] * radius^2 */
      else if( iter <= iterlim * 0.1 )
         scores[i] = SCIPrandomGetReal(randnumgen, 0.0, 1.0) * rext * rext;

      /* use [0,1] * radius * lambda */
      else if( iter <= iterlim * 0.4 )
         scores[i] = SCIPrandomGetReal(randnumgen, 0.0, 1.0) * rext * lambdas[elemtype];

      /* use [-1,1] *  lambda / radius */
      else if( iter <= iterlim * 0.8 )
         scores[i] = SCIPrandomGetReal(randnumgen, -1.0, 1.0) * rext * lambdas[elemtype];

      /* use a random order */
      else
         scores[i] = SCIPrandomGetReal(randnumgen, 0.0, 1.0);
   }
}

/** tries to find a column with negative reduced cost by using a greedy packing heuristic */
static
SCIP_RETCODE solvePricingHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PRICERDATA*      pricerdata,         /**< pricer data */
   SCIP_Real*            lambdas,            /**< dual multipliers for pattern constraints */
   SCIP_Real             timelim,            /**< time limit */
   int                   iterlim,            /**< iteration limit */
   SCIP_Bool*            addedvar            /**< pointer to store whether a variable with negative reduced cost has been added */
   )
{
   SCIP_Real* scores;
   SCIP_Real* rexts;
   SCIP_Real* xs;
   SCIP_Real* ys;
   SCIP_Bool* ispacked;
   int* demands;
   int* elements;
   SCIP_Real width;
   SCIP_Real height;
   SCIP_Real timestart;
   SCIP_Real bestredcosts;
   SCIP_Real bestvol;
   int nelements;
   int ntypes;
   int npacked;
   int niters;
   int t;

   assert(pricerdata != NULL);
   assert(addedvar != NULL);

   *addedvar = FALSE;
   niters = 0;
   timestart = SCIPgetTotalTime(scip);
   bestredcosts = 0.0;
   bestvol = 0.0;

   /* get problem data */
   rexts = SCIPprobdataGetRexts(probdata);
   demands = SCIPprobdataGetDemands(probdata);
   width = SCIPprobdataGetWidth(probdata);
   height = SCIPprobdataGetHeight(probdata);
   ntypes = SCIPprobdataGetNTypes(probdata);

   /* compute the total number of elements that need to be considered */
   nelements = 0;
   for( t = 0; t < ntypes; ++t )
      nelements += getNCircles(scip, rexts[t], demands[t], width, height, lambdas[t]);

   /* allocate enough memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &elements, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ispacked, nelements) );

   /* create entry for each element */
   nelements = 0;
   for( t = 0; t < ntypes; ++t )
   {
      int i;
      int n;

      n = getNCircles(scip, rexts[t], demands[t], width, height, lambdas[t]);

      for( i = 0; i < n; ++i )
      {
         elements[nelements] = t;
         ++nelements;
      }
   }

   /* main loop */
   while( niters < iterlim
      && SCIPgetTotalTime(scip) - timestart <= timelim
      && !SCIPisStopped(scip) )
   {
      SCIP_Real redcosts = 1.0;
      SCIP_Real vol = 0.0;
      int i;

      /* compute scores */
      computeScores(scip, rexts, pricerdata->randnumgen, elements, nelements, lambdas, scores, niters, iterlim);

      /* sort elements in non-increasing order */
      SCIPsortDownRealInt(scores, elements, nelements);

      /* call heuristic */
      SCIPpackCirclesGreedy(scip, rexts, xs, ys, -1.0, width, height, ispacked, elements, nelements,
         SCIP_PATTERNTYPE_RECTANGULAR, &npacked, niters);

      /* compute reduced costs */
      for( i = 0; i < nelements; ++i )
      {
         if( ispacked[i] )
         {
            redcosts -= lambdas[elements[i]];
            vol += rexts[elements[i]];
         }
      }

      /* add pattern if reduced costs are negative */
      if( SCIPisFeasLT(scip, redcosts, bestredcosts) || (SCIPisGT(scip, vol, bestvol) && SCIPisFeasNegative(scip, redcosts)) )
      {
         SCIPdebugMsg(scip, "pricing heuristic found column with reduced costs %g and volume %g after %d iterations\n", redcosts, vol, niters + 1);

         SCIP_CALL( addVariable(scip, probdata, elements, xs, ys, ispacked, nelements) );
         *addedvar = TRUE;
         bestredcosts = MIN(bestredcosts, redcosts);
         bestvol = MAX(bestvol, vol);
      }

      ++niters;
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &ispacked);
   SCIPfreeBufferArray(scip, &scores);
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);
   SCIPfreeBufferArray(scip, &elements);

   return SCIP_OKAY;
}

/** auxiliary function for solving the pricing problem exactly */
static
SCIP_RETCODE solvePricingMINLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real*            lambdas,            /**< dual multipliers for pattern constraints */
   SCIP_Real             timelim,            /**< time limit */
   SCIP_Longint          nodelim,            /**< node limit */
   SCIP_Bool*            addedvar,           /**< pointer to store whether a variable with negative reduced cost has been added */
   SCIP_STATUS*          solstat,            /**< pointer to store the solution status */
   SCIP_Real*            dualbound           /**< pointer to store the dual bound */
   )
{
   SCIP* subscip;
   SCIP_VAR** x;
   SCIP_VAR** y;
   SCIP_VAR** w;
   SCIP_VAR* quadvars1[6];
   SCIP_VAR* quadvars2[6];
   SCIP_VAR* linvars[2];
   SCIP_Real quadcoefs[6];
   SCIP_Real lincoefs[2];
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   SCIP_Real* rexts;
   SCIP_Real* vols;
   int* types;
   int* demands;
   SCIP_Real width;
   SCIP_Real height;
   int nvars;
   int ntypes;
   int pos;
   int t;
   int i;

   assert(probdata != NULL);
   assert(lambdas != NULL);
   assert(nodelim >= -1L);
   assert(addedvar != NULL);
   assert(solstat != NULL);
   assert(dualbound != NULL);

   *addedvar = FALSE;
   *solstat = SCIP_STATUS_UNKNOWN;
   *dualbound = -SCIPinfinity(scip);

   /* no time left */
   if( timelim <= 0.0 )
      return SCIP_OKAY;

   width = SCIPprobdataGetWidth(probdata);
   height = SCIPprobdataGetHeight(probdata);
   ntypes = SCIPprobdataGetNTypes(probdata);
   rexts = SCIPprobdataGetRexts(probdata);
   demands = SCIPprobdataGetDemands(probdata);
   assert(ntypes > 0);

   /* create a sub-SCIP */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPcreateProbBasic(subscip, "pricing problem") );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* set heuristics to aggressive */
   SCIP_CALL( SCIPsetHeuristics(subscip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/mpec/freq", -1) );

#ifndef SCIP_DEBUG
   SCIPsetMessagehdlrQuiet(subscip, TRUE);
#endif

   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelim) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelim) );

   /* count how many variables are needed */
   nvars = 0;
   for( t = 0; t < ntypes; ++t )
   {
      nvars += getNCircles(scip, rexts[t], demands[t], width, height, lambdas[t]);
      SCIPdebugMsg(scip, "use %d/%d circles of type %d\n", getNCircles(scip, rexts[t], demands[t], width, height, lambdas[t]), demands[t], t);
   }
   assert(nvars > 0);

   /* allocate enough memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &types, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &vols, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &x, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &y, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &w, nvars) );

   /* create variables */
   pos = 0;
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_Real obj = lambdas[t];
      int ncircles = getNCircles(scip, rexts[t], demands[t], width, height, lambdas[t]);
      int k;

      for( k = 0; k < ncircles; ++k )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", pos);
         SCIP_CALL( SCIPcreateVarBasic(subscip, &x[pos], name, rexts[t], width - rexts[t], 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(subscip, x[pos]) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", pos);
         SCIP_CALL( SCIPcreateVarBasic(subscip, &y[pos], name, rexts[t], height - rexts[t], 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(subscip, y[pos]) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "w_%d", pos);
         SCIP_CALL( SCIPcreateVarBasic(subscip, &w[pos], name, 0.0, 1.0, -obj, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(subscip, w[pos]) );

         vols[pos] = SQR(rexts[t]) * M_PI;
         types[pos] = t;
         ++pos;
      }
   }

   /* create non-overlapping constraints */
   for( i = 0; i < nvars; ++i )
   {
      int j;
      int type1;

      type1 = types[i];
      assert(type1 >= 0 && type1 < ntypes);

      for( j = i+1; j < nvars; ++j )
      {
         SCIP_Real c;
         int types2;

         types2 = types[j];
         assert(types2 >= 0 && types2 < ntypes);

         c = (rexts[type1] + rexts[types2]) * (rexts[type1] + rexts[types2]);

         /* linear part */
         linvars[0] = w[i]; lincoefs[0] = -c;
         linvars[1] = w[j]; lincoefs[1] = -c;

         /* quadratic part */
         quadvars1[0] = x[i]; quadvars2[0] = x[i]; quadcoefs[0] =  1.0;
         quadvars1[1] = x[i]; quadvars2[1] = x[j]; quadcoefs[1] = -2.0;
         quadvars1[2] = x[j]; quadvars2[2] = x[j]; quadcoefs[2] =  1.0;
         quadvars1[3] = y[i]; quadvars2[3] = y[i]; quadcoefs[3] =  1.0;
         quadvars1[4] = y[i]; quadvars2[4] = y[j]; quadcoefs[4] = -2.0;
         quadvars1[5] = y[j]; quadvars2[5] = y[j]; quadcoefs[5] =  1.0;

         /* create quadratic constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nonoverlap_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicQuadratic(subscip, &cons, name, 2, linvars, lincoefs, 6, quadvars1, quadvars2,
               quadcoefs, -c, SCIPinfinity(subscip)) );

         /* add and release constraint */
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      }
   }

   /* w_i >= w_{i+1} if type(i) == type(i+1) */
   for( i = 0; i < nvars - 1; ++i )
   {
      int type1;
      int type2;

      type1 = types[i];
      type2 = types[i+1];

      if( type1 != type2 )
         continue;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons_w_%d_%d", i, i+1);

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, name, 0, NULL, NULL,
            0.0, SCIPinfinity(subscip)) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, cons, w[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, cons, w[i+1], -1.0) );

      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }

   /* x_i <= x_{i+1} if type(i) == type(i+1) */
   for( i = 0; i < nvars - 1; ++i )
   {
      int type1;
      int type2;

      type1 = types[i];
      type2 = types[i+1];

      if( type1 != type2 )
         continue;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symmetry_%d_%d", i, i+1);
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, name, 0, NULL, NULL,
            0.0, SCIPinfinity(subscip)) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, cons, x[i], -1.0) );
      SCIP_CALL( SCIPaddCoefLinear(subscip, cons, x[i+1], 1.0) );

      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }

   /* sum_{i} vol_i w_i <= W*H */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "volume");
   SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, name, nvars, w, vols, 0.0, width * height) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* solve the pricing problem */
   SCIPdebugMsg(scip, "----------------------- solve pricing problem -----------------------\n");
   SCIP_CALL( SCIPsolve(subscip) );
   SCIPdebugMsg(scip, "---------------------------------------------------------------------\n");

   /* check solution status */
   *dualbound = SCIPgetDualbound(subscip);
   *solstat = SCIPgetStatus(subscip);

   /* add variable with negative reduced costs */
   SCIP_CALL( extractVariablesMINLP(scip, probdata, subscip, x, y, w, types, nvars, addedvar) );

   /* free variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &w[i]) );
      SCIP_CALL( SCIPreleaseVar(subscip, &y[i]) );
      SCIP_CALL( SCIPreleaseVar(subscip, &x[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(subscip, &w, nvars);
   SCIPfreeBlockMemoryArray(subscip, &y, nvars);
   SCIPfreeBlockMemoryArray(subscip, &x, nvars);
   SCIPfreeBlockMemoryArray(subscip, &vols, nvars);
   SCIPfreeBlockMemoryArray(subscip, &types, nvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/**@} */

/**name Callback methods
 *
 * @{
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeRingpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata->randnumgen == NULL);

   SCIPfreeBlockMemoryNull(scip, &pricerdata);

   return SCIP_OKAY;
}

/** initialization method of variable pricer (called after problem was transformed and pricer is active) */
static
SCIP_DECL_PRICERINIT(pricerInitRingpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->randnumgen == NULL);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &pricerdata->randnumgen, 0, TRUE) );

   return SCIP_OKAY;
}

/** deinitialization method of variable pricer (called before transformed problem is freed and pricer is active) */
static
SCIP_DECL_PRICEREXIT(pricerExitRingpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->randnumgen != NULL);

   SCIPfreeRandom(scip, &pricerdata->randnumgen);

   return SCIP_OKAY;
}

/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostRingpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;
   SCIP_PROBDATA* probdata;
   SCIP_CONS** conss;
   SCIP_Real* lambdas;
   SCIP_STATUS solstat;
   SCIP_Real redcostslb;
   SCIP_Real nlptimelimit;
   SCIP_Real heurtimelimit;
   SCIP_Real totaltilim;
   SCIP_Bool success;
   int t;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* switch to price-and-price algorithm when dual bound has become invalid */
   *result = SCIPprobdataIsDualboundInvalid(probdata) ? SCIP_SUCCESS : SCIP_DIDNOTRUN;

   /* only run pricer in the root node */
   if( SCIPgetDepth(scip) > 0 )
   {
      SCIPprobdataInvalidateDualbound(scip, probdata);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   conss = SCIPprobdataGetPatternConss(probdata);
   assert(conss != NULL);

   /* collect dual multipliers */
   SCIP_CALL( SCIPallocBufferArray(scip, &lambdas, SCIPprobdataGetNTypes(probdata)) );
   for( t = 0; t < SCIPprobdataGetNTypes(probdata); ++t )
   {
      assert(conss[t] != NULL);
      assert( !strncmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[t])), "linear", 6) );

      lambdas[t] = SCIPgetDualsolLinear(scip, conss[t]);
      SCIPdebugMsg(scip, "lambda_%d = %g\n", t, lambdas[t]);
   }

   /* collect working limits */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &totaltilim) );

   /* solve pricing problem with heuristic */
   heurtimelimit = MIN(pricerdata->heurtilim, totaltilim - SCIPgetSolvingTime(scip)); /*lint !e666*/
   pricerdata->timeleft += SCIPgetSolvingTime(scip);
   SCIP_CALL( solvePricingHeuristic(scip, probdata, pricerdata, lambdas, heurtimelimit, pricerdata->heuriterlim, &success) );
   pricerdata->timeleft -= SCIPgetSolvingTime(scip);

   if( success )
   {
      *result = SCIP_SUCCESS;
   }
   /* solve pricing problem as MINLP if heuristic was not successful and dual bound is still valid */
   else if ( !SCIPprobdataIsDualboundInvalid(probdata) )
   {
      nlptimelimit = MIN3(pricerdata->timeleft, pricerdata->nlptilim, totaltilim - SCIPgetSolvingTime(scip)); /*lint !e666*/
      pricerdata->timeleft += SCIPgetSolvingTime(scip);
      SCIP_CALL( solvePricingMINLP(scip, probdata, lambdas, nlptimelimit, pricerdata->nlpnodelim, &success, &solstat,
            &redcostslb) );
      pricerdata->timeleft -= SCIPgetSolvingTime(scip);
      redcostslb += 1.0;
      SCIPdebugMsg(scip, "result of pricing MINLP: addedvar=%u soltat=%d\n", success, solstat);

      /* check whether pricing problem could be solved to optimality */
      if( SCIPisFeasGE(scip, redcostslb, 0.0) )
      {
         *lowerbound = SCIPgetLPObjval(scip);
         SCIPinfoMessage(scip, NULL, "+++++++++++++ LP(master) = ceil(%g) = %g\n", *lowerbound, SCIPfeasCeil(scip, *lowerbound));
      }
      else
      {
         /* compute Farley's bound */
         *lowerbound = SCIPgetLPObjval(scip) / (1.0 - redcostslb);
         SCIPinfoMessage(scip, NULL, "+++++++++++++ Farley's bound = ceil(%g/%g) = %g\n", SCIPgetLPObjval(scip), 1.0 - redcostslb,
            SCIPfeasCeil(scip, *lowerbound));
      }
      *lowerbound = SCIPfeasCeil(scip, *lowerbound);

      /* updates dual bound that is stored in the problem data */
      SCIPprobdataUpdateDualbound(scip, probdata, *lowerbound);

      /* MINLP found an improving column or pricing problem could have been solved to optimality */
      if( success || solstat == SCIP_STATUS_OPTIMAL || SCIPisFeasGE(scip, redcostslb, 0.0) )
         *result = SCIP_SUCCESS;
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &lambdas);

   return SCIP_OKAY;
}

/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasRingpacking)
{  /*lint --e{715}*/

   /* farkas pricing should not happen */
   SCIPABORT();

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates the ringpacking variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerRpa(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create ringpacking variable pricer data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );
   BMSclearMemory(pricerdata);

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostRingpacking, pricerFarkasRingpacking, pricerdata) );

   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeRingpacking) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitRingpacking) );
   SCIP_CALL( SCIPsetPricerExit(scip, pricer, pricerExitRingpacking) );

   /* variable pricer parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/pricing/nlptilim",
         "time limit for each pricing NLP",
         &pricerdata->nlptilim, FALSE, DEFAULT_PRICING_NLPTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "ringpacking/pricing/nlpnodelim",
         "node limit for each pricing NLP",
         &pricerdata->nlpnodelim, FALSE, DEFAULT_PRICING_NLPNODELIM, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/pricing/heurtilim",
         "time limit for each heuristic pricing",
         &pricerdata->heurtilim, FALSE, DEFAULT_PRICING_HEURTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "ringpacking/pricing/heuriterlim",
         "iteration limit for each heuristic pricing",
         &pricerdata->heuriterlim, FALSE, DEFAULT_PRICING_HEURITERLIM, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/pricing/totaltilim",
         "total time limit for all pricing NLPs and heuristic calls",
         &pricerdata->timeleft, FALSE, DEFAULT_PRICING_TOTALTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerRpaActivate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/**@} */
