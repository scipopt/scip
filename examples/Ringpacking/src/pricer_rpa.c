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

/**@file   pricer_ringpacking.c
 * @brief  Ringpacking variable pricer
 * @author Benjamin Mueller
 *
 * This file implements the variable pricer which check if variables exist with negative reduced cost. See
 * @ref PRICER for more details.
 *
 * @page PRICER Pricing new variables
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
};


/**@name Local methods
 *
 * @{
 */

/** returns an upper bound on the denity for n equal circles in a square (holds also for rectangles); this is a result
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
      (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", types[i]);
      (void) strcat(name, strtmp);

      SCIP_CALL( SCIPpatternAddElement(pattern, types[i], xs[i], ys[i]) );
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
      assert(types[i] >= 0 && types[i] < SCIPprobdataGetNTypes(probdata));
      SCIP_CALL( SCIPaddCoefLinear(scip, conss[types[i]], var, 1.0) );
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
   SCIP_Bool*            success             /**< pointer to store if we could add at least one variable with reduced costs */
   )
{
   SCIP_SOL* sol;
   SCIP_Real* xs;
   SCIP_Real* ys;
   int* selectedtypes;
   int nselected;
   int i;

   assert(SCIPgetNSols(subscip) > 0);
   assert(success != NULL);

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
   SCIP_CALL( addVariable(scip, probdata, selectedtypes, xs, ys, nselected) );

   /* free memory */
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);
   SCIPfreeBufferArray(scip, &selectedtypes);

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
   assert(timelim >= 0.0);
   assert(nodelim >= -1L);
   assert(addedvar != NULL);
   assert(solstat != NULL);
   assert(dualbound != NULL);

   *addedvar = FALSE;
   *solstat = SCIP_STATUS_UNKNOWN;
   *dualbound = -SCIPinfinity(scip);

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
   SCIPfreeBlockMemoryNull(scip, &pricerdata);

   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitRingpacking)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolRingpacking)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostRingpacking)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_CONS** conss;
   SCIP_Real* lambdas;
   SCIP_STATUS solstat;
   SCIP_Real redcostslb;
   SCIP_Bool success;
   int t;

   *result = SCIP_SUCCESS;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* only run pricer in the root node */
   if( SCIPgetDepth(scip) > 0 )
   {
      SCIPprobdataInvalidateDualbound(scip, probdata);
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

   /* TODO add parameter for time and node limit */
   SCIP_CALL( solvePricingMINLP(scip, probdata, lambdas, 100.0, 10000L, &success, &solstat, &redcostslb) );
   redcostslb += 1.0;
   SCIPdebugMsg(scip, "result of pricing MINLP: addedvar=%u soltat=%d\n", success, solstat);

   /* compute Farley's bound */
   if( SCIPisFeasGE(scip, redcostslb, 0.0) )
   {
      *lowerbound = SCIPgetLPObjval(scip);
      SCIPinfoMessage(scip, NULL, "+++++++++++++ LP(master) = ceil(%g) = %g\n", *lowerbound, SCIPfeasCeil(scip, *lowerbound));
   }
   else
   {
      *lowerbound = SCIPgetLPObjval(scip) / (1.0 - redcostslb);
      SCIPinfoMessage(scip, NULL, "+++++++++++++ Farley's bound = ceil(%g/%g) = %g\n", SCIPgetLPObjval(scip), 1.0 - redcostslb,
         SCIPfeasCeil(scip, *lowerbound));
   }
   *lowerbound = SCIPfeasCeil(scip, *lowerbound);

   /* updates dual bound that is stored in the problem data */
   SCIPprobdataUpdateDualbound(scip, probdata, *lowerbound);

   /* invalidate dual bound if no variable has been added and the pricing problem has not been solved to optimality */
   if( !success && solstat != SCIP_STATUS_OPTIMAL )
   {
      assert(SCIPisFeasLE(scip, redcostslb, 0.0));
      SCIPprobdataInvalidateDualbound(scip, probdata);
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
SCIP_RETCODE SCIPincludePricerRingpacking(
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
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolRingpacking) );

   /* variable pricer parameters */

   return SCIP_OKAY;
}

/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerRingpackingActivate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/**@} */
