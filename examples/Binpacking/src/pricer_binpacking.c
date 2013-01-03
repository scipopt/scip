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

/**@file   pricer_binpacking.c
 * @brief  Binpacking variable pricer
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the variable pricer which check if variables exist with negative reduced cost. See
 * @ref PRICER for more details.
 *
 * @page PRICER Pricing new variables
 *
 * The task of the pricer is to search for new variables with negative reduced costs. For this, the following integer
 * program is solved:
 *
 *  \f[
 *  \begin{array}[t]{rll}
 *       \max & \displaystyle \sum_{i=1}^n (\lambda_S)_i y^\star_i\\
 *        & \\
 *        subject \ to & \displaystyle \sum_{i=0}^n (\lambda_S)_i s_i \leq \kappa \\
 *        & \\
 *        & (\lambda_S)_i \in \{0,1\} & \quad \forall i \in \{ 1, \dots , n \} \\
 *  \end{array}
 * \f]
 *
 * where \f$ (\lambda_S)_i \f$ for \f$i\in\{1,\dots,n\}\f$ are binary variables and \f$y^\star_i\f$ given by the dual
 * solution of the restricted master problem. See the \ref PROBLEM "problem description" for more details.
 *
 * To solve the above integer program, we create a new SCIP instance within SCIP and use the usual functions to create
 * variables and constraints. Besides, we need the current dual solutions to all set covering constraints (each stands
 * for one item) which are the objective coefficients of the binary variables. Therefore, we use the function
 * SCIPgetDualsolSetppc() which returns the dual solutions for the given set covering constraint.
 *
 * Since we also want to generate new variables during search, we have to care that we do not generate variables over
 * and over again. For example, if we branched or fixed a certain packing to zero, we have to make sure that we do not
 * generate the corresponding variables at that node again. For this, we have to add constraints forbidding to generate
 * variables which are locally fixed to zero. See the function addFixedVarsConss() for more details. While using the
 * \ref BRANCHING "Ryan/Foster branching", we also have to ensure that these branching decisions are respected. This is
 * realized within the function addBranchingDecisionConss().
 *
 * @note In case of this binpacking example, the master LP should not get infeasible after branching, because of the way
 *       branching is performed. Therefore, the Farkas pricing is not implemented.
 *       1. In case of Ryan/Foster branching, the two items are selected in a way such that the sum of the LP values of
 *          all columns/packings containing both items is fractional. Hence, it exists at least one column/packing which
 *          contains both items and also at least one column/packing for each item containing this but not the other
 *          item. That means, branching in the "same" direction stays LP feasible since there exists at least one
 *          column/packing with both items and branching in the "differ" direction stays LP feasible since there exists
 *          at least one column/packing containing one item, but not the other.
 *       2. In case of variable branching, we only branch on fractional variables. If a variable is fixed to one, there
 *          is no issue.  If a variable is fixed to zero, then we know that for each item which is part of that
 *          column/packing, there exists at least one other column/packing containing this particular item due to the
 *          covering constraints.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_knapsack.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/scipdefplugins.h"

#include "cons_samediff.h"
#include "pricer_binpacking.h"
#include "probdata_binpacking.h"
#include "vardata_binpacking.h"

/**@name Pricer properties
 *
 * @{
 */

#define PRICER_NAME            "binpacking"
#define PRICER_DESC            "pricer for binpacking tours"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

/**@} */


/*
 * Data structures
 */

/** @brief Variable pricer data used in the \ref pricer_binpacking.c "pricer" */
struct SCIP_PricerData
{
   SCIP_CONSHDLR*        conshdlr;           /**< comstraint handler for "same" and "diff" constraints */
   SCIP_CONS**           conss;              /**< set covering constraints for the items */
   SCIP_Longint*         weights;            /**< weight of the items */
   int*                  ids;                /**< array of item ids */
   int                   nitems;             /**< number of items to be packed */
   SCIP_Longint          capacity;           /**< capacity of the bins */
};



/**@name Local methods
 *
 * @{
 */

/** add branching decisions constraints to the sub SCIP */
static
SCIP_RETCODE addBranchingDecisionConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array of the subscuip oder variables */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for branching data */
   )
{
   SCIP_CONS** conss;
   SCIP_CONS* cons;
   int nconss;
   int id1;
   int id2;
   CONSTYPE type;

   SCIP_Real vbdcoef;
   SCIP_Real lhs;
   SCIP_Real rhs;

   int c;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( conshdlr != NULL );

   /* collect all branching decision constraints */
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   /* loop over all branching decision constraints and apply the branching decision if the corresponding constraint is
    * active
    */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];

      /* ignore constraints which are not active since these are not laying on the current active path of the search
       * tree
       */
      if( !SCIPconsIsActive(cons) )
         continue;

      /* collect the two item ids and the branching type (SAME or DIFFER) on which the constraint branched */
      id1 = SCIPgetItemid1Samediff(scip, cons);
      id2 = SCIPgetItemid2Samediff(scip, cons);
      type = SCIPgetTypeSamediff(scip, cons);

      SCIPdebugMessage("create varbound for %s(%d,%d)\n", type == SAME ? "same" : "diff",
         SCIPprobdataGetIds(SCIPgetProbData(scip))[id1], SCIPprobdataGetIds(SCIPgetProbData(scip))[id2]);

      /* depending on the branching type select the correct left and right hand side for the linear constraint which
       * enforces this branching decision in the pricing problem MIP
       */
      if( type == SAME )
      {
         lhs = 0.0;
         rhs = 0.0;
         vbdcoef = -1.0;
      }
      else if( type == DIFFER )
      {
         lhs = -SCIPinfinity(scip);
         rhs = 1.0;
         vbdcoef = 1.0;
      }
      else
      {
         SCIPerrorMessage("unknow constraint type <%d>\n, type");
         return SCIP_INVALIDDATA;
      }

      /* add linear (in that case a variable bound) constraint to pricing MIP depending on the branching type:
       *
       * - branching type SAME:  x1 = x2 <=> x1 - x2 = 0 <=> 0 <= x1 - x2 <= 0
       *
       * - branching type DIFFER:  x1 - x2 <= 1 <=> -inf <= x1 - x2 <= 1
       *
       */
      SCIP_CALL( SCIPcreateConsBasicVarbound(subscip, &cons, SCIPconsGetName(conss[c]),
            vars[id1], vars[id2], vbdcoef, lhs, rhs) );
      
      SCIPdebug( SCIPprintCons(subscip, cons, NULL) );

      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }

   return SCIP_OKAY;
}

/** avoid to generate columns which are fixed to zero; therefore add for each variable which is fixed to zero a
 *  corresponding logicor constraint to forbid this column
 *
 * @note variable which are fixed locally to zero should not be generated again by the pricing MIP
 */
static
SCIP_RETCODE addFixedVarsConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array of the subscuip */
   SCIP_CONS**           conss,              /**< array of setppc constraint for each item one */
   int                   nitems              /**< number of items */
   )
{
   SCIP_VAR** origvars;
   int norigvars;

   SCIP_CONS* cons;
   int* consids;
   int nconsids;
   int consid;
   int nvars;

   SCIP_VAR** logicorvars;
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   SCIP_Bool needed;
   int nlogicorvars;

   int v;
   int c;
   int o;

   /* collect all variable which are currently existing */
   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);

   /* loop over all these variables and check if they are fixed to zero */
   for( v = 0; v < norigvars; ++v )
   {
      assert(SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_BINARY);

      /* if the upper bound is smaller than 0.5 if follows due to the integrality that the binary variable is fixed to zero */
      if( SCIPvarGetUbLocal(origvars[v]) < 0.5 )
      {
         SCIPdebugMessage("variable <%s> glb=[%.15g,%.15g] loc=[%.15g,%.15g] is fixed to zero\n",
            SCIPvarGetName(origvars[v]), SCIPvarGetLbGlobal(origvars[v]), SCIPvarGetUbGlobal(origvars[v]),
            SCIPvarGetLbLocal(origvars[v]), SCIPvarGetUbLocal(origvars[v]) );

         /* coolect the constraints/items the variable belongs to */
         vardata = SCIPvarGetData(origvars[v]);
         nconsids = SCIPvardataGetNConsids(vardata);
         consids = SCIPvardataGetConsids(vardata);
         needed = TRUE;

         SCIP_CALL( SCIPallocBufferArray(subscip, &logicorvars, nitems) );
         nlogicorvars = 0;
         consid = consids[0];
         nvars = 0;

         /* loop over these items and create a linear (logicor) constraint which forbids this item combination in the
          * pricing problem; thereby check if this item combination is already forbidden
          */
         for( c = 0, o = 0; o < nitems && needed; ++o )
         {
            assert(o <= consid);
            cons = conss[o];

            if( SCIPconsIsEnabled(cons) )
            {
               assert( SCIPgetNFixedonesSetppc(scip, cons) == 0 );

               var = vars[nvars];
               nvars++;
               assert(var != NULL);

               if( o == consid )
               {
                  SCIP_CALL( SCIPgetNegatedVar(subscip, var, &var) );
               }

               logicorvars[nlogicorvars] = var;
               nlogicorvars++;
            }
            else if( o == consid )
               needed = FALSE;

            if( o == consid )
            {
               c++;
               if ( c == nconsids )
                  consid = nitems + 100;
               else
               {
                  assert(consid < consids[c]);
                  consid = consids[c];
               }
            }
         }

         if( needed )
         {
            SCIP_CALL( SCIPcreateConsBasicLogicor(subscip, &cons, SCIPvarGetName(origvars[v]), nlogicorvars, logicorvars) );
            SCIP_CALL( SCIPsetConsInitial(subscip, cons, FALSE) );

            SCIP_CALL( SCIPaddCons(subscip, cons) );
            SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         }

         SCIPfreeBufferArray(subscip, &logicorvars);
      }
   }

   return SCIP_OKAY;
}

/** initializes the pricing problem for the given capacity */
static
SCIP_RETCODE initPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricer data */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars                /**< variable array for the items */
   )
{
   SCIP_CONS** conss;
   SCIP_Longint* vals;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Real dual;

   int nitems;
   int nvars;
   int c;

   assert( SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM );
   assert(pricerdata != NULL);

   nitems = pricerdata->nitems;
   conss = pricerdata->conss;
   weights = pricerdata->weights;
   capacity = pricerdata->capacity;
   nvars = 0;

   SCIP_CALL( SCIPallocBufferArray(subscip, &vals, nitems) );

   /* create for each order, which is not assigned yet, a variable with objective coefficient */
   for( c = 0; c < nitems; ++c )
   {
      cons = conss[c];

      /* check if each constraint is setppc constraint */
      assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(cons) ), "setppc", 6) );

      /* constraints which are (locally) disabled/redundant are not of
       * interest since the corresponding job is assigned to a packing
       */
      if( !SCIPconsIsEnabled(cons) )
         continue;

      if( SCIPgetNFixedonesSetppc(scip, cons) == 1 )
      {
         /* disable constraint locally */
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         continue;
      }

      /* dual value in original SCIP */
      dual = SCIPgetDualsolSetppc(scip, cons);
      
      SCIP_CALL( SCIPcreateVarBasic(subscip, &var, SCIPconsGetName(cons), 0.0, 1.0, dual, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(subscip, var) );

      vals[nvars] = weights[c];
      vars[nvars] = var;
      nvars++;

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(subscip, &var) );
   }

   /* create capacity constraint */
   SCIP_CALL( SCIPcreateConsBasicKnapsack(subscip, &cons, "capacity", nvars, vars, vals,
         capacity) );
   
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* add constraint of the branching decisions */
   SCIP_CALL( addBranchingDecisionConss(scip, subscip, vars, pricerdata->conshdlr) );

   /* avoid to generate columns which are fixed to zero */
   SCIP_CALL( addFixedVarsConss(scip, subscip, vars, conss, nitems) );

   SCIPfreeBufferArray(subscip, &vals);

   return SCIP_OKAY;
}

/**@} */

/**name Callback methods
 *
 * @{
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeBinpacking)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   if( pricerdata != NULL)
   {
      /* free memory */
      SCIPfreeMemoryArrayNull(scip, &pricerdata->conss);
      SCIPfreeMemoryArrayNull(scip, &pricerdata->weights);
      SCIPfreeMemoryArrayNull(scip, &pricerdata->ids);

      SCIPfreeMemory(scip, &pricerdata);
   }

   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitBinpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS* cons;
   int c;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get transformed constraints */
   for( c = 0; c < pricerdata->nitems; ++c )
   {
      cons = pricerdata->conss[c];

      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->conss[c]) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->conss[c]) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, pricerdata->conss[c]) );
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolBinpacking)
{
   SCIP_PRICERDATA* pricerdata;
   int c;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get release constraints */
   for( c = 0; c < pricerdata->nitems; ++c )
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->conss[c])) );
   }

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostBinpacking)
{  /*lint --e{715}*/
   SCIP* subscip;
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   int* ids;
   SCIP_Bool addvar;

   SCIP_SOL** sols;
   int nsols;
   int s;

   int nitems;
   SCIP_Longint capacity;

   SCIP_Real timelimit;
   SCIP_Real memorylimit;

   assert(scip != NULL);
   assert(pricer != NULL);

   (*result) = SCIP_DIDNOTRUN;

   /* get the pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   capacity = pricerdata->capacity;
   conss = pricerdata->conss;
   ids = pricerdata->ids;
   nitems = pricerdata->nitems;

   /* get the remaining time and memory limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* create problem in sub SCIP */
   SCIP_CALL( SCIPcreateProbBasic(subscip, "pricing") );
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* set time and memory limit */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   SCIP_CALL( SCIPallocMemoryArray(subscip, &vars, nitems) );

   /* initialization local pricing problem */
   SCIP_CALL( initPricing(scip, pricerdata, subscip, vars) );

   SCIPdebugMessage("solve pricer problem\n");

   /* solve sub SCIP */
   SCIP_CALL( SCIPsolve(subscip) );

   sols = SCIPgetSols(subscip);
   nsols = SCIPgetNSols(subscip);
   addvar = FALSE;

   /* loop over all solutions and create the corresponding column to master if the reduced cost are negative for master,
    * that is the objective value i greater than 1.0
    */
   for( s = 0; s < nsols; ++s )
   {
      SCIP_Bool feasible;
      SCIP_SOL* sol;

      /* the soultion should be sorted w.r.t. the objective function value */
      assert(s == 0 || SCIPisFeasGE(subscip, SCIPgetSolOrigObj(subscip, sols[s-1]), SCIPgetSolOrigObj(subscip, sols[s])));

      sol = sols[s];
      assert(sol != NULL);

      /* check if solution is feasible in original sub SCIP */
      SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, FALSE, FALSE ) );

      if( !feasible )
      {
         SCIPwarningMessage(scip, "solution in pricing problem (capacity <%d>) is infeasible\n", capacity);
         continue;
      }

      /* check if the solution has a value greater than 1.0 */
      if( SCIPisFeasGT(subscip, SCIPgetSolOrigObj(subscip, sol), 1.0) )
      {
         SCIP_VAR* var;
         SCIP_VARDATA* vardata;
         int* consids;
         char strtmp[SCIP_MAXSTRLEN];
         char name[SCIP_MAXSTRLEN];
         int nconss;
         int o;
         int v;

         SCIPdebug( SCIP_CALL( SCIPprintSol(subscip, sol, NULL, FALSE) ) );

         nconss = 0;
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "items");

         SCIP_CALL( SCIPallocBufferArray(scip, &consids, nitems) );

         /* check which variables are fixed -> which item belongs to this packing */
         for( o = 0, v = 0; o < nitems; ++o )
         {
            if( !SCIPconsIsEnabled(conss[o]) )
               continue;

            assert(SCIPgetNFixedonesSetppc(scip, conss[o]) == 0);

            if( SCIPgetSolVal(subscip, sol, vars[v]) > 0.5 )
            {
               (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", ids[o]);
               strcat(name, strtmp);

               consids[nconss] = o;
               nconss++;
            }
            else
               assert( SCIPisFeasEQ(subscip, SCIPgetSolVal(subscip, sol, vars[v]), 0.0) );

            v++;
         }

         SCIP_CALL( SCIPvardataCreateBinpacking(scip, &vardata, consids, nconss) );

         /* create variable for a new column with objective function coefficient 0.0 */
         SCIP_CALL( SCIPcreateVarBinpacking(scip, &var, name, 1.0, FALSE, TRUE, vardata) );

         /* add the new variable to the pricer store */
         SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
         addvar = TRUE;

         /* change the upper bound of the binary variable to lazy since the upper bound is already enforced due to
          * the objective function the set covering constraint; The reason for doing is that, is to avoid the bound
          * of x <= 1 in the LP relaxation since this bound constraint would produce a dual variable which might have
          * a positive reduced cost
          */
         SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

         /* check which variable are fixed -> which orders belong to this packing */
         for( v = 0; v < nconss; ++v )
         {
            assert(SCIPconsIsEnabled(conss[consids[v]]));
            SCIP_CALL( SCIPaddCoefSetppc(scip, conss[consids[v]], var) );
         }

         SCIPdebug(SCIPprintVar(scip, var, NULL) );
         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         SCIPfreeBufferArray(scip, &consids);
      }
      else
         break;
   }

   /* free pricer MIP */
   SCIPfreeMemoryArray(subscip, &vars);

   if( addvar || SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      (*result) = SCIP_SUCCESS;

   /* free sub SCIP */
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasBinpacking)
{  /*lint --e{715}*/

   /** @note In case of this binpacking example, the master LP should not get infeasible after branching, because of the
    *        way branching is performed. Therefore, the Farkas pricing is not implemented.
    *        1. In case of Ryan/Foster branching, the two items are selected in a way such that the sum of the LP values
    *           of all columns/packings containing both items is fractional. Hence, it exists at least one
    *           column/packing which contains both items and also at least one column/packing for each item containing
    *           this but not the other item. That means, branching in the "same" direction stays LP feasible since there
    *           exists at least one column/packing with both items and branching in the "differ" direction stays LP
    *           feasible since there exists at least one column/packing containing one item, but not the other.
    *        2. In case of variable branching, we only branch on fractional variables. If a variable is fixed to one,
    *           there is no issue.  If a variable is fixed to zero, then we know that for each item which is part of
    *           that column/packing, there exists at least one other column/packing containing this particular item due
    *           to the covering constraints.
    */
   SCIPwarningMessage(scip, "Current master LP is infeasible, but Farkas pricing was not implemented\n");
   SCIPABORT();

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates the binpacking variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerBinpacking(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create binpacking variable pricer data */
   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );

   pricerdata->conshdlr = SCIPfindConshdlr(scip, "samediff");
   assert(pricerdata->conshdlr != NULL);

   pricerdata->conss = NULL;
   pricerdata->weights = NULL;
   pricerdata->ids = NULL;
   pricerdata->nitems = 0;
   pricerdata->capacity = 0;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostBinpacking, pricerFarkasBinpacking, pricerdata) );

   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeBinpacking) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitBinpacking) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolBinpacking) );

   /* add binpacking variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}


/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerBinpackingActivate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< set covering constraints for the items */
   SCIP_Longint*         weights,            /**< weight of the items */
   int*                  ids,                /**< array of item ids */
   int                   nitems,             /**< number of items to be packed */
   SCIP_Longint          capacity            /**< capacity of the bins */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int c;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(weights != NULL);
   assert(nitems > 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* copy arrays */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &pricerdata->conss, conss, nitems) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &pricerdata->weights, weights, nitems) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &pricerdata->ids, ids, nitems) );

   pricerdata->nitems = nitems;
   pricerdata->capacity = capacity;

   SCIPdebugMessage("   nitems: %d capacity: %"SCIP_LONGINT_FORMAT"  \n", nitems, capacity);
   SCIPdebugMessage("      # profits    weights   x  \n");   /* capture constraints */

   /* capture all constraints */
   for( c = 0; c < nitems; ++c )
   {
      SCIP_CALL( SCIPcaptureCons(scip, conss[c]) );
      SCIPdebugPrintf("%4d %3"SCIP_LONGINT_FORMAT"\n", c, weights[c]);
   }

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/**@} */
