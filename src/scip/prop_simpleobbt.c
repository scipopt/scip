/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    prop_simpleobbt.c
 * @ingroup PROPAGATORS
 * @brief   simple optimization-based bound tightening propagator
 * @author  Stefan Weltge
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_simpleobbt.h"
#include "scip/prop_genvbounds.h"

#define PROP_NAME                       "simpleobbt"
#define PROP_DESC                       "simple optimization-based bound tightening propagator"
#define PROP_TIMING                     SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY                -1000000      /**< propagator priority */
#define PROP_FREQ                           0      /**< propagator frequency */
#define PROP_DELAY                       TRUE      /**< should propagation method be delayed, if other propagators
                                                    *   found reductions? */
#define DEFAULT_DUALFEASTOL              1e-9      /**< feasibility tolerance for reduced costs used in obbt; this value
                                                    *   is used if SCIP's dual feastol is greater */
#define DEFAULT_RELAXBOUNDS              TRUE      /**< should bounds be relaxed before solving the corresponding OBBT
                                                    *   LPs in order to gain more genbvounds? */
#define GENVBOUND_PROP_NAME             "genvbounds"


/*
 * Data structures
 */

/** bound data */
struct Bound
{
   SCIP_VAR*             var;                /**< variable */
   SCIP_Real             newval;             /**< stores a probably tighter value for this bound */
   SCIP_BOUNDTYPE        boundtype;          /**< type of bound */
   unsigned int          found:1;            /**< stores whether a probably tighter value for this bound was found */
};
typedef struct Bound BOUND;

/** propagator data */
struct SCIP_PropData
{
   BOUND**               bounds;             /**< array of interesting bounds */
   SCIP_ROW*             cutoffrow;          /**< pointer to current objective cutoff row */
   SCIP_PROP*            genvboundprop;      /**< pointer to genvbound propagator */
   SCIP_Longint          lastnode;           /**< number of last node where obbt was performed */
   SCIP_Real             dualfeastol;        /**< feasibility tolerance for reduced costs used in obbt; this value is
                                              *   used if SCIP's dual feastol is greater */
   SCIP_Bool             relaxbounds;        /**< should bounds be relaxed before solving the corresponding OBBT LPs in
                                              *   order to gain more genbvounds? */
   int                   nbounds;            /**< length of interesting bounds array */
};


/*
 * Local methods
 */

/** solves the LP and handles errors */
static
SCIP_RETCODE solveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlimit,            /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            error,              /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            optimal             /**< was the LP solved to optimalilty? */
   )
{
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(itlimit == -1 || itlimit >= 0);
   assert(error != NULL);
   assert(optimal != NULL);

   *optimal = FALSE;
   *error = FALSE;

   retcode = SCIPsolveDiveLP(scip, itlimit, error);
   lpsolstat = SCIPgetLPSolstat(scip);

   /* an error should not kill the overall solving process */
   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "   error while solving LP in obbt propagator; LP solve terminated with code <%d>\n", retcode);
      SCIPwarningMessage(scip, "   this does not affect the remaining solution procedure --> continue\n");

      *error = TRUE;

      return SCIP_OKAY;
   }

   if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      assert(!*error);
      *optimal = TRUE;
   }
#ifdef SCIP_DEBUG
   else
   {
      switch( lpsolstat )
      {
      case SCIP_LPSOLSTAT_ITERLIMIT:
         SCIPdebugMessage("   reached lp iteration limit\n");
         break;
      case SCIP_LPSOLSTAT_TIMELIMIT:
         SCIPdebugMessage("   reached time limit while solving lp\n");
         break;
      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         SCIPdebugMessage("   lp was unbounded\n");
         break;
      case SCIP_LPSOLSTAT_NOTSOLVED:
         SCIPdebugMessage("   lp was not solved\n");
         break;
      case SCIP_LPSOLSTAT_ERROR:
         SCIPdebugMessage("   an error occured during solving lp\n");
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
      case SCIP_LPSOLSTAT_OBJLIMIT:
      case SCIP_LPSOLSTAT_OPTIMAL: /* should not appear because it is handled earlier */
      default:
         SCIPdebugMessage("   received an unexpected solstat during solving lp: %d\n", lpsolstat);
      }
   }
#endif

   return SCIP_OKAY;
}

/** adds the objective cutoff to the LP; must be in diving mode */
static
SCIP_RETCODE addObjCutoff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_ROW* row;
   SCIP_VAR** vars;
   char rowname[SCIP_MAXSTRLEN];

   int nvars;
   int i;

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(propdata != NULL);
   assert(propdata->cutoffrow == NULL);

   if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
   {
      SCIPdebugMessage("no objective cutoff since there is no cutoff bound\n");
      return SCIP_OKAY;
   }

   SCIPdebugMessage("create objective cutoff and add it to the LP\n");

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create objective cutoff row; set local flag to FALSE since primal cutoff is globally valid */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "obbt_objcutoff");
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, rowname, -SCIPinfinity(scip), SCIPgetCutoffbound(scip), FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], SCIPvarGetObj(vars[i])) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* add row to the LP */
   SCIP_CALL( SCIPaddRowDive(scip, row) );

   propdata->cutoffrow = row;
   assert(SCIProwIsInLP(propdata->cutoffrow));

   return SCIP_OKAY;
}

/** determines, whether a variable is already locally fixed */
static
SCIP_Bool varIsFixedLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to check */
   )
{
   return SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
}

/** applies possible bound changes that were found */
static
SCIP_RETCODE applyBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   int i;

   assert(scip != NULL);
   assert(!SCIPinDive(scip));
   assert(propdata != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND);

   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;                          /* shortcut to the current bound */
      SCIP_Bool infeas;                      /* stores wether a tightening approach forced an infeasibilty */
      SCIP_Bool tightened;                   /* stores wether a tightening approach was successful */

      bound = propdata->bounds[i];

      if( bound->found )
      {
         if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, bound->var, bound->newval, FALSE, &infeas, &tightened) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarUb(scip, bound->var, bound->newval, FALSE, &infeas, &tightened) );
         }

         /* handle information about the success */
         if( infeas )
         {
            *result = SCIP_CUTOFF;
            SCIPdebugMessage("cut off\n");
            break;
         }

         if( tightened )
         {
            *result = SCIP_REDUCEDDOM;
         }
      }
   }

   return SCIP_OKAY;
}

/** tries to tighten a bound in diving mode  */
static
SCIP_RETCODE tightenBoundDive(
   SCIP*                 scip,               /**< SCIP data structure */
   BOUND*                bound,              /**< bound that could be tightened */
   SCIP_Real             newval              /**< new bound value */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(bound != NULL);

   /* get old bounds */
   lb = SCIPgetVarLbDive(scip, bound->var);
   ub = SCIPgetVarUbDive(scip, bound->var);

   /* round bounds new value if variable is integral */
   if( SCIPvarIsIntegral(bound->var) )
   {
      newval = bound->boundtype == SCIP_BOUNDTYPE_LOWER ? SCIPceil(scip, newval) : SCIPfloor(scip, newval);
   }

   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisLbBetter(scip, newval, lb, ub) )
   {
      SCIP_CALL( SCIPchgVarLbDive(scip, bound->var, newval) );
   }
   else if( bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisUbBetter(scip, newval, lb, ub) )
   {
      SCIP_CALL( SCIPchgVarUbDive(scip, bound->var, newval) );
   }

   return SCIP_OKAY;
}

/** creates a genvbound if the dual LP solution provides such information
 *
 *  Consider the problem
 *
 *     min { +/- x_i : obj * x <= z, lb <= Ax <= ub, l <= x <= u },
 *
 *  where z is the current cutoff bound. Let (mu, nu, gamma, alpha, beta) >= 0 be the optimal solution of the dual of
 *  problem (P), where the variables correspond to the primal inequalities in the following way:
 *
 *           Ax >=  lb    <->   mu
 *          -Ax >= -ub    <->   nu
 *     -obj * x >=  -z    <->   gamma
 *            x >=   l    <->   alpha
 *           -x >=  -u    <->   beta
 *
 *  Fixing these multipliers, by weak duality, we obtain the inequality
 *
 *     +/- x_i >= lb*mu - ub*nu - z*gamma + l*alpha - u*beta
 *
 *  that holds for all primal feasible points x with objective value at least z. Setting
 *
 *     c = lb*mu - ub*nu, redcost_k = alpha_k - beta_k
 *
 *  we obtain the inequality
 *
 *     +/- x_i >= sum ( redcost_k * x_k ) + (-gamma) * cutoff_bound + c,
 *
 *  that holds for all primal feasible points with objective value at least cutoff_bound. Therefore, the latter
 *  inequality can be added as a generalized variable bound.
 */
static
SCIP_RETCODE createGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   BOUND*                bound               /**< bound of x_i */
   )
{
   assert(scip != NULL);
   assert(bound != NULL);
   assert(propdata != NULL);
   assert(propdata->genvboundprop != NULL);

   /* make sure we are in dive mode having an optimal LP solution */
   assert(SCIPinDive(scip));
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* only genvbounds created in the root node are globally valid */
   assert(SCIPgetDepth(scip) == 0);

   SCIPdebugMessage("      try to create a genvbound for <%s>...\n", SCIPvarGetName(bound->var));

   /* a genvbound with a multiplier for x_i would not help us */
   if( SCIPisZero(scip, SCIPgetVarRedcost(scip, bound->var)) )
   {
      SCIP_VAR** vars;                          /* global variables array */
      SCIP_VAR** genvboundvars;                 /* genvbound variables array */

      SCIP_VAR* xi;                             /* variable x_i */

      SCIP_Real* genvboundcoefs;                /* genvbound coefficients array */

      SCIP_Real gamma_dual;                     /* dual multiplier of objective cutoff */

      int k;                                    /* variable for indexing global variables array */
      int ncoefs;                               /* number of nonzero coefficients in genvbound */
      int nvars;                                /* number of global variables */

      /* set x_i */
      xi = bound->var;

      /* get variable data */
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

      /* count nonzero coefficients in genvbound */
      ncoefs = 0;
      for( k = 0; k < nvars; k++ )
      {
         if( !SCIPisZero(scip, SCIPgetVarRedcost(scip, vars[k])) )
         {
            assert(vars[k] != xi);
            ncoefs++;
         }
      }

      /* get dual multiplier for the objective cutoff (set to zero if there is no) */
      if( propdata->cutoffrow == NULL )
      {
         gamma_dual = 0.0;
      }
      else
      {
         assert(!SCIPisInfinity(scip, SCIPgetCutoffbound(scip)));

         /* note that the objective cutoff is of the form
          *    -inf <= obj * x <= cutoff_bound
          * but we want the positive dual multiplier!
          */
         gamma_dual = -SCIProwGetDualsol(propdata->cutoffrow);
      }

      /* we need at least one nonzero coefficient or a nonzero dual multiplier for the objective cutoff */
      if( ncoefs > 0 || !SCIPisZero(scip, gamma_dual) )
      {
         SCIP_Real c;                           /* helper variable to calculate constant term in genvbound */
         int idx;                               /* variable for indexing genvbound's coefficients array */

         /* there should be no coefficient for x_i */
         assert(SCIPisZero(scip, SCIPgetVarRedcost(scip, xi)));

         /* allocate memory for storing the genvbounds right-hand side variables and coefficients */
         SCIP_CALL( SCIPallocMemoryArray(scip, &(genvboundvars), ncoefs) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(genvboundcoefs), ncoefs) );

         /* set c = lb*mu - ub*nu - z*gamma + l*alpha - u*beta */
         c = SCIPgetLPObjval(scip);

         /* subtract ( - z * gamma ) from c */
         c += SCIPgetCutoffbound(scip) * gamma_dual;

         /* subtract ( l*alpha - u*beta ) from c and set the coefficients of the variables */
         idx = 0;
         for( k = 0; k < nvars; k++ )
         {
            SCIP_VAR* xk;
            SCIP_Real redcost;

            xk = vars[k];
            redcost = SCIPgetVarRedcost(scip, xk);

            if( !SCIPisZero(scip, redcost) )
            {
               assert(xk != xi);

               /* store coefficients */
               assert(idx < ncoefs);
               genvboundvars[idx] = xk;
               genvboundcoefs[idx] = redcost;
               idx++;

               /* if redcost > 0, then redcost = alpha_k, otherwise redcost = - beta_k */
               c -= redcost > 0 ? redcost * SCIPgetVarLbDive(scip, xk) : redcost * SCIPgetVarUbDive(scip, xk);
            }
         }

         /* add genvbound */
         SCIPdebugMessage("         adding genvbound\n");
         SCIP_CALL( SCIPgenVBoundAdd(scip, propdata->genvboundprop, genvboundvars, xi, genvboundcoefs, ncoefs,
               !SCIPisPositive(scip, gamma_dual) ? 0.0 : -gamma_dual, c, bound->boundtype) );

         /* free arrays */
         SCIPfreeMemoryArray(scip, &genvboundcoefs);
         SCIPfreeMemoryArray(scip, &genvboundvars);
      }
      else
      {
         SCIPdebugMessage("         trivial genvbound, skipping\n");
      }
   }
   else
   {
      SCIPdebugMessage("         found multiplier for <%s>: %g, skipping\n",
         SCIPvarGetName(bound->var), SCIPgetVarRedcost(scip, bound->var));
   }

   return SCIP_OKAY;
}

/** tries to find tighter values for bounds and stores them in the bound data structure */
static
SCIP_RETCODE findNewBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */

   int i;
   int nvars;                                /* number of the problems variables */

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(propdata != NULL);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* set objective coefficients to zero */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], 0.0) );
   }

   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;
      SCIP_VAR* var;

      SCIP_Real oldbound;
      SCIP_Bool error;
      SCIP_Bool optimal;                     /* was the LP solved to optimalilty? */

      oldbound = 0.0;
      bound = propdata->bounds[i];
      var = bound->var;

      SCIPdebugMessage("   applying obbt on %s bound of <%s> (local bounds: [%f,%f])\n",
         bound->boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

      /* set objective coefficient to +/- 1 (note that we minimize) */
      SCIP_CALL( SCIPchgVarObjDive(scip, var, (bound->boundtype == SCIP_BOUNDTYPE_LOWER) ? 1.0 : -1.0 ) );

      /* relax bound */
      SCIPdebugMessage("relaxbounds: %d\n", propdata->relaxbounds);
      if( propdata->relaxbounds )
      {
         if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
         {
            oldbound = SCIPgetVarLbDive(scip, var);
            SCIP_CALL( SCIPchgVarLbDive(scip, var, oldbound - MAX(0.1 * oldbound, 0.1)) );
            SCIPdebugMessage("relaxed lower bound from %f to %f\n", oldbound, SCIPgetVarLbDive(scip, var));
         }
         else
         {
            oldbound = SCIPgetVarUbDive(scip, var);
            SCIP_CALL( SCIPchgVarUbDive(scip, var, oldbound + MAX(0.1 * oldbound, 0.1)) );
            SCIPdebugMessage("relaxed upper bound from %f to %f\n", oldbound, SCIPgetVarUbDive(scip, var));
         }
      }

      /* solve LP */
      SCIP_CALL( solveLP(scip, -1, &error, &optimal) );
      SCIPdebugMessage("solved lp\n");

      /* stop this procedure if an error occured */
      if( error )
      {
         SCIPdebugMessage("aborting since lp error\n");
         return SCIP_OKAY;
      }

      if( optimal )
      {
         /* store this value in the bound data structure */
         bound->newval = SCIPvarGetLPSol(var);
         bound->found = TRUE;

         SCIPdebugMessage("      LP value: %f\n", bound->newval);

         /* in root node we may want to create a genvbound (independent of tightening success) */
         if( SCIPgetDepth(scip) == 0 && propdata->genvboundprop != NULL )
         {
            SCIP_CALL( createGenVBound(scip, propdata, bound) );
         }

         /* restore bound */
         if( propdata->relaxbounds )
         {
            if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
            {
               SCIP_CALL( SCIPchgVarLbDive(scip, var, oldbound) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarUbDive(scip, var, oldbound) );
            }
         }

         /* try to tighten bound in dive mode */
         SCIP_CALL( tightenBoundDive(scip, bound, bound->newval) );
      }

      /* set objective coefficient back to zero */
      SCIP_CALL( SCIPchgVarObjDive(scip, var, 0.0 ) );
   }

   return SCIP_OKAY;
}

/** main function of obbt */
static
SCIP_RETCODE applyObbt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint          itlimit,            /**< LP iteration limit (-1: no limit) */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_Real olddualfeastol;
   SCIP_Bool error;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   *result = SCIP_DIDNOTFIND;
   olddualfeastol = SCIPdualfeastol(scip);
   error = FALSE;

   /* start diving */
   SCIP_CALL( SCIPstartDive(scip) );

   /* try to (re-)solve root LP */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIP_Bool optimal;
      optimal = FALSE;

      SCIPdebugMessage("try to (re-)solve root lp...\n");
      SCIP_CALL( solveLP(scip, -1, &error, &optimal) );
   }

   if( !error )
   {
      /* set dual feastol */
      if( propdata->dualfeastol < olddualfeastol )
      {
         SCIP_CALL( SCIPchgDualfeastol(scip, propdata->dualfeastol) );
      }

      /* add objective cutoff */
      SCIP_CALL( addObjCutoff(scip, propdata) );

      /* try to find new bounds and store them in the bound data structure */
      SCIP_CALL( findNewBounds(scip, propdata) );
   }

   if( error )
   {
      SCIPdebugMessage("skipping obbt since an error occured in (re-)solving the LP\n");
   }

   /* reset dual feastol */
   SCIP_CALL( SCIPchgDualfeastol(scip, olddualfeastol) );

   /* end diving */
   SCIP_CALL( SCIPendDive(scip) );

   /* release cutoff row if there is one */
   if( propdata->cutoffrow != NULL )
   {
      assert(!SCIProwIsInLP(propdata->cutoffrow));
      SCIP_CALL( SCIPreleaseRow(scip, &(propdata->cutoffrow)) );
   }

   /* apply buffered bound changes */
   SCIP_CALL( applyBoundChgs(scip, propdata, result) );

   return SCIP_OKAY;
}

/** determines whether a variable is interesting */
static
SCIP_Bool varIsInteresting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Bool             nlrows              /**< is the variable contained in a nonlinear row? */
   )
{
   assert(SCIPgetDepth(scip) == 0);

   return !SCIPvarIsBinary(var) && !varIsFixedLocal(scip, var) && nlrows;
}

/** initializes interesting bounds */
static
SCIP_RETCODE initBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   SCIP_Bool* nlarray;                       /* array that stores for each variable wether it is contained in a
                                              * nonlinear constraint appears */

   int bdidx;                                /* bound index inside propdata->bounds */
   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* count nonlinearities */
   SCIP_CALL( SCIPallocBufferArray(scip, &nlarray, nvars) );
   if( SCIPisNLPConstructed(scip) )
   {
      assert(SCIPgetNNLPVars(scip) == nvars);
      SCIP_CALL( SCIPgetNLPVarsNonlinearRows(scip, nlarray) );
   }
   else
   {
      BMSclearMemoryArray(nlarray, nvars);
   }

   /* allocate interesting bounds array */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->bounds), 2 * nvars) );

   /* get all interesting variables and their bounds */
   bdidx = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varIsInteresting(scip, vars[i], nlarray[i]) )
      {
         BOUND** bdaddress;

         /* create lower bound */
         bdaddress = &(propdata->bounds[bdidx]);
         SCIP_CALL( SCIPallocMemory(scip, bdaddress) );
         propdata->bounds[bdidx]->boundtype = SCIP_BOUNDTYPE_LOWER;
         propdata->bounds[bdidx]->var = vars[i];
         propdata->bounds[bdidx]->found = FALSE;
         propdata->bounds[bdidx]->newval = 0.0;
         bdidx++;

         /* create upper bound */
         bdaddress = &(propdata->bounds[bdidx]);
         SCIP_CALL( SCIPallocMemory(scip, bdaddress) );
         propdata->bounds[bdidx]->boundtype = SCIP_BOUNDTYPE_UPPER;
         propdata->bounds[bdidx]->var = vars[i];
         propdata->bounds[bdidx]->found = FALSE;
         propdata->bounds[bdidx]->newval = 0.0;
         bdidx++;
      }
   }

   /* free memory for buffering nonlinearities */
   assert(nlarray != NULL);
   SCIPfreeBufferArray(scip, &nlarray);

   /* set number of interesting bounds */
   propdata->nbounds = bdidx;

   /* resize propdata->bounds array */
   if( propdata->nbounds > 0 )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(propdata->bounds), propdata->nbounds) );
   }
   else
   {
      assert(propdata->nbounds == 0);
      SCIPfreeMemoryArray(scip, &(propdata->bounds));
   }

   SCIPdebugMessage("problem has %d/%d interesting bounds\n", propdata->nbounds, 2 * nvars);
   SCIPstatisticMessage("problem has %d/%d interesting bounds\n", propdata->nbounds, 2 * nvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolSimpleObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->bounds = NULL;
   propdata->nbounds = -1;
   propdata->cutoffrow = NULL;
   propdata->lastnode = -1;

   /* if genvbounds propagator is not available, we cannot create genvbounds */
   propdata->genvboundprop = SCIPfindProp(scip, GENVBOUND_PROP_NAME);

   SCIPdebugMessage("creating genvbounds: %s\n", propdata->genvboundprop != NULL ? "true" : "false");

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSimpleObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   /* do not run in: presolving, repropagation, probing mode */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip) || SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* only run for nonlinear problems, i.e., if NLP is constructed */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMessage("NLP not constructed, skipping obbt\n");
      SCIPstatisticMessage("NLP not constructed\n");
      return SCIP_OKAY;
   }

   /* only run if LP all columns are in the LP, i.e., the LP is a relaxation; e.g., do not run if pricers are active
    * since pricing is not performed in diving mode
    */
   if( !SCIPallColsInLP(scip) )
   {
      SCIPdebugMessage("not all columns in LP, skipping obbt\n");
      SCIPstatisticMessage("not all columns in LP\n");
      return SCIP_OKAY;
   }

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* ensure that bounds are initialized */
   if( propdata->nbounds == -1 )
   {
      /* bounds must be initialized at root node */
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( initBounds(scip, propdata) );
      }
      else
      {
         assert(!SCIPinProbing(scip));
         SCIPpropSetFreq(prop, -1);
         return SCIP_OKAY;
      }
   }
   assert(propdata->nbounds >= 0);

   /* disable obbt if there are no interesting bounds */
   if( propdata->nbounds == 0 )
   {
      SCIPdebugMessage("there are no interesting bounds, disabling obbt\n");
      SCIPstatisticMessage("there are no interesting bounds, disabling obbt\n");
      SCIPpropSetFreq(prop, -1);

      return SCIP_OKAY;
   }

   /* only run once in a node */
   if( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode )
   {
      return SCIP_OKAY;
   }

   SCIPdebugMessage("applying obbt for problem <%s> at depth %d\n", SCIPgetProbName(scip), SCIPgetDepth(scip));

   /* apply obbt */
   SCIP_CALL( applyObbt(scip, propdata, -1, result) );
   assert(*result != SCIP_DIDNOTRUN);

   /* set current node as last node */
   propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropSimpleObbt)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolSimpleObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* free memory allocated for the bounds */
   if( propdata->nbounds > 0 )
   {
      assert(propdata->bounds != NULL);

      /* free bounds */
      for( i = propdata->nbounds - 1; i >= 0; i-- )
      {
         SCIPfreeMemory(scip, &(propdata->bounds[i]));
      }
      SCIPfreeMemoryArray(scip, &(propdata->bounds));
   }

   propdata->nbounds = -1;

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSimpleObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the obbt propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropSimpleObbt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create obbt propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecSimpleObbt, propdata) );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSimpleObbt) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolSimpleObbt) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolSimpleObbt) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropSimpleObbt) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/relaxbounds",
         "should bounds be relaxed before solving the corresponding OBBT LPs in order to gain more genbvounds?",
         &propdata->relaxbounds, TRUE, DEFAULT_RELAXBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/dualfeastol",
         "feasibility tolerance for reduced costs used in obbt; this value is used if SCIP's dual feastol is greater",
         &propdata->dualfeastol, FALSE, DEFAULT_DUALFEASTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
