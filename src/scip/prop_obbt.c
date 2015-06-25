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

/**@file    prop_obbt.c
 * @ingroup PROPAGATORS
 * @brief   optimization-based bound tightening propagator
 * @author  Stefan Weltge
 * @author  Benjamin Mueller
 */

/**@todo if bound tightenings of other propagators are the reason for lpsolstat != SCIP_LPSOLSTAT_OPTIMAL, resolve LP */
/**@todo only run more than once in root node if primal bound improved or many cuts were added to the LP */
/**@todo filter bounds of a variable already if SCIPisLbBetter()/SCIPisUbBetter() would return FALSE */
/**@todo improve warmstarting of LP solving */
/**@todo include bound value (finite/infinite) into getScore() function */
/**@todo use unbounded ray in filtering */
/**@todo do we want to run if the LP is unbounded, maybe for infinite variable bounds? */
/**@todo add first filter round in direction of objective function */
/**@todo implement conflict resolving callback by calling public method of genvbounds propagator, since the reason are
 *       exactly the variable bounds with nonnegative reduced costs stored in the right-hand side of the generated
 *       generalized variable bound (however, this only makes sense if we run locally)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_obbt.h"
#include "scip/prop_genvbounds.h"
#include "scip/debug.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_abspower.h"
#include "scip/cons_bivariate.h"

#define PROP_NAME                       "obbt"
#define PROP_DESC                       "optimization-based bound tightening propagator"
#define PROP_TIMING                     SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY                -1000000      /**< propagator priority */
#define PROP_FREQ                           0      /**< propagator frequency */
#define PROP_DELAY                       TRUE      /**< should propagation method be delayed, if other propagators
                                                    *   found reductions? */

#define DEFAULT_CREATE_GENVBOUNDS        TRUE      /**< should obbt try to provide genvbounds if possible? */
#define DEFAULT_FILTERING_NORM           TRUE      /**< should coefficients in filtering be normalized w.r.t. the
                                                    *   domains sizes? */
#define DEFAULT_APPLY_FILTERROUNDS      FALSE      /**< try to filter bounds in so-called filter rounds by solving
                                                    *   auxiliary LPs? */
#define DEFAULT_APPLY_TRIVIALFITLERING   TRUE      /**< should obbt try to use the LP solution to filter some bounds? */
#define DEFAULT_GENVBDSDURINGFILTER      TRUE      /**< try to genrate genvbounds during trivial and aggressive filtering? */
#define DEFAULT_DUALFEASTOL              1e-9      /**< feasibility tolerance for reduced costs used in obbt; this value
                                                    *   is used if SCIP's dual feastol is greater */
#define DEFAULT_CONDITIONLIMIT           -1.0      /**< maximum condition limit used in LP solver (-1.0: no limit) */
#define DEFAULT_BOUNDSTREPS             0.001      /**< minimal relative improve for strengthening bounds */
#define DEFAULT_FILTERING_MIN               2      /**< minimal number of filtered bounds to apply another filter
                                                    *   round */
#define DEFAULT_ITLIMITFACTOR            10.0      /**< multiple of root node LP iterations used as total LP iteration
                                                    *   limit for obbt (<= 0: no limit ) */
#define DEFAULT_MINITLIMIT              5000L      /**< minimum LP iteration limit */
#define DEFAULT_ONLYNONCONVEXVARS       FALSE      /**< only apply obbt on non-convex variables */
#define DEFAULT_TIGHTINTBOUNDSPROBING    TRUE      /**< should bounds of integral variables be tightened during
                                                    *   the probing mode? */
#define DEFAULT_TIGHTCONTBOUNDSPROBING  FALSE      /**< should bounds of continuous variables be tightened during
                                                    *   the probing mode? */
#define DEFAULT_ORDERINGALGO                1      /**< which type of ordering algorithm should we use?
                                                    *   (0: no, 1: greedy, 2: greedy reverse) */
#define OBBT_SCOREBASE                      5      /**< base that is used to calculate a bounds score value */
#define GENVBOUND_PROP_NAME             "genvbounds"
#define INTERVALINFTY                   1E+43      /**< value for infinity in interval operations */

#define DEFAULT_SEPARATESOL             FALSE      /**< should the obbt LP solution be separated? note that that by
                                                    *   separating solution OBBT will apply all bound tightenings
                                                    *   immediatly */
#define DEFAULT_SEPAMINITER                 0      /**< minimum number of iteration spend to separate an obbt LP solution */
#define DEFAULT_SEPAMAXITER                10      /**< maximum number of iteration spend to separate an obbt LP solution */
#define DEFAULT_GENVBDSDURINGSEPA        TRUE      /**< try to create genvbounds during separation process? */
#define DEFAULT_PROPAGATEFREQ               0      /**< trigger a propagation round after that many bound tightenings
                                                    *   (0: no propagation) */

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/*
 * Data structures
 */

/** bound data */
struct Bound
{
   SCIP_VAR*             var;                /**< variable */
   SCIP_Real             newval;             /**< stores a probably tighter value for this bound */
   SCIP_BOUNDTYPE        boundtype;          /**< type of bound */
   unsigned int          score;              /**< score value that is used to group bounds */
   unsigned int          filtered:1;         /**< thrown out during pre-filtering step */
   unsigned int          found:1;            /**< stores whether a probably tighter value for this bound was found */
   unsigned int          done:1;             /**< has this bound been processed already? */
   unsigned int          nonconvex:1;        /**< is this bound affecting a nonconvex term? */
   int                   index;              /**< unique index */
};
typedef struct Bound BOUND;


/** propagator data */
struct SCIP_PropData
{
   BOUND**               bounds;             /**< array of interesting bounds */
   SCIP_ROW*             cutoffrow;          /**< pointer to current objective cutoff row */
   SCIP_PROP*            genvboundprop;      /**< pointer to genvbound propagator */
   SCIP_Longint          lastnode;           /**< number of last node where obbt was performed */
   SCIP_Longint          npropagatedomreds;  /**< number of domain reductions found during propagation */
   SCIP_Longint          nprobingiterations; /**< number of LP iterations during the probing mode */
   SCIP_Longint          nfilterlpiters;     /**< number of LP iterations spend for filtering */
   SCIP_Longint          minitlimit;         /**< minimum LP iteration limit */
   SCIP_Real             dualfeastol;        /**< feasibility tolerance for reduced costs used in obbt; this value is
                                              *   used if SCIP's dual feastol is greater */
   SCIP_Real             conditionlimit;     /**< maximum condition limit used in LP solver (-1.0: no limit) */
   SCIP_Real             boundstreps;        /**< minimal relative improve for strengthening bounds */
   SCIP_Real             itlimitfactor;      /**< LP iteration limit for obbt will be this factor times total LP
                                              *   iterations in root node */
   SCIP_Bool             applyfilterrounds;  /**< apply filter rounds? */
   SCIP_Bool             applytrivialfilter; /**< should obbt try to use the LP solution to filter some bounds? */
   SCIP_Bool             genvbdsduringfilter;/**< should we try to generate genvbounds during trivial and aggressive
                                              *   filtering? */
   SCIP_Bool             genvbdsduringsepa;  /**< try to create genvbounds during separation process? */
   SCIP_Bool             creategenvbounds;   /**< should obbt try to provide genvbounds if possible? */
   SCIP_Bool             normalize;          /**< should coefficients in filtering be normalized w.r.t. the domains
                                              *   sizes? */
   SCIP_Bool             onlynonconvexvars;  /**< only apply obbt on non-convex variables */
   SCIP_Bool             tightintboundsprobing; /**< should bounds of integral variables be tightened during
                                              *   the probing mode? */
   SCIP_Bool             tightcontboundsprobing;/**< should bounds of continuous variables be tightened during
                                              *   the probing mode? */
   SCIP_Bool             separatesol;        /**< should the obbt LP solution be separated? note that that by
                                              *   separating solution OBBT will apply all bound tightenings
                                              *   immediatly */
   int                   orderingalgo;       /**< which type of ordering algorithm should we use?
                                              *   (0: no, 1: greedy, 2: greedy reverse) */
   int                   nbounds;            /**< length of interesting bounds array */
   int                   nminfilter;         /**< minimal number of filtered bounds to apply another filter round */
   int                   nfiltered;          /**< number of filtered bounds by solving auxiliary variables */
   int                   ntrivialfiltered;   /**< number of filtered bounds because the LP value was equal to the bound */
   int                   nsolvedbounds;      /**< number of solved bounds during the loop in applyObbt() */
   int                   ngenvboundsprobing; /**< number of non-trivial genvbounds generated and added during obbt */
   int                   ngenvboundsaggrfil; /**< number of non-trivial genvbounds found during aggressive filtering */
   int                   ngenvboundstrivfil; /**< number of non-trivial genvbounds found during trivial filtering */
   int                   lastidx;            /**< index to store the last undone and unfiltered bound */
   int                   sepaminiter;        /**< minimum number of iteration spend to separate an obbt LP solution */
   int                   sepamaxiter;        /**< maximum number of iteration spend to separate an obbt LP solution */
   int                   propagatefreq;      /**< trigger a propagation round after that many bound tightenings
                                              *   (0: no propagation) */
   int                   propagatecounter;   /**< number of bound tightenings since the last propagation round */
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

   retcode = SCIPsolveProbingLP(scip, itlimit, error, NULL);

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

/** adds the objective cutoff to the LP; must be in probing mode */
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
   assert(SCIPinProbing(scip));
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
   SCIP_CALL( SCIPaddRowProbing(scip, row) );

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

/** sets objective to minimize or maximize a single variable */
static
SCIP_RETCODE setObjProbing(
   SCIP*                 scip,
   SCIP_PROPDATA*        propdata,
   BOUND*                bound,
   SCIP_Real             coef
   )
{
#ifdef SCIP_DEBUG
   SCIP_VAR** vars;
   int nvars;
   int counter;
   int i;
#endif

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( bound != NULL );

   /* set the objective for bound->var */
   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, coef) );
   }
   else
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, -coef) );
   }

#ifdef SCIP_DEBUG
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   counter = 0;

   for( i = 0; i < nvars; ++i )
   {
      if( SCIPgetVarObjProbing(scip, vars[i]) != 0.0 )
         ++counter;
   }

   assert((counter == 0 && coef == 0.0) || (counter == 1 && coef != 0.0));
#endif

   return SCIP_OKAY;
}

/** determines whether variable should be included in the right-hand side of the generalized variable bound */
static
SCIP_Bool includeVarGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to check */
   )
{
   SCIP_Real redcost;

   assert(scip != NULL);
   assert(var != NULL);

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return FALSE;

   redcost = SCIPgetVarRedcost(scip, var);
   assert(redcost != SCIP_INVALID); /*lint !e777 */

   if( redcost == SCIP_INVALID ) /*lint !e777 */
      return FALSE;

   if( redcost < SCIPdualfeastol(scip) && redcost > -SCIPdualfeastol(scip) )
      return FALSE;

   return TRUE;
}

/** returns number of LP iterations left (-1: no limit ) */
static
int getIterationsLeft(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          nolditerations,     /**< iterations count at the beginning of the corresponding function */
   SCIP_Longint          itlimit             /**< LP iteration limit (-1: no limit) */
   )
{
   SCIP_Longint itsleft;

   assert(scip != NULL);
   assert(nolditerations >= 0);
   assert(itlimit == -1 || itlimit >= 0);

   if( itlimit == -1 )
   {
      SCIPdebugMessage("iterations left: unlimited\n");
      return -1;
   }
   else
   {
      itsleft = itlimit - ( SCIPgetNLPIterations(scip) - nolditerations );
      itsleft = MAX(itsleft, 0);
      itsleft = MIN(itsleft, INT_MAX);

      SCIPdebugMessage("iterations left: %d\n", (int) itsleft);
      return (int) itsleft;
   }
}

/** returns the objective coefficient for a variable's bound that will be chosen during filtering */
static
SCIP_Real getFilterCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_VAR*             var,                /**< variable */
   SCIP_BOUNDTYPE        boundtype           /**< boundtype to be filtered? */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(var != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   /* this function should not be called for fixed variables */
   assert(!varIsFixedLocal(scip, var));

   /* infinite bounds will not be reached */
   if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisInfinity(scip, -lb) )
      return 0.0;
   if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisInfinity(scip, ub) )
      return 0.0;

   if( propdata->normalize )
   {
      /* if the length of the domain is too large then the coefficient should be set to +/- 1.0 */
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisInfinity(scip, ub) )
         return 1.0;
      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisInfinity(scip, -lb) )
         return -1.0;

      /* otherwise the coefficient is +/- 1.0 / ( ub - lb ) */
      return boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 / (ub - lb) : -1.0 / (ub - lb);
   }
   else
   {
      return boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 : -1.0;
   }
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
   BOUND*                bound,              /**< bound of x_i */
   SCIP_Bool*            found               /**< pointer to store if we have found a non-trivial genvbound */
   )
{
   assert(scip != NULL);
   assert(bound != NULL);
   assert(propdata != NULL);
   assert(propdata->genvboundprop != NULL);
   assert(found != NULL);

   *found = FALSE;

   /* make sure we are in probing mode having an optimal LP solution */
   assert(SCIPinProbing(scip));

   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* only genvbounds created in the root node are globally valid
    *
    * note: depth changes to one if we use the probing mode to solve the obbt LPs
    */
   assert(SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1));

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
         if( includeVarGenVBound(scip, vars[k]) )
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
         SCIP_Bool addgenvbound;                /* if everything is fine with the redcosts and the bounds, add the genvbound */
         SCIP_Real c;                           /* helper variable to calculate constant term in genvbound */
         int idx;                               /* variable for indexing genvbound's coefficients array */

         /* add the bound if the bool is still TRUE after the loop */
         addgenvbound = TRUE;

         /* there should be no coefficient for x_i */
         assert(SCIPisZero(scip, SCIPgetVarRedcost(scip, xi)));

         /* allocate memory for storing the genvbounds right-hand side variables and coefficients */
         SCIP_CALL( SCIPallocBufferArray(scip, &(genvboundvars), ncoefs) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(genvboundcoefs), ncoefs) );

         /* set c = lb*mu - ub*nu - z*gamma + l*alpha - u*beta */
         c = SCIPgetLPObjval(scip);

         /* subtract ( - z * gamma ) from c */
         c += SCIPgetCutoffbound(scip) * gamma_dual;

         /* subtract ( l*alpha - u*beta ) from c and set the coefficients of the variables */
         idx = 0;
         for( k = 0; k < nvars; k++ )
         {
            SCIP_VAR* xk;

            xk = vars[k];

            if( includeVarGenVBound(scip, xk) )
            {
               SCIP_Real redcost;

               redcost = SCIPgetVarRedcost(scip, xk);

               assert(redcost != SCIP_INVALID); /*lint !e777 */
               assert(xk != xi);

               /* in this case dont add a genvbound */
               if( ( (redcost > SCIPdualfeastol(scip))  && SCIPisInfinity(scip, -SCIPvarGetLbLocal(xk)) ) ||
                  ( (redcost < -SCIPdualfeastol(scip))  && SCIPisInfinity(scip, SCIPvarGetUbLocal(xk)) ) )
               {
                  addgenvbound = FALSE;
                  break;
               }

               /* store coefficients */
               assert(idx < ncoefs);
               genvboundvars[idx] = xk;
               genvboundcoefs[idx] = redcost;
               idx++;

               /* if redcost > 0, then redcost = alpha_k, otherwise redcost = - beta_k */
               assert(redcost <= 0 || !SCIPisInfinity(scip, -SCIPvarGetLbLocal(xk)));
               assert(redcost >= 0 || !SCIPisInfinity(scip, SCIPvarGetUbLocal(xk)));
               c -= redcost > 0 ? redcost * SCIPvarGetLbLocal(xk) : redcost * SCIPvarGetUbLocal(xk);
            }
         }

         assert(!addgenvbound || idx == ncoefs);

         /* add genvbound */
         if( addgenvbound && !SCIPisInfinity(scip, -c) )
         {
            SCIPdebugMessage("         adding genvbound\n");
            SCIP_CALL( SCIPgenVBoundAdd(scip, propdata->genvboundprop, genvboundvars, xi, genvboundcoefs, ncoefs,
                  !SCIPisPositive(scip, gamma_dual) ? 0.0 : -gamma_dual, c, bound->boundtype) );

            *found = TRUE;
         }

         /* free arrays */
         SCIPfreeBufferArray(scip, &genvboundcoefs);
         SCIPfreeBufferArray(scip, &genvboundvars);
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

/** exchange a bound which has been processed and updates the last undone and unfiltered bound index
 *  NOTE: this method has to be called after filtering or processing a bound
 */
static
void exchangeBounds(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   i                   /**< bound that was filtered or processed */
   )
{
   assert(i >= 0 && i < propdata->nbounds);
   assert(propdata->lastidx >= 0 && propdata->lastidx < propdata->nbounds);

   /* exchange the bounds */
   if( propdata->lastidx != i )
   {
      BOUND* tmp;

      tmp = propdata->bounds[i];
      propdata->bounds[i] = propdata->bounds[propdata->lastidx];
      propdata->bounds[propdata->lastidx] = tmp;
   }

   propdata->lastidx -= 1;
}

/** trying to filter some bounds using the existing LP solution */
static
SCIP_RETCODE filterExistingLP(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   int*                  nfiltered,          /**< how many bounds were filtered this round? */
   BOUND*                currbound           /**< bound for which OBBT LP was solved (Note: might be NULL) */
   )
{
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(nfiltered != NULL);

   *nfiltered = 0;

   /* only apply filtering if an LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("can't filter using existing lp solution since it was not solved to optimality\n");
      return SCIP_OKAY;
   }

   /* check if a bound is tight */
   for( i = propdata->nbounds - 1; i >= 0; --i )
   {
      BOUND* bound;                          /* shortcut for current bound */

      SCIP_Real solval;                      /* the variables value in the current solution */
      SCIP_Real boundval;                    /* current local bound for the variable */

      bound = propdata->bounds[i];
      if( bound->filtered || bound->done )
         continue;

      boundval = bound->boundtype == SCIP_BOUNDTYPE_UPPER ?
         SCIPvarGetUbLocal(bound->var) : SCIPvarGetLbLocal(bound->var);
      solval = SCIPvarGetLPSol(bound->var);

      /* bound is tight; since this holds for all fixed variables, those are filtered here automatically */
      if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGE(scip, solval, boundval))
         || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLE(scip, solval, boundval)) )
      {
         SCIP_BASESTAT basestat;

         /* mark bound as filtered */
         bound->filtered = TRUE;
         SCIPdebugMessage("trivial filtered var: %s boundval=%e solval=%e\n", SCIPvarGetName(bound->var), boundval, solval);

         /* get the basis status of the variable */
         basestat = SCIPcolGetBasisStatus(SCIPvarGetCol(bound->var));

         /* solve corresponding OBBT LP and try to generate a nontrivial genvbound */
         if( propdata->genvbdsduringfilter && currbound != NULL && basestat == SCIP_BASESTAT_BASIC )
         {
#ifndef NDEBUG
            int j;
#endif
            SCIP_Bool optimal;
            SCIP_Bool error;

            /* set objective coefficient of the bound */
            SCIP_CALL( SCIPchgVarObjProbing(scip, currbound->var, 0.0) );
            SCIP_CALL( setObjProbing(scip, propdata, bound, 1.0) );

#ifndef NDEBUG
            for( j = 0; j < SCIPgetNVars(scip); ++j )
            {
               SCIP_VAR* var;

               var = SCIPgetVars(scip)[j];
               assert(var != NULL);
               assert(SCIPisZero(scip, SCIPgetVarObjProbing(scip, var)) || var == bound->var);
            }
#endif

            /* solve the OBBT LP */
            propdata->nprobingiterations -= SCIPgetNLPIterations(scip);
            SCIP_CALL( solveLP(scip, -1, &error, &optimal) );
            propdata->nprobingiterations += SCIPgetNLPIterations(scip);
            assert(propdata->nprobingiterations >= 0);

            /* try to generate a genvbound if we have solved the OBBT LP */
            if( optimal && propdata->genvboundprop != NULL
                  && (SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1)) )
            {
               SCIP_Bool found;

               assert(!error);
               SCIP_CALL( createGenVBound(scip, propdata, bound, &found) );

               if( found )
               {
                  propdata->ngenvboundstrivfil += 1;
                  SCIPdebugMessage("found genvbound during trivial filtering\n");
               }
            }

            /* restore objective function */
            SCIP_CALL( setObjProbing(scip, propdata, bound, 0.0) );
            SCIP_CALL( setObjProbing(scip, propdata, currbound, 1.0) );
         }

         /* exchange bound i with propdata->bounds[propdata->lastidx] */
         if( propdata->lastidx >= 0 )
            exchangeBounds(propdata, i);

         /* increase number of filtered variables */
         (*nfiltered)++;
      }
   }

   return SCIP_OKAY;
}

/** enforces one round of filtering */
static
SCIP_RETCODE filterRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   int                   itlimit,            /**< LP iteration limit (-1: no limit) */
   int*                  nfiltered,          /**< how many bounds were filtered this round */
   SCIP_Real*            objcoefs,           /**< array to store the nontrivial objective coefficients */
   int*                  objcoefsinds,       /**< array to store bound indices for which their corresponding variables
                                               *  has a nontrivial objective coefficient */
   int                   nobjcoefs           /**< number of nontrivial objective coefficients */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   SCIP_Bool error;
   SCIP_Bool optimal;

   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);
   assert(nfiltered != NULL);
   assert(objcoefs != NULL);
   assert(objcoefsinds != NULL);
   assert(nobjcoefs >= 0);

   *nfiltered = 0;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* solve LP */
   propdata->nfilterlpiters -= (int) SCIPgetNLPIterations(scip);
   SCIP_CALL( solveLP(scip, itlimit, &error, &optimal) );
   propdata->nfilterlpiters += (int) SCIPgetNLPIterations(scip);
   assert(propdata->nfilterlpiters >= 0);

   if( !optimal )
   {
      SCIPdebugMessage("skipping filter round since the LP was not solved to optimality\n");
      return SCIP_OKAY;
   }

   assert(!error);

   /* check if a bound is tight */
   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;                          /* shortcut for current bound */

      SCIP_Real solval;                      /* the variables value in the current solution */
      SCIP_Real boundval;                    /* current local bound for the variable */

      bound = propdata->bounds[i];

      /* if bound is filtered it was handled already before */
      if( bound->filtered )
         continue;

      boundval = bound->boundtype == SCIP_BOUNDTYPE_UPPER ?
         SCIPvarGetUbLocal(bound->var) : SCIPvarGetLbLocal(bound->var);
      solval = SCIPvarGetLPSol(bound->var);

      /* bound is tight */
      if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGE(scip, solval, boundval))
         || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLE(scip, solval, boundval)) )
      {
         SCIP_Real objcoef;
         SCIP_BASESTAT basestat;

         /* mark bound as filtered */
         bound->filtered = TRUE;

         /* get the basis status of the variable */
         basestat = SCIPcolGetBasisStatus(SCIPvarGetCol(bound->var));

         /* increase number of filtered variables */
         (*nfiltered)++;

         /* solve corresponding OBBT LP and try to generate a nontrivial genvbound */
         if( propdata->genvbdsduringfilter && basestat == SCIP_BASESTAT_BASIC )
         {
            int j;

            /* set all objective coefficients to zero */
            for( j = 0; j < nobjcoefs; ++j )
            {
               BOUND* filterbound;

               filterbound = propdata->bounds[ objcoefsinds[j] ];
               assert(filterbound != NULL);

               SCIP_CALL( SCIPchgVarObjProbing(scip, filterbound->var, 0.0) );
            }

#ifndef NDEBUG
            for( j = 0; j < nvars; ++j )
               assert(SCIPisZero(scip, SCIPgetVarObjProbing(scip, vars[j])));
#endif

            /* set objective coefficient of the bound */
            SCIP_CALL( setObjProbing(scip, propdata, bound, 1.0) );

            /* solve the OBBT LP */
            propdata->nfilterlpiters -= (int) SCIPgetNLPIterations(scip);
            SCIP_CALL( solveLP(scip, -1, &error, &optimal) );
            propdata->nfilterlpiters += (int) SCIPgetNLPIterations(scip);
            assert(propdata->nfilterlpiters >= 0);

            /* try to generate a genvbound if we have solved the OBBT LP */
            if( optimal && propdata->genvboundprop != NULL
                  && (SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1)) )
            {
               SCIP_Bool found;

               assert(!error);
               SCIP_CALL( createGenVBound(scip, propdata, bound, &found) );

               if( found )
               {
                  propdata->ngenvboundsaggrfil += 1;
                  SCIPdebugMessage("found genvbound during aggressive filtering\n");
               }

            }

            /* restore objective function */
            for( j = 0; j < nobjcoefs; ++j )
            {
               BOUND* filterbound;

               filterbound = propdata->bounds[ objcoefsinds[j] ];
               assert(filterbound != NULL);

               /* NOTE: only restore coefficients of nonfiltered bounds */
               if( !filterbound->filtered )
               {
                  assert(!SCIPisZero(scip, objcoefs[j]));
                  SCIP_CALL( SCIPchgVarObjProbing(scip, propdata->bounds[ objcoefsinds[j] ]->var, objcoefs[j]) );
               }
            }
         }

         /* get the corresponding variable's objective coefficient */
         objcoef = SCIPgetVarObjProbing(scip, bound->var);

         /* change objective coefficient if it was set up for this bound */
          if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisNegative(scip, objcoef))
             || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisPositive(scip, objcoef)) )
          {
             SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, 0.0) );
          }
      }
   }

   return SCIP_OKAY;
}

/** filter some bounds that are not improvable by solving auxiliary LPs */
static
SCIP_RETCODE filterBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint          itlimit             /**< LP iteration limit (-1: no limit) */
   )
{
   SCIP_VAR** vars;
   SCIP_Longint nolditerations;
   SCIP_Real* objcoefs;               /* array to store the nontrivial objective coefficients */
   int* objcoefsinds;                 /* array to store bound indices for which the corresponding variable
                                       * has a nontrivial objective coefficient */
   int nobjcoefs;                     /* number of nontrivial objective coefficients */
   int nleftiterations;
   int i;
   int nfiltered;
   int ntotalfiltered;
   int nvars;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   ntotalfiltered = 0;
   nolditerations = SCIPgetNLPIterations(scip);
   nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIPdebugMessage("start filter rounds\n");

   SCIP_CALL( SCIPallocBufferArray(scip, &objcoefs, propdata->nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objcoefsinds, propdata->nbounds) );
   nobjcoefs = 0;

   /*
    * 1.) Try first to filter lower bounds of interesting variables, whose bounds are not already filtered
    */

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[i], 0.0) );
   }

   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->boundtype == SCIP_BOUNDTYPE_LOWER && !propdata->bounds[i]->filtered
            && !propdata->bounds[i]->done )
      {
         SCIP_Real objcoef;

         objcoef = getFilterCoef(scip, propdata, propdata->bounds[i]->var, SCIP_BOUNDTYPE_LOWER);

         if( !SCIPisZero(scip, objcoef) )
         {
            SCIP_CALL( SCIPchgVarObjProbing(scip, propdata->bounds[i]->var, objcoef) );

            /* store nontrivial objective coefficients */
            objcoefs[nobjcoefs] = objcoef;
            objcoefsinds[nobjcoefs] = i;
            ++nobjcoefs;
         }
      }
   }

   do
   {
      SCIPdebugMessage("doing a lower bounds round\n");
      SCIP_CALL( filterRound(scip, propdata, nleftiterations, &nfiltered, objcoefs, objcoefsinds, nobjcoefs) );
      ntotalfiltered += nfiltered;
      SCIPdebugMessage("filtered %d more bounds in lower bounds round\n", nfiltered);

      /* update iterations left */
      nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   }
   while( nfiltered >= propdata->nminfilter && ( nleftiterations == -1 ||  nleftiterations > 0 ) );

   /*
    * 2.) Now try to filter the remaining upper bounds of interesting variables, whose bounds are not already filtered
    */

   /* set all objective coefficients to zero */
   for( i = 0; i < nobjcoefs; i++ )
   {
      BOUND* bound;

      assert(objcoefsinds[i] >= 0 && objcoefsinds[i] < propdata->nbounds);
      bound = propdata->bounds[ objcoefsinds[i] ];
      assert(bound != NULL);
      SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, 0.0) );
   }

   /* reset number of nontrivial objective coefficients */
   nobjcoefs = 0;

#ifndef NDEBUG
   for( i = 0; i < nvars; ++i )
      assert(SCIPisZero(scip, SCIPgetVarObjProbing(scip, vars[i])));
#endif

   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->boundtype == SCIP_BOUNDTYPE_UPPER && !propdata->bounds[i]->filtered )
      {
         SCIP_Real objcoef;

         objcoef = getFilterCoef(scip, propdata, propdata->bounds[i]->var, SCIP_BOUNDTYPE_UPPER);

         if( !SCIPisZero(scip, objcoef) )
         {
            SCIP_CALL( SCIPchgVarObjProbing(scip, propdata->bounds[i]->var, objcoef) );

            /* store nontrivial objective coefficients */
            objcoefs[nobjcoefs] = objcoef;
            objcoefsinds[nobjcoefs] = i;
            ++nobjcoefs;
         }
      }
   }

   do
   {
      SCIPdebugMessage("doing an upper bounds round\n");
      SCIP_CALL( filterRound(scip, propdata, nleftiterations, &nfiltered, objcoefs, objcoefsinds, nobjcoefs) );
      SCIPdebugMessage("filtered %d more bounds in upper bounds round\n", nfiltered);
      ntotalfiltered += nfiltered;
      /* update iterations left */
      nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   }
   while( nfiltered >= propdata->nminfilter && ( nleftiterations == -1 ||  nleftiterations > 0 ) );

   SCIPdebugMessage("filtered %d this round\n", ntotalfiltered);
   propdata->nfiltered += ntotalfiltered;

   /* free array */
   SCIPfreeBufferArray(scip, &objcoefsinds);
   SCIPfreeBufferArray(scip, &objcoefs);

   return SCIP_OKAY;
}

/** applies possible bound changes that were found */
static
SCIP_RETCODE applyBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
#ifdef SCIP_DEBUG
   int ntightened;                           /* stores the number of successful bound changes */
#endif
   int i;

   assert(scip != NULL);
   assert(!SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND);

   SCIPdebug( ntightened = 0 );

   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;                          /* shortcut to the current bound */
      SCIP_Bool infeas;                      /* stores wether a tightening approach forced an infeasibilty */
      SCIP_Bool tightened;                   /* stores wether a tightening approach was successful */

      bound = propdata->bounds[i];

      if( bound->found )
      {
         SCIPdebug( double oldbound = (bound->boundtype == SCIP_BOUNDTYPE_LOWER)
            ? SCIPvarGetLbLocal(bound->var)
            : SCIPvarGetUbLocal(bound->var) );

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
            SCIPdebug( SCIPdebugMessage("tightended: %s old: %e new: %e\n" , SCIPvarGetName(bound->var), oldbound,
                  bound->newval) );
            *result = SCIP_REDUCEDDOM;
            SCIPdebug( ntightened++ );
         }
      }
   }

   SCIPdebug( SCIPdebugMessage("tightened bounds: %d\n", ntightened) );

   return SCIP_OKAY;
}

/** tries to tighten a bound in probing mode  */
static
SCIP_RETCODE tightenBoundProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   BOUND*                bound,              /**< bound that could be tightened */
   SCIP_Real             newval,             /**< new bound value */
   SCIP_Bool*            tightened           /**< was tightening successful? */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(bound != NULL);
   assert(tightened != NULL);

   *tightened = FALSE;

   /* get old bounds */
   lb = SCIPvarGetLbLocal(bound->var);
   ub = SCIPvarGetUbLocal(bound->var);

   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      /* round bounds new value if variable is integral */
      if( SCIPvarIsIntegral(bound->var) )
         newval = SCIPceil(scip, newval);

      /* ensure that we give consistent bounds to the LP solver */
      if( newval > ub )
         newval = ub;

      /* tighten if really better */
      if( SCIPisLbBetter(scip, newval, lb, ub) )
      {
         SCIP_CALL( SCIPchgVarLbProbing(scip, bound->var, newval) );
         *tightened = TRUE;
      }
   }
   else
   {
      /* round bounds new value if variable is integral */
      if( SCIPvarIsIntegral(bound->var) )
         newval = SCIPfloor(scip, newval);

      /* ensure that we give consistent bounds to the LP solver */
      if( newval < lb )
         newval = lb;

      /* tighten if really better */
      if( SCIPisUbBetter(scip, newval, lb, ub) )
      {
         SCIP_CALL( SCIPchgVarUbProbing(scip, bound->var, newval) );
         *tightened = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** comparison method for two bounds w.r.t. their scores */
static
SCIP_DECL_SORTPTRCOMP(compBoundsScore)
{
   BOUND* bound1 = (BOUND*) elem1;
   BOUND* bound2 = (BOUND*) elem2;

   return bound1->score == bound2->score ? 0 : ( bound1->score > bound2->score ? 1 : -1 );
}

/** comparison method for two bounds w.r.t. their boundtype */
static
SCIP_DECL_SORTPTRCOMP(compBoundsBoundtype)
{
   int diff;
   BOUND* bound1 = (BOUND*) elem1;
   BOUND* bound2 = (BOUND*) elem2;

   /* prioritize undone bounds */
   diff = (!bound1->done ? 1 : 0) - (!bound2->done ? 1 : 0);
   if( diff != 0 )
      return diff;

   /* prioritize unfiltered bounds */
   diff = (!bound1->filtered ? 1 : 0) - (!bound2->filtered ? 1 : 0);
   if( diff != 0 )
      return diff;

   diff = (bound1->boundtype == SCIP_BOUNDTYPE_LOWER ? 1 : 0) - (bound2->boundtype == SCIP_BOUNDTYPE_LOWER ? 1 : 0);

   if( diff == 0 )
      return (bound1->score == bound2->score) ? 0 : (bound1->score > bound2->score ? 1 : -1);
   else
      return diff;
}

/** sort the propdata->bounds array with their distance or their boundtype key */
static
SCIP_RETCODE sortBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   SCIPdebugMessage("sort bounds\n");
   SCIPsortDownPtr((void**) propdata->bounds, compBoundsBoundtype, propdata->nbounds);

   return SCIP_OKAY;
}

/** evaluates a bound for the current LP solution */
static
SCIP_Real evalBound(
   SCIP*                 scip,
   BOUND*                bound
   )
{
   assert(scip != NULL);
   assert(bound != NULL);

   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
      return REALABS( SCIPvarGetLPSol(bound->var) - SCIPvarGetLbLocal(bound->var) );
   else
      return REALABS( SCIPvarGetUbLocal(bound->var) - SCIPvarGetLPSol(bound->var) );
}

/** returns the index of the next undone and unfiltered bound with the smallest distance */
static
int nextBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Bool             convexphase         /**< consider only convex variables? */
   )
{
   SCIP_Real bestval;
   int bestidx;
   int k;

   assert(scip != NULL);
   assert(propdata != NULL);

   bestidx = -1;
   bestval = SCIPinfinity(scip);

   for( k = 0; k <= propdata->lastidx; ++k )
   {
      BOUND* tmpbound;
      tmpbound = propdata->bounds[k];

      assert(tmpbound != NULL);

      if( !tmpbound->filtered && !tmpbound->done && (tmpbound->nonconvex == !convexphase) )
      {
         SCIP_Real boundval;

         /* return the next bound which is not done or unfiltered yet */
         if( propdata->orderingalgo == 0 )
            return k;

         boundval = evalBound(scip, tmpbound);

         /* negate boundval if we use the reverse greedy algorithm */
         boundval = (propdata->orderingalgo == 2) ? -1.0 * boundval : boundval;

         if( bestidx == -1 || boundval < bestval )
         {
            bestidx = k;
            bestval = boundval;
         }
      }
   }

   return bestidx;
}

/** try to separate the solution of the last OBBT LP in order to learn better variable bounds; we apply additional
 *  separation rounds as long as the routine finds better bounds; because of dual degeneracy we apply a minimum number of
 *  separation rounds
 */
static
SCIP_RETCODE applySeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   BOUND*                currbound,          /**< current bound */
   SCIP_Longint*         nleftiterations,    /**< number of left iterations (-1 for no limit) */
   SCIP_Bool*            success             /**< pointer to store if we have found a better bound */
   )
{
   SCIP_Bool inroot;
   int i;

   assert(nleftiterations != NULL);
   assert(success != NULL);
   assert(SCIPinProbing(scip));

   *success = FALSE;

   /* check if we are originally in the root node */
   inroot = SCIPgetDepth(scip) == 1;

   for( i = 0; i <= propdata->sepamaxiter; ++i )
   {
      SCIP_Longint nlpiter;
      SCIP_Real oldval;
      SCIP_Bool cutoff;
      SCIP_Bool delayed;
      SCIP_Bool error;
      SCIP_Bool optimal;
      SCIP_Bool tightened;

      oldval = SCIPvarGetLPSol(currbound->var);

      /* find and store cuts to separate the current LP solution */
      SCIP_CALL( SCIPseparateSol(scip, NULL, inroot, FALSE, &delayed, &cutoff) );
      SCIPdebugMessage("applySeparation() - ncuts = %d\n", SCIPgetNCuts(scip));

      /* leave if we did not found any cut */
      if( SCIPgetNCuts(scip) == 0 )
         break;

      /* apply cuts and resolve LP */
      SCIP_CALL( SCIPapplyCutsProbing(scip, &cutoff) );
      assert(SCIPgetNCuts(scip) == 0);
      nlpiter = SCIPgetNLPIterations(scip);
      SCIP_CALL( solveLP(scip, (int) *nleftiterations, &error, &optimal) );
      nlpiter = SCIPgetNLPIterations(scip) - nlpiter;
      SCIPdebugMessage("applySeparation() - optimal=%u error=%u lpiter=%" SCIP_LONGINT_FORMAT "\n", optimal, error, nlpiter);
      SCIPdebugMessage("oldval = %e newval = %e\n", oldval, SCIPvarGetLPSol(currbound->var));

      /* leave if we did not solve the LP to optimality or an error occured */
      if( error || !optimal )
         break;

      /* try to generate a genvbound */
      if( inroot && propdata->genvboundprop != NULL && propdata->genvbdsduringsepa )
      {
         SCIP_Bool found;
         SCIP_CALL( createGenVBound(scip, propdata, currbound, &found) );
         propdata->ngenvboundsprobing += found ? 1 : 0;
      }

      /* try to tight the variable bound */
      tightened = FALSE;
      if( !SCIPisEQ(scip, oldval, SCIPvarGetLPSol(currbound->var)) )
      {
         SCIP_CALL( tightenBoundProbing(scip, currbound, SCIPvarGetLPSol(currbound->var), &tightened) );
         SCIPdebugMessage("apply separation - tightened=%u oldval=%e newval=%e\n", tightened, oldval,
            SCIPvarGetLPSol(currbound->var));

         *success |= tightened;
      }

      /* leave the separation if we did not tighten the bound and proceed at least propdata->sepaminiter iterations */
      if( !tightened && i >= propdata->sepaminiter )
         break;
   }

   return SCIP_OKAY;
}

/** finds new variable bounds until no iterations left or all bounds have been checked */
static
SCIP_RETCODE findNewBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint*         nleftiterations,    /**< pointer to store the number of left iterations */
   SCIP_Bool             convexphase         /**< consider only convex variables? */
   )
{
   SCIP_Longint nolditerations;
   SCIP_Bool iterationsleft;
   BOUND* currbound;
   SCIP_Longint itlimit;
   int nextboundidx;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(nleftiterations != NULL);

   /* update the number of left iterations */
   nolditerations = SCIPgetNLPIterations(scip);
   itlimit = *nleftiterations;
   assert(*nleftiterations == getIterationsLeft(scip, nolditerations, itlimit));
   iterationsleft = (*nleftiterations == -1) || (*nleftiterations > 0);

   /* To improve the performance we sort the bound in such a way that the undone and
    * unfiltered bounds are at the end of propdata->bounds. We calculate and update
    * the position of the last unfiltered and undone bound in propdata->lastidx
    */
   if( !convexphase )
   {
      /* sort bounds */
      SCIP_CALL( sortBounds(scip, propdata) );

      /* if the first bound is filtered or done then there is no bound left */
      if( propdata->bounds[0]->done || propdata->bounds[0]->filtered )
      {
         SCIPdebugMessage("no unprocessed/unfiltered bound left\n");
         return SCIP_OKAY;
      }

      /* compute the last undone and unfiltered node */
      propdata->lastidx = 0;
      while( propdata->lastidx < propdata->nbounds - 1 && !propdata->bounds[propdata->lastidx]->done &&
            !propdata->bounds[propdata->lastidx]->filtered )
         ++propdata->lastidx;

      SCIPdebugMessage("lastidx = %d\n", propdata->lastidx);
   }

   /* find the first unprocessed bound */
   nextboundidx = nextBound(scip, propdata, convexphase);

   /* skip if there is no bound left */
   if( nextboundidx == -1 )
   {
      SCIPdebugMessage("no unprocessed/unfiltered bound left\n");
      return SCIP_OKAY;
   }

   currbound = propdata->bounds[nextboundidx];
   assert(!currbound->done && !currbound->filtered);

   /* main loop */
   while( iterationsleft &&  !SCIPisStopped(scip) )
   {
      SCIP_Bool optimal;
      SCIP_Bool error;
      int nfiltered;

      assert(currbound != NULL);
      assert(currbound->done == FALSE);
      assert(currbound->filtered == FALSE);

      /* do not visit currbound more than once */
      currbound->done = TRUE;
      exchangeBounds(propdata, nextboundidx);

      /* set objective for curr */
      SCIP_CALL( setObjProbing(scip, propdata, currbound, 1.0) );

      SCIPdebugMessage("before solving      Boundtype: %d , LB: %e , UB: %e\n",
         currbound->boundtype == SCIP_BOUNDTYPE_LOWER, SCIPvarGetLbLocal(currbound->var),
         SCIPvarGetUbLocal(currbound->var) );
      SCIPdebugMessage("before solving      var <%s>, LP value: %f\n",
         SCIPvarGetName(currbound->var), SCIPvarGetLPSol(currbound->var));

      SCIPdebugMessage("probing iterations before solve: %lld \n", SCIPgetNLPIterations(scip));

      propdata->nprobingiterations -= SCIPgetNLPIterations(scip);

      /* now solve the LP */
      SCIP_CALL( solveLP(scip, (int) *nleftiterations, &error, &optimal) );

      propdata->nprobingiterations += SCIPgetNLPIterations(scip);
      propdata->nsolvedbounds++;

      SCIPdebugMessage("probing iterations after solve: %lld \n", SCIPgetNLPIterations(scip));
      SCIPdebugMessage("OPT: %u ERROR: %u\n" , optimal, error);
      SCIPdebugMessage("after solving      Boundtype: %d , LB: %e , UB: %e\n",
         currbound->boundtype == SCIP_BOUNDTYPE_LOWER, SCIPvarGetLbLocal(currbound->var),
         SCIPvarGetUbLocal(currbound->var) );
      SCIPdebugMessage("after solving      var <%s>, LP value: %f\n",
         SCIPvarGetName(currbound->var), SCIPvarGetLPSol(currbound->var));

      /* update nleftiterations */
      *nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
      iterationsleft = (*nleftiterations == -1) || (*nleftiterations > 0);

      if( error )
      {
         SCIPdebugMessage("ERROR during LP solving\n");

         /* set the objective of currbound to zero to null the whole objective; otherwise the objective is wrong when
          * we call findNewBounds() for the convex phase
          */
         SCIP_CALL( SCIPchgVarObjProbing(scip, currbound->var, 0.0) );

         return SCIP_OKAY;
      }

      if( optimal )
      {
         SCIP_Bool success;

         currbound->newval = SCIPvarGetLPSol(currbound->var);
         currbound->found = TRUE;

         /* in root node we may want to create a genvbound (independent of tightening success) */
         if( (SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1))
               && propdata->genvboundprop != NULL )
         {
            SCIP_Bool found;

            SCIP_CALL( createGenVBound(scip, propdata, currbound, &found) );

            if( found )
               propdata->ngenvboundsprobing += 1;
         }

         /* try to tighten bound in probing mode */
         success = FALSE;
         if( propdata->tightintboundsprobing && SCIPvarIsIntegral(currbound->var) )
         {
            SCIPdebugMessage("tightening bound %s = %e bounds: [%e, %e]\n", SCIPvarGetName(currbound->var),
                currbound->newval, SCIPvarGetLbLocal(currbound->var), SCIPvarGetUbLocal(currbound->var) );
            SCIP_CALL( tightenBoundProbing(scip, currbound, currbound->newval, &success) );
            SCIPdebugMessage("tightening bound %s\n", success ? "successful" : "not successful");
         }
         else if( propdata->tightcontboundsprobing && !SCIPvarIsIntegral(currbound->var) )
         {
            SCIPdebugMessage("tightening bound %s = %e bounds: [%e, %e]\n", SCIPvarGetName(currbound->var),
               currbound->newval, SCIPvarGetLbLocal(currbound->var), SCIPvarGetUbLocal(currbound->var) );
            SCIP_CALL( tightenBoundProbing(scip, currbound, currbound->newval, &success) );
            SCIPdebugMessage("tightening bound %s\n", success ? "successful" : "not successful");
         }

         /* separate current OBBT LP solution */
         if( iterationsleft && propdata->separatesol )
         {
            propdata->nprobingiterations -= SCIPgetNLPIterations(scip);
            SCIP_CALL( applySeparation(scip, propdata, currbound, nleftiterations, &success) );
            propdata->nprobingiterations += SCIPgetNLPIterations(scip);

            /* remember best solution value after solving additional separations LPs */
            if( success )
            {
#ifndef NDEBUG
               SCIP_Real newval = SCIPvarGetLPSol(currbound->var);

               /* round new bound if the variable is integral */
               if( SCIPvarIsIntegral(currbound->var) )
                  newval = currbound->boundtype == SCIP_BOUNDTYPE_LOWER ?
                     SCIPceil(scip, newval) : SCIPfloor(scip, newval);

               assert((currbound->boundtype == SCIP_BOUNDTYPE_LOWER &&
                     SCIPisGT(scip, newval, currbound->newval))
                  || (currbound->boundtype == SCIP_BOUNDTYPE_UPPER &&
                     SCIPisLT(scip, newval, currbound->newval)));
#endif

               currbound->newval = SCIPvarGetLPSol(currbound->var);
            }
         }

         /* filter bound candidates by using the current LP solution */
         if( propdata->applytrivialfilter )
         {
            SCIP_CALL( filterExistingLP(scip, propdata, &nfiltered, currbound) );
            SCIPdebugMessage("filtered %d bounds via inspecting present LP solution\n", nfiltered);
            propdata->ntrivialfiltered += nfiltered;
         }

         propdata->propagatecounter += success ? 1 : 0;

         /* propagate if we have found enough bound tightenings */
         if( propdata->propagatefreq != 0 && propdata->propagatecounter >= propdata->propagatefreq )
         {
            SCIP_Longint ndomredsfound;
            SCIP_Bool cutoff;

            SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, &ndomredsfound) );
            SCIPdebugMessage("propagation - cutoff %u  ndomreds %" SCIP_LONGINT_FORMAT "\n", cutoff, ndomredsfound);

            propdata->npropagatedomreds += ndomredsfound;
            propdata->propagatecounter = 0;
         }
      }

      /* set objective to zero */
      SCIP_CALL( setObjProbing(scip, propdata, currbound, 0.0) );

      /* find the first unprocessed bound */
      nextboundidx = nextBound(scip, propdata, convexphase);

      /* check if there is no unprocessed and unfiltered node left */
      if( nextboundidx == -1 )
      {
         SCIPdebugMessage("NO unvisited/unfiltered bound left!\n");
         break;
      }

      currbound = propdata->bounds[nextboundidx];
      assert(!currbound->done && !currbound->filtered);
   }

   if( iterationsleft )
   {
      SCIPdebugMessage("still iterations left: %" SCIP_LONGINT_FORMAT "\n", *nleftiterations);
   }
   else
   {
      SCIPdebugMessage("no iterations left\n");
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
   SCIP_VAR** vars;
   SCIP_Real* oldlbs;
   SCIP_Real* oldubs;
   SCIP_Longint lastnpropagatedomreds;
   SCIP_Longint nleftiterations;
   SCIP_Real oldconditionlimit;
   SCIP_Real oldboundstreps;
   SCIP_Real olddualfeastol;
   SCIP_Bool hasconditionlimit;
   SCIP_Bool continuenode;
   SCIP_Bool boundleft;
   int nfiltered;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   SCIPdebugMessage("apply obbt\n");

   oldlbs = NULL;
   oldubs = NULL;
   lastnpropagatedomreds = propdata->npropagatedomreds;
   nleftiterations = itlimit;
   continuenode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode;
   propdata->lastidx = -1;
   boundleft = FALSE;
   *result = SCIP_DIDNOTFIND;

   /* store old variable bounds if we use propagation during obbt */
   if( propdata->propagatefreq > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &oldlbs, propdata->nbounds) );
      SCIP_CALL( SCIPallocBufferArray(scip, &oldubs, propdata->nbounds) );
   }

   /* reset bound data structure flags; fixed variables are marked as filtered */
   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound = propdata->bounds[i];
      bound->found = FALSE;

      /* store old variable bounds */
      if( oldlbs != NULL && oldubs != NULL )
      {
         oldlbs[bound->index] = SCIPvarGetLbLocal(bound->var);
         oldubs[bound->index] = SCIPvarGetUbLocal(bound->var);
      }

      /* reset 'done' and 'filtered' flag in a new B&B node */
      if( !continuenode )
      {
         bound->done = FALSE;
         bound->filtered = FALSE;
      }

      /* mark fixed variables as filtered */
      bound->filtered |= varIsFixedLocal(scip, bound->var);

      /* check for an unprocessed bound */
      if( !bound->filtered && !bound->done )
         boundleft = TRUE;
   }

   /* no bound left to check */
   if( !boundleft )
      goto TERMINATE;

   /* filter variables via inspecting present LP solution */
   if( propdata->applytrivialfilter && !continuenode )
   {
      SCIP_CALL( filterExistingLP(scip, propdata, &nfiltered, NULL) );
      SCIPdebugMessage("filtered %d bounds via inspecting present LP solution\n", nfiltered);
      propdata->ntrivialfiltered += nfiltered;
   }

   /* store old dualfeasibletol */
   olddualfeastol = SCIPdualfeastol(scip);

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );
   SCIPdebugMessage("start probing\n");

   /* tighten dual feastol */
   if( propdata->dualfeastol < olddualfeastol )
   {
      SCIP_CALL( SCIPchgDualfeastol(scip, propdata->dualfeastol) );
   }

   /* tighten condition limit */
   hasconditionlimit = (SCIPgetRealParam(scip, "lp/conditionlimit", &oldconditionlimit) == SCIP_OKAY);
   if( !hasconditionlimit )
   {
      SCIPwarningMessage(scip, "obbt propagator could not set condition limit in LP solver - running without\n");
   }
   else if( propdata->conditionlimit > 0.0 && (oldconditionlimit < 0.0 || propdata->conditionlimit < oldconditionlimit) )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "lp/conditionlimit", propdata->conditionlimit) );
   }

   /* tighten relative bound improvement limit */
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/boundstreps", &oldboundstreps) );
   if( !SCIPisEQ(scip, oldboundstreps, propdata->boundstreps) )
   {
     SCIP_CALL( SCIPsetRealParam(scip, "numerics/boundstreps", propdata->boundstreps) );
   }

   /* add objective cutoff */
   SCIP_CALL( addObjCutoff(scip, propdata) );

   /* apply filtering */
   if( propdata->applyfilterrounds )
   {
      SCIP_CALL( filterBounds(scip, propdata, nleftiterations) );
   }

   /* set objective coefficients to zero */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[i], 0.0) );
   }

   /* find new bounds for the variables */
   SCIP_CALL( findNewBounds(scip, propdata, &nleftiterations, FALSE) );

   if( nleftiterations > 0 || itlimit < 0 )
   {
      SCIP_CALL( findNewBounds(scip, propdata, &nleftiterations, TRUE) );
   }

   /* reset dual feastol and condition limit */
   SCIP_CALL( SCIPchgDualfeastol(scip, olddualfeastol) );
   if( hasconditionlimit )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "lp/conditionlimit", oldconditionlimit) );
   }

   /* update bound->newval if we have learned additional bound tightenings during SCIPpropagateProbing() */
   if( oldlbs != NULL && oldubs != NULL && propdata->npropagatedomreds - lastnpropagatedomreds > 0 )
   {
      assert(propdata->propagatefreq > 0);
      for( i = 0; i < propdata->nbounds; ++i )
      {
         BOUND* bound = propdata->bounds[i];

         /* it might be the case that a bound found by the additional propagation is better than the bound found after solving an OBBT
          * LP
          */
         if( bound->found )
         {
            if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
               bound->newval = MAX(bound->newval, SCIPvarGetLbLocal(bound->var)); /*lint !e666*/
            else
               bound->newval = MIN(bound->newval, SCIPvarGetUbLocal(bound->var)); /*lint !e666*/
         }
         else
         {
            SCIP_Real oldlb;
            SCIP_Real oldub;

            oldlb = oldlbs[bound->index];
            oldub = oldubs[bound->index];

            if( bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisLbBetter(scip, SCIPvarGetLbLocal(bound->var), oldlb, oldub) )
            {
               SCIPdebugMessage("tighter lower bound due to propagation: %d - %e -> %e\n", i, oldlb, SCIPvarGetLbLocal(bound->var));
               bound->newval = SCIPvarGetLbLocal(bound->var);
               bound->found = TRUE;
            }

            if( bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisUbBetter(scip, SCIPvarGetUbLocal(bound->var), oldlb, oldub) )
            {
               SCIPdebugMessage("tighter upper bound due to propagation: %d - %e -> %e\n", i, oldub, SCIPvarGetUbLocal(bound->var));
               bound->newval = SCIPvarGetUbLocal(bound->var);
               bound->found = TRUE;
            }
         }
      }
   }

   /* reset relative bound improvement limit */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/boundstreps", oldboundstreps) );

   /* end probing */
   SCIP_CALL( SCIPendProbing(scip) );
   SCIPdebugMessage("end probing!\n");

   /* release cutoff row if there is one */
   if( propdata->cutoffrow != NULL )
   {
      assert(!SCIProwIsInLP(propdata->cutoffrow));
      SCIP_CALL( SCIPreleaseRow(scip, &(propdata->cutoffrow)) );
   }

   /* apply buffered bound changes */
   SCIP_CALL( applyBoundChgs(scip, propdata, result) );

TERMINATE:
   SCIPfreeBufferArrayNull(scip, &oldubs);
   SCIPfreeBufferArrayNull(scip, &oldlbs);

   return SCIP_OKAY;
}

/** computes the score of a bound */
static
unsigned int getScore(
   SCIP*                 scip,               /**< SCIP data structure */
   BOUND*                bound,              /**< pointer of bound */
   int                   nlcount,            /**< number of nonlinear constraints containing the bounds variable */
   int                   maxnlcount          /**< maximal number of nonlinear constraints a variable appears in */
   )
{
   unsigned int score;                       /* score to be computed */

   assert(scip != NULL);
   assert(bound != NULL);
   assert(nlcount >= 0);
   assert(maxnlcount >= nlcount);

   /* score = ( nlcount * ( BASE - 1 ) / maxnlcount ) * BASE^2 + vartype * BASE + boundtype */
   score = (unsigned int) ( nlcount > 0 ? (OBBT_SCOREBASE * nlcount * ( OBBT_SCOREBASE - 1 )) / maxnlcount : 0 );
   switch( SCIPvarGetType(bound->var) )
   {
   case SCIP_VARTYPE_INTEGER:
      score += 1;
      break;
   case SCIP_VARTYPE_IMPLINT:
      score += 2;
      break;
   case SCIP_VARTYPE_CONTINUOUS:
      score += 3;
      break;
   case SCIP_VARTYPE_BINARY:
      score += 4;
      break;
   default:
      break;
   }

   score *= OBBT_SCOREBASE;
   if( bound->boundtype == SCIP_BOUNDTYPE_UPPER )
      score += 1;

   return score;
}

/** count the variables which appear in non-convex term of nlrow  */
static
SCIP_RETCODE countNLRowVarsNonConvexity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcounts,           /**< store the number each variable appears in a
                                              *   non-convex term */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   int t;
   int nexprtreevars;
   SCIP_VAR** exprtreevars;
   SCIP_EXPRTREE* exprtree;

   assert(scip != NULL);
   assert(nlcounts != NULL);
   assert(nlrow != NULL);

   /* go through all quadratic terms */
   for( t = SCIPnlrowGetNQuadElems(nlrow) - 1; t >= 0; --t )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* bilinvar1;
      SCIP_VAR* bilinvar2;

      /* get quadratic term */
      quadelem = &SCIPnlrowGetQuadElems(nlrow)[t];

      /* get involved variables */
      bilinvar1 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      bilinvar2 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];

      assert(bilinvar1 != NULL);
      assert(bilinvar2 != NULL);

      /* we have a non-convex square term */
      if( bilinvar1 == bilinvar2 && !(quadelem->coef >= 0 ? SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)) : SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow))) )
      {
         ++nlcounts[SCIPvarGetProbindex(bilinvar1)];
         ++nlcounts[SCIPvarGetProbindex(bilinvar2)];
      }

      /* bilinear terms are in general non-convex */
      if( bilinvar1 != bilinvar2 )
      {
         ++nlcounts[SCIPvarGetProbindex(bilinvar1)];
         ++nlcounts[SCIPvarGetProbindex(bilinvar2)];
      }
   }

   exprtree = SCIPnlrowGetExprtree(nlrow);
   if( exprtree != NULL )
   {
      nexprtreevars = SCIPexprtreeGetNVars(exprtree);
      exprtreevars = SCIPexprtreeGetVars(exprtree);

      /* assume that the expression tree represents a non-convex constraint */
      for( t = 0; t < nexprtreevars; ++t)
      {
         SCIP_VAR* var;
         var = exprtreevars[t];
         assert(var != NULL);

         ++nlcounts[SCIPvarGetProbindex(var)];
      }
   }

   return SCIP_OKAY;
}

/** count how often each variable appears in a non-convex term */
static
SCIP_RETCODE getNLPVarsNonConvexity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcounts            /**< store the number each variable appears in a
                                              *   non-convex term */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nvars;
   int nconss;
   int i;

   assert(scip != NULL);
   assert(nlcounts != NULL);

   nvars = SCIPgetNVars(scip);
   BMSclearMemoryArray(nlcounts, nvars);

   /* quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL )
   {

      /*SCIPdebugMessage("cons_quadratic is there!\n");*/
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMessage("nconss(quadratic) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         /* only check the nlrow if the constraint is not convex */
         if( SCIPisConvexQuadratic(scip, conss[i]) == FALSE )
         {
            SCIP_NLROW* nlrow;
            SCIP_CALL( SCIPgetNlRowQuadratic(scip, conss[i], &nlrow) );
            assert(nlrow != NULL);

            SCIP_CALL( countNLRowVarsNonConvexity(scip, nlcounts, nlrow) );
         }
      }
   }

   /* nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   if( conshdlr != NULL )
   {
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMessage("nconss(nonlinear) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_EXPRCURV curvature;
         SCIP_CALL( SCIPgetCurvatureNonlinear(scip, conss[i], TRUE, &curvature) );

         /* only check the nlrow if the constraint is not convex */
         if(  curvature != SCIP_EXPRCURV_CONVEX )
         {
            SCIP_NLROW* nlrow;
            SCIP_CALL( SCIPgetNlRowNonlinear(scip, conss[i], &nlrow) );
            assert(nlrow != NULL);

            SCIP_CALL( countNLRowVarsNonConvexity(scip, nlcounts, nlrow) );
         }
      }
   }

   /* bivariate constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "bivariate");
   if( conshdlr != NULL )
   {
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMessage("nconss(bivariate) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_EXPRCURV curvature;
         SCIP_INTERVAL* varbounds;
         SCIP_EXPRTREE* exprtree;
         int j;

         exprtree = SCIPgetExprtreeBivariate(scip, conss[i]);
         if( exprtree != NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &varbounds, SCIPexprtreeGetNVars(exprtree)) );
            for( j = 0; j < SCIPexprtreeGetNVars(exprtree); ++j )
            {
               SCIP_VAR* var;
               var = SCIPexprtreeGetVars(exprtree)[j];

               SCIPintervalSetBounds(&varbounds[j],
                  -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))),    /*lint !e666*/
                  +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))) );  /*lint !e666*/
            }

            SCIP_CALL( SCIPexprtreeCheckCurvature(exprtree, SCIPinfinity(scip), varbounds, &curvature, NULL) );

            /* increase counter for all variables in the expression tree if the constraint is non-convex */
            if( curvature != SCIP_EXPRCURV_CONVEX )
            {
               for( j = 0; j < SCIPexprtreeGetNVars(exprtree); ++j )
               {
                  SCIP_VAR* var;
                  var = SCIPexprtreeGetVars(exprtree)[j];

                  ++nlcounts[SCIPvarGetProbindex(var)];
               }
            }
         }
      }
   }

   /* abspower constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "abspower");
   if( conshdlr != NULL )
   {
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMessage("nconss(abspower) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         /* constraint is non-convex in general */
         SCIP_NLROW* nlrow;
         SCIP_CALL( SCIPgetNlRowAbspower(scip, conss[i], &nlrow) );
         assert(nlrow != NULL);

         SCIP_CALL( countNLRowVarsNonConvexity(scip, nlcounts, nlrow) );
      }
   }

   return SCIP_OKAY;
}


/** determines whether a variable is interesting */
static
SCIP_Bool varIsInteresting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   int                   nlcount             /**< number of nonlinear constraints containing the variable
                                               *  or number of non-convex terms containing the variable
                                               * (depends on propdata->onlynonconvexvars)  */
   )
{
   assert(SCIPgetDepth(scip) == 0);

   return !SCIPvarIsBinary(var) && SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && nlcount > 0
      && !varIsFixedLocal(scip, var);
}

/** initializes interesting bounds */
static
SCIP_RETCODE initBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   int* nlcount;                             /* array that stores in how many nonlinearities each variable appears */
   int* nccount;                             /* array that stores in how many nonconvexities each variable appears */

   int bdidx;                                /* bound index inside propdata->bounds */
   int maxnlcount;                           /* maximal number of nonlinear constraints a variable appears in */
   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(SCIPisNLPConstructed(scip));

   SCIPdebugMessage("initialize bounds\n");

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* count nonlinearities */
   assert(SCIPgetNNLPVars(scip) == nvars);

   SCIP_CALL( SCIPallocBufferArray(scip, &nlcount, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nccount, nvars) );

   SCIP_CALL( SCIPgetNLPVarsNonlinearity(scip, nlcount) );
   SCIP_CALL( getNLPVarsNonConvexity(scip, nccount) );

   maxnlcount = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( maxnlcount < nlcount[i] )
         maxnlcount = nlcount[i];
   }

   /* allocate interesting bounds array */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->bounds), 2 * nvars) );

   /* get all interesting variables and their bounds */
   bdidx = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varIsInteresting(scip, vars[i], (propdata->onlynonconvexvars ? nccount[i] : nlcount[i])) )
      {
         BOUND** bdaddress;

         /* create lower bound */
         bdaddress = &(propdata->bounds[bdidx]);
         SCIP_CALL( SCIPallocMemory(scip, bdaddress) );
         propdata->bounds[bdidx]->boundtype = SCIP_BOUNDTYPE_LOWER;
         propdata->bounds[bdidx]->var = vars[i];
         propdata->bounds[bdidx]->found = FALSE;
         propdata->bounds[bdidx]->filtered = FALSE;
         propdata->bounds[bdidx]->newval = 0.0;
         propdata->bounds[bdidx]->score = getScore(scip, propdata->bounds[bdidx], nlcount[i], maxnlcount);
         propdata->bounds[bdidx]->done = FALSE;
         propdata->bounds[bdidx]->nonconvex = (nccount[i] > 0);
         propdata->bounds[bdidx]->index = bdidx;
         bdidx++;

         /* create upper bound */
         bdaddress = &(propdata->bounds[bdidx]);
         SCIP_CALL( SCIPallocMemory(scip, bdaddress) );
         propdata->bounds[bdidx]->boundtype = SCIP_BOUNDTYPE_UPPER;
         propdata->bounds[bdidx]->var = vars[i];
         propdata->bounds[bdidx]->found = FALSE;
         propdata->bounds[bdidx]->filtered = FALSE;
         propdata->bounds[bdidx]->newval = 0.0;
         propdata->bounds[bdidx]->score = getScore(scip, propdata->bounds[bdidx], nlcount[i], maxnlcount);
         propdata->bounds[bdidx]->done = FALSE;
         propdata->bounds[bdidx]->nonconvex = (nccount[i] > 0);
         propdata->bounds[bdidx]->index = bdidx;
         bdidx++;
      }
   }

   /* free memory for buffering nonlinearities */
   assert(nlcount != NULL);
   assert(nccount != NULL);
   SCIPfreeBufferArray(scip, &nlcount);
   SCIPfreeBufferArray(scip, &nccount);

   /* set number of interesting bounds */
   propdata->nbounds = bdidx;

   /*  propdata->bounds array if empty */
   if( propdata->nbounds <= 0 )
   {
      assert(propdata->nbounds == 0);
      SCIPfreeMemoryArray(scip, &(propdata->bounds));
   }

   SCIPdebugMessage("problem has %d/%d interesting bounds\n", propdata->nbounds, 2 * nvars);

   if( propdata->nbounds > 0 )
   {
      /* sort bounds according to decreasing score; although this initial order will be overruled by the distance
       * criterion later, gives a more well-defined starting situation for OBBT and might help to reduce solver
       * variability
       */
      SCIPsortDownPtr((void**) propdata->bounds, compBoundsScore, propdata->nbounds);
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolObbt)
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
   propdata->genvboundprop = propdata->creategenvbounds ? SCIPfindProp(scip, GENVBOUND_PROP_NAME) : NULL;

   SCIPdebugMessage("creating genvbounds: %s\n", propdata->genvboundprop != NULL ? "true" : "false");

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Longint itlimit;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   /* do not run in: presolving, repropagation, probing mode, if no objective propagation is allowed  */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip) || SCIPinProbing(scip) || !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* only run for nonlinear problems, i.e., if NLP is constructed */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMessage("NLP not constructed, skipping obbt\n");
      return SCIP_OKAY;
   }

   /* only run if LP all columns are in the LP, i.e., the LP is a relaxation; e.g., do not run if pricers are active
    * since pricing is not performed in probing mode
    */
   if( !SCIPallColsInLP(scip) )
   {
      SCIPdebugMessage("not all columns in LP, skipping obbt\n");
      return SCIP_OKAY;
   }

   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

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
         return SCIP_OKAY;
      }
   }
   assert(propdata->nbounds >= 0);

   /* do not run if there are no interesting bounds */
   /**@todo disable */
   if( propdata->nbounds <= 0 )
   {
      SCIPdebugMessage("there are no interesting bounds\n");
      return SCIP_OKAY;
   }

   /* only run once in a node != root */
   if( SCIPgetDepth(scip) > 0 && SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode )
   {
      return SCIP_OKAY;
   }

   SCIPdebugMessage("applying obbt for problem <%s> at depth %d\n", SCIPgetProbName(scip), SCIPgetDepth(scip));

   /* without an optimal LP solution we don't want to run; this may be because propagators with higher priority have
    * already found reductions or numerical troubles occured during LP solving
    */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPdebugMessage("aborting since no optimal LP solution is at hand\n");
      return SCIP_OKAY;
   }

   /* compute iteration limit */
   if( propdata->itlimitfactor > 0.0 )
      itlimit = (SCIP_Longint) MAX(propdata->itlimitfactor * SCIPgetNRootLPIterations(scip),
         propdata->minitlimit); /*lint !e666*/
   else
      itlimit = -1;

   /* apply obbt */
   SCIP_CALL( applyObbt(scip, propdata, itlimit, result) );
   assert(*result != SCIP_DIDNOTRUN);

   /* set current node as last node */
   propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropObbt)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* note that because we reset filtered flags to false at each call to obbt, the same bound may be filtered multiple
    * times
    */
   SCIPstatisticMessage("DIVE-LP: %" SCIP_LONGINT_FORMAT "  NFILTERED: %d NTRIVIALFILTERED: %d NSOLVED: %d "
      "FILTER-LP: %" SCIP_LONGINT_FORMAT " NGENVB(dive): %d NGENVB(aggr.): %d NGENVB(triv.) %d\n",
      propdata->nprobingiterations, propdata->nfiltered, propdata->ntrivialfiltered, propdata->nsolvedbounds,
      propdata->nfilterlpiters, propdata->ngenvboundsprobing, propdata->ngenvboundsaggrfil, propdata->ngenvboundstrivfil);

   /* free memory allocated for the bounds */
   if( propdata->nbounds > 0 )
   {
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
SCIP_DECL_PROPFREE(propFreeObbt)
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
SCIP_RETCODE SCIPincludePropObbt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create obbt propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );

   /* initialize statistic variables */
   propdata->nprobingiterations = 0;
   propdata->nfiltered = 0;
   propdata->ntrivialfiltered = 0;
   propdata->nsolvedbounds = 0;
   propdata->ngenvboundsprobing = 0;
   propdata->ngenvboundsaggrfil = 0;
   propdata->ngenvboundstrivfil = 0;
   propdata->nfilterlpiters = 0;
   propdata->lastidx = -1;
   propdata->propagatecounter = 0;
   propdata->npropagatedomreds = 0;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecObbt, propdata) );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeObbt) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolObbt) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolObbt) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropObbt) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/creategenvbounds",
         "should obbt try to provide genvbounds if possible?",
         &propdata->creategenvbounds, TRUE, DEFAULT_CREATE_GENVBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/normalize",
         "should coefficients in filtering be normalized w.r.t. the domains sizes?",
         &propdata->normalize, TRUE, DEFAULT_FILTERING_NORM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/applyfilterrounds",
         "try to filter bounds in so-called filter rounds by solving auxiliary LPs?",
         &propdata->applyfilterrounds, TRUE, DEFAULT_APPLY_FILTERROUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/applytrivialfilter",
         "try to filter bounds with the LP solution after each solve?",
         &propdata->applytrivialfilter, TRUE, DEFAULT_APPLY_TRIVIALFITLERING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/genvbdsduringfilter",
         "should we try to generate genvbounds during trivial and aggressive filtering?",
         &propdata->genvbdsduringfilter, TRUE, DEFAULT_GENVBDSDURINGFILTER, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/genvbdsduringsepa",
         "try to create genvbounds during separation process?",
         &propdata->genvbdsduringsepa, TRUE, DEFAULT_GENVBDSDURINGSEPA, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/minfilter",
         "minimal number of filtered bounds to apply another filter round",
         &propdata->nminfilter, TRUE, DEFAULT_FILTERING_MIN, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/itlimitfactor",
         "multiple of root node LP iterations used as total LP iteration limit for obbt (<= 0: no limit )",
         &propdata->itlimitfactor, FALSE, DEFAULT_ITLIMITFACTOR, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "propagating/" PROP_NAME "/minitlimit",
         "minimum LP iteration limit",
         &propdata->minitlimit, FALSE, DEFAULT_MINITLIMIT, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/dualfeastol",
         "feasibility tolerance for reduced costs used in obbt; this value is used if SCIP's dual feastol is greater",
         &propdata->dualfeastol, FALSE, DEFAULT_DUALFEASTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/conditionlimit",
         "maximum condition limit used in LP solver (-1.0: no limit)",
         &propdata->conditionlimit, FALSE, DEFAULT_CONDITIONLIMIT, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/boundstreps",
         "minimal relative improve for strengthening bounds",
         &propdata->boundstreps, FALSE, DEFAULT_BOUNDSTREPS, 0.0, 1.0, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/onlynonconvexvars",
         "only apply obbt on non-convex variables",
         &propdata->onlynonconvexvars, TRUE, DEFAULT_ONLYNONCONVEXVARS, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/tightintboundsprobing",
         "should integral bounds be tightened during the probing mode?",
         &propdata->tightintboundsprobing, TRUE, DEFAULT_TIGHTINTBOUNDSPROBING, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/tightcontboundsprobing",
         "should continuous bounds be tightened during the probing mode?",
         &propdata->tightcontboundsprobing, TRUE, DEFAULT_TIGHTCONTBOUNDSPROBING, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/orderingalgo",
        "select the type of ordering algorithm which should be used (0: no special ordering, 1: greedy, 2: greedy reverse)",
        &propdata->orderingalgo, TRUE, DEFAULT_ORDERINGALGO, 0, 2, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/separatesol",
         "should the obbt LP solution be separated?",
         &propdata->separatesol, TRUE, DEFAULT_SEPARATESOL, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/sepaminiter",
        "minimum number of iteration spend to separate an obbt LP solution",
        &propdata->sepaminiter, TRUE, DEFAULT_SEPAMINITER, 0, INT_MAX, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/sepamaxiter",
        "maximum number of iteration spend to separate an obbt LP solution",
        &propdata->sepamaxiter, TRUE, DEFAULT_SEPAMAXITER, 0, INT_MAX, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/propagatefreq",
        "trigger a propagation round after that many bound tightenings (0: no propagation)",
        &propdata->propagatefreq, TRUE, DEFAULT_PROPAGATEFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
