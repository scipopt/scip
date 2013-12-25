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

/**@file    prop_obbt.c
 * @ingroup PROPAGATORS
 * @brief   optimization-based bound tightening propagator
 * @author  Stefan Weltge
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
#define DEFAULT_DUALFEASTOL              1e-9      /**< feasibility tolerance for reduced costs used in obbt; this value
                                                    *   is used if SCIP's dual feastol is greater */
#define DEFAULT_FILTERING_MIN               2      /**< minimal number of filtered bounds to apply another filter
                                                    *   round */
#define DEFAULT_ITLIMITFACTOR             5.0      /**< multiple of root node LP iterations used as total LP iteration
                                                    *   limit for obbt (<= 0: no limit ) */
#define DEFAULT_MAXLOOKAHEAD                3      /**< maximal number of bounds evaluated without success per group
                                                    *   (-1: no limit) */

#define OBBT_SCOREBASE                      5      /**< base that is used to calculate a bounds score value */
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
   unsigned int          score;              /**< score value that is used to group bounds */
   unsigned int          filtered:1;         /**< thrown out during pre-filtering step */
   unsigned int          found:1;            /**< stores whether a probably tighter value for this bound was found */
};
typedef struct Bound BOUND;

/** bound group */
struct BoundGroup
{
   int                   firstbdindex;       /**< index of the first bound of this group in propdata->bounds array */
   int                   nbounds;            /**< number of bounds in this group */
};
typedef struct BoundGroup BOUNDGROUP;

/** propagator data */
struct SCIP_PropData
{
   BOUND**               bounds;             /**< array of interesting bounds */
   BOUNDGROUP*           boundgroups;        /**< array of bound groups */
   SCIP_ROW*             cutoffrow;          /**< pointer to current objective cutoff row */
   SCIP_PROP*            genvboundprop;      /**< pointer to genvbound propagator */
   SCIP_Longint          lastnode;           /**< number of last node where obbt was performed */
   SCIP_Real             dualfeastol;        /**< feasibility tolerance for reduced costs used in obbt; this value is
                                              *   used if SCIP's dual feastol is greater */
   SCIP_Real             itlimitfactor;      /**< LP iteration limit for obbt will be this factor times total LP
                                              *   iterations in root node */
   SCIP_Bool             applyfilterrounds;  /**< apply filter rounds? */
   SCIP_Bool             creategenvbounds;   /**< should obbt try to provide genvbounds if possible? */
   SCIP_Bool             normalize;          /**< should coefficients in filtering be normalized w.r.t. the domains
                                              *   sizes? */
   int                   maxlookahead;       /**< maximal number of bounds evaluated without success per group
                                              *   (-1: no limit) */
   int                   nbounds;            /**< length of interesting bounds array */
   int                   nboundgroups;       /**< length of boundgroups array */
   int                   nminfilter;         /**< minimal number of filtered bounds to apply another filter round */
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

/** trying to filter some bounds using the existing LP solution */
static
SCIP_RETCODE filterExistingLP(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   int*                  nfiltered           /**< how many bounds were filtered this round? */
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
   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;                          /* shortcut for current bound */

      SCIP_Real solval;                      /* the variables value in the current solution */
      SCIP_Real boundval;                    /* current local bound for the variable */

      bound = propdata->bounds[i];
      if( bound->filtered )
         continue;

      boundval = bound->boundtype == SCIP_BOUNDTYPE_UPPER ?
         SCIPvarGetUbLocal(bound->var) : SCIPvarGetLbLocal(bound->var);
      solval = SCIPvarGetLPSol(bound->var);

      /* bound is tight; since this holds for all fixed variables, those are filtered here automatically */
      if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGE(scip, solval, boundval))
         || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLE(scip, solval, boundval)) )
      {
         /* mark bound as filtered */
         bound->filtered = TRUE;

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
   int*                  nfiltered           /**< how many bounds were filtered this round */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   SCIP_Bool error;
   SCIP_Bool optimal;

   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);
   assert(nfiltered != NULL);

   *nfiltered = 0;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* solve LP */
   SCIP_CALL( solveLP(scip, itlimit, &error, &optimal) );

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

         /* mark bound as filtered */
         bound->filtered = TRUE;

         /* increase number of filtered variables */
         (*nfiltered)++;

         /* get the corresponding variable's objective coeffient */
         objcoef = SCIPgetVarObjDive(scip, bound->var);

         /* change objective coefficient if it was set up for this bound */
          if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisNegative(scip, objcoef))
             || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisPositive(scip, objcoef)) )
          {
             SCIP_CALL( SCIPchgVarObjDive(scip, bound->var, 0.0) );
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
   int nleftiterations;
   int i;
   int nfiltered;
   int nvars;

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   nolditerations = SCIPgetNLPIterations(scip);
   nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIPdebugMessage("start filter rounds\n");

   /*
    * 1.) Try first to filter lower bounds of interesting variables, whose bounds are not already filtered
    */

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], 0.0) );
   }

   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->boundtype == SCIP_BOUNDTYPE_LOWER && !propdata->bounds[i]->filtered )
      {
         SCIP_CALL( SCIPchgVarObjDive(scip, propdata->bounds[i]->var,
               getFilterCoef(scip, propdata, propdata->bounds[i]->var, SCIP_BOUNDTYPE_LOWER)) );
      }
   }

   do
   {
      SCIPdebugMessage("doing a lower bounds round\n");
      SCIP_CALL( filterRound(scip, propdata, nleftiterations, &nfiltered) );
      SCIPdebugMessage("filtered %d more bounds in lower bounds round\n", nfiltered);

      /* update iterations left */
      nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   }
   while( nfiltered >= propdata->nminfilter && ( nleftiterations == -1 ||  nleftiterations > 0 ) );

   /*
    * 2.) Now try to filter the remaining upper bounds of interesting variables, whose bounds are not already filtered
    */

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], 0.0) );
   }

   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->boundtype == SCIP_BOUNDTYPE_UPPER && !propdata->bounds[i]->filtered )
      {
         SCIP_CALL( SCIPchgVarObjDive(scip, propdata->bounds[i]->var,
               getFilterCoef(scip, propdata, propdata->bounds[i]->var, SCIP_BOUNDTYPE_UPPER)) );
      }
   }

   do
   {
      SCIPdebugMessage("doing an upper bounds round\n");
      SCIP_CALL( filterRound(scip, propdata, nleftiterations, &nfiltered) );
      SCIPdebugMessage("filtered %d more bounds in upper bounds round\n", nfiltered);

      /* update iterations left */
      nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   }
   while( nfiltered >= propdata->nminfilter && ( nleftiterations == -1 ||  nleftiterations > 0 ) );

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
   assert(!SCIPinDive(scip));
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
            SCIPdebug( ntightened++ );
         }
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("tightened %d bounds\n", ntightened);
#endif

   return SCIP_OKAY;
}

/** tries to tighten a bound in diving mode  */
static
SCIP_RETCODE tightenBoundDive(
   SCIP*                 scip,               /**< SCIP data structure */
   BOUND*                bound,              /**< bound that could be tightened */
   SCIP_Real             newval,             /**< new bound value */
   SCIP_Bool*            tightened           /**< was tightening successful? */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(bound != NULL);
   assert(tightened != NULL);

   *tightened = FALSE;

   /* get old bounds */
   lb = SCIPgetVarLbDive(scip, bound->var);
   ub = SCIPgetVarUbDive(scip, bound->var);

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
         SCIP_CALL( SCIPchgVarLbDive(scip, bound->var, newval) );
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
         SCIP_CALL( SCIPchgVarUbDive(scip, bound->var, newval) );
         *tightened = TRUE;
      }
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
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint          itlimit             /**< LP iteration limit (-1: no limit) */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   int* nextindices;                         /* next relative index to try per bound group */

   SCIP_Bool boundsleft;
   SCIP_Bool iterationsleft;
   SCIP_Longint nolditerations;
   int nleftiterations;
   int i;
   int obbtround;
   int nvars;                                /* number of the problems variables */

   assert(scip != NULL);
   assert(SCIPinDive(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   SCIPdebugMessage("solve LPs...\n");

   nolditerations = SCIPgetNLPIterations(scip);
   nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* set objective coefficients to zero */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjDive(scip, vars[i], 0.0) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &nextindices, propdata->nboundgroups) );
   BMSclearMemoryArray(nextindices, propdata->nboundgroups);

   boundsleft = TRUE;
   iterationsleft = nleftiterations == -1 || nleftiterations > 0;

   /* as long as bounds are left and maxlookahead is not exceeded, iterate over all groups */
   for( obbtround = 1; boundsleft && (propdata->maxlookahead == -1 || obbtround <= propdata->maxlookahead)
           && iterationsleft; obbtround++ )
   {
      SCIPdebugMessage("obbt round %d\n", obbtround);

      boundsleft = FALSE;

      for( i = 0; i < propdata->nboundgroups && iterationsleft; i++ )
      {
         BOUNDGROUP* bdgroup;
         SCIP_Bool lastboundsuccessful;

         bdgroup = &(propdata->boundgroups[i]);
         lastboundsuccessful = TRUE;

         while( lastboundsuccessful && nextindices[i] < bdgroup->nbounds && iterationsleft )
         {
            BOUND* bound;
            SCIP_VAR* var;

            SCIP_Bool error;
            SCIP_Bool optimal;                     /* was the LP solved to optimalilty? */

            bound = propdata->bounds[bdgroup->firstbdindex + nextindices[i]];
            var = bound->var;
            nextindices[i]++;

            SCIPdebugMessage("   applying obbt on %s bound of <%s> (local bounds: [%f,%f])\n",
               bound->boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(var),
               SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

            /* ignore filtered bounds */
            if( bound->filtered )
            {
               SCIPdebugMessage("      bound already filtered\n");
               continue;
            }

            assert(!varIsFixedLocal(scip, var));

            lastboundsuccessful = FALSE;

            /* set objective coefficient to +/- 1 (note that we minimize) */
            SCIP_CALL( SCIPchgVarObjDive(scip, var, (bound->boundtype == SCIP_BOUNDTYPE_LOWER) ? 1.0 : -1.0 ) );

            /* solve LP */
            SCIP_CALL( solveLP(scip, nleftiterations, &error, &optimal) );

            /* update iterations left */
            nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
            iterationsleft = nleftiterations == -1 || nleftiterations > 0;

            /* stop this procedure if an error occured */
            if( error )
            {
               SCIPfreeBufferArray(scip, &nextindices);
               return SCIP_OKAY;
            }

            /* this variable should also be ignored in further calls to filterExistingLP() */
            bound->filtered = TRUE;

            if( optimal )
            {
               int nfiltered;

               /* try to filter further bounds using the lp solution */
               nfiltered = 0;
               SCIP_CALL( filterExistingLP(scip, propdata, &nfiltered));
               if( nfiltered > 0 )
               {
                  SCIPdebugMessage("      filtered %d more bounds using existing lp solution\n", nfiltered);
               }

               /* store this value in the bound data structure */
               bound->newval = SCIPvarGetLPSol(var);
               bound->found = TRUE;

               SCIPdebugMessage("      var <%s>, LP value: %f\n", SCIPvarGetName(var), bound->newval);

               /* in root node we may want to create a genvbound (independent of tightening success) */
               if( SCIPgetDepth(scip) == 0 && propdata->genvboundprop != NULL )
               {
                  SCIP_CALL( createGenVBound(scip, propdata, bound) );
               }

               /* try to tighten bound in dive mode */
               SCIP_CALL( tightenBoundDive(scip, bound, bound->newval, &lastboundsuccessful) );
            }

            /* set objective coefficient back to zero */
            SCIP_CALL( SCIPchgVarObjDive(scip, var, 0.0 ) );
         }

         if( nextindices[i] < bdgroup->nbounds )
            boundsleft = TRUE;
      }
   }

   SCIPfreeBufferArray(scip, &nextindices);
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
   SCIP_Longint nolditerations;
   SCIP_Real olddualfeastol;
   int nfiltered;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   *result = SCIP_DIDNOTFIND;
   nolditerations = SCIPgetNLPIterations(scip);
   olddualfeastol = SCIPdualfeastol(scip);
   nfiltered = 0;

   /* reset bound data structure flags; fixed variables are marked as filtered */
   for( i = 0; i < propdata->nbounds; i++ )
   {
      propdata->bounds[i]->filtered = varIsFixedLocal(scip, propdata->bounds[i]->var);
      propdata->bounds[i]->found = FALSE;
   }

   /* filter variables via inspecting present LP solution */
   SCIP_CALL( filterExistingLP(scip, propdata, &nfiltered) );
   SCIPdebugMessage("filtered %d bounds via inspecting present LP solution\n", nfiltered);

   /* start diving */
   SCIP_CALL( SCIPstartDive(scip) );

   /* set dual feastol */
   if( propdata->dualfeastol < olddualfeastol )
   {
      SCIP_CALL( SCIPchgDualfeastol(scip, propdata->dualfeastol) );
   }

   /* add objective cutoff */
   SCIP_CALL( addObjCutoff(scip, propdata) );

   /* apply filtering */
   if( propdata->applyfilterrounds )
   {
      SCIP_CALL( filterBounds(scip, propdata, itlimit) );
   }

   /**@todo maybe endDive, startDive, addObjCutoff to restore old LP basis information here */

   if( itlimit > 0 )
   {
      itlimit = itlimit - ( SCIPgetNLPIterations(scip) - nolditerations );
      itlimit = MAX(itlimit, 0);
   }

   /* try to find new bounds and store them in the bound data structure */
   SCIP_CALL( findNewBounds(scip, propdata, itlimit) );

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

/** comparison method for two bounds w.r.t. their scores */
static
SCIP_DECL_SORTPTRCOMP(compBounds)
{
   BOUND* bound1 = (BOUND*) elem1;
   BOUND* bound2 = (BOUND*) elem2;

   return bound1->score == bound2->score ? 0 : ( bound1->score > bound2->score ? 1 : -1 );
}

#ifdef SCIP_DEBUG
/** prints groups of variables */
static
void printGroups(
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   int i;

   assert(propdata != NULL);
   assert(propdata->nbounds > 0);

   SCIPdebugPrintf("groups={\n");

   for( i = 0; i < propdata->nboundgroups; i++ )
   {
      int j;

      SCIPdebugPrintf("  {\n");

      for( j = 0; j < propdata->boundgroups[i].nbounds; j++ )
      {
         BOUND* bound;

         bound = propdata->bounds[propdata->boundgroups[i].firstbdindex + j];

         SCIPdebugPrintf("      %s bound of <%s>, scoreval=%u\n", bound->boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper",
            SCIPvarGetName(bound->var), bound->score);
      }

      SCIPdebugPrintf("  }\n");
   }

   SCIPdebugPrintf("}\n");
}
#endif

/** creates groups for the bounds array */
static
SCIP_RETCODE createGroups(
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   unsigned int oldscoreval;                 /* old scores value */

   int i;
   int j;

   assert(propdata != NULL);
   assert(propdata->nbounds > 0);

   /* count boundgroups */
   propdata->nboundgroups = 0;
   oldscoreval = propdata->bounds[0]->score + 1;
   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->score != oldscoreval )
      {
         propdata->nboundgroups++;
         oldscoreval = propdata->bounds[i]->score;
      }
   }

   /* allocate bound groups array */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->boundgroups), propdata->nboundgroups) );

   oldscoreval = propdata->bounds[0]->score + 1;
   j = 0;
   for( i = 0; i < propdata->nbounds; i++ )
   {
      /* if a new score value appears, create a new group */
      if( propdata->bounds[i]->score != oldscoreval )
      {
         BOUNDGROUP* newgroup;

         newgroup = &(propdata->boundgroups[j]);
         newgroup->firstbdindex = i;
         newgroup->nbounds = 1;

         j++;

         oldscoreval = propdata->bounds[i]->score;
      }
      /* otherwise the group has another member */
      else
      {
         assert(j >= 1);
         propdata->boundgroups[j - 1].nbounds++;
      }
   }
   assert(j == propdata->nboundgroups);

   return SCIP_OKAY;
}

/** determines whether a variable is interesting */
static
SCIP_Bool varIsInteresting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   int                   nlcount             /**< number of nonlinear constraints containing the variable */
   )
{
   assert(SCIPgetDepth(scip) == 0);

   return !SCIPvarIsBinary(var) && !varIsFixedLocal(scip, var) && nlcount > 0;
}

/** initializes interesting bounds */
static
SCIP_RETCODE initBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   int* nlcount;                             /* array that stores in how many nonlinear constraints each variable
                                              * appears */

   int bdidx;                                /* bound index inside propdata->bounds */
   int maxnlcount;                           /* maximal number of nonlinear constraints a variable appears in */
   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);

   SCIPdebugMessage("initialize bounds\n");

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* count nonlinearities */
   SCIP_CALL( SCIPallocBufferArray(scip, &nlcount, nvars) );
   maxnlcount = 0;
   if( SCIPisNLPConstructed(scip) )
   {
      assert(SCIPgetNNLPVars(scip) == nvars);
      SCIP_CALL( SCIPgetNLPVarsNonlinearity(scip, nlcount) );

      for( i = 0; i < nvars; i++ )
         if( maxnlcount < nlcount[i] )
            maxnlcount = nlcount[i];
   }
   else
   {
      BMSclearMemoryArray(nlcount, nvars);
   }

   /* allocate interesting bounds array */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(propdata->bounds), 2 * nvars) );

   /* get all interesting variables and their bounds */
   bdidx = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varIsInteresting(scip, vars[i], nlcount[i]) )
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
         bdidx++;
      }
   }

   /* free memory for buffering nonlinearities */
   assert(nlcount != NULL);
   SCIPfreeBufferArray(scip, &nlcount);

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

   if( propdata->nbounds > 0 )
   {
      /* sort bounds according to decreasing score */
      SCIPsortDownPtr((void**) propdata->bounds, compBounds, propdata->nbounds);

      /* create groups */
      SCIP_CALL( createGroups(propdata) );
      SCIPdebug( printGroups(propdata) );
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
   propdata->boundgroups = NULL;
   propdata->nboundgroups = -1;
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
      return SCIP_OKAY;
   }

   /* only run if LP all columns are in the LP, i.e., the LP is a relaxation; e.g., do not run if pricers are active
    * since pricing is not performed in diving mode
    */
   if( !SCIPallColsInLP(scip) )
   {
      SCIPdebugMessage("not all columns in LP, skipping obbt\n");
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

   /* apply obbt */
   SCIP_CALL( applyObbt(scip, propdata, propdata->itlimitfactor > 0.0 ?
      (SCIP_Longint) (propdata->itlimitfactor * SCIPgetNRootLPIterations(scip)) : -1, result) );
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

   /* free memory allocated for the bounds */
   if( propdata->nbounds > 0 )
   {
      assert(propdata->bounds != NULL);
      assert(propdata->boundgroups != NULL);

      /* free bound groups */
      SCIPfreeMemoryArray(scip, &(propdata->boundgroups));

      /* free bounds */
      for( i = propdata->nbounds - 1; i >= 0; i-- )
      {
         SCIPfreeMemory(scip, &(propdata->bounds[i]));
      }
      SCIPfreeMemoryArray(scip, &(propdata->bounds));
   }

   propdata->nbounds = -1;
   propdata->nboundgroups = -1;

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

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecObbt, propdata) );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeObbt) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolObbt) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolObbt) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropObbt) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/creategenvbounds",
         "should obbt try to provide genvbounds if possible?",
         &propdata->creategenvbounds, TRUE, DEFAULT_CREATE_GENVBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/normalize",
         "should coefficients in filtering be normalized w.r.t. the domains sizes?",
         &propdata->normalize, TRUE, DEFAULT_FILTERING_NORM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/"PROP_NAME"/applyfilterrounds",
         "try to filter bounds in so-called filter rounds by solving auxiliary LPs?",
         &propdata->applyfilterrounds, TRUE, DEFAULT_APPLY_FILTERROUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/"PROP_NAME"/minfilter",
         "minimal number of filtered bounds to apply another filter round",
         &propdata->nminfilter, TRUE, DEFAULT_FILTERING_MIN, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/"PROP_NAME"/maxlookahead",
         "maximal number of bounds evaluated without success per group (-1: no limit)",
         &propdata->maxlookahead, FALSE, DEFAULT_MAXLOOKAHEAD, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/itlimitfactor",
         "multiple of root node LP iterations used as total LP iteration limit for obbt (<= 0: no limit )",
         &propdata->itlimitfactor, FALSE, DEFAULT_ITLIMITFACTOR, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/"PROP_NAME"/dualfeastol",
         "feasibility tolerance for reduced costs used in obbt; this value is used if SCIP's dual feastol is greater",
         &propdata->dualfeastol, FALSE, DEFAULT_DUALFEASTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
