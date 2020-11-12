#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_erf.h"
#include "scip/cons_expr_exp.h"
#include "scip/cons_expr_log.h"
#include "scip/cons_expr_abs.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_entropy.h"
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_cos.h"
#include "scip/cons_expr_nlhdlr_bilinear.h"
#include "scip/cons_expr_nlhdlr_convex.h"
#include "scip/cons_expr_nlhdlr_default.h"
#include "scip/cons_expr_nlhdlr_perspective.h"
#include "scip/cons_expr_nlhdlr_quadratic.h"
#include "scip/cons_expr_nlhdlr_quotient.h"
#include "scip/cons_expr_nlhdlr_soc.h"
#include "scip/cons_expr_iterator.h"
#include "scip/cons_expr_rowprep.h"



/** ensures that a block memory array has at least a given size
 *
 *  if cursize is 0, then *array1 can be NULL
 */
#define ENSUREBLOCKMEMORYARRAYSIZE(scip, array1, cursize, minsize)      \
   do {                                                                 \
      int __newsize;                                                    \
      assert((scip)  != NULL);                                          \
      if( (cursize) >= (minsize) )                                      \
         break;                                                         \
      __newsize = SCIPcalcMemGrowSize(scip, minsize);                   \
      assert(__newsize >= (minsize));                                   \
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(array1), cursize, __newsize) ); \
      (cursize) = __newsize;                                            \
   } while( FALSE )

/*
 * Local methods
 */

/** create and include conshdlr to SCIP and set everything except for expression handlers */
static
SCIP_RETCODE includeConshdlrExprBasic(SCIP* scip);

/** copy expression and nonlinear handlers from sourceconshdlr to (target's) scip consexprhdlr */
static
SCIP_RETCODE copyConshdlrExprExprHdlr(
   SCIP*                 scip,               /**< (target) SCIP data structure */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint expression handler */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   int                i;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLRDATA* sourceconshdlrdata;

   assert(strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   assert(conshdlr != sourceconshdlr);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   sourceconshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert(sourceconshdlrdata != NULL);

   /* copy nonlinear handlers */
   for( i = 0; i < sourceconshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_NLHDLR* sourcenlhdlr;

      /* TODO for now just don't copy disabled nlhdlr, a clean way would probably be to first copy and disable then */
      sourcenlhdlr = sourceconshdlrdata->nlhdlrs[i];
      if( sourcenlhdlr->copyhdlr != NULL && sourcenlhdlr->enabled )
      {
         SCIP_CALL( sourcenlhdlr->copyhdlr(scip, conshdlr, sourceconshdlr, sourcenlhdlr) );
      }
   }

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** compares nonlinear handler by enforcement priority
 *
 * if handlers have same enforcement priority, then compare by detection priority, then by name
 */
static
int nlhdlrEnfoCmp(
   void*                 hdlr1,              /**< first handler */
   void*                 hdlr2               /**< second handler */
)
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(hdlr1 != NULL);
   assert(hdlr2 != NULL);

   h1 = (SCIP_NLHDLR*)hdlr1;
   h2 = (SCIP_NLHDLR*)hdlr2;

   if( h1->enfopriority != h2->enfopriority )
      return (int)(h1->enfopriority - h2->enfopriority);

   if( h1->detectpriority != h2->detectpriority )
      return (int)(h1->detectpriority - h2->detectpriority);

   return strcmp(h1->name, h2->name);
}
#endif

/** tries to automatically convert an expression constraint into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_Bool*            upgraded,           /**< buffer to store whether constraint was upgraded */
   int*                  nupgdconss,         /**< buffer to increase if constraint was upgraded */
   int*                  naddconss           /**< buffer to increase with number of additional constraints created during upgrade */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsModifiable(cons));
   assert(upgraded   != NULL);
   assert(nupgdconss != NULL);
   assert(naddconss  != NULL);

   *upgraded = FALSE;

   nupgdconss_ = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if there are no upgrade methods, we can stop */
   if( conshdlrdata->nexprconsupgrades == 0 )
      return SCIP_OKAY;

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   /* call the upgrading methods */
   SCIPdebugMsg(scip, "upgrading expression constraint <%s> (up to %d upgrade methods): ",
      SCIPconsGetName(cons), conshdlrdata->nexprconsupgrades);
   SCIPdebugPrintCons(scip, cons, NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nexprconsupgrades; ++i )
   {
      if( !conshdlrdata->exprconsupgrades[i]->active )
         continue;

      assert(conshdlrdata->exprconsupgrades[i]->exprconsupgd != NULL);

      SCIP_CALL( conshdlrdata->exprconsupgrades[i]->exprconsupgd(scip, cons, consdata->nvarexprs, &nupgdconss_,
         upgdconss, upgdconsssize) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->exprconsupgrades[i]->exprconsupgd(scip, cons, consdata->nvarexprs, &nupgdconss_,
            upgdconss, upgdconsssize) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         int j;

         SCIPdebugMsg(scip, " -> upgraded to %d constraints:\n", nupgdconss_);

         /* add the upgraded constraints to the problem and forget them */
         for( j = 0; j < nupgdconss_; ++j )
         {
            SCIPdebugMsgPrint(scip, "\t");
            SCIPdebugPrintCons(scip, upgdconss[j], NULL);

            SCIP_CALL( SCIPaddCons(scip, upgdconss[j]) );      /*lint !e613*/
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconss[j]) ); /*lint !e613*/
         }

         /* count the first upgrade constraint as constraint upgrade and the remaining ones as added constraints */
         *nupgdconss += 1;
         *naddconss += nupgdconss_ - 1;
         *upgraded = TRUE;

         /* delete upgraded constraint */
         SCIPdebugMsg(scip, "delete constraint <%s> after upgrade\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );

         break;
      }
   }

   SCIPfreeBufferArray(scip, &upgdconss);

   return SCIP_OKAY;
}

/** interval evaluation of variables as used in bound tightening
 *
 * Returns slightly relaxed local variable bounds of a variable as interval.
 * Does not relax beyond integer values, thus does not relax bounds on integer variables at all.
 */
static
SCIP_DECL_EXPR_INTEVALVAR(intEvalVarBoundTightening)
{
   SCIP_INTERVAL interval;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)intevalvardata;
   assert(conshdlrdata != NULL);

   if( conshdlrdata->globalbounds )
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }
   else
   {
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
   assert(lb <= ub);  /* can SCIP ensure by now that variable bounds are not contradicting? */

   /* implicit integer variables may have non-integer bounds, apparently (run space25a) */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
   {
      lb = EPSROUND(lb, 0.0); /*lint !e835*/
      ub = EPSROUND(ub, 0.0); /*lint !e835*/
   }

   /* integer variables should always have integral bounds in SCIP */
   assert(EPSFRAC(lb, 0.0) == 0.0 || !SCIPvarIsIntegral(var));  /*lint !e835*/
   assert(EPSFRAC(ub, 0.0) == 0.0 || !SCIPvarIsIntegral(var));  /*lint !e835*/

   switch( conshdlrdata->varboundrelax )
   {
      case 'n' : /* no relaxation */
         break;

      case 'a' : /* relax by absolute value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         if( !SCIPisInfinity(scip, -lb) )
         {
            /* reduce lb by epsilon, or to the next integer value, which ever is larger */
            SCIP_Real bnd = floor(lb);
            lb = MAX(bnd, lb - conshdlrdata->varboundrelaxamount);
         }

         if( !SCIPisInfinity(scip, ub) )
         {
            /* increase ub by epsilon, or to the next integer value, which ever is smaller */
            SCIP_Real bnd = ceil(ub);
            ub = MIN(bnd, ub + conshdlrdata->varboundrelaxamount);
         }

         break;
      }

      case 'b' : /* relax always by absolute value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         if( !SCIPisInfinity(scip, -lb) )
            lb -= conshdlrdata->varboundrelaxamount;

         if( !SCIPisInfinity(scip, ub) )
            ub += conshdlrdata->varboundrelaxamount;

         break;
      }

      case 'r' : /* relax by relative value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         /* relax bounds by epsilon*max(1,|bnd|), instead of just epsilon as in case 'a', thus we trust the first log(epsilon) digits
          * however, when domains get small, relaxing can excessively weaken bound tightening, thus do only fraction of |ub-lb| if that is smaller
          * further, do not relax beyond next integer value
          */
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_Real bnd = floor(lb);
            lb = MAX(bnd, lb - MIN(conshdlrdata->varboundrelaxamount * MAX(1.0, REALABS(lb)), 0.001 * REALABS(ub-lb)));  /*lint !e666*/
         }

         if( !SCIPisInfinity(scip, ub) )
         {
            SCIP_Real bnd = ceil(ub);
            ub = MIN(bnd, ub + MIN(conshdlrdata->varboundrelaxamount * MAX(1.0, REALABS(ub)), 0.001 * REALABS(ub-lb)));  /*lint !e666*/
         }

         break;
      }

      default :
      {
         SCIPerrorMessage("Unsupported value '%c' for varboundrelax option.\n");
         SCIPABORT();
         break;
      }
   }

   /* convert SCIPinfinity() to SCIP_INTERVAL_INFINITY */
   lb = -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb);
   ub =  infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, ub);
   assert(lb <= ub);

   SCIPintervalSetBounds(&interval, lb, ub);

   return interval;
}

/** interval evaluation of variables as used in redundancy check
 *
 * Returns local variable bounds of a variable, relaxed by feastol, as interval.
 */
static
SCIP_DECL_EXPR_INTEVALVAR(intEvalVarRedundancyCheck)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL interval;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)intevalvardata;
   assert(conshdlrdata != NULL);

   if( conshdlrdata->globalbounds )
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }
   else
   {
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
   assert(lb <= ub);  /* can SCIP ensure by now that variable bounds are not contradicting? */

   /* relax variable bounds, if there are bounds and variable is not fixed
    * (actually some assert complains if trying SCIPisRelEQ if both bounds are at different infinity)
    */
   if( !(SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub)) && !SCIPisRelEQ(scip, lb, ub) )
   {
      if( !SCIPisInfinity(scip, -lb) )
         lb -= SCIPfeastol(scip);

      if( !SCIPisInfinity(scip, ub) )
         ub += SCIPfeastol(scip);
   }

   /* convert SCIPinfinity() to SCIP_INTERVAL_INFINITY */
   lb = -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb);
   ub =  infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  ub);
   assert(lb <= ub);

   SCIPintervalSetBounds(&interval, lb, ub);

   return interval;
}

/** returns whether intersecting oldinterval with newinterval would provide a properly smaller interval
 *
 * If subsetsufficient is TRUE, then the intersection being smaller than oldinterval is sufficient.
 * If subsetsufficient is FALSE, then we require
 *  - a change from an unbounded interval to a bounded one, or
 *  - or a change from an unfixed (width > epsilon) to a fixed interval, or
 *  - a minimal tightening of one of the interval bounds as defined by SCIPis{Lb,Ub}Better.
 */
static
SCIP_Bool isIntervalBetter(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_Bool               subsetsufficient, /**< whether the intersection being a proper subset of oldinterval is sufficient */
   SCIP_INTERVAL           newinterval,      /**< new interval */
   SCIP_INTERVAL           oldinterval       /**< old interval */
   )
{
   assert(scip != NULL);
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newinterval));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, oldinterval));

   if( subsetsufficient )
      /* oldinterval \cap newinterval < oldinterval iff not oldinterval is subset of newinterval */
      return !SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, oldinterval, newinterval);

   /* check whether lower bound of interval becomes finite */
   if( oldinterval.inf <= -SCIP_INTERVAL_INFINITY && newinterval.inf > -SCIP_INTERVAL_INFINITY )
      return TRUE;

   /* check whether upper bound of interval becomes finite */
   if( oldinterval.sup >=  SCIP_INTERVAL_INFINITY && newinterval.sup >  SCIP_INTERVAL_INFINITY )
      return TRUE;

   /* check whether intersection will have width <= epsilon, if oldinterval doesn't have yet */
   if( !SCIPisEQ(scip, oldinterval.inf, oldinterval.sup) && SCIPisEQ(scip, MAX(oldinterval.inf, newinterval.inf), MIN(oldinterval.sup, newinterval.sup)) ) /*lint !e666*/
      return TRUE;

   /* check whether lower bound on interval will be better by SCIP's quality measures for boundchanges */
   if( SCIPisLbBetter(scip, newinterval.inf, oldinterval.inf, oldinterval.sup) )
      return TRUE;

   /* check whether upper bound on interval will be better by SCIP's quality measures for boundchanges */
   if( SCIPisUbBetter(scip, newinterval.sup, oldinterval.inf, oldinterval.sup) )
      return TRUE;

   return FALSE;
}


/** tightens the bounds of the auxiliary variable associated with an expression (or original variable if being a
 * variable-expression) according to given bounds
 *
 *  The given bounds may very well be the exprs activity (when called from forwardPropExpr), but can also be some
 *  tighter bounds (when called from SCIPtightenExprIntervalNonlinear).
 *
 *  Nothing will happen if SCIP is not in presolve or solve.
 */
static
SCIP_RETCODE tightenAuxVarBounds(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_EXPR*     expr,             /**< expression whose auxvar is to be tightened */
   SCIP_INTERVAL           bounds,           /**< bounds to be used for tightening (must not be empty) */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a cutoff was detected */
   int*                    ntightenings      /**< buffer to add the total number of tightenings, or NULL */
   )
{
   SCIP_VAR* var;
   SCIP_Bool tightenedlb;
   SCIP_Bool tightenedub;
   SCIP_Bool force;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);

   /* the given bounds must not be empty (we could cope, but we shouldn't be called in this situation) */
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bounds));

   *cutoff = FALSE;

   /* do not tighten variable in problem stage (important for unittests)
    * TODO put some kind of #ifdef UNITTEST around this once the unittest are modified to include the .c file (again)?
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   var = SCIPgetExprAuxVarNonlinear(expr);
   if( var == NULL )
      return SCIP_OKAY;

   /* force tightening if conshdlrdata says so or it would mean fixing the variable */
   force = SCIPconshdlrGetData(conshdlr)->forceboundtightening || SCIPisEQ(scip, bounds.inf, bounds.sup);

   /* try to tighten lower bound of (auxiliary) variable */
   SCIP_CALL( SCIPtightenVarLb(scip, var, bounds.inf, force, cutoff, &tightenedlb) );
   if( tightenedlb )
   {
      if( ntightenings != NULL )
         ++*ntightenings;
      SCIPdebugMsg(scip, "tightened lb on auxvar <%s> to %.15g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var));
   }
   if( *cutoff )
   {
      SCIPdebugMsg(scip, "cutoff when tightening lb on auxvar <%s> to %.15g\n", SCIPvarGetName(var), bounds.inf);
      return SCIP_OKAY;
   }

   /* try to tighten upper bound of (auxiliary) variable */
   SCIP_CALL( SCIPtightenVarUb(scip, var, bounds.sup, force, cutoff, &tightenedub) );
   if( tightenedub )
   {
      if( ntightenings != NULL )
         ++*ntightenings;
      SCIPdebugMsg(scip, "tightened ub on auxvar <%s> to %.15g\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var));
   }
   if( *cutoff )
   {
      SCIPdebugMsg(scip, "cutoff when tightening ub on auxvar <%s> to %.15g\n", SCIPvarGetName(var), bounds.sup);
      return SCIP_OKAY;
   }

   /* TODO expr->activity should have been reevaluated now due to boundchange-events, but it used to relax bounds
    * that seems unnecessary and we could easily undo this here, e.g.,
    * if( tightenedlb ) expr->activity.inf = bounds.inf
    */

   return SCIP_OKAY;
}


/** propagate bounds of the expressions in a given expression tree (that is, updates activity intervals)
 *  and tries to tighten the bounds of the auxiliary variables accordingly
 */
static
SCIP_RETCODE forwardPropExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_EXPR*     rootexpr,         /**< expression */
   SCIP_Bool               tightenauxvars,   /**< should the bounds of auxiliary variables be tightened? */
   SCIP_Bool*              infeasible,       /**< buffer to store whether the problem is infeasible (NULL if not needed) */
   int*                    ntightenings      /**< buffer to store the number of auxiliary variable tightenings (NULL if not needed) */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(rootexpr != NULL);

   if( infeasible != NULL )
      *infeasible = FALSE;
   if( ntightenings != NULL )
      *ntightenings = 0;

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   /* if value is valid and empty, then we cannot improve, so do nothing */
   if( rootexpr->activitytag >= conshdlrdata->lastboundrelax && SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rootexpr->activity) )
   {
      SCIPdebugMsg(scip, "stored activity of root expr is empty and valid (activitytag >= lastboundrelax (%u)), skip forwardPropExpr -> cutoff\n", conshdlrdata->lastboundrelax);

      if( infeasible != NULL )
         *infeasible = TRUE;

      rootexpr->activitytag = conshdlrdata->curboundstag;

      return SCIP_OKAY;
   }

   /* if value is up-to-date, then nothing to do */
   if( rootexpr->activitytag == conshdlrdata->curboundstag )
   {
      SCIPdebugMsg(scip, "activitytag of root expr equals curboundstag (%u), skip forwardPropExpr\n", conshdlrdata->curboundstag);

      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rootexpr->activity)); /* handled in previous if() */

      return SCIP_OKAY;
   }

   /* if activity of rootexpr is not used, but expr participated in detect (nenfos >= 0), then we do nothing
    * it seems wrong to be called for such an expression (unless we are in detect at the moment), so I add a SCIPABORT()
    * during detect, we are in some in-between state where we may want to eval activity
    * on exprs that we did not notify about their activity usage
    */
   if( rootexpr->nenfos >= 0 && rootexpr->nactivityusesprop == 0 && rootexpr->nactivityusessepa == 0 && !conshdlrdata->indetect)
   {
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "root expr activity is not used but enfo initialized, skip inteval\n");
#endif
      SCIPABORT();
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexpriterCreate(&it, consexprhdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it);  )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITINGCHILD :
         {
            /* skip child if it has been evaluated already */
            SCIP_EXPR* child;

            child = SCIPexpriterGetChildExprDFS(it);
            if( conshdlrdata->curboundstag == child->activitytag )
            {
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, child->activity) )
               {
                  if( infeasible != NULL )
                     *infeasible = TRUE;
               }

               expr = SCIPexpriterSkipDFS(it);
               continue;
            }

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
         {
            /* we should not have entered this expression if its activity was already uptodate */
            assert(expr->activitytag < conshdlrdata->curboundstag);

            /* for var exprs where varevents are catched, activity is updated immediately when the varbound has been changed
             * so we can assume that the activity is uptodate for all these variables
             * UNLESS we changed the method used to evaluate activity of variable expressions
             *   or we currently use global bounds (varevents are catched for local bound changes only)
             */
            if( expr->exprhdlr == conshdlrdata->exprvarhdlr && SCIPisConsExprExprVarEventCatched(expr) &&
                  expr->activitytag >= conshdlrdata->lastvaractivitymethodchange && !conshdlrdata->globalbounds )
            {
#ifndef NDEBUG
               SCIP_INTERVAL exprhdlrinterval;

               SCIP_CALL( SCIPexprhdlrIntEvalExpr(scip, expr, &exprhdlrinterval, conshdlrdata->intevalvar, conshdlrdata) );
               assert(SCIPisRelEQ(scip, exprhdlrinterval.inf, expr->activity.inf));
               assert(SCIPisRelEQ(scip, exprhdlrinterval.sup, expr->activity.sup));
#endif
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, "skip interval evaluation of expr for var <%s> [%g,%g]\n", SCIPvarGetName(SCIPexprvarGetVar(expr)), expr->activity.inf, expr->activity.sup);
#endif
               expr->activitytag = conshdlrdata->curboundstag;

               break;
            }

            if( expr->activitytag < conshdlrdata->lastboundrelax )
            {
               /* reset activity to entire if invalid, so we can use it as starting point below */
               SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &expr->activity);
            }
            else if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
            {
               /* If already empty, then don't try to compute even better activity.
                * If cons_expr were alone, then we should have noted that we are infeasible
                * so an assert(infeasible == NULL || *infeasible) should work here.
                * However, after reporting a cutoff due to expr->activity being empty,
                * SCIP may wander to a different node and call propagation again.
                * If no bounds in an expr-constraint have been relaxed when switching nodes
                * (so expr->activitytag >= conshdlrdata->lastboundrelax), then
                * we will still have expr->activity being empty, but will have forgotten
                * that we found infeasibility here before (!2221#note_134120).
                * Therefore we just set *infeasibility=TRUE here and stop.
                */
               if( infeasible != NULL )
                  *infeasible = TRUE;
               SCIPdebugMsg(scip, "expr %p already has empty activity -> cutoff\n", (void*)expr);
               break;
            }

            /* if activity of expr is not used, but expr participated in detect (nenfos >= 0), then do nothing */
            if( expr->nenfos >= 0 && expr->nactivityusesprop == 0 && expr->nactivityusessepa == 0 && !conshdlrdata->indetect )
            {
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, "expr %p activity is not used but enfo initialized, skip inteval\n", (void*)expr);
#endif
               break;
            }

#ifdef DEBUG_PROP
            SCIPdebugMsg(scip, "interval evaluation of expr %p ", (void*)expr);
            SCIP_CALL( SCIPprintExpr(scip, consexprhdlr, expr, NULL) );
            SCIPdebugMsgPrint(scip, ", current activity = [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif

            /* run interval eval of nonlinear handlers or expression handler */
            if( expr->nenfos > 0 )
            {
               SCIP_NLHDLR* nlhdlr;
               SCIP_INTERVAL nlhdlrinterval;
               int e;

               /* for expressions with enforcement, nlhdlrs take care of interval evaluation */
               for( e = 0; e < expr->nenfos && !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity); ++e )
               {
                  /* skip nlhdlr if it does not want to participate in activity computation */
                  if( (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY) == 0 )
                     continue;

                  nlhdlr = expr->enfos[e]->nlhdlr;
                  assert(nlhdlr != NULL);

                  /* skip nlhdlr if it does not provide interval evaluation (so it may only provide reverse propagation) */
                  if( !SCIPnlhdlrHasIntEval(nlhdlr) )
                     continue;

                  /* let nlhdlr evaluate current expression */
                  nlhdlrinterval = expr->activity;
                  SCIP_CALL( SCIPcallNlhdlrIntEvalNonlinear(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata,
                     &nlhdlrinterval, conshdlrdata->intevalvar, conshdlrdata) );
#ifdef DEBUG_PROP
                  SCIPdebugMsg(scip, " nlhdlr <%s>::inteval = [%.20g, %.20g]", nlhdlr->name, nlhdlrinterval.inf, nlhdlrinterval.sup);
#endif

                  /* update expr->activity by intersecting with computed activity */
                  SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, nlhdlrinterval);
#ifdef DEBUG_PROP
                  SCIPdebugMsgPrint(scip, " -> new activity: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
               }
            }
            else
            {
               /* for node without enforcement (before or during detect), call the callback of the exprhdlr directly */
               SCIP_INTERVAL exprhdlrinterval = expr->activity;
               SCIP_CALL( SCIPexprhdlrIntEvalExpr(scip, expr, &exprhdlrinterval, conshdlrdata->intevalvar, conshdlrdata) );
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " exprhdlr <%s>::inteval = [%.20g, %.20g]", expr->exprhdlr->name, exprhdlrinterval.inf, exprhdlrinterval.sup);
#endif

               /* update expr->activity by intersecting with computed activity */
               SCIPintervalIntersectEps(&expr->activity, SCIPepsilon(scip), expr->activity, exprhdlrinterval);
#ifdef DEBUG_PROP
               SCIPdebugMsgPrint(scip, " -> new activity: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
            }

            /* if expression is integral, then we try to tighten the interval bounds a bit
             * this should undo the addition of some unnecessary safety added by use of nextafter() in interval arithmetics, e.g., when doing pow()
             * it would be ok to use ceil() and floor(), but for safety we use SCIPceil and SCIPfloor for now
             * do this only if using boundtightening-inteval and not in redundancy check (there we really want to relax all variables)
             * boundtightening-inteval does not relax integer variables, so can omit expressions without children
             * (constants should be ok, too)
             */
            if( expr->isintegral && conshdlrdata->intevalvar == intEvalVarBoundTightening && expr->nchildren > 0 )
            {
               if( expr->activity.inf > -SCIP_INTERVAL_INFINITY )
                  expr->activity.inf = SCIPceil(scip, expr->activity.inf);
               if( expr->activity.sup <  SCIP_INTERVAL_INFINITY )
                  expr->activity.sup = SCIPfloor(scip, expr->activity.sup);
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " applying integrality: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
            }

            /* mark the current node to be infeasible if either the lower/upper bound is above/below +/- SCIPinfinity()
             * TODO this is a problem if dual-presolve fixed a variable to +/- infinity
             */
            if( SCIPisInfinity(scip, expr->activity.inf) || SCIPisInfinity(scip, -expr->activity.sup) )
            {
               SCIPdebugMsg(scip, "cut off due to activity [%g,%g] beyond infinity\n", expr->activity.inf, expr->activity.sup);
               SCIPintervalSetEmpty(&expr->activity);
            }

            /* remember that activity is uptodate now */
            expr->activitytag = conshdlrdata->curboundstag;

            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity) )
            {
               if( infeasible != NULL )
                  *infeasible = TRUE;
            }
            else if( tightenauxvars && expr->auxvar != NULL )
            {
               SCIP_Bool tighteninfeasible;

               SCIP_CALL( tightenAuxVarBounds(scip, consexprhdlr, expr, expr->activity, &tighteninfeasible, ntightenings) );
               if( tighteninfeasible )
               {
                  if( infeasible != NULL )
                     *infeasible = TRUE;
                  SCIPintervalSetEmpty(&expr->activity);
               }
            }

            break;
         }

         default:
            /* you should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** propagates bounds for each sub-expression in the reversepropqueue by starting from the root expressions
 *
 *  the expression will be traversed in breadth first search by using this queue
 *
 *  @note calling this function requires feasible intervals for each sub-expression; this is guaranteed by calling
 *  forwardPropExpr() before calling this function
 *
 *  @note calling this function with *infeasible == TRUE will only empty the queue
 */
static
SCIP_RETCODE reversePropQueue(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_Bool*              infeasible,       /**< buffer to update whether an expression's bounds were propagated to an empty interval */
   int*                    ntightenings      /**< buffer to store the number of (variable) tightenings */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(infeasible != NULL);
   assert(ntightenings != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *ntightenings = 0;

   /* main loop that calls reverse propagation for expressions on the queue
    * when reverseprop finds a tightening for an expression, then that expression is added to the queue (within the reverseprop call)
    */
   while( !SCIPqueueIsEmpty(conshdlrdata->reversepropqueue) && !(*infeasible) )
   {
      SCIP_EXPR* expr;
      SCIP_INTERVAL propbounds;
      int e;

      expr = (SCIP_EXPR*) SCIPqueueRemove(conshdlrdata->reversepropqueue);
      assert(expr != NULL);
      assert(expr->inpropqueue);
      /* mark that the expression is not in the queue anymore */
      expr->inpropqueue = FALSE;

      /* since the expr was in the propagation queue, the propbounds should belong to current propagation and should not be empty
       * (propbounds being entire doesn't make much sense, so assert this for now, too, but that could be removed)
       */
      assert(expr->propboundstag == conshdlrdata->curpropboundstag);
      assert(!SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, expr->propbounds));
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->propbounds));

      /* this intersects propbounds with activity and auxvar bounds
       * I doubt this would be much helpful, since propbounds are already subset of activity and we also propagate
       * auxvar bounds separately, so disabling this for now
       */
#ifdef SCIP_DISABLED_CODE
      propbounds = SCIPgetExprBoundsNonlinear(scip, conshdlr, expr);
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, propbounds) )
      {
         *infeasible = TRUE;
         break;
      }
#else
      propbounds = expr->propbounds;
#endif

      if( expr->nenfos > 0 )
      {
         /* for nodes with enforcement, call reverse propagation callbacks of nlhdlrs */
         for( e = 0; e < expr->nenfos && !*infeasible; ++e )
         {
            SCIP_NLHDLR* nlhdlr;
            int nreds;

            /* skip nlhdlr if it does not want to participate in activity computation */
            if( (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY) == 0 )
               continue;

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* call the reverseprop of the nlhdlr */
#ifdef SCIP_DEBUG
            SCIPdebugMsg(scip, "call reverse propagation for ");
            SCIP_CALL( SCIPprintExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), expr, NULL) );
            SCIPdebugMsgPrint(scip, " in [%g,%g] using nlhdlr <%s>\n", propbounds.inf, propbounds.sup, nlhdlr->name);
#endif

            nreds = 0;
            SCIP_CALL( SCIPcallNlhdlrReversePropNonlinear(scip, conshdlr, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata,
                     propbounds, infeasible, &nreds) );
            assert(nreds >= 0);
            *ntightenings += nreds;
         }
      }
      else
      {
         /* if expr without enforcement (before detect), call reverse propagation callback of exprhdlr directly */
         int nreds = 0;

#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "call reverse propagation for ");
         SCIP_CALL( SCIPprintExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), expr, NULL) );
         SCIPdebugMsgPrint(scip, " in [%g,%g] using exprhdlr <%s>\n", expr->activity.inf, expr->activity.sup, expr->exprhdlr->name);
#endif

         /* if someone added an expr without nlhdlr into the reversepropqueue, then this must be because its enfo hasn't
          * been initialized in detectNlhdlr yet (nenfos < 0)
          */
         assert(expr->nenfos < 0);

         /* call the reverseprop of the exprhdlr */
         SCIP_CALL( SCIPexprhdlrReversePropExpr(scip, conshdlr, expr, propbounds, infeasible, &nreds) );
         assert(nreds >= 0);
         *ntightenings += nreds;
      }
   }

   /* reset inpropqueue for all remaining expr's in queue (can happen in case of early stop due to infeasibility) */
   while( !SCIPqueueIsEmpty(conshdlrdata->reversepropqueue) )
   {
      SCIP_EXPR* expr;

      expr = (SCIP_EXPR*) SCIPqueueRemove(conshdlrdata->reversepropqueue);

      /* mark that the expression is not in the queue anymore */
      expr->inpropqueue = FALSE;
   }

   return SCIP_OKAY;
}

/** calls domain propagation for a given set of constraints
 *
 *  The algorithm alternates calls of forward and reverse propagation.
 *  Forward propagation ensures that activity of expressions is uptodate.
 *  Reverse propagation tries to derive tighter variable bounds by reversing the activity computation, using the constraints
 *  [lhs,rhs] interval as starting point.
 *
 *  the propagation algorithm works as follows:
 *
 *   1.) apply forward propagation (update activities) for all constraints not marked as propagated
 *
 *   2.) if presolve or propauxvars is disabled: collect expressions for which the constraint sides provide tighter bounds
 *       if solve and propauxvars is enabled: collect expressions for which auxvars (including those in root exprs)
 *       provide tighter bounds
 *
 *   3.) apply reverse propagation to all collected expressions; don't explore
 *       sub-expressions which have not changed since the beginning of the propagation loop
 *
 *   4.) if we have found enough tightenings go to 1.) otherwise leave propagation loop
 *
 *  @note after calling forward propagation for a constraint we mark this constraint as propagated; this flag might be
 *  reset during the reverse propagation when we find a bound tightening of a variable expression contained in the
 *  constraint; resetting this flag is done in the EVENTEXEC callback of the event handler
 *
 *  @TODO should we distinguish between expressions where activity information is used for separation and those where not,
 *    e.g., try less to propagate on convex constraints?
 */
static
SCIP_RETCODE propConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds             /**< buffer to add the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff = FALSE;
   SCIP_INTERVAL conssides;
   int ntightenings;
   int roundnr;
   SCIP_EXPRITER* revpropcollectit = NULL;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(*nchgbds >= 0);

   /* no constraints to propagate */
   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->intevalvar == intEvalVarBoundTightening);
   assert(!conshdlrdata->globalbounds);

   *result = SCIP_DIDNOTFIND;
   roundnr = 0;

   /* tightenAuxVarBounds needs to know whether boundtightenings are to be forced */
   conshdlrdata->forceboundtightening = force;

   /* invalidate all propbounds (probably not needed) */
   ++conshdlrdata->curpropboundstag;

   /* create iterator that we will use if we need to look at all auxvars */
   if( conshdlrdata->propauxvars )
   {
      SCIP_CALL( SCIPexpriterCreate(&revpropcollectit, conshdlr, SCIPblkmem(scip)) );
   }

   /* main propagation loop */
   do
   {
      SCIPdebugMsg(scip, "start propagation round %d\n", roundnr);

      assert(SCIPqueueIsEmpty(conshdlrdata->reversepropqueue));

      /* apply forward propagation (update expression activities)
       * and add promising root expressions into queue for reversepropagation
       */
      for( i = 0; i < nconss; ++i )
      {
         consdata = SCIPconsGetData(conss[i]);
         assert(consdata != NULL);

         /* skip deleted, non-active, or propagation-disabled constraints */
         if( SCIPconsIsDeleted(conss[i]) || !SCIPconsIsActive(conss[i]) || !SCIPconsIsPropagationEnabled(conss[i]) )
            continue;

         /* skip already propagated constraints, i.e., constraints where no (original) variable has changed and thus
          * activity didn't change
          */
         if( consdata->ispropagated )
            continue;

         /* update activities in expression */
         SCIPdebugMsg(scip, "call forwardPropExpr() for constraint <%s> (round %d): ", SCIPconsGetName(conss[i]), roundnr);
         SCIPdebugPrintCons(scip, conss[i], NULL);

         ntightenings = 0;
         SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, TRUE, &cutoff, &ntightenings) );
         assert(cutoff || !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, consdata->expr->activity));

         if( cutoff )
         {
            SCIPdebugMsg(scip, " -> cutoff in forwardPropExpr (due to domain error or auxvar tightening) of constraint <%s>\n", SCIPconsGetName(conss[i]));
            *result = SCIP_CUTOFF;
            break;
         }

         /* TODO for a constraint that only has an auxvar for consdata->expr (e.g., convex quadratic), we could also just do the if(TRUE)-branch */
         if( !conshdlrdata->propauxvars || consdata->expr->auxvar == NULL )
         {
            /* check whether constraint sides (relaxed by epsilon) or auxvar bounds provide a tightening
             *   (if we have auxvar (not in presolve), then bounds of the auxvar are initially set to constraint sides,
             *   so taking auxvar bounds is enough)
             */
            if( consdata->expr->auxvar == NULL )
            {
               /* relax sides by SCIPepsilon() and handle infinite sides */
               SCIP_Real lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIP_INTERVAL_INFINITY : consdata->lhs - conshdlrdata->conssiderelaxamount;
               SCIP_Real rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIP_INTERVAL_INFINITY : consdata->rhs + conshdlrdata->conssiderelaxamount;
               SCIPintervalSetBounds(&conssides, lhs, rhs);
            }
            else
            {
               conssides = intEvalVarBoundTightening(scip, consdata->expr->auxvar, (void*)SCIPconshdlrGetData(conshdlr));
            }
            SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, conshdlr, consdata->expr, conssides, &cutoff, &ntightenings) );
         }
         else
         {
            /* check whether bounds of any auxvar used in constraint provides a tightening
             *   (for the root expression, bounds of auxvar are initially set to constraint sides)
             */
            SCIP_EXPR* expr;

            assert(revpropcollectit != NULL);
            SCIP_CALL( SCIPexpriterInit(revpropcollectit, consdata->expr, SCIP_EXPRITER_BFS, FALSE) );
            for( expr = SCIPexpriterGetCurrent(revpropcollectit); !SCIPexpriterIsEnd(revpropcollectit) && !cutoff; expr = SCIPexpriterGetNext(revpropcollectit) )  /*lint !e441*/
            {
               if( expr->auxvar == NULL )
                  continue;

               conssides = intEvalVarBoundTightening(scip, expr->auxvar, (void*)SCIPconshdlrGetData(conshdlr));
               SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, conshdlr, expr, conssides, &cutoff, &ntightenings) );
            }
         }

         if( cutoff )
         {
            SCIPdebugMsg(scip, " -> cutoff after intersect with conssides of constraint <%s>\n", SCIPconsGetName(conss[i]));
            *result = SCIP_CUTOFF;
            break;
         }

         assert(ntightenings >= 0);
         if( ntightenings > 0 )
         {
            *nchgbds += ntightenings;
            *result = SCIP_REDUCEDDOM;
         }

         /* mark constraint as propagated; this will be reset via the event system when we find a variable tightening */
         consdata->ispropagated = TRUE;
      }

      /* apply backward propagation (if cutoff is TRUE, then this call empties the queue) */
      SCIP_CALL( reversePropQueue(scip, conshdlr, &cutoff, &ntightenings) );
      assert(ntightenings >= 0);
      assert(SCIPqueueIsEmpty(conshdlrdata->reversepropqueue));

      if( cutoff )
      {
         SCIPdebugMsg(scip, " -> cutoff\n");
         *result = SCIP_CUTOFF;
         break;
      }

      if( ntightenings > 0 )
      {
         *nchgbds += ntightenings;
         *result = SCIP_REDUCEDDOM;
      }
   }
   while( ntightenings > 0 && ++roundnr < conshdlrdata->maxproprounds );

   if( conshdlrdata->propauxvars )
   {
      SCIPfreeExpriter(&revpropcollectit);
   }

   conshdlrdata->forceboundtightening = FALSE;

   /* invalidate propbounds in all exprs, so noone accidentally uses them outside propagation */
   ++conshdlrdata->curpropboundstag;

   return SCIP_OKAY;
}

/* calls the reverseprop callbacks of all nlhdlrs in all expressions in all constraints using activity as bounds
 *
 * This is meant to propagate any domain restricitions on functions onto variable bounds, if possible.
 *
 * Assumes that activities are still valid and curpropboundstag does not need to be increased.
 * That is, a good place to call this function is immediately after propConss or after forwardPropExpr if outsite propagation.
 */
static
SCIP_RETCODE propExprDomains(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds             /**< buffer to add the number of changed bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_Bool cutoff = FALSE;
   int ntightenings;
   int c;
   int e;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(*nchgbds >= 0);

   assert(SCIPconshdlrGetData(conshdlr)->intevalvar == intEvalVarBoundTightening);
   assert(!SCIPconshdlrGetData(conshdlr)->globalbounds);
   assert(SCIPqueueIsEmpty(SCIPconshdlrGetData(conshdlr)->reversepropqueue));

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( c = 0; c < nconss && !cutoff; ++c )
   {
      /* skip deleted, non-active, or propagation-disabled constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsActive(conss[c]) || !SCIPconsIsPropagationEnabled(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it) && !cutoff; expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         /* call reverseprop for those nlhdlr that participate in this expr's activity computation
          * this will propagate the current activity
          */
         for( e = 0; e < expr->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;
            assert(expr->enfos[e] != NULL);

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);
            if( (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY) == 0 )
               continue;

            SCIPdebugMsg(scip, "propExprDomains calling reverseprop for expression %p [%g,%g]\n", (void*)expr,
                  expr->activity.inf, expr->activity.sup);
            ntightenings = 0;
            SCIP_CALL( SCIPcallNlhdlrReversePropNonlinear(scip, conshdlr, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata,
                     expr->activity, &cutoff, &ntightenings) );

            if( cutoff )
            {
               /* stop everything if we detected infeasibility */
               SCIPdebugMsg(scip, "detect infeasibility for constraint <%s> during reverseprop()\n", SCIPconsGetName(conss[c]));
               *result = SCIP_CUTOFF;
               break;
            }

            assert(ntightenings >= 0);
            if( ntightenings > 0 )
            {
               *nchgbds += ntightenings;
               *result = SCIP_REDUCEDDOM;
            }
         }
      }
   }

   /* apply backward propagation (if cutoff is TRUE, then this call empties the queue) */
   SCIP_CALL( reversePropQueue(scip, conshdlr, &cutoff, &ntightenings) );
   assert(ntightenings >= 0);

   if( cutoff )
   {
      SCIPdebugMsg(scip, " -> cutoff\n");
      *result = SCIP_CUTOFF;
   }
   else if( ntightenings > 0 )
   {
      *nchgbds += ntightenings;
      *result = SCIP_REDUCEDDOM;
   }

   SCIPfreeExpriter(&it);

   /* invalidate propbounds in all exprs, so noone accidentally uses them outside propagation */
   ++SCIPconshdlrGetData(conshdlr)->curpropboundstag;

   return SCIP_OKAY;
}

/** checks constraints for redundancy
 *
 * Checks whether the activity of constraint functions is a subset of the constraint sides (relaxed by feastol).
 * To compute the activity, we use forwardPropExpr(), but relax variable bounds by feastol, because solutions to be checked
 * might violate variable bounds by up to feastol, too.
 * This is the main reason why the redundancy check is not done in propConss(), which relaxes variable bounds by epsilon only.
 *
 * Also removes constraints of the form lhs <= variable <= rhs.
 *
 * @todo it would be sufficient to check constraints for which we know that they are not currently violated by a valid solution
 *
 * @note This could should not run during solving, because the forwardProp takes the bounds of auxiliary variables into account.
 * For the root expression, these bounds are already set to the constraint sides, so that the activity of every expression
 * would appear as if the constraint is redundant.
 */
static
SCIP_RETCODE checkRedundancyConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool*            cutoff,             /**< pointer to store whether infeasibility has been identified */
   int*                  ndelconss,          /**< buffer to add the number of deleted constraints */
   int*                  nchgbds             /**< buffer to add the number of variable bound tightenings */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL sides;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgbds != NULL);

   /* no constraints to check */
   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* increase curboundstag and set lastvaractivitymethodchange
    * we do this here to trigger a reevaluation of all variable bounds, since we will relax variable bounds
    * for the redundancy check differently than for domain propagation
    * we also update lastboundrelax to ensure activites of all expressions are indeed reevaluated
    */
   ++conshdlrdata->curboundstag;
   assert(conshdlrdata->curboundstag > 0);
   conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
   conshdlrdata->intevalvar = intEvalVarRedundancyCheck;

   SCIPdebugMsg(scip, "checking %d constraints for redundancy\n", nconss);

   *cutoff = FALSE;
   for( i = 0; i < nconss; ++i )
   {
      if( !SCIPconsIsActive(conss[i]) || SCIPconsIsDeleted(conss[i]) )
         continue;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* handle constant expressions separately: either the problem is infeasible or the constraint is redundant */
      if( consdata->expr->exprhdlr == SCIPgetExprHdlrValue(conshdlr) )
      {
         SCIP_Real value = SCIPexprvalGetValue(consdata->expr);

         if(  (!SCIPisInfinity(scip, -consdata->lhs) && value < consdata->lhs - SCIPfeastol(scip))
            || (!SCIPisInfinity(scip, consdata->rhs) && value > consdata->rhs + SCIPfeastol(scip)) )
         {
            SCIPdebugMsg(scip, "constant constraint <%s> is infeasible: %g in [%g,%g] ", SCIPconsGetName(conss[i]), value, consdata->lhs, consdata->rhs);
            *cutoff = TRUE;

            goto TERMINATE;
         }

         SCIPdebugMsg(scip, "constant constraint <%s> is redundant: %g in [%g,%g] ", SCIPconsGetName(conss[i]), value, consdata->lhs, consdata->rhs);

         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      /* handle variable expressions separately: tighten variable bounds to constraint sides, then remove constraint (now redundant) */
      if( consdata->expr->exprhdlr == SCIPgetExprHdlrVar(conshdlr) )
      {
         SCIP_VAR* var;
         SCIP_Bool tightened;

         var = SCIPexprvarGetVar(consdata->expr);
         assert(var != NULL);

         SCIPdebugMsg(scip, "variable constraint <%s> can be made redundant: <%s>[%g,%g] in [%g,%g]\n", SCIPconsGetName(conss[i]), SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), consdata->lhs, consdata->rhs);

         /* ensure that variable bounds are within constraint sides */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, consdata->lhs, TRUE, cutoff, &tightened) );

            if( tightened )
               ++*nchgbds;

            if( *cutoff )
               goto TERMINATE;
         }

         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, consdata->rhs, TRUE, cutoff, &tightened) );

            if( tightened )
               ++*nchgbds;

            if( *cutoff )
               goto TERMINATE;
         }

         /* delete the (now) redundant constraint locally */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      /* reevaluate all bounds to remove some possible leftovers that could be in this
       * expression from a reverse propagation in a previous propagation round
       *
       * we relax variable bounds by feastol here, as solutions that are checked later can also violate
       * variable bounds by up to feastol
       * (relaxing fixed variables seems to be too much, but they would be removed by presolve soon anyway)
       */
      SCIPdebugMsg(scip, "call forwardPropExpr() for constraint <%s>: ", SCIPconsGetName(conss[i]));
      SCIPdebugPrintCons(scip, conss[i], NULL);

      SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, FALSE, cutoff, NULL) );
      assert(*cutoff || !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, consdata->expr->activity));

      /* it is unlikely that we detect infeasibility by doing forward propagation */
      if( *cutoff )
      {
         SCIPdebugMsg(scip, " -> cutoff\n");
         goto TERMINATE;
      }

      assert(consdata->expr->activitytag == conshdlrdata->curboundstag);
      activity = consdata->expr->activity;

      /* relax sides by feastol
       * we could accept every solution that violates constraints up to feastol as redundant, so this is the most permissive we can be
       */
      SCIPintervalSetBounds(&sides,
         SCIPisInfinity(scip, -consdata->lhs) ? -SCIP_INTERVAL_INFINITY : consdata->lhs - SCIPfeastol(scip),
         SCIPisInfinity(scip,  consdata->rhs) ?  SCIP_INTERVAL_INFINITY : consdata->rhs + SCIPfeastol(scip));

      if( SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, activity, sides) )
      {
         SCIPdebugMsg(scip, " -> redundant: activity [%g,%g] within sides [%g,%g]\n", activity.inf, activity.sup, consdata->lhs, consdata->rhs);

         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      SCIPdebugMsg(scip, " -> not redundant: activity [%g,%g] not within sides [%g,%g]\n", activity.inf, activity.sup, consdata->lhs, consdata->rhs);
   }

TERMINATE:
   /* make sure all activities are reevaluated again, since we relaxed bounds in a different way */
   ++conshdlrdata->curboundstag;
   conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
   conshdlrdata->intevalvar = intEvalVarBoundTightening;

   return SCIP_OKAY;
}

/** compares enfodata by enforcement priority of nonlinear handler
 *
 * if handlers have same enforcement priority, then compare by detection priority, then by name
 */
static
int enfodataCmp(
   void*                 enfo1,              /**< first enfo data */
   void*                 enfo2               /**< second enfo data */
)
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(enfo1 != NULL);
   assert(enfo2 != NULL);

   h1 = ((EXPRENFO*)enfo1)->nlhdlr;
   h2 = ((EXPRENFO*)enfo2)->nlhdlr;

   assert(h1 != NULL);
   assert(h2 != NULL);

   if( h1->enfopriority != h2->enfopriority )
      return (int)(h1->enfopriority - h2->enfopriority);

   if( h1->detectpriority != h2->detectpriority )
      return (int)(h1->detectpriority - h2->detectpriority);

   return strcmp(h1->name, h2->name);
}

/** install nlhdlrs in one expression */
static
SCIP_RETCODE detectNlhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression for which to run detection routines */
   SCIP_CONS*            cons                /**< constraint for which expr == consdata->expr, otherwise NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD enforcemethodsallowed;
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD enforcemethods;
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD enforcemethodsnew;
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD nlhdlrenforcemethods;
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD nlhdlrparticipating;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
   int enfossize;  /* allocated length of expr->enfos array */
   int h;

   assert(conshdlr != NULL);
   assert(expr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);
   assert(!conshdlrdata->indetect);

   /* there should be no enforcer yet and detection should not even have considered expr yet */
   assert(expr->nenfos < 0);
   assert(expr->enfos == NULL);

   /* check which enforcement methods are required by setting flags in enforcemethods for those that are NOT required
    * - if no auxiliary variable usage, then do not need sepabelow or sepaabove
    * - if auxiliary variable usage, but nobody positively (up) locks expr -> only need to enforce expr >= auxvar -> no need for underestimation
    * - if auxiliary variable usage, but nobody negatively (down) locks expr -> only need to enforce expr <= auxvar -> no need for overestimation
    * - if no one uses activity, then do not need activity methods
    */
   enforcemethods = SCIP_CONSNONLINEAR_EXPRENFO_NONE;
   if( expr->nauxvaruses == 0 )
      enforcemethods |= SCIP_CONSNONLINEAR_EXPRENFO_SEPABOTH;
   else
   {
      if( SCIPgetExprNLocksPosNonlinear(expr) == 0 )  /* no need for underestimation */
         enforcemethods |= SCIP_CONSNONLINEAR_EXPRENFO_SEPABELOW;
      if( SCIPgetExprNLocksNegNonlinear(expr) == 0 )  /* no need for overestimation */
         enforcemethods |= SCIP_CONSNONLINEAR_EXPRENFO_SEPAABOVE;
   }
   if( expr->nactivityusesprop == 0 && expr->nactivityusessepa == 0 )
      enforcemethods |= SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY;

   /* it doesn't make sense to have been called on detectNlhdlr, if the expr isn't used for anything */
   assert(enforcemethods != SCIP_CONSNONLINEAR_EXPRENFO_ALL);

   /* all methods that have not been flagged above are the ones that we want to be handled by nlhdlrs */
   enforcemethodsallowed = ~enforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_ALL;

   expr->nenfos = 0;
   enfossize = 2;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr->enfos, enfossize) );
   conshdlrdata->indetect = TRUE;

   SCIPdebugMsg(scip, "detecting nlhdlrs for %s expression %p (%s); requiring%s%s%s\n",
      cons != NULL ? "root" : "non-root", (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)),
      (enforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_SEPABELOW) != 0 ? "" : " sepabelow",
      (enforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_SEPAABOVE) != 0 ? "" : " sepaabove",
      (enforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY) != 0 ? "" : " activity");

   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
   {
      SCIP_NLHDLR* nlhdlr;

      nlhdlr = conshdlrdata->nlhdlrs[h];
      assert(nlhdlr != NULL);

      /* skip disabled nlhdlrs */
      if( !nlhdlr->enabled )
         continue;

      /* call detect routine of nlhdlr */
      nlhdlrexprdata = NULL;
      enforcemethodsnew = enforcemethods;
      nlhdlrparticipating = SCIP_CONSNONLINEAR_EXPRENFO_NONE;
      conshdlrdata->registerusesactivitysepabelow = FALSE;  /* SCIPregisterExprUsageNonlinear() as called by detect may set this to TRUE */
      conshdlrdata->registerusesactivitysepaabove = FALSE;  /* SCIPregisterExprUsageNonlinear() as called by detect may set this to TRUE */
      SCIP_CALL( SCIPcallNlhdlrDetectNonlinear(scip, conshdlr, nlhdlr, expr, cons, &enforcemethodsnew, &nlhdlrparticipating, &nlhdlrexprdata) );

      /* nlhdlr might have claimed more than needed: clean up sepa flags */
      nlhdlrparticipating &= enforcemethodsallowed;

      /* detection is only allowed to augment to nlhdlrenforcemethods, so previous enforcemethods must still be set */
      assert((enforcemethodsnew & enforcemethods) == enforcemethods);

      /* Because of the previous assert, nlhdlrenforcenew ^ enforcemethods are the methods enforced by this nlhdlr.
       * They are also cleaned up here to ensure that only the needed methods are claimed.
       */
      nlhdlrenforcemethods = (enforcemethodsnew ^ enforcemethods) & enforcemethodsallowed;

      /* nlhdlr needs to participate for the methods it is enforcing */
      assert((nlhdlrparticipating & nlhdlrenforcemethods) == nlhdlrenforcemethods);

      if( nlhdlrparticipating == SCIP_CONSNONLINEAR_EXPRENFO_NONE )
      {
         /* nlhdlr might not have detected anything, or all set flags might have been removed by
          * clean up; in the latter case, we may need to free nlhdlrexprdata */

         /* free nlhdlr exprdata, if there is any and there is a method to free this data */
         if( nlhdlrexprdata != NULL && nlhdlr->freeexprdata != NULL )
         {
            SCIP_CALL( (*nlhdlr->freeexprdata)(scip, nlhdlr, expr, &nlhdlrexprdata) );
            assert(nlhdlrexprdata == NULL);
         }
         /* nlhdlr cannot have added an enforcement method if it doesn't participate (actually redundant due to previous asserts) */
         assert(nlhdlrenforcemethods == SCIP_CONSNONLINEAR_EXPRENFO_NONE);

         SCIPdebugMsg(scip, "nlhdlr <%s> detect unsuccessful\n", SCIPnlhdlrGetName(nlhdlr));

         continue;
      }

      SCIPdebugMsg(scip, "nlhdlr <%s> detect successful; sepabelow: %s, sepaabove: %s, activity: %s\n",
         SCIPnlhdlrGetName(nlhdlr),
         ((nlhdlrenforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_SEPABELOW) != 0) ? "enforcing" : ((nlhdlrparticipating & SCIP_CONSNONLINEAR_EXPRENFO_SEPABELOW) != 0) ? "participating" : "no",
         ((nlhdlrenforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_SEPAABOVE) != 0) ? "enforcing" : ((nlhdlrparticipating & SCIP_CONSNONLINEAR_EXPRENFO_SEPAABOVE) != 0) ? "participating" : "no",
         ((nlhdlrenforcemethods & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY) != 0) ? "enforcing" : ((nlhdlrparticipating & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY) != 0) ? "participating" : "no");

      /* store nlhdlr and its data */
      if( expr->nenfos == enfossize )
      {
         enfossize = SCIPcalcMemGrowSize(scip, expr->nenfos+1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &expr->enfos, expr->nenfos, enfossize) );
      }
      SCIP_CALL( SCIPallocBlockMemory(scip, &expr->enfos[expr->nenfos]) );  /*lint !e866*/
      expr->enfos[expr->nenfos]->nlhdlr = nlhdlr;
      expr->enfos[expr->nenfos]->nlhdlrexprdata = nlhdlrexprdata;
      expr->enfos[expr->nenfos]->nlhdlrparticipation = nlhdlrparticipating;
      expr->enfos[expr->nenfos]->issepainit = FALSE;
      expr->enfos[expr->nenfos]->sepabelowusesactivity = conshdlrdata->registerusesactivitysepabelow;
      expr->enfos[expr->nenfos]->sepaaboveusesactivity = conshdlrdata->registerusesactivitysepaabove;
      expr->nenfos++;

      /* update enforcement flags */
      enforcemethods = enforcemethodsnew;
   }

   conshdlrdata->indetect = FALSE;

   /* stop if an enforcement method is missing but we are already in solving stage
    * (as long as the expression provides its callbacks, the default nlhdlr should have provided all enforcement methods)
    */
   if( enforcemethods != SCIP_CONSNONLINEAR_EXPRENFO_ALL && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIPerrorMessage("no nonlinear handler provided some of the required enforcement methods\n");
      return SCIP_ERROR;
   }

   assert(expr->nenfos > 0);

   /* sort nonlinear handlers by enforcement priority, in decreasing order */
   if( expr->nenfos > 1 )
      SCIPsortDownPtr((void**)expr->enfos, enfodataCmp, expr->nenfos);

   /* resize expr->enfos array to be nenfos long */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &expr->enfos, enfossize, expr->nenfos) );

   return SCIP_OKAY;
}

/** detect nlhdlrs that can handle the expressions */
static
SCIP_RETCODE detectNlhdlrs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   SCIP_EXPRITER* it;
   int i;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE || SCIPgetStage(scip) == SCIP_STAGE_SOLVING);  /* should only be called in presolve or initsolve or consactive */

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, TRUE) );

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetDepth(scip) != 0 )
   {
      /* ensure that activities are recomputed w.r.t. the global variable bounds if CONSACTIVE is called in a local node;
       * for example, this happens if globally valid expression constraints are added during the tree search
       */
      SCIPincrementCurBoundsTagNonlinear(conshdlr, TRUE);
      conshdlrdata->globalbounds = TRUE;
      conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   }

   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL && conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      /* if a constraint is separated, we currently need it to be initial, too
       * this is because INITLP will create the auxiliary variables that are used for any separation
       * TODO we may relax this with a little more programming effort when required, see also TODO in INITLP
       */
      assert((!SCIPconsIsSeparated(conss[i]) && !SCIPconsIsEnforced(conss[i])) || SCIPconsIsInitial(conss[i]));

      /* because of common sub-expressions it might happen that we already detected a nonlinear handler and added it to the expr
       * then we would normally skip to run DETECT again
       * HOWEVER: most likely we have been running DETECT with cons == NULL, which may interest less nlhdlrs
       * thus, if expr is the root expression, we rerun DETECT
       */
      if( consdata->expr->nenfos > 0 )
      {
         SCIP_CALL( freeEnfoData(scip, consdata->expr, FALSE) );
         assert(consdata->expr->nenfos < 0);
      }

      /* if constraint will be enforced, and we are in solve, then ensure auxiliary variable for root expression
       *   this way we can treat the root expression like any other expression when enforcing via separation
       * if constraint will be propagated, then register activity usage of root expression
       * this can trigger a call to forwardPropExpr, for which we better have the indetect flag set
       */
      conshdlrdata->indetect = TRUE;
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, conshdlr, consdata->expr,
         SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && (SCIPconsIsSeparated(conss[i]) || SCIPconsIsEnforced(conss[i])),
         SCIPconsIsPropagated(conss[i]),
         FALSE, FALSE) );
      conshdlrdata->indetect = FALSE;

      /* compute integrality information for all subexpressions */
      SCIP_CALL( SCIPcomputeExprIntegrality(scip, conshdlr, consdata->expr) );

      /* run detectNlhdlr on all expr where required */
      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )  /*lint !e441*/
      {
         /* skip exprs that we already looked at */
         if( expr->nenfos >= 0 )
            continue;

         /* if there is auxvar usage, then someone requires that
          *   auxvar == expr (or auxvar >= expr or auxvar <= expr) or we are at the root expression (expr==consdata->expr)
          *   thus, we need to find nlhdlrs that separate or estimate
          * if there is activity usage, then there is someone requiring that
          *   activity of this expression is updated; this someone would also benefit from better bounds on the activity
          *   of this expression thus, we need to find nlhdlrs that do interval-evaluation
          */
         if( expr->nauxvaruses > 0 || expr->nactivityusesprop > 0 || expr->nactivityusessepa > 0 )
         {
            SCIP_CALL( detectNlhdlr(scip, conshdlr, expr, expr == consdata->expr ? conss[i] : NULL) );

            assert(expr->nenfos >= 0);
         }
         else
         {
            /* remember that we looked at this expression during detectNlhdlrs
             * even though we have not actually run detectNlhdlr, because no nlhdlr showed interest in this expr,
             * in some situations (forwardPropExpr, to be specific) we will have to distinguish between exprs for which
             * we have not initialized enforcement yet (nenfos < 0) and expressions which are just not used in enforcement (nenfos == 0)
             */
            expr->nenfos = 0;
         }
      }

      /* include this constraint into the next propagation round because the added nlhdlr may do find tighter bounds now */
      if( SCIPconsIsPropagated(conss[i]) )
         consdata->ispropagated = FALSE;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetDepth(scip) != 0 )
   {
      /* ensure that the local bounds are used again when reevaluating the expressions later;
       * this is only needed if CONSACTIVE is called in a local node (see begin of this function)
       */
      SCIPincrementCurBoundsTagNonlinear(conshdlr, FALSE);
      conshdlrdata->globalbounds = FALSE;
      conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   }
   else
   {
      /* ensure that all activities (except for var-exprs) are reevaluated since better methods may be available now */
      SCIPincrementCurBoundsTagNonlinear(conshdlr, FALSE);
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** checks for a linear variable that can be increased or decreased without harming feasibility */
static
void findUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int poslock;
   int neglock;
   int i;

   assert(conshdlr != NULL);
   assert(consdata != NULL);

   consdata->linvarincr = NULL;
   consdata->linvardecr = NULL;
   consdata->linvarincrcoef = 0.0;
   consdata->linvardecrcoef = 0.0;

   /* root expression is not a sum -> no unlocked linear variable available */
   if( SCIPexprGetHdlr(consdata->expr) != SCIPgetExprHdlrSum(conshdlr) )
      return;

   for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
   {
      SCIP_EXPR* child;

      child = SCIPexprGetChildren(consdata->expr)[i];
      assert(child != NULL);

      /* check whether the child is a variable expression */
      if( SCIPisExprVar(child) )
      {
         SCIP_VAR* var = SCIPexprvarGetVar(child);
         SCIP_Real coef = SCIPexprsumGetCoefs(consdata->expr)[i];

         if( coef > 0.0 )
         {
            poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         }
         else
         {
            poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         }

         if( SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) - neglock == 0 )
         {
            /* for a*x + f(y) \in [lhs, rhs], we can decrease x without harming other constraints */
            /* if we have already one candidate, then take the one where the loss in the objective function is less */
            if( (consdata->linvardecr == NULL) ||
               (SCIPvarGetObj(consdata->linvardecr) / consdata->linvardecrcoef > SCIPvarGetObj(var) / coef) )
            {
               consdata->linvardecr = var;
               consdata->linvardecrcoef = coef;
            }
         }

         if( SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) - poslock == 0 )
         {
            /* for a*x + f(y) \in [lhs, rhs], we can increase x without harm */
            /* if we have already one candidate, then take the one where the loss in the objective function is less */
            if( (consdata->linvarincr == NULL) ||
               (SCIPvarGetObj(consdata->linvarincr) / consdata->linvarincrcoef > SCIPvarGetObj(var) / coef) )
            {
               consdata->linvarincr = var;
               consdata->linvarincrcoef = coef;
            }
         }
      }
   }

   assert(consdata->linvarincr == NULL || consdata->linvarincrcoef != 0.0);
   assert(consdata->linvardecr == NULL || consdata->linvardecrcoef != 0.0);

#ifdef SCIP_DEBUG
   if( consdata->linvarincr != NULL )
   {
      SCIPdebugMsg(scip, "may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvarincr));
   }
   if( consdata->linvardecr != NULL )
   {
      SCIPdebugMsg(scip, "may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvardecr));
   }
#endif
}

/* forward declaration */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< expression constraint */
   );

/** initializes (pre)solving data of constraints
 *
 * This initializes data in a constraint that is used for separation, propagation, etc, and assumes that expressions will
 * not be modified.
 * In particular, this function
 * - runs the detection method of nlhldrs
 * - looks for unlocked linear variables
 * - checks curvature (if not in presolve)
 * - creates and add row to NLP (if not in presolve)
 *
 * This function can be called in presolve and solve and can be called several times with different sets of constraints,
 * e.g., it should be called in INITSOL and for constraints that are added during solve.
 */
static
SCIP_RETCODE initSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int c;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      /* check for a linear variable that can be increase or decreased without harming feasibility */
      findUnlockedLinearVar(scip, conshdlr, consdata);

      if( SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE || SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         /* call the curvature detection algorithm of the convex nonlinear handler
          * Check only for those curvature that may result in a convex inequality, i.e.,
          * whether f(x) is concave when f(x) >= lhs and/or f(x) is convex when f(x) <= rhs.
          * Also we can assume that we are nonlinear, so do not check for convex if already concave.
          */
         SCIP_Bool success = FALSE;
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPhasConsExprExprCurvature(scip, conshdlr, consdata->expr, SCIP_EXPRCURV_CONCAVE, &success, NULL) );
            if( success )
               consdata->curv = SCIP_EXPRCURV_CONCAVE;
         }
         if( !success && !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPhasConsExprExprCurvature(scip, conshdlr, consdata->expr, SCIP_EXPRCURV_CONVEX, &success, NULL) );
            if( success )
               consdata->curv = SCIP_EXPRCURV_CONVEX;
         }
         SCIPdebugMsg(scip, "root curvature of constraint %s = %d\n", SCIPconsGetName(conss[c]), consdata->curv);

         /* add nlrow representation to NLP, if NLP had been constructed */
         if( SCIPisNLPConstructed(scip) && SCIPconsIsEnabled(conss[c]) )
         {
            if( consdata->nlrow == NULL )
            {
               SCIP_CALL( createNlRow(scip, conss[c]) );
               assert(consdata->nlrow != NULL);
            }
            SCIPnlrowSetCurvature(consdata->nlrow, consdata->curv);
            SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
         }
      }
   }

   /* register non linear handlers */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** deinitializes (pre)solving data of constraints
 *
 * This removes the initialization data created in initSolve().
 *
 * This function can be called in presolve and solve.
 * TODO At the moment, it should not be called for a constraint if there are other constraints
 * that use the same expressions but still require their nlhdlr.
 * We should probably only decrement the auxvar and acitivity usage for the root expr and then
 * proceed as in detectNlhdlrs, i.e., free enfo data only where none is used.
 */
static
SCIP_RETCODE deinitSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL activity;
   SCIP_Bool rootactivityvalid;
   int c;

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   /* call deinitialization callbacks of expression and nonlinear handlers
    * free nonlinear handlers information from expressions
    * remove auxiliary variables and activityusagecounts from expressions
    */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      /* check and remember whether activity in root is valid */
      rootactivityvalid = consdata->expr->activitytag >= SCIPconshdlrGetData(conshdlr)->lastboundrelax;

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         SCIPdebugMsg(scip, "exitsepa and free nonlinear handler data for expression %p\n", (void*)expr);

         /* remove nonlinear handlers in expression and their data and auxiliary variables; reset activityusage count */
         SCIP_CALL( freeEnfoData(scip, expr, TRUE) );

         /* remove quadratic info */
         quadFree(scip, expr);

         if( rootactivityvalid )
         {
            /* ensure activity is valid if consdata->expr activity is valid
             * this is mainly to ensure that we do not leave invalid activities in parts of the expression tree where activity was not used,
             * e.g., an expr's activity was kept uptodate by a nlhdlr, but without using some childs activity
             * so this childs activity would be invalid, which can generate confusion
             */
            SCIP_CALL( SCIPevalExprActivity(scip, conshdlr, expr, &activity, TRUE) );
         }
      }

      if( consdata->nlrow != NULL )
      {
         /* remove row from NLP, if still in solving
          * if we are in exitsolve, the whole NLP will be freed anyway
          */
         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
         {
            SCIP_CALL( SCIPdelNlRow(scip, consdata->nlrow) );
         }

         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      /* forget about linear variables that can be increased or decreased without harming feasibility */
      consdata->linvardecr = NULL;
      consdata->linvarincr = NULL;

      /* forget about curvature */
      consdata->curv = SCIP_EXPRCURV_UNKNOWN;
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** stores all variable expressions into a given constraint */
static
SCIP_RETCODE storeVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_CONSDATA*          consdata          /**< constraint data */
   )
{
   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already */
   if( consdata->varexprs != NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs == NULL);
   assert(consdata->nvarexprs == 0);

   /* create array to store all variable expressions; the number of variable expressions is bounded by SCIPgetNTotalVars() */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNTotalVars(scip)) );

   SCIP_CALL( SCIPgetExprVarExprs(scip, conshdlr, consdata->expr, consdata->varexprs, &(consdata->nvarexprs)) );
   assert(SCIPgetNTotalVars(scip) >= consdata->nvarexprs);

   /* realloc array if there are less variable expression than variables */
   if( SCIPgetNTotalVars(scip) > consdata->nvarexprs )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNTotalVars(scip), consdata->nvarexprs) );
   }

   return SCIP_OKAY;
}

/** frees all variable expression stored in storeVarExprs() */
static
SCIP_RETCODE freeVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSDATA*          consdata          /**< constraint data */
   )
{
   int i;

   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already*/
   if( consdata->varexprs == NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* release variable expressions */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);
      SCIP_CALL( SCIPreleaseExpr(scip, &consdata->varexprs[i]) );
      assert(consdata->varexprs[i] == NULL);
   }

   /* free variable expressions */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->varexprs, consdata->nvarexprs);
   consdata->varexprs = NULL;
   consdata->nvarexprs = 0;

   return SCIP_OKAY;
}

/** forbid multiaggrations of variables that appear nonlinear in constraints */
static
SCIP_RETCODE forbidNonlinearVariablesMultiaggration(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_EXPRITER* it;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   if( !SCIPconshdlrGetData(conshdlr)->forbidmultaggrnlvar )
      return SCIP_OKAY;

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* if root expression is sum, then forbid multiaggregation only for variables that are not in linear terms of sum,
       *   i.e., skip children of sum that are variables
       */
      if( SCIPexprGetHdlr(consdata->expr) == SCIPgetExprHdlrSum(conshdlr) )
      {
         int i;
         SCIP_EXPR* child;
         for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
         {
            child = SCIPexprGetChildren(consdata->expr)[i];

            /* skip variable expression, as they correspond to a linear term */
            if( SCIPexprGetHdlr(child) == SCIPgetExprHdlrVar(conshdlr) )
               continue;

            for( expr = SCIPexpriterRestartDFS(it, child); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
               if( SCIPexprGetHdlr(expr) == SCIPgetExprHdlrVar(conshdlr) )
               {
                  SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, SCIPexprvarGetVar(expr)) );
               }
         }
      }
      else
      {
         for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
            if( SCIPexprGetHdlr(expr) == SCIPgetExprHdlrVar(conshdlr) )
            {
               SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, SCIPexprvarGetVar(expr)) );
            }
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}


/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPevalExpr(scip, SCIPconsGetHdlr(cons), consdata->expr, sol, soltag) );
   activity = SCIPexprGetEvalValue(consdata->expr);

   /* consider constraint as violated if it is undefined in the current point */
   if( activity == SCIP_INVALID ) /*lint !e777*/
   {
      consdata->lhsviol = SCIPinfinity(scip);
      consdata->rhsviol = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   /* compute violations */
   consdata->lhsviol = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : consdata->lhs  - activity;
   consdata->rhsviol = SCIPisInfinity(scip,  consdata->rhs) ? -SCIPinfinity(scip) : activity - consdata->rhs;

   return SCIP_OKAY;
}

/** returns absolute violation of a constraint
 *
 * @note This does not reevaluate the violation, but assumes that @ref computeViolation has been called before.
 */
static
SCIP_Real getConsAbsViolation(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return MAX3(0.0, consdata->lhsviol, consdata->rhsviol);
}

/** computes relative violation of a constraint
 *
 * @note This does not reevaluate the absolute violation, but assumes that @ref computeViolation has been called before.
 */
static
SCIP_RETCODE getConsRelViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            viol,               /**< buffer to store violation */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real scale;

   assert(cons != NULL);
   assert(viol != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *viol = getConsAbsViolation(cons);

   if( conshdlrdata->violscale == 'n' )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, *viol) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( conshdlrdata->violscale == 'a' )
   {
      scale = MAX(1.0, REALABS(SCIPexprGetEvalValue(consdata->expr)));  /*lint !e666*/

      /* consider value of side that is violated for scaling, too */
      if( consdata->lhsviol > 0.0 && REALABS(consdata->lhs) > scale )
      {
         assert(!SCIPisInfinity(scip, -consdata->lhs));
         scale = REALABS(consdata->lhs);
      }
      else if( consdata->rhsviol > 0.0 && REALABS(consdata->rhs) > scale )
      {
         assert(!SCIPisInfinity(scip,  consdata->rhs));
         scale = REALABS(consdata->rhs);
      }

      *viol /= scale;
      return SCIP_OKAY;
   }

   /* if not 'n' or 'a', then it has to be 'g' at the moment */
   assert(conshdlrdata->violscale == 'g');
   if( soltag == 0 || consdata->gradnormsoltag != soltag )
   {
      /* we need the varexprs to conveniently access the gradient */
      SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

      /* update cached value of norm of gradient */
      consdata->gradnorm = 0.0;

      /* compute gradient */
      SCIP_CALL( SCIPevalExprGradient(scip, conshdlr, consdata->expr, sol, soltag) );

      /* gradient evaluation error -> no scaling */
      if( SCIPexprGetDerivative(consdata->expr) != SCIP_INVALID ) /*lint !e777*/
      {
         int i;
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_Real deriv;

            assert(consdata->expr->difftag == consdata->varexprs[i]->difftag);
            deriv = SCIPexprGetDerivative(consdata->varexprs[i]);
            if( deriv == SCIP_INVALID ) /*lint !e777*/
            {
               /* SCIPdebugMsg(scip, "gradient evaluation error for component %d\n", i); */
               consdata->gradnorm = 0.0;
               break;
            }

            consdata->gradnorm += deriv*deriv;
         }
      }
      consdata->gradnorm = sqrt(consdata->gradnorm);
      consdata->gradnormsoltag = soltag;
   }

   *viol /= MAX(1.0, consdata->gradnorm);

   return SCIP_OKAY;
}

/** returns whether constraint is currently violated
 *
 * @note This does not reevaluate the violation, but assumes that @ref computeViolation has been called before.
 */
static
SCIP_Bool isConsViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   return getConsAbsViolation(cons) > SCIPfeastol(scip);
}

/** returns absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(x) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(x) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(x).
 * If f could not be evaluated, then return SCIPinfinity and set both violover and violunder to TRUE.
 *
 * @note This does not reevaluate the violation, but assumes that the expression has been evaluated
 */
static
SCIP_Real getExprAbsOrigViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   )
{
   SCIP_Real auxvarvalue;

   assert(expr != NULL);
   assert(expr->auxvar != NULL);

   if( expr->evalvalue == SCIP_INVALID ) /*lint !e777*/
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = TRUE;
      return SCIPinfinity(scip);
   }

   auxvarvalue = SCIPgetSolVal(scip, sol, expr->auxvar);

   if( SCIPgetExprNLocksNegNonlinear(expr) > 0 && auxvarvalue > expr->evalvalue )
   {
      if( violunder != NULL )
         *violunder = FALSE;
      if( violover != NULL )
         *violover = TRUE;
      return auxvarvalue - expr->evalvalue;
   }

   if( SCIPgetExprNLocksPosNonlinear(expr) > 0 && expr->evalvalue > auxvarvalue )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = FALSE;
      return expr->evalvalue - auxvarvalue;
   }

   if( violunder != NULL )
      *violunder = FALSE;
   if( violover != NULL )
      *violover = FALSE;
   return 0.0;
}

/** returns absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(w) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(w) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(w).
 * If f could not be evaluated, then return SCIPinfinity and set both violover and violunder to TRUE.
 *
 * @note This does not reevaluate the violation, but assumes that f(w) is passed in with auxvalue.
 */
static
SCIP_Real getExprAbsAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   SCIP_Real auxvarvalue;

   assert(expr != NULL);
   assert(expr->auxvar != NULL);

   if( auxvalue == SCIP_INVALID )  /*lint !e777*/
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = TRUE;
      return SCIPinfinity(scip);
   }

   auxvarvalue = SCIPgetSolVal(scip, sol, expr->auxvar);

   if( SCIPgetExprNLocksNegNonlinear(expr) > 0 && auxvarvalue > auxvalue )
   {
      if( violunder != NULL )
         *violunder = FALSE;
      if( violover != NULL )
         *violover = TRUE;
      return auxvarvalue - auxvalue;
   }

   if( SCIPgetExprNLocksPosNonlinear(expr) > 0 && auxvalue > auxvarvalue )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = FALSE;
      return auxvalue - auxvarvalue;
   }

   if( violunder != NULL )
      *violunder = FALSE;
   if( violover != NULL )
      *violover = FALSE;
   return 0.0;
}

/** catch variable events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   int i;

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* check if we have catched variable events already */
   if( consdata->catchedevents )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->intevalvar == intEvalVarBoundTightening);

   SCIPdebugMsg(scip, "catchVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      expr = consdata->varexprs[i];

      assert(expr != NULL);
      assert(SCIPisExprVar(expr));

      SCIP_CALL( SCIPcatchConsExprExprVarEvent(scip, expr, eventhdlr, cons) );

      /* from now on, activity of var-expr will usually be updated in processVarEvent if variable bound is changing
       * since we just registered this eventhdlr, we should make sure that the activity is also uptodate now
       */
      if( expr->activitytag < conshdlrdata->curboundstag )
      {
         SCIP_CALL( SCIPexprhdlrIntEvalExpr(scip, expr, &expr->activity, intEvalVarBoundTightening, conshdlrdata) );
         expr->activitytag = conshdlrdata->curboundstag;
#ifdef DEBUG_PROP
         SCIPdebugMsg(scip, "var-exprhdlr::inteval for var <%s> = [%.20g, %.20g]\n", SCIPvarGetName(SCIPexprvarGetVar(expr)), expr->activity.inf, expr->activity.sup);
#endif
      }
   }

   consdata->catchedevents = TRUE;

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to drop bound change events */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we have catched variable events already */
   if( !consdata->catchedevents )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   SCIPdebugMsg(scip, "dropVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = consdata->nvarexprs - 1; i >= 0; --i )
   {
      assert(consdata->varexprs[i] != NULL);

      SCIP_CALL( SCIPdropConsExprExprVarEvent(scip, consdata->varexprs[i], eventhdlr, cons) );
   }

   consdata->catchedevents = FALSE;

   return SCIP_OKAY;
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* expr;

   eventtype = SCIPeventGetType(event);
   assert(eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED));

   assert(eventdata != NULL);
   expr = (SCIP_EXPR*) eventdata;
   assert(SCIPisExprVar(expr));

   SCIPdebugMsg(scip, "  exec event %#x for variable <%s> (local [%g,%g], global [%g,%g])\n", eventtype,
         SCIPvarGetName(SCIPeventGetVar(event)),
         SCIPvarGetLbLocal(SCIPeventGetVar(event)), SCIPvarGetUbLocal(SCIPeventGetVar(event)),
         SCIPvarGetLbGlobal(SCIPeventGetVar(event)), SCIPvarGetUbGlobal(SCIPeventGetVar(event)));

   /* we only catch varevents for variables in constraints, so there should be constraints */
   assert(SCIPgetConsExprExprVarNConss(expr) > 0);
   conshdlr = SCIPconsGetHdlr(SCIPgetConsExprExprVarConss(expr)[0]);  /*lint !e613*/
   assert(conshdlr != NULL);

   /* notify constraints that use this variable expression (expr) to repropagate and possibly resimplify
    * - propagation can only find something new if a bound was tightened
    * - simplify can only find something new if a var is fixed (or maybe a bound is tightened)
    *   and we look at global changes (that is, we are not looking at boundchanges in probing)
    */
   if( eventtype & (SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED) )
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS** conss;
      int nconss;
      int c;

      nconss = SCIPgetConsExprExprVarNConss(expr);
      conss = SCIPgetConsExprExprVarConss(expr);
      assert(conss != NULL || nconss == 0);

      for( c = 0; c < nconss; ++c )
      {
         assert(conss[c] != NULL);  /*lint !e613*/
         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/

         /* if boundtightening, then mark constraints to be propagated again
          * TODO we could try be more selective here and only trigger a propagation if a relevant bound has changed,
          *   that is, we don't need to repropagate x + ... <= rhs if only the upper bound of x has been tightened
          *   the locks could help if they were available on a per-constraint base, but they aren't (and it may not be worth it)
          */
         if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         {
            consdata->ispropagated = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for propagate\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
         }

         /* if still in presolve (but not probing), then mark constraints to be unsimplified */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !SCIPinProbing(scip) )
         {
            consdata->issimplified = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for simplify\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
         }
      }
   }

   /* update curboundstag, lastboundrelax, and expr activity */
   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* increase tag on bounds */
      ++conshdlrdata->curboundstag;
      assert(conshdlrdata->curboundstag > 0);

      /* remember also if we relaxed bounds now */
      if( eventtype & SCIP_EVENTTYPE_BOUNDRELAXED )
         conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;

      /* update the activity of the var-expr here immediately
       * (we could call expr->activity = intevalvar(var, consdhlr) directly, but then the exprhdlr statistics are not updated)
       */
      SCIP_CALL( SCIPexprhdlrIntEvalExpr(scip, expr, &expr->activity, conshdlrdata->intevalvar, conshdlrdata) );
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "  var-exprhdlr::inteval = [%.20g, %.20g]\n", expr->exprhdlr->name, expr->activity.inf, expr->activity.sup);
#endif
      expr->activitytag = conshdlrdata->curboundstag;
   }

   return SCIP_OKAY;
}

/** propagates variable locks through expression and adds lock to variables */
static
SCIP_RETCODE propagateLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< expression */
   int                   nlockspos,          /**< number of positive locks */
   int                   nlocksneg           /**< number of negative locks */
   )
{
   SCIP_EXPRITER* it;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPRITER_USERDATA ituserdata;

   assert(expr != NULL);

   /* if no locks, then nothing to do, then do nothing */
   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);
   assert(SCIPexpriterGetCurrent(it) == expr); /* iterator should not have moved */

   /* store locks in root node */
   ituserdata.intvals[0] = nlockspos;
   ituserdata.intvals[1] = nlocksneg;
   SCIPexpriterSetCurrentUserData(it, ituserdata);

   while( !SCIPexpriterIsEnd(it) )
   {
      /* collect locks */
      ituserdata = SCIPexpriterGetCurrentUserData(it);
      nlockspos = ituserdata.intvals[0];
      nlocksneg = ituserdata.intvals[1];

      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR:
         {
            if( SCIPexprGetHdlr(expr) == SCIPgetExprHdlrVar(conshdlr) )
            {
               /* if a variable, then also add nlocksneg/nlockspos via SCIPaddVarLocks() */
               SCIP_CALL( SCIPaddVarLocks(scip, SCIPexprvarGetVar(expr), nlocksneg, nlockspos) );
            }

            /* add locks to expression */
            expr->nlockspos += nlockspos;
            expr->nlocksneg += nlocksneg;

            /* add monotonicity information if expression has been locked for the first time */
            if( expr->nlockspos == nlockspos && expr->nlocksneg == nlocksneg && expr->nchildren > 0
               && expr->exprhdlr->monotonicity != NULL )
            {
               int i;

               assert(expr->monotonicity == NULL);
               assert(expr->monotonicitysize == 0);

               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr->monotonicity, expr->nchildren) );
               expr->monotonicitysize = expr->nchildren;

               /* store the monotonicity for each child */
               for( i = 0; i < expr->nchildren; ++i )
               {
                  SCIP_CALL( (*expr->exprhdlr->monotonicity)(scip, conshdlr, expr, i, &expr->monotonicity[i]) );
               }
            }
            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
         {
            /* remove monotonicity information if expression has been unlocked */
            if( expr->nlockspos == 0 && expr->nlocksneg == 0 && expr->monotonicity != NULL )
            {
               assert(expr->monotonicitysize > 0);
               /* keep this assert for checking whether someone changed an expression without updating locks properly */
               assert(expr->monotonicitysize == expr->nchildren);

               SCIPfreeBlockMemoryArray(scip, &expr->monotonicity, expr->monotonicitysize);
               expr->monotonicitysize = 0;
            }
            break;
         }

         case SCIP_EXPRITER_VISITINGCHILD :
         {
            SCIP_MONOTONE monotonicity;

            assert(expr->monotonicity != NULL || expr->exprhdlr->monotonicity == NULL);

            /* get monotonicity of child */
            /* NOTE: the monotonicity stored in an expression might be different from the result obtained by
             * SCIPgetExprMonotonicity
             */
            monotonicity = expr->monotonicity != NULL ? expr->monotonicity[SCIPexpriterGetChildIdxDFS(it)] : SCIP_MONOTONE_UNKNOWN;

            /* compute resulting locks of the child expression */
            switch( monotonicity )
            {
               case SCIP_MONOTONE_INC:
                  ituserdata.intvals[0] = nlockspos;
                  ituserdata.intvals[1] = nlocksneg;
                  break;
               case SCIP_MONOTONE_DEC:
                  ituserdata.intvals[0] = nlocksneg;
                  ituserdata.intvals[1] = nlockspos;
                  break;
               case SCIP_MONOTONE_UNKNOWN:
                  ituserdata.intvals[0] = nlockspos + nlocksneg;
                  ituserdata.intvals[1] = nlockspos + nlocksneg;
                  break;
               case SCIP_MONOTONE_CONST:
                  ituserdata.intvals[0] = 0;
                  ituserdata.intvals[1] = 0;
                  break;
            }
            /* set locks in child expression */
            SCIPexpriterSetChildUserData(it, ituserdata);

            break;
         }

         default :
            /* you should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** main function for adding locks to expressions and variables; locks for an expression constraint are used to update
 *  locks for all sub-expressions and variables; locks of expressions depend on the monotonicity of expressions
 *  w.r.t. their children, e.g., consider the constraint x^2 <= 1 with x in [-2,-1] implies an up-lock for the root
 *  expression (pow) and a down-lock for its child x because x^2 is decreasing on [-2,-1]; since the monotonicity (and thus
 *  the locks) might also depend on variable bounds, the function remembers the computed monotonicity information ofcan
 *  each expression until all locks of an expression have been removed, which implies that updating the monotonicity
 *  information during the next locking of this expression does not break existing locks
 *
 *  @note when modifying the structure of an expression, e.g., during simplification, it is necessary to remove all
 *        locks from an expression and repropagating them after the structural changes have been applied; because of
 *        existing common sub-expressions, it might be necessary to remove the locks of all constraints to ensure
 *        that an expression is unlocked (see canonicalizeConstraints() for an example)
 */
static
SCIP_RETCODE addLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< expression constraint */
   int                   nlockspos,          /**< number of positive rounding locks */
   int                   nlocksneg           /**< number of negative rounding locks */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* no constraint sides -> nothing to lock */
   if( SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, -consdata->lhs) )
      return SCIP_OKAY;

   /* remember locks */
   consdata->nlockspos += nlockspos;
   consdata->nlocksneg += nlocksneg;

   assert(consdata->nlockspos >= 0);
   assert(consdata->nlocksneg >= 0);

   /* compute locks for lock propagation */
   if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlockspos + nlocksneg, nlockspos + nlocksneg));
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlockspos, nlocksneg));
   }
   else
   {
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlocksneg, nlockspos));
   }

   return SCIP_OKAY;
}

/** scales the sides of the constraint l <= sum_i c_i f_i(x) <= r according to the following rules:
 *
 *  let n_+ the number of positive coefficients c_i and n_- be the number of negative coefficients
 *
 *   i. scale by -1 if n_+ < n_-
 *
 *  ii. scale by -1 if n_+ = n_- & r = INF
 */
static
SCIP_RETCODE scaleConsSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_Bool*            changed             /**< buffer to store if the expression of cons changed */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPexprGetHdlr(consdata->expr) == SCIPgetExprHdlrSum(conshdlr) )
   {
      SCIP_Real* coefs;
      SCIP_Real constant;
      int nchildren;
      int counter = 0;

      coefs = SCIPexprsumGetCoefs(consdata->expr);
      constant = SCIPexprsumGetConstant(consdata->expr);
      nchildren = SCIPexprGetNChildren(consdata->expr);

      /* handle special case when constraint is l <= -f(x) <= r and f(x) not a sum: simplfy ensures f is not a sum */
      if( nchildren == 1 && constant == 0.0 && coefs[0] == -1.0 )
      {
         SCIP_EXPR* expr;
         expr = consdata->expr;

         consdata->expr = SCIPexprGetChildren(expr)[0];
         assert(SCIPexprGetHdlr(consdata->expr) != SCIPgetExprHdlrSum(conshdlr));

         SCIPcaptureExpr(consdata->expr);

         SCIPswapReals(&consdata->lhs, &consdata->rhs);
         consdata->lhs = -consdata->lhs;
         consdata->rhs = -consdata->rhs;

         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
         *changed = TRUE;
         return SCIP_OKAY;
      }

      /* compute n_+ - n_i */
      for( i = 0; i < nchildren; ++i )
         counter += coefs[i] > 0 ? 1 : -1;

      if( counter < 0 || (counter == 0 && SCIPisInfinity(scip, consdata->rhs)) )
      {
         SCIP_EXPR* expr;
         SCIP_Real* newcoefs;

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &newcoefs, nchildren) );

         for( i = 0; i < nchildren; ++i )
            newcoefs[i] = -coefs[i];

         /* create a new sum expression */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &expr, nchildren, SCIPexprGetChildren(consdata->expr), newcoefs, -constant) );

         /* replace expression in constraint data and scale sides */
         SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );
         consdata->expr = expr;
         SCIPswapReals(&consdata->lhs, &consdata->rhs);
         consdata->lhs = -consdata->lhs;
         consdata->rhs = -consdata->rhs;

         /* free memory */
         SCIPfreeBufferArray(scip, &newcoefs);

         *changed = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** helper method to decide whether a given expression is product of at least two binary variables */
static
SCIP_Bool isBinaryProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr                /**< product expression */
   )
{
   int nchildren;
   int i;

   assert(expr != NULL);

   /* check whether the expression is a product */
   if( SCIPexprGetHdlr(expr) != SCIPgetExprHdlrProduct(conshdlr) )
      return FALSE;

   nchildren = SCIPexprGetNChildren(expr);

   /* don't consider products with a coefficient != 1 and products with a single child; simplification will take care
    * of this expression later
    */
   if( nchildren <= 1 || SCIPexprprodGetCoef(expr) != 1.0 )
      return FALSE;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_EXPR* child;
      SCIP_VAR* var;
      SCIP_Real ub;
      SCIP_Real lb;

      child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      if( !SCIPisExprVar(child) )
         return FALSE;

      var = SCIPexprvarGetVar(child);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      /* check whether variable is integer and has [0,1] as variable bounds */
      if( !SCIPvarIsIntegral(var) || !SCIPisEQ(scip, lb, 0.0) || !SCIPisEQ(scip, ub, 1.0) )
         return FALSE;
   }

   return TRUE;
}

/** helper method to collect all bilinear binary product terms */
static
SCIP_RETCODE getBilinearBinaryTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   sumexpr,            /**< sum expression */
   SCIP_VAR**            xs,                 /**< array to collect first variable of each bilinear binary product */
   SCIP_VAR**            ys,                 /**< array to collect second variable of each bilinear binary product */
   int*                  childidxs,          /**< array to store the index of the child of each stored bilinear binary product */
   int*                  nterms              /**< pointer to store the total number of bilinear binary terms */
   )
{
   int i;

   assert(sumexpr != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(childidxs != NULL);
   assert(nterms != NULL);

   *nterms = 0;

   for( i = 0; i < SCIPexprGetNChildren(sumexpr); ++i )
   {
      SCIP_EXPR* child;

      child = SCIPexprGetChildren(sumexpr)[i];
      assert(child != NULL);

      if( SCIPexprGetNChildren(child) == 2 && isBinaryProduct(scip, conshdlr, child) )
      {
         SCIP_VAR* x = SCIPexprvarGetVar(SCIPexprGetChildren(child)[0]);
         SCIP_VAR* y = SCIPexprvarGetVar(SCIPexprGetChildren(child)[1]);

         assert(x != NULL);
         assert(y != NULL);

         if( x != y )
         {
            xs[*nterms] = x;
            ys[*nterms] = y;
            childidxs[*nterms] = i;
            ++(*nterms);
         }
      }
   }

   return SCIP_OKAY;
}

/** helper method to reformulate x_i * sum_j c_ij x_j */
static
SCIP_RETCODE reformulateFactorizedBinaryQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR*             facvar,             /**< variable that has been factorized */
   SCIP_VAR**            vars,               /**< variables of sum_j c_ij x_j */
   SCIP_Real*            coefs,              /**< coefficients of sum_j c_ij x_j */
   int                   nvars,              /**< total number of variables in sum_j c_ij x_j */
   SCIP_EXPR**  newexpr,            /**< pointer to store the new expression */
   int*                  naddconss           /**< pointer to update the total number of added constraints (might be NULL) */
   )
{
   SCIP_VAR* auxvar;
   SCIP_CONS* newcons;
   SCIP_Real minact = 0.0;
   SCIP_Real maxact = 0.0;
   SCIP_Bool integral = TRUE;
   char name [SCIP_MAXSTRLEN];
   int i;

   assert(facvar != NULL);
   assert(vars != NULL);
   assert(nvars > 1);
   assert(newexpr != NULL);

   /* compute minimum and maximum activity of sum_j c_ij x_j */
   /* TODO could compute minact and maxact for facvar=0 and facvar=1 separately, taking implied bounds into account, allowing for possibly tighter big-M's below */
   for( i = 0; i < nvars; ++i )
   {
      minact += MIN(coefs[i], 0.0);
      maxact += MAX(coefs[i], 0.0);
      integral = integral && SCIPisIntegral(scip, coefs[i]);
   }
   assert(minact <= maxact);

   /* create and add auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateVarBasic(scip, &auxvar, name, minact, maxact, 0.0, integral ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );

   /* create and add z - maxact x <= 0 */
   if( !SCIPisZero(scip, maxact) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_1", SCIPconsGetName(cons), SCIPvarGetName(facvar));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &newcons, name, auxvar, facvar, -maxact, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      if( naddconss != NULL )
         ++(*naddconss);
   }

   /* create and add  0 <= z - minact x */
   if( !SCIPisZero(scip, minact) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_2", SCIPconsGetName(cons), SCIPvarGetName(facvar));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &newcons, name, auxvar, facvar, -minact, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      if( naddconss != NULL )
         ++(*naddconss);
   }

   /* create and add minact <= sum_j c_j x_j - z + minact x_i */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_3", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &newcons, name, nvars, vars, coefs, minact, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCoefLinear(scip, newcons, auxvar, -1.0) );
   if( !SCIPisZero(scip, minact) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, facvar, minact) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   if( naddconss != NULL )
      ++(*naddconss);

   /* create and add sum_j c_j x_j - z + maxact x_i <= maxact */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_4", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &newcons, name, nvars, vars, coefs, -SCIPinfinity(scip), maxact) );
   SCIP_CALL( SCIPaddCoefLinear(scip, newcons, auxvar, -1.0) );
   if( !SCIPisZero(scip, maxact) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, facvar, maxact) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   if( naddconss != NULL )
      ++(*naddconss);

   /* create variable expression */
   SCIP_CALL( createExprVar(scip, conshdlr, newexpr, auxvar) );

   /* release auxvar */
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

   return SCIP_OKAY;
}

/** helper method to generate an expression for a sum of product of binary variables; note that the method captures the generated expression */
static
SCIP_RETCODE getFactorizedBinaryQuadraticExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_EXPR*   sumexpr,            /**< sum expression */
   int                   minterms,           /**< minimum number of terms in a the sum of x_i sum_j c_j x_j */
   SCIP_EXPR**  newexpr,            /**< pointer to store the expression that represents the binary quadratic */
   int*                  naddconss           /**< pointer to update the total number of added constraints (might be NULL) */
   )
{
   SCIP_EXPR** exprs = NULL;
   SCIP_VAR** tmpvars = NULL;
   SCIP_VAR** vars = NULL;
   SCIP_VAR** xs = NULL;
   SCIP_VAR** ys = NULL;
   SCIP_Real* exprcoefs = NULL;
   SCIP_Real* tmpcoefs = NULL;
   SCIP_Real* sumcoefs;
   SCIP_Bool* isused  = NULL;
   int* childidxs = NULL;
   int* count = NULL;
   int nchildren;
   int nexprs = 0;
   int nterms;
   int nvars;
   int ntotalvars;
   int i;

   assert(sumexpr != NULL);
   assert(minterms > 1);
   assert(newexpr != NULL);

   *newexpr = NULL;

   /* check whether sumexpr is indeed a sum */
   if( SCIPexprGetHdlr(sumexpr) != SCIPgetExprHdlrSum(conshdlr) )
      return SCIP_OKAY;

   nchildren = SCIPexprGetNChildren(sumexpr);
   sumcoefs = SCIPexprsumGetCoefs(sumexpr);
   nvars = SCIPgetNVars(scip);
   ntotalvars = SCIPgetNTotalVars(scip);

   /* check whether there are enough terms available */
   if( nchildren < minterms )
      return SCIP_OKAY;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &childidxs, nchildren) );

   /* collect all bilinear binary product terms */
   SCIP_CALL( getBilinearBinaryTerms(scip, conshdlr, sumexpr, xs, ys, childidxs, &nterms) );

   /* check whether there are enough terms available */
   if( nterms < minterms )
      goto TERMINATE;

   /* store how often each variable appears in a bilinear binary product */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, SCIPgetVars(scip), nvars) ); /*lint !e666*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &count, ntotalvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &isused, nchildren) );

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprcoefs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, MIN(nterms, nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoefs, MIN(nterms, nvars)) );

   for( i = 0; i < nterms; ++i )
   {
      int xidx;
      int yidx;

      assert(xs[i] != NULL);
      assert(ys[i] != NULL);

      xidx = SCIPvarGetIndex(xs[i]);
      assert(xidx < ntotalvars);
      yidx = SCIPvarGetIndex(ys[i]);
      assert(yidx < ntotalvars);

      ++count[xidx];
      ++count[yidx];

      SCIPdebugMsg(scip, "increase counter for %s to %d\n", SCIPvarGetName(xs[i]), count[xidx]);
      SCIPdebugMsg(scip, "increase counter for %s to %d\n", SCIPvarGetName(ys[i]), count[yidx]);
   }

   /* sort variables; don't change order of count array because it depends on problem indices */
   {
      int* tmpcount;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpcount, count, nvars) );
      SCIPsortDownIntPtr(tmpcount, (void**)vars, nvars);
      SCIPfreeBufferArray(scip, &tmpcount);
   }

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* facvar = vars[i];
      int ntmpvars = 0;
      int j;

      /* skip candidate if there are not enough terms left */
      if( count[SCIPvarGetIndex(vars[i])] < minterms )
         continue;

      SCIPdebugMsg(scip, "consider facvar = %s with count = %d\n", SCIPvarGetName(facvar), count[SCIPvarGetIndex(vars[i])]);

      /* collect variables for x_i * sum_j c_ij x_j */
      for( j = 0; j < nterms; ++j )
      {
         int childidx = childidxs[j];
         assert(childidx >= 0 && childidx < nchildren);

         if( !isused[childidx] && (xs[j] == facvar || ys[j] == facvar) )
         {
            SCIP_Real coef;
            int xidx;
            int yidx;

            coef = sumcoefs[childidx];
            assert(coef != 0.0);

            /* collect corresponding variable */
            tmpvars[ntmpvars] = (xs[j] == facvar) ? ys[j] : xs[j];
            tmpcoefs[ntmpvars] = coef;
            ++ntmpvars;

            /* update counters */
            xidx = SCIPvarGetIndex(xs[j]);
            assert(xidx < ntotalvars);
            yidx = SCIPvarGetIndex(ys[j]);
            assert(yidx < ntotalvars);
            --count[xidx];
            --count[yidx];
            assert(count[xidx] >= 0);
            assert(count[yidx] >= 0);

            /* mark term to be used */
            isused[childidx] = TRUE;
         }
      }
      assert(ntmpvars >= minterms);
      assert(SCIPvarGetIndex(facvar) < ntotalvars);
      assert(count[SCIPvarGetIndex(facvar)] == 0); /* facvar should not appear in any other bilinear term */

      /* create required constraints and store the generated expression */
      SCIP_CALL( reformulateFactorizedBinaryQuadratic(scip, conshdlr, cons, facvar, tmpvars, tmpcoefs, ntmpvars, &exprs[nexprs], naddconss) );
      exprcoefs[nexprs] = 1.0;
      ++nexprs;
   }

   /* factorization was only successful if at least one expression has been generated */
   if( nexprs > 0 )
   {
      int nexprsold = nexprs;

      /* add all children of the sum that have not been used */
      for( i = 0; i < nchildren; ++i )
      {
         if( !isused[i] )
         {
            exprs[nexprs] = SCIPexprGetChildren(sumexpr)[i];
            exprcoefs[nexprs] = sumcoefs[i];
            ++nexprs;
         }
      }

      /* create a new sum expression */
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, newexpr, nexprs, exprs, exprcoefs, SCIPexprsumGetConstant(sumexpr)) );

      /* release all expressions that have been generated by reformulateFactorizedBinaryQuadratic() */
      for( i = 0; i < nexprsold; ++i )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &exprs[i]) );
      }
   }

TERMINATE:
   /* free memory */
   SCIPfreeBufferArrayNull(scip, &tmpcoefs);
   SCIPfreeBufferArrayNull(scip, &tmpvars);
   SCIPfreeBufferArrayNull(scip, &exprcoefs);
   SCIPfreeBufferArrayNull(scip, &exprs);
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &isused);
   SCIPfreeBufferArrayNull(scip, &count);
   SCIPfreeBufferArray(scip, &childidxs);
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);

   return SCIP_OKAY;
}

/** helper method to create an AND constraint or varbound constraints for a given binary product expression */
static
SCIP_RETCODE getBinaryProductExprDo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   prodexpr,           /**< product expression */
   SCIP_EXPR**  newexpr,            /**< pointer to store the expression that represents the product */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   SCIP_Bool             empathy4and         /**< whether to use an AND constraint, if possible */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS* cons;
   SCIP_Real* coefs;
   SCIP_VAR* w;
   char name[SCIP_MAXSTRLEN];
   int nchildren;
   int i;

   assert(conshdlr != NULL);
   assert(prodexpr != NULL);
   assert(newexpr != NULL);

   nchildren = SCIPexprGetNChildren(prodexpr);
   assert(nchildren >= 2);

   /* memory to store the variables of the variable expressions (+1 for w) */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nchildren + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nchildren + 1) );

   /* prepare the names of the variable and the constraints */
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform");
   for( i = 0; i < nchildren; ++i )
   {
      vars[i] = SCIPexprvarGetVar(SCIPexprGetChildren(prodexpr)[i]);
      coefs[i] = 1.0;
      assert(vars[i] != NULL);
      (void) strcat(name, "_");
      (void) strcat(name, SCIPvarGetName(vars[i]));
   }

   /* create and add variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIPdebugMsg(scip, "  created auxiliary variable %s\n", name);

   /* use variable bound constraints if it is a bilinear product and there is no empathy for an AND constraint */
   if( nchildren == 2 && !empathy4and )
   {
      SCIP_VAR* x = vars[0];
      SCIP_VAR* y = vars[1];

      assert(x != NULL);
      assert(y != NULL);
      assert(x != y);

      /* create and add x - w >= 0 */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_1", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, name, x, w, -1.0, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* create and add y - w >= 0 */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_2", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, name, y, w, -1.0, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* create and add x + y - w <= 1 */
      vars[2] = w;
      coefs[2] = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_3", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 3, vars, coefs, -SCIPinfinity(scip), 1.0) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* update number of added constraints */
      if( naddconss != NULL )
         *naddconss += 3;
   }
   else
   {
      /* create, add, and release AND constraint */
      SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, name, w, nchildren, vars) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIPdebugMsg(scip, "  create AND constraint\n");

      /* update number of added constraints */
      if( naddconss != NULL )
         *naddconss += 1;
   }

   /* create variable expression */
   SCIP_CALL( createExprVar(scip, conshdlr, newexpr, w) );

   /* release created variable */
   SCIP_CALL( SCIPreleaseVar(scip, &w) );

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** helper method to generate an expression for the product of binary variables; note that the method captures the generated expression */
static
SCIP_RETCODE getBinaryProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         exprmap,            /**< map to remember generated variables for visited product expressions */
   SCIP_EXPR*   prodexpr,           /**< product expression */
   SCIP_EXPR**  newexpr,            /**< pointer to store the expression that represents the product */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to update the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nchildren;

   assert(prodexpr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;

   /* only consider products of binary variables */
   if( !isBinaryProduct(scip, conshdlr, prodexpr) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   nchildren = SCIPexprGetNChildren(prodexpr);
   assert(nchildren >= 2);

   /* check whether there is already an expression that represents the product */
   if( SCIPhashmapExists(exprmap, (void*)prodexpr) )
   {
      *newexpr = (SCIP_EXPR*) SCIPhashmapGetImage(exprmap, (void*)prodexpr);
      assert(*newexpr != NULL);

      /* capture expression */
      SCIPcaptureExpr(*newexpr);
   }
   else
   {
      SCIPdebugMsg(scip, "  product expression %p has been considered for the first time\n", (void*)prodexpr);

      if( nchildren == 2 )
      {
         SCIP_CLIQUE** xcliques;
         SCIP_VAR* x;
         SCIP_VAR* y;
         SCIP_Bool found_clique = FALSE;
         int c;

         /* get variables from the product expression */
         x = SCIPexprvarGetVar(SCIPexprGetChildren(prodexpr)[0]);
         assert(x != NULL);
         y = SCIPexprvarGetVar(SCIPexprGetChildren(prodexpr)[1]);
         assert(y != NULL);
         assert(x != y);

         /* first try to find a clique containing both variables */
         xcliques = SCIPvarGetCliques(x, TRUE);

         /* look in cliques containing x */
         for( c = 0; c < SCIPvarGetNCliques(x, TRUE); ++c )
         {
            if( SCIPcliqueHasVar(xcliques[c], y, TRUE) ) /* x + y <= 1 => x*y = 0 */
            {
               /* create zero value expression */
               SCIP_CALL( SCIPcreateExprValue(scip, conshdlr, newexpr, 0.0) );

               if( nchgcoefs != NULL )
                  *nchgcoefs += 1;

               found_clique = TRUE;
               break;
            }

            if( SCIPcliqueHasVar(xcliques[c], y, FALSE) ) /* x + (1-y) <= 1 => x*y = x */
            {
               /* create variable expression for x */
               SCIP_CALL( createExprVar(scip, conshdlr, newexpr, x) );

               if( nchgcoefs != NULL )
                  *nchgcoefs += 2;

               found_clique = TRUE;
               break;
            }
         }

         if( !found_clique )
         {
            xcliques = SCIPvarGetCliques(x, FALSE);

            /* look in cliques containing complement of x */
            for( c = 0; c < SCIPvarGetNCliques(x, FALSE); ++c )
            {
               if( SCIPcliqueHasVar(xcliques[c], y, TRUE) ) /* (1-x) + y <= 1 => x*y = y */
               {
                  /* create variable expression for y */
                  SCIP_CALL( createExprVar(scip, conshdlr, newexpr, y) );

                  if( nchgcoefs != NULL )
                     *nchgcoefs += 1;

                  found_clique = TRUE;
                  break;
               }

               if( SCIPcliqueHasVar(xcliques[c], y, FALSE) ) /* (1-x) + (1-y) <= 1 => x*y = x + y - 1 */
               {
                  /* create sum expression */
                  SCIP_EXPR* sum_children[2];
                  SCIP_Real sum_coefs[2];
                  SCIP_CALL( createExprVar(scip, conshdlr, &sum_children[0], x) );
                  SCIP_CALL( createExprVar(scip, conshdlr, &sum_children[1], y) );
                  sum_coefs[0] = 1.0;
                  sum_coefs[1] = 1.0;
                  SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, newexpr, 2, sum_children, sum_coefs, -1.0) );

                  SCIP_CALL( SCIPreleaseExpr(scip, &sum_children[0]) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &sum_children[1]) );

                  if( nchgcoefs != NULL )
                     *nchgcoefs += 3;

                  found_clique = TRUE;
                  break;
               }
            }
         }

         /* if the variables are not in a clique, do standard linearization */
         if( !found_clique )
         {
            SCIP_CALL( getBinaryProductExprDo(scip, conshdlr, prodexpr, newexpr, naddconss,
               conshdlrdata->reformbinprodsand) );
         }
      }
      else
      {
         /* linearize binary product using an AND constraint because nchildren > 2 */
         SCIP_CALL( getBinaryProductExprDo(scip, conshdlr, prodexpr, newexpr, naddconss,
            conshdlrdata->reformbinprodsand) );
      }

      /* hash variable expression */
      SCIP_CALL( SCIPhashmapInsert(exprmap, (void*)prodexpr, *newexpr) );
   }

   return SCIP_OKAY;
}

/** helper function to replace binary products in a given expression constraints */
static
SCIP_RETCODE replaceBinaryProducts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_HASHMAP*         exprmap,            /**< map to remember generated variables for visited product expressions */
   SCIP_EXPRITER* it,               /**< expression iterator */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to update the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* expr;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(exprmap != NULL);
   assert(it != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   SCIPdebugMsg(scip, "  check constraint %s\n", SCIPconsGetName(cons));

   for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
   {
      SCIP_EXPR* newexpr = NULL;
      SCIP_EXPR* childexpr;
      int childexpridx;

      childexpridx = SCIPexpriterGetChildIdxDFS(it);
      assert(childexpridx >= 0 && childexpridx < SCIPexprGetNChildren(expr));
      childexpr = SCIPexpriterGetChildExprDFS(it);
      assert(childexpr != NULL);

      /* try to factorize variables in a sum expression that contains several products of binary variables */
      if( conshdlrdata->reformbinprodsfac > 1 )
      {
         SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, cons, childexpr,
            conshdlrdata->reformbinprodsfac, &newexpr, naddconss) );
      }

      /* try to create an expression that represents a product of binary variables */
      if( newexpr == NULL )
      {
         SCIP_CALL( getBinaryProductExpr(scip, conshdlr, exprmap, childexpr, &newexpr, naddconss, nchgcoefs) );
      }

      if( newexpr != NULL )
      {
         assert(naddconss == NULL || *naddconss > 0 || nchgcoefs == NULL || *nchgcoefs > 0);

         /* replace product expression */
         SCIP_CALL( SCIPreplaceExprChild(scip, expr, childexpridx, newexpr) );

         /* note that the expression has been captured by getBinaryProductExpr and SCIPreplaceExprChild */
         SCIP_CALL( SCIPreleaseExpr(scip, &newexpr) );

         /* mark the constraint to not be simplied anymore */
         consdata->issimplified = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** reformulates products of binary variables during presolving in the following way:
 *
 * Let sum_{i,j} Q_ij x_i x_j be a subexpression that only contains binary variables. Each term x_i x_j is
 * reformulated with the help of an extra (implicit integer) variable z_ij in {0,1}:
 *
 *    z_ij <= x_i, z_ij <= x_j, x_i + x_j - z_ij <= 1
 *
 * Before reformulating x_i x_j in this way, it is checked whether there is a clique that contains x_i and x_j. These
 * cliques allows for a better reformulation. There are four cases:
 *
 *    1. x_i + x_j <= 1 implies that x_i x_j = 0
 *
 *    2. x_i + (1 - x_j) <= 1 implies x_i x_j = x_i
 *
 *    3. (1 - x_i) + x_j <= 1 implies x_i x_j = x_j
 *
 *    4. (1 - x_i) + (1 - x_j) <= 1 implies x_i x_j = x_i + x_j - 1
 *
 * The reformulation using z_ij or the cliques is implemented in getBinaryProductExpr().
 *
 * Introducing too many extra variables and constraints can have a negative impact on the performance (e.g., due to
 * slow probing). For this reason, it is checked in getFactorizedBinaryQuadraticExpr() whether sum_{i,j} Q_ij x_i x_j
 * contains large (>= reformbinprodsfac parameter) lower sums of the form x_i sum_{j} Q_ij x_j. Such a lower sum is
 * reformulated with only one extra variable w_i:
 *
 *    maxact := sum_j max{0, Q_ij}, minact := sum_j min{0, Q_ij}
 *    minact x_i <= w_i, w_i <= maxact x_i
 *    minact <= sum_j Q_ij x_j - w_i + minact x_i
 *    maxact >= sum_j Q_ij x_j - w_i + maxact x_i
 *
 * We mark w_i to be implicit integer if all Q_ij are integer. After each replacment of a lower sum, it
 * is checked whether there are enough terms left to factorize other binary variables. Lower sums with a larger number
 * of terms are prioritized.
 */
static
SCIP_RETCODE presolveBinaryProducts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< expression constraints */
   int                   nconss,             /**< total number of expression constraints */
   int*                  naddconss,          /**< pointer to store the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to store the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_HASHMAP* exprmap;
   SCIP_EXPRITER* it;
   int c;

   assert(conshdlr != NULL);

   /* no expression constraints or binary variables -> skip */
   if( nconss == 0 || SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;
   assert(conss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create expression hash map */
   SCIP_CALL( SCIPhashmapCreate(&exprmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* create expression iterator */
   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

   SCIPdebugMsg(scip, "call presolveBinaryProducts()\n");

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_EXPR* newexpr = NULL;

      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* try to reformulate the root expression */
      if( conshdlrdata->reformbinprodsfac > 1 )
      {
         SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, conss[c], consdata->expr,
            conshdlrdata->reformbinprodsfac, &newexpr, naddconss) );
      }

      /* release the root node if another expression has been found */
      if( newexpr != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );
         consdata->expr = newexpr;

         /* mark constraint to be not simplified anymore */
         consdata->issimplified = FALSE;
      }

      /* replace each product of binary variables separately */
      SCIP_CALL( replaceBinaryProducts(scip, conshdlr, conss[c], exprmap, it, naddconss, nchgcoefs) );
   }

   /* free memory */
   SCIPhashmapFree(&exprmap);
   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** simplifies expressions and replaces common subexpressions for a set of constraints
 * @todo put the constant to the constraint sides
 */
static
SCIP_RETCODE canonicalizeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< total number of constraints */
   SCIP_PRESOLTIMING     presoltiming,       /**< presolve timing (SCIP_PRESOLTIMING_ALWAYS if not in presolving) */
   SCIP_Bool*            infeasible,         /**< buffer to store whether infeasibility has been detected */
   int*                  ndelconss,          /**< counter to add number of deleted constraints, or NULL */
   int*                  naddconss,          /**< counter to add number of added constraints, or NULL */
   int*                  nchgcoefs           /**< counter to add number of changed coefficients, or NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int* nlockspos;
   int* nlocksneg;
   SCIP_Bool havechange;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* update number of canonicalize calls */
   ++(conshdlrdata->ncanonicalizecalls);

   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->canonicalizetime) );

   *infeasible = FALSE;

   /* set havechange to TRUE in the first call of canonicalize; otherwise we might not replace common subexpressions */
   havechange = conshdlrdata->ncanonicalizecalls == 1;

   /* free nonlinear handlers information from expressions */  /* TODO can skip this in first presolve round */
   SCIP_CALL( deinitSolve(scip, conshdlr, conss, nconss) );

   /* allocate memory for storing locks of each constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &nlockspos, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksneg, nconss) );

   /* unlock all constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* remember locks */
      nlockspos[i] = consdata->nlockspos;
      nlocksneg[i] = consdata->nlocksneg;

      /* remove locks */
      SCIP_CALL( addLocks(scip, conss[i], -consdata->nlockspos, -consdata->nlocksneg) );
      assert(consdata->nlockspos == 0);
      assert(consdata->nlocksneg == 0);
   }

#ifndef NDEBUG
   /* check whether all locks of each expression have been removed */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_EXPR* expr;
      SCIP_EXPRITER* it;

      SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      SCIP_CALL( SCIPexpriterInit(it, consdata->expr, SCIP_EXPRITER_RTOPOLOGIC, TRUE) );
      for( expr = consdata->expr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         assert(expr != NULL);
         assert(expr->nlocksneg == 0);
         assert(expr->nlockspos == 0);
      }
      SCIPfreeExpriter(&it);
   }
#endif

   /* reformulate products of binary variables */
   if( conshdlrdata->reformbinprods && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING
      && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      int tmpnaddconss = 0;
      int tmpnchgcoefs = 0;

      /* call this function before simplification because expressions might not be simplified after reformulating
       * binary products; the detection of some nonlinear handlers might assume that expressions are simplified
       */
      SCIP_CALL( presolveBinaryProducts(scip, conshdlr, conss, nconss, &tmpnaddconss, &tmpnchgcoefs) );

      /* update counters */
      if( naddconss != NULL )
         *naddconss = tmpnaddconss;
      if( nchgcoefs != NULL )
         *nchgcoefs = tmpnchgcoefs;

      /* check whether at least one expression has changed */
      if( tmpnaddconss + tmpnchgcoefs > 0 )
         havechange = TRUE;
   }

   for( i = 0; i < nconss; ++i )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* call simplify for each expression */
      if( !consdata->issimplified && consdata->expr != NULL )
      {
         SCIP_EXPR* simplified;
         SCIP_Bool changed;

         changed = FALSE;
         SCIP_CALL( SCIPsimplifyExpr(scip, conshdlr, consdata->expr, &simplified, &changed, infeasible) );
         consdata->issimplified = TRUE;

         if( changed )
            havechange = TRUE;

         /* If root expression changed, then we need to take care updating the locks as well (the consdata is the one holding consdata->expr "as a child").
          * If root expression did not change, some subexpression may still have changed, but the locks were taking care of in the corresponding SCIPreplaceExprChild() call.
          */
         if( simplified != consdata->expr )
         {
            assert(changed);

            /* release old expression */
            SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );

            /* store simplified expression */
            consdata->expr = simplified;
         }
         else
         {
            /* The simplify captures simplified in any case, also if nothing has changed.
             * Therefore, we have to release it here.
             */
            SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
         }

         if( *infeasible )
            break;

         /* scale constraint sides */
         SCIP_CALL( scaleConsSides(scip, conshdlr, conss[i], &changed) );

         if( changed )
            havechange = TRUE;

         /* handle constant root expression; either the problem is infeasible or the constraint is redundant */
         if( consdata->expr->exprhdlr == SCIPgetExprHdlrValue(conshdlr) )
         {
            SCIP_Real value = SCIPexprvalGetValue(consdata->expr);
            if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasNegative(scip, value - consdata->lhs)) ||
                (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasPositive(scip, value - consdata->rhs)) )
            {
               SCIPdebugMsg(scip, "<%s> with constant expression found infeasible\n", SCIPconsGetName(conss[i]));
               SCIPdebugPrintCons(scip, conss[i], NULL);
               *infeasible = TRUE;
               break;
            }
            else
            {
               SCIP_CALL( addLocks(scip, conss[i], nlockspos[i], nlocksneg[i]) );
               SCIP_CALL( SCIPdelCons(scip, conss[i]) );
               if( ndelconss != NULL )
                  ++*ndelconss;
               havechange = TRUE;
            }
         }
      }
   }

   /* replace common subexpressions */
   if( havechange && !*infeasible )
   {
      SCIP_CONS** consssorted;

      SCIP_CALL( replaceCommonSubexpressions(scip, conss, nconss) );

      /* FIXME: this is a dirty hack for updating the variable expressions stored inside an expression which might have
       * been changed after simplification; now we completely recollect all variable expression and variable events
       */

      /* Each variable stores the constraints for which it catched varbound events sorted by the constraint index.
       * Thus, for performance reasons, it is better to call dropVarEvents in descending order of constraint index.
       */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consssorted, conss, nconss) );
      SCIPsortPtr((void**)consssorted, SCIPcompareIndexConsNonlinear, nconss);

      for( i = nconss-1; i >= 0; --i )
      {
         assert(i == 0 || SCIPcompareIndexConsNonlinear((void*)consssorted[i-1], (void*)consssorted[i]) < 0);
         if( SCIPconsIsDeleted(consssorted[i]) )
            continue;

         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
         SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(consssorted[i])) );
      }
      for( i = 0; i < nconss; ++i )
      {
         if( SCIPconsIsDeleted(consssorted[i]) )
            continue;

         SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(consssorted[i])) );
         SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
      }

      SCIPfreeBufferArray(scip, &consssorted);

      /* forbid multiaggregation for nonlinear variables again (in case new variables appeared now)
       * a multiaggregation of a nonlinear variable can yield to a large increase in expressions due to
       * expanding terms in simplify, e.g. ,(sum_i x_i)^2, so we just forbid these
       */
      SCIP_CALL( forbidNonlinearVariablesMultiaggration(scip, conshdlr, conss, nconss) );
   }

   /* restore locks */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsDeleted(conss[i]) )
         continue;

      SCIP_CALL( addLocks(scip, conss[i], nlockspos[i], nlocksneg[i]) );
   }

   /* run nlhdlr detect if in presolving stage (that is, not in exitpre)
    * TODO can we skip this in presoltiming fast?
    */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !*infeasible )
   {
      /* reset one of the number of detections counter to count only current presolving round */
      for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
         conshdlrdata->nlhdlrs[i]->ndetectionslast = 0;

      SCIP_CALL( initSolve(scip, conshdlr, conss, nconss) );
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &nlocksneg);
   SCIPfreeBufferArray(scip, &nlockspos);

   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->canonicalizetime) );

   return SCIP_OKAY;
}

/** given a cons_expr expression, creates an equivalent classic (nlpi-) expression */
static
SCIP_RETCODE makeClassicExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   sourceexpr,         /**< expression to convert */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to created expression */
   SCIP_EXPR**  varexprs,           /**< variable expressions that might occur in expr, their position in this array determines the varidx */
   int                   nvarexprs           /**< number of variable expressions */
   )
{
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_EXPR** children = NULL;
   int nchildren;
   int c;

   assert(scip != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);

   exprhdlr = SCIPexprGetHdlr(sourceexpr);
   nchildren = SCIPexprGetNChildren(sourceexpr);

   /* collect children expressions from children, if any */
   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );
      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( makeClassicExpr(scip, SCIPexprGetChildren(sourceexpr)[c], &children[c], varexprs, nvarexprs) );
         assert(children[c] != NULL);
      }
   }

   /* create target expression */
   if( strcmp(SCIPexprhdlrGetName(exprhdlr), "var") == 0 )
   {
      int varidx;

      /* find variable expression in varexprs array
       * the position in the array determines the index of the variable in the classic expression
       * TODO if varexprs are sorted, then can do this more efficient
       */
      for( varidx = 0; varidx < nvarexprs; ++varidx )
         if( varexprs[varidx] == sourceexpr )
            break;
      assert(varidx < nvarexprs);

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_VARIDX, varidx) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "val") == 0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_CONST, SCIPexprvalGetValue(sourceexpr)) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "sum") == 0 )
   {
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), targetexpr, nchildren, children, SCIPexprsumGetCoefs(sourceexpr), SCIPexprsumGetConstant(sourceexpr)) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "pow") == 0 )
   {
      SCIP_Real exponent;

      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);

      exponent = SCIPexprpowGetExponent(sourceexpr);
      if( EPSISINT(exponent, 0.0) )  /*lint !e835*/
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_INTPOWER, *children, (int)exponent) );
      }
      else
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_REALPOWER, *children, exponent) );
      }
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "signpower") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_SIGNPOWER, *children,
         SCIPexprpowGetExponent(sourceexpr)) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "prod") == 0 )
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, SCIPexprprodGetCoef(sourceexpr), nchildren, NULL, NULL) );
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), targetexpr, nchildren, children, 1, &monomial, 0.0, FALSE) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "abs") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_ABS, children[0]) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "exp") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_EXP, children[0]) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "log") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_LOG, children[0]) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "sin") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_SIN, children[0]) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "cos") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_COS, children[0]) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "entropy") == 0 )
   {
      SCIP_EXPR* childcopy;
      SCIP_Real minusone = -1.0;

      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &childcopy, children[0]) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &childcopy, SCIP_EXPR_LOG, childcopy) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_MUL, children[0], childcopy) );
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), targetexpr, 1, targetexpr, &minusone, 0.0) );
   }
   else
   {
      SCIPerrorMessage("unsupported expression handler <%s>, cannot convert to classical expression\n", SCIPexprhdlrGetName(exprhdlr));
      return SCIP_ERROR;
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
}

/** create a nonlinear row representation of an expr constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< expression constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* classicexpr = NULL;
   SCIP_VAR** nlvars = NULL;
   int nnlvars = 0;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* better curvature info will be set in INITSOL just before nlrow is added to NLP */
   SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
         0, NULL, NULL, 0, NULL, 0, NULL, NULL, consdata->lhs, consdata->rhs, SCIP_EXPRCURV_UNKNOWN) );

   if( consdata->expr == NULL )
      return SCIP_OKAY;

   if( SCIPexprGetHdlr(consdata->expr) == conshdlrdata->exprsumhdlr )
   {
      /* if root is a sum, then split into linear, quadratic, and expression */
      SCIP_EXPR* child;
      SCIP_Real* coefs;

      /* constant term of sum */
      SCIP_CALL( SCIPchgNlRowConstant(scip, consdata->nlrow, SCIPexprsumGetConstant(consdata->expr)) );

      coefs = SCIPexprsumGetCoefs(consdata->expr);

      for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
      {
         child = SCIPexprGetChildren(consdata->expr)[i];

         if( SCIPisExprVar(child) )
         {
            /* linear term */
            SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(child), coefs[i]) );
         }
         else if( SCIPexprGetHdlr(child) == conshdlrdata->exprpowhdlr &&
            SCIPexprpowGetExponent(child) == 2.0 &&
            SCIPisExprVar(SCIPexprGetChildren(child)[0]) )
         {
            /* square term  */
            SCIP_QUADELEM quadelem;

            quadelem.idx1 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(child)[0]));
            if( quadelem.idx1 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(child)[0])) );
               quadelem.idx1 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }
            quadelem.idx2 = quadelem.idx1;
            quadelem.coef = coefs[i];

            SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
         }
         else if( SCIPexprGetHdlr(child) == conshdlrdata->exprprodhdlr &&
            SCIPexprGetNChildren(child) == 2 &&
            SCIPisExprVar(SCIPexprGetChildren(child)[0]) &&
            SCIPisExprVar(SCIPexprGetChildren(child)[1]) )
         {
            /* bilinear term */
            SCIP_QUADELEM quadelem;

            quadelem.idx1 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(child)[0]));
            if( quadelem.idx1 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(child)[0])) );
               quadelem.idx1 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }

            quadelem.idx2 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(child)[1]));
            if( quadelem.idx2 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(child)[1])) );
               quadelem.idx2 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }

            quadelem.coef = coefs[i];

            SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
         }
         else
         {
            /* general nonlinear term */
            SCIP_EXPR* classicchild;

            /* make classic expression of child i */
            SCIP_CALL( makeClassicExpr(scip, child, &classicchild, consdata->varexprs, consdata->nvarexprs) );

            /* create or extend classicexpr */
            if( classicexpr == NULL )
            {
               SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &classicexpr, 1, &classicchild, coefs + i, 0.0) );
            }
            else
            {
               SCIP_CALL( SCIPexprAddToLinear(SCIPblkmem(scip), classicexpr, 1, &coefs[i], &classicchild, 0.0) );
            }
         }
      }

      if( classicexpr != NULL )
      {
         /* reindex variables in classicexpr so that only used variables are left */
         int* varsusage;
         int* reindexvars;

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &nlvars, consdata->nvarexprs) );
         SCIP_CALL( SCIPallocBufferArray(scip, &reindexvars, consdata->nvarexprs) );
         SCIP_CALL( SCIPallocClearBufferArray(scip, &varsusage, consdata->nvarexprs) );

         /* get count how often variables are used in expr */
         SCIPexprGetVarsUsage(classicexpr, varsusage);

         /* sort out unused variables and collect and reindex remaining variables */
         nnlvars = 0;
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            if( varsusage[i] == 0 )
            {
               reindexvars[i] = -1;
            }
            else
            {
               reindexvars[i] = nnlvars;
               nlvars[nnlvars] = SCIPexprvarGetVar(consdata->varexprs[i]);
               ++nnlvars;
            }
         }

         SCIPexprReindexVars(classicexpr, reindexvars);

         SCIPfreeBufferArray(scip, &varsusage);
         SCIPfreeBufferArray(scip, &reindexvars);
      }
   }
   else if( SCIPexprGetHdlr(consdata->expr) == conshdlrdata->exprpowhdlr &&
      SCIPexprpowGetExponent(consdata->expr) == 2.0 &&
      SCIPisExprVar(SCIPexprGetChildren(consdata->expr)[0]) )
   {
      /* if root is a x^2, then set the quadratic part of the nlrow */
      SCIP_QUADELEM quadelem;

      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(consdata->expr)[0])) );
      quadelem.idx1 = 0;
      quadelem.idx2 = 0;
      quadelem.coef = 1.0;

      SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
   }
   else if( SCIPexprGetHdlr(consdata->expr) == conshdlrdata->exprprodhdlr &&
      SCIPexprGetNChildren(consdata->expr) == 2 &&
      SCIPisExprVar(SCIPexprGetChildren(consdata->expr)[0]) &&
      SCIPisExprVar(SCIPexprGetChildren(consdata->expr)[1]) )
   {
      /* if root is a bilinear term x*y, then set the quadratic part of the nlrow */
      SCIP_QUADELEM quadelem;

      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(consdata->expr)[0])) );
      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPexprvarGetVar(SCIPexprGetChildren(consdata->expr)[1])) );

      quadelem.idx1 = 0;
      quadelem.idx2 = 1;
      quadelem.coef = 1.0;

      SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, quadelem) );
   }
   else
   {
      /* make classic expression */
      SCIP_CALL( makeClassicExpr(scip, consdata->expr, &classicexpr, consdata->varexprs, consdata->nvarexprs) );

      /* collect variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &nlvars, consdata->nvarexprs) );

      nnlvars = consdata->nvarexprs;
      for( i = 0; i < consdata->nvarexprs; ++i )
         nlvars[i] = SCIPexprvarGetVar(consdata->varexprs[i]);
   }
   assert((classicexpr != NULL) == (nlvars != NULL));

   if( classicexpr != NULL )
   {
      /* make classic expression tree */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, classicexpr, nnlvars, 0, NULL) );

      /* set variables in expression tree */
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, nnlvars, nlvars) );
      SCIPfreeBufferArray(scip, &nlvars);

      /* add expression tree in nlrow (this will make a copy) */
      SCIP_CALL( SCIPsetNlRowExprtree(scip, consdata->nlrow, exprtree) );

      /* free exprtree */
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** compares nonlinear handler by detection priority
 *
 * if handlers have same detection priority, then compare by name
 */
static
int nlhdlrCmp(
   void*                 hdlr1,              /**< first handler */
   void*                 hdlr2               /**< second handler */
   )
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(hdlr1 != NULL);
   assert(hdlr2 != NULL);

   h1 = (SCIP_NLHDLR*)hdlr1;
   h2 = (SCIP_NLHDLR*)hdlr2;

   if( h1->detectpriority != h2->detectpriority )
      return (int)(h1->detectpriority - h2->detectpriority);

   return strcmp(h1->name, h2->name);
}

/** creates auxiliary variable for a given expression
 *
 * @note for a variable expression it does nothing
 * @note this function can only be called in stage SCIP_STAGE_SOLVING
 */
static
SCIP_RETCODE createAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr                /**< expression */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VARTYPE vartype;
   SCIP_INTERVAL activity;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(expr->nauxvaruses > 0);

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
   {
      SCIPerrorMessage("it is not possible to create auxiliary variables during stage=%d\n", SCIPgetStage(scip));
      return SCIP_INVALIDCALL;
   }

   /* if we already have auxvar, then do nothing */
   if( expr->auxvar != NULL )
      return SCIP_OKAY;

   /* if expression is a variable-expression, then do nothing */
   if( expr->exprhdlr == SCIPgetExprHdlrVar(conshdlr) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);

   /* it doesn't harm much to have an auxvar for a constant, as this can be handled well by the default hdlr,
    * but it usually indicates a missing simplify
    * if we find situations where we need to have an auxvar for a constant, then remove this assert
    */
   assert(expr->exprhdlr != SCIPgetExprHdlrValue(conshdlr));

   /* create and capture auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "auxvar_%s_%d", expr->exprhdlr->name, conshdlrdata->auxvarid);
   ++conshdlrdata->auxvarid;

   /* type of auxiliary variable depends on integrality information of the expression */
   vartype = SCIPexprIsIntegral(expr) ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS;

   /* get activity of expression to initialize variable bounds */
   SCIP_CALL( SCIPevalExprActivity(scip, conshdlr, expr, &activity, TRUE) );
   /* we cannot handle a domain error here at the moment, but it seems unlikely that it could occur
    * if it appear, then we could change code to handle this properly, but for now we just ensure that we continue correctly
    * and abort in debug mode only
    */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity) )
   {
      SCIPABORT();
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &activity);
   }

   SCIP_CALL( SCIPcreateVarBasic(scip, &expr->auxvar, name, MAX( -SCIPinfinity(scip), activity.inf ),
      MIN( SCIPinfinity(scip), activity.sup ), 0.0, vartype) ); /*lint !e666*/

   /* mark the auxiliary variable to be added for the relaxation only
    * this prevents SCIP to create linear constraints from cuts or conflicts that contain auxiliary variables,
    * or to copy the variable to a subscip
    */
   SCIPvarMarkRelaxationOnly(expr->auxvar);

   SCIP_CALL( SCIPaddVar(scip, expr->auxvar) );

   SCIPdebugMsg(scip, "added auxiliary variable <%s> [%g,%g] for expression %p\n", SCIPvarGetName(expr->auxvar), SCIPvarGetLbGlobal(expr->auxvar), SCIPvarGetUbGlobal(expr->auxvar), (void*)expr);

   /* add variable locks in both directions
    * TODO should be sufficient to lock only according to expr->nlockspos/neg,
    *   but then we need to also update the auxvars locks when the expr locks change
    */
   SCIP_CALL( SCIPaddVarLocks(scip, expr->auxvar, 1, 1) );

#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      /* store debug solution value of auxiliary variable
       * assumes that expression has been evaluated in debug solution before
       */
      SCIP_CALL( SCIPdebugAddSolVal(scip, expr->auxvar, SCIPexprGetEvalValue(expr)) );
   }
#endif

   return SCIP_OKAY;
}

/** initializes separation for constraint
 *
 * - ensures that activities are uptodate in all expressions
 * - creates auxiliary variables where required
 * - calls propExprDomains to possibly tighten auxvar bounds
 * - calls separation initialization callback of nlhdlrs
 */
static
SCIP_RETCODE initSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether the problem is infeasible or not */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_INTERVAL activity;
   SCIP_RESULT result;
   int nreductions = 0;
   int c, e;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* start with new propbounds (just to be sure, should not be needed) */
   ++conshdlrdata->curpropboundstag;

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   /* first ensure activities are uptodate and create auxvars */
   *infeasible = FALSE;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_SOL* debugsol;

         SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

         if( debugsol != NULL ) /* it can be compiled WITH_DEBUG_SOLUTION, but still no solution given */
         {
            /* evaluate expression in debug solution, so we can set the solution value of created auxiliary variables
             * in createAuxVar()
             */
            SCIP_CALL( SCIPevalExpr(scip, conshdlr, consdata->expr, debugsol, 0) );
         }
      }
#endif

      /* ensure we have a valid activity for auxvars and propExprDomains() call below */
      SCIP_CALL( SCIPevalExprActivity(scip, conshdlr, consdata->expr, &activity, TRUE) );

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         if( expr->nauxvaruses > 0 )
         {
            SCIP_CALL( createAuxVar(scip, conshdlr, expr) );
         }
      }

      if( consdata->expr->auxvar != NULL )
      {
         SCIPdebugMsg(scip, "tighten auxvar <%s> bounds using constraint sides [%g,%g]\n",
               SCIPvarGetName(consdata->expr->auxvar), consdata->lhs, consdata->rhs);
         /* change the bounds of the auxiliary variable of the root node to [lhs,rhs] */
         SCIP_CALL( SCIPtightenVarLb(scip, consdata->expr->auxvar, consdata->lhs, TRUE, infeasible, NULL) );
         if( *infeasible )
         {
            SCIPdebugMsg(scip, "infeasibility detected while tightening auxvar lb (%g) using lhs of constraint (%g)\n",
               SCIPvarGetLbLocal(consdata->expr->auxvar), consdata->lhs);
            break;
         }

         SCIP_CALL( SCIPtightenVarUb(scip, consdata->expr->auxvar, consdata->rhs, TRUE, infeasible, NULL) );
         if( *infeasible )
         {
            SCIPdebugMsg(scip, "infeasibility detected while tightening auxvar ub (%g) using rhs of constraint (%g)\n",
               SCIPvarGetUbLocal(consdata->expr->auxvar), consdata->rhs);
            break;
         }
      }
   }

   /* now run a special version of reverseprop to ensure that important bound information (like function domains) is stored in bounds of auxvars,
    * since sometimes they cannot be recovered from activity evaluation even after some rounds of domain propagation
    * (e.g., log(x*y), which becomes log(w), w=x*y
    *  log(w) implies w >= 0, but we may not be able to derive bounds on x and y such that w >= 0 is ensured)
    */
   SCIP_CALL( propExprDomains(scip, conshdlr, conss, nconss, &result, &nreductions) );
   if( result == SCIP_CUTOFF )
      *infeasible = TRUE;

   /* now call initsepa of nlhdlrs
    * TODO skip if !SCIPconsIsInitial(conss[c]) ?
    *   but at the moment, initSepa() is called from INITLP anyway, so we have SCIPconsIsInitial(conss[c]) anyway
    */
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   for( c = 0; c < nconss && !*infeasible; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it) && !*infeasible; expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         if( expr->nauxvaruses == 0 )
            continue;

         for( e = 0; e < expr->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;
            SCIP_Bool underestimate;
            SCIP_Bool overestimate;
            assert(expr->enfos[e] != NULL);

            /* skip if initsepa was already called, e.g., because this expression is also part of a constraint
             * which participated in a previous initSepa() call
             */
            if( expr->enfos[e]->issepainit )
               continue;

            /* only call initsepa if it will actually separate */
            if( (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPABOTH) == 0 )
               continue;

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* only init sepa if there is an initsepa callback */
            if( !SCIPnlhdlrHasInitSepa(nlhdlr) )
               continue;

            /* check whether expression needs to be under- or overestimated */
            overestimate = SCIPgetExprNLocksNegNonlinear(expr) > 0;
            underestimate = SCIPgetExprNLocksPosNonlinear(expr) > 0;
            assert(underestimate || overestimate);

            SCIPdebugMsg(scip, "initsepa under=%u over=%u for expression %p\n", underestimate, overestimate, (void*)expr);

            /* call the separation initialization callback of the nonlinear handler */
            SCIP_CALL( SCIPcallNlhdlrInitSepaNonlinear(scip, conshdlr, conss[c], nlhdlr, expr,
               expr->enfos[e]->nlhdlrexprdata, overestimate, underestimate, infeasible) );
            expr->enfos[e]->issepainit = TRUE;

            if( *infeasible )
            {
               /* stop everything if we detected infeasibility */
               SCIPdebugMsg(scip, "detect infeasibility for constraint %s during initsepa()\n", SCIPconsGetName(conss[c]));
               break;
            }
         }
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** gets weight of variable when splitting violation score onto several variables in an expression */
static
SCIP_Real getViolSplitWeight(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_VAR*               var,              /**< variable */
   SCIP_SOL*               sol               /**< current solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   switch( conshdlrdata->branchviolsplit )
   {
      case 'e' :  /* evenly: everyone gets the same score */
         return 1.0;

      case 'm' :  /* midness of solution: 0.5 if in middle of domain, 0.05 if close to lower or upper bound */
      {
         SCIP_Real weight;
         weight = MIN(SCIPgetSolVal(scip, sol, var) - SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var) - SCIPgetSolVal(scip, sol, var)) / (SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var)); /*lint !e666*/
         return MAX(0.05, weight);
      }

      case 'd' :  /* domain width */
         return SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);

      case 'l' :  /* logarithmic domain width: log-scale if width is below 0.1 or above 10, otherwise actual width */
      {
         SCIP_Real width = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);
         assert(width > 0.0);
         if( width > 10.0 )
            return 10.0*log10(width);
         if( width < 0.1 )
            return 0.1/(-log10(width));
         return width;
      }

      default :
         SCIPerrorMessage("invalid value for parameter constraints/expr/branching/violsplit");
         SCIPABORT();
         return SCIP_INVALID;
   }
}

/** adds violation-branching score to a set of expressions, thereby distributing the score
 *
 * Each expression must either be a variable expression or have an aux-variable.
 *
 * If unbounded variables are present, each unbounded var gets an even score.
 * If no unbounded variables, then parameter constraints/expr/branching/violsplit decides weight for each var.
 */
static
void addConsExprExprsViolScore(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_EXPR**    exprs,            /**< expressions where to add branching score */
   int                     nexprs,           /**< number of expressions */
   SCIP_Real               violscore,        /**< violation-branching score to add to expression */
   SCIP_SOL*               sol,              /**< current solution */
   SCIP_Bool*              success           /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_VAR* var;
   SCIP_Real weight;
   SCIP_Real weightsum = 0.0; /* sum of weights over all candidates with bounded domain */
   int nunbounded = 0;  /* number of candidates with unbounded domain */
   int i;

   assert(exprs != NULL);
   assert(nexprs > 0);
   assert(success != NULL);

   if( nexprs == 1 )
   {
      SCIPaddExprViolScoreNonlinear(scip, conshdlr, exprs[0], violscore);
      *success = TRUE;
      return;
   }

   for( i = 0; i < nexprs; ++i )
   {
      var = SCIPgetExprAuxVarNonlinear(exprs[i]);
      assert(var != NULL);

      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
         ++nunbounded;
      else if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         weightsum += getViolSplitWeight(scip, conshdlr, var, sol);
   }

   *success = FALSE;
   for( i = 0; i < nexprs; ++i )
   {
      var = SCIPgetExprAuxVarNonlinear(exprs[i]);
      assert(var != NULL);

      if( nunbounded > 0 )
      {
         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIPaddExprViolScoreNonlinear(scip, conshdlr, exprs[i], violscore / nunbounded);
            *success = TRUE;
         }
      }
      else if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
      {
         assert(weightsum > 0.0);

         weight = getViolSplitWeight(scip, conshdlr, var, sol);
         SCIPaddExprViolScoreNonlinear(scip, conshdlr, exprs[i], violscore * weight / weightsum);
         SCIPdebugMsg(scip, "add score %g (%g%% of %g) to <%s>[%g,%g]\n", violscore * weight / weightsum,
            100*weight / weightsum, violscore,
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         *success = TRUE;
      }
      else
      {
         SCIPdebugMsg(scip, "skip score for fixed variable <%s>[%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      }
   }
}

/** adds violation-branching score to children of expression for given auxiliary variables
 *
 * Iterates over the successors of expr to find expressions that are associated with one of the given auxiliary variables.
 * Adds violatoin-branching scores to all found exprs by means of addConsExprExprsViolScore().
 *
 * @note This method may modify the given auxvars array by means of sorting.
 */
static
SCIP_RETCODE addConsExprExprViolScoresAuxVars(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_EXPR*     expr,             /**< expression where to start searching */
   SCIP_Real               violscore,        /**< violation score to add to expression */
   SCIP_VAR**              auxvars,          /**< auxiliary variables for which to find expression */
   int                     nauxvars,         /**< number of auxiliary variables */
   SCIP_SOL*               sol,              /**< current solution (NULL for the LP solution) */
   SCIP_Bool*              success           /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_EXPRITER* it;
   SCIP_VAR* auxvar;
   SCIP_EXPR** exprs;
   int nexprs;
   int pos;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvars != NULL);
   assert(success != NULL);

   /* sort variables to make lookup below faster */
   SCIPsortPtr((void**)auxvars, SCIPvarComp, nauxvars);

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_BFS, FALSE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nauxvars) );
   nexprs = 0;

   for( expr = SCIPexpriterGetNext(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )  /*lint !e441*/
   {
      auxvar = SCIPgetExprAuxVarNonlinear(expr);
      if( auxvar == NULL )
         continue;

      /* if auxvar of expr is contained in auxvars array, add branching score to expr */
      if( SCIPsortedvecFindPtr((void**)auxvars, SCIPvarComp, auxvar, nauxvars, &pos) )
      {
         assert(auxvars[pos] == auxvar);

         SCIPdebugMsg(scip, "adding branchingscore for expr %p with auxvar <%s>\n", expr, SCIPvarGetName(auxvar));
         exprs[nexprs++] = expr;

         if( nexprs == nauxvars )
            break;
      }
   }

   SCIPfreeExpriter(&it);

   if( nexprs > 0 )
   {
      SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, conshdlr, exprs, nexprs, violscore, sol, success) );
   }
   else
      *success = FALSE;

   SCIPfreeBufferArray(scip, &exprs);

   return SCIP_OKAY;
}

/** registers all unfixed variables in violated constraints as branching candidates */
static
SCIP_RETCODE registerBranchingCandidatesAllUnfixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int c;
   int i;

   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nnotify != NULL);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* consider only violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      /* register all variables that have not been fixed yet */
      assert(consdata->varexprs != NULL);
      for( i = 0; i < consdata->nvarexprs; ++i )
      {
         var = SCIPexprvarGetVar(consdata->varexprs[i]);
         assert(var != NULL);

         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, getConsAbsViolation(conss[c]), SCIP_INVALID) );
            ++(*nnotify);
         }
      }
   }

   return SCIP_OKAY;
}

/** registers all variables in violated constraints with branching scores as external branching candidates */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            success             /**< buffer to store whether at least one branching candidate was added */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it = NULL;
   int c;

   assert(conshdlr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( SCIPgetBranchAuxNonlinear(scip, conshdlr) )
   {
      SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   }

   /* register external branching candidates */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->varexprs != NULL);

      /* consider only violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      if( !SCIPgetBranchAuxNonlinear(scip, conshdlr) )
      {
         int i;

         /* if not branching on auxvars, then violation-branching scores will have been added to original variables
          * only, so we can loop over variable expressions
          */
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_Real violscore;
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_VAR* var;

            violscore = SCIPgetExprViolScoreNonlinear(conshdlr, consdata->varexprs[i]);

            /* skip variable expressions that do not have a violation score */
            if( violscore == 0.0 )
               continue;

            var = SCIPexprvarGetVar(consdata->varexprs[i]);
            assert(var != NULL);

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* consider variable for branching if it has not been fixed yet */
            if( !SCIPisEQ(scip, lb, ub) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " add variable <%s>[%g,%g] as extern branching candidate "\
                        "with score %g\n", SCIPvarGetName(var), lb, ub, violscore); )
               SCIP_CALL( SCIPaddExternBranchCand(scip, var, violscore, SCIP_INVALID) );
               *success = TRUE;
            }
            else
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
            }

            /* invalidate violscore-tag, so that we do not register variables that appear in multiple constraints
             * several times as external branching candidate, see SCIPgetExprViolScoreNonlinear()
             */
            consdata->varexprs[i]->violscoretag = 0;
         }
      }
      else
      {
         SCIP_EXPR* expr;
         SCIP_VAR* var;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real violscore;

         for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
         {
            violscore = SCIPgetExprViolScoreNonlinear(conshdlr, expr);
            if( violscore == 0.0 )
               continue;

            /* if some nlhdlr added a branching score for this expression, then it considered this expression as a
             * variable, so this expression should either be an original variable or have an auxiliary variable
             */
            var = SCIPgetExprAuxVarNonlinear(expr);
            assert(var != NULL);

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* consider variable for branching if it has not been fixed yet */
            if( !SCIPisEQ(scip, lb, ub) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " add variable <%s>[%g,%g] as extern branching candidate "\
                        "with score %g\n", SCIPvarGetName(var), lb, ub, violscore); )

               SCIP_CALL( SCIPaddExternBranchCand(scip, var, violscore, SCIP_INVALID) );
               *success = TRUE;
            }
            else
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
            }
         }
      }
   }

   if( SCIPgetBranchAuxNonlinear(scip, conshdlr) )
      SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** collect branching candidates from violated constraints
 *
 * Fills array with expressions that serve as branching candidates.
 * Collects those expressions that have a branching score assigned and stores the score in the auxviol field of the
 * branching candidate.
 *
 * If branching on aux-variables is allowed, then iterate through expressions of violated constraints, otherwise iterate
 * through variable-expressions only.
 */
static
SCIP_RETCODE collectBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Real             maxrelconsviol,     /**< maximal scaled constraint violation */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   BRANCHCAND*           cands,              /**< array where to store candidates, must be at least SCIPgetNVars() long */
   int*                  ncands              /**< number of candidates found */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it = NULL;
   int c;
   int attempt;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cands != NULL);
   assert(ncands != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetBranchAuxNonlinear(scip, conshdlr) )
   {
      SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   }

   *ncands = 0;
   for( attempt = 0; attempt < 2; ++attempt )
   {
      /* collect branching candidates from violated constraints
       * in the first attempt, consider only constraints with large violation
       * in the second attempt, consider all remaining violated constraints
       */
      for( c = 0; c < nconss; ++c )
      {
         SCIP_Real consviol;

         assert(conss != NULL && conss[c] != NULL);

         /* consider only violated constraints */
         if( !isConsViolated(scip, conss[c]) )
            continue;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         assert(consdata->varexprs != NULL);

         SCIP_CALL( getConsRelViolation(scip, conss[c], &consviol, sol, soltag) );

         if( attempt == 0 && consviol < conshdlrdata->branchhighviolfactor * maxrelconsviol )
            continue;
         else if( attempt == 1 && consviol >= conshdlrdata->branchhighviolfactor * maxrelconsviol )
            continue;

         if( !SCIPgetBranchAuxNonlinear(scip, conshdlr) )
         {
            int i;

            /* if not branching on auxvars, then violation-branching scores will be available for original variables
             * only, so we can loop over variable expressions
             * unfortunately, we don't know anymore which constraint contributed the violation-branching score to the
             * variable, therefore we invalidate the score of a variable after processing it.
             */
            for( i = 0; i < consdata->nvarexprs; ++i )
            {
               SCIP_Real lb;
               SCIP_Real ub;

               /* skip variable expressions that do not have a valid violation score */
               if( conshdlrdata->enforound != consdata->varexprs[i]->violscoretag )
                  continue;

               var = SCIPexprvarGetVar(consdata->varexprs[i]);
               assert(var != NULL);

               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               /* skip already fixed variable */
               if( SCIPisEQ(scip, lb, ub) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
                  continue;
               }

               assert(*ncands + 1 < SCIPgetNVars(scip));
               cands[*ncands].expr = consdata->varexprs[i];
               cands[*ncands].auxviol = SCIPgetExprViolScoreNonlinear(conshdlr, consdata->varexprs[i]);
               ++(*ncands);

               /* invalidate violscore-tag, so that we do not register variables that appear in multiple constraints
                * several times as external branching candidate */
               consdata->varexprs[i]->violscoretag = 0;
            }
         }
         else
         {
            SCIP_EXPR* expr;
            SCIP_Real lb;
            SCIP_Real ub;

            for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
            {
               if( expr->violscoretag != conshdlrdata->enforound )
                  continue;

               /* if some nlhdlr added a branching score for this expression, then it considered this expression as
                * variables, so this expression should either be an original variable or have an auxiliary variable
                */
               var = SCIPgetExprAuxVarNonlinear(expr);
               assert(var != NULL);

               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               /* skip already fixed variable */
               if( SCIPisEQ(scip, lb, ub) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
                  continue;
               }

               assert(*ncands + 1 < SCIPgetNVars(scip));
               cands[*ncands].expr = expr;
               cands[*ncands].auxviol = SCIPgetExprViolScoreNonlinear(conshdlr, expr);
               ++(*ncands);
            }
         }
      }

      /* if we have branching candidates, then we don't need another attempt */
      if( *ncands > 0 )
         break;
   }

   if( SCIPgetBranchAuxNonlinear(scip, conshdlr) )
      SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** computes a branching score for a variable that reflects how important branching on this variable would be for
 * improving the dual bound from the LP relaxation
 *
 * Assume the Lagrangian for the current LP is something of the form
 *   L(x,z,lambda) = c'x + sum_i lambda_i (a_i'x - z_i + b_i) + ...
 * where x are the original variables, z the auxiliary variables,
 * and a_i'x - z_i + b_i <= 0 are the rows of the LP.
 *
 * Assume that a_i'x + b_i <= z_i was derived from some nonlinear constraint f(x) <= z and drop index i.
 * If we could have used not only an estimator, but the actual function f(x), then this would
 * have contributed lambda*(f(x) - z) to the Lagrangian function (though the value of z would be different).
 * Using a lot of handwaving, we claim that
 *   lambda_i * (f(x) - a_i'x + b_i)
 * is a value that can be used to quantity how much improving the estimator a'x + b <= z could change the dual bound.
 * If an estimator depended on local bounds, then it could be improved by branching.
 * We use row-is-local as proxy for estimator-depending-on-lower-bounds.
 *
 * To score a variable, we then sum the values lambda_i * (f(x) - a_i'x + b_i) for all rows in which the variable appears.
 * To scale, we divide by the LP objective value (if >1).
 *
 * TODO if we branch only on original variables, we neglect here estimators that are build on auxiliary variables
 *     these are affected by the bounds on original variables indirectly (through forward-propagation)
 * TODO if we branch also on auxiliary variables, then separating z from the x-variables in the row a'x+b <= z should happen
 *     in effect, we should go from the row to the expression for which it was generated and consider only variables that
 *     would also be branching candidates
 */
static
SCIP_Real getDualBranchscore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraints handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_COL* col;
   SCIP_ROW** rows;
   int nrows;
   int r;
   SCIP_Real dualscore;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(var != NULL);

   /* if LP not solved, then the dual branching score is not available */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return 0.0;

   /* if var is not in the LP, then the dual branching score is not available */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0.0;

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
      return 0.0;

   nrows = SCIPcolGetNLPNonz(col);  /* TODO there is a big warning on when not to use this method; is the check for SCIPcolIsInLP sufficient? */
   rows = SCIPcolGetRows(col);

   /* SCIPinfoMessage(scip, enfologfile, " dualscoring <%s>\n", SCIPvarGetName(var)); */

   /* aggregate duals from all rows from consexpr with non-zero dual
    * TODO: this is a quick-and-dirty implementation, and not used by default
    *   in the long run, this should be either removed or replaced by a proper implementation
    */
   dualscore = 0.0;
   for( r = 0; r < nrows; ++r )
   {
      SCIP_Real estimategap;
      const char* estimategapstr;

      /* rows from cuts that may be replaced by tighter ones after branching are the interesting ones
       * these would typically be local, unless they are created at the root node
       * so not check for local now, but trust that estimators that do not improve after branching will have an estimategap of 0
      if( !SCIProwIsLocal(rows[r]) )
         continue;
       */
      if( SCIProwGetOriginConshdlr(rows[r]) != conshdlr )
         continue;
      if( SCIPisZero(scip, SCIProwGetDualsol(rows[r])) )
         continue;

      estimategapstr = strstr(SCIProwGetName(rows[r]), "_estimategap=");
      if( estimategapstr == NULL ) /* gap not stored, maybe because it was 0 */
         continue;
      estimategap = atof(estimategapstr + 13);
      assert(estimategap >= 0.0);
      if( !SCIPisFinite(estimategap) || SCIPisHugeValue(scip, estimategap) )
         estimategap = SCIPgetHugeValue(scip);

      /* SCIPinfoMessage(scip, enfologfile, "  row <%s> contributes %g*|%g|: ", SCIProwGetName(rows[r]), estimategap, SCIProwGetDualsol(rows[r]));
      SCIP_CALL( SCIPprintRow(scip, rows[r], enfologfile) ); */

      dualscore += estimategap * REALABS(SCIProwGetDualsol(rows[r]));
   }

   /* divide by optimal value of LP for scaling */
   dualscore /= MAX(1.0, REALABS(SCIPgetLPObjval(scip)));  /*lint !e666*/

   return dualscore;
}

/** computes branching scores (including weighted score) for a set of candidates
 *
 * For each candidate in the array, compute and store the various branching scores (violation, pseudo-costs, vartype, domainwidth).
 * For pseudo-costs, it's possible that the score is not available, in which case cands[c].pscost will be set to SCIP_INVALID.
 *
 * For each score, compute the maximum over all candidates.
 *
 * Then compute for each candidate a "weighted" score using the weights as specified by parameters
 * and the scores as previously computed, but scale each score to be in [0,1], i.e., divide each score by the maximum
 * score all candidate.
 * Further divide by the sum of all weights where a score was available (even if the score was 0).
 *
 * For example:
 * - Let variable x have violation-score 10.0 and pseudo-cost-score 5.0.
 * - Let variable y have violation-score 12.0 but no pseudo-cost-score (because it hasn't yet been branched on sufficiently often).
 * - Assuming violation is weighted by 2.0 and pseudo-costs are weighted by 3.0.
 * - Then the weighted scores for x will be (2.0 * 10.0/12.0 + 3.0 * 5.0/5.0) / (2.0 + 3.0) = 0.9333.
 *   The weighted score for y will be (2.0 * 12.0/12.0) / 2.0 = 1.0.
 */
static
void scoreBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BRANCHCAND*           cands,              /**< branching candidates */
   int                   ncands,             /**< number of candidates */
   SCIP_SOL*             sol                 /**< solution to enforce (NULL for the LP solution) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   BRANCHCAND maxscore;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cands != NULL);
   assert(ncands > 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initialize counts to 0 */
   memset(&maxscore, 0, sizeof(BRANCHCAND));

   for( c = 0; c < ncands; ++c )
   {
      if( conshdlrdata->branchviolweight > 0.0 )
      {
         /* cands[c].auxviol was set in collectBranchingCandidates, so only update maxscore here */
         maxscore.auxviol = MAX(maxscore.auxviol, cands[c].auxviol);
      }

      if( conshdlrdata->branchdomainweight > 0.0 )
      {
         SCIP_Real domainwidth;
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         /* get domain width, taking infinity at 1e20 on purpose */
         domainwidth = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);

         /* domain-score is going to be log(2*infinity / domainwidth) if domain width >= 1
          * and log(2 * infinity *  MAX(epsilon, domainwidth)) for domain width < 1
          * the idea is to penalize very large and very small domains
          */
         if( domainwidth >= 1.0 )
            cands[c].domain = log10(2 * SCIPinfinity(scip) / domainwidth);
         else
            cands[c].domain = log10(2 * SCIPinfinity(scip) * MAX(SCIPepsilon(scip), domainwidth));  /*lint !e666*/

         maxscore.domain = MAX(cands[c].domain, maxscore.domain);
      }
      else
         cands[c].domain = 0.0;

      if( conshdlrdata->branchdualweight > 0.0 )
      {
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         cands[c].dual = getDualBranchscore(scip, conshdlr, var);
         maxscore.dual = MAX(cands[c].dual, maxscore.dual);
      }

      if( conshdlrdata->branchpscostweight > 0.0 && SCIPgetNObjVars(scip) > 0 )
      {
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
            cands[c].pscost = SCIP_INVALID;
         else
         {
            SCIP_Real brpoint;
            SCIP_Real pscostdown;
            SCIP_Real pscostup;
            char strategy;

            /* decide how to compute pseudo-cost scores
             * this should be consistent with the way how pseudo-costs are updated in the core, which is decided by
             * branching/lpgainnormalize for continuous variables and move in LP-value for non-continuous variables
             */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               strategy = conshdlrdata->branchpscostupdatestrategy;
            else
               strategy = 'l';

            brpoint = SCIPgetBranchingPoint(scip, var, SCIP_INVALID);

            /* branch_relpscost deems pscosts as reliable, if the pseudo-count is at least something between 1 and 4
             * or it uses some statistical tests involving SCIPisVarPscostRelerrorReliable
             * For here, I use a simple #counts >= branchpscostreliable.
             * TODO use SCIPgetVarPseudocostCount() instead?
             */
            if( SCIPgetVarPseudocostCountCurrentRun(scip, var, SCIP_BRANCHDIR_DOWNWARDS) >= conshdlrdata->branchpscostreliable )
            {
               switch( strategy )
               {
                  case 's' :
                     pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPvarGetUbLocal(var) - SCIPadjustedVarLb(scip, var, brpoint)));
                     break;
                  case 'd' :
                     pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPadjustedVarUb(scip, var, brpoint) - SCIPvarGetLbLocal(var)));
                     break;
                  case 'l' :
                     if( SCIPisInfinity(scip, SCIPgetSolVal(scip, sol, var)) )
                        pscostdown = SCIP_INVALID;
                     else if( SCIPgetSolVal(scip, sol, var) <= SCIPadjustedVarUb(scip, var, brpoint) )
                        pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0);
                     else
                        pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPgetSolVal(scip, NULL, var) - SCIPadjustedVarUb(scip, var, brpoint)));
                     break;
                  default :
                     SCIPerrorMessage("pscost update strategy %c unknown\n", strategy);
                     pscostdown = SCIP_INVALID;
               }
            }
            else
               pscostdown = SCIP_INVALID;

            if( SCIPgetVarPseudocostCountCurrentRun(scip, var, SCIP_BRANCHDIR_UPWARDS) >= conshdlrdata->branchpscostreliable )
            {
               switch( strategy )
               {
                  case 's' :
                     pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPadjustedVarUb(scip, var, brpoint) - SCIPvarGetLbLocal(var));
                     break;
                  case 'd' :
                     pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPvarGetUbLocal(var) - SCIPadjustedVarLb(scip, var, brpoint));
                     break;
                  case 'l' :
                     if( SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, var)) )
                        pscostup = SCIP_INVALID;
                     else if( SCIPgetSolVal(scip, NULL, var) >= SCIPadjustedVarLb(scip, var, brpoint) )
                        pscostup = SCIPgetVarPseudocostVal(scip, var, 0.0);
                     else
                        pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPadjustedVarLb(scip, var, brpoint) - SCIPgetSolVal(scip, NULL, var) );
                     break;
                  default :
                     SCIPerrorMessage("pscost update strategy %c unknown\n", strategy);
                     pscostup = SCIP_INVALID;
               }
            }
            else
               pscostup = SCIP_INVALID;

            if( pscostdown == SCIP_INVALID && pscostup == SCIP_INVALID )  /*lint !e777*/
               cands[c].pscost = SCIP_INVALID;
            else if( pscostdown == SCIP_INVALID )  /*lint !e777*/
               cands[c].pscost = pscostup;
            else if( pscostup == SCIP_INVALID )  /*lint !e777*/
               cands[c].pscost = pscostdown;
            else
               cands[c].pscost = SCIPgetBranchScore(scip, NULL, pscostdown, pscostup);  /* pass NULL for var to avoid multiplication with branch-factor */
         }

         if( cands[c].pscost != SCIP_INVALID )  /*lint !e777*/
            maxscore.pscost = MAX(cands[c].pscost, maxscore.pscost);
      }

      if( conshdlrdata->branchvartypeweight > 0.0 )
      {
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         switch( SCIPvarGetType(var) )
         {
            case SCIP_VARTYPE_BINARY :
               cands[c].vartype = 1.0;
               break;
            case SCIP_VARTYPE_INTEGER :
               cands[c].vartype = 0.1;
               break;
            case SCIP_VARTYPE_IMPLINT :
               cands[c].vartype = 0.01;
               break;
            case SCIP_VARTYPE_CONTINUOUS :
            default:
               cands[c].vartype = 0.0;
         }
         maxscore.vartype = MAX(cands[c].vartype, maxscore.vartype);
      }
   }

   /* now computed a weighted score for each candidate from the single scores
    * the single scores are scaled to be in [0,1] for this
    */
   for( c = 0; c < ncands; ++c )
   {
      SCIP_Real weightsum;

      ENFOLOG(
         SCIP_VAR* var;
         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         SCIPinfoMessage(scip, enfologfile, " scoring <%8s>[%7.1g,%7.1g]:(", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         )

      cands[c].weighted = 0.0;
      weightsum = 0.0;

      if( maxscore.auxviol > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchviolweight * cands[c].auxviol / maxscore.auxviol;
         weightsum += conshdlrdata->branchviolweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(viol)", conshdlrdata->branchviolweight, cands[c].auxviol / maxscore.auxviol); )
      }

      if( maxscore.domain > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchdomainweight * cands[c].domain / maxscore.domain;
         weightsum += conshdlrdata->branchdomainweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(domain)", conshdlrdata->branchdomainweight, cands[c].domain / maxscore.domain); )
      }

      if( maxscore.dual > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchdualweight * cands[c].dual / maxscore.dual;
         weightsum += conshdlrdata->branchdualweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(dual)", conshdlrdata->branchdualweight, cands[c].dual / maxscore.dual); )
      }

      /* use pseudo-costs, if we have some for at least half the candidates */
      if( maxscore.pscost > 0.0 )
      {
         if( cands[c].pscost != SCIP_INVALID )  /*lint !e777*/
         {
            cands[c].weighted += conshdlrdata->branchpscostweight * cands[c].pscost / maxscore.pscost;
            weightsum += conshdlrdata->branchpscostweight;

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(pscost)", conshdlrdata->branchpscostweight, cands[c].pscost / maxscore.pscost); )
         }
         else
         {
            /* do not add pscostscore, if not available, also do not add into weightsum */
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " +0.0*    n/a(pscost)"); )
         }
      }

      if( maxscore.vartype > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchvartypeweight * cands[c].vartype / maxscore.vartype;
         weightsum += conshdlrdata->branchvartypeweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%6.2g(vartype)", conshdlrdata->branchvartypeweight, cands[c].vartype / maxscore.vartype); )
      }
      assert(weightsum > 0.0);  /* we should have got at least one valid score */
      cands[c].weighted /= weightsum;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " ) / %g = %g\n", weightsum, cands[c].weighted); )
   }
}

/** compare two branching candidates by their weighted score
 *
 * if weighted score is equal, use variable index of (aux)var
 */
static
SCIP_DECL_SORTINDCOMP(branchcandCompare)
{
   BRANCHCAND* cands = (BRANCHCAND*)dataptr;

   if( cands[ind1].weighted != cands[ind2].weighted )  /*lint !e777*/
      return cands[ind1].weighted < cands[ind2].weighted ? -1 : 1;
   else
      return SCIPvarGetIndex(SCIPgetExprAuxVarNonlinear(cands[ind1].expr)) - SCIPvarGetIndex(SCIPgetExprAuxVarNonlinear(cands[ind2].expr));
}

/** do branching or register branching candidates */
static
SCIP_RETCODE branching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Real             maxrelconsviol,     /**< maximal scaled constraint violation */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_RESULT*          result              /**< pointer to store the result of branching */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   BRANCHCAND* cands;
   int ncands;
   SCIP_VAR* var;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

   assert(conshdlr != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->branchexternal )
   {
      /* just register branching candidates as external */
      SCIP_Bool success;

      SCIP_CALL( registerBranchingCandidates(scip, conshdlr, conss, nconss, &success) );
      if( success )
         *result = SCIP_INFEASIBLE;

      return SCIP_OKAY;
   }

   /* collect branching candidates and their auxviol-score */
   SCIP_CALL( SCIPallocBufferArray(scip, &cands, SCIPgetNVars(scip)) );
   SCIP_CALL( collectBranchingCandidates(scip, conshdlr, conss, nconss, maxrelconsviol, sol, soltag, cands, &ncands) );

   /* if no unfixed branching candidate in all violated constraint, then it's probably numerics that prevented us to separate or decide a cutoff
    * we will return here and let the fallbacks in consEnfo() decide how to proceed
    */
   if( ncands == 0 )
      goto TERMINATE;

   if( ncands > 1 )
   {
      /* if there are more than one candidate, then compute scores and select */
      int* perm;
      int c;
      int left;
      int right;
      SCIP_Real threshold;

      /* compute additional scores on branching candidates and weighted score */
      scoreBranchingCandidates(scip, conshdlr, cands, ncands, sol);

      /* sort candidates by weighted score */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, ncands) );
      SCIPsortDown(perm, branchcandCompare, (void*)cands, ncands);

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %d branching candidates <%s>(%g)...<%s>(%g)\n", ncands,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[0]].expr)), cands[perm[0]].weighted,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[ncands - 1]].expr)), cands[perm[ncands - 1]].weighted); )

      /* binary search to find first low-scored (score below branchhighscorefactor * maximal-score)  candidate */
      left = 0;
      right = ncands - 1;
      threshold = conshdlrdata->branchhighscorefactor * cands[perm[0]].weighted;
      while( left < right )
      {
         int mid = (left + right) / 2;
         if( cands[perm[mid]].weighted >= threshold )
            left = mid + 1;
         else
            right = mid;
      }
      assert(left <= ncands);

      if( left < ncands )
      {
         if( cands[perm[left]].weighted >= threshold )
         {
            assert(left + 1 == ncands || cands[perm[left + 1]].weighted < threshold);
            ncands = left + 1;
         }
         else
         {
            assert(cands[perm[left]].weighted < threshold);
            ncands = left;
         }
      }
      assert(ncands > 0);

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %d branching candidates <%s>(%g)...<%s>(%g) after removing low scores\n", ncands,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[0]].expr)), cands[perm[0]].weighted,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[ncands - 1]].expr)), cands[perm[ncands - 1]].weighted); )

      if( ncands > 1 )
      {
         /* choose at random from candidates 0..ncands-1 */
         if( conshdlrdata->branchrandnumgen == NULL )
         {
            SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->branchrandnumgen, BRANCH_RANDNUMINITSEED, TRUE) );
         }
         c = SCIPrandomGetInt(conshdlrdata->branchrandnumgen, 0, ncands - 1);
         var = SCIPgetExprAuxVarNonlinear(cands[perm[c]].expr);
      }
      else
         var = SCIPgetExprAuxVarNonlinear(cands[perm[0]].expr);

      SCIPfreeBufferArray(scip, &perm);
   }
   else
   {
      var = SCIPgetExprAuxVarNonlinear(cands[0].expr);
   }
   assert(var != NULL);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " branching on variable <%s>[%g,%g]\n", SCIPvarGetName(var),
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)); )

   SCIP_CALL( SCIPbranchVarVal(scip, var, SCIPgetBranchingPoint(scip, var, SCIP_INVALID), &downchild, &eqchild,
            &upchild) );
   if( downchild != NULL || eqchild != NULL || upchild != NULL )
      *result = SCIP_BRANCHED;
   else
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      *result = SCIP_REDUCEDDOM;

 TERMINATE:
   SCIPfreeBufferArray(scip, &cands);

   return SCIP_OKAY;
}

/** call enforcement or estimate callback of nonlinear handler
 *
 * Calls the enforcement callback, if available.
 * Otherwise, calls the estimate callback, if available, and constructs a cut from the estimator.
 *
 * If cut is weak, but estimator is not tight, tries to add branching candidates.
 */
static
SCIP_RETCODE enforceExprNlhdlr(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_NLHDLR* nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, /**< nonlinear handler data of expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_Bool             separated,          /**< whether another nonlinear handler already added a cut for this expression */
   SCIP_Bool             allowweakcuts,      /**< whether we allow for weak cuts */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement (and not just separation) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   assert(result != NULL);

   /* call enforcement callback of the nlhdlr */
   SCIP_CALL( SCIPcallNlhdlrEnfoNonlinear(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
            allowweakcuts, separated, inenforcement, result) );

   /* if it was not running (e.g., because it was not available) or did not find anything, then try with estimator callback */
   if( *result != SCIP_DIDNOTRUN && *result != SCIP_DIDNOTFIND )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    sepa of nlhdlr %s succeeded with result %d\n",
               SCIPnlhdlrGetName(nlhdlr), *result); )
      return SCIP_OKAY;
   }
   else
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    sepa of nlhdlr <%s> did not succeed with result %d\n", SCIPnlhdlrGetName(nlhdlr), *result); )
   }

   *result = SCIP_DIDNOTFIND;

   /* now call the estimator callback of the nlhdlr */
   if( SCIPnlhdlrHasEstimate(nlhdlr) )
   {
      SCIP_VAR* auxvar;
      SCIP_Bool sepasuccess = FALSE;
      SCIP_Bool branchscoresuccess = FALSE;
      SCIP_PTRARRAY* rowpreps;
      int minidx;
      int maxidx;
      int r;
      SCIP_ROWPREP* rowprep;

      SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );

      auxvar = SCIPgetExprAuxVarNonlinear(expr);
      assert(auxvar != NULL);

      SCIP_CALL( SCIPcallNlhdlrEstimateNonlinear(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
               SCIPgetSolVal(scip, sol, auxvar), rowpreps, &sepasuccess, inenforcement, &branchscoresuccess) );

      minidx = SCIPgetPtrarrayMinIdx(scip, rowpreps);
      maxidx = SCIPgetPtrarrayMaxIdx(scip, rowpreps);

      assert((sepasuccess && minidx <= maxidx) || (!sepasuccess && minidx > maxidx));

      if( !sepasuccess )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s failed\n",
                  SCIPnlhdlrGetName(nlhdlr)); )
      }

      for( r = minidx; r <= maxidx; ++r )
      {
         rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, r);

         assert(rowprep != NULL);

         /* complete estimator to cut */
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0) );

         /* add the cut and/or branching scores */
         SCIP_CALL( SCIPprocessRowprepNonlinear(scip, conshdlr, nlhdlr, cons, expr, rowprep, overestimate, auxvar,
               auxvalue, allowweakcuts, branchscoresuccess, inenforcement, sol, result) );

         SCIPfreeRowprep(scip, &rowprep);
      }

      SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   }

   return SCIP_OKAY;
}

/** tries to enforce violation in an expression by separation, bound tightening, or finding a branching candidate
 *
 * if not inenforcement, then we should be called by consSepa, and thus only try separation
 */
static
SCIP_RETCODE enforceExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraints handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Bool             allowweakcuts,      /**< whether we allow weak cuts */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement (and not just separation) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_Real origviol;
   SCIP_Bool underestimate;
   SCIP_Bool overestimate;
   SCIP_Real auxviol;
   SCIP_Bool auxunderestimate;
   SCIP_Bool auxoverestimate;
   SCIP_RESULT hdlrresult;
   int e;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr->auxvar != NULL);  /* there must be a variable attached to the expression in order to construct a cut here */
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* make sure that this expression has been evaluated */
   SCIP_CALL( SCIPevalExpr(scip, conshdlr, expr, sol, soltag) );

   /* decide whether under- or overestimate is required and get amount of violation */
   origviol = getExprAbsOrigViolation(scip, expr, sol, &underestimate, &overestimate);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* no sufficient violation w.r.t. the original variables -> skip expression */
   if( !overestimate && !underestimate )
   {
      return SCIP_OKAY;
   }

   /* check aux-violation w.r.t. each nonlinear handlers and try to enforce when there is a decent violation */
   for( e = 0; e < expr->nenfos; ++e )
   {
      SCIP_NLHDLR* nlhdlr;

      /* skip nlhdlr that do not want to participate in any separation */
      if( (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPABOTH) == 0 )
         continue;

      nlhdlr = expr->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
      SCIP_CALL( SCIPcallNlhdlrEvalauxNonlinear(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &expr->enfos[e]->auxvalue, sol) );
      ENFOLOG(
         SCIPinfoMessage(scip, enfologfile, "  expr ");
         SCIPprintExpr(scip, conshdlr, expr, enfologfile);
         SCIPinfoMessage(scip, enfologfile, " (%p): evalvalue %.15g auxvarvalue %.15g [%.15g,%.15g], nlhdlr <%s> " \
            "auxvalue: %.15g\n", (void*)expr, expr->evalvalue, SCIPgetSolVal(scip, sol, expr->auxvar),
            expr->activity.inf, expr->activity.sup, nlhdlr->name, expr->enfos[e]->auxvalue);
      )

      /* TODO if expr is root of constraint (consdata->expr == expr),
       * then compare auxvalue with constraint sides instead of auxvarvalue, as the former is what actually matters
       * that is, if auxvalue is good enough for the constraint to be satisfied, but when looking at evalvalue we see
       * the the constraint is violated, then some of the auxvars that nlhdlr uses is not having a good enough value,
       * so we should enforce in these auxiliaries first
       * if changing this here, we must also adapt analyzeViolation
       */

      auxviol = getExprAbsAuxViolation(scip, expr, expr->enfos[e]->auxvalue, sol, &auxunderestimate, &auxoverestimate);
      assert(auxviol >= 0.0);

      /* if aux-violation is much smaller than orig-violation, then better enforce further down in the expression first */
      if( !SCIPisInfinity(scip, auxviol) && auxviol < conshdlrdata->enfoauxviolfactor * origviol )  /*lint !e777*/
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   skip enforce using nlhdlr <%s> for expr %p (%s) with " \
                  "auxviolation %g << origviolation %g under:%d over:%d\n", nlhdlr->name, (void*)expr,
                  expr->exprhdlr->name, auxviol, origviol, underestimate, overestimate); )

         /* TODO expr->lastenforced = conshdlrdata->enforound;  ??? */
         continue;
      }

      /* if aux-violation is small (below feastol) and we look only for strong cuts, then it's unlikely to give a strong cut, so skip it */
      if( !allowweakcuts && auxviol < SCIPfeastol(scip) )  /*lint !e777*/
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   skip enforce using nlhdlr <%s> for expr %p (%s) with tiny " \
                  "auxviolation %g under:%d over:%d\n", nlhdlr->name, (void*)expr, expr->exprhdlr->name, auxviol,
                  underestimate, overestimate); )

         /* TODO expr->lastenforced = conshdlrdata->enforound;  ??? */
         continue;
      }

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   enforce using nlhdlr <%s> for expr %p (%s) with auxviolation " \
               "%g origviolation %g under:%d over:%d weak:%d\n", nlhdlr->name, (void*)expr, expr->exprhdlr->name,
               auxviol, origviol, underestimate, overestimate, allowweakcuts); )

      /* if we want to overestimate and violation w.r.t. auxiliary variables is also present on this side and nlhdlr
       * wants to be called for separation on this side, then call separation of nlhdlr
       */
      if( overestimate && auxoverestimate && (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPAABOVE) != 0 )  /*lint !e777*/
      {
         /* call the separation or estimation callback of the nonlinear handler for overestimation */
         hdlrresult = SCIP_DIDNOTFIND;
         SCIP_CALL( enforceExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, sol,
            expr->enfos[e]->auxvalue, TRUE, *result == SCIP_SEPARATED, allowweakcuts, inenforcement, &hdlrresult) );

         if( hdlrresult == SCIP_CUTOFF )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    found a cutoff -> stop separation\n"); )
            *result = SCIP_CUTOFF;
            expr->lastenforced = conshdlrdata->enforound;
            break;
         }

         if( hdlrresult == SCIP_SEPARATED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by cut\n", nlhdlr->name); )
            *result = SCIP_SEPARATED;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we give other nlhdlr another chance? (also #3070) */
            break;
         }

         if( hdlrresult == SCIP_REDUCEDDOM )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by boundchange\n", nlhdlr->name); )
            *result = SCIP_REDUCEDDOM;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_BRANCHED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> added branching candidate\n", nlhdlr->name); )
            assert(inenforcement);

            /* separation and domain reduction takes precedence over branching */
            assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED);
            if( *result == SCIP_DIDNOTFIND )
               *result = SCIP_BRANCHED;
            expr->lastenforced = conshdlrdata->enforound;
         }
      }

      /* if we want to underestimate and violation w.r.t. auxiliary variables is also present on this side and nlhdlr
       * wants to be called for separation on this side, then call separation of nlhdlr
       */
      if( underestimate && auxunderestimate && (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPABELOW) != 0 )  /*lint !e777*/
      {
         /* call the separation or estimation callback of the nonlinear handler for underestimation */
         hdlrresult = SCIP_DIDNOTFIND;
         SCIP_CALL( enforceExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, sol,
            expr->enfos[e]->auxvalue, FALSE, *result == SCIP_SEPARATED, allowweakcuts, inenforcement, &hdlrresult) );

         if( hdlrresult == SCIP_CUTOFF )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    found a cutoff -> stop separation\n"); )
            *result = SCIP_CUTOFF;
            expr->lastenforced = conshdlrdata->enforound;
            break;
         }

         if( hdlrresult == SCIP_SEPARATED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by cut\n", nlhdlr->name); )
            *result = SCIP_SEPARATED;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we give other nlhdlr another chance? (also #3070) */
            break;
         }

         if( hdlrresult == SCIP_REDUCEDDOM )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by boundchange\n", nlhdlr->name); )
            *result = SCIP_REDUCEDDOM;
            expr->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_BRANCHED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> added branching candidate\n", nlhdlr->name); )
            assert(inenforcement);

            /* separation takes precedence over branching */
            assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED);
            if( *result == SCIP_DIDNOTFIND )
               *result = SCIP_BRANCHED;
            expr->lastenforced = conshdlrdata->enforound;
         }
      }
   }

   return SCIP_OKAY;
}

/** helper function to enforce a single constraint */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_EXPRITER* it,               /**< expression iterator that we can just use here */
   SCIP_Bool             allowweakcuts,      /**< whether to allow weak cuts in this round */
   SCIP_Bool             inenforcement,      /**< whether to we are in enforcement, and not just separation */
   SCIP_RESULT*          result,             /**< pointer to update with result of the enforcing call */
   SCIP_Bool*            success             /**< buffer to store whether some enforcement took place */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* expr;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(it != NULL);
   assert(result != NULL);
   assert(success != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr->nenfos >= 0);

   *success = FALSE;

   if( inenforcement && !consdata->ispropagated )
   {
      /* If there are boundchanges that haven't been propagated to activities yet, then do this now and update bounds of
       * auxiliary variables, since some nlhdlr/exprhdlr may look at auxvar bounds or activities
       * (TODO: nlhdlr tells us now whether they do and so we could skip).
       * For now, update bounds of auxiliary variables only if called from enforcement, since updating auxvar bounds in
       * separation doesn't seem to be right (it would be ok if the boundchange cuts off the current LP solution by a
       * nice amount, but if not, we may just add a boundchange that doesn't change the dual bound much and could
       * confuse the stalling check for how long to do separation).
       */
      SCIP_Bool infeasible;
      int ntightenings;

      SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, inenforcement, &infeasible, &ntightenings) );
      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      /* if we tightened an auxvar bound, we better communicate that */
      if( ntightenings > 0 )
         *result = SCIP_REDUCEDDOM;
   }

   for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
   {
      SCIP_RESULT resultexpr;

      /* we can only enforce if there is an auxvar to compare with */
      if( expr->auxvar == NULL )
         continue;

      assert(expr->lastenforced <= conshdlrdata->enforound);
      if( expr->lastenforced == conshdlrdata->enforound )
      {
         ENFOLOG(
            SCIPinfoMessage(scip, enfologfile, "  skip expr ");
            SCIPprintExpr(scip, conshdlr, expr, enfologfile);
            SCIPinfoMessage(scip, enfologfile, " as already enforced in this enforound\n");
         )
         *success = TRUE;
         continue;
      }

      SCIP_CALL( enforceExpr(scip, conshdlr, cons, expr, sol, soltag, allowweakcuts, inenforcement, &resultexpr) );

      /* if not enforced, then we must not have found a cutoff, cut, domain reduction, or branchscore */
      assert((expr->lastenforced == conshdlrdata->enforound) == (resultexpr != SCIP_DIDNOTFIND));
      if( expr->lastenforced == conshdlrdata->enforound )
         *success = TRUE;

      if( resultexpr == SCIP_CUTOFF )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if( resultexpr == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;

      if( resultexpr == SCIP_REDUCEDDOM && *result != SCIP_SEPARATED )
         *result = SCIP_REDUCEDDOM;

      if( resultexpr == SCIP_BRANCHED && *result != SCIP_SEPARATED && *result != SCIP_REDUCEDDOM )
         *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** try to separate violated constraints and, if in enforcement, register branching scores
 *
 * Sets result to
 * - SCIP_DIDNOTFIND, if nothing of the below has been done
 * - SCIP_CUTOFF, if node can be cutoff,
 * - SCIP_SEPARATED, if a cut has been added,
 * - SCIP_REDUCEDDOM, if a domain reduction has been found,
 * - SCIP_BRANCHED, if branching has been done,
 * - SCIP_REDUCEDDOM, if a variable got fixed (in an attempt to branch on it),
 * - SCIP_INFEASIBLE, if external branching candidates were registered
 */
static
SCIP_RETCODE enforceConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, and not just separation */
   SCIP_Real             maxrelconsviol,     /**< largest scaled violation among all violated expr-constraints, only used if in enforcement */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRITER* it;
   SCIP_Bool consenforced;  /* whether any expression in constraint could be enforced */
   int c;

   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* increase tag to tell whether branching scores in expression belong to this sweep
    * and which expressions have already been enforced in this sweep
    * (we also want to distinguish sepa rounds, so this need to be here and not in consEnfo)
    */
   ++(conshdlrdata->enforound);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, TRUE) );

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      /* skip constraints that are not enabled or deleted */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      /* skip constraints that have separation disabled if we are only in separation */
      if( !inenforcement && !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;

      /* skip non-violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      ENFOLOG(
      {
         SCIP_CONSDATA* consdata;
         int i;
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         SCIPinfoMessage(scip, enfologfile, " constraint ");
         SCIP_CALL( SCIPprintCons(scip, conss[c], enfologfile) );
         SCIPinfoMessage(scip, enfologfile, "\n with viol %g and point\n", getConsAbsViolation(conss[c]));
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_VAR* var;
            var = SCIPexprvarGetVar(consdata->varexprs[i]);
            SCIPinfoMessage(scip, enfologfile, "  %-10s = %15g bounds: [%15g,%15g]\n", SCIPvarGetName(var),
                  SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
      })

      SCIP_CALL( enforceConstraint(scip, conshdlr, conss[c], sol, soltag, it, FALSE, inenforcement, result, &consenforced) );

      if( *result == SCIP_CUTOFF )
         break;

      if( !consenforced && inenforcement )
      {
         SCIP_Real viol;

         SCIP_CALL( getConsRelViolation(scip, conss[c], &viol, sol, soltag) );
         if( viol > conshdlrdata->weakcutminviolfactor * maxrelconsviol )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " constraint <%s> could not be enforced, try again with weak "\
                     "cuts allowed\n", SCIPconsGetName(conss[c])); )

            SCIP_CALL( enforceConstraint(scip, conshdlr, conss[c], sol, soltag, it, TRUE, inenforcement, result, &consenforced) );

            if( consenforced )
               ++conshdlrdata->nweaksepa;  /* TODO maybe this should not be counted per constraint, but per enforcement round? */

            if( *result == SCIP_CUTOFF )
               break;
         }
      }
   }

   SCIPfreeExpriter(&it);

   ENFOLOG( if( enfologfile != NULL ) fflush( enfologfile); )

   /* if having branching scores, then propagate them from expressions with children to variable expressions */
   if( *result == SCIP_BRANCHED )
   {
      /* having result set to branched here means only that we have branching candidates, we still need to do the actual
       * branching
       */
      SCIP_CALL( branching(scip, conshdlr, conss, nconss, maxrelconsviol, sol, soltag, result) );

      /* branching should either have branched: result == SCIP_BRANCHED,
       * or fixed a variable: result == SCIP_REDUCEDDOM,
       * or have registered external branching candidates: result == SCIP_INFEASIBLE,
       * or have not done anything: result == SCIP_DIDNOTFIND
       */
      assert(*result == SCIP_BRANCHED || *result == SCIP_REDUCEDDOM || *result == SCIP_INFEASIBLE || *result == SCIP_DIDNOTFIND);
   }

   ENFOLOG( if( enfologfile != NULL ) fflush( enfologfile); )

   return SCIP_OKAY;
}

/** collect (and print (if debugging enfo)) information on violation in expressions
 *
 * assumes that constraint violations have been computed
 */
static
SCIP_RETCODE analyzeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Real*            maxabsconsviol,     /**< buffer to store maximal absolute violation of constraints */
   SCIP_Real*            maxrelconsviol,     /**< buffer to store maximal relative violation of constraints */
   SCIP_Real*            minauxviol,         /**< buffer to store minimal (nonzero) violation of auxiliaries */
   SCIP_Real*            maxauxviol,         /**< buffer to store maximal violation of auxiliaries (violation in "extended formulation") */
   SCIP_Real*            maxvarboundviol     /**< buffer to store maximal violation of variable bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_Real v;
   int c;

   assert(conss != NULL || nconss == 0);
   assert(maxabsconsviol != NULL);
   assert(maxrelconsviol != NULL);
   assert(maxauxviol != NULL);
   assert(maxvarboundviol != NULL);

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   *maxabsconsviol = 0.0;
   *maxrelconsviol = 0.0;
   *minauxviol = SCIPinfinity(scip);
   *maxauxviol = 0.0;
   *maxvarboundviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* skip constraints that are not enabled, deleted, or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      v = getConsAbsViolation(conss[c]);
      *maxabsconsviol = MAX(*maxabsconsviol, v);

      /* skip non-violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      SCIP_CALL( getConsRelViolation(scip, conss[c], &v, sol, soltag) );
      *maxrelconsviol = MAX(*maxrelconsviol, v);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         SCIP_Real auxvarvalue;
         SCIP_Real auxvarlb;
         SCIP_Real auxvarub;
         SCIP_Bool violunder;
         SCIP_Bool violover;
         SCIP_Real origviol;
         SCIP_Real auxviol;
         int e;

         if( expr->auxvar == NULL )
         {
            /* check violation of variable bounds of original variable */
            if( SCIPisExprVar(expr) )
            {
               SCIP_VAR* var;
               var = SCIPexprvarGetVar(expr);
               auxvarvalue = SCIPgetSolVal(scip, sol, var);
               auxvarlb = SCIPvarGetLbLocal(var);
               auxvarub = SCIPvarGetUbLocal(var);

               origviol = 0.0;
               if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
                  origviol = auxvarlb - auxvarvalue;
               else if( auxvarub < auxvarvalue && !SCIPisInfinity(scip, auxvarub) )
                  origviol = auxvarvalue - auxvarub;
               if( origviol <= 0.0 )
                  continue;

               *maxvarboundviol = MAX(*maxvarboundviol, origviol);

               ENFOLOG(
               SCIPinfoMessage(scip, enfologfile, "var <%s>[%.15g,%.15g] = %.15g", SCIPvarGetName(var), auxvarlb, auxvarub, auxvarvalue);
               if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
                  SCIPinfoMessage(scip, enfologfile, " var >= lb violated by %g", auxvarlb - auxvarvalue);
               if( auxvarub < auxvarvalue && !SCIPisInfinity(scip,  auxvarub) )
                  SCIPinfoMessage(scip, enfologfile, " var <= ub violated by %g", auxvarvalue - auxvarub);
               SCIPinfoMessage(scip, enfologfile, "\n");
               )
            }

            continue;
         }

         auxvarvalue = SCIPgetSolVal(scip, sol, expr->auxvar);
         auxvarlb = SCIPvarGetLbLocal(expr->auxvar);
         auxvarub = SCIPvarGetUbLocal(expr->auxvar);

         /* check violation of variable bounds of auxiliary variable */
         if( auxvarlb - auxvarvalue > *maxvarboundviol && !SCIPisInfinity(scip, -auxvarlb) )
            *maxvarboundviol = auxvarlb - auxvarvalue;
         else if( auxvarvalue - auxvarub > *maxvarboundviol && !SCIPisInfinity(scip,  auxvarub) )
            *maxvarboundviol = auxvarvalue - auxvarub;

         origviol = getExprAbsOrigViolation(scip, expr, sol, &violunder, &violover);

         ENFOLOG(
         if( origviol > 0.0 || auxvarlb > auxvarvalue || auxvarub < auxvarvalue )
         {
            SCIPinfoMessage(scip, enfologfile, "expr ");
            SCIP_CALL( SCIPprintExpr(scip, conshdlr, expr, enfologfile) );
            SCIPinfoMessage(scip, enfologfile, " (%p)[%.15g,%.15g] = %.15g\n", (void*)expr, expr->activity.inf, expr->activity.sup, expr->evalvalue);

            SCIPinfoMessage(scip, enfologfile, "  auxvar <%s>[%.15g,%.15g] = %.15g", SCIPvarGetName(expr->auxvar), auxvarlb, auxvarub, auxvarvalue);
            if( origviol > 0.0 )
               SCIPinfoMessage(scip, enfologfile, " auxvar %s expr violated by %g", violunder ? ">=" : "<=", origviol);
            if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
               SCIPinfoMessage(scip, enfologfile, " auxvar >= auxvar's lb violated by %g", auxvarlb - auxvarvalue);
            if( auxvarub < auxvarvalue && !SCIPisInfinity(scip,  auxvarub) )
               SCIPinfoMessage(scip, enfologfile, " auxvar <= auxvar's ub violated by %g", auxvarvalue - auxvarub);
            SCIPinfoMessage(scip, enfologfile, "\n");
         }
         )

         /* no violation w.r.t. the original variables -> skip expression */
         if( origviol == 0.0 )
            continue;

         /* TODO remove? origviol shouldn't be mixed up with auxviol */
         *maxauxviol = MAX(*maxauxviol, origviol);  /*lint !e666*/
         *minauxviol = MIN(*minauxviol, origviol);  /*lint !e666*/

         /* compute aux-violation for each nonlinear handlers */
         for( e = 0; e < expr->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;

            /* eval in auxvars is only defined for nlhdrs that separate; there might not even be auxvars otherwise */
            if( (expr->enfos[e]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPABOTH) == 0 )
               continue;

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
            SCIP_CALL( SCIPcallNlhdlrEvalauxNonlinear(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &expr->enfos[e]->auxvalue, sol) );

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "  nlhdlr <%s> = %.15g", nlhdlr->name, expr->enfos[e]->auxvalue); )

            auxviol = getExprAbsAuxViolation(scip, expr, expr->enfos[e]->auxvalue, sol, &violunder, &violover);

            if( auxviol > 0.0 )  /*lint !e777*/
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " auxvar %s nlhdlr-expr violated by %g", violover ? "<=" : ">=", auxviol); )
               *maxauxviol = MAX(*maxauxviol, auxviol);
               *minauxviol = MIN(*minauxviol, auxviol);
            }
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "\n"); )
         }
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
} /*lint !e715*/

/** enforcement of constraints called by enfolp and enforelax */
static
SCIP_RETCODE consEnfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real maxabsconsviol;
   SCIP_Real maxrelconsviol;
   SCIP_Real minauxviol;
   SCIP_Real maxauxviol;
   SCIP_Real maxvarboundviol;
   unsigned int soltag;
   int nnotify;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   soltag = ++conshdlrdata->exprlastsoltag;

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_FEASIBLE )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: all expr-constraints feasible, skip enforcing\n",
               SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )
      return SCIP_OKAY;
   }

   SCIP_CALL( analyzeViolation(scip, conshdlr, conss, nconss, sol, soltag, &maxabsconsviol, &maxrelconsviol,
            &minauxviol, &maxauxviol, &maxvarboundviol) );

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: enforcing constraints with max conssviol=%e (rel=%e), "\
            "auxviolations in %g..%g, variable bounds violated by at most %g\n",
            SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), maxabsconsviol, maxrelconsviol, minauxviol, maxauxviol,
            maxvarboundviol); )

   assert(maxvarboundviol <= SCIPgetLPFeastol(scip));

   /* try to propagate */
   if( conshdlrdata->propinenforce )
   {
      SCIP_RESULT propresult;
      int nchgbds = 0;

      SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds) );

      if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
      {
         *result = propresult;
         return SCIP_OKAY;
      }
   }

   /* tighten the LP tolerance if violation in variables bounds is larger than aux-violation (max |expr - auxvar| over
    * all violated expr/auxvar in violated constraints)
    */
   if( conshdlrdata->tightenlpfeastol && maxvarboundviol > maxauxviol && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) &&
         sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxvarboundviol / 2.0, SCIPgetLPFeastol(scip) / 2.0)));  /*lint !e666*/
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " variable bound violation %g larger than auxiliary violation %g, "\
               "reducing LP feastol to %g\n", maxvarboundviol, maxauxviol, SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   SCIP_CALL( enforceConstraints(scip, conshdlr, conss, nconss, sol, soltag, TRUE, maxrelconsviol, result) );

   if( *result == SCIP_CUTOFF || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED ||
         *result == SCIP_INFEASIBLE )
      return SCIP_OKAY;

   assert(*result == SCIP_DIDNOTFIND);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " could not enforce violation %g in regular ways, LP feastol=%g, "\
            "becoming desperate now...\n", maxabsconsviol, SCIPgetLPFeastol(scip)); )

   if( conshdlrdata->tightenlpfeastol && SCIPisPositive(scip, maxvarboundviol) && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxvarboundviol / 2.0, SCIPgetLPFeastol(scip) / 2.0)));  /*lint !e666*/
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " variable bounds are violated by more than eps, reduced LP "\
               "feasibility tolerance to %g\n", SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   if( conshdlrdata->tightenlpfeastol && SCIPisPositive(scip, maxauxviol) && SCIPisPositive(scip,
            SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      /* try whether tighten the LP feasibility tolerance could help
       * maybe it is just some cut that hasn't been taken into account sufficiently
       * in the next enforcement round, we would then also allow even weaker cuts, as we want a minimal cut violation of LP's feastol
       * unfortunately, we do not know the current LP solution primal infeasibility, so sometimes this just repeats without effect
       * until the LP feastol reaches epsilon
       */
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxauxviol / 2.0, SCIPgetLPFeastol(scip) / 10.0)));  /*lint !e666*/
      ++conshdlrdata->ndesperatetightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " reduced LP feasibility tolerance to %g and hope\n", SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   /* try to propagate, if not tried above TODO(?) allow to disable this as well */
   if( !conshdlrdata->propinenforce )
   {
      SCIP_RESULT propresult;
      int nchgbds = 0;

      SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds) );

      if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
      {
         *result = propresult;
         return SCIP_OKAY;
      }
   }

   /* could not find branching candidates even when looking at minimal violated (>eps) expressions
    * now look if we find any unfixed variable that we could still branch on
    */
   SCIP_CALL( registerBranchingCandidatesAllUnfixed(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify > 0 )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " registered %d unfixed variables as branching candidates\n", nnotify); )
      ++conshdlrdata->ndesperatebranch;

      *result = SCIP_INFEASIBLE;  /* enforceConstraints may have changed it to SCIP_DIDNOTFIND */

      return SCIP_OKAY;
   }

   /* if everything is fixed in violated constraints, then let's cut off the node
    * - bound tightening with all vars fixed should prove cutoff, but interval arithmetic overestimates and so the
    *   result may not be conclusive (when constraint violations are small)
    * - if tightenlpfeastol=FALSE, then the LP solution that we try to enforce here may just not be within bounds
    *   sufficiently (see st_e40)
    * - but if the LP solution is really within bounds and since variables are fixed, cutting off the node is actually
    *   not "desperate", but a pretty obvious thing to do
    */
   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " enforcement with max. violation %g failed; cutting off node\n", maxabsconsviol); )
   *result = SCIP_CUTOFF;

   /* it's only "desperate" if the LP solution does not coincide with variable fixings (should we use something tighter than epsilon here?) */
   if( !SCIPisZero(scip, maxvarboundviol) )
      ++conshdlrdata->ndesperatecutoff;

   return SCIP_OKAY;
}

/** separation for all violated constraints to be used by SEPA callbacks */
static
SCIP_RETCODE consSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   unsigned int soltag;
   SCIP_Bool haveviol = FALSE;
   int c;

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   soltag = ++conshdlrdata->exprlastsoltag;

   /* compute violations */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);

      /* skip constraints that are not enabled, deleted, or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
         haveviol = TRUE;
   }

   /* if none of our constraints are violated, don't attempt separation */
   if( !haveviol )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: skip separation of non-violated constraints\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )
      return SCIP_OKAY;
   }

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: separation\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )

   /* call separation */
   SCIP_CALL( enforceConstraints(scip, conshdlr, conss, nconss, sol, soltag, FALSE, SCIP_INVALID, result) );

   return SCIP_OKAY;
}

/** Given a solution where every expression constraint is either feasible or can be made feasible by
 *  moving a linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
 *
 *  The method assumes that this is always possible and that not all constraints are feasible already.
 */
static
SCIP_RETCODE proposeFeasibleSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to process */
   SCIP_Bool*            success             /**< buffer to store whether we succeeded to construct a solution that satisfies all provided constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* newsol;
   int c;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(success != NULL);

   *success = FALSE;

   /* don't propose new solutions if not in presolve or solving */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );
   SCIPdebugMsg(scip, "attempt to make solution from <%s> feasible by shifting linear variable\n",
      sol != NULL ? (SCIPsolGetHeur(sol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(sol)) : "tree") : "LP");

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      SCIP_Real viol = 0.0;
      SCIP_Real delta;
      SCIP_Real gap;

      assert(consdata != NULL);

      /* get absolute violation and sign */
      if( consdata->lhsviol > SCIPfeastol(scip) )
         viol = consdata->lhsviol; /* lhs - activity */
      else if( consdata->rhsviol > SCIPfeastol(scip) )
         viol = -consdata->rhsviol; /* rhs - activity */
      else
         continue; /* constraint is satisfied */

      if( consdata->linvarincr != NULL &&
         ((viol > 0.0 && consdata->linvarincrcoef > 0.0) || (viol < 0.0 && consdata->linvarincrcoef < 0.0)) )
      {
         SCIP_VAR* var = consdata->linvarincr;

         /* compute how much we would like to increase var */
         delta = viol / consdata->linvarincrcoef;
         assert(delta > 0.0);

         /* if var has an upper bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            gap = SCIPvarGetUbGlobal(var) - SCIPgetSolVal(scip, newsol, var);
            delta = MIN(MAX(0.0, gap), delta);
         }
         if( SCIPisPositive(scip, delta) )
         {
            /* if variable is integral, round delta up so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPceil(scip, delta);

            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy lhs-violation %g of cons <%s>\n",
               SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));  /*lint !e613*/

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->linvarincrcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvardecr != NULL &&
         ((viol > 0.0 && consdata->linvardecrcoef < 0.0) || (viol < 0.0 && consdata->linvardecrcoef > 0.0)) )
      {
         SCIP_VAR* var = consdata->linvardecr;

         /* compute how much we would like to decrease var */
         delta = viol / consdata->linvardecrcoef;
         assert(delta < 0.0);

         /* if var has a lower bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
         {
            gap = SCIPgetSolVal(scip, newsol, var) - SCIPvarGetLbGlobal(var);
            delta = MAX(MIN(0.0, gap), delta);
         }
         if( SCIPisNegative(scip, delta) )
         {
            /* if variable is integral, round delta down so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPfloor(scip, delta);
            SCIP_CALL( SCIPincSolVal(scip, newsol, consdata->linvardecr, delta) );
            /*lint --e{613} */
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy rhs-violation %g of cons <%s>\n",
               SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->linvardecrcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so probably we could not make constraint feasible due to variable bounds, thus give up */
      break;
   }

   /* if we have a solution that should satisfy all quadratic constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic
    */
   if( c == nconss && (SCIPisInfinity(scip, SCIPgetUpperbound(scip)) || SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip))) )
   {
      SCIPdebugMsg(scip, "pass solution with objective val %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      assert(conshdlrdata->trysolheur != NULL);
      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );

      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** merges constraints that have the same root expression */
static
SCIP_RETCODE presolMergeConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            success             /**< pointer to store whether at least one constraint could be deleted */
   )
{
   SCIP_HASHMAP* expr2cons;
   SCIP_Bool* updatelocks;
   int* nlockspos;
   int* nlocksneg;
   int c;

   assert(success != NULL);

   *success = FALSE;

   /* not enough constraints available */
   if( nconss <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapCreate(&expr2cons, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &updatelocks, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlockspos, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksneg, nconss) );

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      /* ignore deleted constraints */
      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add expression to the hash map if not seen so far */
      if( !SCIPhashmapExists(expr2cons, (void*)consdata->expr) )
      {
         SCIP_CALL( SCIPhashmapInsertInt(expr2cons, (void*)consdata->expr, c) );
      }
      else
      {
         SCIP_CONSDATA* imgconsdata;
         int idx;

         idx = SCIPhashmapGetImageInt(expr2cons, (void*)consdata->expr);
         assert(idx >= 0 && idx < nconss);

         imgconsdata = SCIPconsGetData(conss[idx]);
         assert(imgconsdata != NULL);
         assert(imgconsdata->expr == consdata->expr);

         SCIPdebugMsg(scip, "merge constraint %g <= %s <= %g with %g <= %s <= %g\n", consdata->lhs,
            SCIPconsGetName(conss[c]), consdata->rhs, imgconsdata->lhs, SCIPconsGetName(conss[idx]), imgconsdata->rhs);

         /* check whether locks need to be updated */
         if( !updatelocks[idx] && ((SCIPisInfinity(scip, -imgconsdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs))
            || (SCIPisInfinity(scip, imgconsdata->rhs) && !SCIPisInfinity(scip, consdata->rhs))) )
         {
            nlockspos[idx] = imgconsdata->nlockspos;
            nlocksneg[idx] = imgconsdata->nlocksneg;
            SCIP_CALL( addLocks(scip, conss[idx], -imgconsdata->nlockspos, -imgconsdata->nlocksneg) );
            updatelocks[idx] = TRUE;
         }

         /* update constraint sides */
         imgconsdata->lhs = MAX(imgconsdata->lhs, consdata->lhs);
         imgconsdata->rhs = MIN(imgconsdata->rhs, consdata->rhs);

         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         *success = TRUE;
      }
   }

   /* restore locks of updated constraints */
   if( *success )
   {
      for( c = 0; c < nconss; ++c )
      {
         if( updatelocks[c] )
         {
            SCIP_CALL( addLocks(scip, conss[c], nlockspos[c], nlocksneg[c]) );
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &nlocksneg);
   SCIPfreeBufferArray(scip, &nlockspos);
   SCIPfreeBufferArray(scip, &updatelocks);
   SCIPhashmapFree(&expr2cons);

   return SCIP_OKAY;
}

/** print statistics for nonlinear handlers */
static
void printNlhdlrStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPinfoMessage(scip, file, "Nlhdlrs            : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "Detects", "EnfoCalls", "#IntEval", "PropCalls", "DetectAll", "Separated", "Cutoffs", "DomReds", "BranchScor", "DetectTime", "EnfoTime", "PropTime", "IntEvalTi");

   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_NLHDLR* nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      /* skip disabled nlhdlr */
      if( !nlhdlr->enabled )
         continue;

      SCIPinfoMessage(scip, file, "  %-17s:", nlhdlr->name);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndetectionslast);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nenfocalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nintevalcalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->npropcalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndetections);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nseparated);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ncutoffs);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndomreds);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nbranchscores);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->detecttime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->enfotime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->proptime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->intevaltime));
      SCIPinfoMessage(scip, file, "\n");
   }
}

/** print statistics for constraint handlers */
static
void printConshdlrStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPinfoMessage(scip, file, "Enforce            : %10s %10s %10s %10s %10s %10s\n", "WeakSepa", "TightenLP", "DespTghtLP", "DespBranch", "DespCutoff", "ForceLP");
   SCIPinfoMessage(scip, file, "  consexpr%-9s:", "");
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nweaksepa);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ntightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatetightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatebranch);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatecutoff);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nforcelp);
   SCIPinfoMessage(scip, file, "\n");
   SCIPinfoMessage(scip, file, "Presolve           : %10s\n", "CanonTime");
   SCIPinfoMessage(scip, file, "  consexpr%-9s:", "");
   SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, conshdlrdata->canonicalizetime));
   SCIPinfoMessage(scip, file, "\n");
}


/*
 * vertex polyhedral separation
 */

/** builds LP used to compute facets of the convex envelope of vertex-polyhedral functions */
static
SCIP_RETCODE buildVertexPolyhedralSeparationLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of (unfixed) variables in vertex-polyhedral functions */
   SCIP_LPI**            lp                  /**< pointer to store created LP */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   int* beg;
   int* ind;
   unsigned int nnonz;
   unsigned int ncols;
   unsigned int nrows;
   unsigned int i;
   unsigned int k;

   assert(scip != NULL);
   assert(lp != NULL);
   assert(nvars > 0);
   assert(nvars <= SCIP_MAXVERTEXPOLYDIM);

   SCIPdebugMsg(scip, "Building LP for computing facets of convex envelope of vertex-polyhedral function\n");

   /* create lpi to store the LP */
   SCIP_CALL( SCIPlpiCreate(lp, SCIPgetMessagehdlr(scip), "facet finding LP", SCIP_OBJSEN_MINIMIZE) );

   nrows = (unsigned int)nvars + 1;
   ncols = POWEROFTWO((unsigned int)nvars);
   nnonz = (ncols * (nrows + 1)) / 2;

   /* allocate necessary memory; set obj, lb, and ub to zero */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );

   /* calculate nonzero entries in the LP */
   for( i = 0, k = 0; i < ncols; ++i )
   {
      int row;
      unsigned int a;

      /* an upper bound of 1.0 is implied by the last row, but I presume that LP solvers prefer unbounded variables */
      ub[i] = SCIPlpiInfinity(*lp);

      SCIPdebugMsg(scip, "col %i starts at position %d\n", i, k);
      beg[i] = (int)k;
      row = 0;

      /* iterate through the bit representation of i */
      a = 1;
      while( a <= i )
      {
         if( (a & i) != 0 )
         {
            val[k] = 1.0;
            ind[k] = row;

            SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", row, i, k);

            ++k;
         }

         a <<= 1;
         ++row;
         assert(0 <= row && row <= SCIP_MAXVERTEXPOLYDIM);
         assert(POWEROFTWO(row) == a);
      }

      /* put 1 as a coefficient for sum_{i} \lambda_i = 1 row (last row) */
      val[k] = 1.0;
      ind[k] = (int)nrows - 1;
      ++k;
      SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", nrows - 1, i, k);
   }
   assert(k == nnonz);

   /* load all data into LP interface
    * we can assume nrows (=nvars+1) <= ncols (=2^nvars), so we can pass lb as dummy lhs and rhs
    */
   assert(nrows <= ncols);
   SCIP_CALL( SCIPlpiLoadColLP(*lp, SCIP_OBJSEN_MINIMIZE,
      (int)ncols, obj, lb, ub, NULL,
      (int)nrows, lb, lb, NULL,
      (int)nnonz, beg, ind, val) );

   /* for the last row, we can set the rhs to 1.0 already */
   ind[0] = (int)nrows - 1;
   val[0] = 1.0;
   SCIP_CALL( SCIPlpiChgSides(*lp, 1, ind, val, val) );

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** the given facet might not be a valid under(over)estimator, because of numerics and bad fixings; we compute \f$
 * \max_{v \in V} f(v) - (\alpha v + \beta) \f$ (\f$\max_{v \in V} \alpha v + \beta - f(v) \f$) where \f$ V \f$ is the
 * set of vertices of the domain
 */
static
SCIP_Real computeVertexPolyhedralMaxFacetError(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             overestimate,       /**< whether we check for an over or underestimator */
   SCIP_Real*            funvals,            /**< array containing the evaluation of the function at all corners, length: 2^nvars */
   SCIP_Real*            box,                /**< box for which facet was computed, length: 2*nallvars */
   int                   nallvars,           /**< number of all variables */
   int                   nvars,              /**< number of unfixed variables */
   int*                  nonfixedpos,        /**< indices of unfixed variables, length: nvars */
   SCIP_Real*            facetcoefs,         /**< current facet candidate's coefficients, length: nallvars */
   SCIP_Real             facetconstant       /**< current facet candidate's constant, length: nallvars */
   )
{
   SCIP_Real maxerror;
   SCIP_Real facetval;
   SCIP_Real funval;
   SCIP_Real error;
   unsigned int i;
   unsigned int ncorners;
   unsigned int prev;

   assert(scip != NULL);
   assert(funvals != NULL);
   assert(box != NULL);
   assert(nonfixedpos != NULL);
   assert(facetcoefs != NULL);

   ncorners = POWEROFTWO(nvars);
   maxerror = 0.0;

   /* check the origin (all variables at lower bound) */
   facetval = facetconstant;
   for( i = 0; i < (unsigned int) nallvars; ++i )
      facetval += facetcoefs[i] * box[2*i];

   /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
   funval = funvals[0];
   if( overestimate )
      error = funval - facetval;
   else
      error = facetval - funval;

   /* update maximum error */
   maxerror = MAX(error, maxerror);

   prev = 0;
   for( i = 1; i < ncorners; ++i )
   {
      unsigned int gray;
      unsigned int diff;
      unsigned int pos;
      int origpos;

      gray = i ^ (i >> 1);
      diff = gray ^ prev;

      /* compute position of unique 1 of diff */
      pos = 0;
      while( (diff >>= 1) != 0 )
         ++pos;
      assert(pos < (unsigned int)nvars);

      origpos = nonfixedpos[pos];

      if( gray > prev )
         facetval += facetcoefs[origpos] * (box[2*origpos+1] - box[2*origpos]);
      else
         facetval -= facetcoefs[origpos] * (box[2*origpos+1] - box[2*origpos]);

      /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
      funval = funvals[gray];
      if( overestimate )
         error = funval - facetval;
      else
         error = facetval - funval;

      /* update  maximum error */
      maxerror = MAX(error, maxerror);

      prev = gray;
   }

   SCIPdebugMsg(scip, "maximum error of facet: %2.8e\n", maxerror);

   return maxerror;
}

/** computes a facet of the convex or concave envelope of a vertex polyhedral function using by solving an LP */
static
SCIP_RETCODE computeVertexPolyhedralFacetLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   int*                  nonfixedpos,        /**< indices of nonfixed variables */
   SCIP_Real*            funvals,            /**< values of function in all corner points (w.r.t. nonfixed variables) */
   int                   nvars,              /**< number of nonfixed variables */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an zero'ed array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LPI* lp;
   SCIP_Real* aux; /* used to transform x^* and then to store LP solution */
   int* inds;
   int ncols;
   int nrows;
   int i;
   SCIP_Real facetvalue;
   SCIP_Real mindomwidth;
   SCIP_RETCODE lpsolveretcode;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(xstar != NULL);
   assert(box != NULL);
   assert(nonfixedpos != NULL);
   assert(funvals != NULL);
   assert(nvars <= SCIP_MAXVERTEXPOLYDIM);
   assert(success != NULL);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->vp_randnumgen == NULL && conshdlrdata->vp_maxperturb > 0.0 )
   {
      SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->vp_randnumgen, VERTEXPOLY_RANDNUMINITSEED, TRUE) );
   }

   /* construct an LP for this size, if not having one already */
   if( conshdlrdata->vp_lp[nvars] == NULL )
   {
      SCIP_CALL( buildVertexPolyhedralSeparationLP(scip, nvars, &conshdlrdata->vp_lp[nvars]) );
   }
   lp = conshdlrdata->vp_lp[nvars];
   assert(lp != NULL);

   /* get number of cols and rows of separation lp */
   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );

   /* number of columns should equal the number of corners = 2^nvars */
   assert(ncols == (int)POWEROFTWO(nvars));

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &aux, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );

   /*
    * set up the described LP on the transformed space
    */

   for( i = 0; i < ncols; ++i )
      inds[i] = i;

   /* compute T^-1(x^*), i.e. T^-1(x^*)_i = (x^*_i - lb_i)/(ub_i - lb_i) */
   mindomwidth = 2*SCIPinfinity(scip);
   for( i = 0; i < nrows-1; ++i )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;
      int varpos;

      assert(i < nvars);

      varpos = nonfixedpos[i];
      lb = box[2 * varpos];
      ub = box[2 * varpos + 1];
      solval = xstar[varpos];

      if( ub - lb < mindomwidth )
         mindomwidth = ub - lb;

      /* explicitly handle solution which violate bounds of variables (this can happen because of tolerances) */
      if( solval <= lb )
         aux[i] = 0.0;
      else if( solval >= ub )
         aux[i] = 1.0;
      else
         aux[i] = (solval - lb) / (ub - lb);

      /* perturb point to hopefully obtain a facet of the convex envelope */
      if( conshdlrdata->vp_maxperturb > 0.0 )
      {
         assert(conshdlrdata->vp_randnumgen != NULL);

         if( aux[i] == 1.0 )
            aux[i] -= SCIPrandomGetReal(conshdlrdata->vp_randnumgen, 0.0, conshdlrdata->vp_maxperturb);
         else if( aux[i] == 0.0 )
            aux[i] += SCIPrandomGetReal(conshdlrdata->vp_randnumgen, 0.0, conshdlrdata->vp_maxperturb);
         else
         {
            SCIP_Real perturbation;

            perturbation = MIN( aux[i], 1.0 - aux[i] ) / 2.0;
            perturbation = MIN( perturbation, conshdlrdata->vp_maxperturb );
            aux[i] += SCIPrandomGetReal(conshdlrdata->vp_randnumgen, -perturbation, perturbation);
         }
         assert(0.0 < aux[i] && aux[i] < 1.0);
      }

      SCIPdebugMsg(scip, "LP row %d in [%e, %e]\n", i, aux[i], aux[i]);
   }

   /* update LP */
   SCIP_CALL( SCIPlpiChgObj(lp, ncols, inds, funvals) );
   SCIP_CALL( SCIPlpiChgSides(lp, nrows-1, inds, aux, aux) );
   SCIP_CALL( SCIPlpiChgObjsen(lp, overestimate ? SCIP_OBJSEN_MAXIMIZE : SCIP_OBJSEN_MINIMIZE) );

   /* we can stop the LP solve if will not meet the target value anyway, but only if xstar hasn't been perturbed */
   if( conshdlrdata->vp_maxperturb == 0.0 && !SCIPisInfinity(scip, REALABS(targetvalue)) )
   {
      SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_OBJLIM, targetvalue) );
   }
   /* set an iteration limit so we do not run forever */
   SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_LPITLIM, 100*ncols) );
   /* since we work with the dual of the LP, primal feastol determines how much we want the computed facet to be the best possible one */
   SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_FEASTOL, SCIPfeastol(scip)) );
   /* since we work with the dual of the LP, dual feastol determines validity of the facet
    * if some ub-lb is small, we need higher accuracy, since below we divide coefs by ub-lb (we moved and scaled the box)
    * thus, we set the dual feastol to be between SCIPepsilon and SCIPfeastol
    */
   SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, MIN(SCIPfeastol(scip), MAX(SCIPepsilon(scip), mindomwidth * SCIPfeastol(scip)))) ); /*lint !e666*/

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_LPINFO, 1) );
#endif
   /* SCIP_CALL( SCIPlpiWriteLP(lp, "lp.lp") ); */

   /*
    * solve the LP and store the resulting facet for the transformed space
    */
   if( conshdlrdata->vp_dualsimplex )
   {
      lpsolveretcode = SCIPlpiSolveDual(lp);
   }
   else
   {
      lpsolveretcode = SCIPlpiSolvePrimal(lp);
   }
   if( lpsolveretcode == SCIP_LPERROR )
   {
      SCIPdebugMsg(scip, "LP error, aborting.\n");
      goto CLEANUP;
   }
   SCIP_CALL( lpsolveretcode );

   /* any dual feasible solution should provide a valid estimator (and a dual optimal one a facet) */
   if( !SCIPlpiIsDualFeasible(lp) )
   {
      SCIPdebugMsg(scip, "LP not solved to dual feasibility, aborting.\n");
      goto CLEANUP;
   }

   /* get dual solution (facet of convex envelope); again, we have to be careful since the LP can have more rows and
    * columns than needed, in particular, \bar \beta is the last dual multiplier
    */
   SCIP_CALL( SCIPlpiGetSol(lp, NULL, NULL, aux, NULL, NULL) );

   for( i = 0; i < nvars; ++i )
      facetcoefs[nonfixedpos[i]] = aux[i];
   /* last dual multiplier is the constant */
   *facetconstant = aux[nrows - 1];


#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "facet for the transformed problem: ");
   for( i = 0; i < nallvars; ++i )
   {
      SCIPdebugMsgPrint(scip, "%3.4e * x%d + ", facetcoefs[i], i);
   }
   SCIPdebugMsgPrint(scip, "%3.4e\n", *facetconstant);
#endif

   /*
    * transform the facet to original space and compute value at x^*, i.e., alpha x + beta
    */

   SCIPdebugMsg(scip, "facet in orig. space: ");

   facetvalue = 0.0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      int varpos;

      varpos = nonfixedpos[i];
      lb = box[2 * varpos];
      ub = box[2 * varpos + 1];
      assert(!SCIPisEQ(scip, lb, ub));

      /* alpha_i := alpha_bar_i / (ub_i - lb_i) */
      facetcoefs[varpos] = facetcoefs[varpos] / (ub - lb);

      /* beta = beta_bar - sum_i alpha_i * lb_i */
      *facetconstant -= facetcoefs[varpos] * lb;

      /* evaluate */
      facetvalue += facetcoefs[varpos] * xstar[varpos];

      SCIPdebugMsgPrint(scip, "%3.4e * x%d + ", facetcoefs[varpos], varpos);
   }
   SCIPdebugMsgPrint(scip, "%3.4e ", *facetconstant);

   /* add beta to the facetvalue: at this point in the code, facetvalue = g(x^*) */
   facetvalue += *facetconstant;

   SCIPdebugMsgPrint(scip, "has value %g, target = %g\n", facetvalue, targetvalue);

    /* if overestimate, then we want facetvalue < targetvalue
    * if underestimate, then we want facetvalue > targetvalue
    * if none holds, give up
    * so maybe here we should check against the minimal violation
    */
   if( overestimate == (facetvalue > targetvalue) )
   {
      SCIPdebugMsg(scip, "missed the target, facetvalue %g targetvalue %g, overestimate=%d\n", facetvalue, targetvalue, overestimate);
      goto CLEANUP;
   }

   /* if we made it until here, then we have a nice facet */
   *success = TRUE;

CLEANUP:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &aux);

   return SCIP_OKAY;
}

/** computes a facet of the convex or concave envelope of a univariant vertex polyhedral function
 *
 * In other words, compute the line that passes through two given points.
 */
static
SCIP_RETCODE computeVertexPolyhedralFacetUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             left,               /**< left coordinate */
   SCIP_Real             right,              /**< right coordinate */
   SCIP_Real             funleft,            /**< value of function in left coordinate */
   SCIP_Real             funright,           /**< value of function in right coordinate */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoef,          /**< buffer to store coefficient of facet defining inequality */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{
   assert(scip != NULL);
   assert(SCIPisLE(scip, left, right));
   assert(!SCIPisInfinity(scip, -left));
   assert(!SCIPisInfinity(scip, right));
   assert(SCIPisFinite(funleft) && funleft != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(funright) && funright != SCIP_INVALID);  /*lint !e777*/
   assert(success != NULL);
   assert(facetcoef != NULL);
   assert(facetconstant != NULL);

   *facetcoef = (funright - funleft) / (right - left);
   *facetconstant = funleft - *facetcoef * left;

   *success = TRUE;

   return SCIP_OKAY;
}

/** computes a facet of the convex or concave envelope of a bivariate vertex polyhedral function */
static
SCIP_RETCODE computeVertexPolyhedralFacetBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real             p1[2],              /**< first vertex of box */
   SCIP_Real             p2[2],              /**< second vertex of box */
   SCIP_Real             p3[2],              /**< third vertex of box */
   SCIP_Real             p4[2],              /**< forth vertex of box */
   SCIP_Real             p1val,              /**< value in p1 */
   SCIP_Real             p2val,              /**< value in p2 */
   SCIP_Real             p3val,              /**< value in p3 */
   SCIP_Real             p4val,              /**< value in p4 */
   SCIP_Real             xstar[2],           /**< point to be separated */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least 2 */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{
   SCIP_Real alpha, beta, gamma_, delta;
   SCIP_Real xstarval, candxstarval = 0.0;
   int leaveout;

   assert(scip != NULL);
   assert(success != NULL);
   assert(SCIPisFinite(p1val) && p1val != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(p2val) && p2val != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(p3val) && p3val != SCIP_INVALID);  /*lint !e777*/
   assert(SCIPisFinite(p4val) && p4val != SCIP_INVALID);  /*lint !e777*/
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   /* if we want an underestimator, flip f(x,y), i.e., do as if we compute an overestimator for -f(x,y) */
   if( !overestimate )
   {
      p1val = -p1val;
      p2val = -p2val;
      p3val = -p3val;
      p4val = -p4val;
      targetvalue = -targetvalue;
   }

   SCIPdebugMsg(scip, "p1 = (%g, %g), f(p1) = %g\n", p1[0], p1[1], p1val);
   SCIPdebugMsg(scip, "p2 = (%g, %g), f(p2) = %g\n", p2[0], p2[1], p2val);
   SCIPdebugMsg(scip, "p3 = (%g, %g), f(p3) = %g\n", p3[0], p3[1], p3val);
   SCIPdebugMsg(scip, "p4 = (%g, %g), f(p4) = %g\n", p4[0], p4[1], p4val);

   /* Compute coefficients alpha, beta, gamma (>0), delta such that
    *   alpha*x + beta*y + gamma*z = delta
    * is satisfied by at least three of the corner points (p1,f(p1)), ..., (p4,f(p4)) and
    * the fourth corner point lies below this hyperplane.
    * Since we assume that f is vertex-polyhedral, we then know that all points (x,y,f(x,y)) are below this hyperplane, i.e.,
    *    alpha*x + beta*y - delta <= -gamma * f(x,y),
    * or, equivalently,
    *   -alpha/gamma*x - beta/gamma*y + delta/gamma >= f(x,y).
    */
   for( leaveout = 1; leaveout <= 4; ++leaveout )
   {
      switch( leaveout)
      {
         case 1 :
            /* get hyperplane through p2, p3, p4 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p1, then go to next candidate */
            if( alpha * p1[0] + beta * p1[1] + gamma_ * p1val - delta > 0.0 )
               continue;
            break;

         case 2 :
            /* get hyperplane through p1, p3, p4 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p2, then go to next candidate */
            if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val - delta > 0.0 )
               continue;
            break;

         case 3 :
            /* get hyperplane through p1, p2, p4 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p3, then go to next candidate */
            if( alpha * p3[0] + beta * p3[1] + gamma_ * p3val - delta > 0.0 )
               continue;
            break;

         case 4 :
            /* get hyperplane through p1, p2, p3 */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p4, then stop */
            if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val - delta > 0.0 )
               continue;
            break;

         default: /* only for lint */
            alpha = SCIP_INVALID;
            beta = SCIP_INVALID;
            gamma_ =  SCIP_INVALID;
            delta = SCIP_INVALID;
            break;
      }

      /* check if bad luck: should not happen if numerics are fine */
      if( SCIPisZero(scip, gamma_) )
         continue;
      assert(!SCIPisNegative(scip, gamma_));

      /* if coefficients become tiny because division by gamma makes them < SCIPepsilon(scip), then skip, too */
      if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
         ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta/gamma_)) )
         continue;

      SCIPdebugMsg(scip, "alpha = %g, beta = %g, gamma = %g, delta = %g\n", alpha, beta, gamma_, delta);

      /* value of hyperplane candidate in xstar */
      xstarval = -alpha/gamma_ * xstar[0] -beta/gamma_ * xstar[1] + delta/gamma_;

      /* if reaching target and first or better than previous candidate, then update */
      if( xstarval <= targetvalue && (!*success || xstarval < candxstarval) )
      {
         /* flip hyperplane */
         if( !overestimate )
            gamma_ = -gamma_;

         facetcoefs[0] = -alpha/gamma_;
         facetcoefs[1] = -beta/gamma_;
         *facetconstant = delta/gamma_;

         *success = TRUE;
         candxstarval = xstarval;
      }
   }

   return SCIP_OKAY;
}

/** hash key retrieval function for bilinear term entries */
static
SCIP_DECL_HASHGETKEY(bilinearTermsGetHashkey)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int idx;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userptr;
   assert(conshdlrdata != NULL);

   idx = ((int)(size_t)elem) - 1;
   assert(idx >= 0 && idx < conshdlrdata->nbilinterms);

   return (void*)&conshdlrdata->bilinterms[idx];
}

/** returns TRUE iff the bilinear term entries are equal */
static
SCIP_DECL_HASHKEYEQ(bilinearTermsIsHashkeyEq)
{  /*lint --e{715}*/
   SCIP_CONSNONLINEAR_BILINTERM* entry1;
   SCIP_CONSNONLINEAR_BILINTERM* entry2;

   /* get corresponding entries */
   entry1 = (SCIP_CONSNONLINEAR_BILINTERM*)key1;
   entry2 = (SCIP_CONSNONLINEAR_BILINTERM*)key2;
   assert(entry1->x != NULL && entry1->y != NULL);
   assert(entry2->x != NULL && entry2->y != NULL);
   assert(SCIPvarCompare(entry1->x, entry1->y) < 1);
   assert(SCIPvarCompare(entry2->x, entry2->y) < 1);

   return entry1->x == entry2->x && entry1->y == entry2->y;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(bilinearTermsGetHashkeyVal)
{  /*lint --e{715}*/
   SCIP_CONSNONLINEAR_BILINTERM* entry;

   entry = (SCIP_CONSNONLINEAR_BILINTERM*)key;
   assert(entry->x != NULL && entry->y != NULL);
   assert(SCIPvarCompare(entry->x, entry->y) < 1);

   return SCIPhashTwo(SCIPvarGetIndex(entry->x), SCIPvarGetIndex(entry->y));
}

/** ensures the auxexprs array in the bilinear term is large enough to store n entries */
static
SCIP_RETCODE ensureBilinTermSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSNONLINEAR_BILINTERM* term,            /**< bilinear term */
   int                   n                   /**< number of entries to store */
   )
{
   int newsize;

   assert(term != NULL);

   /* check whether array is large enough */
   if( n <= term->sauxexprs )
      return SCIP_OKAY;

   /* compute new size */
   newsize = SCIPcalcMemGrowSize(scip, n);
   assert(n <= newsize);

   /* realloc array */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &term->aux.exprs, term->sauxexprs, newsize) );
   term->sauxexprs = newsize;

   return SCIP_OKAY;
}

/** compare two auxiliary expressions
 *
 *  Compares auxiliary variables, followed by coefficients and then constants.
 */
static
SCIP_DECL_SORTPTRCOMP(auxexprComp)
{
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr1 = (SCIP_CONSNONLINEAR_AUXEXPR*)elem1;
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr2 = (SCIP_CONSNONLINEAR_AUXEXPR*)elem2;
   int compvars;
   int i;

   /* compare the auxiliary variables */
   compvars = SCIPvarCompare(auxexpr1->auxvar, auxexpr2->auxvar); /* TODO can one of these be NULL? */

   if( compvars != 0 )
      return compvars;

   /* compare the coefficients and constants */
   for( i = 0; i < 3; ++i )
   {
      if( auxexpr1->coefs[i] != auxexpr2->coefs[i] ) /*lint !e777*/
         return auxexpr1->coefs[i] < auxexpr2->coefs[i] ? -1 : 1;
   }

   return auxexpr1->cst < auxexpr2->cst ? -1 : auxexpr1->cst == auxexpr2->cst ? 0 : 1; /*lint !e777*/
}

/* add an auxiliary expression to a bilinear term */
static
SCIP_RETCODE bilinTermAddAuxExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< expr constraint handler data */
   SCIP_CONSNONLINEAR_BILINTERM* term,            /**< bilinear term */
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr,           /**< auxiliary expression to add */
   SCIP_Bool*            added               /**< pointer to store whether auxexpr has been added */
   )
{
   SCIP_Bool found;
   int pos;
   int i;

   *added = FALSE;

   /* check if auxexpr has already been added to term */
   if( term->nauxexprs == 0 )
   {
      found = FALSE;
      pos = 0;
   }
   else
   {
      found = SCIPsortedvecFindPtr((void**)term->aux.exprs, auxexprComp, auxexpr, term->nauxexprs, &pos);
   }

   if( !found )
   {
      if( term->nauxexprs >= conshdlrdata->bilinmaxnauxexprs )
         return SCIP_OKAY;

      SCIP_CALL( ensureBilinTermSize(scip, term, term->nauxexprs + 1) );
      assert(term->sauxexprs >= term->nauxexprs + 1);

      /* insert expression at the correct position */
      for( i = term->nauxexprs; i > pos; --i )
      {
         term->aux.exprs[i] = term->aux.exprs[i-1];
      }
      term->aux.exprs[pos] = auxexpr;
      ++(term->nauxexprs);
      *added = TRUE;
   }
   else
   {
      term->aux.exprs[pos]->underestimate += auxexpr->underestimate;
      term->aux.exprs[pos]->overestimate += auxexpr->overestimate;
   }

   return SCIP_OKAY;
}

/** resizes an array of bilinear terms */
static
SCIP_RETCODE bilinearTermsResize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int                   reqsize             /**< required size */
   )
{
   int newsize;

   assert(conshdlrdata != NULL);

   /* check whether array is large enough */
   if( reqsize <= conshdlrdata->bilintermssize )
      return SCIP_OKAY;

   /* compute new size */
   newsize = SCIPcalcMemGrowSize(scip, reqsize);
   assert(reqsize <= newsize);

   /* realloc array */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->bilinterms, conshdlrdata->bilintermssize,
      newsize) );
   conshdlrdata->bilintermssize = newsize;

   return SCIP_OKAY;
}

/** iterates through all expressions of all expression constraints and adds the corresponding bilinear terms to the
 *  hash table
 */
static
SCIP_RETCODE bilinearTermsInsertAll(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< expression constraints */
   int                   nconss              /**< total number of expression constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRITER* it;
   SCIP_EXPRHDLR* producthdlr;
   SCIP_EXPRHDLR* powhdlr;
   int c;

   assert(conss != NULL || nconss == 0);

   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether the bilinear terms have been stored already */
   if( conshdlrdata->bilinterms != NULL )
      return SCIP_OKAY;

   /* create and initialize iterator */
   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

   /* get product and pow expression handlers */
   producthdlr = SCIPgetExprHdlrProduct(conshdlr);
   powhdlr = SCIPgetExprHdlrPower(conshdlr);

   /* iterate through all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_EXPR* expr;

      assert(conss != NULL && conss[c] != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* iterate through all expressions */
      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         SCIP_EXPR** children = SCIPexprGetChildren(expr);
         SCIP_VAR* x = NULL;
         SCIP_VAR* y = NULL;

         /* check whether the expression is of the form f(..)^2 */
         if( SCIPexprGetHdlr(expr) == powhdlr && SCIPexprpowGetExponent(expr) == 2.0 )
         {
            x = SCIPgetExprAuxVarNonlinear(children[0]);
            y = x;
         }
         /* check whether the expression is of the form f(..) * g(..) */
         else if( SCIPexprGetHdlr(expr) == producthdlr && SCIPexprGetNChildren(expr) == 2 )
         {
            x = SCIPgetExprAuxVarNonlinear(children[0]);
            y = SCIPgetExprAuxVarNonlinear(children[1]);
         }

         /* add variables to the hash table */
         if( x != NULL && y != NULL )
         {
            SCIP_CALL( SCIPinsertBilinearTermExistingNonlinear(scip, conshdlr, x, y, SCIPgetExprAuxVarNonlinear(expr),
               SCIPgetExprNLocksPosNonlinear(expr), SCIPgetExprNLocksNegNonlinear(expr)) );
         }
      }
   }

   /* release iterator */
   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** store x, y and the locks in a new bilinear term */
static
SCIP_RETCODE bilinearTermsInsertEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expr constraint handler */
   SCIP_VAR*             x,                  /**< the first variable */
   SCIP_VAR*             y,                  /**< the second variable */
   int                   nlockspos,          /**< number of positive locks of the bilinear term */
   int                   nlocksneg,          /**< number of negative locks of the bilinear term */
   int*                  idx,                /**< pointer to store the position of the term in bilinterms array */
   SCIP_Bool             existing            /**< whether the term exists explicitly in the problem */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM* term;

   assert(conshdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(nlockspos >= 0);
   assert(nlocksneg >= 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *idx = SCIPgetBilinTermIdxNonlinear(conshdlr, x, y);

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
   }
   assert(SCIPvarCompare(x, y) < 1);

   /* update or create the term */
   if( *idx >= 0 )
   { /* the term has already been added */
      assert(conshdlrdata->bilinterms[*idx].x == x);
      assert(conshdlrdata->bilinterms[*idx].y == y);

      /* get term and add locks */
      term = &conshdlrdata->bilinterms[*idx];
      assert(existing <= term->existing); /* implicit terms are added after existing ones */
      term->nlockspos += nlockspos;
      term->nlocksneg += nlocksneg;
   }
   else
   { /* this is the first time we encounter this product */
      /* ensure size of bilinterms array */
      SCIP_CALL( bilinearTermsResize(scip, conshdlrdata, conshdlrdata->nbilinterms + 1) );

      *idx = conshdlrdata->nbilinterms;

      /* get term and set values in the created bilinear term */
      term = &conshdlrdata->bilinterms[*idx];
      assert(term != NULL);
      term->x = x;
      term->y = y;
      term->nauxexprs = 0;
      term->sauxexprs = 0;
      term->nlockspos = nlockspos;
      term->nlocksneg = nlocksneg;
      term->existing = existing;
      if( existing )
         term->aux.var = NULL;
      else
         term->aux.exprs = NULL;

      /* increase the total number of bilinear terms */
      ++(conshdlrdata->nbilinterms);

      /* save to the hashtable */
      if( conshdlrdata->bilinhashtable == NULL )
      {
         SCIP_CALL( SCIPhashtableCreate(&conshdlrdata->bilinhashtable, SCIPblkmem(scip), conshdlrdata->nbilinterms,
               bilinearTermsGetHashkey, bilinearTermsIsHashkeyEq, bilinearTermsGetHashkeyVal,
               (void*)conshdlrdata) );
      }
      assert(conshdlrdata->bilinhashtable != NULL);

      /* insert the index of the bilinear term into the hash table; note that the index of the i-th element is (i+1)
       * because zero can not be inserted into hash table
       */
      SCIP_CALL( SCIPhashtableInsert(conshdlrdata->bilinhashtable, (void*)(size_t)(*idx + 1)) );/*lint !e571 !e776*/

      /* capture product variables */
      SCIP_CALL( SCIPcaptureVar(scip, x) );
      SCIP_CALL( SCIPcaptureVar(scip, y) );
   }

   return SCIP_OKAY;
}

/** frees array of bilinear terms and hash table */
static
SCIP_RETCODE bilinearTermsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int i;
   int j;

   assert(conshdlrdata != NULL);

   /* check whether bilinear terms have been stored */
   if( conshdlrdata->bilinterms == NULL )
   {
      assert(conshdlrdata->bilinterms == NULL);
      assert(conshdlrdata->nbilinterms == 0);
      assert(conshdlrdata->bilintermssize == 0);

      return SCIP_OKAY;
   }

   /* release variables */
   for( i = 0; i < conshdlrdata->nbilinterms; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].y) );
      SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].x) );

      for( j = 0; j < conshdlrdata->bilinterms[i].nauxexprs; ++j )
      {
         if( conshdlrdata->bilinterms[i].aux.exprs[j]->auxvar != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].aux.exprs[j]->auxvar) );
         }
         SCIPfreeBlockMemory(scip, &(conshdlrdata->bilinterms[i].aux.exprs[j])); /*lint !e866*/
      }

      if( conshdlrdata->bilinterms[i].nauxexprs > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->bilinterms[i].aux.exprs), conshdlrdata->bilinterms[i].sauxexprs);
         continue;
      }

      /* the rest is for simple terms with a single auxvar */

      /* it might be that there is a bilinear term without a corresponding auxiliary variable */
      if( conshdlrdata->bilinterms[i].aux.var != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].aux.var) );
      }
   }

   /* free hash table */
   if( conshdlrdata->bilinhashtable != NULL )
   {
      SCIPhashtableFree(&conshdlrdata->bilinhashtable);
   }

   /* free bilinterms array; reset counters */
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bilinterms, conshdlrdata->bilintermssize);
   conshdlrdata->nbilinterms = 0;
   conshdlrdata->bilintermssize = 0;

   return SCIP_OKAY;
}

/** returns whether the variable of a given variable expression is a candidate for presolSingleLockedVars(), i.e.,
 *  the variable is only contained in a single expression constraint, has no objective coefficient, has finite
 *  variable bounds, and is not binary
 */
static
SCIP_Bool isSingleLockedCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_VAR* var;

   assert(SCIPisExprVar(expr));

   var = SCIPexprvarGetVar(expr);
   assert(var != NULL);

   return SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) == SCIPgetExprNLocksNegNonlinear(expr)
      && SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) == SCIPgetExprNLocksPosNonlinear(expr)
      && SCIPgetConsExprExprVarNConss(expr) == 1 && SCIPisZero(scip, SCIPvarGetObj(var))
      && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var))
      && SCIPvarGetType(var) != SCIP_VARTYPE_BINARY
      && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
}

/** removes all variable expressions that are contained in a given expression from a hash map */
static
SCIP_RETCODE removeSingleLockedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_EXPRITER* it,               /**< expression iterator */
   SCIP_HASHMAP*         exprcands           /**< map to hash variable expressions */
   )
{
   SCIP_EXPR* e;

   for( e = SCIPexpriterRestartDFS(it, expr); !SCIPexpriterIsEnd(it); e = SCIPexpriterGetNext(it) ) /*lint !e441*/
   {
      if( SCIPisExprVar(e) && SCIPhashmapExists(exprcands, (void*)e) )
      {
         SCIP_CALL( SCIPhashmapRemove(exprcands, (void*)e) );
      }
   }

   return SCIP_OKAY;
}

/* presolving method to check if there is a single linear continuous variable that can be made implicit integer */
static
SCIP_RETCODE presolveImplint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< expression constraints */
   int                   nconss,             /**< total number of expression constraints */
   int*                  nchgvartypes,       /**< pointer to update the total number of changed variable types */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   SCIP_EXPRHDLR* sumhdlr;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nchgvartypes != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* nothing can be done if there are no integer variables available */
   if( SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   /* get sum expression handler */
   sumhdlr = SCIPgetExprHdlrSum(conshdlr);

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_EXPR** children;
      int nchildren;

      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      children = SCIPexprGetChildren(consdata->expr);
      nchildren = SCIPexprGetNChildren(consdata->expr);

      /* the constraint must be an equality constraint with an integer constraint side; also, the root expression
       * needs to be a sum expression with at least two children
       */
      if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) && SCIPisIntegral(scip, consdata->lhs)
         && nchildren > 1 && SCIPexprGetHdlr(consdata->expr) == sumhdlr
         && SCIPisIntegral(scip, SCIPexprsumGetConstant(consdata->expr)) )
      {
         SCIP_Real* coefs;
         SCIP_VAR* cand = NULL;
         SCIP_Bool fail = FALSE;
         int i;

         coefs = SCIPexprsumGetCoefs(consdata->expr);

         /* find candidate variable and check whether all coefficients are integral */
         for( i = 0; i < nchildren; ++i )
         {
            /* check coefficient */
            if( !SCIPisIntegral(scip, coefs[i]) )
            {
               fail = TRUE;
               break;
            }

            if( !SCIPexprIsIntegral(children[i]) )
            {
               /* the child must be a variable expression and the first non-integral expression */
               if( cand != NULL || !SCIPisExprVar(children[i]) || !SCIPisEQ(scip, REALABS(coefs[i]), 1.0) )
               {
                  fail = TRUE;
                  break;
               }

               /* store candidate variable */
               cand = SCIPexprvarGetVar(children[i]);
            }
         }

         if( !fail && cand != NULL )
         {
            SCIPdebugMsg(scip, "make variable <%s> implicit integer due to constraint <%s>\n",
               SCIPvarGetName(cand), SCIPconsGetName(conss[c]));

            /* change variable type */
            SCIP_CALL( SCIPchgVarType(scip, cand, SCIP_VARTYPE_IMPLINT, infeasible) );

            if( *infeasible )
               return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** presolving method to fix a variable x_i to one of its bounds if the variable is only contained in a single
 *  expression contraint g(x) <= rhs (>= lhs) if g is concave (convex) in x_i;  if a continuous variable has bounds
 *  [0,1], then the variable type is changed to be binary; otherwise a bound disjunction constraint is added
 *
 *  @todo the same reduction can be applied if g(x) is not concave, but monotone in x_i for g(x) <= rhs
 *  @todo extend this to cases where a variable can appear in a monomial with an exponent, essentially relax
 *    g(x) to sum_i [a_i,b_i] x^{p_i} for a single variable x and try to conclude montonicity or convexity/concavity
 *    on this (probably have one or two flags per variable and update this whenever another x^{p_i} is found)
 */
static
SCIP_RETCODE presolSingleLockedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   int*                  nchgvartypes,       /**< pointer to store the total number of changed variable types */
   int*                  naddconss,          /**< pointer to store the total number of added constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   SCIP_EXPR** singlelocked;
   SCIP_EXPRHDLR* prodhdlr;
   SCIP_EXPRHDLR* powhdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* exprcands;
   SCIP_Bool hasbounddisj;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int nsinglelocked = 0;
   int i;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(nchgvartypes != NULL);
   assert(naddconss != NULL);
   assert(infeasible != NULL);

   *nchgvartypes = 0;
   *naddconss = 0;
   *infeasible = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* only consider constraints with one finite side */
   if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
      return SCIP_OKAY;

   /* only consider sum expressions */
   if( consdata->expr == NULL || SCIPexprGetHdlr(consdata->expr) != SCIPgetExprHdlrSum(conshdlr) )
      return SCIP_OKAY;

   /* remember which side is finite */
   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* get product and power handlers */
   prodhdlr = SCIPgetExprHdlrProduct(conshdlr);
   powhdlr = SCIPgetExprHdlrPower(conshdlr);

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(&exprcands, SCIPblkmem(scip), consdata->nvarexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &singlelocked, consdata->nvarexprs) );

   /* check all variable expressions for single locked variables */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);

      if( isSingleLockedCand(scip, consdata->varexprs[i]) )
      {
         SCIP_CALL( SCIPhashmapInsert(exprcands, (void*)consdata->varexprs[i], NULL) );
         singlelocked[nsinglelocked++] = consdata->varexprs[i];
      }
   }
   SCIPdebugMsg(scip, "found %d single locked variables for constraint %s\n", nsinglelocked, SCIPconsGetName(cons));

   if( nsinglelocked > 0 )
   {
      SCIP_EXPR** children;
      SCIP_EXPRITER* it;
      int nchildren;

      children = SCIPexprGetChildren(consdata->expr);
      nchildren = SCIPexprGetNChildren(consdata->expr);

      /* create iterator */
      SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
      SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

      for( i = 0; i < nchildren; ++i )
      {
         SCIP_EXPR* child;
         SCIP_Real coef;

         child = children[i];
         assert(child != NULL);
         coef = SCIPexprsumGetCoefs(consdata->expr)[i];

         /* ignore linear terms */
         if( SCIPisExprVar(child) )
            continue;

         /* consider products prod_j f_j(x); ignore f_j(x) if it is a single variable, otherwise iterate through the
          * expression that represents f_j and remove each variable expression from exprcands
          */
         else if( SCIPexprGetHdlr(child) == prodhdlr )
         {
            int j;

            for( j = 0; j < SCIPexprGetNChildren(child); ++j )
            {
               SCIP_EXPR* grandchild = SCIPexprGetChildren(child)[j];

               if( !SCIPisExprVar(grandchild) )
               {
                  /* mark all variable expressions that are contained in the expression */
                  SCIP_CALL( removeSingleLockedVars(scip, grandchild, it, exprcands) );
               }
            }
         }
         /* fixing a variable x to one of its bounds is only valid for ... +x^p >= lhs or ... -x^p <= rhs if p = 2k
          * for an integer k > 1
          */
         else if( SCIPexprGetHdlr(child) == powhdlr )
         {
            SCIP_EXPR* grandchild = SCIPexprGetChildren(child)[0];
            SCIP_Real exponent = SCIPexprpowGetExponent(child);
            SCIP_Bool valid;

            /* check for even integral exponent */
            valid = exponent > 1.0 && fmod(exponent, 2.0) == 0.0;

            if( !valid || !SCIPisExprVar(grandchild) || (hasrhs && coef > 0.0) || (haslhs && coef < 0.0) )
            {
               /* mark all variable expressions that are contained in the expression */
               SCIP_CALL( removeSingleLockedVars(scip, grandchild, it, exprcands) );
            }
         }
         /* all other cases cannot be handled */
         else
         {
            /* mark all variable expressions that are contained in the expression */
            SCIP_CALL( removeSingleLockedVars(scip, child, it, exprcands) );
         }
      }

      /* free expression iterator */
      SCIPfreeExpriter(&it);
   }

   /* check whether the bound disjunction constraint handler is available */
   hasbounddisj = SCIPfindConshdlr(scip, "bounddisjunction") != NULL;

   /* fix variable to one of its bounds by either changing its variable type or adding a disjunction constraint */
   for( i = 0; i < nsinglelocked; ++i )
   {
      /* only consider expressions that are still contained in the exprcands map */
      if( SCIPhashmapExists(exprcands, (void*)singlelocked[i]) )
      {
         SCIP_CONS* newcons;
         SCIP_VAR* vars[2];
         SCIP_BOUNDTYPE boundtypes[2];
         SCIP_Real bounds[2];
         char name[SCIP_MAXSTRLEN];
         SCIP_VAR* var;

         var = SCIPexprvarGetVar(singlelocked[i]);
         assert(var != NULL);
         SCIPdebugMsg(scip, "found single locked variable %s in [%g,%g] that can be fixed to one of its bounds\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));

         /* try to change the variable type to binary */
         if( conshdlrdata->checkvarlocks == 't' && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0) && SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0) )
         {
            assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            ++(*nchgvartypes);

            if( *infeasible )
            {
               SCIPdebugMsg(scip, "detect infeasibility after changing variable type of <%s>\n", SCIPvarGetName(var));
               break;
            }
         }
         /* add bound disjunction constraint if bounds of the variable are finite */
         else if( hasbounddisj && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            vars[0] = var;
            vars[1] = var;
            boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
            boundtypes[1] = SCIP_BOUNDTYPE_UPPER;
            bounds[0] = SCIPvarGetUbGlobal(var);
            bounds[1] = SCIPvarGetLbGlobal(var);

            SCIPdebugMsg(scip, "add bound disjunction constraint for %s\n", SCIPvarGetName(var));

            /* create, add, and release bound disjunction constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quadvarbnddisj_%s", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &newcons, name, 2, vars, boundtypes, bounds, TRUE, TRUE,
               TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            ++(*naddconss);
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &singlelocked);
   SCIPhashmapFree(&exprcands);

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExpr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(valid != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* create basic data of constraint handler and include it to scip */
   SCIP_CALL( includeConshdlrExprBasic(scip) );

   /* copy expression and nonlinear handlers */
   SCIP_CALL( copyConshdlrExprExprHdlr(scip, conshdlr, valid) );

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_NLHDLR* nlhdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize);

   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      if( nlhdlr->freehdlrdata != NULL )
      {
         SCIP_CALL( (*nlhdlr->freehdlrdata)(scip, nlhdlr, &nlhdlr->data) );
      }

      /* free clocks */
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->detecttime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->enfotime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->proptime) );
      SCIP_CALL( SCIPfreeClock(scip, &nlhdlr->intevaltime) );

      SCIPfreeMemory(scip, &nlhdlr->name);
      SCIPfreeMemoryNull(scip, &nlhdlr->desc);

      SCIPfreeMemory(scip, &nlhdlr);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->nlhdlrs, conshdlrdata->nlhdlrssize);
   conshdlrdata->nlhdlrssize = 0;

   /* free upgrade functions */
   for( i = 0; i < conshdlrdata->nexprconsupgrades; ++i )
   {
      assert(conshdlrdata->exprconsupgrades[i] != NULL);
      SCIPfreeBlockMemory(scip, &conshdlrdata->exprconsupgrades[i]);  /*lint !e866*/
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->exprconsupgrades, conshdlrdata->exprconsupgradessize);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->canonicalizetime) );

   SCIPqueueFree(&conshdlrdata->reversepropqueue);

   assert(conshdlrdata->vp_randnumgen == NULL);
#ifndef NDEBUG
   for( i = 0; i <= SCIP_MAXVERTEXPOLYDIM; ++i )
      assert(conshdlrdata->vp_lp[i] == NULL);
#endif

   assert(conshdlrdata->branchrandnumgen == NULL);

   assert(SCIPhashmapGetNElements(conshdlrdata->var2expr) == 0);
   SCIPhashmapFree(&conshdlrdata->var2expr);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_NLHDLR* nlhdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* make sure current activity tags in expressions are invalid, because we start catching variable events only now */
   conshdlrdata->lastboundrelax = ++conshdlrdata->curboundstag;
   /* set to 1 so it is larger than initial value of lastenforound in exprs */
   conshdlrdata->enforound = 1;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(conss[i])) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
   }

   /* sort nonlinear handlers by detection priority, in decreasing order */
   if( conshdlrdata->nnlhdlrs > 1 )
      SCIPsortDownPtr((void**)conshdlrdata->nlhdlrs, nlhdlrCmp, conshdlrdata->nnlhdlrs);

   /* get heuristics for later use */
   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   /* reset statistics in expression handlers */
   for( i = 0; i < conshdlrdata->nexprhdlrs; ++i )
   {
      exprhdlr = conshdlrdata->exprhdlrs[i];
      assert(exprhdlr != NULL);

      exprhdlr->nestimatecalls = 0;
      exprhdlr->nintevalcalls = 0;
      exprhdlr->npropcalls = 0;
      exprhdlr->ncutsfound = 0;
      exprhdlr->ncutoffs = 0;
      exprhdlr->ndomreds = 0;
      exprhdlr->nbranchscores = 0;
      exprhdlr->nsimplifycalls = 0;
      exprhdlr->nsimplified = 0;

      SCIP_CALL( SCIPresetClock(scip, exprhdlr->estimatetime) );
      SCIP_CALL( SCIPresetClock(scip, exprhdlr->proptime) );
      SCIP_CALL( SCIPresetClock(scip, exprhdlr->intevaltime) );
      SCIP_CALL( SCIPresetClock(scip, exprhdlr->simplifytime) );
   }

   /* reset statistics in nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      nlhdlr->nenfocalls = 0;
      nlhdlr->nintevalcalls = 0;
      nlhdlr->npropcalls = 0;
      nlhdlr->nseparated = 0;
      nlhdlr->ncutoffs = 0;
      nlhdlr->ndomreds = 0;
      nlhdlr->nbranchscores = 0;
      nlhdlr->ndetections = 0;
      nlhdlr->ndetectionslast = 0;

      SCIP_CALL( SCIPresetClock(scip, nlhdlr->detecttime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->enfotime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->proptime) );
      SCIP_CALL( SCIPresetClock(scip, nlhdlr->intevaltime) );

      if( nlhdlr->init != NULL )
      {
         SCIP_CALL( (*nlhdlr->init)(scip, nlhdlr) );
      }
   }

   /* reset statistics in constraint handler */
   conshdlrdata->nweaksepa = 0;
   conshdlrdata->ntightenlp = 0;
   conshdlrdata->ndesperatebranch = 0;
   conshdlrdata->ndesperatecutoff = 0;
   conshdlrdata->ndesperatetightenlp = 0;
   conshdlrdata->nforcelp = 0;
   SCIP_CALL( SCIPresetClock(scip, conshdlrdata->canonicalizetime) );

#ifdef ENFOLOGFILE
   ENFOLOG( enfologfile = fopen(ENFOLOGFILE, "w"); )
#endif

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS** consssorted;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( nconss > 0 )
   {
      /* for better performance of dropVarEvents, we sort by index, descending */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consssorted, conss, nconss) );
      SCIPsortDownPtr((void**)consssorted, SCIPcompareIndexConsNonlinear, nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
         SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(consssorted[i])) );
      }

      SCIPfreeBufferArray(scip, &consssorted);
   }

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   if( conshdlrdata->vp_randnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->vp_randnumgen);

   /* free LPs used to construct facets of envelops of vertex-polyhedral functions */
   for( i = 0; i <= SCIP_MAXVERTEXPOLYDIM; ++i )
   {
      if( conshdlrdata->vp_lp[i] != NULL )
      {
         SCIP_CALL( SCIPlpiFree(&conshdlrdata->vp_lp[i]) );
      }
   }

   if( conshdlrdata->branchrandnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->branchrandnumgen);

   /* deinitialize nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_NLHDLR* nlhdlr;

      nlhdlr = conshdlrdata->nlhdlrs[i];
      if( nlhdlr->exit != NULL )
      {
         SCIP_CALL( (*nlhdlr->exit)(scip, nlhdlr) );
      }
   }

   ENFOLOG(
      if( enfologfile != NULL )
      {
         fclose(enfologfile);
         enfologfile = NULL;
      }
   )

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreExpr)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreExpr)
{  /*lint --e{715}*/
   SCIP_Bool infeasible;

   if( nconss == 0 )
      return SCIP_OKAY;

   /* skip some extra work if already known to be infeasible */
   if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
      return SCIP_OKAY;

   /* simplify constraints and replace common subexpressions */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, nconss, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );

   /* currently SCIP does not offer to communicate this,
    * but at the moment this can only become true if canonicalizeConstraints called detectNlhdlrs (which it doesn't do in EXITPRESOLVE stage)
    * or if a constraint expression became constant
    */
   assert(!infeasible);

   /* tell SCIP that we have something nonlinear */
   SCIPenableNLP(scip);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* skip a number of initializations if we have solved already
    * if infeasibility was found by our boundtightening, then curvature check may also fail as some exprhdlr (e.g., pow)
    * assumes nonempty activities in expressions
    */
   if( SCIPgetStatus(scip) != SCIP_STATUS_OPTIMAL && SCIPgetStatus(scip) != SCIP_STATUS_INFEASIBLE &&
      SCIPgetStatus(scip) != SCIP_STATUS_UNBOUNDED && SCIPgetStatus(scip) != SCIP_STATUS_INFORUNBD )
   {
      int i;

      /* reset one of the number of detections counter to count only current round */
      for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
         conshdlrdata->nlhdlrs[i]->ndetectionslast = 0;

      SCIP_CALL( initSolve(scip, conshdlr, conss, nconss) );
   }

   if( conshdlrdata->branchpscostweight > 0.0 )
   {
      SCIP_CALL( SCIPgetCharParam(scip, "branching/lpgainnormalize", &(conshdlrdata->branchpscostupdatestrategy)) );
      if( strchr("lds", conshdlrdata->branchpscostupdatestrategy) == NULL )
      {
         SCIPerrorMessage("branching/lpgainnormalize strategy %c unknown\n", conshdlrdata->branchpscostupdatestrategy);
         SCIPABORT();
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( deinitSolve(scip, conshdlr, conss, nconss) );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free hash table for bilinear terms */
   SCIP_CALL( bilinearTermsFree(scip, conshdlrdata) );

   /* reset flag to allow another call of presolSingleLockedVars() after a restart */
   conshdlrdata->checkedvarlocks = FALSE;

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExpr)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->expr != NULL);

   /* constraint locks should have been removed */
   assert((*consdata)->nlockspos == 0);
   assert((*consdata)->nlocksneg == 0);

   /* free variable expressions */
   SCIP_CALL( freeVarExprs(scip, *consdata) );

   SCIP_CALL( SCIPreleaseExpr(scip, &(*consdata)->expr) );

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransExpr)
{  /*lint --e{715}*/
   SCIP_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( copyExpr(scip, scip, conshdlr, sourcedata->expr, &targetexpr, transformVar, NULL) );
   assert(targetexpr != NULL);  /* copyExpr cannot fail if source and target scip are the same */

   /* create transformed cons (only captures targetexpr, no need to copy again) */
   SCIP_CALL( createConsExpr(scip, conshdlr, targetcons, SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs, FALSE,
      SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
      SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
      SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
      SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons)) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpExpr)
{
   /* create auxiliary variables and call separation initialization callbacks of the expression handlers
    * TODO if we ever want to allow constraints that are separated but not initial, then we need to call initSepa also
    *   during SEPALP, ENFOLP, etc, whenever a constraint may be separated the first time
    *   for now, there is an assert in detectNlhdlrs to require initial if separated
    */
   SCIP_CALL( initSepa(scip, conshdlr, conss, nconss, infeasible) );

   /* collect all bilinear terms for which an auxvar is present
    * TODO this will only do something for the first call of initlp after initsol, because it cannot handle
    * addition (and removal?) of constraints during solve
    * this is typically the majority of constraints, but the method should be made more flexible
    */
   SCIP_CALL( bilinearTermsInsertAll(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consSepa(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consSepa(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consEnfo(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxExpr)
{  /*lint --e{715}*/
   SCIP_CALL( consEnfo(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   SCIP_RESULT propresult;
   unsigned int soltag;
   int nchgbds;
   int nnotify;
   int c;

   soltag = ++conshdlrdata->exprlastsoltag;

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], NULL, soltag) );

      if( isConsViolated(scip, conss[c]) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* try to propagate
    * TODO obey propinenfo parameter, but we need something to recognize cutoff
    */
   nchgbds = 0;
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds) );

   if( (propresult == SCIP_CUTOFF) || (propresult == SCIP_REDUCEDDOM) )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* register all unfixed variables in all violated constraints as branching candidates */
   SCIP_CALL( registerBranchingCandidatesAllUnfixed(scip, conshdlr, conss, nconss, &nnotify) );
   if( nnotify > 0 )
   {
      SCIPdebugMsg(scip, "registered %d external branching candidates\n", nnotify);

      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "could not find branching candidates, forcing to solve LP\n");
   *result = SCIP_SOLVELP;
   ++conshdlrdata->nforcelp;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   SCIP_Bool          maypropfeasible;
   unsigned int soltag;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;
   soltag = ++(conshdlrdata->exprlastsoltag);
   maxviol = 0.0;
   maypropfeasible = conshdlrdata->trysolheur != NULL && SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED
      && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING;

   /* check nonlinear constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
      {
         *result = SCIP_INFEASIBLE;
         maxviol = MAX(maxviol, getConsAbsViolation(conss[c]));  /*lint !e666*/

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* print reason for infeasibility */
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( consdata->lhsviol > SCIPfeastol(scip) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhsviol);
            }
            if( consdata->rhsviol > SCIPfeastol(scip) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", consdata->rhsviol);
            }
         }
         else if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible && !completely )
         {
            /* if we don't want to pass to subnlp heuristic and don't need to print reasons, then can stop checking here */
            return SCIP_OKAY;
         }

         /* do not try to shift linear variables if violation is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, getConsAbsViolation(conss[c])) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            if( consdata->lhsviol > SCIPfeastol(scip) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such variable, then we cannot get feasible
                */
               if( !(consdata->linvarincr != NULL && consdata->linvarincrcoef > 0.0) &&
                  !(consdata->linvardecr != NULL && consdata->linvardecrcoef < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(consdata->rhsviol > SCIPfeastol(scip));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such variable, then we cannot get feasible
                */
               if( !(consdata->linvarincr != NULL && consdata->linvarincrcoef < 0.0) &&
                  !(consdata->linvardecr != NULL && consdata->linvardecrcoef > 0.0) )
                  maypropfeasible = FALSE;
            }
         }
      }
   }

   if( *result == SCIP_INFEASIBLE && maypropfeasible )
   {
      SCIP_Bool success;

      SCIP_CALL( proposeFeasibleSolution(scip, conshdlr, conss, nconss, sol, &success) );

      /* do not pass solution to NLP heuristic if we made it feasible this way */
      if( success )
         return SCIP_OKAY;
   }

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL && sol != NULL && !SCIPisInfinity(scip, maxviol) )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropExpr)
{  /*lint --e{715}*/
   int nchgbds = 0;

   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, SCIPgetDepth(scip) == 0, result, &nchgbds) );
   assert(nchgbds >= 0);

   /* TODO would it make sense to check for redundant constraints? */

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible;
   int c;

   *result = SCIP_DIDNOTFIND;

   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* simplify constraints and replace common subexpressions, reinit nlhdlrs */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, nconss, presoltiming, &infeasible, ndelconss, naddconss, nchgcoefs) );
   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* merge constraints with the same root expression */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
   {
      SCIP_Bool success;

      SCIP_CALL( presolMergeConss(scip, conss, nconss, &success) );
      if( success )
         *result = SCIP_SUCCESS;
   }

   /* propagate constraints */
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, (presoltiming & (SCIP_PRESOLTIMING_MEDIUM | SCIP_PRESOLTIMING_EXHAUSTIVE)) != 0, result, nchgbds) );
   if( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   /* propagate function domains (TODO integrate with simplify?) */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) || nrounds == 0 )
   {
      SCIP_RESULT localresult;
      SCIP_CALL( propExprDomains(scip, conshdlr, conss, nconss, &localresult, nchgbds) );
      if( localresult == SCIP_CUTOFF )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( localresult == SCIP_REDUCEDDOM )
         *result = SCIP_REDUCEDDOM;
   }

   /* check for redundant constraints, remove constraints that are a value expression */
   SCIP_CALL( checkRedundancyConss(scip, conshdlr, conss, nconss, &infeasible, ndelconss, nchgbds) );
   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* try to upgrade constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_Bool upgraded;

      /* skip inactive and deleted constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsActive(conss[c]) )
         continue;

      SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );  /*lint !e794*/
   }

   /* try to change continuous variables that appear linearly to be implicit integer */
   if( (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
   {
      SCIP_CALL( presolveImplint(scip, conshdlr, conss, nconss, nchgvartypes, &infeasible) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, "presolveImplint() detected infeasibility\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* fix variables that are contained in only one expression constraint to their upper or lower bounds, if possible */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && SCIPisPresolveFinished(scip)
      && !conshdlrdata->checkedvarlocks && conshdlrdata->checkvarlocks != 'd' )
   {
      /* run this presolving technique only once because we don't want to generate identical bound disjunction
       * constraints multiple times
       */
      conshdlrdata->checkedvarlocks = TRUE;

      for( c = 0; c < nconss; ++c )
      {
         int tmpnchgvartypes = 0;
         int tmpnaddconss = 0;

         SCIP_CALL( presolSingleLockedVars(scip, conshdlr, conss[c], &tmpnchgvartypes, &tmpnaddconss, &infeasible) );
         SCIPdebugMsg(scip, "presolSingleLockedVars() for %s: nchgvartypes=%d naddconss=%d infeas=%u\n",
            SCIPconsGetName(conss[c]), tmpnchgvartypes, tmpnaddconss, infeasible);

         if( infeasible )
         {
            SCIPdebugMsg(scip, "presolSingleLockedVars() detected infeasibility\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         (*nchgvartypes) += tmpnchgvartypes;
         (*naddconss) += tmpnaddconss;
      }
   }

   if( *ndelconss > 0 || *nchgbds > 0 || *nupgdconss > 0 || *naddconss > 0 || *nchgvartypes > 0 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropExpr NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool reinitsolve = FALSE;

   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->expr == NULL )
      return SCIP_OKAY;

   /* check whether we need to initSolve again because
    * - we have enfo initialized (nenfos >= 0)
    * - and locks appeared (going from zero to nonzero) or disappeared (going from nonzero to zero) now
    */
   if( consdata->expr->nenfos >= 0 )
   {
      if( (consdata->nlockspos == 0) != (nlockspos == 0) )
         reinitsolve = TRUE;
      if( (consdata->nlocksneg == 0) != (nlocksneg == 0) )
         reinitsolve = TRUE;
   }

   if( reinitsolve )
   {
      SCIP_CALL( deinitSolve(scip, conshdlr, &cons, 1) );
   }

   /* add locks */
   SCIP_CALL( addLocks(scip, cons, nlockspos, nlocksneg) );

   if( reinitsolve )
   {
      SCIP_CALL( initSolve(scip, conshdlr, &cons, 1) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions */
   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );
   }

   /* simplify root expression if the constraint has been added after presolving */
   if( SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE )
   {
      if( !consdata->issimplified )
      {
         SCIP_EXPR* simplified;
         SCIP_Bool changed;

         /* simplify constraint */
         SCIP_CALL( SCIPsimplifyExpr(scip, conshdlr, consdata->expr, &simplified, &changed, &infeasible) );
         SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );
         assert(simplified != NULL);
         consdata->expr = simplified;
         consdata->issimplified = TRUE;
      }
   }

   /* add manually locks to constraints that are not checked for feasibility */
   if( !SCIPconsIsChecked(cons) )
   {
      assert(consdata->nlockspos == 0);
      assert(consdata->nlocksneg == 0);

      SCIP_CALL( addLocks(scip, cons, 1, 0) );
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_INITPRESOLVE && !infeasible )
   {
      SCIP_CALL( initSolve(scip, conshdlr, &cons, 1) );
   }

   /* TODO deal with infeasibility */
   assert(!infeasible);

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) < SCIP_STAGE_EXITSOLVE )
   {
      SCIP_CALL( deinitSolve(scip, conshdlr, &cons, 1) );
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(cons)) );
   }

   /* remove locks that have been added in consActiveExpr() */
   if( !SCIPconsIsChecked(cons) )
   {
      SCIP_CALL( addLocks(scip, cons, -1, 0) );

      assert(SCIPconsGetData(cons)->nlockspos == 0);
      assert(SCIPconsGetData(cons)->nlocksneg == 0);
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsExpr NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExpr)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged constraints */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print expression */
   if( consdata->expr != NULL )
   {
      SCIP_CALL( SCIPprintExpr(scip, conshdlr, consdata->expr, file) );
   }
   else
   {
      SCIPinfoMessage(scip, file, "0");
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyExpr)
{  /*lint --e{715}*/
   COPY_MAPVAR_DATA mapvardata;
   SCIP_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   assert(cons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   mapvardata.varmap = varmap;
   mapvardata.consmap = consmap;
   mapvardata.global = global;
   mapvardata.valid = TRUE; /* hope the best */

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( copyExpr(sourcescip, scip, sourceconshdlr, sourcedata->expr, &targetexpr, copyVar, &mapvardata) );

   if( targetexpr == NULL )
   {
      *cons = NULL;
      *valid = FALSE;

      return SCIP_OKAY;
   }

   /* validity depends only on the SCIPgetVarCopy() returns from copyVar, which are accumulated in mapvardata.valid */
   *valid = mapvardata.valid;

   /* create copy (only capture targetexpr, no need to copy again) */
   SCIP_CALL( createConsExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), cons, name != NULL ? name : SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs, FALSE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseExpr)
{  /*lint --e{715}*/
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   SCIP_EXPR* consexprtree;

   SCIPdebugMsg(scip, "cons_expr::consparse parsing %s\n",str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = FALSE;

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, (char**)&endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the beginning of the expression and not a left-hand-side */
         lhs = -SCIPinfinity(scip);
      }
      else
      {
         /* it was indeed a left-hand-side, so continue parsing after it */
         str = endptr + 2;

         /* ignore whitespace */
         while( isspace((unsigned char)*str) )
            ++str;
      }
   }

   debugParse("str should start at beginning of expr: %s\n", str); /*lint !e506 !e681*/

   /* parse expression: so far we did not allocate memory, so can just return in case of readerror */
   SCIP_CALL( SCIPparseExpr(scip, conshdlr, str, &str, &consexprtree) );

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
         return SCIP_OKAY;
      }
      *success = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, (char**)&endptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
            return SCIP_OKAY;
      }
   }

   /* create constraint */
   SCIP_CALL( createConsExpr(scip, conshdlr, cons, name,
      consexprtree, lhs, rhs, TRUE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   assert(*cons != NULL);

   SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );

   debugParse("created expression constraint: <%s>\n", SCIPconsGetName(*cons)); /*lint !e506 !e681*/

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions if not done so far */
   SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

   /* check whether array is too small in order to store all variables */
   if( varssize < consdata->nvarexprs )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      vars[i] = SCIPexprvarGetVar(consdata->varexprs[i]);
      assert(vars[i] != NULL);
   }

   *success = TRUE;

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions if not done so far */
   SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

   *nvars = consdata->nvarexprs;
   *success = TRUE;

   return SCIP_OKAY;
}


/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0 /* TODO? */
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsExpr NULL
#endif


/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputExpr)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* print statistics for nonlinear handlers */
   printNlhdlrStatistics(scip, conshdlr, file);

   /* print statistics for constraint handler */
   printConshdlrStatistics(scip, conshdlr, file);

   return SCIP_OKAY;
}

/** returns whether we are ok to branch on auxiliary variables
 *
 * Currently returns whether depth of node in B&B tree is at least value of constraints/expr/branching/aux parameter.
 */
SCIP_Bool SCIPgetBranchAuxNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->branchauxmindepth <= SCIPgetDepth(scip);
}

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_VAR* SCIPgetExprAuxVarNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return SCIPisExprVar(expr) ? SCIPexprvarGetVar(expr) : expr->auxvar;
}

// TODO this should become a cons_nonlinear-specific version that prints ownerdata
/** prints structure of an expression a la Maple's dismantle */
SCIP_RETCODE SCIPdismantleExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_EXPR*     expr              /**< expression to dismantle */
   )
{
   SCIP_EXPRITER* it;
   int depth = -1;

   SCIP_CALL( SCIPexpriterCreate(&it, SCIPfindConshdlr(scip, CONSHDLR_NAME), SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);

   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR:
         {
            int nspaces;
            const char* type;

            ++depth;
            nspaces = 3 * depth;
            type = SCIPexprhdlrGetName(SCIPexprGetHdlr(expr));

            /* use depth of expression to align output */
            SCIPinfoMessage(scip, file, "%*s[%s]: ", nspaces, "", type);

            if( strcmp(type, "var") == 0 )
            {
               SCIP_VAR* var;

               var = SCIPexprvarGetVar(expr);
               SCIPinfoMessage(scip, file, "%s in [%g, %g]", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
                  SCIPvarGetUbLocal(var));
            }
            else if(strcmp(type, "sum") == 0)
               SCIPinfoMessage(scip, file, "%g", SCIPexprsumGetConstant(expr));
            else if(strcmp(type, "prod") == 0)
               SCIPinfoMessage(scip, file, "%g", SCIPexprprodGetCoef(expr));
            else if(strcmp(type, "val") == 0)
               SCIPinfoMessage(scip, file, "%g", SCIPexprvalGetValue(expr));
            else if(strcmp(type, "pow") == 0 || strcmp(type, "signpower") == 0)
               SCIPinfoMessage(scip, file, "%g", SCIPexprpowGetExponent(expr));

            /* print nl handlers associated to expr */
            if(expr->nenfos > 0 )
            {
               int i;
               SCIPinfoMessage(scip, file, "   {");

               for( i = 0; i < expr->nenfos; ++i )
               {
                  SCIPinfoMessage(scip, file, "%s:", expr->enfos[i]->nlhdlr->name);
                  if( expr->enfos[i]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_ACTIVITY )
                     SCIPinfoMessage(scip, file, "a");
                  if( expr->enfos[i]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPABELOW )
                     SCIPinfoMessage(scip, file, "u");
                  if( expr->enfos[i]->nlhdlrparticipation & SCIP_CONSNONLINEAR_EXPRENFO_SEPAABOVE )
                     SCIPinfoMessage(scip, file, "o");
                  if( i < expr->nenfos-1 )
                     SCIPinfoMessage(scip, file, ", ");
               }

               SCIPinfoMessage(scip, file, "}");
            }

            /* print aux var associated to expr */
            if( expr->auxvar != NULL )
               SCIPinfoMessage(scip, file, "  (%s in [%g, %g])", SCIPvarGetName(expr->auxvar),
                     SCIPvarGetLbLocal(expr->auxvar), SCIPvarGetUbLocal(expr->auxvar));
            SCIPinfoMessage(scip, file, "\n");

            break;
         }

         case SCIP_EXPRITER_VISITINGCHILD:
         {
            int nspaces;
            const char* type;

            nspaces = 3 * depth;
            type = SCIPexprhdlrGetName(SCIPexprGetHdlr(expr));

            if( strcmp(type, "sum") == 0 )
            {
               SCIPinfoMessage(scip, file, "%*s   ", nspaces, "");
               SCIPinfoMessage(scip, file, "[coef]: %g\n", SCIPexprsumGetCoefs(expr)[SCIPexpriterGetChildIdxDFS(it)]);
            }

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR:
         {
            --depth;
            break;
         }

         default:
            /* shouldn't be here */
            SCIPABORT();
            break;
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** returns the partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @note expression must belong to a constraint
 */
SCIP_Real SCIPgetExprPartialDiffNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression which has been used in the last SCIPevalExprGradient() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(consexprhdlr), CONSHDLR_NAME) == 0);
   assert(expr != NULL);
   assert(var != NULL);

   /* return 0.0 for value expression */
   if( expr->exprhdlr == SCIPgetExprHdlrValue(consexprhdlr) )
   {
      assert(expr->derivative == 0.0);
      return 0.0;
   }

   /* check if an error occurred during the last SCIPevalExprGradient() call */
   if( expr->derivative == SCIP_INVALID ) /*lint !e777*/
      return SCIP_INVALID;

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   /* use variable to expressions mapping which is stored in the constraint handler data */
   assert(SCIPhashmapExists(conshdlrdata->var2expr, var));

   varexpr = (SCIP_EXPR*)SCIPhashmapGetImage(conshdlrdata->var2expr, var);
   assert(varexpr != NULL);
   assert(SCIPisExprVar(varexpr));

   /* use difftag to decide whether the variable belongs to the expression */
   return (expr->difftag != varexpr->difftag) ? 0.0 : varexpr->derivative;
}

/** returns the var's coordinate of Hu partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @note expression must belong to a constraint
 */
SCIP_Real SCIPgetExprPartialDiffGradientDirNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< root expression of constraint used in the last SCIPevalExprHessianDir() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(consexprhdlr), CONSHDLR_NAME) == 0);
   assert(expr != NULL);
   assert(var != NULL);

   /* return 0.0 for value expression */
   if( expr->exprhdlr == SCIPgetExprHdlrValue(consexprhdlr) )
      return 0.0;

   /* check if an error occurred during the last SCIPevalExprHessianDir() call */
   if( expr->bardot == SCIP_INVALID ) /*lint !e777*/
      return SCIP_INVALID;

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   /* use variable to expressions mapping which is stored in the constraint handler data;
    * if this fails it means that we are asking for the var's component of H*u for a var
    * that doesn't appear in any nonlinear constraint, so maybe we can also just return 0.0
    */
   assert(SCIPhashmapExists(conshdlrdata->var2expr, var));

   varexpr = (SCIP_EXPR*)SCIPhashmapGetImage(conshdlrdata->var2expr, var);
   assert(varexpr != NULL);
   assert(SCIPisExprVar(varexpr));

   /* use difftag to decide whether the variable belongs to the expression */
   return (expr->difftag != varexpr->difftag) ? 0.0 : varexpr->bardot;
}

/** returns the activity of the expression
 *
 * The caller needs to make sure that the activity is valid.
 * For expression and nonlinear handlers, this is made sure when the following callbacks are called:
 * - interval evaluation (intervals for children only)
 * - reverse propagation
 * - estimate and enforce (for exprs where activity usage was signaled during nlhdlr detect)
 */
SCIP_INTERVAL SCIPexprGetActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*     expr              /**< expression */
   )
{
#ifndef NDEBUG
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);

   /* check whether activity is valid */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert(expr->activitytag >= conshdlrdata->lastboundrelax);
#endif

   return expr->activity;
}

/** returns the tag associated with the activity of the expression
 *
 * Can be compared with SCIPgetCurBoundsTagNonlinear() and SCIPgetLastBoundRelaxTagNonlinear()
 * to check whether the activity currently stored in this expression is current and valid, respectively.
 */
unsigned int SCIPexprGetActivityTag(
   SCIP_EXPR*     expr              /**< expression */
   );

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is not valid (some bound was relaxed since last evaluation).
 * If validsufficient is set to FALSE, then it will also reevaluate activity if a bound tightening was happening
 * since last evaluation.
 */
SCIP_RETCODE SCIPevalExprActivity(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler, or NULL */
   SCIP_EXPR*     expr,             /**< expression */
   SCIP_INTERVAL*          activity,         /**< interval where to store expression */
   SCIP_Bool               validsufficient   /**< whether any valid activity is sufficient */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(activity != NULL);

   if( consexprhdlr == NULL )
      consexprhdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   if( expr->activitytag < conshdlrdata->lastboundrelax ||
      (!validsufficient && expr->activitytag < conshdlrdata->curboundstag) )
   {
      /* update activity of expression */
      SCIP_CALL( forwardPropExpr(scip, consexprhdlr, expr, FALSE, NULL, NULL) );

      assert(expr->activitytag == conshdlrdata->curboundstag);
   }

   *activity = expr->activity;

   return SCIP_OKAY;
}

/** returns bounds on the expression
 *
 * This gives an intersection of bounds from
 * - activity calculation (\ref SCIPexprGetActivity), if valid,
 * - auxiliary variable, if present,
 * - stored by \ref SCIPtightenExprIntervalNonlinear during domain propagation
 *
 * @note The returned interval can be empty!
 */
SCIP_INTERVAL SCIPgetExprBoundsNonlinear(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_EXPR*     expr              /**< expression */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL bounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* SCIPdebugMsg(scip, "get bounds expr %p:", expr); */

   /* start with propbounds if they belong to current propagation */
   if( expr->propboundstag == conshdlrdata->curpropboundstag )
   {
      bounds = expr->propbounds;
      /* SCIPdebugMsgPrint(scip, " propbounds [%.15g,%.15g]", expr->propbounds.inf, expr->propbounds.sup); */
   }
   else
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &bounds);

   if( expr->activitytag >= conshdlrdata->lastboundrelax )
   {
      /* apply propbounds to expr activity, but ensure it's not-empty if very close disjoint intervals */
      /* SCIPdebugMsgPrint(scip, " activity [%.15g,%.15g]", expr->activity.inf, expr->activity.sup); */
      SCIPintervalIntersectEps(&bounds, SCIPepsilon(scip), expr->activity, bounds);
   }

   if( expr->auxvar != NULL )
   {
      /* apply auxiliary variable bounds to bounds */
      SCIP_INTERVAL auxvarbounds;

      auxvarbounds = conshdlrdata->intevalvar(scip, expr->auxvar, conshdlrdata);
      /* SCIPdebugMsgPrint(scip, " auxvar [%.15g,%.15g]", auxvarbounds.inf, auxvarbounds.sup); */
      SCIPintervalIntersectEps(&bounds, SCIPepsilon(scip), bounds, auxvarbounds);
   }

   /* SCIPdebugMsgPrint(scip, " -> [%.15g,%.15g]\n", bounds.inf, bounds.sup); */

   return bounds;
}

/** informs the expression about new bounds that can be used for reverse-propagation and to tighten bounds of
 * corresponding (auxiliary) variable (if any)
 *
 * @attention this function should only be called during domain propagation in cons_expr
 */
SCIP_RETCODE SCIPtightenExprIntervalNonlinear(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_EXPR*     expr,             /**< expression to be tightened */
   SCIP_INTERVAL           newbounds,        /**< new bounds for the expression */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a cutoff was detected */
   int*                    ntightenings      /**< buffer to add the total number of tightenings, or NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);

   /* the code below assumes that current activity is valid
    * if it turns out that we cannot ensure that, then we should change code
    */
   assert(expr->activitytag >= SCIPconshdlrGetData(conshdlr)->lastboundrelax || SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, expr->activity));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr->activity));

   *cutoff = FALSE;

#ifdef DEBUG_PROP
   SCIPdebugMsg(scip, "Trying to tighten bounds of expr ");
   SCIP_CALL( SCIPprintExpr(scip, SCIPfindConshdlr(scip, CONSHDLR_NAME), expr, NULL) );
   SCIPdebugMsgPrint(scip, " with activity [%.15g,%.15g] to [%.15g,%.15g] (force=%d)\n", expr->activity.inf, expr->activity.sup, newbounds.inf, newbounds.sup, SCIPconshdlrGetData(conshdlr)->forceboundtightening);
#endif

   if( expr->isintegral )
   {
      /* apply integrality to new bounds
       * it should be ok to use normal ceil() and floor(), but for safety, we use SCIPceil and SCIPfloor for now
       */
      if( newbounds.inf > -SCIP_INTERVAL_INFINITY )
         newbounds.inf = SCIPceil(scip, newbounds.inf);
      if( newbounds.sup <  SCIP_INTERVAL_INFINITY )
         newbounds.sup = SCIPfloor(scip, newbounds.sup);
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, " applied integrality: [%.15g,%.15g]\n", newbounds.inf, newbounds.sup);
#endif
   }

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      SCIPdebugMsg(scip, " cut off due to new bounds being empty\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* treat the new bounds as empty if either the lower/upper bound is above/below +/- SCIPinfinity() */
   if( SCIPisInfinity(scip, newbounds.inf) || SCIPisInfinity(scip, -newbounds.sup) )
   {
      SCIPdebugMsg(scip, " cut off due to new bounds being beyond infinity\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* tighten newbounds w.r.t. existing expr->propbounds or activity */
   if( expr->propboundstag == conshdlrdata->curpropboundstag )
   {
      /* if already having propbounds in expr, then tighten newbounds by propbounds */
      SCIPintervalIntersectEps(&newbounds, SCIPepsilon(scip), expr->propbounds, newbounds);
   }
   else
   {
      /* first time we have propbounds for expr in this propagation rounds:
       * intersect with activity (though don't let it become empty if very close intervals)
       */
      SCIPintervalIntersectEps(&newbounds, SCIPepsilon(scip), expr->activity, newbounds);
   }
#ifdef DEBUG_PROP
   SCIPdebugMsg(scip, " applied %s: [%.20g,%.20g]\n", expr->propboundstag == conshdlrdata->curpropboundstag ? "previous propbounds" : "activity", newbounds.inf, newbounds.sup, (void*)expr);
#endif

   /* check if the new bounds lead to an empty interval */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      SCIPdebugMsg(scip, " cut off due to empty intersection with previous propbounds or activity\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if expr is not constant or variable, then store newbounds in expr->propbounds
    * - for constant, the intersection with activity should have been sufficient to determine infeasibilty
    * - for variable, the tightenAuxVarBounds call below should be suffient to have to new bounds acknowledged
    */
   if( expr->nchildren > 0 )
   {
      expr->propbounds = newbounds;
      expr->propboundstag = conshdlrdata->curpropboundstag;
   }

   /* if updated propbounds do not allow a sufficient tightening, then do not consider adding to queue for reverse
    * propagation or update of auxvar bounds
    * TODO? if we first had a considerable tightening and then only get small tightenings under the same
    *   curpropboundstag, then these will still be considered as isIntervalBetter, since we compare with activity here and
    *   not with the propbounds as set in the beginning; I'm not sure, though, that comparing always with previous
    *   propbounds would be better, since a number of small updates to propbounds could eventually lead to a considerable
    *   one or should we not even update propbounds to newbounds if the update is small?
    */
   if( !isIntervalBetter(scip, conshdlrdata->forceboundtightening, newbounds, expr->activity) )
   {
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, " new bounds [%g,%g] for expr %p not sufficiently tighter than activity -- not adding to propqueue or tightening auxvar\n", newbounds.inf, newbounds.sup, (void*)expr);
#endif
      return SCIP_OKAY;
   }

   if( expr->nchildren > 0 && !expr->inpropqueue && (expr->nactivityusesprop > 0 || expr->nactivityusessepa > 0 || expr->nenfos < 0) )
   {
      /* add expression to propagation queue if not there yet and not var or constant and
       * if it should have a nlhdlr with a reverseprop callback or nlhdlrs are not initialized yet (nenfos < 0)
       */
#ifdef DEBUG_PROP
         SCIPdebugMsg(scip, " insert expr <%p> (%s) into reversepropqueue\n", (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));
#endif
         SCIP_CALL( SCIPqueueInsert(conshdlrdata->reversepropqueue, expr) );
         expr->inpropqueue = TRUE;
   }

   /* update bounds on variable or auxiliary variable */
   SCIP_CALL( tightenAuxVarBounds(scip, conshdlr, expr, newbounds, cutoff, ntightenings) );

   return SCIP_OKAY;
}

/** mark constraints that include this expression to be propagated again
 *
 * This can be used by, e.g., nlhdlrs, to trigger a new propagation of constraints without
 * a change of variable bounds, e.g., because new information on the expression is available
 * that could potentially lead to tighter expression activity values.
 *
 * Note, that this call marks also constraints for propagation which only share some variable
 * with this expression.
 */
SCIP_RETCODE SCIPmarkExprPropagateNonlinear(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*     expr              /**< expression to propagate again */
   )
{
   SCIP_EXPRITER* it;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )  /*lint !e441*/
   {
      if( !SCIPisExprVar(expr) )
         continue;

      conss = SCIPgetConsExprExprVarConss(expr);
      nconss = SCIPgetConsExprExprVarNConss(expr);

      for( c = 0; c < nconss; ++c )
      {
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         consdata->ispropagated = FALSE;
      }
   }

   SCIPfreeExpriter(&it);

   SCIPincrementCurBoundsTagNonlinear(conshdlr, FALSE);

   return SCIP_OKAY;
}

/** increments the curboundstag and resets lastboundrelax in constraint handler data
 *
 * @note This method is not intended for normal use.
 *   These tags are maintained by the event handler for variable bound change events.
 *   This method is used by some unittests.
 */
void SCIPincrementCurBoundsTagNonlinear(
   SCIP_CONSHDLR*          conshdlr,         /**< expression constraint handler */
   SCIP_Bool               boundrelax        /**< indicates whether a bound was relaxed, i.e., lastboundrelax should be set too */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ++conshdlrdata->curboundstag;
   assert(conshdlrdata->curboundstag > 0);

   if( boundrelax )
      conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
}

/** adds violation-branching score to an expression
 *
 * Adds a score to the expression-specific violation-branching score, thereby marking it as branching candidate.
 * The expression must either be a variable expression or have an aux-variable.
 * In the latter case, branching on auxiliary variables must have been enabled.
 * In case of doubt, use SCIPaddExprsViolScoreNonlinear(). Roughly, the difference between these functions is that the current
 * function adds the violscore to the expression directly, while SCIPaddExprsViolScoreNonlinear() will split the
 * violation score among all the given expressions according to constraints/expr/branching/violsplit. See
 * SCIPaddExprsViolScoreNonlinear() for more details.
 */
void SCIPaddExprViolScoreNonlinear(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_EXPR*     expr,             /**< expression where to add branching score */
   SCIP_Real               violscore         /**< violation score to add to expression */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(violscore >= 0.0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if not allowing to branch on auxvars, then expr must be a var-expr */
   assert(SCIPgetBranchAuxNonlinear(scip, conshdlr) || expr->exprhdlr == conshdlrdata->exprvarhdlr);
   /* if allowing to branch on auxvars, then expr must be a var-expr or have an auxvar */
   assert(!SCIPgetBranchAuxNonlinear(scip, conshdlr) || (expr->exprhdlr == conshdlrdata->exprvarhdlr || expr->auxvar != NULL));

   /* reset branching score if we are in a different enfo round */
   if( expr->violscoretag != conshdlrdata->enforound )
   {
      expr->violscoresum = violscore;
      expr->violscoremax = violscore;
      expr->nviolscores = 1;
      expr->violscoretag = conshdlrdata->enforound;
      return;
   }

   expr->violscoresum += violscore;
   if( violscore > expr->violscoremax )
      expr->violscoremax = violscore;
   ++(expr->nviolscores);
}

/** adds violation-branching score to a set of expressions, distributing the score among all the expressions.
 *
 * Each expression must either be a variable expression or have an aux-variable.
 * If branching on aux-variables is disabled, then the violation branching score will be distributed among all among the
 * variables present in exprs
 */
SCIP_RETCODE SCIPaddExprsViolScoreNonlinear(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          conshdlr,         /**< expr constraint handler */
   SCIP_EXPR**    exprs,            /**< expressions where to add branching score */
   int                     nexprs,           /**< number of expressions */
   SCIP_Real               violscore,        /**< violation-branching score to add to expression */
   SCIP_SOL*               sol,              /**< current solution */
   SCIP_Bool*              success           /**< buffer to store whether at least one branchscore was added */
   )
{
   /* distribute violation as branching score to original variables in children of expr that are marked in branchcand */
   SCIP_EXPRITER* it;
   SCIP_EXPR** varexprs;
   SCIP_EXPR* e;
   int nvars;
   int varssize;
   int i;

   assert(exprs != NULL || nexprs == 0);
   assert(success != NULL);

   if( nexprs == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* if allowing to branch on auxiliary variables, then call internal addConsExprExprsViolScore immediately */
   if( SCIPgetBranchAuxNonlinear(scip, conshdlr) )
   {
      addConsExprExprsViolScore(scip, conshdlr, exprs, nexprs, violscore, sol, success);
      return SCIP_OKAY;
   }

   /* if not allowing to branch on aux vars, then create new array containing var expressions that exprs depend on */
   nvars = 0;
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, varssize) );

   SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( i = 0; i < nexprs; ++i )
   {
      for( e = SCIPexpriterRestartDFS(it, exprs[i]); !SCIPexpriterIsEnd(it); e = SCIPexpriterGetNext(it) ) /*lint !e441*/
      {
         assert(e != NULL);

         if( SCIPisExprVar(e) )
         {
            /* add variable expression to vars array */
            if( varssize == nvars )
            {
               varssize = SCIPcalcMemGrowSize(scip, nvars + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, &varexprs, varssize) );
            }
            assert(varssize > nvars);

            varexprs[nvars++] = e;
         }
      }
   }

   SCIPfreeExpriter(&it);

   addConsExprExprsViolScore(scip, conshdlr, varexprs, nvars, violscore, sol, success);

   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}

/** gives violation-branching score stored in expression, or 0.0 if no valid score has been stored */
SCIP_Real SCIPgetExprViolScoreNonlinear(
   SCIP_CONSHDLR*          conshdlr,         /**< constraint handler */
   SCIP_EXPR*     expr              /**< expression */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(expr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->enforound != expr->violscoretag )
      return 0.0;

   if( expr->nviolscores == 0 )
      return 0.0;

   switch( conshdlrdata->branchscoreagg )
   {
      case 'a' :
         /* average */
         return expr->violscoresum / expr->nviolscores;

      case 'm' :
         /* maximum */
         return expr->violscoremax;

      case 's' :
         /* sum */
         return expr->violscoresum;

      default:
         SCIPerrorMessage("Invalid value %c for branchscoreagg parameter\n", conshdlrdata->branchscoreagg);
         SCIPABORT();
         return SCIP_INVALID;
   }
}

/** returns the number of positive rounding locks of an expression */
int SCIPgetExprNLocksPosNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nlockspos;
}

/** returns the number of negative rounding locks of an expression */
int SCIPgetExprNLocksNegNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nlocksneg;
}

/** number of nonlinear handlers whose activity computation and propagation methods depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
unsigned int SCIPgetExprNPropUsesActivityNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nactivityusesprop;
}

/** number of nonlinear handlers whose separation methods (estimate or enforcement) depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
unsigned int SCIPgetExprNSepaUsesActivityNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nactivityusessepa;
}

/** number of nonlinear handlers whose separation methods (estimate or enforcement) use auxiliary variable of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
unsigned int SCIPgetExprNAuxvarUsesNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nauxvaruses;
}

/** method to be called by a nlhdlr during NLHDLRDETECT to notify expression that it will be used
 *
 * - if useauxvar is enabled, then ensures that an auxiliary variable will be created in INITLP
 * - if useactivityforprop or useactivityforsepa{below,above} is enabled, then ensured that activity will be updated for expr
 * - if useactivityforprop is enabled, then increments the count returned by \ref SCIPgetExprNPropUsesActivityNonlinear
 * - if useactivityforsepa{below,above} is enabled, then increments the count returned by \ref SCIPgetExprNSepaUsesActivityNonlinear
 *   and also increments this count for all variables in the expression.
 *
 * The distinction into useactivityforprop and useactivityforsepa{below,above} is to recognize variables which domain influences
 * under/overestimators. Domain propagation routines (like OBBT) may invest more work for these variables.
 * The distinction into useactivityforsepabelow and useactivityforsepaabove is to recognize whether a nlhdlr that called this method
 * will use activity of expr in enfomethod sepabelow or enfomethod sepaabove.
 */
SCIP_RETCODE SCIPregisterExprUsageNonlinear(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,         /**< expression constraint handler */
   SCIP_EXPR*   expr,             /**< expression */
   SCIP_Bool             useauxvar,        /**< whether an auxiliary variable will be used for estimate or cut generation */
   SCIP_Bool             useactivityforprop, /**< whether activity of expr will be used by domain propagation or activity calculation (inteval) */
   SCIP_Bool             useactivityforsepabelow, /**< whether activity of expr will be used by underestimation */
   SCIP_Bool             useactivityforsepaabove  /**< whether activity of expr will be used by overestimation */
   )
{
   assert(conshdlr != NULL);
   assert(expr != NULL);

   /* do not store auxvar request for variable expressions */
   if( useauxvar && SCIPisExprVar(expr) )
      useauxvar = FALSE;

   if( expr->nenfos >= 0 &&
      ( (expr->nactivityusesprop == 0 && expr->nactivityusessepa == 0 && (useactivityforprop || useactivityforsepabelow || useactivityforsepaabove)) ||
        (expr->nauxvaruses == 0 && useauxvar)
      ) )
   {
      /* if we already have ran detect of nlhdlrs on expr (nenfos >= 0), then we need to rerun detection if
       * we require additional enforcement methods, that is,
       * - activity of expr was not used before but will be used now, or
       * - auxiliary variable of expr was not required before but will be used now
       */
      SCIP_CALL( freeEnfoData(scip, expr, FALSE) );
   }

   if( useauxvar )
      ++(expr->nauxvaruses);

   if( useactivityforprop )
      ++(expr->nactivityusesprop);

   if( useactivityforsepabelow || useactivityforsepaabove )
      ++(expr->nactivityusessepa);

   /* remember that SCIPregisterExprUsageNonlinear() has been called with useactivityforsepa{below,above}=TRUE; this
    * information is used in detectNlhdlr()
    */
   if( useactivityforsepabelow )
      SCIPconshdlrGetData(conshdlr)->registerusesactivitysepabelow = TRUE;
   if( useactivityforsepaabove )
      SCIPconshdlrGetData(conshdlr)->registerusesactivitysepaabove = TRUE;

   if( useactivityforprop )
   {
      /* if activity will be used for propagation, then make sure there is a valid activity
       * this way, we can do a reversepropcall after detectNlhdlr
       */
      SCIP_INTERVAL activity;
      SCIP_CALL( SCIPevalExprActivity(scip, conshdlr, expr, &activity, FALSE) );
   }

   /* increase the nactivityusedsepa counter for all variables used in the given expression */
   if(( useactivityforsepabelow || useactivityforsepaabove) && SCIPexprGetNChildren(expr) > 0 )
   {
      SCIP_EXPRITER* it;

      /* create and initialize iterator */
      SCIP_CALL( SCIPexpriterCreate(&it, conshdlr, SCIPblkmem(scip)) );
      SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

      for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
         if( SCIPisExprVar(expr) )
            ++(expr->nactivityusessepa);

      /* free iterator */
      SCIPfreeExpriter(&it);
   }

   return SCIP_OKAY;
}

/** computes absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(x) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(x) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(x).
 *
 * If necessary, f is evaluated in the given solution. If that fails (domain error),
 * then viol is set to SCIPinfinity and both violover and violunder are set to TRUE.
 */
SCIP_RETCODE SCIPgetExprAbsOrigViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution */
   unsigned int          soltag,             /**< tag of solution */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* make sure expression has been evaluated */
   SCIP_CALL( SCIPevalExpr(scip, conshdlr, expr, sol, soltag) );

   /* get violation from internal method */
   *viol = getExprAbsOrigViolation(scip, expr, sol, violunder, violover);

   return SCIP_OKAY;
}

/** computes absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z <= f(w) and sets violover to TRUE.
 * If there are positive locks, then return the violation of z >= f(w) and sets violunder to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z == f(w).
 *
 * If the given value of f(w) is SCIP_INVALID, then viol is set to SCIPinfinity and
 * both violover and violunder are set to TRUE.
 */
SCIP_RETCODE SCIPgetExprAbsAuxViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* get violation from internal method */
   *viol = getExprAbsAuxViolation(scip, expr, auxvalue, sol, violunder, violover);

   return SCIP_OKAY;
}

/** computes relative violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * Taking the absolute violation from SCIPgetExprAbsAuxViolationNonlinear, this function returns
 * the absolute violation divided by max(1,|f(w)|).
 *
 * If the given value of f(w) is SCIP_INVALID, then viol is set to SCIPinfinity and
 * both violover and violunder are set to TRUE.
 */
SCIP_RETCODE SCIPgetExprRelAuxViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* get violation from internal method */
   *viol = getExprAbsAuxViolation(scip, expr, auxvalue, sol, violunder, violover);

   if( !SCIPisInfinity(scip, *viol) )
   {
      assert(auxvalue != SCIP_INVALID);  /*lint !e777*/
      *viol /= MAX(1.0, REALABS(auxvalue));  /*lint !e666*/
   }

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** gets tag indicating current local variable bounds */
unsigned int SCIPgetCurBoundsTagNonlinear(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);

   return conshdlrdata->curboundstag;
}

/** gets the curboundstag at the last time where variable bounds were relaxed */
unsigned int SCIPgetLastBoundRelaxTagNonlinear(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);

   return conshdlrdata->lastboundrelax;
}

/** returns the hashmap that is internally used to map variables to their corresponding variable expressions */
SCIP_HASHMAP* SCIPgetVarExprHashmapNonlinear(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   assert(consexprhdlr != NULL);

   return SCIPconshdlrGetData(consexprhdlr)->var2expr;
}

/** notifies conshdlr that a variable expression is to be freed
 *
 * the conshdlr will then update its var2expr hashmap
 *
 * @note To be called only by var-exprhdlr.
 * @note Temporary method that will be replaced by ownerdata-free
 */
SCIP_RETCODE SCIPnotifyExprVarFreedNonlinear(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_EXPR*        varexpr         /**< variable expression to be freed */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(varexpr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   var = SCIPexprvarGetVar(varexpr);
   assert(var != NULL);

   /* if no variable-expression stored for var hashmap, then the var hasn't been used in any constraint, so do nothing
    * if variable-expression stored for var is different, then also do nothing
    */
   if( SCIPhashmapGetImage(conshdlrdata->var2expr, var) != (void*)varexpr )
      return SCIP_OKAY;

   /* remove var -> varexpr map from hashmap */
   SCIP_CALL( SCIPhashmapRemove(conshdlrdata->var2expr, var) );

   return SCIP_OKAY;
}

/** collects all bilinear terms for a given set of constraints
 *
 * @note This method should only be used for unit tests that depend on SCIPgetBilinTermsNonlinear(),
 *       SCIPgetBilinTermNonlinear() or SCIPgetBilinTermIdxNonlinear().
 */
SCIP_RETCODE SCIPcollectBilinTermsNonlinear(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_CONS**                conss,          /**< expression constraints */
   int                        nconss          /**< total number of expression constraints */
   )
{
   assert(consexprhdlr != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( bilinearTermsInsertAll(scip, consexprhdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** returns the total number of bilinear terms that are contained in all expression constraints
 *
 *  @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
int SCIPgetNBilinTermsNonlinear(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nbilinterms;
}

/** returns all bilinear terms that are contained in all expression constraints
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @note The value of the auxiliary variable of a bilinear term might be NULL, which indicates that the term does not have an auxiliary variable.
 */
SCIP_CONSNONLINEAR_BILINTERM* SCIPgetBilinTermsNonlinear(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(consexprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->bilinterms;
}

/** returns the index of the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns -1 if the variables do not appear bilinearly.
 */
int SCIPgetBilinTermIdxNonlinear(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM entry;
   int idx;

   assert(consexprhdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->bilinhashtable == NULL )
   {
      return -1;
   }

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
   }
   assert(SCIPvarCompare(x, y) < 1);

   /* use a new entry to find the image in the bilinear hash table */
   entry.x = x;
   entry.y = y;
   idx = (int)(size_t)SCIPhashtableRetrieve(conshdlrdata->bilinhashtable, (void*)&entry) - 1;
   assert(idx >= -1 && idx < conshdlrdata->nbilinterms);
   assert(idx < 0 || conshdlrdata->bilinterms[idx].x == x);
   assert(idx < 0 || conshdlrdata->bilinterms[idx].y == y);

   return idx;
}

/** returns the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns NULL if the variables do not appear bilinearly.
 */
SCIP_CONSNONLINEAR_BILINTERM* SCIPgetBilinTermNonlinear(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int idx;

   assert(consexprhdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);

   conshdlrdata = SCIPconshdlrGetData(consexprhdlr);
   assert(conshdlrdata != NULL);

   idx = SCIPgetBilinTermIdxNonlinear(consexprhdlr, x, y);
   assert(idx >= -1 && idx < conshdlrdata->nbilinterms);

   if( idx >= 0 )
   {
      return &conshdlrdata->bilinterms[idx];
   }

   return NULL;
}

/** evaluates an auxiliary expression for a bilinear term */
SCIP_Real SCIPevalBilinAuxExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable of the bilinear term */
   SCIP_VAR*             y,                  /**< second variable of the bilinear term */
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr,           /**< auxiliary expression */
   SCIP_SOL*             sol                 /**< solution at which to evaluate (can be NULL) */
   )
{
   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(auxexpr != NULL);
   assert(auxexpr->auxvar != NULL);

   return auxexpr->cst + auxexpr->coefs[0] * SCIPgetSolVal(scip, sol, auxexpr->auxvar) +
          auxexpr->coefs[1] * SCIPgetSolVal(scip, sol, x) + auxexpr->coefs[2] * SCIPgetSolVal(scip, sol, y);
}

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_RETCODE SCIPinsertBilinearTermExistingNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   int                   nlockspos,          /**< number of positive expression locks */
   int                   nlocksneg           /**< number of negative expression locks */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM* term;
   int idx;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( bilinearTermsInsertEntry(scip, conshdlr, x, y, nlockspos, nlocksneg, &idx, TRUE) );

   term = &conshdlrdata->bilinterms[idx];
   assert(term != NULL);

   /* store the auxiliary variable */
   assert(term->nauxexprs == 0); /* existing terms should be added before implicit terms */
   term->aux.var = auxvar;

   /* capture auxiliary variable */
   if( auxvar != NULL )
   {
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   return SCIP_OKAY;
}

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_RETCODE SCIPinsertBilinearTermImplicitNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   SCIP_Real             coefx,              /**< coefficient of x in the auxiliary expression */
   SCIP_Real             coefy,              /**< coefficient of y in the auxiliary expression */
   SCIP_Real             coefaux,            /**< coefficient of the auxiliary variable in the auxiliary expression */
   SCIP_Real             cst,                /**< constant of the auxiliary expression */
   SCIP_Bool             overestimate        /**< whether the auxiliary expression overestimates the bilinear product */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM* term;
   int idx;
   int nlockspos;
   int nlocksneg;
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr;
   SCIP_Bool added;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nlockspos = overestimate ? 1 : 0;
   nlocksneg = overestimate ? 0 : 1;

   SCIP_CALL( bilinearTermsInsertEntry(scip, conshdlr, x, y, nlockspos, nlocksneg, &idx, FALSE) );

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
      SCIPswapReals(&coefx, &coefy);
   }
   assert(SCIPvarCompare(x, y) < 1);

   term = &conshdlrdata->bilinterms[idx];
   assert(term != NULL);
   assert(SCIPvarCompare(term->x, term->y) < 1);

   if( term->existing && term->nauxexprs == 0 && term->aux.var != NULL )
   {
      SCIP_CONSNONLINEAR_AUXEXPR* auxvarexpr;
      /* this is the case where we are adding an implicitly defined relation for a product that has already
       * been explicitly defined; convert auxvar into an auxexpr */

      /* nothing to do if we aren't allowed to add more than one auxexpr per term */
      if( conshdlrdata->bilinmaxnauxexprs <= 1 )
         return SCIP_OKAY;

      SCIP_CALL( SCIPallocBlockMemory(scip, &auxvarexpr) );
      auxvarexpr->cst = 0.0;
      auxvarexpr->coefs[0] = 1.0;
      auxvarexpr->coefs[1] = 0.0;
      auxvarexpr->coefs[2] = 0.0;
      auxvarexpr->auxvar = term->aux.var;
      auxvarexpr->underestimate = term->nlocksneg > 0;
      auxvarexpr->overestimate = term->nlockspos > 0;

      /* before we were working with term->aux.var; now aux.var has been saved and aux.exprs can be initialised to NULL */
      term->aux.exprs = NULL;

      SCIP_CALL( bilinTermAddAuxExpr(scip, conshdlrdata, term, auxvarexpr, &added) );

      /* since there were no auxexprs before and we've already checked for bilinmaxnauxexprs, auxvarexpr should always be added */
      assert(added);
   }

   /* create and add auxexpr */
   SCIP_CALL( SCIPallocBlockMemory(scip, &auxexpr) );
   auxexpr->underestimate = !overestimate;
   auxexpr->overestimate = overestimate;
   auxexpr->auxvar = auxvar;
   auxexpr->coefs[0] = coefaux;
   auxexpr->coefs[1] = coefx;
   auxexpr->coefs[2] = coefy;
   auxexpr->cst = cst;
   SCIP_CALL( bilinTermAddAuxExpr(scip, conshdlrdata, term, auxexpr, &added) );

   if( !added )
   {
      SCIPfreeBlockMemory(scip, &auxexpr);
   }
   else if( auxvar != NULL )
   { /* capture auxiliary variable */
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   return SCIP_OKAY;
}

/* returns the number of enforcements for an expression */
int SCIPgetExprNEnfosNonlinear(
   SCIP_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);
   return expr->nenfos;
}

/** returns the data for one of the enforcements of an expression */
void SCIPgetExprEnfoDataNonlinear(
   SCIP_EXPR*   expr,                         /**< expression */
   int                   idx,                          /**< position of enforcement in enfos array */
   SCIP_NLHDLR** nlhdlr,                      /**< buffer to store nlhldr */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,      /**< buffer to store nlhdlr data for expression, or NULL */
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD* nlhdlrparticipation, /**< buffer to store methods where nonlinear handler participates, or NULL */
   SCIP_Bool*            sepabelowusesactivity,        /**< buffer to store whether sepabelow uses activity of some expression, or NULL */
   SCIP_Bool*            sepaaboveusesactivity,        /**< buffer to store whether sepaabove uses activity of some expression, or NULL */
   SCIP_Real*            auxvalue                      /**< buffer to store current auxvalue, or NULL */
   )
{
   assert(expr != NULL);
   assert(idx >= 0);
   assert(idx < expr->nenfos);
   assert(expr->enfos[idx] != NULL);
   assert(nlhdlr != NULL);

   *nlhdlr = expr->enfos[idx]->nlhdlr;

   if( nlhdlrexprdata != NULL )
      *nlhdlrexprdata = expr->enfos[idx]->nlhdlrexprdata;

   if( nlhdlrparticipation != NULL )
      *nlhdlrparticipation = expr->enfos[idx]->nlhdlrparticipation;

   if( sepabelowusesactivity != NULL )
      *sepabelowusesactivity = expr->enfos[idx]->sepabelowusesactivity;

   if( sepaaboveusesactivity != NULL )
      *sepaaboveusesactivity = expr->enfos[idx]->sepaaboveusesactivity;

   if( auxvalue != NULL )
      *auxvalue = expr->enfos[idx]->auxvalue;
}

/** sets the auxiliary value of expression for one of the enforcements of an expression */
void SCIPsetExprEnfoAuxValueNonlinear(
   SCIP_EXPR*   expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_Real             auxvalue            /**< the new value of auxval */
)
{
   assert(expr != NULL);
   assert(idx >= 0);
   assert(idx < expr->nenfos);
   assert(expr->enfos[idx] != NULL);

   expr->enfos[idx]->auxvalue = auxvalue;
}

/** create and include conshdlr to SCIP and set everything except for expression handlers */
static
SCIP_RETCODE includeConshdlrExprBasic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create expr constraint handler data */
   SCIP_CALL( SCIPallocClearMemory(scip, &conshdlrdata) );
   conshdlrdata->intevalvar = intEvalVarBoundTightening;
   conshdlrdata->exprlastsoltag = 1;
   conshdlrdata->curboundstag = 1;
   conshdlrdata->lastboundrelax = 1;
   conshdlrdata->curpropboundstag = 1;
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->canonicalizetime) );
   SCIP_CALL( SCIPqueueCreate(&conshdlrdata->reversepropqueue, 100, 2.0) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->var2expr, SCIPblkmem(scip), 100) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyExpr,
         consFreeExpr, consInitExpr, consExitExpr,
         consInitpreExpr, consExitpreExpr, consInitsolExpr, consExitsolExpr,
         consDeleteExpr, consTransExpr, consInitlpExpr,
         consSepalpExpr, consSepasolExpr, consEnfolpExpr, consEnforelaxExpr, consEnfopsExpr, consCheckExpr,
         consPropExpr, consPresolExpr, consRespropExpr, consLockExpr,
         consActiveExpr, consDeactiveExpr,
         consEnableExpr, consDisableExpr, consDelvarsExpr,
         consPrintExpr, consCopyExpr, consParseExpr,
         consGetVarsExpr, consGetNVarsExpr, consGetDiveBdChgsExpr, conshdlrdata) );

   /* add expr constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a set of constraints within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/propauxvars",
         "whether to check bounds of all auxiliary variable to seed reverse propagation",
         &conshdlrdata->propauxvars, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/varboundrelax",
         "strategy on how to relax variable bounds during bound tightening: relax (n)ot, relax by (a)bsolute value, relax always by a(b)solute value, relax by (r)relative value",
         &conshdlrdata->varboundrelax, TRUE, 'r', "nabr", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/varboundrelaxamount",
         "by how much to relax variable bounds during bound tightening if strategy 'a', 'b', or 'r'",
         &conshdlrdata->varboundrelaxamount, TRUE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/conssiderelaxamount",
         "by how much to relax constraint sides during bound tightening",
         &conshdlrdata->conssiderelaxamount, TRUE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/vpmaxperturb",
         "maximal relative perturbation of reference point when computing facet of envelope of vertex-polyhedral function (dim>2)",
         &conshdlrdata->vp_maxperturb, TRUE, VERTEXPOLY_MAXPERTURBATION, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/vpadjfacetthresh",
         "adjust computed facet of envelope of vertex-polyhedral function up to a violation of this value times LP feasibility tolerance",
         &conshdlrdata->vp_adjfacetthreshold, TRUE, VERTEXPOLY_ADJUSTFACETFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/vpdualsimplex",
         "whether to use dual simplex instead of primal simplex for LP that computes facet of vertex-polyhedral function",
         &conshdlrdata->vp_dualsimplex, TRUE, VERTEXPOLY_USEDUALSIMPLEX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/bilinmaxnauxexprs",
           "maximal number of auxiliary expressions per bilinear term",
           &conshdlrdata->bilinmaxnauxexprs, FALSE, BILIN_MAXNAUXEXPRS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprods",
         "whether to reformulate products of binary variables during presolving",
         &conshdlrdata->reformbinprods, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsand",
         "whether to use the AND constraint handler for reformulating binary products",
         &conshdlrdata->reformbinprodsand, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsfac",
         "minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled)",
         &conshdlrdata->reformbinprodsfac, FALSE, 50, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forbidmultaggrnlvar",
         "whether to forbid multiaggregation of nonlinear variables",
         &conshdlrdata->forbidmultaggrnlvar, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/tightenlpfeastol",
         "whether to tighten LP feasibility tolerance during enforcement, if it seems useful",
         &conshdlrdata->tightenlpfeastol, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/propinenforce",
         "whether to (re)run propagation in enforcement",
         &conshdlrdata->propinenforce, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/weakcutthreshold",
         "threshold for when to regard a cut from an estimator as weak (lower values allow more weak cuts)",
         &conshdlrdata->weakcutthreshold, TRUE, 0.2, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/strongcutmaxcoef",
         "\"strong\" cuts will be scaled to have their maximal coef in [1/strongcutmaxcoef,strongcutmaxcoef]",
         &conshdlrdata->strongcutmaxcoef, TRUE, 1000.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/strongcutefficacy",
         "consider efficacy requirement when deciding whether a cut is \"strong\"",
         &conshdlrdata->strongcutefficacy, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forcestrongcut",
         "whether to force \"strong\" cuts in enforcement",
         &conshdlrdata->forcestrongcut, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/enfoauxviolfactor",
         "an expression will be enforced if the \"auxiliary\" violation is at least this factor times the \"original\" violation",
         &conshdlrdata->enfoauxviolfactor, TRUE, 0.01, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/weakcutminviolfactor",
         "retry enfo of constraint with weak cuts if violation is least this factor of maximal violated constraints",
         &conshdlrdata->weakcutminviolfactor, TRUE, 0.5, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/violscale",
         "method how to scale violations to make them comparable (not used for feasibility check): (n)one, (a)ctivity and side, norm of (g)radient",
         &conshdlrdata->violscale, TRUE, 'n', "nag", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/checkvarlocks",
         "whether variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction)",
         &conshdlrdata->checkvarlocks, TRUE, 't', "bdt", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/branching/aux",
         "from which depth on in the tree to allow branching on auxiliary variables (variables added for extended formulation)",
         &conshdlrdata->branchauxmindepth, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branching/external",
         "whether to use external branching candidates and branching rules for branching",
         &conshdlrdata->branchexternal, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/highviolfactor",
         "consider a constraint highly violated if its violation is >= this factor * maximal violation among all constraints",
         &conshdlrdata->branchhighviolfactor, FALSE, 0.0, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/highscorefactor",
         "consider a variable branching score high if its branching score >= this factor * maximal branching score among all variables",
         &conshdlrdata->branchhighscorefactor, FALSE, 0.9, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/violweight",
         "weight by how much to consider the violation assigned to a variable for its branching score",
         &conshdlrdata->branchviolweight, FALSE, 1.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/dualweight",
         "weight by how much to consider the dual values of rows that contain a variable for its branching score",
         &conshdlrdata->branchdualweight, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/pscostweight",
         "weight by how much to consider the pseudo cost of a variable for its branching score",
         &conshdlrdata->branchpscostweight, FALSE, 1.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/domainweight",
         "weight by how much to consider the domain width in branching score",
         &conshdlrdata->branchdomainweight, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/vartypeweight",
         "weight by how much to consider variable type (continuous: 0, binary: 1, integer: 0.1, impl-integer: 0.01) in branching score",
         &conshdlrdata->branchvartypeweight, FALSE, 0.5, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branching/scoreagg",
         "how to aggregate several branching scores given for the same expression: 'a'verage, 'm'aximum, 's'um",
         &conshdlrdata->branchscoreagg, TRUE, 's', "ams", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branching/violsplit",
         "method used to split violation in expression onto variables: 'e'venly, 'm'idness of solution, 'd'omain width, 'l'ogarithmic domain width",
         &conshdlrdata->branchviolsplit, TRUE, 'm', "emdl", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/pscostreliable",
         "minimum pseudo-cost update count required to consider pseudo-costs reliable",
         &conshdlrdata->branchpscostreliable, FALSE, 2.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   /* include handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, CONSHDLR_NAME "_boundchange",
         "signals a bound change to an expression constraint", processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   /* include table for statistics */
   assert(SCIPfindTable(scip, TABLE_NAME_EXPR) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_EXPR, TABLE_DESC_EXPR, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputExpr,
         NULL, TABLE_POSITION_EXPR, TABLE_EARLIEST_STAGE_EXPR) );

   return SCIP_OKAY;
}

/** creates the handler for expr constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( includeConshdlrExprBasic(scip) );

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* include and remember handler for variable expression */
   SCIP_CALL( SCIPincludeExprHdlrVar(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "var") == 0);
   conshdlrdata->exprvarhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for constant value expression */
   SCIP_CALL( SCIPincludeExprHdlrValue(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "val") == 0);
   conshdlrdata->exprvalhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for sum expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "sum") == 0);
   conshdlrdata->exprsumhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for product expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "prod") == 0);
   conshdlrdata->exprprodhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for exponential expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrExp(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "exp") == 0);
   conshdlrdata->exprexphdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for logarithmic expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrLog(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "log") == 0);
   conshdlrdata->exprloghdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for absolute expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrAbs(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "abs") == 0);

   /* include handler for power expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrPow(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "pow") == 0);
   conshdlrdata->exprpowhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for signed power expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSignpower(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "signpower") == 0);
   conshdlrdata->exprsignpowhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for entropy expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrEntropy(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "entropy") == 0);

   /* include handler for sine expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSin(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "sin") == 0);

   /* include handler for cosine expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrCos(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "cos") == 0);

   /* include handler for error function expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrErf(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "erf") == 0);

   /* include default nonlinear handler */
   SCIP_CALL( SCIPincludeConsExprNlhdlrDefault(scip, conshdlr) );

   /* include nonlinear handler for quadratics */
   SCIP_CALL( SCIPincludeConsExprNlhdlrQuadratic(scip, conshdlr) );

   /* include nonlinear handler for convex expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrConvex(scip, conshdlr) );

   /* include nonlinear handler for concave expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrConcave(scip, conshdlr) );

   /* include nonlinear handler for bilinear expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrBilinear(scip, conshdlr) );

   /* include nonlinear handler for SOC constraints */
   SCIP_CALL( SCIPincludeConsExprNlhdlrSoc(scip, conshdlr) );

   /* include nonlinear handler for use of perspective formulations */
   SCIP_CALL( SCIPincludeConsExprNlhdlrPerspective(scip, conshdlr) );

   /* include nonlinear handler for quotient expressions */
   SCIP_CALL( SCIPincludeConsExprNlhdlrQuotient(scip, conshdlr) );

   return SCIP_OKAY;
}

/** includes an expression constraint upgrade method into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsUpgradeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*exprconsupgd)),  /**< method to call for upgrading expression constraint, or NULL */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   CONSUPGRADE* exprconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;

   assert(conshdlrname != NULL );

   /* ignore empty upgrade functions */
   if( exprconsupgd == NULL )
      return SCIP_OKAY;

   /* find the expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether upgrade method exists already */
   for( i = conshdlrdata->nexprconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->exprconsupgrades[i]->exprconsupgd == exprconsupgd )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade method %p for constraint handler <%s>.\n", exprconsupgd, conshdlrname); /*lint !e611*/
#endif
         return SCIP_OKAY;
      }
   }

   /* create an expression constraint upgrade data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &exprconsupgrade) );
   exprconsupgrade->exprconsupgd = exprconsupgd;
   exprconsupgrade->priority   = priority;
   exprconsupgrade->active     = active;

   /* insert expression constraint upgrade method into constraint handler data */
   assert(conshdlrdata->nexprconsupgrades <= conshdlrdata->exprconsupgradessize);
   if( conshdlrdata->nexprconsupgrades+1 > conshdlrdata->exprconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nexprconsupgrades+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->exprconsupgrades, conshdlrdata->nexprconsupgrades, newsize) );
      conshdlrdata->exprconsupgradessize = newsize;
   }
   assert(conshdlrdata->nexprconsupgrades+1 <= conshdlrdata->exprconsupgradessize);

   for( i = conshdlrdata->nexprconsupgrades; i > 0 && conshdlrdata->exprconsupgrades[i-1]->priority < exprconsupgrade->priority; --i )
      conshdlrdata->exprconsupgrades[i] = conshdlrdata->exprconsupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nexprconsupgrades);
   conshdlrdata->exprconsupgrades[i] = exprconsupgrade;
   conshdlrdata->nexprconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/" CONSHDLR_NAME "/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable expression upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &exprconsupgrade->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a expr constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsNonlinear() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;

   assert(expr != NULL);

   /* find the expr constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("expr constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( createConsExpr(scip, conshdlr, cons, name, expr, lhs, rhs, TRUE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   return SCIP_OKAY;
}

/** creates and captures a expr constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name, expr, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic expression constraint */
SCIP_RETCODE SCIPcreateConsQuadraticNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* expr;

   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadterms == 0 || (quadvars1 != NULL && quadvars2 != NULL && quadcoefs != NULL));

   /* get expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* create quadratic expression */
   SCIP_CALL( SCIPcreateExprQuadratic(scip, conshdlr, &expr, nlinvars, linvars, lincoefs, nquadterms,
      quadvars1, quadvars2, quadcoefs) );
   assert(expr != NULL);

   /* create expression constraint */
   SCIP_CALL( createConsExpr(scip, conshdlr, cons, name, expr, lhs, rhs, TRUE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* release quadratic expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** returns the expression of the given expression constraint */
SCIP_EXPR* SCIPgetExprConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not expression\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->expr;
}

/** returns the root curvature of the given expression constraint
 *
 * @note The curvature information are computed during CONSINITSOL.
 */
SCIP_EXPRCURV SCIPgetCurvatureConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not expression\n");
      SCIPABORT();
      return SCIP_EXPRCURV_UNKNOWN;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->curv;
}

/** returns representation of the expression of the given expression constraint as quadratic form, if possible
 *
 * Only sets *isquadratic to TRUE if the whole expression is quadratic (in the non-extended formulation) and non-linear.
 * That is, the expr in each SCIP_QUADEXPR_QUADTERM will be a variable expressions and
 * \ref SCIPexprvarGetVar() can be used to retrieve the variable.
 */
SCIP_RETCODE SCIPcheckQuadraticConsNonlinear(
   SCIP*                    scip,               /**< SCIP data structure */
   SCIP_CONS*               cons,               /**< constraint data */
   SCIP_Bool*               isquadratic         /**< buffer to store whether constraint is quadratic */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(isquadratic != NULL);

   /* check whether constraint expression is quadratic in extended formulation */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, SCIPgetExprConsNonlinear(scip, cons), isquadratic) );

   /* if not quadratic in non-extended formulation, then do indicate quadratic */
   if( *isquadratic )
      *isquadratic = SCIPexprAreQuadraticExprsVariables(SCIPgetExprConsNonlinear(scip, cons));

   return SCIP_OKAY;
}

/** gets the expr constraint as a nonlinear row representation. */
SCIP_RETCODE SCIPgetNlRowConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** gets the left hand side of an expression constraint */
SCIP_Real SCIPgetLhsConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets the right hand side of an expression constraint */
SCIP_Real SCIPgetRhsConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** adds coef * var to expression constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPaddLinearTermConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             coef,               /**< coefficient */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPaddLinearTermConsNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   if( coef == 0.0 )
      return SCIP_OKAY;

   /* we should not have collected additional data for it
    * if some of these asserts fail, we may have to remove it and add some code to keep information uptodate
    */
   assert(consdata->nvarexprs == 0);
   assert(consdata->varexprs == NULL);
   assert(!consdata->catchedevents);

   SCIP_CALL( createExprVar(scip, conshdlr, &varexpr, var) );

   /* append to sum, if consdata->expr is sum and not used anywhere else */
   if( SCIPexprGetNUses(consdata->expr) == 1 && SCIPexprGetHdlr(consdata->expr) == SCIPgetExprHdlrSum(conshdlr) )
   {
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, consdata->expr, varexpr, coef) );
   }
   else
   {
      /* create new expression = 1 * consdata->expr + coef * var */
      SCIP_EXPR* children[2] = { consdata->expr, varexpr };
      SCIP_Real coefs[2] = { 1.0, coef };

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &consdata->expr, 2, children, coefs, 0.0) );

      /* release old root expr */
      SCIP_CALL( SCIPreleaseExpr(scip, &children[0]) );
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );

   /* not sure we care about any of these flags for original constraints */
   consdata->issimplified = FALSE;
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** gets absolute violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * If this value is at most SCIPfeastol(scip), the constraint would be considered feasible.
 */
SCIP_RETCODE SCIPgetAbsViolationConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   )
{
   assert(cons != NULL);
   assert(viol != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, 0) );
   *viol = getConsAbsViolation(cons);

   return SCIP_OKAY;
}

/** gets scaled violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * The scaling that is applied to the absolute violation of the constraint
 * depends on the setting of parameter constraints/expr/violscale.
 */
SCIP_RETCODE SCIPgetRelViolationConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   )
{
   assert(cons != NULL);
   assert(viol != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, 0) );
   SCIP_CALL( getConsRelViolation(scip, cons, viol, sol, 0) );

   return SCIP_OKAY;
}

/** gives the unique index of an expression constraint
 *
 * Each expression constraint gets an index assigned when it is created.
 * This index never changes and is unique among all expression constraints
 * within the same SCIP instance.
 * Thus, it can be used to sort a set of expression constraints.
 */
int SCIPgetIndexConsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->consindex;
}

/** compares two expression constraints by its index
 *
 * Usable as compare operator in array sort functions.
 */
int SCIPcompareIndexConsNonlinear(
   void*                 cons1,
   void*                 cons2
   )
{
   return SCIPgetIndexConsNonlinear((SCIP_CONS*)cons1) - SCIPgetIndexConsNonlinear((SCIP_CONS*)cons2);
}

/** returns an equivalent linear constraint if possible */
SCIP_RETCODE SCIPgetLinearConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_CONS**           lincons             /**< buffer to store linear constraint data */
   )
{
   SCIP_EXPRHDLR* sumhdlr;
   SCIP_EXPRHDLR* varhdlr;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   SCIP_VAR** vars;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(lincons != NULL);

   *lincons = NULL;
   expr = SCIPgetExprConsNonlinear(scip, cons);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   sumhdlr = SCIPgetExprHdlrSum(conshdlr);
   assert(sumhdlr != NULL);
   varhdlr = SCIPgetExprHdlrVar(conshdlr);
   assert(varhdlr != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* not a linear constraint if the root expression is not a sum */
   if( expr == NULL || expr->exprhdlr != sumhdlr )
      return SCIP_OKAY;

   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];

      /* at least one child is not a variable -> not a linear constraint */
      if( child->exprhdlr != varhdlr )
         return SCIP_OKAY;
   }

   /* collect all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPexprGetNChildren(expr)) );
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];

      assert(child->exprhdlr == varhdlr);
      vars[i] = SCIPexprvarGetVar(child);
   }

   /* consider constant part of the sum expression */
   lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : (consdata->lhs - SCIPexprsumGetConstant(expr));
   rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : (consdata->rhs - SCIPexprsumGetConstant(expr));

   SCIP_CALL( SCIPcreateConsLinear(scip, lincons, SCIPconsGetName(cons),
         SCIPexprGetNChildren(expr), vars, SCIPexprsumGetCoefs(expr),
         lhs, rhs,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

   /* free memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** returns a variable that appears linearly that may be decreased without making any other constraint infeasible */
SCIP_RETCODE SCIPgetLinvarMayDecreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check for a linear variable that can be increased or decreased without harming feasibility */
   findUnlockedLinearVar(scip, conshdlr, consdata);

   *var = consdata->linvardecr;
   *coef = consdata->linvardecrcoef;

   return SCIP_OKAY;
}

/** returns a variable that appears linearly that may be increased without making any other constraint infeasible */
SCIP_RETCODE SCIPgetLinvarMayIncreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check for a linear variable that can be increased or decreased without harming feasibility */
   findUnlockedLinearVar(scip, conshdlr, consdata);

   *var = consdata->linvarincr;
   *coef = consdata->linvarincrcoef;

   return SCIP_OKAY;
}

/** detects nonlinear handlers that can handle the expressions and creates needed auxiliary variables
 *
 *  @note this method is only used for testing purposes
 */
SCIP_RETCODE SCIPdetectNlhdlrsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   )
{
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( detectNlhdlrs(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** add the cut and maybe report branchscores */
SCIP_RETCODE SCIPprocessRowprepNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_NLHDLR* nlhdlr,             /**< nonlinear handler which provided the estimator */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_ROWPREP*         rowprep,            /**< cut to be added */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             allowweakcuts,      /**< whether we should only look for "strong" cuts, or anything that separates is fine */
   SCIP_Bool             branchscoresuccess, /**< whether the branching score callback of the estimator was successful */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, or only in separation */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result */
)
{
   SCIP_Real cutviol;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real auxvarvalue;
   SCIP_Bool sepasuccess;
   SCIP_Real estimateval = SCIP_INVALID;
   SCIP_Real mincutviolation;

   /* decide on minimal violation of cut */
   if( sol == NULL )
      mincutviolation = SCIPgetLPFeastol(scip);  /* we enforce an LP solution */
   else
      mincutviolation = SCIPfeastol(scip);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sepasuccess = TRUE;

   cutviol = SCIPgetRowprepViolation(scip, rowprep, sol, NULL);
   if( cutviol > 0.0 )
   {
      auxvarvalue = SCIPgetSolVal(scip, sol, auxvar);

      /* check whether cut is weak (if f(x) not defined, then it's never weak) */
      if( !allowweakcuts && auxvalue != SCIP_INVALID )  /*lint !e777*/
      {
         /* let the estimator be c'x-b, the auxvar is z (=auxvarvalue), and the expression is f(x) (=auxvalue)
          * then if we are underestimating and since the cut is violated, we should have z <= c'x-b <= f(x)
          * cutviol is c'x-b - z, so estimator value is c'x-b = z + cutviol
          * if the estimator value (c'x-b) is too close to z (auxvarvalue), when compared to f(x) (auxvalue),
          * then let's call this a weak cut that is, it's a weak cut if c'x-b <= z + weakcutthreshold * (f(x)-z)
          *   <->   c'x-b - z <= weakcutthreshold * (f(x)-z)
          *
          * if we are overestimating, we have z >= c'x-b >= f(x)
          * cutviol is z - (c'x-b), so estimator value is c'x-b = z - cutviol
          * it's weak if c'x-b >= f(x) + (1-weakcutthreshold) * (z - f(x))
          *   <->   c'x-b - z >= weakcutthreshold * (f(x)-z)
          *
          * when linearizing convex expressions, then we should have c'x-b = f(x), so they would never be weak
          */
         if( (!overestimate && ( cutviol <= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) ||
             ( overestimate && (-cutviol >= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded, but cut is too "\
                           "weak: auxvarvalue %g estimateval %g auxvalue %g (over %d)\n",
                                     SCIPnlhdlrGetName(nlhdlr), auxvarvalue,
                                     auxvarvalue + (overestimate ? -cutviol : cutviol), auxvalue, overestimate); )
            sepasuccess = FALSE;
         }
      }

      /* save estimator value for later, see long comment above why this gives the value for c'x-b */
      estimateval = auxvarvalue + (!overestimate ? cutviol : -cutviol);
   }
   else
   {
      sepasuccess = FALSE;
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded, but cut does not "\
                     "separate\n", SCIPnlhdlrGetName(nlhdlr)); )
   }

   /* clean up estimator */
   if( sepasuccess )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded: auxvarvalue %g "\
                     "estimateval %g auxvalue %g (over %d)\n    ", SCIPnlhdlrGetName(nlhdlr), auxvarvalue,
                               auxvarvalue + (overestimate ? -cutviol : cutviol), auxvalue, overestimate);
                       SCIPprintRowprep(scip, rowprep, enfologfile); )

      /* if not allowweakcuts, then do not attempt to get cuts more violated by scaling them up,
       * instead, may even scale them down, that is, scale so that max coef is close to 1
       */
      if( !allowweakcuts )
      {
         SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIP_CONSNONLINEAR_CUTMAXRANGE,
                                        conshdlrdata->strongcutmaxcoef, &sepasuccess) );

         if( !sepasuccess )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup cut failed due to bad numerics\n"); )
         }
         else
         {
            cutviol = SCIPgetRowprepViolation(scip, rowprep, sol, &sepasuccess);
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup succeeded, violation = %g and %sreliable, "\
                           "min requ viol = %g\n", cutviol, sepasuccess ? "" : "not ", mincutviolation); )
            if( sepasuccess )
               sepasuccess = cutviol > mincutviolation;
         }

         if( sepasuccess && auxvalue != SCIP_INVALID ) /*lint !e777*/
         {
            /* check whether cut is weak now
             * auxvar z may now have a coefficient due to scaling (down) in cleanup - take this into account when
             * reconstructing estimateval from cutviol (TODO improve or remove?)
             */
            SCIP_Real auxvarcoef = 0.0;
            int i;

            /* get absolute value of coef of auxvar in row - this makes the whole check here more expensive than
             * it should be...
             */
            for( i = 0; i < rowprep->nvars; ++i )
            {
               if( rowprep->vars[i] == auxvar )
               {
                  auxvarcoef = REALABS(rowprep->coefs[i]);
                  break;
               }
            }

            if( auxvarcoef == 0.0 ||
                (!overestimate && ( cutviol / auxvarcoef <= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) ||
                ( overestimate && (-cutviol / auxvarcoef >= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) )  /*lint !e644*/
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut is too weak after cleanup: auxvarvalue %g "\
                              "estimateval %g auxvalue %g (over %d)\n", auxvalue, auxvarvalue,
                                        auxvarvalue + (overestimate ? -cutviol : cutviol) / auxvarcoef, auxvalue, overestimate); )
               sepasuccess = FALSE;
            }
         }
      }
      else
      {
         /* TODO if violations are really tiny, then maybe handle special (decrease LP feastol, for example) */

         /* if estimate didn't report branchscores explicitly, then consider branching on those children for
          * which the following cleanup changes coefficients (we had/have this in cons_expr_sum this way)
          */
         if( !branchscoresuccess )
            rowprep->recordmodifications = TRUE;

         SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSNONLINEAR_CUTMAXRANGE, mincutviolation, &cutviol,
                                       &sepasuccess) );

         if( !sepasuccess )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup failed, %d coefs modified, cutviol %g\n",
                                     rowprep->nmodifiedvars, cutviol); )
         }

         /* if cleanup left us with a useless cut, then consider branching on variables for which coef were
          * changed
          */
         if( !sepasuccess && !branchscoresuccess && rowprep->nmodifiedvars > 0 )
         {
            SCIP_Real violscore;

#ifdef BRSCORE_ABSVIOL
            violscore = getExprAbsAuxViolation(scip, expr, auxvalue, sol, NULL, NULL);
#else
            SCIP_CALL( SCIPgetExprRelAuxViolationNonlinear(scip, conshdlr, expr, auxvalue, sol, &violscore, NULL,
                                                          NULL) );
#endif
            SCIP_CALL( addConsExprExprViolScoresAuxVars(scip, conshdlr, expr, violscore, rowprep->modifiedvars,
                                                        rowprep->nmodifiedvars, sol, &branchscoresuccess) );

            /* addConsExprExprBranchScoresAuxVars can fail if the only vars for which the coef was changed
             * - were fixed,
             * - are this expr's auxvar (I don't think it makes sense to branch on that one (would it?)), or
             * - if a variable in the rowprep is not in expr (can happen with indicator added by perspective)
             * the first case came up again in #3085 and I don't see how to exclude this in the assert,
             * so I'm disabling the assert for now
             */
            /* assert(branchscoresuccess || (rowprep->nmodifiedvars == 1 && rowprep->modifiedvars[0] == auxvar) ||
                  strcmp(SCIPnlhdlrGetName(nlhdlr), "perspective")==0); */
         }
      }
   }

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( sepasuccess )
   {
      SCIP_ROW* row;

      if( conshdlrdata->branchdualweight > 0.0 )
      {
         /* store remaining gap |f(x)-estimateval| in row name, which could be used in getDualBranchscore
          * skip if gap is zero
          */
         if( auxvalue == SCIP_INVALID )  /*lint !e777*/
            strcat(rowprep->name, "_estimategap=inf");
         else if( !SCIPisEQ(scip, auxvalue, estimateval) )
         {
            char gap[40];
            sprintf(gap, "_estimategap=%g", REALABS(auxvalue - estimateval));
            strcat(rowprep->name, gap);
         }
      }

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

      if( !allowweakcuts && conshdlrdata->strongcutefficacy && !SCIPisCutEfficacious(scip, sol, row) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut efficacy %g is too low (minefficacy=%g)\n",
                                  SCIPgetCutEfficacy(scip, sol, row), SCIPgetSepaMinEfficacy(scip)); )
      }
      else
      {
         SCIP_Bool infeasible;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    adding cut ");
                          SCIP_CALL( SCIPprintRow(scip, row, enfologfile) ); )

         /* I take !allowweakcuts as equivalent for having a strong cut (we usually have allowweakcuts=TRUE only
          * if we haven't found strong cuts before)
          */
         SCIP_CALL( SCIPaddRow(scip, row, conshdlrdata->forcestrongcut && !allowweakcuts && inenforcement,
                               &infeasible) );

#ifdef SCIP_CONSEXPR_ROWNOTREMOVABLE
         /* mark row as not removable from LP for current node, if in enforcement (this can prevent some cycling) */
         if( inenforcement )
            SCIPmarkRowNotRemovableLocal(scip, row);
#endif

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            ++nlhdlr->ncutoffs;
         }
         else
         {
            *result = SCIP_SEPARATED;
            ++nlhdlr->nseparated;
         }
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   else if( branchscoresuccess )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    separation with estimate of nlhdlr %s failed, but "\
                     "branching candidates added\n", SCIPnlhdlrGetName(nlhdlr)); )

      /* well, not branched, but addConsExprExprViolScoresAuxVars() added scores to (aux)variables and that makes the
       * expressions eligible for branching candidate, see enforceConstraints() and branching()
       */
      *result = SCIP_BRANCHED;
   }
   else
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    separation with estimate of nlhdlr %s failed and no "\
                     "branching candidates%s\n", SCIPnlhdlrGetName(nlhdlr), (allowweakcuts && inenforcement) ?
                                                                                    " (!)" : ""); )
   }

   return SCIP_OKAY;
}

/** creates the nonlinearity handler and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrNonlinear(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,           /**< expression constraint handler */
   SCIP_NLHDLR**      nlhdlr,             /**< buffer where to store nonlinear handler */
   const char*                 name,               /**< name of nonlinear handler (must not be NULL) */
   const char*                 desc,               /**< description of nonlinear handler (can be NULL) */
   int                         detectpriority,     /**< detection priority of nonlinear handler */
   int                         enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_NLHDLRDETECT((*detect)),  /**< structure detection callback of nonlinear handler */
   SCIP_DECL_NLHDLREVALAUX((*evalaux)), /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_NLHDLRDATA*   data                /**< data of nonlinear handler (can be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nlhdlr != NULL);
   assert(name != NULL);
   assert(detect != NULL);
   assert(evalaux != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, nlhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*nlhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*nlhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*nlhdlr)->detectpriority = detectpriority;
   (*nlhdlr)->enfopriority = enfopriority;
   (*nlhdlr)->data = data;
   (*nlhdlr)->detect = detect;
   (*nlhdlr)->evalaux = evalaux;

   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->detecttime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->enfotime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->proptime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->intevaltime) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/expr/nlhdlr/%s/enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname, "should this nonlinear handler be used",
      &(*nlhdlr)->enabled, FALSE, TRUE, NULL, NULL) );

   ENSUREBLOCKMEMORYARRAYSIZE(scip, conshdlrdata->nlhdlrs, conshdlrdata->nlhdlrssize, conshdlrdata->nnlhdlrs+1);

   conshdlrdata->nlhdlrs[conshdlrdata->nnlhdlrs] = *nlhdlr;
   ++conshdlrdata->nnlhdlrs;

   /* sort nonlinear handlers by detection priority, in decreasing order
    * will happen in INIT, so only do when called late
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_INIT && conshdlrdata->nnlhdlrs > 1 )
      SCIPsortDownPtr((void**)conshdlrdata->nlhdlrs, nlhdlrCmp, conshdlrdata->nnlhdlrs);

   return SCIP_OKAY;
}

/** set the copy handler callback of a nonlinear handler */
void SCIPnlhdlrSetCopyHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRCOPYHDLR((*copy)) /**< copy callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->copyhdlr = copy;
}

/** set the nonlinear handler callback to free the nonlinear handler data */
void SCIPnlhdlrSetFreeHdlrData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->freehdlrdata = freehdlrdata;
}

/** set the expression handler callback to free expression specific data of nonlinear handler */
void SCIPnlhdlrSetFreeExprData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback (can be NULL if data does not need to be freed) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->freeexprdata = freeexprdata;
}

/** set the initialization and deinitialization callback of a nonlinear handler */
void SCIPnlhdlrSetInitExit(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRINIT((*init)),   /**< initialization callback (can be NULL) */
   SCIP_DECL_NLHDLREXIT((*exit_))    /**< deinitialization callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->init = init;
   nlhdlr->exit = exit_;
}

/** set the propagation callbacks of a nonlinear handler */
void SCIPnlhdlrSetProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRINTEVAL((*inteval)), /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->inteval = inteval;
   nlhdlr->reverseprop = reverseprop;
}

/** set the separation callbacks of a nonlinear handler */
void SCIPnlhdlrSetSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_NLHDLRENFO((*enfo)),         /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_NLHDLRESTIMATE((*estimate)), /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_NLHDLREXITSEPA((*exitsepa))  /**< separation deinitialization callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);
   assert(enfo != NULL || estimate != NULL);

   nlhdlr->initsepa = initsepa;
   nlhdlr->enfo = enfo;
   nlhdlr->estimate = estimate;
   nlhdlr->exitsepa = exitsepa;
}

/** gives name of nonlinear handler */
const char* SCIPnlhdlrGetName(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->name;
}

/** gives description of nonlinear handler, can be NULL */
const char* SCIPnlhdlrGetDesc(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->desc;
}

/** gives detection priority of nonlinear handler */
int SCIPnlhdlrGetDetectPriority(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->detectpriority;
}

/** gives enforcement priority of nonlinear handler */
int SCIPnlhdlrGetEnfoPriority(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->enfopriority;
}

/** returns a nonlinear handler of a given name (or NULL if not found) */
SCIP_NLHDLR* SCIPfindNlhdlrNonlinear(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of nonlinear handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int h;

   assert(conshdlr != NULL);
   assert(name != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
      if( strcmp(SCIPnlhdlrGetName(conshdlrdata->nlhdlrs[h]), name) == 0 )
         return conshdlrdata->nlhdlrs[h];

   return NULL;
}

/** gives handler data of nonlinear handler */
SCIP_NLHDLRDATA* SCIPnlhdlrGetData(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->data;
}

/** gives nonlinear handler expression data
 *
 * @return NULL if expr has not been detected by nlhdlr or nlhdlr did not store data
 */
SCIP_NLHDLREXPRDATA* SCIPgetNlhdlrExprDataNonlinear(
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_EXPR*        expr           /**< expression */
   )
{
   int e;

   for( e = 0; e < expr->nenfos; ++e )
   {
      if( expr->enfos[e]->nlhdlr == nlhdlr )
         return expr->enfos[e]->nlhdlrexprdata;
   }

   return NULL;
}

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_Bool SCIPnlhdlrHasIntEval(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->inteval != NULL;
}

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_Bool SCIPnlhdlrHasReverseProp(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->reverseprop != NULL;
}

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_Bool SCIPnlhdlrHasInitSepa(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->initsepa != NULL;
}

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_Bool SCIPnlhdlrHasExitSepa(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->exitsepa != NULL;
}

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_Bool SCIPnlhdlrHasEnfo(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->enfo != NULL;
}

/** returns whether nonlinear handler implements the estimator callback */
SCIP_Bool SCIPnlhdlrHasEstimate(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->estimate != NULL;
}

/** call the detect callback of a nonlinear handler */
SCIP_DECL_NLHDLRDETECT(SCIPcallNlhdlrDetectNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->detect != NULL);
   assert(nlhdlr->detecttime != NULL);
   assert(participating != NULL);

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->detecttime) );
   SCIP_CALL( nlhdlr->detect(scip, conshdlr, nlhdlr, expr, cons, enforcing, participating, nlhdlrexprdata) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->detecttime) );

   if( *participating != SCIP_CONSNONLINEAR_EXPRENFO_NONE )
   {
      ++nlhdlr->ndetections;
      ++nlhdlr->ndetectionslast;
   }

   return SCIP_OKAY;
}

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLREVALAUX(SCIPcallNlhdlrEvalauxNonlinear)
{
   assert(nlhdlr != NULL);
   assert(nlhdlr->evalaux != NULL);

   SCIP_CALL( nlhdlr->evalaux(scip, nlhdlr, expr, nlhdlrexprdata, auxvalue, sol) );

   return SCIP_OKAY;
}

/** calls the interval evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLRINTEVAL(SCIPcallNlhdlrIntEvalNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->intevaltime != NULL);

   if( nlhdlr->inteval != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->intevaltime) );
      SCIP_CALL( nlhdlr->inteval(scip, nlhdlr, expr, nlhdlrexprdata, interval, intevalvar, intevalvardata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->intevaltime) );

      ++nlhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** calls the reverse propagation callback of a nonlinear handler */
SCIP_DECL_NLHDLRREVERSEPROP(SCIPcallNlhdlrReversePropNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->proptime != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   if( nlhdlr->reverseprop == NULL )
   {
      *infeasible = FALSE;
      *nreductions = 0;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->proptime) );
   SCIP_CALL( nlhdlr->reverseprop(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, bounds, infeasible, nreductions) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->proptime) );

   /* update statistics */
   nlhdlr->ndomreds += *nreductions;
   if( *infeasible )
      ++nlhdlr->ncutoffs;
   ++nlhdlr->npropcalls;

   return SCIP_OKAY;
}

/** calls the separation initialization callback of a nonlinear handler */
SCIP_DECL_NLHDLRINITSEPA(SCIPcallNlhdlrInitSepaNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(infeasible != NULL);

   if( nlhdlr->initsepa == NULL )
   {
      *infeasible = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->initsepa(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, overestimate, underestimate, infeasible) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   ++nlhdlr->nenfocalls;
   if( *infeasible )
      ++nlhdlr->ncutoffs;

   return SCIP_OKAY;
}

/** calls the separation deinitialization callback of a nonlinear handler */
SCIP_DECL_NLHDLREXITSEPA(SCIPcallNlhdlrExitSepaNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);

   if( nlhdlr->exitsepa != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
      SCIP_CALL( nlhdlr->exitsepa(scip, nlhdlr, expr, nlhdlrexprdata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );
   }

   return SCIP_OKAY;
}

/** calls the enforcement callback of a nonlinear handler */
SCIP_DECL_NLHDLRENFO(SCIPcallNlhdlrEnfoNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(result != NULL);

   if( nlhdlr->enfo == NULL )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPcallNlhdlrEvalauxNonlinear(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      assert(auxvalue == auxvaluetest);  /* we should get EXACTLY the same value from calling evalaux with the same solution as before */  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->enfo(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, allowweakcuts, separated, addbranchscores, result) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;
   switch( *result )
   {
      case SCIP_SEPARATED :
         ++nlhdlr->nseparated;
         break;
      case SCIP_BRANCHED:
         ++nlhdlr->nbranchscores;
         break;
      case SCIP_CUTOFF:
         ++nlhdlr->ncutoffs;
         break;
      case SCIP_REDUCEDDOM:
         ++nlhdlr->ndomreds;
         break;
      default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** calls the estimator callback of a nonlinear handler */
SCIP_DECL_NLHDLRESTIMATE(SCIPcallNlhdlrEstimateNonlinear)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(success != NULL);
   assert(addedbranchscores != NULL);

   if( nlhdlr->estimate == NULL )
   {
      *success = FALSE;
      *addedbranchscores = FALSE;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPcallNlhdlrEvalauxNonlinear(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      assert(auxvalue == auxvaluetest);  /* we should get EXACTLY the same value from calling evalaux with the same solution as before */  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->estimate(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, targetvalue,
         rowpreps, success, addbranchscores, addedbranchscores) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;

   return SCIP_OKAY;
}

/* Quadratic expression functions */


/** evaluates quadratic term in a solution w.r.t. auxiliary variables
 *
 * \note This assumes that for every expr used in the quadratic data, an auxiliary variable is available.
 */
SCIP_Real SCIPevalExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_QUADEXPR* quaddata,         /**< quadratic form */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   )
{
   SCIP_Real auxvalue;
   int i;

   assert(scip != NULL);
   assert(quaddata != NULL);

   auxvalue = quaddata->constant;

   for( i = 0; i < quaddata->nlinexprs; ++i ) /* linear exprs */
      auxvalue += quaddata->lincoefs[i] * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(quaddata->linexprs[i]));

   for( i = 0; i < quaddata->nquadexprs; ++i ) /* quadratic terms */
   {
      SCIP_QUADEXPR_QUADTERM* quadexprterm;
      SCIP_Real solval;

      quadexprterm = &quaddata->quadexprterms[i];
      assert(SCIPgetExprAuxVarNonlinear(quadexprterm->expr) != NULL);

      solval = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(quadexprterm->expr));
      auxvalue += (quadexprterm->lincoef + quadexprterm->sqrcoef * solval) * solval;
   }

   for( i = 0; i < quaddata->nbilinexprterms; ++i ) /* bilinear terms */
   {
      SCIP_QUADEXPR_BILINTERM* bilinexprterm;

      bilinexprterm = &quaddata->bilinexprterms[i];
      assert(SCIPgetExprAuxVarNonlinear(bilinexprterm->expr1) != NULL);
      assert(SCIPgetExprAuxVarNonlinear(bilinexprterm->expr2) != NULL);
      auxvalue += bilinexprterm->coef *
         SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(bilinexprterm->expr1)) *
         SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(bilinexprterm->expr2));
   }

   return auxvalue;
}

/* computes a facet of the convex or concave envelope of a vertex polyhedral function
 * see (doxygen-)comment of this function in cons_expr.h
 * (this is by intention not a doxygen comment)
 */
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedralNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_DECL_VERTEXPOLYFUN((*function)),     /**< pointer to vertex polyhedral function */
   void*                 fundata,            /**< data for function evaluation (can be NULL) */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
)
{
   SCIP_Real* corner;
   SCIP_Real* funvals;
   int* nonfixedpos;
   SCIP_Real maxfaceterror;
   int nvars; /* number of nonfixed variables */
   unsigned int ncorners;
   unsigned int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(function != NULL);
   assert(xstar != NULL);
   assert(box != NULL);
   assert(success != NULL);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   /* identify fixed variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &nonfixedpos, nallvars) );
   nvars = 0;
   for( j = 0; j < nallvars; ++j )
   {
      if( SCIPisRelEQ(scip, box[2 * j], box[2 * j + 1]) )
         continue;
      nonfixedpos[nvars] = j;
      nvars++;
   }

   /* if all variables are fixed, then we could provide something trivial, but that wouldn't be the job of separation
    * if too many variables are not fixed, then we do nothing currently
    */
   if( nvars == 0 || nvars > SCIP_MAXVERTEXPOLYDIM )
   {
      SCIPwarningMessage(scip, "SCIPcomputeFacetVertexPolyhedralNonlinear() called with %d nonfixed variables. Must be between [1,%d].\n", nvars, SCIP_MAXVERTEXPOLYDIM);
      SCIPfreeBufferArray(scip, &nonfixedpos);
      return SCIP_OKAY;
   }

   /* compute f(v^i) for each corner v^i of [l,u] */
   ncorners = POWEROFTWO(nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &funvals, ncorners) );
   SCIP_CALL( SCIPallocBufferArray(scip, &corner, nallvars) );
   for( j = 0; j < nallvars; ++j )
   {
      if( SCIPisRelEQ(scip, box[2 * j], box[2 * j + 1]) )
         corner[j] = (box[2 * j] + box[2 * j + 1]) / 2.0;
   }
   for( i = 0; i < ncorners; ++i )
   {
      SCIPdebugMsg(scip, "corner %d: ", i);
      for( j = 0; j < nvars; ++j )
      {
         int varpos = nonfixedpos[j];
         /* if j'th bit of row index i is set, then take upper bound on var j, otherwise lower bound var j
          * we check this by shifting i for j positions to the right and checking whether the last bit is set
          */
         if( (i >> j) & 0x1 )
            corner[varpos] = box[2 * varpos + 1]; /* ub of var */
         else
            corner[varpos] = box[2 * varpos ]; /* lb of var */
         SCIPdebugMsgPrint(scip, "%g, ", corner[varpos]);
         assert(!SCIPisInfinity(scip, REALABS(corner[varpos])));
      }

      funvals[i] = function(corner, nallvars, fundata);

      SCIPdebugMsgPrint(scip, "obj = %e\n", funvals[i]);

      if( funvals[i] == SCIP_INVALID || SCIPisInfinity(scip, REALABS(funvals[i])) )  /*lint !e777*/
      {
         SCIPdebugMsg(scip, "cannot compute underestimator; function value at corner is too large %g\n", funvals[i]);
         goto CLEANUP;
      }
   }

   /* clear coefs array; below we only fill in coefs for nonfixed variables */
   BMSclearMemoryArray(facetcoefs, nallvars);

   if( nvars == 1 )
   {
      SCIP_CALL( computeVertexPolyhedralFacetUnivariate(scip, box[2 * nonfixedpos[0]], box[2 * nonfixedpos[0] + 1], funvals[0], funvals[1], success, &facetcoefs[nonfixedpos[0]], facetconstant) );

      /* check whether target has been missed */
      if( *success && overestimate == (*facetconstant + facetcoefs[nonfixedpos[0]] * xstar[nonfixedpos[0]] > targetvalue) )
      {
         SCIPdebugMsg(scip, "computed secant, but missed target %g (facetvalue=%g, overestimate=%d)\n", targetvalue, *facetconstant + facetcoefs[nonfixedpos[0]] * xstar[nonfixedpos[0]], overestimate);
         *success = FALSE;
      }
   }
   else if( nvars == 2 )
   {
      int idx1 = nonfixedpos[0];
      int idx2 = nonfixedpos[1];
      double p1[2] = { box[2*idx1],   box[2*idx2]   }; /* corner 0: 0>>0 & 0x1 = 0, 0>>1 & 0x1 = 0 */
      double p2[2] = { box[2*idx1+1], box[2*idx2]   }; /* corner 1: 1>>0 & 0x1 = 1, 1>>1 & 0x1 = 0 */
      double p3[2] = { box[2*idx1],   box[2*idx2+1] }; /* corner 2: 2>>0 & 0x1 = 0, 2>>1 & 0x1 = 1 */
      double p4[2] = { box[2*idx1+1], box[2*idx2+1] }; /* corner 3: 3>>0 & 0x1 = 1, 3>>1 & 0x1 = 1 */
      double xstar2[2] = { xstar[idx1], xstar[idx2] };
      double coefs[2] = { 0.0, 0.0 };

      SCIP_CALL( computeVertexPolyhedralFacetBivariate(scip, overestimate, p1, p2, p3, p4, funvals[0], funvals[1], funvals[2], funvals[3], xstar2, targetvalue, success, coefs, facetconstant) );

      facetcoefs[idx1] = coefs[0];
      facetcoefs[idx2] = coefs[1];
   }
   else
   {
      SCIP_CALL( computeVertexPolyhedralFacetLP(scip, conshdlr, overestimate, xstar, box, nallvars, nonfixedpos, funvals, nvars, targetvalue, success, facetcoefs, facetconstant) );
   }
   if( !*success )
   {
      SCIPdebugMsg(scip, "no success computing facet, %d vars\n", nvars);
      goto CLEANUP;
   }

   /*
    *  check and adjust facet with the algorithm of Rikun et al.
    */

   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, overestimate, funvals, box, nallvars, nvars, nonfixedpos, facetcoefs, *facetconstant);

   /* adjust constant part of the facet by maxerror to make it a valid over/underestimator (not facet though) */
   if( maxfaceterror > 0.0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real midval;
      SCIP_Real feastol;

      feastol = SCIPgetStage(scip) == SCIP_STAGE_SOLVING ? SCIPgetLPFeastol(scip) : SCIPfeastol(scip);

      /* evaluate function in middle point to get some idea for a scaling */
      for( j = 0; j < nvars; ++j )
         corner[nonfixedpos[j]] = (box[2 * nonfixedpos[j]] + box[2 * nonfixedpos[j] + 1]) / 2.0;
      midval = function(corner, nallvars, fundata);
      if( midval == SCIP_INVALID )  /*lint !e777*/
         midval = 1.0;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* there seem to be numerical problems if the error is too large; in this case we reject the facet */
      if( maxfaceterror > conshdlrdata->vp_adjfacetthreshold * feastol * fabs(midval) )
      {
         SCIPdebugMsg(scip, "ignoring facet due to instability, it cuts off a vertex by %g (midval=%g).\n", maxfaceterror, midval);
         *success = FALSE;
         goto CLEANUP;
      }

      SCIPdebugMsg(scip, "maximum facet error %g (midval=%g), adjust constant to make cut valid!\n", maxfaceterror, midval);

      if( overestimate )
         *facetconstant += maxfaceterror;
      else
         *facetconstant -= maxfaceterror;
   }

   /* if we made it until here, then we have a nice facet */
   assert(*success);

CLEANUP:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &funvals);
   SCIPfreeBufferArray(scip, &corner);
   SCIPfreeBufferArray(scip, &nonfixedpos);

   return SCIP_OKAY;
}

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
SCIP_RETCODE SCIPcomputeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   )
{
   assert(scip != NULL);
   assert(alpha != NULL);
   assert(beta  != NULL);
   assert(gamma_ != NULL);
   assert(delta != NULL);

   *alpha  = -b3*c2 + a3*(-b2+c2) + a2*(b3-c3) + b2*c3;
   *beta   = -(-b3*c1 + a3*(-b1+c1) + a1*(b3-c3) + b1*c3);
   *gamma_ = -a2*b1 + a1*b2 + a2*c1 - b2*c1 - a1*c2 + b1*c2;
   *delta  = -a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3;

   /* SCIPdebugMsg(scip, "alpha: %g beta: %g gamma: %g delta: %g\n", *alpha, *beta, *gamma_, *delta); */

   if( SCIPisInfinity(scip, REALABS(*gamma_ * a3)) ||
      SCIPisInfinity(scip, REALABS(*gamma_ * b3)) ||
      SCIPisInfinity(scip, REALABS(*gamma_ * c3)) )
   {
      SCIPdebugMsg(scip, "activity above SCIP infinity\n");
      *delta  = 0.0;
      *alpha  = 0.0;
      *beta   = 0.0;
      *gamma_ = 0.0;
      return SCIP_OKAY;
   }

   /* check if hyperplane contains all three points (necessary because of numerical troubles) */
   if( !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3) ||
      !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3) ||
      !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3) )
   {
      SCIP_Real m[9];
      SCIP_Real rhs[3];
      SCIP_Real x[3];
      SCIP_Bool success;

      /*
      SCIPdebugMsg(scip, "a = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", a1, a2, a3, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3, SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3));
      SCIPdebugMsg(scip, "b = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", b1, b2, b3, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3, SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3));
      SCIPdebugMsg(scip, "c = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", c1, c2, c3, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3, SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3));
      */

      /* initialize matrix column-wise */
      m[0] = a1;
      m[1] = b1;
      m[2] = c1;
      m[3] = a2;
      m[4] = b2;
      m[5] = c2;
      m[6] = a3;
      m[7] = b3;
      m[8] = c3;

      rhs[0] = 1.0;
      rhs[1] = 1.0;
      rhs[2] = 1.0;

      SCIPdebugMsg(scip, "numerical troubles - try to solve the linear system via an LU factorization\n");

      /* solve the linear problem */
      SCIP_CALL( SCIPsolveLinearProb(3, m, rhs, x, &success) );
      /* assert(success); */

      *delta  = rhs[0];
      *alpha  = x[0];
      *beta   = x[1];
      *gamma_ = x[2];

      /* set all coefficients to zero if one of the points is not contained in the hyperplane; this ensures that we do
       * not add a cut to SCIP and that all assertions are trivially fulfilled
       */
      if( !success || !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3) ||
         !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3) ||
         !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3) ) /*lint !e774*/
      {
         SCIPdebugMsg(scip, "could not resolve numerical difficulties\n");
         *delta  = 0.0;
         *alpha  = 0.0;
         *beta   = 0.0;
         *gamma_ = 0.0;
      }
   }

   if( *gamma_ < 0.0 )
   {
      *alpha  = -*alpha;
      *beta   = -*beta;
      *gamma_ = -*gamma_;
      *delta  = -*delta;
   }

   return SCIP_OKAY;
}
