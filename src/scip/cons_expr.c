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

/*
 * Local methods
 */


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
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_CONST, SCIPgetValueExprValue(sourceexpr)) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "sum") == 0 )
   {
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), targetexpr, nchildren, children, SCIPgetCoefsExprSum(sourceexpr), SCIPgetConstantExprSum(sourceexpr)) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "pow") == 0 )
   {
      SCIP_Real exponent;

      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);

      exponent = SCIPgetExponentExprPow(sourceexpr);
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
         SCIPgetExponentExprPow(sourceexpr)) );
   }
   else if( strcmp(SCIPexprhdlrGetName(exprhdlr), "prod") == 0 )
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, SCIPgetCoefExprProduct(sourceexpr), nchildren, NULL, NULL) );
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
      SCIP_CALL( SCIPchgNlRowConstant(scip, consdata->nlrow, SCIPgetConstantExprSum(consdata->expr)) );

      coefs = SCIPgetCoefsExprSum(consdata->expr);

      for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
      {
         child = SCIPexprGetChildren(consdata->expr)[i];

         if( SCIPisExprVar(child) )
         {
            /* linear term */
            SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(child), coefs[i]) );
         }
         else if( SCIPexprGetHdlr(child) == conshdlrdata->exprpowhdlr &&
            SCIPgetExponentExprPow(child) == 2.0 &&
            SCIPisExprVar(SCIPexprGetChildren(child)[0]) )
         {
            /* square term  */
            SCIP_QUADELEM quadelem;

            quadelem.idx1 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(child)[0]));
            if( quadelem.idx1 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(child)[0])) );
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

            quadelem.idx1 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(child)[0]));
            if( quadelem.idx1 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(child)[0])) );
               quadelem.idx1 = SCIPnlrowGetNQuadVars(consdata->nlrow)-1;
            }

            quadelem.idx2 = SCIPnlrowSearchQuadVar(consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(child)[1]));
            if( quadelem.idx2 == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(child)[1])) );
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
               nlvars[nnlvars] = SCIPgetVarExprVar(consdata->varexprs[i]);
               ++nnlvars;
            }
         }

         SCIPexprReindexVars(classicexpr, reindexvars);

         SCIPfreeBufferArray(scip, &varsusage);
         SCIPfreeBufferArray(scip, &reindexvars);
      }
   }
   else if( SCIPexprGetHdlr(consdata->expr) == conshdlrdata->exprpowhdlr &&
      SCIPgetExponentExprPow(consdata->expr) == 2.0 &&
      SCIPisExprVar(SCIPexprGetChildren(consdata->expr)[0]) )
   {
      /* if root is a x^2, then set the quadratic part of the nlrow */
      SCIP_QUADELEM quadelem;

      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(consdata->expr)[0])) );
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

      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(consdata->expr)[0])) );
      SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(SCIPexprGetChildren(consdata->expr)[1])) );

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
         nlvars[i] = SCIPgetVarExprVar(consdata->varexprs[i]);
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
            if( (expr->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
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
            SCIP_CALL( SCIPnlhdlrInitsepa(scip, conshdlr, conss[c], nlhdlr, expr,
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
         var = SCIPgetVarExprVar(consdata->varexprs[i]);
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

            var = SCIPgetVarExprVar(consdata->varexprs[i]);
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

               var = SCIPgetVarExprVar(consdata->varexprs[i]);
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
   SCIP_CALL( SCIPnlhdlrEnfo(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
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

      SCIP_CALL( SCIPnlhdlrEstimate(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
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
      if( (expr->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
         continue;

      nlhdlr = expr->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &expr->enfos[e]->auxvalue, sol) );
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
      if( overestimate && auxoverestimate && (expr->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0 )  /*lint !e777*/
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
      if( underestimate && auxunderestimate && (expr->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABELOW) != 0 )  /*lint !e777*/
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
            var = SCIPgetVarExprVar(consdata->varexprs[i]);
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
               var = SCIPgetVarExprVar(expr);
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
            if( (expr->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
               continue;

            nlhdlr = expr->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
            SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, expr->enfos[e]->nlhdlrexprdata, &expr->enfos[e]->auxvalue, sol) );

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
         if( SCIPexprGetHdlr(expr) == powhdlr && SCIPgetExponentExprPow(expr) == 2.0 )
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

/*
 * Callback methods of constraint handler
 */

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

/** creates the handler for expr constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;


   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

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
