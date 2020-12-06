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
