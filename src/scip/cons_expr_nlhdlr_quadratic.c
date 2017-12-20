/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_quadratic.c
 * @brief  nonlinear handler to handle quadratic expressions
 * @author Felipe Serrano
 *
 * Some definitions:
 * - a BILINEXPRTERM is a product of two expressions
 * - a SCIP_QUADEXPRTERM stores an expression expr that is known to appear in a nonlinear, quadratic term, that is
 *   expr^2 or expr * other_expr. It stores its sqrcoef (that can be 0), its linear coef and all the bilinear expression
 *   terms in which expr appears.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#define SCIP_PRIVATE_ROWPREP

#include "scip/cons_quadratic.h" /* for SCIP_ROWPREP */
#include "scip/cons_expr_nlhdlr_quadratic.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_product.h"

#include "nlpi/nlpi_ipopt.h" /* for LAPACK */

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "quadratic"
#define NLHDLR_DESC         "handler for quadratic expressions"
#define NLHDLR_PRIORITY     100

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   int                   nlinexprs;          /**< number of expressions that appear linearly */
   int                   linexprssize;       /**< size of linexprs and lincoefs arrays */
   SCIP_CONSEXPR_EXPR**  linexprs;           /**< expressions that appear linearly */
   SCIP_Real*            lincoefs;           /**< coefficients of expressions that appear linearly */

   int                   nquadexprs;         /**< number of expressions in quadratic terms */
   int                   quadexprssize;      /**< size of quadexprterms array */
   SCIP_QUADEXPRTERM*    quadexprterms;      /**< array with quadratic expression terms */

   int                   nbilinexprterms;    /**< number of bilinear expressions terms */
   int                   bilinexprtermssize; /**< size of bilinexprterms array */
   SCIP_BILINEXPRTERM*   bilinexprterms;     /**< bilinear expression terms array */

   SCIP_EXPRCURV         curvature;          /**< curvature of the quadratic representation of the expression */
};

/*
 * static methods
 */


/** frees nlhdlrexprdata structure */
static
void freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata /**< nlhdlr expression data */
   )
{
   int i;

   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->linexprs), nlhdlrexprdata->linexprssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->lincoefs), nlhdlrexprdata->linexprssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->bilinexprterms), nlhdlrexprdata->bilinexprtermssize);

   for( i = 0; i < nlhdlrexprdata->nquadexprs; ++i )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->quadexprterms[i].adjbilin),
            nlhdlrexprdata->quadexprterms[i].adjbilinsize);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->quadexprterms), nlhdlrexprdata->quadexprssize);
}

/** ensures, that linear vars and coefs arrays can store at least num entries */
static
SCIP_RETCODE nlhdlrexprdataEnsureLinearVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlinexprs <= nlhdlrexprdata->linexprssize);

   if( num > nlhdlrexprdata->linexprssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->linexprs,  nlhdlrexprdata->linexprssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->lincoefs, nlhdlrexprdata->linexprssize, newsize) );
      nlhdlrexprdata->linexprssize = newsize;
   }
   assert(num <= nlhdlrexprdata->linexprssize);

   return SCIP_OKAY;
}

/** ensures, that quadratic variable terms array can store at least num entries */
static
SCIP_RETCODE nlhdlrexprdataEnsureQuadVarTermsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nquadexprs <= nlhdlrexprdata->quadexprssize);

   if( num > nlhdlrexprdata->quadexprssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->quadexprterms, nlhdlrexprdata->quadexprssize, newsize) );
      nlhdlrexprdata->quadexprssize = newsize;
   }
   assert(num <= nlhdlrexprdata->quadexprssize);

   return SCIP_OKAY;
}

/** ensures, that adjacency array can store at least num entries */
static
SCIP_RETCODE nlhdlrexprdataEnsureAdjBilinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_QUADEXPRTERM*    quadexprterm,       /**< quadratic expression term */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(quadexprterm != NULL);
   assert(quadexprterm->nadjbilin <= quadexprterm->adjbilinsize);

   if( num > quadexprterm->adjbilinsize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &quadexprterm->adjbilin, quadexprterm->adjbilinsize, newsize) );
      quadexprterm->adjbilinsize = newsize;
   }
   assert(num <= quadexprterm->adjbilinsize);

   return SCIP_OKAY;
}

/** ensures, that bilinear term arrays can store at least num entries */
static
SCIP_RETCODE nlhdlrexprdataEnsureBilinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nbilinexprterms <= nlhdlrexprdata->bilinexprtermssize);

   if( num > nlhdlrexprdata->bilinexprtermssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->bilinexprterms, nlhdlrexprdata->bilinexprtermssize, newsize) );
      nlhdlrexprdata->bilinexprtermssize = newsize;
   }
   assert(num <= nlhdlrexprdata->bilinexprtermssize);

   return SCIP_OKAY;
}

/** if expr is in seenexpr then sets seen to TRUE, otherwise, inserts expr to seenexpr */
/* TODO: better name for this function? or nice way to factorize the code? */
static
SCIP_RETCODE checkProperQuadratic(
   SCIP_CONSEXPR_EXPR*   expr,               /**< the expression */
   SCIP_HASHMAP*         seenexpr,           /**< hash map */
   SCIP_Bool*            seen                /**< buffer to store whether expr is in hashmap */
   )
{
   if( SCIPhashmapExists(seenexpr, (void *)expr) )
      *seen = TRUE;
   else
   {
      SCIP_CALL( SCIPhashmapInsert(seenexpr, (void *)expr, NULL) );
   }

   return SCIP_OKAY;
}

/** Checks the curvature of the quadratic function, x^T Q x + b^T x stored in nlhdlrexprdata; for this, it builds the
 * matrix Q and computes its eigenvalues using LAPACK; if Q is
 * - semidefinite positive -> provided is set to sepaunder
 * - semidefinite negative -> provided is set to sepaover
 * - otherwise -> provided is set to none
 */
/* TODO: make more simple test; like diagonal entries don't change sign, etc */
static
SCIP_RETCODE checkCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata /**< nlhdlr expression data */
   )
{
   SCIP_HASHMAP* expr2matrix;
   double* matrix;
   double* alleigval;
   int nvars;
   int nn;
   int n;
   int i;

   nlhdlrexprdata->curvature = SCIP_EXPRCURV_UNKNOWN;

   n  = nlhdlrexprdata->nquadexprs;
   nn = n * n;

   /* do not check curvature if nn is too large */
   if( nn < 0 || (unsigned) (int) nn > UINT_MAX / sizeof(SCIP_Real) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "nlhdr_quadratic - number of quadratic variables is too large (%d) to check the curvature; will not handle this expression\n", n);

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &alleigval, n) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &matrix, nn) );

   SCIP_CALL( SCIPhashmapCreate(&expr2matrix, SCIPblkmem(scip), n) );

   /* fill matrix's diagonal */
   nvars = 0;
   for( i = 0; i < n; ++i )
   {
      SCIP_QUADEXPRTERM quadexprterm;

      quadexprterm = nlhdlrexprdata->quadexprterms[i];

      assert(!SCIPhashmapExists(expr2matrix, (void*)quadexprterm.expr));

      if( quadexprterm.sqrcoef == 0.0 )
      {
         assert(quadexprterm.nadjbilin > 0);
         /* SCIPdebugMsg(scip, "var <%s> appears in bilinear term but is not squared --> indefinite quadratic\n", SCIPvarGetName(quadexprterm.var)); */
         goto CLEANUP;
      }

      matrix[nvars * n + nvars] = quadexprterm.sqrcoef;

      /* remember row of variable in matrix */
      SCIP_CALL( SCIPhashmapInsert(expr2matrix, (void *)quadexprterm.expr, (void *)(size_t)nvars) );
      nvars++;
   }

   /* fill matrix's upper-diagonal */
   for( i = 0; i < nlhdlrexprdata->nbilinexprterms; ++i )
   {
      SCIP_BILINEXPRTERM bilinexprterm;
      int col;
      int row;

      bilinexprterm = nlhdlrexprdata->bilinexprterms[i];

      assert(SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr1));
      assert(SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr2));
      row = (int)(size_t)SCIPhashmapGetImage(expr2matrix, bilinexprterm.expr1);
      col = (int)(size_t)SCIPhashmapGetImage(expr2matrix, bilinexprterm.expr2);

      assert(row != col);

      if( row < col )
         matrix[row * n + col] = bilinexprterm.coef / 2.0;
      else
         matrix[col * n + row] = bilinexprterm.coef / 2.0;
   }

   /* compute eigenvalues */
   if( LapackDsyev(FALSE, n, matrix, alleigval) != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Failed to compute eigenvalues of quadratic coefficient matrix --> don't know curvature\n");
      goto CLEANUP;
   }

   /* check convexity */
   if( !SCIPisNegative(scip, alleigval[0]) )
   {
      nlhdlrexprdata->curvature = SCIP_EXPRCURV_CONVEX;
   }
   else if( !SCIPisPositive(scip, alleigval[n-1]) )
   {
      nlhdlrexprdata->curvature = SCIP_EXPRCURV_CONCAVE;
   }

CLEANUP:
   SCIPhashmapFree(&expr2matrix);
   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &alleigval);

   return SCIP_OKAY;
}


/** add expression expr to quadratic terms: this means several things depending on what is known about expr
 * - if is the first time seeing this expr -> creates new quadratic term
 * - if it has been seen linearly before -> removes it from the linear exprs and creates a new quadratic term
 * - if it has been seen quadratically before, then
 *    - if is the first time seeing the expr quadratically (i.e sqrcoef * expr^2) -> add sqrcoef to existing quad term
 *    - expr is being seen in a bilinear term (expr * other_expr) -> add bilinear information to quadexprterm
 */
static
SCIP_RETCODE addExprToQuadexprterms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to add to quadratic terms */
   SCIP_Real             sqrcoef,            /**< coefficient of variable in quadartic term */
   int                   bilinexprtermidx,   /**< index of bilin term where expr was seen: -1 if not coming from bilin term */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data (where the quadratic terms live) */
   SCIP_HASHMAP*         expridx             /**< map between expressoins and its position in linear vars (if positive) or
                                              * quadratic vars (if negative) */
   )
{
   SCIP_QUADEXPRTERM* quadexprterm;

   assert(expr != NULL);

   if( SCIPhashmapExists(expridx, expr) ) /* expr has been seen before */
   {
      int idx;

      idx = (int)(size_t)SCIPhashmapGetImage(expridx, expr);

      if( idx >= 0 ) /* expr has been seen before, but only as a linearly -> new quadexprterm */
      {
         assert(expr == nlhdlrexprdata->linexprs[idx]);

         SCIP_CALL( nlhdlrexprdataEnsureQuadVarTermsSize(scip, nlhdlrexprdata, nlhdlrexprdata->nquadexprs + 1) );

         quadexprterm = &nlhdlrexprdata->quadexprterms[nlhdlrexprdata->nquadexprs];

         quadexprterm->expr = expr;
         quadexprterm->sqrcoef = sqrcoef;
         quadexprterm->lincoef = nlhdlrexprdata->lincoefs[idx];
         quadexprterm->nadjbilin = 0;
         quadexprterm->adjbilinsize = 0;
         quadexprterm->adjbilin = NULL;

         /* expr now appears quadratically --> remove it from nlhdlrexprdata->linexprs */
         nlhdlrexprdata->nlinexprs--;
         if( idx < nlhdlrexprdata->nlinexprs )
         {
            nlhdlrexprdata->linexprs[idx] = nlhdlrexprdata->linexprs[nlhdlrexprdata->nlinexprs];
            nlhdlrexprdata->lincoefs[idx] = nlhdlrexprdata->lincoefs[nlhdlrexprdata->nlinexprs];
            SCIP_CALL( SCIPhashmapSetImage(expridx, (void *)nlhdlrexprdata->linexprs[idx], (void *)(size_t)idx) );
         }

         /* update index of expr */
         SCIP_CALL( SCIPhashmapSetImage(expridx, (void *)expr, (void *)(size_t)(-nlhdlrexprdata->nquadexprs - 1)) );

         /* update number of quadexprs */
         nlhdlrexprdata->nquadexprs++;
      }
      else /* expr has been seen before quadratically (expr^2 or expr * other_expr) */
      {
         quadexprterm = &nlhdlrexprdata->quadexprterms[-idx-1];

         if( bilinexprtermidx >= 0 ) /* quadratic expression seen again in a bilinear term; store which */
         {
            SCIP_CALL( nlhdlrexprdataEnsureAdjBilinSize(scip, quadexprterm, quadexprterm->nadjbilin + 1) );

            quadexprterm->adjbilin[quadexprterm->nadjbilin] = bilinexprtermidx;
            quadexprterm->nadjbilin++;
         }
         else /* first time seing quadratic expression in expr^2 from */
         {
            assert(bilinexprtermidx == -1);
            assert(quadexprterm->expr == expr);
            assert(quadexprterm->sqrcoef == 0.0);

            quadexprterm->sqrcoef = sqrcoef;
         }
      }
   }
   else /* first time seeing expression; it appears in a quadratically (expr^2 or expr * other_expr) -> new quadexprterm */
   {
      SCIP_CALL( nlhdlrexprdataEnsureQuadVarTermsSize(scip, nlhdlrexprdata, nlhdlrexprdata->nquadexprs + 1) );

      quadexprterm = &nlhdlrexprdata->quadexprterms[nlhdlrexprdata->nquadexprs];

      quadexprterm->expr = expr;
      quadexprterm->lincoef = 0.0;
      quadexprterm->nadjbilin = 0;
      quadexprterm->adjbilinsize = 0;
      quadexprterm->adjbilin = NULL;

      if( bilinexprtermidx >= 0 ) /* seen from a bilinear term; store which */
      {
         SCIP_CALL( nlhdlrexprdataEnsureAdjBilinSize(scip, quadexprterm, 1) );

         quadexprterm->adjbilin[quadexprterm->nadjbilin] = bilinexprtermidx;
         quadexprterm->nadjbilin++;
         quadexprterm->sqrcoef = 0.0;
      }
      else /* seen from a quadratic term */
      {
         quadexprterm->sqrcoef = sqrcoef;
      }

      /* add expr to expridx */
      SCIP_CALL( SCIPhashmapInsert(expridx, (void *)expr, (void *)(size_t)(-nlhdlrexprdata->nquadexprs - 1)) );

      /* update number of quadvars */
      nlhdlrexprdata->nquadexprs++;
   }
   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataQuadratic)
{  /*lint --e{715}*/
   freeNlhdlrExprData(scip, *nlhdlrexprdata);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree
 *
 * A term is quadratic if:
 * - It is a product expression of two expressions
 * - It is power expression of an expression with exponent 2.0
 *
 * A proper quadratic expression (i.e the only quadratic expressions that can be handled by this nlhdlr) is a sum
 * expression such that there is at least one (aux) variable that appears at least twice, in a quadratic terms and
 * somewhere else. In addition, it has to be convex or concave.
 * For example: x^2 + y^2 is not a proper quadratic expression; x^2 + x is proper quadratic expression;
 * x^2 + x * y is not a proper quadratic expression because it is not convex nor concave.
 *
 * @note:
 * - the expression needs to be simplified (in particular, it is assumed to be sorted)
 * - common subexpressions are also assumed to have been identified.
 *
 * Sorted implies that:
 *  - expr < expr^2: bases are the same, but exponent 1 < 2
 *  - expr < expr * other_expr: u*v < w holds if and only if v < w (OR8), but here w = u < v, since expr comes before
 *  other_expr in the product
 *  - expr < other_expr * expr: u*v < w holds if and only if v < w (OR8), but here v = w
 *
 * It also implies that
 *  - expr^2 < expr * other_expr
 *  - other_expr * expr < expr^2
 *
 * It also implies that x^-2 < x^-1, but since, so far, we do not interpret x^-2 as (x^-1)^2, it is not a problem.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(detectHdlrQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlexprdata;
   SCIP_HASHMAP*  expridx;
   SCIP_HASHMAP*  seenexpr;
   SCIP_Bool properquadratic;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* don't check if enforcement is already ensured */
   if( *enforcedbelow && *enforcedabove )
      return SCIP_OKAY;

   /* if it is not a sum of at least two terms, it cannot be a proper quadratic expressions */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(expr) < 2 )
      return SCIP_OKAY;

   /* check if expression is a proper quadratic expression */
   properquadratic = FALSE;
   SCIP_CALL( SCIPhashmapCreate(&seenexpr, SCIPblkmem(scip), 2*SCIPgetConsExprExprNChildren(expr)) );
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_CONSEXPR_EXPR* child;

      child = SCIPgetConsExprExprChildren(expr)[c];

      assert(child != NULL);

      if( strcmp("pow", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 &&
            SCIPgetConsExprExprPowExponent(child) == 2.0 ) /* quadratic term */
      {
         SCIP_CALL( checkProperQuadratic(SCIPgetConsExprExprChildren(child)[0], seenexpr, &properquadratic) );
      }
      else if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrProduct(conshdlr) &&
            SCIPgetConsExprExprNChildren(child) == 2 ) /* bilinear term */
      {
         SCIP_CALL( checkProperQuadratic(SCIPgetConsExprExprChildren(child)[0], seenexpr, &properquadratic) );
         SCIP_CALL( checkProperQuadratic(SCIPgetConsExprExprChildren(child)[1], seenexpr, &properquadratic) );
      }
      else
      {
         SCIP_CALL( checkProperQuadratic(child, seenexpr, &properquadratic) );
      }

      if( properquadratic )
         break;
   }
   SCIPhashmapFree(&seenexpr);

   if( ! properquadratic )
      return SCIP_OKAY;

   /* expridx maps expressions to indices; if index > 0, it is its index in the linexprs array, otherwise -index-1 is
    * its index in the quadexprterms array
    */
   SCIP_CALL( SCIPhashmapCreate(&expridx, SCIPblkmem(scip), SCIPgetConsExprExprNChildren(expr)) );

   /* sets everything to 0; nlexprdata->nquadexprs, etc */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   nlexprdata = *nlhdlrexprdata;

   /* for every term of the expr */
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_CONSEXPR_EXPR* child;
      SCIP_Real coef;

      child = SCIPgetConsExprExprChildren(expr)[c];
      coef = SCIPgetConsExprExprSumCoefs(expr)[c];

      assert(child != NULL);
      assert(! SCIPisZero(scip, coef)); /* TODO maybe this should be only coef != 0.0, since the original problem might have bad numerics */

      if( strcmp("pow", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 &&
            SCIPgetConsExprExprPowExponent(child) == 2.0 ) /* quadratic term */
      {
         assert(SCIPgetConsExprExprNChildren(child) == 1);

         SCIP_CALL( addExprToQuadexprterms(scip, SCIPgetConsExprExprChildren(child)[0], coef, -1, nlexprdata, expridx) );
      }
      else if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrProduct(conshdlr) &&
            SCIPgetConsExprExprNChildren(child) == 2 ) /* bilinear term */
      {
         SCIP_BILINEXPRTERM* bilinexprterm;
         SCIP_CONSEXPR_EXPR* expr1;
         SCIP_CONSEXPR_EXPR* expr2;

         assert(SCIPgetConsExprExprProductCoef(child) == 1.0);

         expr1 = SCIPgetConsExprExprChildren(child)[0];
         expr2 = SCIPgetConsExprExprChildren(child)[1];
         assert(expr1 != NULL && expr2 != NULL);

         SCIP_CALL( nlhdlrexprdataEnsureBilinSize(scip, nlexprdata, nlexprdata->nbilinexprterms + 1) );

         bilinexprterm = &nlexprdata->bilinexprterms[nlexprdata->nbilinexprterms];

         bilinexprterm->coef = coef;
         bilinexprterm->expr1 = expr1;
         bilinexprterm->expr2 = expr2;

         /* expression involved in a bilinear term that are not in a quadexprterm -> needs to be added to a quadexprterm
          * and removed from nlexprdata->linexprs
          */
         SCIP_CALL( addExprToQuadexprterms(scip, expr1, 0.0, nlexprdata->nbilinexprterms, nlexprdata, expridx) );
         SCIP_CALL( addExprToQuadexprterms(scip, expr2, 0.0, nlexprdata->nbilinexprterms, nlexprdata, expridx) );

         nlexprdata->nbilinexprterms++;
      }
      else /* linear term */
      {
         SCIP_CALL( nlhdlrexprdataEnsureLinearVarsSize(scip, nlexprdata, nlexprdata->nlinexprs + 1) );

         /* store its index in linexprs in case we see this expr again later in a bilinear or quadratic term */
         SCIP_CALL( SCIPhashmapInsert(expridx, child, (void *)(size_t)nlexprdata->nlinexprs) );
         nlexprdata->linexprs[nlexprdata->nlinexprs] = child;
         nlexprdata->lincoefs[nlexprdata->nlinexprs] = coef;
         nlexprdata->nlinexprs++;
      }
   }
   SCIPhashmapFree(&expridx);

   /* check curvature of quadratic function stored in nlexprdata */
   SCIP_CALL( checkCurvature(scip, nlexprdata) );

   if( nlexprdata->curvature == SCIP_EXPRCURV_CONVEX )
   {
      /* we will estimate the expression from below, that is handle expr <= auxvar */
      *enforcedbelow = TRUE;
      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   }
   else if( nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE )
   {
      /* we will estimate the expression from above, that is handle expr >= auxvar */
      *enforcedabove = TRUE;
      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   }
   else
   {
      /* we can't handle this expression, free data */
      SCIP_CALL( nlhdlrfreeExprDataQuadratic(scip, nlhdlr, nlhdlrexprdata) );
      return SCIP_OKAY;
   }

   {
      int i;

      /* create auxiliary variables for all expressions stored in nlhdlrexprdata */
      for( i = 0; i < nlexprdata->nlinexprs; ++i ) /* expressions appearing linearly */
      {
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexprdata->linexprs[i], NULL) );
      }
      for( i = 0; i < nlexprdata->nquadexprs; ++i ) /* quadratic terms */
      {
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexprdata->quadexprterms[i].expr, NULL) );
      }
      for( i = 0; i < nlexprdata->nbilinexprterms; ++i ) /* bilinear terms */
      {
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexprdata->bilinexprterms[i].expr1, NULL) );
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexprdata->bilinexprterms[i].expr2, NULL) );
      }
   }

   return SCIP_OKAY;
}

/** nonlinear handler separation callback */
static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrsepaHdlrQuadratic)
{  /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Real coef2;
   SCIP_Bool success;
   int j;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr));
   assert(result != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX || nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE);

   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );

   /* check violation side; one must use expression stored in nlhdlrexprdata!
    * TODO: maybe have a flag to know whether the expression is quadratic in the original variables
    */
   {
      SCIP_Real activity;
      SCIP_Real side;
      int i;

      activity = SCIPgetConsExprExprSumConstant(expr); /* TODO: is this okay or should the constant be stored at the moment of creation? */
      for( i = 0; i < nlhdlrexprdata->nlinexprs; ++i ) /* linear exprs */
         activity += nlhdlrexprdata->lincoefs[i] *
            SCIPgetSolVal(scip, sol, SCIPgetConsExprExprLinearizationVar(nlhdlrexprdata->linexprs[i]));

      for( i = 0; i < nlhdlrexprdata->nquadexprs; ++i ) /* quadratic terms */
      {
         SCIP_QUADEXPRTERM quadexprterm;
         SCIP_Real solval;

         quadexprterm = nlhdlrexprdata->quadexprterms[i];
         solval = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprLinearizationVar(quadexprterm.expr));
         activity += (quadexprterm.lincoef + quadexprterm.sqrcoef * solval) * solval;
      }
      for( i = 0; i < nlhdlrexprdata->nbilinexprterms; ++i ) /* bilinear terms */
      {
         SCIP_BILINEXPRTERM bilinexprterm;

         bilinexprterm = nlhdlrexprdata->bilinexprterms[i];
         activity += bilinexprterm.coef *
            SCIPgetSolVal(scip, sol, SCIPgetConsExprExprLinearizationVar(bilinexprterm.expr1)) *
            SCIPgetSolVal(scip, sol, SCIPgetConsExprExprLinearizationVar(bilinexprterm.expr2));
      }

      side = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprLinearizationVar(expr));

      SCIPdebugMsg(scip, "Activity = %g (act of expr is %g), side = %g, curvature %s\n", activity,
            SCIPgetConsExprExprValue(expr), side, nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX ? "convex" :
            "concave");

      if( activity > side && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX )
         rowprep->sidetype = SCIP_SIDETYPE_RIGHT;
      else if( activity < side && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE )
         rowprep->sidetype = SCIP_SIDETYPE_LEFT;
      else
         goto CLEANUP;
   }

   /*
    * compute cut: quadfun(sol) + \nabla quadfun(sol) (x - sol) - auxvar(expr)
    */

   /* constant */
   SCIPaddRowprepConstant(rowprep, SCIPgetConsExprExprSumConstant(expr));

   /* handle purely linear variables */
   for( j = 0; j < nlhdlrexprdata->nlinexprs; ++j )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprLinearizationVar(nlhdlrexprdata->linexprs[j]),
               nlhdlrexprdata->lincoefs[j]) );
   }

   /* quadratic variables */
   success = TRUE;
   for( j = 0; j < nlhdlrexprdata->nquadexprs && success; ++j )
   {
      int k;
      SCIP_VAR* var;
      var = SCIPgetConsExprExprLinearizationVar(nlhdlrexprdata->quadexprterms[j].expr);

      /* initialize coefficients to linear coefficients of quadratic variables */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, nlhdlrexprdata->quadexprterms[j].lincoef) );

      /* add linearization of square term */
      coef = 0.0;
      constant = 0.0;
      SCIPaddSquareLinearization(scip, nlhdlrexprdata->quadexprterms[j].sqrcoef, SCIPgetSolVal(scip, sol, var),
         nlhdlrexprdata->quadexprterms[j].nadjbilin == 0 && SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS, &coef, &constant, &success);

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
      SCIPaddRowprepConstant(rowprep, constant);

      /* add linearization of bilinear terms that have var as first variable */
      for( k = 0; k < nlhdlrexprdata->quadexprterms[j].nadjbilin; ++k )
      {
         SCIP_BILINEXPRTERM* bilinexprterm;
         SCIP_VAR* var2;

         bilinexprterm = &nlhdlrexprdata->bilinexprterms[nlhdlrexprdata->quadexprterms[j].adjbilin[k]];
         if( SCIPgetConsExprExprLinearizationVar(bilinexprterm->expr1) != var )
            continue;

         var2 = SCIPgetConsExprExprLinearizationVar(bilinexprterm->expr2);
         assert(var2 != NULL);
         assert(var2 != var);

         coef = 0.0;
         coef2 = 0.0;
         constant = 0.0;
         SCIPaddBilinLinearization(scip, bilinexprterm->coef, SCIPgetSolVal(scip, sol, var), SCIPgetSolVal(scip, sol,
                  var2), &coef, &coef2, &constant, &success);
         if( success )
         {
            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var2, coef2) );
            SCIPaddRowprepConstant(rowprep, constant);
         }
         else
            break;
      }
   }
   if( !success )
      goto CLEANUP;

   /* add auxiliary variable */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprLinearizationVar(expr), -1.0) );

   /* check build cut and check violation */
   {
      SCIP_Real coefrange;
      SCIP_Real viol;
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      viol = SCIPgetRowprepViolation(scip, rowprep, sol);

      if( viol <= 0.0 )
         goto CLEANUP;

      SCIPmergeRowprepTerms(scip, rowprep);

      /* improve coefficients */
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, 1e7, minviolation, &coefrange, &viol) );
      success = coefrange <= 1e7; /* magic number should maybe be given in as argument? */

      if( !success || viol < minviolation )
         goto CLEANUP;

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, conshdlr) );

      SCIP_CALL( SCIPaddRow(scip, row, TRUE, &infeasible) );
      (*ncuts)++;

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;
   }

CLEANUP:
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** nonlinear handler copy callback
 *
 * the method includes the nonlinear handler into a expression constraint handler
 *
 * This method is usually called when doing a copy of an expression constraint handler.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrcopyHdlrQuadratic)
{  /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrQuadratic(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes quadratic nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY,
            detectHdlrQuadratic, NULL) );

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrcopyHdlrQuadratic);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataQuadratic);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, nlhdlrsepaHdlrQuadratic, NULL);

   return SCIP_OKAY;
}
