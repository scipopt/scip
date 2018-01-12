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

   SCIP_INTERVAL         linactivity;        /**< activity of linear part */

   /* activities of quadratic parts as defined in nlhdlrIntervarQuadratic */
   SCIP_Real             minquadfiniteact;   /**< minimum activity of quadratic part where only terms with finite min
                                               activity contribute */
   SCIP_Real             maxquadfiniteact;   /**< maximum activity of quadratic part, where only terms with finite max
                                               activity contribute */
   int                   nneginfinityquadact;/**< number of quadratic terms contributing -infinity to activity */
   int                   nposinfinityquadact;/**< number of quadratic terms contributing +infinity to activity */
   SCIP_INTERVAL*        quadactivities;     /**< activity of each quadratic function as defined in nlhdlrIntevalQuadratic */
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

   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->quadactivities), nlhdlrexprdata->nquadexprs);

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

         /* when seen in a bilinear term -> store it */
         if( bilinexprtermidx >= 0 )
         {
            SCIP_CALL( nlhdlrexprdataEnsureAdjBilinSize(scip, quadexprterm, quadexprterm->nadjbilin + 1) );

            quadexprterm->adjbilin[quadexprterm->nadjbilin] = bilinexprtermidx;
            quadexprterm->nadjbilin++;
         }

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

/** creates auxiliary variable when necessary */
static
SCIP_RETCODE createAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expr conshdlr */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to add to quadratic terms */
   SCIP_Bool*            originalvar         /**< set it to false when expression is not var */
   )
{
   if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrVar(conshdlr) )
      return SCIP_OKAY;

   *originalvar = FALSE;
   SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, NULL) );

   return SCIP_OKAY;
}

/** solves a quadratic equation \f$ a expr^2 + b expr \in rhs \f$ (with b an interval) and reduces bounds on expr or
 * deduces infeasibility if possible; expr is quadexpr.expr
 */
static
SCIP_RETCODE propagateBoundsQuadExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_QUADEXPRTERM     quadexpr,           /**< quadratic expression to propagate */
   SCIP_INTERVAL         b,                  /**< interval acting as linear coefficient */
   SCIP_INTERVAL         rhs,                /**< interval acting as rhs */
   SCIP_QUEUE*           reversepropqueue,   /**< queue used in reverse prop, pass to SCIPtightenConsExprExprInterval */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions,        /**< buffer to store the number of interval reductions */
   SCIP_Bool             force               /**< to force tightening */
   )
{
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   /* compute solution of a*x^2 + b*x \in rhs */
   if( quadexpr.sqrcoef == 0.0 && SCIPintervalGetInf(b) == 0.0 && SCIPintervalGetSup(b) == 0.0 )
   {
      /* relatively easy case: 0.0 \in rhs, thus check if infeasible or just redundant */
      if( SCIPintervalGetInf(rhs) > 0.0 || SCIPintervalGetSup(rhs) < 0.0 )
         SCIPintervalSetEmpty(&newrange);
      else
         return SCIP_OKAY;
   }
   else if( SCIPintervalGetInf(SCIPgetConsExprExprInterval(quadexpr.expr)) >= 0.0 )
   {
      SCIP_INTERVAL a;

      /* need only positive solutions */
      SCIPintervalSet(&a, quadexpr.sqrcoef);
      SCIPintervalSolveUnivariateQuadExpressionPositive(SCIP_INTERVAL_INFINITY, &newrange, a, b, rhs);
   }
   else if( SCIPintervalGetSup(SCIPgetConsExprExprInterval(quadexpr.expr)) <= 0.0 )
   {
      /* need only negative solutions */
      SCIP_INTERVAL a;
      SCIP_INTERVAL tmp;
      SCIPintervalSet(&a, quadexpr.sqrcoef);
      SCIPintervalSetBounds(&tmp, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));
      SCIPintervalSolveUnivariateQuadExpressionPositive(SCIP_INTERVAL_INFINITY, &tmp, a, tmp, rhs);

      SCIPintervalSetBounds(&newrange, -SCIPintervalGetSup(tmp), -SCIPintervalGetInf(tmp));
   }
   else
   {
      /* need both positive and negative solution */
      SCIP_INTERVAL a;
      SCIPintervalSet(&a, quadexpr.sqrcoef);
      SCIPintervalSolveUnivariateQuadExpression(SCIP_INTERVAL_INFINITY, &newrange, a, b, rhs);
   }

   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, quadexpr.expr, newrange, force, reversepropqueue, infeasible,
            nreductions) );

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

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "Nlhdlr quadratic detecting expr %p aka", (void*)expr);
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif
   SCIPdebugMsg(scip, "checking if expr %p is a proper quadratic\n", (void*)expr);
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

   SCIPdebugMsg(scip, "expr %p is proper quadratic: checking convexity\n", (void*)expr);

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
      assert(coef != 0.0);

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

   /* every detected quadratic expression will be handled since we can propagate */
   *enforcedbelow = FALSE;
   *enforcedabove = FALSE;
   *success = TRUE;
   *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_INTEVAL | SCIP_CONSEXPR_EXPRENFO_REVERSEPROP;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlexprdata->quadactivities, nlexprdata->nquadexprs) );

   /* check if we can do something more: check curvature of quadratic function stored in nlexprdata */
   SCIP_CALL( checkCurvature(scip, nlexprdata) );

   if( nlexprdata->curvature == SCIP_EXPRCURV_CONVEX )
   {
      SCIPdebugMsg(scip, "expr %p is convex when replacing factors of bilinear terms, bases of squares and every other term by their aux vars\n",
            (void*)expr);

      /* we will estimate the expression from below, that is handle expr <= auxvar */
      *enforcedbelow = TRUE;
      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   }
   else if( nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE )
   {
      SCIPdebugMsg(scip, "expr %p is concave when replacing factors of bilinear terms, bases of squares and every other term by their aux vars\n",
            (void*)expr);

      /* we will estimate the expression from above, that is handle expr >= auxvar */
      *enforcedabove = TRUE;
      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   }
   else
   {
      /* we cannot do more with this quadratic function */
      return SCIP_OKAY;
   }

   /* quadratic expression is concave/convex -> create aux vars for all expressions stored in nlhdlrexprdata */
   {
      int i;
      SCIP_Bool originalvars = TRUE;

      for( i = 0; i < nlexprdata->nlinexprs; ++i ) /* expressions appearing linearly */
      {
         SCIP_CALL( createAuxVar(scip, conshdlr, nlexprdata->linexprs[i], &originalvars) );
      }
      for( i = 0; i < nlexprdata->nquadexprs; ++i ) /* quadratic terms */
      {
         SCIP_CALL( createAuxVar(scip, conshdlr, nlexprdata->quadexprterms[i].expr, &originalvars) );
      }
      for( i = 0; i < nlexprdata->nbilinexprterms; ++i ) /* bilinear terms */
      {
         SCIP_CALL( createAuxVar(scip, conshdlr, nlexprdata->bilinexprterms[i].expr1, &originalvars) );
         SCIP_CALL( createAuxVar(scip, conshdlr, nlexprdata->bilinexprterms[i].expr2, &originalvars) );
      }

      if( originalvars )
      {
         SCIPsetConsExprExprCurvature(expr, nlexprdata->curvature);
         SCIPdebugMsg(scip, "expr is %s in the original variables\n", nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE ? "concave" : "convex");
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
            SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(nlhdlrexprdata->linexprs[i]));

      for( i = 0; i < nlhdlrexprdata->nquadexprs; ++i ) /* quadratic terms */
      {
         SCIP_QUADEXPRTERM quadexprterm;
         SCIP_Real solval;

         quadexprterm = nlhdlrexprdata->quadexprterms[i];
         solval = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(quadexprterm.expr));
         activity += (quadexprterm.lincoef + quadexprterm.sqrcoef * solval) * solval;
      }
      for( i = 0; i < nlhdlrexprdata->nbilinexprterms; ++i ) /* bilinear terms */
      {
         SCIP_BILINEXPRTERM bilinexprterm;

         bilinexprterm = nlhdlrexprdata->bilinexprterms[i];
         activity += bilinexprterm.coef *
            SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(bilinexprterm.expr1)) *
            SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(bilinexprterm.expr2));
      }

      side = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));

      SCIPdebugMsg(scip, "Activity = %g (act of expr is %g), side = %g, curvature %s\n", activity,
            SCIPgetConsExprExprValue(expr), side, nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX ? "convex" :
            "concave");

      if( activity - side > minviolation && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX )
         rowprep->sidetype = SCIP_SIDETYPE_RIGHT;
      else if( minviolation < side - activity && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE )
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
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprAuxVar(nlhdlrexprdata->linexprs[j]),
               nlhdlrexprdata->lincoefs[j]) );
   }

   /* quadratic variables */
   success = TRUE;
   for( j = 0; j < nlhdlrexprdata->nquadexprs && success; ++j )
   {
      int k;
      SCIP_VAR* var;
      var = SCIPgetConsExprExprAuxVar(nlhdlrexprdata->quadexprterms[j].expr);

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
         if( SCIPgetConsExprExprAuxVar(bilinexprterm->expr1) != var )
            continue;

         var2 = SCIPgetConsExprExprAuxVar(bilinexprterm->expr2);
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
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprAuxVar(expr), -1.0) );

   /* check build cut and check violation */
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      SCIPmergeRowprepTerms(scip, rowprep);

      /* improve coefficients */
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, minviolation, NULL, &success) );

      if( !success )
         goto CLEANUP;

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, conshdlr) );

      SCIP_CALL( SCIPaddRow(scip, row, TRUE, &infeasible) );
      (*ncuts)++;

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

CLEANUP:
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** nonlinear handler forward propagation callback
 * This method should solve the problem
 * max/min quad expression over box constraints
 * However, this problem is difficult so we are satisfied with a proxy.
 * Interval arithmetic suffices when no variable appears twice, however this is seldom the case, so we try
 * to take care of the dependency problem to some extent:
 * 1. partition the quadratic expression as sum of quadratic functions
 * \sum_l q_l
 * where q_l = a_l expr_l^2 + \sum_{i \in P_l} b_il expr_i expr_l + c_l expr_l
 * 2. build interval quadratic functions, i.e, a x^2 + b x where b is an interval as
 * a_l expr_l^2 + [\sum_{i \in P_l} b_il expr_i + c_l] expr_l
 * 3. compute \min and \max { a x^2 + b x : x \in [x] } for each interval quadratic, i.e.
 * \min and \max a_l expr_l^2 + [\sum_{i \in P_l} b_il expr_i + c_l] expr_l : expr_l \in [expr_l]
 *
 * In particular, P_l = \{i : expr_l expr_i is a bilinear expr\}. Note that the
 * order matters, that is in P_l, expr_l is the the first expression.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalQuadratic)
{ /*lint --e{715}*/

   SCIP_INTERVAL quadactivity;

   assert(scip != NULL);
   assert(expr != NULL);

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nquadexprs != 0);

   SCIPdebugMsg(scip, "Interval evaluation of quadratic expr\n");
#ifdef DEBUG_PROP
   SCIP_CALL( SCIPdismantleConsExprExpr(scip, expr) );
#endif

   /*
    * compute activity of linear part
    */
   {
      int i;

      SCIPdebugMsg(scip, "Computing activity of linear part\n");

      SCIPintervalSet(&nlhdlrexprdata->linactivity, SCIPgetConsExprExprSumConstant(expr));
      for( i = 0; i < nlhdlrexprdata->nlinexprs; ++i )
      {
         SCIP_INTERVAL linterminterval;

         SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &linterminterval,
               SCIPgetConsExprExprInterval(nlhdlrexprdata->linexprs[i]), nlhdlrexprdata->lincoefs[i]);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &nlhdlrexprdata->linactivity, nlhdlrexprdata->linactivity, linterminterval);
      }

      SCIPdebugMsg(scip, "Activity of linear part is [%g, %g]\n", nlhdlrexprdata->linactivity.inf,
            nlhdlrexprdata->linactivity.sup);
   }

   /*
    * compute activity of quadratic part
    */
   nlhdlrexprdata->nneginfinityquadact = 0;
   nlhdlrexprdata->nposinfinityquadact = 0;
   nlhdlrexprdata->minquadfiniteact = 0.0;
   nlhdlrexprdata->maxquadfiniteact = 0.0;
   SCIPintervalSet(&quadactivity, 0.0);
   {
      SCIP_BILINEXPRTERM* bilinterms;
      int i;

      bilinterms = nlhdlrexprdata->bilinexprterms;
      for( i = 0; i < nlhdlrexprdata->nquadexprs; ++i )
      {
         int j;
         SCIP_INTERVAL b;
         SCIP_QUADEXPRTERM quadexpr;
         SCIP_Real quadub;
         SCIP_Real quadlb;

         /* b = [c_l] */
         quadexpr = nlhdlrexprdata->quadexprterms[i];
         SCIPintervalSet(&b, quadexpr.lincoef);
         for( j = 0; j < quadexpr.nadjbilin; ++j )
         {
            SCIP_BILINEXPRTERM bilinterm;
            SCIP_INTERVAL bterm;

            bilinterm = bilinterms[quadexpr.adjbilin[j]];
            if( bilinterm.expr1 != quadexpr.expr )
               continue;

            /* b += [b_jl * expr_j] for j \in P_l */
            SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, SCIPgetConsExprExprInterval(bilinterm.expr2),
                  bilinterm.coef);
            SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);
         }

         /* TODO: under which assumptions do we know that we just need to compute min or max? */
         quadub = SCIPintervalQuadUpperBound(SCIP_INTERVAL_INFINITY, quadexpr.sqrcoef, b,
               SCIPgetConsExprExprInterval(quadexpr.expr));

         /* TODO: implement SCIPintervalQuadLowerBound */
         {
            SCIP_INTERVAL minusb;
            SCIPintervalSetBounds(&minusb, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));

            quadlb = -SCIPintervalQuadUpperBound(SCIP_INTERVAL_INFINITY, -quadexpr.sqrcoef, minusb,
                  SCIPgetConsExprExprInterval(quadexpr.expr));
         }

         SCIPintervalSetBounds(&nlhdlrexprdata->quadactivities[i], quadlb, quadub);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &quadactivity, quadactivity, nlhdlrexprdata->quadactivities[i]);

         /* get number of +/-infinity contributions and compute finite activity */
         if( quadlb <= -SCIP_INTERVAL_INFINITY )
            nlhdlrexprdata->nneginfinityquadact++;
         else
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();

            nlhdlrexprdata->minquadfiniteact += quadlb;

            SCIPintervalSetRoundingMode(roundmode);
         }
         if( quadub >= SCIP_INTERVAL_INFINITY )
            nlhdlrexprdata->nneginfinityquadact++;
         else
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();

            nlhdlrexprdata->maxquadfiniteact += quadub;

            SCIPintervalSetRoundingMode(roundmode);
         }
      }

      SCIPdebugMsg(scip, "Activity of quadratic part is [%g, %g]\n", quadactivity.inf, quadactivity.sup);
   }

   /* interval evaluation is linear activity + quadactivity */
   SCIPintervalAdd(SCIP_INTERVAL_INFINITY, interval, nlhdlrexprdata->linactivity,  quadactivity);

   return SCIP_OKAY;
}


/** nonlinear handler reverse propagation callback
 * @note: the implemented technique is a proxy for solving the OBBT problem min/max{ x_i : quad expr \in [quad expr] }
 * and as such can be improved.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - reversepropqueue : expression queue in reverse propagation, to be passed on to SCIPtightenConsExprExprInterval
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 *  - force : force tightening even if it is below the bound strengthening tolerance
 */
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropQuadratic)
{ /*lint --e{715}*/
   SCIP_INTERVAL rhs;
   SCIP_INTERVAL quadactivity;

   SCIPdebugMsg(scip, "Reverse propagation of quadratic expr\n");

   assert(scip != NULL);
   assert(expr != NULL);
   assert(reversepropqueue != NULL);
   assert(infeasible != NULL);

   /* not possible to conclude finite bounds if the interval of the expression is [-inf,inf] */
   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, SCIPgetConsExprExprInterval(expr)) )
      return SCIP_OKAY;

   /* propagate linear part in rhs = expr's interval - quadratic activity; first, reconstruct the quadratic activity */
   SCIPintervalSetBounds(&quadactivity,
         nlhdlrexprdata->nneginfinityquadact > 0 ? -SCIP_INTERVAL_INFINITY : nlhdlrexprdata->minquadfiniteact,
         nlhdlrexprdata->nposinfinityquadact > 0 ?  SCIP_INTERVAL_INFINITY : nlhdlrexprdata->maxquadfiniteact);

   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, SCIPgetConsExprExprInterval(expr), quadactivity);
   SCIP_CALL( SCIPreverseConsExprExprPropagateWeightedSum(scip, nlhdlrexprdata->nlinexprs,
            nlhdlrexprdata->linexprs, nlhdlrexprdata->lincoefs, SCIPgetConsExprExprSumConstant(expr),
            rhs, reversepropqueue, infeasible, nreductions, force) );

   /* stop if we find infeasibility */
   if( *infeasible )
      return SCIP_OKAY;

   /* propagate quadratic part in expr's interval - linear activity:
    * linear activity was computed in INTEVAL
    * One way of achieving this is by, for each expression expr_i, write the quadratic expression as
    * a_i expr^2_i + expr_i ( \sum_{j \in J_i} b_ij expr_j + c_i ) + quadratic expression in expr_k for k \neq i
    * then compute the interval b = [\sum_{j \in J_i} b_ij expr_j + c_i], where J_i are all the indices j such that the
    * bilinear expression expr_i expr_j appears, and use some technique (like the one in nlhdlrIntevalQuadratic), to
    * evaluate the activity rest_i = [quadratic expression in expr_k for k \neq i].
    * Then, solve a_i expr_i^2 + b expr_i = [expr] - rest_i =: rhs_i.
    * However, this might be expensive, specially computing rest_i. Hence, we implement a simpler version, namely,
    * we use the same partition as in nlhdlrIntevalQuadratic for the bilinear terms. This way,
    * b = [\sum_{j \in P_i} b_ij expr_j + c_i], where P_i is defined as in nlhdlrIntevalQuadratic, all the indices j
    * such that expr_i expr_j appears in that order, and rest_i = sum_{k \neq i} [\min q_k, \max q_k] where
    * q_k = a_k expr_k^2 + [\sum_{j \in P_k} b_jk expr_j + c_k] expr_k. The intervals [\min q_k, \max q_k] were
    * already computed in nlhdlrIntevalQuadratic, so we just reuse them.
    *
    * TODO: in cons_quadratic there seems to be a further technique that tries, when propagating expr_i, to borrow a
    * bilinear term expr_k expr_i when the quadratic function for expr_k is simple enough.
    *
    * TODO: handle simple cases
    * TODO: identify early when there is nothing to be gain
    */
   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, SCIPgetConsExprExprInterval(expr), nlhdlrexprdata->linactivity);
   {
      SCIP_BILINEXPRTERM* bilinterms;
      int i;

      bilinterms = nlhdlrexprdata->bilinexprterms;
      for( i = 0; i < nlhdlrexprdata->nquadexprs; ++i )
      {
         int j;
         SCIP_INTERVAL b;
         SCIP_INTERVAL rhs_i;
         SCIP_INTERVAL rest_i;
         SCIP_QUADEXPRTERM quadexpr;

         /* b = [c_l] */
         quadexpr = nlhdlrexprdata->quadexprterms[i];
         SCIPintervalSet(&b, quadexpr.lincoef);
         for( j = 0; j < quadexpr.nadjbilin; ++j )
         {
            SCIP_BILINEXPRTERM bilinterm;
            SCIP_INTERVAL bterm;

            bilinterm = bilinterms[quadexpr.adjbilin[j]];
            if( bilinterm.expr1 != quadexpr.expr )
               continue;

            /* b += [b_jl * expr_j] for j \in P_l */
            SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, SCIPgetConsExprExprInterval(bilinterm.expr2),
                  bilinterm.coef);
            SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);
         }

         /* rhs_i = rhs - rest_i.
          * to compute rest_i = [\sum_{k \neq i} q_k] we just have to substract
          * the activity of q_i from quadactivity; however, care must be taken about infinities;
          * if [q_i].sup = +infinity and there is = 1 contributing +infinity -> rest_i.sup = maxquadfiniteact
          * if [q_i].sup = +infinity and there is > 1 contributing +infinity -> rest_i.sup = +infinity
          * if [q_i].sup = finite and there is > 0 contributing +infinity -> rest_i.sup = +infinity
          * if [q_i].sup = finite and there is = 0 contributing +infinity -> rest_i.sup = maxquadfiniteact - [q_i].sup
          *
          * the same holds when replacing sup with inf, + with - and max(quadfiniteact) with min(...)
          */
         /* compute rest_i.sup */
         if( SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]) < SCIP_INTERVAL_INFINITY &&
               nlhdlrexprdata->nposinfinityquadact == 0 )
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeUpwards();
            rest_i.sup = nlhdlrexprdata->maxquadfiniteact - SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]);

            SCIPintervalSetRoundingMode(roundmode);
         }
         else if( SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]) >= SCIP_INTERVAL_INFINITY &&
               nlhdlrexprdata->nposinfinityquadact == 1 )
            rest_i.sup = nlhdlrexprdata->maxquadfiniteact;
         else
            rest_i.sup = SCIP_INTERVAL_INFINITY;

         /* compute rest_i.inf */
         if( SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]) > -SCIP_INTERVAL_INFINITY &&
               nlhdlrexprdata->nneginfinityquadact == 0 )
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();
            rest_i.inf = nlhdlrexprdata->minquadfiniteact - SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]);

            SCIPintervalSetRoundingMode(roundmode);
         }
         else if( SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]) <= -SCIP_INTERVAL_INFINITY &&
               nlhdlrexprdata->nneginfinityquadact == 1 )
            rest_i.inf = nlhdlrexprdata->minquadfiniteact;
         else
            rest_i.inf = -SCIP_INTERVAL_INFINITY;

         /* compute rhs_i */
         SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs_i, rhs, rest_i);

         /* solve a_i expr_i^2 + b expr_i = rhs_i */
         if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, rhs_i) )
            continue;

         SCIP_CALL( propagateBoundsQuadExpr(scip, quadexpr, b, rhs_i, reversepropqueue, infeasible, nreductions, force) );

         /* stop if we find infeasibility */
         if( *infeasible )
            return SCIP_OKAY;
      }
   }

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
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalQuadratic, nlhdlrReversepropQuadratic);

   return SCIP_OKAY;
}
