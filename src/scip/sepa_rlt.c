/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_rlt.c
 * @brief  RLT separator
 * @author Fabian Wegscheider
 *
 * @todo implement the possibility to add extra auxiliary variables for RLT (like in DOI 10.1080/10556788.2014.916287)
 * @todo add RLT cuts for the product of equality constraints
 * @todo implement dynamic addition of RLT cuts during branching (see DOI 10.1007/s10898-012-9874-7)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "scip/sepa_rlt.h"
#include "scip/cons_expr.h"
#include "cons_linear.h"
#include "cons_knapsack.h"
#include "cons_varbound.h"
#include "cons_setppc.h"
#include "cons_expr_iterator.h"
#include "cons_expr_pow.h"
#include "cons_expr_product.h"
#include "cons_expr_sum.h"
#include "cons_expr_var.h"
#include "struct_cons_expr.h"
#include "scip/struct_scip.h"
#include "set.h"


#define SEPA_NAME              "rlt"
#define SEPA_DESC              "rlt separator"
#define SEPA_PRIORITY                10 /**< priority for separation */
#define SEPA_FREQ                     0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define SEPA_MAXBOUNDDIST           1.0 /**< maximal relative distance from the current node's dual bound to primal bound
+                                        *   compared to best node's dual bound for applying separation.*/
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXLINEXPRS          10 /**< default value for parameter maxlinexprs */
#define DEFAULT_MAXUNKNOWNTERMS      -1 /**< default value for parameter maxunknownterms */
#define DEFAULT_MAXUSEDVARS          -1 /**< default value for parameter maxusedvars */
#define DEFAULT_MAXNONZEROPROP      0.0 /**< default value for parameter maxnonzeroprop */
#define DEFAULT_MAXNCUTS             -1 /**< default value for parameter maxncuts */
#define DEFAULT_MAXROUNDS             1 /**< default value for parameter maxrounds */
#define DEFAULT_MAXROUNDSROOT        10 /**< default value for parameter maxroundsroot */
#define DEFAULT_ONLYEQROWS        FALSE /**< default value for parameter eqrowsfirst */
#define DEFAULT_ONLYCONTROWS      FALSE /**< default value for parameter eqrowsfirst */
#define DEFAULT_ONLYINITIAL        TRUE /**< default value for parameter onlyinitial */
#define DEFAULT_USEINSUBSCIP      FALSE /**< default value for parameter useinsubscip */

#define MAXVARBOUND                1e+5 /**< maximum allowed variable bound for computing an RLT-cut */

/*
 * Data structures
 */

/** data object to compare constraint easier */
struct HashData
{
   SCIP_ROW**            rows;               /**< pointer to the rows defined by the same three variables */
   SCIP_VAR**            vars;               /**< constraint variables used for hash comparison */
   int                   nvars;              /**< number of variables */
   int                   nrows;              /**< number of rows */
};
typedef struct HashData HASHDATA;


/** separator data */
struct SCIP_SepaData
{
   SCIP_CONSHDLR*        conshdlr;           /**< expression constraint handler */
   SCIP_VAR**            varssorted;         /**< variables that occur in bilinear terms sorted by priority */

   SCIP_CONSEXPR_EXPR*** linexprs;           /**< arrays of linearisation expressions for each bilinear term */
   int*                  nlinexprs;          /**< number of linearisation expressions for each bilinear term */
   int*                  slinexprs;          /**< sizes of linexprs arrays for each bilinear term */

   SCIP_Bool*            isimplicit;         /**< whether bilinterms[i] is an implicit bilinear relation */
   SCIP_CONSEXPR_EXPR**  bilinterms;         /**< bilinear terms */
   int*                  bilinlockspos;      /**< positive locks of the bilinear terms */
   int*                  bilinlocksneg;      /**< negative locks of the bilinear terms */
   SCIP_HASHMAP*         bilinvarsmap;       /**< map for accessing the exprs and linearization variables of each bilinear term */
   SCIP_VAR***           varbilinvars;       /**< arrays of vars appearing in a bilinear term together with xj for each xj from varssorted */
   int*                  nvarbilinvars;      /**< number of vars for each element of varbilinvars */
   int*                  varpriorities;      /**< priorities of the variables in varssorted */
   int                   maxvarindex;        /**< maximum variable index when creating bilinvarsmap */
   int                   nbilinterms;        /**< total number of bilinear terms */
   int                   nbilinvars;         /**< total number of variables occurring in bilinear terms */
   int                   currentnunknown;    /**< number of unknown terms in current row (not printed) */
   SCIP_Bool             iscreated;          /**< indicates whether the sepadata has been initialized yet */
   SCIP_Bool             isinitialround;     /**< indicates that this is the first round and initial rows are used */

   /* parameters */
   int                   maxlinexprs;        /**< maximum number of linearisation expressions per bilinear term */
   SCIP_Real             maxnonzeroprop;     /**< maximum acceptable proportion of known bilinear terms to non-zeroes */
   int                   maxunknownterms;    /**< maximum number of unknown bilinear terms a row can have to be used */
   int                   maxusedvars;        /**< maximum number of variables that will be used to compute rlt cuts */
   int                   maxncuts;           /**< maximum number of cuts that will be added per round */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   SCIP_Bool             onlyeqrows;         /**< indicates wether only equality rows should be used for rlt cuts */
   SCIP_Bool             onlycontrows;       /**< indicates wether only continuous rows should be used for rlt cuts */
   SCIP_Bool             onlyinitial;        /**< indicates whether only initial rows should be uswed for rlt cuts */
   SCIP_Bool             useinsubscip;       /**< indicates whether the seperator should also be used in sub-scips */
};

/** projected LP data structure */
struct ProjLP
{
   SCIP_Real** coefs;       /* arrays of coefficients for each row */
   SCIP_VAR*** vars;        /* arrays of variables for each row */
   SCIP_Real*  lhss;        /* row left hand sides */
   SCIP_Real*  rhss;        /* row right hand sides */
   SCIP_Real*  consts;      /* row constants */
   int*        nNonz;       /* number of nonzeros in each row */
};
typedef struct ProjLP PROJLP;

/*
 * Local methods
 */

/** returns TRUE iff both keys are equal; two constraint arrays are equal if they have the same pointer */
static
SCIP_DECL_HASHKEYEQ(hashdataKeyEqConss)
{
   HASHDATA* hashdata1;
   HASHDATA* hashdata2;
   int v;

   hashdata1 = (HASHDATA*)key1;
   hashdata2 = (HASHDATA*)key2;

   /* check data structure */
   assert(hashdata1->nvars == 3);
   assert(hashdata2->nvars == 3);
   assert(hashdata1->rows != NULL || hashdata2->rows != NULL);

   for( v = 2; v >= 0; --v )
   {
      /* tests if variables are equal */
      if( hashdata1->vars[v] != hashdata2->vars[v] )
         return FALSE;

      assert(SCIPvarCompare(hashdata1->vars[v], hashdata2->vars[v]) == 0);
   }

   /* a hashdata object is only equal if it has the same constraint array pointer */
   if( hashdata1->rows == NULL || hashdata2->rows == NULL || hashdata1->rows == hashdata2->rows )
      return TRUE;
   else
      return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashdataKeyValConss)
{  /*lint --e{715}*/
   HASHDATA* hashdata;
   int minidx, mididx, maxidx;

   hashdata = (HASHDATA*)key;
   assert(hashdata != NULL);
   assert(hashdata->vars != NULL);
   assert(hashdata->nvars == 3);

   minidx = SCIPvarGetIndex(hashdata->vars[0]);
   mididx = SCIPvarGetIndex(hashdata->vars[1]);
   maxidx = SCIPvarGetIndex(hashdata->vars[2]);

   /* vars should already be sorted by index */
   assert(minidx < mididx && mididx < maxidx);

   return SCIPhashTwo(SCIPcombineTwoInt(minidx, mididx), maxidx);
}


/* helper method to free the separation data */
static
SCIP_RETCODE freeSepaData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separation data */
   )
{  /*lint --e{715}*/
   int i, j;

   assert(sepadata->iscreated);
   assert(sepadata->bilinvarsmap != NULL);

   /* release auxiliary variables that were captured for rlt */
   for( i = 0; i < sepadata->nbilinterms; ++i )
   {
      assert(sepadata->bilinterms[i] != NULL);
      assert(sepadata->linexprs[i] != NULL);

      for( j = 0; j < sepadata->nlinexprs[i]; ++j )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(sepadata->linexprs[i][j])) );
      }
      SCIPfreeBlockMemoryArray(scip, &(sepadata->linexprs[i]), sepadata->slinexprs[i]);

      if( sepadata->isimplicit[i] )
      {  /* the separator sets the locks for implicit products, so they have to be removed here */
         sepadata->bilinterms[i]->nlockspos = 0;
         sepadata->bilinterms[i]->nlocksneg = 0;
      }
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(sepadata->bilinterms[i])) );
   }

   /* release bilinvars that were captured for rlt */
   for( i = 0; i < sepadata->nbilinvars; ++i )
   {
      assert(sepadata->varssorted[i] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &(sepadata->varssorted[i])) );
   }

   /* free arrays */
   for( i = 0; i < sepadata->nbilinvars; ++i )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &sepadata->varbilinvars[i], sepadata->nvarbilinvars[i]);
   }

   SCIPfreeBlockMemoryArray(scip, &sepadata->isimplicit, sepadata->nbilinterms);
   SCIPfreeBlockMemoryArray(scip, &sepadata->bilinterms, sepadata->nbilinterms);
   SCIPfreeBlockMemoryArray(scip, &sepadata->bilinlockspos, sepadata->nbilinterms);
   SCIPfreeBlockMemoryArray(scip, &sepadata->bilinlocksneg, sepadata->nbilinterms);

   SCIPfreeBlockMemoryArray(scip, &sepadata->linexprs, sepadata->nbilinterms);
   SCIPfreeBlockMemoryArray(scip, &sepadata->nlinexprs, sepadata->nbilinterms);
   SCIPfreeBlockMemoryArray(scip, &sepadata->slinexprs, sepadata->nbilinterms);

   SCIPfreeBlockMemoryArray(scip, &sepadata->varpriorities, sepadata->nbilinvars);
   SCIPfreeBlockMemoryArray(scip, &sepadata->nvarbilinvars, sepadata->nbilinvars);
   SCIPfreeBlockMemoryArray(scip, &sepadata->varbilinvars, sepadata->nbilinvars);
   SCIPfreeBlockMemoryArray(scip, &sepadata->varssorted, sepadata->nbilinvars);
   /* TODO check order of freeing */

   /* free the hashmap */
   SCIPhashmapFree(&sepadata->bilinvarsmap);

   sepadata->iscreated = FALSE;

   return SCIP_OKAY;
}

/** creates and returns rows of initial linear constraints */
static
SCIP_RETCODE getInitialRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW***           rows,               /**< buffer to store the rows */
   int*                  nrows               /**< buffer to store the number of linear rows */
)
{
   SCIP_CONS** conss;
   SCIP_CONSHDLR* linhdlr;
   SCIP_CONSHDLR* knpsckhdlr;
   SCIP_CONSHDLR* varbndhdlr;
   SCIP_CONSHDLR* setppchdlr;
   int nconss;
   int i;

   assert(rows != NULL);
   assert(nrows != NULL);

   linhdlr = SCIPfindConshdlr(scip, "linear");
   knpsckhdlr = SCIPfindConshdlr(scip, "knapsack");
   varbndhdlr = SCIPfindConshdlr(scip, "varbound");
   setppchdlr = SCIPfindConshdlr(scip, "setppc");
   /* TODO any other handlers? */

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   *nrows = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, rows, nconss) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_ROW *row;

      if( SCIPconsGetHdlr(conss[i]) == linhdlr )
      {
         row = SCIPgetRowLinear(scip, conss[i]);
         SCIPdebugMsg(scip, "linear constraint found\n");
      }
      else if( SCIPconsGetHdlr(conss[i]) == knpsckhdlr )
      {
         row = SCIPgetRowKnapsack(scip, conss[i]);
         SCIPdebugMsg(scip, "knapsack constraint found\n");
      }
      else if( SCIPconsGetHdlr(conss[i]) == varbndhdlr )
      {
         row = SCIPgetRowVarbound(scip, conss[i]);
         SCIPdebugMsg(scip, "varbound constraint found\n");
      }
      else if( SCIPconsGetHdlr(conss[i]) == setppchdlr )
      {
         row = SCIPgetRowSetppc(scip, conss[i]);
         SCIPdebugMsg(scip, "setppc constraint found\n");
      }
      else
      {
         continue;
      }

      if( row != NULL)
      {
         (*rows)[*nrows] = row;
         ++*nrows;
      }
   }

   return SCIP_OKAY;
}


/** store the bilinear product and all related information in sepadata */
static
SCIP_RETCODE addBilinProduct(
   SCIP*               scip,       /**< SCIP data structure */
   SCIP_SEPADATA*      sepadata,   /**< separator data */
   SCIP_CONSEXPR_EXPR* expr,       /**< product expression */
   SCIP_VAR*           x,          /**< a product variable with smaller index */
   SCIP_VAR*           y,          /**< a product variable with larger index */
   SCIP_CONSEXPR_EXPR* linexpr,    /**< linearisation expression */
   SCIP_Bool           isimplicit, /**< whether the product is implicit */
   SCIP_Bool           lhs,        /**< is the product bounded from below? */
   SCIP_Bool           rhs,        /**< is the product bounded from above? */
   SCIP_HASHMAP*       varmap      /**< map containing vars that have already been encountered in products */
)
{
   int mapidx, xidx, yidx;
   int xpos, ypos;
   SCIP_CONSEXPR_EXPR* prodexprs[2];
   SCIP_Bool new = FALSE;
   int termpos, linpos;
   SCIP_Bool found;

   assert(lhs || rhs);

   /* the variables should be given in the correct order */
   if( SCIPvarComp(x, y) > 0 )
      SCIPswapPointers((void**)&x, (void**)&y);

   xidx = SCIPvarGetIndex(x);
   yidx = SCIPvarGetIndex(y);

   mapidx = xidx * sepadata->maxvarindex + yidx;

   termpos = SCIPhashmapGetImageInt(sepadata->bilinvarsmap, (void*)(size_t) mapidx);

   /* if the bilinear expression is not given, either create it or get it from the array */
   if( expr == NULL )
   {
      if( termpos == INT_MAX )
      {
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, sepadata->conshdlr, &(prodexprs[0]), x) );
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, sepadata->conshdlr, &(prodexprs[1]), y) );
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, sepadata->conshdlr, &expr, 2, prodexprs, 1.0) );

         SCIPreleaseConsExprExpr(scip, &(prodexprs[1]));
         SCIPreleaseConsExprExpr(scip, &(prodexprs[0]));

         new = TRUE;
      }
      else
         expr = sepadata->bilinterms[termpos];
   }

   if( termpos == INT_MAX )
   {
      /* store variables if it's the first time they are found in a bilinear term */
      if( !SCIPhashmapExists(varmap, (void*)(size_t) xidx) )
      {
         SCIP_CALL( SCIPhashmapInsertInt(varmap, (void*)(size_t) xidx, sepadata->nbilinvars) ); /*lint !e571*/
         sepadata->varssorted[sepadata->nbilinvars] = x;
         SCIP_CALL( SCIPcaptureVar(scip, x) );
         xpos = sepadata->nbilinvars;
         ++sepadata->nbilinvars;
      }
      else
      {
         xpos = SCIPhashmapGetImageInt(varmap, (void*)(size_t) xidx);
      }
      if( sepadata->nvarbilinvars[xpos] == 0 )
         SCIPallocBlockMemoryArray(scip, &sepadata->varbilinvars[xpos], SCIPgetNVars(scip));
      sepadata->varbilinvars[xpos][sepadata->nvarbilinvars[xpos]] = y;
      ++sepadata->nvarbilinvars[xpos];

      if( !SCIPhashmapExists(varmap, (void*)(size_t) yidx) )
      {
         SCIP_CALL( SCIPhashmapInsertInt(varmap, (void*)(size_t) yidx, sepadata->nbilinvars) ); /*lint !e571*/
         sepadata->varssorted[sepadata->nbilinvars] = y;
         SCIP_CALL( SCIPcaptureVar(scip, y) );
         ypos = sepadata->nbilinvars;
         ++sepadata->nbilinvars;
      }
      else
      {
         ypos = SCIPhashmapGetImageInt(varmap, (void*)(size_t) yidx);
      }
      if( sepadata->nvarbilinvars[ypos] == 0 )
         SCIPallocBlockMemoryArray(scip, &sepadata->varbilinvars[ypos], SCIPgetNVars(scip));
      if( xidx != yidx )
      {
         sepadata->varbilinvars[ypos][sepadata->nvarbilinvars[ypos]] = x;
         ++sepadata->nvarbilinvars[ypos];
      }

      /* insert the position of the linearization expression into auxvar hashmap */
      SCIP_CALL( SCIPhashmapInsertInt(sepadata->bilinvarsmap, (void*)(size_t) mapidx,
                                      sepadata->nbilinterms) ); /*lint !e571*/

      /* add variables and exprs to bilin-arrays and capture them */
      SCIPinfoMessage(scip, NULL, "\nat pos %d adding bilinear term ", sepadata->nbilinterms);
      SCIPprintConsExprExpr(scip, sepadata->conshdlr, expr, NULL);
      sepadata->bilinterms[sepadata->nbilinterms] = expr;
      SCIPcaptureConsExprExpr(expr);
      sepadata->isimplicit[sepadata->nbilinterms] = isimplicit;

      sepadata->nlinexprs[sepadata->nbilinterms] = 0;
      sepadata->slinexprs[sepadata->nbilinterms] = 0;
      sepadata->linexprs[sepadata->nbilinterms] = NULL;
      sepadata->bilinlockspos[sepadata->nbilinterms] = SCIPgetConsExprExprNLocksPos(expr);
      sepadata->bilinlocksneg[sepadata->nbilinterms] = SCIPgetConsExprExprNLocksNeg(expr);

      termpos = sepadata->nbilinterms;
      ++sepadata->nbilinterms;
   }

   /* search for linexpr in the expressions already found for this product */
   if( sepadata->linexprs[termpos] == NULL )
   {
      found = FALSE;
      linpos = 0;
   }
   else
      found = SCIPsortedvecFindPtr((void**)sepadata->linexprs[termpos], SCIPexprsComp, linexpr, sepadata->nlinexprs[termpos], &linpos);

   if( found )
   {
      SCIPinfoMessage(scip, NULL, "\nthis linearisation has already been added");
      assert(SCIPcompareConsExprExprs(sepadata->linexprs[termpos][linpos], linexpr) == 0);
   }
   else
   {  /* linexpr has not been added yet, add it here */
      SCIPinfoMessage(scip, NULL, "\nnew linearisation for termpos %d", termpos);

      SCIPprintConsExprExpr(scip, sepadata->conshdlr, sepadata->bilinterms[sepadata->nbilinterms-1], NULL);
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &(sepadata->linexprs[termpos]), &(sepadata->slinexprs[termpos]),
         sepadata->nlinexprs[termpos] + 1) );

      /* insert expression at the correct position */
      for( int i = sepadata->nlinexprs[termpos]; i > linpos; --i )
      {
         sepadata->linexprs[termpos][i] = sepadata->linexprs[termpos][i-1];
      }
      assert(linexpr != NULL);
      SCIPinfoMessage(scip, NULL, "\nadding in position %d", sepadata->nlinexprs[termpos]);
      sepadata->linexprs[termpos][sepadata->nlinexprs[termpos]] = linexpr;
      ++sepadata->nlinexprs[termpos];
   }

   /* if it is an implicit relation, need to increase the numbers of locks */
   if( isimplicit )
   {
      if( lhs )
         sepadata->bilinlocksneg[termpos]++;
      if( rhs )
         sepadata->bilinlockspos[termpos]++;
   }

   /* add locks to priorities of both variables */
   sepadata->varpriorities[SCIPhashmapGetImageInt(varmap, (void*)(size_t) xidx)] += sepadata->bilinlockspos[termpos] + sepadata->bilinlocksneg[termpos]; /*lint !e571*/
   sepadata->varpriorities[SCIPhashmapGetImageInt(varmap, (void*)(size_t) yidx)] += sepadata->bilinlockspos[termpos] + sepadata->bilinlocksneg[termpos]; /*lint !e571*/

   /* if expr was created here, it was captured then */
   if( new )
      SCIPreleaseConsExprExpr(scip, &expr);

   return SCIP_OKAY;
}


static
SCIP_RETCODE createLinearisation(
   SCIP*                scip,
   SCIP_SEPADATA*       sepadata,
   SCIP_VAR*            w,
   SCIP_VAR*            x,
   SCIP_VAR*            y,
   SCIP_Real*           coefs,
   SCIP_Real            cst,
   SCIP_CONSEXPR_EXPR** linexpr
   )
{
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONSEXPR_EXPR* linchildren[3];
   SCIP_Bool changed, infeasible;

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, sepadata->conshdlr, &linchildren[0], w) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, sepadata->conshdlr, &linchildren[1], x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, sepadata->conshdlr, &linchildren[2], y) );

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, sepadata->conshdlr, linexpr, 3, linchildren, coefs, cst) );

   SCIPinfoMessage(scip, NULL, "\nLinear expression:");
   SCIPprintConsExprExpr(scip, sepadata->conshdlr, *linexpr, NULL);

   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, sepadata->conshdlr, *linexpr, &simplified, &changed, &infeasible) );

   SCIPreleaseConsExprExpr(scip, linexpr);
   *linexpr = simplified;

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &linchildren[2]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &linchildren[1]) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &linchildren[0]) );

   SCIPinfoMessage(scip, NULL, "\nadding to term %s and %s a linearisation: ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPprintConsExprExpr(scip, sepadata->conshdlr, *linexpr, NULL);

   return SCIP_OKAY;
}

/** extract all possible bilinear products from two implied relations
 *
 * Given two linear rows having the same three variables, at least one of which is binary
 * (i.e. implied relations), detect all possible bilinear products encoded in them. A
 * bilinear product can be described by relations that can be written in the following form:
 *
 * a_1w + b_1x + c_1y <=/>= d_1,
 * a_2w + b_2x + c_2y <=/>= d_2,
 *
 * where x is binary, b_1*b_2 < 0 and a_1*a_2 > 0. It is possible for different variables
 * to assume the roles of w, x and y here, which produces different products.
 *
 */
static
SCIP_RETCODE extractImplrelProducts(
   SCIP*          scip,      /**< SCIP data structure */
   SCIP_SEPADATA* sepadata,  /**< separator data */
   SCIP_ROW*      impl_rel1, /**< the first implied relation */
   SCIP_ROW*      impl_rel2, /**< the second implied relation */
   SCIP_HASHMAP*  varmap     /**< variable map */
   )
{
   SCIP_Real* coefs1;
   SCIP_Real* coefs2;
   SCIP_COL** cols1;
   SCIP_COL** cols2;
   int xpos, wpos, ypos, i;
   SCIP_Real sign1, sign2;
   SCIP_Real lhs1, rhs1, lhs2, rhs2;
   SCIP_Real mult;
   SCIP_CONSEXPR_EXPR* bilinexpr;
   SCIP_CONSEXPR_EXPR* linexpr;
   SCIP_Real lincoefs[3];
   SCIP_VAR* w;
   SCIP_VAR* x;
   SCIP_VAR* y;

   /* get row coefficients and columns */
   coefs1 = SCIProwGetVals(impl_rel1);
   coefs2 = SCIProwGetVals(impl_rel2);
   cols1 = SCIProwGetCols(impl_rel1);
   cols2 = SCIProwGetCols(impl_rel2);

   /* an implied relation must have three variables */
   assert(SCIProwGetNNonz(impl_rel1) == 3);
   assert(SCIProwGetNNonz(impl_rel2) == 3);

   /* both rows have the same variables and the same number of columns */
   assert(SCIPcolGetIndex(cols1[0]) == SCIPcolGetIndex(cols2[0]));
   assert(SCIPcolGetIndex(cols1[1]) == SCIPcolGetIndex(cols2[1]));
   assert(SCIPcolGetIndex(cols1[2]) == SCIPcolGetIndex(cols2[2]));

   /* at least one variable must be binary */
   assert(SCIPvarGetType(SCIPcolGetVar(cols1[0])) == SCIP_VARTYPE_BINARY ||
              SCIPvarGetType(SCIPcolGetVar(cols1[1])) == SCIP_VARTYPE_BINARY ||
              SCIPvarGetType(SCIPcolGetVar(cols1[2])) == SCIP_VARTYPE_BINARY);

   /* look for suitable binary variables */
   for( xpos = 0; xpos < 3; ++xpos )
   {
      if( SCIPvarGetType(SCIPcolGetVar(cols1[xpos])) != SCIP_VARTYPE_BINARY )
         continue;

      SCIPdebugMsg(scip, "binary var = %s\n", SCIPvarGetName(SCIPcolGetVar(cols1[xpos])));

      /* we flip the 1st row so that coef of x is positive */
      sign1 = coefs1[xpos] > 0 ? 1.0 : -1.0;

      /* the indicator var's coefficients in the two rows should have different signs */
      if( coefs1[xpos]*coefs2[xpos] > 0 )
      { /* take 2 combinations of lhs and rhs */
         sign2 = -sign1;
      }
      else
      { /* take 2 lhs and 2 rhs, keep the second row as it is */
         sign2 = sign1;
      }

      /* get lhs and rhs depending on whether we flip the rows */
      if( sign1 > 0 )
      {
         lhs1 = SCIProwGetLhs(impl_rel1);
         rhs1 = SCIProwGetRhs(impl_rel1);
      }
      else
      {
         lhs1 = -SCIProwGetRhs(impl_rel1);
         rhs1 = -SCIProwGetLhs(impl_rel1);
      }

      if( sign2 > 0 )
      {
         lhs2 = SCIProwGetLhs(impl_rel2);
         rhs2 = SCIProwGetRhs(impl_rel2);
      }
      else
      {
         lhs2 = -SCIProwGetRhs(impl_rel2);
         rhs2 = -SCIProwGetLhs(impl_rel2);
      }
      SCIPdebugMsg(scip, "\nsigns: %f, %f", sign1, sign2);
      /* from here on, we consider only the flipped (by multiplying by signi) rows */

      /* try two different combinations of w and y */
      for( i = 1; i <= 2; ++i )
      {
         wpos = (xpos + i) % 3;

         /* coefficients of w should have the same sign in the two imol_rels */
         if( coefs1[wpos]*sign1*coefs2[wpos]*sign2 < 0 )
            continue;

         ypos = (xpos - i + 3) % 3;

         /* all conditions satisfied, we can extract the product */

         w = SCIPcolGetVar(cols1[wpos]);
         x = SCIPcolGetVar(cols1[xpos]);
         y = SCIPcolGetVar(cols1[ypos]);

         /* do lhs */
         if( lhs1 != -SCIPinfinity(scip) && lhs2 != -SCIPinfinity(scip) )
         {
            mult = 1/(coefs2[wpos]*sign2*coefs1[ypos]*sign1 - coefs2[ypos]*sign2*coefs1[wpos]*sign1);
            assert(mult != 0.0); /* TODO with some tolerance? */

            SCIPinfoMessage(scip, NULL, "\nLhs, found suitable implied rels (w,x,y): %f%s + %f%s + %f%s %s %f\n",
                            sign1*coefs1[wpos], SCIPvarGetName(w), sign1*coefs1[xpos], SCIPvarGetName(x),
                            sign1*coefs1[ypos], SCIPvarGetName(y), sign1 > 0 ? ">=" : "<=", lhs1);

            SCIPinfoMessage(scip, NULL, "\nand %f%s + %f%s + %f%s %s %f\n",
                            sign2*coefs2[wpos], SCIPvarGetName(w), sign2*coefs2[xpos], SCIPvarGetName(x),
                            sign2*coefs2[ypos], SCIPvarGetName(y), sign2 > 0 ? ">=" : "<=", lhs2);

            lincoefs[0] = sign2*coefs2[wpos]*sign1*coefs1[wpos]*mult;
            lincoefs[1] = (sign2*coefs2[xpos]*sign1*coefs1[wpos] - lhs2*sign1*coefs1[wpos] + lhs1*sign2*coefs2[wpos] )*mult;
            lincoefs[2] = sign2*coefs2[wpos]*sign1*coefs1[ypos]*mult;

            SCIPinfoMessage(scip, NULL, "\nproduct: %s%s %s %f%s + %f%s + %f%s + %f", SCIPvarGetName(x), SCIPvarGetName(y),
               sign2*coefs2[wpos]*mult > 0 ? "<=" : ">=", lincoefs[0], SCIPvarGetName(w), lincoefs[1], SCIPvarGetName(x), lincoefs[2],
               SCIPvarGetName(y), -coefs2[wpos]*lhs1*mult);

            SCIP_CALL( createLinearisation(scip, sepadata, w, x, y, lincoefs, -sign2*coefs2[wpos]*lhs1*mult, &linexpr) );
            SCIP_CALL( addBilinProduct(scip, sepadata, NULL, x, y, linexpr, TRUE, sign2*coefs2[wpos]*mult > 0, sign2*coefs2[wpos]*mult < 0, varmap) );
         }

         /* do rhs */
         if( rhs1 != SCIPinfinity(scip) && rhs2 != SCIPinfinity(scip) )
         {
            mult = 1/(coefs1[wpos]*sign1*coefs2[ypos]*sign2 - coefs1[ypos]*sign1*coefs2[wpos]*sign2);

            SCIPinfoMessage(scip, NULL, "\nRhs, found suitable implied rels (w,x,y): %f%s + %f%s + %f%s %s %f\n",
                            sign1*coefs1[wpos], SCIPvarGetName(w), sign1*coefs1[xpos], SCIPvarGetName(x),
                            sign1*coefs1[ypos], SCIPvarGetName(y), sign1 > 0 ? "<=" : ">=", rhs1);

            SCIPinfoMessage(scip, NULL, "\nand %f%s + %f%s + %f%s %s %f\n",
                            sign2*coefs2[wpos], SCIPvarGetName(w), sign2*coefs2[xpos], SCIPvarGetName(x),
                            sign2*coefs2[ypos], SCIPvarGetName(y), sign2 > 0 ? "<=" : ">=", rhs2);

            lincoefs[0] = sign1*coefs1[wpos]*sign2*coefs2[wpos]*mult;
            lincoefs[1] = (sign1*coefs1[xpos]*sign2*coefs2[wpos] - rhs1*sign2*coefs2[wpos] + rhs2*sign1*coefs1[wpos] )*mult;
            lincoefs[2] = sign1*coefs1[wpos]*sign2*coefs2[ypos]*mult;

            SCIPinfoMessage(scip, NULL, "\nproduct: %s%s %s %f%s + %f%s + %f%s + %f", SCIPvarGetName(x), SCIPvarGetName(y),
               sign1*coefs1[wpos]*mult > 0 ? ">=" : "<=", lincoefs[0], SCIPvarGetName(w), lincoefs[1], SCIPvarGetName(x), lincoefs[2],
               SCIPvarGetName(y), -coefs1[wpos]*rhs2*mult);

            SCIP_CALL( createLinearisation(scip, sepadata, w, x, y, lincoefs, -sign1*coefs1[wpos]*rhs2*mult, &linexpr) );
            SCIP_CALL( addBilinProduct(scip, sepadata, bilinexpr, x, y, linexpr, TRUE, sign1*coefs1[wpos]*mult < 0, sign1*coefs1[wpos]*mult > 0, varmap) );
         }
      }
   }

   return SCIP_OKAY;
}


/* detect bilinear products encoded in linear constraints */
static
SCIP_RETCODE detectHiddenProducts(
   SCIP*          scip,       /**< SCIP data structure */
   SCIP_SEPADATA* sepadata,   /**< separation data */
   SCIP_HASHMAP*  varmap      /**< variable map */
   )
{
   int v, r, r2, i, nrows;
   SCIP_ROW** prob_rows;
   SCIP_HASHTABLE* hashdatatable;
   HASHDATA hashdata;
   HASHDATA* foundhashdata;
   SCIP_COL** cols;

   /* get the rows, depending on settings */
   if( sepadata->isinitialround || sepadata->onlyinitial )
   {
      SCIP_CALL( getInitialRows(scip, &prob_rows, &nrows) );
   }
   else
   {
      SCIP_CALL( SCIPgetLPRowsData(scip, &prob_rows, &nrows) );
   }

   /* create a table of implied relations */
   SCIP_CALL( SCIPhashtableCreate(&hashdatatable, SCIPblkmem(scip), nrows, SCIPhashGetKeyStandard,
      hashdataKeyEqConss, hashdataKeyValConss, (void*) scip) );

   /* look for implied relations */
   for( r = 0; r < nrows; ++r ) /* go through rows */
   {
      assert(prob_rows[r] != NULL);

      cols = SCIProwGetCols(prob_rows[r]);
      SCIPinfoMessage(scip, NULL, "\nrows %s:", SCIProwGetName(prob_rows[r]));
      for( v = 0; v < SCIProwGetNNonz(prob_rows[r]); ++v )
         SCIPinfoMessage(scip, NULL, "%s(%d) ", SCIPvarGetName(SCIPcolGetVar(cols[v])), SCIPcolGetIndex(cols[v]));

      /* look for rows with 3 variables */
      if( SCIProwGetNNonz(prob_rows[r]) != 3 )
         continue;

      cols = SCIProwGetCols(prob_rows[r]);

      /* an implied relation contains at least one binary variable */
      if( SCIPvarGetType(SCIPcolGetVar(cols[0])) != SCIP_VARTYPE_BINARY
            && SCIPvarGetType(SCIPcolGetVar(cols[1])) != SCIP_VARTYPE_BINARY
            && SCIPvarGetType(SCIPcolGetVar(cols[2])) != SCIP_VARTYPE_BINARY)
         continue;

      /* fill in hashdata */
      hashdata.nvars = 3;
      hashdata.rows = NULL;
      SCIP_CALL( SCIPallocBufferArray(scip, &(hashdata.vars), 3) );
      for( v = 0; v < 3; ++v )
         hashdata.vars[v] = SCIPcolGetVar(cols[v]);

      /* get the element corresponsing to the three variables */
      foundhashdata = (HASHDATA*)SCIPhashtableRetrieve(hashdatatable, &hashdata);

      if( foundhashdata != NULL )
      {
         /* if element exists, update it */
         foundhashdata->rows[foundhashdata->nrows] = prob_rows[r];
         ++foundhashdata->nrows;
         SCIPfreeBufferArray(scip, &(hashdata.vars));
      }
      else
      {
         /* create an element for the combination of three variables */
         SCIP_CALL( SCIPallocBuffer(scip, &foundhashdata) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(foundhashdata->rows), nrows) );

         foundhashdata->rows[0] = prob_rows[r];
         foundhashdata->nvars = 3;
         foundhashdata->nrows = 1;
         foundhashdata->vars = hashdata.vars;

         SCIPhashtableInsert(hashdatatable, (void*)foundhashdata);
      }
   }


   SCIPinfoMessage(scip, NULL, "\n");
   SCIPdebugMsg(scip, "Implied relations table:\n");
   /* go through all sets of three variables */
   for( i = 0; i < SCIPhashtableGetNEntries(hashdatatable); ++i )
   {
      foundhashdata = (HASHDATA*)SCIPhashtableGetEntry(hashdatatable, i);
      if( foundhashdata == NULL )
         continue;

      SCIPdebugMsg(scip, "(%s, %s, %s): ", SCIPvarGetName(foundhashdata->vars[0]),
         SCIPvarGetName(foundhashdata->vars[1]), SCIPvarGetName(foundhashdata->vars[2]));

      /* go through implied relations for the corresponsing three variables */
      for( r = 0; r < foundhashdata->nrows; ++r )
      {
         SCIPinfoMessage(scip, NULL, "%s; ", SCIProwGetName(foundhashdata->rows[r]));

         /* go through the remaining rows in this entry */
         for( r2 = r + 1; r2 < foundhashdata->nrows; ++r2 )
         {
            SCIP_CALL( extractImplrelProducts(scip, sepadata, foundhashdata->rows[r], foundhashdata->rows[r2], varmap) );
         }

         /* TODO go through the implied bounds of the binary variable */

         /* TODO for each var that can be yij, check global bound */



      }

      SCIPfreeBufferArray(scip, &(foundhashdata->vars));
      SCIPfreeBufferArray(scip, &(foundhashdata->rows));
      SCIPfreeBuffer(scip, &foundhashdata);
   }

   SCIPhashtableFree(&hashdatatable);

   if( sepadata->isinitialround || sepadata->onlyinitial )
      SCIPfreeBufferArray(scip, &prob_rows);
   /* no need to free prob_rows in the other case since SCIPgetLPRowsData does not allocate it */

   return SCIP_OKAY;
}


/** helper method to create separation data */
static
SCIP_RETCODE createSepaData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separation data */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPR_EXPR* linexpr;
   SCIP_HASHMAP* varmap;
   SCIP_VAR* x;
   SCIP_VAR* y;
   int xidx;
   int yidx;
   int i;
   int nconss;
   int nvars;
   int maxidx = 0;

   assert(sepadata != NULL);

   sepadata->nbilinvars = 0;
   sepadata->nbilinterms = 0;

   conss = SCIPconshdlrGetConss(sepadata->conshdlr);
   nconss = SCIPconshdlrGetNConss(sepadata->conshdlr);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* create variable map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* create iterator */
   SCIP_CALL( SCIPexpriteratorCreate(&it, sepadata->conshdlr, SCIPblkmem(scip)) );

   /* create the empty map for bilinear terms */
   SCIP_CALL( SCIPhashmapCreate(&sepadata->bilinvarsmap, SCIPblkmem(scip), nvars) );

   /* allocate memory for arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->varssorted, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &sepadata->varbilinvars, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &sepadata->nvarbilinvars, nvars) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &sepadata->varpriorities, nvars) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->linexprs, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->nlinexprs, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->slinexprs, nvars) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->bilinterms, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->bilinlockspos, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->bilinlocksneg, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sepadata->isimplicit, nvars) );

   /* find maximum variable index */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
      maxidx = MAX(maxidx, SCIPvarGetIndex(vars[i]));  /*lint !e666*/
   sepadata->maxvarindex = maxidx; /* TODO don't we need +1 here? */

   SCIP_CALL( detectHiddenProducts(scip, sepadata, varmap) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSEXPR_EXPR* expr;

      SCIP_CALL(SCIPexpriteratorInit(it, SCIPgetExprConsExpr(scip, conss[i]), SCIP_CONSEXPRITERATOR_DFS, TRUE));
      SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR);
      expr = SCIPexpriteratorGetCurrent(it);

      /* collect bilinear terms */
      while( !SCIPexpriteratorIsEnd(it) ) /*lint !e441*/
      {
         switch( SCIPexpriteratorGetStageDFS(it) )
         {
            case SCIP_CONSEXPRITERATOR_ENTEREXPR:
            {
               SCIP_VAR* auxvar;
               int mapidx;

               assert(expr != NULL);

               auxvar = SCIPgetConsExprExprAuxVar(expr);

               /* no linearization variable available */
               if( auxvar == NULL )
               {
                  SCIPdebugMsg(scip, "auxvar = NULL, ignoring the expression ");
#ifdef SCIP_DEBUG
                  SCIPprintConsExprExpr(scip, sepadata->conshdlr, expr, NULL);
#endif
                  break;
               }

               x = NULL;
               y = NULL;

               /* test if expression is quadratic */
               if( SCIPgetConsExprExprHdlr(expr) == SCIPfindConsExprExprHdlr(sepadata->conshdlr, "pow")
                  && SCIPgetConsExprExprPowExponent(expr) == 2.0 )
               {
                  /* if only initial rows are requested, skip products of non-variable expressions */
                  if( sepadata->onlyinitial && !SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(expr)[0]) )
                     break;

                  x = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]);
                  y = x;
               }
                  /* test if expression is bilinear */
               else if( SCIPgetConsExprExprHdlr(expr) == SCIPfindConsExprExprHdlr(sepadata->conshdlr, "prod")
                  && SCIPgetConsExprExprNChildren(expr) == 2 )
               {
                  /* if only initial rows are requested, skip products of non-variable expressions */
                  if( sepadata->onlyinitial && (!SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(expr)[0])
                                                || !SCIPisConsExprExprVar(SCIPgetConsExprExprChildren(expr)[1])) )
                     break;

                  x = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]);
                  y = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[1]);
               }

               /* if children don't have linearization variables, there's nothing to do */
               if( x != NULL && y != NULL )
               {
                  /* switch variables if necessary */
                  if( SCIPvarComp(x, y) > 0 )
                     SCIPswapPointers((void**)&x, (void**)&y);

                  SCIPdebugMsg(scip, "found a product of vars %s and %s\n", SCIPvarGetName(x), SCIPvarGetName(y));

                  assert(auxvar != NULL);

                  xidx = SCIPvarGetIndex(x);
                  yidx = SCIPvarGetIndex(y);

                  /* compute unique index of the bilinear term */
                  mapidx = xidx * sepadata->maxvarindex + yidx;

                  /*  */
                  if( !SCIPhashmapExists(sepadata->bilinvarsmap, (void*)(size_t) mapidx) )
                  {
                     SCIP_CALL( SCIPcreateConsExprExprVar(scip, sepadata->conshdlr, &linexpr, auxvar) );
                     SCIP_CALL( addBilinProduct(scip, sepadata, expr, x, y, linexpr, FALSE, SCIPgetConsExprExprNLocksNeg(expr) > 0,
                        SCIPgetConsExprExprNLocksPos(expr) > 0, varmap) );
                  }
                  else
                  {
                     expr = SCIPexpriteratorSkipDFS(it);
                     continue;
                  }
               }

               break;
            }

            default:
               SCIPABORT();
               break;
         }

         expr = SCIPexpriteratorGetNext(it);
      }
   }

   /* reallocate arrays to fit actually sizes */
   if( sepadata->nbilinvars < nvars )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->varssorted, nvars, sepadata->nbilinvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->varbilinvars, nvars, sepadata->nbilinvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->nvarbilinvars, nvars, sepadata->nbilinvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->varpriorities, nvars, sepadata->nbilinvars) );
   }

   if( sepadata->nbilinterms < nvars )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->linexprs, nvars, sepadata->nbilinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->nlinexprs, nvars, sepadata->nbilinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->slinexprs, nvars, sepadata->nbilinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->bilinterms, nvars, sepadata->nbilinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->bilinlockspos, nvars, sepadata->nbilinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->bilinlocksneg, nvars, sepadata->nbilinterms) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->isimplicit, nvars, sepadata->nbilinterms) );
   }

   for( i = 0; i < sepadata->nbilinvars; ++i )
   {
      if( sepadata->nvarbilinvars[i] > 0 && sepadata->nvarbilinvars[i] < nvars )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &sepadata->varbilinvars[i], nvars, sepadata->nvarbilinvars[i]) );
      }
   }

   /* sort maxnumber of variables according to their occurrences */
   SCIPselectDownIntIntPtrPtr(sepadata->varpriorities, sepadata->nvarbilinvars, (void**) sepadata->varssorted,
      (void**) sepadata->varbilinvars, sepadata->maxusedvars, sepadata->nbilinvars);

   SCIPexpriteratorFree(&it);
   SCIPhashmapFree(&varmap);

   sepadata->iscreated = TRUE;
   sepadata->isinitialround = TRUE;

   if( sepadata->nbilinterms > 0 )
      SCIPinfoMessage(scip, NULL, "\nFound bilinear terms\n");
   else
      SCIPinfoMessage(scip, NULL, "\nNo bilinear terms");

   return SCIP_OKAY;
}


/** helper method to get the position of linearization terms of a bilinear term xy
 *
 *  @return -1 if no linearization variable exists
 */
static
int getBilinPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separation data */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y                   /**< second variable */
   )
{
   int idx;
   int img;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(x != NULL);
   assert(y != NULL);

   /* it seems that x or y have been added after initsol -> no linearization variable available */
   if( SCIPvarGetIndex(x) > sepadata->maxvarindex || SCIPvarGetIndex(y) > sepadata->maxvarindex )
   {
      return -1;
   }

   /* switch variables if necessary */
   if( x != y && SCIPvarComp(x, y) > 0 )
      SCIPswapPointers((void**) &x, (void**) &y);

   /* compute unique index of the bilinear term */
   idx = SCIPvarGetIndex(x) * sepadata->maxvarindex + SCIPvarGetIndex(y);

   if( SCIPhashmapExists(sepadata->bilinvarsmap, (void*)(size_t) idx) )
   {
      img = SCIPhashmapGetImageInt(sepadata->bilinvarsmap, (void*)(size_t) idx); /*lint !e571*/
      return img;
   }

   return -1;
}

/** tests if a row contains too many unknown bilinear terms w.r.t. the parameters */
static
SCIP_RETCODE isAcceptableRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separation data */
   SCIP_ROW*             row,                /**< the row to be tested */
   SCIP_VAR*             var,                /**< the variable that is to be multiplied with row */
   int                   nlocks,             /**< the number of locks of the variable */
   SCIP_Bool*            acceptable          /**< buffer to store the result */
   )
{
   SCIP_VAR* linvar;
   int i;
   int nterms = 0;
   int linpos;

   assert(row != NULL);
   assert(var != NULL);

   /* test if the ratio of non-zeroes and known terms of this variable is ok */
   if( SCIProwGetNNonz(row) * sepadata->maxnonzeroprop > nlocks )
   {
      *acceptable = FALSE;
      return SCIP_OKAY;
   }

   for( i = 0; (i < SCIProwGetNNonz(row)) && (sepadata->maxunknownterms >= 0 || nterms <= sepadata->maxunknownterms); ++i )
   {
      linpos = getBilinPos(scip, sepadata, var, SCIPcolGetVar(SCIProwGetCols(row)[i]) );

      if( linpos == -1 )
         ++nterms;
   }

   sepadata->currentnunknown = nterms;

   *acceptable = nterms <= sepadata->maxunknownterms;

   return SCIP_OKAY;
}

/* add a linearisation of term coef*colvar*var to cut
 *
 * adds the linear term involving colvar to cut and updates coefvar and finalside
 */
static
SCIP_RETCODE addRltTerm(
   SCIP*          scip,         /**< SCIP data structure */
   SCIP_SEPADATA* sepadata,     /**< separator data */
   SCIP_SOL*      sol,          /**< the point to be separated (can be NULL) */
   SCIP_ROW*      cut,          /**< cut to which the term is to be added */
   SCIP_VAR*      var,          /**< multiplier variable */
   SCIP_VAR*      colvar,       /**< row variable to be multiplied */
   SCIP_Real      coef,         /**< coefficient of the bilinear term */
   SCIP_Bool      uselb,        /**< whether we multiply with (var - lb) or (ub - var) */
   SCIP_Bool      uselhs,       /**< whether to create a cut for the lhs or rhs */
   SCIP_Bool      local,        /**< whether local or global cuts should be computed */
   SCIP_Bool      computeEqCut, /**< whether conditions are fulfilled to compute equality cuts */
   SCIP_Real*     coefvar,      /**< coefficient of var */
   SCIP_Real*     finalside,    /**< buffer to store the left or right hand side of cut */
   SCIP_Bool*     success       /**< buffer to store whether cut was created successfully */
   )
{
   SCIP_VAR* auxvar;
   SCIP_Real lbvar, ubvar;
   SCIP_Real refpointvar;
   SCIP_Real signfactor, boundfactor;
   SCIP_Real coefauxvar, coefcolvar;
   int auxpos;

   auxpos = getBilinPos(scip, sepadata, var, colvar);
   auxvar = auxpos == -1 ? NULL : SCIPgetConsExprExprVarVar(sepadata->linexprs[auxpos][0]);

   lbvar = local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
   ubvar = local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

   refpointvar = MAX(lbvar, MIN(ubvar, SCIPgetSolVal(scip, sol, var))); /*lint !e666*/

   signfactor = (uselb ? 1.0 : -1.0);
   boundfactor = (uselb ? -lbvar : ubvar);

   coefauxvar = coef * signfactor;
   coefcolvar = coef * boundfactor;

   /* if the auxiliary variable for this term exists, simply add it to the cut with the previous coefficient */
   if( auxvar != NULL )
   {
      SCIPdebugMsg(scip, "auxvar for %s and %s found, will be added to cut\n", SCIPvarGetName(colvar), SCIPvarGetName(var));
      assert(!SCIPisInfinity(scip, REALABS(coefauxvar)));
      SCIP_CALL( SCIPaddVarToRow(scip, cut, auxvar, coefauxvar) );
   }

   /* otherwise, use the McCormick estimator in place of the bilinear term */
   else if( colvar != var )
   {
      SCIP_Bool found_clique = FALSE;
      SCIP_Real lbcolvar = local ? SCIPvarGetLbLocal(colvar) : SCIPvarGetLbGlobal(colvar);
      SCIP_Real ubcolvar = local ? SCIPvarGetUbLocal(colvar) : SCIPvarGetUbGlobal(colvar);
      SCIP_Real refpointcolvar = MAX(lbcolvar, MIN(ubcolvar, SCIPgetSolVal(scip, sol, colvar))); /*lint !e666*/

      assert(!computeEqCut);

      if( REALABS(lbcolvar) > MAXVARBOUND || REALABS(ubcolvar) > MAXVARBOUND )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, "auxvar for %s and %s not found, will linearise the product\n", SCIPvarGetName(colvar), SCIPvarGetName(var));

      /* if both variables are binary. check if they are contained together in some clique */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY &&  SCIPvarGetType(colvar) == SCIP_VARTYPE_BINARY )
      {
         int c;
         SCIP_CLIQUE** varcliques;

         varcliques = SCIPvarGetCliques(var, TRUE);

         /* look through cliques containing var */
         for( c = 0; c < SCIPvarGetNCliques(var, TRUE); ++c )
         {
            if( SCIPcliqueHasVar(varcliques[c], colvar, TRUE) ) /* var + colvar <= 1 => var*colvar = 0 */
            {
               /* product is zero, add nothing */
               found_clique = TRUE;
               break;
            }

            if( SCIPcliqueHasVar(varcliques[c], colvar, FALSE) ) /* var + (1-colvar) <= 1 => var*colvar = var */
            {
               *coefvar += coefauxvar;
               found_clique = TRUE;
               break;
            }
         }

         if(!found_clique)
         {
            varcliques = SCIPvarGetCliques(var, FALSE);

            /* look through cliques containing complement of var */
            for( c = 0; c < SCIPvarGetNCliques(var, FALSE); ++c )
            {
               if( SCIPcliqueHasVar(varcliques[c], colvar, TRUE) ) /* (1-var) + colvar <= 1 => var*colvar = colvar */
               {
                  coefcolvar += coefauxvar;
                  found_clique = TRUE;
                  break;
               }

               if( SCIPcliqueHasVar(varcliques[c], colvar, FALSE) ) /* (1-var) + (1-colvar) <= 1 => var*colvar = var + colvar - 1 */
               {
                  *coefvar += coefauxvar;
                  coefcolvar += coefauxvar;
                  *finalside -= coefauxvar;
                  found_clique = TRUE;
                  break;
               }
            }
         }
      }

      if( !found_clique )
      {
         SCIPdebugMsg(scip, "clique for %s and %s not found or at least one of them is not binary, will use McCormick\n", SCIPvarGetName(colvar), SCIPvarGetName(var));
         SCIPaddBilinMcCormick(scip, coefauxvar, lbvar, ubvar, refpointvar, lbcolvar,
                                ubcolvar, refpointcolvar, uselhs, coefvar, &coefcolvar, finalside, success);
         if( !*success )
            return SCIP_OKAY;
      }
   }

   /* or, if it's a quadratic term, use a secant for overestimation and a gradient for underestimation */
   else
   {
      SCIPdebugMsg(scip, "auxvar for %s^2 not found, will use gradient and secant estimators\n", SCIPvarGetName(colvar));

      assert(!computeEqCut);

      /* for a binary var, var^2 = var */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         *coefvar += coefauxvar;
      }
      else
      {
         /* depending on over-/underestimation and the sign of the column variable, compute secant or tangent */
         if( (uselhs && coefauxvar > 0.0) || (!uselhs && coefauxvar < 0.0) )
            SCIPaddSquareSecant(scip, coefauxvar, lbvar, ubvar, refpointvar, coefvar, finalside, success);
         else
            SCIPaddSquareLinearization(scip, coefauxvar, refpointvar, SCIPvarIsIntegral(var), coefvar, finalside, success);

         if( !*success )
            return SCIP_OKAY;
      }
   }

   /* add the linear term for this column */
   if( colvar != var )
   {
      assert(!SCIPisInfinity(scip, REALABS(coefcolvar)));
      SCIP_CALL( SCIPaddVarToRow(scip, cut, colvar, coefcolvar) );
   }
   else
      *coefvar += coefcolvar;
}

/** creates the RLT-cuts formed by multiplying a given row with (x - lb) or (ub - x)
 *
 * in detail:
 * -The row is multiplied either with (x - lb(x)) or with (ub(x) - x), depending on parameter uselb.
 * -The cut is computed either for lhs or rhs, depending on parameter uselhs.
 * -Terms for which no auxiliary variable exists are replaced by either McCormick, secants, or linearization cuts
 */
static
SCIP_RETCODE computeRltCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separation data */
   SCIP_ROW**            cut,                /**< buffer to store the cut */
   SCIP_ROW*             row,                /**< the row that is used for the rlt cuts */
   SCIP_SOL*             sol,                /**< the point to be separated (can be NULL) */
   SCIP_VAR*             var,                /**< the variable that is used for the rlt cuts */
   SCIP_Bool*            success,            /**< buffer to store whether cut was created successfully */
   SCIP_Bool             uselb,              /**< whether we multiply with (var - lb) or (ub - var) */
   SCIP_Bool             uselhs,             /**< whether to create a cut for the lhs or rhs */
   SCIP_Bool             local,              /**< whether local or global cuts should be computed */
   SCIP_Bool             computeEqCut        /**< whether conditions are fulfilled to compute equality cuts */
   )
{
   SCIP_Real signfactor;
   SCIP_Real boundfactor;
   SCIP_Real lbvar;
   SCIP_Real ubvar;
   SCIP_Real coefvar;
   SCIP_Real constside;
   SCIP_Real finalside;
   int i;
   char cutname[SCIP_MAXSTRLEN];

   /* create cut name */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s%d_%d%s%s", "rlt_cut", SCIProwGetIndex(row), SCIPvarGetIndex(var), uselb ? "l" : "r", uselhs ? "l" : "r");

   assert(sepadata != NULL);
   assert(cut != NULL);
   assert(row != NULL);
   assert(var != NULL);
   assert(success != NULL);
   assert(!computeEqCut || SCIPisEQ(scip, SCIProwGetLhs(row), SCIProwGetRhs(row)));

   *cut = NULL;

   /* get data for given variable */
   lbvar = local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
   ubvar = local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
   constside = uselhs ? SCIProwGetLhs(row) : SCIProwGetRhs(row);

   /* if the bounds are too large or the respective side is infinity, skip this cut */
   if( REALABS(lbvar) > MAXVARBOUND || REALABS(ubvar) > MAXVARBOUND || SCIPisInfinity(scip, REALABS(constside)) )
   {
      SCIPdebugMsg(scip, "cut generation for row %s, %s and variable %s with its %s %f not possible\n",
         SCIProwGetName(row), uselhs ? "lhs" : "rhs", SCIPvarGetName(var),
         uselb ? "lower bound" : "upper bound", uselb ? lbvar : ubvar);

      if( REALABS(lbvar) > MAXVARBOUND )
         SCIPdebugMsg(scip, " because of lower bound\n");
      if( REALABS(ubvar) > MAXVARBOUND )
         SCIPdebugMsg(scip, " because of upper bound\n");
      if( SCIPisInfinity(scip, REALABS(constside)) )
         SCIPdebugMsg(scip, " because of side %f\n", constside);

      *success = FALSE;
      return SCIP_OKAY;
   }

   /* initialize some factors needed for computation */
   coefvar = 0.0;
   finalside = 0.0;
   signfactor = (uselb ? 1.0 : -1.0);
   boundfactor = (uselb ? -lbvar : ubvar);

   *success = TRUE;

   /* create an empty row which we then fill with variables step by step */
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, cut, sepa, cutname, -SCIPinfinity(scip), SCIPinfinity(scip),
         TRUE, FALSE, FALSE) );

   /* iterate over all variables in the row and add the corresponding terms to the cuts */
   for( i = 0; i < SCIProwGetNNonz(row); ++i )
   {
      SCIP_VAR* colvar;
      colvar = SCIPcolGetVar(SCIProwGetCols(row)[i]);
      addRltTerm(scip, sepadata, sol, *cut, var, colvar, SCIProwGetVals(row)[i], uselb, uselhs, local, computeEqCut,
         &coefvar, &finalside, success);
   }

   if( REALABS(finalside) > MAXVARBOUND )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* multiply (x-lb) or (ub -x) with the lhs and rhs of the row */
   coefvar += signfactor * (SCIProwGetConstant(row) - constside);
   finalside = boundfactor * (constside - SCIProwGetConstant(row)) - finalside;

   /* set the coefficient of var and the constant side */
   assert(!SCIPisInfinity(scip, REALABS(coefvar)));
   SCIP_CALL( SCIPaddVarToRow(scip, *cut, var, coefvar) );

   assert(!SCIPisInfinity(scip, REALABS(finalside)));
   if( uselhs || computeEqCut )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, *cut, finalside) );
   }
   if( !uselhs || computeEqCut )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, *cut, finalside) );
   }

   SCIPdebugMsg(scip, "cut was generated successfully:\n");
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, *cut, NULL) );
#endif

   return SCIP_OKAY;
}

/** creates the RLT cuts formed by multiplying a given projected row with (x - lb) or (ub - x)
 *
 * in detail:
 * -The row is multiplied either with (x - lb(x)) or with (ub(x) - x), depending on parameter uselb.
 * -The cut is computed either for lhs or rhs, depending on parameter uselhs.
 * -Terms for which no auxiliary variable exists are replaced by either McCormick, secants, or linearization cuts
 */
static
SCIP_RETCODE computeProjRltCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separation data */
   SCIP_ROW**            cut,                /**< buffer to store the cut */
   PROJLP*               projlp,             /**< projected lp */
   int                   idx,                /**< index of the row that is used for the rlt cut */
   SCIP_SOL*             sol,                /**< the point to be separated (can be NULL) */
   SCIP_VAR*             var,                /**< the variable that is used for the rlt cuts */
   SCIP_Bool*            success,            /**< buffer to store whether cut was created successfully */
   SCIP_Bool             uselb,              /**< whether we multiply with (var - lb) or (ub - var) */
   SCIP_Bool             uselhs,             /**< whether to create a cut for the lhs or rhs */
   SCIP_Bool             local,              /**< whether local or global cuts should be computed */
   SCIP_Bool             computeEqCut        /**< whether conditions are fulfilled to compute equality cuts */
)
{
   SCIP_Real signfactor;
   SCIP_Real boundfactor;
   SCIP_Real lbvar;
   SCIP_Real ubvar;
   SCIP_Real coefvar;
   SCIP_Real constside;
   SCIP_Real finalside;
   int i;

   assert(sepadata != NULL);
   assert(cut != NULL);
   assert(projlp != NULL);
   assert(var != NULL);
   assert(success != NULL);
   assert(!computeEqCut || SCIPisEQ(scip, projlp->lhss[idx], projlp->rhss[idx]));

   *cut = NULL;

   /* get data for given variable */
   lbvar = local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
   ubvar = local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
   constside = uselhs ? projlp->lhss[idx] : projlp->rhss[idx];

   /* if the bounds are too large or the respective side is infinity, skip this cut */
   if( REALABS(lbvar) > MAXVARBOUND || REALABS(ubvar) > MAXVARBOUND || SCIPisInfinity(scip, REALABS(constside)) )
   {
      SCIPdebugMsg(scip, "cut generation for projected row %d, %s and variable %s with its %s %f not possible\n",
                   idx, uselhs ? "lhs" : "rhs", SCIPvarGetName(var),
                   uselb ? "lower bound" : "upper bound", uselb ? lbvar : ubvar);

      *success = FALSE;
      return SCIP_OKAY;
   }

   /* initialize some factors needed for computation */
   coefvar = 0.0;
   finalside = 0.0;
   signfactor = (uselb ? 1.0 : -1.0);
   boundfactor = (uselb ? -lbvar : ubvar);

   *success = TRUE;

   /* create an empty row which we then fill with variables step by step */
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, cut, sepa, "rlt_cut", -SCIPinfinity(scip), SCIPinfinity(scip),
                                     TRUE, FALSE, FALSE) );

   /* iterate over all variables in the row and add the corresponding terms to the cuts */
   for( i = 0; i < projlp->nNonz[idx]; ++i )
   {
      addRltTerm(scip, sepadata, sol, *cut, var, projlp->vars[idx][i], projlp->coefs[idx][i], uselb, uselhs,
         local, computeEqCut, &coefvar, &finalside, success);
   }

   if( REALABS(finalside) > MAXVARBOUND )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* multiply (x-lb) or (ub -x) with the lhs and rhs of the row */
   coefvar += signfactor * (projlp->consts[idx] - constside);
   finalside = boundfactor * (constside - projlp->consts[idx]) - finalside;

   /* set the coefficient of var and the constant side */
   assert(!SCIPisInfinity(scip, REALABS(coefvar)));
   SCIP_CALL( SCIPaddVarToRow(scip, *cut, var, coefvar) );

   assert(!SCIPisInfinity(scip, REALABS(finalside)));
   if( uselhs || computeEqCut )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, *cut, finalside) );
   }
   if( !uselhs || computeEqCut )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, *cut, finalside) );
   }

   SCIPdebugMsg(scip, "projected cut was generated successfully:\n");
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, *cut, NULL) );
#endif

   return SCIP_OKAY;
}

/** creates the projected problem
 *
 *  All variables that are at their bounds at the current solution are added
 *  to left and/or right hand sides as constant values.
 */
static
SCIP_RETCODE createProjLP(
   SCIP*            scip,       /**< SCIP data structure */
   SCIP_ROW**       rows,       /**< problem rows */
   int              nrows,      /**< number of rows */
   SCIP_SOL*        sol,        /**< the point to be separated (can be NULL) */
   PROJLP**         projlp,     /**< the projected problem data structure */
   SCIP_Bool        local       /**< are local cuts allowed? */
   )
{
   SCIP_COL** cols;
   int i, v;
   SCIP_VAR* var;
   SCIP_Real val, vlb, vub;

   assert(scip != NULL);
   assert(rows != NULL);
   assert(projlp != NULL);

   SCIP_CALL( SCIPallocBuffer(scip, projlp) );

   SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->coefs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->vars, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->nNonz, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->lhss, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->rhss, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->consts, nrows) );

   for( i = 0; i < nrows; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->coefs[i], SCIProwGetNNonz(rows[i])) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(*projlp)->vars[i], SCIProwGetNNonz(rows[i])) );
      (*projlp)->nNonz[i] = 0;
      (*projlp)->lhss[i] = SCIProwGetLhs(rows[i]);
      (*projlp)->rhss[i] = SCIProwGetRhs(rows[i]);
      (*projlp)->consts[i] = SCIProwGetConstant(rows[i]);

      cols = SCIProwGetCols(rows[i]);
      for( v = 0; v < SCIProwGetNNonz(rows[i]); ++v )
      {
         var = SCIPcolGetVar(cols[v]);
         val = SCIPgetSolVal(scip, sol, var);
         vlb = local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
         vub = local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
         if( vlb == val || vub == val )
         {
            /* add var as a constant to row of projlp */
            if( !SCIPisInfinity(scip, -(*projlp)->lhss[i]) )
               (*projlp)->lhss[i] -= SCIProwGetVals(rows[i])[v]*val;
            if( !SCIPisInfinity(scip, (*projlp)->rhss[i]) )
               (*projlp)->rhss[i] -= SCIProwGetVals(rows[i])[v]*val;
         }
         else
         {
            /* add the entry to row of projlp */
            (*projlp)->coefs[i][(*projlp)->nNonz[i]] = SCIProwGetVals(rows[i])[v];
            (*projlp)->vars[i][(*projlp)->nNonz[i]] = var;
            ++(*projlp)->nNonz[i];
         }
      }
   }

   return SCIP_OKAY;
}

/* prints the projected LP */
static
void printProjLP(
   SCIP*   scip,        /**< SCIP data structure */
   PROJLP* projlp,      /**< the projected LP */
   int     nrows,       /**< number of rows in projlp */
   FILE*   file         /**< output file (or NULL for standard output) */
   )
{
   int i,j;

   assert(projlp != NULL);

   for( i = 0; i < nrows; ++i )
   {
      SCIPinfoMessage(scip, file, "\nproj_row[%d]: ", i);
      if( projlp->lhss[i] != -SCIPinfinity(scip) )
         SCIPinfoMessage(scip, file, "%.15g <= ", projlp->lhss[i]);
      for( j = 0; j < projlp->nNonz[i]; ++j )
      {
         if( j == 0 )
         {
            if( projlp->coefs[i][j] < 0 )
               SCIPinfoMessage(scip, file, "-");
         }
         else
         {
            if( projlp->coefs[i][j] < 0 )
               SCIPinfoMessage(scip, file, " - ");
            else
               SCIPinfoMessage(scip, file, " + ");
         }

         if( projlp->coefs[i][j] != 1.0 )
            SCIPinfoMessage(scip, file, "%.15g*", REALABS(projlp->coefs[i][j]));
         SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(projlp->vars[i][j]));
      }
      if( projlp->consts[i] > 0 )
         SCIPinfoMessage(scip, file, " + %.15g", projlp->consts[i]);
      else if( projlp->consts[i] < 0 )
         SCIPinfoMessage(scip, file, " - %.15g", REALABS(projlp->consts[i]));

      if( projlp->rhss[i] != SCIPinfinity(scip) )
         SCIPinfoMessage(scip, file, " <= %.15g", projlp->rhss[i]);
   }
   SCIPinfoMessage(scip, file, "\n");
}

/** frees the projected LP
 */
static
void freeProjLP(
   SCIP*    scip,   /**< SCIP data structure */
   PROJLP** projlp, /**< the projected LP */
   int      nrows   /**< number of rows in projlp */
   )
{
   int i;

   for( i = 0; i < nrows; ++i )
   {
      SCIPfreeBufferArray(scip, &(*projlp)->vars[i]);
      SCIPfreeBufferArray(scip, &(*projlp)->coefs[i]);
   }

   SCIPfreeBufferArray(scip, &(*projlp)->consts);
   SCIPfreeBufferArray(scip, &(*projlp)->rhss);
   SCIPfreeBufferArray(scip, &(*projlp)->lhss);
   SCIPfreeBufferArray(scip, &(*projlp)->nNonz);
   SCIPfreeBufferArray(scip, &(*projlp)->vars);
   SCIPfreeBufferArray(scip, &(*projlp)->coefs);
   SCIPfreeBuffer(scip, projlp);
}

/* mark a row for rlt cut selection
 *
 * depending on the sign of value and row inequality type, set the mark to:
 * 1 - cuts for axy < aw case,
 * 2 - cuts for axy > aw case,
 * 3 - cuts for both cases
 */
static
void addRowMark(
   int           ridx,         /**< row index */
   SCIP_Real     value,        /**< ai*(w - xy) */
   int*          row_idcs,     /**< sparse array with indices of marked rows */
   int*          row_marks,    /**< sparse array to store the marks */
   int*          nmarked       /**< number of marked rows */
)
{
   int newmark;
   int pos;
   SCIP_Bool exists;

   assert(value != 0.0);

   if( value > 0.0)
      newmark = 1;
   else newmark = 2;

   /* find row idx in row_idcs */
   exists = SCIPsortedvecFindInt(row_idcs, ridx, *nmarked, &pos);

   if( exists )
   {
      /* we found the row index: update the mark at pos1 */
      if( (newmark == 1 && row_marks[pos] == 2) || (newmark == 2 && row_marks[pos] == 1) )
      {
         row_marks[pos] = 3;
      }
   }
   else /* the given row index does not yet exist in row_idcs */
   {
      int i;

      /* insert variable at the correct position */
      for( i = *nmarked; i > pos; --i )
      {
         row_idcs[i] = row_idcs[i-1];
         row_marks[i] = row_marks[i-1];
      }
      row_idcs[pos] = ridx;
      row_marks[pos] = newmark;
      (*nmarked)++;
   }
}

/* mark all rows that should be multiplied by xj */
static
SCIP_RETCODE markRowsXj(
   SCIP*          scip,       /**< SCIP data structure */
   SCIP_SEPADATA* sepadata,   /**< separator data */
   SCIP_CONSHDLR* conshdlr,   /**< constraint handler */
   SCIP_SOL*      sol,        /**< point to be separated (can be NULL) */
   int            j,          /**< index of the multiplier variable in sepadata */
   SCIP_Bool      local,      /**< are local cuts allowed? */
   int*           row_marks,  /**< sparse array storing the row marks */
   int*           row_idcs,   /**< sparse array storing the marked row positions */
   int*           nmarked     /**< number of marked rows */
)
{
   int i, idx, img, ncolrows, r, ridx;
   SCIP_VAR* xi;
   SCIP_VAR* xj;
   SCIP_Real vlb, vub, val;
   SCIP_Real prod_viol, a;
   SCIP_COL* coli;
   SCIP_Real* colvals;
   SCIP_ROW** colrows;
   SCIP_Real auxval;

   *nmarked = 0;

   xj = sepadata->varssorted[j];
   assert(xj != NULL);

   val = SCIPgetSolVal(scip, sol, xj);
   vlb = local ? SCIPvarGetLbLocal(xj) : SCIPvarGetLbGlobal(xj);
   vub = local ? SCIPvarGetUbLocal(xj) : SCIPvarGetUbGlobal(xj);

   if( vlb == val || vub == val )
   {
      /* we don't want to multiply by variables that are at bound */
      SCIPdebugMsg(scip, "Rejected multiplier %s in [%f,%f] because it is at bound (current value %f)\n", SCIPvarGetName(xj), vlb, vub, val);
      return SCIP_OKAY;
   }

   /* for each var which appears in a bilinear product together with xj, mark rows */
   for( i = 0; i < sepadata->nvarbilinvars[j]; ++i )
   {
      xi = sepadata->varbilinvars[j][i];

      val = SCIPgetSolVal(scip, sol, xi);
      vlb = local ? SCIPvarGetLbLocal(xi) : SCIPvarGetLbGlobal(xi);
      vub = local ? SCIPvarGetUbLocal(xi) : SCIPvarGetUbGlobal(xi);

      if( vlb == val || vub == val ) /* we aren't interested in products with variables that are at bound */
         break;

      /* find the bilinear product */
      if( SCIPvarComp(xj, xi) < 0 )
         idx = SCIPvarGetIndex(xj) * sepadata->maxvarindex + SCIPvarGetIndex(xi);
      else
         idx = SCIPvarGetIndex(xi) * sepadata->maxvarindex + SCIPvarGetIndex(xj);
      assert( SCIPhashmapExists(sepadata->bilinvarsmap, (void*)(size_t)idx) );
      img = SCIPhashmapGetImageInt(sepadata->bilinvarsmap, (void*)(size_t) idx);
      SCIPevalConsExprExpr(scip, conshdlr, sepadata->bilinterms[img], sol, 0);
      auxval = SCIPevalConsExprExpr(scip, sepadata->conshdlr, sepadata->linexprs[img][0], sol, 0);
      prod_viol = auxval - SCIPgetConsExprExprValue(sepadata->bilinterms[img]);
      SCIPdebugMsg(scip, "aux val = %f, prod val = %f, prod viol = %f\n", auxval,
         SCIPgetConsExprExprValue(sepadata->bilinterms[img]), prod_viol);

      /* we are interested only in violated product relations */
      if( (SCIPgetConsExprExprNLocksNeg(sepadata->bilinterms[img]) == 0 || SCIPisFeasLE(scip, prod_viol, 0.0)) &&
         (SCIPgetConsExprExprNLocksPos(sepadata->bilinterms[img]) == 0 && SCIPisFeasGE(scip, prod_viol, 0.0)) )
      {
         SCIPdebugMsg(scip, "the product for vars %s, %s is not violated\n", SCIPvarGetName(xj), SCIPvarGetName(xi));
         continue;
      }

      /* get the column of xi */
      coli = SCIPvarGetCol(xi);
      colvals = SCIPcolGetVals(coli);
      ncolrows = SCIPcolGetNNonz(coli);

      SCIPdebugMsg(scip, NULL, "marking rows for xj, xi = %s, %s\n", SCIPvarGetName(xj), SCIPvarGetName(xi));

      /* mark the rows */
      for( r = 0; r < ncolrows; ++r )
      {
         a = colvals[r];
         if( a == 0.0 )
            continue;

         SCIPdebugMsg(scip, "Marking row %d\n", ridx);
         addRowMark(ridx, a*prod_viol, row_idcs, row_marks, nmarked);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE separateRltCuts(
   SCIP*          scip,         /**< SCIP data structure */
   SCIP_SEPA*     sepa,         /**< separator */
   SCIP_SEPADATA* sepadata,     /**< separator data */
   SCIP_CONSHDLR* conshdlr,     /**< constraint handler */
   SCIP_SOL*      sol,          /**< the point to be separated (can be NULL) */
   SCIP_HASHMAP*  row_to_pos,   /**< hashmap linking row indices to positions in array */
   PROJLP*        projlp,       /**< the projected LP */
   SCIP_ROW**     rows,         /**< problem rows */
   int            nrows,        /**< number of problem rows */
   SCIP_Bool      allowlocal,   /**< are local cuts allowed? */
   int*           ncuts,        /**< buffer to store the number of generated cuts */
   SCIP_RESULT*   result        /**< buffer to store whether separation was successful */
)
{
   int j, r, k, nmarked;
   SCIP_VAR* xj;
   int* row_marks;
   int* row_idcs;
   SCIP_ROW* cut;
   SCIP_Bool uselb[4] = {TRUE, TRUE, FALSE, FALSE};
   SCIP_Bool uselhs[4] = {TRUE, FALSE, TRUE, FALSE};
   SCIP_Bool success, infeasible, accepted, buildeqcut, iseqrow;

   assert(projlp != NULL);

   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &row_marks, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row_idcs, nrows) );

   /* loop through all variables that appear in bilinear products */
   for( j = 0; j < sepadata->nbilinvars && (sepadata->maxusedvars < 0 || j < sepadata->maxusedvars); ++j )
   {
      xj = sepadata->varssorted[j];

      /* mark all rows for multiplier xj */
      SCIP_CALL( markRowsXj(scip, sepadata, conshdlr, sol, j, allowlocal, row_marks, row_idcs, &nmarked) );

      /* generate the projected cut and if it is violated, generate the actual cut */
      for( r = 0; r < nmarked; ++r )
      {
         int pos;
         SCIP_ROW* row;

         assert(row_marks[r] != 0);

         if( !SCIPhashmapExists(row_to_pos, (void*)(size_t)row_idcs[r]) )
         {
            row_marks[r] = 0;
            continue; /* if row index is not in row_to_pos, it means that storeSuitableRows decided to ignore this row */
         }

         pos = SCIPhashmapGetImageInt(row_to_pos, (void*)(size_t)row_idcs[r]);
         row = rows[pos];
         assert(SCIProwGetIndex(row) == row_idcs[r]);

         /* check whether this row and var fulfill the conditions */
         /* for now this is disabled */
         /* TODO decide what to do with this */
#if 0
         SCIP_CALL( isAcceptableRow(scip, sepadata, row, xj, sepadata->varpriorities[j], &accepted) );

         if( !accepted )
         {
            SCIPdebugMsg(scip, "rejected row %s for variable %s\n", SCIProwGetName(row), SCIPvarGetName(xj));
            row_marks[r] = 0;
            continue;
         }
#endif

         SCIPdebugMsg(scip, "accepted row %s for variable %s\n", SCIProwGetName(rows[r]), SCIPvarGetName(xj));
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPprintRow(scip, rows[r], NULL) );
#endif
         iseqrow = SCIPisEQ(scip, SCIProwGetLhs(row), SCIProwGetRhs(row));

         /* if all terms are known and it is an equality row, compute equality cuts */
         buildeqcut = (sepadata->currentnunknown == 0 && iseqrow);

         /* go over all suitable combinations of sides and bounds and compute the respective cuts */
         for( k = 0; k < 4; ++k )
         {
            /* if equality cuts are possible, lhs and rhs cuts are equal so skip rhs */
            if( buildeqcut )
            {
               if( k % 2 == 1 )
                  continue;
            }
            else
            {
               if( row_marks[r] == 1 && uselb[k] == uselhs[k] )
                  continue;

               if( row_marks[r] == 2 && uselb[k] != uselhs[k] )
                  continue;
            }

            success = TRUE;

            SCIPdebugMsg(scip, "uselb = %d, uselhs = %d\n", uselb[k], uselhs[k]);

            /* if no variables are left in the projected row, the RLT cut will not be violated */
            if( projlp->nNonz[pos] == 0 )
               continue;

            /* compute the rlt cut for a projected row first */
            SCIP_CALL( computeProjRltCut(scip, sepa, sepadata, &cut, projlp, pos, sol, xj, &success, uselb[k], uselhs[k],
                                      allowlocal, buildeqcut) );

            /* if the projected cut is not violated, set success to FALSE */
            if( cut != NULL )
               SCIPdebugMsg(scip, "proj cut viol = %f\n", SCIPgetRowFeasibility(scip, cut));
            if( cut != NULL && !SCIPisFeasLT(scip, SCIPgetRowFeasibility(scip, cut), 0.0) )
            {
               SCIPdebugMsg(scip, "projected cut is not violated, feasibility = %f\n", SCIPgetRowFeasibility(scip, cut));
               success = FALSE;
            }

            /* release the projected cut */
            if( cut != NULL )
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );

            /* if the projected cut was generated successfully and is violated, generate the actual cut */
            if( success )
               SCIP_CALL( computeRltCuts(scip, sepa, sepadata, &cut, row, sol, xj, &success, uselb[k], uselhs[k],
                  allowlocal, buildeqcut) );

            /* if the cut was created successfully and is violated, it is added to SCIP */
            if( success )
            {
               if( SCIPisFeasLT(scip, SCIPgetRowFeasibility(scip, cut), 0.0) )
               {
                  /* add the row to SCIP; equality cuts are forced to be added to the LP */
                  SCIP_CALL(SCIPaddRow(scip, cut, FALSE, &infeasible));
                  ++*ncuts;

                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, "CUTOFF! At least one of the cuts revealed infeasibility!\n");
                     *result = SCIP_CUTOFF;
                  } else
                  {
                     SCIPdebugMsg(scip, "SEPARATED: added cut to scip\n");
                     *result = SCIP_SEPARATED;
                  }
               }
               else
                  SCIPdebugMsg(scip, "the cut was created successfully, but not accepted by scip\n");
            } else
               SCIPdebugMsg(scip, "the generation of the cut failed\n");

            /* release the cut */
            if( cut != NULL)
            {
               SCIP_CALL(SCIPreleaseRow(scip, &cut));
            }

            if( (sepadata->maxncuts >= 0 && *ncuts >= sepadata->maxncuts) || *result == SCIP_CUTOFF )
            {
               SCIPdebugMsg(scip, "exit separator because we found enough cuts or a cutoff -> skip\n");
               SCIPdebugMsg(scip, "maxncuts = %d, ncuts = %d\n", sepadata->maxncuts, *ncuts);
               SCIPdebugMsg(scip, "result = %d\n", *result);
               /* entries of row_marks must be set to 0 before the array is freed */
               for( int r1 = r; r1 < nmarked; ++r1 )
                  row_marks[r1] = 0;
               goto TERMINATE;
            }
         }
         /* clear row_mark since it will be used for the next multiplier */
         row_marks[r] = 0;
      }
   }

   SCIPdebugMsg(scip, "exit separator because cut calculation is finished\n");

   TERMINATE:
   SCIPfreeCleanBufferArray(scip, &row_marks);
   SCIPfreeBufferArray(scip, &row_idcs);
   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyRlt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of separator */
   SCIP_CALL( SCIPincludeSepaRlt(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeRlt)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* free separator data */
   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolRlt)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   if( sepadata->iscreated )
   {
      SCIP_CALL( freeSepaData(scip, sepadata) );
   }

   return SCIP_OKAY;
}

/** creates an array of rows suitable for RLT cut generation */
static
void storeSuitableRows(
   SCIP*          scip,        /**< SCIP data structure */
   SCIP_SEPA*     sepa,        /**< separator */
   SCIP_SEPADATA* sepadata,    /**< separator data */
   SCIP_ROW**     prob_rows,   /**< problem rows */
   SCIP_ROW**     rows,        /**< an array to be filled with suitable rows */
   int*           nrows,       /**< buffer to store the number of suitable rows */
   SCIP_HASHMAP*  row_to_pos,  /**< hashmap linking row indices to positions in rows */
   SCIP_Bool      allowlocal   /**< are local rows allowed? */
   )
{
   int new_nrows, r, j;
   SCIP_Bool iseqrow;
   SCIP_COL** cols;
   SCIP_Bool iscontrow;

   new_nrows = 0;

   for( r = 0; r < *nrows; ++r )
   {
      iseqrow = SCIPisEQ(scip, SCIProwGetLhs(prob_rows[r]), SCIProwGetRhs(prob_rows[r]));

      /* if equality rows are requested, only those can be used */
      if( sepadata->onlyeqrows && !iseqrow )
         continue;

      /* if global cuts are requested, only globally valid rows can be used */
      if( !allowlocal && SCIProwIsLocal(prob_rows[r]))
         continue;

      /* if continuous rows are requested, only those can be used */
      if( sepadata->onlycontrows )
      {
         cols = SCIProwGetCols(prob_rows[r]);
         iscontrow = TRUE;

         /* check row for integral variables */
         for( j = 0; j < SCIProwGetNNonz(prob_rows[r]); ++j )
         {
            if( SCIPcolIsIntegral(cols[j]) )
            {
               iscontrow = FALSE;
               break;
            }
         }

         if( !iscontrow )
            continue;
      }

      /* don't try to use rows that have been generated by the RLT separator
       *
       * @TODO check whether name for McCormick cuts changes
       */
      if( SCIProwGetOriginSepa(prob_rows[r]) == sepa || strcmp(SCIProwGetName(prob_rows[r]), "mccormick") == 0 )
         continue;

      /* if we are here, the row has passed all checks and should be added to rows */
      rows[new_nrows] = prob_rows[r];
      SCIPhashmapSetImageInt(row_to_pos, (void*)(size_t)SCIProwGetIndex(prob_rows[r]), new_nrows);
      ++new_nrows;
   }

   *nrows = new_nrows;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpRlt)
{  /*lint --e{715}*/
   SCIP_ROW** prob_rows;
   SCIP_ROW** rows;
   SCIP_SEPADATA* sepadata;
   int ncalls;
   int depth;
   int ncuts;
   int nrows;
   SCIP_HASHMAP* row_to_pos;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   sepadata = SCIPsepaGetData(sepa);

   *result = SCIP_DIDNOTRUN;

   if( sepadata->maxncuts == 0 )
   {
      SCIPdebugMsg(scip, "exit separator because maxncuts is set to 0\n");
      return SCIP_OKAY;
   }

   /* don't run in a sub-SCIP or in probing */
   if( SCIPgetSubscipDepth(scip) > 0 && !sepadata->useinsubscip )
   {
      SCIPdebugMsg(scip, "exit separator because in sub-SCIP\n");
      return SCIP_OKAY;
   }

   /* don't run in a sub-SCIP or in probing */
   if( SCIPinProbing(scip) )
   {
      SCIPdebugMsg(scip, "exit separator because in or probing\n");
      return SCIP_OKAY;
   }

   /* only call separator a given number of times at each node */
   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
        || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
   {
      SCIPdebugMsg(scip, "exit separator because round limit for this node is reached\n");
      return SCIP_OKAY;
   }

   /* if this is called for the first time, create the sepadata and start the initial separation round */
   if( !sepadata->iscreated )
   {
      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( createSepaData(scip, sepadata) );
   }

   /* no bilinear terms available -> skip */
   if( sepadata->nbilinvars == 0 )
   {
      SCIPdebugMsg(scip, "exit separator because there are no known bilinear terms\n");
      return SCIP_OKAY;
   }

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
   {
      SCIPdebugMsg(scip, "exit separator because we are too close to terminating\n");
      return SCIP_OKAY;
   }

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMsg(scip, "exit separator because there is no LP solution at hand\n");
      return SCIP_OKAY;
   }

   /* get the rows, depending on settings */
   if( sepadata->isinitialround || sepadata->onlyinitial )
   {
      SCIP_CALL( getInitialRows(scip, &prob_rows, &nrows) );
   }
   else
   {
      SCIP_CALL( SCIPgetLPRowsData(scip, &prob_rows, &nrows) );
   }

   /* save the suitable rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &rows, nrows) );
   SCIP_CALL( SCIPhashmapCreate(&row_to_pos, SCIPblkmem(scip), nrows) );

   storeSuitableRows(scip, sepa, sepadata, prob_rows, rows, &nrows, row_to_pos, allowlocal);

   if( sepadata->isinitialround || sepadata->onlyinitial )
   {
      SCIPfreeBufferArray(scip, &prob_rows);
      sepadata->isinitialround = FALSE;
   } /* no need to free memory in the other case since SCIPgetLPRowsData does not allocate it */

   if( nrows == 0 ) /* no suitable rows found, free memory and exit */
   {
      SCIPhashmapFree(&row_to_pos);
      SCIPfreeBufferArray(scip, &rows);
      return SCIP_OKAY;
   }
   else /* suitable rows have been found */
   {
      SCIPreallocBufferArray(scip, &rows, nrows);
   }

   /* create the projected problem */
   PROJLP* projlp;
   createProjLP(scip, rows, nrows, NULL, &projlp, allowlocal);

#ifdef SCIP_DEBUG
   printProjLP(scip, projlp, nrows, NULL);
#endif

   /* separate the cuts */
   separateRltCuts(scip, sepa, sepadata, sepadata->conshdlr, NULL, row_to_pos, projlp, rows, nrows, allowlocal, &ncuts, result);

   /* free the projected problem */
   freeProjLP(scip, &projlp, nrows);

   SCIPhashmapFree(&row_to_pos);
   SCIPfreeBufferArray(scip, &rows);

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the RLT separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaRlt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create RLT separator data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpRlt, NULL, sepadata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyRlt) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeRlt) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolRlt) );

   /* add RLT separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/maxlinexprs",
      "maximal number of linearisation expressions per bilinear term (-1: unlimited)",
      &sepadata->maxlinexprs, FALSE, DEFAULT_MAXLINEXPRS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxncuts",
         "maximal number of rlt-cuts that are added per round (-1: unlimited)",
         &sepadata->maxncuts, FALSE, DEFAULT_MAXNCUTS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxunknownterms",
         "maximal number of unknown bilinear terms a row is still used with (-1: unlimited)",
         &sepadata->maxunknownterms, FALSE, DEFAULT_MAXUNKNOWNTERMS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxusedvars",
         "maximal number of variables used to compute rlt cuts (-1: unlimited)",
         &sepadata->maxusedvars, FALSE, DEFAULT_MAXUSEDVARS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/maxnonzeroprop",
         "maximal proportion of known bilinear terms of a variable to non-zeroes of a row that is accepted",
         &sepadata->maxnonzeroprop, FALSE, DEFAULT_MAXNONZEROPROP, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "separating/" SEPA_NAME "/maxrounds",
      "maximal number of eccuts separation rounds per node (-1: unlimited)",
      &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "separating/" SEPA_NAME "/maxroundsroot",
      "maximal number of eccuts separation rounds in the root node (-1: unlimited)",
      &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "separating/" SEPA_NAME "/onlyeqrows",
      "if set to true, only equality rows are used for rlt cuts",
      &sepadata->onlyeqrows, FALSE, DEFAULT_ONLYEQROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "separating/" SEPA_NAME "/onlycontrows",
      "if set to true, only continuous rows are used for rlt cuts",
      &sepadata->onlycontrows, FALSE, DEFAULT_ONLYCONTROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "separating/" SEPA_NAME "/onlyinitial",
      "if set to true, only initial constraints are used",
      &sepadata->onlyinitial, FALSE, DEFAULT_ONLYINITIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "separating/" SEPA_NAME "/useinsubscip",
      "if set to true, rlt is also used in sub-scips",
      &sepadata->useinsubscip, FALSE, DEFAULT_USEINSUBSCIP, NULL, NULL) );

   return SCIP_OKAY;
}
