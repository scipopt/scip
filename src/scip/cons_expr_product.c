/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_product.c
 * @brief  product expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 *
 * Implementation of the product expression, representing a product of expressions
 * and a constant, i.e., coef * prod_i x_i.
 *
 * @todo initsepaProduct
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_product.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_pow.h"

#include "scip/pub_misc.h"

#define PRODUCT_PRECEDENCE  50000
#define PRODUCT_HASHKEY     SCIPcalcFibHash(54949.0)

#define ADJUSTFACETTOL             1e-6 /**< adjust resulting facets in checkRikun() up to a violation of this value */
#define USEDUALSIMPLEX             TRUE /**< use dual or primal simplex algorithm? */

#define MAXMULTILINSEPALPSIZE      14   /**< maximum size of the multilinear separation LP */

#define DEFAULT_RANDSEED           101  /**< initial random seed */
#define MAXPERTURBATION            1e-3 /**< maximum perturbation */

/** first values for 2^n */
static const int poweroftwo[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 };

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

/** macro to activate/deactivate debugging information of simplify method */
#ifdef SIMPLIFY_DEBUG
#define debugSimplify                   printf
#else
#define debugSimplify                   while( FALSE ) printf
#endif

/*
 * Data structures
 */

struct SCIP_ConsExpr_ExprData
{
   SCIP_Real             coefficient;        /**< coefficient */

   SCIP_ROW*             row;                /**< row created during initLP() */
};

struct SCIP_ConsExpr_ExprHdlrData
{
   SCIP_LPI*             multilinearseparationlp; /**< lp to separate product expressions */
   int                   lpsize;             /**< number of rows - 1 of multilinearseparationlp */

   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
};

/** node for linked list of expressions */
struct exprnode
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< expression in node */
   struct exprnode*      next;               /**< next node */
};

typedef struct exprnode EXPRNODE;


/*
 * Local methods
 */

/*  methods for handling linked list of expressions */
/** inserts newnode at beginning of list */
static
void insertFirstList(
   EXPRNODE*             newnode,            /**< node to insert */
   EXPRNODE**            list                /**< list */
   )
{
   assert(list != NULL);
   assert(newnode != NULL);

   newnode->next = *list;
   *list = newnode;
}

/** removes first element of list and returns it */
static
EXPRNODE* listPopFirst(
   EXPRNODE**            list                /**< list */
   )
{
   EXPRNODE* first;

   assert(list != NULL);

   if( *list == NULL )
      return NULL;

   first = *list;
   *list = (*list)->next;
   first->next = NULL;

   return first;
}

/** returns length of list */
static
int listLength(
   EXPRNODE*             list                /**< list */
   )
{
   int length;

   if( list == NULL )
      return 0;

   length = 1;
   while( (list=list->next) != NULL )
      ++length;

   return length;
}

/** creates expression node and capture expression */
static
SCIP_RETCODE createExprNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression stored at node */
   EXPRNODE**            newnode             /**< pointer to store node */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, newnode) );

   (*newnode)->expr = expr;
   (*newnode)->next = NULL;
   SCIPcaptureConsExprExpr(expr);

   return SCIP_OKAY;
}

/** creates expression list from expressions */
static
SCIP_RETCODE createExprlistFromExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  exprs,              /**< expressions stored in list */
   int                   nexprs,             /**< number of expressions */
   EXPRNODE**            list                /**< pointer to store list */
   )
{
   int i;

   assert(*list == NULL);
   assert(nexprs > 0);

   debugSimplify("building expr list from %d expressions\n", nexprs); /*lint !e506 !e681*/
   for( i = nexprs - 1; i >= 0; --i )
   {
      EXPRNODE* newnode;

      SCIP_CALL( createExprNode(scip, exprs[i], &newnode) );
      insertFirstList(newnode, list);
   }
   assert(nexprs > 1 || (*list)->next == NULL);

   return SCIP_OKAY;
}

/** frees expression node and release expressions */
static
SCIP_RETCODE freeExprNode(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE**            node                /**< node to be freed */
   )
{
   assert(node != NULL && *node != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*node)->expr) );
   SCIPfreeBuffer(scip, node);

   return SCIP_OKAY;
}

/** frees an expression list */
static
SCIP_RETCODE freeExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE**            exprlist            /**< list */
   )
{
   EXPRNODE* current;

   if( *exprlist == NULL )
      return SCIP_OKAY;

   current = *exprlist;
   while( current != NULL )
   {
      EXPRNODE* tofree;

      tofree = current;
      current = current->next;
      SCIP_CALL( freeExprNode(scip, &tofree) );
   }
   assert(current == NULL);
   *exprlist = NULL;

   return SCIP_OKAY;
}

/** builds LP used to compute facets of the convex envelope of multilinear functions, see @ref separatePointProduct() */
static
SCIP_RETCODE buildMultilinearSeparationLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   lpsize,             /**< size of the LP */
   SCIP_LPI**            lp                  /**< pointer to store created LP */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   int* beg;
   int* ind;
   int nnonz;
   int ncols;
   int nrows;
   int i;
   int k;

   assert(scip != NULL);
   assert(lp != NULL);
   assert(0 < lpsize && lpsize < MAXMULTILINSEPALPSIZE);

   SCIPdebugMsg(scip, "Building LP for computing facets of convex envelope of multilinear terms\n");

   /* create lpi to store the LP */
   SCIP_CALL( SCIPlpiCreate(lp, SCIPgetMessagehdlr(scip), "edge concave LP", SCIP_OBJSEN_MINIMIZE) );

   nrows = lpsize + 1;
   ncols = poweroftwo[nrows - 1];
   nnonz = (ncols * (nrows + 1)) / 2;
   k = 0;

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );

   /* calculate nonzero entries in the LP; set obj, lb, and ub to zero */
   for( i = 0; i < ncols; ++i )
   {
      int row;
      int a;

      obj[i] = 0.0;
      lb[i] = 0.0;
      ub[i] = 0.0;

      SCIPdebugMsg(scip, "col %i starts at position %d\n", i, k);
      beg[i] = k;
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

         a <<= 1; /*lint !e701*/
         ++row;
         assert(0 <= row && row < MAXMULTILINSEPALPSIZE);
         assert(poweroftwo[row] == a);
      }

      /* put 1 as a coefficient for sum_{i} \lambda_i = 1 row (last row) */
      val[k] = 1.0;
      ind[k] = nrows - 1;
      ++k;
      SCIPdebugMsg(scip, " val[%d][%d] = 1 (position  %d)\n", nrows - 1, i, k);
   }
   assert(k == nnonz);

   /* add all columns to the LP interface; CPLEX needs the row to exist before adding columns, so we create the rows with
    * dummy sides; note that the assert is not needed once somebody fixes the LPI
    * FIXME: was the LPI fixed???
    */
   assert(nrows <= ncols);
   SCIP_CALL( SCIPlpiAddRows(*lp, nrows, obj, obj, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(*lp, ncols, obj, lb, ub, NULL, nnonz, beg, ind, val) );

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** evaluates multilinear function at corner of the domain \f$ [lb,ub] \f$ encoded by binary expansion of \f$ k \f$ */
static
SCIP_Real evalCorner(
   SCIP_VAR**            vars,               /**< vars in multilinear function */
   int                   nvars,              /**< number of variables in multilinear function */
   SCIP_Real             constant,           /**< constant of multilinear function */
   int                   k                   /**< k-th corner */
   )
{
   SCIP_Real val;
   int i;

   assert(k >= 0 && k < poweroftwo[nvars]);

   val = constant;
   for( i = 0; i < nvars; ++i )
      val *= ((poweroftwo[i]) & k) == 0 ? SCIPvarGetLbLocal(vars[i]) : SCIPvarGetUbLocal(vars[i]);

   return val;
}

/** the given facet might not be a valid under(over)estimator, because of numerics and bad fixings; we compute \f$
 * \max_{v \in V} f(v) - (\alpha v + \beta) \f$ (\f$\max_{v \in V} \alpha v + \beta - f(v) \f$) where \f$ V \f$ are the
 * vertices of the domain, see separatePointProduct
 */
static
SCIP_Real computeMaxFacetError(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            funvals,            /**< array containing the evaluation of the function at all corners */
   SCIP_VAR**            vars,               /**< variables representing \f$ x_i \f$ */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             overestimate,       /**< whether we check for an over or underestimator */
   SCIP_Real             midval,             /**< coefficient representing \f$ mid(y) \f$ */
   SCIP_INTERVAL         fixedinterval,      /**< interval evaluation of the fixed variables, \f$ \Pi_j y_j \f$ */
   SCIP_Real*            facet               /**< current facet candidate (array of size nvars + 1 */
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
   assert(facet != NULL);

   ncorners = (unsigned int) poweroftwo[nvars];
   maxerror = 0.0;

   /* check the origin */
   facetval = facet[nvars];
   for( i = 0; i < (unsigned int) nvars; ++i )
      facetval += facet[i] * SCIPvarGetLbLocal(vars[i]);

   /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
   funval = funvals[0] / midval;
   if( overestimate )
   {
      funval *= funval > 0 ? fixedinterval.sup : fixedinterval.inf;
      error = funval - facetval;
   }
   else
   {
      funval *= funval > 0 ? fixedinterval.inf : fixedinterval.sup;
      error = facetval - funval;
   }

   /* update maximum error */
   maxerror = MAX(error, maxerror);

   prev = 0;
   for( i = 1; i < ncorners; ++i )
   {
      unsigned int gray;
      unsigned int diff;
      unsigned int pos;

      gray = i ^ (i >> 1);
      diff = gray ^ prev;

      /* compute position of unique 1 of diff */
      pos = 0;
      while( (diff >>= 1) != 0 )
         ++pos;

      if( gray > prev )
         facetval += facet[pos] * (SCIPvarGetUbLocal(vars[pos]) - SCIPvarGetLbLocal(vars[pos]));
      else
         facetval -= facet[pos] * (SCIPvarGetUbLocal(vars[pos]) - SCIPvarGetLbLocal(vars[pos]));


      /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
      funval = funvals[gray] / midval;
      if( overestimate )
      {
         funval *= funval > 0 ? fixedinterval.sup : fixedinterval.inf;
         error = funval - facetval;
      }
      else
      {
         funval *= funval > 0 ? fixedinterval.inf : fixedinterval.sup;
         error = facetval - funval;
      }

      /* update  maximum error */
      maxerror = MAX(error, maxerror);

      prev = gray;
   }

   SCIPdebugMsg(scip, "maximum error of facet: %2.8e\n", maxerror);

   return maxerror;
}

/** computes a facet of convex/concave envelope of \f$ w = c \Pi_i x_i \f$ that tries to separate \f$ (w^*, x^*) \f$;
 * checks that the facet is valid for \f$ w = \frac{c}{mid(y)} \Pi_i x_i \Pi_j y_j \f$, where the \f$ y_j \f$ are
 * (almost) fixed, see technical details of @ref separatePointProduct()
 */
static
SCIP_RETCODE computeFacet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for perturbation */
   SCIP_LPI*             lp,                 /**< lp used to compute facet */
   SCIP_SOL*             sol,                /**< solution representing \f$ (w^*, x^*) \f$ */
   SCIP_VAR**            vars,               /**< variables representing \f$ x_i \f$ */
   int                   nvars,              /**< number of variables */
   SCIP_VAR*             auxvar,             /**< variable representing \f$ w \f$ */
   SCIP_Real             coefficient,        /**< coefficient representing \f$ c \f$ */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real             midval,             /**< coefficient representing \f$ mid(y) \f$ */
   SCIP_INTERVAL         fixedinterval,      /**< interval evaluation of the fixed variables, \f$ \Pi_j y_j \f$ */
   SCIP_Real*            violation,          /**< buffer to store the violation between facet and sol */
   SCIP_Real*            facet               /**< buffer to store the facet; facet[nvars] holds facet's constant */
   )
{
   SCIP_Real* funvals;
   SCIP_Real* aux; /* used for settings sides and getting the dual solution */
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real maxfaceterror; /* stores violation between facet and function */
   int* inds;
   int ncorners;
   int ncols;
   int nrows;
   int i;

   assert(scip != NULL);
   assert(randnumgen != NULL);
   assert(lp != NULL);
   assert(sol != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(auxvar != NULL);
   assert(violation != NULL);
   assert(facet != NULL);

   /* get number of cols and rows of separation lp */
   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &funvals, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aux, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );

   /*
    * 1. compute f(v^i) for each corner v^i of [l,u]
    * 2. set up the described LP on the transformed space
    */

   /* get numbre of corners: 2^nvars */
   ncorners = poweroftwo[nvars];

   for( i = 0; i < ncols; ++i )
   {
      funvals[i] = i < ncorners ? evalCorner(vars, nvars, coefficient, i) : 0.0;
      inds[i] = i;

      /* update bounds; variables that are not in the LP are get fixed to 0 */
      lb[i] = 0.0;
      ub[i] = i < ncorners ? 1.0 : 0.0;

      SCIPdebugMsg(scip, "bounds of LP col %d = [%e, %e]; obj = %e\n", i, lb[i], ub[i], funvals[i]);
   }

   /* compute T^-1(x^*), i.e. T^-1(x^*)_i = (x^*_i - lb_i)/(ub_i - lb_i) */
   for( i = 0; i < nrows; ++i )
   {
      if( i < nvars )
      {
         SCIP_Real solval;

         assert(vars[i] != NULL);
         solval = SCIPgetSolVal(scip, sol, vars[i]);

         /* explicitely handle solution which violate bounds of variables (this can happen because of tolerances) */
         if( solval < SCIPvarGetLbLocal(vars[i]) )
            aux[i] = 0.0;
         else if( solval > SCIPvarGetUbLocal(vars[i]) )
            aux[i] = 1.0;
         else
            aux[i] = (SCIPgetSolVal(scip, sol, vars[i]) - SCIPvarGetLbLocal(vars[i])) /
               (SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]));

         /* perturb point to hopefuly obtain a facet of the convex envelope */
         if( aux[i] == 1.0 )
            aux[i] -= SCIPrandomGetReal(randnumgen, 0.0, MAXPERTURBATION);
         else if( aux[i] == 0.0 )
            aux[i] += SCIPrandomGetReal(randnumgen, 0.0, MAXPERTURBATION);
         else
         {
            SCIP_Real perturbation;

            perturbation = MIN( aux[i], 1.0 - aux[i] ) / 2.0;
            perturbation = MIN( perturbation, MAXPERTURBATION );
            aux[i] += SCIPrandomGetReal(randnumgen, -perturbation, perturbation);
         }
         assert(0.0 < aux[i] && aux[i] < 1.0);
      }
      else
      {
         /* constraints between nvars and nrows - 2 should be 0 == 0; last row corresponds to sum_{j} \lambda_j = 1 */
         aux[i] = (i == nrows - 1) ? 1.0 : 0.0;
      }

      SCIPdebugMsg(scip, "LP row %d in [%e, %e]\n", i, aux[i], aux[i]);
   }

   /* update LP */
   SCIP_CALL( SCIPlpiChgObj(lp, ncols, inds, funvals) );
   SCIP_CALL( SCIPlpiChgBounds(lp, ncols, inds, lb, ub) );
   SCIP_CALL( SCIPlpiChgSides(lp, nrows, inds, aux, aux) );
   SCIP_CALL( SCIPlpiChgObjsen(lp, overestimate ? SCIP_OBJSEN_MAXIMIZE : SCIP_OBJSEN_MINIMIZE) );
   SCIP_CALL( SCIPlpiWriteLP(lp, "lp.lp") );

   /* free memory used to update the LP */
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &inds);

   /*
    * 3. solve the LP and store the resulting facet for the transformed space
    */
   if( USEDUALSIMPLEX ) /*lint !e774 !e506*/
   {
      SCIP_CALL( SCIPlpiSolveDual(lp) );
   }
   else
   {
      SCIP_CALL( SCIPlpiSolvePrimal(lp) );
   }
   /* @todo: check solution status */

   /* get dual solution (facet of convex envelope); again, we have to be careful since the LP can have more rows and
    * columns than needed, in particular, \bar \beta is the last dual multiplier
    */
   SCIP_CALL( SCIPlpiGetSol(lp, NULL, NULL, aux, NULL, NULL) );

   for( i = 0; i < nvars; ++i )
   {
      facet[i] = aux[i];
   }
   /* last dual multiplier is the constant */
   facet[nvars] = aux[nrows - 1];


#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "facet for the transformed problem: ");
   for( i = 0; i < nvars; ++i )
   {
      SCIPdebugMsgPrint(scip, "%3.4e * %s + ", facet[i], SCIPvarGetName(vars[i]));
   }
   SCIPdebugMsgPrint(scip, "%3.4e\n", facet[nvars]);
#endif

   /*
    *  4. transform the facet to original space and compute violation at x^*, i.e., \alpha x + \beta - w when
    *  underestimating, w - \alpha x - \beta when overestimating
    */

   SCIPdebugMsg(scip, "facet in orig. space: ");

   *violation = 0.0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real varlb;
      SCIP_Real varub;

      varlb = SCIPvarGetLbLocal(vars[i]);
      varub = SCIPvarGetUbLocal(vars[i]);
      assert(!SCIPisEQ(scip, varlb, varub));

      /* \alpha_i := \bar \alpha_i / (ub_i - lb_i) */
      facet[i] = facet[i] / (varub - varlb);

      /* \beta = \bar \beta - \sum_i \alpha_i * lb_i */
      facet[nvars] -= facet[i] * varlb;

      /* evaluate */
      *violation += facet[i] * SCIPgetSolVal(scip, sol, vars[i]);

      SCIPdebugMsgPrint(scip, "%3.4e * %s + ", facet[i], SCIPvarGetName(vars[i]));
   }
   SCIPdebugMsgPrint(scip, "%3.4e ", facet[nvars]);

   /* add \beta to the violation: at this point in the code, violation = g(x^*) */
   *violation += facet[nvars];

   /* compute actual violation */
   if( overestimate )
      *violation = SCIPgetSolVal(scip, sol, auxvar) - *violation;
   else
      *violation = *violation - SCIPgetSolVal(scip, sol, auxvar);

   SCIPdebugMsgPrint(scip, "has a violation of %g\n", violation);

   /* if cut doesn't separate x^* (i.e. violation <= 0) there is no point in going on, since we only weaking the cut */
   if( SCIPisLE(scip, *violation, 0.0) )
   {
      /* free memory and return */
      SCIPfreeBufferArray(scip, &aux);
      SCIPfreeBufferArray(scip, &funvals);

      return SCIP_OKAY;
   }

   /*
    *  5. check and adjust facet with the algorithm of Rikun et al.
    */

   maxfaceterror = computeMaxFacetError(scip, funvals, vars, nvars, overestimate, midval, fixedinterval, facet);

   /* adjust constant part of the facet by maxerror to make it a valid over/underestimator (not facet though) */
   if( maxfaceterror > 0 )
   {
      /* there seem to be numerical problems if the error is too large; in this case we reject the facet */
      if( maxfaceterror > ADJUSTFACETTOL )
      {
         SCIPdebugMsg(scip, "ignoring facet due to instability, it cuts off a vertex by %g.\n", maxfaceterror);
         *violation = -1.0;
      }

      if( overestimate )
         facet[nvars] += maxfaceterror;
      else
         facet[nvars] -= maxfaceterror;

      /* update violation */
      *violation -= maxfaceterror;
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &aux);
   SCIPfreeBufferArray(scip, &funvals);

   return SCIP_OKAY;
}

/* helper functions for simplifying expressions */

/** creates a product expression with the elements of exprlist as its children */
static
SCIP_RETCODE createExprProductFromExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             exprlist,           /**< list containing the children of expr */
   SCIP_Real             coef,               /**< coef of expr */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the product expression */
   )
{
   int i;
   int nchildren;
   SCIP_CONSEXPR_EXPR** children;

   /* asserts SP8 */
   assert(coef == 1.0);
   nchildren = listLength(exprlist);

   SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      children[i] = exprlist->expr;
      exprlist = exprlist->next;
   }

   assert(exprlist == NULL);

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, SCIPfindConshdlr(scip, "expr"), expr, nchildren, children, coef) );

   SCIPfreeBufferArray(scip, &children);

   return SCIP_OKAY;
}

/** simplifies a factor of a product expression: base, so that it is a valid children of a simplified product expr
 * @note: in contrast to other simplify methods, this does *not* return a simplified expression.
 * Instead, the method is intended to be called only when simplifying a product expression,
 * Since in general, base is not a simplified child of a product expression, this method returns
 * a list of expressions L, such that (prod L) = baset *and* each expression in L
 * is a valid child of a simplified product expression.
 */
static
SCIP_RETCODE simplifyFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   factor,             /**< expression to be simplified */
   SCIP_Real*            simplifiedcoef,     /**< coefficient of parent product expression */
   EXPRNODE**            simplifiedfactor    /**< pointer to store the resulting expression node/list of nodes */
   )
{
   const char* factortype;

   assert(simplifiedfactor != NULL);
   assert(*simplifiedfactor == NULL);
   assert(factor != NULL);

   factortype = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(factor));

   /* enforces SP7 */
   if( strcmp(factortype, "val") == 0 )
   {
      *simplifiedcoef *= SCIPgetConsExprExprValueValue(factor);
      return SCIP_OKAY;
   }

   /* enforces SP2 */
   if( strcmp(factortype, "prod") == 0 )
   {
      /* assert SP8 */
      assert(SCIPgetConsExprExprProductCoef(factor) == 1.0);
      debugSimplify("[simplifyFactor] seeing a product: include its children\n"); /*lint !e506 !e681*/

      SCIP_CALL( createExprlistFromExprs(scip, SCIPgetConsExprExprChildren(factor),
               SCIPgetConsExprExprNChildren(factor), simplifiedfactor) );

      return SCIP_OKAY;
   }

   /* the given (simplified) expression `factor`, can be a child of a simplified product */
   assert(strcmp(factortype, "prod") != 0);
   assert(strcmp(factortype, "val") != 0);
   SCIP_CALL( createExprNode(scip, factor, simplifiedfactor) );

   return SCIP_OKAY;
}

/** merges tomerge into finalchildren
 * Both, tomerge and finalchildren contain expressions that could be the children of a simplified product
 * (except for SP8 and SP10 which are enforced later).
 * However, the concatenation of both lists will not in general yield a simplified product expression,
 * because both SP4 and SP5 could be violated. So the purpose of this method is to enforce SP4 and SP5.
 * In the process of enforcing SP4, it could happen that SP2. Since enforcing SP2
 * could generate further violations, we remove the affected children from finalchildren
 * and include them in unsimplifiedchildren for further processing.
 * @note: if tomerge has more than one element, then they are the children of a simplified product expression
 */
static
SCIP_RETCODE mergeProductExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             tomerge,            /**< list to merge */
   EXPRNODE**            finalchildren,      /**< pointer to store the result of merge between tomerge and *finalchildren */
   EXPRNODE**            unsimplifiedchildren/**< the list of children that should go to the product expression; they are
                                                  unsimplified when seen as children of a simplified product */
   )
{
   EXPRNODE* tomergenode;
   EXPRNODE* current;
   EXPRNODE* previous;

   if( tomerge == NULL )
      return SCIP_OKAY;

   if( *finalchildren == NULL )
   {
      *finalchildren = tomerge;
      return SCIP_OKAY;
   }

   tomergenode = tomerge;
   current = *finalchildren;
   previous = NULL;

   while( tomergenode != NULL && current != NULL )
   {
      int compareres;
      EXPRNODE* aux;
      SCIP_CONSEXPR_EXPR* base1;
      SCIP_CONSEXPR_EXPR* base2;
      SCIP_Real expo1;
      SCIP_Real expo2;

      /* assert invariants */
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "val") != 0);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "val") != 0);
      assert(previous == NULL || previous->next == current);

      /* in general the base of an expression is itself if type(expr) != pow, otherwise it is child of pow */
      /* TODO: better documentation
       *       clean code */
      if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "pow") == 0 )
      {
         base1 = SCIPgetConsExprExprChildren(current->expr)[0];
         expo1 = SCIPgetConsExprExprPowExponent(current->expr);
      }
      else
      {
         base1 = current->expr;
         expo1 = 1.0;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "pow") == 0 )
      {
         base2 = SCIPgetConsExprExprChildren(tomergenode->expr)[0];
         expo2 = SCIPgetConsExprExprPowExponent(tomergenode->expr);
      }
      else
      {
         base2 = tomergenode->expr;
         expo2 = 1.0;
      }

      /* if both bases are the same: have to build simplifiy(base^(expo1 + expo2)) */
      if( SCIPcompareConsExprExprs(base1, base2) == 0 )
      {
         SCIP_CONSEXPR_EXPR* power;
         SCIP_CONSEXPR_EXPR* simplifiedpower;

         SCIP_CALL( SCIPcreateConsExprExprPow(scip, SCIPfindConshdlr(scip, "expr"), &power, base1, expo1 + expo2) );
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, power, &simplifiedpower) ); /* FIXME: call simplifyPow */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &power) );

         /* replace tomergenode's expression with simplifiedpower */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &tomergenode->expr) );
         tomergenode->expr = simplifiedpower;
         /* move tomergenode to unsimplifiedchildren */
         aux = tomergenode;
         tomergenode = tomergenode->next;
         insertFirstList(aux, unsimplifiedchildren);

         /* destroy current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = listPopFirst(finalchildren);
            assert(aux == current);
            current = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            aux = current;
            current = current->next;
            previous->next = current;
         }
         SCIP_CALL( freeExprNode(scip, &aux) );

         /* continue */
         continue;
      }

      /* bases are not the same, then expressions cannot be the same */
      compareres = SCIPcompareConsExprExprs(current->expr, tomergenode->expr);
      if( compareres == -1 )
      {
         /* current < tomergenode => move current */
         previous = current;
         current = current->next;
      }
      else
      {
         assert(compareres == 1);

         /* insert: if current is the first node, then insert at beginning; otherwise, insert between previous and current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = tomergenode;
            tomergenode = tomergenode->next;
            insertFirstList(aux, finalchildren);
            previous = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            /* extract */
            aux = tomergenode;
            tomergenode = tomergenode->next;
            /* insert */
            previous->next = aux;
            aux->next = current;
            previous = aux;
         }
      }
   }

   /* if all nodes of tomerge were merged, we are done */
   if( tomergenode == NULL )
      return SCIP_OKAY;

   assert(current == NULL);

   /* if all nodes of finalchildren were cancelled by nodes of tomerge (ie, transfered to unsimplifiedchildren),
    * then the rest of tomerge is finalchildren */
   if( *finalchildren == NULL )
   {
      assert(previous == NULL);
      *finalchildren = tomergenode;
      return SCIP_OKAY;
   }

   /* there are still nodes of tomerge unmerged; these nodes are larger than finalchildren, so append at end */
   assert(previous != NULL && previous->next == NULL);
   previous->next = tomergenode;

   return SCIP_OKAY;
}

static
SCIP_RETCODE createData(
   SCIP*                    scip,            /**< SCIP data structure */
   SCIP_CONSEXPR_EXPRDATA** exprdata,        /**< pointer where to store expression data */
   SCIP_Real                coefficient      /**< coefficient of product */
   )
{
   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, exprdata) );

   (*exprdata)->coefficient  = coefficient;
   (*exprdata)->row          = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** simplifies a product expression
 *
 * Summary: we first build a list of expressions (called finalchildren) which will be the children of the simplified product
 * and then we process this list in order to enforce SP8 and SP10
 * Description: In order to build finalchildren, we first build list of unsimplified children (called unsimplifiedchildren)
 * with the children of the product. Each node of the list is manipulated (see simplifyFactor) in order to satisfy
 * SP2 and SP7 as follows
 * SP7: if the node's expression is a value, multiply the value to the products's coef
 * SP2: if the node's expression is a product, then build a list with the child's children
 * Then, we merge the built list (or the simplified node) into finalchildren. While merging, nodes from finalchildren
 * can go back to unsimplifiedchildren for further processing (see mergeProductExprlist for more details)
 * After building finalchildren, we create the simplified product out of it, taking care that SP8 and SP10 are satisfied
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyProduct)
{
   EXPRNODE* unsimplifiedchildren;
   EXPRNODE* finalchildren;
   SCIP_Real simplifiedcoef;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") == 0);

   /* set up list of current children (when looking at each of them individually, they are simplified, but as
    * children of a product expression they might be unsimplified) */
   unsimplifiedchildren = NULL;
   SCIP_CALL( createExprlistFromExprs(scip, SCIPgetConsExprExprChildren(expr), SCIPgetConsExprExprNChildren(expr),
            &unsimplifiedchildren) );

   /* while there are still children to process */
   finalchildren  = NULL;
   simplifiedcoef = SCIPgetConsExprExprProductCoef(expr);
   while( unsimplifiedchildren != NULL )
   {
      EXPRNODE* tomerge;
      EXPRNODE* first;

      /* if the simplified coefficient is 0, we can return value 0 */
      if( simplifiedcoef == 0.0 )
      {
         SCIP_CALL( freeExprlist(scip, &finalchildren) );
         SCIP_CALL( freeExprlist(scip, &unsimplifiedchildren) );
         assert(finalchildren == NULL);
         break;
      }

      first = listPopFirst(&unsimplifiedchildren);
      assert(first != NULL);

      /* enforces SP2 and SP7 */
      tomerge = NULL;
      SCIP_CALL( simplifyFactor(scip, first->expr, &simplifiedcoef, &tomerge) );

      /* enforces SP4 and SP5
       * note: merge frees (or uses) the nodes of the tomerge list */
      SCIP_CALL( mergeProductExprlist(scip, tomerge, &finalchildren, &unsimplifiedchildren) );

      /* free first */
      SCIP_CALL( freeExprlist(scip, &first) );
   }

   /* build product expression from finalchildren and post-simplify */
   debugSimplify("[simplifyProduct] finalchildren has length %d\n", listLength(finalchildren));

   /* enforces SP10: if list is empty, return value */
   if( finalchildren == NULL )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr, simplifiedcoef) );
   }
   /* enforces SP10: if finalchildren has only one expr with coef 1.0, return that expr */
   else if( finalchildren->next == NULL && simplifiedcoef == 1.0 )
   {
      *simplifiedexpr = finalchildren->expr;
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }
   /* enforces SP10: if finalchildren has only one expr and coef != 1.0, return (sum 0 coef expr) */
   else if( finalchildren->next == NULL )
   {
      SCIP_CONSEXPR_EXPR* aux;

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), &aux,
               1, &(finalchildren->expr), &simplifiedcoef, 0.0) );

      /* simplifying here is necessary, the product could have sums as children e.g., (prod 2 (sum 1 <x>))
       * -> (sum 0 2 (sum 1 <x>)) and that needs to be simplified to (sum 0 2 <x>)
       */
      SCIP_CALL( SCIPsimplifyConsExprExpr(scip, aux, simplifiedexpr) ); /*FIXME: how to call simplifySum ? */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
   }
   /* enforces SP8: if simplifiedcoef != 1.0, transform it into a sum with the (simplified) product as child */
   else if( simplifiedcoef != 1.0 )
   {
      SCIP_CONSEXPR_EXPR* aux;

      SCIP_CALL( createExprProductFromExprlist(scip, finalchildren, 1.0, &aux) );
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr,
               1, &aux, &simplifiedcoef, 0.0) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
   }
   /* build product expression from list */
   else
   {
      SCIP_CALL( createExprProductFromExprlist(scip, finalchildren, simplifiedcoef, simplifiedexpr) );
   }

   /* free memory */
   SCIP_CALL( freeExprlist(scip, &finalchildren) );
   assert(finalchildren == NULL);

   assert(*simplifiedexpr != NULL);
   return SCIP_OKAY;
}

/** the order of two product expressions, u and v, is a lexicographical order on the factors.
 *  Starting from the *last*, we find the first child where they differ, say, the i-th.
 *  Then u < v <=> u_i < v_i.
 *  If there is no such children and they have different number of children, then u < v <=> nchildren(u) < nchildren(v)
 *  If all children are the same and they have the same number of childre, u < v <=> coeff(u) < coeff(v)
 *  Otherwise, they are the same.
 *  Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc
 *  Example: y * z < x * y * z
 */
static
SCIP_DECL_CONSEXPR_EXPRCMP(compareProduct)
{  /*lint --e{715}*/
   int compareresult;
   int i;
   int j;
   int nchildren1;
   int nchildren2;
   SCIP_CONSEXPR_EXPR** children1;
   SCIP_CONSEXPR_EXPR** children2;

   nchildren1 = SCIPgetConsExprExprNChildren(expr1);
   nchildren2 = SCIPgetConsExprExprNChildren(expr2);
   children1 = SCIPgetConsExprExprChildren(expr1);
   children2 = SCIPgetConsExprExprChildren(expr2);

   for( i = nchildren1 - 1, j = nchildren2 - 1; i >= 0 && j >= 0; --i, --j )
   {
      compareresult = SCIPcompareConsExprExprs(children1[i], children2[j]);
      if( compareresult != 0 )
         return compareresult;
      /* expressions are equal, continue */
   }

   /* all children of one expression are children of the other expression, use number of children as a tie-breaker */
   if( i < j )
   {
      assert(i == -1);
      /* expr1 has less elements, hence expr1 < expr2 */
      return -1;
   }
   if( i > j )
   {
      assert(j == -1);
      /* expr1 has more elements, hence expr1 > expr2 */
      return 1;
   }

   /* everything is equal, use coefficient as tie-breaker */
   assert(i == -1 && j == -1);
   if( SCIPgetConsExprExprProductCoef(expr1) < SCIPgetConsExprExprProductCoef(expr2) )
      return -1;
   if( SCIPgetConsExprExprProductCoef(expr1) > SCIPgetConsExprExprProductCoef(expr2) )
      return 1;

   /* they are equal */
   return 0;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrProduct)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** expression handler free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEHDLR(freehdlrProduct)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(exprhdlr != NULL);
   assert(exprhdlrdata != NULL);
   assert(*exprhdlrdata != NULL);

   /* free lp to separate product expressions */
   if( (*exprhdlrdata)->lpsize > 0 )
   {
      SCIP_CALL( SCIPlpiFree(&((*exprhdlrdata)->multilinearseparationlp)) );
      (*exprhdlrdata)->lpsize = 0;
   }
   assert((*exprhdlrdata)->lpsize == 0);
   assert((*exprhdlrdata)->multilinearseparationlp == NULL);

   /* free random number generator */
   SCIPrandomFree(&(*exprhdlrdata)->randnumgen);

   SCIPfreeBlockMemory(scip, exprhdlrdata);
   assert(*exprhdlrdata == NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataProduct)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* sourceexprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = SCIPgetConsExprExprData(sourceexpr);
   assert(sourceexprdata != NULL);

   SCIP_CALL( createData(targetscip, targetexprdata, sourceexprdata->coefficient) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemory(scip, &exprdata);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print opening parenthesis, if necessary */
         if( PRODUCT_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr) )
         {
            SCIPinfoMessage(scip, file, "(");
         }

         /* print coefficient, if not one */
         if( exprdata->coefficient != 1.0 )
         {
            if( exprdata->coefficient < 0.0 && PRODUCT_PRECEDENCE > SCIPgetConsExprExprWalkParentPrecedence(expr) )
            {
               SCIPinfoMessage(scip, file, "(%g)", exprdata->coefficient);
            }
            else
            {
               SCIPinfoMessage(scip, file, "%g", exprdata->coefficient);
            }
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         int childidx = SCIPgetConsExprExprWalkCurrentChild(expr);

         /* print multiplication sign, if not first factor */
         if( exprdata->coefficient != 1.0 || childidx > 0 )
         {
            SCIPinfoMessage(scip, file, "*");
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      {
         break;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* print closing parenthesis, if necessary */
         if( PRODUCT_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr) )
         {
            SCIPinfoMessage(scip, file, ")");
         }
         break;
      }
   }

   return SCIP_OKAY;
}

/** product hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *hashkey = PRODUCT_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(exprdata->coefficient);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      unsigned int childhash;

      assert(SCIPhashmapExists(expr2key, (void*)SCIPgetConsExprExprChildren(expr)[c]));
      childhash = (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[c]);

      *hashkey ^= childhash;
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalProduct)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real childval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *val = exprdata->coefficient;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr) && (*val != 0.0); ++c )
   {
      childval = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]);
      assert(childval != SCIP_INVALID); /*lint !e777*/

      *val *= childval;
   }

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffProduct)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));

   child = SCIPgetConsExprExprChildren(expr)[idx];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);
   assert(SCIPgetConsExprExprValue(child) != SCIP_INVALID);

   if( !SCIPisZero(scip, SCIPgetConsExprExprValue(child)) )
      *val = SCIPgetConsExprExprValue(expr) / SCIPgetConsExprExprValue(child);
   else
   {
      int i;

      *val = SCIPgetConsExprExprData(expr)->coefficient;
      for( i = 0; i < SCIPgetConsExprExprNChildren(expr) && (*val != 0.0); ++i )
      {
         if( i == idx )
            continue;

         *val *= SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[i]);
      }
   }

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPintervalSet(interval, exprdata->coefficient);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_INTERVAL childinterval;

      childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[c]);
      assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

      /* multiply childinterval with the so far computed interval */
      SCIPintervalMul(SCIPinfinity(scip), interval, *interval, childinterval);
   }

   return SCIP_OKAY;
}

/** separates a multilinear constraint of the form \f$ f(x) := a \Pi_{i = 1}^n x_i = w \f$ where \f$ x_i \f$ are the
 * auxiliary variables of the children and \f$ w \f$ is the auxiliary variable of expr. If \f$ f(x^*) > w^* \f$, then we
 * look for an affine underestimator of \f$ f(x) \f$ which separates \f$ (x^*, w^*) \f$ from the feasible region, i.e.
 * \f$ g(x) := \alpha^T x + \beta \le f(x) = w \f$ for all \f$ x \f$ in the domain, such that \f$ \alpha x^* > w^* \f$.
 *
 * Since \f$ f(x) \f$ is componentwise linear, its convex envelope is piecewise linear and its value can be computed by
 * finding the largest affine underestimator. Furthermore, it implies that \f$ g \f$ is an underestimator if and only if
 * \f$ g(v^i) \leq f(v^i), \forall i \f$, where \f$ \{ v^i \}_{i = 1}^{2^n} \subseteq \mathbb R^n \f$ are the vertices
 * of the domain of \f$ x \f$, \f$ [\ell,u] \f$. Hence, we can compute a linear underestimator by solving the following
 * LP (we don't necessarily get a facet of the convex envelope, see Technical detail below):
 *
 * \f[
 *              \max \, \alpha^T x^* + \beta
 * \f]
 * \f[
 *     s.t. \; \alpha^T v^i + \beta \le f(v^i), \, \forall i = 1, \ldots, 2^n
 * \f]
 *
 * In principle, one would need to update the LP whenever the domain changes. However, \f$ [\ell,u] = T([0, 1]^n) \f$,
 * where \f$ T \f$ is an affine linear invertible transformation given by \f$ T(y)_i = (u_i - \ell_i) y_i + \ell_i \f$.
 * Working with the change of variables \f$ x = T(y) \f$ allows us to keep the constraints of the LP, even it the domain
 * changes. Indeed, after the change of variables, the problem is: find an affine underestimator \f$ g \f$ such that \f$
 * g(T(y)) \le f(T(y)) \f$, for all \f$ y \in [0, 1]^n \f$. Now \f$ f(T(y)) \f$ is componentwise affine, but still
 * satisfies that \f$ g \f$ is a valid underestimator if and only if \f$ g(T(u)) \leq f(T(u)), \forall u \in \{0, 1\}^n
 * \f$. So we now look for \f$ \bar g(y) := g(T(y)) = g(((u_i - \ell_i) y_i + \ell_i)_i) = \bar \alpha^T y + \bar \beta
 * \f$, where \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i \f$ and \f$ \bar \beta = \sum_i \alpha_i \ell_i + \beta \f$. So
 * we find \f$ \bar g \f$ by solving the LP:
 *
 * \f[
 *              \max \, \bar \alpha^T T^{-1}(x^*) + \bar \beta
 * \f]
 * \f[
 *     s.t. \; \bar \alpha^T u + \bar \beta \le f(T(u)), \, \forall u \in \{0, 1\}^n
 * \f]
 *
 * and recover \f$ g \f$ by solving \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i, \bar \beta = \sum_i \alpha_i \ell_i +
 * \beta \f$. Notice that \f$ f(T(u^i)) = f(v^i) \f$ so the right hand side doesn't change after the change of variables
 *
 * Furthermore, the LP has more constraints than variables, so we solve its dual:
 *
 * \f[
 *              \min \, \sum_i \lambda_i f(v^i)
 * \f]
 * \f[
 *     s.t. \; \sum_i \lambda_i u^i = T^{-1}(x^*)
 * \f]
 * \f[
 *             \sum_i \lambda_i = 1
 * \f]
 * \f[
 *             \forall i, \, \lambda_i \geq 0
 * \f]
 *
 * In case we violate the constraint in the other direction, i.e., \f$ f(x^*) < w^* \f$, we do exactly the same, but
 * looking for an overestimator instead of an underestimator. This means, for the dual, that we have to maximize instead
 * of minimize.
 *
 * #### Technical and implementation details
 * -# The special case \f$ n = 2 \f$ is handled separately
 * -# \f$ U \f$ has exponentially many variables, so we only apply this separator for \f$ n \leq 10 \f$
 * -# We store a unique LP containing \f$ U = [u^1 | u^2 | \cdots | u^{2^n}] \f$, and \f$ U \f$ is build in such a way
 * that its submatrices consisting of the first \f$ k \f$ rows and first \f$ 2^k \f$ columns contains all the vectors in
 * \f$ \{0, 1\}^k \f$. This way, the same matrix can be used to separate a multilinear constraint with only \f$ k \f$
 * variables just by fixing \f$ \lambda_i = 0, i > 2^k \f$ to 0. The \f$ n + 1 \f$-th row is the row representing the
 * constraint \f$ \sum_i \lambda_i = 1 \f$, where \f$ n \f$ is the minimum between 10 and the maximum number of products
 * among all product expressions.
 * -# If the bounds are not finite, there is no underestimator. Also, \f$ T^{-1}(x^*) \f$ must be in the domain,
 * otherwise the dual is infeasible
 * -# After a facet is computed, we check whether it is a valid facet (i.e. we check \f$ \alpha^T v + \beta \le f(v) \f$
 * for every vertex \f$ v \f$). If we find a violation of at most ADJUSTFACETTOL, then we weaken \f$ \beta \f$ by this
 * amount, otherwise, we discard the cut.
 * -# If a variable is fixed within tolerances, we replace it with its value and compute the facet of the remaining
 * expression. Note that since we are checking the cut for validity, this will never produce wrong result.
 * -# In every iteration we set _all_ \f$ 2^n \f$ bounds and objective values. The reason is that different products
 * have different number of children, so we might need to fix/unfix more variables
 * -# If \f$ x^* \f$ is in the boundary of the domain, then the LP has infinitely many solutions, some of which might
 * have very bad numerical properties. For this reason, we perturb \f$ x^* \f$ to be in the interior of the region.
 * Furthermore, for some interior points, there might also be infinitely many solutions (e.g. for \f$ x y \f$ in \f$
 * [0,1]^2 \f$ any point \f$ (x^*, y^*) \f$ such that \f$ y^* = 1 - x^* \f$ has infinitely many solutions). For this
 * reason, we perturb any given \f$ x^* \f$. The idea is to try to get a facet of the convex/concave envelope. This only
 * happens when the solution has \f$ n + 1 \f$ non zero \f$ \lambda \f$'s (i.e. the primal has a unique solution)
 * -# We need to compute \f$ f(v^i) \f$ for every vertex of \f$ [\ell,u] \f$. A vertex is encoded by a number between 0
 * and \f$ 2^n - 1 \f$, via its binary representation (0 bit is lower bound, 1 bit is upper bound), so we can compute
 * all these values by iterating between 0 and \f$ 2^n - 1 \f$.
 * @note This could be accelerated: Since \f$ f(v^i) = 0 \f$, whenever a component of the vertex is 0, we can sort the
 * variables so that the first \f$ p \f$ have a 0.0 bound. This way, if \f$ i_0 \f$ is the first vertex such that \f$
 * f(v^{i_0}) = 0 \f$, then we can loop increasing by \f$ 2^p \f$ in every iteration.
 * -# to check that the computed cut is valid we do the following: when there are no fixed variables, we use a gray code
 * to loop over the vertices of the box domain in order to evaluate the underestimator. When there are fixed variables,
 * we compute the underestimator for a different function: suppose \f$ f(x,y) = c \Pi_i x_i \Pi_j y_j \f$ is the
 * function we are underestimating and we fix the \f$ y \f$ variables at their mid point. Let \f$ y_m = \Pi_j mid(y_j)
 * \f$, then we are actually computing an underestimator for \f$ h(x) = f(x, mid(y)) = c y_m \Pi_i x_i \f$. Let \f$
 * \alpha x + \beta \f$ be the underestimator and \f$ [y_l, y_u] \f$ be the interval of \f$ \Pi_j y_j \f$. To ensure the
 * validity of the underestimator, we check whether \f$ \alpha x + \beta \le f(x,y) \f$ for every vertex of the domains
 * \f$ [x] \f$ and \f$ [y] \f$.  Given a vertex of \f$ x^i \in [x] \f$ we just need to pick the vertex of \f$ [y] \f$
 * that produces the worst value for \f$ f(x^i, y) \f$. However, since \f$ y \f$ does not appear in the underestimator,
 * we don't care about the vertex itself, but of the value. This value has to be \f$ y_l \f$ or \f$ y_u \f$. So we
 * actually just check whether \f$ \alpha x^i + \beta \le \min{ \frac{h(x^i) y_l}{y_m}, \frac{h(x^i) y_u}{y_m} } \f$
 *
 * @todo the solution is a facet if all variables of the primal have positive reduced costs (i.e. the solution is
 * unique). In the dual, this means that there are \f$ n + 1 \f$ variables with positive value. Can we use this or some
 * other information to handle any of both cases (point in the boundary or point in the intersection of polytopes
 * defining different pieces of the convex envelope)? In the case where the point is in the boundary, can we use that
 * information to maybe solve another to find a facet? How do the polytopes defining the pieces where the convex
 * envelope is linear looks like, i.e, given a point in the interior of a facet of the domain, does the midpoint of the
 * segment joining \f$ x^* \f$ with the center of the domain, always belongs to the interior of one of those polytopes?
 */
static
SCIP_RETCODE separatePointProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< product expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_ROW**            cut                 /**< pointer to store the row */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* var;
   SCIP_Real violation;
   SCIP_Bool overestimate;
   SCIP_Bool success;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") == 0);
   assert(cut != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);

   *cut = NULL;

   /* compute violation of the expression by evaluating auxiliary variables */
   violation = exprdata->coefficient;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];

      /* value expressions should have been removed during simplification */
      assert(SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrValue(conshdlr));

      var = SCIPgetConsExprExprLinearizationVar(child);
      assert(var != NULL);

      violation *= SCIPgetSolVal(scip, sol, var);
   }
   violation -= SCIPgetSolVal(scip, sol, auxvar);

   /* no violation in this sub-expression */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   overestimate = SCIPisLT(scip, violation, 0.0);
   success = FALSE;

   /* debug output: prints expression we are trying to separate, bounds of variables and point */
#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "separating product with %d variables: will try to separate violated point (%g) by an %s\n",
         SCIPgetConsExprExprNChildren(expr), violation, overestimate ? "overestimator": "underestimator");
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];
      var = SCIPgetConsExprExprLinearizationVar(child);
      assert(var != NULL);
      SCIPdebugMsg(scip, "var: %s = %g in [%g, %g]\n", SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

      if( SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) || SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
      {
         SCIPdebugMsg(scip, "unbounded factor related to\n");
         SCIP_CALL( SCIPdismantleConsExprExpr(scip, child) );
      }
   }
   SCIPdebugMsg(scip, "The product should be equal to auxvar: %s = %g in [%g, %g]\n", SCIPvarGetName(auxvar),
         SCIPgetSolVal(scip, sol, auxvar), SCIPvarGetLbLocal(auxvar), SCIPvarGetUbLocal(auxvar));
#endif

   /* bilinear term */
   if( SCIPgetConsExprExprNChildren(expr) == 2 )
   {
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_Real lincoefx;
      SCIP_Real lincoefy;
      SCIP_Real linconstant;
      SCIP_Real refpointx;
      SCIP_Real refpointy;

      /* collect first variable */
      child = SCIPgetConsExprExprChildren(expr)[0];
      x = SCIPgetConsExprExprLinearizationVar(child);
      assert(x != NULL);

      /* collect second variable */
      child = SCIPgetConsExprExprChildren(expr)[1];
      y = SCIPgetConsExprExprLinearizationVar(child);
      assert(y != NULL);

      lincoefx = 0.0;
      lincoefy = 0.0;
      linconstant = 0.0;
      refpointx = SCIPgetSolVal(scip, sol, x);
      refpointy = SCIPgetSolVal(scip, sol, y);
      success = TRUE;

      /* adjust the reference points */
      refpointx = SCIPisLT(scip, refpointx, SCIPvarGetLbLocal(x)) ? SCIPvarGetLbLocal(x) : refpointx;
      refpointx = SCIPisGT(scip, refpointx, SCIPvarGetUbLocal(x)) ? SCIPvarGetUbLocal(x) : refpointx;
      refpointy = SCIPisLT(scip, refpointy, SCIPvarGetLbLocal(y)) ? SCIPvarGetLbLocal(y) : refpointy;
      refpointy = SCIPisGT(scip, refpointy, SCIPvarGetUbLocal(y)) ? SCIPvarGetUbLocal(y) : refpointy;
      assert(SCIPisLE(scip, refpointx, SCIPvarGetUbLocal(x)) && SCIPisGE(scip, refpointx, SCIPvarGetLbLocal(x)));
      assert(SCIPisLE(scip, refpointy, SCIPvarGetUbLocal(y)) && SCIPisGE(scip, refpointy, SCIPvarGetLbLocal(y)));

      SCIPaddBilinMcCormick(scip, exprdata->coefficient, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), refpointx,
            SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y), refpointy, overestimate, &lincoefx, &lincoefy, &linconstant,
            &success);

      /* @todo allow lhs/rhs of +/- infinity? */
      if( success && !SCIPisInfinity(scip, REALABS(linconstant)) )
      {
         /* McCormicks are only valid locally */
         SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "mccormick", 0, NULL, NULL, -SCIPinfinity(scip),
                  SCIPinfinity(scip), TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *cut, x, lincoefx) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cut, y, lincoefy) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );

         if( overestimate )
         {
            SCIP_CALL( SCIPchgRowLhs(scip, *cut, -linconstant) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowRhs(scip, *cut, -linconstant) );
         }
      }

      return SCIP_OKAY;
   }
   else
   {
      /* general case */
      SCIP_CONSEXPR_EXPRHDLRDATA* exprhdlrdata;
      SCIP_INTERVAL fixedinterval;
      SCIP_Real* facet;
      SCIP_VAR** vars;
      SCIP_Real midval;
      int nvars;
      int pos;
      int i;

      /* prepare data to compute a facet of the envelope */
      nvars = SCIPgetConsExprExprNChildren(expr);

      /* allocate necessary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &facet, nvars + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

      /* store variables: we only store variables that are not fixed; compute the product of the fixed mid point of the
       * fix variables along with their interval product
       */
      pos = 0;
      midval = 1.0;
      SCIPintervalSet(&fixedinterval, 1.0);
      for( i = 0; i < nvars; ++i )
      {
         var = SCIPgetConsExprExprLinearizationVar(SCIPgetConsExprExprChildren(expr)[i]);

         if( SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) || SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
         {
            SCIPdebugMsg(scip, "a factor is unbounded, no cut is possible\n");
            goto CLEANUP;
         }

         if( !SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)) )
         {
            vars[pos] = var;
            ++pos;
         }
         else
         {
            SCIP_INTERVAL bounds;

            midval *= (SCIPvarGetUbLocal(var) + SCIPvarGetLbLocal(var)) / 2.0;
            SCIPintervalSetBounds(&bounds, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
            SCIPintervalMul(SCIPinfinity(scip), &fixedinterval, fixedinterval, bounds);
         }
      }
      nvars = pos; /* update number of vars */
      assert(nvars >= 0);

      /* nothing to do if all variables are fixed */
      /* TODO if a single variable is unfix, should we add the linear equality w = constant * the_unfix_var ?? */
      if( nvars == 0 )
         goto CLEANUP;

      /* currently, we can't separate if there are too many variables (we could if we split the product in smaller
       * products, but this is not so simple with the current design) */
      if( nvars >= MAXMULTILINSEPALPSIZE )
         goto CLEANUP;

      exprhdlrdata = SCIPgetConsExprExprHdlrData(SCIPgetConsExprExprHdlrProduct(conshdlr));

      /* if LP is not constructed, build it now */
      if( exprhdlrdata->lpsize == 0 )
      {
         assert(exprhdlrdata->multilinearseparationlp == NULL);

         SCIP_CALL( buildMultilinearSeparationLP(scip, nvars, &(exprhdlrdata->multilinearseparationlp)) );
         exprhdlrdata->lpsize = nvars;
      }
      /* if LP is not large enough, rebuild it */
      else if( exprhdlrdata->lpsize < nvars )
      {
         assert(exprhdlrdata->multilinearseparationlp != NULL);
         SCIP_CALL( SCIPlpiFree(&(exprhdlrdata->multilinearseparationlp)) );

         SCIP_CALL( buildMultilinearSeparationLP(scip, nvars, &(exprhdlrdata->multilinearseparationlp)) );
         exprhdlrdata->lpsize = nvars;
      }
      assert(exprhdlrdata->multilinearseparationlp != NULL);
      assert(exprhdlrdata->lpsize >= nvars);

      SCIPdebugMsg(scip, "computing multilinear cut with %d variables, coef %g, midval %g\n", nvars, midval * exprdata->coefficient, midval );

      /* compute facet of fixmid * coefficient * \Pi_i vars[i] */
      SCIP_CALL( computeFacet(scip, exprhdlrdata->randnumgen, exprhdlrdata->multilinearseparationlp, sol, vars, nvars,
               auxvar, midval * exprdata->coefficient, overestimate, midval, fixedinterval, &violation, facet) );

      /* if we can't separate the point, return */
      if( SCIPisLE(scip, violation, 0.0) )
         goto CLEANUP;

      /* build cut; multilinear cuts are only valid locally */
      SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "multilinear", 0, NULL, NULL, -SCIPinfinity(scip),
               SCIPinfinity(scip), TRUE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddVarsToRow(scip, *cut, nvars, vars, facet) );
      SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );

      if( overestimate )
      {
         SCIP_CALL( SCIPchgRowLhs(scip, *cut, -facet[nvars]) );
      }
      else
      {
         SCIP_CALL( SCIPchgRowRhs(scip, *cut, -facet[nvars]) );
      }

CLEANUP:
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &facet);
   }

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaProduct)
{
   SCIP_ROW* cut;
   SCIP_Bool infeasible;

   *result = SCIP_DIDNOTFIND;
   *ncuts = 0;

   /* try to find a cut */
   SCIP_CALL( separatePointProduct(scip, conshdlr, expr, sol, &cut) );

   if( cut == NULL )
      return SCIP_OKAY;

   /* check whether its violation and numerical properties are ok (and maybe improve) */
   SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cut, sol, minviolation) );

   if( cut == NULL )
      return SCIP_OKAY;

   assert(-SCIPgetRowSolFeasibility(scip, cut, sol) >= minviolation);

   /* add cut */
   SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
   *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   *ncuts += 1;

#ifdef SCIP_DEBUG
   if( *result == infeasible )
   {
      SCIPdebugMsg(scip, "add cut makes node infeasible!\n");
   }
   else
   {
      SCIPdebugMsg(scip, "add cut with violation %e\n", -SCIPgetRowSolFeasibility(scip, cut, sol));
   }
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   /* release cut */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_INTERVAL childbounds;
   int i;
   int j;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) > 0);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   *nreductions = 0;
   *infeasible = FALSE;

   /* too expensive (runtime here is quadratic in number of children) */
   if( SCIPgetConsExprExprNChildren(expr) > 10 )
      return SCIP_OKAY;

   /* not possible to learn bounds if expression interval is unbounded in both directions */
   if( SCIPintervalIsEntire(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)) )
      return SCIP_OKAY;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   /* f = const * prod_k c_k => c_i = f / (const * prod_{j:j!=i} c_j ) */
   for( i = 0; i < SCIPgetConsExprExprNChildren(expr) && !(*infeasible); ++i )
   {
      SCIPintervalSet(&childbounds, exprdata->coefficient);

      /* compute prod_{j:j!=i} c_j */
      for( j = 0; j < SCIPgetConsExprExprNChildren(expr); ++j )
      {
         if( i == j )
            continue;

         SCIPintervalMul(SCIPinfinity(scip), &childbounds, childbounds,
               SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[j]));

         /* if there is 0.0 in the product, then later division will hardly give useful bounds, so give up for this i */
         if( childbounds.inf <= 0.0 && childbounds.sup >= 0.0 )
            break;
      }

      /* if the previous for finish, not because of the break */
      if( j == SCIPgetConsExprExprNChildren(expr) )
      {
         /* f / (const * prod_{j:j!=i} c_j) */
         SCIPintervalDiv(SCIPinfinity(scip), &childbounds, SCIPgetConsExprExprInterval(expr), childbounds);

         /* try to tighten the bounds of the expression */
         SCIP_CALL( SCIPtightenConsExprExprInterval(scip, SCIPgetConsExprExprChildren(expr)[i], childbounds, force,
               infeasible, nreductions) );
      }
   }

   return SCIP_OKAY;
}

/** creates the handler for product expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLRDATA* exprhdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* allocate expression handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &exprhdlrdata) );

   /* initialize all data in exprhdlrdata to 0/NULL */
   BMSclearMemory(exprhdlrdata);

   /* create/initialize random number generator */
   SCIP_CALL( SCIPrandomCreate(&exprhdlrdata->randnumgen, SCIPblkmem(scip),
         SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED)) );

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "prod", "product of children",
            PRODUCT_PRECEDENCE, evalProduct, exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrProduct, freehdlrProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataProduct, freedataProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffProduct) );

   return SCIP_OKAY;
}

/** creates a product expression */
SCIP_RETCODE SCIPcreateConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real             coefficient         /**< constant coefficient of product */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, coefficient) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrProduct(consexprhdlr), exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** gets the constant coefficient of a product expression */
SCIP_Real SCIPgetConsExprExprProductCoef(
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->coefficient;
}

/** appends an expression to a product expression */
SCIP_RETCODE SCIPappendConsExprExprProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< product expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
   )
{
   assert(expr != NULL);

   SCIP_CALL( SCIPappendConsExprExpr(scip, expr, child) );

   return SCIP_OKAY;
}
