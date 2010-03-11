/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpioracle.c,v 1.1 2010/03/11 11:22:32 bzfviger Exp $"

/**@file    nlpioracle.c
 * @brief   implementation of NLPI oracle interface
 * @author  Stefan Vigerske
 * 
 * @todo jacobi evaluation should be sparse
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpioracle.h"
#ifdef WITH_NL
#include "expression.h"
#include "exprinterpret.h"
#endif
#include "scip/pub_misc.h"

#include <string.h> /* for strlen */

/**@name NLPI Oracle data structures */
/**@{ */

struct SCIP_NlpiOracle
{
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SCIP_Real             infinity;           /**< value for infinity */
   
   int                   nvars;              /**< number of variables */
   SCIP_Real*            varlbs;             /**< array with variable lower bounds */
   SCIP_Real*            varubs;             /**< array with variable upper bounds */
   char**                varnames;           /**< array with variable names */
   
   int                   nconss;             /**< number of constraints */
   SCIP_Real*            conslhss;           /**< array with left-hand sides of constraints */
   SCIP_Real*            consrhss;           /**< array with right-hand sides of constraints */
   int*                  conslinlens;        /**< linlens[.] gives the length of linear part, or 0 if no linear part in this constraint */
   int**                 conslinidxs;        /**< linidxs[.] gives the variable indices in linear part, or NULL if no linear part in this constraint */
   SCIP_Real**           conslincoefs;       /**< lincoefs[.] gives the variable coefficients in linear part, or NULL if no linear part in this constraint */
   int*                  consquadlens;       /**< quadlens[.] gives length of quadratic part, or 0 if no quadratic part in this constraint */
   int**                 consquadrows;       /**< quadrows[.] gives row indices for quadratic part, or NULL if no quadratic part in this constraint */
   int**                 consquadcols;       /**< quadcols[.] gives column indices for quadratic part, or NULL if no quadratic part in this constraint */
   SCIP_Real**           consquadvals;       /**< quadvals[.] gives coefficients in quadratic part, or NULL if no quadratic part in this constraint */
   int**                 consexprvaridxs;    /**< exprvaridx[.] gives indices of variables from expression in NLP, or NULL if no nonquadratic part in this constraint */
   SCIP_EXPRTREE**       consexprtrees;      /**< exprtrees[.] gives nonquadratic part, or NULL if no nonquadratic part in this constraint */
   char**                consnames;          /**< array with constraint names */

   SCIP_Real             objconstant;        /**< constant part of objective */
   int                   objnlin;            /**< number of linear variable coefficients in objective */ 
   int*                  objlinidxs;         /**< array with indices of linear variables in objective, or NULL if no linear part */
   SCIP_Real*            objlinvals;         /**< array with coefficients of linear variables in objective, or NULL if no linear part */
   int                   objquadlen;         /**< length of quadratic part of objective, or 0 if no quadratic part in objective */
   int*                  objquadrows;        /**< array with row indices in quadratic part of objective, or NULL if no quadratic part in objective */
   int*                  objquadcols;        /**< array with column indices in quadratic part in objective, or NULL if no quadratic part in objective */ 
   SCIP_Real*            objquadvals;        /**< array with coefficients in quadratic part in objective, or NULL if no quadratic part in objective */
   int*                  objexprvaridxs;     /**< exprvaridxs[.] gives indices of variables from expression in NLP, or NULL if no nonquadratic part in objective */
   SCIP_EXPRTREE*        objexprtree;        /**< expression tree of nonquadratic part in objective, or NULL if no nonquadratic part in objective */

   int*                  jacoffsets;         /**< rowwise jacobi sparsity pattern: constraint offsets in jaccols */
   int*                  jaccols;            /**< rowwise jacobi sparsity pattern: indices of variables appearing in constraints */
   
   int*                  heslagoffsets;      /**< rowwise sparsity pattern of hessian matrix of Lagrangian: row offsets in heslagcol */
   int*                  heslagcols;         /**< rowwise sparsity pattern of hessian matrix of Lagrangian: column indices; sorted for each row */
   
   int*                  vardegrees;         /**< array with maximal degree of variable over objective and all constraints */
#ifdef WITH_NL
   SCIP_EXPRINT*         exprinterpreter;    /**< interpreter for expression trees: evaluation and derivatives */
#endif
};

/**@} */

/**@name Local functions */
/**@{ */

/** Invalidates the sparsity pattern of the Jacobian.
 *  Should be called when constraints are added or deleted.
 */
static
void invalidateJacobiSparsity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   if( oracle->jacoffsets == NULL )
   { /* nothing to do */
      assert(oracle->jaccols == NULL);
      return;
   }
   
   assert(oracle->jaccols != NULL);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->jaccols,    oracle->jacoffsets[oracle->nconss]);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->jacoffsets, oracle->nconss + 1);
}

/** Invalidates the sparsity pattern of the Hessian of the Lagragian.
 *  Should be called when the objective is set or constraints are added or deleted.
 */
static
void invalidateHessianLagSparsity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   if( oracle->heslagoffsets == NULL )
   { /* nothing to do */
      assert(oracle->heslagcols == NULL);
      return;
   }
   
   assert(oracle->heslagcols != NULL);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->heslagcols,    oracle->heslagoffsets[oracle->nvars]);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->heslagoffsets, oracle->nvars + 1);
}

/** Moves one constraint.
 * The place where it moves to need to be empty (all NULL) but allocated.
 */
static
SCIP_RETCODE moveConstraint(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to store NLPIORACLE data structure */
   int                   fromidx,            /**< index of constraint to move */
   int                   toidx               /**< index of place where to move constraint to */
   )
{
   assert(oracle != NULL);

   assert(0 <= fromidx);
   assert(0 <= toidx);
   assert(fromidx < oracle->nconss);
   assert(toidx   < oracle->nconss);
   
   assert(oracle->conslinlens     == NULL || oracle->conslinlens    [toidx] == 0);
   assert(oracle->consquadlens    == NULL || oracle->consquadlens   [toidx] == 0);
   assert(oracle->consexprvaridxs == NULL || oracle->consexprvaridxs[toidx] == 0);
   assert(oracle->consnames       == NULL || oracle->consnames      [toidx] == NULL);
   
   oracle->conslhss[toidx] = oracle->conslhss[fromidx];
   oracle->consrhss[toidx] = oracle->consrhss[fromidx];
   oracle->conslhss[fromidx] = -oracle->infinity;
   oracle->consrhss[fromidx] =  oracle->infinity;

   if( oracle->conslinlens )
   {
      oracle->conslinlens [toidx] = oracle->conslinlens [fromidx];
      oracle->conslinidxs [toidx] = oracle->conslinidxs [fromidx];
      oracle->conslincoefs[toidx] = oracle->conslincoefs[fromidx];
      
      oracle->conslinlens [fromidx] = 0;
      oracle->conslinidxs [fromidx] = NULL;
      oracle->conslincoefs[fromidx] = NULL;
   }

   if( oracle->consquadlens )
   {
      oracle->consquadlens[toidx] = oracle->consquadlens[fromidx];
      oracle->consquadrows[toidx] = oracle->consquadrows[fromidx];
      oracle->consquadcols[toidx] = oracle->consquadcols[fromidx];
      oracle->consquadvals[toidx] = oracle->consquadvals[fromidx];
      
      oracle->consquadlens[fromidx] = 0;
      oracle->consquadrows[fromidx] = NULL;
      oracle->consquadcols[fromidx] = NULL;
      oracle->consquadvals[fromidx] = NULL;
   }

   if( oracle->consexprtrees )
   {
      oracle->consexprtrees  [toidx] = oracle->consexprtrees  [fromidx];
      oracle->consexprvaridxs[toidx] = oracle->consexprvaridxs[fromidx];
      
      oracle->consexprtrees  [fromidx] = NULL;
      oracle->consexprvaridxs[fromidx] = NULL;
   }

   if( oracle->consnames )
   {
      oracle->consnames[toidx] = oracle->consnames[fromidx];
      
      oracle->consnames[fromidx] = NULL;
   }
   
   return SCIP_OKAY;
}

/** Frees one constraint.
 */
static
SCIP_RETCODE freeConstraint(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to store NLPIORACLE data structure */
   int                   considx             /**< index of constraint to free */
   )
{
   assert(oracle != NULL);

   assert(considx >= 0);
   assert(considx < oracle->nconss);
   
   oracle->conslhss[considx] = -oracle->infinity;
   oracle->consrhss[considx] =  oracle->infinity;
   
   if( oracle->conslinlens )
   {
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conslinidxs [considx], oracle->conslinlens[considx]);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conslincoefs[considx], oracle->conslinlens[considx]);
      oracle->conslinlens[considx] = 0;
   }
   
   if( oracle->consquadlens )
   {
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadrows[considx], oracle->consquadlens[considx]);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadcols[considx], oracle->consquadlens[considx]);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadvals[considx], oracle->consquadlens[considx]);
      oracle->consquadlens[considx] = 0;
   }
   
   if( oracle->consexprtrees && oracle->consexprtrees[considx] != NULL )
   {
#ifdef WITH_NL
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consexprvaridxs[considx], SCIPexprtreeGetNVars(oracle->consexprtrees[considx]));
      SCIP_CALL( SCIPexprtreeFree(&oracle->consexprtrees[considx]) );
#endif
   }

   if( oracle->consnames && oracle->consnames[considx] != NULL )
   {
      BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->consnames[considx], strlen(oracle->consnames[considx])+1);
   }
   
   return SCIP_OKAY;
}

/** Frees all constraints.
 */
static
SCIP_RETCODE freeConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   int i;
   
   assert(oracle != NULL);
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_CALL( freeConstraint(oracle, i) );
   }

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conslhss,        oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consrhss,        oracle->nconss);
      
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conslinlens,     oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conslinidxs,     oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conslincoefs,    oracle->nconss);

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadlens,    oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadrows,    oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadcols,    oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadvals,    oracle->nconss);

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consexprtrees,   oracle->nconss);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consexprvaridxs, oracle->nconss);
   
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consnames,       oracle->nconss);
      
   oracle->nconss = 0;
   
   return SCIP_OKAY;
}


/** Moves one variable.
 * The place where it moves to need to be empty (all NULL) but allocated.
 * Note that this function does not update the variable indices in the constraints!
 */
static
SCIP_RETCODE moveVariable(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to store NLPIORACLE data structure */
   int                   fromidx,            /**< index of variable to move */
   int                   toidx               /**< index of place where to move variable to */
   )
{
   assert(oracle != NULL);

   assert(0 <= fromidx);
   assert(0 <= toidx);
   assert(fromidx < oracle->nconss);
   assert(toidx   < oracle->nconss);

   assert(oracle->varlbs[toidx] <= -oracle->infinity);
   assert(oracle->varubs[toidx] >=  oracle->infinity);
   assert(oracle->varnames == NULL || oracle->varnames[toidx] == NULL);
   assert(oracle->vardegrees[toidx] == -1);
   
   oracle->varlbs[toidx] = oracle->varlbs[fromidx];
   oracle->varubs[toidx] = oracle->varubs[fromidx];
   
   oracle->varlbs[fromidx] = -oracle->infinity;
   oracle->varubs[fromidx] =  oracle->infinity;
   
   oracle->vardegrees[toidx]   = oracle->vardegrees[fromidx];
   oracle->vardegrees[fromidx] = -1;
   
   if( oracle->varnames != NULL )
   {
      oracle->varnames[toidx]   = oracle->varnames[fromidx];
      oracle->varnames[fromidx] = NULL;
   }
      
   return SCIP_OKAY;
}

/** Frees one variable.
 */
static
SCIP_RETCODE freeVariable(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to store NLPIORACLE data structure */
   int                   varidx              /**< index of variable to free */
   )
{
   assert(oracle != NULL);

   assert(varidx >= 0);
   assert(varidx < oracle->nvars);

   oracle->varlbs[varidx] = -oracle->infinity;
   oracle->varubs[varidx] =  oracle->infinity;

   oracle->vardegrees[varidx] = -1;
   
   if( oracle->varnames && oracle->varnames[varidx] != NULL )
   {
      BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varnames[varidx], strlen(oracle->varnames[varidx])+1);
   }

   return SCIP_OKAY;
}

/** Frees all variables.
 */
static
SCIP_RETCODE freeVariables(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   int i;
   
   assert(oracle != NULL);
   
   for( i = 0; i < oracle->nvars; ++i )
   {
      SCIP_CALL( freeVariable(oracle, i) );
   }
   
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varlbs,     oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varubs,     oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varnames,   oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->vardegrees, oracle->nvars);
   
   oracle->nvars = 0;
   
   return SCIP_OKAY;
}

/** applies a mapping of indices to one array of indices */
static
void mapIndices(
   int*                  indexmap,           /**< mapping from old variable indices to new indices */
   int                   nindices,           /**< number of indices in indices1 and indices2 */
   int*                  indices             /**< array of indices to adjust */
   )
{
   assert(indexmap != NULL);
   assert(nindices == 0 || indices != NULL);

   for( ; nindices ; --nindices, ++indices )
   {
      assert(indexmap[*indices] >= 0);
      *indices = indexmap[*indices];
   }
}

/** applies a mapping of indices to two array of indices (of the same length) */
static
void mapIndices2(
   int*                  indexmap,           /**< mapping from old variable indices to new indices */
   int                   nindices,           /**< number of indices in indices1 and indices2 */
   int*                  indices1,           /**< first array of indices to adjust */
   int*                  indices2            /**< second array of indices to adjust */
   )
{
   assert(indexmap != NULL);
   assert(nindices == 0 || indices1 != NULL);
   assert(nindices == 0 || indices2 != NULL);

   for( ; nindices ; --nindices, ++indices1, ++indices2 )
   {
      assert(indexmap[*indices1] >= 0);
      assert(indexmap[*indices2] >= 0);
      *indices1 = indexmap[*indices1];
      *indices2 = indexmap[*indices2];
   }
}

/** computes the value of a function */
static
SCIP_RETCODE evalFunctionValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nlin,               /**< length of linear part */
   const int*            lininds,            /**< indices of linear variables */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables */
   int                   quadlen,            /**< length of quadratic part matrix */
   int*                  quadrows,           /**< indices of rows in quadratic part matrix */
   int*                  quadcols,           /**< indices of columns in quadratic part matrix */
   SCIP_Real*            quadvals,           /**< coefficients in quadratic part matrix */
   int*                  exprvaridxs,        /**< indices of variables in nonquadratic part */
   SCIP_EXPRTREE*        exprtree,           /**< nonquadratic part */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Real*            val                 /**< pointer to store function value */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(val != NULL);
   
   *val = 0.0;
   
   if( nlin != 0 )
   {
      assert(lininds != NULL);
      assert(linvals != NULL);
      assert(x != NULL);
      
      for( ; nlin > 0; --nlin, ++lininds, ++linvals )
         *val += *linvals * x[*lininds]; 
   }
   
   if( quadlen != 0 )
   {
      assert(quadrows != NULL);
      assert(quadcols != NULL);
      assert(quadvals != NULL);
      assert(x != NULL);
      
      for( ; quadlen > 0; --quadlen, ++quadrows, ++quadcols, ++quadvals )
         *val += *quadvals * x[*quadrows] * x[*quadcols];
   }
   
   if( exprtree != NULL )
   {
#ifdef WITH_NL
      SCIP_Real* xx;
      int        i;
      SCIP_Real  nlval;
      int        nvars;
      
      nvars = SCIPexprtreeGetNVars(exprtree);
      
      if( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) == NULL )
         return SCIP_NOMEMORY;
      for (i = 0; i < nvars; ++i)
      {
         assert(exprvaridxs[i] >= 0);
         xx[i] = x[exprvaridxs[i]];
      }
      
      SCIP_CALL( SCIPexprintEval(oracle->exprinterpreter, exprtree, xx, NULL, &nlval) );
      if( nlval != nlval || ABS(nlval) >= -oracle->infinity )
         *val = nlval;
      else
         *val += nlval;
      
      BMSfreeBlockMemoryArray(oracle->blkmem, &xx, nvars);
#else
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
#endif
   }

   return SCIP_OKAY;
}

/** computes the value and gradient of a function */
static
SCIP_RETCODE evalFunctionGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nlin,               /**< length of linear part */
   const int*            lininds,            /**< indices of linear variables */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables */
   int                   quadlen,            /**< length of quadratic part matrix */
   int*                  quadrows,           /**< indices of rows in quadratic part matrix */
   int*                  quadcols,           /**< indices of columns in quadratic part matrix */
   SCIP_Real*            quadvals,           /**< coefficients in quadratic part matrix */
   int*                  exprvaridxs,        /**< indices of variables in nonquadratic part */
   SCIP_EXPRTREE*        exprtree,           /**< nonquadratic part */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */
   SCIP_Real*            val,                /**< pointer to store function value */
   SCIP_Real*            grad                /**< pointer to store function gradient */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(val != NULL);
   assert(grad != NULL);
   
   *val = 0.0;
   BMSclearMemoryArray(grad, oracle->nvars);
   
   if( nlin != 0 )
   {
      assert(lininds != NULL);
      assert(linvals != NULL);
      assert(x != NULL);
      
      for( ; nlin > 0; --nlin, ++lininds, ++linvals )
      {
         *val += *linvals * x[*lininds];
         assert(grad[*lininds] == 0.0);   /* we do not like duplicate indices */
         grad[*lininds] = *linvals;
      }
   }
   
   if( quadlen != 0 )
   {
      SCIP_Real tmp;
      assert(quadrows != NULL);
      assert(quadcols != NULL);
      assert(quadvals != NULL);
      assert(x != NULL);
      
      for( ; quadlen > 0; --quadlen, ++quadrows, ++quadcols, ++quadvals )
      {
         tmp = *quadvals * x[*quadrows];
         *val += tmp * x[*quadcols];
         grad[*quadcols] += tmp;
         grad[*quadrows] += *quadvals * x[*quadcols];
      }
   }

   if( exprtree != NULL )
   {
#ifdef WITH_NL
      SCIP_Real* xx = NULL;
      SCIP_Real* g;
      int        i;
      SCIP_Real  nlval;
      int        nvars;
      
      nvars = SCIPexprtreeGetNVars(exprtree);

      if( BMSallocBlockMemoryArray(oracle->blkmem, &g, nvars) == NULL )
         return SCIP_NOMEMORY;

      if (newx)
      {
         if( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) == NULL )
            return SCIP_NOMEMORY;
         for (i = 0; i < nvars; ++i)
         {
            assert(exprvaridxs[i] >= 0);
            xx[i] = x[exprvaridxs[i]];
         }
      }
      
      SCIP_CALL( SCIPexprintGradDense(oracle->exprinterpreter, exprtree, xx, newx, NULL, &nlval, g) );
      if( nlval != nlval || ABS(nlval) >= oracle->infinity )
      {
         SCIPfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
         SCIPfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
         SCIPdebugMessage("gradient evaluation yield invalid function value %g\n", nlval);
         return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
      }
      else
      {
         *val += nlval;
         for (i = 0; i < nvars; ++i)
            if (g[i] != g[i])
            {
               SCIPdebugMessage("gradient evaluation yield invalid gradient value %g\n", g[i]);
               SCIPfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
               SCIPfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
               return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
            }
            else
               grad[exprvaridxs[i]] += g[i];
      }
      
      SCIPfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
      SCIPfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
#else
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
#endif
   }

   return SCIP_OKAY;
}

/** calculate memory size for dynamically allocated arrays (copied from scip/set.c) */
static
int calcGrowSize(
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
   size = 4;
   while( size < num )
      size = (int)(1.2 * size + 4);

   return size;
}

/** increases the size of an an array of integers so that it can contain at least one more element than in its current size */
static
SCIP_RETCODE increaseIntArraySize(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int**                 intarray,           /**< array of integers */
   int*                  len                 /**< length of array (may be increased) */
   )
{
   int oldsize;
   
   assert(intarray != NULL);
   assert(len != NULL);

   oldsize = *len;
   *len = calcGrowSize(*len + 1);

   if( *intarray == NULL )
   {
      if( BMSallocBlockMemoryArray(blkmem, intarray, *len) == NULL )
         return SCIP_NOMEMORY;
   }
   else
   {
      if( BMSreallocBlockMemoryArray(blkmem, intarray, oldsize, *len) == NULL )
         return SCIP_NOMEMORY;
   }
   
   return SCIP_OKAY;
}

/** collects nonzeros entries in colnz and increases the nzcount given indices of quadratic terms */
static
SCIP_RETCODE hessLagSparsitySetNzFlagForQuad(
   SCIP_NLPIORACLE*      oracle,     /**< NLPI oracle */
   int**                 colnz,              /**< indices of nonzero entries for each column */
   int*                  collen,             /**< space allocated to store indices of nonzeros for each column */
   int*                  colnnz,             /**< number of nonzero entries for each column */
   int*                  nzcount,            /**< counter for total number of nonzeros; should be increased whenever some colnnz is increased */
   int*                  rowidx,             /**< row indices */
   int*                  colidx,             /**< column indices */
   int                   length              /**< length of quadratic part */
   )
{
   int pos;

   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(rowidx != NULL);
   assert(colidx != NULL);
   assert(length >= 0);

   for( ; length > 0; --length, ++rowidx, ++colidx )
   {
      assert(*colidx <= *rowidx);
      
      if( colnz[*rowidx] == NULL || !SCIPsortedvecFindInt(colnz[*rowidx], *colidx, colnnz[*rowidx], &pos) )
      {
         if( colnnz[*rowidx] >= collen[*rowidx] )
         { /* allocate more space so that we can insert one more element */
            SCIP_CALL( increaseIntArraySize(oracle->blkmem, &colnz[*rowidx], &collen[*rowidx]) );
         }
         assert(collen[*rowidx] > colnnz[*rowidx]);
         
         SCIPsortedvecInsertInt(colnz[*rowidx], *colidx, &colnnz[*rowidx]);
         ++(*nzcount);
      }
   }
   
   return SCIP_OKAY;
}

#ifdef WITH_NL
static
SCIP_RETCODE hessLagSparsitySetNzFlagForExprtree(
   SCIP_NLPIORACLE*      oracle,     /**< NLPI oracle */
   int**                 colnz,      /**< indices of nonzero entries for each column */
   int*                  collen,     /**< space allocated to store indices of nonzeros for each column */
   int*                  colnnz,     /**< number of nonzero entries for each column */
   int*                  nzcount,    /**< counter for total number of nonzeros; should be increased when nzflag is set to 1 the first time */
   int*                  exprvaridx, /**< indices of variables from expression tree in NLP */
   SCIP_EXPRTREE*        exprtree,   /**< expression tree */
   int                   dim         /**< dimension of matrix */
)
{
   SCIP_Real*  x;
   SCIP_Bool*  hesnz;
   int         i, j, n, nn, row, col, pos;
   
   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(exprvaridx != NULL);
   assert(exprtree != NULL);
   assert(dim >= 0);
   
   n  = SCIPexprtreeGetNVars(exprtree);
   nn = n*n;
   
   if( BMSallocBlockMemoryArray(oracle->blkmem, &x, n) == NULL )
      return SCIP_NOMEMORY;
   if( BMSallocBlockMemoryArray(oracle->blkmem, &hesnz, nn) == NULL )
      return SCIP_NOMEMORY;
   
   for (i = 0; i < n; ++i)
      x[i] = 2.0; /* hope that this value does not make much trouble for the evaluation routines */
   
   SCIP_CALL( SCIPexprintHessianSparsityDense(oracle->exprinterpreter, exprtree, x, NULL, hesnz) );
   
   for (i = 0; i < n; ++i) /* rows */
      for (j = 0; j <= i; ++j) /* cols */
      {
         if (!hesnz[i*n + j])
            continue;
         
         row = MAX(exprvaridx[i], exprvaridx[j]);
         col = MIN(exprvaridx[i], exprvaridx[j]);
         
         assert(row < dim);
         assert(col <= row);

         if( colnz[row] == NULL || !SCIPsortedvecFindInt(colnz[row], col, colnnz[row], &pos) )
         {
            if( colnnz[row] >= collen[row] )
            { /* allocate more space so that we can insert one more element */
               SCIP_CALL( increaseIntArraySize(oracle->blkmem, &colnz[row], &collen[row]) );
            }
            assert(collen[row] > colnnz[row]);
            
            SCIPsortedvecInsertInt(colnz[row], col, &colnnz[row]);
            ++(*nzcount);
         }
      }
   
   BMSfreeBlockMemoryArray(oracle->blkmem, &x, n);
   BMSfreeBlockMemoryArray(oracle->blkmem, &hesnz, nn);
   
   return SCIP_OKAY;
}
#endif

/** adds quadratic part of a constraint into hessian structure */
static
SCIP_RETCODE hessLagAddQuad(
   SCIP_Real             weight,             /**< weight of quadratic part */
   int*                  rowidx,             /**< row indices in matrix of quadratic part */
   int*                  colidx,             /**< column indices in matrix of quadratic part */
   SCIP_Real*            coeff,              /**< coefficients in matrix of quadratic part */
   int                   length,             /**< number of coefficients */
   int*                  hesoffset,          /**< row offsets in sparse matrix that is to be filled */ 
   int*                  hescol,             /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values              /**< buffer for values of sparse matrix that is to be filled */
   )
{
   int idx;

   assert(rowidx != NULL);
   assert(colidx != NULL);
   assert(length >= 0);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   for( ; length > 0; --length, ++rowidx, ++colidx, ++coeff )
   {
      assert(*colidx <= *rowidx);
      if( !SCIPsortedvecFindInt(&hescol[hesoffset[*rowidx]], *colidx, hesoffset[*rowidx+1] - hesoffset[*rowidx], &idx) )
      {
         SCIPerrorMessage("Could not find entry in hessian sparsity\n");
         return SCIP_ERROR;
      }
      values[hesoffset[*rowidx] + idx] += weight * ((*colidx == *rowidx) ? 2 * *coeff : *coeff);
   }

   return SCIP_OKAY;
}

#ifdef WITH_NL
static
SCIP_RETCODE hessLagAddExprtree(
   SCIP_NLPIORACLE*      oracle,     /**< oracle */
   SCIP_Real             weight,     /**< weight of quadratic part */
   const SCIP_Real*      x,          /**< point for which hessian should be returned */
   SCIP_Bool             new_x,      /**< whether point has been evaluated before */
   int*                  exprvaridx, /**< NLP indices for variables in expression tree */
   SCIP_EXPRTREE*        exprtree,   /**< expression tree */
   int*                  hesoffset,  /**< row offsets in sparse matrix that is to be filled */ 
   int*                  hescol,     /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values      /**< buffer for values of sparse matrix that is to be filled */
)
{
   SCIP_Real* xx = NULL;
   SCIP_Real* h;
   SCIP_Real* hh;
   int        i, j, n, nn, row, col, idx;
   SCIP_Real  val;
   
   assert(oracle != NULL);
   assert(x != NULL || new_x == FALSE);
   assert(exprvaridx != NULL);
   assert(exprtree != NULL);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   n = SCIPexprtreeGetNVars(exprtree);
   nn = n*n;

   if( BMSallocBlockMemoryArray(oracle->blkmem, &h, nn) == NULL )
      return SCIP_NOMEMORY;

   if (new_x)
   {
      if( BMSallocBlockMemoryArray(oracle->blkmem, &xx, n) == NULL )
         return SCIP_NOMEMORY;
      for (i = 0; i < n; ++i)
      {
         assert(exprvaridx[i] >= 0);
         xx[i] = x[exprvaridx[i]];
      }
   }
   
   SCIP_CALL( SCIPexprintHessianDense(oracle->exprinterpreter, exprtree, xx, new_x, NULL, &val, h) );
   if (val != val)
   {
      SCIPdebugMessage("hessian evaluation yield invalid function value %g\n", val);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, n);
      BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
      return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
   }

   hh = h;
   for (i = 0; i < n; ++i) /* rows */
   {
      for (j = 0; j <= i; ++j, ++hh) /* cols */
      {
         if (!*hh)
            continue;
         
         if (*hh != *hh)
         {
            SCIPdebugMessage("hessian evaluation yield invalid hessian value %g\n", *hh);
            BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, n);
            BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
            return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
         }
         
         row = MAX(exprvaridx[i], exprvaridx[j]);
         col = MIN(exprvaridx[i], exprvaridx[j]);
         
         if (!SCIPsortedvecFindInt(&hescol[hesoffset[row]], col, hesoffset[row+1] - hesoffset[row], &idx))
         {
            SCIPerrorMessage("Could not find entry (%d, %d) in hessian sparsity\n", row, col);
            BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, n);
            BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
            return SCIP_ERROR;
         }
         values[hesoffset[row] + idx] += weight * *hh;
      }
      hh += (n-j);
   }
   
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, n);
   BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);

   return SCIP_OKAY;
}
#endif

/** prints a function */ 
static
SCIP_RETCODE printFunction(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file,               /**< file to print to, has to be not NULL */
   int                   nlin,               /**< length of linear part */
   const int*            lininds,            /**< indices of linear variables */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables */
   int                   quadlen,            /**< number of indices in quadratic part matrix */
   int*                  quadrows,            /**< indices of rows in quadratic part matrix */
   int*                  quadcols,            /**< indices of columns in quadratic part matrix */
   SCIP_Real*            quadvals,            /**< coefficients in quadratic part matrix */
   SCIP_EXPRTREE*        exprtree            /**< nonquadratic part */
   )
{  /*lint --e{715}*/
   int i;
   int j;

   assert(oracle != NULL);
   assert(file != NULL);
   
   for( i = 0; i < nlin; ++i )
   {
      SCIPmessageFPrintInfo(file, "%+g*", linvals[i]);
      if( oracle->varnames != NULL && oracle->varnames[lininds[i]] != NULL )
         SCIPmessageFPrintInfo(file, oracle->varnames[lininds[i]]);
      else
         SCIPmessageFPrintInfo(file, "x%d", lininds[i]);
   }
   
   if( quadlen != 0 )
   {
      j = 0;
      for( i = 0; i < quadlen; ++i )
      {
         SCIPmessageFPrintInfo(file, "%+g*", quadvals[j]);
         if( oracle->varnames != NULL && oracle->varnames[quadrows[i]] != NULL )
            SCIPmessageFPrintInfo(file, oracle->varnames[quadrows[i]]);
         else
            SCIPmessageFPrintInfo(file, "x%d", quadrows[i]);
         SCIPmessageFPrintInfo(file, "*");
         if( oracle->varnames != NULL && oracle->varnames[quadcols[i]] != NULL )
            SCIPmessageFPrintInfo(file, oracle->varnames[quadcols[i]]);
         else
            SCIPmessageFPrintInfo(file, "x%d", quadcols[i]);
      }
   }

   if (exprtree)
   {
      fprintf(file, " +");
#ifdef WITH_NL
      SCIPexprtreePrint(file, exprtree);
#endif
   }

   return SCIP_OKAY;
}

/**@} */

/**@name public function */
/**@{ */

/** creates an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(blkmem != NULL);
   assert(oracle != NULL);
   
   if( BMSallocMemory(oracle) != NULL )
      return SCIP_NOMEMORY;
   BMSclearMemory(*oracle);
   
   (*oracle)->blkmem   = blkmem;
   (*oracle)->infinity = SCIP_DEFAULT_INFINITY;
   
#ifdef WITH_NL
   SCIPmessageFPrintInfo(NULL, "Oracle initializes expression interpreter %s\n", SCIPexprintGetName());
   SCIP_CALL( SCIPexprintCreate(blkmem, &(*oracle)->exprinterpreter) );
#endif

   return SCIP_OKAY;
}

/** frees an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle  != NULL);
   assert(*oracle != NULL);
   
   invalidateJacobiSparsity(*oracle);
   invalidateHessianLagSparsity(*oracle);
   
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objlinidxs,  (*oracle)->objnlin);
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objlinvals,  (*oracle)->objnlin);
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objquadrows, (*oracle)->objquadlen);
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objquadcols, (*oracle)->objquadlen);
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objquadvals, (*oracle)->objquadlen);
   if( (*oracle)->objexprtree != NULL )
   {
#ifdef WITH_NL
      BMSfreeBlockMemoryArray((*oracle)->blkmem, &(*oracle)->objexprvaridxs, SCIPexprtreeGetNVars((*oracle)->objexprtree));
      SCIP_CALL( SCIPexprtreeFree(&(*oracle)->objexprtree) );
#endif
   }

   SCIP_CALL( freeConstraints(*oracle) );
   SCIP_CALL( freeVariables(*oracle) );

#ifdef WITH_NL
   SCIP_CALL( SCIPexprintFree(&(*oracle)->exprinterpreter) );
#endif

   BMSfreeMemory(oracle);
   
   return SCIP_OKAY;
}

/** sets the value for infinity */
SCIP_RETCODE SCIPnlpiOracleSetInfinity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             infinity            /**< value to use for infinity */
   )
{
   assert(oracle != NULL);
   
   if( infinity < 0.0 )
      return SCIP_PARAMETERWRONGVAL;
   
   oracle->infinity = infinity;
   
   return SCIP_OKAY;
}

/** gets the value for infinity */
SCIP_Real SCIPnlpiOracleGetInfinity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->infinity;
}

/** adds variables */
SCIP_RETCODE SCIPnlpiOracleAddVars(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to add */
   const SCIP_Real*      lbs,                /**< array with lower bounds of new variables, or NULL if all -infinity */
   const SCIP_Real*      ubs,                /**< array with upper bounds of new variables, or NULL if all +infinity */
   const char**          varnames            /**< array with names of new variables, or NULL if no names should be stored */
   )
{
   int i;

   assert(oracle != NULL);
   
   if( nvars == 0 )
      return SCIP_OKAY;
   
   assert(nvars > 0);
   
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varlbs, oracle->nvars, oracle->nvars + nvars) == NULL )
      return SCIP_NOMEMORY;
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varubs, oracle->nvars, oracle->nvars + nvars) == NULL )
      return SCIP_NOMEMORY;
   
   if( lbs != NULL )
   {
      BMScopyMemoryArray(&((oracle->varlbs)[oracle->nvars]), lbs, nvars);
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varlbs[oracle->nvars+i] = -oracle->infinity;
   
   if( ubs != NULL )
   {
      BMScopyMemoryArray(&((oracle->varubs)[oracle->nvars]), ubs, nvars);
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varubs[oracle->nvars+i] =  oracle->infinity;
   
   if( oracle->varnames != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varnames, oracle->nvars, oracle->nvars + nvars) == NULL )
         return SCIP_NOMEMORY;
      if( varnames != NULL )
      {
         for( i = 0; i < nvars; ++i )
         {
            if( varnames[i] != NULL )
            {
               if( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->varnames)[oracle->nvars+i]), varnames[i], strlen(varnames[i])+1) == NULL )
                  return SCIP_NOMEMORY;
            }
            else
               oracle->varnames[oracle->nvars+i] = NULL;
         }
      }
      else
      {
         BMSclearMemoryArray(&((oracle->varnames)[oracle->nvars]), nvars);
      }
   }
   else if( varnames != NULL )
   {
      if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->varnames), oracle->nvars + nvars) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray(oracle->varnames, oracle->nvars);
      for( i = 0; i < nvars; ++i )
      {
         if( varnames[i] != NULL )
         {
            if( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->varnames)[oracle->nvars+i]), varnames[i], strlen(varnames[i])+1) == NULL )
               return SCIP_NOMEMORY;
         }
         else
            oracle->varnames[oracle->nvars+i] = NULL;
      }
   }
   
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->vardegrees), oracle->nvars, oracle->nvars + nvars) == NULL )
      return SCIP_NOMEMORY;
   BMSclearMemoryArray(&((oracle->vardegrees)[oracle->nvars]), nvars);

   /* @TODO update sparsity pattern by extending heslagoffsets */
   invalidateHessianLagSparsity(oracle);

   oracle->nvars += nvars;

   return SCIP_OKAY;
}

/** adds constraints 
 * 
 *  linear coefficients: row(=constraint) oriented matrix;
 *  quadratic coefficiens: row oriented matrix for each constraint
 */
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   const int*            nquadrows,          /**< NULL if no quadratic parts, otherwise nquadrows[.] gives the number of columns in the matrix 
                                              *   of the quadratic part, or 0 if no quadratic part in this constraint */
   int* const*           quadrowidxs,        /**< NULL if no quadratic parts, otherwise quadrowidx[.] gives the indices of variables 
                                              *   for which a quadratic part is specified, or NULL if no quadratic part in this constraint */
   int* const*           quadoffsets,        /**< NULL if no quadratic parts, otherwise quadoffsets[.] gives start index of each rows quadratic coefficients 
                                              *   in quadidxs[.] and quadvals[.] (quadoffsets[.][nvars] gives length of quadidxs[.] and quadvals[.]), 
                                              *   or NULL if no quadratic part in this constraint */
   int* const*           quadidxs,           /**< NULL if no quadratic parts, otherwise quadidxs[.] gives column indices for quadratic part, 
                                              *   or NULL if no quadratic part in this constraint */
   SCIP_Real* const*     quadvals,           /**< NULL if no quadratic parts, otherwise quadvals[.] gives coefficients in quadratic part, 
                                              *   or NULL if no quadratic part in this constraint */
   int* const*           exprvaridxs,        /**< NULL if no nonquadratic parts, otherwise epxrvaridxs[.] maps variable indices in expression tree to indices in nlp */
   SCIP_EXPRTREE* const* exprtrees,          /**< NULL if no nonquadratic parts, otherwise exprtrees[.] gives nonquadratic part, 
                                              *   or NULL if no nonquadratic part in this constraint */
   const char**          consnames           /**< names of new constraints, or NULL if no names should be stored */
   )
{  /*lint --e{715}*/
   SCIP_Bool addednlcon;  /* whether a nonlinear constraint was added */
   int i, j;

   assert(oracle != NULL);
   
   if( nconss == 0 )
      return SCIP_OKAY;
   
   assert(nconss > 0);

   addednlcon = FALSE;

   invalidateJacobiSparsity(oracle); /* @TODO we could also update (extend) the sparsity pattern */

   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslhss, oracle->nconss, oracle->nconss + nconss) == NULL )
      return SCIP_NOMEMORY;
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consrhss, oracle->nconss, oracle->nconss + nconss) == NULL )
      return SCIP_NOMEMORY;

   if( lhss != NULL )
   {
      BMScopyMemoryArray(&((oracle->conslhss)[oracle->nconss]), lhss, nconss);
   }
   else
      for( i = 0; i < nconss; ++i )
         oracle->conslhss[oracle->nconss+i] = -oracle->infinity;
   
   if( rhss != NULL )
   {
      BMScopyMemoryArray(&((oracle->consrhss)[oracle->nconss]), rhss, nconss);
   }
   else
      for( i = 0; i < nconss; ++i )
         oracle->consrhss[oracle->nconss+i] =  oracle->infinity;

   if( nlininds != NULL )
   {
      assert(lininds != NULL);
      assert(linvals != NULL);

      if( oracle->conslinlens != NULL || oracle->nconss == 0 )
      {
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
      }
      else
      { /* had no linear parts so far */
         if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->conslinlens,  oracle->nconss);
         BMSclearMemoryArray(oracle->conslinidxs,  oracle->nconss);
         BMSclearMemoryArray(oracle->conslincoefs, oracle->nconss);
      }
      
      for( i = 0; i < nconss; ++i )
      {
         oracle->conslinlens[oracle->nconss+i] = nlininds[i];
         if( oracle->conslinlens[oracle->nconss+i] )
         {
            if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs [oracle->nconss+i], lininds[i], oracle->conslinlens[oracle->nconss+i]) == NULL )
               return SCIP_NOMEMORY;
            if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs[oracle->nconss+i], linvals[i], oracle->conslinlens[oracle->nconss+i]) == NULL )
               return SCIP_NOMEMORY;
            SCIPsortIntReal(oracle->conslinidxs[oracle->nconss+i], oracle->conslincoefs[oracle->nconss+i], oracle->conslinlens[oracle->nconss+i]);
            for( j = 0; j < nlininds[i]; ++j )
            {
               assert(lininds[i][j] >= 0);
               assert(lininds[i][j] <= oracle->nvars);
               oracle->vardegrees[lininds[i][j]] = MAX(1, oracle->vardegrees[lininds[i][j]]);
            }
         }
         else
         {
            oracle->conslinidxs [oracle->nconss+i] = NULL;
            oracle->conslincoefs[oracle->nconss+i] = NULL;
         }
      }
   }
   else if( oracle->conslinlens != NULL )
   { /* no new linear parts, but we had linear parts before */
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray(&oracle->conslinlens [oracle->nconss], nconss);
      BMSclearMemoryArray(&oracle->conslinidxs [oracle->nconss], nconss);
      BMSclearMemoryArray(&oracle->conslincoefs[oracle->nconss], nconss);
   }
   
   if( quadoffsets != NULL )
   {
      assert(nquadrows != NULL);
      assert(quadrowidxs != NULL);
      assert(quadidxs != NULL);
      assert(quadvals != NULL);
      if( (oracle->consquadlens == NULL) && (oracle->nconss != 0) )
      {
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens), oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadrows), oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadcols), oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadvals), oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->consquadlens, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadrows, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadcols, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadvals, oracle->nconss);
      }
      else
      {
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens), oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadrows), oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadcols), oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadvals), oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
      }
      for( i = 0; i < nconss; ++i )
      {
         if( quadoffsets[i] != NULL && nquadrows[i] != 0 )
         {
            int k;

            assert(quadrowidxs[i] != NULL);
            assert(quadidxs[i] != NULL);
            assert(quadvals[i] != NULL);
            oracle->consquadlens[oracle->nconss + i] = quadoffsets[i][nquadrows[i]];
            if( BMSallocBlockMemoryArray(oracle->blkmem, &((oracle->consquadrows)[oracle->nconss + i]), quadoffsets[i][nquadrows[i]]) == NULL )
               return SCIP_NOMEMORY;
            if( BMSallocBlockMemoryArray(oracle->blkmem, &((oracle->consquadcols)[oracle->nconss + i]), quadoffsets[i][nquadrows[i]]) == NULL )
               return SCIP_NOMEMORY;
            if( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->consquadvals)[oracle->nconss + i]), quadvals[i], quadoffsets[i][nquadrows[i]]) == NULL )
               return SCIP_NOMEMORY;
            
            k = quadoffsets[i][0];
            for( j = 0; j < nquadrows[i]; ++j )
            {
               if( quadoffsets[i][j] == quadoffsets[i][j+1] )
                  continue;
               for( ; k < quadoffsets[i][j+1]; ++k )
               {
                  assert(quadidxs[i][k] < nquadrows[i]);
                  if( quadrowidxs[i][j] > quadrowidxs[i][quadidxs[i][k]] )
                  {
                     oracle->consquadrows[oracle->nconss + i][k] = quadrowidxs[i][j];
                     oracle->consquadcols[oracle->nconss + i][k] = quadrowidxs[i][quadidxs[i][k]];
                  }
                  else
                  {
                     oracle->consquadrows[oracle->nconss + i][k] = quadrowidxs[i][quadidxs[i][k]];
                     oracle->consquadcols[oracle->nconss + i][k] = quadrowidxs[i][j];
                  }
                  oracle->vardegrees[quadrowidxs[i][quadidxs[i][k]]] = MAX(2, oracle->vardegrees[quadrowidxs[i][quadidxs[i][k]]]);
               }
               oracle->vardegrees[quadrowidxs[i][j]] = MAX(2, oracle->vardegrees[quadrowidxs[i][j]]);
            }
            assert(k == oracle->consquadlens[oracle->nconss + i]);
            
            /* sort quadrows and quadcols */
            SCIPsortIntIntReal(oracle->consquadrows[oracle->nconss+i], oracle->consquadcols[oracle->nconss+i], oracle->consquadvals[oracle->nconss+i], oracle->consquadlens[oracle->nconss+i]);
            j = 0;
            k = 0;
            while( j < oracle->consquadlens[oracle->nconss+i] )
            {
               while( k < oracle->consquadlens[oracle->nconss+i] && oracle->consquadrows[oracle->nconss+i][j] == oracle->consquadrows[oracle->nconss+i][k] )
                  ++k;
               SCIPsortIntReal(&oracle->consquadcols[oracle->nconss+i][j], &oracle->consquadvals[oracle->nconss+i][j], k-j);
               j = k;
            }
            
            addednlcon = TRUE;
         }
         else
         {
            oracle->consquadlens[oracle->nconss + i] = 0;
            oracle->consquadrows[oracle->nconss + i] = NULL;
            oracle->consquadcols[oracle->nconss + i] = NULL;
            oracle->consquadvals[oracle->nconss + i] = NULL;
         }
      }
   }
   else if( oracle->consquadlens != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens), oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadrows), oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadcols), oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadvals), oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray(&((oracle->consquadlens)[oracle->nconss]), nconss);
      BMSclearMemoryArray(&((oracle->consquadrows)[oracle->nconss]), nconss);
      BMSclearMemoryArray(&((oracle->consquadcols)[oracle->nconss]), nconss);
      BMSclearMemoryArray(&((oracle->consquadvals)[oracle->nconss]), nconss);
   }
   
   if( exprtrees != NULL )
   {
      if( oracle->consexprtrees == NULL && oracle->nconss != 0 )
      {
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprtrees),   oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->consexprtrees,   oracle->nconss);
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprvaridxs), oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->consexprvaridxs, oracle->nconss);
      }
      else
      {
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprtrees),   oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprvaridxs), oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
      }
      for( i = 0; i < nconss; ++i )
         if( exprtrees[i] != NULL )
         {
#ifdef WITH_NL
            assert(oracle->exprinterpreter != NULL);
            assert(SCIPexprtreeHasVarsAsIndex(exprtrees[i]));
            
            SCIP_CALL( SCIPexprtreeCopy(&oracle->consexprtrees[oracle->nconss + i], exprtrees[i]) );
            
            SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->consexprtrees[oracle->nconss + i]) );
            
            if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs[oracle->nconss + i], exprvaridxs[i], SCIPexprtreeGetNVars(exprtrees[i])) == NULL )
               return SCIP_NOMEMORY;
            for (j = 0; j < SCIPexprtreeGetNVars(exprtrees[i]); ++j)
            {
               assert(exprvaridxs[i][j] >= 0);
               assert(exprvaridxs[i][j] <  oracle->nvars);
               oracle->vardegrees[exprvaridxs[i][j]] = INT_MAX; /* @TODO could try to be more clever, maybe use getMaxDegree function in exprtree */
            }
#else
            SCIPerrorMessage("nonquadratic functions not supported in NLPI yet.\n");
            return SCIP_ERROR;
#endif
         }
         else
         {
            oracle->consexprtrees[oracle->nconss + i] = NULL;
            oracle->consexprvaridxs[oracle->nconss + i] = NULL;
         }
   }
   else if( oracle->consexprtrees != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprtrees),   oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray(&((oracle->consexprtrees)[oracle->nconss]), nconss);
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprvaridxs), oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray(&((oracle->consexprvaridxs)[oracle->nconss]), nconss);
   }

   if( consnames != NULL )
   {
      if( oracle->consnames == NULL && oracle->nconss != 0 )
      {
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consnames), oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->consnames, oracle->nconss);
      }
      else
      {
         if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consnames), oracle->nconss, oracle->nconss + nconss) == NULL )
            return SCIP_NOMEMORY;
      }
      for( i = 0; i < nconss; ++i )
         if( consnames[i] != NULL )
         {
            if( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->consnames)[oracle->nconss + i]), consnames[i], strlen(consnames[i])+1) == NULL )
               return SCIP_NOMEMORY;
         }
         else
            oracle->consnames[oracle->nconss + i] = NULL;
   }
   else if( oracle->consnames != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consnames), oracle->nconss, oracle->nconss + nconss) == NULL )
         return SCIP_NOMEMORY;
      BMSclearMemoryArray(&((oracle->consnames)[oracle->nconss]), nconss);
   }
   
   oracle->nconss += nconss;

   if( addednlcon == TRUE )
      invalidateHessianLagSparsity(oracle);
   
   return SCIP_OKAY;
}

/** sets or overwrites objective, a minization problem is expected
 * 
 *  May change sparsity pattern.
 */
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real       constant,           /**< constant part of objective */
   int                   nlin,               /**< number of linear variable coefficients */ 
   const int*            lininds,            /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables, or NULL if no linear part */
   int                   nquadrows,          /**< number of columns in the matrix of the quadratic part */
   const int*            quadrowidxs,        /**< indices of variables appearing in quadratic part, or NULL if no quadratic part */
   const int*            quadoffsets,        /**< start index of each rows quadratic coefficients in quadidxs and quadvals 
                                              *   (quadoffsets[.][nvars] gives length of quadidxs and quadvals), or NULL if no quadratic part */
   const int*            quadidxs,           /**< column indices in quadratic part, or NULL if no quadratic part */ 
   const SCIP_Real*      quadvals,           /**< coefficients in quadratic part, or NULL if no quadratic part */
   const int*            exprvaridxs,        /**< maps variable indices in expression tree to indices in nlp, or NULL if no nonquadratic part */
   const SCIP_EXPRTREE*  exprtree            /**< expression tree of nonquadratic part, or NULL if no nonquadratic part */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   
   /* clear previous objective */
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objlinidxs,  oracle->objnlin);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objlinvals,  oracle->objnlin);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objquadrows, oracle->objquadlen);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objquadcols, oracle->objquadlen);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objquadvals, oracle->objquadlen);
   if( oracle->objexprtree != NULL )
   {
#ifdef WITH_NL
      BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->objexprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree));
      SCIP_CALL( SCIPexprtreeFree(&(*oracle)->objexprtree) );
#endif
      /* @TODO this does not clear the vardegrees's */
   }
   
   if( nquadrows != 0 || oracle->objquadlen != 0 || exprtree != NULL || oracle->objexprtree != NULL )
      invalidateHessianLagSparsity(oracle);
   
   oracle->objconstant = constant;
   
   if( nlin != 0 )
   {
      int i;

      assert(lininds != NULL);
      assert(linvals != NULL);

      if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objlinidxs, lininds, nlin) == NULL )
         return SCIP_NOMEMORY;
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objlinvals, linvals, nlin) == NULL )
         return SCIP_NOMEMORY;
      oracle->objnlin = nlin;
      for( i = 0; i < nlin; ++i )
         oracle->vardegrees[lininds[i]] = MAX(1, oracle->vardegrees[lininds[i]]);
      
      SCIPsortIntReal(oracle->objlinidxs, oracle->objlinvals, nlin);
   }
   else
      oracle->objnlin = 0;
   
   if( nquadrows != 0 && quadoffsets != NULL )
   {
      int j;
      int k;

      assert(quadrowidxs != NULL);
      assert(quadidxs != NULL);
      assert(quadvals != NULL);
      
      oracle->objquadlen = quadoffsets[nquadrows];
      if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->objquadrows, oracle->objquadlen) == NULL )
         return SCIP_NOMEMORY;
      if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->objquadcols, oracle->objquadlen) == NULL )
         return SCIP_NOMEMORY;
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objquadvals, quadvals, oracle->objquadlen) == NULL )
         return SCIP_NOMEMORY;
      
      k = quadoffsets[0];
      for( j = 0; j < nquadrows; ++j )
      {
         if( quadoffsets[j] == quadoffsets[j+1] )
            continue;
         for( ; k < quadoffsets[j+1]; ++k )
         {
            assert(quadidxs[k] < nquadrows);
            if( quadrowidxs[j] > quadrowidxs[quadidxs[k]] )
            {
               oracle->objquadrows[k] = quadrowidxs[j];
               oracle->objquadcols[k] = quadrowidxs[quadidxs[k]];
            }
            else
            {
               oracle->objquadrows[k] = quadrowidxs[quadidxs[k]];
               oracle->objquadcols[k] = quadrowidxs[j];
            }
            oracle->vardegrees[quadrowidxs[quadidxs[k]]] = MAX(2, oracle->vardegrees[quadrowidxs[quadidxs[k]]]);
         }
         oracle->vardegrees[quadrowidxs[j]] = MAX(2, oracle->vardegrees[quadrowidxs[j]]);
      }
      assert(k == oracle->objquadlen);

      /* sort quadrows and quadcols */
      SCIPsortIntIntReal(oracle->objquadrows, oracle->objquadcols, oracle->objquadvals, oracle->objquadlen);
      j = 0;
      k = 0;
      while( j < oracle->objquadlen )
      {
         while( k < oracle->objquadlen && oracle->objquadrows[j] == oracle->objquadrows[k] )
            ++k;
         SCIPsortIntReal(&oracle->objquadcols[j], &oracle->objquadvals[j], k-j);
         j = k;
      }
   }
   else
      oracle->objquadlen = 0;
   
   if( exprtree != NULL )
   {
#ifdef WITH_NL
      int j;
      
      assert(oracle->exprinterpreter != NULL);
      assert(SCIPexprtreeHasVarsAsIndex((SCIP_EXPRTREE*)exprtree));
      
      SCIP_CALL( SCIPexprtreeCopy(&oracle->objexprtree, (SCIP_EXPRTREE*)exprtree) );
      
      SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->objexprtree) );
      
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objexprvaridxs, exprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree)) == NULL )
         return SCIP_NOMEMORY;
      for (j = 0; j < SCIPexprtreeGetNVars(oracle->objexprtree); ++j)
      {
         assert(exprvaridxs[j] >= 0);
         assert(exprvaridxs[j] <  oracle->nvars);
         oracle->vardegrees[exprvaridxs[j]] = INT_MAX; /* @TODO could try to be more clever, maybe use getMaxDegree function in exprtree */
      }
#else
      SCIPerrorMessage("nonquadratic functions not supported in NLPI yet\n");
      return SCIP_ERROR;
#endif
   }
   
   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ubs                 /**< new upper bounds, or NULL if all should be +infty */
   )
{
   int i;

   assert(oracle != NULL);
   assert(indices != NULL || nvars == 0);
   
   for( i = 0; i < nvars; ++i )
   {
      assert(indices != NULL);
      assert(indices[i] >= 0);
      assert(indices[i] < oracle->nvars);
      
      oracle->varlbs[indices[i]] = (lbs != NULL ? lbs[i] : -oracle->infinity);
      oracle->varubs[indices[i]] = (ubs != NULL ? ubs[i] :  oracle->infinity);
      
      if( oracle->varlbs[indices[i]] > oracle->varubs[indices[i]] )
      { /* inconsistent bounds; let's assume it's due to rounding and make them equal */
         assert(EPSLE(oracle->varlbs[indices[i]], oracle->varubs[indices[i]], SCIP_DEFAULT_EPSILON) == TRUE);
         oracle->varlbs[indices[i]] = oracle->varubs[indices[i]];
      }
   }
   
   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_RETCODE SCIPnlpiOracleChgConsBounds(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             nconss,             /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lhss,               /**< new left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhss                /**< new right-hand sides, or NULL if all should be +infty */
   )
{
   int i;
   
   assert(oracle != NULL);
   assert(indices != NULL || nconss == 0);
   
   for( i = 0; i < nconss; ++i )
   {
      assert(indices != NULL);
      assert(indices[i] >= 0);
      assert(indices[i] < oracle->nconss);
      
      oracle->conslhss[indices[i]] = (lhss != NULL ? lhss[i] : -oracle->infinity);
      oracle->consrhss[indices[i]] = (rhss != NULL ? rhss[i] :  oracle->infinity);
   }
   
   return SCIP_OKAY;
}

/** deletes a set of variables */
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< deletion status of vars in input (1 if var should be deleted, 0 if not); 
                                              *   new position of var in output (-1 if var was deleted) */
   )
{  /*lint --e{715}*/
   int c, j, k;
   int lastgood; /* index of the last variable that should be kept */ 
   
   assert(oracle != NULL);

   invalidateJacobiSparsity(oracle);
   invalidateHessianLagSparsity(oracle);
   
   lastgood = oracle->nvars - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1)
      --lastgood;
   if( lastgood < 0 )
   { /* all variables should be deleted */
      assert(oracle->nconss == 0); /* we could relax this by checking that all constraints are constant */
      assert(oracle->objquadlen == 0);
#ifdef WITH_NL
      assert(oracle->objexprtree == NULL || SCIPexprtreeGetNVars(oracle->objexprtree) == 0);
#endif
      if( oracle->objnlin > 0 )
      {
         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->objlinidxs, oracle->objnlin);
         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->objlinvals, oracle->objnlin);
         oracle->objnlin = 0;
      }
      for( c = 0; c < oracle->nvars; ++c )
         delstats[c] = -1;
      SCIP_CALL( freeVariables(oracle) );
      return SCIP_OKAY;
   }
   
   /* delete variables at the end */
   for( c = oracle->nvars - 1; c > lastgood; --c )
   {
      SCIP_CALL( freeVariable(oracle, c) );
      delstats[c] = -1;
   }

   /* go through variables from the beginning on
    * if variable should be deleted, free it and move lastgood variable to this position
    * then update lastgood */
   for( c = 0; c <= lastgood; ++c )
   {
      if( delstats[c] == 0 )
      { /* variable should not be deleted and is kept on position c */
         delstats[c] = c;
         continue;
      }
      assert(delstats[c] == 1); /* variable should be deleted */
      
      SCIP_CALL( freeVariable(oracle, c) );
      delstats[c] = -1;
      
      /* move constraint at position lastgood to position c */
      SCIP_CALL( moveVariable(oracle, lastgood, c) );
      delstats[lastgood] = c; /* mark that lastgood variable is now at position c */
      
      /* move lastgood forward, delete variables on the way */
      --lastgood;
      while( lastgood > c && delstats[lastgood] == 1)
      {
         SCIP_CALL( freeVariable(oracle, lastgood) );
         delstats[lastgood] = -1;
         --lastgood;
      }
   }
   assert(c == lastgood);

   mapIndices(delstats, oracle->objnlin, oracle->objlinidxs); /* TODO delete entries for deleted variables */
   SCIPsortIntReal(oracle->objlinidxs, oracle->objlinvals, oracle->objnlin);
   
   mapIndices2(delstats, oracle->objquadlen, oracle->objquadrows, oracle->objquadcols);
   /* sort quadrows and quadcols */
   SCIPsortIntIntReal(oracle->objquadrows, oracle->objquadcols, oracle->objquadvals, oracle->objquadlen);
   j = 0;
   k = 0;
   while( j < oracle->objquadlen )
   {
      while( k < oracle->objquadlen && oracle->objquadrows[j] == oracle->objquadrows[k] )
         ++k;
      SCIPsortIntReal(&oracle->objquadcols[j], &oracle->objquadvals[j], k-j);
      j = k;
   }

#ifdef WITH_NL
   if( oracle->objexprtree )
      mapIndices(delstats, SCIPexprtreeGetNVars(oracle->objexprtree), oracle->objexprvaridxs);
#endif
   
   for( c = 0; c < oracle->nconss; ++c )
   {
      if( oracle->conslinlens )
      {
         mapIndices(delstats, oracle->conslinlens[c], oracle->conslinidxs[c]);
         SCIPsortIntReal(oracle->conslinidxs[c], oracle->conslincoefs[c], oracle->conslinlens[c]);
      }
      
      if( oracle->consquadlens )
      {
         mapIndices2(delstats, oracle->consquadlens[c], oracle->consquadrows[c], oracle->consquadcols[c]);
         /* sort quadrows and quadcols */
         SCIPsortIntIntReal(oracle->consquadrows[c], oracle->consquadcols[c], oracle->consquadvals[c], oracle->consquadlens[c]);
         j = 0;
         k = 0;
         while( j < oracle->consquadlens[c] )
         {
            while( k < oracle->consquadlens[c] && oracle->consquadrows[c][j] == oracle->consquadrows[c][k] )
               ++k;
            SCIPsortIntReal(&oracle->consquadcols[c][j], &oracle->consquadvals[c][j], k-j);
            j = k;
         }
      }
      
#ifdef WITH_NL
      if( oracle->consexprtrees )
         mapIndices(delstats, SCIPexprtreeGetNVars(oracle->consexprtrees[c]), oracle->consexprvaridxs[c]);
#endif
   }
   
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varlbs, oracle->nvars, lastgood+1) == NULL )
      return SCIP_NOMEMORY;
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varubs, oracle->nvars, lastgood+1) == NULL )
      return SCIP_NOMEMORY;
   
   if( oracle->varnames != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varnames, oracle->nvars, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
   }

   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->vardegrees, oracle->nvars, lastgood+1) == NULL )
      return SCIP_NOMEMORY;

   oracle->nvars = lastgood+1;

   return SCIP_OKAY;
}

/** deletes a set of constraints */
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not); 
                                              *   new position of row in output (-1 if row was deleted) */
   )
{  /*lint --e{715}*/
   int c;
   int lastgood; /* index of the last constraint that should be kept */ 
   
   assert(oracle != NULL);

   invalidateJacobiSparsity(oracle);
   invalidateHessianLagSparsity(oracle);
   
   lastgood = oracle->nconss - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1)
      --lastgood;
   if( lastgood < 0 )
   { /* all constraints should be deleted */
      for( c = 0; c < oracle->nconss; ++c )
         delstats[c] = -1;
      SCIP_CALL( freeConstraints(oracle) );
      return SCIP_OKAY;
   }
   
   /* delete constraints at the end */
   for( c = oracle->nconss - 1; c > lastgood; --c )
   {
      SCIP_CALL( freeConstraint(oracle, c) );
      delstats[c] = -1;
   }

   /* go through constraint from the beginning on
    * if constraint should be deleted, free it and move lastgood constraint to this position
    * then update lastgood */
   for( c = 0; c <= lastgood; ++c )
   {
      if( delstats[c] == 0 )
      { /* constraint should not be deleted and is kept on position c */
         delstats[c] = c;
         continue;
      }
      assert(delstats[c] == 1); /* constraint should be deleted */
      
      SCIP_CALL( freeConstraint(oracle, c) );
      delstats[c] = -1;
      
      /* move constraint at position lastgood to position c */
      SCIP_CALL( moveConstraint(oracle, lastgood, c) );
      delstats[lastgood] = c; /* mark that lastgood constraint is now at position c */
      
      /* move lastgood forward, delete constraints on the way */
      --lastgood;
      while( lastgood > c && delstats[lastgood] == 1)
      {
         SCIP_CALL( freeConstraint(oracle, lastgood) );
         delstats[lastgood] = -1;
         --lastgood;
      }
   }
   assert(c == lastgood);
   
   /* free memory */
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslhss, oracle->nconss, lastgood+1) == NULL )
      return SCIP_NOMEMORY;
   if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consrhss, oracle->nconss, lastgood+1) == NULL )
      return SCIP_NOMEMORY;

   if( oracle->conslinlens != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
   }

   if( oracle->consquadlens != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadlens, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadrows, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadcols, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadvals, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
   }
   
   if( oracle->consexprvaridxs != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consexprtrees,   oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
   }
   
   if( oracle->consnames != NULL )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consnames, oracle->nconss, lastgood+1) == NULL )
         return SCIP_NOMEMORY;
   }

   oracle->nconss = lastgood+1;

   return SCIP_OKAY;
}

/** changes linear coefficients in one constraint or objective
 */
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidxs,            /**< array with indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoefs            /**< array with new coefficients of variables */
   )
{  /*lint --e{715}*/
   SCIP_Bool addednew;
   int       i;
   
   int**       linidxs;
   SCIP_Real** lincoefs;
   int*        linlen;
   
   assert(oracle != NULL);
   assert(varidxs != NULL || nentries == 0);
   assert(newcoefs != NULL || nentries == 0);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   
   if( nentries == 0 )
      return SCIP_OKAY;
   
   addednew = FALSE;
   
   if( considx == -1 )
   {
      linidxs  = &oracle->objlinidxs;
      lincoefs = &oracle->objlinvals;
      linlen   = &oracle->objnlin;
   }
   else
   {
      if( oracle->conslinlens == NULL )
      { /* first time we have linear coefficients in a constraint */
         if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->conslinlens,  oracle->nconss);
         BMSclearMemoryArray(oracle->conslinidxs,  oracle->nconss);
         BMSclearMemoryArray(oracle->conslincoefs, oracle->nconss);
      }
      
      linidxs  = &oracle->conslinidxs[considx];
      lincoefs = &oracle->conslincoefs[considx];
      linlen   = &oracle->conslinlens[considx];
   }
   
   if( *linlen == 0 )
   { /* first time we have linear coefficients in this constraint (or objective) */
      assert(*linidxs == NULL);
      assert(*lincoefs == NULL);
      
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, linidxs,  varidxs,  nentries) == NULL )
         return SCIP_NOMEMORY;
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, lincoefs, newcoefs, nentries) == NULL )
         return SCIP_NOMEMORY;
      SCIPsortIntReal(*linidxs, *lincoefs, nentries);
            
      addednew = TRUE;
   }
   else
   {
      int pos;
      int len = *linlen;
      
      for( i = 0; i < nentries; ++i )
      {
         /* we assume that indices are not repeating in varidxs (!) */
         if( SCIPsortedvecFindInt(*linidxs, varidxs[i], *linlen, &pos) )
         {
            (*lincoefs)[pos] = newcoefs[i];
         }
         else
         {
            if( !addednew )
            { /* first coefficient that is added new, realloc memory */
               int newsize = *linlen + (nentries-i); /* new size for arrays (upper bound) */
               if( BMSreallocBlockMemoryArray(oracle->blkmem, linidxs,  *linlen, newsize) == NULL )
                  return SCIP_NOMEMORY;
               if( BMSreallocBlockMemoryArray(oracle->blkmem, lincoefs, *linlen, newsize) == NULL )
                  return SCIP_NOMEMORY;
            }
            /* append new entries */
            (*linidxs)[len]  = varidxs[i];
            (*lincoefs)[len] = newcoefs[i];
            ++len;
            addednew = TRUE;
            /* increase degree of variable to 1 */
            oracle->vardegrees[varidxs[i]] = MAX(1, oracle->vardegrees[varidxs[i]]);
         }
      }
      
      if( addednew )
      {
         /* shrink to actual needed size */
         if( BMSreallocBlockMemoryArray(oracle->blkmem, linidxs , *linlen + (nentries-i), len) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, lincoefs, *linlen + (nentries-i), len) == NULL )
            return SCIP_NOMEMORY;
         *linlen = len;
         
         SCIPsortIntReal(*linidxs, *lincoefs, len);
      }
   }
   
   if( addednew )
      invalidateJacobiSparsity(oracle);
   
   assert(*linlen == 0 || (*linidxs)[0] >= 0);
   assert(*linlen == 0 || (*linidxs)[*linlen-1] < oracle->nvars);
   
   return SCIP_OKAY;
}

/** changes (or adds) coefficients in the quadratic part of one constraint or objective
 */
SCIP_RETCODE SCIPnlpiOracleChgQuadCoefs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where quadratic coefficients should be changed, or -1 for objective */
   const int             nentries,           /**< number of coefficients to change */
   const int*            rowidxs,            /**< array with row indices of quadratic matrix entries for which new values are provided */
   const int*            colidxs,            /**< array with column indices of quadratic matrix entries for which new values are provided */
   SCIP_Real*            newcoefs            /**< new quadratic coefficients */ 
   )
{  /*lint --e{715}*/
   SCIP_Bool addednew;
   int       i;
   
   int**       quadrowidxs;
   int**       quadcolidxs;
   SCIP_Real** quadcoefs;
   int*        quadlen;
   
   assert(oracle != NULL);
   assert(rowidxs != NULL || nentries == 0);
   assert(colidxs != NULL || nentries == 0);
   assert(newcoefs != NULL || nentries == 0);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   
   if( nentries == 0 )
      return SCIP_OKAY;
   
   addednew = FALSE;
   
   if( considx == -1 )
   {
      quadrowidxs = &oracle->objquadrows;
      quadcolidxs = &oracle->objquadcols;
      quadcoefs   = &oracle->objquadvals;
      quadlen     = &oracle->objquadlen;
   }
   else
   {
      if( oracle->consquadlens == NULL )
      { /* first time we have quadratic coefficients in a constraint */
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens), oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadrows), oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadcols), oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         if( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadvals), oracle->nconss) == NULL )
            return SCIP_NOMEMORY;
         BMSclearMemoryArray(oracle->consquadlens, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadrows, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadcols, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadvals, oracle->nconss);
      }
      
      quadrowidxs = &oracle->consquadrows[considx];
      quadcolidxs = &oracle->consquadcols[considx];
      quadcoefs   = &oracle->consquadvals[considx];
      quadlen     = &oracle->consquadlens[considx];
   }
   
   if( *quadlen == 0 )
   { /* first time we have quadratic coefficients in this constraint (or objective) */
      assert(*quadrowidxs == NULL);
      assert(*quadcolidxs == NULL);
      assert(*quadcoefs   == NULL);
      
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, quadrowidxs, rowidxs,  nentries) == NULL )
         return SCIP_NOMEMORY;
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, quadcolidxs, colidxs,  nentries) == NULL )
         return SCIP_NOMEMORY;
      if( BMSduplicateBlockMemoryArray(oracle->blkmem, quadcoefs,   newcoefs, nentries) == NULL )
         return SCIP_NOMEMORY;
            
      addednew = TRUE;
   }
   else
   {
      int rowpos;
      int colpos;
      int len = *quadlen;
      
      for( i = 0; i < nentries; ++i )
      {
         assert(rowidxs[i] >= 0);
         assert(colidxs[i] >= 0);
         assert(rowidxs[i] < oracle->nvars);
         assert(colidxs[i] < oracle->nvars);
         
         /* we assume that index pairs are not repeating */
         /* see if we already have an entry for (rowidxs[i], colidxs[i]) */
         if( SCIPsortedvecFindInt(*quadrowidxs, rowidxs[i], *quadlen, &rowpos) )
         {
            int nextrowpos;
            if( !SCIPsortedvecFindInt(&((*quadrowidxs)[rowpos+1]), rowidxs[i]+1, *quadlen-rowpos, &nextrowpos) )
               nextrowpos = *quadlen;
            if( !SCIPsortedvecFindInt(&((*quadcolidxs)[rowpos]), colidxs[i], nextrowpos-rowpos, &colpos) )
               colpos = -1;
         }
         else
            colpos = -1;
         
         if( colpos >= 0 )
         {
            (*quadcoefs)[colpos] = newcoefs[i];
         }
         else
         {
            if( !addednew )
            { /* first coefficient that is added new, realloc memory */
               int newsize = *quadlen + (nentries-i); /* new size for arrays (upper bound) */
               if( BMSreallocBlockMemoryArray(oracle->blkmem, quadrowidxs, *quadlen, newsize) == NULL )
                  return SCIP_NOMEMORY;
               if( BMSreallocBlockMemoryArray(oracle->blkmem, quadcolidxs, *quadlen, newsize) == NULL )
                  return SCIP_NOMEMORY;
               if( BMSreallocBlockMemoryArray(oracle->blkmem, quadcoefs,   *quadlen, newsize) == NULL )
                  return SCIP_NOMEMORY;
            }
            /* append new entries */
            (*quadrowidxs)[len] = rowidxs[i];
            (*quadcolidxs)[len] = colidxs[i];
            (*quadcoefs)  [len] = newcoefs[i];
            ++len;
            addednew = TRUE;
            /* increase degree of variable to 1 */
            oracle->vardegrees[rowidxs[i]] = MAX(2, oracle->vardegrees[rowidxs[i]]);
            oracle->vardegrees[colidxs[i]] = MAX(2, oracle->vardegrees[colidxs[i]]);
         }
      }
      
      if( addednew )
      {
         /* shrink to actual needed size */
         if( BMSreallocBlockMemoryArray(oracle->blkmem, quadrowidxs, *quadlen + (nentries-i), len) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, quadcolidxs, *quadlen + (nentries-i), len) == NULL )
            return SCIP_NOMEMORY;
         if( BMSreallocBlockMemoryArray(oracle->blkmem, quadcoefs,   *quadlen + (nentries-i), len) == NULL )
            return SCIP_NOMEMORY;
         *quadlen = len;
      }
   }
   
   if( addednew )
   {
      int j, k;
      
      invalidateJacobiSparsity(oracle);
      
      /* sort quadrows and quadcols */
      SCIPsortIntIntReal(*quadrowidxs, *quadcolidxs, *quadcoefs, *quadlen);
      j = 0;
      k = 0;
      while( j < *quadlen )
      {
         while( k < *quadlen && (*quadrowidxs)[j] == (*quadrowidxs)[k] )
            ++k;
         SCIPsortIntReal(&((*quadcolidxs)[j]), &((*quadcoefs)[j]), k-j);
         j = k;
      }
   }
   
   assert(*quadlen == 0 || (*quadrowidxs)[0] >= 0);
   assert(*quadlen == 0 || (*quadrowidxs)[*quadlen-1] < oracle->nvars);
   assert(*quadlen == 0 || (*quadcolidxs)[0] >= 0);
   assert(*quadlen == 0 || (*quadcolidxs)[*quadlen-1] < oracle->nvars);
   
   return SCIP_OKAY;
}

/** gives the current number of variables */
int SCIPnlpiOracleGetNVars(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->nvars;
}

/** gives the current number of constraints */
int SCIPnlpiOracleGetNConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->nconss;
}

/** gives the variables lower bounds */
const SCIP_Real* SCIPnlpiOracleGetVarLbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->varlbs;
}

/** gives the variables upper bounds */
const SCIP_Real* SCIPnlpiOracleGetVarUbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->varubs;
}

/** Gives maximum degree of a variable w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
int SCIPnlpiOracleGetVarDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx
   )
{
   assert(oracle != NULL);
   assert(varidx >= 0);
   assert(varidx < oracle->nvars);
   
   return oracle->vardegrees[varidx];
}

/** Gives maximum degree of all variables w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
int* SCIPnlpiOracleGetVarDegrees(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->vardegrees;
}

/** gives the constraints left-hand sides */
const SCIP_Real* SCIPnlpiOracleGetConstraintLhss(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->conslhss;
}

/** gives the constraints right-hand sides */
const SCIP_Real* SCIPnlpiOracleGetConstraintRhss(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->consrhss;
}

/** Gives maximum degree of a constraints.
 *  The degree of a constraint is the maximal degree of all summands which appear in it, and is infinity for nonpolynomial terms.
 */ 
int SCIPnlpiOracleGetConstraintDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< index of constraint for which the degree is requested */
   )
{
   assert(oracle != NULL);
   assert(considx >= 0);
   assert(considx < oracle->nconss);
   
   if( oracle->consexprtrees && oracle->consexprtrees[considx] )
      return INT_MAX;  /* @TODO something more clever, e.g., use exprtreeGetMaxDegree*/
   
   if( oracle->consquadlens && oracle->consquadlens[considx] )
      return 2;

   if( oracle->conslinlens && oracle->conslinlens[considx] )
      return 1;

   return 0;
}

/** evaluates the objective function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */  
   )
{
   assert(oracle != NULL);
   
   SCIP_CALL( evalFunctionValue(oracle,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals,
      oracle->objexprvaridxs, oracle->objexprtree,
      x, objval) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** evaluates one constraint function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */  
   )
{
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( evalFunctionValue(oracle,
      (oracle->conslinlens   != NULL ? oracle->conslinlens [considx] : 0),
      (oracle->conslinlens   != NULL ? oracle->conslinidxs [considx] : NULL),
      (oracle->conslinlens   != NULL ? oracle->conslincoefs[considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadlens[considx] : 0),
      (oracle->consquadlens  != NULL ? oracle->consquadrows[considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadcols[considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadvals[considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprtrees[considx] : NULL),
      x, conval) );
   
   return SCIP_OKAY;
}

/** evaluates all constraint functions in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< buffer to store constraint values */  
   )
{
   int i;
   
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(convals != NULL);

   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_CALL( SCIPnlpiOracleEvalConstraintValue(oracle, i, x, &convals[i]) );
   }
   
   return SCIP_OKAY;
}

/** computes the objective gradient in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Real*            objgrad             /**< pointer to store (dense) objective gradient */  
   )
{
   assert(oracle != NULL);
   
   SCIP_CALL( evalFunctionGradient(oracle,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals,
      oracle->objexprvaridxs, oracle->objexprtree,
      x, newx, objval, objgrad) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** computes a constraints gradient in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */ 
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */  
   )
{
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( evalFunctionGradient(oracle,
      (oracle->conslinlens   != NULL ? oracle->conslinlens [considx] : 0),
      (oracle->conslinlens   != NULL ? oracle->conslinidxs [considx] : NULL),
      (oracle->conslinlens   != NULL ? oracle->conslincoefs[considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadlens[considx] : 0),
      (oracle->consquadlens  != NULL ? oracle->consquadrows[considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadcols[considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadvals[considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprtrees[considx] : NULL),
      x, newx, conval, congrad) );
   
   return SCIP_OKAY;
}

/** Gets sparsity pattern (rowwise) of Jacobian matrix.
 * 
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting constraints destroys the sparsity structure and make another call to this function necessary. 
 */
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[nconss] gives length of col, can be NULL */
   )
{
   SCIP_Bool* nzflag;
   int nnz;
   int maxnnz;
   int i;
   int j;
   
   assert(oracle != NULL);
   
   if( oracle->jacoffsets != NULL )
   {
      assert(oracle->jaccols != NULL);
      if( offset != NULL )
         *offset = oracle->jacoffsets;
      if( col != NULL )
         *col = oracle->jaccols;
      return SCIP_OKAY;
   }

   if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->jacoffsets, oracle->nconss + 1) == NULL )
      return SCIP_NOMEMORY;

   maxnnz = MIN(oracle->nvars, 10) * oracle->nconss;  /* initial guess */
   if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, maxnnz) == NULL )
      return SCIP_NOMEMORY;
   
   if( maxnnz == 0 )
   {  /* no variables */
      BMSclearMemoryArray(oracle->jacoffsets, oracle->nconss + 1);
      if( offset != NULL )
         *offset = oracle->jacoffsets;
      if( col != NULL )
         *col = oracle->jaccols;
      return SCIP_OKAY;
   }
   nnz = 0;
   
   if( BMSallocBlockMemoryArray(oracle->blkmem, &nzflag, oracle->nvars) == NULL )
      return SCIP_NOMEMORY;

   for( i = 0; i < oracle->nconss; ++i )
   {
      oracle->jacoffsets[i] = nnz;
      
      if( oracle->conslinlens != NULL && (oracle->consquadlens == NULL || oracle->consquadlens[i] == 0) && (oracle->consexprtrees == NULL || oracle->consexprtrees[i] == NULL) )
      { /* linear constraint: since we just want to copy the conslinvals at EvalJacobian, we need to copy conslinidxs here too */
         int nz = oracle->conslinlens[i];
         if( nz > 0 )
         {
            if( nnz + nz > maxnnz )
            {
               int oldsize = maxnnz;
               maxnnz = MAX(2*maxnnz, nnz+nz);
               if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, oldsize, maxnnz) == NULL )
                  return SCIP_OKAY;
            }
            assert(maxnnz >= nnz + nz);
            BMScopyMemoryArray(&oracle->jaccols[nnz], oracle->conslinidxs[i], nz);
            nnz += nz;
         }
         continue;
      }
      
      /* check which variables appear in constraint i */
      BMSclearMemoryArray(nzflag, oracle->nvars);
      
      if( oracle->conslinlens != NULL )
      {
         assert(oracle->conslinidxs);
         for( j = 0; j < oracle->conslinlens[i]; ++j )
            nzflag[oracle->conslinidxs[i][j]] = TRUE;
      }
      
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
      {
         assert(oracle->consquadrows    != NULL );
         assert(oracle->consquadrows[i] != NULL );
         assert(oracle->consquadcols    != NULL );
         assert(oracle->consquadcols[i] != NULL );
         for( j = 0; j < oracle->consquadlens[i]; ++j )
         {
            nzflag[oracle->consquadrows[i][j]] = TRUE;
            nzflag[oracle->consquadcols[i][j]] = TRUE;
         }
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         assert(oracle->consexprvaridxs != NULL );
         assert(oracle->consexprvaridxs[i] != NULL );
#ifdef WITH_NL
         for (j = 0; j < SCIPexprtreeGetNVars(oracle->consexprtrees[i]); ++j)
            nzflag[oracle->consexprvaridxs[i][j]] = TRUE;
#else
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
#endif
      }
      
      /* store variables indices in jaccols */
      for( j = 0; j < oracle->nvars; ++j )
      {
         if( nzflag[j] == FALSE )
            continue;
         
         if( nnz >= maxnnz )
         {
            int oldsize = maxnnz;
            maxnnz = MAX(nnz, 2*maxnnz);
            if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, oldsize, maxnnz) == NULL )
               return SCIP_NOMEMORY;
         }
         
         oracle->jaccols[nnz] = j;
         ++nnz;
      }
   }
   
   oracle->jacoffsets[oracle->nconss] = nnz;
   
   /* shrink jaccols array to nnz */
   if( nnz < maxnnz )
   {
      if( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, maxnnz, nnz) == NULL )
         return SCIP_NOMEMORY;
   }
   
   BMSfreeBlockMemoryArray(oracle->blkmem, &nzflag, oracle->nvars);

   if( offset != NULL )
      *offset = oracle->jacoffsets;
   if( col != NULL )
      *col = oracle->jaccols;

   return SCIP_OKAY;
}

/** Evaluates the Jacobi matrix in a given point.
 * 
 *  The values in the Jacobi matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity.
 *  The user need to call SCIPnlpiOracleGetJacobianSparsity at least ones before using this function. 
 */
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real*            convals,            /**< pointer to store constraint values, can be NULL */ 
   SCIP_Real*            jacobi              /**< pointer to store sparse jacobian values */  
   )
{
   SCIP_Real* grad;
   SCIP_Real dummy;
   int i;
   int j;
   int k;
   int l;
   
   assert(oracle != NULL);
   assert(jacobi != NULL);
   
   assert(oracle->jacoffsets != NULL);
   assert(oracle->jaccols != NULL);
   
   if( BMSallocBlockMemoryArray(oracle->blkmem, &grad, oracle->nvars) == NULL )
      return SCIP_NOMEMORY;
   
   j = oracle->jacoffsets[0];
   k = 0;
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->conslinlens != NULL && (oracle->consquadlens == NULL || oracle->consquadlens[i] == 0) && (oracle->consexprtrees == NULL || oracle->consexprtrees[i] == NULL) )
      { /* linear constraints should be easy */
         l = oracle->jacoffsets[i+1] - j;
         assert(l == oracle->conslinlens[i]);
         if( l > 0 )
         {
            BMScopyMemoryArray(&jacobi[k], oracle->conslincoefs[i], l);
            j += l; 
            k += l;
         }
      }
      else
      { /* @TODO do this sparse too */
         SCIP_RETCODE retcode = SCIPnlpiOracleEvalConstraintGradient(oracle, i, x, newx, (convals ? &convals[i] : &dummy), grad);
         if( retcode != SCIP_OKAY )
         {
            BMSfreeBlockMemoryArray(oracle->blkmem, &grad, oracle->nvars);
            return retcode;
         }
      
         for( ; j < oracle->jacoffsets[i+1]; ++j, ++k )
            jacobi[k] = grad[oracle->jaccols[j]];
      }
   }
   
   BMSfreeBlockMemoryArray(oracle->blkmem, &grad, oracle->nvars);
   
   return SCIP_OKAY;
}

/** Gets sparsity pattern of the Hessian matrix of the Lagrangian.
 * 
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 *  Only elements of the lower left triangle and the diagonal are counted.
 */
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[nconss] gives length of col, can be NULL */
   )
{
   int** colnz;   /** nonzeros in Hessian corresponding to one column */
   int*  collen;  /** collen[i] is length of array colnz[i] */
   int*  colnnz;  /** colnnz[i] is number of entries in colnz[i] (<= collen[i]) */ 
   int nnz;
   int i;
   int j;
   int cnt;
   
   assert(oracle != NULL);
   
   if( oracle->heslagoffsets != NULL )
   {
      assert(oracle->heslagcols != NULL);
      if( offset != NULL )
         *offset = oracle->heslagoffsets;
      if( col != NULL )
         *col = oracle->heslagcols;
      return SCIP_OKAY;
   }

   if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->heslagoffsets, oracle->nvars + 1) == NULL )
      return SCIP_NOMEMORY;
   nnz = 0;
   
   if( BMSallocBlockMemoryArray(oracle->blkmem, &colnz,  oracle->nvars) == NULL )
      return SCIP_NOMEMORY;
   if( BMSallocBlockMemoryArray(oracle->blkmem, &collen, oracle->nvars) == NULL )
      return SCIP_NOMEMORY;
   if( BMSallocBlockMemoryArray(oracle->blkmem, &colnnz, oracle->nvars) == NULL )
      return SCIP_NOMEMORY;
   BMSclearMemoryArray(colnz,  oracle->nvars);
   BMSclearMemoryArray(collen, oracle->nvars);
   BMSclearMemoryArray(colnnz, oracle->nvars);
   nnz = 0;
   
   if( oracle->objquadlen != 0 )
   {
      SCIP_CALL( hessLagSparsitySetNzFlagForQuad(oracle, colnz, collen, colnnz, &nnz, oracle->objquadrows, oracle->objquadcols, oracle->objquadlen) );
   }

   if( oracle->objexprtree != NULL )
   {
#ifdef WITH_NL
      SCIP_CALL( hessLagSparsitySetNzFlagForExprtree(oracle, colnz, collen, colnnz, &nnz, oracle->objexprvaridxs, oracle->objexprtree, oracle->nvars) );
#else
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
#endif
   }

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
      {
         SCIP_CALL( hessLagSparsitySetNzFlagForQuad(oracle, colnz, collen, colnnz, &nnz, oracle->consquadrows[i], oracle->consquadcols[i], oracle->consquadlens[i]) );
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
#ifdef WITH_NL
         SCIP_CALL( hessLagSparsitySetNzFlagForExprtree(oracle, colnz, collen, colnnz, &nnz, oracle->consexprvaridxs[i], oracle->consexprtrees[i], oracle->nvars) );
#else
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
#endif
      }
   }
   
   if( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->heslagcols, nnz) == NULL )
      return SCIP_NOMEMORY;
   
   /* set hessian sparsity from colnz, colnnz */
   cnt = 0;
   for( i = 0; i < oracle->nvars; ++i )
   {
      oracle->heslagoffsets[i] = cnt;
      for( j = 0; j < colnnz[i]; ++j )
      {
         assert(cnt < nnz);
         oracle->heslagcols[cnt++] = colnz[i][j];
      }
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &colnz[i], collen[i]);
      collen[i] = 0;
   }
   oracle->heslagoffsets[oracle->nvars] = cnt;
   assert(cnt == nnz);
   
   BMSfreeBlockMemoryArray(oracle->blkmem, &colnz,  oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &colnnz, oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &collen, oracle->nvars);
   
   if( offset != NULL )
      *offset = oracle->heslagoffsets;
   if( col != NULL )
      *col = oracle->heslagcols;

   return SCIP_OKAY;
}

/** Evaluates the Hessian matrix of the Lagrangian in a given point.
 * 
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity.
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity at least ones before using this function. 
 *  Only elements of the lower left triangle and the diagonal are computed.
 */
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambda,             /**< weights (Lagrangian multipliers) for the constraints */ 
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */  
   )
{  /*lint --e{715}*/
   int i;

   assert(oracle != NULL);
   assert(x != NULL);
   assert(lambda != NULL);
   assert(hessian != NULL);

   assert(oracle->heslagoffsets != NULL);
   assert(oracle->heslagcols != NULL);
   
   for( i = oracle->heslagoffsets[oracle->nvars] - 1; i >= 0; --i )
      hessian[i] = 0.0;
   
   if( objfactor != 0.0 )
   {
      if( oracle->objquadlen != 0 )
      {
         assert(oracle->objquadrows != NULL);
         assert(oracle->objquadcols != NULL);
         assert(oracle->objquadvals != NULL);
         SCIP_CALL( hessLagAddQuad(objfactor, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals, 
               oracle->objquadlen, oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
      
      if( oracle->objexprtree != NULL )
      {
         assert(oracle->objexprvaridxs != NULL );
#ifdef WITH_NL
         SCIP_CALL( hessLagAddExprtree(oracle, objfactor, x, newx, oracle->objexprvaridxs, oracle->objexprtree, oracle->heslagoffsets, oracle->heslagcols, hessian) );
#else
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
#endif
      }
   }
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( lambda[i] == 0.0 )
         continue;
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
      {
         assert(oracle->consquadrows != NULL);
         assert(oracle->consquadrows[i] != NULL);
         assert(oracle->consquadcols != NULL);
         assert(oracle->consquadcols[i] != NULL);
         assert(oracle->consquadvals != NULL);
         assert(oracle->consquadvals[i] != NULL);
         SCIP_CALL( hessLagAddQuad(lambda[i], oracle->consquadrows[i], oracle->consquadcols[i], oracle->consquadvals[i], 
               oracle->consquadlens[i], oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         assert(oracle->consexprvaridxs != NULL);
         assert(oracle->consexprvaridxs[i] != NULL);
#ifdef WITH_NL
         SCIP_CALL( hessLagAddExprtree(oracle, lambda[i], x, newx, oracle->consexprvaridxs[i], oracle->consexprtrees[i], oracle->heslagoffsets, oracle->heslagcols, hessian) );
#else
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
#endif
      }
   }
   
   return SCIP_OKAY;
}

/** prints the problem to a file. */
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;

   assert(oracle != NULL);
   
   if( file == NULL )
      file = stdout;
   
   SCIPmessageFPrintInfo(file, "NLPI Oracle: %d variables and %d constraints\n", oracle->nvars, oracle->nconss);
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%10s", oracle->varnames[i]);
      else
         SCIPmessageFPrintInfo(file, "x%09d", i);
      SCIPmessageFPrintInfo(file, ": [%8g, %8g]\n", oracle->varlbs[i], oracle->varubs[i]);
   }
   
   SCIPmessageFPrintInfo(file, "objective: ");
   SCIP_CALL( printFunction(oracle, file,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals,
      oracle->objexprtree) );
   if( oracle->objconstant != 0.0 )
      SCIPmessageFPrintInfo(file, "%+g", oracle->objconstant);
   SCIPmessageFPrintInfo(file, "\n");
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%10s", oracle->consnames[i]);
      else
         SCIPmessageFPrintInfo(file, "con%07d", i);
      
      SCIPmessageFPrintInfo(file, ": ");
      if( oracle->conslhss[i] > -oracle->infinity && oracle->consrhss[i] < oracle->infinity && oracle->conslhss[i] != oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, "%g <= ", oracle->conslhss[i]);
      
      SCIP_CALL( printFunction(oracle, file,
         (oracle->conslinlens   != NULL ? oracle->conslinlens [i] : 0),
         (oracle->conslinlens   != NULL ? oracle->conslinidxs [i] : NULL),
         (oracle->conslinlens   != NULL ? oracle->conslincoefs[i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadlens[i] : 0),
         (oracle->consquadlens  != NULL ? oracle->consquadrows[i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadcols[i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadvals[i] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL)) );
      
      if( oracle->conslhss[i] == oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, " = %g", oracle->consrhss[i]);
      else if( oracle->consrhss[i] <  oracle->infinity )
         SCIPmessageFPrintInfo(file, " <= %g", oracle->consrhss[i]);
      else if( oracle->conslhss[i] > -oracle->infinity )
         SCIPmessageFPrintInfo(file, " >= %g", oracle->conslhss[i]);
      
      SCIPmessageFPrintInfo(file, "\n");
   }

   return SCIP_OKAY;
}

/** prints the problem to a file in GAMS format */
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;

   assert(oracle != NULL);
   
   if( file == NULL )
      file = stdout;
   
   SCIPmessageFPrintInfo(file, "$offlisting\n");
   SCIPmessageFPrintInfo(file, "* NLPI Oracle Problem\n");
   SCIPmessageFPrintInfo(file, "Variables ");
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%s, ", oracle->varnames[i]);
      else
         SCIPmessageFPrintInfo(file, "x%d, ", i);
   }
   SCIPmessageFPrintInfo(file, "OBJVAR;\n\n");
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varlbs[i] == oracle->varubs[i] )
      {
         if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
            SCIPmessageFPrintInfo(file, "%s", oracle->varnames[i]);
         else
            SCIPmessageFPrintInfo(file, "x%d", i);
         SCIPmessageFPrintInfo(file, ".fx = %g;\t", oracle->varlbs[i]);
      }
      else
      {
         if( oracle->varlbs[i] > -oracle->infinity )
         {
            if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
               SCIPmessageFPrintInfo(file, "%s", oracle->varnames[i]);
            else
               SCIPmessageFPrintInfo(file, "x%d", i);
            SCIPmessageFPrintInfo(file, ".lo = %g;\t", oracle->varlbs[i]);
         }
         if( oracle->varubs[i] <  oracle->infinity )
         {
            if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
               SCIPmessageFPrintInfo(file, "%s", oracle->varnames[i]);
            else
               SCIPmessageFPrintInfo(file, "x%d", i);
            SCIPmessageFPrintInfo(file, ".up = %g;\t", oracle->varubs[i]);
         }
      }
      if( initval != NULL )
      {
         if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
            SCIPmessageFPrintInfo(file, "%s", oracle->varnames[i]);
         else
            SCIPmessageFPrintInfo(file, "x%d", i);
         SCIPmessageFPrintInfo(file, ".l = %g;\t", initval[i]);
      }
      SCIPmessageFPrintInfo(file, "\n");
   }
   SCIPmessageFPrintInfo(file, "\n");
   
   SCIPmessageFPrintInfo(file, "Equations ");
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%s, ", oracle->consnames[i]);
      else
         SCIPmessageFPrintInfo(file, "e%d, ", i);

      if( oracle->conslhss[i] > -oracle->infinity && oracle->consrhss[i] < oracle->infinity && oracle->conslhss[i] != oracle->consrhss[i] )
      { /* ranged row: add second constraint */
         if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
            SCIPmessageFPrintInfo(file, "%s_RNG, ", oracle->consnames[i]);
         else
            SCIPmessageFPrintInfo(file, "e%d_RNG, ", i);
      }
   }
   SCIPmessageFPrintInfo(file, "OBJ;\n\n");
   
   SCIPmessageFPrintInfo(file, "OBJ.. OBJVAR =E= ");
   SCIP_CALL( printFunction(oracle, file,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals,
      oracle->objexprtree) );
   if( oracle->objconstant != 0.0 )
      SCIPmessageFPrintInfo(file, "%+g", oracle->objconstant);
   SCIPmessageFPrintInfo(file, ";\n");
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%s", oracle->consnames[i]);
      else
         SCIPmessageFPrintInfo(file, "e%d", i);
      SCIPmessageFPrintInfo(file, ".. ");
      
      SCIP_CALL( printFunction(oracle, file,
         (oracle->conslinlens   != NULL ? oracle->conslinlens [i] : 0),
         (oracle->conslinlens   != NULL ? oracle->conslinidxs [i] : NULL),
         (oracle->conslinlens   != NULL ? oracle->conslincoefs[i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadlens[i] : 0),
         (oracle->consquadlens  != NULL ? oracle->consquadrows[i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadcols[i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadvals[i] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL)) );

      if( oracle->conslhss[i] == oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, " =E= %g", oracle->consrhss[i]);
      else if( oracle->consrhss[i] <  oracle->infinity )
         SCIPmessageFPrintInfo(file, " =L= %g", oracle->consrhss[i]);
      else if( oracle->conslhss[i] > -oracle->infinity )
         SCIPmessageFPrintInfo(file, " =G= %g", oracle->conslhss[i]);
      else
         SCIPmessageFPrintInfo(file, " =N= 0");
      SCIPmessageFPrintInfo(file, ";\n");

      if( oracle->conslhss[i] > -oracle->infinity && oracle->consrhss[i] < oracle->infinity && oracle->conslhss[i] != oracle->consrhss[i] )
      {
         if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
            SCIPmessageFPrintInfo(file, "%s", oracle->consnames[i]);
         else
            SCIPmessageFPrintInfo(file, "e%d", i);
         SCIPmessageFPrintInfo(file, "_RNG.. ");
        
         SCIP_CALL( printFunction(oracle, file,
            (oracle->conslinlens   != NULL ? oracle->conslinlens [i] : 0),
            (oracle->conslinlens   != NULL ? oracle->conslinidxs [i] : NULL),
            (oracle->conslinlens   != NULL ? oracle->conslincoefs[i] : NULL),
            (oracle->consquadlens  != NULL ? oracle->consquadlens[i] : 0),
            (oracle->consquadlens  != NULL ? oracle->consquadrows[i] : NULL),
            (oracle->consquadlens  != NULL ? oracle->consquadcols[i] : NULL),
            (oracle->consquadlens  != NULL ? oracle->consquadvals[i] : NULL),
            (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL)) );
         
         SCIPmessageFPrintInfo(file, " =G= %g;\n", oracle->conslhss[i]);
      }
   }
   
   SCIPmessageFPrintInfo(file, "Model m / all /;\n");
   SCIPmessageFPrintInfo(file, "option limrow = 0;\n");
   SCIPmessageFPrintInfo(file, "option limcol = 0;\n");
   SCIPmessageFPrintInfo(file, "Solve m minimizing OBJVAR using NLP;\n");

   return SCIP_OKAY;
}

/**@} */
