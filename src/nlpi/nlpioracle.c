/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpioracle.c,v 1.30 2010/11/04 07:36:03 bzfviger Exp $"

/**@file    nlpioracle.c
 * @brief   implementation of NLPI oracle interface
 * @author  Stefan Vigerske
 * 
 * @todo jacobi evaluation should be sparse
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpioracle.h"
#include "nlpi/pub_expr.h"
#include "nlpi/exprinterpret.h"
#include "scip/pub_misc.h"

#include <string.h> /* for strlen */

/**@name NLPI Oracle data structures */
/**@{ */

struct SCIP_NlpiOracle
{
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SCIP_Real             infinity;           /**< value for infinity */
   char*                 name;               /**< name of problem */
   
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
   SCIP_QUADELEM**       consquadelems;      /**< quadelems[.] gives elements of quadratic part, or NULL if no quadratic part in this constraint */
   int**                 consexprvaridxs;    /**< exprvaridx[.] gives indices of variables from expression in NLP, or NULL if no nonquadratic part in this constraint */
   SCIP_EXPRTREE**       consexprtrees;      /**< exprtrees[.] gives nonquadratic part, or NULL if no nonquadratic part in this constraint */
   char**                consnames;          /**< array with constraint names */

   SCIP_Real             objconstant;        /**< constant part of objective */
   int                   objnlin;            /**< number of linear variable coefficients in objective */ 
   int*                  objlinidxs;         /**< array with indices of linear variables in objective, or NULL if no linear part */
   SCIP_Real*            objlinvals;         /**< array with coefficients of linear variables in objective, or NULL if no linear part */
   int                   objquadlen;         /**< length of quadratic part of objective, or 0 if no quadratic part in objective */
   SCIP_QUADELEM*        objquadelems;       /**< array with elements of quadratic part in objective, or NULL if no quadratic part in objective */
   int*                  objexprvaridxs;     /**< exprvaridxs[.] gives indices of variables from expression in NLP, or NULL if no nonquadratic part in objective */
   SCIP_EXPRTREE*        objexprtree;        /**< expression tree of nonquadratic part in objective, or NULL if no nonquadratic part in objective */

   int*                  jacoffsets;         /**< rowwise jacobi sparsity pattern: constraint offsets in jaccols */
   int*                  jaccols;            /**< rowwise jacobi sparsity pattern: indices of variables appearing in constraints */
   
   int*                  heslagoffsets;      /**< rowwise sparsity pattern of hessian matrix of Lagrangian: row offsets in heslagcol */
   int*                  heslagcols;         /**< rowwise sparsity pattern of hessian matrix of Lagrangian: column indices; sorted for each row */
   
   int*                  vardegrees;         /**< array with maximal degree of variable over objective and all constraints */
   SCIP_Bool             vardegreesuptodate; /**< whether the variable degrees are up to date */
   
   SCIP_EXPRINT*         exprinterpreter;    /**< interpreter for expression trees: evaluation and derivatives */
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
      oracle->consquadlens [toidx] = oracle->consquadlens [fromidx];
      oracle->consquadelems[toidx] = oracle->consquadelems[fromidx];
      
      oracle->consquadlens [fromidx] = 0;
      oracle->consquadelems[fromidx] = NULL;
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
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadelems[considx], oracle->consquadlens[considx]);
      oracle->consquadlens[considx] = 0;
   }
   
   if( oracle->consexprtrees && oracle->consexprtrees[considx] != NULL )
   {
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consexprvaridxs[considx], SCIPexprtreeGetNVars(oracle->consexprtrees[considx]));
      SCIP_CALL( SCIPexprtreeFree(&oracle->consexprtrees[considx]) );
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
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consquadelems,   oracle->nconss);

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
   assert(!oracle->vardegreesuptodate || oracle->vardegrees[toidx] == -1);
   
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
   
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->varlbs,     oracle->nvars);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->varubs,     oracle->nvars);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->varnames,   oracle->nvars);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->vardegrees, oracle->nvars);
   
   oracle->nvars = 0;
   oracle->vardegreesuptodate = TRUE; 
   
   return SCIP_OKAY;
}

/** Updates the degrees of all variables. */
static
void updateVariableDegrees(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   int i;
   int j;
   
   assert(oracle != NULL);
   assert(oracle->nvars == 0 || oracle->vardegrees != NULL);
   
   if( oracle->vardegreesuptodate || oracle->nvars == 0 )
      return;
   
   /* assume all variables do not appear in NLP */
   BMSclearMemoryArray(oracle->vardegrees, oracle->nvars);

   if( oracle->objlinidxs != NULL )
      for( j = 0; j < oracle->objnlin; ++j )
         oracle->vardegrees[oracle->objlinidxs[j]] = 1;
   
   if( oracle->conslinidxs != NULL )
      for( i = 0; i < oracle->nconss; ++i )
         if( oracle->conslinidxs[i] != NULL )
            for( j = 0; j < oracle->conslinlens[i]; ++j )
               oracle->vardegrees[oracle->conslinidxs[i][j]] = 1;
   
   if( oracle->objquadelems != NULL )
      for( j = 0; j < oracle->objquadlen; ++j )
      {
         oracle->vardegrees[oracle->objquadelems[j].idx1] = 2;
         oracle->vardegrees[oracle->objquadelems[j].idx2] = 2;
      }
   
   if( oracle->consquadelems != NULL )
      for( i = 0; i < oracle->nconss; ++i )
         if( oracle->consquadelems[i] != NULL )
            for( j = 0; j < oracle->consquadlens[i]; ++j )
            {
               oracle->vardegrees[oracle->consquadelems[i][j].idx1] = 2;
               oracle->vardegrees[oracle->consquadelems[i][j].idx2] = 2;
            }
   
   /* @todo could use exprtreeGetDegree to get actual degree of a variable in tree */
   if( oracle->objexprvaridxs != NULL )
      for( j = 0; j < SCIPexprtreeGetNVars(oracle->objexprtree); ++j )
         oracle->vardegrees[oracle->objexprvaridxs[j]] = INT_MAX;
   
   if( oracle->consexprvaridxs != NULL )
      for( i = 0; i < oracle->nconss; ++i )
         if( oracle->consexprvaridxs[i] != NULL )
            for( j = 0; j < SCIPexprtreeGetNVars(oracle->consexprtrees[i]); ++j )
               oracle->vardegrees[oracle->consexprvaridxs[i][j]] = INT_MAX;
   
   oracle->vardegreesuptodate = TRUE;
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

/** removes entries with index -1 (marked as deleted) from array of linear elements
 * assumes that array is sorted by index, i.e., all -1 are at the beginning
 * shortens memory accordingly */
static
SCIP_RETCODE clearDeletedLinearElements(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int**                 linidxs,            /**< variable indices */
   SCIP_Real**           coefs,              /**< variable coefficients */
   int*                  nidxs               /**< number of indices */
   )
{
   int i;
   int offset;
   
   assert(blkmem   != NULL);
   assert(linidxs  != NULL);
   assert(*linidxs != NULL);
   assert(coefs    != NULL);
   assert(*coefs   != NULL);
   assert(nidxs    != NULL);
   assert(*nidxs   > 0);
   
   for( offset = 0; offset < *nidxs; ++offset )
      if( (*linidxs)[offset] >= 0 )
         break;

   /* nothing was deleted */
   if( offset == 0 )
      return SCIP_OKAY;

   /* all elements were deleted */
   if( offset == *nidxs )
   {
      BMSfreeBlockMemoryArray(blkmem, linidxs, *nidxs);
      BMSfreeBlockMemoryArray(blkmem, coefs,   *nidxs);
      *nidxs = 0;
      return SCIP_OKAY;
   }
   
   /* some elements were deleted */
   for( i = 0; i < *nidxs - offset; ++i )
   {
      (*linidxs)[i] = (*linidxs)[i+offset];
      (*coefs)[i]   = (*coefs)  [i+offset];
   }
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, linidxs, *nidxs, *nidxs - offset) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, coefs,   *nidxs, *nidxs - offset) );
   *nidxs -= offset;
   
   return SCIP_OKAY;
}

/** removes entries with index pair (-1,-1) (marked as deleted) from array of quadratic elements
 * assumes that array is sorted, i.e., all deleted elements are at the beginning
 * moves and reallocs memory accordingly */
static
SCIP_RETCODE clearDeletedQuadElements(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_QUADELEM**       quadelems,          /**< quadratic elements */
   int*                  nquadelems          /**< number of quadratic elements */
   )
{
   int i;
   int offset;
   
   assert(blkmem      != NULL);
   assert(quadelems   != NULL);
   assert(*quadelems  != NULL);
   assert(nquadelems  != NULL);
   assert(*nquadelems > 0);
   
   for( offset = 0; offset < *nquadelems; ++offset )
   {
      /* either both variables are marked as deleted or none of them */
      assert(((*quadelems)[offset].idx1 >= 0) == ((*quadelems)[offset].idx2 >= 0));
      if( (*quadelems)[offset].idx1 >= 0 )
         break;
   }

   /* nothing was deleted */
   if( offset == 0 )
      return SCIP_OKAY;

   /* all elements were deleted */
   if( offset == *nquadelems )
   {
      BMSfreeBlockMemoryArray(blkmem, quadelems, *nquadelems);
      *nquadelems = 0;
      return SCIP_OKAY;
   }
   
   /* some elements were deleted */
   for( i = 0; i < *nquadelems - offset; ++i )
      (*quadelems)[i] = (*quadelems)[i+offset];
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, quadelems, *nquadelems, *nquadelems - offset) );
   *nquadelems -= offset;
   
   return SCIP_OKAY;
}

#if 0
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
#endif

/** applies a mapping of indices to an array of quadratic elements */
static
void mapIndicesQuad(
   int*                  indexmap,           /**< mapping from old variable indices to new indices */
   int                   nelems,             /**< number of quadratic elements */
   SCIP_QUADELEM*        elems               /**< array of quadratic elements to adjust */
   )
{
   assert(indexmap != NULL);
   assert(nelems == 0 || elems != NULL);

   for( ; nelems ; --nelems, ++elems )
   {
      assert(indexmap[elems->idx1] >= 0);
      assert(indexmap[elems->idx2] >= 0);
      elems->idx1 = indexmap[elems->idx1];
      elems->idx2 = indexmap[elems->idx2];
      /* swap indices if not idx1 <= idx2 */
      if( elems->idx1 > elems->idx2 )
      {
         int tmp = elems->idx2;
         elems->idx2 = elems->idx1;
         elems->idx1 = tmp;
      }
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
   const SCIP_QUADELEM*  quadelems,          /**< elements of quadratic part matrix */
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
      assert(quadelems != NULL);
      assert(x != NULL);
      
      for( ; quadlen > 0; --quadlen, ++quadelems )
         *val += quadelems->coef * x[quadelems->idx1] * x[quadelems->idx2];
   }
   
   if( exprtree != NULL )
   {
      SCIP_Real* xx;
      int        i;
      SCIP_Real  nlval;
      int        nvars;
      
      nvars = SCIPexprtreeGetNVars(exprtree);
      
      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) );
      for (i = 0; i < nvars; ++i)
      {
         assert(exprvaridxs[i] >= 0);
         xx[i] = x[exprvaridxs[i]];  /*lint !e613*/
      }
      
      SCIP_CALL( SCIPexprintEval(oracle->exprinterpreter, exprtree, xx, &nlval) );
      if( nlval != nlval || ABS(nlval) >= oracle->infinity )  /*lint !e777*/
         *val = nlval;
      else
         *val += nlval;
      
      BMSfreeBlockMemoryArray(oracle->blkmem, &xx, nvars);
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
   const SCIP_QUADELEM*  quadelems,          /**< elements of quadratic part matrix */
   int*                  exprvaridxs,        /**< indices of variables in nonquadratic part */
   SCIP_EXPRTREE*        exprtree,           /**< nonquadratic part */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
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
      assert(quadelems != NULL);
      assert(x != NULL);
      
      for( ; quadlen > 0; --quadlen, ++quadelems )
      {
         tmp = quadelems->coef * x[quadelems->idx1];
         *val += tmp * x[quadelems->idx2];
         grad[quadelems->idx2] += tmp;
         grad[quadelems->idx1] += quadelems->coef * x[quadelems->idx2];
      }
   }

   if( exprtree != NULL )
   {
      SCIP_Real* xx;
      SCIP_Real* g;
      int        i;
      SCIP_Real  nlval;
      int        nvars;
      
      xx = NULL;
      nvars = SCIPexprtreeGetNVars(exprtree);

      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &g, nvars) );

      if( isnewx )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) );
         for (i = 0; i < nvars; ++i)
         {
            assert(exprvaridxs[i] >= 0);
            xx[i] = x[exprvaridxs[i]];  /*lint !e613*/
         }
      }
      
      SCIPdebugMessage("eval gradient of ");
      SCIPdebug( SCIPexprtreePrint(exprtree, NULL, NULL, NULL) );
      SCIPdebug( printf("\nx ="); for( i = 0; i < nvars; ++i) printf(" %g", xx[i]); printf("\n"); )
      
      SCIP_CALL( SCIPexprintGrad(oracle->exprinterpreter, exprtree, xx, isnewx, &nlval, g) );
      
      SCIPdebug( printf("g ="); for( i = 0; i < nvars; ++i) printf(" %g", g[i]); printf("\n"); )

      if( nlval != nlval || ABS(nlval) >= oracle->infinity )  /*lint !e777*/
      {
         BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
         BMSfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
         SCIPdebugMessage("gradient evaluation yield invalid function value %g\n", nlval);
         return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
      }
      else
      {
         *val += nlval;
         for( i = 0; i < nvars; ++i )
            if( g[i] != g[i] )  /*lint !e777*/
            {
               SCIPdebugMessage("gradient evaluation yield invalid gradient value %g\n", g[i]);
               BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
               BMSfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
               return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
            }
            else
            {
               grad[exprvaridxs[i]] += g[i];
            }
      }
      
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
      BMSfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
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
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, intarray, *len) );
   }
   else
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, intarray, oldsize, *len) );
   }
   
   return SCIP_OKAY;
}

/** collects nonzeros entries in colnz and increases the nzcount given indices of quadratic terms */
static
SCIP_RETCODE hessLagSparsitySetNzFlagForQuad(
   SCIP_NLPIORACLE*      oracle,             /**< NLPI oracle */
   int**                 colnz,              /**< indices of nonzero entries for each column */
   int*                  collen,             /**< space allocated to store indices of nonzeros for each column */
   int*                  colnnz,             /**< number of nonzero entries for each column */
   int*                  nzcount,            /**< counter for total number of nonzeros; should be increased whenever some colnnz is increased */
   int                   length,             /**< length of quadratic part */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements */
   )
{
   int pos;

   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(quadelems != NULL);
   assert(length >= 0);

   for( ; length > 0; --length, ++quadelems )
   {
      assert(quadelems->idx1 <= quadelems->idx2);
      
      if( colnz[quadelems->idx2] == NULL || !SCIPsortedvecFindInt(colnz[quadelems->idx2], quadelems->idx1, colnnz[quadelems->idx2], &pos) )
      {
         if( colnnz[quadelems->idx2] >= collen[quadelems->idx2] )
         { /* allocate more space so that we can insert one more element */
            SCIP_CALL( increaseIntArraySize(oracle->blkmem, &colnz[quadelems->idx2], &collen[quadelems->idx2]) );
         }
         assert(collen[quadelems->idx2] > colnnz[quadelems->idx2]);
         
         SCIPsortedvecInsertInt(colnz[quadelems->idx2], quadelems->idx1, &colnnz[quadelems->idx2]);
         ++(*nzcount);
      }
   }
   
   return SCIP_OKAY;
}

/** collects indices of nonzero entries in the lower-left part of the hessian matrix of an expression
 * adds the indices to a given set of indices, avoiding duplicates */
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
   int         i;
   int         j;
   int         nvars;
   int         nn;
   int         row;
   int         col;
   int         pos;
   
   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(exprvaridx != NULL);
   assert(exprtree != NULL);
   assert(dim >= 0);
   
   nvars = SCIPexprtreeGetNVars(exprtree);
   nn = nvars * nvars;
   
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &x,     nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &hesnz, nn) );
   
   for( i = 0; i < nvars; ++i )
      x[i] = 2.0; /* hope that this value does not make much trouble for the evaluation routines */
   
   SCIP_CALL( SCIPexprintHessianSparsityDense(oracle->exprinterpreter, exprtree, x, hesnz) );
   
   for( i = 0; i < nvars; ++i ) /* rows */
      for( j = 0; j <= i; ++j ) /* cols */
      {
         if( !hesnz[i*nvars + j] )
            continue;
         
         row = MAX(exprvaridx[i], exprvaridx[j]);
         col = MIN(exprvaridx[i], exprvaridx[j]);
         
         assert(row <  dim);
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
   
   BMSfreeBlockMemoryArray(oracle->blkmem, &x, nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &hesnz, nn);
   
   return SCIP_OKAY;
}

/** adds quadratic part into hessian structure */
static
SCIP_RETCODE hessLagAddQuad(
   SCIP_Real             weight,             /**< weight of quadratic part */
   int                   length,             /**< number of elements in matrix of quadratic part */
   SCIP_QUADELEM*        quadelems,          /**< elements in matrix of quadratic part */
   int*                  hesoffset,          /**< row offsets in sparse matrix that is to be filled */ 
   int*                  hescol,             /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values              /**< buffer for values of sparse matrix that is to be filled */
   )
{
   int idx;

   assert(length >= 0);
   assert(quadelems != NULL);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   for( ; length > 0; --length, ++quadelems )
   {
      assert(quadelems->idx1 <= quadelems->idx2);
      if( !SCIPsortedvecFindInt(&hescol[hesoffset[quadelems->idx2]], quadelems->idx1, hesoffset[quadelems->idx2 + 1] - hesoffset[quadelems->idx2], &idx) )
      {
         SCIPerrorMessage("Could not find entry in hessian sparsity\n");
         return SCIP_ERROR;
      }
      values[hesoffset[quadelems->idx2] + idx] += weight * ((quadelems->idx1 == quadelems->idx2) ? 2 * quadelems->coef : quadelems->coef);
   }

   return SCIP_OKAY;
}

/** adds hessian of an expression into hessian structure */
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
   SCIP_Real* xx;
   SCIP_Real* h;
   SCIP_Real* hh;
   int        i;
   int        j;
   int        nvars;
   int        nn;
   int        row;
   int        col;
   int        idx;
   SCIP_Real  val;
   
   xx = NULL;

   assert(oracle != NULL);
   assert(x != NULL || new_x == FALSE);
   assert(exprvaridx != NULL);
   assert(exprtree != NULL);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   nvars = SCIPexprtreeGetNVars(exprtree);
   nn = nvars * nvars;

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &h, nn) );

   if( new_x )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) );
      for( i = 0; i < nvars; ++i )
      {
         assert(exprvaridx[i] >= 0);
         xx[i] = x[exprvaridx[i]];  /*lint !e613*/
      }
   }
   
   SCIP_CALL( SCIPexprintHessianDense(oracle->exprinterpreter, exprtree, xx, new_x, &val, h) );
   if( val != val )  /*lint !e777*/
   {
      SCIPdebugMessage("hessian evaluation yield invalid function value %g\n", val);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
      BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
      return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
   }

   hh = h;
   for( i = 0; i < nvars; ++i ) /* rows */
   {
      for( j = 0; j <= i; ++j, ++hh ) /* cols */
      {
         if( !*hh )
            continue;
         
         if( *hh != *hh )  /*lint !e777*/
         {
            SCIPdebugMessage("hessian evaluation yield invalid hessian value %g\n", *hh);
            BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
            BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
            return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
         }
         
         row = MAX(exprvaridx[i], exprvaridx[j]);
         col = MIN(exprvaridx[i], exprvaridx[j]);
         
         if( !SCIPsortedvecFindInt(&hescol[hesoffset[row]], col, hesoffset[row+1] - hesoffset[row], &idx) )
         {
            SCIPerrorMessage("Could not find entry (%d, %d) in hessian sparsity\n", row, col);
            BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
            BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
            return SCIP_ERROR;
         }
         
         values[hesoffset[row] + idx] += weight * *hh;
      }
      hh += nvars - j;
   }
   
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);

   return SCIP_OKAY;
}

/** prints a name, if available, makes sure it has not more than 64 characters, and adds a unique prefix if the longnames flag is set */
static
void printName(
   char*                 buffer,             /**< buffer to print to, has to be not NULL */
   char**                names,              /**< array of names, or NULL */
   int                   idx,                /**< index which name we want to print */
   char                  prefix,             /**< a letter (typically 'x' or 'e') to distinguish variable and equation names, if names[idx] is not available */
   const char*           suffix,             /**< a suffer to add to the name, or NULL */
   SCIP_Bool             longnames           /**< whether prefixes for long names should be added */
)
{
   if( longnames )
   {
      if( names != NULL && names[idx] != NULL )
         sprintf(buffer, "%c%05d%.*s%s", prefix, idx, suffix ? (int)(57-strlen(suffix)) : 57, names[idx], suffix ? suffix : "");
      else
         sprintf(buffer, "%c%05d", prefix, idx);
   }
   else
   {
      if( names != NULL && names[idx] != NULL )
      {
         assert(strlen(names[idx]) + (suffix ? strlen(suffix) : 0) <= 64);
         sprintf(buffer, "%s%s", names[idx], suffix ? suffix : "");
      }
      else
         sprintf(buffer, "%c%d%s", prefix, idx, suffix ? suffix : "");
   }
}

/** prints a function */ 
static
SCIP_RETCODE printFunction(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file,               /**< file to print to, has to be not NULL */
   int                   nlin,               /**< length of linear part */
   const int*            lininds,            /**< indices of linear variables */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables */
   int                   quadlen,            /**< number of indices in quadratic part matrix */
   const SCIP_QUADELEM*  quadelems,          /**< elements in quadratic part matrix */
   SCIP_EXPRTREE*        exprtree,           /**< nonquadratic part */
   int*                  exprvaridxs,        /**< indices of variables in nonquadratic part */
   SCIP_Bool             longvarnames,       /**< whether variable names need to be shorten to 64 characters */
   SCIP_Bool             longequnames        /**< whether equation names need to be shorten to 64 characters */
   )
{  /*lint --e{715}*/
   int i;
   int j;
   char namebuf[70];

   assert(oracle != NULL);
   assert(file != NULL);
   
   for( i = 0; i < nlin; ++i )
   {
      printName(namebuf, oracle->varnames, lininds[i], 'x', NULL, longvarnames);
      SCIPmessageFPrintInfo(file, "%+.20g*%s", linvals[i], namebuf);
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(file, "\n");
   }
   
   if( quadlen != 0 )
   {
      j = 0;
      for( i = 0; i < quadlen; ++i )
      {
         printName(namebuf, oracle->varnames, quadelems[i].idx1, 'x', NULL, longvarnames);
         SCIPmessageFPrintInfo(file, "%+.20g*%s", quadelems[j].coef, namebuf);
         printName(namebuf, oracle->varnames, quadelems[i].idx2, 'x', NULL, longvarnames);
         SCIPmessageFPrintInfo(file, "*%s", namebuf);
         if( i % 10 == 9 )
            SCIPmessageFPrintInfo(file, "\n");
      }
   }

   if( exprtree != NULL )
   {
      char** varnames;
      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &varnames, SCIPexprtreeGetNVars(exprtree)) );

      /* setup variable names */
      for( i = 0; i < SCIPexprtreeGetNVars(exprtree); ++i )
      {
         assert(exprvaridxs[i] < 1e+20);
         BMSallocBlockMemoryArray(oracle->blkmem, &varnames[i], 70);
         printName(varnames[i], oracle->varnames, exprvaridxs[i], 'x', NULL, longvarnames);
      }

      SCIPmessageFPrintInfo(file, " +");
      SCIPexprtreePrint(exprtree, file, (const char**)varnames, NULL);

      for( i = 0; i < SCIPexprtreeGetNVars(exprtree); ++i )
      {
         BMSfreeBlockMemoryArray(oracle->blkmem, &varnames[i], 70);
      }
      BMSfreeBlockMemoryArray(oracle->blkmem, &varnames, SCIPexprtreeGetNVars(exprtree));
   }

   return SCIP_OKAY;
}

/** returns whether an expression is contains nonsmooth operands (min, max, abs, ...) */
static
SCIP_Bool exprIsNonSmooth(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int i;
   
   assert(expr != NULL);
   assert(SCIPexprGetChildren(expr) != NULL || SCIPexprGetNChildren(expr) == 0);

   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      if( exprIsNonSmooth(SCIPexprGetChildren(expr)[i]) )
         return TRUE;
   }

   switch( SCIPexprGetOperator(expr) )
   {
      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      case SCIP_EXPR_ABS:
      case SCIP_EXPR_SIGN:
      case SCIP_EXPR_SIGNPOWER:
         return TRUE;

      default: ;
   } /*lint !e788*/

   return FALSE;
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
   
   SCIP_ALLOC( BMSallocMemory(oracle) );
   BMSclearMemory(*oracle);
   
   (*oracle)->blkmem   = blkmem;
   (*oracle)->infinity = SCIP_DEFAULT_INFINITY;
   (*oracle)->vardegreesuptodate = TRUE;
   
   SCIPdebugMessage("Oracle initializes expression interpreter %s\n", SCIPexprintGetName());
   SCIP_CALL( SCIPexprintCreate(blkmem, &(*oracle)->exprinterpreter) );

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
   
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objlinidxs,   (*oracle)->objnlin);
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objlinvals,   (*oracle)->objnlin);
   BMSfreeBlockMemoryArrayNull((*oracle)->blkmem, &(*oracle)->objquadelems, (*oracle)->objquadlen);
   if( (*oracle)->objexprtree != NULL )
   {
      BMSfreeBlockMemoryArray((*oracle)->blkmem, &(*oracle)->objexprvaridxs, SCIPexprtreeGetNVars((*oracle)->objexprtree));
      SCIP_CALL( SCIPexprtreeFree(&(*oracle)->objexprtree) );
   }

   SCIP_CALL( freeConstraints(*oracle) );
   SCIP_CALL( freeVariables(*oracle) );

   SCIP_CALL( SCIPexprintFree(&(*oracle)->exprinterpreter) );
   
   if( (*oracle)->name != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleSetProblemName(*oracle, NULL) );
   }

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

/** sets the problem name (used for printing) */
SCIP_RETCODE SCIPnlpiOracleSetProblemName(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const char*           name                /**< name of problem */
   )
{
   assert(oracle != NULL);

   if( oracle->name != NULL )
   {
      BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->name, strlen(oracle->name)+1);
   }

   if( name != NULL )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->name, name, strlen(name)+1) );
   }
   
   return SCIP_OKAY;
}

/** gets the problem name, or NULL if none set */
const char* SCIPnlpiOracleGetProblemName(
   SCIP_NLPIORACLE*     oracle               /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->name;
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
   
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varlbs, oracle->nvars, oracle->nvars + nvars) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varubs, oracle->nvars, oracle->nvars + nvars) );
   
   if( lbs != NULL )
   {
      BMScopyMemoryArray(&((oracle->varlbs)[oracle->nvars]), lbs, nvars);  /*lint !e866*/
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varlbs[oracle->nvars+i] = -oracle->infinity;
   
   if( ubs != NULL )
   {
      BMScopyMemoryArray(&((oracle->varubs)[oracle->nvars]), ubs, nvars);  /*lint !e866*/
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varubs[oracle->nvars+i] =  oracle->infinity;
   
   if( oracle->varnames != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varnames, oracle->nvars, oracle->nvars + nvars) );
      if( varnames != NULL )
      {
         for( i = 0; i < nvars; ++i )
         {
            if( varnames[i] != NULL )
            {
               SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->varnames)[oracle->nvars+i]), varnames[i], strlen(varnames[i])+1) );  /*lint !e866*/
            }
            else
               oracle->varnames[oracle->nvars+i] = NULL;
         }
      }
      else
      {
         BMSclearMemoryArray(&((oracle->varnames)[oracle->nvars]), nvars);  /*lint !e866*/
      }
   }
   else if( varnames != NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->varnames), oracle->nvars + nvars) );
      BMSclearMemoryArray(oracle->varnames, oracle->nvars);
      for( i = 0; i < nvars; ++i )
      {
         if( varnames[i] != NULL )
         {
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->varnames)[oracle->nvars+i]), varnames[i], strlen(varnames[i])+1) );  /*lint !e866*/
         }
         else
            oracle->varnames[oracle->nvars+i] = NULL;
      }
   }
   
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->vardegrees), oracle->nvars, oracle->nvars + nvars) );  /*lint !e866*/
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
   const int*            nquadelems,         /**< number of elements in matrix of quadratic part for each constraint,
                                              * may be NULL in case of no quadratic part in any constraint */
   SCIP_QUADELEM* const* quadelems,          /**< quadratic elements specifying quadratic part for each constraint, entry of array may be NULL in case of no quadratic part,
                                              * may be NULL in case of no quadratic part in any constraint */
   int* const*           exprvaridxs,        /**< NULL if no nonquadratic parts, otherwise epxrvaridxs[.] maps variable indices in expression tree to indices in nlp */
   SCIP_EXPRTREE* const* exprtrees,          /**< NULL if no nonquadratic parts, otherwise exprtrees[.] gives nonquadratic part, 
                                              *   or NULL if no nonquadratic part in this constraint */
   const char**          consnames           /**< names of new constraints, or NULL if no names should be stored */
   )
{  /*lint --e{715}*/
   SCIP_Bool addednlcon;  /* whether a nonlinear constraint was added */
   int       i;
   int       j;
   int       newlen;

   assert(oracle != NULL);
   
   if( nconss == 0 )
      return SCIP_OKAY;
   
   assert(nconss > 0);

   addednlcon = FALSE;

   invalidateJacobiSparsity(oracle); /* @TODO we could also update (extend) the sparsity pattern */

   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslhss, oracle->nconss, oracle->nconss + nconss) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consrhss, oracle->nconss, oracle->nconss + nconss) );

   if( lhss != NULL )
   {
      BMScopyMemoryArray(&((oracle->conslhss)[oracle->nconss]), lhss, nconss);  /*lint !e866*/
   }
   else
      for( i = 0; i < nconss; ++i )
         oracle->conslhss[oracle->nconss+i] = -oracle->infinity;
   
   if( rhss != NULL )
   {
      BMScopyMemoryArray(&((oracle->consrhss)[oracle->nconss]), rhss, nconss);  /*lint !e866*/
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
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss, oracle->nconss + nconss) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss, oracle->nconss + nconss) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss, oracle->nconss + nconss) );
      }
      else
      { /* had no linear parts so far */
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss + nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss + nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->conslinlens,  oracle->nconss);
         BMSclearMemoryArray(oracle->conslinidxs,  oracle->nconss);
         BMSclearMemoryArray(oracle->conslincoefs, oracle->nconss);
      }
      
      for( i = 0; i < nconss; ++i )
      {
         oracle->conslinlens[oracle->nconss+i] = nlininds[i];
         if( oracle->conslinlens[oracle->nconss+i] )
         {
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs [oracle->nconss+i], lininds[i], oracle->conslinlens[oracle->nconss+i]) );  /*lint !e866*/
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs[oracle->nconss+i], linvals[i], oracle->conslinlens[oracle->nconss+i]) );  /*lint !e866*/
            SCIPsortIntReal(oracle->conslinidxs[oracle->nconss+i], oracle->conslincoefs[oracle->nconss+i], oracle->conslinlens[oracle->nconss+i]);
            for( j = 0; j < nlininds[i]; ++j )
            {
               assert(lininds[i][j] >= 0);
               assert(lininds[i][j] <= oracle->nvars);
               assert(j == 0 || lininds[i][j-1] != lininds[i][j]); /* in gradient evaluation, we do not like duplicate indices */
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
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss, oracle->nconss + nconss) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss, oracle->nconss + nconss) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss, oracle->nconss + nconss) );
      BMSclearMemoryArray(&oracle->conslinlens [oracle->nconss], nconss);  /*lint !e866*/
      BMSclearMemoryArray(&oracle->conslinidxs [oracle->nconss], nconss);  /*lint !e866*/
      BMSclearMemoryArray(&oracle->conslincoefs[oracle->nconss], nconss);  /*lint !e866*/
   }
   
   if( nquadelems != NULL )
   {
      assert(quadelems != NULL);
      if( (oracle->consquadlens == NULL) && (oracle->nconss != 0) )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens),  oracle->nconss + nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadelems), oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consquadlens,  oracle->nconss);
         BMSclearMemoryArray(oracle->consquadelems, oracle->nconss);
      }
      else
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens),  oracle->nconss, oracle->nconss + nconss) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadelems), oracle->nconss, oracle->nconss + nconss) );
      }
      for( i = 0; i < nconss; ++i )
      {
         if( quadelems[i] != NULL && nquadelems[i] != 0 )
         {
            oracle->consquadlens[oracle->nconss + i] = nquadelems[i];
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->consquadelems[oracle->nconss + i], quadelems[i], nquadelems[i]) );  /*lint !e866*/
            
            /* sort and squeeze quadratic part */
            SCIPquadelemSort(oracle->consquadelems[oracle->nconss+i], nquadelems[i]);
            SCIPquadelemSqueeze(oracle->consquadelems[oracle->nconss+i], nquadelems[i], &newlen);
            assert(newlen >= 0);
            assert(newlen <= nquadelems[i]);
            
            if( newlen == 0 )
            {
               /* quadratic elements turned out to be all 0.0 */
               BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->consquadelems[oracle->nconss+i], nquadelems[i]);  /*lint !e866*/
               assert(oracle->consquadelems[oracle->nconss+i] == NULL);
               oracle->consquadlens[oracle->nconss+i] = 0;
            }
            else
            {
               if( newlen < nquadelems[i] )
               {
                  SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadelems[oracle->nconss+i], nquadelems[i], newlen) );  /*lint !e866*/
                  oracle->consquadlens[oracle->nconss+i] = newlen;
               }
               /* update variable degress for variables in quadratic part to keep uptodate */
               if( oracle->vardegreesuptodate )
                  for( j = 0; j < oracle->consquadlens[oracle->nconss+i]; ++j )
                  {
                     oracle->vardegrees[oracle->consquadelems[oracle->nconss+i][j].idx1] = MAX(2, oracle->vardegrees[oracle->consquadelems[oracle->nconss+i][j].idx1]);
                     oracle->vardegrees[oracle->consquadelems[oracle->nconss+i][j].idx2] = MAX(2, oracle->vardegrees[oracle->consquadelems[oracle->nconss+i][j].idx2]);
                  }

               addednlcon = TRUE;
            }
         }
         else
         {
            oracle->consquadlens [oracle->nconss + i] = 0;
            oracle->consquadelems[oracle->nconss + i] = NULL;
         }
      }
   }
   else if( oracle->consquadlens != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens),  oracle->nconss, oracle->nconss + nconss) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadelems), oracle->nconss, oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consquadlens) [oracle->nconss]), nconss);  /*lint !e866*/
      BMSclearMemoryArray(&((oracle->consquadelems)[oracle->nconss]), nconss);  /*lint !e866*/
   }
   
   if( exprtrees != NULL )
   {
      if( oracle->consexprtrees == NULL && oracle->nconss != 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprtrees),   oracle->nconss + nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprvaridxs), oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consexprtrees,   oracle->nconss);
         BMSclearMemoryArray(oracle->consexprvaridxs, oracle->nconss);
      }
      else
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprtrees),   oracle->nconss, oracle->nconss + nconss) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprvaridxs), oracle->nconss, oracle->nconss + nconss) );
      }
      for( i = 0; i < nconss; ++i )
         if( exprtrees[i] != NULL )
         {
            assert(oracle->exprinterpreter != NULL);
            
            SCIP_CALL( SCIPexprtreeCopy(oracle->blkmem, &oracle->consexprtrees[oracle->nconss + i], exprtrees[i]) );
            
            SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->consexprtrees[oracle->nconss + i]) );
            
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs[oracle->nconss + i], exprvaridxs[i], SCIPexprtreeGetNVars(exprtrees[i])) );  /*lint !e866*/
            /* update variable degrees of variables to keep them up to date */
            if( oracle->vardegreesuptodate )
               for (j = 0; j < SCIPexprtreeGetNVars(exprtrees[i]); ++j)
               {
                  assert(exprvaridxs[i][j] >= 0);
                  assert(exprvaridxs[i][j] <  oracle->nvars);
                  oracle->vardegrees[exprvaridxs[i][j]] = INT_MAX; /* @TODO could try to be more clever, maybe use getMaxDegree function in exprtree */
               }
         }
         else
         {
            oracle->consexprtrees[oracle->nconss + i] = NULL;
            oracle->consexprvaridxs[oracle->nconss + i] = NULL;
         }
   }
   else if( oracle->consexprtrees != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprtrees),   oracle->nconss, oracle->nconss + nconss) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consexprvaridxs), oracle->nconss, oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consexprtrees)  [oracle->nconss]), nconss);  /*lint !e866*/
      BMSclearMemoryArray(&((oracle->consexprvaridxs)[oracle->nconss]), nconss);  /*lint !e866*/
   }

   if( consnames != NULL )
   {
      if( oracle->consnames == NULL && oracle->nconss != 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consnames), oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consnames, oracle->nconss);
      }
      else
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consnames), oracle->nconss, oracle->nconss + nconss) );
      }
      for( i = 0; i < nconss; ++i )
         if( consnames[i] != NULL )
         {
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &((oracle->consnames)[oracle->nconss + i]), consnames[i], strlen(consnames[i])+1) );  /*lint !e866*/
         }
         else
            oracle->consnames[oracle->nconss + i] = NULL;
   }
   else if( oracle->consnames != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &(oracle->consnames), oracle->nconss, oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consnames)[oracle->nconss]), nconss);  /*lint !e866*/
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
   int                   nquadelems,         /**< number of entries in matrix of quadratic part */
   const SCIP_QUADELEM*  quadelems,          /**< entries in matrix of quadratic part, may be NULL in case of no quadratic part */
   const int*            exprvaridxs,        /**< maps variable indices in expression tree to indices in nlp, or NULL if no nonquadratic part */
   const SCIP_EXPRTREE*  exprtree            /**< expression tree of nonquadratic part, or NULL if no nonquadratic part */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   
   /* clear previous objective */
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objlinidxs,   oracle->objnlin);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objlinvals,   oracle->objnlin);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objquadelems, oracle->objquadlen);
   if( oracle->objexprtree != NULL )
   {
      BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->objexprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree));
      SCIP_CALL( SCIPexprtreeFree(&oracle->objexprtree) );
      
      oracle->vardegreesuptodate = FALSE;
   }
   
   if( nquadelems != 0 || oracle->objquadlen != 0 || exprtree != NULL || oracle->objexprtree != NULL )
      invalidateHessianLagSparsity(oracle);
   
   oracle->objconstant = constant;
   
   if( nlin != 0 )
   {
      int i;

      assert(lininds != NULL);
      assert(linvals != NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objlinidxs, lininds, nlin) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objlinvals, linvals, nlin) );
      oracle->objnlin = nlin;
      if( oracle->vardegreesuptodate )
         for( i = 0; i < nlin; ++i )
            oracle->vardegrees[lininds[i]] = MAX(1, oracle->vardegrees[lininds[i]]);
      
      SCIPsortIntReal(oracle->objlinidxs, oracle->objlinvals, nlin);
   }
   else
      oracle->objnlin = 0;
   
   oracle->objquadlen = nquadelems;
   if( nquadelems != 0 )
   {
      int j;
      int newlen;

      assert(quadelems != NULL);
      
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objquadelems, quadelems, oracle->objquadlen) );

      /* sort and squeeze quadratic elements */
      SCIPquadelemSort(oracle->objquadelems, oracle->objquadlen);
      SCIPquadelemSqueeze(oracle->objquadelems, oracle->objquadlen, &newlen);
      assert(newlen >= 0);
      assert(newlen <= nquadelems);
      
      if( newlen == 0 )
      {
         /* quadratic elements turned out to be all 0.0 */
         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->objquadelems, nquadelems);
         assert(oracle->objquadelems == NULL);
         oracle->objquadlen = 0;
      }
      else if( newlen < nquadelems )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->objquadelems, nquadelems, newlen) );
         oracle->objquadlen = newlen;
      }

      /* update variable degrees for variables in quadratic part to keep them up to date */
      if( oracle->vardegreesuptodate )
         for( j = 0; j < oracle->objquadlen; ++j )
         {
            oracle->vardegrees[oracle->objquadelems[j].idx1] = MAX(2, oracle->vardegrees[oracle->objquadelems[j].idx1]);
            oracle->vardegrees[oracle->objquadelems[j].idx2] = MAX(2, oracle->vardegrees[oracle->objquadelems[j].idx2]);
         }
   }
   
   if( exprtree != NULL )
   {
      int j;
      
      assert(oracle->exprinterpreter != NULL);
      
      SCIP_CALL( SCIPexprtreeCopy(oracle->blkmem, &oracle->objexprtree, (SCIP_EXPRTREE*)exprtree) );
      
      SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->objexprtree) );
      
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objexprvaridxs, exprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree)) );
      if( oracle->vardegreesuptodate )
         for (j = 0; j < SCIPexprtreeGetNVars(oracle->objexprtree); ++j)
         {
            assert(exprvaridxs[j] >= 0);
            assert(exprvaridxs[j] <  oracle->nvars);
            oracle->vardegrees[exprvaridxs[j]] = INT_MAX; /* @TODO could try to be more clever, maybe use getMaxDegree function in exprtree */
         }
   }
   
   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to change bounds */
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
SCIP_RETCODE SCIPnlpiOracleChgConsSides(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to change bounds */
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
   int c;
   int lastgood; /* index of the last variable that should be kept */ 
   
   assert(oracle != NULL);

   invalidateJacobiSparsity(oracle);
   invalidateHessianLagSparsity(oracle);
   
   lastgood = oracle->nvars - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1 )
      --lastgood;
   if( lastgood < 0 )
   { /* all variables should be deleted */
      assert(oracle->nconss == 0); /* we could relax this by checking that all constraints are constant */
      assert(oracle->objquadlen == 0);
      assert(oracle->objexprtree == NULL || SCIPexprtreeGetNVars(oracle->objexprtree) == 0);
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

   /* update indices in linear part, sort indices, and then clear elements that are marked as deleted */
   mapIndices(delstats, oracle->objnlin, oracle->objlinidxs);
   SCIPsortIntReal(oracle->objlinidxs, oracle->objlinvals, oracle->objnlin);
   SCIP_CALL( clearDeletedLinearElements(oracle->blkmem, &oracle->objlinidxs, &oracle->objlinvals, &oracle->objnlin) );
   
   /* update indices in quadratic part, sort elements, and then clear elements that are marked as deleted */
   mapIndicesQuad(delstats, oracle->objquadlen, oracle->objquadelems);
   SCIPquadelemSort(oracle->objquadelems, oracle->objquadlen);
   SCIP_CALL( clearDeletedQuadElements(oracle->blkmem, &oracle->objquadelems, &oracle->objquadlen) );

   if( oracle->objexprtree )
   {
      mapIndices(delstats, SCIPexprtreeGetNVars(oracle->objexprtree), oracle->objexprvaridxs);
      /* assert that all variables from this expression have been deleted */
      assert(SCIPexprtreeGetNVars(oracle->objexprtree) == 0 || oracle->objexprvaridxs[SCIPexprtreeGetNVars(oracle->objexprtree)-1] == -1);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->objexprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree));
      SCIP_CALL( SCIPexprtreeFree(&oracle->objexprtree) );
   }
   
   for( c = 0; c < oracle->nconss; ++c )
   {
      if( oracle->conslinlens )
      {
         mapIndices(delstats, oracle->conslinlens[c], oracle->conslinidxs[c]);
         SCIPsortIntReal(oracle->conslinidxs[c], oracle->conslincoefs[c], oracle->conslinlens[c]);
         SCIP_CALL( clearDeletedLinearElements(oracle->blkmem, &oracle->conslinidxs[c], &oracle->conslincoefs[c], &oracle->conslinlens[c]) );
      }
      
      if( oracle->consquadlens )
      {
         mapIndicesQuad(delstats, oracle->consquadlens[c], oracle->consquadelems[c]);
         SCIPquadelemSort(oracle->consquadelems[c], oracle->consquadlens[c]);
         SCIP_CALL( clearDeletedQuadElements(oracle->blkmem, &oracle->consquadelems[c], &oracle->consquadlens[c]) );
      }
      
      if( oracle->consexprtrees )
      {
         mapIndices(delstats, SCIPexprtreeGetNVars(oracle->consexprtrees[c]), oracle->consexprvaridxs[c]);
         /* assert that all variables from this expression have been deleted */
         assert(SCIPexprtreeGetNVars(oracle->consexprtrees[c]) == 0 || oracle->consexprvaridxs[c][SCIPexprtreeGetNVars(oracle->consexprtrees[c])-1] == -1);
         BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->consexprvaridxs[c], SCIPexprtreeGetNVars(oracle->consexprtrees[c]));
         SCIP_CALL( SCIPexprtreeFree(&oracle->consexprtrees[c]) );
      }
   }
   
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varlbs, oracle->nvars, lastgood+1) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varubs, oracle->nvars, lastgood+1) );
   
   if( oracle->varnames != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varnames, oracle->nvars, lastgood+1) );
   }

   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->vardegrees, oracle->nvars, lastgood+1) );

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
   oracle->vardegreesuptodate = FALSE;
   
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
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslhss, oracle->nconss, lastgood+1) );
   SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consrhss, oracle->nconss, lastgood+1) );

   if( oracle->conslinlens != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss, lastgood+1) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss, lastgood+1) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss, lastgood+1) );
   }

   if( oracle->consquadlens != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadlens,  oracle->nconss, lastgood+1) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consquadelems, oracle->nconss, lastgood+1) );
   }
   
   if( oracle->consexprvaridxs != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs, oracle->nconss, lastgood+1) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consexprtrees,   oracle->nconss, lastgood+1) );
   }
   
   if( oracle->consnames != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->consnames, oracle->nconss, lastgood+1) );
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
   
   SCIPdebugMessage("change %d linear coefficients in cons %d\n", nentries, considx);

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
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinlens,  oracle->nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslinidxs,  oracle->nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->conslincoefs, oracle->nconss) );
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
      
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, linidxs,  varidxs,  nentries) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, lincoefs, newcoefs, nentries) );
      SCIPsortIntReal(*linidxs, *lincoefs, nentries);

      addednew = TRUE;
   }
   else
   {
      int pos;
      int len;

      len = *linlen;
      
      for( i = 0; i < nentries; ++i )
      {
         /* we assume that indices are not repeating in varidxs (!) */
         if( SCIPsortedvecFindInt(*linidxs, varidxs[i], *linlen, &pos) )  /*lint !e613*/
         {
            (*lincoefs)[pos] = newcoefs[i];  /*lint !e613*/
            SCIPdebugMessage("replace coefficient of var %d at pos %d by %g\n", varidxs[i], pos, newcoefs[i]);
         }
         else if( newcoefs[i] != 0.0 )  /*lint !e613*/
         {
            if( !addednew )
            { /* first coefficient that is added new, realloc memory */
               int newsize = *linlen + nentries; /* new size for arrays (upper bound) */
               SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, linidxs,  *linlen, newsize) );
               SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, lincoefs, *linlen, newsize) );
            }
            else
               assert(len < *linlen);
            /* append new entries */
            (*linidxs)[len]  = varidxs[i];  /*lint !e613*/
            (*lincoefs)[len] = newcoefs[i]; /*lint !e613*/
            ++len;
            addednew = TRUE;
            /* increase degree of variable to 1 */
            oracle->vardegrees[varidxs[i]] = MAX(1, oracle->vardegrees[varidxs[i]]);  /*lint !e613*/
            SCIPdebugMessage("added coefficient of var %d at pos %d, value %g\n", varidxs[i], len-1, newcoefs[i]);
         }
      }
      assert(*linlen == len || addednew);

      if( addednew )
      {
         /* shrink to actual needed size */
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, linidxs , *linlen + nentries, len) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, lincoefs, *linlen + nentries, len) );
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
   int                   nquadelems,         /**< number of entries in quadratic constraint to change */
   const SCIP_QUADELEM*  quadelems           /**< new elements in quadratic matrix (replacing already existing ones or adding new ones) */
   )
{  /*lint --e{715}*/
   SCIP_Bool addednew;
   int       i;
   
   SCIP_QUADELEM** myquadelems;
   int*            myquadlen;
   
   assert(oracle != NULL);
   assert(quadelems != NULL || nquadelems == 0);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   
   if( nquadelems == 0 )
      return SCIP_OKAY;
   
   addednew = FALSE;
   
   if( considx == -1 )
   {
      myquadelems = &oracle->objquadelems;
      myquadlen   = &oracle->objquadlen;
   }
   else
   {
      if( oracle->consquadlens == NULL )
      { /* first time we have quadratic coefficients in a constraint */
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadlens),  oracle->nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &(oracle->consquadelems), oracle->nconss) );
         BMSclearMemoryArray(oracle->consquadlens,  oracle->nconss);
         BMSclearMemoryArray(oracle->consquadelems, oracle->nconss);
      }
      
      myquadelems = &oracle->consquadelems[considx];
      myquadlen   = &oracle->consquadlens [considx];
   }
   
   if( *myquadlen == 0 )
   { /* first time we have quadratic coefficients in this constraint (or objective) */
      assert(*myquadelems == NULL);
      
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, myquadelems, quadelems, nquadelems) );
            
      addednew = TRUE;
   }
   else
   {
      int pos;
      int len;

      len = *myquadlen;
      
      for( i = 0; i < nquadelems; ++i )
      {
         assert(quadelems[i].idx1 >= 0);  /*lint !e613*/
         assert(quadelems[i].idx2 >= 0);  /*lint !e613*/
         assert(quadelems[i].idx1 < oracle->nvars);  /*lint !e613*/
         assert(quadelems[i].idx2 < oracle->nvars);  /*lint !e613*/
         
         /* if we already have an entry for quadelems[i], then just replace the coefficient, otherwise append new entry */
         if( SCIPquadelemSortedFind(*myquadelems, quadelems[i].idx1, quadelems[i].idx2, *myquadlen, &pos) )  /*lint !e613*/
         {
            /* we can assume that index pairs are not repeating, since we squeezed them when creating the array at first time and changing coefficients should not create duplicate entries */
            assert(!SCIPquadelemSortedFind(&(*myquadelems)[pos+1], quadelems[i].idx1, quadelems[i].idx2, *myquadlen-pos-1, NULL));  /*lint !e613*/
            
            (*myquadelems)[pos].coef = quadelems[i].coef;  /*lint !e613*/
         }
         else
         {
            if( !addednew )
            { /* first coefficient that is added new, realloc memory */
               int newsize = *myquadlen + nquadelems; /* new size for arrays (upper bound) */
               SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, myquadelems, *myquadlen, newsize) );
            }
            /* append new entry */
            (*myquadelems)[len] = quadelems[i];  /*lint !e613*/
            ++len;
            addednew = TRUE;
            /* increase degree of variables to 2 */
            oracle->vardegrees[quadelems[i].idx1] = MAX(2, oracle->vardegrees[quadelems[i].idx1]);  /*lint !e613*/
            oracle->vardegrees[quadelems[i].idx2] = MAX(2, oracle->vardegrees[quadelems[i].idx2]);  /*lint !e613*/
         }
      }
      
      if( addednew )
      {
         /* shrink to actual needed size */
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, myquadelems, *myquadlen + nquadelems, len) );
         *myquadlen = len;
      }
   }
   
   if( addednew )
   {
      invalidateJacobiSparsity(oracle);
      invalidateHessianLagSparsity(oracle);
      SCIPquadelemSort(*myquadelems, *myquadlen);
   }
   
   return SCIP_OKAY;
}

/** replaces expression tree of one constraint or objective
 */
SCIP_RETCODE SCIPnlpiOracleChgExprtree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where expression tree should be changed, or -1 for objective */
   const int*            exprvaridxs,        /**< problem indices of variables in expression tree */
   const SCIP_EXPRTREE*  exprtree            /**< new expression tree, or NULL */
   )
{
   int j;

   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   assert((exprvaridxs != NULL) == (exprtree != NULL));

   invalidateHessianLagSparsity(oracle);
   oracle->vardegreesuptodate = FALSE;
   
   /* free previous expression tree */
   if( considx == -1 )
   {
      if( oracle->objexprtree != NULL )
      {
         assert(oracle->objexprvaridxs != NULL);

         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->objexprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree));
         SCIP_CALL( SCIPexprtreeFree(&oracle->objexprtree) );
      }
      else if( exprtree == NULL )
      {
         /* nothing to do */
         return SCIP_OKAY;
      }
   }
   else
   {
      assert(considx >= 0);
      assert(considx <  oracle->nconss);
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[considx] != NULL )
      {
         assert(oracle->consexprvaridxs != NULL);
         assert(oracle->consexprvaridxs[considx] != NULL);

         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs[considx], SCIPexprtreeGetNVars(oracle->consexprtrees[considx]));
         SCIP_CALL( SCIPexprtreeFree(&oracle->consexprtrees[considx]) );
      }
   }

   invalidateJacobiSparsity(oracle);

   /* if user did not want to set new tree, then we are done */
   if( exprtree == NULL )
      return SCIP_OKAY;

   assert(oracle->exprinterpreter != NULL);

   /* install new expression tree */
   if( considx >= 0 )
   {
      if( oracle->consexprtrees == NULL )
      {
         assert(oracle->consexprvaridxs == NULL);

         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->consexprtrees,   oracle->nconss) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs, oracle->nconss) );
         BMSclearMemoryArray(oracle->consexprtrees,   oracle->nconss);
         BMSclearMemoryArray(oracle->consexprvaridxs, oracle->nconss);
      }

      SCIP_CALL( SCIPexprtreeCopy(oracle->blkmem, &oracle->consexprtrees[considx], (SCIP_EXPRTREE*)exprtree) );
      SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->consexprtrees[considx]) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->consexprvaridxs[considx], exprvaridxs, SCIPexprtreeGetNVars((SCIP_EXPRTREE*)exprtree)) );
   }
   else
   {
      SCIP_CALL( SCIPexprtreeCopy(oracle->blkmem, &oracle->objexprtree, (SCIP_EXPRTREE*)exprtree) );
      SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->objexprtree) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->objexprvaridxs, exprvaridxs, SCIPexprtreeGetNVars(oracle->objexprtree)) );
   }

   /* increase variable degree to keep them up to date */
   if( oracle->vardegreesuptodate )
      for( j = 0; j < SCIPexprtreeGetNVars((SCIP_EXPRTREE*)exprtree); ++j )
      {
         assert(exprvaridxs[j] >= 0);
         assert(exprvaridxs[j] <  oracle->nvars);
         oracle->vardegrees[exprvaridxs[j]] = INT_MAX; /* @TODO could try to be more clever, maybe use getMaxDegree function in exprtree */
      }

   return SCIP_OKAY;
}

/** changes one parameter of expression tree of one constraint or objective
 */
SCIP_RETCODE SCIPnlpiOracleChgExprParam(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where parameter should be changed in expression tree, or -1 for objective */
   int                   paramidx,           /**< index of parameter */
   SCIP_Real             paramval            /**< new value of parameter */
   )
{
   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   assert(paramidx >= 0);
   assert(considx >= 0  || oracle->objexprtree != NULL);
   assert(considx >= 0  || paramidx < SCIPexprtreeGetNParams(oracle->objexprtree));
   assert(considx == -1 || oracle->consexprtrees != NULL);
   assert(considx == -1 || oracle->consexprtrees[considx] != NULL);
   assert(considx == -1 || paramidx < SCIPexprtreeGetNParams(oracle->consexprtrees[considx]));

   SCIPexprtreeSetParamVal(considx >= 0 ? oracle->consexprtrees[considx] : oracle->objexprtree, paramidx, paramval);

   return SCIP_OKAY;
}

/** changes the constant value in the objective function
 */
SCIP_RETCODE SCIPnlpiOracleChgObjConstant(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             objconstant         /**< new value for objective constant */
   )
{
   assert(oracle != NULL);

   oracle->objconstant = objconstant;
   
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

/** gives the variables names, or NULL if not set */
char** SCIPnlpiOracleGetVarNames(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->varnames;
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
   
   updateVariableDegrees(oracle);
   
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

   updateVariableDegrees(oracle);

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

/** gives the constraint names, or NULL if not set */
char** SCIPnlpiOracleGetConstraintNames(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   
   return oracle->consnames;
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

/** Gives maximum degree over all constraints and the objective (or over all variables, resp.).
 * Thus, if this function returns 0, then the objective and all constraints are constant.
 * If it returns 1, then the problem in linear.
 * If it returns 2, then its a QP, QCP, or QCQP.
 * And if it returns > 2, then it is an NLP. */
int SCIPnlpiOracleGetMaxDegree(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   int i;
   int maxdegree;
   
   assert(oracle != NULL);

   updateVariableDegrees(oracle);

   maxdegree = 0;
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->vardegrees[i] > maxdegree )
         maxdegree = oracle->vardegrees[i];
   }
   
   return maxdegree;
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
      oracle->objquadlen, oracle->objquadelems,
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
      (oracle->conslinlens   != NULL ? oracle->conslinlens    [considx] : 0),
      (oracle->conslinlens   != NULL ? oracle->conslinidxs    [considx] : NULL),
      (oracle->conslinlens   != NULL ? oracle->conslincoefs   [considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadlens   [considx] : 0),
      (oracle->consquadlens  != NULL ? oracle->consquadelems  [considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprtrees  [considx] : NULL),
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
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Real*            objgrad             /**< pointer to store (dense) objective gradient */  
   )
{
   assert(oracle != NULL);
   
   SCIP_CALL( evalFunctionGradient(oracle,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadelems,
      oracle->objexprvaridxs, oracle->objexprtree,
      x, isnewx, objval, objgrad) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** computes a constraints gradient in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */  
   )
{
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( evalFunctionGradient(oracle,
      (oracle->conslinlens   != NULL ? oracle->conslinlens    [considx] : 0),
      (oracle->conslinlens   != NULL ? oracle->conslinidxs    [considx] : NULL),
      (oracle->conslinlens   != NULL ? oracle->conslincoefs   [considx] : NULL),
      (oracle->consquadlens  != NULL ? oracle->consquadlens   [considx] : 0),
      (oracle->consquadlens  != NULL ? oracle->consquadelems  [considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[considx] : NULL),
      (oracle->consexprtrees != NULL ? oracle->consexprtrees  [considx] : NULL),
      x, isnewx, conval, congrad) );
   
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

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->jacoffsets, oracle->nconss + 1) );

   maxnnz = MIN(oracle->nvars, 10) * oracle->nconss;  /* initial guess */
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, maxnnz) );
   
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
   
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &nzflag, oracle->nvars) );

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
               SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, oldsize, maxnnz) );
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
         assert(oracle->consquadelems    != NULL );
         assert(oracle->consquadelems[i] != NULL );
         for( j = 0; j < oracle->consquadlens[i]; ++j )
         {
            nzflag[oracle->consquadelems[i][j].idx1] = TRUE;
            nzflag[oracle->consquadelems[i][j].idx2] = TRUE;
         }
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         assert(oracle->consexprvaridxs != NULL );
         assert(oracle->consexprvaridxs[i] != NULL );
         for (j = 0; j < SCIPexprtreeGetNVars(oracle->consexprtrees[i]); ++j)
            nzflag[oracle->consexprvaridxs[i][j]] = TRUE;
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
            SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, oldsize, maxnnz) );
         }
         
         oracle->jaccols[nnz] = j;
         ++nnz;
      }
   }
   
   oracle->jacoffsets[oracle->nconss] = nnz;
   
   /* shrink jaccols array to nnz */
   if( nnz < maxnnz )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, maxnnz, nnz) );
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
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
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
   
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &grad, oracle->nvars) );
   
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
         SCIP_RETCODE retcode = SCIPnlpiOracleEvalConstraintGradient(oracle, i, x, isnewx, (convals ? &convals[i] : &dummy), grad);
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
   int   nnz;
   int   i;
   int   j;
   int   cnt;
   
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

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->heslagoffsets, oracle->nvars + 1) );
   
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &colnz,  oracle->nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &collen, oracle->nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &colnnz, oracle->nvars) );
   BMSclearMemoryArray(colnz,  oracle->nvars);
   BMSclearMemoryArray(collen, oracle->nvars);
   BMSclearMemoryArray(colnnz, oracle->nvars);
   nnz = 0;
   
   if( oracle->objquadlen != 0 )
   {
      SCIP_CALL( hessLagSparsitySetNzFlagForQuad(oracle, colnz, collen, colnnz, &nnz, oracle->objquadlen, oracle->objquadelems) );
   }

   if( oracle->objexprtree != NULL )
   {
      SCIP_CALL( hessLagSparsitySetNzFlagForExprtree(oracle, colnz, collen, colnnz, &nnz, oracle->objexprvaridxs, oracle->objexprtree, oracle->nvars) );
   }

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
      {
         SCIP_CALL( hessLagSparsitySetNzFlagForQuad(oracle, colnz, collen, colnnz, &nnz, oracle->consquadlens[i], oracle->consquadelems[i]) );
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         SCIP_CALL( hessLagSparsitySetNzFlagForExprtree(oracle, colnz, collen, colnnz, &nnz, oracle->consexprvaridxs[i], oracle->consexprtrees[i], oracle->nvars) );
      }
   }
   
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->heslagcols, nnz) );
   
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
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
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
         assert(oracle->objquadelems != NULL);
         SCIP_CALL( hessLagAddQuad(objfactor, oracle->objquadlen, oracle->objquadelems, oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
      
      if( oracle->objexprtree != NULL )
      {
         assert(oracle->objexprvaridxs != NULL );
         SCIP_CALL( hessLagAddExprtree(oracle, objfactor, x, isnewx, oracle->objexprvaridxs, oracle->objexprtree, oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
   }
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( lambda[i] == 0.0 )
         continue;
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
      {
         assert(oracle->consquadelems    != NULL);
         assert(oracle->consquadelems[i] != NULL);
         SCIP_CALL( hessLagAddQuad(lambda[i], oracle->consquadlens[i], oracle->consquadelems[i], oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         assert(oracle->consexprvaridxs != NULL);
         assert(oracle->consexprvaridxs[i] != NULL);
         SCIP_CALL( hessLagAddExprtree(oracle, lambda[i], x, isnewx, oracle->consexprvaridxs[i], oracle->consexprtrees[i], oracle->heslagoffsets, oracle->heslagcols, hessian) );
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
   
   SCIPmessageFPrintInfo(file, "NLPI Oracle %s: %d variables and %d constraints\n", oracle->name ? oracle->name : "", oracle->nvars, oracle->nconss);
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%10s", oracle->varnames[i]);
      else
         SCIPmessageFPrintInfo(file, "x%09d", i);
      SCIPmessageFPrintInfo(file, ": [%8g, %8g]", oracle->varlbs[i], oracle->varubs[i]);
      if( oracle->vardegreesuptodate )
         SCIPmessageFPrintInfo(file, "\t degree: %d", oracle->vardegrees[i]);
      SCIPmessageFPrintInfo(file, "\n");
   }
   
   SCIPmessageFPrintInfo(file, "objective: ");
   SCIP_CALL( printFunction(oracle, file,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadelems,
      oracle->objexprtree, oracle->objexprvaridxs,
      FALSE, FALSE) );
   if( oracle->objconstant != 0.0 )
      SCIPmessageFPrintInfo(file, "%+.20g", oracle->objconstant);
   SCIPmessageFPrintInfo(file, "\n");
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
         SCIPmessageFPrintInfo(file, "%10s", oracle->consnames[i]);
      else
         SCIPmessageFPrintInfo(file, "con%07d", i);
      
      SCIPmessageFPrintInfo(file, ": ");
      if( oracle->conslhss[i] > -oracle->infinity && oracle->consrhss[i] < oracle->infinity && oracle->conslhss[i] != oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, "%.20g <= ", oracle->conslhss[i]);
      
      SCIP_CALL( printFunction(oracle, file,
         (oracle->conslinlens   != NULL ? oracle->conslinlens  [i] : 0),
         (oracle->conslinlens   != NULL ? oracle->conslinidxs  [i] : NULL),
         (oracle->conslinlens   != NULL ? oracle->conslincoefs [i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadlens [i] : 0),
         (oracle->consquadlens  != NULL ? oracle->consquadelems[i] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[i] : NULL),
         FALSE, FALSE) );
      
      if( oracle->conslhss[i] == oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, " = %.20g", oracle->consrhss[i]);
      else if( oracle->consrhss[i] <  oracle->infinity )
         SCIPmessageFPrintInfo(file, " <= %.20g", oracle->consrhss[i]);
      else if( oracle->conslhss[i] > -oracle->infinity )
         SCIPmessageFPrintInfo(file, " >= %.20g", oracle->conslhss[i]);
      
      SCIPmessageFPrintInfo(file, "\n");
   }

   return SCIP_OKAY;
}

/** prints the problem to a file in GAMS format
 * If there are variable (equation, resp.) names with more than 9 characters, then variable (equation, resp.) names are prefixed with an unique identifier.
 * This is to make it easier to identify variables solution output in the listing file.
 * Names with more than 64 characters are shorten to 64 letters due to GAMS limits.
 */
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;
   int nllevel; /* level of nonlinearity of problem: linear = 0, quadratic, smooth nonlinear, nonsmooth */
   static const char* nllevelname[4] = { "LP", "QCP", "NLP", "DNLP" };
   const char* problemname;
   char namebuf[70];
   SCIP_Bool havelongvarnames;
   SCIP_Bool havelongequnames;

   assert(oracle != NULL);
   
   if( file == NULL )
      file = stdout;
   
   nllevel = 0;

   havelongvarnames = FALSE;
   for( i = 0; i < oracle->nvars; ++i )
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL && strlen(oracle->varnames[i]) > 9 )
      {
         havelongvarnames = TRUE;
         break;
      }

   havelongequnames = FALSE;
   for( i = 0; i < oracle->nconss; ++i )
      if( oracle->consnames != NULL && oracle->consnames[i] != NULL && strlen(oracle->consnames[i]) > 9 )
      {
         havelongequnames = TRUE;
         break;
      }

   SCIPmessageFPrintInfo(file, "$offlisting\n");
   SCIPmessageFPrintInfo(file, "$offdigit\n");
   SCIPmessageFPrintInfo(file, "* NLPI Oracle Problem %s\n", oracle->name ? oracle->name : "");
   SCIPmessageFPrintInfo(file, "Variables ");
   for( i = 0; i < oracle->nvars; ++i )
   {
      printName(namebuf, oracle->varnames, i, 'x', NULL, havelongvarnames);
      SCIPmessageFPrintInfo(file, "%s, ", namebuf);
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(file, "\n");
   }
   SCIPmessageFPrintInfo(file, "NLPIORACLEOBJVAR;\n\n");
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varlbs[i] == oracle->varubs[i] )
      {
         printName(namebuf, oracle->varnames, i, 'x', NULL, havelongvarnames);
         SCIPmessageFPrintInfo(file, "%s.fx = %.20g;\t", namebuf, oracle->varlbs[i]);
      }
      else
      {
         if( oracle->varlbs[i] > -oracle->infinity )
         {
            printName(namebuf, oracle->varnames, i, 'x', NULL, havelongvarnames);
            SCIPmessageFPrintInfo(file, "%s.lo = %.20g;\t", namebuf, oracle->varlbs[i]);
         }
         if( oracle->varubs[i] <  oracle->infinity )
         {
            printName(namebuf, oracle->varnames, i, 'x', NULL, havelongvarnames);
            SCIPmessageFPrintInfo(file, "%s.up = %.20g;\t", namebuf, oracle->varubs[i]);
         }
      }
      if( initval != NULL )
      {
         printName(namebuf, oracle->varnames, i, 'x', NULL, havelongvarnames);
         SCIPmessageFPrintInfo(file, "%s.l = %.20g;\t", namebuf, initval[i]);
      }
      SCIPmessageFPrintInfo(file, "\n");
   }
   SCIPmessageFPrintInfo(file, "\n");
   
   SCIPmessageFPrintInfo(file, "Equations ");
   for( i = 0; i < oracle->nconss; ++i )
   {
      printName(namebuf, oracle->consnames, i, 'e', NULL, havelongequnames);
      SCIPmessageFPrintInfo(file, "%s, ", namebuf);

      if( oracle->conslhss[i] > -oracle->infinity && oracle->consrhss[i] < oracle->infinity && oracle->conslhss[i] != oracle->consrhss[i] )
      { /* ranged row: add second constraint */
         printName(namebuf, oracle->consnames, i, 'e', "_RNG", havelongequnames);
         SCIPmessageFPrintInfo(file, "%s, ", namebuf);
      }
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(file, "\n");
   }
   SCIPmessageFPrintInfo(file, "NLPIORACLEOBJ;\n\n");
   
   SCIPmessageFPrintInfo(file, "NLPIORACLEOBJ.. NLPIORACLEOBJVAR =E= ");
   SCIP_CALL( printFunction(oracle, file,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadelems,
      oracle->objexprtree, oracle->objexprvaridxs,
      havelongvarnames, havelongequnames) );
   if( oracle->objconstant != 0.0 )
      SCIPmessageFPrintInfo(file, "%+.20g", oracle->objconstant);
   SCIPmessageFPrintInfo(file, ";\n");
   
   for( i = 0; i < oracle->nconss; ++i )
   {
      printName(namebuf, oracle->consnames, i, 'e', NULL, havelongequnames);
      SCIPmessageFPrintInfo(file, "%s.. ", namebuf);
      
      SCIP_CALL( printFunction(oracle, file,
         (oracle->conslinlens   != NULL ? oracle->conslinlens  [i] : 0),
         (oracle->conslinlens   != NULL ? oracle->conslinidxs  [i] : NULL),
         (oracle->conslinlens   != NULL ? oracle->conslincoefs [i] : NULL),
         (oracle->consquadlens  != NULL ? oracle->consquadlens [i] : 0),
         (oracle->consquadlens  != NULL ? oracle->consquadelems[i] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[i] : NULL),
         havelongvarnames, havelongequnames) );

      if( oracle->conslhss[i] == oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, " =E= %.20g", oracle->consrhss[i]);
      else if( oracle->consrhss[i] <  oracle->infinity )
         SCIPmessageFPrintInfo(file, " =L= %.20g", oracle->consrhss[i]);
      else if( oracle->conslhss[i] > -oracle->infinity )
         SCIPmessageFPrintInfo(file, " =G= %.20g", oracle->conslhss[i]);
      else
         SCIPmessageFPrintInfo(file, " =N= 0");
      SCIPmessageFPrintInfo(file, ";\n");

      if( oracle->conslhss[i] > -oracle->infinity && oracle->consrhss[i] < oracle->infinity && oracle->conslhss[i] != oracle->consrhss[i] )
      {
         printName(namebuf, oracle->consnames, i, 'e', "_RNG", havelongequnames);
         SCIPmessageFPrintInfo(file, "%s.. ", namebuf);
        
         SCIP_CALL( printFunction(oracle, file,
            (oracle->conslinlens   != NULL ? oracle->conslinlens  [i] : 0),
            (oracle->conslinlens   != NULL ? oracle->conslinidxs  [i] : NULL),
            (oracle->conslinlens   != NULL ? oracle->conslincoefs [i] : NULL),
            (oracle->consquadlens  != NULL ? oracle->consquadlens [i] : 0),
            (oracle->consquadlens  != NULL ? oracle->consquadelems[i] : NULL),
            (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL),
            (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[i] : NULL),
            havelongvarnames, havelongequnames) );
         
         SCIPmessageFPrintInfo(file, " =G= %.20g;\n", oracle->conslhss[i]);
      }
      
      if( nllevel <= 0 && oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
         nllevel = 1;
      if( nllevel <= 1 && oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
         nllevel = 2;
      if( nllevel <= 2 && oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL && exprIsNonSmooth(SCIPexprtreeGetRoot(oracle->consexprtrees[i])) )
         nllevel = 3;
   }
   
   problemname = oracle->name ? oracle->name : "m";
   
   SCIPmessageFPrintInfo(file, "Model %s / all /;\n", problemname);
   SCIPmessageFPrintInfo(file, "option limrow = 0;\n");
   SCIPmessageFPrintInfo(file, "option limcol = 0;\n");
   SCIPmessageFPrintInfo(file, "Solve %s minimizing NLPIORACLEOBJVAR using %s;\n", problemname, nllevelname[nllevel]);

   return SCIP_OKAY;
}

/**@} */
