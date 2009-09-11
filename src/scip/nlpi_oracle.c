/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpi_oracle.c,v 1.17 2009/09/11 16:02:48 bzfwinkm Exp $"

/**@file    nlpi_oracle.c
 * @ingroup NLPIS
 * @brief   implementation of NLPI oracle interface
 * @author  Stefan Vigerske
 * 
 * @todo jacobi evaluation should be sparse
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_oracle.h"

#include <string.h> /* for strlen */

struct SCIP_NlpiOracle
{
   int                   nvars;              /**< number of variables */
   SCIP_Real*            varlbs;             /**< array with variable lower bounds */
   SCIP_Real*            varubs;             /**< array with variable upper bounds */
   char**                varnames;           /**< array with variable names */
   
   int                   nconss;             /**< number of constraints */
   SCIP_Real*            conslhss;           /**< array with left-hand sides of constraints */
   SCIP_Real*            consrhss;           /**< array with right-hand sides of constraints */
   int*                  conslinoffsets;     /**< array with start index of each constraints linear coefficients in lininds and linvals 
                                              *   (length: nconss + 1, linoffsets[nconss] gives length of lininds and linvals) */
   int*                  conslininds;        /**< array with variable indices in linear part */
   SCIP_Real*            conslinvals;        /**< array with variable coefficients in linear part */ 
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
   
   int*                  heslagoffsets;      /**< rowwise sparsity pattern of hessian matrix of Lagrangian: row offsets in laghescol */
   int*                  heslagcols;         /**< rowwise sparsity pattern of hessian matrix of Lagrangian: column indices; sorted for each row */
   
   int*                  vardegrees;         /**< array with maximal degree of variable over objective and all constraints */
};

/** Invalidates the sparsity pattern of the Jacobian.
 *  Should be called when constraints are added or deleted.
 */
static
void SCIPnlpiOracleInvalidateJacobiSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(scip != NULL);
   assert(oracle != NULL);
   
   if( oracle->jacoffsets == NULL )
   { /* nothing to do */
      assert(oracle->jaccols == NULL);
      return;
   }
   
   assert(oracle->jaccols != NULL);
   SCIPfreeMemoryArray(scip, &oracle->jacoffsets);
   SCIPfreeMemoryArray(scip, &oracle->jaccols);
}

/** Invalidates the sparsity pattern of the Hessian of the Lagragian.
 *  Should be called when the objective is set or constraints are added or deleted.
 */
static
void SCIPnlpiOracleInvalidateHessianLagSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( oracle->heslagoffsets == NULL )
   { /* nothing to do */
      assert(oracle->heslagcols == NULL);
      return;
   }
   
   assert(oracle->heslagcols != NULL);
   SCIPfreeMemoryArray(scip, &oracle->heslagoffsets);
   SCIPfreeMemoryArray(scip, &oracle->heslagcols);
}

/** creates an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleCreate(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(scip != NULL);
   assert(oracle != NULL);
   
   SCIP_CALL( SCIPallocMemory(scip, oracle) );
   
   (*oracle)->nvars = 0;
   (*oracle)->varlbs = NULL;
   (*oracle)->varubs = NULL;
   (*oracle)->varnames = NULL;

   (*oracle)->nconss = 0;
   (*oracle)->conslhss = NULL;
   (*oracle)->consrhss = NULL;
   (*oracle)->conslinoffsets = NULL;
   (*oracle)->conslininds = NULL;
   (*oracle)->conslinvals = NULL;
   (*oracle)->consquadlens = NULL;
   (*oracle)->consquadrows = NULL;
   (*oracle)->consquadcols = NULL;
   (*oracle)->consquadvals = NULL;
   (*oracle)->consexprvaridxs = NULL;
   (*oracle)->consexprtrees = NULL;
   (*oracle)->consnames = NULL;

   (*oracle)->objconstant = 0.0;
   (*oracle)->objnlin = 0;
   (*oracle)->objlinidxs = NULL;
   (*oracle)->objlinvals = NULL;
   (*oracle)->objquadlen = 0;
   (*oracle)->objquadrows = NULL;
   (*oracle)->objquadcols = NULL;
   (*oracle)->objquadvals = NULL;
   (*oracle)->objexprvaridxs = NULL;
   (*oracle)->objexprtree = NULL;
   
   (*oracle)->jacoffsets = NULL;
   (*oracle)->jaccols = NULL;
   
   (*oracle)->heslagoffsets = NULL;
   (*oracle)->heslagcols = NULL;
   
   (*oracle)->vardegrees = NULL;
   
   return SCIP_OKAY;
}

/** initializes NLPI oracle */
SCIP_RETCODE SCIPnlpiOracleInit(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< NLPIORACLE data structure */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(oracle != NULL);

   return SCIP_OKAY;
}

/** frees an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(oracle != NULL);
   assert(*oracle != NULL);
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->varlbs);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->varubs);
   if( (*oracle)->varnames != NULL )
   {
      for( i = 0; i < (*oracle)->nvars; ++i )
      {
        SCIPfreeMemoryArrayNull(scip, &(*oracle)->varnames[i]);
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->varnames);
   }
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conslhss);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->consrhss);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conslinoffsets);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conslininds);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conslinvals);
   
   if( (*oracle)->consquadlens != NULL )
   {
      assert((*oracle)->consquadrows != NULL);
      assert((*oracle)->consquadcols != NULL);
      assert((*oracle)->consquadvals != NULL);
      for( i = 0; i < (*oracle)->nconss; ++i )
      {
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->consquadrows[i]);
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->consquadcols[i]);
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->consquadvals[i]);
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->consquadlens);
      SCIPfreeMemoryArray(scip, &(*oracle)->consquadrows);
      SCIPfreeMemoryArray(scip, &(*oracle)->consquadcols);
      SCIPfreeMemoryArray(scip, &(*oracle)->consquadvals);
   }
   
   if( (*oracle)->consexprtrees != NULL )
   {
      for( i = 0; i < (*oracle)->nconss; ++i )
      {
         if( (*oracle)->consexprtrees[i] != NULL )
         {
            SCIPfreeMemoryArrayNull(scip, &(*oracle)->consexprvaridxs[i]);
         }
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->consexprtrees);
      SCIPfreeMemoryArray(scip, &(*oracle)->consexprvaridxs);
   }
   
   if( (*oracle)->consnames != NULL )
   {
      for( i = 0; i < (*oracle)->nconss; ++i )
      {
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->consnames[i]);
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->consexprtrees);
   }
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objlinidxs);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objlinvals);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objquadrows);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objquadcols);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objquadvals);
   if( (*oracle)->objexprtree != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*oracle)->objexprvaridxs);
   }

   SCIPnlpiOracleInvalidateJacobiSparsity(scip, *oracle);
   SCIPnlpiOracleInvalidateHessianLagSparsity(scip, *oracle);
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->vardegrees);

   SCIPfreeMemory(scip, oracle);
   
   return SCIP_OKAY;
}

/** adds variables */
SCIP_RETCODE SCIPnlpiOracleAddVars(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to add */
   const SCIP_Real*      lbs,                /**< array with lower bounds of new variables, or NULL if all -infinity */
   const SCIP_Real*      ubs,                /**< array with upper bounds of new variables, or NULL if all +infinity */
   const char**          varnames            /**< array with names of new variables, or NULL if no names should be stored */
   )
{
   int i;

   assert(scip != NULL);
   assert(oracle != NULL);
   
   if( nvars == 0 )
      return SCIP_OKAY;
   
   assert(nvars > 0);
   
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->varlbs, oracle->nvars + nvars) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->varubs, oracle->nvars + nvars) );
   
   if( lbs != NULL )
   {
      BMScopyMemoryArray(&((oracle->varlbs)[oracle->nvars]), lbs, nvars);
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varlbs[oracle->nvars+i] = -SCIPinfinity(scip);
   
   if( ubs != NULL )
   {
      BMScopyMemoryArray(&((oracle->varubs)[oracle->nvars]), ubs, nvars);
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varubs[oracle->nvars+i] =  SCIPinfinity(scip);
   
   if( oracle->varnames != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->varnames, oracle->nvars + nvars) );
      if( varnames != NULL )
      {
         for( i = 0; i < nvars; ++i )
         {
            if( varnames[i] != NULL )
            {
               SCIP_CALL( SCIPduplicateMemoryArray(scip, &((oracle->varnames)[oracle->nvars+i]), varnames[i], strlen(varnames[i])+1) );
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
      SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->varnames), oracle->nvars + nvars) );
      BMSclearMemoryArray(oracle->varnames, oracle->nvars);
      for( i = 0; i < nvars; ++i )
      {
         if( varnames[i] != NULL )
         {
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &((oracle->varnames)[oracle->nvars+i]), varnames[i], strlen(varnames[i])+1) );
         }
         else
            oracle->varnames[oracle->nvars+i] = NULL;
      }
   }
   
   SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->vardegrees), oracle->nvars + nvars) );
   BMSclearMemoryArray(&((oracle->vardegrees)[oracle->nvars]), nvars);

   /* @TODO update sparsity pattern by extending heslagoffsets */
   SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);

   oracle->nvars += nvars;

   return SCIP_OKAY;
}

/** adds constraints 
 * 
 *  linear coefficients: row(=constraint) oriented matrix;
 *  quadratic coefficiens: row oriented matrix for each constraint
 */
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            linoffsets,         /**< array with start indices of each constraints linear coefficients in linind and linval 
                                              *   (length: nconss + 1, linoffset[nconss] gives length of linind and linval), or NULL if no linear part */
   const int*            lininds,            /**< array with variable indices in linear part, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< array with variable coefficients in linear part, of NULL if no linear part */ 
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
   int i;
   int oldnnz;   

   assert(scip != NULL);
   assert(oracle != NULL);
   
   if( nconss == 0 )
      return SCIP_OKAY;
   
   assert(nconss > 0);

   addednlcon = FALSE;

   SCIPnlpiOracleInvalidateJacobiSparsity(scip, oracle); /* @TODO we could also update (extend) the sparsity pattern */

   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conslhss, oracle->nconss + nconss) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->consrhss, oracle->nconss + nconss) );

   if( lhss != NULL )
   {
      BMScopyMemoryArray(&((oracle->conslhss)[oracle->nconss]), lhss, nconss);
   }
   else
      for( i = 0; i < nconss; ++i )
         oracle->conslhss[oracle->nconss+i] = -SCIPinfinity(scip);
   
   if( rhss != NULL )
   {
      BMScopyMemoryArray(&((oracle->consrhss)[oracle->nconss]), rhss, nconss);
   }
   else
      for( i = 0; i < nconss; ++i )
         oracle->consrhss[oracle->nconss+i] =  SCIPinfinity(scip);

   /* old number of nonzeros */
   oldnnz = (oracle->nconss != 0 ? oracle->conslinoffsets[oracle->nconss] : 0);

   SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->conslinoffsets), oracle->nconss + nconss + 1) );
   if( linoffsets != NULL )
   {
      assert(lininds != NULL);
      assert(linvals != NULL);
      
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->conslininds), oldnnz + linoffsets[nconss]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->conslinvals), oldnnz + linoffsets[nconss]) );
   
      for( i = 0; i <= nconss; ++i )
         oracle->conslinoffsets[oracle->nconss + i] = oldnnz + linoffsets[i];
      BMScopyMemoryArray(&((oracle->conslininds)[oldnnz]), lininds, linoffsets[nconss]);
      BMScopyMemoryArray(&((oracle->conslinvals)[oldnnz]), linvals, linoffsets[nconss]);
      
      for( i = 0; i < linoffsets[nconss]; ++i )
         oracle->vardegrees[lininds[i]] = MAX(1, oracle->vardegrees[lininds[i]]);
   }
   else if( oracle->conslinoffsets != NULL )
   {
      for( i = 0; i <= nconss; ++i )
         oracle->conslinoffsets[oracle->nconss + i] = oldnnz;
   }
   
   if( quadoffsets != NULL )
   {
      assert(nquadrows != NULL);
      assert(quadrowidxs != NULL);
      assert(quadidxs != NULL);
      assert(quadvals != NULL);
      if( (oracle->consquadlens == NULL) && (oracle->nconss != 0) )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consquadlens), oracle->nconss + nconss) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consquadrows), oracle->nconss + nconss) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consquadcols), oracle->nconss + nconss) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consquadvals), oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consquadlens, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadrows, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadcols, oracle->nconss);
         BMSclearMemoryArray(oracle->consquadvals, oracle->nconss);
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadlens), oracle->nconss + nconss) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadrows), oracle->nconss + nconss) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadcols), oracle->nconss + nconss) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadvals), oracle->nconss + nconss) );
      }
      for( i = 0; i < nconss; ++i )
      {
         if( quadoffsets[i] != NULL && nquadrows[i] != 0 )
         {
            int j;
            int k;

            assert(quadrowidxs[i] != NULL);
            assert(quadidxs[i] != NULL);
            assert(quadvals[i] != NULL);
            oracle->consquadlens[oracle->nconss + i] = quadoffsets[i][nquadrows[i]];
            SCIP_CALL( SCIPallocMemoryArray(scip, &((oracle->consquadrows)[oracle->nconss + i]), quadoffsets[i][nquadrows[i]]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &((oracle->consquadcols)[oracle->nconss + i]), quadoffsets[i][nquadrows[i]]) );
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &((oracle->consquadvals)[oracle->nconss + i]), quadvals[i], quadoffsets[i][nquadrows[i]]) );
            
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
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadlens), oracle->nconss + nconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadrows), oracle->nconss + nconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadcols), oracle->nconss + nconss) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consquadvals), oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consquadlens)[oracle->nconss]), nconss);
      BMSclearMemoryArray(&((oracle->consquadrows)[oracle->nconss]), nconss);
      BMSclearMemoryArray(&((oracle->consquadcols)[oracle->nconss]), nconss);
      BMSclearMemoryArray(&((oracle->consquadvals)[oracle->nconss]), nconss);
   }
   
   if( exprtrees != NULL )
   {
      if( oracle->consexprtrees == NULL && oracle->nconss != 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consexprtrees),   oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consexprtrees,   oracle->nconss);
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consexprvaridxs), oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consexprvaridxs, oracle->nconss);
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consexprtrees),   oracle->nconss + nconss) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consexprvaridxs), oracle->nconss + nconss) );
      }
      for( i = 0; i < nconss; ++i )
         if( exprtrees[i] != NULL )
         {
            SCIPerrorMessage("nonquadratic functions not supported in NLPI yet.\n");
            return SCIP_ERROR;
         }
         else
         {
            oracle->consexprtrees[oracle->nconss + i] = NULL;
            oracle->consexprvaridxs[oracle->nconss + i] = NULL;
         }
   }
   else if( oracle->consexprtrees != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consexprtrees),   oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consexprtrees)[oracle->nconss]), nconss);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consexprvaridxs), oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consexprvaridxs)[oracle->nconss]), nconss);
   }

   if( consnames != NULL )
   {
      if( oracle->consnames == NULL && oracle->nconss != 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(oracle->consnames), oracle->nconss + nconss) );
         BMSclearMemoryArray(oracle->consnames, oracle->nconss);
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consnames), oracle->nconss + nconss) );
      }
      for( i = 0; i < nconss; ++i )
         if( consnames[i] != NULL )
         {
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &((oracle->consnames)[oracle->nconss + i]), consnames[i], strlen(consnames[i])+1) );
         }
         else
            oracle->consnames[oracle->nconss + i] = NULL;
   }
   else if( oracle->consnames != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(oracle->consnames), oracle->nconss + nconss) );
      BMSclearMemoryArray(&((oracle->consnames)[oracle->nconss]), nconss);
   }
   
   oracle->nconss += nconss;

   if( addednlcon == TRUE )
      SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);
   
   return SCIP_OKAY;
}

/** sets or overwrites objective, a minization problem is expected
 * 
 *  May change sparsity pattern.
 */
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP*                 scip,               /**< pointer to SCIP */
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
   assert(scip != NULL);
   assert(oracle != NULL);
   
   /* clear previous objective */
   SCIPfreeMemoryArrayNull(scip, &oracle->objlinidxs);
   SCIPfreeMemoryArrayNull(scip, &oracle->objlinvals);
   SCIPfreeMemoryArrayNull(scip, &oracle->objquadrows);
   SCIPfreeMemoryArrayNull(scip, &oracle->objquadcols);
   SCIPfreeMemoryArrayNull(scip, &oracle->objquadvals);
   if( oracle->objexprtree != NULL )
   {
      SCIPfreeMemoryArray(scip, &oracle->objexprvaridxs);
      /* @TODO this does not clear the vardegrees's */
   }
   
   if( nquadrows != 0 || oracle->objquadlen != 0 || exprtree != NULL || oracle->objexprtree != NULL )
      SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);
   
   oracle->objconstant = constant;
   
   if( nlin != 0 )
   {
      int i;

      assert(lininds != NULL);
      assert(linvals != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->objlinidxs, lininds, nlin) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->objlinvals, linvals, nlin) );
      oracle->objnlin = nlin;
      for( i = 0; i < nlin; ++i )
         oracle->vardegrees[lininds[i]] = MAX(1, oracle->vardegrees[lininds[i]]);
   }
   
   if( nquadrows != 0 && quadoffsets != NULL )
   {
      int j;
      int k;

      assert(quadrowidxs != NULL);
      assert(quadidxs != NULL);
      assert(quadvals != NULL);
      
      oracle->objquadlen = quadoffsets[nquadrows];
      SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->objquadrows, oracle->objquadlen) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->objquadcols, oracle->objquadlen) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->objquadvals, quadvals, oracle->objquadlen) );
      
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
   }
   else
      oracle->objquadlen = 0;
   
   if( exprtree != NULL )
   {
      SCIPerrorMessage("nonquadratic functions not supported in NLPI yet\n");
      return SCIP_ERROR;
   }
   
   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ubs                 /**< new upper bounds, or NULL if all should be +infty */
   )
{
   int i;

   assert(scip != NULL);
   assert(oracle != NULL);
   assert(indices != NULL || nvars == 0);
   
   for( i = 0; i < nvars; ++i )
   {
      assert(indices != NULL);
      assert(indices[i] >= 0);
      assert(indices[i] < oracle->nvars);
      
      oracle->varlbs[indices[i]] = (lbs != NULL ? lbs[i] : -SCIPinfinity(scip));
      oracle->varubs[indices[i]] = (ubs != NULL ? ubs[i] : SCIPinfinity(scip));
      
      if( oracle->varlbs[indices[i]] > oracle->varubs[indices[i]] )
      { /* inconsistent bounds; let's assume it's due to rounding and make them equal */
         assert(SCIPisLE(scip, oracle->varlbs[indices[i]], oracle->varubs[indices[i]]) == TRUE);
         oracle->varlbs[indices[i]] = oracle->varubs[indices[i]];
      }
   }
   
   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_RETCODE SCIPnlpiOracleChgConsBounds(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             nconss,             /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lhss,               /**< new left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhss                /**< new right-hand sides, or NULL if all should be +infty */
   )
{
   int i;
   
   assert(scip != NULL);
   assert(oracle != NULL);
   assert(indices != NULL || nconss == 0);
   
   for( i = 0; i < nconss; ++i )
   {
      assert(indices != NULL);
      assert(indices[i] >= 0);
      assert(indices[i] < oracle->nconss);
      
      oracle->conslhss[indices[i]] = (lhss != NULL ? lhss[i] : -SCIPinfinity(scip));
      oracle->consrhss[indices[i]] = (rhss != NULL ? rhss[i] : SCIPinfinity(scip));
   }
   
   return SCIP_OKAY;
}

/** deletes a set of variables */
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< deletion status of vars in input (1 if var should be deleted, 0 if not); 
                                              *   new position of var in output (-1 if var was deleted) */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(oracle != NULL);
   
   /* TODO */
   SCIPerrorMessage("SCIPnlpiOracleDelVarSet is not implemented\n");
   
   return SCIP_ERROR;
}

/** deletes a set of constraints */
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not); 
                                              *   new position of row in output (-1 if row was deleted) */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(oracle != NULL);

   SCIPnlpiOracleInvalidateJacobiSparsity(scip, oracle);
   SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);

   /* TODO */
   SCIPerrorMessage("SCIPnlpiOracleDelConsSet is not implemented\n");
   
   return SCIP_ERROR;
}

/** changes linear coefficients in one constraint or objective
 *  @return error if coefficient did not exist before
 */
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidxs,            /**< array with indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoeffs           /**< array with new coefficients of variables */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(oracle != NULL);
   
   /* @TODO */
   SCIPerrorMessage("SCIPnlpiOracleChgLinearCoefs is not implemented\n");
   
   return SCIP_ERROR;
}
  
/** changes coefficients in the quadratic part of one constraint or objective
 *  @return error if coefficient did not exist before
 */
SCIP_RETCODE SCIPnlpiOracleChgQuadCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where quadratic coefficients should be changed, or -1 for objective */
   const int             nentries,           /**< number of coefficients to change */
   const int*            rowoffset,          /**< row offset containing modified indices */
   const int*            colidx,             /**< columns containint modified indices to the corresponding row offset */
   SCIP_Real*            newcoeff            /**< new quadratic coefficients */ 
   )
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(oracle != NULL);
   
   /* @TODO */
   SCIPerrorMessage("SCIPnlpiOracleChgQuadCoefs is not implemented\n");
   
   return SCIP_ERROR;
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

   if( oracle->conslinoffsets && oracle->conslinoffsets[considx] )
      return 1;

   return 0;
}

/** computes the value of a function */
static
SCIP_RETCODE SCIPnlpiOracleEvalFunctionValue(
   SCIP*                 scip,               /**< pointer to SCIP */
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
   assert(scip != NULL);
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
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** computes the value and gradient of a function */
static
SCIP_RETCODE SCIPnlpiOracleEvalFunctionGradient(
   SCIP*                 scip,               /**< pointer to SCIP */
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
   assert(scip != NULL);
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
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** evaluates the objective function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */  
   )
{
   assert(scip != NULL);
   assert(oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionValue(scip, oracle,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals,
      oracle->objexprvaridxs, oracle->objexprtree,
      x, objval) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** evaluates one constraint function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */  
   )
{
   assert(scip != NULL);
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionValue(scip, oracle,
         (oracle->conslinoffsets != NULL ? (oracle->conslinoffsets[considx+1] - oracle->conslinoffsets[considx]) : 0),
         (oracle->conslinoffsets != NULL ? &oracle->conslininds[oracle->conslinoffsets[considx]] : NULL),
         (oracle->conslinoffsets != NULL ? &oracle->conslinvals[oracle->conslinoffsets[considx]] : NULL),
         (oracle->consquadlens != NULL ? oracle->consquadlens[considx] : 0),
         (oracle->consquadlens != NULL ? oracle->consquadrows[considx] : NULL),
         (oracle->consquadlens != NULL ? oracle->consquadcols[considx] : NULL),
         (oracle->consquadlens != NULL ? oracle->consquadvals[considx] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprvaridxs[considx] : NULL),
         (oracle->consexprtrees != NULL ? oracle->consexprtrees[considx] : NULL),
         x, conval) );
   
   return SCIP_OKAY;
}

/** evaluates all constraint functions in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< buffer to store constraint values */  
   )
{
   int i;
   
   assert(scip != NULL);
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(convals != NULL);

   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_CALL( SCIPnlpiOracleEvalConstraintValue(scip, oracle, i, x, &convals[i]) );
   }
   
   return SCIP_OKAY;
}

/** computes the objective gradient in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Real*            objgrad             /**< pointer to store (dense) objective gradient */  
   )
{
   assert(scip != NULL);
   assert(oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionGradient(scip, oracle,
      oracle->objnlin, oracle->objlinidxs, oracle->objlinvals,
      oracle->objquadlen, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals,
      oracle->objexprvaridxs, oracle->objexprtree,
      x, newx, objval, objgrad) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** computes a constraints gradient in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether the function has not been evaluated for this point before */ 
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */  
   )
{
   assert(scip != NULL);
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionGradient(scip, oracle,
         (oracle->conslinoffsets != NULL ? (oracle->conslinoffsets[considx+1] - oracle->conslinoffsets[considx]) : 0),
         (oracle->conslinoffsets != NULL ? &oracle->conslininds[oracle->conslinoffsets[considx]] : NULL),
         (oracle->conslinoffsets != NULL ? &oracle->conslinvals[oracle->conslinoffsets[considx]] : NULL),
         (oracle->consquadlens != NULL ? oracle->consquadlens[considx] : 0),
         (oracle->consquadlens != NULL ? oracle->consquadrows[considx] : NULL),
         (oracle->consquadlens != NULL ? oracle->consquadcols[considx] : NULL),
         (oracle->consquadlens != NULL ? oracle->consquadvals[considx] : NULL),
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
   SCIP*                 scip,               /**< pointer to SCIP */
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
   
   assert(scip != NULL);
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->jacoffsets, oracle->nconss + 1) );

   maxnnz = MIN(oracle->nvars, 10) * oracle->nconss;  /* initial guess */
   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->jaccols, maxnnz) );
   
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
   
   SCIP_CALL( SCIPallocBufferArray(scip, &nzflag, oracle->nvars) );

   for( i = 0; i < oracle->nconss; ++i )
   {
      oracle->jacoffsets[i] = nnz;
      
      if( oracle->conslinoffsets != NULL && (oracle->consquadlens == NULL || oracle->consquadlens[i] == 0) && (oracle->consexprtrees == NULL || oracle->consexprtrees[i] == NULL) )
      { /* linear constraint: since we just want to copy the conslinvals at EvalJacobian, we need to copy conslininds here too (it can be unsorted, and sorting would confuse ChgLinearCoef later) */
         int nz = oracle->conslinoffsets[i+1] - oracle->conslinoffsets[i];
         if( nnz + nz > maxnnz )
         {
           maxnnz *= 2;
           SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->jaccols, maxnnz) );
         }
         BMScopyMemoryArray(&oracle->jaccols[nnz], &oracle->conslininds[oracle->conslinoffsets[i]], nz);
         nnz += nz;
         continue;
      }
      
      /* check which variables appear in constraint i */
      BMSclearMemoryArray(nzflag, oracle->nvars);
      
      if( oracle->conslinoffsets != NULL )
      {
         assert(oracle->conslininds);
         for( j = oracle->conslinoffsets[i]; j < oracle->conslinoffsets[i+1]; ++j )
            nzflag[oracle->conslininds[j]] = TRUE;
      }
      
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
      {
         assert(oracle->consquadrows != NULL );
         assert(oracle->consquadrows[i] != NULL );
         assert(oracle->consquadcols != NULL );
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
         
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
      }
      
      /* store variables indices in jaccols */
      for( j = 0; j < oracle->nvars; ++j )
      {
         if( nzflag[j] == FALSE )
            continue;
         
         if( nnz >= maxnnz )
         {
            maxnnz *= 2;
            SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->jaccols, maxnnz) );
         }
         
         oracle->jaccols[nnz] = j;
         ++nnz;
      }
   }
   
   oracle->jacoffsets[oracle->nconss] = nnz;
   
   if( nnz < maxnnz )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->jaccols, nnz) );
   }
   
   SCIPfreeBufferArray(scip, &nzflag);

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
   SCIP*                 scip,               /**< pointer to SCIP */
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
   
   assert(scip != NULL);
   assert(oracle != NULL);
   assert(jacobi != NULL);
   
   assert(oracle->jacoffsets != NULL);
   assert(oracle->jaccols != NULL);
   
   SCIP_CALL( SCIPallocBufferArray(scip, &grad, oracle->nvars) );
   
   j = oracle->jacoffsets[0];
   k = 0;
   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->conslinoffsets != NULL && (oracle->consquadlens == NULL || oracle->consquadlens[i] == 0) && (oracle->consexprtrees == NULL || oracle->consexprtrees[i] == NULL) )
      { /* linear constraints should be easy */
         l = oracle->jacoffsets[i+1] - j;
         assert(l == oracle->conslinoffsets[i+1] - oracle->conslinoffsets[i]);
         BMScopyMemoryArray(&jacobi[k], &oracle->conslinvals[oracle->conslinoffsets[i]], l);
         j += l; 
         k += l;
      }
      else
      { /* @TODO do this sparse too */
         SCIP_RETCODE retcode = SCIPnlpiOracleEvalConstraintGradient(scip, oracle, i, x, newx, (convals ? &convals[i] : &dummy), grad);
         if( retcode != SCIP_OKAY )
         {
            SCIPfreeBufferArray(scip, &grad);
            return retcode;
         }
      
         for( ; j < oracle->jacoffsets[i+1]; ++j, ++k )
            jacobi[k] = grad[oracle->jaccols[j]];
      }
   }
   
   SCIPfreeBufferArray(scip, &grad);

   return SCIP_OKAY;
}

/** sets the nzflags and increases the nzcount given indices of quadratic terms */
static
void SCIPnlpiOracleHessLagSparsitySetNzFlagForQuad(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_Bool*            nzflag,             /**< where to set flags */
   int*                  nzcount,            /**< counter for total number of nonzeros; should be increased when nzflag is set to 1 the first time */
   int*                  rowidx,             /**< row indices */
   int*                  colidx,             /**< column indices */
   int                   length              /**< length of quadratic part */
   )
{
   int i;

   assert(scip != NULL);
   assert(nzflag != NULL);
   assert(nzcount != NULL);
   assert(rowidx != NULL);
   assert(colidx != NULL);
   assert(length >= 0);

   for( ; length > 0; --length, ++rowidx, ++colidx )
   {
      assert(*colidx <= *rowidx);
      
      i = (*rowidx * (*rowidx+1)) / 2 + *colidx;
      if( nzflag[i] == FALSE )
      {
         nzflag[i] = TRUE;
         ++(*nzcount);
      }
   }
}

/** Gets sparsity pattern of the Hessian matrix of the Lagrangian.
 * 
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 *  Only elements of the lower left triangle and the diagonal are counted.
 */
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[nconss] gives length of col, can be NULL */
   )
{
   SCIP_Bool* nzflag;
   SCIP_Bool* nz;
   int nnz;
   int i;
   int j;
   int cnt;
   
   assert(scip  != NULL);
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->heslagoffsets, oracle->nvars + 1) );
   nnz = 0;
   
   /* @TODO should be improved */
   SCIP_CALL( SCIPallocBufferArray(scip, &nzflag, (oracle->nvars * (oracle->nvars + 1)) / 2) );
   BMSclearMemoryArray(nzflag, (oracle->nvars * (oracle->nvars + 1)) / 2);
   
   if( oracle->objquadlen != 0 )
      SCIPnlpiOracleHessLagSparsitySetNzFlagForQuad(scip, nzflag, &nnz, oracle->objquadrows, oracle->objquadcols, oracle->objquadlen);

   if( oracle->objexprtree != NULL )
   {
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
   }

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->consquadlens != NULL && oracle->consquadlens[i] != 0 )
         SCIPnlpiOracleHessLagSparsitySetNzFlagForQuad(scip, nzflag, &nnz, oracle->consquadrows[i], oracle->consquadcols[i], oracle->consquadlens[i]);
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
      }
   }
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->heslagcols, nnz) );
   
   /* set hessian sparsity from nzflag */
   cnt = 0;
   nz = nzflag;
   for( i = 0; i < oracle->nvars; ++i )
   {
      oracle->heslagoffsets[i] = cnt;
      for( j = 0; j <= i; ++j, ++nz )
         if( *nz == TRUE )
         {
            assert(cnt < nnz);
            oracle->heslagcols[cnt++] = j;
         }
   }
   oracle->heslagoffsets[oracle->nvars] = cnt;
   assert(cnt == nnz);
   
   SCIPfreeBufferArray(scip, &nzflag);
   
   if( offset != NULL )
      *offset = oracle->heslagoffsets;
   if( col != NULL )
      *col = oracle->heslagcols;

   return SCIP_OKAY;
}

/** adds quadratic part of a constraint into hessian structure */
static
SCIP_RETCODE SCIPnlpiOracleHessLagAddQuad(
   SCIP*                 scip,               /**< pointer to SCIP */
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

   assert(scip != NULL);
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

/** Evaluates the Hessian matrix of the Lagrangian in a given point.
 * 
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity.
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity at least ones before using this function. 
 *  Only elements of the lower left triangle and the diagonal are computed.
 */
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             newx,               /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambda,             /**< weights (Lagrangian multipliers) for the constraints */ 
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */  
   )
{  /*lint --e{715}*/
   int i;

   assert(scip != NULL);
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
         SCIP_CALL( SCIPnlpiOracleHessLagAddQuad(scip, objfactor, oracle->objquadrows, oracle->objquadcols, oracle->objquadvals, 
               oracle->objquadlen, oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
      
      if( oracle->objexprtree != NULL )
      {
         assert(oracle->objexprvaridxs != NULL );
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
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
         SCIP_CALL( SCIPnlpiOracleHessLagAddQuad(scip, lambda[i], oracle->consquadrows[i], oracle->consquadcols[i], oracle->consquadvals[i], 
               oracle->consquadlens[i], oracle->heslagoffsets, oracle->heslagcols, hessian) );
      }
      
      if( oracle->consexprtrees != NULL && oracle->consexprtrees[i] != NULL )
      {
         assert(oracle->consexprvaridxs != NULL);
         assert(oracle->consexprvaridxs[i] != NULL);
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
      }
   }
   
   return SCIP_OKAY;
}

/** prints a function */ 
static
SCIP_RETCODE SCIPnlpiOraclePrintFunction(
   SCIP*                 scip,               /**< pointer to SCIP */
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

   assert(scip != NULL);
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
      
   return SCIP_OKAY;
}

/** prints the problem to a file. */
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;

   assert(scip != NULL);
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
   SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
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
      if( !SCIPisInfinity(scip, -oracle->conslhss[i]) && !SCIPisInfinity(scip, oracle->consrhss[i]) && oracle->conslhss[i] != oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, "%g <= ", oracle->conslhss[i]);
      
      SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
            (oracle->conslinoffsets != NULL ? (oracle->conslinoffsets[i+1] - oracle->conslinoffsets[i]) : 0),
            (oracle->conslinoffsets != NULL ? &oracle->conslininds[oracle->conslinoffsets[i]] : NULL),
            (oracle->conslinoffsets != NULL ? &oracle->conslinvals[oracle->conslinoffsets[i]] : NULL),
            (oracle->consquadlens != NULL ? oracle->consquadlens[i] : 0),
            (oracle->consquadlens != NULL ? oracle->consquadrows[i] : NULL),
            (oracle->consquadlens != NULL ? oracle->consquadcols[i] : NULL),
            (oracle->consquadlens != NULL ? oracle->consquadvals[i] : NULL),
            (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL)) );
      
      if( oracle->conslhss[i] == oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, " = %g", oracle->consrhss[i]);
      else if( !SCIPisInfinity(scip, oracle->consrhss[i]) )
         SCIPmessageFPrintInfo(file, " <= %g", oracle->consrhss[i]);
      else if( !SCIPisInfinity(scip, -oracle->conslhss[i]) )
         SCIPmessageFPrintInfo(file, " >= %g", oracle->conslhss[i]);
      
      SCIPmessageFPrintInfo(file, "\n");
   }

   return SCIP_OKAY;
}

/** prints the problem to a file in GAMS format */
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;

   assert(scip != NULL);
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
         if( !SCIPisInfinity(scip, -oracle->varlbs[i]) )
         {
            if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
               SCIPmessageFPrintInfo(file, "%s", oracle->varnames[i]);
            else
               SCIPmessageFPrintInfo(file, "x%d", i);
            SCIPmessageFPrintInfo(file, ".lo = %g;\t", oracle->varlbs[i]);
         }
         if( !SCIPisInfinity(scip, oracle->varubs[i]) )
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

      if( !SCIPisInfinity(scip, -oracle->conslhss[i]) && !SCIPisInfinity(scip, oracle->consrhss[i]) && oracle->conslhss[i] != oracle->consrhss[i] )
      { /* ranged row: add second constraint */
         if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
            SCIPmessageFPrintInfo(file, "%s_RNG, ", oracle->consnames[i]);
         else
            SCIPmessageFPrintInfo(file, "e%d_RNG, ", i);
      }
   }
   SCIPmessageFPrintInfo(file, "OBJ;\n\n");
   
   SCIPmessageFPrintInfo(file, "OBJ.. OBJVAR =E= ");
   SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
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
      
      SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
            (oracle->conslinoffsets != NULL ? (oracle->conslinoffsets[i+1] - oracle->conslinoffsets[i]) : 0),
            (oracle->conslinoffsets != NULL ? &oracle->conslininds[oracle->conslinoffsets[i]] : NULL),
            (oracle->conslinoffsets != NULL ? &oracle->conslinvals[oracle->conslinoffsets[i]] : NULL),
            (oracle->consquadlens != NULL ? oracle->consquadlens[i] : 0),
            (oracle->consquadlens != NULL ? oracle->consquadrows[i] : NULL),
            (oracle->consquadlens != NULL ? oracle->consquadcols[i] : NULL),
            (oracle->consquadlens != NULL ? oracle->consquadvals[i] : NULL),
            (oracle->consexprtrees != NULL ? oracle->consexprtrees[i] : NULL)) );

      if( oracle->conslhss[i] == oracle->consrhss[i] )
         SCIPmessageFPrintInfo(file, " =E= %g", oracle->consrhss[i]);
      else if( !SCIPisInfinity(scip, oracle->consrhss[i]) )
         SCIPmessageFPrintInfo(file, " =L= %g", oracle->consrhss[i]);
      else if( !SCIPisInfinity(scip, -oracle->conslhss[i]) )
         SCIPmessageFPrintInfo(file, " =G= %g", oracle->conslhss[i]);
      else
         SCIPmessageFPrintInfo(file, " =N= 0");
      SCIPmessageFPrintInfo(file, ";\n");

      if( !SCIPisInfinity(scip, -oracle->conslhss[i]) && !SCIPisInfinity(scip, oracle->consrhss[i]) && oracle->conslhss[i] != oracle->consrhss[i] )
      {
         if( oracle->consnames != NULL && oracle->consnames[i] != NULL )
            SCIPmessageFPrintInfo(file, "%s", oracle->consnames[i]);
         else
            SCIPmessageFPrintInfo(file, "e%d", i);
         SCIPmessageFPrintInfo(file, "_RNG.. ");
        
         SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
               (oracle->conslinoffsets != NULL ? (oracle->conslinoffsets[i+1] - oracle->conslinoffsets[i]) : 0),
               (oracle->conslinoffsets != NULL ? &oracle->conslininds[oracle->conslinoffsets[i]] : NULL),
               (oracle->conslinoffsets != NULL ? &oracle->conslinvals[oracle->conslinoffsets[i]] : NULL),
               (oracle->consquadlens != NULL ? oracle->consquadlens[i] : 0),
               (oracle->consquadlens != NULL ? oracle->consquadrows[i] : NULL),
               (oracle->consquadlens != NULL ? oracle->consquadcols[i] : NULL),
               (oracle->consquadlens != NULL ? oracle->consquadvals[i] : NULL),
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
