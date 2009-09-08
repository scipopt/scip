/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpi_oracle.c,v 1.14 2009/09/08 19:35:59 bzfviger Exp $"

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
   SCIP_Real*            varlb;              /**< variable lower bounds */
   SCIP_Real*            varub;              /**< variable upper bounds */
   char**                varname;            /**< variable names */
   
   int                   ncons;              /**< number of constraints */
   SCIP_Real*            conlhs;             /**<  left-hand side of constraints */
   SCIP_Real*            conrhs;             /**< right-hand side of constraints */
   int*                  conlinoffset;       /**< start index of each constraints linear coefficients in linind and linval (length: ncons + 1, linoffset[ncons] gives length of linind and linval) */
   int*                  conlinind;          /**< variable indices in linear part */
   SCIP_Real*            conlinval;          /**< variable coefficients in linear part */ 
   int*                  conquadlen;         /**< quadlen[.] gives length of quadratic part, or 0 if no quadratic part in this constraint */
   int**                 conquadrow;         /**< quadrow[.] gives row    indices for quadratic part, or NULL if no quadratic part in this constraint */
   int**                 conquadcol;         /**< quadcol[.] gives column indices for quadratic part, or NULL if no quadratic part in this constraint */
   SCIP_Real**           conquadval;         /**< quadval[.] gives coefficients in quadratic part, or NULL if no quadratic part in this constraint */
   int**                 conexprvaridx;      /**< exprvaridx[.] gives indices of variables from expression in NLP, or NULL if no nonquadratic part in this constraint */
   SCIP_EXPRTREE**       conexprtree;        /**< exprtree[.] gives nonquadratic part, or NULL if no nonquadratic part in this constraint */
   char**                conname;            /**< constraint names */

   SCIP_Real             objconstant;        /**< constant part of objective */
   int                   objnlin;            /**< number of linear variable coefficients in objective */ 
   int*                  objlinind;          /**< indices of linear variables in objective, or NULL if no linear part */
   SCIP_Real*            objlinval;          /**< coefficients of linear variables in objective, or NULL if no linear part */
   int                   objquadlen;         /**< length of quadratic part of objective, or 0 if no quadratic part in objective */
   int*                  objquadrow;         /**< row    indices in quadratic part of objective, or NULL if no quadratic part in objective */
   int*                  objquadcol;         /**< column indices in quadratic part in objective, or NULL if no quadratic part in objective */ 
   SCIP_Real*            objquadval;         /**< coefficients in quadratic part in objective, or NULL if no quadratic part in objective */
   int*                  objexprvaridx;      /**< exprvaridx[.] gives indices of variables from expression in NLP, or NULL if no nonquadratic part in objective */
   SCIP_EXPRTREE*        objexprtree;        /**< expression tree of nonquadratic part in objective, or NULL if no nonquadratic part in objective */

   int*                  jacoffset;          /**< rowwise jacobi sparsity pattern: constraint offsets in jaccol */
   int*                  jaccol;             /**< rowwise jacobi sparsity pattern: indices of variables appearing in constraints */
   
   int*                  heslagoffset;       /**< rowwise sparsity pattern of hessian matrix of Lagrangian: row offsets in laghescol */
   int*                  heslagcol;          /**< rowwise sparsity pattern of hessian matrix of Lagrangian: column indices; sorted for each row */
   
   int*                  vardegree;          /**< maximal degree of variable over objective and all constraints */
};

/** Invalidates the sparsity pattern of the Jacobian.
 * Should be called when constraints are added or deleted.
 */
static
void SCIPnlpiOracleInvalidateJacobiSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( oracle->jacoffset == NULL )
   { /* nothing to do */
      assert(oracle->jaccol == NULL);
      return;
   }
   
   assert(oracle->jaccol != NULL);
   SCIPfreeMemoryArray(scip, &oracle->jacoffset);
   SCIPfreeMemoryArray(scip, &oracle->jaccol);
}

/** Invalidates the sparsity pattern of the Hessian of the Lagragian.
 * Should be called when the objective is set or constraints are added or deleted.
 */
static
void SCIPnlpiOracleInvalidateHessianLagSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( oracle->heslagoffset == NULL )
   { /* nothing to do */
      assert(oracle->heslagcol == NULL);
      return;
   }
   
   assert(oracle->heslagcol != NULL);
   SCIPfreeMemoryArray(scip, &oracle->heslagoffset);
   SCIPfreeMemoryArray(scip, &oracle->heslagcol);
}

/** creates an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleCreate(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   SCIP_CALL( SCIPallocMemory(scip, oracle) );
   
   (*oracle)->nvars         = 0;
   (*oracle)->varlb         = NULL;
   (*oracle)->varub         = NULL;
   (*oracle)->varname       = NULL;

   (*oracle)->ncons         = 0;
   (*oracle)->conlhs        = NULL;
   (*oracle)->conrhs        = NULL;
   (*oracle)->conlinoffset  = NULL;
   (*oracle)->conlinind     = NULL;
   (*oracle)->conlinval     = NULL;
   (*oracle)->conquadlen    = NULL;
   (*oracle)->conquadrow    = NULL;
   (*oracle)->conquadcol    = NULL;
   (*oracle)->conquadval    = NULL;
   (*oracle)->conexprvaridx = NULL;
   (*oracle)->conexprtree   = NULL;
   (*oracle)->conname       = NULL;

   (*oracle)->objconstant   = 0.0;
   (*oracle)->objnlin       = 0;
   (*oracle)->objlinind     = NULL;
   (*oracle)->objlinval     = NULL;
   (*oracle)->objquadlen    = 0;
   (*oracle)->objquadrow    = NULL;
   (*oracle)->objquadcol    = NULL;
   (*oracle)->objquadval    = NULL;
   (*oracle)->objexprvaridx = NULL;
   (*oracle)->objexprtree   = NULL;
   
   (*oracle)->jacoffset     = NULL;
   (*oracle)->jaccol        = NULL;
   
   (*oracle)->heslagoffset  = NULL;
   (*oracle)->heslagcol     = NULL;
   
   (*oracle)->vardegree     = NULL;
   
   return SCIP_OKAY;
}

/** initializes NLPI oracle */
SCIP_RETCODE SCIPnlpiOracleInit(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle              /**< NLPIORACLE data structure */
)
{  /*lint --e{715}*/
   assert(scip != NULL);

   return SCIP_OKAY;
}

/** frees an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
)
{
   int i;
   assert(scip    != NULL);
   assert(oracle  != NULL);
   assert(*oracle != NULL);
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->varlb);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->varub);
   if( (*oracle)->varname != NULL )
   {
      for( i = 0; i < (*oracle)->nvars; ++i )
      {
        SCIPfreeMemoryArrayNull(scip, &(*oracle)->varname[i]);
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->varname);
   }
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conlhs);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conrhs);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conlinoffset);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conlinind);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->conlinval);
   
   if( (*oracle)->conquadlen != NULL )
   {
      assert((*oracle)->conquadrow != NULL);
      assert((*oracle)->conquadcol != NULL);
      assert((*oracle)->conquadval != NULL);
      for( i = 0; i < (*oracle)->ncons; ++i )
      {
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->conquadrow[i]);
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->conquadcol[i]);
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->conquadval[i]);
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->conquadlen);
      SCIPfreeMemoryArray(scip, &(*oracle)->conquadrow);
      SCIPfreeMemoryArray(scip, &(*oracle)->conquadcol);
      SCIPfreeMemoryArray(scip, &(*oracle)->conquadval);
   }
   
   if( (*oracle)->conexprtree != NULL )
   {
      for( i = 0; i < (*oracle)->ncons; ++i )
      {
         if( (*oracle)->conexprtree[i] != NULL )
         {
            SCIPfreeMemoryArrayNull(scip, &(*oracle)->conexprvaridx[i]);
         }
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->conexprtree);
      SCIPfreeMemoryArray(scip, &(*oracle)->conexprvaridx);
   }
   
   if( (*oracle)->conname != NULL )
   {
      for( i = 0; i < (*oracle)->ncons; ++i )
      {
         SCIPfreeMemoryArrayNull(scip, &(*oracle)->conname[i]);
      }
      SCIPfreeMemoryArray(scip, &(*oracle)->conexprtree);
   }
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objlinind);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objlinval);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objquadrow);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objquadcol);
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->objquadval);
   if( (*oracle)->objexprtree != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*oracle)->objexprvaridx);
   }

   SCIPnlpiOracleInvalidateJacobiSparsity(scip, *oracle);
   SCIPnlpiOracleInvalidateHessianLagSparsity(scip, *oracle);
   
   SCIPfreeMemoryArrayNull(scip, &(*oracle)->vardegree);

   SCIPfreeMemory(scip, oracle);
   
   return SCIP_OKAY;
}

/** adds variables */
SCIP_RETCODE SCIPnlpiOracleAddVars(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to add */
   const SCIP_Real*      lb,                 /**< lower bounds of new variables, or NULL if all -infty */
   const SCIP_Real*      ub,                 /**< upper bounds of new variables, or NULL if all +infty */
   const char**          varnames            /**< names of new variables, or NULL if no names should be stored */
)
{
   int i;
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( nvars == 0 )
      return SCIP_OKAY;
   
   assert(nvars > 0);
   
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->varlb, oracle->nvars + nvars) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->varub, oracle->nvars + nvars) );
   
   if( lb != NULL )
      BMScopyMemoryArray(oracle->varlb + oracle->nvars, lb, nvars);
   else
      for( i = 0; i < nvars; ++i )
         oracle->varlb[oracle->nvars+i] = -SCIPinfinity(scip);
   
   if( ub != NULL )
      BMScopyMemoryArray(oracle->varub + oracle->nvars, ub, nvars);
   else
      for( i = 0; i < nvars; ++i )
         oracle->varub[oracle->nvars+i] =  SCIPinfinity(scip);
   
   if( oracle->varname != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->varname, oracle->nvars + nvars) );
      if( varnames != NULL )
      {
         for( i = 0; i < nvars; ++i )
         {
            if( varnames[i] != NULL )
               SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->varname[oracle->nvars+i], varnames[i], strlen(varnames[i])+1) );
            else
               oracle->varname[oracle->nvars+i] = NULL;
         }
      }
      else
      {
         BMSclearMemoryArray(oracle->varname + oracle->nvars, nvars);
      }
   }
   else if( varnames != NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->varname, oracle->nvars + nvars) );
      BMSclearMemoryArray(oracle->varname, oracle->nvars);
      for( i = 0; i < nvars; ++i )
      {
         if( varnames[i] != NULL )
         {
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->varname[oracle->nvars+i], varnames[i], strlen(varnames[i])+1) );
         }
         else
            oracle->varname[oracle->nvars+i] = NULL;
      }
   }
   
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->vardegree, oracle->nvars + nvars) );
   BMSclearMemoryArray(oracle->vardegree + oracle->nvars, nvars);

   /* @TODO update sparsity pattern by extending heslagoffset */
   SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);

   oracle->nvars += nvars;

   return SCIP_OKAY;
}

/** adds constraints 
 * 
 * linear coefficients: row(=constraint) oriented matrix;
 * quadratic coefficiens: row oriented matrix for each constraint
 */
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   ncons,              /**< number of constraints to add */
   const SCIP_Real*      lhs,                /**<  left-hand side of constraints, or NULL if all -infty */
   const SCIP_Real*      rhs,                /**< right-hand side of constraints, or NULL if all +infty */
   const int*            linoffset,          /**< start index of each constraints linear coefficients in linind and linval (length: ncons + 1, linoffset[ncons] gives length of linind and linval), or NULL if no linear part */
   const int*            linind,             /**< variable indices in linear part, or NULL if no linear part */
   const SCIP_Real*      linval,             /**< variable coefficients in linear part, of NULL if no linear part */ 
   const int*            nquadrows,          /**< NULL if no quadratic parts, otherwise nquadrows[.] gives the number of columns in the matrix of the quadratic part, or 0 if no quadratic part in this constraint */
   int* const*           quadrowidx,         /**< NULL if no quadratic parts, otherwise quadrowidx[.] gives the indices of variables for which a quadratic part is specified, or NULL if no quadratic part in this constraint */
   int* const*           quadoffset,         /**< NULL if no quadratic parts, otherwise quadoffset[.] gives start index of each rows quadratic coefficients in quadind[.] and quadval[.] (quadoffset[.][nvars] gives length of quadind[.] and quadval[.]), or NULL if no quadratic part in this constraint */
   int* const*           quadind,            /**< NULL if no quadratic parts, otherwise quadind[.] gives column indices for quadratic part, or NULL if no quadratic part in this constraint */
   SCIP_Real* const*     quadval,            /**< NULL if no quadratic parts, otherwise quadval[.] gives coefficients in quadratic part, or NULL if no quadratic part in this constraint */
   int* const*           exprvaridx,         /**< NULL if no nonquadratic parts, otherwise epxrvaridx[.] maps variable indices in expression tree to indices in nlp */
   SCIP_EXPRTREE* const* exprtree,           /**< NULL if no nonquadratic parts, otherwise exprtree[.] gives nonquadratic part, or NULL if no nonquadratic part in this constraint */
   const char**          connames            /**< names of new constraints, or NULL if no names should be stored */
   )
{  /*lint --e{715}*/
   int        i, oldnnz;
   SCIP_Bool  addednlcon = FALSE;  /* whether a nonlinear constraint was added */
   
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( ncons == 0 )
      return SCIP_OKAY;
   
   assert(ncons > 0);

   SCIPnlpiOracleInvalidateJacobiSparsity(scip, oracle); /* @TODO we could also update (extend) the sparsity pattern */

   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conlhs, oracle->ncons + ncons) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conrhs, oracle->ncons + ncons) );

   if( lhs != NULL )
   {
      BMScopyMemoryArray(oracle->conlhs + oracle->ncons, lhs, ncons);
   }
   else
      for( i = 0; i < ncons; ++i )
         oracle->conlhs[oracle->ncons+i] = -SCIPinfinity(scip);
   
   if( rhs != NULL )
   {
     BMScopyMemoryArray(oracle->conrhs + oracle->ncons, rhs, ncons);
   }
   else
      for( i = 0; i < ncons; ++i )
         oracle->conrhs[oracle->ncons+i] =  SCIPinfinity(scip);

   /* old number of nonzeros */
   oldnnz = oracle->ncons != 0 ? oracle->conlinoffset[oracle->ncons] : 0;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conlinoffset, oracle->ncons + ncons + 1) );
   if( linoffset != NULL )
   {
      assert(linind != NULL);
      assert(linval != NULL);
      
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conlinind, oldnnz + linoffset[ncons]) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conlinval, oldnnz + linoffset[ncons]) );
   
      for( i = 0; i <= ncons; ++i )
         oracle->conlinoffset[oracle->ncons + i] = oldnnz + linoffset[i];
      BMScopyMemoryArray(oracle->conlinind + oldnnz, linind, linoffset[ncons]);
      BMScopyMemoryArray(oracle->conlinval + oldnnz, linval, linoffset[ncons]);
      
      for( i = 0; i < linoffset[ncons]; ++i )
         oracle->vardegree[linind[i]] = MAX(1, oracle->vardegree[linind[i]]);
   }
   else if( oracle->conlinoffset != NULL )
   {
      for( i = 0; i <= ncons; ++i )
         oracle->conlinoffset[oracle->ncons + i] = oldnnz;
   }
   
   if( quadoffset != NULL )
   {
      assert(nquadrows  != NULL);
      assert(quadrowidx != NULL);
      assert(quadind    != NULL);
      assert(quadval    != NULL);
      if( (oracle->conquadlen == NULL) && (oracle->ncons != 0) )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conquadlen, oracle->ncons + ncons) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conquadrow, oracle->ncons + ncons) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conquadcol, oracle->ncons + ncons) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conquadval, oracle->ncons + ncons) );
         BMSclearMemoryArray(oracle->conquadlen, oracle->ncons);
         BMSclearMemoryArray(oracle->conquadrow, oracle->ncons);
         BMSclearMemoryArray(oracle->conquadcol, oracle->ncons);
         BMSclearMemoryArray(oracle->conquadval, oracle->ncons);
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadlen, oracle->ncons + ncons) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadrow, oracle->ncons + ncons) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadcol, oracle->ncons + ncons) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadval, oracle->ncons + ncons) );
      }
      for( i = 0; i < ncons; ++i )
      {
         if( quadoffset[i] != NULL && nquadrows[i] != 0 )
         {
            int j, k;
            assert(quadrowidx[i] != NULL);
            assert(quadind[i]    != NULL);
            assert(quadval[i]    != NULL);
            oracle->conquadlen[oracle->ncons + i] = quadoffset[i][nquadrows[i]];
            SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conquadrow[oracle->ncons + i], quadoffset[i][nquadrows[i]]) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conquadcol[oracle->ncons + i], quadoffset[i][nquadrows[i]]) );
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->conquadval[oracle->ncons + i], quadval[i], quadoffset[i][nquadrows[i]]) );
            
            k = quadoffset[i][0];
            for( j = 0; j < nquadrows[i]; ++j )
            {
               if( quadoffset[i][j] == quadoffset[i][j+1] )
                  continue;
               for( ; k < quadoffset[i][j+1]; ++k )
               {
                  assert(quadind[i][k] < nquadrows[i]);
                  if( quadrowidx[i][j] > quadrowidx[i][quadind[i][k]] )
                  {
                     oracle->conquadrow[oracle->ncons + i][k] = quadrowidx[i][j];
                     oracle->conquadcol[oracle->ncons + i][k] = quadrowidx[i][quadind[i][k]];
                  }
                  else
                  {
                     oracle->conquadrow[oracle->ncons + i][k] = quadrowidx[i][quadind[i][k]];
                     oracle->conquadcol[oracle->ncons + i][k] = quadrowidx[i][j];
                  }
                  oracle->vardegree[quadrowidx[i][quadind[i][k]]] = MAX(2, oracle->vardegree[quadrowidx[i][quadind[i][k]]]);
               }
               oracle->vardegree[quadrowidx[i][j]] = MAX(2, oracle->vardegree[quadrowidx[i][j]]);
            }
            assert(k == oracle->conquadlen[oracle->ncons + i]);
            
            addednlcon = TRUE;
         }
         else
         {
            oracle->conquadlen[oracle->ncons + i] = 0;
            oracle->conquadrow[oracle->ncons + i] = NULL;
            oracle->conquadcol[oracle->ncons + i] = NULL;
            oracle->conquadval[oracle->ncons + i] = NULL;
         }
      }
   }
   else if( oracle->conquadlen != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadlen, oracle->ncons + ncons) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadrow, oracle->ncons + ncons) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadcol, oracle->ncons + ncons) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conquadval, oracle->ncons + ncons) );
      BMSclearMemoryArray(oracle->conquadlen + oracle->ncons, ncons);
      BMSclearMemoryArray(oracle->conquadrow + oracle->ncons, ncons);
      BMSclearMemoryArray(oracle->conquadcol + oracle->ncons, ncons);
      BMSclearMemoryArray(oracle->conquadval + oracle->ncons, ncons);
   }
   
   if( exprtree != NULL )
   {
      if( oracle->conexprtree == NULL && oracle->ncons != 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conexprtree,   oracle->ncons + ncons) );
         BMSclearMemoryArray(oracle->conexprtree,   oracle->ncons);
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conexprvaridx, oracle->ncons + ncons) );
         BMSclearMemoryArray(oracle->conexprvaridx, oracle->ncons);
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conexprtree,   oracle->ncons + ncons) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conexprvaridx, oracle->ncons + ncons) );
      }
      for( i = 0; i < ncons; ++i )
         if( exprtree[i] != NULL )
         {
            SCIPerrorMessage("nonquadratic functions not supported in NLPI yet.\n");
            return SCIP_ERROR;
         }
         else
         {
            oracle->conexprtree[oracle->ncons + i]   = NULL;
            oracle->conexprvaridx[oracle->ncons + i] = NULL;
         }
   }
   else if( oracle->conexprtree != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conexprtree,   oracle->ncons + ncons) );
      BMSclearMemoryArray(oracle->conexprtree   + oracle->ncons, ncons);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conexprvaridx, oracle->ncons + ncons) );
      BMSclearMemoryArray(oracle->conexprvaridx + oracle->ncons, ncons);
   }

   if( connames != NULL )
   {
      if( oracle->conname == NULL && oracle->ncons != 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->conname, oracle->ncons + ncons) );
         BMSclearMemoryArray(oracle->conname, oracle->ncons);
      }
      else
         SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conname, oracle->ncons + ncons) );
      for( i = 0; i < ncons; ++i )
         if( connames[i] != NULL )
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->conname[oracle->ncons + i], connames[i], sizeof(connames[i])+1) );
         else
            oracle->conname[oracle->ncons + i] = NULL;
   }
   else if( oracle->conname != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->conname, oracle->ncons + ncons) );
      BMSclearMemoryArray(oracle->conname + oracle->ncons, ncons);
   }
   
   oracle->ncons += ncons;

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
   const int*            linind,             /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linval,             /**< coefficients of linear variables, or NULL if no linear part */
   int                   nquadrows,          /**< number of columns in the matrix of the quadratic part */
   const int*            quadrowidx,         /**< indices of variables appearing in quadratic part, or NULL if no quadratic part */
   const int*            quadoffset,         /**< start index of each rows quadratic coefficients in quadind and quadval (quadoffset[.][nvars] gives length of quadind and quadval), or NULL if no quadratic part */
   const int*            quadind,            /**< column indices in quadratic part, or NULL if no quadratic part */ 
   const SCIP_Real*      quadval,            /**< coefficients in quadratic part, or NULL if no quadratic part */
   const int*            exprvaridx,         /**< maps variable indices in expression tree to indices in nlp, or NULL if no nonquadratic part */
   const SCIP_EXPRTREE*  exprtree            /**< expression tree of nonquadratic part, or NULL if no nonquadratic part */
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   /* clear previous objective */
   SCIPfreeMemoryArrayNull(scip, &oracle->objlinind);
   SCIPfreeMemoryArrayNull(scip, &oracle->objlinval);
   SCIPfreeMemoryArrayNull(scip, &oracle->objquadrow);
   SCIPfreeMemoryArrayNull(scip, &oracle->objquadcol);
   SCIPfreeMemoryArrayNull(scip, &oracle->objquadval);
   if( oracle->objexprtree != NULL )
   {
      SCIPfreeMemoryArray(scip, &oracle->objexprvaridx);
      /* @TODO this does not clear the vardegree's */
   }
   
   if( nquadrows != 0 || oracle->objquadlen != 0 || exprtree != NULL || oracle->objexprtree != NULL )
      SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);
   
   oracle->objconstant = constant;
   
   if( nlin != 0 )
   {
      int i;
      assert(linind != NULL);
      assert(linval != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->objlinind, linind, nlin) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->objlinval, linval, nlin) );
      oracle->objnlin = nlin;
      for( i = 0; i < nlin; ++i )
         oracle->vardegree[linind[i]] = MAX(1, oracle->vardegree[linind[i]]);
   }
   
   if( nquadrows != 0 && quadoffset != NULL )
   {
      int j, k;
      assert(quadrowidx != NULL);
      assert(quadind    != NULL);
      assert(quadval    != NULL);
      
      oracle->objquadlen = quadoffset[nquadrows];
      SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->objquadrow, oracle->objquadlen) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->objquadcol, oracle->objquadlen) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &oracle->objquadval, quadval, oracle->objquadlen) );
      
      k = quadoffset[0];
      for( j = 0; j < nquadrows; ++j )
      {
         if( quadoffset[j] == quadoffset[j+1] )
            continue;
         for( ; k < quadoffset[j+1]; ++k )
         {
            assert(quadind[k] < nquadrows);
            if( quadrowidx[j] > quadrowidx[quadind[k]] )
            {
               oracle->objquadrow[k] = quadrowidx[j];
               oracle->objquadcol[k] = quadrowidx[quadind[k]];
            }
            else
            {
               oracle->objquadrow[k] = quadrowidx[quadind[k]];
               oracle->objquadcol[k] = quadrowidx[j];
            }
            oracle->vardegree[quadrowidx[quadind[k]]] = MAX(2, oracle->vardegree[quadrowidx[quadind[k]]]);
         }
         oracle->vardegree[quadrowidx[j]] = MAX(2, oracle->vardegree[quadrowidx[j]]);
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
   const SCIP_Real*      lb,                 /**< new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ub                  /**< new upper bounds, or NULL if all should be +infty */
)
{
   int i;

   assert(scip    != NULL);
   assert(oracle  != NULL);
   assert(indices != NULL || nvars == 0);
   
   for( i = 0; i < nvars; ++i )
   {
      assert(indices    != NULL); /* for lint */
      assert(indices[i] >= 0);
      assert(indices[i] <  oracle->nvars);
      
      oracle->varlb[indices[i]] = lb != NULL ? lb[i] : -SCIPinfinity(scip);
      oracle->varub[indices[i]] = ub != NULL ? ub[i] :  SCIPinfinity(scip);
      
      if( oracle->varlb[indices[i]] > oracle->varub[indices[i]] )
      { /* inconsistent bounds; let's assume it's due to rounding and make them equal */
         assert(SCIPisLE(scip, oracle->varlb[indices[i]], oracle->varub[indices[i]]) == TRUE);
         oracle->varlb[indices[i]] = oracle->varub[indices[i]];
      }
   }
   
   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_RETCODE SCIPnlpiOracleChgConsBounds(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             ncons,              /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lhs,                /**< new  left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhs                 /**< new right-hand sides, or NULL if all should be +infty */
)
{
   int i;
   
   assert(scip    != NULL);
   assert(oracle  != NULL);
   assert(indices != NULL || ncons == 0);
   
   for( i = 0; i < ncons; ++i )
   {
      assert(indices    != NULL); /* for lint */
      assert(indices[i] >= 0);
      assert(indices[i] <  oracle->ncons);
      
      oracle->conlhs[indices[i]] = lhs != NULL ? lhs[i] : -SCIPinfinity(scip);
      oracle->conrhs[indices[i]] = rhs != NULL ? rhs[i] :  SCIPinfinity(scip);
   }
   
   return SCIP_OKAY;
}

/** deletes a set of variables */
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstat             /**< deletion status of vars in input (1 if var should be deleted, 0 if not); new position of var in output (-1 if var was deleted) */
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   /* TODO */
   SCIPerrorMessage("SCIPnlpiOracleDelVarSet is not implemented\n");
   
   return SCIP_ERROR;
}

/** deletes a set of constraints */
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstat             /**< deletion status of rows in input (1 if row should be deleted, 0 if not); new position of row in output (-1 if row was deleted) */
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
   assert(oracle != NULL);

   SCIPnlpiOracleInvalidateJacobiSparsity(scip, oracle);
   SCIPnlpiOracleInvalidateHessianLagSparsity(scip, oracle);

   /* TODO */
   SCIPerrorMessage("SCIPnlpiOracleDelConsSet is not implemented\n");
   
   return SCIP_ERROR;
}

/** changes linear coefficients in one constraint or objective
 * @return error if coefficient did not exist before
 */
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   cons,               /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidx,             /**< indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoeff            /**< new coefficients of variables */
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   /* @TODO */
   SCIPerrorMessage("SCIPnlpiOracleChgLinearCoefs is not implemented\n");
   
   return SCIP_ERROR;
}
  
/** changes coefficients in the quadratic part of one constraint or objective
 * @return error if coefficient did not exist before
 */
SCIP_RETCODE SCIPnlpiOracleChgQuadCoefs(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   cons,               /**< index of constraint where quadratic coefficients should be changed, or -1 for objective */
   const int             nentries,           /**< number of coefficients to change */
   const int*            rowoffset,          /**< row offset containing modified indices */
   const int*            colidx,             /**< columns containint modified indices to the corresponding row offset */
   SCIP_Real*            newcoeff            /**< new quadratic coefficients */ 
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
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
   
   return oracle->ncons;
}

/** gives the variables lower bounds */
const SCIP_Real* SCIPnlpiOracleGetVarLb(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
)
{
   assert(oracle != NULL);
   
   return oracle->varlb;
}

/** gives the variables upper bounds */
const SCIP_Real* SCIPnlpiOracleGetVarUb(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
)
{
   assert(oracle != NULL);
   
   return oracle->varub;
}

/** Gives maximum degree of a variable w.r.t. objective and all constraints.
 * The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
int SCIPnlpiOracleGetVarDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx
)
{
   assert(oracle != NULL);
   assert(varidx >= 0);
   assert(varidx < oracle->nvars);
   
   return oracle->vardegree[varidx];
}

/** Gives maximum degree of all variables w.r.t. objective and all constraints.
 * The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
int* SCIPnlpiOracleGetVarsDegree(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
)
{
   assert(oracle != NULL);

   return oracle->vardegree;
}

/** gives the constraints left-hand sides */
const SCIP_Real* SCIPnlpiOracleGetConstraintsLhs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
)
{
   assert(oracle != NULL);
   
   return oracle->conlhs;
}

/** gives the constraints right-hand sides */
const SCIP_Real* SCIPnlpiOracleGetConstraintsRhs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
)
{
   assert(oracle != NULL);
   
   return oracle->conrhs;
}

/** Gives maximum degree of a constraints.
 * The degree of a constraint is the maximal degree of all summands which appear in it, and is infinity for nonpolynomial terms.
 */ 
int SCIPnlpiOracleGetConstraintDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   conidx              /**< index of constraint for which the degree is requested */
)
{
   assert(oracle != NULL);
   assert(conidx >= 0);
   assert(conidx < oracle->ncons);
   
   if( oracle->conexprtree && oracle->conexprtree[conidx] )
      return INT_MAX;  /* @TODO something more clever, e.g., use exprtreeGetMaxDegree*/
   
   if( oracle->conquadlen && oracle->conquadlen[conidx] )
      return 2;

   if( oracle->conlinoffset && oracle->conlinoffset[conidx] )
      return 1;

   return 0;
}

/** computes the value of a function */
static
SCIP_RETCODE SCIPnlpiOracleEvalFunctionValue(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nlin,               /**< length of linear part */
   const int*            linind,             /**< indices of linear variables */
   const SCIP_Real*      linval,             /**< coefficients of linear variables */
   int                   quadlen,            /**< length of quadratic part matrix */
   int*                  quadrow,            /**< indices of    rows in quadratic part matrix */
   int*                  quadcol,            /**< indices of columns in quadratic part matrix */
   SCIP_Real*            quadval,            /**< coefficients in quadratic part matrix */
   int*                  exprvaridx,         /**< indices of variables in nonquadratic part */
   SCIP_EXPRTREE*        exprtree,           /**< nonquadratic part */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Real*            val                 /**< buffer to store function value */
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
   assert(oracle != NULL);
   assert(x      != NULL || oracle->nvars == 0);
   assert(val    != NULL);
   
   *val = 0.0;
   
   if( nlin != 0 )
   {
      assert(linind != NULL);
      assert(linval != NULL);
      assert(x != NULL); /* for lint */
      
      for( ; nlin; --nlin, ++linind, ++linval )
         *val += *linval * x[*linind]; 
   }
   
   if( quadlen != 0 )
   {
      assert(quadrow != NULL);
      assert(quadcol != NULL);
      assert(quadval != NULL);
      assert(x       != NULL); /* for lint */
      
      for( ; quadlen; --quadlen, ++quadrow, ++quadcol, ++quadval )
         *val += *quadval * x[*quadrow] * x[*quadcol];
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
   const int*            linind,             /**< indices of linear variables */
   const SCIP_Real*      linval,             /**< coefficients of linear variables */
   int                   quadlen,            /**< length of quadratic part matrix */
   int*                  quadrow,            /**< indices of    rows in quadratic part matrix */
   int*                  quadcol,            /**< indices of columns in quadratic part matrix */
   SCIP_Real*            quadval,            /**< coefficients in quadratic part matrix */
   int*                  exprvaridx,         /**< indices of variables in nonquadratic part */
   SCIP_EXPRTREE*        exprtree,           /**< nonquadratic part */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Bool             new_x,              /**< indicates whether the function has not been evaluated for this point before */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            grad                /**< buffer to store function gradient */
)
{  /*lint --e{715}*/
   assert(scip   != NULL);
   assert(oracle != NULL);
   assert(x      != NULL || oracle->nvars == 0);
   assert(val    != NULL);
   assert(grad   != NULL);
   
   *val = 0.0;
   BMSclearMemoryArray(grad, oracle->nvars);
   
   if( nlin != 0 )
   {
      assert(linind != NULL);
      assert(linval != NULL);
      assert(x      != NULL); /* for lint */
      
      for( ; nlin; --nlin, ++linind, ++linval )
      {
         *val += *linval * x[*linind];
         assert(grad[*linind] == 0.0);   /* we do not like duplicate indices */
         grad[*linind] = *linval;
      }
   }
   
   if( quadlen != 0 )
   {
      SCIP_Real tmp;
      assert(quadrow != NULL);
      assert(quadcol != NULL);
      assert(quadval != NULL);
      assert(x       != NULL); /* for lint */
      
      for( ; quadlen; --quadlen, ++quadrow, ++quadcol, ++quadval )
      {
         tmp = *quadval * x[*quadrow];
         *val += tmp * x[*quadcol];
         grad[*quadcol] += tmp;
         grad[*quadrow] += *quadval * x[*quadcol];
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
   SCIP_Real*            objval              /**< pointer to buffer to store objective value */  
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionValue(scip, oracle,
      oracle->objnlin, oracle->objlinind, oracle->objlinval,
      oracle->objquadlen, oracle->objquadrow, oracle->objquadcol, oracle->objquadval,
      oracle->objexprvaridx, oracle->objexprtree,
      x, objval) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** evaluates one constraint function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             conidx,             /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to buffer to store constraint value */  
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   assert(x      != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionValue(scip, oracle,
      oracle->conlinoffset != NULL ? (oracle->conlinoffset[conidx+1] - oracle->conlinoffset[conidx]) : 0,
      oracle->conlinoffset != NULL ? &oracle->conlinind[oracle->conlinoffset[conidx]] : NULL,
      oracle->conlinoffset != NULL ? &oracle->conlinval[oracle->conlinoffset[conidx]] : NULL,
      oracle->conquadlen   != NULL ? oracle->conquadlen[conidx]    : 0,
      oracle->conquadlen   != NULL ? oracle->conquadrow[conidx]    : NULL,
      oracle->conquadlen   != NULL ? oracle->conquadcol[conidx]    : NULL,
      oracle->conquadlen   != NULL ? oracle->conquadval[conidx]    : NULL,
      oracle->conexprtree  != NULL ? oracle->conexprvaridx[conidx] : NULL,
      oracle->conexprtree  != NULL ? oracle->conexprtree[conidx]   : NULL,
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
   
   assert(scip    != NULL);
   assert(oracle  != NULL);
   assert(x       != NULL || oracle->nvars == 0);
   assert(convals != NULL);

   for( i = 0; i < oracle->ncons; ++i )
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
   SCIP_Bool             new_x,              /**< indicates whether the function has not been evaluated for this point before */
   SCIP_Real*            objval,             /**< pointer to buffer to store objective value */
   SCIP_Real*            objgrad             /**< pointer to buffer to store (dense) objective gradient */  
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionGradient(scip, oracle,
      oracle->objnlin, oracle->objlinind, oracle->objlinval,
      oracle->objquadlen, oracle->objquadrow, oracle->objquadcol, oracle->objquadval,
      oracle->objexprvaridx, oracle->objexprtree,
      x, new_x, objval, objgrad) );
   
   *objval += oracle->objconstant;
   
   return SCIP_OKAY;
}

/** computes a constraints gradient in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             conidx,             /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             new_x,              /**< indicates whether the function has not been evaluated for this point before */ 
   SCIP_Real*            conval,             /**< pointer to buffer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to buffer to store (dense) constraint gradient */  
)
{
   assert(scip   != NULL);
   assert(oracle != NULL);
   assert(x      != NULL || oracle->nvars == 0);
   assert(conval != NULL);
   
   SCIP_CALL( SCIPnlpiOracleEvalFunctionGradient(scip, oracle,
      oracle->conlinoffset != NULL ? (oracle->conlinoffset[conidx+1] - oracle->conlinoffset[conidx]) : 0,
      oracle->conlinoffset != NULL ? &oracle->conlinind[oracle->conlinoffset[conidx]] : NULL,
      oracle->conlinoffset != NULL ? &oracle->conlinval[oracle->conlinoffset[conidx]] : NULL,
      oracle->conquadlen   != NULL ? oracle->conquadlen[conidx]    : 0,
      oracle->conquadlen   != NULL ? oracle->conquadrow[conidx]    : NULL,
      oracle->conquadlen   != NULL ? oracle->conquadcol[conidx]    : NULL,
      oracle->conquadlen   != NULL ? oracle->conquadval[conidx]    : NULL,
      oracle->conexprtree  != NULL ? oracle->conexprvaridx[conidx] : NULL,
      oracle->conexprtree  != NULL ? oracle->conexprtree[conidx]   : NULL,
      x, new_x, conval, congrad) );
   
   return SCIP_OKAY;
}

/** Gets sparsity pattern (rowwise) of Jacobian matrix.
 * 
 * Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 * Adding or deleting constraints destroys the sparsity structure and make another call to this function necessary. 
 */
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[ncons] gives length of col, can be NULL */
)
{
   int nnz, maxnnz, i, j;
   SCIP_Bool* nzflag;
   
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( oracle->jacoffset != NULL )
   {
      assert(oracle->jaccol != NULL);
      if( offset != NULL )
         *offset = oracle->jacoffset;
      if( col    != NULL )
         *col    = oracle->jaccol;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->jacoffset, oracle->ncons + 1) );

   maxnnz = MIN(oracle->nvars, 10) * oracle->ncons;  /* initial guess */
   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->jaccol, maxnnz) );
   
   if( maxnnz == 0 )
   {  /* no variables */
      BMSclearMemoryArray(oracle->jacoffset, oracle->ncons + 1);
      if( offset != NULL )
         *offset = oracle->jacoffset;
      if( col    != NULL )
         *col    = oracle->jaccol;
      return SCIP_OKAY;
   }
   nnz = 0;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &nzflag, oracle->nvars) );

   for( i = 0; i < oracle->ncons; ++i )
   {
      oracle->jacoffset[i] = nnz;
      
      if( oracle->conlinoffset != NULL && (oracle->conquadlen == NULL || oracle->conquadlen[i] == 0) && (oracle->conexprtree == NULL || oracle->conexprtree[i] == NULL) )
      { /* linear constraint: since we just want to copy the conlinval at EvalJacobian, we need to copy conlinind here too (it can be unsorted, and sorting would confuse ChgLinearCoef later) */
         int nz = oracle->conlinoffset[i+1] - oracle->conlinoffset[i];
         if( nnz + nz > maxnnz )
         {
           maxnnz *= 2;
           SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->jaccol, maxnnz) );
         }
         BMScopyMemoryArray(&oracle->jaccol[nnz], &oracle->conlinind[oracle->conlinoffset[i]], nz);
         nnz += nz;
         continue;
      }
      
      /* check which variables appear in constraint i */
      BMSclearMemoryArray(nzflag, oracle->nvars);
      
      if( oracle->conlinoffset != NULL )
      {
         assert(oracle->conlinind);
         for( j = oracle->conlinoffset[i]; j < oracle->conlinoffset[i+1]; ++j )
            nzflag[oracle->conlinind[j]] = TRUE;
      }
      
      if( oracle->conquadlen != NULL && oracle->conquadlen[i] != 0 )
      {
         assert(oracle->conquadrow    != NULL );
         assert(oracle->conquadrow[i] != NULL );
         assert(oracle->conquadcol    != NULL );
         assert(oracle->conquadcol[i] != NULL );
         for( j = 0; j < oracle->conquadlen[i]; ++j )
         {
            nzflag[oracle->conquadrow[i][j]] = TRUE;
            nzflag[oracle->conquadcol[i][j]] = TRUE;
         }
      }
      
      if( oracle->conexprtree != NULL && oracle->conexprtree[i] != NULL )
      {
         assert(oracle->conexprvaridx    != NULL );
         assert(oracle->conexprvaridx[i] != NULL );
         
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
      }
      
      /* store variables indices in jaccol */
      for( j = 0; j < oracle->nvars; ++j )
      {
         if( nzflag[j] == FALSE )
            continue;
         
         if( nnz >= maxnnz )
         {
            maxnnz *= 2;
            SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->jaccol, maxnnz) );
         }
         
         oracle->jaccol[nnz] = j;
         ++nnz;
      }
   }
   
   oracle->jacoffset[oracle->ncons] = nnz;
   
   if( nnz < maxnnz )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &oracle->jaccol, nnz) );
   }
   
   SCIPfreeBufferArray(scip, &nzflag);

   if( offset != NULL )
      *offset = oracle->jacoffset;
   if( col    != NULL )
      *col    = oracle->jaccol;

   return SCIP_OKAY;
}

/** Evaluates the Jacobi matrix in a given point.
 * 
 * The values in the Jacobi matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity.
 * The user need to call SCIPnlpiOracleGetJacobianSparsity at least ones before using this function. 
 */
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             new_x,              /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real*            convals,            /**< pointer to buffer to store constraint values, can be NULL */ 
   SCIP_Real*            jacobi              /**< pointer to buffer to store sparse jacobian values */  
)
{
   SCIP_Real* grad;
   int i, j, k, l;
   SCIP_Real dummy;
   
   assert(scip   != NULL);
   assert(oracle != NULL);
   assert(jacobi != NULL);
   
   assert(oracle->jacoffset != NULL);
   assert(oracle->jaccol    != NULL);
   
   SCIP_CALL( SCIPallocBufferArray(scip, &grad, oracle->nvars) );
   
   j = oracle->jacoffset[0];
   k = 0;
   for( i = 0; i < oracle->ncons; ++i )
   {
      if( oracle->conlinoffset != NULL && (oracle->conquadlen == NULL || oracle->conquadlen[i] == 0) && (oracle->conexprtree == NULL || oracle->conexprtree[i] == NULL) )
      { /* linear constraints should be easy */
         l = oracle->jacoffset[i+1] - j;
         assert(l == oracle->conlinoffset[i+1] - oracle->conlinoffset[i]);
         BMScopyMemoryArray(&jacobi[k], &oracle->conlinval[oracle->conlinoffset[i]], l);
         j += l; k += l;
      }
      else
      { /* @TODO do this sparse too */
         SCIP_RETCODE retcode = SCIPnlpiOracleEvalConstraintGradient(scip, oracle, i, x, new_x, convals ? &convals[i] : &dummy, grad);
         if( retcode != SCIP_OKAY )
         {
            SCIPfreeBufferArray(scip, &grad);
            return retcode;
         }
      
         for( ; j < oracle->jacoffset[i+1]; ++j, ++k )
            jacobi[k] = grad[oracle->jaccol[j]];
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
   assert(scip    != NULL);
   assert(nzflag  != NULL);
   assert(nzcount != NULL);
   assert(rowidx  != NULL);
   assert(colidx  != NULL);
   assert(length  >= 0);

   for( ; length; --length, ++rowidx, ++colidx )
   {
      assert(*colidx <= *rowidx);
      
      i = (*rowidx * (*rowidx+1)) / 2 + *colidx;
      if( nzflag[i] == FALSE )
      {
         nzflag[i] = TRUE;
         ++*nzcount;
      }
   }
}

/** Gets sparsity pattern of the Hessian matrix of the Lagrangian.
 * 
 * Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 * Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 * Only elements of the lower left triangle and the diagonal are counted.
 */
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[ncons] gives length of col, can be NULL */
)
{
   int nnz, i, j, cnt;
   SCIP_Bool* nzflag;
   SCIP_Bool* nz;
   
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( oracle->heslagoffset != NULL )
   {
      assert(oracle->heslagcol != NULL);
      if( offset != NULL )
         *offset = oracle->heslagoffset;
      if( col != NULL )
         *col    = oracle->heslagcol;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->heslagoffset, oracle->nvars + 1) );
   nnz = 0;
   
   /* @TODO should be improved */
   SCIP_CALL( SCIPallocBufferArray(scip, &nzflag, (oracle->nvars * (oracle->nvars + 1)) / 2) );
   BMSclearMemoryArray(nzflag, (oracle->nvars * (oracle->nvars + 1)) / 2);
   
   if( oracle->objquadlen != 0 )
      SCIPnlpiOracleHessLagSparsitySetNzFlagForQuad(scip, nzflag, &nnz, oracle->objquadrow, oracle->objquadcol, oracle->objquadlen);

   if( oracle->objexprtree != NULL )
   {
      SCIPerrorMessage("expression tree support not available\n");
      return SCIP_ERROR;
   }

   for( i = 0; i < oracle->ncons; ++i )
   {
      if( oracle->conquadlen != NULL && oracle->conquadlen[i] != 0 )
         SCIPnlpiOracleHessLagSparsitySetNzFlagForQuad(scip, nzflag, &nnz, oracle->conquadrow[i], oracle->conquadcol[i], oracle->conquadlen[i]);
      
      if( oracle->conexprtree != NULL && oracle->conexprtree[i] != NULL )
      {
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
      }
   }
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &oracle->heslagcol, nnz) );
   
   /* set hessian sparsity from nzflag */
   cnt = 0;
   nz = nzflag;
   for( i = 0; i < oracle->nvars; ++i )
   {
      oracle->heslagoffset[i] = cnt;
      for( j = 0; j <= i; ++j, ++nz )
         if( *nz == TRUE )
         {
            assert(cnt < nnz);
            oracle->heslagcol[cnt++] = j;
         }
   }
   oracle->heslagoffset[oracle->nvars] = cnt;
   assert(cnt == nnz);
   
   SCIPfreeBufferArray(scip, &nzflag);
   
   if( offset != NULL )
      *offset = oracle->heslagoffset;
   if( col    != NULL )
      *col    = oracle->heslagcol;

   return SCIP_OKAY;
}

/** adds quadratic part of a constraint into hessian structure */
static
SCIP_RETCODE SCIPnlpiOracleHessLagAddQuad(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_Real             weight,             /**< weight of quadratic part */
   int*                  rowidx,             /**< row    indices in matrix of quadratic part */
   int*                  colidx,             /**< column indices in matrix of quadratic part */
   SCIP_Real*            coeff,              /**< coefficients in matrix of quadratic part */
   int                   length,             /**< number of coefficients */
   int*                  hesoffset,          /**< row offsets in sparse matrix that is to be filled */ 
   int*                  hescol,             /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values              /**< buffer for values of sparse matrix that is to be filled */
)
{
   int idx;
   assert(scip      != NULL);
   assert(rowidx    != NULL);
   assert(colidx    != NULL);
   assert(length    >= 0);
   assert(hesoffset != NULL);
   assert(hescol    != NULL);
   assert(values    != NULL);

   for( ; length; --length, ++rowidx, ++colidx, ++coeff )
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
 * The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity.
 * The user must call SCIPnlpiOracleGetHessianLagSparsity at least ones before using this function. 
 * Only elements of the lower left triangle and the diagonal are computed.
 */
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP*                 scip,               /**< pointer to SCIP */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             new_x,              /**< indicates whether some function has not been evaluated for this point before */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambda,             /**< weights (Lagrangian multipliers) for the constraints */ 
   SCIP_Real*            hessian             /**< pointer to buffer to store sparse hessian values */  
)
{  /*lint --e{715}*/
   int i;
   assert(scip    != NULL);
   assert(oracle  != NULL);
   assert(x       != NULL);
   assert(lambda  != NULL);
   assert(hessian != NULL);

   assert(oracle->heslagoffset != NULL);
   assert(oracle->heslagcol    != NULL);
   
   for( i = oracle->heslagoffset[oracle->nvars]-1; i >= 0; --i )
      hessian[i] = 0.0;
   
   if( objfactor != 0.0 )
   {
      if( oracle->objquadlen != 0 )
      {
         assert(oracle->objquadrow != NULL);
         assert(oracle->objquadcol != NULL);
         assert(oracle->objquadval != NULL);
         SCIP_CALL( SCIPnlpiOracleHessLagAddQuad(scip, objfactor, oracle->objquadrow, oracle->objquadcol, oracle->objquadval, oracle->objquadlen, oracle->heslagoffset, oracle->heslagcol, hessian) );
      }
      
      if( oracle->objexprtree != NULL )
      {
         assert(oracle->objexprvaridx != NULL );
         SCIPerrorMessage("expression tree support not available\n");
         return SCIP_ERROR;
      }
   }
   
   for( i = 0; i < oracle->ncons; ++i )
   {
      if( lambda[i] == 0.0 )
         continue;
      if( oracle->conquadlen != NULL && oracle->conquadlen[i] != 0 )
      {
         assert(oracle->conquadrow    != NULL);
         assert(oracle->conquadrow[i] != NULL);
         assert(oracle->conquadcol    != NULL);
         assert(oracle->conquadcol[i] != NULL);
         assert(oracle->conquadval    != NULL);
         assert(oracle->conquadval[i] != NULL);
         SCIP_CALL( SCIPnlpiOracleHessLagAddQuad(scip, lambda[i], oracle->conquadrow[i], oracle->conquadcol[i], oracle->conquadval[i], oracle->conquadlen[i], oracle->heslagoffset, oracle->heslagcol, hessian) );
      }
      
      if( oracle->conexprtree != NULL && oracle->conexprtree[i] != NULL )
      {
         assert(oracle->conexprvaridx    != NULL);
         assert(oracle->conexprvaridx[i] != NULL);
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
   const int*            linind,             /**< indices of linear variables */
   const SCIP_Real*      linval,             /**< coefficients of linear variables */
   int                   quadlen,            /**< number of indices  in quadratic part matrix */
   int*                  quadrow,            /**< indices of rows    in quadratic part matrix */
   int*                  quadcol,            /**< indices of columns in quadratic part matrix */
   SCIP_Real*            quadval,            /**< coefficients in quadratic part matrix */
   SCIP_EXPRTREE*        exprtree            /**< nonquadratic part */
)
{  /*lint --e{715}*/
   int i, j;
   assert(scip   != NULL);
   assert(oracle != NULL);
   assert(file   != NULL);
   
   for( i = 0; i < nlin; ++i )
   {
      fprintf(file, "%+g*", linval[i]);
      if( oracle->varname != NULL && oracle->varname[linind[i]] != NULL )
         fprintf(file, oracle->varname[linind[i]]);
      else
         fprintf(file, "x%d", linind[i]);
   }
   
   if( quadlen != 0 )
   {
      j = 0;
      for( i = 0; i < quadlen; ++i )
      {
         fprintf(file, "%+g*", quadval[j]);
         if( oracle->varname != NULL && oracle->varname[quadrow[i]] != NULL )
            fprintf(file, oracle->varname[quadrow[i]]);
         else
            fprintf(file, "x%d", quadrow[i]);
         fprintf(file, "*");
         if( oracle->varname != NULL && oracle->varname[quadcol[i]] != NULL )
            fprintf(file, oracle->varname[quadcol[i]]);
         else
            fprintf(file, "x%d", quadcol[i]);
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
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( file == NULL )
      file = stdout;
   
   fprintf(file, "NLPI Oracle: %d variables and %d constraints\n", oracle->nvars, oracle->ncons);
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varname != NULL && oracle->varname[i] != NULL )
         fprintf(file, "%10s", oracle->varname[i]);
      else
         fprintf(file, "x%09d", i);
      fprintf(file, ": [%8g, %8g]\n", oracle->varlb[i], oracle->varub[i]);
   }
   
   fprintf(file, "objective: ");
   SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
      oracle->objnlin, oracle->objlinind, oracle->objlinval,
      oracle->objquadlen, oracle->objquadrow, oracle->objquadcol, oracle->objquadval,
      oracle->objexprtree) );
   if( oracle->objconstant != 0.0 )
      fprintf(file, "%+g", oracle->objconstant);
   fprintf(file, "\n");
   
   for( i = 0; i < oracle->ncons; ++i )
   {
      if( oracle->conname != NULL && oracle->conname[i] != NULL )
         fprintf(file, "%10s", oracle->conname[i]);
      else
         fprintf(file, "con%07d", i);
      
      fprintf(file, ": ");
      if( !SCIPisInfinity(scip, -oracle->conlhs[i]) && !SCIPisInfinity(scip, oracle->conrhs[i]) && oracle->conlhs[i] != oracle->conrhs[i] )
         fprintf(file, "%g <= ", oracle->conlhs[i]);
      
      SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
         oracle->conlinoffset != NULL ? (oracle->conlinoffset[i+1] - oracle->conlinoffset[i]) : 0,
         oracle->conlinoffset != NULL ? &oracle->conlinind[oracle->conlinoffset[i]] : NULL,
         oracle->conlinoffset != NULL ? &oracle->conlinval[oracle->conlinoffset[i]] : NULL,
         oracle->conquadlen   != NULL ? oracle->conquadlen[i]  : 0,
         oracle->conquadlen   != NULL ? oracle->conquadrow[i]  : NULL,
         oracle->conquadlen   != NULL ? oracle->conquadcol[i]  : NULL,
         oracle->conquadlen   != NULL ? oracle->conquadval[i]  : NULL,
         oracle->conexprtree  != NULL ? oracle->conexprtree[i] : NULL) );
      
      if( oracle->conlhs[i] == oracle->conrhs[i] )
         fprintf(file, " = %g", oracle->conrhs[i]);
      else if( !SCIPisInfinity(scip, oracle->conrhs[i]) )
         fprintf(file, " <= %g", oracle->conrhs[i]);
      else if( !SCIPisInfinity(scip, -oracle->conlhs[i]) )
         fprintf(file, " >= %g", oracle->conlhs[i]);
      
      fprintf(file, "\n");
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
   assert(scip   != NULL);
   assert(oracle != NULL);
   
   if( file == NULL )
      file = stdout;
   
   fprintf(file, "$offlisting\n");
   fprintf(file, "* NLPI Oracle Problem\n");
   fprintf(file, "Variables ");
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varname != NULL && oracle->varname[i] != NULL )
         fprintf(file, "%s, ", oracle->varname[i]);
      else
         fprintf(file, "x%d, ", i);
   }
   fprintf(file, "OBJVAR;\n\n");
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varlb[i] == oracle->varub[i] )
      {
         if( oracle->varname != NULL && oracle->varname[i] != NULL )
            fprintf(file, "%s", oracle->varname[i]);
         else
            fprintf(file, "x%d", i);
         fprintf(file, ".fx = %g;\t", oracle->varlb[i]);
      }
      else
      {
         if( !SCIPisInfinity(scip, -oracle->varlb[i]) )
         {
            if( oracle->varname != NULL && oracle->varname[i] != NULL )
               fprintf(file, "%s", oracle->varname[i]);
            else
               fprintf(file, "x%d", i);
            fprintf(file, ".lo = %g;\t", oracle->varlb[i]);
         }
         if( !SCIPisInfinity(scip, oracle->varub[i]) )
         {
            if( oracle->varname != NULL && oracle->varname[i] != NULL )
               fprintf(file, "%s", oracle->varname[i]);
            else
               fprintf(file, "x%d", i);
            fprintf(file, ".up = %g;\t", oracle->varub[i]);
         }
      }
      if( initval != NULL )
      {
         if( oracle->varname != NULL && oracle->varname[i] != NULL )
            fprintf(file, "%s", oracle->varname[i]);
         else
            fprintf(file, "x%d", i);
         fprintf(file, ".l = %g;\t", initval[i]);
      }
      fprintf(file, "\n");
   }
   fprintf(file, "\n");
   
   fprintf(file, "Equations ");
   for( i = 0; i < oracle->ncons; ++i )
   {
      if( oracle->conname != NULL && oracle->conname[i] != NULL )
         fprintf(file, "%s, ", oracle->conname[i]);
      else
         fprintf(file, "e%d, ", i);

      if( !SCIPisInfinity(scip, -oracle->conlhs[i]) && !SCIPisInfinity(scip, oracle->conrhs[i]) && oracle->conlhs[i] != oracle->conrhs[i] )
      { /* ranged row: add second constraint */
         if( oracle->conname != NULL && oracle->conname[i] != NULL )
            fprintf(file, "%s_RNG, ", oracle->conname[i]);
         else
            fprintf(file, "e%d_RNG, ", i);
      }
   }
   fprintf(file, "OBJ;\n\n");
   
   fprintf(file, "OBJ.. OBJVAR =E= ");
   SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
      oracle->objnlin, oracle->objlinind, oracle->objlinval,
      oracle->objquadlen, oracle->objquadrow, oracle->objquadcol, oracle->objquadval,
      oracle->objexprtree) );
   if( oracle->objconstant != 0.0 )
      fprintf(file, "%+g", oracle->objconstant);
   fprintf(file, ";\n");
   
   for( i = 0; i < oracle->ncons; ++i )
   {
      if( oracle->conname != NULL && oracle->conname[i] != NULL )
         fprintf(file, "%s", oracle->conname[i]);
      else
         fprintf(file, "e%d", i);
      fprintf(file, ".. ");
      
      SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
         oracle->conlinoffset != NULL ? (oracle->conlinoffset[i+1] - oracle->conlinoffset[i]) : 0,
         oracle->conlinoffset != NULL ? &oracle->conlinind[oracle->conlinoffset[i]] : NULL,
         oracle->conlinoffset != NULL ? &oracle->conlinval[oracle->conlinoffset[i]] : NULL,
         oracle->conquadlen   != NULL ? oracle->conquadlen[i]  : 0,
         oracle->conquadlen   != NULL ? oracle->conquadrow[i]  : NULL,
         oracle->conquadlen   != NULL ? oracle->conquadcol[i]  : NULL,
         oracle->conquadlen   != NULL ? oracle->conquadval[i]  : NULL,
         oracle->conexprtree  != NULL ? oracle->conexprtree[i] : NULL) );

      if( oracle->conlhs[i] == oracle->conrhs[i] )
         fprintf(file, " =E= %g", oracle->conrhs[i]);
      else if( !SCIPisInfinity(scip, oracle->conrhs[i]) )
         fprintf(file, " =L= %g", oracle->conrhs[i]);
      else if( !SCIPisInfinity(scip, -oracle->conlhs[i]) )
         fprintf(file, " =G= %g", oracle->conlhs[i]);
      else
         fprintf(file, " =N= 0");
      fprintf(file, ";\n");

      if( !SCIPisInfinity(scip, -oracle->conlhs[i]) && !SCIPisInfinity(scip, oracle->conrhs[i]) && oracle->conlhs[i] != oracle->conrhs[i] )
      {
         if( oracle->conname != NULL && oracle->conname[i] != NULL )
            fprintf(file, "%s", oracle->conname[i]);
         else
            fprintf(file, "e%d", i);
         fprintf(file, "_RNG.. ");
        
         SCIP_CALL( SCIPnlpiOraclePrintFunction(scip, oracle, file,
            oracle->conlinoffset != NULL ? (oracle->conlinoffset[i+1] - oracle->conlinoffset[i]) : 0,
            oracle->conlinoffset != NULL ? &oracle->conlinind[oracle->conlinoffset[i]] : NULL,
            oracle->conlinoffset != NULL ? &oracle->conlinval[oracle->conlinoffset[i]] : NULL,
            oracle->conquadlen   != NULL ? oracle->conquadlen[i]  : 0,
            oracle->conquadlen   != NULL ? oracle->conquadrow[i]  : NULL,
            oracle->conquadlen   != NULL ? oracle->conquadcol[i]  : NULL,
            oracle->conquadlen   != NULL ? oracle->conquadval[i]  : NULL,
            oracle->conexprtree  != NULL ? oracle->conexprtree[i] : NULL) );
         
         fprintf(file, " =G= %g;\n", oracle->conlhs[i]);
      }
   }
   
   fprintf(file, "Model m / all /;\n");
   fprintf(file, "option limrow = 0;\n");
   fprintf(file, "option limcol = 0;\n");
   fprintf(file, "Solve m minimizing OBJVAR using NLP;\n");

   return SCIP_OKAY;
}
