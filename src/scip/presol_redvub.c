/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_redvub.c
 * @brief  remove redundant variable upper bound constraints
 * @author Dieter Weninger
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"

#include "presol_redvub.h"

#define PRESOL_NAME            "redvub"
#define PRESOL_DESC            "remove redundant vubs"
#define PRESOL_PRIORITY         24000000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               FALSE     /**< should presolver be delayed, if other presolvers found reductions? */


/*
 * Local methods
 */

/** is the constraint a vub? */
static
SCIP_RETCODE isVubCons(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix instance */
   int                   row,                /**< row index */
   SCIP_Bool*            isvub               /**< flag indicating if constraint is a vub */
   )
{
   int* rowpnt;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < SCIPmatrixGetNRows(matrix));
   assert(isvub != NULL);

   *isvub = FALSE;

   if( SCIPmatrixGetRowNNonzs(matrix,row) == 2 )
   {
      SCIP_VARTYPE type1;
      SCIP_VARTYPE type2;
      int idx1;
      int idx2;
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real lhs;
      SCIP_Real rhs;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix,row);
      idx1 = *rowpnt;
      var1 = SCIPmatrixGetVar(matrix,idx1);
      type1 = SCIPvarGetType(var1);
      rowpnt++;
      idx2 = *rowpnt;
      var2 = SCIPmatrixGetVar(matrix,idx2);
      type2 = SCIPvarGetType(var2);
      lhs = SCIPmatrixGetRowLhs(matrix,row);
      rhs = SCIPmatrixGetRowRhs(matrix,row);

      if( ((type1 == SCIP_VARTYPE_CONTINUOUS && type2 == SCIP_VARTYPE_BINARY) ||
            (type2 == SCIP_VARTYPE_CONTINUOUS && type1 == SCIP_VARTYPE_BINARY)) &&
         (SCIPisInfinity(scip,-lhs) || SCIPisInfinity(scip,rhs)) )
      {
         *isvub = TRUE;
      }
   }

   return SCIP_OKAY;
}

/**< detect vub's on the same continuous variable */
static
SCIP_RETCODE detectParallelVubs(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   int                   nvubs,              /**< number of vub's */
   int*                  vubs,               /**< row indices of the vub's */
   int*                  parallelvubs        /**< parallel vub's concerning the continuous variable */
   )
{
   int i;
   int j;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real* scale;
   SCIP_Real* scaledx;
   SCIP_Real* scaledy;
   SCIP_Real* scaledc;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vubs != NULL);
   assert(parallelvubs != NULL);

   /* we assume the following form: c <= ax + by <= inf */

   assert(nvubs >= 2);

   SCIP_CALL( SCIPallocBufferArray(scip, &scale, nvubs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scaledx, nvubs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scaledy, nvubs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scaledc, nvubs) );

   for( i = 0; i < nvubs; i++ )
   {
      assert(SCIPmatrixGetRowNNonzs(matrix,vubs[i]) == 2);

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix,vubs[i]);
      rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix,vubs[i]);
      valpnt = SCIPmatrixGetRowValPtr(matrix,vubs[i]);

      for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         int col;
         SCIP_VAR* var;

         col = *rowpnt;
         var = SCIPmatrixGetVar(matrix,col);

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            scale[i] = *valpnt;
            scaledx[i] = 1.0;
         }
         else
         {
            scaledy[i] = *valpnt;
         }
      }

      scaledc[i] = SCIPmatrixGetRowLhs(matrix,vubs[i]);
   }

   for( i = 0; i < nvubs; i++ )
   {
      scaledy[i] /= scale[i];
      scaledc[i] /= scale[i];
   }

   for( i = 0; i < nvubs; i++ )
      parallelvubs[i] = i;

   for( i = 0; i < nvubs; i++ )
   {
      for( j = i+1; j < nvubs; j++ )
      {
         /* currently we only treat the case with the same coef on y and lhs=0 */
         if( SCIPisEQ(scip,scaledy[i],scaledy[j]) && SCIPisEQ(scip,scaledc[i],scaledc[j]) &&
            SCIPisEQ(scip,scaledc[i],0) &&
            ((SCIPisLT(scip,scale[i],0) && SCIPisLT(scip,scale[j],0)) || (SCIPisGT(scip,scale[i],0) && SCIPisGT(scip,scale[j],0))) )
         {
            parallelvubs[j] = parallelvubs[i];
         }
      }
   }

   SCIPfreeBufferArray(scip, &scaledc);
   SCIPfreeBufferArray(scip, &scaledy);
   SCIPfreeBufferArray(scip, &scaledx);
   SCIPfreeBufferArray(scip, &scale);

   return SCIP_OKAY;
}

/**< detect variables which should be substituted by other ones */
static
SCIP_RETCODE substitutevar(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix */
   int*                  parallelvubs,       /**< parallel vub's concerning the same continuous variable */
   int                   fill,               /**< number of parallel vub's */
   int*                  nvarsub,            /**< number of substituted variables */
   SCIP_Bool*            isvartosub,         /**< flag array if variable could be substituted */
   SCIP_VAR**            subvars,            /**< pointers to the variables by which the substitution should be done */
   int*                  ndeletecons,        /**< number of deleteable constraints */
   SCIP_Bool*            deletecons          /**< flags which constraints could be deleted */
   )
{
   int i;
   int holdcol;
   int* substitute;
   int maxsupport;
   int* rowpnt;
   int* rowend;

   SCIP_CALL( SCIPallocBufferArray(scip, &substitute, fill) );
   for( i = 0; i < fill; i++ )
      substitute[i] = -1;

   holdcol = -1;

   /* detect binary variable which should not be substituted */
   maxsupport = 0;
   for( i = 0; i < fill; i++ )
   {
      rowpnt = SCIPmatrixGetRowIdxPtr(matrix,parallelvubs[i]);
      rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix,parallelvubs[i]);

      for( ; (rowpnt < rowend); rowpnt++ )
      {
         int col;
         SCIP_VAR* var;

         col = *rowpnt;
         var = SCIPmatrixGetVar(matrix,col);

         if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         {
            substitute[i] = col;
            if( SCIPmatrixGetColNNonzs(matrix,col) > maxsupport &&
               SCIPmatrixGetColNDownlocks(matrix,col) == 1 &&
               SCIPisGE(scip,SCIPvarGetObj(var),0) )
            {
               maxsupport = SCIPmatrixGetColNNonzs(matrix,col);
               holdcol = col;
            }
         }
      }
   }

   /* substitute all binary variables by the hold binary variable
    * and delete the corresponding vub
    */
   assert(holdcol > -1);
   for( i = 0; i < fill; i++ )
   {
      assert(substitute[i] > -1);

      if( substitute[i] != holdcol &&
         SCIPmatrixGetColNDownlocks(matrix,substitute[i]) == 1 &&
         SCIPisGE(scip,SCIPvarGetObj(SCIPmatrixGetVar(matrix,substitute[i])),0) )
      {
         assert(isvartosub[substitute[i]] == FALSE);
         isvartosub[substitute[i]] = TRUE;
         subvars[substitute[i]] = SCIPmatrixGetVar(matrix,holdcol);
         (*nvarsub)++;

         assert(deletecons[parallelvubs[i]] == FALSE);
         deletecons[parallelvubs[i]] = TRUE;
         (*ndeletecons)++;
      }
   }

   SCIPfreeBufferArray(scip, &substitute);

   return SCIP_OKAY;
}

/**< get substitutable variables and redundant constraints */
static
SCIP_RETCODE getSubDelcons(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix */
   int*                  nvarsub,            /**< number of redundant variables */
   SCIP_Bool*            isvartosub,         /**< flags indicating which variables could be substituted */
   SCIP_VAR**            subvars,            /**< pointers to the variables by which the substitution should be done */
   int*                  ndeletecons,        /**< number of redundant constraints */
   SCIP_Bool*            deletecons          /**< flags which constraints could be deleted */
   )
{
   int c;
   int* colpnt;
   int* colend;
   int* vubcons;
   int nvubcons;
   int* pclass;
   int pc;
   int pclassstart;
   int fill;
   int* pcons;
   int ncols;
   int nrows;

   ncols = SCIPmatrixGetNColumns(matrix);
   nrows = SCIPmatrixGetNRows(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &vubcons, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pclass, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcons, nrows) );

   for( c = 0; c < ncols; c++ )
   {
      SCIP_VAR* var;

      var = SCIPmatrixGetVar(matrix,c);

      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         continue;

      /* search vubs per variable */
      nvubcons = 0;
      colpnt = SCIPmatrixGetColIdxPtr(matrix,c);
      colend = colpnt + SCIPmatrixGetColNNonzs(matrix,c);
      for( ; (colpnt < colend); colpnt++ )
      {
         SCIP_Bool isvub;

         SCIP_CALL( isVubCons(scip,matrix,*colpnt,&isvub) );
         if( isvub )
            vubcons[nvubcons++] = *colpnt;
      }

      if( nvubcons < 2 )
         continue;

      /* detect vubs working on the same continuous variable */
      SCIP_CALL( detectParallelVubs(scip,matrix,nvubcons,vubcons,pclass) );
      SCIPsortIntInt(pclass, vubcons, nvubcons);

      pc = 0;
      while( pc < nvubcons )
      {
         fill = 0;
         pclassstart = pclass[pc];
         while( pc < nvubcons && pclassstart == pclass[pc] )
         {
            pcons[fill++] = vubcons[pc];
            pc++;
         }

         if( fill > 1 )
         {
            SCIP_CALL( substitutevar(scip,matrix,pcons,
                  fill,nvarsub,isvartosub,subvars,ndeletecons,deletecons) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &pcons);
   SCIPfreeBufferArray(scip, &pclass);
   SCIPfreeBufferArray(scip, &vubcons);

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecRedvub)
{  /*lint --e{715}*/
   SCIPMILPMATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPgetNVars(scip) == 0 || SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( !complete )
      SCIPdebugMessage("Warning: milp matrix incomplete!\n");

   if( initialized )
   {
      int nvarsub;
      SCIP_Bool* isvartosub;
      int ndeletecons;
      SCIP_Bool* deletecons;
      SCIP_VAR** subvars;
      int ncols;
      int nrows;

      ncols = SCIPmatrixGetNColumns(matrix);
      nrows = SCIPmatrixGetNRows(matrix);
      nvarsub = 0;
      ndeletecons = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &isvartosub, ncols) );
      BMSclearMemoryArray(isvartosub, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &deletecons, ncols) );
      BMSclearMemoryArray(deletecons, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, ncols) );

      SCIP_CALL( getSubDelcons(scip,matrix,&nvarsub,isvartosub,subvars,&ndeletecons,deletecons) );

      if( nvarsub > 0 )
      {
         int v;
         for( v = 0; v < ncols; v++ )
         {
            if( isvartosub[v] )
            {
               SCIP_Bool infeasible;
               SCIP_Bool redundant;
               SCIP_Bool aggregated;

               /* substitute/aggregate variable */
               SCIP_CALL( SCIPaggregateVars(scip, SCIPmatrixGetVar(matrix,v), subvars[v], 1.0, -1.0,
                     0, &infeasible, &redundant, &aggregated) );

               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible aggregation\n");
                  return SCIP_OKAY;
               }

               if( aggregated )
                  (*naggrvars)++;
            }
         }
      }

      if( ndeletecons > 0 )
      {
         int r;
         for( r = 0; r < nrows; r++ )
         {
            if( deletecons[r] )
            {
               SCIP_CONS* cons;

               /* remove redundant constraint */
               cons = SCIPmatrixGetCons(matrix,r);
               SCIP_CALL( SCIPdelCons(scip, cons) );

               (*ndelconss)++;
            }
         }
      }

      SCIPfreeBufferArray(scip, &subvars);
      SCIPfreeBufferArray(scip, &deletecons);
      SCIPfreeBufferArray(scip, &isvartosub);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the redvub presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolRedvub(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_DELAY, presolExecRedvub, NULL) );

   return SCIP_OKAY;
}
