/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_implint.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  Presolver that detects implicit integer variables
 * @author Rolf van der Hulst
 */

/* TODO: explore integer to implicit integer conversion */
/* TODO: fix to work in MINLP context */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_implint.h"


#include "scip/pub_matrix.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_network.h"
#include "scip/pub_var.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_presol.h"
#include "scip/scip_pricer.h"
#include "scip/scip_probing.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"


#define PRESOL_NAME            "implint"
#define PRESOL_DESC            "detects implicit integer variables"
#define PRESOL_PRIORITY               100 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Data structures
 */

/* TODO: fill in the necessary presolver data */

/** presolver data */
struct SCIP_PresolData
{
   int temp; //TODO fix to stop complaining compiler
};


/*
 * Local methods
 */
typedef enum
{
   VAR_CLASS_INTEGER_CANDIDATE = 0,          /**< The variable has integrality constraints and could potentially
                                               *< be an implicit integer. */
   VAR_CLASS_INTEGER_FIXED = 1,              /**< All other variables with integrality constraints*/
   VAR_CLASS_CONTINUOUS_CANDIDATE = 2,       /**< The variable is continuous and could be an implicit integer */
   VAR_CLASS_CONTINUOUS_FIXED = 3            /**< All other continuous variables */
} VAR_CLASS;

/**
 * Struct that contains information about the blocks/components of the submatrix given by the continuous columns
 */
typedef struct{
   int nmatrixrows;                          /**< Number of rows in the matrix for the linear part of the problem */
   int nmatrixcols;                          /**< Number of columns in the matrix for the linear part of the problem */

   SCIP_VARTYPE * coltype;                   /**< SCIP_VARTYPE of the associated column */

   int* rowcomponent;                        /**< Maps a row to the index of the component it belongs to */
   int* colcomponent;                        /**< Maps a column to the index of the component it belongs to */

   int* componentrows;                       /**< Flattened array of array of rows that are in a given component. */
   int* componentcols;                       /**< Flattened array of array of columns that are in a given component. */
   int* componentrowstart;                   /**< The index of componentrows where the given component starts. */
   int* componentcolstart;                   /**< The index of componentcols where the given component starts. */
   int* ncomponentrows;                      /**< The number of rows in the given component. */
   int* ncomponentcols;                      /**< The number of columns in the given component */

   int ncomponents;
} MATRIX_COMPONENTS;

static
SCIP_RETCODE createMatrixComponents(
   SCIP* scip,
   SCIP_MATRIX* matrix,
   MATRIX_COMPONENTS** pmatrixcomponents
)
{
   SCIP_CALL( SCIPallocBlockMemory(scip, pmatrixcomponents) );
   MATRIX_COMPONENTS * comp = *pmatrixcomponents;

   int nrows = SCIPmatrixGetNRows(matrix);
   int ncols = SCIPmatrixGetNColumns(matrix);

   comp->nmatrixrows = nrows;
   comp->nmatrixcols = ncols;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->coltype, ncols) );
   for( int i = 0; i < ncols; ++i )
   {
      SCIP_VAR * var = SCIPmatrixGetVar(matrix,i);
      comp->coltype[i] = SCIPvarGetType(var);
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->rowcomponent, nrows) );
   for( int i = 0; i < nrows; ++i )
   {
      comp->rowcomponent[i] = -1;
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->colcomponent, ncols) );
   for( int i = 0; i < ncols; ++i )
   {
      comp->colcomponent[i] = -1;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->componentrows, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->componentcols, ncols) );
   //There will be at most ncols components
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->componentrowstart, ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->componentcolstart, ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->ncomponentrows, ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comp->ncomponentcols, ncols) );

   comp->ncomponents = 0;

   return SCIP_OKAY;
}

static
void freeMatrixInfo(
   SCIP* scip,
   MATRIX_COMPONENTS** pmatrixcomponents
)
{
   MATRIX_COMPONENTS* comp = *pmatrixcomponents;
   SCIPfreeBlockMemoryArray(scip, &comp->ncomponentcols, comp->nmatrixcols);
   SCIPfreeBlockMemoryArray(scip, &comp->ncomponentrows, comp->nmatrixcols);
   SCIPfreeBlockMemoryArray(scip, &comp->componentcolstart, comp->nmatrixcols);
   SCIPfreeBlockMemoryArray(scip, &comp->componentrowstart, comp->nmatrixcols);
   SCIPfreeBlockMemoryArray(scip, &comp->componentcols, comp->nmatrixcols);
   SCIPfreeBlockMemoryArray(scip, &comp->componentrows, comp->nmatrixrows);
   SCIPfreeBlockMemoryArray(scip, &comp->colcomponent, comp->nmatrixcols);
   SCIPfreeBlockMemoryArray(scip, &comp->rowcomponent, comp->nmatrixrows);
   SCIPfreeBlockMemoryArray(scip, &comp->coltype, comp->nmatrixcols);

   SCIPfreeBlockMemory(scip, pmatrixcomponents);
}

static
SCIP_RETCODE computeContinuousComponents(
   SCIP * scip,
   SCIP_MATRIX* matrix,
   MATRIX_COMPONENTS* comp
)
{
   int currentcomponent = 0;
   int rowindex = 0;
   int colindex = 0;
   /* We let rows and columns share an index by mapping column i to index nrows + i*/
   int* dfsstack = NULL;
   int ndfsstack = 0;
   SCIP_CALL(SCIPallocBufferArray(scip, &dfsstack, comp->nmatrixcols + comp->nmatrixrows));

   for( int i = 0; i < comp->nmatrixcols; ++i )
   {
      /* Check if earlier DFS already found column */
      if( comp->colcomponent[i] != -1 || comp->coltype[i] != SCIP_VARTYPE_CONTINUOUS){
         continue;
      }
      /* If not, then we have a new component, and perform a DFS to find all connected columns and rows */
      comp->componentrowstart[currentcomponent] = rowindex;
      comp->componentcolstart[currentcomponent] = colindex;

      assert(ndfsstack == 0);
      /* Push column to DFS search */
      dfsstack[ndfsstack] = comp->nmatrixrows + i;
      ++ndfsstack;
      comp->colcomponent[i] = currentcomponent;

      while( ndfsstack != 0 ){
         --ndfsstack;
         int index = dfsstack[ndfsstack];
         /* process column or row, adding new connected rows/columns to the dfs stack */
         if( index >= comp->nmatrixrows ){
            int column = index - comp->nmatrixrows;
            assert(comp->coltype[column] == SCIP_VARTYPE_CONTINUOUS);
            comp->componentcols[colindex] = column;
            ++colindex;

            /* Push connected rows to the dfs stack */
            int colnnonzs = SCIPmatrixGetColNNonzs(matrix, column);
            int* colrows = SCIPmatrixGetColIdxPtr(matrix, column);
            for( int j = 0; j < colnnonzs; ++j )
            {
               int row = colrows[j];
               if( comp->rowcomponent[row] == -1 ){
                  dfsstack[ndfsstack] = row;
                  ++ndfsstack;
                  comp->rowcomponent[row] = currentcomponent;
               }
            }
         }
         else
         {
            int row = index;
            comp->componentrows[rowindex] = row;
            ++rowindex;

            /* Push any connected continuous columns on the dfs search stack */
            int rownnonzs = SCIPmatrixGetRowNNonzs(matrix,row);
            int* rowcols = SCIPmatrixGetRowIdxPtr(matrix,row);
            for( int j = 0; j < rownnonzs; ++j )
            {
               int col = rowcols[j];
               if( comp->colcomponent[col] == -1 && comp->coltype[col] == SCIP_VARTYPE_CONTINUOUS){
                  comp->colcomponent[col] = currentcomponent;
                  dfsstack[ndfsstack] = comp->nmatrixrows + col;
                  ++ndfsstack;
               }
            }
         }
      }

      /* We are done with DFS for this component, save the sizes */
      comp->ncomponentrows[currentcomponent] = rowindex - comp->componentrowstart[currentcomponent];
      comp->ncomponentcols[currentcomponent] = colindex - comp->componentcolstart[currentcomponent];
      ++currentcomponent;
      ++comp->ncomponents;
      assert(currentcomponent == comp->ncomponents);
   }
   SCIPfreeBufferArray(scip,&dfsstack);
   return SCIP_OKAY;
}

typedef struct{
   SCIP_Bool* rowintegral;                   /**< Are all the non-continuous column entries and lhs rhs integral? */
   SCIP_Bool* rowequality;                   /**< Is the row an equality? */
   SCIP_Bool* rowbadnumerics;                /**< Does the row contain large entries that make numerics difficult? */
   int* rownnonz;                            /**< Number of nonzeros in the column */
   int* rowncontinuous;                      /**< The number of those nonzeros that are in continuous columns */
   int* rowncontinuouspmone;                 /**< The number of continuous columns +-1 entries */

} MATRIX_STATISTICS;

static
SCIP_RETCODE computeMatrixStatistics(
   SCIP * scip,
   SCIP_MATRIX* matrix,
   MATRIX_STATISTICS** pstats
)
{
   SCIP_CALL( SCIPallocBuffer(scip,pstats) );
   MATRIX_STATISTICS* stats = *pstats;

   int nrows = SCIPmatrixGetNRows(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip,&stats->rowintegral,nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip,&stats->rowequality,nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip,&stats->rowbadnumerics,nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip,&stats->rownnonz,nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip,&stats->rowncontinuous,nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip,&stats->rowncontinuouspmone,nrows) );


   for( int i = 0; i < nrows; ++i )
   {
      double lhs = SCIPmatrixGetRowLhs(matrix,i);
      double rhs = SCIPmatrixGetRowRhs(matrix,i);
      int * cols = SCIPmatrixGetRowIdxPtr(matrix,i);
      double * vals = SCIPmatrixGetRowValPtr(matrix,i);
      int nnonz = SCIPmatrixGetRowNNonzs(matrix,i);
      stats->rownnonz[i] = nnonz;
      stats->rowequality[i] = SCIPisFeasEQ(scip,lhs,rhs) && !( SCIPisInfinity(scip,-lhs) || SCIPisInfinity(scip, rhs) );

      SCIP_Bool integral = ( SCIPisInfinity(scip,-lhs) || SCIPisIntegral(scip,lhs)) &&
                           ( SCIPisInfinity(scip,rhs) || SCIPisIntegral(scip,rhs));
      SCIP_Bool badnumerics = FALSE;

      int ncontinuous = 0;
      int ncontinuouspmone = 0;
      for( int j = 0; j < nnonz; ++j )
      {
         SCIP_Bool continuous = SCIPvarGetType(SCIPmatrixGetVar(matrix,cols[j])) == SCIP_VARTYPE_CONTINUOUS;
         double value = vals[j];
         if(continuous){
            ++ncontinuous;
         }
         if(continuous && ABS(value) == 1.0){
            ++ncontinuouspmone;
         }
         if(!continuous){
            integral = integral && SCIPisIntegral(scip,value);
         }
         if(ABS(value) > 1e7){
            badnumerics = TRUE;
         }
      }

      stats->rowncontinuous[i] = ncontinuous;
      stats->rowncontinuouspmone[i] = ncontinuouspmone;
      stats->rowintegral[i] = integral;
      stats->rowbadnumerics[i] = badnumerics;
   }


   return SCIP_OKAY;
}

static
void freeMatrixStatistics(
   SCIP* scip,
   MATRIX_STATISTICS** pstats
)
{
   MATRIX_STATISTICS* stats= *pstats;
   SCIPfreeBufferArray(scip,&stats->rowncontinuouspmone);
   SCIPfreeBufferArray(scip,&stats->rowncontinuous);
   SCIPfreeBufferArray(scip,&stats->rownnonz);
   SCIPfreeBufferArray(scip,&stats->rowequality);
   SCIPfreeBufferArray(scip,&stats->rowintegral);
   SCIPfreeBufferArray(scip,&stats->rowbadnumerics);

   SCIPfreeBuffer(scip,pstats);
}

static
SCIP_RETCODE findImpliedIntegers(
   SCIP * scip,
   SCIP_MATRIX* matrix,
   MATRIX_COMPONENTS* comp,
   MATRIX_STATISTICS* stats,
   int* nchgvartypes
)
{
   //TODO: some checks to prevent expensive memory initialization if not necessary (e.g. there must be some candidates)
   SCIP_NETMATDEC * dec = NULL;
   SCIP_CALL(SCIPnetmatdecCreate(SCIPblkmem(scip),&dec,comp->nmatrixrows,comp->nmatrixcols));

   SCIP_NETMATDEC * transdec = NULL;
   SCIP_CALL(SCIPnetmatdecCreate(SCIPblkmem(scip),&transdec,comp->nmatrixcols,comp->nmatrixrows));

   int planarcomponents = 0;
   int goodcomponents = 0;
   int nbadnumerics = 0;
   int nbadintegrality = 0;
   int nnonnetwork = 0;

   /* Because the rows may also contain non-continuous columns, we need to remove these from the array that we
   * pass to the network matrix decomposition method. We use these working arrays for this purpose. */
   double* tempValArray;
   int* tempIdxArray;
   SCIP_CALL(SCIPallocBufferArray(scip,&tempValArray,comp->nmatrixcols));
   SCIP_CALL(SCIPallocBufferArray(scip,&tempIdxArray,comp->nmatrixcols));


   for( int component = 0; component < comp->ncomponents; ++component )
   {
      int startrow = comp->componentrowstart[component];
      int nrows = comp->ncomponentrows[component];
      SCIP_Bool componentokay = TRUE;
      for( int i = startrow; i < startrow + nrows; ++i )
      {
         int row = comp->componentrows[i];
         if(stats->rowncontinuous[row] != stats->rowncontinuouspmone[row]){
            componentokay = FALSE;
            ++nbadintegrality;
            break;
         }
         if(!stats->rowintegral[row]){
            componentokay = FALSE;
            ++nbadintegrality;
            break;
         }
         if(stats->rowbadnumerics[row]){
            componentokay = FALSE;
            ++nbadnumerics;
            break;
         }
      }
      if(!componentokay){
         continue;
      }
      int startcol = comp->componentcolstart[component];
      int ncols = comp->ncomponentcols[component];

      /* Check if the component is a network matrix */
      SCIP_Bool componentnetwork = TRUE;

      /* We use the row-wise algorithm only if the number of columns is much larger than the number of rows.
       * Generally, the column-wise algorithm will be faster, but in these extreme cases, the row algorithm is faster.
       * TODO: test/tune this parameter
       */
      //TODO for ~50-150 seems to be about even on neos-...-inde and neos-...-isar, test further
      if(nrows * 20 < ncols){
         printf("Row addition: %f \n",(double) ncols / (double) nrows);
         for( int i = startrow; i < startrow + nrows && componentnetwork; ++i )
         {
            int row = comp->componentrows[i];
            int nrownnoz = SCIPmatrixGetRowNNonzs(matrix,row);
            int* rowcols = SCIPmatrixGetRowIdxPtr(matrix,row);
            double* rowvals = SCIPmatrixGetRowValPtr(matrix,row);
            int ncontnonz = 0;
            for( int j = 0; j < nrownnoz; ++j )
            {
               int col = rowcols[j];
               if(SCIPvarGetType(SCIPmatrixGetVar(matrix,col)) == SCIP_VARTYPE_CONTINUOUS)
               {
                  tempIdxArray[ncontnonz] = col;
                  tempValArray[ncontnonz] = rowvals[j];
                  ++ncontnonz;
                  assert(ABS(rowvals[j]) == 1.0);
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddRow(dec,row,tempIdxArray,tempValArray,ncontnonz,&componentnetwork) );
         }
      }else{
         for( int i = startcol; i < startcol + ncols && componentnetwork; ++i )
         {
            int col = comp->componentcols[i];
            int ncolnnonz = SCIPmatrixGetColNNonzs(matrix,col);
            int* colrows = SCIPmatrixGetColIdxPtr(matrix,col);
            double* colvals = SCIPmatrixGetColValPtr(matrix,col);
            SCIP_CALL( SCIPnetmatdecTryAddCol(dec,col,colrows,colvals,ncolnnonz,&componentnetwork) );
         }
      }

      if( !componentnetwork )
      {
         SCIPnetmatdecRemoveComponent(dec,&comp->componentrows[startrow], nrows, &comp->componentcols[startcol], ncols);
         ++nnonnetwork;
      }

      SCIP_Bool componenttransnetwork = TRUE;

      /* For the transposed matrix, the situation is exactly reversed because the row/column algorithms are swapped */
      //TODO: tune parameter
      if(nrows < ncols * 20){
         for( int i = startrow; i < startrow + nrows && componenttransnetwork ; ++i )
         {
            int row = comp->componentrows[i];
            int nrownnoz = SCIPmatrixGetRowNNonzs(matrix,row);
            int* rowcols = SCIPmatrixGetRowIdxPtr(matrix,row);
            double* rowvals = SCIPmatrixGetRowValPtr(matrix,row);
            int ncontnonz = 0;
            for( int j = 0; j < nrownnoz; ++j )
            {
               int col = rowcols[j];
               if(SCIPvarGetType(SCIPmatrixGetVar(matrix,col)) == SCIP_VARTYPE_CONTINUOUS)
               {
                  tempIdxArray[ncontnonz] = col;
                  tempValArray[ncontnonz] = rowvals[j];
                  ++ncontnonz;
                  assert(ABS(rowvals[j]) == 1.0);
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddCol(transdec,row,tempIdxArray,tempValArray,ncontnonz,
                                              &componenttransnetwork) );
         }
      }else{
         for( int i = startcol; i < startcol + ncols && componenttransnetwork; ++i )
         {
            int col = comp->componentcols[i];
            int ncolnnonz = SCIPmatrixGetColNNonzs(matrix,col);
            int* colrows = SCIPmatrixGetColIdxPtr(matrix,col);
            double* colvals = SCIPmatrixGetColValPtr(matrix,col);
            SCIP_CALL( SCIPnetmatdecTryAddRow(transdec,col,colrows,colvals,ncolnnonz,&componenttransnetwork) );
         }
      }

      if( !componenttransnetwork )
      {
         SCIPnetmatdecRemoveComponent(transdec,&comp->componentcols[startcol],
                                      ncols, &comp->componentrows[startrow], nrows);
      }

      if( !componentnetwork && !componenttransnetwork){
         ++nnonnetwork;
         continue;
      }
      ++goodcomponents;
      if(componentnetwork && componenttransnetwork){
         ++planarcomponents;
      }
      for( int i = startcol; i < startcol + ncols; ++i )
      {
         int col = comp->componentcols[i];
         SCIP_VAR * var = SCIPmatrixGetVar(matrix,col);
         SCIP_Bool infeasible = FALSE;
         SCIP_CALL(SCIPchgVarType(scip,var,SCIP_VARTYPE_IMPLINT,&infeasible));
         (*nchgvartypes)++;
         //TODO: inform of changed variable types
         assert(!infeasible);
      }
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH,NULL,
                   "implied integer components: %d (%d planar) / %d (disqualified: %d by integrality, %d by numerics, %d not network) \n",
                   goodcomponents,planarcomponents,comp->ncomponents,nbadintegrality,nbadnumerics,nnonnetwork);

   SCIPfreeBufferArray(scip,&tempIdxArray);
   SCIPfreeBufferArray(scip,&tempValArray);

   SCIPnetmatdecFree(&dec);

   return SCIP_OKAY;
}
/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyImplint NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRESOLFREE(presolFreeImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolFreeImplint NULL
#endif


/** initialization method of presolver (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRESOLINIT(presolInitImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitImplint NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitImplint NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreImplint NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitpreImplint NULL
#endif


/** execution method of presolver */

static
SCIP_DECL_PRESOLEXEC(presolExecImplint)
{
   *result = SCIP_DIDNOTRUN;

   //TODO: re-check these conditions again
   //Disable implicit integer detection if we are probing or in NLP context
   if(( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING ) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip))
   {
      return SCIP_OKAY;
   }
   //Since implied integer detection relies on rows not being changed, we disable it for branch-and-price applications
   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
   {
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   double starttime = SCIPgetSolvingTime(scip);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                   "   (%.1fs) implied integer detection started\n", starttime);

   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   SCIP_MATRIX* matrix = NULL;
   SCIP_Bool onlyifcomplete = TRUE;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, onlyifcomplete, &initialized, &complete, &infeasible,
                              naddconss, ndelconss, nchgcoefs, nchgbds, nfixedvars) );
   /*If infeasibility was detected during matrix creation, we return. */
   if( infeasible )
   {
      if( initialized )
      {
         SCIPmatrixFree(scip, &matrix);
      }
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /*For now, we only work on pure MILP's TODO; use uplocks/downlocks */
   if( !( initialized && complete ))
   {
      if( initialized )
      {
         SCIPmatrixFree(scip, &matrix);
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "   (%.1fs) implied integer detection stopped because problem is not an MILP\n",
                      SCIPgetSolvingTime(scip));
      return SCIP_OKAY;
   }

   int beforechanged = *nchgvartypes;
   MATRIX_COMPONENTS* comp = NULL;
   MATRIX_STATISTICS* stats = NULL;
   SCIP_CALL( createMatrixComponents(scip, matrix, &comp) );
   SCIP_CALL( computeMatrixStatistics(scip, matrix, &stats) );
   SCIP_CALL( computeContinuousComponents(scip, matrix, comp) );
   SCIP_CALL( findImpliedIntegers(scip,matrix,comp, stats, nchgvartypes) );
   int afterchanged = *nchgvartypes;


   double endtime = SCIPgetSolvingTime(scip);
   if(afterchanged == beforechanged){
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "   (%.1fs) no implied integers detected (time: %.2fs)\n", endtime,endtime-starttime);
      *result = SCIP_DIDNOTFIND;
   }else{
      *result = SCIP_SUCCESS;
   }
   freeMatrixStatistics(scip,&stats);
   freeMatrixInfo(scip, &comp);
   SCIPmatrixFree(scip, &matrix);
   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the implint presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolImplint(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create implint presolver data */
   presoldata = NULL;
   /* TODO: (optional) create presolver specific data here */

   presol = NULL;

   /* use SCIPincludePresolBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
                                     PRESOL_TIMING, presolExecImplint, presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplint) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeImplint) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitImplint) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitImplint) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreImplint) );
   SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreImplint) );

   /* add implint presolver parameters */
   /* TODO: (optional) add presolver specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
