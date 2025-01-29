/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_presol.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"

#define PRESOL_NAME             "implint"
#define PRESOL_DESC             "detects implicit integer variables"
#define PRESOL_PRIORITY         100 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS        0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_COLUMNROWRATIO  50.0
#define DEFAULT_NUMERICSLIMIT   1e6

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Real             columnrowratio;     /**< Use the network row addition algorithm when the column to row ratio
                                               * becomes larger than this threshold. Otherwise, use column addition. */
   SCIP_Real             numericslimit;      /**< A row that contains variables with coefficients that are greater in
                                                * absolute value than this limit is not considered for
                                                * implied integrality detection. */
};


/** Struct that contains information about the blocks/components of the submatrix given by the continuous columns.
 *
 * Note that currently, the matrix represents exactly the SCIP_MATRIX created by SCIPmatrixCreate(), but this may change
 * once MINLP problems are also accounted for.
 * @todo extend the plugin to also work for MINLP problems. This changes the computed matrix.
 */
struct MatrixComponents
{
   int nmatrixrows;                          /**< Number of rows in the matrix for the linear part of the problem */
   int nmatrixcols;                          /**< Number of columns in the matrix for the linear part of the problem */

   SCIP_VARTYPE* coltype;                    /**< SCIP_VARTYPE of the associated column */

   int* rowcomponent;                        /**< Maps a row to the index of the component it belongs to */
   int* colcomponent;                        /**< Maps a column to the index of the component it belongs to */

   int* componentrows;                       /**< Flattened array of arrays of rows that are in a given component. */
   int* componentcols;                       /**< Flattened array of arrays of columns that are in a given component. */
   int* componentrowend;                     /**< The index of componentrows where the given component ends. */
   int* componentcolend;                     /**< The index of componentcols where the given component ends. */
   int ncomponents;                          /**< The number of components. */
};
typedef struct MatrixComponents MATRIX_COMPONENTS;

/** A temporary data structure that stores some statistics/data on the rows and columns.
 *
 * This is freed again after implied integer detection is finished.
 */
struct MatrixStatistics
{
   SCIP_Bool* rowintegral;                   /**< Are all row entries of non-continuous columns and the row sides integral? */
   SCIP_Bool* rowequality;                   /**< Is the row an equality? */
   SCIP_Bool* rowbadnumerics;                /**< Does the row contain large entries that make numerics difficult? */
   int* rownnonz;                            /**< Number of nonzeros in the row */
   int* rowncontinuous;                      /**< The number of those nonzeros that are in continuous columns */
   int* rowncontinuouspmone;                 /**< The number of +-1 entries in continuous columns */
   SCIP_Bool* colintegralbounds;             /**< Does the column have integral bounds? */
};
typedef struct MatrixStatistics MATRIX_STATISTICS;

/** Creates the matrix components data structure */
static
SCIP_RETCODE createMatrixComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< The constraint matrix */
   MATRIX_COMPONENTS**   pmatrixcomponents   /**< Pointer to create the matrix components data structure */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, pmatrixcomponents) );
   MATRIX_COMPONENTS* comp = *pmatrixcomponents;

   int nrows = SCIPmatrixGetNRows(matrix);
   int ncols = SCIPmatrixGetNColumns(matrix);

   comp->nmatrixrows = nrows;
   comp->nmatrixcols = ncols;

   SCIP_CALL( SCIPallocBufferArray(scip, &comp->coltype, ncols) );
   for( int i = 0; i < ncols; ++i )
   {
      SCIP_VAR* var = SCIPmatrixGetVar(matrix,i);
      comp->coltype[i] = SCIPvarGetType(var);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &comp->rowcomponent, nrows) );
   for( int i = 0; i < nrows; ++i )
   {
      comp->rowcomponent[i] = -1;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->colcomponent, ncols) );
   for( int i = 0; i < ncols; ++i )
   {
      comp->colcomponent[i] = -1;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentrows, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentcols, ncols) );
   /* There will be at most ncols components */
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentrowend, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentcolend, ncols) );

   comp->ncomponents = 0;

   return SCIP_OKAY;
}

/** Frees the matrix components data structure */
static
void freeMatrixComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   MATRIX_COMPONENTS**   pmatrixcomponents   /**< Pointer to the allocated matrix components data structure */
   )
{
   MATRIX_COMPONENTS* comp = *pmatrixcomponents;

   /* Make sure to free in reverse */
   SCIPfreeBufferArray(scip, &comp->componentcolend);
   SCIPfreeBufferArray(scip, &comp->componentrowend);
   SCIPfreeBufferArray(scip, &comp->componentcols);
   SCIPfreeBufferArray(scip, &comp->componentrows);
   SCIPfreeBufferArray(scip, &comp->colcomponent);
   SCIPfreeBufferArray(scip, &comp->rowcomponent);
   SCIPfreeBufferArray(scip, &comp->coltype);

   SCIPfreeBuffer(scip, pmatrixcomponents);
}

/** Finds the representative of an element in the disjoint set datastructure.
 *
 * Afterwards compresses the path to speed up subsequent queries.
 */
static
int disjointSetFind(
   int*                  disjointset,        /**< The array storing the disjoint set representatives */
   int                   ind                 /**< The index to find */
   )
{
   assert(disjointset != NULL);

   int current = ind;
   int next;
   /* traverse down tree */
   while( (next = disjointset[current]) >= 0 )
   {
      current = next;
   }
   int root = current;

   /* compress indices along path */
   current = ind;
   while( (next = disjointset[current]) >= 0 )
   {
      disjointset[current] = root;
      current = next;
   }

   return root;
}

/** Merges two sets/elements into one set. Returns the index of the merged element.
 *
 * The provided elements to be merged must be representative (i.e. returned by disjointSetFind()).
 */
static
int disjointSetMerge(
   int*                  disjointset,        /**< The array storing the disjoint set representatives */
   int                   first,              /**< The first index to merge */
   int                   second              /**< The second index to merge */
   )
{
   assert(disjointset);
   assert(disjointset[first] <= -1);
   assert(disjointset[second] <= -1);
   assert(first != second); /* We cannot merge a node into itself */

   /* The rank is stored as a negative number: we decrement it making the negative number larger.
    * The rank is an upper bound on the height of the tree. We want the new root to be the one with 'largest' rank,
    * so smallest number. This way, we ensure that the tree remains shallow. If they are equal, we decrement.
    */
   int firstRank = disjointset[first];
   int secondRank = disjointset[second];
   if( firstRank > secondRank )
   {
      SCIPswapInts(&first, &second);
   }
   /* first becomes representative */
   disjointset[second] = first;
   if( firstRank == secondRank )
   {
      --disjointset[first];
   }

   return first;
}

/** Computes the continuous connected components, i.e. the connected components of the submatrix given by all
 * continuous columns.
 */
static
SCIP_RETCODE computeContinuousComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< The constraint matrix to compute the components for */
   MATRIX_COMPONENTS*    comp                /**< The connected components data structure to store the components in */
   )
{
   /* We let rows and columns share an index by mapping row index i to artificial column index i + nmatrixcols */
   int* disjointset = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &disjointset, comp->nmatrixcols + comp->nmatrixrows) );
   /* First n entries belong to columns, last entries to rows */
   for( int i = 0; i < comp->nmatrixcols + comp->nmatrixrows; ++i )
   {
      disjointset[i] = -1;
   }

   for( int col = 0; col < comp->nmatrixcols; ++col )
   {
      if( comp->coltype[col] != SCIP_VARTYPE_CONTINUOUS )
      {
         continue;
      }
      int colnnonzs = SCIPmatrixGetColNNonzs(matrix, col);
      int* colrows = SCIPmatrixGetColIdxPtr(matrix, col);

      int colrep = disjointSetFind(disjointset, col);
      for( int i = 0; i < colnnonzs; ++i )
      {
         int colrow = colrows[i];
         int ind = colrow + comp->nmatrixcols;
         int rowrep = disjointSetFind(disjointset, ind);
         if( colrep != rowrep )
         {
            colrep = disjointSetMerge(disjointset, colrep, rowrep);
         }
      }
   }

   /* Now, fill in the relevant data. */
   int* representativecomponent;
   SCIP_CALL( SCIPallocBufferArray(scip, &representativecomponent, comp->nmatrixcols + comp->nmatrixrows) );

   for( int i = 0; i < comp->nmatrixcols + comp->nmatrixrows; ++i )
   {
      representativecomponent[i] = -1;
   }
   comp->ncomponents = 0;
   for( int col = 0; col < comp->nmatrixcols; ++col )
   {
      if( comp->coltype[col] != SCIP_VARTYPE_CONTINUOUS )
      {
         continue;
      }
      int colroot = disjointSetFind(disjointset, col);
      int component = representativecomponent[colroot];
      if( component < 0 )
      {
         assert(component == -1);
         /* add new component */
         component = comp->ncomponents;
         representativecomponent[colroot] = component;
         comp->componentcolend[component] = 0;
         comp->componentrowend[component] = 0;
         ++comp->ncomponents;
      }
      comp->colcomponent[col] = component;
      ++comp->componentcolend[component];
   }
   for( int row = 0; row < comp->nmatrixrows; ++row )
   {
      int rowroot = disjointSetFind(disjointset, row + comp->nmatrixcols);
      int component = representativecomponent[rowroot];
      if( component < 0 )
      {
         assert(component == -1);
         /* Any rows that have roots that we have not seen yet are rows that have no continuous columns
          * We can safely skip these for finding the continuous connected components
          */
         continue;
      }
      comp->rowcomponent[row] = component;
      ++comp->componentrowend[component];
   }
   if( comp->ncomponents >= 1 )
   {
      for( int i = 1; i < comp->ncomponents; ++i )
      {
         comp->componentrowend[i] += comp->componentrowend[i-1];
         comp->componentcolend[i] += comp->componentcolend[i-1];
      }
      int * componentnextrowindex;
      int * componentnextcolindex;
      SCIP_CALL( SCIPallocBufferArray(scip, &componentnextrowindex, comp->ncomponents) );
      SCIP_CALL( SCIPallocBufferArray(scip, &componentnextcolindex, comp->ncomponents) );
      componentnextrowindex[0] = 0;
      componentnextcolindex[0] = 0;
      for( int i = 1; i < comp->ncomponents; ++i )
      {
         componentnextcolindex[i] = comp->componentcolend[i-1];
         componentnextrowindex[i] = comp->componentrowend[i-1];
      }

      for( int col = 0; col < comp->nmatrixcols; ++col )
      {
         int component = comp->colcomponent[col];
         if( component < 0 )
         {
            assert(component == -1);
            continue;
         }
         int ind = componentnextcolindex[component];
         comp->componentcols[ind] = col;
         ++componentnextcolindex[component];
      }
      for( int row = 0; row < comp->nmatrixrows; ++row )
      {
         int component = comp->rowcomponent[row];
         if( component < 0 )
         {
            assert(component == -1);
            continue;
         }
         int ind = componentnextrowindex[component];
         comp->componentrows[ind] = row;
         ++componentnextrowindex[component];
      }

#ifndef NDEBUG
      for( int i = 0; i < comp->ncomponents; ++i )
      {
         assert(componentnextrowindex[i] == comp->componentrowend[i]);
         assert(componentnextcolindex[i] == comp->componentcolend[i]);
      }
#endif

      SCIPfreeBufferArray(scip, &componentnextcolindex);
      SCIPfreeBufferArray(scip, &componentnextrowindex);
   }

   SCIPfreeBufferArray(scip, &representativecomponent);
   SCIPfreeBufferArray(scip, &disjointset);

   return SCIP_OKAY;
}

/** Creates the matrix statistics data structure */
static
SCIP_RETCODE computeMatrixStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< The constraint matrix to compute the statistics for */
   MATRIX_COMPONENTS*    comp,               /**< Datastructure that contains the components of the matrix */
   MATRIX_STATISTICS**   pstats,             /**< Pointer to allocate the statistics data structure at */
   SCIP_Real             numericslimit       /**< The limit beyond which we consider integrality of coefficients
                                              *   to be unreliable */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, pstats) );
   MATRIX_STATISTICS* stats = *pstats;

   int nrows = SCIPmatrixGetNRows(matrix);
   int ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowintegral, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowequality, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowbadnumerics, nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rownnonz, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowncontinuous, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowncontinuouspmone, nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &stats->colintegralbounds, ncols) );

   for( int i = 0; i < nrows; ++i )
   {
      SCIP_Real lhs = SCIPmatrixGetRowLhs(matrix, i);
      SCIP_Real rhs = SCIPmatrixGetRowRhs(matrix, i);
      int* cols = SCIPmatrixGetRowIdxPtr(matrix, i);
      SCIP_Real* vals = SCIPmatrixGetRowValPtr(matrix, i);
      int nnonz = SCIPmatrixGetRowNNonzs(matrix, i);
      stats->rownnonz[i] = nnonz;
      stats->rowequality[i] = !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && SCIPisEQ(scip, lhs, rhs);

      SCIP_Bool integral = ( SCIPisInfinity(scip, -lhs) || SCIPisIntegral(scip, lhs) )
                        && ( SCIPisInfinity(scip, rhs) || SCIPisIntegral(scip, rhs) );
      SCIP_Bool badnumerics = FALSE;

      int ncontinuous = 0;
      int ncontinuouspmone = 0;
      for( int j = 0; j < nnonz; ++j )
      {
         SCIP_Bool continuous = comp->coltype[cols[j]] == SCIP_VARTYPE_CONTINUOUS;
         SCIP_Real value = vals[j];
         if( continuous )
         {
            ++ncontinuous;
            if( SCIPisEQ(scip, ABS(value), 1.0) )
            {
               ++ncontinuouspmone;
            }
         }
         else
         {
            integral = integral && SCIPisIntegral(scip, value);
         }
         if( ABS(value) > numericslimit )
         {
            badnumerics = TRUE;
         }
      }

      stats->rowncontinuous[i] = ncontinuous;
      stats->rowncontinuouspmone[i] = ncontinuouspmone;
      stats->rowintegral[i] = integral;
      stats->rowbadnumerics[i] = badnumerics;
   }

   for( int i = 0; i < ncols; ++i )
   {
      SCIP_Real lb = SCIPmatrixGetColLb(matrix, i);
      SCIP_Real ub = SCIPmatrixGetColUb(matrix, i);
      stats->colintegralbounds[i] = ( SCIPisInfinity(scip, -lb) || SCIPisIntegral(scip, lb) )
                                 && ( SCIPisInfinity(scip, ub) || SCIPisIntegral(scip, ub) );

      /* Check that integer variables have integer bounds, as expected. */
      assert(comp->coltype[i] == SCIP_VARTYPE_CONTINUOUS || stats->colintegralbounds[i]);
   }

   return SCIP_OKAY;
}

/** Frees the matrix statistics data structure */
static
void freeMatrixStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   MATRIX_STATISTICS**   pstats              /**< Pointer to the statistics data structure to be freed */
   )
{
   MATRIX_STATISTICS* stats= *pstats;

   /* Make sure, for performance, that these frees occur in reverse */
   SCIPfreeBufferArray(scip, &stats->colintegralbounds);
   SCIPfreeBufferArray(scip, &stats->rowncontinuouspmone);
   SCIPfreeBufferArray(scip, &stats->rowncontinuous);
   SCIPfreeBufferArray(scip, &stats->rownnonz);
   SCIPfreeBufferArray(scip, &stats->rowbadnumerics);
   SCIPfreeBufferArray(scip, &stats->rowequality);
   SCIPfreeBufferArray(scip, &stats->rowintegral);

   SCIPfreeBuffer(scip, pstats);
}

/** Given the continuous components and statistics on the matrix, detect components that have implied integer variables
 * by checking if the component is a (transposed) network matrix and if all the bounds/sides/coefficients are integral.
 *
 * For every component, we detect if the associated matrix is either a network matrix or a transposed network matrix
 * (or both, in which case it represents a planar graph).
 * We choose to check if it is a (transposed) network matrix either in a row-wise or in a column-wise fashion,
 * depending on the size of the component.
 * Finally, every variable that is in a network matrix or transposed network matrix is changed to an implied integer.
 */
static
SCIP_RETCODE findImpliedIntegers(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< The data belonging to the presolver */
   SCIP_MATRIX*          matrix,             /**< The constraint matrix to compute implied integers for */
   MATRIX_COMPONENTS*    comp,               /**< The continuous connected components of the matrix */
   MATRIX_STATISTICS*    stats,              /**< Statistics of the matrix */
   int*                  nchgvartypes        /**< Pointer to count the number of changed variable types */
   )
{
   /* TODO: some checks to prevent expensive memory initialization if not necessary (e.g. there must be some candidates) */
   SCIP_NETMATDEC* dec = NULL;
   SCIP_CALL( SCIPnetmatdecCreate(SCIPblkmem(scip), &dec, comp->nmatrixrows, comp->nmatrixcols) );

   SCIP_NETMATDEC* transdec = NULL;
   SCIP_CALL( SCIPnetmatdecCreate(SCIPblkmem(scip), &transdec, comp->nmatrixcols, comp->nmatrixrows) );

   int ngoodcomponents = 0;
   int nbadnumerics = 0;
   int nbadintegrality = 0;
   int nnonnetwork = 0;

   /* Because the rows may also contain non-continuous columns, we need to remove these from the array that we
    * pass to the network matrix decomposition method. We use these working arrays for this purpose.
    */
   SCIP_Real* tempValArray;
   int* tempIdxArray;
   SCIP_CALL( SCIPallocBufferArray(scip, &tempValArray, comp->nmatrixcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tempIdxArray, comp->nmatrixcols) );

   for( int component = 0; component < comp->ncomponents; ++component )
   {
      int startrow = (component == 0) ? 0 : comp->componentrowend[component-1];
      int nrows = comp->componentrowend[component] - startrow;
      SCIP_Bool componentokay = TRUE;
      for( int i = startrow; i < startrow + nrows; ++i )
      {
         int row = comp->componentrows[i];
         if( stats->rowncontinuous[row] != stats->rowncontinuouspmone[row] || !stats->rowintegral[row] )
         {
            componentokay = FALSE;
            ++nbadintegrality;
            break;
         }
         if( stats->rowbadnumerics[row] )
         {
            componentokay = FALSE;
            ++nbadnumerics;
            break;
         }
      }
      int startcol = ( component == 0 ) ? 0 : comp->componentcolend[component-1];
      int ncols = comp->componentcolend[component] - startcol;

      for( int i = startcol; i < startcol + ncols && componentokay; ++i )
      {
         int col = comp->componentcols[i];
         if( !stats->colintegralbounds[col] )
         {
            componentokay = FALSE;
            ++nbadintegrality;
            break;
         }
      }

      if(!componentokay)
      {
         continue;
      }

      /* Check if the component is a network matrix */
      SCIP_Bool componentnetwork = TRUE;

      /* We use the row-wise algorithm only if the number of columns is much larger than the number of rows.
       * Generally, the column-wise algorithm will be faster, but in these extreme cases, the row algorithm is faster.
       * Only very few instances should use the row-wise algorithm.
       */
      if( nrows * presoldata->columnrowratio < ncols )
      {
         for( int i = startrow; i < startrow + nrows && componentnetwork; ++i )
         {
            int row = comp->componentrows[i];
            int nrownnoz = SCIPmatrixGetRowNNonzs(matrix, row);
            int* rowcols = SCIPmatrixGetRowIdxPtr(matrix, row);
            SCIP_Real* rowvals = SCIPmatrixGetRowValPtr(matrix, row);
            int ncontnonz = 0;
            for( int j = 0; j < nrownnoz; ++j )
            {
               int col = rowcols[j];
               if( comp->coltype[col] == SCIP_VARTYPE_CONTINUOUS )
               {
                  tempIdxArray[ncontnonz] = col;
                  tempValArray[ncontnonz] = rowvals[j];
                  ++ncontnonz;
                  assert(SCIPisEQ(scip, ABS(rowvals[j]), 1.0));
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddRow(dec, row, tempIdxArray, tempValArray, ncontnonz, &componentnetwork) );
         }
      }
      else
      {
         for( int i = startcol; i < startcol + ncols && componentnetwork; ++i )
         {
            int col = comp->componentcols[i];
            int ncolnnonz = SCIPmatrixGetColNNonzs(matrix, col);
            int* colrows = SCIPmatrixGetColIdxPtr(matrix, col);
            SCIP_Real* colvals = SCIPmatrixGetColValPtr(matrix, col);
            SCIP_CALL( SCIPnetmatdecTryAddCol(dec, col, colrows, colvals, ncolnnonz, &componentnetwork) );
         }
      }

      if( !componentnetwork )
      {
         SCIPnetmatdecRemoveComponent(dec, &comp->componentrows[startrow], nrows, &comp->componentcols[startcol], ncols);
      }

      SCIP_Bool componenttransnetwork = TRUE;

      /* For the transposed matrix, the situation is exactly reversed because the row/column algorithms are swapped
       * We only run transposed network detection if network detection failed
       */
      if( !componentnetwork )
      {
         if( nrows <= ncols * presoldata->columnrowratio )
         {
            for( int i = startrow; i < startrow + nrows && componenttransnetwork; ++i )
            {
               int row = comp->componentrows[i];
               int nrownnoz = SCIPmatrixGetRowNNonzs(matrix, row);
               int* rowcols = SCIPmatrixGetRowIdxPtr(matrix, row);
               SCIP_Real* rowvals = SCIPmatrixGetRowValPtr(matrix, row);
               int ncontnonz = 0;
               for( int j = 0; j < nrownnoz; ++j )
               {
                  int col = rowcols[j];
                  if( comp->coltype[col] == SCIP_VARTYPE_CONTINUOUS )
                  {
                     tempIdxArray[ncontnonz] = col;
                     tempValArray[ncontnonz] = rowvals[j];
                     ++ncontnonz;
                     assert(SCIPisEQ(scip, ABS(rowvals[j]), 1.0));
                  }
               }

               SCIP_CALL( SCIPnetmatdecTryAddCol(transdec, row, tempIdxArray, tempValArray, ncontnonz,
                                                &componenttransnetwork) );
            }
         }
         else
         {
            for( int i = startcol; i < startcol + ncols && componenttransnetwork; ++i )
            {
               int col = comp->componentcols[i];
               int ncolnnonz = SCIPmatrixGetColNNonzs(matrix, col);
               int* colrows = SCIPmatrixGetColIdxPtr(matrix, col);
               SCIP_Real* colvals = SCIPmatrixGetColValPtr(matrix, col);
               SCIP_CALL( SCIPnetmatdecTryAddRow(transdec, col, colrows, colvals, ncolnnonz, &componenttransnetwork) );
            }
         }
      }

      if( !componenttransnetwork )
      {
         SCIPnetmatdecRemoveComponent(transdec, &comp->componentcols[startcol], ncols,
                                      &comp->componentrows[startrow], nrows);
      }

      if( !componentnetwork && !componenttransnetwork)
      {
         ++nnonnetwork;
         continue;
      }
      ++ngoodcomponents;

      for( int i = startcol; i < startcol + ncols; ++i )
      {
         int col = comp->componentcols[i];
         SCIP_VAR* var = SCIPmatrixGetVar(matrix, col);
         SCIP_Bool infeasible = FALSE;
         SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_VARIMPLTYPE_WEAK, &infeasible) );
         (*nchgvartypes)++;
         assert(!infeasible);
      }
   }
   if( *nchgvartypes == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "No implied integers detected \n");
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
                      "%d implied integers in %d / %d components (disqualified: %d by integrality, %d by numerics, %d not network) \n",
                      *nchgvartypes, ngoodcomponents, comp->ncomponents, nbadintegrality, nbadnumerics, nnonnetwork);
   }

   SCIPfreeBufferArray(scip, &tempIdxArray);
   SCIPfreeBufferArray(scip, &tempValArray);

   SCIPnetmatdecFree(&transdec);
   SCIPnetmatdecFree(&dec);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for presolver plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyImplint)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolImplint(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeImplint)
{
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecImplint)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   /* TODO: re-check these conditions again */
   /* Disable implicit integer detection if we are probing or in NLP context */
   if( ( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING ) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
   {
      return SCIP_OKAY;
   }
   /* Since implied integer detection relies on rows being static, we disable it for branch-and-price applications*/
   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
   {
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* Exit early if there are no candidates variables to upgrade */
   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);
   if( SCIPgetNContVars(scip) == 0 )
   {
      return SCIP_OKAY;
   }

   SCIP_Real starttime = SCIPgetSolvingTime(scip);
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
   if( !( initialized && complete ) )
   {
      if( initialized )
      {
         SCIPmatrixFree(scip, &matrix);
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
                      "   (%.1fs) implied integer detection stopped because problem is not an MILP\n",
                      SCIPgetSolvingTime(scip));
      return SCIP_OKAY;
   }

   int beforechanged = *nchgvartypes;
   MATRIX_COMPONENTS* comp = NULL;
   MATRIX_STATISTICS* stats = NULL;
   SCIP_CALL( createMatrixComponents(scip, matrix, &comp) );
   SCIP_CALL( computeMatrixStatistics(scip, matrix, comp, &stats, presoldata->numericslimit) );
   SCIP_CALL( computeContinuousComponents(scip, matrix, comp) );
   SCIP_CALL( findImpliedIntegers(scip, presoldata, matrix, comp, stats, nchgvartypes) );
   int afterchanged = *nchgvartypes;

   SCIP_Real endtime = SCIPgetSolvingTime(scip);
   if( afterchanged == beforechanged )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "   (%.1fs) no implied integers detected (time: %.2fs)\n", endtime, endtime - starttime);
      *result = SCIP_DIDNOTFIND;
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "   (%.1fs) %d implied integers detected (time: %.2fs)\n",endtime,*nchgvartypes,endtime-starttime);

      *result = SCIP_SUCCESS;
   }
   freeMatrixStatistics(scip, &stats);
   freeMatrixComponents(scip, &comp);
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
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include implint presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
                                     PRESOL_TIMING, presolExecImplint, presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplint) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeImplint) );

   SCIP_CALL( SCIPaddRealParam(scip,
                               "presolving/implint/columnrowratio",
                               "use the network row addition algorithm when the column to row ratio becomes larger than this threshold",
                               &presoldata->columnrowratio, TRUE, DEFAULT_COLUMNROWRATIO, 0.0, 1e98, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
                               "presolving/implint/numericslimit",
                               "a row that contains variables with coefficients that are greater in absolute value than this limit is not considered for implied integrality detection",
                               &presoldata->numericslimit, TRUE, DEFAULT_NUMERICSLIMIT, 1.0, 1e98, NULL, NULL) );
   return SCIP_OKAY;
}
