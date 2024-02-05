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

/**@file   symmetry_orbitopal.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling orbitopal symmetries
 * @author Jasper van Doornmalen
 *
 * This implements orbitopal reducion, which generalizes full orbitope propagation to work for non-binary variable
 * domains, and is dynamified. See cons_orbitope.c for the variant for binary variables, both the static and partially
 * dynamic variant.
 * Just as in orbital reduction (cf. symmetry_orbital.c), the variable order is chosen as the variables
 * branched on from the root node to the focus node.
 *
 * See Section 4.2, Example 12 and Section 5.1 in [vD,H]:@n
 * J. van Doornmalen, C. Hojny, "A Unified Framework for Symmetry Handling", preprint, 2023,
 * https://doi.org/10.48550/arXiv.2211.01295.
 *
 * Orbitopal reduction can be used to handle symmetries of the following type.
 * If the variables can be arranged in a matrix and the symmetries correspond to all column permutations of this matrix,
 * then these symmetries are called orbitopal.
 * Symmetry is completely handled by enforcing that the columns are lexicographically decreasing.
 * If a reduction on a variable is applied, and if this variable is high up in the variable matrix, then this has a
 * relatively large impact on the lexicographic ordering. Moreover, the ordering of the columns also matter.
 * Dynamification allows for choosing a custom ordering of a subset of rows and a permutation of the columns.
 * For every node, we maintain the ordered subset of rows and the column order.
 * The root node assumes no rows and an arbitrary column order (we choose the identity).
 * If towards a new node it is branched on a variable, that appears in a row which is not included in the subset of
 * rows for the current node, then the row set of the new children is the ordered row set of the current node, appended
 * by this new row.
 * For the column order, if at the current node columns are symmetrically equivalent, then these can be permuted for
 * the sake of symmetry handling. In the implementation, we only swap the column containing the branching variable
 * with a symmetrically equivalent column elsewhere. We use one of the following heuristics:
 *
 * - None: Keep the column-order as-is.
 * - First: Permute such that the column containing the branching variable is in the symmetrically equivalent column
 *   with minimal index.
 * - Last: The same, but to the symmetrically equivalent column with maximal index.
 * - Centre: The same, but to the symmetrically equivalent column closest to to the middlemost column among all columns.
 * - Median: The same, but to the median of all symmetrically equivalent columns. (This is the default)
 *
 * Since the dynamic row and column ordering rule for a branch-and-bound tree node depends on the decisions made up to
 * that branch-and-bound tree node, we compute and store the row and column order for the branch-and-bound tree children
 * at the moment of branching. This is done by the eventhandler in this file.
 * Instead of storing those, we could have chosen to reconstruct this row and column ordering to save memory.
 * However, we cannot reliably reconstruct this order from the branch-and-bound tree itself,
 * because the row and column ordering depends on symmetrical equivalence of columns in the orbitope matrix,
 * and because SCIP can change the tree structure during solving that may re-write historic variable bound changes
 * (for instance when global variable bound changes are found, or when the root node is moved down the tree to become
 * the new effective root node).
 * We are not concerned about storing the row and column ordering, since we only store the mutations with its parent.
 * These are usually at most one column swap and usually at most one additional row.
 *
 * @todo Exploiting packing-partitioning structures
 * @todo for packing-partitioning rows, use the FIRST column swap heuristic.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/symmetry_orbitopal.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/struct_var.h"
#include "scip/type_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/struct_scip.h"
#include "scip/struct_mem.h"
#include "scip/struct_tree.h"
#include "scip/symmetry.h"
#include "scip/debug.h"
#include <string.h>
#include <symmetry/type_symmetry.h>


/* symmetry handler properties */
#define SYMHDLR_NAME           "orbitopalreduction"

/* orbitopal symmetry handler properties */
#define EVENTHDLR_NAME         "symmetry_orbitopal_eventhdlr"
#define EVENTHDLR_DESC         "event handler for maintaining the branch-and-bound tree"
#define DEFAULT_COLUMNORDERING SCIP_COLUMNORDERING_MEDIAN /**< the column ordering variant */

/*
 * Data structures
 */

/** orbitopal symmetry handling data for a single orbitope */
struct OrbitopeData
{
   SCIP_VAR**            vars;               /**< orbitope variable array representing orbitope matrix row-wise */
   int                   nrows;              /**< number of rows */
   int                   ncols;              /**< number of columns */
   int                   nbranchrows;        /**< number of rows containing variables potentially used for branching */
   SCIP_HASHMAP*         rowindexmap;        /**< map of variables to row index in orbitope matrix */
   SCIP_HASHMAP*         colindexmap;        /**< map of variables to column index in orbitope matrix */
#ifndef NDEBUG
   SCIP_Longint          lastnodenumber;     /**< the last node number where the row and column order is computed */
   int                   dbghash;            /**< a hash for the column order in the last iteration */
#endif
   SCIP_HASHTABLE*       nodeinfos;          /**< symmetry handling information per branch-and-bound tree node */
   SCIP_COLUMNORDERING   columnordering;     /**< policy for the column ordering */
   SCIP_ROWORDERING      rowordering;        /**< policy for the row ordering */
};
typedef struct OrbitopeData ORBITOPEDATA; /**< orbitopal symmetry handling data for a single orbitope */

/** wrapper for all orbitopes in orbitopal symmetry handling data */
struct SCIP_OrbitopalReductionData
{
   SCIP_COLUMNORDERING   defaultcolumnordering; /**< default policy for the column ordering */
   SCIP_EVENTHDLR*       eventhdlr;          /**< pointer to the event handler for managing the branching tree */
   ORBITOPEDATA**        orbitopes;          /**< array of pointers to orbitope data structs */
   int                   norbitopes;         /**< number of orbitope data structs in array */
   int                   maxnorbitopes;      /**< allocated orbitopes array size */
   SCIP_CONSHDLR*        conshdlr_nonlinear; /**< nonlinear constraint handler,
                                              * is used to determine if a variable is a branching variable */
   SCIP_Bool             conshdlr_nonlinear_checked; /**< nonlinear constraint handler is already added? */
   int                   nred;               /**< total number of reductions */
   int                   ncutoff;            /**< total number of cutoffs */
};

/*
 * Local methods
 */

/** gets whether a variable type is a branchrow-type */
static
SCIP_Bool vartypeIsBranchRowType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< pointer to the dynamic orbitopal reduction data */
   SCIP_VARTYPE          vartype             /**< var type */
)
{
   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( orbireddata->conshdlr_nonlinear_checked );

   switch (vartype)
   {
   case SCIP_VARTYPE_BINARY:
   case SCIP_VARTYPE_INTEGER:
      return TRUE;
   case SCIP_VARTYPE_CONTINUOUS:
   case SCIP_VARTYPE_IMPLINT:
      /* potential branching variables if nonlinear constraints exist */
      assert( orbireddata->conshdlr_nonlinear_checked );
      return orbireddata->conshdlr_nonlinear == NULL ? FALSE :
         SCIPconshdlrGetNActiveConss(orbireddata->conshdlr_nonlinear) > 0;
   default:
      SCIPerrorMessage("unknown vartype\n");
      SCIPABORT();
      /* resolve compiler warning: no asserts in optimized mode */
      return FALSE;
   }
}


/** container for column index permutations */
struct ColSwap
{
   int                   from;               /**< from which column ID */
   int                   to;                 /**< to which column ID */
};
typedef struct ColSwap COLSWAP;

/** information stored for branch-and-bound nodes */
struct BnbNodeInfo
{
   SCIP_Longint          nodenumber;         /**< node number of the branch-and-bound tree node */
   COLSWAP*              colswaps;           /**< list containing column swaps by node branching decisions */
   int                   ncolswaps;          /**< number of elements in colswaps. ncolswaps == 0 <=> colswaps == NULL */
   int*                  rows;               /**< list of new variable rows by node branching decisions */
   int                   nrows;              /**< number of new variable rows added. nrows == 0 <=> rows == NULL */
};
typedef struct BnbNodeInfo BNBNODEINFO;

/** hash key for virtual branch and bound nodeinfo struct */
static
SCIP_DECL_HASHGETKEY(hashGetKeyBnbnodeinfo)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both node numbers are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqBnbnodeinfo)
{  /*lint --e{715}*/
   BNBNODEINFO* nodeinfo1 = (BNBNODEINFO*) key1;
   BNBNODEINFO* nodeinfo2 = (BNBNODEINFO*) key2;
   return nodeinfo1->nodenumber == nodeinfo2->nodenumber;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValBnbnodeinfo)
{  /*lint --e{715}*/
   BNBNODEINFO* nodeinfo = (BNBNODEINFO*) key;
   return (unsigned int) nodeinfo->nodenumber;
}


/** tests if two columns are symmetrically equivalent
 *
 * We test if the columns with index col1 and col2 have elementwise the same bounds.
 * If all symmetry-compatible reductions are applied, then it suffices to check only as many rows as are selected
 * for orbitopal reduction. However, to be resilient to reductions that are not symmetry-compatible,
 * we test all variables in the columns.
 */
static
SCIP_Bool testColumnsAreSymmetricallyEquivalent(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope information */
   int                   col1,               /**< first column to compare */
   int                   col2                /**< second column to compare */
   )
{
   int i;
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( col1 >= 0 );
   assert( col1 < orbidata->ncols );
   assert( col2 >= 0 );
   assert( col2 < orbidata->ncols );

   /* @todo test only for the selected rows (see function description) */
   for (i = 0; i < orbidata->nrows; ++i)
   {
      var1 = orbidata->vars[i * orbidata->ncols + col1];
      var2 = orbidata->vars[i * orbidata->ncols + col2];

      /* if variable bounds differ: columns c and origcolid are not the same */
      if (
         (! SCIPsymEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetLbLocal(var2)))
         ||
         (! SCIPsymEQ(scip, SCIPvarGetUbLocal(var1), SCIPvarGetUbLocal(var2)))
      )
         return FALSE;
   }

   /* loop terminated, so columns are equal */
   return TRUE;
}

/** updates the column order with a bound change
 *
 *  When it is branched on a variable in a column, update the column order for the children of the focusnode.
 *  Symmetrically equivalent columns, that is the columns where the variables have elementwise the same domain,
 *  at the focusnode at the moment of branching can be permuted.
 *  In this function, we select such a permutation, based on the column containing the branching variable(s).
 *  In all cases, we swap the column containing the branching variable with a symmetrically equivalent column,
 *  and the @param columnordering specifies if we prefer it to be the leftmost, rightmost, centermost symmetrically
 *  equivalent column, or the median column among the symmetrically equivalent columns.
 *
 *  The column ordering is determined and stored at the moment of branching.
 */
static
SCIP_RETCODE updateColumnOrderWhenBranchingOnColumn(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   int*                  colorder,           /**< array to populate with column order, of size colorder */
   int*                  colorderinv,        /**< inverse array of the column order, of size colorder */
   SCIP_VAR*             var,                /**< variable that we branch on */
   COLSWAP*              thiscolswap         /**< the colswap to populate */
   )
{
   int origcolid;
   int swaporigcolid;
   int c;
   int ncols;
   int* origequalcolids;
   int norigequalcolids;
   int middlecolumn = 0;
   int positionorigcolidincolorder;
   int positionswaporigcolidincolorder;

#ifndef NDEBUG
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int i;
   int nrows;
#endif

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( colorder != NULL );
   assert( colorderinv != NULL );
   assert( var != NULL );
   assert( thiscolswap != NULL );
   assert( orbidata->ncols > 0 );
   assert( orbidata->nrows > 0 );

   ncols = orbidata->ncols;
   assert( ncols > 0 );
#ifndef NDEBUG
   nrows = orbidata->nrows > 0;
   assert( nrows > 0 );
#endif

   /* do not apply a column swap if no column permutations are applied */
   if ( orbidata->columnordering == SCIP_COLUMNORDERING_NONE )
   {
      thiscolswap->from = 0;
      thiscolswap->to = 0;
      return SCIP_OKAY;
   }

   /* only variables from the orbitope matrix are of interest */
   origcolid = SCIPhashmapGetImageInt(orbidata->colindexmap, (void*) var);
   if ( origcolid == INT_MAX )
   {
      thiscolswap->from = 0;
      thiscolswap->to = 0;
      return SCIP_OKAY;
   }
   assert( origcolid >= 0 );
   assert( origcolid < ncols );

   /* policy: swap with identical column that is closest to the center in relabeled order */
   /* @todo other policies: If the variable is in a ppc-row, then select the minimal/second to minimal to branch on */
   swaporigcolid = origcolid;

   switch (orbidata->columnordering)
   {
   case SCIP_COLUMNORDERING_CENTRE:
      /* CENTRE uses the same code as FIRST and LAST, but requires that middlecolumn = ncols / 2 is set */
      middlecolumn = ncols / 2;
      /*lint -fallthrough*/
   case SCIP_COLUMNORDERING_FIRST:
   case SCIP_COLUMNORDERING_LAST:
      /* for each column, test column ordering condition, then if the column is identical to the var-column */
      for (c = 0; c < ncols; ++c)
      {
         /* origcolid is not interesting */
         if ( c == origcolid )
            continue;

         /* test if c is a better choice than swaporigcolid,
          * otherwise continue to next iteration through CONDITIONFAIL
          */
         switch (orbidata->columnordering)
         {
            case SCIP_COLUMNORDERING_FIRST:
               /* only swap with c if c is earlier in column order than swaporigcolid */
               if ( colorderinv[c] >= colorderinv[swaporigcolid] )
                  goto CONDITIONFAIL;
               break;
            case SCIP_COLUMNORDERING_LAST:
               /* only swap with c if c is later in column order than swaporigcolid */
               if ( colorderinv[c] <= colorderinv[swaporigcolid] )
                  goto CONDITIONFAIL;
               break;
            case SCIP_COLUMNORDERING_CENTRE:
               /* if the column is not more central than swaporigcolid, ignore */
               if ( ABS(colorderinv[c] - middlecolumn) >=
                     ABS(colorderinv[swaporigcolid] - middlecolumn) )
                  goto CONDITIONFAIL;
               break;
            default:
               return SCIP_ERROR;
         }

         /* test: are c and origcolid the same columns w.r.t. the variable domain restrictions? */
         if ( !testColumnsAreSymmetricallyEquivalent(scip, orbidata, c, origcolid) )
            continue;

         /* the variable domain reductions in c and origcolid are the same */
         swaporigcolid = c;

      CONDITIONFAIL:
         ;  /* no-op for going to the next iteration */
      }

      /* end switch */
      break;

   case SCIP_COLUMNORDERING_MEDIAN:
      /* collect columns identical to the var-column, then select column satisfying ordering condition */
      norigequalcolids = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &origequalcolids, ncols) );

      /* collect equal columns */
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("Detect columns identical to original column %d: ", origcolid);
#endif
      for (c = 0; c < ncols; ++c)
      {
         /* column origcolid is always equal to itself */
         if ( c == origcolid )
         {
            origequalcolids[norigequalcolids++] = c;
#ifdef SCIP_MORE_DEBUG
            SCIPdebugPrintf("%d ", c);
#endif
            assert( norigequalcolids <= ncols );
            continue;
         }

         /* test: are c and origcolid the same columns w.r.t. the variable domain restrictions? */
         if ( !testColumnsAreSymmetricallyEquivalent(scip, orbidata, c, origcolid) )
            continue;

         /* the variable domain reductions in c and origcolid are the same */
         origequalcolids[norigequalcolids++] = c;
#ifdef SCIP_MORE_DEBUG
         SCIPdebugPrintf("%d ", c);
#endif
         assert( norigequalcolids <= ncols );
      }
#ifdef SCIP_MORE_DEBUG
      SCIPdebugPrintf("\n");
#endif

      /* we should have found origcolid, at least */
      assert( norigequalcolids >= 1 );

      /* from origequalcolids, select the column satisfying the column ordering policy */

      /* get median column; since colorder maps origcolids to our ordering,
       * we need to use colorderinv as the argument. */
      /* @todo computing the median is O(n) by repeated median-of-medians (CLRS, Ch9), but let's just sort things */
      SCIPsortInd(origequalcolids, SCIPsortArgsortInt, colorderinv, norigequalcolids);
      /* get the median, that is swaporigcolid */
      swaporigcolid = origequalcolids[norigequalcolids / 2];

      SCIPfreeBufferArray(scip, &origequalcolids);

      /* end switch */
      break;

   case SCIP_COLUMNORDERING_NONE:
      /* already handled earlier in this function */
   default:
      /* unknown column ordering variant */
      return SCIP_ERROR;
   }

   thiscolswap->from = swaporigcolid;
   thiscolswap->to = origcolid;

   /* if we do not replace origcolid */
   if ( swaporigcolid == origcolid )
      return SCIP_OKAY;

#ifndef NDEBUG
   /* swapped columns should be equivalent */
   for (i = 0; i < nrows; ++i)
   {
      var1 = orbidata->vars[i * ncols + swaporigcolid];
      var2 = orbidata->vars[i * ncols + origcolid];
      assert( SCIPsymEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetLbLocal(var2)) );
      assert( SCIPsymEQ(scip, SCIPvarGetUbLocal(var1), SCIPvarGetUbLocal(var2)) );
   }
#endif

   /* now swap the permuted column indices of swaporigcolid and origcolid */

   /* at which column is origcolid? */
   positionorigcolidincolorder = colorderinv[origcolid];
   assert( positionorigcolidincolorder >= 0 );
   assert( positionorigcolidincolorder < ncols );
   assert( colorder[positionorigcolidincolorder] == origcolid );

   /* at which column is swaporigcolid? */
   positionswaporigcolidincolorder = colorderinv[swaporigcolid];
   assert( positionswaporigcolidincolorder >= 0 );
   assert( positionswaporigcolidincolorder < ncols );
   assert( colorder[positionswaporigcolidincolorder] == swaporigcolid );

   SCIPdebugMessage("Orbitope %p: Swapping column %d (at position %d) with column %d (at position %d)\n",
      (void*) orbidata, origcolid, positionorigcolidincolorder, swaporigcolid, positionswaporigcolidincolorder);

   /* swap them, also keep track of the inverses */
   colorder[positionswaporigcolidincolorder] = origcolid;
   colorder[positionorigcolidincolorder] = swaporigcolid;
   colorderinv[origcolid] = positionswaporigcolidincolorder;
   colorderinv[swaporigcolid] = positionorigcolidincolorder;

   return SCIP_OKAY;
}


/** yields entry at index in array, or returns entry if array is NULL */
static
int getArrayEntryOrIndex(
   int*                  arr,                /**< array */
   int                   idx                 /**< index */
)
{
   assert( idx >= 0 );
   if ( arr == NULL )
      return idx;
   return arr[idx];
}

/** frees the row order */
static
void freeRowOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   int**                 roworder            /**< roworder array that is initialized with the roworder in the dynamic
                                               *  case, and NULL in the static case */
)
{
   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( roworder != NULL );

   if ( orbidata->rowordering == SCIP_ROWORDERING_NONE )
   {
      assert( *roworder == NULL );
      return;
   }

   assert( *roworder != NULL );
   assert( orbidata->rowordering == SCIP_ROWORDERING_BRANCHING );
   SCIPfreeBlockMemoryArray(scip, roworder, orbidata->nrows);

   return;
}

/** gets the row order at the node
 *
 *  this is NULL (i.e., the identity map) in the static (none) setting.
 *  this is an array of size orbidata->nrows in the dynamic (branching) setting.
 *
 *  The row order is given in the order of the variables that is branched on.
 *  @todo combine with variant of cons_orbitope.c
 */
static
SCIP_RETCODE getRowOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   SCIP_NODE*            node,               /**< node for which the row order should be detected */
   int**                 roworder,           /**< array to populate with row order */
   int*                  nselrows            /**< pointer to populate with the number of rows part of the row order */
   )
{
   int i;
   int j;
   BNBNODEINFO* ancestornodeinfo;
   BNBNODEINFO tmpnodeinfo;  /* used for lookups in hash table */

   assert( orbidata != NULL );
   assert( orbidata->nrows > 0 );
   assert( orbidata->ncols > 0 );
   assert( node != NULL );
   assert( roworder != NULL );
   assert( nselrows != NULL );

   if ( orbidata->rowordering == SCIP_ROWORDERING_NONE )
   {
      *roworder = NULL;
      *nselrows = orbidata->nrows;
      return SCIP_OKAY;
   }

   assert( orbidata->rowordering == SCIP_ROWORDERING_BRANCHING );

   /* allocate number of rows */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, roworder, orbidata->nrows) );

   *nselrows = 0;

   /* get the present row order up to this node (excluding the node itself) */
   node = SCIPnodeGetParent(node);
   while (node != NULL)
   {
      /* retrieve the nodeinfo of this ancestor node */
      tmpnodeinfo.nodenumber = SCIPnodeGetNumber(node);
      ancestornodeinfo = (BNBNODEINFO*) SCIPhashtableRetrieve(orbidata->nodeinfos, (void*) &tmpnodeinfo);
      if ( ancestornodeinfo != NULL )
      {
         assert( ancestornodeinfo->nrows >= 0 );
         for (i = ancestornodeinfo->nrows - 1; i >= 0; --i)
         {
            (*roworder)[(*nselrows)++] = ancestornodeinfo->rows[i];
#ifndef NDEBUG
            {
               /* check if this row is not featured earlier */
               for (j = 0; j < (*nselrows) - 1; ++j)
               {
                  assert( ancestornodeinfo->rows[i] != (*roworder)[j] );
               }
            }
#endif
         }
      }

      node = SCIPnodeGetParent(node);
   }

   /* row order is in reverse order now, so reverse the array */
   for (i = 0; i < (*nselrows) / 2; ++i)
   {
      /* swap entry i with nselrows - 1 - i */
      j = (*roworder)[i];
      (*roworder)[i] = (*roworder)[(*nselrows) - 1 - i];
      (*roworder)[(*nselrows) - 1 - i] = j;
   }

   return SCIP_OKAY;
}


/** gets rooted path up to node and populates column ordering array */
static
SCIP_RETCODE populateRootedPathColumnOrder(
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   SCIP_NODE*            node,               /**< node considered */
   SCIP_NODE**           rootedpath,         /**< array to populate with the rooted path, must be sufficiently long */
   int*                  colorder,           /**< array to populate with the column order, must be nvars long */
   int*                  colorderinv         /**< array to populate with the inverse column order, must be nvars long */
   )
{
   int i;
   int j;
   int depth;
   BNBNODEINFO* ancestornodeinfo;
   BNBNODEINFO tmpnodeinfo;
   COLSWAP* thiscolswap;

   assert( orbidata != NULL );
   assert( node != NULL );
   assert( rootedpath != NULL );
   assert( colorder != NULL );
   assert( colorderinv != NULL );

   depth = SCIPnodeGetDepth(node);
   i = depth;
   while ( node != NULL )
   {
      assert( SCIPnodeGetDepth(node) == i );
      rootedpath[i--] = node;
      node = SCIPnodeGetParent(node);
   }
   assert( i == -1 );

   for (i = 0; i <= depth; ++i)
   {
      node = rootedpath[i];

      assert( SCIPnodeGetDepth(node) == i );

      /* get the node info of that node */
      tmpnodeinfo.nodenumber = SCIPnodeGetNumber(node);
      ancestornodeinfo = (BNBNODEINFO*) SCIPhashtableRetrieve(orbidata->nodeinfos, (void*) &tmpnodeinfo);

      /* skip nodes that do not imply any row or column swaps */
      if ( ancestornodeinfo == NULL )
         continue;

      /* ncolswaps == 0 iff colswaps == NULL */
      assert( (ancestornodeinfo->ncolswaps == 0) != (ancestornodeinfo->colswaps != NULL) );

      for (j = 0; j < ancestornodeinfo->ncolswaps; ++j)
      {
         int positionfromincolorder;
         int positiontoincolorder;

         thiscolswap = &ancestornodeinfo->colswaps[j];
         assert( thiscolswap->from != thiscolswap->to );  /* there are no trivial swaps in the list */
         assert( thiscolswap->from >= 0 && thiscolswap->from < orbidata->ncols );
         assert( thiscolswap->to >= 0 && thiscolswap->to < orbidata->ncols );

         /* at which column is origcolid? */
         positionfromincolorder = colorderinv[thiscolswap->from];
         assert( positionfromincolorder >= 0 );
         assert( positionfromincolorder < orbidata->ncols );
         assert( colorder[positionfromincolorder] == thiscolswap->from );

         /* at which column is swaporigcolid? */
         positiontoincolorder = colorderinv[thiscolswap->to];
         assert( positiontoincolorder >= 0 );
         assert( positiontoincolorder < orbidata->ncols );
         assert( colorder[positiontoincolorder] == thiscolswap->to );

         /* swap them, also keep track of the inverses */
         colorder[positiontoincolorder] = thiscolswap->from;
         colorder[positionfromincolorder] = thiscolswap->to;
         colorderinv[thiscolswap->from] = positiontoincolorder;
         colorderinv[thiscolswap->to] = positionfromincolorder;
      }
   }

   return SCIP_OKAY;
}

/** at branching decisions, maintains the column swap and potential new rows in the orbitope */
static
SCIP_DECL_EVENTEXEC(eventExecNodeBranched)
{
   ORBITOPEDATA* orbidata;
   SCIP_NODE* node;
   SCIP_NODE* eventnode;
   SCIP_NODE** children;
   SCIP_NODE* childnode;
   SCIP_DOMCHG* domchg;
   SCIP_BOUNDCHG* boundchg;
   SCIP_VAR* var;
   SCIP_VAR** branchvars;
   int maxnbranchvars;
   int nbranchvars;
   int nboundchgs;
   int nchildren;
   int i;
   int j;
   int c;
   int rowid;
   BNBNODEINFO* newnodeinfo;
   SCIP_NODE** rootedpath;

   assert( eventdata != NULL );
   assert( !SCIPinProbing(scip) );

   eventnode = SCIPeventGetNode(event);
   assert( SCIPgetFocusNode(scip) == eventnode );

   orbidata = (ORBITOPEDATA*) eventdata;
   assert( orbidata != NULL );
   assert( orbidata->nrows > 0 );
   assert( orbidata->ncols > 0 );
   assert( orbidata->vars != NULL );
   assert( orbidata->colindexmap != NULL );
   assert( orbidata->rowindexmap != NULL );

   SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );

   /* arrays used within the loop */
   maxnbranchvars = 1;  /* it's a good guess that there's one branching variable, because that's likely the number */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, maxnbranchvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rootedpath, SCIPnodeGetDepth(eventnode)) );

   /* get all variables branched upon (check all branches) */
   nbranchvars = 0;
   for (c = 0; c < nchildren; ++c)
   {
      childnode = children[c];
      domchg = SCIPnodeGetDomchg(childnode);

      /* loop through all bound changes */
      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
      for (i = 0; i < nboundchgs; ++i)
      {
         /* get bound change info */
         boundchg = SCIPdomchgGetBoundchg(domchg, i);
         assert( boundchg != NULL );

         /* branching decisions have to be in the beginning of the bound change array */
         if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
            break;

         /* get corresponding branching variable */
         var = SCIPboundchgGetVar(boundchg);

         /* only variables from the orbitope matrix are of interest */
         if ( ! SCIPhashmapExists(orbidata->rowindexmap, (void*) var) )
            continue;

         /* skip variables that are already stored */
         if ( nbranchvars > 0 )
         {
            for (j = 0; j < nbranchvars; ++j)
            {
               if ( branchvars[j] == var )
                  break;
            }
            /* if the loop above is stopped with `break`, `j < nbranchvars` is not satisfied.
               * then, go to the next iteration
               */
            if ( j < nbranchvars )
               continue;
         }

         /* the variable is not already in the array, so store it */
         if ( nbranchvars >= maxnbranchvars )
         {
            assert( nbranchvars == maxnbranchvars );
            assert( maxnbranchvars > 0 );
            maxnbranchvars = SCIPcalcMemGrowSize(scip, maxnbranchvars + 1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &branchvars, maxnbranchvars) );
         }
         assert( nbranchvars < maxnbranchvars );
         branchvars[nbranchvars++] = var;
      }
   }

   /* skip orbitopes whose variable matrices do not contain any branching variable */
   if ( nbranchvars <= 0 )
      goto FREE;

   SCIP_CALL( SCIPallocBlockMemory(scip, &newnodeinfo) );
   newnodeinfo->nodenumber = SCIPnodeGetNumber(eventnode);
   newnodeinfo->colswaps = NULL;
   newnodeinfo->ncolswaps = 0;
   newnodeinfo->rows = NULL;
   newnodeinfo->nrows = 0;

   /* store data about row ordering */
   if ( orbidata->rowordering != SCIP_ROWORDERING_NONE )
   {
      int* roworder;
      int nselrows;

      assert( orbidata->nrows > 0 );
      assert( orbidata->rowordering == SCIP_ROWORDERING_BRANCHING );

      /* get the present row order up to this node */
      SCIP_CALL( getRowOrder(scip, orbidata, eventnode, &roworder, &nselrows) );
      assert( roworder != NULL );

      /* extend the row fixings with the steps from this node */
      for (i = 0; i < nbranchvars; ++i)
      {
         var = branchvars[i];

         assert( SCIPhashmapExists(orbidata->rowindexmap, (void*) var) ); /* otherwise was not added to branchvars */
         rowid = SCIPhashmapGetImageInt(orbidata->rowindexmap, (void*) var);
         assert( rowid >= 0 );
         assert( rowid < orbidata->nrows );

         /* avoid adding row to row order twice */
         if ( nselrows > 0 )
         {
            for (j = 0; j < nselrows; ++j)
            {
               if ( rowid == roworder[j] )
                  break;
            }
            if ( j < nselrows )  /* if the loop is interrupted */
               continue;
         }

         /* if we end up here, the row index does not appear for any ancestor or the present row order */

         /* append rowid to present roworder */
         roworder[nselrows++] = rowid;

         /* mark that this row index is the new one in the node */
         if ( newnodeinfo->rows == NULL )
         {
            assert( newnodeinfo->nrows == 0 );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newnodeinfo->rows, newnodeinfo->nrows + 1) );
         }
         else
         {
            /* reallocate with linear increments, because we expect 1 branching variable most of the time */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &newnodeinfo->rows,
               newnodeinfo->nrows, newnodeinfo->nrows + 1) );
         }
         newnodeinfo->rows[newnodeinfo->nrows++] = rowid;
      }

      freeRowOrder(scip, orbidata, &roworder);
   }

   /* store data about column ordering */
   if ( orbidata->columnordering != SCIP_COLUMNORDERING_NONE )
   {
      int* colorder;
      int* colorderinv;
      COLSWAP* thiscolswap;
      COLSWAP tmpcolswap;

      assert( orbidata->ncols > 0 );
      SCIP_CALL( SCIPallocBufferArray(scip, &colorder, orbidata->ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colorderinv, orbidata->ncols) );

      /* populate colorder with standard ordering */
      for (i = 0; i < orbidata->ncols; ++i)
         colorder[i] = i;

      /* introduce inverse column ordering */
      for (i = 0; i < orbidata->ncols; ++i)
         colorderinv[i] = i;

      /* get the rooted path
      *
      * We want to iterate through the bound changes in the order of the rooted path to this node.
      */
      node = SCIPnodeGetParent(eventnode);
      if ( node != NULL )
      {
         SCIP_CALL( populateRootedPathColumnOrder(orbidata, node, rootedpath, colorder, colorderinv) );
      }

      /* get the swap for this node */
      for (i = 0; i < nbranchvars; ++i)
      {
         SCIP_CALL( updateColumnOrderWhenBranchingOnColumn(scip, orbidata, colorder,
            colorderinv, branchvars[i], &tmpcolswap) );

         /* skip trivial swaps of columns */
         if ( tmpcolswap.from == tmpcolswap.to )
            continue;

         /* mark that this row index is the new one in the node */
         if ( newnodeinfo->rows == NULL )
         {
            assert( newnodeinfo->nrows == 0 );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newnodeinfo->colswaps, newnodeinfo->ncolswaps + 1) );
         }
         else
         {
            /* reallocate with linear increments, because we expect 1 branching variable most of the time */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &newnodeinfo->colswaps, newnodeinfo->ncolswaps,
               newnodeinfo->ncolswaps + 1) );
         }
         thiscolswap = &(newnodeinfo->colswaps[newnodeinfo->ncolswaps++]);
         thiscolswap->from = tmpcolswap.from;
         thiscolswap->to = tmpcolswap.to;
      }

      SCIPfreeBufferArray(scip, &colorder);
      SCIPfreeBufferArray(scip, &colorderinv);
   }

   /* store updates of row/column order or free memory if no change applied */
   if ( newnodeinfo->nrows > 0 || newnodeinfo->ncolswaps > 0 )
   {
      SCIP_CALL( SCIPhashtableSafeInsert(orbidata->nodeinfos, newnodeinfo) );
   }
   else
   {
      SCIPfreeBlockMemory(scip, &newnodeinfo);
   }

FREE:
   SCIPfreeBufferArray(scip, &rootedpath);
   SCIPfreeBufferArray(scip, &branchvars);

   return SCIP_OKAY;
} /*lint !e715*/


/** at branching decisions, maintains the column swap and potential new rows in the orbitope */
static
SCIP_DECL_EVENTEXEC(eventExec)
{
   switch (SCIPeventGetType(event))
   {
   case SCIP_EVENTTYPE_NODEBRANCHED:
      return eventExecNodeBranched(scip, eventhdlr, event, eventdata);
   default:
      SCIPerrorMessage("Eventhandler " EVENTHDLR_NAME " is called with an unsupported eventtype.\n");
      return SCIP_ERROR;
   }
}


/** returns whether a row contains potential branching variables */
static
SCIP_Bool rowIsBranchRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< pointer to the dynamic orbitopal reduction data */
   ORBITOPEDATA*         orbidata,           /**< symmetry handling data for orbitopal structure */
   int                   rowid               /**< row id for which to check */
   )
{
   SCIP_VAR* var;
#ifndef NDEBUG
   int c;
#endif

   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( orbidata != NULL );
   assert( orbidata->nrows > 0 );
   assert( orbidata->ncols > 0 );
   assert( rowid >= 0 );
   assert( rowid < orbidata->nrows );
   assert( orbidata->vars != NULL );
   assert( orbidata->vars[rowid * orbidata->ncols] );  /* variable in first column must be set */

   /* get the first variable from the row */
   var = orbidata->vars[rowid * orbidata->ncols];

   /* debugging: the variable types in a row should all be the same */
#ifndef NDEBUG
   for (c = 1; c < orbidata->ncols; ++c)
   {
      /* the actual vartypes can be different,
       * for example when an INTEGER vartype turns into BINARY due to bound changes
       */
      assert( vartypeIsBranchRowType(scip, orbireddata, SCIPvarGetType(var)) ==
         vartypeIsBranchRowType(scip, orbireddata, SCIPvarGetType(orbidata->vars[rowid * orbidata->ncols + c])) );
   }
#endif

   return vartypeIsBranchRowType(scip, orbireddata, SCIPvarGetType(var));
}


/** frees orbitope data */
static
SCIP_RETCODE freeOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< pointer to the dynamic orbitopal reduction data */
   ORBITOPEDATA**        orbidata            /**< pointer to orbitope data */
   )
{
   BNBNODEINFO* nodeinfo;
   int i;
   int nentries;
   int nelem;

   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( orbidata != NULL );
   assert( *orbidata != NULL );
   assert( (*orbidata)->vars != NULL );
   assert( (*orbidata)->nrows > 0 );
   assert( (*orbidata)->ncols > 0 );
   assert( (*orbidata)->nrows * (*orbidata)->ncols > 0 );
   assert( SCIPisTransformed(scip) );

   /* free data if orbitopal reduction is dynamic */
   if ( (*orbidata)->columnordering != SCIP_COLUMNORDERING_NONE || (*orbidata)->rowordering != SCIP_ROWORDERING_NONE )
   {
      /* drop event */
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED, orbireddata->eventhdlr,
            (SCIP_EVENTDATA*) *orbidata, -1 ) );

      /* free nodeinfos */
      nentries = SCIPhashtableGetNEntries((*orbidata)->nodeinfos);
      for (i = 0; i < nentries; ++i)
      {
         /* @todo in principle, can deal with memory sparsity by first getting all nodeinfos,
          * then sorting by address and free them in descending order
          */
         nodeinfo = (BNBNODEINFO*) (SCIPhashtableGetEntry((*orbidata)->nodeinfos, i));
         if ( nodeinfo == NULL )
            continue;

         assert( nodeinfo != NULL );
         assert( nodeinfo->nrows > 0 || nodeinfo->ncolswaps > 0 );

         assert( (nodeinfo->ncolswaps == 0) != (nodeinfo->colswaps != NULL) );
         SCIPfreeBlockMemoryArrayNull(scip, &(nodeinfo->colswaps), nodeinfo->ncolswaps);

         assert( (nodeinfo->nrows == 0) != (nodeinfo->rows != NULL) );
         SCIPfreeBlockMemoryArrayNull(scip, &(nodeinfo->rows), nodeinfo->nrows);

         SCIPfreeBlockMemory(scip, &nodeinfo);
      }
      SCIPhashtableFree(&((*orbidata)->nodeinfos));
   }

   /* free index lookup hashsets */
   SCIPhashmapFree(&((*orbidata)->colindexmap));
   SCIPhashmapFree(&((*orbidata)->rowindexmap));

   /* free and release vars */
   nelem = (*orbidata)->nrows * (*orbidata)->ncols;
   assert( nelem > 0 );
   for (i = 0; i < nelem; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*orbidata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArray(scip, &((*orbidata)->vars), (*orbidata)->nrows * (*orbidata)->ncols); /*lint !e647*/

   SCIPfreeBlockMemory(scip, orbidata);

   return SCIP_OKAY;
}


/** adds an orbitope to the orbitopal reduction data */
static
SCIP_RETCODE addOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< pointer to the dynamic orbitopal reduction data */
   SCIP_ROWORDERING      rowordering,        /**< specifies how rows of orbitope are ordered */
   SCIP_COLUMNORDERING   colordering,        /**< specifies how columnss of orbitope are ordered */
   SCIP_VAR**            vars,               /**< variables array, must have size nrows * ncols */
   int                   nrows,              /**< number of rows in orbitope */
   int                   ncols,              /**< number of columns in orbitope */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   )
{
   ORBITOPEDATA* orbidata;
   SCIP_VAR* var;
   int i;
   int rowid;
   int colid;
   int nelem;

   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( orbireddata->eventhdlr != NULL );
   assert( vars != NULL );
   assert( nrows >= 0 );
   assert( ncols >= 0 );

   nelem = nrows * ncols;
   assert( nelem >= 0 );

   /* prevent trivial case with empty orbitope */
   if ( nelem == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   *success = TRUE;

   SCIP_CALL( SCIPallocBlockMemory(scip, &orbidata) );

   orbidata->nrows = nrows;
   orbidata->ncols = ncols;
   orbidata->columnordering = colordering;
   orbidata->rowordering = rowordering;

#ifndef NDEBUG
   orbidata->lastnodenumber = -1;
   orbidata->dbghash = 0;
#endif

   /* variable array enumerates the orbitope matrix row-wise */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &orbidata->vars, nelem) );

   /* create row and column index lookup maps */
   SCIP_CALL( SCIPhashmapCreate(&orbidata->rowindexmap, SCIPblkmem(scip), nrows) );
   SCIP_CALL( SCIPhashmapCreate(&orbidata->colindexmap, SCIPblkmem(scip), ncols) );

   SCIPdebugMessage("Orbitope variables for (%dx%d) orbitope with orbidata %p\n", nrows, ncols, (void*) orbidata);

   /* populate variable array defining orbitope matrix for orbitope data */
   for (i = 0, rowid = 0, colid = 0; i < nelem; ++i, ++colid)
   {
      if ( colid == ncols )
      {
         colid = 0;
         ++rowid;
      }
      assert( nrows > 0 );
      assert( ncols > 0 );
      assert( rowid == i / ncols );
      assert( colid == i % ncols );

      var = vars[i];
      assert( var != NULL );
      assert( SCIPvarIsTransformed(var) );

      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );
      SCIP_CALL( SCIPcaptureVar(scip, var) );

      orbidata->vars[i] = var;

      /* variables cannot be repeated in the variable matrix */
      assert( ! SCIPhashmapExists(orbidata->rowindexmap, var) );
      SCIP_CALL( SCIPhashmapInsertInt(orbidata->rowindexmap, var, rowid) );

      assert( ! SCIPhashmapExists(orbidata->colindexmap, var) );
      SCIP_CALL( SCIPhashmapInsertInt(orbidata->colindexmap, var, colid) );

      SCIPdebugMessage("%4d %4d -> %s\n", rowid, colid, var->name);
   }

   /* count number of branchable rows in orbitope */
   orbidata->nbranchrows = 0;
   /* @todo at getRowData: If nselrows == nbranchrows, append the non-branch rows (like before) */
   for (i = 0; i < nrows; ++i)
   {
      if ( rowIsBranchRow(scip, orbireddata, orbidata, i) )
         ++orbidata->nbranchrows;
   }

   /* cannot add orbitope when already branching */
   assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING ? SCIPgetNNodes(scip) == 0 : TRUE );

   /* possibly create data needed for dynamic orbitopal reduction */
   if ( orbidata->columnordering != SCIP_COLUMNORDERING_NONE || orbidata->rowordering != SCIP_ROWORDERING_NONE )
   {
      /* add the event to store the row and column updates of nodes in the branch-and-bound tree */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED, orbireddata->eventhdlr,
            (SCIP_EVENTDATA*) orbidata, NULL) );

      /* nodeinfos: every node that implies a column swap is represented
       *
       * Assuming at most one branching on every variable implying a column swap, initial hashtable size nelem.
       * In case that there are many more rows than columns, we do not expect too many column swaps.
       */
      SCIP_CALL( SCIPhashtableCreate(&orbidata->nodeinfos, scip->mem->probmem, MIN(16 * ncols + 64, nelem),
            hashGetKeyBnbnodeinfo, hashKeyEqBnbnodeinfo, hashKeyValBnbnodeinfo, NULL) );
   }

   /* resize orbitope array if needed */
   assert( orbireddata->norbitopes >= 0 );
   assert( (orbireddata->norbitopes == 0) == (orbireddata->orbitopes == NULL) );
   assert( orbireddata->norbitopes <= orbireddata->maxnorbitopes );
   if ( orbireddata->norbitopes == orbireddata->maxnorbitopes )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, orbireddata->norbitopes + 1);
      assert( newsize >= 0 );

      if ( orbireddata->norbitopes == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &orbireddata->orbitopes, newsize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &orbireddata->orbitopes, orbireddata->norbitopes, newsize) );
      }

      orbireddata->maxnorbitopes = newsize;
   }
   assert( orbireddata->orbitopes != NULL );
   assert( orbireddata->norbitopes < orbireddata->maxnorbitopes );

   /* add orbitope to orbitopal reduction data */
   assert( orbireddata->norbitopes < orbireddata->maxnorbitopes );
   orbireddata->orbitopes[orbireddata->norbitopes++] = orbidata;

   SCIPdebugMsg(scip, "Added orbitope for orbitopal reduction of size %d by %d\n", nrows, ncols);

   return SCIP_OKAY;
}


/** frees the column order */
static
void freeColumnOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   int**                 colorder,           /**< colorder array that is initialized with the colorder in the dynamic
                                               *  case, of size ncols, and NULL in the static case */
   int**                 colorderinv         /**< array with the inverse column order, of size ncols */
)
{
   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( colorder != NULL );
   assert( colorderinv != NULL );

   if ( orbidata->columnordering == SCIP_COLUMNORDERING_NONE )
   {
      assert( *colorder == NULL );
      assert( *colorderinv == NULL );
      return;
   }
   assert( *colorder != NULL );
   assert( *colorderinv != NULL );

   SCIPfreeBlockMemoryArray(scip, colorder, orbidata->ncols);
   SCIPfreeBlockMemoryArray(scip, colorderinv, orbidata->ncols);
}


/** gets the column order at the node
 *
 *  The column order is (deterministically) dynamically decided based on the policy for column ordering.
 */
static
SCIP_RETCODE getColumnOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   SCIP_NODE*            eventnode,          /**< node where this should be determined at */
   int**                 colorder,           /**< array to populate with column order, of size ncols */
   int**                 colorderinv         /**< array to populate with inverse column order, of size ncols */
   )
{
   SCIP_NODE* node;
   SCIP_NODE** rootedpath;
   int i;
   int depth;
   int ncols;

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( eventnode != NULL );
   assert( colorder != NULL );
   assert( colorderinv != NULL );

   if ( orbidata->columnordering == SCIP_COLUMNORDERING_NONE )
   {
      *colorder = NULL;
      *colorderinv = NULL;
      return SCIP_OKAY;
   }
   ncols = orbidata->ncols;
   assert( ncols > 0 );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, colorder, ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, colorderinv, ncols) );

   /* populate colorder with standard ordering */
   for (i = 0; i < ncols; ++i)
      (*colorder)[i] = i;

   /* introduce inverse column ordering */
   for (i = 0; i < ncols; ++i)
      (*colorderinv)[i] = i;

   /* get the rooted path
    *
    * We want to iterate through the bound changes in the order of the rooted path to this node.
    */
   node = SCIPnodeGetParent(eventnode);
   if ( node != NULL )
   {
      depth = SCIPnodeGetDepth(node);
      SCIP_CALL( SCIPallocBufferArray(scip, &rootedpath, depth + 1) );
      SCIP_CALL( populateRootedPathColumnOrder(orbidata, node, rootedpath, *colorder, *colorderinv) );
      SCIPfreeBufferArray(scip, &rootedpath);
   }

   return SCIP_OKAY;
}


#ifndef NDEBUG
/** checks if the columns of the matrix are lexicographically decreasing, using the specified row and column ordering */
static
void assertIsOrbitopeMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   int*                  roworder,           /**< array with the row order */
   int*                  colorder,           /**< array with the column order */
   SCIP_Real*            matrix,             /**< a matrix */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of cols of matrix */
   int*                  infinitesimal,      /**< array specifying where the infinitesimals are at */
   SCIP_Bool             addinfinitesimals   /**< whether infinitesimals are added (TRUE) or subtracted (FALSE) */
   )
{
   int rowid;
   int colid;
   int idx;
   int origrowid;
   int origcolid;
   int origidx;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( matrix != NULL );
   assert( orbidata != NULL );
   assert( orbidata->vars != NULL );
   assert( nrows >= 0 );
   assert( nrows <= orbidata->nrows );
   assert( ncols >= 0 );
   assert( ncols <= orbidata->ncols );
   assert( infinitesimal != NULL );

   /* respect variable bounds */
   for (rowid = 0; rowid < nrows; ++rowid)
   {
      origrowid = getArrayEntryOrIndex(roworder, rowid);
      for (colid = 0; colid < ncols; ++colid)
      {
         origcolid = getArrayEntryOrIndex(colorder, colid);
         idx = rowid * ncols + colid;
         origidx = origrowid * ncols + origcolid;
         var = orbidata->vars[origidx];
         assert( SCIPsymGE(scip, matrix[idx], SCIPvarGetLbLocal(var)) );
         assert( SCIPsymLE(scip, matrix[idx], SCIPvarGetUbLocal(var)) );
      }
   }

   /* is orbitope */
   for (colid = 0; colid < ncols - 1; ++colid)
   {
      /* compare column colid with colid + 1 */
      for (rowid = 0; rowid < nrows; ++rowid)
      {
         /* entry is >= entry to the right */
         assert( SCIPsymGE(scip, matrix[rowid * ncols + colid], matrix[rowid * ncols + colid + 1]) );

         if ( SCIPsymGT(scip, matrix[rowid * ncols + colid], matrix[rowid * ncols + colid + 1]) )
         {
            /* critical row */
            break;
         }
         else
         {
            /* check for infinitesimal values
             * If infinitesimals are added (lexminface case), then if the left column has a +epsilon,
             * it does not matter whether the right column has +epsilon or not, then the left column is >,
             * due to the axioms x + epsilon > x + epsilon and x + epsilon > x.
             * Analogously, x > x - epsilon and x - epsilon > x - epsilon.
             */
            assert( SCIPsymEQ(scip, matrix[rowid * ncols + colid], matrix[rowid * ncols + colid + 1]) );
            if ( addinfinitesimals
               ? (infinitesimal[colid] == rowid) /* left has +epsilon term */
               : (infinitesimal[colid + 1] == rowid) /* right has -epsilon term */
            )
            {
               /* critical row */
               break;
            }
         }
      }
   }
}
#endif

#ifndef NDEBUG
/** to test if arrays are the same, generates some hash for an array of integers */
static
int debugGetArrayHash(
   int*                  array,              /** array */
   int                   len                 /** array length */
   )
{
   int i;
   unsigned int hash = 0;

   assert( array != NULL );
   assert( len >= 0 );

   for (i = 0; i < len; i++)
   {
      hash ^= (unsigned int) (array[i]);
      hash = (hash << 1) ^ (hash >> 1);
   }

   return (int) hash;
}
#endif

#ifdef SCIP_MORE_DEBUG
/** prints nrows × ncols matrix of floats with 2 decimals */
static
void debugPrintMatrix(
   SCIP_Real*            matrix,             /** matrix, encoded as array enumerating the elements row-wise */
   int                   nrows,              /** number of rows */
   int                   ncols               /** number of rows */
   )
{
   int row;
   int col;

   assert( matrix != NULL );
   assert( nrows >= 0 );
   assert( ncols >= 0 );

   for (row = 0; row < nrows; ++row)
   {
      SCIPdebugPrintf("[");
      for (col = 0; col < ncols; ++col)
      {
         SCIPdebugPrintf(" %+10.2f", matrix[row * ncols + col]);
      }
      SCIPdebugPrintf(" ]\n");
   }
}
#endif


/** gets the column order at the node */
static
SCIP_RETCODE propagateStaticOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   int*                  roworder,           /**< array with the row order (or NULL if identity function is used) */
   int                   nselrows,           /**< number of selected rows */
   int*                  colorder,           /**< array with the column order (or NULL if identity function is used) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nfixedvars          /**< pointer to counter of number of variable domain reductions */
   )
{
   /* @todo also make "nselcols" to allow for colorders smaller than orbidata->ncols */
   SCIP_Real* lexminface = NULL;
   int* lexminepsrow = NULL;
   SCIP_Real* lexmaxface = NULL;
   int* lexmaxepsrow = NULL;
   int nelem;
   int rowid;
   int colid;
   int ncols;
   int origrowid;
   int origcolid;
   int origidx;
   int i;
   int lastunfixed;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool iseq;
   SCIP_Bool success;
   SCIP_VAR* var;
   SCIP_VAR* othervar;

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( orbidata->vars != NULL );
   assert( nselrows >= 0 );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   *infeasible = FALSE;

   assert( orbidata->nrows > 0 );
   assert( orbidata->nrows >= nselrows );
   ncols = orbidata->ncols;
   assert( ncols > 1 );

   /* nothing to propagate without any rows */
   if ( nselrows <= 0 )
      return SCIP_OKAY;

#ifdef SCIP_MORE_DEBUG
   /* print matrix for debugging purposes */
   {
      int k;
      int r;
      SCIP_VAR* thisvar;
      SCIPdebugMessage("Start propagating static orbitope: \n");
      SCIPdebugPrintf(">");
      for (k = 0; k < ncols; ++k)
      {
         SCIPdebugPrintf("%12d ", colorder[k]);
      }
      SCIPdebugPrintf("< (IDs)\n");

      for (r = 0; r < nselrows; ++r)
      {
         SCIPdebugPrintf("[");
         for (k = 0; k < ncols; ++k)
         {
            thisvar = orbidata->vars[roworder[r] * ncols + colorder[k]];
            SCIPdebugPrintf("%4s %+1.2f,%+1.2f ", SCIPvarGetName(thisvar),
               SCIPvarGetLbLocal(thisvar), SCIPvarGetUbLocal(thisvar));
         }
         SCIPdebugPrintf("] (row %d)\n", roworder[r]);
      }
   }
#endif

   nelem = nselrows * ncols;

   /* compute lexmin face */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexminface, nelem) );

   /* array to store for each column at which row we add an infinitesimal value, initially at none (-1) */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexminepsrow, ncols) );
   for (colid = 0; colid < ncols; ++colid)
      lexminepsrow[colid] = -1;

   /* last column takes the minimally possible values. */
   colid = ncols - 1;
   origcolid = getArrayEntryOrIndex(colorder, colid);
   for (rowid = 0, i = colid; rowid < nselrows; ++rowid, i += ncols)
   {
      origrowid = getArrayEntryOrIndex(roworder, rowid);
      origidx = origrowid * ncols + origcolid;
      var = orbidata->vars[origidx];

      assert( i == rowid * ncols + colid );
      assert( var != NULL );

      lexminface[i] = SCIPvarGetLbLocal(var);
   }
   /* all previous columns: One-column replacement algorithm */
   for (colid = ncols - 2; colid >= 0; --colid)
   {
      /* "rowid" of the last unfixed entry whose domain allows for larger values than the current chosen.
       * If there is none, -1. */
      lastunfixed = -1;
      /* whether column "colid" is the same as column "colid + 1" up (but excluding) to "rowid" */
      iseq = TRUE;

      origcolid = getArrayEntryOrIndex(colorder, colid);
      for (rowid = 0, i = colid; rowid < nselrows; ++rowid, i += ncols)
      {
         origrowid = getArrayEntryOrIndex(roworder, rowid);
         origidx = origrowid * ncols + origcolid;
         assert( i == rowid * ncols + colid );

         /* the entry one to the right is not the first column */
         assert( (i + 1) % ncols > 0 );

         var = orbidata->vars[origidx];
         assert( var != NULL );

         if ( iseq )
         {
            /* equality holds up to this row
             * Compare to the entry value on the column immediately right.
             * The value we choose on the left must be at least this.
             * 2 Options:
             * Option 1: The upper bound is smaller. Then we're in an infeasible situation. Resolve as described below.
             * Option 2: The upper bound is greater or equal.
             */
            ub = SCIPvarGetUbLocal(var);

            /* compare to the value in the column right of it */
            if ( SCIPsymLT(scip, ub, lexminface[i + 1]) ||
               ( lexminepsrow[colid + 1] == rowid && SCIPsymEQ(scip, ub, lexminface[i + 1]) ) )
            {
               /* value of this column can only be strictly smaller than the value in the column to its right
                * This may not be possible.
                * Try to repair: Go back to the last row with "room" left, and make the value minimally larger.
                */
               if ( lastunfixed >= 0 )
               {
                  /* repair: return to the last row with "room", and increase the lexmin-value at that row. */
                  assert( SCIPsymEQ(scip, lexminface[lastunfixed * ncols + colid],
                     lexminface[lastunfixed * ncols + colid + 1]) );
                  othervar = orbidata->vars[getArrayEntryOrIndex(roworder, lastunfixed) * ncols + origcolid];
                  switch (SCIPvarGetType(othervar))
                  {
                  case SCIP_VARTYPE_BINARY:
                  case SCIP_VARTYPE_IMPLINT:
                  case SCIP_VARTYPE_INTEGER:
                     /* discrete type with unit steps: Add one to the bound. */
                     /* @todo @question Are variable bounds for SCIP_VARTYPE_IMPLINT always integral? */
                     assert( SCIPisIntegral(scip, lexminface[lastunfixed * ncols + colid]) );
                     lexminface[lastunfixed * ncols + colid] += 1.0;
                     assert( SCIPisIntegral(scip, lexminface[lastunfixed * ncols + colid]) );
                     assert( SCIPsymLE(scip, lexminface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     break;
                  case SCIP_VARTYPE_CONTINUOUS:
                     /* continuous type, so add an infinitesimal value to the bound */
                     assert( SCIPsymLE(scip, lexminface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     assert( lexminepsrow[colid] == -1 );
                     lexminepsrow[colid] = lastunfixed;
                     break;
                  default:
                     return SCIP_ERROR;
                  }
                  /* now row "lastunfixed" is greater. Restart from here. */
                  iseq = FALSE;
                  rowid = lastunfixed; /* the next iteration considers "lastunfixed + 1" */
                  i = rowid * ncols + colid;
                  continue;
               }
               else
               {
                  /* cannot repair. It is infeasible. */
                  *infeasible = TRUE;
                  SCIPdebugMessage("Cannot repair infeasibility for column %d (original: %d), min\n", colid, origcolid);
                  goto FREE;
               }
            }
            else
            {
               assert( SCIPsymGE(scip, ub, lexminface[i + 1]) );
               lb = SCIPvarGetLbLocal(var);
               assert( SCIPsymLE(scip, lb, ub) );
               lexminface[i] = MAX(lexminface[i + 1], lb);
               assert( SCIPsymGE(scip, lexminface[i], lexminface[i + 1]) );

               /* are we still equal? */
               if ( SCIPsymGT(scip, lexminface[i], lexminface[i + 1]) )
                  iseq = FALSE;
               else if ( lexminepsrow[colid + 1] == rowid )
               {
                  assert( SCIPsymEQ(scip, lexminface[i], lexminface[i + 1]) );
                  assert( SCIPvarGetType(orbidata->vars[getArrayEntryOrIndex(roworder, rowid) * ncols + origcolid])
                     == SCIP_VARTYPE_CONTINUOUS );
                  assert( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS );
                  /* right column (colid+1) has value x + epsilon, left column (colid) has value x, now
                   * must also become  x + epsilon in order to be larger or equal
                   * by axioms, we can squeeze infinitesimals between one other; epsilon > epsilon.
                   */
                  iseq = FALSE;
                  assert( lexminepsrow[colid] == -1 );
                  lexminepsrow[colid] = rowid;
               }

               /* is there room left to increase this variable further? */
               switch (SCIPvarGetType(var))
               {
               case SCIP_VARTYPE_BINARY:
               case SCIP_VARTYPE_IMPLINT:
               case SCIP_VARTYPE_INTEGER:
                  /* discrete type with unit steps: Add one to the bound. */
                  /* @todo @question Are variable bounds for SCIP_VARTYPE_IMPLINT always integral? */
                  /* @todo in principle, this can be made more tight using the hole-lists... */
                  assert( SCIPisIntegral(scip, lexminface[i]) );
                  if ( SCIPsymLE(scip, lexminface[i] + 1.0, ub) )
                     lastunfixed = rowid;
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous type: if we can add an infinitesimal value to the current lexminface[i] value,
                   * mark row as 'lastunfixed'
                   */
                  if ( SCIPsymLT(scip, lexminface[i], ub) )
                     lastunfixed = rowid;
                  break;
               default:
                  return SCIP_ERROR;
               }
            }
         }
         else
         {
            /* there had been a row before which breaks the equality-condition, choose minimally possible value */
            lexminface[i] = SCIPvarGetLbLocal(var);
         }
      }
   }

#ifndef NDEBUG
   /* sanity checks */
   assertIsOrbitopeMatrix(scip, orbidata, roworder, colorder, lexminface, nselrows, ncols, lexminepsrow, TRUE);
#endif

   /* compute lexmax face */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxface, nelem) );

   /* array to store for each column at which row we subtract an infinitesimal value, initially at none (-1) */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxepsrow, ncols) );
   for (colid = 0; colid < ncols; ++colid)
      lexmaxepsrow[colid] = -1;

   /* first column, fill all unfixed entries with maximally possible values */
   colid = 0;
   origcolid = getArrayEntryOrIndex(colorder, colid);
   for (rowid = 0, i = colid; rowid < nselrows; ++rowid, i += ncols)
   {
      origrowid = getArrayEntryOrIndex(roworder, rowid);
      origidx = origrowid * ncols + origcolid;
      var = orbidata->vars[origidx];

      assert( i == rowid * ncols + colid );
      assert( var != NULL );

      lexmaxface[i] = SCIPvarGetUbLocal(var);
   }
   /* all next columns: One-column replacement algorithm */
   for (colid = 1; colid < ncols; ++colid)
   {
      /* "rowid" of the last unfixed entry whose domain allows for smaller values than the current chosen.
       * If there is none, -1. */
      lastunfixed = -1;
      /* whether column "colid" is the same as column "colid - 1" up (but excluding) to "rowid" */
      iseq = TRUE;

      origcolid = getArrayEntryOrIndex(colorder, colid);
      for (rowid = 0, i = colid; rowid < nselrows; ++rowid, i += ncols)
      {
         origrowid = getArrayEntryOrIndex(roworder, rowid);
         origidx = origrowid * ncols + origcolid;
         assert( i == rowid * ncols + colid );

         /* the entry one to the left is not the last column, i.e., this column cannot be the first column */
         assert( i % ncols > 0 );

         var = orbidata->vars[origidx];
         assert( var != NULL );

         if ( iseq )
         {
            /* equality holds up to this row
             * Compare to the entry value on the column immediately left.
             * The value we choose on the right must be at most this.
             * 2 Options:
             * Option 1: The lower bound is larger. Then we're in an infeasible situation. Resolve as described below.
             * Option 2: The lower bound is smaller or equal.
             */
            lb = SCIPvarGetLbLocal(var);

            /* compare to the value in the column left of it */
            if ( SCIPsymGT(scip, lb, lexmaxface[i - 1]) ||
               ( lexmaxepsrow[colid - 1] == rowid && SCIPsymEQ(scip, lb, lexmaxface[i - 1]) ) )
            {
               /* value of this column can only be strictly larger than the value in the column to its left
                * This may not be possible.
                * Try to repair: Go back to the last row with "room" left, and make the value minimally smaller.
                */
               if ( lastunfixed >= 0 )
               {
                  /* repair: return to the last row with "room", and decrease the lexmax-value at that row. */
                  assert( SCIPsymEQ(scip, lexmaxface[lastunfixed * ncols + colid],
                     lexmaxface[lastunfixed * ncols + colid - 1]) );
                  othervar = orbidata->vars[getArrayEntryOrIndex(roworder, lastunfixed) * ncols + origcolid];
                  switch (SCIPvarGetType(othervar))
                  {
                  case SCIP_VARTYPE_BINARY:
                  case SCIP_VARTYPE_IMPLINT:
                  case SCIP_VARTYPE_INTEGER:
                     /* discrete type with unit steps: Remove one from the lexmax-value. */
                     /* @todo @question Are variable bounds for SCIP_VARTYPE_IMPLINT always integral? */
                     assert( SCIPisIntegral(scip, lexmaxface[lastunfixed * ncols + colid]) );
                     lexmaxface[lastunfixed * ncols + colid] -= 1.0;
                     assert( SCIPisIntegral(scip, lexmaxface[lastunfixed * ncols + colid]) );
                     assert( SCIPsymGE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetLbLocal(othervar)) );
                     assert( SCIPsymLE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     break;
                  case SCIP_VARTYPE_CONTINUOUS:
                     /* continuous type, so subtract an infinitesimal value to the bound */
                     assert( SCIPsymGE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetLbLocal(othervar)) );
                     assert( SCIPsymLE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     assert( lexmaxepsrow[colid] == -1 );
                     lexmaxepsrow[colid] = lastunfixed;
                     break;
                  default:
                     return SCIP_ERROR;
                  }
                  /* now row "lastunfixed" is greater. Restart from here. */
                  iseq = FALSE;
                  rowid = lastunfixed; /* the next iteration considers "lastunfixed + 1" */
                  i = rowid * ncols + colid;
                  continue;
               }
               else
               {
                  /* cannot repair. It is infeasible. */
                  *infeasible = TRUE;
                  SCIPdebugMessage("Cannot repair infeasibility for column %d (original: %d), max\n", colid, origcolid);
                  goto FREE;
               }
            }
            else
            {
               assert( SCIPsymLE(scip, lb, lexmaxface[i - 1]) );
               ub = SCIPvarGetUbLocal(var);
               assert( SCIPsymLE(scip, lb, ub) );
               lexmaxface[i] = MIN(lexmaxface[i - 1], ub);
               assert( SCIPsymGE(scip, lexmaxface[i - 1], lexmaxface[i]) );

               /* are we still equal? */
               if ( SCIPsymGT(scip, lexmaxface[i - 1], lexmaxface[i]) )
                  iseq = FALSE;
               else if ( lexmaxepsrow[colid - 1] == rowid )
               {
                  assert( SCIPsymEQ(scip, lexmaxface[i - 1], lexmaxface[i]) );
                  assert( SCIPvarGetType(orbidata->vars[getArrayEntryOrIndex(roworder, rowid) * ncols + origcolid])
                     == SCIP_VARTYPE_CONTINUOUS );
                  assert( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS );
                  /* left column (colid-1) has value x - epsilon, right column (colid) has value x, now
                   * must also become  x - epsilon in order to be larger or equal
                   * by axioms, we can squeeze infinitesimals between one other; epsilon > epsilon.
                   */
                  iseq = FALSE;
                  assert( lexmaxepsrow[colid] == -1 );
                  lexmaxepsrow[colid] = rowid;
               }

               /* is there room left to decrease this variable further? */
               switch (SCIPvarGetType(var))
               {
               case SCIP_VARTYPE_BINARY:
               case SCIP_VARTYPE_IMPLINT:
               case SCIP_VARTYPE_INTEGER:
                  /* discrete type with unit steps: Remove one from the lexmax-value. */
                  /* @todo @question Are variable bounds for SCIP_VARTYPE_IMPLINT always integral? */
                  /* @todo in principle, this can be made more tight using the hole-lists... */
                  assert( SCIPisIntegral(scip, lexmaxface[i]) );
                  if ( SCIPsymGE(scip, lexmaxface[i] - 1.0, lb) )
                     lastunfixed = rowid;
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous type: if we can subtract an infinitesimal value to the current lexmaxface[i] value,
                   * mark row as 'lastunfixed'
                   */
                  if ( SCIPsymGT(scip, lexmaxface[i], lb) )
                     lastunfixed = rowid;
                  break;
               default:
                  return SCIP_ERROR;
               }
            }
         }
         else
         {
            /* there had been a row before which breaks the equality-condition, choose maximally possible value */
            lexmaxface[i] = SCIPvarGetUbLocal(var);
         }
      }
   }

#ifndef NDEBUG
   /* sanity checks */
   assertIsOrbitopeMatrix(scip, orbidata, roworder, colorder, lexmaxface, nselrows, ncols, lexmaxepsrow, FALSE);
#endif

#ifdef SCIP_MORE_DEBUG
   /* show lexmin and lexmax face */
   SCIPdebugMessage("Lex min face\n");
   debugPrintMatrix(lexminface, nselrows, ncols);
   SCIPdebugMessage("Lex max face\n");
   debugPrintMatrix(lexmaxface, nselrows, ncols);
#endif

   /* compare the two column-wise and apply domain reductions */
   for (colid = 0; colid < ncols; ++colid)
   {
      for (rowid = 0, i = colid; rowid < nselrows; ++rowid, i += ncols)
      {
         assert( i == rowid * ncols + colid );

         /* get var */
         origrowid = getArrayEntryOrIndex(roworder, rowid);
         origcolid = getArrayEntryOrIndex(colorder, colid);
         origidx = origrowid * ncols + origcolid;
         var = orbidata->vars[origidx];

         if ( SCIPsymEQ(scip, lexminface[i], lexmaxface[i]) )
         {
            /* tighten LB and UB to same value (i.e. fixing) */
            SCIP_CALL( SCIPtightenVarLb(scip, var, lexminface[i], FALSE, infeasible, &success) );
            if ( success )
            {
               SCIPdebugMessage("Fixing variable LB %12s (%3d,%3d) to %5.2f\n", var->name, rowid, colid, lexminface[i]);
               *nfixedvars += 1;
            }
            else
            {
               SCIPdebugMessage("Fixing variable LB %12s (%3d,%3d) to %5.2f (no success)\n", var->name, rowid, colid,
                 lexminface[i]);
            }
            if ( *infeasible )
            {
               SCIPdebugMessage("Detected infeasibility fixing variable %12s (%3d,%3d) to %5.2f\n",
                  var->name, rowid, colid, lexminface[i]);
               goto FREE;
            }

            SCIP_CALL( SCIPtightenVarUb(scip, var, lexminface[i], FALSE, infeasible, &success) );
            if ( success )
            {
               SCIPdebugMessage("Fixing variable UB %12s (%3d,%3d) to %5.2f\n", var->name, rowid, colid, lexminface[i]);
               *nfixedvars += 1;
            }
            else
            {
               SCIPdebugMessage("Fixing variable UB %12s (%3d,%3d) to %5.2f (no success)\n", var->name, rowid, colid,
                 lexminface[i]);
            }
            if ( *infeasible )
            {
               SCIPdebugMessage("Detected infeasibility fixing variable %12s (%3d,%3d) to %5.2f\n",
                  var->name, rowid, colid, lexminface[i]);
               goto FREE;
            }
         }
         else
         {
            /* This is the row index where the min- and max-face have a different value for this column entry.
             * Update the lower bound and upper bound */

            /* lower bound, based on lexminface */
            SCIP_CALL( SCIPtightenVarLb(scip, var, lexminface[i], FALSE, infeasible, &success) );
            if ( success )
            {
               SCIPdebugMessage("Restricting variable LB %12s (%3d,%3d) to %5.2f\n", var->name, rowid, colid,
                  lexminface[i]);
               *nfixedvars += 1;
            }
            else
            {
               SCIPdebugMessage("Restricting variable LB %12s (%3d,%3d) to %5.2f (no success)\n", var->name,
                 rowid, colid, lexminface[i]);
            }
            if ( *infeasible )
            {
               SCIPdebugMessage("Detected infeasibility restricting variable LB %12s (%3d,%3d) to %5.2f\n",
                  var->name, rowid, colid, lexminface[i]);
               goto FREE;
            }

            /* upper bound, based on lexmaxface */
            SCIP_CALL( SCIPtightenVarUb(scip, var, lexmaxface[i], FALSE, infeasible, &success) );
            if ( success )
            {
               SCIPdebugMessage("Restricting variable UB %12s (%3d,%3d) to %5.2f\n", var->name, rowid, colid,
                  lexmaxface[i]);
               *nfixedvars += 1;
            }
            else
            {
               SCIPdebugMessage("Restricting variable UB %12s (%3d,%3d) to %5.2f (no success)\n", var->name,
                 rowid, colid, lexmaxface[i]);
            }
            if ( *infeasible )
            {
               SCIPdebugMessage("Detected infeasibility restricting variable UB %12s (%3d,%3d) to %5.2f\n",
                  var->name, rowid, colid, lexmaxface[i]);
               goto FREE;
            }
            break;
         }
      }
   }

FREE:
   SCIPfreeBufferArrayNull(scip, &lexmaxepsrow);
   SCIPfreeBufferArrayNull(scip, &lexmaxface);
   SCIPfreeBufferArrayNull(scip, &lexminepsrow);
   SCIPfreeBufferArrayNull(scip, &lexminface);

   return SCIP_OKAY;
}


/** propagation method for a single orbitope matrix */
static
SCIP_RETCODE propagateOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nfixedvars          /**< pointer to store the number of found domain reductions */
   )
{
   SCIP_NODE* focusnode;
   int* roworder;
   int nselrows;
   int* colorder;
   int* colorderinv;

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   *nfixedvars = 0;
   *infeasible = FALSE;

   assert( orbidata->ncols > 0 );
   assert( orbidata->nrows > 0 );

   focusnode = SCIPgetFocusNode(scip);
   assert( focusnode != NULL );

   /* get row order */
   SCIP_CALL( getRowOrder(scip, orbidata, focusnode, &roworder, &nselrows) );
   assert( nselrows >= 0 );
   assert( nselrows <= orbidata->nrows );
   if ( nselrows == 0 )
      goto FREEROWS;

   /* get column order */
   SCIP_CALL( getColumnOrder(scip, orbidata, focusnode, &colorder, &colorderinv) );

#ifndef NDEBUG
   /* DEBUG: if propagation is repeated in the same node, the same column order and row order is needed */
   /* @todo: performance: move roworder and colorder to orbidata, then re-use */
   {
      int colhash = (colorder == NULL) ? 1 : debugGetArrayHash(colorder, orbidata->ncols);
      int rowhash = (roworder == NULL) ? 0 : debugGetArrayHash(roworder, nselrows);
      int hash = colhash ^ rowhash;

#ifdef SCIP_DEBUG
      SCIPdebugPrintf("Col hash %32d; Row hash %32d; Hash %32d\n", colhash, rowhash, hash);
      {
         SCIP_NODE* tmpnode;
         tmpnode = SCIPgetFocusNode(scip);
         while ( tmpnode != NULL )
         {
            int nbranchings, nconsprop, nprop;
            SCIPnodeGetDomchg(tmpnode);
            SCIPnodeGetNDomchg(tmpnode, &nbranchings, &nconsprop, &nprop);
            SCIPdebugPrintf("  node %lld: (%d, %d, %d) \n", tmpnode->number, nbranchings, nconsprop, nprop);
            tmpnode = SCIPnodeGetParent(tmpnode);
         }
      }
#endif

      assert( SCIPgetCurrentNode(scip) == SCIPgetFocusNode(scip) );  /* no probing */
      if ( orbidata->lastnodenumber == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
      {
         assert( orbidata->dbghash == hash );
      }
      orbidata->dbghash = hash;
   }
   orbidata->lastnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
#endif

   SCIP_CALL( propagateStaticOrbitope(scip, orbidata, roworder, nselrows, colorder, infeasible, nfixedvars) );

   freeColumnOrder(scip, orbidata, &colorder, &colorderinv);
FREEROWS:
   freeRowOrder(scip, orbidata, &roworder);

#ifdef SCIP_MORE_DEBUG
   SCIPdebugPrintf("\n\n");
#endif

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** gets the number of reductions */
SCIP_RETCODE SCIPorbitopalReductionGetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< orbitopal reduction data structure */
   int*                  nred,               /**< total number of reductions applied */
   int*                  ncutoff             /**< total number of cutoffs applied */
   )
{
   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( nred != NULL );

   *nred = orbireddata->nred;
   *ncutoff = orbireddata->ncutoff;

   return SCIP_OKAY;
}


/** prints orbitopal reduction data */
SCIP_RETCODE SCIPorbitopalReductionPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata        /**< orbitopal reduction data structure */
   )
{
   int i;

   assert( scip != NULL );
   assert( orbireddata != NULL );

   if ( orbireddata->norbitopes == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   orbitopal reduction:       no components\n");
      return SCIP_OKAY;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   orbitopal reduction:     %4d components: ", orbireddata->norbitopes);
   for (i = 0; i < orbireddata->norbitopes; ++i)
   {
      if ( i > 0 )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", ");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "%dx%d", orbireddata->orbitopes[i]->nrows, orbireddata->orbitopes[i]->ncols);
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");

   return SCIP_OKAY;
}


/** propagates orbitopal reduction */
SCIP_RETCODE SCIPorbitopalReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< orbitopal reduction data structure */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nred,               /**< pointer to store the number of domain reductions */
   SCIP_Bool*            didrun              /**< a global pointer maintaining if any symmetry propagator has run
                                              *   only set this to TRUE when a reduction is found, never set to FALSE */
   )
{
   ORBITOPEDATA* orbidata;
   int c;
   int thisfixedvars;

   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( (orbireddata->norbitopes == 0) == (orbireddata->orbitopes == NULL) );
   assert( infeasible != NULL );
   assert( nred != NULL );

   *infeasible = FALSE;
   *nred = 0;

   /* @todo Can the following be removed? */
   /* @todo shouldn't this be SCIPallowWeakDualReds, since we do not regard the objective */
   if ( ! SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   /* cannot do anything during probing
    * @todo can find deductions for the probing node, maybe?
    */
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* propagate all orbitopes */
   for (c = 0; c < orbireddata->norbitopes; ++c)
   {
      orbidata = orbireddata->orbitopes[c];
      assert( orbidata != NULL );

      SCIP_CALL( propagateOrbitope(scip, orbidata, infeasible, &thisfixedvars) );
      SCIPdebugMessage("Found %d reductions during orbitopal reduction for orbitope %d\n", thisfixedvars, c);
      *nred += thisfixedvars;

      /* a symmetry propagator has ran, so set didrun to TRUE */
      *didrun = TRUE;

      /* stop if we find infeasibility in one of the methods */
      if ( *infeasible )
      {
         SCIPdebugMessage("Detected infeasibility during orbitopal reduction for orbitope %d\n", c);
         break;
      }
   }

   orbireddata->nred += *nred;
   if ( *infeasible )
      ++orbireddata->ncutoff;

   return SCIP_OKAY;
}

/** adds orbitopal component to orbitopal symmetry handler */
SCIP_RETCODE SCIPorbitopalReductionAddOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< orbitopal reduction data structure */
   SCIP_ROWORDERING      rowordering,        /**< specifies how rows of orbitope are ordered */
   SCIP_COLUMNORDERING   colordering,        /**< specifies how columnss of orbitope are ordered */
   SCIP_VAR**            vars,               /**< matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows */
   int                   ncols,              /**< number of columns */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   )
{
   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );

   /* dynamic symmetry reductions cannot be performed on original problem */
   assert( SCIPisTransformed(scip) );

   /* if this is the first time adding an orbitope, check if the nonlinear conshlr exists */
   if ( !orbireddata->conshdlr_nonlinear_checked )
   {
      orbireddata->conshdlr_nonlinear = SCIPfindConshdlr(scip, "nonlinear");
      orbireddata->conshdlr_nonlinear_checked = TRUE;
   }

   /* create orbitope data */
   SCIP_CALL( addOrbitope(scip, orbireddata, rowordering, colordering, vars, nrows, ncols, success) );

   return SCIP_OKAY;
}


/** resets orbitopal reduction data structure (clears all orbitopes) */
SCIP_RETCODE SCIPorbitopalReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata        /**< pointer to orbitopal reduction structure to populate */
   )
{
   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( orbireddata->norbitopes >= 0 );
   assert( (orbireddata->norbitopes == 0) == (orbireddata->orbitopes == NULL) );
   assert( orbireddata->norbitopes <= orbireddata->maxnorbitopes );
   assert( orbireddata->eventhdlr != NULL );

   /* free orbitopes that are added */
   while (orbireddata->norbitopes > 0)
   {
      SCIP_CALL( freeOrbitope(scip, orbireddata, &(orbireddata->orbitopes[--orbireddata->norbitopes])) );
   }
   assert( orbireddata->norbitopes == 0 );
   SCIPfreeBlockMemoryArrayNull(scip, &orbireddata->orbitopes, orbireddata->maxnorbitopes);
   orbireddata->orbitopes = NULL;
   orbireddata->maxnorbitopes = 0;

   return SCIP_OKAY;
}


/** frees orbitopal reduction data */
SCIP_RETCODE SCIPorbitopalReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA** orbireddata       /**< pointer to orbitopal reduction structure to populate */
   )
{
   assert( scip != NULL );
   assert( orbireddata != NULL );
   assert( *orbireddata != NULL );

   SCIP_CALL( SCIPorbitopalReductionReset(scip, *orbireddata) );

   SCIPfreeBlockMemory(scip, orbireddata);
   return SCIP_OKAY;
}


/** initializes structures needed for orbitopal reduction
 *
 *  This is only done exactly once.
 */
SCIP_RETCODE SCIPincludeOrbitopalReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA** orbireddata       /**< pointer to orbitopal reduction structure to populate */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert( scip != NULL );
   assert( orbireddata != NULL );

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeOrbitopalReduction", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* create orbitope handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, orbireddata) );

   /* default column ordering param */
   SCIP_CALL( SCIPaddIntParam(scip, "propagating/symmetry/" SYMHDLR_NAME "/columnordering",
         "The column ordering variant, respects enum SCIP_ColumnOrdering.",
         (int*) &(*orbireddata)->defaultcolumnordering, TRUE, DEFAULT_COLUMNORDERING, 0, 4,
         NULL, NULL) ); /*lint !e641*/

   /* initialize event handler. */
   assert( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL );
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExec, NULL) );
   assert( eventhdlr != NULL );
   (*orbireddata)->eventhdlr = eventhdlr;

   /* orbitopes array */
   (*orbireddata)->orbitopes = NULL;
   (*orbireddata)->norbitopes = 0;
   (*orbireddata)->maxnorbitopes = 0;

   /* conshdlr nonlinear */
   (*orbireddata)->conshdlr_nonlinear = NULL;
   (*orbireddata)->conshdlr_nonlinear_checked = FALSE;

   /* counter of total number of reductions and cutoffs */
   (*orbireddata)->nred = 0;
   (*orbireddata)->ncutoff = 0;

   return SCIP_OKAY;
}


/** returns the default column ordering */
SCIP_COLUMNORDERING SCIPorbitopalReductionGetDefaultColumnOrdering(
   SCIP_ORBITOPALREDDATA* orbireddata        /**< pointer to orbitopal reduction structure to populate */
   )
{
   assert( orbireddata != NULL );

   return orbireddata->defaultcolumnordering;
}
