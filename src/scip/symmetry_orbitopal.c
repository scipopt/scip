/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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
 * @author Christopher Hojny
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
#include <ctype.h>
#include <string.h>
#include <symmetry/type_symmetry.h>

#include <memory.h>

/* constraint handler properties */
#define EVENTHDLR_NAME         "symmetry_orbitopal_eventhdlr"
#define EVENTHDLR_DESC         "event handler for maintaining the branch-and-bound tree"
#define DEFAULT_COLUMNORDERING SCIP_COLUMNORDERING_MEDIAN /**< the column ordering variant */

/*
 * Data structures
 */

/* local declarations for structs used in this file */
struct SuborbitopeData;
typedef struct SuborbitopeData SUBORBITOPEDATA;
struct OrbitopeData;
typedef struct OrbitopeData ORBITOPEDATA;

/** data for the sub-orbitopes */
struct SuborbitopeData
{
   ORBITOPEDATA*         masterorbitopedata; /**< the master orbitope */
   SCIP_VAR**            sovars;             /**< matrix of this suborbitope */
   SCIP_HASHMAP*         rowindexmap;        /**< map of variables to row index in suborbitope matrix */
   SCIP_HASHMAP*         colindexmap;        /**< map of variables to column index in suborbitope matrix */
   int                   ncols;              /**< number of columns */
#ifndef NDEBUG
   int                   lastnodenumber;     /**< the last node number where the row and column order is computed */
   int                   dbghash;            /**< a hash for the column order in the last iteration */
#endif
   SCIP_HASHTABLE*       nodeinfos;          /**< constraint information per branch-and-bound tree node */
   SCIP_ORBITOPECOLUMNORDERING columnordering; /**< policy for the column ordering */
};

/** constraint data for orbitope constraints */
struct OrbitopeData
{
   SCIP_VAR**            vars;               /**< matrix of variables on which the symmetry acts */
   int                   nrows;              /**< number of rows */
   int                   ncols;              /**< number of columns */
   int                   nbranchrows;        /**< number of rows containing variables potentially used for branching */
   SUBORBITOPEDATA*      suborbitopes;       /**< array of the suborbitopes */
   int                   nsuborbitopes;      /**< number of elements in suborbitope array */
};

/** constraint handler data */
struct SCIP_OrbitopalFixingData
{
   SCIP_ORBITOPECOLUMNORDERING defaultcolumnordering; /**< default policy for the column ordering */
   SCIP_EVENTHDLR*       eventhdlr;          /**< pointer to the event handler for managing the branching tree */
   ORBITOPEDATA**        orbitopes;          /**< array of pointers to orbitope datas */
   int                   norbitopes;         /**< number of orbitope datas in array */
   int                   maxnorbitopes;      /**< allocated orbitopes array size */
   SCIP_CONSHDLR*        conshdlr_nonlinear; /**< nonlinear constraint handler, used for determining if branchingvar */
   SCIP_Bool             conshdlr_nonlinear_checked; /**< nonlinear constraint handler is already added? */
};

/*
 * Local methods
 */

/** helper functions for LT, GE, LE, GE, EQ, that do take infinity-values into account */
static
SCIP_Bool EQ(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return TRUE;
   if ( inf1 != inf2 )
      return FALSE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return TRUE;
   if ( minf1 != minf2 )
      return FALSE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisEQ(scip, val1, val2);
}

static
SCIP_Bool LE(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return TRUE;
   if ( !inf1 && inf2 )
      return TRUE;
   if ( inf1 && !inf2 )
      return FALSE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return TRUE;
   if ( !minf1 && minf2 )
      return FALSE;
   if ( minf1 && !minf2 )
      return TRUE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisLE(scip, val1, val2);
}

static
SCIP_Bool GE(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return TRUE;
   if ( !inf1 && inf2 )
      return FALSE;
   if ( inf1 && !inf2 )
      return TRUE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return TRUE;
   if ( !minf1 && minf2 )
      return TRUE;
   if ( minf1 && !minf2 )
      return FALSE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisGE(scip, val1, val2);
}

static
SCIP_Bool LT(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return FALSE;
   if ( !inf1 && inf2 )
      return TRUE;
   if ( inf1 && !inf2 )
      return FALSE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return FALSE;
   if ( !minf1 && minf2 )
      return FALSE;
   if ( minf1 && !minf2 )
      return TRUE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisLT(scip, val1, val2);
}

static
SCIP_Bool GT(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return FALSE;
   if ( !inf1 && inf2 )
      return FALSE;
   if ( inf1 && !inf2 )
      return TRUE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return FALSE;
   if ( !minf1 && minf2 )
      return TRUE;
   if ( minf1 && !minf2 )
      return FALSE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisGT(scip, val1, val2);
}


/** get whether a variable type is a branchrow-type */
static
SCIP_Bool vartypeIsBranchRowType(
   SCIP_ORBITOPALFIXINGDATA* orbifixdata,    /**< pointer to the dynamic orbitopal fixing data */
   SCIP_VARTYPE          vartype             /**< var type */
)
{
   switch (vartype)
   {
   case SCIP_VARTYPE_BINARY:
   case SCIP_VARTYPE_INTEGER:
      return TRUE;
   case SCIP_VARTYPE_CONTINUOUS:
   case SCIP_VARTYPE_IMPLINT:
      /* potential branching variables if nonlinear constraints exist */
      assert( orbifixdata->conshdlr_nonlinear_checked );
      return orbifixdata->conshdlr_nonlinear == NULL ? FALSE : 
         SCIPconshdlrGetNActiveConss(orbifixdata->conshdlr_nonlinear) > 0;
   default:
      assert( FALSE );
      /* resolve compiler warning: no asserts in optimized mode */
      return FALSE;
   }
}


/** container for column index permutations */
struct colswap
{
   int                   from;               /**< from which column ID */
   int                   to;                 /**< to which column ID */
};
typedef struct colswap COLSWAP;

/** information stored for branch-and-bound nodes */
struct bnbnodeinfo
{
   SCIP_Longint          nodenumber;         /**< node number of the branch-and-bound tree node */
   COLSWAP*              colswaps;           /**< list containing column swaps by node branching decisions */
   int                   ncolswaps;          /**< number of elements in colswaps. ncolswaps == 0 <=> colswaps == NULL */
   int*                  rows;               /**< list of new variable rows by node branching decisions */
   int                   nrows;              /**< number of new variable rows added. nrows == 0 <=> rows == NULL */
};
typedef struct bnbnodeinfo BNBNODEINFO;

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


/** test if two columns are symmetrically equivalent */
static
SCIP_Bool testColumnsAreSymmetricalyEquivalent(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< suborbitope information */
   int                   col1,               /**< first column to compare */
   int                   col2                /**< second column to compare */
)
{
   int i;
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   /* test: are col1 and col2 the same columns w.r.t. the variable domain restrictions?
    * Performance: only need to check the columns in roworder.
    * Note: if no branching has happened on a row, then the complete symmetry handling assumption
    * yields that then the columns are identical only if the columns restricted to
    * the branched rows are identical.
    * In particular, this could mean that these row variables are domain-reduced or not.
    */
   /* the performance improvement assumes the symmetry assumption.
    * This assumption is not guaranteed, so a different result could be found.
    *
    * It could go wrong in the following way: Suppose branch-and-bound path root->b1->b2->b3, branching on a
    * variable from row r1, r2, r3 respectively, in column c1, c2, c3. Suppose that entry (r3,c2) has a
    * different value. Then, on the branch b1->b2, nselrows does not contain row r3, so c2 and c3 are
    * identical. But, on branch b2->b3, when "replaying" the branching steps, it will regard c2 and c3 as
    * different on the branch b1->b2. A different choice could be made.
    *
    * For now: Check all rows.
    */
   for (i = 0; i < sodata->masterorbitopedata->nrows; ++i)
   {
      var1 = sodata->sovars[i * sodata->ncols + col1];
      var2 = sodata->sovars[i * sodata->ncols + col2];
      /* if variable bounds differ: columns c and origcolid are not the same */
      if (
         (! EQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetLbLocal(var2)))
         ||
         (! EQ(scip, SCIPvarGetUbLocal(var1), SCIPvarGetUbLocal(var2)))
      )
         return FALSE;
   }

   /* loop terminated, so columns are equal */
   return TRUE;
}

/** update the column order with a bound change */
static
SCIP_RETCODE updateColumnOrderWhenBranchingOnColumn(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< suborbitope data */
   int*                  roworder,           /**< array with the row order, of size nselrows */
   int                   nselrows,           /**< number of rows (required to be positive) */
   int*                  colorder,           /**< array to populate with column order, of size colorder */
   int*                  colorderinv,        /**< inverse array of the column order, of size colorder */
   SCIP_VAR*             var,                /**< variable that we branch on */
   COLSWAP*              thiscolswap         /**< the colswap to populate */
)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int origcolid;
   int swaporigcolid;
   int c;
   int i;
   int nrows;
   int ncols;
   int* origequalcolids;
   int norigequalcolids;
   int middlecolumn;
   int positionorigcolidincolorder;
   int positionswaporigcolidincolorder;

   assert( scip != NULL );
   assert( sodata != NULL );
   assert( roworder != NULL );
   assert( nselrows > 0 );
   assert( colorder != NULL );
   assert( sodata->ncols > 0 );
   assert( sodata->masterorbitopedata != NULL );

   ncols = sodata->ncols;
   assert( ncols > 0 );
   nrows = sodata->masterorbitopedata->nrows;
   assert( nrows > 0 );

   if ( sodata->columnordering == SCIP_COLUMNORDERING_NONE )
   {
      thiscolswap->from = 0;
      thiscolswap->to = 0;
      return SCIP_OKAY;
   }

   /* only variables from the orbitope matrix are of interest */
   if ( ! SCIPhashmapExists(sodata->colindexmap, (void*) var) )
   {
      thiscolswap->from = 0;
      thiscolswap->to = 0;
      return SCIP_OKAY;
   }

   origcolid = (int) (size_t) SCIPhashmapGetImage(sodata->colindexmap, (void*) var);
   assert( origcolid >= 0 );
   assert( origcolid < ncols );

   /* policy: swap with identical column that is closest to the center in relabeled order */
   /* @todo other policies: If the variable is in a ppc-row, then select the minimal/second to minimal to branch on */
   swaporigcolid = origcolid;

   /* we often use ncols / 2 */
   middlecolumn = ncols / 2;

   switch (sodata->columnordering)
   {
   case SCIP_COLUMNORDERING_FIRST:
   case SCIP_COLUMNORDERING_LAST:
   case SCIP_COLUMNORDERING_CENTRE:
      /* to each column, test column ordering condition, then if the column is identical to the var-column */
      for (c = 0; c < ncols; ++c)
      {
         /* origcolid is not interesting */
         if ( c == origcolid )
            continue;

         /* test if c is a better choice than swaporigcolid,
            * otherwise continue to next iteration through CONDITIONFAIL
            */
         switch (sodata->columnordering)
         {
            case SCIP_COLUMNORDERING_FIRST:
               /* only swap with c if c is earlier in column order than swaporigcolid */
               if ( colorder[c] >= colorder[swaporigcolid] )
                  goto CONDITIONFAIL;
               break;
            case SCIP_COLUMNORDERING_LAST:
               /* only swap with c if c is later in column order than swaporigcolid */
               if ( colorder[c] <= colorder[swaporigcolid] )
                  goto CONDITIONFAIL;
               break;
            case SCIP_COLUMNORDERING_CENTRE:
               /* if the column is not more central than swaporigcolid, ignore */
               if ( ABS(colorder[c] - middlecolumn) >=
                     ABS(colorder[swaporigcolid] - middlecolumn) )
                  goto CONDITIONFAIL;
               break;
            default:
               return SCIP_ERROR;
         }

         /* test: are c and origcolid the same columns w.r.t. the variable domain restrictions? */
         if ( !testColumnsAreSymmetricalyEquivalent(scip, sodata, c, origcolid) )
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
      SCIPallocBufferArray(scip, &origequalcolids, ncols);

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
         if ( !testColumnsAreSymmetricalyEquivalent(scip, sodata, c, origcolid) )
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
      var1 = sodata->sovars[i * ncols + swaporigcolid];
      var2 = sodata->sovars[i * ncols + origcolid];
      assert( EQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetLbLocal(var2)) );
      assert( EQ(scip, SCIPvarGetUbLocal(var1), SCIPvarGetUbLocal(var2)) );
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

   SCIPdebugMessage("Swapping %d (at %d) with %d (at %d)\n", origcolid, positionorigcolidincolorder,
      swaporigcolid, positionswaporigcolidincolorder);

   /* swap them, also keep track of the inverses */
   colorder[positionswaporigcolidincolorder] = origcolid;
   colorder[positionorigcolidincolorder] = swaporigcolid;
   colorderinv[origcolid] = positionswaporigcolidincolorder;
   colorderinv[swaporigcolid] = positionorigcolidincolorder;

   return SCIP_OKAY;
}


/** get the row order at the node
 *
 * The row order is given in the order of the variables that is branched on.
 * @todo combine with variant of cons_orbitope.c
 */
static
SCIP_RETCODE getRowOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< data for suborbitope */
   SCIP_NODE*            node,               /**< node to which the row order should be detected */
   int*                  roworder,           /**< array to populate with row order */
   int*                  nselrows            /**< pointer to number of rows part of the row order */
)
{
   int i;
   int j;
   BNBNODEINFO* ancestornodeinfo;
   BNBNODEINFO tmpnodeinfo;  /* used for lookups in hash table */

   assert( sodata != NULL );
   assert( sodata->masterorbitopedata != NULL );

   assert( sodata->masterorbitopedata->nrows > 0 );
   assert( sodata->ncols > 0 );

   *nselrows = 0;

   /* get the present row order up to this node (excluding the node itself) */
   assert( node != NULL );
   node = SCIPnodeGetParent(node);
   while (node != NULL)
   {
      /* retrieve the nodeinfo of this ancestor node */
      tmpnodeinfo.nodenumber = SCIPnodeGetNumber(node);
      ancestornodeinfo = (BNBNODEINFO*) SCIPhashtableRetrieve(sodata->nodeinfos, (void*) &tmpnodeinfo);
      if ( ancestornodeinfo != NULL )
      {
         assert( ancestornodeinfo->nrows >= 0 );
         for (i = ancestornodeinfo->nrows - 1; i >= 0; --i)
         {
            roworder[(*nselrows)++] = ancestornodeinfo->rows[i];
#ifndef NDEBUG
            {
               /* check if this row is not featured earlier */
               for (j = 0; j < (*nselrows) - 1; ++j)
               {
                  assert( ancestornodeinfo->rows[i] != roworder[j] );
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
      j = roworder[i];
      roworder[i] = roworder[(*nselrows) - 1 - i];
      roworder[(*nselrows) - 1 - i] = j;
   }

   return SCIP_OKAY;
}


/** gets rooted path up to node and populates column ordering array */
static
SCIP_RETCODE populateRootedPathColumnOrder(
   SUBORBITOPEDATA*      sodata,             /**< suborbitope data*/
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

   assert( node != NULL );

   depth = SCIPnodeGetDepth(node);
   i = depth;
   while ( node )
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
      ancestornodeinfo = (BNBNODEINFO*) SCIPhashtableRetrieve(sodata->nodeinfos, (void*) &tmpnodeinfo);
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
         assert( thiscolswap->from >= 0 && thiscolswap->from < sodata->ncols );
         assert( thiscolswap->to >= 0 && thiscolswap->to < sodata->ncols );

         /* at which column is origcolid? */
         positionfromincolorder = colorderinv[thiscolswap->from];
         assert( positionfromincolorder >= 0 );
         assert( positionfromincolorder < sodata->ncols );
         assert( colorder[positionfromincolorder] == thiscolswap->from );

         /* at which column is swaporigcolid? */
         positiontoincolorder = colorderinv[thiscolswap->to];
         assert( positiontoincolorder >= 0 );
         assert( positiontoincolorder < sodata->ncols );
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

/** at branching decisions, maintain the column swap and potential new rows in the orbitope */
static
SCIP_DECL_EVENTEXEC(eventExecNodeBranched)
{
   ORBITOPEDATA* consdata;
   SUBORBITOPEDATA* sodata;
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
   int s;
   int nrows;
   int rowid;
   /* data for row ordering */
   BNBNODEINFO* newnodeinfo;
   int* roworder;
   int nselrows;
   /* data for column ordering */
   int* colorder;
   int* colorderinv;
   SCIP_NODE** rootedpath;
   COLSWAP* thiscolswap;
   COLSWAP tmpcolswap;

   assert( eventdata != NULL );
   assert( !SCIPinProbing(scip) );

   eventnode = SCIPeventGetNode(event);
   assert( SCIPgetFocusNode(scip) == eventnode );

   consdata = (ORBITOPEDATA*) eventdata;
   assert( consdata != NULL );

   SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );

   /* arrays used within the loop */
   maxnbranchvars = 1;  /* it's a good guess that there's one branching variable, because that's likely the number */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, maxnbranchvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rootedpath, SCIPnodeGetDepth(eventnode)) );

   /* an orbitope could be split in multiple sub-orbitopes; handle each suborbitope separately */
   for (s = 0; s < consdata->nsuborbitopes; ++s)
   {
      sodata = &(consdata->suborbitopes[s]);
      assert( sodata != NULL);
      assert( sodata->masterorbitopedata != NULL );
      assert( sodata->masterorbitopedata == consdata );
      nrows = consdata->nrows;

      /* get all variables branched upon (check both branches) */
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
            if ( ! SCIPhashmapExists(sodata->rowindexmap, (void*) var) )
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

      /* skip suborbitopes whose variable matrix do not contain any branching variable */
      if ( nbranchvars <= 0 )
         continue;

      SCIP_CALL( SCIPallocBlockMemory(scip, &newnodeinfo) );
      newnodeinfo->nodenumber = SCIPnodeGetNumber(eventnode);
      newnodeinfo->colswaps = NULL;
      newnodeinfo->ncolswaps = 0;
      newnodeinfo->rows = NULL;
      newnodeinfo->nrows = 0;

      /* store data about row ordering */
      assert( nrows > 0 );
      SCIP_CALL( SCIPallocBufferArray(scip, &roworder, nrows) );
      nselrows = 0;

      /* get the present row order up to this node */
      SCIP_CALL( getRowOrder(scip, sodata, eventnode, roworder, &nselrows) );

      /* extend the row fixings with the steps from this node */
      for (i = 0; i < nbranchvars; ++i)
      {
         var = branchvars[i];

         assert( SCIPhashmapExists(sodata->rowindexmap, (void*) var) ); /* otherwise was not added to branchvars */
         rowid = (int) (size_t) SCIPhashmapGetImage(sodata->rowindexmap, (void*) var);
         assert( rowid >= 0 );
         assert( rowid < sodata->masterorbitopedata->nrows );

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
            SCIPallocBlockMemoryArray(scip, &newnodeinfo->rows, newnodeinfo->nrows + 1);
         }
         else
         {
            /* reallocate with linear increments, because we expect we expect 1 branching variable most of the time */
            SCIPreallocBlockMemoryArray(scip, &newnodeinfo->rows, newnodeinfo->nrows, newnodeinfo->nrows + 1);
         }
         newnodeinfo->rows[newnodeinfo->nrows++] = rowid;
      }

      /* store data about column ordering */
      assert( consdata->ncols > 0 );
      SCIP_CALL( SCIPallocBufferArray(scip, &colorder, consdata->ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colorderinv, consdata->ncols) );

      /* populate colorder with standard ordering */
      for (i = 0; i < consdata->ncols; ++i)
         colorder[i] = i;

      /* introduce inverse column ordering */
      for (i = 0; i < consdata->ncols; ++i)
         colorderinv[i] = i;

      /* get the rooted path
      *
      * We want to iterate through the bound changes in the order of the rooted path to this node.
      */
      node = SCIPnodeGetParent(eventnode);
      if ( node != NULL )
      {
         SCIP_CALL( populateRootedPathColumnOrder(sodata, node, rootedpath, colorder, colorderinv) );
      }

      /* get the swap for this node */
      for (i = 0; i < nbranchvars; ++i)
      {
         SCIP_CALL( updateColumnOrderWhenBranchingOnColumn(scip, sodata, roworder, nselrows, colorder,
            colorderinv, branchvars[i], &tmpcolswap) );
         /* skip trivial swaps of columns */
         if ( tmpcolswap.from == tmpcolswap.to )
            continue;

         /* mark that this row index is the new one in the node */
         if ( newnodeinfo->rows == NULL )
         {
            assert( newnodeinfo->nrows == 0 );
            SCIPallocBlockMemoryArray(scip, &newnodeinfo->colswaps, newnodeinfo->ncolswaps + 1);
         }
         else
         {
            /* reallocate with linear increments, because we expect we expect 1 branching variable most of the time */
            SCIPreallocBlockMemoryArray(scip, &newnodeinfo->colswaps, newnodeinfo->ncolswaps, newnodeinfo->ncolswaps + 1);
         }
         thiscolswap = &(newnodeinfo->colswaps[newnodeinfo->ncolswaps++]);
         thiscolswap->from = tmpcolswap.from;
         thiscolswap->to = tmpcolswap.to;
      }

      SCIPfreeBufferArray(scip, &colorderinv);
      SCIPfreeBufferArray(scip, &colorder);

      SCIPfreeBufferArray(scip, &roworder);

      /* store updates of row/column order or free memory if no change applied */
      if ( newnodeinfo->nrows > 0 || newnodeinfo->ncolswaps > 0 )
      {
         SCIP_CALL( SCIPhashtableSafeInsert(sodata->nodeinfos, newnodeinfo) );
      }
      else
      {
         SCIPfreeBlockMemory(scip, &newnodeinfo);
      }
   }

   SCIPfreeBufferArray(scip, &rootedpath);
   SCIPfreeBufferArray(scip, &branchvars);

   return SCIP_OKAY;
}


/** at branching decisions, maintain the column swap and potential new rows in the orbitope */
static
SCIP_DECL_EVENTEXEC(eventExec)
{
   switch (SCIPeventGetType(event))
   {
   case SCIP_EVENTTYPE_NODEBRANCHED:
      return eventExecNodeBranched(scip, eventhdlr, event, eventdata);
   default:
      return SCIP_ERROR;
   }
}


/** returns whether a row contains potential branching variables */
static
SCIP_Bool rowIsBranchRow(
   SCIP_ORBITOPALFIXINGDATA* orbifixdata,    /**< pointer to the dynamic orbitopal fixing data */
   ORBITOPEDATA*         orbidata,           /**< dynamic orbitope constraint data */
   int                   rowid               /**< row id for which to check */
)
{
   SCIP_VAR* var;
#ifndef NDEBUG
   int c;
#endif

   assert( orbidata != NULL );
   assert( orbidata->nrows > 0 );
   assert( orbidata->ncols > 0 );
   assert( rowid >= 0 );
   assert( rowid < orbidata->nrows );
   assert( orbidata->vars != NULL );
   assert( orbidata->vars[rowid * orbidata->ncols + 0] );

   /* get the first variable from the row */
   var = orbidata->vars[rowid * orbidata->ncols];

   /* debugging: the variable types in a row should all be the same */
#ifndef NDEBUG
   for (c = 1; c < orbidata->ncols; ++c)
   {
      /* the actual vartypes can be different,
       * for example when an INTEGER vartype turns into BINARY due to bound changes
       */
      assert( vartypeIsBranchRowType(orbifixdata, SCIPvarGetType(var)) ==
         vartypeIsBranchRowType(orbifixdata, SCIPvarGetType(orbidata->vars[rowid * orbidata->ncols + c])) );
   }
#endif

   return vartypeIsBranchRowType(orbifixdata, SCIPvarGetType(var));
}


/** clear suborbitope data */
static
SCIP_RETCODE clearSuborbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata              /**< pointer to suborbitope data to clear */
)
{
   BNBNODEINFO* nodeinfo;
   int nentries;
   int i;

   assert( sodata != NULL );
   assert( sodata->masterorbitopedata != NULL );

   nentries = SCIPhashtableGetNEntries(sodata->nodeinfos);
   for (i = 0; i < nentries; ++i)
   {
      /* @todo in principle, can deal with memory sparsity by first getting all nodeinfos,
      * then sorting by address and free them in descending order
      */
      nodeinfo = (BNBNODEINFO*) (SCIPhashtableGetEntry(sodata->nodeinfos, i));
      if ( nodeinfo == NULL )
         continue;

      assert( nodeinfo != NULL );
      assert( nodeinfo->nrows > 0 || nodeinfo->ncolswaps > 0 );

      assert( (nodeinfo->ncolswaps == 0) != (nodeinfo->colswaps != NULL) );
      if ( nodeinfo->colswaps != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &(nodeinfo->colswaps), nodeinfo->ncolswaps);
      }

      assert( (nodeinfo->nrows == 0) != (nodeinfo->rows != NULL) );
      if ( nodeinfo->rows != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &(nodeinfo->rows), nodeinfo->nrows);
      }

      SCIPfreeBlockMemory(scip, &nodeinfo);
   }
   SCIPhashtableFree(&(sodata->nodeinfos));

   SCIPhashmapFree(&(sodata->colindexmap));
   SCIPhashmapFree(&(sodata->rowindexmap));

   SCIPfreeBlockMemoryArrayNull(scip, &(sodata->sovars), sodata->masterorbitopedata->nrows * sodata->ncols);

   return SCIP_OKAY;
}


/** frees orbitope constraint data */
static
SCIP_RETCODE freeOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA* orbifixdata,    /**< pointer to the dynamic orbitopal fixing data */
   ORBITOPEDATA**        orbidata            /**< pointer to orbitope constraint data */
   )
{
   SUBORBITOPEDATA* sodata;
   int s;

   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( orbidata != NULL );
   assert( *orbidata != NULL );
   assert( (*orbidata)->vars != NULL );
   assert( (*orbidata)->nrows > 0 );
   assert( (*orbidata)->ncols > 0 );
   assert( (*orbidata)->nrows * (*orbidata)->ncols > 0 );
   assert( (*orbidata)->nsuborbitopes >= 0 );
   assert( (*orbidata)->suborbitopes != NULL || (*orbidata)->nsuborbitopes == 0 );
   assert( SCIPisTransformed(scip) );

   /* drop event */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED, orbifixdata->eventhdlr,
      (SCIP_EVENTDATA*) *orbidata, -1 ) );

   /* clear all the suborbitopes */
   if ( (*orbidata)->nsuborbitopes > 0 )
   {
      for (s = (*orbidata)->nsuborbitopes - 1; s >= 0; --s)
      {
         sodata = &((*orbidata)->suborbitopes[s]);
         assert( sodata != NULL );
         assert( sodata->masterorbitopedata == (*orbidata) );
         SCIP_CALL( clearSuborbitope(scip, sodata) );
      }

      SCIPfreeBlockMemoryArray(scip, &((*orbidata)->suborbitopes), (*orbidata)->nsuborbitopes);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*orbidata)->vars), (*orbidata)->nrows * (*orbidata)->ncols);

   SCIPfreeBlockMemory(scip, orbidata);

   return SCIP_OKAY;
}


/** initialize suborbitope data */
static
SCIP_RETCODE computeSuborbitopesAndInitialize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA*  orbifixdata,        /**< pointer to the dynamic orbitopal fixing data */
   ORBITOPEDATA*         orbidata            /**< pointer to constraint data */
)
{
   SUBORBITOPEDATA* sodata;
   SCIP_VAR* var;
   int i;
   int j;
   int c;
   int c_;
   int* columnsubset;
   int* thissubset;
   int nthissubset;
   int ncolumnsubset;
   int nidenticalcolumns;
   int nelem;
   int nrows;
   int ncols;
   int rowid;
   int colid;

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( orbidata->suborbitopes == NULL );
   assert( orbidata->nsuborbitopes == 0 );
   assert( orbidata->nrows > 0 );
   assert( orbidata->ncols > 0 );

   /* for each column, an index that indicates which class of columns it is part of. */
   SCIPallocBufferArray(scip, &columnsubset, orbidata->ncols);

   /* initialize all columns to be in their own class */
   for (c = 0; c < orbidata->ncols; ++c)
      columnsubset[c] = -1;

   /* determine which columns are equal and which are not */
   ncolumnsubset = 0;

   for (c = 0; c < orbidata->ncols; ++c)
   {
      /* skip if already in a class */
      if ( columnsubset[c] >= 0 )
         continue;

      /* to which columns is c equal? */
      columnsubset[c] = ncolumnsubset;
      nidenticalcolumns = 0;
      for (c_ = c + 1; c_ < orbidata->ncols; ++c_)
      {
         /* skip columns already contained in a different subset
          * if c_ is already put in a subset of columns, then it's not the same as c because we would've added c, then
          */
         if ( columnsubset[c_] >= 0 )
            continue;

         /* are c_ and c identical? */
         for (i = 0; i < orbidata->nrows; ++i)
         {
            if ( ! EQ(scip,
               SCIPvarGetLbLocal(orbidata->vars[i * orbidata->ncols + c]),
               SCIPvarGetLbLocal(orbidata->vars[i * orbidata->ncols + c_])
            ))
               break;
            if ( ! EQ(scip,
               SCIPvarGetUbLocal(orbidata->vars[i * orbidata->ncols + c]),
               SCIPvarGetUbLocal(orbidata->vars[i * orbidata->ncols + c_])
            ))
               break;
         }
         if ( i < orbidata->nrows )  /* loop is stopped by break */
            continue;

         /* columns c and c_ are identical */
         columnsubset[c_] = ncolumnsubset;
         ++nidenticalcolumns;
      }

      /* if c is alone, columnsubset to -1; otherwise, prepare for the next class */
      if ( nidenticalcolumns == 0 )
         columnsubset[c] = -1;
      else
         ++ncolumnsubset;
   }

   orbidata->nsuborbitopes = ncolumnsubset;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &orbidata->suborbitopes, orbidata->nsuborbitopes) );

   SCIP_CALL( SCIPallocBufferArray(scip, &thissubset, orbidata->ncols) );

   for (i = 0; i < orbidata->nsuborbitopes; ++i)
   {
      /* determine all columns that are labeled with i  */
      nthissubset = 0;
      for (c = 0; c < orbidata->ncols; ++c)
      {
         if ( columnsubset[c] == i )
            thissubset[nthissubset++] = c;
      }
      assert( nthissubset > 1 );

      /* get the suborbitope data that we are going to populate */
      sodata = &(orbidata->suborbitopes[i]);

      nrows = orbidata->nrows;
      ncols = nthissubset;
      nelem = nrows * ncols;

      sodata->masterorbitopedata = orbidata;
      sodata->ncols = ncols;
#ifndef NDEBUG
      sodata->lastnodenumber = -1;
      sodata->dbghash = -1;
#endif

      /* @todo allow more dynamic changes per suborbitope */
      sodata->columnordering = orbifixdata->defaultcolumnordering;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sodata->sovars, nelem) );
      SCIP_CALL( SCIPhashmapCreate(&sodata->rowindexmap, SCIPblkmem(scip), nrows) );
      SCIP_CALL( SCIPhashmapCreate(&sodata->colindexmap, SCIPblkmem(scip), ncols) );

      /* nodeinfos: every node number is mapped to the struct
       *
       * Hashtables of initial size 1 seem to cover most cases as there is usually just one branching variable.
       */
      SCIP_CALL( SCIPhashtableCreate(&sodata->nodeinfos, scip->mem->probmem, 1,
         hashGetKeyBnbnodeinfo, hashKeyEqBnbnodeinfo, hashKeyValBnbnodeinfo, NULL) );

      SCIPdebugMessage("Orbitope variables for (%dx%d) suborbitope for consdata %p and subortibtope %d\n",
         nrows, ncols, (void*) orbidata, i);
      for (j = 0, rowid = 0, colid = 0; j < nelem; ++j, ++colid)
      {
         if ( colid == ncols )
         {
            colid = 0;
            ++rowid;
         }
         assert( rowid == j / ncols );
         assert( colid == j % ncols );

         var = orbidata->vars[rowid * orbidata->ncols + thissubset[colid]];
         assert( var != NULL );

         /* variables cannot be repeated in the variable matrix */
         assert( ! SCIPhashmapExists(sodata->rowindexmap, var) );
         SCIP_CALL( SCIPhashmapInsert(sodata->rowindexmap, var, (void*) (size_t) rowid) );

         assert( ! SCIPhashmapExists(sodata->colindexmap, var) );
         SCIP_CALL( SCIPhashmapInsert(sodata->colindexmap, var, (void*) (size_t) colid) );

         sodata->sovars[j] = var;

         SCIPdebugMessage("%4d %4d -> %s\n", rowid, colid, var->name);
      }
   }

   SCIPfreeBufferArray(scip, &thissubset);
   SCIPfreeBufferArray(scip, &columnsubset);

   return SCIP_OKAY;
}


/** adds an orbitope to the orbital fixing data */
static
SCIP_RETCODE addOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA* orbifixdata,    /**< pointer to the dynamic orbitopal fixing data */
   SCIP_VAR**            vars,               /**< variables array, must have size nrows * ncols */
   int                   nrows,              /**< number of rows in orbitope */
   int                   ncols               /**< number of columns in orbitope */
   )
{
   ORBITOPEDATA* orbidata;
   SCIP_VAR* var;
   int i;
   int rowid;
   int colid;
   int nelem;

   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( orbifixdata->eventhdlr != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );

   nelem = nrows * ncols;
   assert( nelem > 0 );

   SCIP_CALL( SCIPallocBlockMemory(scip, &orbidata) );

   orbidata->nrows = nrows;
   orbidata->ncols = ncols;
   orbidata->nbranchrows = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &orbidata->vars, nelem) );

   SCIPdebugMessage("Orbitope variables for (%dx%d) orbitope with consdata %p\n", nrows, ncols, (void*) orbidata);
   for (i = 0, rowid = 0, colid = 0; i < nelem; ++i, ++colid)
   {
      if ( colid == ncols )
      {
         colid = 0;
         ++rowid;
      }
      assert( rowid == i / ncols );
      assert( colid == i % ncols );

      var = vars[i];
      assert( var != NULL );

      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var ) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );

      orbidata->vars[i] = var;

      SCIPdebugMessage("%4d %4d -> %s\n", rowid, colid, var->name);
   }

   orbidata->suborbitopes = NULL;
   orbidata->nsuborbitopes = 0;

   /* @todo at getRowData: If nselrows == nbranchrows, append the non-branch rows (like before) */
   for (i = 0; i < nrows; ++i)
   {
      if ( rowIsBranchRow(orbifixdata, orbidata, i) )
         ++orbidata->nbranchrows;
   }

   /* constraint handler data */
   assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING ? SCIPgetNNodes(scip) == 0 : TRUE );

   /* add the event to store the row and column updates of nodes in the branch-and-bound tree */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED, orbifixdata->eventhdlr, 
      (SCIP_EVENTDATA*) orbidata, NULL) );

   /* construct all suborbitope stuff */
   SCIP_CALL( computeSuborbitopesAndInitialize(scip, orbifixdata, orbidata) );

   /* resize orbitope array if needed */
   assert( orbifixdata->norbitopes >= 0 );
   assert( (orbifixdata->norbitopes == 0) == (orbifixdata->orbitopes == NULL) );
   assert( orbifixdata->norbitopes <= orbifixdata->maxnorbitopes );
   if ( orbifixdata->norbitopes == orbifixdata->maxnorbitopes )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, orbifixdata->norbitopes + 1);
      assert( newsize >= 0 );

      if ( orbifixdata->norbitopes == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &orbifixdata->orbitopes, newsize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &orbifixdata->orbitopes, orbifixdata->norbitopes, newsize) );
      }

      orbifixdata->maxnorbitopes = newsize;
   }

   /* add orbitope to orbitopal fixing data */
   assert( orbifixdata->norbitopes < orbifixdata->maxnorbitopes );
   orbifixdata->orbitopes[orbifixdata->norbitopes++] = orbidata;

   SCIPdebugMsg(scip, "Added orbitope for orbitopal fixing of size %d by %d\n", nrows, ncols);

   return SCIP_OKAY;
}

/** get the column order at the node
 *
 * The column order is (deterministically) dynamically decided based on the policy for column ordering
 */
static
SCIP_RETCODE getColumnOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< suborbitope data */
   SCIP_NODE*            eventnode,          /**< node where this should be determined at */
   int*                  roworder,           /**< array with the row order, of size nselrows */
   int                   nselrows,           /**< number of rows (required to be positive) */
   int*                  colorder,           /**< array to populate with column order, of size ncols */
   int*                  colorderinv         /**< array to populate with inverse column order, of size ncols */
)
{
   SCIP_NODE* node;
   SCIP_NODE** rootedpath;
   int i;
   int depth;
   int ncols;

   assert( scip != NULL );
   assert( sodata != NULL );
   assert( nselrows > 0 );
   ncols = sodata->ncols;

   /* populate colorder with standard ordering */
   for (i = 0; i < ncols; ++i)
      colorder[i] = i;

   /* introduce inverse column ordering */
   for (i = 0; i < ncols; ++i)
      colorderinv[i] = i;

   /* get the rooted path
    *
    * We want to iterate through the bound changes in the order of the rooted path to this node.
    */
   node = SCIPnodeGetParent(eventnode);
   if ( node != NULL )
   {
      depth = SCIPnodeGetDepth(node);
      SCIP_CALL( SCIPallocBufferArray(scip, &rootedpath, depth + 1) );
      SCIP_CALL( populateRootedPathColumnOrder(sodata, node, rootedpath, colorder, colorderinv) );
      SCIPfreeBufferArray(scip, &rootedpath);
   }

   return SCIP_OKAY;
}


#ifndef NDEBUG
/** checks if a matrix is contained in an orbitope, i.e., columns are sorted lexicographically non-increasingly */
static
void assertIsOrbitopeMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< constraint data */
   int*                  roworder,           /**< array with the row order */
   int*                  colorder,           /**< array with the column order */
   SCIP_Real*            matrix,             /**< a matrix */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of cols of matrix */
   int*                  infinitesimal,      /**< array specifying where the infinitesimals are at */
   SCIP_Bool             addinfinitesimals   /**< whether infinitesimals are added or subtracted */
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
   assert( roworder != NULL );
   assert( colorder != NULL );
   assert( matrix != NULL );
   assert( sodata != NULL );
   assert( sodata->masterorbitopedata != NULL );
   assert( sodata->sovars != NULL );
   assert( sodata->ncols > 0 );

   /* respect variable bounds */
   for (rowid = 0; rowid < nrows; ++rowid)
   {
      origrowid = roworder[rowid];
      for (colid = 0; colid < ncols; ++colid)
      {
         origcolid = colorder[colid];
         idx = rowid * ncols + colid;
         origidx = origrowid * ncols + origcolid;
         var = sodata->sovars[origidx];
         assert( GE(scip, matrix[idx], SCIPvarGetLbLocal(var)) );
         assert( LE(scip, matrix[idx], SCIPvarGetUbLocal(var)) );
      }
   }

   /* is orbitope */
   for (colid = 0; colid < ncols - 1; ++colid)
   {
      /* compare column colid with colid + 1 */
      for (rowid = 0; rowid < nrows; ++rowid)
      {
         /* entry is >= entry to the right */
         assert( GE(scip, matrix[rowid * ncols + colid], matrix[rowid * ncols + colid + 1]) );

         if ( GT(scip, matrix[rowid * ncols + colid], matrix[rowid * ncols + colid + 1]) )
         {
            /* critical row */
            break;
         }
         else
         {
            /* check for infinitisimal values
             * If infinitesimals are added (lexminface case), then if the left column has a +epsilon, 
             * it does not matter whether the right column has +epsilon or not, then the left column is >, 
             * due to the axioms x + epsilon > x + epsilon and x + epsilon > x.
             * Analogously, x > x - epsilon and x - epsilon > x - epsilon. 
             */
            assert( EQ(scip, matrix[rowid * ncols + colid], matrix[rowid * ncols + colid + 1]) );
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
/** to test if arrays are the same, generate some hash for an array of integers */
static
int debugGetArrayHash(
   int*                  array,              /** array */
   int                   len                 /** array length */
)
{
   int i;
   unsigned int hash = 0;

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
   SCIP_Real*            matrix,             /** Matrix */
   int                   nrows,              /** Number of rows */
   int                   ncols               /** Number of rows */
)
{
   int row;
   int col;

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


/** get the column order at the node */
static
SCIP_RETCODE propagateStaticOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< suborbitope data */
   int*                  roworder,           /**< array with the row order */
   int                   nselrows,           /**< number of selected rows */
   int*                  colorder,           /**< array with the column order */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nfixedvars          /**< pointer to counter of number of variable domain reductions */
)
{
   /* @todo also make "nselcols" to allow for colorders smaller than consdata->ncols */
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
   assert( sodata != NULL );
   assert( sodata->sovars != NULL );
   assert( roworder != NULL );
   assert( nselrows >= 0 );
   assert( colorder != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   *infeasible = FALSE;

   assert( sodata->masterorbitopedata != NULL );
   assert( sodata->masterorbitopedata->nrows > 0 );
   assert( sodata->masterorbitopedata->nrows >= nselrows );
   ncols = sodata->ncols;
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
            thisvar = sodata->sovars[roworder[r] * ncols + colorder[k]];
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
   origcolid = colorder[colid];
   for (rowid = 0; rowid < nselrows; ++rowid)
   {
      origrowid = roworder[rowid];
      origidx = origrowid * ncols + origcolid;
      var = sodata->sovars[origidx];
      i = rowid * ncols + colid;

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

      origcolid = colorder[colid];
      for (rowid = 0; rowid < nselrows; ++rowid)
      {
         origrowid = roworder[rowid];
         origidx = origrowid * ncols + origcolid;
         i = rowid * ncols + colid;

         /* the entry one to the right is not the first column */
         assert( (i + 1) % ncols > 0 );

         var = sodata->sovars[origidx];
         assert( var != NULL );

         if ( iseq )
         {
            /* equality holds up to this row
             * Compare to the entry value on the column immediately right.
             * The value we choose on the left must be at least this.
             * 2 Options:
             * Option 1: The upper bound is smaller. Then we're in an infeasible situation. Fix as described below.
             * Option 2: The upper bound is greater or equal.
             */
            ub = SCIPvarGetUbLocal(var);

            /* compare to the value in the column right of it */
            if ( LT(scip, ub, lexminface[i + 1]) || 
               ( lexminepsrow[colid + 1] == rowid && EQ(scip, ub, lexminface[i + 1]) ) )
            {
               /* value of this column can only be strictly smaller than the value in the column to its right
                * This may not be possible.
                * Try to repair: Go back to the last row with "room" left, and make the value minimally larger.
                */
               if ( lastunfixed >= 0 )
               {
                  /* repair: return to the last row with "room", and increase the lexmin-value at that row. */
                  assert( EQ(scip, lexminface[lastunfixed * ncols + colid],
                     lexminface[lastunfixed * ncols + colid + 1]) );
                  othervar = sodata->sovars[roworder[lastunfixed] * ncols + origcolid];
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
                     assert( LE(scip, lexminface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     break;
                  case SCIP_VARTYPE_CONTINUOUS:
                     /* continuous type, so add an infinitisimal value to the bound */
                     assert( LE(scip, lexminface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     assert( lexminepsrow[colid] == -1 );
                     lexminepsrow[colid] = lastunfixed;
                     break;
                  default:
                     return SCIP_ERROR;
                  }
                  /* now row "lastunfixed" is greater. Restart from here. */
                  iseq = FALSE;
                  rowid = lastunfixed; /* the next iteration considers "lastunfixed + 1" */
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
               assert( GE(scip, ub, lexminface[i + 1]) );
               lb = SCIPvarGetLbLocal(var);
               assert( LE(scip, lb, ub) );
               lexminface[i] = MAX(lexminface[i + 1], lb);
               assert( GE(scip, lexminface[i], lexminface[i + 1]) );

               /* are we still equal? */
               if ( GT(scip, lexminface[i], lexminface[i + 1]) )
                  iseq = FALSE;
               else if ( lexminepsrow[colid + 1] == rowid )
               {
                  assert( EQ(scip, lexminface[i], lexminface[i + 1]) );
                  assert( SCIPvarGetType(sodata->sovars[roworder[lastunfixed] * ncols + origcolid]) 
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
                  if ( LE(scip, lexminface[i] + 1.0, ub) )
                     lastunfixed = rowid;
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous type: if we can add an infinitesimal value to the current lexminface[i] value, 
                   * mark row as 'lastunfixed'
                   */
                  if ( LT(scip, lexminface[i], ub) )
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
   assertIsOrbitopeMatrix(scip, sodata, roworder, colorder, lexminface, nselrows, ncols, lexminepsrow, TRUE);
#endif

   /* compute lexmax face */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxface, nelem) );

   /* array to store for each column at which row we subtract an infinitesimal value, initially at none (-1) */
   SCIP_CALL( SCIPallocBufferArray(scip, &lexmaxepsrow, ncols) );
   for (colid = 0; colid < ncols; ++colid)
      lexmaxepsrow[colid] = -1;

   /* first column, fill all unfixed entries with maximally possible values */
   colid = 0;
   origcolid = colorder[colid];
   for (rowid = 0; rowid < nselrows; ++rowid)
   {
      origrowid = roworder[rowid];
      origidx = origrowid * ncols + origcolid;
      var = sodata->sovars[origidx];
      i = rowid * ncols + colid;

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

      origcolid = colorder[colid];
      for (rowid = 0; rowid < nselrows; ++rowid)
      {
         origrowid = roworder[rowid];
         origidx = origrowid * ncols + origcolid;
         i = rowid * ncols + colid;

         /* the entry one to the left is not the last column */
         assert( i % ncols > 0 );

         var = sodata->sovars[origidx];
         assert( var != NULL );

         if ( iseq )
         {
            /* equality holds up to this row
             * Compare to the entry value on the column immediately left.
             * The value we choose on the right must be at most this.
             * 2 Options:
             * Option 1: The lower bound is larger. Then we're in an infeasible situation. Fix as described below.
             * Option 2: The lower bound is smaller or equal.
             */
            lb = SCIPvarGetLbLocal(var);

            /* compare to the value in the column left of it */
            if ( GT(scip, lb, lexmaxface[i - 1]) ||
               ( lexmaxepsrow[colid - 1] == rowid && EQ(scip, lb, lexmaxface[i - 1]) ) )
            {
               /* value of this column can only be strictly larger than the value in the column to its left
                * This may not be possible.
                * Try to repair: Go back to the last row with "room" left, and make the value minimally smaller.
                */
               if ( lastunfixed >= 0 )
               {
                  /* repair: return to the last row with "room", and increase the lexmax-value at that row. */
                  assert( EQ(scip, lexmaxface[lastunfixed * ncols + colid],
                     lexmaxface[lastunfixed * ncols + colid - 1]) );
                  othervar = sodata->sovars[roworder[lastunfixed] * ncols + origcolid];
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
                     assert( GE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetLbLocal(othervar)) );
                     assert( LE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     break;
                  case SCIP_VARTYPE_CONTINUOUS:
                     /* continuous type, so subtract an infinitisimal value to the bound */
                     assert( GE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetLbLocal(othervar)) );
                     assert( LE(scip, lexmaxface[lastunfixed * ncols + colid], SCIPvarGetUbLocal(othervar)) );
                     assert( lexmaxepsrow[colid] == -1 );
                     lexmaxepsrow[colid] = lastunfixed;
                     break;
                  default:
                     return SCIP_ERROR;
                  }
                  /* now row "lastunfixed" is greater. Restart from here. */
                  iseq = FALSE;
                  rowid = lastunfixed; /* the next iteration considers "lastunfixed + 1" */
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
               assert( LE(scip, lb, lexmaxface[i - 1]) );
               ub = SCIPvarGetUbLocal(var);
               assert( LE(scip, lb, ub) );
               lexmaxface[i] = MIN(lexmaxface[i - 1], ub);
               assert( GE(scip, lexmaxface[i - 1], lexmaxface[i]) );

               /* are we still equal? */
               if ( GT(scip, lexmaxface[i - 1], lexmaxface[i]) )
                  iseq = FALSE;
               else if ( lexmaxepsrow[colid - 1] == rowid )
               {
                  assert( EQ(scip, lexmaxface[i - 1], lexmaxface[i]) );
                  assert( SCIPvarGetType(sodata->sovars[roworder[lastunfixed] * ncols + origcolid]) 
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
                  if ( GE(scip, lexmaxface[i] - 1.0, lb) )
                     lastunfixed = rowid;
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous type: if we can subtract an infinitesimal value to the current lexmaxface[i] value, 
                   * mark row as 'lastunfixed'
                   */
                  if ( GT(scip, lexmaxface[i], lb) )
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
   assertIsOrbitopeMatrix(scip, sodata, roworder, colorder, lexmaxface, nselrows, ncols, lexmaxepsrow, FALSE);
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
      for (rowid = 0; rowid < nselrows; ++rowid)
      {
         i = rowid * ncols + colid;

         /* get var */
         origrowid = roworder[rowid];
         origcolid = colorder[colid];
         origidx = origrowid * ncols + origcolid;
         var = sodata->sovars[origidx];

         if ( EQ(scip, lexminface[i], lexmaxface[i]) )
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


/** propagation method for a suborbitope */
static
SCIP_RETCODE propagateSuborbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SUBORBITOPEDATA*      sodata,             /**< suborbitope data */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nfixedvars          /**< pointer to store the number of found domain reductions */
   )
{
   SCIP_NODE* focusnode;
   int* roworder;
   int nselrows;
   int* colorder;
   int* colorderinv;
   int nrows;
   int ncols;

   assert( scip != NULL );
   assert( sodata != NULL );
   assert( sodata->masterorbitopedata != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   ncols = sodata->ncols;
   nrows = sodata->masterorbitopedata->nrows;
   assert( ncols > 0 );
   assert( nrows > 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &roworder, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colorder, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colorderinv, ncols) );

   focusnode = SCIPgetFocusNode(scip);
   assert( focusnode != NULL );

   SCIP_CALL( getRowOrder(scip, sodata, focusnode, roworder, &nselrows) );
   assert( nselrows >= 0 );
   assert( nselrows <= nrows );
   if ( nselrows == 0 )
      goto FREE;

   SCIP_CALL( getColumnOrder(scip, sodata, focusnode, roworder, nselrows, colorder, colorderinv) );

#ifndef NDEBUG
   /* DEBUG: if propagation is repeated in the same node, the same column order and row order is needed */
   /* @todo: performance: move roworder and colorder to consdata, then re-use */
   {
      int colhash = debugGetArrayHash(colorder, ncols);
      int rowhash = debugGetArrayHash(roworder, nselrows);
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
      if ( sodata->lastnodenumber == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
      {
         assert( sodata->dbghash == hash );
      }
      sodata->dbghash = hash;
   }
   sodata->lastnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
#endif

   SCIP_CALL( propagateStaticOrbitope(scip, sodata, roworder, nselrows, colorder, infeasible, nfixedvars) );

   FREE:
   SCIPfreeBufferArray(scip, &colorderinv);
   SCIPfreeBufferArray(scip, &colorder);
   SCIPfreeBufferArray(scip, &roworder);

#ifdef SCIP_MORE_DEBUG
   SCIPdebugPrintf("\n\n");
#endif

   return SCIP_OKAY;
}


/** propagation method for a single orbitope constraint */
static
SCIP_RETCODE propagateOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   ORBITOPEDATA*         orbidata,           /**< orbitope data*/
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nfixedvars          /**< pointer to store the number of found domain reductions */
   )
{
   SUBORBITOPEDATA* sodata;
   int s;

   assert( scip != NULL );
   assert( orbidata != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );
   assert( orbidata->nsuborbitopes > 0 );

   *nfixedvars = 0;
   *infeasible = FALSE;

   for (s = 0; s < orbidata->nsuborbitopes; ++s)
   {
      sodata = &(orbidata->suborbitopes[s]);
      assert( sodata->masterorbitopedata == orbidata );

      SCIP_CALL( propagateSuborbitope(scip, sodata, infeasible, nfixedvars) );
      if ( *infeasible )
         break;
   }

   return SCIP_OKAY;
}


/*
 * Interface methods
 */

/** propagates orbitopal fixing */
SCIP_RETCODE SCIPorbitopalFixingPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA*  orbifixdata,        /**< orbitopal fixing data structure */
   SCIP_Bool*            infeasible,         /**< whether infeasibility is found */
   int*                  nred,               /**< number of domain reductions */
   SCIP_Bool*            didrun              /**< whether propagator actually ran */
   )
{
   ORBITOPEDATA* orbidata;
   int c;
   int thisfixedvars;

   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( (orbifixdata->norbitopes == 0) == (orbifixdata->orbitopes == NULL) );
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

   /* propagate all useful constraints */
   for (c = 0; c < orbifixdata->norbitopes; ++c)
   {
      orbidata = orbifixdata->orbitopes[c];
      assert( orbidata != NULL );

      SCIP_CALL( propagateOrbitope(scip, orbidata, infeasible, &thisfixedvars) );
      SCIPdebugMessage("Found %d reductions during orbitopal fixing for orbitope %d\n", thisfixedvars, c);
      *nred += thisfixedvars;
      *didrun = TRUE;

      /* stop if we find infeasibility in one of the constraints */
      if ( *infeasible )
      {
         SCIPdebugMessage("Detected infeasibility during orbitopal fixing for orbitope %d\n", c);
         break;
      }
   }

   return SCIP_OKAY;
}

/** marks presence of orbitopal symmetries component for orbitopal fixing */
SCIP_RETCODE SCIPorbitopalFixingAddOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA*  orbifixdata,        /**< orbitopal fixing data structure */
   SCIP_VAR**            vars,               /**< matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows */
   int                   ncols               /**< number of columns */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );

   /* dynamic symmetry reductions cannot be performed on original problem */
   assert( SCIPisTransformed(scip) );

   /* if this is the first time adding an orbitope, check if the nonlinear conshlr exists */
   if ( !orbifixdata->conshdlr_nonlinear_checked )
   {
      orbifixdata->conshdlr_nonlinear = SCIPfindConshdlr(scip, "nonlinear");
      orbifixdata->conshdlr_nonlinear_checked = TRUE;
   }

   /* create constraint data */
   SCIP_CALL( addOrbitope(scip, orbifixdata, vars, nrows, ncols) );

   return SCIP_OKAY;
}


/** resets orbitopal fixing data structure (clears all orbitopes) */
SCIP_RETCODE SCIPorbitopalFixingReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA*  orbifixdata         /**< pointer to orbitopal fixing structure to populate */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( orbifixdata->norbitopes >= 0 );
   assert( (orbifixdata->norbitopes == 0) == (orbifixdata->orbitopes == NULL) );
   assert( orbifixdata->norbitopes <= orbifixdata->maxnorbitopes );
   assert( orbifixdata->eventhdlr != NULL );

   /* free orbitopes that are added */
   while (orbifixdata->norbitopes > 0)
   {
      SCIP_CALL( freeOrbitope(scip, orbifixdata, &(orbifixdata->orbitopes[--orbifixdata->norbitopes])) );
   }
   assert( orbifixdata->norbitopes == 0 );
   SCIPfreeBlockMemoryArray(scip, &orbifixdata->orbitopes, orbifixdata->maxnorbitopes);
   orbifixdata->orbitopes = NULL;
   orbifixdata->maxnorbitopes = 0;

   return SCIP_OKAY;
}


/** free orbitopal fixing data */
SCIP_RETCODE SCIPorbitopalFixingFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA** orbifixdata    /**< pointer to orbitopal fixing structure to populate */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( *orbifixdata != NULL );

   SCIP_CALL( SCIPorbitopalFixingReset(scip, *orbifixdata) );

   SCIPfreeBlockMemory(scip, orbifixdata);
   return SCIP_OKAY;
}


/** initializes structures needed for orbitopal fixing
 * This is only done exactly once.
 */
SCIP_RETCODE SCIPorbitopalFixingInclude(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALFIXINGDATA** orbifixdata    /**< pointer to orbitopal fixing structure to populate */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert( scip != NULL );
   assert( orbifixdata != NULL );

   assert( SCIPgetStage(scip) == SCIP_STAGE_INIT );

   /* create orbitope constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, orbifixdata) );

   /* default column ordering param */
   SCIP_CALL( SCIPaddIntParam(scip, "propagating/symmetry/orbitopalfixing/columnordering",
         "The column ordering variant, respects enum SCIP_OrbitopeColumnOrdering.",
         (int*) &(*orbifixdata)->defaultcolumnordering, TRUE, DEFAULT_COLUMNORDERING, 0, 4, NULL, NULL) );

   /* initialize event handler. */
   assert( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL );
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExec, NULL) );
   assert( eventhdlr != NULL );
   (*orbifixdata)->eventhdlr = eventhdlr;

   /* orbitopes array */
   (*orbifixdata)->orbitopes = NULL;
   (*orbifixdata)->norbitopes = 0;
   (*orbifixdata)->maxnorbitopes = 0;

   /* conshdlr nonlinear */
   (*orbifixdata)->conshdlr_nonlinear = NULL;
   (*orbifixdata)->conshdlr_nonlinear_checked = FALSE;

   return SCIP_OKAY;
}
