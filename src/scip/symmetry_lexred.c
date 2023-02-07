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

/**@file   symmetry_lexred.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries by dynamic lexicographic ordering reduction
 * @author Jasper van Doornmalen
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/symmetry_lexred.h"
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
#include "scip/event_shadowtree.h"
#include <ctype.h>
#include <string.h>
#include <symmetry/type_symmetry.h>

#include <memory.h>


/*
 * Data structures
 */


/** data per permutation for lexicographic reduction propagator */
struct LexicographicReductionPermutationData
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   int*                  perm;               /**< permutation associated to the symresack */
   int*                  invperm;            /**< inverse permutation */
};
typedef struct LexicographicReductionPermutationData LEXDATA;


/** data for dynamic lexicographic reduction propagator */
struct SCIP_LexicographicReductionData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr;/**< eventhandler for the shadow tree data structure */
   SCIP_HASHMAP*         symvarmap;          /**< map of variables to indices in permvars array */
   int                   nsymvars;           /**< number of variables in symvarmap */
   LEXDATA**             lexdatas;           /**< array of pointers to individual LEXDATA's */
   int                   nlexdatas;          /**< number of datas in array */
   int                   maxnlexdatas;       /**< allocated datas array size */
};


/** to store branch-and-bound tree paths, (depth, index)-information per variable in rooted path */
struct nodedepthbranchindex
{
   int nodedepth;                            /**< depth of var domain change */
   int index;                                /**< index of var domain change on node at depth */
};
typedef struct nodedepthbranchindex NODEDEPTHBRANCHINDEX;


/** auxiliary struct to pass branch-and-bound tree information to sort function */
struct vararraynodedepthbranchindex
{
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices; /**< pointer to branch-and-bound tree information */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata;  /**< pointer to global data for lexicogrpahic order propagator */
   SCIP_VAR** vars;                          /**< pointer to variable array */
};
typedef struct vararraynodedepthbranchindex VARARRAYNODEDEPTHBRANCHINDEX;


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


/** frees dynamic lexicographic order propagator data */
static
SCIP_RETCODE lexdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA**             lexdata             /**< pointer to individual lexicographic reduction propagator datas */
   )
{
   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( (*lexdata) != NULL );

   if ( (*lexdata)->nvars > 0 )
   {
      assert( (*lexdata)->invperm != NULL );
      assert( (*lexdata)->perm != NULL );
      assert( (*lexdata)->vars != NULL );

      SCIPfreeBlockMemoryArray(scip, &(*lexdata)->invperm, (*lexdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*lexdata)->perm, (*lexdata)->nvars);
      SCIPfreeBlockMemoryArray(scip, &(*lexdata)->vars, (*lexdata)->nvars);
   }
#ifndef NDEBUG
   else
   {
      assert( (*lexdata)->nvars == 0 );
      assert( (*lexdata)->invperm == NULL );
      assert( (*lexdata)->perm == NULL );
      assert( (*lexdata)->vars == NULL );
   }
#endif
   SCIPfreeBlockMemory(scip, lexdata);

   return SCIP_OKAY;
}


/** creates dynamic lexicographic order propagator data
 *
 * fixed points are removed.
 */
static
SCIP_RETCODE lexdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   LEXDATA**             lexdata,            /**< pointer to store the data for this permutation */
   SCIP_VAR*const*       vars,               /**< input variables of the lexicographic reduction propagator */
   int                   nvars,              /**< input number of variables of the lexicographic reduction propagator */
   int*                  perm                /**< input permutation of the lexicographic reduction propagator */
   )
{
   SCIP_VAR* var;
   int* indexcorrection;
   int naffectedvariables;
   int i;
   int j;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( vars != NULL );
   assert( nvars >= 0 );
   assert( perm != NULL );
   assert( SCIPisTransformed(scip) );
   assert( masterdata->shadowtreeeventhdlr != NULL );

   /* initialize the data structures */
   SCIP_CALL( SCIPallocBlockMemory(scip, lexdata) );

   /* remove fixed points */
   naffectedvariables = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &indexcorrection, nvars) );
   for (i = 0; i < nvars; ++i)
   {
      if ( perm[i] == i )
         indexcorrection[i] = -1; /* fixed points get an entry < 0 in the indexcorrection array */
      else
         indexcorrection[i] = naffectedvariables++;
   }

   /* if this is a trivial symresack, do nothing */
   if ( naffectedvariables <= 0 )
   {
      assert( naffectedvariables == 0 );
      SCIPfreeBufferArray(scip, &indexcorrection);

      (*lexdata)->nvars = 0;
      (*lexdata)->vars = NULL;
      (*lexdata)->perm = NULL;
      (*lexdata)->invperm = NULL;
      return SCIP_OKAY;
   }

   /* initialize variable arrays */
   (*lexdata)->nvars = naffectedvariables;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->vars, (*lexdata)->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->perm, (*lexdata)->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->invperm, (*lexdata)->nvars) );

   /* determine the vars and perm */
   for (j = 0; j < nvars; ++j)
   {
      i = indexcorrection[j];
      if ( i < 0 )
         continue;

      /* j is the original index, i is the relabeled index */
      (*lexdata)->vars[i] = vars[j];
      (*lexdata)->perm[i] = indexcorrection[perm[j]];
      assert( perm[j] != j );
      assert( (*lexdata)->perm[i] != i );
      assert( (*lexdata)->perm[i] >= 0 );
      assert( (*lexdata)->perm[i] < (*lexdata)->nvars );
   }

   /* determine invperm */
   for (i = 0; i < (*lexdata)->nvars; ++i)
      (*lexdata)->invperm[(*lexdata)->perm[i]] = i;
   SCIPfreeBufferArray(scip, &indexcorrection);

   /* get transformed variables, if we are in the transformed problem */
   /* Make sure that all variables cannot be multiaggregated */
   for (i = 0; i < (*lexdata)->nvars; ++i)
   {
      assert( SCIPvarIsTransformed((*lexdata)->vars[i]) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*lexdata)->vars[i]) );
   }

   /* create hashmap for all variables */
   assert( (masterdata->symvarmap == NULL) == (masterdata->nsymvars == 0) );
   if ( masterdata->symvarmap == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&masterdata->symvarmap, SCIPblkmem(scip), (*lexdata)->nvars) );
   }
   assert( masterdata->symvarmap != NULL );
   for (i = 0; i < (*lexdata)->nvars; ++i)
   {
      assert( vars[i] != NULL );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*lexdata)->vars[i], &var) );

      /* variable is not transformed yet */
      if ( var == NULL )
         continue;

      /* var already added to hashmap */
      if ( SCIPhashmapExists(masterdata->symvarmap, (void*) var) )
         continue;

      SCIP_CALL( SCIPhashmapInsertInt(masterdata->symvarmap, (void*) var, masterdata->nsymvars++) );
   }

   return SCIP_OKAY;
}


/** comparator for sorting array of NODEDEPTHBRANCHINDEX by depth, and then by index */
static
SCIP_DECL_SORTINDCOMP(sortbynodedepthbranchindices)
{
   /* unpack the dataptr */
   VARARRAYNODEDEPTHBRANCHINDEX* vararraynodedepthbranchindices;
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices;
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata;
   SCIP_VAR** vars;
   NODEDEPTHBRANCHINDEX* index1;
   NODEDEPTHBRANCHINDEX* index2;

   vararraynodedepthbranchindices = (VARARRAYNODEDEPTHBRANCHINDEX*) dataptr;
   nodedepthbranchindices = vararraynodedepthbranchindices->nodedepthbranchindices;
   masterdata = vararraynodedepthbranchindices->masterdata;
   vars = vararraynodedepthbranchindices->vars;

   /* comparing the same element is an identity operation */
   if ( ind1 == ind2 )
      return 0;
   assert( ind1 != ind2 );

   /* sort by depth, then by index, in increasing order */
   index1 = &(nodedepthbranchindices[SCIPhashmapGetImageInt(masterdata->symvarmap, vars[ind1])]);
   index2 = &(nodedepthbranchindices[SCIPhashmapGetImageInt(masterdata->symvarmap, vars[ind2])]);

   if ( index1->nodedepth < index2->nodedepth )
      return -1;
   if ( index1->nodedepth > index2->nodedepth )
      return 1;
   assert( index1->index != index2->index );

   /* depth is the same, sort by index */
   if ( index1->index < index2->index )
      return -1;
   if ( index1->index > index2->index )
      return 1;

   /* this may not happen, since all elements must be different */
   assert( index1->index == index2->index );
   assert( FALSE );
   return 0;
}


/** get the variable ordering based on the branching decisions at the node */
static
SCIP_RETCODE getVarOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array with (depth, index)-information per variable in
                                                   * rooted path to focus node */
   int                   nvarstotal,         /**< length of that array */
   int*                  varorder,           /**< array to populate with variable order */
   int*                  nselvars,           /**< pointer to number of variables in the ordering */
   int                   maxnselvars         /**< maximal size of the number of selected variables */
)
{
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR* var;
   int varindex;
   int i;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( lexdata != NULL );
   assert( nodedepthbranchindices != NULL );
   assert( nvarstotal >= 0 );
   assert( varorder != NULL );
   assert( nselvars != NULL );

   vars = lexdata->vars;
   assert( vars != NULL );
   nvars = lexdata->nvars;
   assert( nvars >= 0 );

   /* first collect every variable that was branched on */
   *nselvars = 0;

   for (i = 0; i < nvars; ++i)
   {
      var = vars[i];
      assert( var != NULL );

      assert( SCIPhashmapExists(masterdata->symvarmap, (void*) var) );
      varindex = SCIPhashmapGetImageInt(masterdata->symvarmap, var);
      assert( varindex >= 0 );
      assert( varindex < masterdata->nsymvars );

      assert( nodedepthbranchindices[varindex].nodedepth >= 0 );

      /* skip variables that have not been used for branching */
      if ( nodedepthbranchindices[varindex].nodedepth <= 0 )
         continue;

      /* add index i of branching variable */
      assert( i >= 0 );
      assert( i < nvars );
      varorder[(*nselvars)++] = i;
   }

   if ( *nselvars > 1 )
   {
      /* Sort the first n elements of varorder by depth, then by index, as indicated by nodedepthbranchindices. */
      VARARRAYNODEDEPTHBRANCHINDEX vararraynodedepthbranchindices;
      vararraynodedepthbranchindices.nodedepthbranchindices = nodedepthbranchindices;
      vararraynodedepthbranchindices.masterdata = masterdata;
      vararraynodedepthbranchindices.vars = vars;
      SCIPsortInd(varorder, sortbynodedepthbranchindices, (void*) &vararraynodedepthbranchindices, *nselvars);
   }

   return SCIP_OKAY;
}


/** check if the static symresack with a certain variable ordering is feasible in the hypothetical scenario where
 * variables with indices i and j are fixed to fixvalue (i.e., peeking) */
static
SCIP_RETCODE peekStaticSymresackIsFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   int*                  varorder,           /**< array populated with variable order */
   int                   nselvars,           /**< number of variables in the ordering */
   int                   fixi,               /**< variable index of left fixed column */
   int                   fixj,               /**< variable index of right fixed column */
   int                   fixrow,             /**< row index in var ordering, from which to start peeking */
   SCIP_Real             fixvalue,           /**< value on which variables i and j are fixed */
   SCIP_Bool*            peekfeasible        /**< whether this is feasible or not */
)
{
   SCIP_Real* peeklbs;
   SCIP_Real* peekubs;
   SCIP_Bool* peekboundspopulated;
   int row;
   int i;
   int j;
   SCIP_VAR* vari;
   SCIP_VAR* varj;
   SCIP_Real lbi;
   SCIP_Real lbj;
   SCIP_Real ubi;
   SCIP_Real ubj;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( lexdata->vars != NULL );
   assert( lexdata->nvars >= 0 );
   assert( nselvars <= lexdata->nvars );
   assert( varorder != NULL );
   assert( nselvars > 0 );
   assert( fixi >= 0 );
   assert( fixi < lexdata->nvars );
   assert( fixj < lexdata->nvars );
   assert( fixi != fixj );
   assert( fixrow >= 0 );
   assert( fixrow < nselvars );
   assert( peekfeasible != NULL );
   assert( varorder[fixrow] == fixi );
   assert( lexdata->invperm[varorder[fixrow]] == fixj );
   assert( lexdata->perm[fixj] == fixi );
   assert( fixi == varorder[fixrow] );
   assert( fixj == lexdata->invperm[varorder[fixrow]] );

   *peekfeasible = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &peeklbs, lexdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &peekubs, lexdata->nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &peekboundspopulated, lexdata->nvars) );
#ifndef NDEBUG
   /* it must be clean */
   for (i = 0; i < lexdata->nvars; ++i)
   {
      assert( peekboundspopulated[i] == FALSE );
   }
#endif

   peeklbs[fixi] = fixvalue;
   peeklbs[fixj] = fixvalue;
   peekubs[fixi] = fixvalue;
   peekubs[fixj] = fixvalue;
   peekboundspopulated[fixi] = TRUE;
   peekboundspopulated[fixj] = TRUE;

   for (row = fixrow + 1; row < nselvars; ++row)
   {
      /* get left and right column indices */
      i = varorder[row];
      j = lexdata->invperm[i];
      assert( i == lexdata->perm[j] );

      /* no fixed points */
      assert( i != j );

      /* data checking */
      assert( i >= 0 );
      assert( i < lexdata->nvars );
      assert( j >= 0 );
      assert( j < lexdata->nvars );

      /* receive variables */
      vari = lexdata->vars[i];
      varj = lexdata->vars[j];
      assert( vari != varj );
      assert( vari != NULL );
      assert( varj != NULL );

      /* receive bounds */
      if ( peekboundspopulated[i] )
      {
         lbi = peeklbs[i];
         ubi = peekubs[i];
      }
      else
      {
         lbi = peeklbs[i] = SCIPvarGetLbLocal(vari);
         ubi = peekubs[i] = SCIPvarGetUbLocal(vari);
         peekboundspopulated[i] = TRUE;
      }
      assert( LE(scip, lbi, ubi) );

      if ( peekboundspopulated[j] )
      {
         lbj = peeklbs[j];
         ubj = peekubs[j];
      }
      else
      {
         lbj = peeklbs[j] = SCIPvarGetLbLocal(varj);
         ubj = peekubs[j] = SCIPvarGetUbLocal(varj);
         peekboundspopulated[j] = TRUE;
      }
      assert( LE(scip, lbj, ubj) );

      /* propagate that vari >= varj */

      /* if the maximal value of vari is smaller than the minimal value of varj, then vari >= varj can never hold. */
      if ( LT(scip, ubi, lbj) )
      {
         *peekfeasible = FALSE;
         SCIPdebugMessage("PEEK: Detected infeasibility at row (%3d): upper bound of %12s (%5.2f) "
            "is smaller than lower bound of %12s (%5.2f)\n",
            row, SCIPvarGetName(vari), ubi, SCIPvarGetName(varj), lbj);
         break;
      }

      /* propagate lower bound for vari */
      if ( LT(scip, lbi, lbj) )  /* lbi < lbj */
      {
         lbi = peeklbs[i] = lbj;
         assert( LE(scip, lbi, ubi) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* propagate upper bound for varj */
      if ( LT(scip, ubi, ubj) )  /* ubi < ubj */
      {
         ubj = peekubs[j] = ubi;
         assert( LE(scip, lbj, ubj) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* if there exists a solution with vari > varj, the symresack is feasible */
      if ( GT(scip, ubi, lbj) )
         break;
   }

   if (row >= nselvars)
      row = nselvars - 1;

   /* clean the peekboundspopulated array by undoing the above loop and the assignments at row fixrow */
   for (; row >= fixrow; --row)
   {
      /* left and right column indices */
      i = varorder[row];
      j = lexdata->invperm[i];
      assert( i == lexdata->perm[j] );

      /* no fixed points */
      assert( i != j );

      /* data checking */
      assert( i >= 0 );
      assert( i < lexdata->nvars );
      assert( j >= 0 );
      assert( j < lexdata->nvars );

      peekboundspopulated[i] = FALSE;
      peekboundspopulated[j] = FALSE;
   }

#ifndef NDEBUG
   /* it must be clean */
   for (i = 0; i < lexdata->nvars; ++i)
   {
      assert( peekboundspopulated[i] == FALSE );
   }
#endif

   SCIPfreeCleanBufferArray(scip, &peekboundspopulated);
   SCIPfreeBufferArray(scip, &peekubs);
   SCIPfreeBufferArray(scip, &peeklbs);

   return SCIP_OKAY;
}


/** propagate static symresack with specified variable ordering */
static
SCIP_RETCODE propagateStaticSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   int*                  varorder,           /**< array populated with variable order */
   int                   nselvars,           /**< number of variables in the ordering */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
)
{
   int row;
   int i;
   int j;
   SCIP_VAR* vari;
   SCIP_VAR* varj;

   SCIP_Real lbi;
   SCIP_Real ubi;
   SCIP_Real lbj;
   SCIP_Real ubj;
   SCIP_Bool success;
   SCIP_Bool peekfeasible;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( varorder != NULL );
   assert( nselvars >= 0 );
   assert( infeasible != NULL );
   assert( !*infeasible );
   assert( nreductions != NULL );
   assert( *nreductions >= 0 );

   /* avoid trivial cases */
   if ( nselvars <= 0 )
      return SCIP_OKAY;

   /* iterate over the symresack entrywise. We see this as two columns, with the left vector being the variable
    * ordering, and the right column the permuted variables of this var ordering. */
   for (row = 0; row < nselvars; ++row)
   {
      /* left and right column indices */
      i = varorder[row];
      j = lexdata->invperm[i];
      assert( i == lexdata->perm[j] );

      /* no fixed points */
      assert( i != j );

      /* data checking */
      assert( i >= 0 );
      assert( i <= lexdata->nvars );
      assert( j >= 0 );
      assert( j <= lexdata->nvars );

      /* receive variables */
      vari = lexdata->vars[i];
      varj = lexdata->vars[j];
      assert( vari != varj );
      assert( vari != NULL );
      assert( varj != NULL );

      /* receive bounds */
      lbi = SCIPvarGetLbLocal(vari);
      lbj = SCIPvarGetLbLocal(varj);
      ubi = SCIPvarGetUbLocal(vari);
      ubj = SCIPvarGetUbLocal(varj);
      assert( LE(scip, lbi, ubi) );
      assert( LE(scip, lbj, ubj) );

      /* propagate that vari >= varj */

      /* if the maximal value of vari is smaller than the minimal value of varj, then vari >= varj can never hold */
      if ( LT(scip, ubi, lbj) )
      {
         *infeasible = TRUE;
         SCIPdebugMessage("Detected infeasibility at row (%3d): upper bound of %12s (%5.2f) "
            "is smaller than lower bound of %12s (%5.2f)\n",
            row, SCIPvarGetName(vari), ubi, SCIPvarGetName(varj), lbj);
         goto FREE;
      }

      /* propagate lower bound for vari */
      if ( LT(scip, lbi, lbj) )
      {
         SCIP_CALL( SCIPtightenVarLb(scip, vari, lbj, FALSE, infeasible, &success) );
         if ( success )
         {
            SCIPdebugMessage("Restricting variable LB %12s (%3d) to %5.2f\n", SCIPvarGetName(vari), row, lbj);
            *nreductions += 1;
         }
         else
         {
            SCIPdebugMessage("Restricting variable LB %12s (%3d) to %5.2f (no success)\n",
               SCIPvarGetName(vari), row, lbj);
         }
         if ( *infeasible )
         {
            SCIPdebugMessage("Detected infeasibility restricting variable LB %12s (%3d) to %5.2f\n",
               SCIPvarGetName(vari), row, lbj);
            goto FREE;
         }
         /* for later reference, update this bound change */
         lbi = lbj;
         assert( LE(scip, lbi, ubi) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* propagate upper bound for varj */
      if ( LT(scip, ubi, ubj) )
      {
         SCIP_CALL( SCIPtightenVarUb(scip, varj, ubi, FALSE, infeasible, &success) );
         if ( success )
         {
            SCIPdebugMessage("Restricting variable UB %12s (%3d) to %5.2f\n", SCIPvarGetName(varj), row, ubi);
            *nreductions += 1;
         }
         else
         {
            SCIPdebugMessage("Restricting variable UB %12s (%3d) to %5.2f (no success)\n",
               SCIPvarGetName(varj), row, ubi);
         }
         if ( *infeasible )
         {
            SCIPdebugMessage("Detected infeasibility restricting variable UB %12s (%3d) to %5.2f\n",
               SCIPvarGetName(varj), row, ubi);
            goto FREE;
         }
         /* for later reference, update this bound change */
         ubj = ubi;
         assert( LE(scip, lbj, ubj) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* terminate if there exists a solution being lexicographically strictly larger than its image */
      if ( GT(scip, ubi, lbj) )
         break;
   }

   if ( row < nselvars )
   {
      /* If we get here, the previous loop is broken at row "row", which allowed for choosing vari > varj.
       * Question: Is vari == varj permitted? If not, tighten even further.
       */

      /* REMARK: For a symresack (without knowing it's an orbisack actually), we do the feasibility check twice.
       * If it is an orbisack, instead, then it suffices to do it only once. Make distinction here.
       */

      /* Case 1: If vari is minimal (lbi).
       * Then, propagation of lbi = vari >= varj can yield two situations:
       *   Option 1: varj can take a value < lbi. Then no further fixings can be detected.
       *   Option 2: varj gets fixed to lbi. Then, we must check if feasibility is found, still.
       *     If it turns out infeasible, then we know vari cannot take value lbi, so we can increase the lower bound.
       */
      assert( LE(scip, lbj, lbi) );  /* this must be the case after reductions in the for-loop */
      if ( EQ(scip, lbj, lbi) )
      {
         /* this is Option 2: varj gets fixed to lbi by propagation. */
         SCIP_CALL( peekStaticSymresackIsFeasible(scip, lexdata, varorder, nselvars, i, j, row, lbi, &peekfeasible) );
         if ( !peekfeasible )
         {
            /* vari cannot take value lbi, so we increase the lower bound. */
            switch ( SCIPvarGetType(vari) )
            {
               case SCIP_VARTYPE_BINARY:
               case SCIP_VARTYPE_IMPLINT:
               case SCIP_VARTYPE_INTEGER:
                  /* discrete variable type: increase lower bound by 1. */
                  assert( SCIPisIntegral(scip, lbi) );
                  SCIP_CALL( SCIPtightenVarLb(scip, vari, lbi + 1.0, FALSE, infeasible, &success) );
                  if ( success )
                     *nreductions += 1;
                  if ( *infeasible )
                     goto FREE;
                  lbi = lbi + 1.0;
                  assert( LE(scip, lbi, ubi) );
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous variable type: act as if we increase the variable by a very little bit.
                   * that is only possible if we're able to increase the variable bound by a bit. */
                  if ( EQ(scip, lbi, ubi) )
                  {
                     *infeasible = TRUE;
                     goto FREE;
                  }
                  break;
               default:
                  return SCIP_ERROR;
            }
         }
      }

      /* Case 2: varj is maximal (ubj).
       * Then, propagation of vari >= varj = ubj can yield two situatiosn:
       *   Option 1: vari can take a value > ubj. Then, no further fixings can be detected.
       *   Option 2: vari gets fixed to ubj. Then, we must check if feasibility is found, still.
       *     If it turns out infeasible, then we know varj cannot take value ubj, so we can decrease the upper bound.
       */
      assert( GE(scip, ubi, ubj) );  /* this must be the case after reductions in the for-loop */
      if ( EQ(scip, ubi, ubj) )
      {
         /* this is Option 2: vari gets fixed to ubj by propagation. */
         SCIP_CALL( peekStaticSymresackIsFeasible(scip, lexdata, varorder, nselvars, i, j, row, ubj, &peekfeasible) );
         if ( !peekfeasible )
         {
            /* varj cannot take value ubj, so we increase the upper bound. */
            switch ( SCIPvarGetType(varj) )
            {
               case SCIP_VARTYPE_BINARY:
               case SCIP_VARTYPE_IMPLINT:
               case SCIP_VARTYPE_INTEGER:
                  /* discrete variable type: decrease upper bound by 1. */
                  assert( SCIPisIntegral(scip, ubj) );
                  SCIP_CALL( SCIPtightenVarUb(scip, varj, ubj - 1.0, FALSE, infeasible, &success) );
                  if ( success )
                     *nreductions += 1;
                  if ( *infeasible )
                     goto FREE;
                  ubj = ubj - 1.0;
                  assert( LE(scip, lbj, ubj) );
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous variable type: act as if we decrease the variable by a very little bit.
                   * that is only possible if we're able to decrease the variable bound by a bit. */
                  if ( EQ(scip, lbj, ubj) )
                  {
                     *infeasible = TRUE;
                     goto FREE;
                  }
                  break;
               default:
                  return SCIP_ERROR;
            }
         }
      }
   }

   FREE:
   ;

   return SCIP_OKAY;
}


/** propagation method for a dynamic symresack */
static
SCIP_RETCODE propagateSymresackDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array with (depth, index)-information per variable in
                                                   * rooted path to focus node */
   int                   nvarstotal,         /**< length of that array */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
   )
{
   int* varorder;
   int nvarorder;
   int nvars;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( lexdata != NULL );

   nvars = lexdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &varorder, nvars) );

   SCIP_CALL( getVarOrder(scip, masterdata, lexdata, nodedepthbranchindices, nvarstotal,
      varorder, &nvarorder, nvars) );
   assert( nvarorder >= 0 );
   assert( nvarorder <= nvars );

   /* if none of the variables of the symresack is in the varorder, then no propagation is needed */
   if ( nvarorder == 0 )
      goto FREE;

   /* propagate static symresack with this variable order */
   SCIP_CALL( propagateStaticSymresack(scip, lexdata, varorder, nvarorder, infeasible, nreductions) );

   FREE:
   SCIPfreeBufferArray(scip, &varorder);

   return SCIP_OKAY;
}


/** propagation method for applying dynamic lexicographic reduction for a single permutation */
static
SCIP_RETCODE propagateLexicographicReductionPerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array with (depth, index)-information per variable in
                                                   * rooted path to focus node */
   int                   nvarstotal,         /**< length of that array */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   *nreductions = 0;
   *infeasible = FALSE;

   SCIP_CALL( propagateSymresackDynamic(scip, masterdata, lexdata,
      nodedepthbranchindices, nvarstotal, infeasible, nreductions) );

   return SCIP_OKAY;
}


/** assuming nodedepthbranchindices is initially clean, populate array with information of first variable change */
static
SCIP_RETCODE shadowtreeFillNodeDepthBranchIndices(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array to populate */
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to shadow tree structure */
   SCIP_NODE*            focusnode           /**< focusnode to which the rooted path is evaluated */
)
{
   SCIP_SHADOWNODE* shadownode;
   SCIP_SHADOWNODE* shadowchild;
   int shadowdepth;
   SCIP_VAR* var;
   int varindex;
   int nlevelvars;
   int c;
   int i;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( masterdata->symvarmap != NULL );
   assert( masterdata->nsymvars >= 0 );
   assert( nodedepthbranchindices != NULL );
   assert( shadowtree != NULL );
   assert( focusnode != NULL );

   shadownode = SCIPshadowtreeGetShadowNode(shadowtree, focusnode);
   assert( shadownode != NULL );
   shadowdepth = SCIPnodeGetDepth(focusnode);

   /* we start looking one level lower, because we consider the branching decisions each time. */
   shadownode = shadownode->parent;

   /* now, walk from the leaf to the root. Each time look at all the children of the node considered,
    * and save the variable depth and index in the branching array. It is important to consider all children each time,
    * because we need to comply with the instance where in different branches it is branched on different variables.
    * This has to be consistent. */
   while (shadownode != NULL)
   {
      assert( shadowdepth > 0 );
      nlevelvars = 0;
      for (c = 0; c < shadownode->nchildren; ++c)
      {
         shadowchild = shadownode->children[c];
         /* get the variables branched on, for each of the children (that's likely 1 variable each time) */
         for (i = 0; i < shadowchild->nbranchingdecisions; ++i)
         {
            var = shadowchild->branchingdecisions[i].var;
            /* ignore variables not relevant for symresacks */
            if ( !SCIPhashmapExists(masterdata->symvarmap, (void*) var) )
               continue;
            assert( SCIPhashmapExists(masterdata->symvarmap, (void*) var) );

            varindex = SCIPhashmapGetImageInt(masterdata->symvarmap, var);
            assert( varindex >= 0 );
            assert( varindex < masterdata->nsymvars );

            /* var already in other child at this level? Continue */
            if ( nodedepthbranchindices[varindex].nodedepth == shadowdepth )
               continue;

            /* the variable is either not seen (nodedepth == 0), or it is at a higher level (nodedepth > shadowdepth) */
            assert( nodedepthbranchindices[varindex].nodedepth == 0 ||
               nodedepthbranchindices[varindex].nodedepth > shadowdepth);

            /* var is not seen on this level yet. Update */
            nodedepthbranchindices[varindex].nodedepth = shadowdepth;
            nodedepthbranchindices[varindex].index = nlevelvars++;
         }
      }

      /* prepare for the next iteration */
      shadownode = shadownode->parent;
      --shadowdepth;
   }
   /* In the last iteration, we handled the branching decisions at the root node, so shadowdepth must have value 0. */
   assert( shadowdepth == 0 );

   return SCIP_OKAY;
}


/** cleans nodedepthbranchindices array */
static
SCIP_RETCODE shadowtreeUndoNodeDepthBranchIndices(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array to populate */
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to shadow tree structure */
   SCIP_NODE*            focusnode           /**< focusnode to which the rooted path is evaluated */
)
{
   /* Do the same as above, but each time we see a variable we NULL all the fields, cleaning the object */
   SCIP_SHADOWNODE* shadownode;
   SCIP_SHADOWNODE* shadowchild;
   int shadowdepth;
   SCIP_VAR* var;
   int varindex;
   int c;
   int i;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( masterdata->symvarmap != NULL );
   assert( masterdata->nsymvars >= 0 );
   assert( nodedepthbranchindices != NULL );
   assert( shadowtree != NULL );
   assert( focusnode != NULL );

   shadownode = SCIPshadowtreeGetShadowNode(shadowtree, focusnode);
   shadowdepth = SCIPnodeGetDepth(focusnode);

   /* we start looking one level lower, because we consider the branching decisions each time. */
   shadownode = shadownode->parent;

   /* now, walk from the leaf to the root. Each time look at all the children of the node considered,
    * and save the variable depth and index in the branching array. It is important to consider all children each time,
    * because we need to comply with the instance where in different branches it is branched on different variables.
    * This has to be consistent. */
   while (shadownode != NULL)
   {
      assert( shadowdepth > 0 );
      for (c = 0; c < shadownode->nchildren; ++c)
      {
         shadowchild = shadownode->children[c];
         /* get the variables branched on, for each of the children (that's likely 1 variable each time) */
         for (i = 0; i < shadowchild->nbranchingdecisions; ++i)
         {
            var = shadowchild->branchingdecisions[i].var;
            /* ignore variables not relevant for symresacks */
            if ( !SCIPhashmapExists(masterdata->symvarmap, (void*) var) )
               continue;
            assert( SCIPhashmapExists(masterdata->symvarmap, (void*) var) );

            varindex = SCIPhashmapGetImageInt(masterdata->symvarmap, var);
            assert( varindex >= 0 );
            assert( varindex < masterdata->nsymvars );

            /* reset */
            nodedepthbranchindices[varindex].index = 0;
            nodedepthbranchindices[varindex].nodedepth = 0;
         }
      }

      /* prepare for the next iteration */
      shadownode = shadownode->parent;
      --shadowdepth;
   }
   /* In the last iteration, we handled the branching decisions at the root node, so shadowdepth must have value 0. */
   assert( shadowdepth == 0 );

   return SCIP_OKAY;
}


/*
 * Interface methods
 */

/** apply lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   SCIP_Bool*            infeasible,         /**< whether infeasibility is found */
   int*                  nred,               /**< number of domain reductions */
   SCIP_Bool*            didrun              /**< whether propagator actually ran */
   )
{
   int nlocalred; 
   int p;
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices;
   SCIP_SHADOWTREE* shadowtree;
   SCIP_NODE* focusnode;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( (masterdata->lexdatas == NULL) == (masterdata->nlexdatas == 0) );
   assert( masterdata->nlexdatas >= 0 );
   assert( masterdata->nlexdatas <= masterdata->maxnlexdatas );
   assert( infeasible != NULL );
   assert( nred != NULL );

   *infeasible = FALSE;
   *nred = 0;

   /* early termination */
   if ( masterdata->nlexdatas == 0 )
      return SCIP_OKAY;

   /* compute the variable ordering based on the branching decisions of the shadowtree */
   assert( masterdata->shadowtreeeventhdlr != NULL );
   shadowtree = SCIPgetShadowTree(masterdata->shadowtreeeventhdlr);
   focusnode = SCIPgetFocusNode(scip);

   /* fill the node-depth-branch-indices structure
    * this is an array that maps every variable index to (depth, index) = (0, 0) if the variable is not branched on,
    * or (depth, index) is the depth (branching at root node: depth 1) and variable index when it was branched thereon.
    * For individual dynamic symresacks, we use this ordering the following way:
    *  1. Choose those variables that have (depth, index) with depth > 0 (these)
    *  2. Sort by depth, then by index, in increasing order.
    */
   SCIPallocCleanBufferArray(scip, &nodedepthbranchindices, masterdata->nsymvars);
   SCIP_CALL( shadowtreeFillNodeDepthBranchIndices(scip, masterdata, nodedepthbranchindices, shadowtree, focusnode) );
   /* ... Do everything using this nodedepthbranchindices structure */

   /* propagate all dynamic lexicographic order constraints */
   for (p = 0; p < masterdata->nlexdatas; ++p)
   {
      assert( masterdata->lexdatas[p] != NULL );

      SCIP_CALL( propagateLexicographicReductionPerm(scip, masterdata, masterdata->lexdatas[p], 
         nodedepthbranchindices, masterdata->nsymvars, infeasible, &nlocalred) );

      /* keep track of the total number of fixed vars */
      *nred += nlocalred;
      *didrun = TRUE;

      /* stop if we find infeasibility */
      if ( *infeasible )
         break;
   }

   /* clean the node-depth-branch-indices structure */
   SCIP_CALL( shadowtreeUndoNodeDepthBranchIndices(scip, masterdata, nodedepthbranchindices, shadowtree, focusnode) );
   SCIPfreeCleanBufferArray(scip, &nodedepthbranchindices);
   
   return SCIP_OKAY;
}


/** adds permutation for lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionAddPermutation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int*                  perm                /**< permutation */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( (masterdata->lexdatas == NULL) == (masterdata->nlexdatas == 0) );
   assert( masterdata->nlexdatas >= 0 );
   assert( masterdata->nlexdatas <= masterdata->maxnlexdatas );
   assert( masterdata->shadowtreeeventhdlr != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( perm != NULL );
   assert( SCIPisTransformed(scip) );

   /* resize component array if needed */
   if ( masterdata->nlexdatas == masterdata->maxnlexdatas )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, masterdata->nlexdatas + 1);
      assert( newsize >= 0 );

      if ( masterdata->nlexdatas == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &masterdata->lexdatas, newsize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &masterdata->lexdatas, 
            masterdata->maxnlexdatas, newsize) );
      }

      masterdata->maxnlexdatas = newsize;
   }

   /* prepare lexdatas */
   SCIP_CALL( lexdataCreate(scip, masterdata, &masterdata->lexdatas[masterdata->nlexdatas++], 
      permvars, npermvars, perm) );

   return SCIP_OKAY;
}


/** resets lexicographic reduction propagation (removes all permutations) */
SCIP_RETCODE SCIPlexicographicReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA* masterdata/**< pointer to global data for lexicogrpahic order propagator */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( (masterdata->lexdatas == NULL) == (masterdata->nlexdatas == 0) );
   assert( masterdata->nlexdatas >= 0 );
   assert( masterdata->nlexdatas <= masterdata->maxnlexdatas );
   assert( masterdata->shadowtreeeventhdlr != NULL );

   while ( masterdata->nlexdatas > 0 )
   {
      SCIP_CALL( lexdataFree(scip, &(masterdata->lexdatas[--masterdata->nlexdatas])) );
   }
   assert( masterdata->nlexdatas == 0 );

   SCIPfreeBlockMemoryArray(scip, &masterdata->lexdatas, masterdata->maxnlexdatas);
   masterdata->lexdatas = NULL;
   masterdata->maxnlexdatas = 0;

   assert( masterdata->nsymvars >= 0 );
   assert( (masterdata->symvarmap == NULL) == (masterdata->nsymvars == 0) );
   if ( masterdata->symvarmap != NULL )
   {
      SCIPhashmapFree(&masterdata->symvarmap);
      masterdata->symvarmap = NULL;
      masterdata->nsymvars = 0;
   }

   return SCIP_OKAY;
}


/** free lexicographic reduction data */
SCIP_RETCODE SCIPlexicographicReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA** masterdata/**< pointer to global data for lexicogrpahic order propagator */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( *masterdata != NULL );

   SCIP_CALL( SCIPlexicographicReductionReset(scip, *masterdata) );
   assert( (*masterdata)->lexdatas == NULL );
   assert( (*masterdata)->symvarmap == NULL );

   SCIPfreeBlockMemory(scip, masterdata);
   return SCIP_OKAY;
}


/** initializes structures needed for lexicographic reduction propagation
 * This is only done exactly once.
 */
SCIP_RETCODE SCIPlexicographicReductionInclude(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXICOGRAPHICREDUCTIONDATA** masterdata,/**< pointer to global data for lexicogrpahic order propagator */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( shadowtreeeventhdlr != NULL );

   assert( SCIPgetStage(scip) == SCIP_STAGE_INIT );

   SCIP_CALL( SCIPallocBlockMemory(scip, masterdata) );

   (*masterdata)->shadowtreeeventhdlr = shadowtreeeventhdlr;
   (*masterdata)->symvarmap = NULL;
   (*masterdata)->nsymvars = 0;
   (*masterdata)->lexdatas = NULL;
   (*masterdata)->nlexdatas = 0;
   (*masterdata)->maxnlexdatas = 0;

   return SCIP_OKAY;
}
