/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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
#include "scip/debug.h"
#include "scip/struct_scip.h"
#include "scip/struct_mem.h"
#include "scip/struct_tree.h"
#include "scip/symmetry.h"
#include "scip/event_shadowtree.h"
#include <ctype.h>
#include <string.h>
#include <memory.h>


/*
 * Data structures
 */


/** data per permutation for lexicographic reduction propagator */
struct LexRedPermData
{
   SCIP_VAR**            vars;               /**< variables affected by permutation */
   int                   nvars;              /**< number of variables */
   int*                  perm;               /**< permutation for lexicographic reduction */
   int*                  invperm;            /**< inverse permutation */
   SCIP_HASHMAP*         varmap;             /**< map of variables to indices in vars array */
};
typedef struct LexRedPermData LEXDATA;


/** data for dynamic lexicographic reduction propagator */
struct SCIP_LexRedData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr;/**< eventhandler for the shadow tree data structure */
   SCIP_HASHMAP*         symvarmap;          /**< map of variables affected by some permutation handled by a LEXDATA */
   int                   nsymvars;           /**< number of variables in symvarmap */
   LEXDATA**             lexdatas;           /**< array of pointers to individual LEXDATA's */
   int                   nlexdatas;          /**< number of datas in array */
   int                   maxnlexdatas;       /**< allocated datas array size */
   int                   nred;               /**< total number of reductions */
};


/** to store branch-and-bound tree paths, (depth, index)-information per variable in rooted path */
struct NodeDepthBranchIndex
{
   int                   nodedepth;          /**< depth of var domain change */
   int                   index;              /**< index of var domain change on node at depth */
};
typedef struct NodeDepthBranchIndex NODEDEPTHBRANCHINDEX;


/** auxiliary struct to pass branch-and-bound tree information to sort function */
struct VarArrayNodeDepthBranchIndex
{
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices; /**< pointer to branch-and-bound tree information */
   SCIP_LEXREDDATA*      masterdata;             /**< pointer to global data for lexicographic reduction propagator */
   SCIP_VAR**            vars;                   /**< pointer to variable array */
};
typedef struct VarArrayNodeDepthBranchIndex VARARRAYNODEDEPTHBRANCHINDEX;


/*
 * Local methods
 */

/** frees dynamic lexicographic reduction propagator data */
static
SCIP_RETCODE lexdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA**             lexdata             /**< pointer to individual lexicographic reduction propagator datas */
   )
{
   int i;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( (*lexdata) != NULL );

   if ( (*lexdata)->nvars > 0 )
   {
      assert( (*lexdata)->invperm != NULL );
      assert( (*lexdata)->perm != NULL );
      assert( (*lexdata)->vars != NULL );

      /* free hashmap */
      SCIPhashmapFree(&((*lexdata)->varmap));

      /* release variables */
      for (i = 0; i < (*lexdata)->nvars; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*lexdata)->vars[i]) );
      }

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
      assert( (*lexdata)->varmap == NULL );
   }
#endif
   SCIPfreeBlockMemory(scip, lexdata);

   return SCIP_OKAY;
}


/** creates dynamic lexicographic reduction propagator data
 *
 *  Fixed points are removed.
 */
static
SCIP_RETCODE lexdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   LEXDATA**             lexdata,            /**< pointer to store data for permutation to be added */
   SCIP_VAR*const*       vars,               /**< input variables of the lexicographic reduction propagator */
   int                   nvars,              /**< input number of variables of the lexicographic reduction propagator */
   int*                  perm,               /**< input permutation of the lexicographic reduction propagator */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   )
{
   SCIP_VAR* var;
   int* indexcorrection;
   int naffectedvariables;
   int i;
   int j;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( lexdata != NULL );
   assert( vars != NULL );
   assert( nvars >= 0 );
   assert( perm != NULL );
   assert( success != NULL );
   assert( SCIPisTransformed(scip) );
   assert( masterdata->shadowtreeeventhdlr != NULL );

   *success = TRUE;

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

   /* do nothing if reductions would be trivial */
   if ( naffectedvariables <= 0 )
   {
      assert( naffectedvariables == 0 );
      SCIPfreeBufferArray(scip, &indexcorrection);

      *success = FALSE;
      SCIPfreeBlockMemory(scip, lexdata);
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

   /* make sure that we deal with transformed variables and that variables cannot be multi-aggregated, and capture */
   for (i = 0; i < (*lexdata)->nvars; ++i)
   {
      assert( SCIPvarIsTransformed((*lexdata)->vars[i]) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*lexdata)->vars[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, (*lexdata)->vars[i]) );
   }

   /* create hashmap for all variables, both globally and just for this lexdata */
   assert( (masterdata->symvarmap == NULL) == (masterdata->nsymvars == 0) );
   if ( masterdata->symvarmap == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&masterdata->symvarmap, SCIPblkmem(scip), (*lexdata)->nvars) );
   }
   assert( masterdata->symvarmap != NULL );

   SCIP_CALL( SCIPhashmapCreate(&(*lexdata)->varmap, SCIPblkmem(scip), (*lexdata)->nvars) );
   assert( (*lexdata)->varmap != NULL );

   for (i = 0; i < (*lexdata)->nvars; ++i)
   {
      var = (*lexdata)->vars[i];
      assert( var != NULL );
      assert( SCIPvarIsTransformed(var) );

      /* the hashmap in lexdata maps to the index in the vars array */
      SCIP_CALL( SCIPhashmapInsertInt((*lexdata)->varmap, (void*) var, i) );

      /* var already added to hashmap */
      if ( SCIPhashmapExists(masterdata->symvarmap, (void*) var) )
         continue;

      /* the hashmap in symvarmap maps to a unique index */
      SCIP_CALL( SCIPhashmapInsertInt(masterdata->symvarmap, (void*) var, masterdata->nsymvars++) );
   }

   return SCIP_OKAY;
}


/** comparator used in the getVarOrder() function, for sorting an array of NODEDEPTHBRANCHINDEX by depth, then by index
 *
 *  @warning this function is only meant to be used in the getVarOrder() function
 *
 *  @pre datapointer is populated with a VARARRAYNODEDEPTHBRANCHINDEX pointer
 *  @pre the comparator is only called on entries with set (depth, index)-information
 *  @pre the (depth, index)-informations are all different
 *
 *  result:
 *    0: the same index is compared, so the (depth, index)-informations are the same
 *   -1: the depth of ind1 is smaller than the depth of ind2, or it's equal and the index is smaller
 *    1: the depth of ind2 is smaller than the depth of ind1, or it's equal and the index is smaller
 */
static
SCIP_DECL_SORTINDCOMP(sortbynodedepthbranchindices)
{
   /* unpack the dataptr */
   VARARRAYNODEDEPTHBRANCHINDEX* vararraynodedepthbranchindices;
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices;
   SCIP_LEXREDDATA* masterdata;
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

   return 0;
}


/** gets the variable ordering based on the branching decisions at the node */
static
SCIP_RETCODE getVarOrder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array with (depth, index)-information per variable in
                                                  *   rooted path to focus node */
   int                   nvarstotal,         /**< length of that array */
   SCIP_VAR**            branchvars,         /**< array populated with variables branched on in the symvarmap hashset */
   int                   nbranchvars,        /**< number of elements in branchvars array */
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
   assert( branchvars != NULL );
   assert( nbranchvars >= 0 );
   assert( varorder != NULL );
   assert( nselvars != NULL );

   vars = lexdata->vars;
   assert( vars != NULL );
   nvars = lexdata->nvars;
   assert( nvars >= 0 );

   /* first collect every variable that was branched on */
   *nselvars = 0;

   if ( nvars < nbranchvars )
   {
      /* for permutations with small support, just check all support entries */
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
         assert( (*nselvars) < maxnselvars );
         varorder[(*nselvars)++] = i;
      }
   }
   else
   {
      /* for permutations where the support is larger than the number of branched vars, check for the branched vars */
      for (i = 0; i < nbranchvars; ++i)
      {
         var = branchvars[i];
         assert( var != NULL );

#ifndef NDEBUG
         /* debugging: test if it is indeed a branched variable! */
         varindex = SCIPhashmapGetImageInt(masterdata->symvarmap, var);
         assert( varindex >= 0 );
         assert( varindex < masterdata->nsymvars );
         assert( nodedepthbranchindices[varindex].nodedepth > 0 );
#endif

         /* get the variable index in the lexdata->vars array */
         varindex = SCIPhashmapGetImageInt(lexdata->varmap, (void*) var);
         assert( varindex == INT_MAX || (varindex >= 0 && varindex < lexdata->nvars) );

         /* skip variables that are not permuted by the permutation */
         if ( varindex == INT_MAX )
            continue;
         assert( lexdata->vars[varindex] == var );

         /* add index varindex of the branching variable */
         varorder[(*nselvars)++] = varindex;
      }
   }

   if ( *nselvars > 1 )
   {
      /* sort the first n elements of varorder by depth, then by index, as indicated by nodedepthbranchindices. */
      VARARRAYNODEDEPTHBRANCHINDEX vararraynodedepthbranchindices;
      vararraynodedepthbranchindices.nodedepthbranchindices = nodedepthbranchindices;
      vararraynodedepthbranchindices.masterdata = masterdata;
      vararraynodedepthbranchindices.vars = vars;
      SCIPsortInd(varorder, sortbynodedepthbranchindices, (void*) &vararraynodedepthbranchindices, *nselvars);
   }

   return SCIP_OKAY;
}


/** checks if the static lexred with a certain variable ordering is feasible in the hypothetical scenario where
 *  variables with indices i and j are fixed to fixvalue (i.e., peeking)
 */
static
SCIP_RETCODE peekStaticLexredIsFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   int*                  varorder,           /**< array populated with variable order */
   int                   nselvars,           /**< number of variables in the ordering */
   int                   fixi,               /**< variable index of left fixed column */
   int                   fixj,               /**< variable index of right fixed column */
   int                   fixrow,             /**< row index in var ordering, from which to start peeking */
   SCIP_Real             fixvalue,           /**< value on which variables i and j are fixed */
   SCIP_Bool*            peekfeasible        /**< pointer to store whether this is feasible or not */
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
         lbi = SCIPvarGetLbLocal(vari);
         ubi = SCIPvarGetUbLocal(vari);
         peeklbs[i] = lbi;
         peekubs[i] = ubi;
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
         lbj = SCIPvarGetLbLocal(varj);
         ubj = SCIPvarGetUbLocal(varj);
         peeklbs[j] = lbj;
         peekubs[j] = ubj;
         peekboundspopulated[j] = TRUE;
      }
      assert( LE(scip, lbj, ubj) );

      /* propagate that vari >= varj */

      /* vari >= varj can never hold if the maximal value of vari is smaller than the minimal value of varj */
      if ( LT(scip, ubi, lbj) )
      {
         *peekfeasible = FALSE;
         SCIPdebugMessage("PEEK: Detected infeasibility at row (%3d): upper bound of %12s (%5.2f) "
            "is smaller than lower bound of %12s (%5.2f)\n",
            row, SCIPvarGetName(vari), ubi, SCIPvarGetName(varj), lbj);
         break;
      }

      /* propagate lower bound for vari */
      if ( LT(scip, lbi, lbj) )
      {
         lbi = lbj;
         peeklbs[i] = lbj;
         assert( LE(scip, lbi, ubi) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* propagate upper bound for varj */
      if ( LT(scip, ubi, ubj) )
      {
         ubj = ubi;
         peekubs[j] = ubi;
         assert( LE(scip, lbj, ubj) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* if there exists a solution with vari > varj, reductions are feasible w.r.t. lexred */
      if ( GT(scip, ubi, lbj) )
         break;
   }

   if ( row >= nselvars )
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


/** propagates static lexicographic reduction with specified variable ordering */
static
SCIP_RETCODE propagateStaticLexred(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   int*                  varorder,           /**< array populated with variable order */
   int                   nselvars,           /**< number of variables in the ordering */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
)
{
   int row;
   int i = -1;
   int j = -1;
   SCIP_VAR* vari = NULL;
   SCIP_VAR* varj = NULL;

   SCIP_Real lbi = 0.0;
   SCIP_Real ubi = 0.0;
   SCIP_Real lbj = 0.0;
   SCIP_Real ubj = 0.0;
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

   /* iterate over the variable array entrywise
    *
    * We see this as two columns, with the left vector being the variable ordering,
    * and the right column the permuted variables of this var ordering.
    */
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
         return SCIP_OKAY;
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
            return SCIP_OKAY;
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
            return SCIP_OKAY;
         }
         /* for later reference, update this bound change */
         ubj = ubi;
         assert( LE(scip, lbj, ubj) );  /* otherwise we returned in the `if ( LT(scip, ubi, lbj) )` block */
      }

      /* terminate if there exists a solution being lexicographically strictly larger than its image */
      if ( GT(scip, ubi, lbj) )
         break;
   }
   assert( i >= 0 );
   assert( j >= 0 );
   assert( vari != NULL );
   assert( varj != NULL );

   if ( row < nselvars )
   {
      /* The previous loop is broken at row "row", which allows for choosing vari > varj.
       *
       * Now check if vari == varj is permitted, and if not, tighten the domain further.
       *
       * @todo we peek twice if vari and varj are unfixed
       * But, if the subcycle only contains var1 and var2, a single peek suffices.
       * This is similar to orbisack and symresack propagation.
       *
       * Case 1: vari is minimal (lbi).
       * Then, propagation of lbi = vari >= varj can yield two situations:
       *   Option 1: varj can take a value < lbi. Then no further reductions can be detected.
       *   Option 2: varj gets fixed to lbi. Then, we must check if feasibility is found, still.
       *     If it turns out infeasible, then we know vari cannot take value lbi, so we can increase the lower bound.
       */

      assert( LE(scip, lbj, lbi) );  /* this must be the case after reductions in the for-loop */
      if ( EQ(scip, lbj, lbi) )
      {
         /* this is Option 2: varj gets fixed to lbi by propagation. */
         SCIP_CALL( peekStaticLexredIsFeasible(scip, lexdata, varorder, nselvars, i, j,
            row, lbi, &peekfeasible) );
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
                     return SCIP_OKAY;
                  lbi = lbi + 1.0;
                  assert( LE(scip, lbi, ubi) );
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous variable type: act as if we increase the variable by a very little bit.
                   * That is only possible if we're able to increase the variable bound by a bit.
                   */
                  if ( EQ(scip, lbi, ubi) )
                  {
                     *infeasible = TRUE;
                     return SCIP_OKAY;
                  }
                  break;
               default:
                  SCIPerrorMessage("unsupported variable type encountered at the lexicographic reduction propagator\n");
                  return SCIP_ERROR;
            }
         }
      }

      /* Case 2: varj is maximal (ubj).
       * Then, propagation of vari >= varj = ubj can yield two situatiosn:
       *   Option 1: vari can take a value > ubj. Then, no further reductions can be detected.
       *   Option 2: vari gets fixed to ubj. Then, we must check if feasibility is found, still.
       *     If it turns out infeasible, then we know varj cannot take value ubj, so we can decrease the upper bound.
       */
      assert( GE(scip, ubi, ubj) );  /* this must be the case after reductions in the for-loop */
      if ( EQ(scip, ubi, ubj) )
      {
         /* this is Option 2: vari gets fixed to ubj by propagation. */
         SCIP_CALL( peekStaticLexredIsFeasible(scip, lexdata, varorder, nselvars, i, j, row, ubj, &peekfeasible) );
         if ( !peekfeasible )
         {
            /* varj cannot take value ubj, so we decrease the upper bound. */
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
                     return SCIP_OKAY;
                  ubj = ubj - 1.0;
                  assert( LE(scip, lbj, ubj) );
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous variable type: act as if we decrease the variable by a very little bit.
                   * that is only possible if we're able to decrease the variable bound by a bit. */
                  if ( EQ(scip, lbj, ubj) )
                  {
                     *infeasible = TRUE;
                     return SCIP_OKAY;
                  }
                  break;
               default:
                  return SCIP_ERROR;
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** propagation method for a dynamic lexicographic reduction */
static
SCIP_RETCODE propagateLexredDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array with (depth, index)-information per variable in
                                                  *   rooted path to focus node */
   int                   nvarstotal,         /**< length of nodedepthbranchindices array */
   SCIP_VAR**            branchvars,         /**< array populated with variables branched on */
   int                   nbranchvars,        /**< number of elements in branchvars array */
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
   assert( nodedepthbranchindices != NULL );
   assert( nvarstotal >= 0 );
   assert( branchvars != NULL );
   assert( nbranchvars >= 0 );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   nvars = lexdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &varorder, nvars) );

   SCIP_CALL( getVarOrder(scip, masterdata, lexdata, nodedepthbranchindices, nvarstotal, branchvars, nbranchvars,
      varorder, &nvarorder, nvars) );
   assert( nvarorder >= 0 );
   assert( nvarorder <= nvars );

   /* possibly propagate the constraint with this variable order */
   if ( nvarorder > 0 )
   {
      SCIP_CALL( propagateStaticLexred(scip, lexdata, varorder, nvarorder, infeasible, nreductions) );
   }
   SCIPfreeBufferArray(scip, &varorder);

   return SCIP_OKAY;
}


/** propagation method for applying dynamic lexicographic reduction for a single permutation */
static
SCIP_RETCODE propagateLexicographicReductionPerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array with (depth, index)-information per variable in
                                                  *   rooted path to focus node */
   int                   nvarstotal,         /**< length of that array */
   SCIP_VAR**            branchvars,         /**< array populated with variables branched on */
   int                   nbranchvars,        /**< number of elements in branchvars array */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( lexdata != NULL );
   assert( nodedepthbranchindices != NULL );
   assert( nvarstotal >= 0 );
   assert( branchvars != NULL );
   assert( nbranchvars >= 0 );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   *nreductions = 0;
   *infeasible = FALSE;

   SCIP_CALL( propagateLexredDynamic(scip, masterdata, lexdata,
      nodedepthbranchindices, nvarstotal, branchvars, nbranchvars, infeasible, nreductions) );

   return SCIP_OKAY;
}


/** populates array with information of first variable change
 *  @pre assuming nodedepthbranchindices is initially clean
 */
static
SCIP_RETCODE shadowtreeFillNodeDepthBranchIndices(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array to populate */
   SCIP_VAR**            branchvars,         /**< array to populate with variables branched on */
   int*                  nbranchvars,        /**< number of elements in branchvars array */
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
   assert( branchvars != NULL );
   assert( nbranchvars != NULL );
   assert( shadowtree != NULL );
   assert( focusnode != NULL );

   shadownode = SCIPshadowTreeGetShadowNode(shadowtree, focusnode);
   assert( shadownode != NULL );
   shadowdepth = SCIPnodeGetDepth(focusnode);

   /* branchvars array is initially empty */
   *nbranchvars = 0;

   /* We start looking one level lower, because we consider the branching decisions each time. */
   shadownode = shadownode->parent;

   /* Now, walk from the leaf to the root. Each time look at all the children of the node considered,
    * and save the variable depth and index in the branching array. It is important to consider all children each time,
    * because we need to comply with the instance where in different branches it is branched on different variables.
    * This has to be consistent.
    */
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
            assert( SCIPvarIsTransformed(var) );

            varindex = SCIPhashmapGetImageInt(masterdata->symvarmap, var);

            /* ignore variables that are irrelevant for lexicographic reduction */
            if ( varindex == INT_MAX )
               continue;

            assert( varindex >= 0 );
            assert( varindex < masterdata->nsymvars );

            /* var already in other child at this level? Continue */
            if ( nodedepthbranchindices[varindex].nodedepth == shadowdepth )
               continue;

            /* the variable is either not seen (nodedepth == 0), or it is at a higher level (nodedepth > shadowdepth) */
            assert( nodedepthbranchindices[varindex].nodedepth == 0 ||
               nodedepthbranchindices[varindex].nodedepth > shadowdepth);

            if ( nodedepthbranchindices[varindex].nodedepth == 0 )
            {
               /* variable is not featured in branchvars, yet */
               assert( *nbranchvars < masterdata->nsymvars );
               branchvars[(*nbranchvars)++] = var;
            }

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
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices, /**< array populated by nodedepthbranchindices to clean */
   SCIP_VAR**            branchvars,         /**< array populated with variables branched on */
   int*                  nbranchvars,        /**< number of elements in branchvars array */
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to shadow tree structure */
   SCIP_NODE*            focusnode           /**< focusnode to which the rooted path is evaluated */
)
{
   /* undo the operations from shadowtreeFillNodeDepthBranchIndices, which makes nodedepthbranchindices clean */
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
   assert( branchvars != NULL );
   assert( nbranchvars != NULL );
   assert( *nbranchvars >= 0 );
   assert( *nbranchvars <= masterdata->nsymvars );
   assert( shadowtree != NULL );
   assert( focusnode != NULL );

   shadownode = SCIPshadowTreeGetShadowNode(shadowtree, focusnode);
   shadowdepth = SCIPnodeGetDepth(focusnode);

   /* clean nbranchvars array */
   while ( *nbranchvars > 0 )
      branchvars[--(*nbranchvars)] = NULL;
   assert( *nbranchvars == 0 );

   /* we start looking one level lower, because we consider the branching decisions each time */
   shadownode = shadownode->parent;

   /* now, walk from the leaf to the root. Each time look at all the children of the node considered,
    * and save the variable depth and index in the branching array. It is important to consider all children each time,
    * because we need to comply with the instance where in different branches it is branched on different variables.
    * This has to be consistent.
    */
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
            /* ignore variables not relevant for lexicographic reduction */
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


/** prints lexicographic reduction propagation data */
SCIP_RETCODE SCIPlexicographicReductionGetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   int*                  nred                /**< total number of reductions applied */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( nred != NULL );

   *nred = masterdata->nred;

   return SCIP_OKAY;
}

/** prints lexicographic reduction propagation data */
SCIP_RETCODE SCIPlexicographicReductionPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata          /**< pointer to global data for lexicographic reduction propagator */
   )
{
   int i;

   assert( scip != NULL );
   assert( masterdata != NULL );

   if ( masterdata->nlexdatas == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   lexicographic reduction:   no permutations\n");
      return SCIP_OKAY;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   lexicographic reduction: %4d permutations with support sizes ",
      masterdata->nlexdatas);

   for (i = 0; i < masterdata->nlexdatas; ++i)
   {
      if ( i > 0 )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", ");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%d", masterdata->lexdatas[i]->nvars);
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");

   return SCIP_OKAY;
}


/** applies lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nred,               /**< pointer to store the number of domain reductions */
   SCIP_Bool*            didrun              /**< a global pointer maintaining if any symmetry propagator has run
                                              *   only set this to TRUE when a reduction is found, never set to FALSE */
   )
{
   int nlocalred;
   int p;
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices;
   SCIP_VAR** branchvars;
   int nbranchvars;
   SCIP_SHADOWTREE* shadowtree;
   SCIP_NODE* focusnode;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( (masterdata->lexdatas == NULL) == (masterdata->nlexdatas == 0) );
   assert( masterdata->nlexdatas >= 0 );
   assert( masterdata->nlexdatas <= masterdata->maxnlexdatas );
   assert( infeasible != NULL );
   assert( nred != NULL );
   assert( didrun != NULL );

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
    *
    * this is an array that maps every variable index to (depth, index) = (0, 0) if the variable is not branched on,
    * or (depth, index) is the depth (branching at root node: depth 1) and variable index when it was branched thereon.
    * For individual dynamic lexicographic reductions, we use this ordering the following way:
    *  1. Choose those variables that have (depth, index) with depth > 0 (these)
    *  2. Sort by depth, then by index, in increasing order.
    */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &nodedepthbranchindices, masterdata->nsymvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &branchvars, masterdata->nsymvars) );
   SCIP_CALL( shadowtreeFillNodeDepthBranchIndices(scip, masterdata, nodedepthbranchindices,
      branchvars, &nbranchvars, shadowtree, focusnode) );
   /* ... Do everything using this nodedepthbranchindices structure */

   if ( nbranchvars > 0 )
   {
      /* apply lexicographic reduction propagator to all lexdata objects */
      for (p = 0; p < masterdata->nlexdatas; ++p)
      {
         assert( masterdata->lexdatas[p] != NULL );

         SCIP_CALL( propagateLexicographicReductionPerm(scip, masterdata, masterdata->lexdatas[p],
            nodedepthbranchindices, masterdata->nsymvars, branchvars, nbranchvars, infeasible, &nlocalred) );

         /* keep track of the total number of fixed vars */
         *nred += nlocalred;

         /* a symmetry propagator has ran, so set didrun to TRUE */
         *didrun = TRUE;

         /* stop if we find infeasibility */
         if ( *infeasible )
            break;
      }

      /* maintain total number of reductions made */
      masterdata->nred += *nred;
   }

   /* clean the node-depth-branch-indices structure */
   SCIP_CALL( shadowtreeUndoNodeDepthBranchIndices(scip, masterdata, nodedepthbranchindices,
      branchvars, &nbranchvars, shadowtree, focusnode) );
   assert( nbranchvars == 0 );
   SCIPfreeCleanBufferArray(scip, &branchvars);
   SCIPfreeCleanBufferArray(scip, &nodedepthbranchindices);

   return SCIP_OKAY;
}


/** adds permutation for lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionAddPermutation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int*                  perm,               /**< permutation */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
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
   SCIP_CALL( lexdataCreate(scip, masterdata, &masterdata->lexdatas[masterdata->nlexdatas],
      permvars, npermvars, perm, success) );

   /* if not successfully added, undo increasing the counter */
   if ( *success )
      ++masterdata->nlexdatas;

   return SCIP_OKAY;
}


/** resets lexicographic reduction propagation (removes all permutations) */
SCIP_RETCODE SCIPlexicographicReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata          /**< pointer to global data for lexicographic reduction propagator */
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

   SCIPfreeBlockMemoryArrayNull(scip, &masterdata->lexdatas, masterdata->maxnlexdatas);
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


/** frees lexicographic reduction data */
SCIP_RETCODE SCIPlexicographicReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA**     masterdata          /**< pointer to global data for lexicographic reduction propagator */
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
 *
 *
 *  This is only done exactly once.
 */
SCIP_RETCODE SCIPincludeLexicographicReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA**     masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( shadowtreeeventhdlr != NULL );

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeLexicographicReduction", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPallocBlockMemory(scip, masterdata) );

   (*masterdata)->shadowtreeeventhdlr = shadowtreeeventhdlr;
   (*masterdata)->symvarmap = NULL;
   (*masterdata)->nsymvars = 0;
   (*masterdata)->lexdatas = NULL;
   (*masterdata)->nlexdatas = 0;
   (*masterdata)->maxnlexdatas = 0;
   (*masterdata)->nred = 0;

   return SCIP_OKAY;
}
