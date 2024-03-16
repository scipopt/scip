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

/**@file   symmetry_lexred.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries by dynamic lexicographic ordering reduction
 * @author Jasper van Doornmalen
 *
 * This implements lexicographic reduction, which generalizes symresack propagation to work for non-binary variable
 * domains, and is dynamified. Whereas static lexicographic reduction propagates that a vector \f$x\f$ is
 * lexicographically larger than its permuted counterpart (i.e., \f$x \succeq \gamma(x)\f$ with \f$\succeq\f$ being
 * standard elementwise lexicographic comparison), the dynamic variant determines the variable vector ordering
 * dynamically. Just as in orbital reduction (cf. symmetry_orbital.c), the variable order is chosen as the variables
 * branched on from the root node to the focus node.
 * Thus, in node \f$\beta\f$, we propagate \f$\sigma_\beta(x) \succeq \sigma_\beta(\gamma(x))\f$,
 * where \f$\sigma_\beta(\cdot)\f$ permutes and restricts the variable vector such that it corresponds to the branched
 * variables in the order from the rooted path to \f$\beta\f$.
 *
 * See Section 4.1 and Example 11 in [vD,H]:@n
 * J. van Doornmalen, C. Hojny, "A Unified Framework for Symmetry Handling", preprint, 2023,
 * https://doi.org/10.48550/arXiv.2211.01295.
 *
 * For dynamic lexicographic reduction, it is crucial that the vectors \f$\sigma_\beta(x)\f$ are the branching
 * variables up to node \f$\beta\f$ in the given order. Since SCIP can change the tree structure during solving
 * (re-writing history), we store the original branching decisions at the moment they are made. See event_shadowtree.c .
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
#include <string.h>


/*
 * Data structures
 */


/** data per permutation for lexicographic reduction propagator */
struct LexRedPermData
{
   SCIP_Bool             isdynamic;          /**< whether permutation shall be propagated with dynamic variable order */
   SCIP_VAR**            vars;               /**< variables affected by permutation */
   int                   nvars;              /**< number of variables */
   int*                  perm;               /**< permutation for lexicographic reduction */
   int*                  invperm;            /**< inverse permutation */
   SCIP_HASHMAP*         varmap;             /**< map of variables to indices in vars array */
   SYM_SYMTYPE           symtype;            /**< type of symmetries in perm */
   SCIP_Real*            vardomaincenter;    /**< array of centers of variable domains */
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
   int                   ncutoff;            /**< total number of cutoffs */
   SCIP_Bool             hasdynamicperm;     /**< whether there exists a permutation that is treated dynamically */
   SCIP_Bool             treewarninggiven;   /**< whether a warning is given for missing nodes in shadowtree */
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
   SCIP_Bool issigned;
   int permlen;
   int i;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( (*lexdata) != NULL );

   if ( (*lexdata)->symtype == SYM_SYMTYPE_SIGNPERM )
   {
      issigned = TRUE;
      permlen = 2 * (*lexdata)->nvars;
   }
   else
   {
      issigned = FALSE;
      permlen = (*lexdata)->nvars;
   }

   if ( (*lexdata)->nvars > 0 )
   {
      assert( (*lexdata)->invperm != NULL );
      assert( (*lexdata)->perm != NULL );
      assert( (*lexdata)->vars != NULL );

      /* free hashmap */
      if ( (*lexdata)->isdynamic )
      {
         assert( (*lexdata)->varmap != NULL );
         SCIPhashmapFree(&((*lexdata)->varmap));
      }

      /* release variables */
      for (i = 0; i < (*lexdata)->nvars; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*lexdata)->vars[i]) );
      }

      SCIPfreeBlockMemoryArray(scip, &(*lexdata)->invperm, permlen);
      SCIPfreeBlockMemoryArray(scip, &(*lexdata)->perm, permlen);
      SCIPfreeBlockMemoryArray(scip, &(*lexdata)->vars, (*lexdata)->nvars);

      if ( issigned )
      {
         SCIPfreeBlockMemoryArray(scip, &(*lexdata)->vardomaincenter, (*lexdata)->nvars);
      }
      else
      {
         assert( (*lexdata)->vardomaincenter == NULL );
      }
   }
#ifndef NDEBUG
   else
   {
      assert( (*lexdata)->nvars == 0 );
      assert( (*lexdata)->invperm == NULL );
      assert( (*lexdata)->vardomaincenter == NULL );
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
   SYM_SYMTYPE           symtype,            /**< type of symmetries in perm */
   SCIP_Real*            permvardomaincenter, /**< array containing center point for each variable domain */
   SCIP_Bool             usedynamicorder,    /**< whether a dynamic variable order shall be used */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   )
{
   SCIP_VAR* var;
   SCIP_Bool issignedperm;
   int* indexcorrection;
   int naffectedvariables;
   int permlen;
   int i;
   int j;

   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( lexdata != NULL );
   assert( vars != NULL );
   assert( nvars >= 0 );
   assert( perm != NULL );
   assert( symtype == SYM_SYMTYPE_PERM || permvardomaincenter != NULL );
   assert( success != NULL );
   assert( SCIPisTransformed(scip) );
   assert( masterdata->shadowtreeeventhdlr != NULL );

   *success = TRUE;
   issignedperm = symtype == SYM_SYMTYPE_PERM ? FALSE : TRUE;

   /* initialize the data structures */
   SCIP_CALL( SCIPallocBlockMemory(scip, lexdata) );
   (*lexdata)->symtype = symtype;

   (*lexdata)->isdynamic = usedynamicorder;

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

   /* require that the shadowtree is active if dynamic propagation is used */
   if ( usedynamicorder )
   {
      masterdata->hasdynamicperm = TRUE;

      SCIP_CALL( SCIPactivateShadowTree(scip, masterdata->shadowtreeeventhdlr) );
   }

   /* initialize variable arrays */
   (*lexdata)->nvars = naffectedvariables;
   permlen = issignedperm ? 2 * (*lexdata)->nvars : (*lexdata)->nvars;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->vars, (*lexdata)->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->perm, permlen) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->invperm, permlen) );
   if ( issignedperm )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexdata)->vardomaincenter, (*lexdata)->nvars) );
   }
   else
      (*lexdata)->vardomaincenter = NULL;

   /* determine the vars, perm, and centers */
   for (j = 0; j < nvars; ++j)
   {
      i = indexcorrection[j];
      if ( i < 0 )
         continue;

      /* j is the original index, i is the relabeled index */
      (*lexdata)->vars[i] = vars[j];

      if ( issignedperm )
      {
         if ( perm[j] >= nvars )
         {
            (*lexdata)->perm[i] = indexcorrection[perm[j] - nvars] + (*lexdata)->nvars;
            (*lexdata)->perm[i + (*lexdata)->nvars] = indexcorrection[perm[j] - nvars];
            assert( (*lexdata)->nvars <= (*lexdata)->perm[i] && (*lexdata)->perm[i] < 2 * (*lexdata)->nvars );
         }
         else
         {
            (*lexdata)->perm[i] = indexcorrection[perm[j]];
            (*lexdata)->perm[i + (*lexdata)->nvars] = indexcorrection[perm[j]] + (*lexdata)->nvars;
         }
      }
      else
         (*lexdata)->perm[i] = indexcorrection[perm[j]];

      if ( issignedperm )
         (*lexdata)->vardomaincenter[i] = permvardomaincenter[j];

      assert( perm[j] != j );
      assert( (*lexdata)->perm[i] != i );
      assert( (*lexdata)->perm[i] >= 0 );
      assert( (*lexdata)->perm[i] < permlen );
   }

   /* determine invperm */
   for (i = 0; i < (*lexdata)->nvars; ++i)
   {
      if ( (*lexdata)->perm[i] >= (*lexdata)->nvars )
      {
         assert( issignedperm );

         (*lexdata)->invperm[(*lexdata)->perm[i]] = i;
         (*lexdata)->invperm[(*lexdata)->perm[i] - (*lexdata)->nvars] = i + (*lexdata)->nvars;
      }
      else
      {
         (*lexdata)->invperm[(*lexdata)->perm[i]] = i;

         if ( issignedperm )
            (*lexdata)->invperm[(*lexdata)->perm[i] + (*lexdata)->nvars] = i + (*lexdata)->nvars;
      }
   }
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
   if ( usedynamicorder )
   {
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
   }
   else
      (*lexdata)->varmap = NULL;

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


/** return the index of a variable at a specific position of a variable order */
static
int varOrderGetIndex(
   int*                  varorder,           /**< variable order (or NULL) */
   int                   pos                 /**< position for which index is returned */
   )
{
   assert( pos >= 0 );

   if ( varorder == NULL )
      return pos;
   return varorder[pos];
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


/** gerts local variable bounds or reads bound from peek data */
static
SCIP_RETCODE getVarBounds(
   SCIP_VAR*             var1,               /**< first variable in comparison */
   SCIP_VAR*             var2,               /**< second variable in comparison */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx1,            /**< index of var1 (or NULL) */
   int                   varidx2,            /**< index of var2 (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset,          /**< whether peek bounds have been set (or NULL) */
   SCIP_Real*            lb1,                /**< pointer to store lower bound of var1 */
   SCIP_Real*            ub1,                /**< pointer to store upper bound of var1 */
   SCIP_Real*            lb2,                /**< pointer to store lower bound of var2 */
   SCIP_Real*            ub2                 /**< pointer to store upper bound of var2 */
   )
{
   assert( var1 != NULL );
   assert( var2 != NULL );
   assert( (!peekmode) || peeklbs != NULL );
   assert( (!peekmode) || peekubs != NULL );
   assert( (!peekmode) || peekbdset != NULL );
   assert( lb1 != NULL );
   assert( ub1 != NULL );
   assert( lb2 != NULL );
   assert( ub2 != NULL );

   if ( peekmode && peekbdset[varidx1] )
   {
      *ub1 = peekubs[varidx1];
      *lb1 = peeklbs[varidx1];
   }
   else
   {
      *ub1 = SCIPvarGetUbLocal(var1);
      *lb1 = SCIPvarGetLbLocal(var1);
   }
   if ( peekmode && peekbdset[varidx2] )
   {
      *ub2 = peekubs[varidx2];
      *lb2 = peeklbs[varidx2];
   }
   else
   {
      *ub2 = SCIPvarGetUbLocal(var2);
      *lb2 = SCIPvarGetLbLocal(var2);
   }

   return SCIP_OKAY;
}

/** returns whether a shifted variable is always smaller than another shifted variable
 *
 *  Shifts are always (var - shift).
 */
static
SCIP_Bool alwaysLTshiftedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable in comparison */
   SCIP_VAR*             var2,               /**< second variable in comparison */
   SCIP_Real             shift1,             /**< shift of var1 */
   SCIP_Real             shift2,             /**< shift of var2 */
   SCIP_Bool             isnegated,          /**< whether shift of var2 is negated */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx1,            /**< index of var1 (or NULL) */
   int                   varidx2,            /**< index of var2 (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   SCIP_Real ub1;
   SCIP_Real ub2;
   SCIP_Real lb1;
   SCIP_Real lb2;

   assert( scip != NULL );
   assert( var1 != NULL );
   assert( var2 != NULL );
   assert( (!peekmode) || peeklbs != NULL );
   assert( (!peekmode) || peekubs != NULL );
   assert( (!peekmode) || peekbdset != NULL );

   SCIP_CALL_ABORT( getVarBounds(var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
         &lb1, &ub1, &lb2, &ub2) );

   /* for negated variables, check: var1 - shift1 < shift2 - var2 */
   if ( isnegated && SCIPisLT(scip, ub1, shift1 + shift2 - ub2) )
      return TRUE;

   /* for non-negated variables, check: var1 - center1 < var2 - center2 */
   if ( (!isnegated) && SCIPisLT(scip, ub1, shift1 - shift2 + lb2) )
      return TRUE;

   return FALSE;
}


/** returns whether a shifted variable can be greater than another shifted variable
 *
 *  Shifts are always (var - shift).
 */
static
SCIP_Bool canGTshiftedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable in comparison */
   SCIP_VAR*             var2,               /**< second variable in comparison */
   SCIP_Real             shift1,             /**< shift of var1 */
   SCIP_Real             shift2,             /**< shift of var2 */
   SCIP_Bool             isnegated,          /**< whether shift of var2 is negated */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx1,            /**< index of var1 (or NULL) */
   int                   varidx2,            /**< index of var2 (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   SCIP_Real ub1;
   SCIP_Real ub2;
   SCIP_Real lb1;
   SCIP_Real lb2;

   assert( scip != NULL );
   assert( var1 != NULL );
   assert( var2 != NULL );
   assert( (!peekmode) || peeklbs != NULL );
   assert( (!peekmode) || peekubs != NULL );
   assert( (!peekmode) || peekbdset != NULL );

   SCIP_CALL_ABORT( getVarBounds(var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
         &lb1, &ub1, &lb2, &ub2) );

   /* for negated variables, check: var1 - shift1 > shift2 - var2 */
   if ( isnegated && SCIPisGT(scip, ub1, shift1 + shift2 - ub2) )
      return TRUE;

   /* for non-negated variables, check: var1 - center1 > var2 - center2 */
   if ( (!isnegated) && SCIPisGT(scip, ub1, shift1 - shift2 + lb2) )
      return TRUE;

   return FALSE;
}


/** propagates lower bound of first variable in relation x >= y for shifted variables */
static
SCIP_RETCODE propagateLowerBoundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable in pair */
   SCIP_VAR*             var2,               /**< second variable in pair */
   SCIP_Real             center1,            /**< center of var1 (original var domain) */
   SCIP_Real             center2,            /**< center of var2 (original var domain) */
   SCIP_Bool             isnegated,          /**< whether var2 is negated by symmetry */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nreductions,        /**< pointer to store number of reductions */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx1,            /**< index of var1 (or NULL) */
   int                   varidx2,            /**< index of var2 (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   SCIP_Real ub1;
   SCIP_Real ub2;
   SCIP_Real lb1;
   SCIP_Real lb2;

   SCIP_Bool tighten = FALSE;
   SCIP_Real newbd;

   assert( (!peekmode) || peeklbs != NULL );
   assert( (!peekmode) || peekubs != NULL );
   assert( (!peekmode) || peekbdset != NULL );

   SCIP_CALL( getVarBounds(var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
         &lb1, &ub1, &lb2, &ub2) );

   /* tighten domain of var1 to ensure that var1 - center1 >= isnegated * (var2 - center2 ) */
   if ( isnegated )
   {
      if ( SCIPisLT(scip, lb1 - center1, center2 - ub2) )
      {
         tighten = TRUE;
         newbd = center1 + center2 - ub2;
      }
   }
   else
   {
      if ( SCIPisLT(scip, lb1 - center1, lb2 - center2) )
      {
         tighten = TRUE;
         newbd = center1 + lb2 - center2;
      }
   }

   if ( tighten )
   {
      /* in peek mode, only store updated bounds */
      if ( peekmode )
      {
         peeklbs[varidx1] = newbd; /*lint !e644*/
         peekubs[varidx1] = ub1;
         peekbdset[varidx1] = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPtightenVarLb(scip, var1, newbd, TRUE, infeasible, &tighten) );
         if ( tighten )
         {
            SCIPdebugMessage("Restricting variable LB %12s to %5.2f\n", SCIPvarGetName(var1), newbd);
            *nreductions += 1;
         }
         else
         {
            SCIPdebugMessage("Restricting variable LB %12s to %5.2f (no success)\n",
               SCIPvarGetName(var1), newbd);
         }
         if ( *infeasible )
         {
            SCIPdebugMessage("Detected infeasibility restricting variable LB %12s to %5.2f\n",
               SCIPvarGetName(var1), newbd);
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** propagates upper bound of second variable in relation x >= y for shifted variables */
static
SCIP_RETCODE propagateUpperBoundSymVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable in pair */
   SCIP_VAR*             var2,               /**< second variable in pair */
   SCIP_Real             center1,            /**< center of var1 (original var domain) */
   SCIP_Real             center2,            /**< center of var2 (original var domain) */
   SCIP_Bool             isnegated,          /**< whether var2 is negated by symmetry */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nreductions,        /**< pointer to store number of reductions */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx1,            /**< index of var1 (or NULL) */
   int                   varidx2,            /**< index of var2 (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   SCIP_Real ub1;
   SCIP_Real ub2;
   SCIP_Real lb1;
   SCIP_Real lb2;

   SCIP_Bool tighten = FALSE;
   SCIP_Real newbd;

   assert( scip != NULL );
   assert( var1 != NULL );
   assert( var2 != NULL );
   assert( infeasible != NULL );
   assert( nreductions != NULL );
   assert( (!peekmode) || peeklbs != NULL );
   assert( (!peekmode) || peekubs != NULL );
   assert( (!peekmode) || peekbdset != NULL );

   SCIP_CALL( getVarBounds(var1, var2, peekmode, varidx1, varidx2, peeklbs, peekubs, peekbdset,
         &lb1, &ub1, &lb2, &ub2) );

   /* tighten domain of var2 to ensure that var1 - center1 >= isnegated * (var2 - center2 ) */
   if ( isnegated )
   {
      if ( SCIPisLT(scip, ub1 - center1, center2 - lb2) )
      {
         tighten = TRUE;
         newbd = center1 + center2 - ub1;
      }
   }
   else
   {
      if ( SCIPisLT(scip, ub1 - center1, ub2 - center2) )
      {
         tighten = TRUE;
         newbd = center2 - center1 + ub1;
      }
   }

   if ( tighten )
   {
      /* in peek mode, only store updated bounds */
      if ( peekmode )
      {
         if ( isnegated )
         {
            peeklbs[varidx2] = newbd; /*lint !e644*/
            peekubs[varidx2] = ub2;
            peekbdset[varidx2] = TRUE;
         }
         else
         {
            peeklbs[varidx2] = lb2;
            peekubs[varidx2] = newbd;
            peekbdset[varidx2] = TRUE;
         }
      }
      else
      {
         if ( isnegated )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var2, newbd, TRUE, infeasible, &tighten) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var2, newbd, TRUE, infeasible, &tighten) );
         }

         if ( tighten )
         {
            SCIPdebugMessage("Restricting variable %sB %12s to %5.2f\n",
               isnegated ? "L" : "U", SCIPvarGetName(var2), newbd);
            *nreductions += 1;
         }
         else
         {
            SCIPdebugMessage("Restricting variable %sB %12s to %5.2f (no success)\n",
               isnegated ? "L" : "U", SCIPvarGetName(var2), newbd);
         }
         if ( *infeasible )
         {
            SCIPdebugMessage("Detected infeasibility restricting variable %sB %12s to %5.2f\n",
               isnegated ? "L" : "U", SCIPvarGetName(var2), newbd);
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** propagates x - c >=  c - x */
static
SCIP_RETCODE propagateSelfReflectionVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             center,             /**< center of var (original var domain) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nreductions,        /**< pointer to store number of reductions */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx,             /**< index of var (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   SCIP_Real ub1;
   SCIP_Real ub2;
   SCIP_Real lb1;
   SCIP_Real lb2;
   SCIP_Bool tighten = FALSE;

   assert( scip != NULL );
   assert( var != NULL );
   assert( infeasible != NULL );
   assert( nreductions != NULL );
   assert( (!peekmode) || peeklbs != NULL );
   assert( (!peekmode) || peekubs != NULL );
   assert( (!peekmode) || peekbdset != NULL );

   SCIP_CALL( getVarBounds(var, var, peekmode, varidx, varidx, peeklbs, peekubs, peekbdset,
         &lb1, &ub1, &lb2, &ub2) );

   if ( SCIPisLT(scip, ub1, center) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   else if ( SCIPisLT(scip, lb1, center) )
   {
      if ( peekmode )
      {
         peeklbs[varidx] = center;
         peekubs[varidx] = ub1;
         peekbdset[varidx] = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPtightenVarLb(scip, var, center, TRUE, infeasible, &tighten) );
         if ( tighten )
         {
            SCIPdebugMessage("Restricting variable LB %12s to %5.2f\n", SCIPvarGetName(var), center);
            *nreductions += 1;
         }
         else
         {
            SCIPdebugMessage("Restricting variable LB %12s to %5.2f (no success)\n",
               SCIPvarGetName(var), center);
         }
         if ( *infeasible )
         {
            SCIPdebugMessage("Detected infeasibility restricting variable LB %12s to %5.2f\n",
               SCIPvarGetName(var), center);
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** propagates lexicographic order for one pair of symmetric variables */
static
SCIP_RETCODE propagateVariablePair(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable in pair */
   SCIP_VAR*             var2,               /**< second variable in pair */
   SCIP_Real             center1,            /**< center of var1 (original var domain) */
   SCIP_Real             center2,            /**< center of var2 (original var domain) */
   SCIP_Bool             isnegated,          /**< whether var2 is negated by symmetry */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nreductions,        /**< pointer to store number of reductions */
   SCIP_Bool             peekmode,           /**< whether function is called in peek mode */
   int                   varidx1,            /**< index of var1 (or NULL) */
   int                   varidx2,            /**< index of var2 (or NULL) */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   assert( scip != NULL );
   assert( var1 != NULL );
   assert( var2 != NULL );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   /* perform lexicographic comparison: var1 - center1 >= +/- (var2 - center2)  */
   if ( alwaysLTshiftedVars(scip, var1, var2, center1, center2, isnegated, peekmode,
         varidx1, varidx2, peeklbs, peekubs, peekbdset) )
   {
#ifdef SCIP_DEBUG
      SCIP_Real ub1;
      SCIP_Real ub2;
      SCIP_Real lb1;
      SCIP_Real lb2;

      /* get bounds of shifted (and possibly negated) variables */
      ub1 = SCIPvarGetUbLocal(var1);
      lb1 = SCIPvarGetLbLocal(var1);
      ub2 = SCIPvarGetUbLocal(var2);
      lb2 = SCIPvarGetLbLocal(var2);

      SCIPdebugMessage("Detected infeasibility: x < y for "
         "x: lb=%5.2f, ub=%5.2f, shift=%5.2f "
         "y: lb=%5.2f, ub=%5.2f, shift=%5.2f negated=%u\n",
         lb1, ub1, center1, lb2, ub2, center2, isnegated);
#endif

      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* for signed permutations, a variable might be mapped to itself */
   if ( var1 == var2 )
   {
      SCIP_CALL( propagateSelfReflectionVar(scip, var1, center1, infeasible, nreductions, peekmode, varidx1,
            peeklbs, peekubs, peekbdset) );
   }
   else
   {
      SCIP_CALL( propagateLowerBoundVar(scip, var1, var2, center1, center2, isnegated, infeasible, nreductions, peekmode,
            varidx1, varidx2, peeklbs, peekubs, peekbdset) );
      if ( *infeasible )
         return SCIP_OKAY;

      SCIP_CALL( propagateUpperBoundSymVar(scip, var1, var2, center1, center2, isnegated, infeasible, nreductions, peekmode,
            varidx1, varidx2, peeklbs, peekubs, peekbdset) );
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
   int*                  varorder,           /**< array populated with variable order (or NULL for static propagation) */
   int                   nselvars,           /**< number of variables in the ordering */
   int                   fixi,               /**< variable index of left fixed column */
   int                   fixj,               /**< variable index of right fixed column */
   int                   fixrow,             /**< row index in var ordering, from which to start peeking */
   SCIP_Real             fixvaluei,          /**< value on which variables i is fixed */
   SCIP_Real             fixvaluej,          /**< value on which variables j is fixed */
   SCIP_Bool*            peekfeasible,       /**< pointer to store whether this is feasible or not */
   SCIP_Real*            peeklbs,            /**< lower bounds of variables in peek mode (or NULL) */
   SCIP_Real*            peekubs,            /**< upper bounds of variables in peek mode (or NULL) */
   SCIP_Bool*            peekbdset           /**< whether peek bounds have been set (or NULL) */
   )
{
   int row;
   int i;
   int j;
   SCIP_VAR* vari;
   SCIP_VAR* varj;
   SCIP_Real centeri;
   SCIP_Real centerj;
   SCIP_Bool issigned;
   SCIP_Bool isnegated;
   SCIP_Bool infeasible = FALSE;
   int nreductions;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( lexdata->vars != NULL );
   assert( lexdata->nvars >= 0 );
   assert( nselvars <= lexdata->nvars );
   assert( nselvars > 0 );
   assert( fixi >= 0 );
   assert( fixi < lexdata->nvars );
   assert( fixj < lexdata->nvars );
   assert( fixi != fixj || lexdata->symtype == SYM_SYMTYPE_SIGNPERM );
   assert( fixi != fixj || fixvaluei == fixvaluej ); /*lint !e777*/
   assert( fixrow >= 0 );
   assert( fixrow < nselvars );
   assert( peekfeasible != NULL );
   assert( fixi == varOrderGetIndex(varorder, fixrow) );
   assert( fixj == (lexdata->invperm[varOrderGetIndex(varorder, fixrow)] % lexdata->nvars) );
   assert( fixi == (lexdata->perm[fixj] % lexdata->nvars) );

   *peekfeasible = TRUE;
   issigned = lexdata->symtype == SYM_SYMTYPE_SIGNPERM ? TRUE : FALSE;
   assert( (!issigned) || lexdata->vardomaincenter != NULL );

   /* intialize peekbdset */
   for (i = 0; i < lexdata->nvars; ++i)
      peekbdset[i] = FALSE;

   peeklbs[fixi] = fixvaluei;
   peeklbs[fixj] = fixvaluej;
   peekubs[fixi] = fixvaluei;
   peekubs[fixj] = fixvaluej;
   peekbdset[fixi] = TRUE;
   peekbdset[fixj] = TRUE;

   for (row = fixrow + 1; row < nselvars; ++row)
   {
      /* get left and right column indices */
      i = varOrderGetIndex(varorder, row);
      j = lexdata->invperm[i];
      assert( i == lexdata->perm[j] );

      /* no fixed points */
      assert( i != j );

      assert( 0 <= i && i < lexdata->nvars );
      assert( 0 <= j && j < 2 * lexdata->nvars );
      assert( issigned || j < lexdata->nvars );

      vari = lexdata->vars[i];
      if ( j >= lexdata->nvars )
      {
         assert( lexdata->symtype == SYM_SYMTYPE_SIGNPERM );
         j = j - lexdata->nvars;
         varj = lexdata->vars[j];
         isnegated = TRUE;
      }
      else
      {
         varj = lexdata->vars[j];
         isnegated = FALSE;
      }

      if ( issigned )
      {
         assert( lexdata->vardomaincenter != NULL );
         centeri = lexdata->vardomaincenter[i];
         centerj = lexdata->vardomaincenter[j];
      }
      else
      {
         centeri = 0.0;
         centerj = 0.0;
      }

      /* propagate that vari >= varj */

      /* vari >= varj can never hold if the maximal value of vari is smaller than the minimal value of varj */
      if ( alwaysLTshiftedVars(scip, vari, varj, centeri, centerj, isnegated, TRUE, i, j, peeklbs, peekubs, peekbdset) )
      {
         *peekfeasible = FALSE;
         SCIPdebugMessage("PEEK: Detected infeasibility at row %3d.\n", row);
         break;
      }

      SCIP_CALL( propagateLowerBoundVar(scip, vari, varj, centeri, centerj, isnegated, &infeasible, &nreductions, TRUE,
            i, j, peeklbs, peekubs, peekbdset) );

      SCIP_CALL( propagateUpperBoundSymVar(scip, vari, varj, centeri, centerj, isnegated, &infeasible, &nreductions, TRUE,
            i, j, peeklbs, peekubs, peekbdset) );

      /* if there exists a solution with vari > varj, reductions are feasible w.r.t. lexred */
      if ( canGTshiftedVars(scip, vari, varj, centeri, centerj, isnegated, TRUE,
            i, j, peeklbs, peekubs, peekbdset) )
         break;
   }

   return SCIP_OKAY;
}

/** propagates static lexicographic reduction with specified variable ordering */
static
SCIP_RETCODE propagateStaticLexred(
   SCIP*                 scip,               /**< SCIP data structure */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   int*                  varorder,           /**< array populated with variable order (or NULL if static propagation) */
   int                   nselvars,           /**< number of variables in the ordering */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
   )
{ /*lint --e{771}*/
   int row;
   int i = -1;
   int j = -1;
   SCIP_VAR* vari = NULL;
   SCIP_VAR* varj = NULL;
   SCIP_Real centeri;
   SCIP_Real centerj;
   SCIP_Bool success;
   SCIP_Bool issigned;
   SCIP_Bool isnegated;

   assert( scip != NULL );
   assert( lexdata != NULL );
   assert( nselvars >= 0 );
   assert( infeasible != NULL );
   assert( !*infeasible );
   assert( nreductions != NULL );
   assert( *nreductions >= 0 );

   /* avoid trivial cases */
   if ( nselvars <= 0 )
      return SCIP_OKAY;

   issigned = lexdata->symtype == SYM_SYMTYPE_SIGNPERM ? TRUE : FALSE;
   assert( (!issigned) || lexdata->vardomaincenter != NULL );

   /* iterate over the variable array entrywise
    *
    * We see this as two columns, with the left vector being the variable ordering,
    * and the right column the permuted variables of this var ordering.
    */
   for (row = 0; row < nselvars; ++row)
   {
      /* left and right column indices */
      i = varOrderGetIndex(varorder, row);
      j = lexdata->invperm[i];
      assert( i == lexdata->perm[j] );

      /* no fixed points */
      assert( i != j );

      assert( 0 <= i && i < lexdata->nvars );
      assert( 0 <= j && j < 2 * lexdata->nvars );
      assert( issigned || j < lexdata->nvars );

      vari = lexdata->vars[i];
      if ( j >= lexdata->nvars )
      {
         assert( issigned );
         j = j - lexdata->nvars;
         varj = lexdata->vars[j];
         isnegated = TRUE;
      }
      else
      {
         varj = lexdata->vars[j];
         isnegated = FALSE;
      }

      if ( issigned )
      {
         assert( lexdata->vardomaincenter != NULL );
         centeri = lexdata->vardomaincenter[i];
         centerj = lexdata->vardomaincenter[j];
      }
      else
      {
         centeri = 0.0;
         centerj = 0.0;
      }

      SCIP_CALL( propagateVariablePair(scip, vari, varj, centeri, centerj, isnegated, infeasible, nreductions, FALSE,
         0, 0, NULL, NULL, NULL) );

      if ( *infeasible )
         return SCIP_OKAY;

      /* terminate if there exists a solution being lexicographically strictly larger than its image */
      if ( canGTshiftedVars(scip, vari, varj, centeri, centerj, isnegated, FALSE,
            0, 0, NULL, NULL, NULL) )
         break;
   }
   assert( vari != NULL );
   assert( varj != NULL );
   assert( 0 <= i && i < lexdata->nvars );
   assert( 0 <= j && j < lexdata->nvars );

   /* if the variables in "row" are fixed to the same value, we might find further propagations */
   if ( row < nselvars )
   {
      SCIP_Real* peeklbs;
      SCIP_Real* peekubs;
      SCIP_Bool* peekbdset;
      SCIP_Real ub1;
      SCIP_Real ub2;
      SCIP_Real lb1;
      SCIP_Real lb2;
      SCIP_Real lbi;
      SCIP_Real ubi;
      SCIP_Real lbj;
      SCIP_Real ubj;
      SCIP_Bool peekfeasible;

      SCIP_CALL( getVarBounds(vari, varj, FALSE, 0, 0, NULL, NULL, NULL, &lb1, &ub1, &lb2, &ub2) );

      /* special treatment if i-th and j-th variable are the same in a signed permutation */
      if ( vari == varj )
      {
         assert( lexdata->symtype == SYM_SYMTYPE_SIGNPERM );
         assert( SCIPsymGE(scip, lb1, lexdata->vardomaincenter[i]) ); /* propagation enforces xi - center >= center - xi */

         /* both variables can only be the same if they are fixed to the domain center */
         if ( SCIPsymGT(scip, lb1, lexdata->vardomaincenter[i]) )
            return SCIP_OKAY;

         SCIP_CALL( SCIPallocBufferArray(scip, &peeklbs, lexdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &peekubs, lexdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &peekbdset, lexdata->nvars) );

         SCIP_CALL( peekStaticLexredIsFeasible(scip, lexdata, varorder, nselvars, i, j,
               row, lexdata->vardomaincenter[i], lexdata->vardomaincenter[i],
               &peekfeasible, peeklbs, peekubs, peekbdset) );
         if ( !peekfeasible )
         {
            /* both variables cannot be the same, so the non-negated variable must be greater than the domain center */
            switch ( SCIPvarGetType(vari) )
            {
            case SCIP_VARTYPE_BINARY:
            case SCIP_VARTYPE_IMPLINT:
            case SCIP_VARTYPE_INTEGER:
               assert( SCIPisIntegral(scip, lb1) );
               SCIP_CALL( SCIPtightenVarLb(scip, vari, lexdata->vardomaincenter[i] + 1.0, TRUE, infeasible, &success) );
               if ( success )
                  *nreductions += 1;
               if ( *infeasible )
                  goto FREEMEMORY;
               lb1 = lexdata->vardomaincenter[i] + 1.0;
               assert( SCIPsymLE(scip, lb1, ub1) );
               break;
            case SCIP_VARTYPE_CONTINUOUS:
               /* continuous variable type: act as if we increase the variable by a very little bit.
                * This is only possible if we're able to increase the variable bound by a bit.
                */
               if ( SCIPsymEQ(scip, lb1, ub1) )
               {
                  *infeasible = TRUE;
                  goto FREEMEMORY;
               }
               break;
            default:
               SCIPerrorMessage("unsupported variable type encountered at the lexicographic reduction propagator\n");
               return SCIP_ERROR;
            }
         }
      }
      else
      {
         /* The previous loop is broken at row "row", which allows for choosing vari > varj.
          *
          * Now check if vari == varj is permitted, and if not, tighten the domain further.
          *
          * @todo We peek twice if vari and varj are unfixed
          * But, if the subcycle only contains var1 and var2, a single peek suffices.
          * This is similar to orbisack and symresack propagation.
          *
          * Case 1: vari is minimal (lbi).
          * Then, propagation of lbi = vari >= varj can yield two situations:
          *   Option 1: varj can take a value < lbi. Then no further reductions can be detected.
          *   Option 2: varj gets fixed to lbi. Then, we must check if feasibility is found, still.
          *     If it turns out infeasible, then we know vari cannot take value lbi, so we can increase the lower bound.
          */
         centeri = 0.0;
         centerj = 0.0;

         if ( lexdata->vardomaincenter != NULL )
         {
            centeri = lexdata->vardomaincenter[i];
            centerj = lexdata->vardomaincenter[j];
         }

         /* translate variable bounds to shifted variable domain and take negation into account */
         lbi = lb1 - centeri;
         ubi = ub1 - centeri;
         if ( isnegated )
         {
            lbj = centerj - ub2;
            ubj = centerj - lb2;
         }
         else
         {
            lbj = lb2 - centerj;
            ubj = ub2 - centerj;
         }

         /* check whether peek is called */
         if ( (!SCIPsymEQ(scip, lbi, lbj)) && (!SCIPsymEQ(scip, ubi, ubj)) )
            return SCIP_OKAY;

         SCIP_CALL( SCIPallocBufferArray(scip, &peeklbs, lexdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &peekubs, lexdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &peekbdset, lexdata->nvars) );

         if ( SCIPsymEQ(scip, lbj, lbi) )
         {
            SCIP_Real fixvalj;

            /* translate lbj back to original variable domain of variable j */
            if ( isnegated )
               fixvalj = centerj - lbj;
            else
               fixvalj = lbj + centerj;

            /* this is Option 2: varj gets fixed to lbi by propagation. */
            SCIP_CALL( peekStaticLexredIsFeasible(scip, lexdata, varorder, nselvars, i, j,
                  row, lbi + centeri, fixvalj, &peekfeasible, peeklbs, peekubs, peekbdset) );
            if ( !peekfeasible )
            {
               /* vari cannot take value lb1, so we increase the lower bound. (do not use lbi as this is a shifted domain bound) */
               switch ( SCIPvarGetType(vari) )
               {
               case SCIP_VARTYPE_BINARY:
               case SCIP_VARTYPE_IMPLINT:
               case SCIP_VARTYPE_INTEGER:
                  /* discrete variable type: increase lower bound by 1. */
                  assert( SCIPisIntegral(scip, lb1) );
                  SCIP_CALL( SCIPtightenVarLb(scip, vari, lb1 + 1.0, TRUE, infeasible, &success) );
                  if ( success )
                     *nreductions += 1;
                  if ( *infeasible )
                     goto FREEMEMORY;
                  lb1 = lb1 + 1.0;
                  assert( SCIPsymLE(scip, lb1, ub1) );
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous variable type: act as if we increase the variable by a very little bit.
                   * That is only possible if we're able to increase the variable bound by a bit.
                   */
                  if ( SCIPsymEQ(scip, lbi, ubi) )
                  {
                     *infeasible = TRUE;
                     goto FREEMEMORY;
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
         assert( SCIPsymGE(scip, ubi, ubj) );  /* this must be the case after reductions in the for-loop */
         if ( SCIPsymEQ(scip, ubi, ubj) )
         {
            SCIP_Real fixvalj;

            /* translate ubj back to original variable domain of variable j */
            if ( isnegated )
               fixvalj = centerj - ubj;
            else
               fixvalj = ubj + centerj;

            /* this is Option 2: vari gets fixed to ubj by propagation. */
            SCIP_CALL( peekStaticLexredIsFeasible(scip, lexdata, varorder, nselvars, i, j,
                  row, ubi + centeri, fixvalj, &peekfeasible, peeklbs, peekubs, peekbdset) );
            if ( !peekfeasible )
            {
               /* varj cannot take value ub2, so we decrease the upper or lower bound. (do not use ubj as this is a shifted domain bound)*/
               switch ( SCIPvarGetType(varj) )
               {
               case SCIP_VARTYPE_BINARY:
               case SCIP_VARTYPE_IMPLINT:
               case SCIP_VARTYPE_INTEGER:
                  /* discrete variable type: decrease upper bound by 1. */
                  if ( isnegated )
                  {
                     assert( SCIPisIntegral(scip, lb2) );
                     SCIP_CALL( SCIPtightenVarUb(scip, varj, lb2 + 1.0, TRUE, infeasible, &success) );
                  }
                  else
                  {
                     assert( SCIPisIntegral(scip, ub2) );
                     SCIP_CALL( SCIPtightenVarUb(scip, varj, ub2 - 1.0, TRUE, infeasible, &success) );
                  }
                  if ( success )
                     *nreductions += 1;
                  if ( *infeasible )
                     goto FREEMEMORY;
                  ubj = ubj - 1.0;
                  assert( SCIPsymLE(scip, lbj, ubj) );
                  break;
               case SCIP_VARTYPE_CONTINUOUS:
                  /* continuous variable type: act as if we decrease the variable by a very little bit.
                   * that is only possible if we're able to decrease the variable bound by a bit. */
                  if ( SCIPsymEQ(scip, lbj, ubj) )
                  {
                     *infeasible = TRUE;
                     goto FREEMEMORY;
                  }
                  break;
               default:
                  return SCIP_ERROR;
               }
            }
         }
      }
   FREEMEMORY:
      SCIPfreeBufferArray(scip, &peekbdset);
      SCIPfreeBufferArray(scip, &peekubs);
      SCIPfreeBufferArray(scip, &peeklbs);
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
   assert( lexdata->isdynamic );
   assert( nodedepthbranchindices != NULL );
   assert( nvarstotal >= 0 );
   assert( branchvars != NULL );
   assert( nbranchvars >= 0 );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   nvars = lexdata->nvars;

   /* get variable order */
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


/** propagation method for a dynamic lexicographic reduction */
static
SCIP_RETCODE propagateLexredStatic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   LEXDATA*              lexdata,            /**< pointer to data for this permutation */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the problem is infeasible */
   int*                  nreductions         /**< pointer to store the number of found domain reductions */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( lexdata != NULL );
   assert( ! lexdata->isdynamic );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   /* skip trivial cases */
   if ( lexdata->nvars == 0 )
      return SCIP_OKAY;

   /* propagate the constraint with this variable order */
   SCIP_CALL( propagateStaticLexred(scip, lexdata, NULL, lexdata->nvars, infeasible, nreductions) );

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
   assert( nodedepthbranchindices != NULL || ! lexdata->isdynamic );
   assert( nvarstotal >= 0 || ! lexdata->isdynamic );
   assert( branchvars != NULL || ! lexdata->isdynamic );
   assert( nbranchvars >= 0 || ! lexdata->isdynamic );
   assert( infeasible != NULL );
   assert( nreductions != NULL );

   *nreductions = 0;
   *infeasible = FALSE;

   if ( lexdata->isdynamic )
   {
      SCIP_CALL( propagateLexredDynamic(scip, masterdata, lexdata,
            nodedepthbranchindices, nvarstotal, branchvars, nbranchvars, infeasible, nreductions) );
   }
   else
   {
      SCIP_CALL( propagateLexredStatic(scip, masterdata, lexdata, infeasible, nreductions) );
   }

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
   SCIP_NODE*            focusnode,          /**< focusnode to which the rooted path is evaluated */
   SCIP_Bool*            inforesolved        /**< pointer to store whether information is successfully resolved */
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
   assert( inforesolved != NULL );

   shadownode = SCIPshadowTreeGetShadowNode(shadowtree, focusnode);

   if ( shadownode == NULL )
   {
      /* the arrays to fill are unchanged, so they remain clean */
      *inforesolved = FALSE;
      if ( !masterdata->treewarninggiven )
      {
         SCIPwarningMessage(scip, "Attempting lexicographic reduction on nodes not existing in the symmetry shadowtree"
            " (and suppressing future warnings)\n");
         masterdata->treewarninggiven = TRUE;
      }
      return SCIP_OKAY;
   }
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

   *inforesolved = TRUE;
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
#ifndef NDEBUG
   int shadowdepth;
#endif
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
#ifndef NDEBUG
   shadowdepth = SCIPnodeGetDepth(focusnode);
#endif

   /* clean nbranchvars array */
   while ( *nbranchvars > 0 )
      branchvars[--(*nbranchvars)] = NULL;
   assert( *nbranchvars == 0 );

   /* we start looking one level lower, because we consider the branching decisions each time */
   /* coverity[dereference] */
   shadownode = shadownode->parent;

   /* now, walk from the leaf to the root. Each time look at all the children of the node considered,
    * and save the variable depth and index in the branching array. It is important to consider all children each time,
    * because we need to comply with the instance where in different branches it is branched on different variables.
    * This has to be consistent.
    */
   while (shadownode != NULL)
   {
#ifndef NDEBUG
      assert( shadowdepth > 0 );
#endif
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
#ifndef NDEBUG
      --shadowdepth;
#endif
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
   int*                  nred,               /**< total number of reductions applied */
   int*                  ncutoff             /**< total number of cutoffs applied */
   )
{
   assert( scip != NULL );
   assert( masterdata != NULL );
   assert( nred != NULL );

   *nred = masterdata->nred;
   *ncutoff = masterdata->ncutoff;

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
   SCIP_SHADOWTREE* shadowtree = NULL;
   SCIP_NODE* focusnode = NULL;
   NODEDEPTHBRANCHINDEX* nodedepthbranchindices = NULL;
   SCIP_VAR** branchvars = NULL;
   int nbranchvars = 0;
   SCIP_Bool inforesolved;

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

   /* compute the variable ordering based on the branching decisions
    * of the shadowtree if there exist dynamic permutations
    */
   if ( masterdata->hasdynamicperm )
   {
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
            branchvars, &nbranchvars, shadowtree, focusnode, &inforesolved) );

      /* if the information cannot be resolved because a node is missing from the shadowtree, do not propagate */
      if ( !inforesolved )
      {
         /* shadowtreeFillNodeDepthBranchIndices keeps the input arrays clean if it terminates early */
         SCIPfreeCleanBufferArray(scip, &branchvars);
         SCIPfreeCleanBufferArray(scip, &nodedepthbranchindices);
         return SCIP_OKAY;
      }
      /* ... Do everything using this nodedepthbranchindices structure */
   }

   if ( nbranchvars > 0 || ! masterdata->hasdynamicperm )
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
      if ( *infeasible )
         ++masterdata->ncutoff;
   }

   /* possibly clean the node-depth-branch-indices structure */
   if ( masterdata->hasdynamicperm )
   {
      assert( shadowtree != NULL );
      assert( focusnode != NULL );
      SCIP_CALL( shadowtreeUndoNodeDepthBranchIndices(scip, masterdata, nodedepthbranchindices,
            branchvars, &nbranchvars, shadowtree, focusnode) );
      assert( nbranchvars == 0 );
      SCIPfreeCleanBufferArray(scip, &branchvars);
      SCIPfreeCleanBufferArray(scip, &nodedepthbranchindices);
   }

   return SCIP_OKAY;
}


/** adds permutation for lexicographic reduction propagation */
SCIP_RETCODE SCIPlexicographicReductionAddPermutation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int*                  perm,               /**< permutation */
   SYM_SYMTYPE           symtype,            /**< type of symmetries in perm */
   SCIP_Real*            permvardomaincenter, /**< array containing center point for each variable domain */
   SCIP_Bool             usedynamicorder,    /**< whether a dynamic variable order shall be used */
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

   if ( symtype != SYM_SYMTYPE_PERM && symtype != SYM_SYMTYPE_SIGNPERM )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   assert( symtype == SYM_SYMTYPE_PERM || permvardomaincenter != NULL );
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
         permvars, npermvars, perm, symtype, permvardomaincenter, usedynamicorder, success) );

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

   masterdata->hasdynamicperm = FALSE;

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

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeLexicographicReduction", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPallocBlockMemory(scip, masterdata) );

   (*masterdata)->shadowtreeeventhdlr = shadowtreeeventhdlr;
   (*masterdata)->symvarmap = NULL;
   (*masterdata)->nsymvars = 0;
   (*masterdata)->lexdatas = NULL;
   (*masterdata)->nlexdatas = 0;
   (*masterdata)->maxnlexdatas = 0;
   (*masterdata)->nred = 0;
   (*masterdata)->ncutoff = 0;
   (*masterdata)->hasdynamicperm = FALSE;
   (*masterdata)->treewarninggiven = FALSE;

   return SCIP_OKAY;
}
