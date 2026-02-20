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

/**@file   sym_rowcol.c
 * @ingroup DEFPLUGINS_SYM
 * @brief  symmetry handler for row and column symmetries
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_linear.h"
#include "scip/cons_orbitope.h"
#include "scip/event_shadowtree.h"
#include "scip/pub_cons.h"
#include "scip/pub_implics.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sym.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_datastructures.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sym.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/sym_rowcol.h"
#include "scip/symmetry.h"
#include "scip/symmetry_lexred.h"
#include "scip/symmetry_orbitopal.h"
#include "scip/type_implics.h"
#include "scip/type_set.h"

/* symmetry handler properties */
#define SYM_NAME              "rowcol"
#define SYM_DESC            "symmetry handler for row and column symmetries"
#define SYM_PRIORITY           -100000       /**< priority of try-add function*/
#define SYM_PROPPRIORITY       -100000       /**< priority of propagation method */
#define SYM_PROPTIMING SCIP_PROPTIMING_BEFORELP /**< timing of propagator method */
#define SYM_PROPFREQ                 1       /**< frequence of propagator method */
#define SYM_DELAYPROP            FALSE       /**< Should propagation method be delayed, if other propagators found reductions? */
#define SYM_PRESOLPRIORITY    -1000000       /**< priority of presolving method */
#define SYM_MAXPRESOLROUNDS          1       /**< maximum number of presolving rounds */

/* default value of parameters */
#define SYM_DEFAULT_DETECTDOUBLELEX TRUE     /**< Should we check whether the components can be handled by double lex matrices? */
#define SYM_DEFAULT_HANDLESIGNEDORBITOPES TRUE /**< Shall signed orbitopes be handled? */
#define SYM_DEFAULT_USEDYNAMICPROP TRUE      /**< Shall dynamic propagation be used for orbitopes? */


/*
 * Data structures
 */

/** symmetry component data */
struct SCIP_SymCompData
{
   SCIP_CONS**           conss;              /**< static symmetry handling constraints */
   int                   nconss;             /**< number of static symmetry handling constraints */
   int                   maxnconss;          /**< maximum number of constraints conss can hold */
   SCIP_LEXREDDATA*      lexreddata;         /**< container for lexicographic reduction propagation; */
   SCIP_Bool             lexredactive;       /**< whether lexicographic reduction is active on this component */
};

/** symmetry handler data */
struct SCIP_SymhdlrData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr; /**< shadow tree event handler */
   SCIP_ORBITOPALREDDATA* orbitopalreddata;  /**< container for orbitopal reduction data */
   SCIP_Bool             detectdoublelex;    /**< Should we check whether the components can be handled by double lex matrices? */
   SCIP_Bool             handlesignedorbitopes;/**< Shall signed orbitopes be handled? */
   SCIP_Bool             usedynamicprop;     /**< Shall dynamic propagation be used for orbitopes? */
};

/*
 * Local methods
 */

/** checks whether a (signed) permutation is an involution */
static
SCIP_RETCODE isPermInvolution(
   int*                  perm,               /**< permutation */
   int                   permlen,            /**< number of original (non-negated) variables in a permutation */
   SCIP_Bool*            isinvolution,       /**< pointer to store whether perm is an involution */
   int*                  ntwocycles          /**< pointer to store number of 2-cycles in an involution */
   )
{
   int v;

   assert( perm != NULL );
   assert( permlen > 0 );
   assert( isinvolution != NULL );
   assert( ntwocycles != NULL );

   *ntwocycles = 0;
   *isinvolution = TRUE;
   for (v = 0; v < permlen && *isinvolution; ++v)
   {
      /* do not handle variables twice */
      if ( perm[v] <= v )
         continue;

      /* detect two cycles */
      if ( perm[perm[v]] == v )
         ++(*ntwocycles);
      else
         *isinvolution = FALSE;
   }

   return SCIP_OKAY;
}

/** checks whether selected permutations define orbitopal symmetries */
static
SCIP_RETCODE detectOrbitopalSymmetries(
   SCIP*                 scip,               /**< SCIP pointer */
   int**                 perms,              /**< array of permutations */
   int*                  selectedperms,      /**< indices of permutations in perm that shall be considered */
   int                   nselectedperms,     /**< number of permutations in selectedperms */
   int                   permlen,            /**< number of variables in a permutation */
   int                   nrows,              /**< number of rows of potential orbitopal symmetries */
   SCIP_Bool*            success,            /**< pointer to store if orbitopal symmetries could be found */
   int****               matrices,           /**< pointer to store matrices of orbitopal symmetries */
   int**                 ncols,              /**< pointer to store number of columns of matrices in matrices */
   int*                  nmatrices           /**< pointer to store the number of matrices in matrices */
   )
{  /*lint --e{771}*/
   SCIP_DISJOINTSET* conncomps;
   SCIP_DISJOINTSET* compcolors;
   int* complastperm;
   int* permstoconsider;
   int* colorbegins;
   int* compidx;
   int* colidx;
   int* varidx;
   int* degrees;
   int* perm;
   int nposdegree = 0;
   int npermstoconsider;
   int colorrepresentative1 = -1;
   int colorrepresentative2 = -1;
   int elemtomove;
   int ncurcols;
   int curcomp1;
   int curcomp2;
   int curdeg1;
   int curdeg2;
   int curcolor1;
   int curcolor2;
   int ncolors;
   int cnt;
   int c;
   int p;
   int v;
   int w;

   assert( scip != NULL );
   assert( perms != NULL );
   assert( selectedperms != NULL );
   assert( nselectedperms >= 0 );
   assert( permlen > 0 );
   assert( nrows > 0 || nselectedperms == 0 );
   assert( success != NULL );
   assert( matrices != NULL );
   assert( ncols != NULL );
   assert( nmatrices != NULL );

   /* initialize data */
   *success = TRUE;
   *matrices = NULL;
   *ncols = NULL;
   *nmatrices = 0;

   /* we have found the empty set of orbitopal symmetries */
   if ( nselectedperms == 0 )
      return SCIP_OKAY;

   /* build a graph to encode potential orbitopal symmetries
    *
    * The 2-cycles of a permutation define a set of edges that need to be added simultaneously. We encode
    * this as a disjoint set data structure to only encode the connected components of the graph. To ensure
    * correctness, we keep track of the degrees of each node and extend a component by a permutation only if
    * all nodes to be extended have the same degree. We also make sure that each connected component is a
    * path. Although this might not detect all orbitopal symmetries, it seems to cover most of the cases when
    * computing symmetries using bliss. We also keep track of which variables are affected by the same
    * permutations by coloring the connected components.
    */
   SCIP_CALL( SCIPcreateDisjointset(scip, &conncomps, permlen) );
   SCIP_CALL( SCIPcreateDisjointset(scip, &compcolors, permlen) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &degrees, permlen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &complastperm, permlen) );
   for (v = 0; v < permlen; ++v)
      complastperm[v] = -1;

   /* try to add edges of permutations to graph */
   for (p = 0; p < nselectedperms; ++p)
   {
      perm = perms[selectedperms[p]];
      curdeg1 = -1;
      curdeg2 = -1;

      /* check whether all potential edges can be added */
      for (v = 0; v < permlen; ++v)
      {
         /* treat each cycle exactly once */
         if ( perm[v] <= v )
            continue;
         w = perm[v];

         curcomp1 = SCIPdisjointsetFind(conncomps, v);
         curcomp2 = SCIPdisjointsetFind(conncomps, w);

         /* an edge is not allowed to connect nodes from the same connected component */
         if ( curcomp1 == curcomp2 )
            break;

         /* permutation p is not allowed to add two edges to the same connected component */
         if ( complastperm[curcomp1] == p || complastperm[curcomp2] == p )
            break;

         /* get colors of nodes */
         curcolor1 = SCIPdisjointsetFind(compcolors, v);
         curcolor2 = SCIPdisjointsetFind(compcolors, w);

         /* an edge is not allowed to connect two nodes with the same color */
         if ( curcolor1 == curcolor2 )
            break;

         if ( curdeg1 == -1 )
         {
            assert( curdeg2 == -1 );

            curdeg1 = degrees[v];
            curdeg2 = degrees[w];
            colorrepresentative1 = curcolor1;
            colorrepresentative2 = curcolor2;

            /* stop, we will generate a vertex with degree 3 */
            if ( curdeg1 == 2 || curdeg2 == 2 )
               break;
         }
         else
         {
            /* check whether nodes have compatible degrees */
            if ( ! ((curdeg1 == degrees[v] && curdeg2 == degrees[w])
                  || (curdeg1 == degrees[w] && curdeg2 == degrees[v])) )
               break;
            assert( colorrepresentative1 >= 0 );
            assert( colorrepresentative2 >= 0 || curdeg2 == -1 );

            /* check whether all components have compatible colors */
            if ( curdeg1 > 0 && curcolor1 != colorrepresentative1 && curcolor2 != colorrepresentative1 )
               break;
            if ( curdeg2 > 0 && curcolor1 != colorrepresentative2 && curcolor2 != colorrepresentative2 )
               break;
         }

         /* store that permutation p extends the connected components */
         complastperm[curcomp1] = p;
         complastperm[curcomp2] = p;
      }

      /* terminate if not all edges can be added */
      if ( v < permlen )
      {
         *success = FALSE;
         goto FREEMEMORY;
      }
      assert( curdeg1 >= 0 && curdeg2 >= 0 );

      /* add edges to graph */
      for (v = 0; v < permlen; ++v)
      {
         /* treat each cycle exactly once */
         if ( perm[v] <= v )
            continue;
         w = perm[v];

#ifndef NDEBUG
         curcomp1 = SCIPdisjointsetFind(conncomps, v);
         curcomp2 = SCIPdisjointsetFind(conncomps, w);
         assert( curcomp1 != curcomp2 );
#endif

         /* add edge */
         SCIPdisjointsetUnion(conncomps, v, w, FALSE);
         ++degrees[v];
         ++degrees[w];

         /* possibly update colors */
         curcolor1 = SCIPdisjointsetFind(compcolors, v);
         curcolor2 = SCIPdisjointsetFind(compcolors, w);

         if ( curcolor1 != curcolor2 )
         {
            /* coverity[negative_returns] */
            SCIPdisjointsetUnion(compcolors, colorrepresentative1, v, TRUE);
            SCIPdisjointsetUnion(compcolors, colorrepresentative1, w, TRUE);
         }
      }
   }

   /* find non-trivial components */
   for (v = 0; v < permlen; ++v)
   {
      if ( degrees[v] > 0 )
         ++nposdegree;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &compidx, nposdegree) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colidx, nposdegree) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varidx, nposdegree) );

   for (v = 0, w = 0; v < permlen; ++v)
   {
      if ( degrees[v] > 0 )
      {
         compidx[w] = SCIPdisjointsetFind(conncomps, v);
         colidx[w] = SCIPdisjointsetFind(compcolors, v);
#ifndef NDEBUG
         if ( w > 0 && compidx[w] == compidx[w-1] )
            assert( colidx[w] == colidx[w-1]);
#endif
         varidx[w++] = v;
      }
   }
   assert( w == nposdegree );

   /* sort variable indices: first by colors, then by components */
   SCIP_CALL( SCIPallocBufferArray(scip, &colorbegins, permlen + 1) );

   SCIPsortIntIntInt(colidx, compidx, varidx, nposdegree);
   w = 0;
   ncolors = 0;
   colorbegins[0] = 0;
   for (v = 1; v < nposdegree; ++v)
   {
      if ( colidx[v] != colidx[w] )
      {
         SCIPsortIntInt(&compidx[w], &varidx[w], v - w);
         colorbegins[++ncolors] = v;
         w = v;
      }
   }
   SCIPsortIntInt(&compidx[w], &varidx[w], nposdegree - w);
   colorbegins[++ncolors] = nposdegree;

   /* try to find the correct order of variable indices per color class */
   SCIP_CALL( SCIPallocBufferArray(scip, &permstoconsider, nselectedperms) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, matrices, ncolors) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, ncols, ncolors) );
   *nmatrices = ncolors;

   for (c = 0; c < ncolors; ++c)
   {
      /* find an element in the first connected component with degree 1 */
      for (v = colorbegins[c]; compidx[v] == compidx[colorbegins[c]]; ++v)
      {
         assert( v < nposdegree ); /* there should be a node of degree 1 */

         if ( degrees[varidx[v]] == 1 )
            break;
      }
      assert( compidx[v] == compidx[colorbegins[c]] );
      elemtomove = varidx[v];

      /* find the permutations affecting the variables in the first connected component */
      npermstoconsider = 0;
      for (p = 0; p < nselectedperms; ++p)
      {
         perm = perms[selectedperms[p]];
         for (v = colorbegins[c]; v < nposdegree && compidx[v] == compidx[colorbegins[c]]; ++v)
         {
            if ( perm[varidx[v]] != varidx[v] )
            {
               permstoconsider[npermstoconsider++] = selectedperms[p];
               break;
            }
         }
      }

      /* allocate memory for matrix */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*matrices)[c], nrows) );
      (*ncols)[c] = npermstoconsider + 1;
      for (p = 0; p < nrows; ++p)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*matrices)[c][p], (*ncols)[c]) );
      }

      /* find a permutation that moves the degree-1 element and iteratively extend this to a matrix */
      assert( degrees[elemtomove] == 1 );

      /* find the first and second column */
      for (p = 0; p < npermstoconsider; ++p)
      {
         perm = perms[permstoconsider[p]];
         if ( perm[elemtomove] != elemtomove )
            break;
      }
      assert( p < npermstoconsider );

      /* elements moved by perm that have degree 1 are in the first column */
      for (v = 0, cnt = 0; v < permlen; ++v)
      {
         if ( perm[v] > v ) /*lint !e771*/
         {
            if ( degrees[v] == 1 )
            {
               (*matrices)[c][cnt][0] = v;
               (*matrices)[c][cnt++][1] = perm[v];
            }
            else
            {
               (*matrices)[c][cnt][0] = perm[v];
               (*matrices)[c][cnt++][1] = v;
            }
         }
      }

      /* if the selected permutation do not form orbitopal symmetries */
      if ( cnt < nrows )
      {
         *success = FALSE;
         for (p = nrows - 1; p >= 0; --p)
         {
            SCIPfreeBufferArray(scip, &(*matrices)[c][p]);
         }
         SCIPfreeBlockMemoryArray(scip, &(*matrices)[c], nrows);
         SCIPfreeBlockMemoryArray(scip, matrices, ncolors);
         SCIPfreeBlockMemoryArray(scip, ncols, ncolors);
         *matrices = NULL;
         *ncols = NULL;
         goto FREEMOREMEMORY;
      }
      assert( cnt == nrows );

      /* remove p from the list of permutations to be considered */
      permstoconsider[p] = permstoconsider[--npermstoconsider];

      ncurcols = 1;
      while ( npermstoconsider > 0 )
      {
         elemtomove = (*matrices)[c][0][ncurcols];

         /* find permutation moving the elemtomove */
         for (p = 0; p < npermstoconsider; ++p)
         {
            perm = perms[permstoconsider[p]];
            if ( perm[elemtomove] != elemtomove )
               break;
         }
         assert( p < npermstoconsider );

         /* extend matrix */
         for (v = 0; v < nrows; ++v)
         {
            assert( perm[(*matrices)[c][v][ncurcols]] != (*matrices)[c][v][ncurcols] );
            (*matrices)[c][v][ncurcols + 1] = perm[(*matrices)[c][v][ncurcols]];
         }
         ++ncurcols;
         permstoconsider[p] = permstoconsider[--npermstoconsider];
      }
   }

 FREEMOREMEMORY:
   SCIPfreeBufferArray(scip, &permstoconsider);
   SCIPfreeBufferArray(scip, &colorbegins);
   SCIPfreeBufferArray(scip, &varidx);
   SCIPfreeBufferArray(scip, &colidx);
   SCIPfreeBufferArray(scip, &compidx);

 FREEMEMORY:
   SCIPfreeBufferArray(scip, &complastperm);
   SCIPfreeBufferArray(scip, &degrees);
   SCIPfreeDisjointset(scip, &compcolors);
   SCIPfreeDisjointset(scip, &conncomps);

   return SCIP_OKAY;
}

/** checks whether two families of orbitopal symmetries define a double lex matrix, and in case of success, generates matrix
 *
 *  The columns of matrix1 will serve as the columns of the matrix to be generated, the columns of matrix2 will
 *  serve as rows.
 */
static
SCIP_RETCODE isDoublelLexSym(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nsymvars,           /**< number of variables on which symmetries act */
   int***                matrices1,          /**< first list of matrices associated with orbitopal symmetries */
   int                   nrows1,             /**< number of rows of first family of matrices */
   int*                  ncols1,             /**< for each matrix in the first family, its number of columns */
   int                   nmatrices1,         /**< number of matrices in the first family */
   int***                matrices2,          /**< second list of matrices associated with orbitopal symmetries */
   int                   nrows2,             /**< number of rows of second family of matrices */
   int*                  ncols2,             /**< for each matrix in the second family, its number of columns */
   int                   nmatrices2,         /**< number of matrices in the second family */
   int***                doublelexmatrix,    /**< pointer to store combined matrix */
   int*                  nrows,              /**< pointer to store number of rows in combined matrix */
   int*                  ncols,              /**< pointer to store number of columns in combined matrix */
   int**                 rowsbegin,          /**< pointer to store the begin positions of a new lex subset of rows */
   int**                 colsbegin,          /**< pointer to store the begin positions of a new lex subset of columns */
   SCIP_Bool*            success             /**< pointer to store whether combined matrix could be generated */
   )
{
   int* idxtomatrix1;
   int* idxtomatrix2;
   int* idxtorow1;
   int* idxtorow2;
   int* idxtocol1;
   int* idxtocol2;
   int* sortvals;
   int elem;
   int mat;
   int col;
   int col2;
   int mat2;
   int cnt;
   int c;
   int d;
   int i;
   int j;

   assert( scip != NULL );
   assert( nsymvars >= 0 );
   assert( matrices1 != NULL );
   assert( nrows1 > 0 );
   assert( ncols1 != NULL );
   assert( nmatrices1 > 0 );
   assert( matrices2 != NULL );
   assert( nrows2 > 0 || nmatrices2 == 0 );
   assert( ncols2 != NULL );
   assert( nmatrices2 >= 0 );
   assert( doublelexmatrix != NULL );
   assert( nrows != NULL );
   assert( ncols != NULL );
   assert( rowsbegin != NULL );
   assert( colsbegin != NULL );
   assert( success != NULL );

   /* initialize data */
   *nrows = nrows1;
   *ncols = nrows2;
   *success = TRUE;

   /* check whether expecteded sizes of matrix match */
   for (j = 0, cnt = 0; j < nmatrices1; ++j)
      cnt += ncols1[j];
   if ( cnt != *ncols )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   for (i = 0, cnt = 0; i < nmatrices2; ++i)
      cnt += ncols2[i];
   if ( cnt != *nrows )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* collect information about entries in matrices */
   SCIP_CALL( SCIPallocBufferArray(scip, &idxtomatrix1, nsymvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &idxtomatrix2, nsymvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &idxtorow1, nsymvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &idxtorow2, nsymvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &idxtocol1, nsymvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &idxtocol2, nsymvars) );

   /* use separate loops for efficiency reasons */
   for (i = 0; i < nsymvars; ++i)
      idxtomatrix1[i] = -1;
   for (i = 0; i < nsymvars; ++i)
      idxtomatrix2[i] = -1;
   for (i = 0; i < nsymvars; ++i)
      idxtorow1[i] = -1;
   for (i = 0; i < nsymvars; ++i)
      idxtorow2[i] = -1;
   for (i = 0; i < nsymvars; ++i)
      idxtocol1[i] = -1;
   for (i = 0; i < nsymvars; ++i)
      idxtocol2[i] = -1;

   for (c = 0; c < nmatrices1; ++c)
   {
      for (i = 0; i < nrows1; ++i)
      {
         for (j = 0; j < ncols1[c]; ++j)
         {
            idxtomatrix1[matrices1[c][i][j]] = c;
            idxtorow1[matrices1[c][i][j]] = i;
            idxtocol1[matrices1[c][i][j]] = j;
         }
      }
   }
   for (c = 0; c < nmatrices2; ++c)
   {
      for (i = 0; i < nrows2; ++i)
      {
         for (j = 0; j < ncols2[c]; ++j)
         {
            idxtomatrix2[matrices2[c][i][j]] = c;
            idxtorow2[matrices2[c][i][j]] = i;
            idxtocol2[matrices2[c][i][j]] = j;
         }
      }
   }

   /* check whether the variables of the two orbitopes coincide */
   for (i = 0; i < nsymvars; ++i)
   {
      if ( (idxtomatrix1[i] == -1) != (idxtomatrix2[i] == -1) )
      {
         *success = FALSE;
         goto FREEINITMEMORY;
      }
   }

   /* Find a big matrix such that the columns of this matrix correspond to the columns of matrices in matrices1
    * and the rows of this matrix correspond to the columns of matrices in matrices2. In total, this leads to
    * a matrix of block shape
    *
    *                     A | B | ... | C
    *                     D | E | ... | F
    *                     . | . |     | .
    *                     G | H | ... | I
    *
    * We start by filling the first column of the big matrix by the first column in matrices1. Sort the column
    * according to the matrices in matrices2.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortvals, MAX(*nrows, *ncols)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, doublelexmatrix, *nrows) );
   for (i = 0; i < *nrows; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*doublelexmatrix)[i], *ncols) );
      (*doublelexmatrix)[i][0] = matrices1[0][i][0];
      sortvals[i] = idxtomatrix2[matrices1[0][i][0]];
   }
   SCIPsortIntPtr(sortvals, (void*) (*doublelexmatrix), *nrows);

   /* check that the first column can be covered by rows of matrices2 */
   cnt = 0;
   for (i = 0; i < nmatrices2; ++i)
   {
      int end;

      end = cnt + ncols2[i] - 1;
      for (j = cnt; j < end; ++j)
      {
         assert( idxtomatrix2[(*doublelexmatrix)[j][0]] == idxtomatrix2[(*doublelexmatrix)[j + 1][0]] );
         if( idxtorow2[(*doublelexmatrix)[j][0]] != idxtorow2[(*doublelexmatrix)[j + 1][0]] )
         {
            *success = FALSE;
            goto FREEMEMORY;
         }
      }
   }

   /* fill first row of big matrix */
   mat = idxtomatrix2[(*doublelexmatrix)[0][0]];
   col = idxtocol2[(*doublelexmatrix)[0][0]];
   cnt = 0;
   for (j = 0; j < *ncols; ++j)
   {
      /* skip the entry that is already contained in the first column */
      if ( matrices2[mat][j][col] == (*doublelexmatrix)[0][0] )
         continue;

      sortvals[cnt++] = idxtomatrix1[matrices2[mat][j][col]];
      (*doublelexmatrix)[0][cnt] = matrices2[mat][j][col];
   }
   assert( cnt == nrows2 - 1);
   SCIPsortIntInt(sortvals, &((*doublelexmatrix)[0][1]), cnt);

   /* fill the remaining entries of the big matrix */
   for (i = 1; i < *nrows; ++i)
   {
      for (j = 1; j < *ncols; ++j)
      {
         /* get the matrices and column/row of the entry */
         mat = idxtomatrix1[(*doublelexmatrix)[0][j]];
         mat2 = idxtomatrix2[(*doublelexmatrix)[i][0]];
         col = idxtocol1[(*doublelexmatrix)[0][j]];
         col2 = idxtocol2[(*doublelexmatrix)[i][0]];

         /* find the unique element in the col column of matrix mat and the row column of matrix mat2 */
         /* @todo improve this by first sorting the columns */
         cnt = 0;
         elem = -1;
         for (c = 0; c < *nrows; ++c)
         {
            for (d = 0; d < *ncols; ++d)
            {
               if ( matrices1[mat][c][col] == matrices2[mat2][d][col2] )
               {
                  ++cnt;
                  elem = matrices1[mat][c][col];
                  break;
               }
            }
         }

         /* stop: the columns do not overlap properly */
         if ( cnt != 1 )
         {
            *success = FALSE;
            goto FREEMEMORY;
         }
         (*doublelexmatrix)[i][j] = elem;
      }
   }

   /* store begin positions of row and column blocks */
   if ( *success )
   {
      assert( nmatrices2 > 0 );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, rowsbegin, nmatrices2 + 1) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, colsbegin, nmatrices1 + 1) );
      (*rowsbegin)[0] = 0;
      (*colsbegin)[0] = 0;
      for (j = 0; j < nmatrices2; ++j)
         (*rowsbegin)[j + 1] = (*rowsbegin)[j] + ncols2[j];
      for (j = 0; j < nmatrices1; ++j)
         (*colsbegin)[j + 1] = (*colsbegin)[j] + ncols1[j];
   }

   /* check whether the rows of doublelexmatrix are covered by rows of matrices1 */
   for (i = 0; i < *nrows; ++i)
   {
      for (c = 0; c < nmatrices1; ++c)
      {
         for (j = (*colsbegin)[c]; j < (*colsbegin)[c + 1] - 1; ++j)
         {
            assert( idxtomatrix1[(*doublelexmatrix)[i][j]] == idxtomatrix1[(*doublelexmatrix)[i][j + 1]] );
            if( idxtorow1[(*doublelexmatrix)[i][j]] != idxtorow1[(*doublelexmatrix)[i][j + 1]] )
            {
               *success = FALSE;
               goto FREEMEMORY;
            }
         }
      }
   }

   /* check whether the columns of doublelexmatrix are covered by rows of matrices1 */
   for (i = 0; i < *ncols; ++i)
   {
      for (c = 0; c < nmatrices2; ++c)
      {
         for (j = (*rowsbegin)[c]; j < (*rowsbegin)[c + 1] - 1; ++j)
         {
            assert( idxtomatrix2[(*doublelexmatrix)[j][i]] == idxtomatrix2[(*doublelexmatrix)[j + 1][i]] );
            if( idxtorow2[(*doublelexmatrix)[j][i]] != idxtorow2[(*doublelexmatrix)[j + 1][i]] )
            {
               *success = FALSE;
               goto FREEMEMORY;
            }
         }
      }
   }

 FREEMEMORY:
   SCIPfreeBufferArray(scip, &sortvals);

   if ( !(*success) )
   {
      for (i = *nrows - 1; i >= 0; --i)
      {
         SCIPfreeBlockMemoryArray(scip, &(*doublelexmatrix)[i], *ncols);
      }
      SCIPfreeBlockMemoryArray(scip, doublelexmatrix, *nrows);
      *doublelexmatrix = NULL;
      *rowsbegin = NULL;
      *colsbegin = NULL;
   }

 FREEINITMEMORY:
   SCIPfreeBufferArray(scip, &idxtocol2);
   SCIPfreeBufferArray(scip, &idxtocol1);
   SCIPfreeBufferArray(scip, &idxtorow2);
   SCIPfreeBufferArray(scip, &idxtorow1);
   SCIPfreeBufferArray(scip, &idxtomatrix2);
   SCIPfreeBufferArray(scip, &idxtomatrix1);

   return SCIP_OKAY;
}

/** detects whether permutations define single or double lex matrices
 *
 *  A single lex matrix is a matrix whose columns can be partitioned into blocks such that the
 *  columns within each block can be permuted arbitrarily. A double lex matrix is a single lex
 *  matrix such that also blocks of rows have the aforementioned property.
 */
static
SCIP_RETCODE SCIPdetectSingleOrDoubleLexMatrices(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_Bool             detectsinglelex,    /**< whether single lex matrices shall be detected */
   int**                 perms,              /**< array of permutations */
   int                   nperms,             /**< number of permutations in perms */
   int                   permlen,            /**< number of variables in a permutation */
   SCIP_Bool*            success,            /**< pointer to store whether structure could be detected */
   SCIP_Bool*            isorbitope,         /**< pointer to store whether detected matrix is orbitopal */
   int***                lexmatrix,          /**< pointer to store single or double lex matrix */
   int*                  nrows,              /**< pointer to store number of rows of lexmatrix */
   int*                  ncols,              /**< pointer to store number of columns of lexmatrix */
   int**                 lexrowsbegin,       /**< pointer to store array indicating begin of new row-lexmatrix */
   int**                 lexcolsbegin,       /**< pointer to store array indicating begin of new col-lexmatrix */
   int*                  nrowmatrices,       /**< pointer to store number of single lex row matrices in rows */
   int*                  ncolmatrices        /**< pointer to store number of single lex column matrices in rows */
   )
{
   int*** matricestype1 = NULL;
   int*** matricestype2 = NULL;
   int* ncolstype1 = NULL;
   int* ncolstype2 = NULL;
   int nmatricestype1 = 0;
   int nmatricestype2 = 0;
   int* permstype1;
   int* permstype2;
   int npermstype1 = 0;
   int npermstype2 = 0;
   int ncycs1 = -1;
   int ncycs2 = -1;
   int tmpncycs;
   int p;
   int i;
   SCIP_Bool isinvolution;

   assert( scip != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( permlen > 0 );
   assert( success != NULL );
   assert( lexmatrix != NULL );
   assert( nrows != NULL );
   assert( ncols != NULL );
   assert( lexrowsbegin != NULL );
   assert( lexcolsbegin != NULL );
   assert( nrowmatrices != NULL );
   assert( ncolmatrices != NULL );

   *success = TRUE;
   *isorbitope = FALSE;
   *nrowmatrices = 0;
   *ncolmatrices = 0;

   /* arrays to store the different types of involutions */
   SCIP_CALL( SCIPallocBufferArray(scip, &permstype1, nperms)  );
   SCIP_CALL( SCIPallocBufferArray(scip, &permstype2, nperms)  );

   /* check whether we can expect lexicographically sorted rows and columns */
   for (p = 0; p < nperms; ++p)
   {
      SCIP_CALL( isPermInvolution(perms[p], permlen, &isinvolution, &tmpncycs) );

      /* terminate if not all permutations are involutions */
      if ( ! isinvolution )
      {
         *success = FALSE;
         goto FREEMEMORY;
      }

      /* store number of cycles or terminate if too many different types of involutions */
      if ( ncycs1 == -1 || ncycs1 == tmpncycs )
      {
         ncycs1 = tmpncycs;
         permstype1[npermstype1++] = p;
      }
      else if ( ncycs2 == -1 || ncycs2 == tmpncycs )
      {
         ncycs2 = tmpncycs;
         permstype2[npermstype2++] = p;
      }
      else
      {
         *success = FALSE;
         goto FREEMEMORY;
      }
   }

   /* for each type, check whether permutations define (disjoint) orbitopal symmetries */
   SCIP_CALL( detectOrbitopalSymmetries(scip, perms, permstype1, npermstype1, permlen, ncycs1, success,
         &matricestype1, &ncolstype1, &nmatricestype1) );
   if ( ! *success )
      goto FREEMEMORY;

   SCIP_CALL( detectOrbitopalSymmetries(scip, perms, permstype2, npermstype2, permlen, ncycs2, success,
         &matricestype2, &ncolstype2, &nmatricestype2) );
   if ( ! *success )
      goto FREEMEMORY;

   /* check whether a double lex matrix is defined */
   *success = FALSE;
   if ( !detectsinglelex && ncycs2 != -1 )
   {
      assert( ncycs1 > 0 );

      SCIP_CALL( isDoublelLexSym(scip, permlen, matricestype1, ncycs1, ncolstype1, nmatricestype1,
            matricestype2, ncycs2, ncolstype2, nmatricestype2,
            lexmatrix, nrows, ncols, lexrowsbegin, lexcolsbegin, success) );

      if ( *success )
      {
         *nrowmatrices = nmatricestype2;
         *ncolmatrices = nmatricestype1;
      }
   }

   /* if no double lex matrix is detected, possibly return orbitope */
   if ( !(*success) && ncycs2 == -1 && nmatricestype1 == 1 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, lexmatrix, ncycs1) );
      for (i = 0; i < ncycs1; ++i)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*lexmatrix)[i], ncolstype1[0]) );
         for (p = 0; p < ncolstype1[0]; ++p)
            (*lexmatrix)[i][p] = matricestype1[0][i][p];
      }
      *nrows = ncycs1;
      *ncols = ncolstype1[0];
      *success = TRUE;
      *isorbitope = TRUE;
   }

 FREEMEMORY:
   for (p = nmatricestype2 - 1; p >= 0; --p)
   {
      for (i = ncycs2 - 1; i >= 0; --i)
      {
         SCIPfreeBlockMemoryArray(scip, &matricestype2[p][i], ncolstype2[p]);
      }
      SCIPfreeBlockMemoryArray(scip, &matricestype2[p], ncycs2);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &matricestype2, nmatricestype2);
   SCIPfreeBlockMemoryArrayNull(scip, &ncolstype2, nmatricestype2);
   for (p = nmatricestype1 - 1; p >= 0; --p)
   {
      for (i = ncycs1 - 1; i >= 0; --i)
      {
         SCIPfreeBlockMemoryArray(scip, &matricestype1[p][i], ncolstype1[p]);
      }
      SCIPfreeBlockMemoryArray(scip, &matricestype1[p], ncycs1);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &matricestype1, nmatricestype1);
   SCIPfreeBlockMemoryArrayNull(scip, &ncolstype1, nmatricestype1);

   SCIPfreeBufferArray(scip, &permstype2);
   SCIPfreeBufferArray(scip, &permstype1);

   return SCIP_OKAY;
}

/** orbitopal reduction */
static
SCIP_RETCODE addOrbitopesDynamic(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   SCIP_VAR**            permvars,           /**< variables the symmetries act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< domain center of variables */
   int                   componentid,        /**< ID of component for which orbitope is added */
   char*                 partialname,        /**< partial name for orbitope constraint */
   int**                 varidxmatrix,       /**< matrix containing variable indices in orbitope matrix */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data needed for orbitopal reduction */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_Bool*            success             /**< pointer to store whether orbitope could be added successfully */
   )
{
   char name[SCIP_MAXSTRLEN];
   int i;
   int j;
   SCIP_Bool ispporbitope;
   SCIP_VAR*** varmatrix;
   SCIP_Bool* pprows;
   int npprows;
   SCIP_ORBITOPETYPE type;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(orbitopalreddata != NULL);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* add linear constraints x_1 >= x_2 >= ... >= x_ncols for single-row orbitopes */
   if( nrows == 1 )
   {
      /* restrict to the packing and partitioning rows */
      SCIP_CONS* cons;
      SCIP_VAR* consvars[2];
      SCIP_Real conscoefs[2] = { -1.0, 1.0 };

      /* for all adjacent column pairs, add linear constraint */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + ncols - 1) );
      for( i = 0; i < ncols - 1; ++i )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_1row_comp_%d_col%d", partialname, componentid, i);

         consvars[0] = permvars[varidxmatrix[0][i]];
         consvars[1] = permvars[varidxmatrix[0][i + 1]];

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, conscoefs, -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         (*orbitopeconss)[(*norbitopeconss)++] = cons;
      }

      *success = TRUE;

      return SCIP_OKAY;
   }

   /* for only 2 columns, the the component can be completely handled by lexicographic reduction */
   if( ncols == 2 )
   {
      int* orbisackperm;

      /* If the component is an orbitope with 2 columns, then there is 1 generator of order 2. */
      orbisackperm = perms[0];

      SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, lexreddata,
            permvars, npermvars, orbisackperm, symtype, permvardomaincenter, TRUE, success) );

      if( *success )
      {
         *lexredactive = TRUE;
         return SCIP_OKAY;
      }
   }

   /* create orbitope variable matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &varmatrix, nrows) );
   for( i = 0; i < nrows; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &varmatrix[i], ncols) );
      for( j = 0; j < ncols; ++j )
         varmatrix[i][j] = permvars[varidxmatrix[i][j]];
   }

   pprows = NULL;
   SCIP_CALL( SCIPisPackingPartitioningOrbitope(scip, varmatrix, nrows, ncols, &pprows, &npprows, &type) );

   /* does it have at least 3 packing-partitioning rows? */
   ispporbitope = npprows >= 3;  /* (use same magic number as cons_orbitope.c) */

   if( ispporbitope ) /* @todo if it's a pporbitope, we do it statically right now. */
   {
      /* restrict to the packing and partitioning rows */
      SCIP_CONS* cons;
      SCIP_VAR*** ppvarsarrayonlypprows;
      int r;

      assert(pprows != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsarrayonlypprows, npprows) );

      r = 0;
      for( i = 0; i < nrows; ++i )
      {
         if( pprows[i] )
         {
            assert(r < npprows);
            ppvarsarrayonlypprows[r++] = varmatrix[i];
         }
      }
      assert(r == npprows);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_pp", partialname);
      SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, ppvarsarrayonlypprows, SCIP_ORBITOPETYPE_PACKING,
            npprows, ncols, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      /* check whether we need to resize */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + 1) );
      (*orbitopeconss)[(*norbitopeconss)++] = cons;
      *success = TRUE;

      SCIPfreeBufferArray(scip, &ppvarsarrayonlypprows);
   }
   else
   {
      /* use orbitopal reduction for component */
      SCIP_COLUMNORDERING columnordering;
      SCIP_VAR** orbitopevarmatrix;
      int nelem;
      int pos = 0;

      /* variable array */
      nelem = nrows * ncols;
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nelem) );
      for( i = 0; i < nrows; ++i )
      {
         for( j = 0; j < ncols; ++j )
            orbitopevarmatrix[pos++] = varmatrix[i][j];
      }

      /* get column ordering */
      columnordering = SCIPorbitopalReductionGetDefaultColumnOrdering(orbitopalreddata);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_full", partialname);
      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
            SCIP_ROWORDERING_BRANCHING, columnordering,
            orbitopevarmatrix, nrows, ncols, success) );
      *success = TRUE;

      SCIPfreeBufferArray(scip, &orbitopevarmatrix);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &pprows, nrows);

   for( i = nrows - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &varmatrix[i]);
   }
   SCIPfreeBufferArray(scip, &varmatrix);

   return SCIP_OKAY;
}

/** returns whether every variable in each row of an orbitope matrix has same center */
static
SCIP_Bool isEquallyCenteredOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_Real*            vardomaincenter,    /**< array containing domain centers of each variable */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   startrow,           /**< row of varidxmatrix in which orbitope begins */
   int                   endrow,             /**< row of varidxmatrix after which orbitope ends */
   int                   startcol,           /**< column of varidxmatrix in which orbitope begins */
   int                   endcol,             /**< column of varidxmatrix after which orbitope ends */
   SCIP_Bool             equalrowcenters     /**< whether rows are centered equally (otherwise, columns) */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(varidxmatrix != NULL);
   assert(0 <= startrow && startrow < endrow);
   assert(0 <= startcol && startcol < endcol);

   /* if we are dealing with unsigned permutations, then all rows have the same center */
   if( vardomaincenter == NULL )
      return TRUE;

   if( equalrowcenters )
   {
      for( i = startrow; i < endrow; ++i )
      {
         for( j = startcol; j < endcol - 1; ++j )
         {
            /* different centers detected */
            if( !SCIPisEQ(scip, vardomaincenter[varidxmatrix[i][j]],
                  vardomaincenter[varidxmatrix[i][j + 1]]) )
               return FALSE;
         }
      }
   }
   else
   {
      for( j = startcol; j < endcol; ++j )
      {
         for( i = startrow; i < endrow - 1; ++i )
         {
            /* different centers detected */
            if( !SCIPisEQ(scip, vardomaincenter[varidxmatrix[i][j]],
                  vardomaincenter[varidxmatrix[i + 1][j]]) )
               return FALSE;
         }
      }
   }

   return TRUE;
}

/** checks whether a proper signed permutation flips a (partial) orbitope column */
static
SCIP_RETCODE hasOrbitopeColumnFlip(
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   startrow,           /**< row of varidxmatrix in which orbitope starts */
   int                   endrow,             /**< row of varidxmatrix after which orbitope ends */
   int                   startcol,           /**< column of varidxmatrix in which orbitope starts */
   int                   endcol,             /**< column of varidxmatrix after which orbitope ends */
   int*                  signedperm,         /**< signed permutation to be checked */
   int                   npermvars,          /**< number of variables symmetries act on */
   SCIP_Bool             transposed,         /**< whether the orbitope is transposed in varidxmatrix */
   int*                  flipablerows,       /**< allocated array to store rows admitting a flip */
   int*                  nflipablerows       /**< pointer to store number of flipable rows */
   )
{
   int flipcolumn = -1;
   int i;
   int j;

   assert(varidxmatrix != NULL);
   assert(0 <= startrow && startrow < endrow);
   assert(0 <= startcol && startcol < endcol);
   assert(signedperm != NULL);
   assert(npermvars >= (endrow - startrow) * (endcol - startcol));
   assert(flipablerows != NULL);
   assert(nflipablerows != NULL);

   *nflipablerows = 0;

   /* iterate over matrix and check for flip */
   for( i = startrow; i < endrow; ++i )
   {
      for( j = startcol; j < endcol; ++j )
      {
         if( signedperm[varidxmatrix[i][j]] == npermvars + varidxmatrix[i][j] )
         {
            if( flipcolumn == -1 )
            {
               flipcolumn = transposed ? i : j;
               flipablerows[(*nflipablerows)++] = transposed ? j : i;
            }
            else if( (transposed && flipcolumn != i) || (!transposed && flipcolumn != j) )
            {
               /* variables from two columns are flipped */
               *nflipablerows = 0;
               return SCIP_OKAY;
            }
            else
               flipablerows[(*nflipablerows)++] = transposed ? j : i;
         }
         else if( signedperm[varidxmatrix[i][j]] != varidxmatrix[i][j] )
         {
            /* variable is not flipped */
            *nflipablerows = 0;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** handles orbitope action by static or dynamic symmetry handling methods */
static
SCIP_RETCODE handleOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< domain centers of permvars */
   int                   componentid,        /**< ID of component to which orbitope is added */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   char*                 partialname,        /**< partial name to be extended by constraints */
   SCIP_Bool             issigned,           /**< whether the first row of the orbitope can be sign-flipped */
   SCIP_Bool             handlestatically,   /**< whether the orbitope shall be handled statically */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool             usedynamicprop,     /**< Shall orbitopes be handled dynamically? */
   SCIP_Bool*            success,            /**< pointer to store whether orbitope could be added successfully */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{
   assert(scip != NULL);
   assert(perms != NULL);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);
   assert(nchgbds != NULL);

   *success = FALSE;

   if( handlestatically )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_Real consvals[2] = {-1.0, 1.0};
      SCIP_VAR* consvars[2];
      SCIP_VAR** orbitopevarmatrix;
      SCIP_CONS* cons;
      int nconss;
      int nelem;
      int pos;
      int i;
      int j;

      /* handle orbitope */
      nelem = nrows * ncols;
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nelem) );
      for( i = 0, pos = 0; i < nrows; ++i )
      {
         for( j = 0; j < ncols; ++j )
            orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
      }

      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata, SCIP_ROWORDERING_NONE,
            SCIP_COLUMNORDERING_NONE, orbitopevarmatrix, nrows, ncols, success) );

      /* compute number of constraints to handle signed part of the orbitope */
      nconss = ncols - 1;
      if( issigned )
         nconss += ncols;

      /* create linear constraints */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + nconss) );

      /* handle symmetries by enforcing sorted variables in first row */
      for( j = 0; j < ncols - 1; ++j )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_%d", partialname, j);
         consvars[0] = orbitopevarmatrix[j];
         consvars[1] = orbitopevarmatrix[j + 1];

         /* enforce constraints to be in LP since this seems to have a positive impact for orbitopes with cont. vars */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals, -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         (*orbitopeconss)[(*norbitopeconss)++] = cons;
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      if( issigned )
      {
         SCIP_VAR* var;
         SCIP_Real bound;

         assert(permvardomaincenter != NULL);

         /* the first row is contained in the upper half of the variable domain */
         for( j = 0; j < ncols; ++j )
         {
            var = orbitopevarmatrix[j];
            bound = permvardomaincenter[varidxmatrix[0][j]];

            if( SCIPisLT(scip, SCIPvarGetLbLocal(var), bound) )
            {
               /* improve lower bound */
               SCIP_CALL( SCIPchgVarLb(scip, var, bound) );
               ++(*nchgbds);
            }
         }
      }
      SCIPfreeBufferArray(scip, &orbitopevarmatrix);

      *success = TRUE;
   }

   /* if symmetries have not been handled yet */
   if( !(*success) )
   {
      /* dynamic propagation */
      if( usedynamicprop )
      {
         SCIP_CALL( addOrbitopesDynamic(scip, symtype, perms, permvars, npermvars, permvardomaincenter, componentid,
               partialname, varidxmatrix, nrows, ncols, orbitopeconss, norbitopeconss, maxnorbitopeconss,
               orbitopalreddata, lexreddata, lexredactive, success) );
      }
      /* static variant only for binary variables */
      else
      {
         SCIP_VAR*** orbitopematrix;
         SCIP_CONS* cons;
         int i;
         int j;
         int nbinrows = 0;

         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, nrows) );
         for( i = 0; i < nrows; ++i )
         {
            /* skip rows without binary variables */
            if( !SCIPvarIsBinary(permvars[varidxmatrix[i][0]]) )
               continue;

            SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[nbinrows], ncols) );
            for( j = 0; j < ncols; ++j )
            {
               assert(SCIPvarIsBinary(permvars[varidxmatrix[i][j]]));
               orbitopematrix[nbinrows][j] = permvars[varidxmatrix[i][j]];
            }
            ++nbinrows;
         }

         if( nbinrows > 0 )
         {
            SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, partialname, orbitopematrix, SCIP_ORBITOPETYPE_FULL,
                  nbinrows, ncols, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                  FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

            /* do not release constraint here - will be done later */
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss, *norbitopeconss + 1) );
            (*orbitopeconss)[(*norbitopeconss)++] = cons;

            *success = TRUE;

            for( i = nbinrows - 1; i >= 0; --i )
            {
               SCIPfreeBufferArray(scip, &orbitopematrix[i]);
            }
         }
         SCIPfreeBufferArray(scip, &orbitopematrix);
      }
   }

   return SCIP_OKAY;
}

/** handles double lex orbitope action by static symmetry handling methods */
static
SCIP_RETCODE handleDoubleLexOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   SCIP_Real*            permvardomaincenter,/**< domain center of variables */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   char*                 partialname,        /**< partial name to be extended by constraints */
   int                   nsignedrows,        /**< the first number of rows that can be sign-flipped */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool*            success,            /**< pointer to store whether orbitope could be added successfully */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real consvals[2] = {-1.0, 1.0};
   SCIP_VAR* consvars[2];
   SCIP_VAR** orbitopevarmatrix;
   SCIP_CONS* cons;
   int nsignedconss;
   int nsortconss;
   int nactiverows;
   int nactrowsprev;
   int nelem;
   int pos;
   int i;
   int j;

   assert(scip != NULL);
   assert(permvars != NULL);
   assert(permvardomaincenter != NULL || nsignedconss == 0);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(0 <= nsignedrows && nsignedrows <= nrows);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);
   assert(nchgbds != NULL);

   /* prepare data for orbitope matrices */
   nelem = nrows * ncols;
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nelem) );

   /* compute number of constraints to handle signed part of the orbitope */
   nsignedconss = 0;
   nsortconss = ncols - 1;
   if( nsignedrows > 0 )
   {
      nactiverows = nsignedrows;
      nactrowsprev = nrows;

      for( j = 0; j < ncols; ++j )
      {
         nsortconss += MAX(nactrowsprev - nactiverows - 1, 0);
         nsignedconss += nactiverows;
         nactrowsprev = nactiverows;

         /* ceil(nactiverows / 2) */
         nactiverows = (int) ((nactiverows + 1) / 2);
      }
      assert(nactiverows >= 1);

      nsortconss += nactiverows - 1;
   }
   else
      nsortconss += nrows - 1;

   /* create linear constraints */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss,
         *norbitopeconss + nsortconss + nsignedconss) );

   /* sort variables in first row */
   for( j = 0; j < ncols - 1; ++j )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_%d", partialname, j);
      consvars[0] = permvars[varidxmatrix[0][j]];
      consvars[1] = permvars[varidxmatrix[0][j + 1]];

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals, -SCIPinfinity(scip), 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      (*orbitopeconss)[(*norbitopeconss)++] = cons;
      SCIP_CALL( SCIPaddCons(scip, cons) );
   }

   /* handle symmetries by enforcing sorted columns via orbitopal reduction */
   for( i = 0, pos = 0; i < nrows; ++i )
   {
      for( j = 0; j < ncols; ++j )
         orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
   }

   SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
         SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
         orbitopevarmatrix, nrows, ncols, success) );

   if( nsignedconss > 0 )
   {
      int k;

      nactiverows = nsignedrows;
      nactrowsprev = nrows;
      for( j = 0; j < ncols; ++j )
      {
         /* ceil(nactiverows / 2) */
         nactiverows = (int) ((nactiverows + 1) / 2);

         /* the second half of active rows can be sorted by linear inequalities */
         for( i = nactiverows; i < nactrowsprev - 1; ++i )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_col_%d_row_%d", partialname, j, i);

            consvars[0] = permvars[varidxmatrix[i][0]];
            consvars[1] = permvars[varidxmatrix[i + 1][0]];

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals,
                  -SCIPinfinity(scip), 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            (*orbitopeconss)[(*norbitopeconss)++] = cons;
            SCIP_CALL( SCIPaddCons(scip, cons) );
         }

         /* we can also sort the second half of rows by orbitopal reduction */
         if( nactrowsprev - nactiverows > 1 )
         {
            for( k = 0, pos = 0; k < ncols; ++k )
            {
               for( i = nactiverows; i < nactrowsprev; ++i )
                  orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][k]];
            }

            SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
                  SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
                  orbitopevarmatrix, ncols, nactrowsprev - nactiverows, success) );
         }
         nactrowsprev = nactiverows;

         /* the first half of the active rows are in the upper part of the variable domain */
         for( i = 0; i < nactiverows; ++i )
         {
            SCIP_Real bound;

            consvars[0] = permvars[varidxmatrix[i][j]];
            bound = permvardomaincenter[varidxmatrix[i][j]];

            /* improve lower bound either by changing the bound or a linear constraint */
            if( SCIPisLT(scip, SCIPvarGetLbLocal(consvars[0]), bound) )
            {
               SCIP_CALL( SCIPchgVarLb(scip, consvars[0], bound) );
               ++(*nchgbds);
            }
         }
      }

      /* within the remaining active rows, the rows can be sorted */
      if( nactiverows > 1 )
      {
         for( i = 0; i < nactiverows - 1; ++i )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sortfirstactive_%d", partialname, i);

            consvars[0] = permvars[varidxmatrix[i][0]];
            consvars[1] = permvars[varidxmatrix[i + 1][0]];

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals,
                  -SCIPinfinity(scip), 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            (*orbitopeconss)[(*norbitopeconss)++] = cons;
            SCIP_CALL( SCIPaddCons(scip, cons) );
         }

         /* also apply orbitopal reduction */
         for( j = 0, pos = 0; j < ncols; ++j )
         {
            for( i = 0; i < nactiverows; ++i )
               orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
         }

         SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
               SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
               orbitopevarmatrix, ncols, nactiverows, success) );
      }
      assert(*norbitopeconss <= *maxnorbitopeconss );
   }
   else
   {
      /* sort first column */
      for( i = 0; i < nrows - 1; ++i )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_sort_col_0_row_%d", partialname, i);

         consvars[0] = permvars[varidxmatrix[i][0]];
         consvars[1] = permvars[varidxmatrix[i + 1][0]];

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, consvals,
               -SCIPinfinity(scip), 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         (*orbitopeconss)[(*norbitopeconss)++] = cons;
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      /* apply orbitopal fixing to row permutations */
      for( j = 0, pos = 0; j < ncols; ++j )
      {
         for( i = 0; i < nrows; ++i )
            orbitopevarmatrix[pos++] = permvars[varidxmatrix[i][j]];
      }
      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, orbitopalreddata,
            SCIP_ROWORDERING_NONE, SCIP_COLUMNORDERING_NONE,
            orbitopevarmatrix, ncols, nrows, success) );
   }

   SCIPfreeBufferArray(scip, &orbitopevarmatrix);

   *success = TRUE;

   return SCIP_OKAY;
}

/** handles double lex matrix */
static
SCIP_RETCODE handleDoublelLexMatrix(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< domain center of variables */
   int                   id,                 /**< ID of component that is handled */
   int**                 varidxmatrix,       /**< matrix containing variable indices of double lex matrix */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   int*                  rowsbegin,          /**< array indicating where a new row block begins */
   int*                  colsbegin,          /**< array indicating where a new column block begins */
   int                   nrowblocks,         /**< number of row blocks */
   int                   ncolblocks,         /**< number of column blocks */
   int**                 signedperms,        /**< array of proper signed permutations */
   int                   nsignedperms,       /**< number of proper signed permutations */
   SCIP_Bool             handlesignedorbitopes,/**< Shall signed orbitopes be handled? */
   SCIP_Bool             usedynamicprop,     /**< Shall dynamic orbitope propagation be used */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool*            success,            /**< pointer to store whether orbitope could be added successfully */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{
   char partialname[SCIP_MAXSTRLEN];
   int** orbitopematrix;
   int* flipableidx;
   int nflipableidx;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_Bool tmpsuccess;
   int maxdim;
   int i;
   int p;
   int j;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(permvars != NULL);
   assert(permvardomaincenter != NULL || symtype != SYM_SYMTYPE_SIGNPERM);
   assert(varidxmatrix != NULL);
   assert(nrows > 0);
   assert(ncols > 0);
   assert(rowsbegin != NULL);
   assert(colsbegin != NULL);
   assert(nrowblocks > 0);
   assert(ncolblocks > 0);
   assert(signedperms != NULL);
   assert(nsignedperms >= 0);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);
   assert(nchgbds != NULL);

   *success = FALSE;

   maxdim = MAX(nrows, ncols);
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, maxdim) );
   for( i = 0; i < maxdim; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[i], maxdim) );
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &flipableidx, maxdim) );

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

   /* detect how double-lex matrices can be handled
    *
    * We use the following strategy:
    *
    * (1) In case of 1 column- and row-block, we check whether one of them can be handled
    *     by signed orbitopes. If this is the case, reorder the rows or columns based
    *     on the signed permutations.
    *
    * (2) If there are multiple column- and row-blocks, just handle classical permutation symmetries.
    */

   if( nrowblocks == 1 && ncolblocks == 1 )
   {
      /* Case 1 */
      SCIP_Bool hascolflip = FALSE;
      SCIP_Bool hasrowflip = FALSE;
      SCIP_Bool canusecolorbitope = FALSE;
      SCIP_Bool canuseroworbitope = FALSE;
      int q;

      /* check whether orbitopes can be used to handle column and row swaps (requires equally centered rows/columns) */
      if( isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix, 0, nrows, 0, ncols, TRUE) )
         canusecolorbitope = TRUE;
      if( isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix, 0, nrows, 0, ncols, FALSE) )
         canuseroworbitope = TRUE;

      nflipableidx = 0;
      if( handlesignedorbitopes && canusecolorbitope )
      {
         /* check whether the signed permutations flip entries within a single column of the orbitope matrix */
         for( q = 0; q < nsignedperms && nflipableidx == 0; ++q )
         {
            SCIP_CALL( hasOrbitopeColumnFlip(varidxmatrix, 0, nrows, 0, ncols,
                  signedperms[q], npermvars, FALSE, flipableidx, &nflipableidx) );
         }

         /* a signed orbitope can be used to handle the column symmetries */
         if( nflipableidx > 0 )
            hascolflip = TRUE;
      }

      if( handlesignedorbitopes && !hascolflip && canuseroworbitope )
      {
         assert(nflipableidx == 0);

         /* check whether the signed permutations flip entries within a single column of transposed orbitope matrix */
         for( q = 0; q < nsignedperms && nflipableidx == 0; ++q )
         {
            SCIP_CALL( hasOrbitopeColumnFlip(varidxmatrix, 0, nrows, 0, ncols,
                  signedperms[q], npermvars, TRUE, flipableidx, &nflipableidx) );
         }

         /* a signed orbitope can be used to handle the row symmetries */
         if( nflipableidx > 0 )
            hasrowflip = TRUE;
      }

      /* handle row and column symmetries by potentially reordered orbitopes */
      if( hascolflip )
      {
         int isigned = 0;
         int iunsigned;

         /* handle column symmetries by reordering the rows of the orbitope */
         iunsigned = nflipableidx;
         for( i = 0; i < nrows; ++i )
         {
            if( isigned < nflipableidx && flipableidx[isigned] == i )
            {
               for( j = 0; j < ncols; ++j )
                  orbitopematrix[isigned][j] = varidxmatrix[i][j];
               ++isigned;
            }
            else
            {
               for( j = 0; j < ncols; ++j )
                  orbitopematrix[iunsigned][j] = varidxmatrix[i][j];
               ++iunsigned;
            }
         }
         assert(isigned == nflipableidx);
         assert(iunsigned == nrows);

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_col_0", id);

         SCIP_CALL( handleDoubleLexOrbitope(scip, permvars, permvardomaincenter, orbitopematrix, nrows, ncols,
               partialname, nflipableidx, orbitopalreddata, orbitopeconss, norbitopeconss, maxnorbitopeconss,
               &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }
      else if( hasrowflip )
      {
         int jsigned = 0;
         int junsigned;

         /* handle row symmetries by reordering the columns of the orbitope */
         junsigned = nflipableidx;
         for( j = 0; j < ncols; ++j )
         {
            if( jsigned < nflipableidx && flipableidx[jsigned] == j )
            {
               for( i = 0; i < nrows; ++i )
                  orbitopematrix[jsigned][i] = varidxmatrix[i][j];
               ++jsigned;
            }
            else
            {
               for( i = 0; i < nrows; ++i )
                  orbitopematrix[junsigned][i] = varidxmatrix[i][j];
               ++junsigned;
            }
         }
         assert(jsigned == nflipableidx);
         assert(junsigned == ncols);

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_row_0", id);

         SCIP_CALL( handleDoubleLexOrbitope(scip, permvars, permvardomaincenter, orbitopematrix, ncols, nrows,
               partialname, nflipableidx, orbitopalreddata, orbitopeconss, norbitopeconss, maxnorbitopeconss,
               &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }
   }

   /* if no symmetries have been handled yet, handle column and row symmetries without signed permutations */
   if( !(*success) )
   {
      /* ensure that we can store orbitope constraints in probdata */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitopeconss, maxnorbitopeconss,
            *norbitopeconss + nrowblocks + ncolblocks) );

      /* handle column symmetries via original column and row ordering */
      for( p = 0; p < ncolblocks; ++p )
      {
         int jj;
         j = 0;

         /* we can only handle the orbitope if all variables in a row have the same domain center */
         if( !isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix, 0, nrows,
               colsbegin[p], colsbegin[p + 1], TRUE) )
            continue;

         /* create the orbitope matrix */
         for( i = 0; i < nrows; ++i )
         {
            for( j = 0, jj = colsbegin[p]; jj < colsbegin[p + 1]; ++j, ++jj )
               orbitopematrix[i][j] = varidxmatrix[i][jj];
         }

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_col_%d", id, p);

         SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, id, orbitopematrix,
               nrows, j, partialname, FALSE, TRUE, lexreddata, lexredactive, orbitopalreddata,
               orbitopeconss, norbitopeconss, maxnorbitopeconss, usedynamicprop, &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }

      /* handle row symmetries via original column and row ordering */
      for( p = 0; p < nrowblocks; ++p )
      {
         int ii;

         /* we can only handle the orbitope if all variables in a row have the same domain center */
         if( !isEquallyCenteredOrbitope(scip, permvardomaincenter, varidxmatrix,
               rowsbegin[p], rowsbegin[p + 1], 0, ncols, FALSE) )
            continue;

         /* create the orbitope matrix */
         for( i = 0, ii = rowsbegin[p]; ii < rowsbegin[p + 1]; ++i, ++ii )
         {
            for( j = 0; j < ncols; ++j )
               orbitopematrix[j][i] = varidxmatrix[ii][j];
         }

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d_doublelex_row_%d", id, p);

         SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, id, orbitopematrix,
               ncols, i, partialname, FALSE, TRUE, lexreddata, lexredactive, orbitopalreddata, orbitopeconss,
               norbitopeconss, maxnorbitopeconss, usedynamicprop, &tmpsuccess, nchgbds) );
         *success = *success || tmpsuccess;
      }
   }

   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &flipableidx);

   for( i = maxdim - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &orbitopematrix[i]);
   }
   SCIPfreeBufferArray(scip, &orbitopematrix);

   return SCIP_OKAY;
}

/** tries to handle symmetries of single lex matrices (orbitopes) or double lex matrices */
static
SCIP_RETCODE tryHandleSingleOrDoubleLexMatrices(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables the permutations act on */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< array of centers of variable domains */
   SCIP_Bool             detectsinglelex,    /**< whether single lex matrices shall be detected */
   SCIP_Bool             handlesignedorbitopes,/**< Shall signed orbitopes be handled? */
   SCIP_Bool             usedynamicprop,     /**< Shall dynamic propagation be used for orbitopes? */
   int                   cidx,               /**< index of component */
   int*                  nchgbds,            /**< pointer to store number of bound changes (or NULL) */
   SCIP_LEXREDDATA*      lexreddata,         /**< data needed for lexicographic reduction */
   SCIP_Bool*            lexredactive,       /**< pointer to store whether lexicographic reduction is active */
   SCIP_ORBITOPALREDDATA* orbitopalreddata,  /**< data for orbitopal reduction */
   SCIP_CONS***          orbitopeconss,      /**< pointer to store array of conss derived from orbitopes */
   int*                  norbitopeconss,     /**< pointer to store number of conss in orbitopeconss */
   int*                  maxnorbitopeconss,  /**< pointer to store maximum number of conss orbitopeconss can hold */
   SCIP_Bool*            success             /**< pointer to store whether symmetries are handled */
   )
{
   int** lexmatrix = NULL;
   int* lexrowsbegin = NULL;
   int* lexcolsbegin = NULL;
   int nrows;
   int ncols;
   int nrowmatrices;
   int ncolmatrices;
   int** selectedperms = NULL;
   int i;
   int nonpermidx = -1;
   SCIP_Real percentageunsigned;
   SCIP_Bool isorbitope;
   SCIP_Bool locsuccess = FALSE;
   SCIP_Bool* isproperperm;
   int nselectedperms = 0;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(nchgbds != NULL);
   assert(lexreddata != NULL);
   assert(lexredactive != NULL);
   assert(orbitopalreddata != NULL);
   assert(orbitopeconss != NULL);
   assert(norbitopeconss != NULL);
   assert(maxnorbitopeconss != NULL);
   assert(success != NULL);

   *nchgbds = 0;
   *success = FALSE;

   /* get proper permutations of component and possibly store index of one non-permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &isproperperm, nperms) );
   if( symtype == SYM_SYMTYPE_PERM )
   {
      selectedperms = perms;
      nselectedperms = nperms;
      for( i = 0; i < nperms; ++i )
         isproperperm[i] = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &selectedperms, nperms) );
      for( i = 0; i < nperms; ++i )
      {
         isproperperm[i] = isProperPerm(symtype, perms[i], npermvars);
         if( isproperperm[i] )
            selectedperms[nselectedperms++] = perms[i];
         else
            nonpermidx = i;
      }
   }
   percentageunsigned = (SCIP_Real) nselectedperms / (SCIP_Real) nperms;

   if( nselectedperms == 0 )
   {
      SCIPfreeBufferArray(scip, &isproperperm);
      if( symtype != SYM_SYMTYPE_PERM )
      {
         SCIPfreeBufferArrayNull(scip, &selectedperms);
      }
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPdetectSingleOrDoubleLexMatrices(scip, detectsinglelex, selectedperms, nselectedperms, npermvars,
         &locsuccess, &isorbitope, &lexmatrix, &nrows, &ncols, &lexrowsbegin, &lexcolsbegin,
         &nrowmatrices, &ncolmatrices) );

   /* possibly handle double lex matrix or orbitope */
   if( locsuccess )
   {
      assert(lexmatrix != NULL);
      assert(nrows > 0);
      assert(ncols > 0);

      /* possibly store non-permutation symmetries (in many cases, there is only one) */
      if( nselectedperms != nperms )
      {
         assert(nselectedperms < nperms);

         if( nselectedperms == nperms - 1 )
         {
            selectedperms[0] = selectedperms[nonpermidx];
            nselectedperms = 1;
         }
         else
         {
            nselectedperms = 0;
            for( i = 0; i < nperms; ++i )
            {
               if( !isproperperm[i] )
                  selectedperms[nselectedperms++] = perms[i];
            }
         }
      }

      if( isorbitope )
      {
         int** orbitopematrix;
         char partialname[SCIP_MAXSTRLEN];

         locsuccess = FALSE;

         /* signed permutations can only handle the orbitope if all variables per row have the same domain center */
         if( symtype != (int) SYM_SYMTYPE_PERM )
         {
            if( !isEquallyCenteredOrbitope(scip, permvardomaincenter, lexmatrix, 0, nrows, 0, ncols, TRUE) )
               goto FREEMEMORY;
         }

         (void) SCIPsnprintf(partialname, SCIP_MAXSTRLEN, "orbitope_component_%d", cidx);

         if( handlesignedorbitopes )
         {
            int* flipablerows;
            int nflipablerows = 0;
            int p;

            SCIP_CALL( SCIPallocBufferArray(scip, &flipablerows, nrows) );

            /* check whether the signed permutations flip entries within a single column of the orbitope matrix
             *
             * It is sufficient to find one such signed permutations, because only one row will incorporate
             * information about the sign change.
             */
            for( p = 0; p < nselectedperms && nflipablerows == 0; ++p )
            {
               SCIP_CALL( hasOrbitopeColumnFlip(lexmatrix, 0, nrows, 0, ncols, selectedperms[p], npermvars, FALSE,
                     flipablerows, &nflipablerows) );
            }

            /* possibly flip rows to be able to handle signed orbitopes */
            if( nflipablerows > 0 )
            {
               int j;
               int isigned = 0;
               int iunsigned;

               iunsigned = nflipablerows;
               SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, nrows) );
               for( i = 0; i < nrows; ++i )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[i], ncols) );

                  if( isigned < nflipablerows && flipablerows[isigned] == i )
                  {
                     for( j = 0; j < ncols; ++j )
                        orbitopematrix[isigned][j] = lexmatrix[i][j];
                     ++isigned;
                  }
                  else
                  {
                     for( j = 0; j < ncols; ++j )
                        orbitopematrix[iunsigned][j] = lexmatrix[i][j];
                     ++iunsigned;
                  }
               }
               assert(isigned == nflipablerows);
               assert(iunsigned == nrows);

               SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, cidx,
                     orbitopematrix, nrows, ncols, partialname, TRUE, TRUE, lexreddata, lexredactive, orbitopalreddata,
                     orbitopeconss, norbitopeconss, maxnorbitopeconss, usedynamicprop, &locsuccess, nchgbds) );

               for( i = nrows - 1; i >= 0; --i )
               {
                  SCIPfreeBufferArray(scip, &orbitopematrix[i]);
               }
               SCIPfreeBufferArray(scip, &orbitopematrix);
            }
            SCIPfreeBufferArray(scip, &flipablerows);
         }

         /* if we have not handled the orbitope yet, handle it as unsigned orbitope and the orbitope is large */
         if( (!locsuccess) && percentageunsigned > 0.8 )
         {
            SCIP_CALL( handleOrbitope(scip, symtype, perms, permvars, npermvars, permvardomaincenter, cidx, lexmatrix,
                  nrows, ncols, partialname, FALSE, FALSE, lexreddata, lexredactive, orbitopalreddata, orbitopeconss,
                  norbitopeconss, maxnorbitopeconss, usedynamicprop, &locsuccess, nchgbds) );
         }
      }
      else
      {
         SCIP_CALL( handleDoublelLexMatrix(scip, symtype, perms, permvars, npermvars, permvardomaincenter, cidx,
               lexmatrix, nrows, ncols, lexrowsbegin, lexcolsbegin, nrowmatrices, ncolmatrices, selectedperms, nselectedperms,
               handlesignedorbitopes, usedynamicprop, lexreddata, lexredactive, orbitopalreddata, orbitopeconss,
               norbitopeconss, maxnorbitopeconss, &locsuccess, nchgbds) );
      }

   FREEMEMORY:
      /* free memory not needed anymore */
      for( i = nrows - 1; i >= 0; --i )
      {
         SCIPfreeBlockMemoryArray(scip, &lexmatrix[i], ncols);
      }
      SCIPfreeBlockMemoryArray(scip, &lexmatrix, nrows);

      if( ncolmatrices > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &lexcolsbegin, ncolmatrices + 1);
      }

      if( nrowmatrices > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &lexrowsbegin, nrowmatrices + 1);
      }
   }
   SCIPfreeBufferArray(scip, &isproperperm);
   if( symtype != SYM_SYMTYPE_PERM )
   {
      SCIPfreeBufferArrayNull(scip, &selectedperms);
   }

   if( locsuccess )
      *success = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of symmetry handler
 */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryaddRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(success != NULL);
   assert(symhdlr != NULL);
   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms >= 0);
   assert(permvars != NULL || npermvars == 0);
   assert(naddedconss != NULL);
   assert(nchgbds != NULL);

   *success = FALSE;
   *nchgbds = 0;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);
   assert(symhdlrdata->shadowtreeeventhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );
   (*symcompdata)->conss = NULL;
   (*symcompdata)->nconss = 0;
   (*symcompdata)->maxnconss = 0;
   (*symcompdata)->lexreddata = NULL;
   (*symcompdata)->lexredactive = FALSE;

   SCIP_CALL( SCIPincludeLexicographicReduction(scip, &(*symcompdata)->lexreddata, symhdlrdata->shadowtreeeventhdlr) );
   assert((*symcompdata)->lexreddata != NULL);

   SCIP_CALL( tryHandleSingleOrDoubleLexMatrices(scip, symtype, perms, nperms, permvars, npermvars, permvardomcenter,
         !symhdlrdata->detectdoublelex, symhdlrdata->handlesignedorbitopes, symhdlrdata->usedynamicprop,
         id, nchgbds, (*symcompdata)->lexreddata, &(*symcompdata)->lexredactive, symhdlrdata->orbitopalreddata,
         &(*symcompdata)->conss, &(*symcompdata)->nconss, &(*symcompdata)->maxnconss, success) );

   if( !(*success) )
   {
      assert(*nchgbds == 0);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, (*symcompdata)->lexreddata) );
      SCIP_CALL( SCIPlexicographicReductionFree(scip, &(*symcompdata)->lexreddata) );

      SCIPfreeBlockMemory(scip, symcompdata);
      *symcompdata = NULL;

      return SCIP_OKAY;
   }

   *naddedconss = (*symcompdata)->nconss;

   return SCIP_OKAY;
}

/** solving process deinitialization method of symmetry handler (called before branch and bound process data is freed) */
static
SCIP_DECL_SYMHDLREXITSOL(symhdlrExitsolRowCol)
{ /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symdata;
   int s;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata != NULL);
      assert(symdata->lexreddata != NULL || !symdata->lexredactive);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
   }
   assert(symhdlrdata->orbitopalreddata != NULL);

   SCIP_CALL( SCIPorbitopalReductionReset(scip, symhdlrdata->orbitopalreddata) );

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symhdlrExitRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(symcomps != NULL || nsymcomps == 0);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata != NULL);
      assert(symdata->lexreddata != NULL || !symdata->lexredactive);

      SCIP_CALL( SCIPlexicographicReductionReset(scip, symdata->lexreddata) );
      SCIP_CALL( SCIPlexicographicReductionFree(scip, &symdata->lexreddata) );

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &symdata->conss[c]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &symdata->conss, symdata->maxnconss);

      SCIPfreeBlockMemory(scip, &symdata);
   }
   assert(symhdlrdata->orbitopalreddata != NULL);

   SCIP_CALL( SCIPorbitopalReductionReset(scip, symhdlrdata->orbitopalreddata) );

   return SCIP_OKAY;
}

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
static
SCIP_DECL_SYMHDLRFREE(symhdlrFreeRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(scip != NULL);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);
   assert(symhdlrdata->orbitopalreddata != NULL);

   SCIP_CALL( SCIPorbitopalReductionFree(scip, &symhdlrdata->orbitopalreddata) );

   SCIPfreeBlockMemory(scip, &symhdlrdata);

   return SCIP_OKAY;
}

/** domain propagation method of symmetry handler */
static
SCIP_DECL_SYMHDLRPROP(symhdlrPropRowCol)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMCOMPDATA* symcompdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool didrun = FALSE;
   SCIP_Bool didrunlocal;
   int nredlocal;
   int nreds = 0;
   int s;

   assert(symhdlr != NULL);
   assert(symcomps != NULL || nsymcomps == 0);
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   /* run orbital reduction */
   SCIP_CALL( SCIPorbitopalReductionPropagate(scip, symhdlrdata->orbitopalreddata,
         &infeasible, &nreds, &didrun) );

   if( infeasible )
      return SCIP_OKAY;

   /* run lexicographic reduction */
   for( s = 0; s < nsymcomps; ++s )
   {
      assert(symcomps[s] != NULL);

      symcompdata = SCIPsymcompGetData(symcomps[s]);
      assert(symcompdata != NULL);

      if( !symcompdata->lexredactive )
         continue;

     SCIP_CALL( SCIPlexicographicReductionPropagate(scip, symcompdata->lexreddata,
            &infeasible, &nredlocal, &didrunlocal) );
      nreds += nredlocal;
      didrun = didrun || didrunlocal;
      if( infeasible )
         return SCIP_OKAY;
   }

   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nreds > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( didrun )
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRESOL(symhdlrPresolRowCol)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(result != NULL);
   assert(symcomps != NULL || nsymcomps == 0);

   *result = nsymcomps > 0 ? SCIP_DIDNOTFIND : SCIP_DIDNOTRUN;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPpresolCons(scip, symdata->conss[c], nrounds, presoltiming, nnewfixedvars, nnewaggrvars,
               nnewchgvartypes, nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs,
               nnewchgsides, nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes, ndelconss, naddconss,
               nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff or unboundedness has been detected */
         if( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n",
               SCIPconsGetName(symdata->conss[c]));
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** symmetry component display method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRINT(symhdlrPrintRowCol)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** include symmetry handler for row and column symmetries */
SCIP_RETCODE SCIPincludeSymhdlrRowCol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;
   SCIP_SYMHDLR* symhdlr;
   SCIP_EVENTHDLR* eventhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &symhdlrdata) );

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, &symhdlr, SYM_NAME, SYM_DESC, SYM_PRIORITY,
         symhdlrTryaddRowCol, symhdlrdata) );

   /* set non-fundamental callback methods */
   SCIP_CALL( SCIPsetSymhdlrFree(scip, symhdlr, symhdlrFreeRowCol) );
   SCIP_CALL( SCIPsetSymhdlrExit(scip, symhdlr, symhdlrExitRowCol) );
   SCIP_CALL( SCIPsetSymhdlrExitsol(scip, symhdlr, symhdlrExitsolRowCol) );
   SCIP_CALL( SCIPsetSymhdlrProp(scip, symhdlr, symhdlrPropRowCol, SYM_PROPFREQ, SYM_DELAYPROP,
         SYM_PROPPRIORITY, SYM_PROPTIMING) );
   SCIP_CALL( SCIPsetSymhdlrPresol(scip, symhdlr, symhdlrPresolRowCol, SYM_MAXPRESOLROUNDS,
         SYM_PRESOLPRIORITY, SCIP_PRESOLTIMING_FAST) );

   /* include shadow tree event handler if it is not included yet */
   eventhdlr = SCIPfindEventhdlr(scip, "event_shadowtree");
   if( eventhdlr == NULL )
   {
      SCIP_CALL( SCIPincludeEventHdlrShadowTree(scip, &symhdlrdata->shadowtreeeventhdlr) );
      assert(symhdlrdata->shadowtreeeventhdlr != NULL);
   }
   else
      symhdlrdata->shadowtreeeventhdlr = eventhdlr;

   SCIP_CALL( SCIPincludeOrbitopalReduction(scip, &symhdlrdata->orbitopalreddata) );
   assert(symhdlrdata->orbitopalreddata != NULL);

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/detectdoublelex",
         "Should we check whether the components can be handled by double lex matrices?",
         &symhdlrdata->detectdoublelex, TRUE, SYM_DEFAULT_DETECTDOUBLELEX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/handlesignedorbitopes",
         "Shall signed orbitopes be handled?",
         &symhdlrdata->handlesignedorbitopes, TRUE, SYM_DEFAULT_HANDLESIGNEDORBITOPES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/usedynamicprop",
         "Shall dynamic orbitope propagation be used?",
         &symhdlrdata->usedynamicprop, TRUE, SYM_DEFAULT_USEDYNAMICPROP, NULL, NULL) );

   return SCIP_OKAY;
}
