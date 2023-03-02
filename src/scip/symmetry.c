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

/**@file   symmetry.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries
 * @author Christopher Hojny
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/symmetry.h"
#include "scip/scip.h"
#include "scip/cons_setppc.h"
#include <scip/misc.h>
#include <symmetry/struct_symmetry.h>
#include <symmetry/type_symmetry.h>


/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 */
SCIP_RETCODE SCIPcomputeOrbitsSym(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< variables considered in a permutation */
   int                   npermvars,          /**< length of a permutation array */
   int**                 perms,              /**< matrix containing in each row a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits             /**< pointer to number of orbits currently stored in orbits */
   )
{
   SCIP_Shortbool* varadded;
   int orbitidx = 0;
   int i;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( npermvars > 0 );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits != NULL );

   /* init data structures*/
   SCIP_CALL( SCIPallocBufferArray(scip, &varadded, npermvars) );

   /* initially, every variable is contained in no orbit */
   for (i = 0; i < npermvars; ++i)
      varadded[i] = FALSE;

   /* find variable orbits */
   *norbits = 0;
   for (i = 0; i < npermvars; ++i)
   {
      int beginorbitidx;
      int j;

      /* skip variable already contained in an orbit of a previous variable */
      if ( varadded[i] )
         continue;

      /* store first variable */
      beginorbitidx = orbitidx;
      orbits[orbitidx++] = i;
      varadded[i] = TRUE;

      /* iterate over variables in curorbit and compute their images */
      j = beginorbitidx;
      while ( j < orbitidx )
      {
         int curelem;
         int image;
         int p;

         curelem = orbits[j];

         for (p = 0; p < nperms; ++p)
         {
            image = perms[p][curelem];

            /* found new element of the orbit of i */
            if ( ! varadded[image] )
            {
               orbits[orbitidx++] = image;
               assert( orbitidx <= npermvars );
               varadded[image] = TRUE;
            }
         }
         ++j;
      }

      /* if the orbit is trivial, reset storage, otherwise store orbit */
      if ( orbitidx <= beginorbitidx + 1 )
         orbitidx = beginorbitidx;
      else
         orbitbegins[(*norbits)++] = beginorbitidx;
   }

   /* store end in "last" orbitbegins entry */
   assert( *norbits < npermvars );
   orbitbegins[*norbits] = orbitidx;

#ifdef SCIP_OUTPUT
   printf("Orbits (total number: %d):\n", *norbits);
   for (i = 0; i < *norbits; ++i)
   {
      int j;

      printf("%d: ", i);
      for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         printf("%d ", orbits[j]);
      printf("\n");
   }
#endif

   /* free memory */
   SCIPfreeBufferArray(scip, &varadded);

   return SCIP_OKAY;
}


/** compute non-trivial orbits of symmetry group using filtered generators
 *
 *  The non-trivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 *
 *  Only permutations that are not inactive (as marked by @p inactiveperms) are used. Thus, one can use this array to
 *  filter out permutations.
 */
SCIP_RETCODE SCIPcomputeOrbitsFilterSym(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< length of a permutation array */
   int**                 permstrans,         /**< transposed matrix containing in each column a
                                              *   permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   SCIP_Shortbool*       inactiveperms,      /**< array to store whether permutations are inactive */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  components,         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins,    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned*             componentblocked,   /**< array to store which symmetry methods have been used on a component
                                              *   using the same bitset information as for misc/usesymmetry */
   int                   ncomponents,        /**< number of components of symmetry group */
   int                   nmovedpermvars      /**< number of variables moved by any permutation in a symmetry component
                                              *   that is handled by orbital fixing */
   )
{
   SCIP_Shortbool* varadded;
   int nvaradded = 0;
   int orbitidx = 0;
   int i;

   assert( scip != NULL );
   assert( permstrans != NULL );
   assert( nperms > 0 );
   assert( npermvars > 0 );
   assert( inactiveperms != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits != NULL );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( vartocomponent != NULL );
   assert( ncomponents > 0 );
   assert( nmovedpermvars > 0 );

   /* init data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &varadded, npermvars) );

   /* initially, every variable is contained in no orbit */
   for (i = 0; i < npermvars; ++i)
      varadded[i] = FALSE;

   /* find variable orbits */
   *norbits = 0;
   for (i = 0; i < npermvars; ++i)
   {
      int beginorbitidx;
      int componentidx;
      int j;

      /* skip unaffected variables and blocked components */
      componentidx = vartocomponent[i];
      if ( componentidx < 0 || componentblocked[componentidx] )
         continue;

      /* skip variable already contained in an orbit of a previous variable */
      if ( varadded[i] )
         continue;

      /* store first variable */
      beginorbitidx = orbitidx;
      orbits[orbitidx++] = i;
      varadded[i] = TRUE;
      ++nvaradded;

      /* iterate over variables in curorbit and compute their images */
      j = beginorbitidx;
      while ( j < orbitidx )
      {
         int* pt;
         int curelem;
         int image;
         int p;

         curelem = orbits[j];

         pt = permstrans[curelem];
         for (p = componentbegins[componentidx]; p < componentbegins[componentidx + 1]; ++p)
         {
            int perm;

            perm = components[p];

            if ( ! inactiveperms[perm] )
            {
               image = pt[perm];
               assert( vartocomponent[image] == componentidx );

               /* found new element of the orbit of i */
               if ( ! varadded[image] )
               {
                  orbits[orbitidx++] = image;
                  assert( orbitidx <= npermvars );
                  varadded[image] = TRUE;
                  ++nvaradded;
               }
            }
         }
         ++j;
      }

      /* if the orbit is trivial, reset storage, otherwise store orbit */
      if ( orbitidx <= beginorbitidx + 1 )
         orbitidx = beginorbitidx;
      else
         orbitbegins[(*norbits)++] = beginorbitidx;

      /* stop if all variables are covered */
      if ( nvaradded >= nmovedpermvars )
         break;
   }

   /* store end in "last" orbitbegins entry */
   assert( *norbits < npermvars );
   orbitbegins[*norbits] = orbitidx;

#ifdef SCIP_OUTPUT
   printf("Orbits (total number: %d):\n", *norbits);
   for (i = 0; i < *norbits; ++i)
   {
      int j;

      printf("%d: ", i);
      for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         printf("%d ", orbits[j]);
      printf("\n");
   }
#endif

   /* free memory */
   SCIPfreeBufferArray(scip, &varadded);

   return SCIP_OKAY;
}

/** Compute orbit of a given variable and store it in @p orbit. The first entry of the orbit will
 *  be the given variable index and the rest is filled with the remaining variables excluding
 *  the ones specified in @p ignoredvars.
 *
 *  @pre orbit is an initialized array of size propdata->npermvars
 *  @pre at least one of @p perms and @p permstrans should not be NULL
 */
SCIP_RETCODE SCIPcomputeOrbitVar(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< number of variables in permvars */
   int**                 perms,              /**< the generators of the permutation group (or NULL) */
   int**                 permstrans,         /**< the transposed matrix of generators (or NULL) */
   int*                  components,         /**< the components of the permutation group */
   int*                  componentbegins,    /**< array containing the starting index of each component */
   SCIP_Shortbool*       ignoredvars,        /**< array indicating which variables should be ignored */
   SCIP_Shortbool*       varfound,           /**< bitmap to mark which variables have been added (or NULL) */
   int                   varidx,             /**< index of variable for which the orbit is requested */
   int                   component,          /**< component that var is in */
   int*                  orbit,              /**< array in which the orbit should be stored */
   int*                  orbitsize           /**< buffer to store the size of the orbit */
   )
{  /*lint --e{571}*/
   SCIP_Shortbool* varadded;
   int* varstotest;
   int nvarstotest;
   int j;
   int p;

   assert( scip != NULL );
   assert( perms != NULL || permstrans != NULL );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( ignoredvars != NULL );
   assert( orbit != NULL );
   assert( orbitsize != NULL );
   assert( 0 <= varidx && varidx < npermvars );
   assert( component >= 0 );
   assert( npermvars > 0 );

   /* init data structures*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &varadded, npermvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &varstotest, npermvars) );

   /* compute and store orbit if it is non-trivial */
   orbit[0] = varidx;
   varstotest[0] = varidx;
   *orbitsize = 1;
   nvarstotest = 1;
   varadded[varidx] = TRUE;

   if ( varfound != NULL )
      varfound[varidx] = TRUE;

   /* iterate over variables in orbit and compute their images */
   j = 0;
   while ( j < nvarstotest )
   {
      int currvar;

      currvar = varstotest[j++];

      for (p = componentbegins[component]; p < componentbegins[component+1]; ++p)
      {
         int image;
         int comp;

         comp = components[p];

         if ( perms != NULL )
            image = perms[comp][currvar]; /*lint !e613*/
         else
            image = permstrans[currvar][comp];

         /* found new element of the orbit of varidx */
         if ( ! varadded[image] )
         {
            varstotest[nvarstotest++] = image;
            varadded[image] = TRUE;

            if ( ! ignoredvars[image] )
            {
               orbit[(*orbitsize)++] = image;

               if ( varfound != NULL )
                  varfound[image] = TRUE;
            }
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &varstotest);
   SCIPfreeBufferArray(scip, &varadded);

   return SCIP_OKAY;
}

/** compute non-trivial orbits of symmetry group
 *
 *  The non-tivial orbits of the group action are stored in the array orbits of length npermvars. This array contains
 *  the indices of variables from the permvars array such that variables that are contained in the same orbit appear
 *  consecutively in the orbits array. The variables of the i-th orbit have indices
 *  orbits[orbitbegins[i]], ... , orbits[orbitbegins[i + 1] - 1].
 *  Note that the description of the orbits ends at orbitbegins[norbits] - 1.
 *
 *  This function is adapted from computeGroupOrbitsFilter().
 */
SCIP_RETCODE SCIPcomputeOrbitsComponentsSym(
   SCIP*                 scip,               /**< SCIP instance */
   int                   npermvars,          /**< length of a permutation array */
   int**                 permstrans,         /**< transposed matrix containing in each column a permutation of the symmetry group */
   int                   nperms,             /**< number of permutations encoded in perms */
   int*                  components,         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins,    /**< array containing in i-th position the first position of component i in components array */
   int*                  vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   int                   ncomponents,        /**< number of components of symmetry group */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int*                  norbits,            /**< pointer to number of orbits currently stored in orbits */
   int*                  varorbitmap         /**< array for storing the orbits for each variable */
   )
{
   SCIP_Shortbool* varadded;
   int orbitidx = 0;
   int i;

   assert( scip != NULL );
   assert( permstrans != NULL );
   assert( nperms > 0 );
   assert( npermvars > 0 );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( vartocomponent != NULL );
   assert( ncomponents > 0 );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits != NULL );
   assert( varorbitmap != NULL );

   /* init data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &varadded, npermvars) );

   /* initially, every variable is contained in no orbit */
   for (i = 0; i < npermvars; ++i)
   {
      varadded[i] = FALSE;
      varorbitmap[i] = -1;
   }

   /* find variable orbits */
   *norbits = 0;
   for (i = 0; i < npermvars; ++i)
   {
      int beginorbitidx;
      int componentidx;
      int j;

      /* skip unaffected variables - note that we also include blocked components */
      componentidx = vartocomponent[i];
      if ( componentidx < 0 )
         continue;

      /* skip variable already contained in an orbit of a previous variable */
      if ( varadded[i] )
         continue;

      /* store first variable */
      beginorbitidx = orbitidx;
      orbits[orbitidx++] = i;
      varadded[i] = TRUE;
      varorbitmap[i] = *norbits;

      /* iterate over variables in curorbit and compute their images */
      j = beginorbitidx;
      while ( j < orbitidx )
      {
         int* pt;
         int curelem;
         int image;
         int p;

         curelem = orbits[j];

         pt = permstrans[curelem];
         for (p = componentbegins[componentidx]; p < componentbegins[componentidx + 1]; ++p)
         {
            int perm;

            perm = components[p];
            image = pt[perm];
            assert( vartocomponent[image] == componentidx );

            /* found new element of the orbit of i */
            if ( ! varadded[image] )
            {
               orbits[orbitidx++] = image;
               assert( orbitidx <= npermvars );
               varadded[image] = TRUE;
               varorbitmap[image] = *norbits;
            }
         }
         ++j;
      }

      /* if the orbit is trivial, reset storage, otherwise store orbit */
      if ( orbitidx <= beginorbitidx + 1 )
      {
         orbitidx = beginorbitidx;
         varorbitmap[i] = -1;
      }
      else
         orbitbegins[(*norbits)++] = beginorbitidx;
   }

   /* store end in "last" orbitbegins entry */
   assert( *norbits < npermvars );
   orbitbegins[*norbits] = orbitidx;

   /* free memory */
   SCIPfreeBufferArray(scip, &varadded);

   return SCIP_OKAY;
}


/** Checks whether a permutation is a composition of 2-cycles and in this case determines the number of overall
 *  2-cycles and binary 2-cycles. It is a composition of 2-cycles iff @p ntwocyclesperm > 0 upon termination.
 */
SCIP_RETCODE SCIPisInvolutionPerm(
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< array of variables perm is acting on */
   int                   nvars,              /**< number of variables */
   int*                  ntwocyclesperm,     /**< pointer to store number of 2-cycles or 0 if perm is not an involution */
   int*                  nbincyclesperm,     /**< pointer to store number of binary cycles */
   SCIP_Bool             earlytermination    /**< whether we terminate early if not all affected variables are binary */
   )
{
   int ntwocycles = 0;
   int i;

   assert( perm != NULL );
   assert( vars != NULL );
   assert( ntwocyclesperm != NULL );
   assert( nbincyclesperm != NULL );

   *ntwocyclesperm = 0;
   *nbincyclesperm = 0;
   for (i = 0; i < nvars; ++i)
   {
      assert( 0 <= perm[i] && perm[i] < nvars );

      /* skip fixed points and avoid treating the same 2-cycle twice */
      if ( perm[i] <= i )
         continue;

      if ( perm[perm[i]] == i )
      {
         if ( SCIPvarIsBinary(vars[i]) && SCIPvarIsBinary(vars[perm[i]]) )
            ++(*nbincyclesperm);
         else if ( earlytermination )
            return SCIP_OKAY;

         ++ntwocycles;
      }
      else
      {
         /* we do not have only 2-cycles */
         return SCIP_OKAY;
      }
   }

   /* at this point the permutation is a composition of 2-cycles */
   *ntwocyclesperm = ntwocycles;

   return SCIP_OKAY;
}


/** determine number of variables affected by symmetry group */
SCIP_RETCODE SCIPdetermineNVarsAffectedSym(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutations */
   int                   nperms,             /**< number of permutations in perms */
   SCIP_VAR**            permvars,           /**< variables corresponding to permutations */
   int                   npermvars,          /**< number of permvars in perms */
   int*                  nvarsaffected       /**< pointer to store number of all affected variables */
   )
{
   SCIP_Shortbool* affected;
   int i;
   int p;

   assert( scip != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( nvarsaffected != NULL );

   *nvarsaffected = 0;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &affected, npermvars) );

   /* iterate over permutations and check which variables are affected by some symmetry */
   for (p = 0; p < nperms; ++p)
   {
      for (i = 0; i < npermvars; ++i)
      {
         if ( affected[i] )
            continue;

         if ( perms[p][i] != i )
         {
            affected[i] = TRUE;
            ++(*nvarsaffected);
         }
      }
   }
   SCIPfreeBufferArray(scip, &affected);

   return SCIP_OKAY;
}


/** Given a matrix with nrows and \#perms + 1 columns whose first nfilledcols columns contain entries of variables, this routine
 *  checks whether the 2-cycles of perm intersect each row of column coltoextend in exactly one position. In this case,
 *  we add one column to the suborbitope of the first nfilledcols columns.
 *
 *  @pre Every non-trivial cycle of perm is a 2-cycle.
 *  @pre perm has nrows many 2-cycles
 */
SCIP_RETCODE SCIPextendSubOrbitope(
   int**                 suborbitope,        /**< matrix containing suborbitope entries */
   int                   nrows,              /**< number of rows of suborbitope */
   int                   nfilledcols,        /**< number of columns of suborbitope which are filled with entries */
   int                   coltoextend,        /**< index of column that should be extended by perm */
   int*                  perm,               /**< permutation */
   SCIP_Bool             leftextension,      /**< whether we extend the suborbitope to the left */
   int**                 nusedelems,         /**< pointer to array storing how often an element was used in the orbitope */
   SCIP_VAR**            permvars,           /**< permutation vars array */
   SCIP_Shortbool*       rowisbinary,        /**< array encoding whether variables in an orbitope row are binary (or NULL) */
   SCIP_Bool*            success,            /**< pointer to store whether extension was successful */
   SCIP_Bool*            infeasible          /**< pointer to store if the number of intersecting cycles is too small */
   )
{
   int nintersections = 0;
   int row;
   int idx1;
   int idx2;

   assert( suborbitope != NULL );
   assert( nrows > 0 );
   assert( nfilledcols > 0 );
   assert( coltoextend >= 0 );
   assert( perm != NULL );
   assert( nusedelems != NULL );
   assert( permvars != NULL );
   assert( success != NULL );
   assert( infeasible != NULL );

   *success = FALSE;
   *infeasible = FALSE;

   /* if we try to extend the first orbitope generator by another one,
    * we can change the order of entries in each row of suborbitope */
   if ( nfilledcols == 2 )
   {
      /* check whether each cycle of perm intersects with a row of suborbitope */
      for (row = 0; row < nrows; ++row)
      {
         idx1 = suborbitope[row][0];
         idx2 = suborbitope[row][1];

         /* if idx1 or idx2 is affected by perm, we can extend the row of the orbitope */
         if ( idx1 != perm[idx1] )
         {
            /* change order of idx1 and idx2 for right extensions */
            if ( ! leftextension )
            {
               suborbitope[row][0] = idx2;
               suborbitope[row][1] = idx1;
            }
            assert( rowisbinary == NULL || rowisbinary[row] == SCIPvarIsBinary(permvars[perm[idx1]]) );

            suborbitope[row][2] = perm[idx1];
            ++nintersections;

            ++(*nusedelems)[idx1];
            ++(*nusedelems)[perm[idx1]];

            /* if an element appears too often in the orbitope matrix */
            if ( (*nusedelems)[idx1] + (*nusedelems)[perm[idx1]] > 3 )
            {
               *infeasible = TRUE;
               break;
            }
         }
         else if ( idx2 != perm[idx2] )
         {
            /* change order of idx1 and idx2 for left extensions */
            if ( leftextension )
            {
               suborbitope[row][0] = idx2;
               suborbitope[row][1] = idx1;
            }
            assert( rowisbinary == NULL || rowisbinary[row] == SCIPvarIsBinary(permvars[perm[idx1]]) );

            suborbitope[row][2] = perm[idx2];
            ++nintersections;

            ++(*nusedelems)[idx2];
            ++(*nusedelems)[perm[idx2]];

            /* if an element appears too often in the orbitope matrix */
            if ( (*nusedelems)[idx2] + (*nusedelems)[perm[idx2]] > 3 )
            {
               *infeasible = TRUE;
               break;
            }
         }
      }
   }
   else
   {
      /* check whether each cycle of perm intersects with a row of suborbitope */
      for (row = 0; row < nrows; ++row)
      {
         idx1 = suborbitope[row][coltoextend];

         /* if idx1 is affected by perm, we can extend the row of the orbitope */
         if ( idx1 != perm[idx1] )
         {
            assert( rowisbinary == NULL || rowisbinary[row] == SCIPvarIsBinary(permvars[perm[idx1]]) );

            suborbitope[row][nfilledcols] = perm[idx1];
            ++nintersections;

            ++(*nusedelems)[idx1];
            ++(*nusedelems)[perm[idx1]];

            /* if an element appears to often in the orbitope matrix */
            if ( (*nusedelems)[idx1] + (*nusedelems)[perm[idx1]] > 3 )
            {
               *infeasible = TRUE;
               break;
            }
         }
      }
   }

   /* if there are too few intersection, this is not an orbitope */
   if ( nintersections > 0 && nintersections < nrows )
      *infeasible = TRUE;
   else if ( nintersections == nrows )
      *success = TRUE;

   return SCIP_OKAY;
}


/** compute components of symmetry group */
SCIP_RETCODE SCIPcomputeComponentsSym(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 perms,              /**< permutation generators as
                                              *   (either nperms x npermvars or npermvars x nperms) matrix */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
   int                   npermvars,          /**< number of variables for permutations */
   SCIP_Bool             transposed,         /**< transposed permutation generators as (npermvars x nperms) matrix */
   int**                 components,         /**< array containing the indices of permutations sorted by components */
   int**                 componentbegins,    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int**                 vartocomponent,     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned**            componentblocked,   /**< array to store which symmetry methods have been used on a component
                                              *   using the same bitset information as for misc/usesymmetry */
   int*                  ncomponents         /**< pointer to store number of components of symmetry group */
   )
{
   SCIP_DISJOINTSET* componentstovar = NULL;
   int* permtovarcomp;
   int* permtocomponent;
   int p;
   int i;
   int idx;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( perms != NULL );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( vartocomponent != NULL );
   assert( componentblocked != NULL );
   assert( ncomponents != NULL );

   if ( nperms <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdisjointsetCreate(&componentstovar, SCIPblkmem(scip), npermvars) );
   *ncomponents = npermvars;

   /* init array that stores for each permutation the representative of its affected variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &permtovarcomp, nperms) );
   for (p = 0; p < nperms; ++p)
      permtovarcomp[p] = -1;

   /* find permutation components and store for each variable an affecting permutation (or -1)  */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, vartocomponent, npermvars) );
   for (i = 0; i < npermvars; ++i)
   {
      (*vartocomponent)[i] = -1;

      for (p = 0; p < nperms; ++p)
      {
         int img;

         img = transposed ? perms[i][p] : perms[p][i];

         /* perm p affects i -> possibly merge var components */
         if ( img != i )
         {
            int component1;
            int component2;
            int representative;

            component1 = SCIPdisjointsetFind(componentstovar, i);
            component2 = SCIPdisjointsetFind(componentstovar, img);
            (*vartocomponent)[i] = p;
            (*vartocomponent)[img] = p;

            /* ensure component1 <= component2 */
            if ( component2 < component1 )
            {
               int swap;

               swap = component1;
               component1 = component2;
               component2 = swap;
            }

            /* init permtovarcomp[p] to component of first moved variable or update the value */
            if ( permtovarcomp[p] == -1 )
            {
               permtovarcomp[p] = component1;
               representative = component1;
            }
            else
            {
               permtovarcomp[p] = SCIPdisjointsetFind(componentstovar, permtovarcomp[p]);
               representative = permtovarcomp[p];
            }

            /* merge both components if they differ */
            if ( component1 != component2 )
            {
               SCIPdisjointsetUnion(componentstovar, component1, component2, TRUE);
               --(*ncomponents);
            }

            /* possibly merge new component and permvartocom[p] and ensure the latter
             * to have the smallest value */
            if ( representative != component1 && representative != component2 )
            {
               if ( representative > component1 )
               {
                  SCIPdisjointsetUnion(componentstovar, component1, representative, TRUE);
                  permtovarcomp[p] = component1;
               }
               else
                  SCIPdisjointsetUnion(componentstovar, representative, component1, TRUE);
               --(*ncomponents);
            }
            else if ( representative > component1 )
            {
               assert( representative == component2 );
               permtovarcomp[p] = component1;
            }
         }
      }

      /* reduce number of components by singletons */
      if ( (*vartocomponent)[i] == -1 )
         --(*ncomponents);
   }
   assert( *ncomponents > 0 );

   /* update permvartocomp array to final variable representatives */
   for (p = 0; p < nperms; ++p)
      permtovarcomp[p] = SCIPdisjointsetFind(componentstovar, permtovarcomp[p]);

   /* init components array by trivial natural order of permutations */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, components, nperms) );
   for (p = 0; p < nperms; ++p)
      (*components)[p] = p;

   /* get correct order of components array */
   SCIPsortIntInt(permtovarcomp, *components, nperms);

   /* determine componentbegins and store components for each permutation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, componentbegins, *ncomponents + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permtocomponent, nperms) );

   (*componentbegins)[0] = 0;
   permtocomponent[(*components)[0]] = 0;
   idx = 0;

   for (p = 1; p < nperms; ++p)
   {
      if ( permtovarcomp[p] > permtovarcomp[p - 1] )
         (*componentbegins)[++idx] = p;

      assert( (*components)[p] >= 0 );
      assert( (*components)[p] < nperms );
      permtocomponent[(*components)[p]] = idx;
   }
   assert( *ncomponents == idx + 1 );
   (*componentbegins)[++idx] = nperms;

   /* determine vartocomponent */
   for (i = 0; i < npermvars; ++i)
   {
      int permidx;
      permidx = (*vartocomponent)[i];
      assert( -1 <= permidx && permidx < nperms );

      if ( permidx != -1 )
      {
         assert( 0 <= permtocomponent[permidx] );
         assert( permtocomponent[permidx] < *ncomponents );

         (*vartocomponent)[i] = permtocomponent[permidx];
      }
   }

   /* init componentblocked */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, componentblocked, *ncomponents) );
   for (i = 0; i < *ncomponents; ++i)
      (*componentblocked)[i] = 0;

   SCIPfreeBufferArray(scip, &permtocomponent);
   SCIPfreeBufferArray(scip, &permtovarcomp);
   SCIPdisjointsetFree(&componentstovar, SCIPblkmem(scip));

#ifdef SCIP_OUTPUT
   printf("number of components: %d\n", propdata->ncomponents);
   for (i = 0; i < *ncomponents; ++i)
   {
      printf("Component %d contains the following permutations:\n\t", i);
      for (p = (*componentbegins)[i]; p < (*componentbegins)[i + 1]; ++p)
      {
         printf("%d, ", (*components)[p]);
      }
      printf("\n");
   }
#endif

   return SCIP_OKAY;
}


/** generate variable matrix for orbitope constraint handler
 *
 * @pre if storelexorder is TRUE, then the permutations define an orbitope
 */
SCIP_RETCODE SCIPgenerateOrbitopeVarsMatrix(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR****          vars,               /**< pointer to matrix of orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder,        /**< permutation to reorder column of orbitopevaridx */
   int*                  nusedelems,         /**< array storing how often an element was used in the orbitope */
   SCIP_Shortbool*       rowisbinary,        /**< array encoding whether a row contains only binary variables (or NULL) */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the potential orbitope is not an orbitope */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing the lexorder (or NULL) */
   int*                  nvarsorder,         /**< pointer to store number of variables in lexorder (or NULL) */
   int*                  maxnvarsorder       /**< pointer to store maximum number of variables in lexorder (or NULL) */
   )
{
   int nfilledcols = 0;
   int curcolumn;
   int i;
   int cnt;
   int nvarsorderold = 0;

   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( orbitopevaridx != NULL );
   assert( columnorder != NULL );
   assert( nusedelems != NULL );
   assert( infeasible != NULL );
   assert( ! storelexorder || lexorder != NULL );
   assert( ! storelexorder || nvarsorder != NULL );
   assert( ! storelexorder || maxnvarsorder != NULL );

   /* possibly store lexicographic order defined by orbitope
    *
    * position (i,j) of orbitope has position nrows * j + i in lexicographic order
    */
   if ( storelexorder )
   {
      assert( *nvarsorder == *maxnvarsorder );
      assert( lexorder != NULL );

      *maxnvarsorder += nrows * ncols;
      nvarsorderold = *nvarsorder;

      if ( *lexorder == NULL )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, lexorder, *maxnvarsorder) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, lexorder, *nvarsorder, *maxnvarsorder) );
      }
   }

   curcolumn = ncols - 1;

   /* start filling vars matrix with the right-most column w.r.t. columnorder */
   while ( curcolumn >= 0 && columnorder[curcolumn] >= 0 && ! *infeasible )
   {
      cnt = 0;
      for (i = 0; i < nrows; ++i)
      {
         /* skip rows containing non-binary variables */
         if ( rowisbinary != NULL && ! rowisbinary[i] )
            continue;

         assert( 0 <= orbitopevaridx[i][curcolumn] && orbitopevaridx[i][curcolumn] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][curcolumn]]) );

         /* elements in first column of orbitope have to appear exactly once in the orbitope */
         if ( nfilledcols == 0 && nusedelems[orbitopevaridx[i][curcolumn]] > 1 )
         {
            *infeasible = TRUE;
            assert( ! storelexorder );
            break;
         }

         if ( storelexorder )
         {
            (*lexorder)[nvarsorderold + nrows * nfilledcols + cnt] = orbitopevaridx[i][curcolumn];
            ++(*nvarsorder);
         }
         (*vars)[cnt++][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
      }
      --curcolumn;
      ++nfilledcols;
   }

   /* There are three possibilities for the structure of columnorder:
    * 1)  [0, 1, -1, -1, ..., -1]
    * 2)  [0, 1, 1, 1, ..., 1]
    * 3)  [0, 1, -1, -1, ...., -1, 1, 1, ..., 1]
    */
   /* Either we are in case 1) or case 3), or all columns should have been added to vars in case 2) */
   assert( curcolumn > 1 || (curcolumn < 0 && nfilledcols == ncols) );

   if ( curcolumn > 1 && ! *infeasible )
   {
      /* add column with columnorder 1 to vars */
      cnt = 0;
      for (i = 0; i < nrows; ++i)
      {
         /* skip rows containing non-binary variables*/
         if ( rowisbinary != NULL && ! rowisbinary[i] )
            continue;

         assert( orbitopevaridx[i][1] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][1]]) );

         if ( storelexorder )
         {
            (*lexorder)[nvarsorderold + nrows * nfilledcols + cnt] = orbitopevaridx[i][1];
            ++(*nvarsorder);
         }
         (*vars)[cnt++][nfilledcols] = permvars[orbitopevaridx[i][1]];
      }
      ++nfilledcols;

      /* add column with columnorder 0 to vars */
      cnt = 0;
      for (i = 0; i < nrows; ++i)
      {
         /* skip rows containing non-binary variables*/
         if ( rowisbinary != NULL && ! rowisbinary[i] )
            continue;

         assert( orbitopevaridx[i][0] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][0]]) );

         if ( storelexorder )
         {
            (*lexorder)[nvarsorderold + nrows * nfilledcols + cnt] = orbitopevaridx[i][0];
            ++(*nvarsorder);
         }
         (*vars)[cnt++][nfilledcols] = permvars[orbitopevaridx[i][0]];
      }
      ++nfilledcols;

      /* add columns with a negative column order to vars */
      if ( nfilledcols < ncols )
      {
         assert( ncols > 2 );

         curcolumn = 2;
         while ( nfilledcols < ncols && ! *infeasible )
         {
            assert( columnorder[curcolumn] < 0 );

            cnt = 0;
            for (i = 0; i < nrows; ++i)
            {
               /* skip rows containing non-binary variables*/
               if ( rowisbinary != NULL && ! rowisbinary[i] )
                  continue;

               assert( orbitopevaridx[i][curcolumn] < npermvars );
               assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][curcolumn]]) );

               /* elements in last column of orbitope have to appear exactly once in the orbitope */
               if ( nfilledcols == ncols - 1 && nusedelems[orbitopevaridx[i][curcolumn]] > 1 )
               {
                  *infeasible = TRUE;
                  assert( ! storelexorder );
                  break;
               }

               if ( storelexorder )
               {
                  (*lexorder)[nvarsorderold + nrows * nfilledcols + cnt] = orbitopevaridx[i][curcolumn];
                  ++(*nvarsorder);
               }
               (*vars)[cnt++][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
            }
            ++curcolumn;
            ++nfilledcols;
         }
      }
   }

   return SCIP_OKAY;
}


/** checks whether an orbitope is a packing or partitioning orbitope; if npprows != NULL,
 *  count how many rows are contained in packing/partitioning constraints
 */
SCIP_RETCODE SCIPisPackingPartitioningOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variable matrix of orbitope constraint */
   int                   nrows,              /**< pointer to number of rows of variable matrix */
   int                   ncols,              /**< number of columns of variable matrix */
   SCIP_Bool**           pprows,             /**< pointer to store which rows are are contained in
                                              *   packing/partitioning constraints or NULL if not needed */
   int*                  npprows,            /**< pointer to store how many rows are contained
                                              *   in packing/partitioning constraints or NULL if not needed */
   SCIP_ORBITOPETYPE*    type                /**< pointer to store type of orbitope constraint after strengthening */
   )
{
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   int nsetppcconss;
   int* covered;
   int nprobvars;
   int* rowidxvar;
   int* rowcoveragesetppc;
   int* rowsinsetppc;
   int ncovered;
   int ncoveredpart;
   int i;
   int j;
   int c;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( type != NULL );

   *type = SCIP_ORBITOPETYPE_FULL;
   if ( npprows != NULL )
      *npprows = 0;

   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
      return SCIP_OKAY;

   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);

   /* we can terminate early if there are not sufficiently many setppc conss
    * (for orbitopes treating a full component, we might allow to remove rows
    * not contained in setppc cons; for this reason we need the second check)
    */
   if ( nsetppcconss == 0 || (nsetppcconss < nrows && npprows == NULL ))
      return SCIP_OKAY;
   assert( setppcconss != NULL );

   /* whether a row is contained in packing/partitioning constraint */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, nrows) );
   ncovered = 0;
   ncoveredpart = 0;

   nprobvars = SCIPgetNTotalVars(scip);

   /* array storing index of orbitope row a variable is contained in */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowidxvar, nprobvars) );

   for (i = 0; i < nprobvars; ++i)
      rowidxvar[i] = -1;

   for (i = 0; i < nrows; ++i)
   {
      for (j = 0; j < ncols; ++j)
      {
         assert( 0 <= SCIPvarGetIndex(vars[i][j]) && SCIPvarGetIndex(vars[i][j]) < nprobvars );
         rowidxvar[SCIPvarGetIndex(vars[i][j])] = i;
      }
   }

   /* storage for number of vars per row that are contained in current setppc cons and
    * labels of rows intersecting with current setppc cons
    */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &rowcoveragesetppc, nrows) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &rowsinsetppc, nrows) );

   /* iterate over set packing and partitioning constraints and check whether the constraint's
    * support is a row r of the orbitope (covered[r] = 2) or contains row r (covered[r] = 1)
    *
    * @todo Check whether we can improve the following loop by using a hash value to check
    * whether the setppccons intersects the orbitope matrix
    */
   for (c = 0; c < nsetppcconss && ncoveredpart < ncols; ++c)
   {
      int nsetppcvars;
      SCIP_VAR** setppcvars;
      int nrowintersect = 0;
      int nvarsinorbitope;

      /* skip covering constraints */
      if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_COVERING )
         continue;

      /* get number of set packing/partitioning variables */
      nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);

      /* constraint does not contain enough variables */
      if ( nsetppcvars < ncols )
         continue;

      setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
      assert( setppcvars != NULL );

      /* upper bound on variables potentially contained in orbitope */
      nvarsinorbitope = nsetppcvars;

      /* for each setppc var, check whether it appears in a row of the orbitope and store
       * for each row the number of such variables; can be terminated early, if less than
       * ncols variables are contained in the orbitope
       */
      for (i = 0; i < nsetppcvars && nvarsinorbitope >= ncols; ++i)
      {
         SCIP_VAR* var;
         int idx;
         int rowidx;

         var = setppcvars[i];
         idx = SCIPvarGetIndex(var);

         assert( 0 <= idx && idx < nprobvars );

         rowidx = rowidxvar[idx];

         /* skip variables not contained in the orbitope */
         if ( rowidx < 0 )
         {
            --nvarsinorbitope;
            continue;
         }

         /* skip variables corresponding to already treated rows */
         if ( covered[rowidx] == 2 || (covered[rowidx] == 1 && (nsetppcvars > ncols || nrowintersect > 1)) )
         {
            --nvarsinorbitope;
            continue;
         }

         /* store information which rows intersect the setppc cons's support */
         if ( rowcoveragesetppc[rowidx] == 0 )
            rowsinsetppc[nrowintersect++] = rowidx;
         ++(rowcoveragesetppc[rowidx]);

         /* we can stop early if not enough variables are left to completely cover one of the rows that
          * intersect the setppc cons
          */
         if ( nsetppcvars - nrowintersect < ncols - 1 )
            break;
      }

      /* store whether rows coincide with set partitioning cons's support or whether
       * row is covered by a set packing/partitioning cons's support
       */
      if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PARTITIONING
           && nrowintersect == 1 && rowcoveragesetppc[rowsinsetppc[0]] == ncols && nsetppcvars == ncols )
      {
         if ( covered[rowsinsetppc[0]] == 1 )
            --ncovered;
         covered[rowsinsetppc[0]] = 2;
         ++ncoveredpart;
         ++ncovered;
      }
      else
      {
         for (i = 0; i < nrowintersect; ++i)
         {
            if ( covered[rowsinsetppc[i]] == 0 && rowcoveragesetppc[rowsinsetppc[i]] >= ncols )
            {
               covered[rowsinsetppc[i]] = 1;
               ++ncovered;
            }
         }
      }

      /* reset data */
      for (i = 0; i < nrowintersect; ++i)
         rowcoveragesetppc[rowsinsetppc[i]] = 0;
   }

   /* check type of orbitope */
   if ( ncovered == nrows )
   {
      if ( ncoveredpart == nrows )
         *type = SCIP_ORBITOPETYPE_PARTITIONING;
      else
         *type = SCIP_ORBITOPETYPE_PACKING;
   }

   if ( npprows != NULL )
      *npprows = ncovered;

   if ( pprows != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, pprows, nrows) );
      for (i = 0; i < nrows; ++i)
         (*pprows)[i] = covered[i] > 0 ? 1 : 0;
   }

   SCIPfreeBufferArray(scip, &rowsinsetppc);
   SCIPfreeCleanBufferArray(scip, &rowcoveragesetppc);
   SCIPfreeBufferArray(scip, &rowidxvar);
   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}

/** checks whether a (signed) permutation is an involution */
static
SCIP_RETCODE isPermIvolution(
   int*                  perm,               /**< permutation */
   int                   permlen,            /**< number of original (non-negated) variables in a permutation */
   SCIP_Bool             issignedperm,       /**< whether permutation is encoded as signed */
   SCIP_Bool*            permissigned,       /**< pointer to store whether perm is a signed permutation */
   SCIP_Bool*            isinvolution,       /**< pointer to store whether perm is an involution */
   int*                  ntwocycles          /**< pointer to store number of 2-cycles in an involution */
   )
{
   int v;

   assert( perm != NULL );
   assert( permlen > 0 );
   assert( permissigned != NULL );
   assert( isinvolution != NULL );
   assert( ntwocycles != NULL );

   *ntwocycles = 0;
   *permissigned = FALSE;
   *isinvolution = TRUE;
   for (v = 0; v < permlen && *isinvolution; ++v)
   {
      /* do not handle variables twice */
      if ( perm[v] <= v )
         continue;

      /* detect signed permutation */
      if ( perm[v] >= permlen )
         *permissigned = TRUE;

      /* detect two cycles */
      if ( perm[perm[v]] == v )
         ++(*ntwocycles);
      else
         *isinvolution = FALSE;
   }

   return SCIP_OKAY;
}

/* checks whether selected permutations define orbitopal symmetrie */
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
{
   SCIP_DISJOINTSET* conncomps;
   SCIP_DISJOINTSET* compcolors;
   int* permstoconsider;
   int* colorbegins;
   int* compidx;
   int* colidx;
   int* varidx;
   int* degrees;
   int* perm;
   int nposdegree = 0;
   int npermstoconsider;
   int colorrepresentative1;
   int colorrepresentative2;
   int elemtomove;
   int ncurcols;
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
   assert( nrows > 0 );
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

   /* try to add edges of permutations to graph */
   for (p = 0; p < nselectedperms; ++p)
   {
      perm = perms[selectedperms[p]];
      curdeg1 = -1;
      curdeg2 = -1;

      for (v = 0; v < permlen; ++v)
      {
         /* treat each cycle exactly once */
         if ( perm[v] <= v )
            continue;
         w = perm[v];

         /* get colors of nodes */
         curcolor1 = SCIPdisjointsetFind(compcolors, v);
         curcolor2 = SCIPdisjointsetFind(compcolors, w);

         /* an edge is not allowed to connect nodes with the same color */
         if ( curcolor1 == curcolor2 )
         {
            *success = FALSE;
            goto FREEMEMORY;
         }

         /* the first edge determines the degrees of the two sets of orbitopal symmetries */
         if ( curdeg1 == -1 )
         {
            assert( curdeg2 == -1 );

            curdeg1 = degrees[v];
            curdeg2 = degrees[w];
            colorrepresentative1 = v;
            colorrepresentative2 = w;

            /* stop, we will generate a vertex with degree 3 */
            if ( curdeg1 == 2 || curdeg2 == 2 )
            {
               *success = FALSE;
               goto FREEMEMORY;
            }
         }

         /* try to add edges and join color classes */
         if ( (curdeg1 == degrees[v] && curdeg2 == degrees[w]) || (curdeg1 == degrees[w] && curdeg2 == degrees[v]) )
         {
            /* check whether colors of different connected components are compatible */
            if ( (curdeg1 > 0 && (curcolor1 != colorrepresentative1 && curcolor1 != colorrepresentative2))
               || (curdeg2 > 0 && (curcolor2 != colorrepresentative1 && curcolor2 != colorrepresentative2)) )
            {
               /* the colors of the connected components to be joinded do not match */
               *success = FALSE;
               goto FREEMEMORY;
            }

            /* add edge for v and w */
            SCIPdisjointsetUnion(conncomps, v, w, FALSE);
            SCIPdisjointsetUnion(compcolors, colorrepresentative1, v, TRUE);
            SCIPdisjointsetUnion(compcolors, colorrepresentative1, w, TRUE);
            ++degrees[v];
            ++degrees[w];
         }
         else
         {
            /* the degrees of the nodes to be connected do not match among the different connected components */
            *success = FALSE;
            goto FREEMEMORY;
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

   /* find the correct order of variable indices per color class */
   SCIP_CALL( SCIPallocBufferArray(scip, &permstoconsider, nselectedperms) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, matrices, ncolors) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, ncols, ncolors) );
   *nmatrices = ncolors;

   for (c = 0; c < ncolors; ++c)
   {
      /* find an element in the first connected component with degree 1 */
      for (v = colorbegins[c]; compidx[v] == compidx[colorbegins[c]]; ++v)
      {
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
         for (v = colorbegins[c]; compidx[v] == compidx[colorbegins[c]]; ++v)
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
         if ( perm[v] > v )
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

   for (c = 0; c < ncolors; ++c)
   {
      for (v = 0; v < nrows; ++v)
      {
         for (w = 0; w < nselectedperms + 1; ++w)
            printf(" %4d", (*matrices)[c][v][w]);
         printf("\n");
      }
      printf("\n");
   }

   SCIPfreeBufferArray(scip, &permstoconsider);
   SCIPfreeBufferArray(scip, &colorbegins);
   SCIPfreeBufferArray(scip, &varidx);
   SCIPfreeBufferArray(scip, &colidx);
   SCIPfreeBufferArray(scip, &compidx);

 FREEMEMORY:
   SCIPfreeBufferArray(scip, &degrees);
   SCIPfreeDisjointset(scip, &compcolors);
   SCIPfreeDisjointset(scip, &conncomps);

   return SCIP_OKAY;
}

/** checks whether to families of orbitopal symmetries defines double lex matrix, in case of success, generate matrix
 *
 *  The columns of matrix1 will serve as the columns of the matrix to be generated, the columns of matrix2 will
 *  serve as rows.
 */
static
SCIP_RETCODE isDoublelLexSym(
   SCIP*                 scip,               /**< SCIP pointer */
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
   SCIP_HASHMAP* idxtomatrix1;
   SCIP_HASHMAP* idxtomatrix2;
   SCIP_HASHMAP* idxtorow1;
   SCIP_HASHMAP* idxtorow2;
   SCIP_HASHMAP* idxtocol1;
   SCIP_HASHMAP* idxtocol2;
   int* sortvals;
   int elem;
   int mat;
   int col;
   int col2;
   int mat2;
   int nidx;
   int cnt;
   int c;
   int d;
   int i;
   int j;

   assert( scip != NULL );
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
   assert (success != NULL );

   /* initialize data */
   *nrows = nrows1;
   *ncols = nrows2;
   rowsbegin = 0;
   colsbegin = 0;
   *success = TRUE;

   /* collect information about entries in matrices */
   nidx = nrows1 * nrows2;
   SCIP_CALL( SCIPhashmapCreate(&idxtomatrix1, SCIPblkmem(scip), nidx) );
   SCIP_CALL( SCIPhashmapCreate(&idxtomatrix2, SCIPblkmem(scip), nidx) );
   SCIP_CALL( SCIPhashmapCreate(&idxtorow1, SCIPblkmem(scip), nidx) );
   SCIP_CALL( SCIPhashmapCreate(&idxtorow2, SCIPblkmem(scip), nidx) );
   SCIP_CALL( SCIPhashmapCreate(&idxtocol1, SCIPblkmem(scip), nidx) );
   SCIP_CALL( SCIPhashmapCreate(&idxtocol2, SCIPblkmem(scip), nidx) );

   for (c = 0; c < nmatrices1; ++c)
   {
      for (i = 0; i < nrows1; ++i)
      {
         for (j = 0; j < ncols1[c]; ++j)
         {
            SCIP_CALL( SCIPhashmapInsertInt(idxtomatrix1, (void*) (size_t) matrices1[c][i][j], c) );
            SCIP_CALL( SCIPhashmapInsertInt(idxtorow1, (void*) (size_t) matrices1[c][i][j], i) );
            SCIP_CALL( SCIPhashmapInsertInt(idxtocol1, (void*) (size_t) matrices1[c][i][j], j) );
         }
      }
   }
   for (c = 0; c < nmatrices2; ++c)
   {
      for (i = 0; i < nrows2; ++i)
      {
         for (j = 0; j < ncols2[c]; ++j)
         {
            SCIP_CALL( SCIPhashmapInsertInt(idxtomatrix2, (void*) (size_t) matrices2[c][i][j], c) );
            SCIP_CALL( SCIPhashmapInsertInt(idxtorow2, (void*) (size_t) matrices2[c][i][j], i) );
            SCIP_CALL( SCIPhashmapInsertInt(idxtocol2, (void*) (size_t) matrices2[c][i][j], j) );
         }
      }
   }

#ifndef NDEBUG
   /* check whether expecteded sizes of matrix match */
   for (j = 0, cnt = 0; j < nmatrices1; ++j)
      cnt += ncols1[j];
   assert( cnt == *ncols );

   for (i = 0, cnt = 0; i < nmatrices2; ++i)
      cnt += ncols2[i];
   assert( cnt == *nrows );
#endif

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
      sortvals[i] = SCIPhashmapGetImageInt(idxtomatrix2, (void*) (size_t) matrices1[0][i][0]);
   }
   SCIPsortIntPtr(sortvals, (void*) (*doublelexmatrix), *nrows);

   /* fill first row of big matrix */
   mat = SCIPhashmapGetImageInt(idxtomatrix2, (void*) (size_t) (*doublelexmatrix)[0][0]);
   col = SCIPhashmapGetImageInt(idxtocol2, (void*) (size_t) (*doublelexmatrix)[0][0]);
   cnt = 0;
   for (j = 0; j < *ncols; ++j)
   {
      /* skip the entry that is already contained in the first column */
      if ( matrices2[mat][j][col] == (*doublelexmatrix)[0][0] )
         continue;

      sortvals[cnt++] = SCIPhashmapGetImageInt(idxtomatrix1, (void*) (size_t) matrices2[mat][j][col]);
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
         mat = SCIPhashmapGetImageInt(idxtomatrix1, (void*) (size_t) (*doublelexmatrix)[0][j]);
         mat2 = SCIPhashmapGetImageInt(idxtomatrix2, (void*) (size_t) (*doublelexmatrix)[i][0]);
         col = SCIPhashmapGetImageInt(idxtocol1, (void*) (size_t) (*doublelexmatrix)[0][j]);
         col2 = SCIPhashmapGetImageInt(idxtocol2, (void*) (size_t) (*doublelexmatrix)[i][0]);

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

   printf("Fill success %d\n", *success);
   for (i = 0; i < *nrows; ++i)
   {
      for (j = 0; j < *ncols; ++j)
         printf(" %4d", (*doublelexmatrix)[i][j]);
      printf("\n");
   }
   printf("\n");

 FREEMEMORY:
   SCIPfreeBufferArray(scip, &sortvals);
   SCIPhashmapFree(&idxtocol2);
   SCIPhashmapFree(&idxtocol1);
   SCIPhashmapFree(&idxtorow2);
   SCIPhashmapFree(&idxtorow1);
   SCIPhashmapFree(&idxtomatrix2);
   SCIPhashmapFree(&idxtomatrix1);

   for (i = *nrows - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*doublelexmatrix)[i], *ncols);
   }
   SCIPfreeBlockMemoryArray(scip, doublelexmatrix, *nrows);

   return SCIP_OKAY;
}

/** tries to handle variable matrices with lex ordered rows and columns */
SCIP_RETCODE tryHandleDoubleLexMatrices(
   SCIP*                 scip,               /**< SCIP pointer */
   int**                 perms,              /**< array of permutations */
   int                   nperms,             /**< number of permutations in perms */
   int                   permlen,            /**< number of original (non-negated) variables in a permutation */
   SCIP_Bool             issignedperm        /**< whether permutations are encoded as signed */
   )
{
   SCIP_Bool success = TRUE;
   SCIP_Bool isinvolution;
   SCIP_Bool permissigned;
   int** doublelexmatrix;
   int* rowsbegin;
   int* colsbegin;
   int nrows;
   int ncols;
   int*** matricestype1;
   int* ncolstype1;
   int nmatricestype1;
   int*** matricestype2;
   int* ncolstype2;
   int nmatricestype2;
   int* permstype1;
   int* permstype2;
   int* signedperms;
   int npermstype1 = 0;
   int npermstype2 = 0;
   int nsignedperms = 0;
   int ncycs1 = -1;
   int ncycs2 = -1;
   int tmpncycs;
   int p;
   int i;

   assert( scip != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( permlen > 0 );

   /* arrays to store the different types of involutions */
   SCIP_CALL( SCIPallocBufferArray(scip, &permstype1, nperms)  );
   SCIP_CALL( SCIPallocBufferArray(scip, &permstype2, nperms)  );
   SCIP_CALL( SCIPallocBufferArray(scip, &signedperms, nperms)  );

   /* check whether we can expect lexicographically sorted rows and columns */
   for (p = 0; p < nperms; ++p)
   {
      SCIP_CALL( isPermIvolution(perms[p], permlen, issignedperm, &permissigned, &isinvolution, &tmpncycs) );

      /* ignore signed permutation for the time being */
      if ( permissigned )
      {
         signedperms[nsignedperms++] = p;
         continue;
      }

      /* terminate if not all permutations are involutions */
      if ( ! isinvolution )
      {
         success = FALSE;
         break;
      }

      /* store number of cycles or termintate if too many different types of involutions */
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
         success = FALSE;
         goto FREEMEMORY;
      }
   }

   /* for each type, check whether permutations define (disjoint) orbitopal symmetries; ignore signed part */
   SCIP_CALL( detectOrbitopalSymmetries(scip, perms, permstype1, npermstype1, permlen, ncycs1, &success,
         &matricestype1, &ncolstype1, &nmatricestype1) );
   if ( ! success )
      goto FREEMEMORY;

   SCIP_CALL( detectOrbitopalSymmetries(scip, perms, permstype2, npermstype2, permlen, ncycs2, &success,
         &matricestype2, &ncolstype2, &nmatricestype2) );
   if ( ! success )
      goto FREEMEMORY;

   /* check whether symmetries of type 2 permute two rows of matrix of type 1 */
   SCIP_CALL( isDoublelLexSym(scip, matricestype1, ncycs1, ncolstype1, nmatricestype1,
         matricestype2, ncycs2, ncolstype2, nmatricestype2,
         &doublelexmatrix, &nrows, &ncols, &rowsbegin, &colsbegin, &success) );

   /* check whether proper signed permutations flip rows or columns */

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

   SCIPfreeBufferArray(scip, &signedperms);
   SCIPfreeBufferArray(scip, &permstype2);
   SCIPfreeBufferArray(scip, &permstype1);

   printf("lex sort detection success: %d\n", success);

   return SCIP_OKAY;
}
