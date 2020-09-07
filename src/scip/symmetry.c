/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
#include "scip/misc.h"


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
   SCIP_Shortbool*       componentblocked,   /**< array to store whether a component is blocked to be considered by
                                              *   further symmetry handling techniques */
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


/** check whether a permutation is a composition of 2-cycles of binary variables and in this case determine the number of 2-cycles */
SCIP_RETCODE SCIPgetPropertiesPerm(
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< array of variables perm is acting on */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            iscompoftwocycles,  /**< pointer to store whether permutation is a composition of 2-cycles */
   int*                  ntwocyclesperm,     /**< pointer to store number of 2-cycles */
   SCIP_Bool*            allvarsbinary       /**< pointer to store whether perm is acting on binary variables only */
   )
{
   int ntwocycles = 0;
   int i;

   assert( perm != NULL );
   assert( vars != NULL );
   assert( iscompoftwocycles != NULL );
   assert( ntwocyclesperm != NULL );
   assert( allvarsbinary != NULL );

   *iscompoftwocycles = FALSE;
   *ntwocyclesperm = 0;
   *allvarsbinary = TRUE;
   for (i = 0; i < nvars; ++i)
   {
      /* skip fixed points and avoid treating the same 2-cycle twice */
      if ( perm[i] <= i )
         continue;

      if ( perm[perm[i]] == i )
      {
         if ( SCIPvarIsBinary(vars[i]) && SCIPvarIsBinary(vars[perm[i]]) )
            ++ntwocycles;
         else
         {
            /* at least one variable is not binary */
            *allvarsbinary = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* we do not have a 2-cycle */
         return SCIP_OKAY;
      }
   }

   /* at this point the permutation is a composition of 2-cycles on binary variables */
   *iscompoftwocycles = TRUE;
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
            suborbitope[row][2] = perm[idx1];
            ++nintersections;

            (*nusedelems)[idx1] += 1;
            (*nusedelems)[perm[idx1]] += 1;

            /* if an element appears too often in the orbitope matrix */
            if ( (*nusedelems)[idx1] > 2 || (*nusedelems)[perm[idx1]] > 2 )
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
            suborbitope[row][2] = perm[idx2];
            ++nintersections;

            (*nusedelems)[idx2] += 1;
            (*nusedelems)[perm[idx2]] += 1;

            /* if an element appears too often in the orbitope matrix */
            if ( (*nusedelems)[idx2] > 2 || (*nusedelems)[perm[idx2]] > 2 )
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
            suborbitope[row][nfilledcols] = perm[idx1];
            ++nintersections;

            (*nusedelems)[idx1] += 1;
            (*nusedelems)[perm[idx1]] += 1;

            /* if an element appears to often in the orbitope matrix */
            if ( (*nusedelems)[idx1] > 2 || (*nusedelems)[perm[idx1]] > 2 )
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
   SCIP_Shortbool**      componentblocked,   /**< array to store whether a component is blocked to be considered by
                                              *   further symmetry handling techniques */
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
      (*componentblocked)[i] = FALSE;

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


/** generate variable matrix for orbitope constraint handler */
SCIP_RETCODE SCIPgenerateOrbitopeVarsMatrix(
   SCIP_VAR****          vars,               /**< pointer to matrix of orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder,        /**< permutation to reorder column of orbitopevaridx */
   int*                  nusedelems,         /**< array storing how often an element was used in the orbitope */
   SCIP_Bool*            infeasible          /**< pointer to store whether the potential orbitope is not an orbitope */
   )
{
   int nfilledcols = 0;
   int curcolumn;
   int i;

   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( orbitopevaridx != NULL );
   assert( columnorder != NULL );
   assert( nusedelems != NULL );
   assert( infeasible != NULL );

   curcolumn = ncols - 1;

   /* start filling vars matrix with the right-most column w.r.t. columnorder */
   while ( curcolumn >= 0 && columnorder[curcolumn] >= 0 )
   {
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][curcolumn] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][curcolumn]]) );

         /* elements in first column of orbitope have to appear exactly once in the orbitope */
         if ( nfilledcols == 0 && nusedelems[orbitopevaridx[i][curcolumn]] > 1 )
         {
            *infeasible = TRUE;
            break;
         }

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
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

   if ( curcolumn > 1 )
   {
      /* add column with columnorder 1 to vars */
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][1] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][1]]) );

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][1]];
      }
      ++nfilledcols;

      /* add column with columnorder 0 to vars */
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][0] < npermvars );
         assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][0]]) );

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][0]];
      }
      ++nfilledcols;

      /* add columns with a negative column order to vars */
      if ( nfilledcols < ncols )
      {
         assert( ncols > 2 );

         curcolumn = 2;
         while ( nfilledcols < ncols )
         {
            assert( columnorder[curcolumn] < 0 );

            for (i = 0; i < nrows; ++i)
            {
               assert( orbitopevaridx[i][curcolumn] < npermvars );
               assert( SCIPvarIsBinary(permvars[orbitopevaridx[i][curcolumn]]) );

               /* elements in last column of orbitope have to appear exactly once in the orbitope */
               if ( nfilledcols == ncols - 1 && nusedelems[orbitopevaridx[i][curcolumn]] > 1 )
               {
                  *infeasible = TRUE;
                  break;
               }

               (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
            }
            ++curcolumn;
            ++nfilledcols;
         }
      }
   }

   return SCIP_OKAY;
}
