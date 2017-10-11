/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_symbreak.cpp
 * @brief  presolver for adding symmetry breaking constraints
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 *
 * This presolver adds the following symmetry breaking constraints:
 *
 * - minimal cover inequalities for symresacks via a constraint handler
 *
 * @note It is important to control the order of presolvers within SCIP in order to avoid contraditions. First, one needs
 * to take care of presolvers that have an effect on symmetry, e.g., presol_domcol. If symmetry is computed on the
 * original formulation, we perform this presolver at the very beginning. Otherwise, we try to compute symmetry as late
 * as possible and then add constraints based on this information.
 *
 * @note Currently, we only allocate memory for pointers to symresack constraints for group generators. If further
 * constraints are considered, we have to reallocate memory.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include <scip/cons_linear.h>
#include <scip/cons_setppc.h>
#include <scip/cons_symresack.h>
#include <scip/misc.h>
#include <scip/prop_probing.h>
#include <scip/presol_symbreak.h>
#include <scip/presol_symmetry.h>
#include <symmetry/type_symmetry.h>

/* presolver properties */
#define PRESOL_NAME            "symbreak"
#define PRESOL_DESC            "presolver for adding symmetry breaking constraints"
#define PRESOL_PRIORITY         -10000000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                 TRUE    /**< should presolver be delayed, if other presolvers found reductions? */
#define PRESOL_TIMING   SCIP_PRESOLTIMING_EXHAUSTIVE   /**< timing for presolving */

/* default parameters */
#define DEFAULT_CONSSADDLP           TRUE    /**< Should the symmetry breaking constraints be added to the LP? */
#define DEFAULT_ADDSYMRESACKS        TRUE    /**< Add inequalities for symresacks for each generator? */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int**                 perms;              /**< list of permutations */
   int                   nperms;             /**< number of permutations in perms */
   int                   npermvars;          /**< number of variables affected by permutations */
   SCIP_VAR**            permvars;           /**< array of variables on which permutations act */
   SCIP_Bool             addedconss;         /**< whether we already added symmetry breaking constraints */
   SCIP_Bool             computedsymmetry;   /**< whether symmetry has been computed already */
   SCIP_Bool             conssaddlp;         /**< Should the symmetry breaking constraints be added to the LP? */
   SCIP_Bool             addsymresacks;      /**< Add symresack constraints for each generator? */
   SCIP_PRESOL*          symmetrypresol;     /**< pointer to symmetry presolver */
   SCIP_Bool             enabled;            /**< run presolver? */
   SCIP_Bool             early;              /**< run presolver as early as possible if symmetry has been detected in initpre() */
   SCIP_CONS**           genconss;           /**< list of generated constraints */
   int                   ngenconss;          /**< number of generated constraints */
   int                   nsymresacks;        /**< number of symresack constraints */
   int                   norbits;            /**< number of non-trivial orbits of permutation group */
   int*                  nvarsinorbits;      /**< number of variables per orbit */
   int**                 orbits;             /**< array containing for each orbit an array with the variable indices in the orbit */
   int                   ncomponents;        /**< number of components of symmetry group */
   int*                  npermsincomponent;  /**< array containing number of permutations per component */
   int**                 components;         /**< array containing for each components the corresponding permutations */
};


/*
 * Local methods
 */

/** compute orbit of a variabe; the method MUST NOT be called if no symmetries were found
 * or symmetries have not been computed yet */
static
SCIP_RETCODE computeOrbitVariable(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata,         /**< memory address of data of symmetry breaking presolver */
   int**                 orbit,              /**< memory address of preliminary orbits array */
   int                   i                   /**< index of variable for which the orbit should be computed */
   )
{
   int** perms;
   int nperms;
   int npermvars;
   int* curorbit;
   int curorbitsize;
   SCIP_Bool* varadded;
   int j;
   int p;
   int curelem;
   int image;
   int norbits;

   assert( orbit != NULL );
   assert( i >= 0 );

   nperms = presoldata->nperms;
   assert( nperms > 0 );

   perms = presoldata->perms;
   npermvars = presoldata->npermvars;

   /* initialize orbit of variable i */
   SCIP_CALL( SCIPallocBufferArray(scip, &curorbit, npermvars - i) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varadded, npermvars - i) ); /* TODO: can we do better by using a map that says whether we have already added a variable to curorbit? */

   for (j = 0; j < npermvars - i; ++j)
      varadded[j] = FALSE;

   curorbit[0] = i;
   curorbitsize = 1;
   varadded[0] = TRUE; /* theoretically, this is varadded[i - i] */

   /* iterate over variables in curorbit and compute their images */
   for (j = 0; j < curorbitsize; ++j)
   {
      curelem = curorbit[j];

      for (p = 0; p < nperms; ++p)
      {
         image = perms[p][curelem];

         /* found new element of the orbit of i */
         if ( ! varadded[image - i] )
         {
            curorbit[curorbitsize++] = image;
            varadded[image - i] = TRUE;

            (*orbit)[image] = i;
         }
      }
   }

   /* orbit is non-trivial -> store it */
   if ( curorbitsize > 1 )
   {
      norbits =  presoldata->norbits;

      if ( norbits == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata)->nvarsinorbits, 1) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata)->orbits, 1) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata)->nvarsinorbits, norbits, norbits + 1) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata)->orbits, norbits, norbits + 1) );
      }

      presoldata->nvarsinorbits[norbits] = curorbitsize;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata)->orbits[norbits], curorbitsize) );
      for (i = 0; i < curorbitsize; ++i)
         presoldata->orbits[norbits][i] = curorbit[i];

      presoldata->norbits = norbits + 1;
   }

   SCIPfreeBufferArray(scip, &varadded);
   SCIPfreeBufferArray(scip, &curorbit);

   return SCIP_OKAY;
}


/** compute orbits of symmetry group; the method MUST NOT be called if no symmetries were found
 * or symmetries have not been computed yet */
static
SCIP_RETCODE computeGroupOrbits(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< pointer to data of symmetry breaking presolver */
   )
{
   int npermvars;
   int* orbit;
   int i;

   assert( scip != NULL );
   assert( presoldata != NULL );

   npermvars = presoldata->npermvars;

   assert( npermvars > 0 );

   /* init data structures*/
   SCIP_CALL( SCIPallocBufferArray(scip, &orbit, npermvars) );

   /* initially, every variable is contained in its own orbit */
   for (i = 0; i < npermvars; ++i)
      orbit[i] = i;

   /* find variable orbits */
   presoldata->norbits = 0;

   for (i = 0; i < npermvars; ++i)
   {
      /* if variable is not contained in an orbit of a previous variable */
      if ( orbit[i] == i )
      {
         /* compute and store orbit */
         SCIP_CALL( computeOrbitVariable(scip, presoldata, &orbit, i) );
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &orbit);

#if 0
   printf("\n\n\nTESTS\n\n");
   printf("Number of orbits:\t\t%d\n", presoldata->norbits);
   for (i = 0; i < presoldata->norbits; ++i)
   {
      int j;
      printf("Orbit %d: Number of variables: %d\n", i, presoldata->nvarsinorbits[i]);
      for (j = 0; j < presoldata->nvarsinorbits[i]; ++j)
      {
         printf("%d ", presoldata->orbits[i][j]);
      }
      printf("\n");
   }
#endif

   return SCIP_OKAY;
}


/** compute components of symmetry group */
static
SCIP_RETCODE computeComponents(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< data of symmetry breaking presolver */
   )
{
   int npermvars;
   int nperms;
   int** perms;
   int i;
   int j;
   int k;
   int npermsincomponent;
   int curcomponent;
   int ncomponents;
   int componentcnt;
   int start;
   int newstart;
   SCIP_Bool startchanged;
   SCIP_DISJOINTSET* componentstoperm;

   assert( presoldata != NULL );

   nperms = presoldata->nperms;

   presoldata->ncomponents = -1;

   if ( nperms <= 0 )
      return SCIP_OKAY;

   /* at this point, we have at least one non-trivial permutation */
   npermvars = presoldata->npermvars;
   perms = presoldata->perms;

   /* init array that assigns each permutation its component of the group */
   SCIP_CALL( SCIPdisjointsetCreate(&componentstoperm, SCIPblkmem(scip), nperms) );
   ncomponents = nperms;

   /* check whether two permutations belong to the same component */
   for (i = 0; i < nperms; ++i)
   {
      for (j = i + 1; j < nperms; ++j)
      {
         int componentI = SCIPdisjointsetFind(componentstoperm, i);
         int componentJ = SCIPdisjointsetFind(componentstoperm, j);
         if ( componentI == componentJ )
            continue;

         /* belong perms[i] and perms[j] to the same component? */
         for (k = 0; k < npermvars; ++k)
         {
            /* both permutations belong to the same component */
            if ( perms[i][k] != k && perms[j][k] != k )
            {
               /* keep the smallest identifier to keep track of where a new component starts */
               if ( componentI < componentJ )
                  SCIPdisjointsetUnion(componentstoperm, i, j, TRUE);
               else
                  SCIPdisjointsetUnion(componentstoperm, j, i, TRUE);

               --ncomponents;
               break;
            }
         }
      }
   }
   assert( ncomponents > 0 );

   /* store data */
   presoldata->ncomponents = ncomponents;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->npermsincomponent, ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->components, ncomponents) );

   /* copy componentstoperm adequatly to components and npermsincomponent*/
   componentcnt = 0;
   start = 0;
   newstart = 0;
   startchanged = TRUE;
   while ( startchanged )
   {
      startchanged = FALSE;
      npermsincomponent = 0;
      curcomponent = SCIPdisjointsetFind(componentstoperm, start);

      /* find number of permutations in current component and detect first perm in another permutation */
      for (i = start; i < nperms; ++i)
      {
         if ( SCIPdisjointsetFind(componentstoperm, i) == curcomponent )
            ++npermsincomponent;
         /* only store other component if it has the same identifier as its node (mark that new component starts) */
         else if ( ! startchanged && SCIPdisjointsetFind(componentstoperm, i) == i )
         {
            newstart = i;
            startchanged = TRUE;
         }
      }

      /* store number of permutations and permutation per component */
      presoldata->npermsincomponent[componentcnt] = npermsincomponent;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->components[componentcnt]), npermsincomponent) );

      j = 0;
      for (i = start; i < nperms; ++i)
      {
         if ( SCIPdisjointsetFind(componentstoperm, i) == curcomponent )
            presoldata->components[componentcnt][j++] = i;
      }

      ++componentcnt;
      start = newstart;
   }

   SCIPdisjointsetFree(&componentstoperm, SCIPblkmem(scip));

#if 0
   printf("\n\n\nTESTS\n\n");
   printf("Number of components:\t\t%d\n", presoldata->ncomponents);
   for (i = 0; i < presoldata->ncomponents; ++i)
   {
      printf("Component %d: Number of perms: %d\n", i, presoldata->npermsincomponent[i]);
      for (j = 0; j < presoldata->npermsincomponent[i]; ++j)
      {
         printf("%d ", presoldata->components[i][j]);
      }
      printf("\n");
   }

   printf("GENERATORS:");
   for (i = 0; i < presoldata->nperms; ++i)
   {
      printf("generator %d:\n", i);
      for (j = 0; j < presoldata->npermvars; ++j)
      {
         if ( presoldata->perms[i][j] != j )
            printf("%d ", presoldata->perms[i][j]);
      }
      printf("\n");
   }
#endif

   return SCIP_OKAY;
}


/** check whether a permutation is a composition of 2-cycles and determines in this case
    the number of 2-cycles */
static
SCIP_RETCODE getPermProperties(
   int*                  perm,               /**< permutation */
   int                   n,                  /**< length of permutation array */
   SCIP_Bool*            iscompoftwocycles,  /**< memory address to store whether permutation is
                                                  a composition of 2-cycles */
   int*                  ntwocyclesperm      /**< memory address to store number of 2-cycles */
   )
{
   int i;
   int ntwocycles = 0;

   assert( perm != NULL );
   assert( iscompoftwocycles != NULL );
   assert( ntwocyclesperm != NULL );

   *iscompoftwocycles = TRUE;
   *ntwocyclesperm = 0;
   for (i = 0; i < n; ++i)
   {
      /* skip fixed points and avoid treating the same 2-cycle twice */
      if ( perm[i] <= i )
         continue;

      if ( perm[perm[i]] == i )
         ++ntwocycles;
      else
         break;
   }

   if ( i < n )
      *iscompoftwocycles = FALSE;
   else
      *ntwocyclesperm = ntwocycles;

   return SCIP_OKAY;
}


/** Given a matrix with nrows and #perms + 1 columns whose first
 *  nfilledcols contain entries of variables, this routine checks
 *  whether the 2-cycles of perm intersect each row of column
 *  coltoextend in exactly one position. In this case, we add one
 *  column to the sub orbitope of the first nfilledcols columns.
 *
 *  @pre Every non-trivial cycle of perm is a 2-cycle.
 *  @pre perm has nrows many 2-cycles
 */
static
SCIP_RETCODE extendSubOrbitope(
   int***                suborbitope,        /**< memory address of matrix containing sub orbitope entries */
   int                   nrows,              /**< number of rows of suborbitope */
   int                   nfilledcols,        /**< number of columns of suborbitope which are filled with entries */
   int                   coltoextend,        /**< index of column that should be extended by perm */
   int*                  perm,               /**< permutation */
   SCIP_Bool             leftextension,      /**< whether we extend the sub orbitope to the left */
   SCIP_Bool*            success,            /**< memory address to store whether extension was successful */
   SCIP_Bool*            infeasible          /**< memory address to store if the number of intersecting cycles is too small */
   )
{
   int row;
   int idx1;
   int idx2;
   int nintersections = 0;

   assert( suborbitope != NULL );
   assert( nfilledcols > 0 );
   assert( perm != NULL );
   assert( success != NULL );

   *success = FALSE;
   *infeasible = FALSE;

   /* if we try to extend the first orbitope generator by another one,
    * we can change the order of entries in each row of suborbitope */
   if ( nfilledcols == 2 )
   {
      /* check whether each cycle of perm intersects with a row of suborbitope */
      for (row = 0; row < nrows; ++row)
      {
         idx1 = (*suborbitope)[row][0];
         idx2 = (*suborbitope)[row][1];

         /* if idx1 or idx2 is affected by perm, we can extend the row of the orbitope */
         if ( idx1 != perm[idx1] )
         {
            /* change order of idx1 and idx2 for right extensions */
            if ( ! leftextension )
            {
               (*suborbitope)[row][0] = idx2;
               (*suborbitope)[row][1] = idx1;
            }
            (*suborbitope)[row][2] = perm[idx1];
            ++nintersections;
         }
         else if ( idx2 != perm[idx2] )
         {
            /* change order of idx1 and idx2 for left extensions */
            if ( leftextension )
            {
               (*suborbitope)[row][0] = idx2;
               (*suborbitope)[row][1] = idx1;
            }
            (*suborbitope)[row][2] = perm[idx2];
            ++nintersections;
         }
      }
   }
   else
   {
      /* check whether each cycle of perm intersects with a row of suborbitope */
      for (row = 0; row < nrows; ++row)
      {
         idx1 = (*suborbitope)[row][coltoextend];

         /* if idx1 is affected by perm, we can extend the row of the orbitope */
         if ( idx1 != perm[idx1] )
         {
            (*suborbitope)[row][nfilledcols] = perm[idx1];
            ++nintersections;
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


/** generate variable matrix for orbitope constraint handler */
static
SCIP_RETCODE generateOrbitopeVarsMatrix(
   SCIP_VAR****          vars,               /**< memory address to store pointer to orbitope variables */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
   SCIP_VAR**            permvars,           /**< superset of variables that are contained in orbitope */
   int                   npermvars,          /**< number of variables in permvars array */
   int**                 orbitopevaridx,     /**< permuted index table of variables in permvars that are contained in orbitope */
   int*                  columnorder         /**< permutation to reorder column of orbitopevaridx */
   )
{
   int i;
   int curcolumn;
   int nfilledcols = 0;

   assert( vars != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( permvars != NULL );
   assert( orbitopevaridx != NULL );
   assert( columnorder != NULL );

   curcolumn = ncols - 1;
   nfilledcols = 0;

   /* start filling vars matrix with the right-most column w.r.t. columnorder */
   while ( curcolumn >= 0 && columnorder[curcolumn] >= 0 )
   {
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][curcolumn] < npermvars );

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
   /* If we are in case 2), all columns should have been added to vars */
   if ( columnorder[curcolumn] == 0 )
      assert( nfilledcols == ncols );
   /* Otherwise, we are in case 1) or 3) and the only remaining non-negative
    *  column orders are 0 and 1. */
   else
   {
      /* add column with columnorder 1 to vars*/
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][1] < npermvars );

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][1]];
      }
      ++nfilledcols;

      /* add column with columnorder 0 to vars*/
      for (i = 0; i < nrows; ++i)
      {
         assert( orbitopevaridx[i][0] < npermvars );

         (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][0]];
      }
      ++nfilledcols;

      /* add columns with a negative column order to vars */
      if ( nfilledcols < ncols )
      {
         curcolumn = 2;
         while ( nfilledcols < ncols )
         {
            if ( columnorder[curcolumn] >= 0 )
               printf("Error\n");
            assert( columnorder[curcolumn] < 0 );

            for (i = 0; i < nrows; ++i)
            {
               assert( orbitopevaridx[i][curcolumn] < npermvars );

               (*vars)[i][nfilledcols] = permvars[orbitopevaridx[i][curcolumn]];
            }
            ++curcolumn;
            ++nfilledcols;
         }
      }
   }

   return SCIP_OKAY;
}


/** check whether components of the symmetry group can be handled completely
 *  by orbitopes */
static
SCIP_RETCODE detectOrbitopes(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< pointer to data of symbreak presolver */
   )
{
   int ncomponents;
   int* npermsincomponent;
   int** components;
   int** perms;
   int npermvars;
   SCIP_VAR** permvars;
   int i;

   assert( scip != NULL );
   assert( presoldata != NULL );

   /* terminate if no symmetry is present or compute components if not done yet */
   if ( presoldata->nperms == 0 )
      return SCIP_OKAY;
   else if ( presoldata->ncomponents == -1 )
   {
      SCIP_CALL( computeComponents(scip, presoldata) );
   }
   assert( presoldata->nperms > 0 );
   assert( presoldata->perms != NULL );
   assert( presoldata->npermvars > 0 );
   assert( presoldata->permvars != NULL );
   assert( presoldata->ncomponents > 0 );
   assert( presoldata->components != NULL );
   assert( presoldata->npermsincomponent != NULL );

   perms = presoldata->perms;
   npermvars = presoldata->npermvars;
   permvars = presoldata->permvars;
   ncomponents = presoldata->ncomponents;
   npermsincomponent = presoldata->npermsincomponent;
   components = presoldata->components;

   /* iterate over components */
   for (i = 0; i < ncomponents; ++i)
   {
      int ntwocyclescomp = -1;
      SCIP_Bool isorbitope = TRUE;
      SCIP_Bool* usedperm;
      int** orbitopevaridx;
      int* columnorder;
      int j;
      int row;
      int nfilledcols;
      int nusedperms;
      SCIP_VAR*** vars;
      int coltoextend;

      /* get properties of permutations */
      for (j = 0; j < npermsincomponent[i]; ++j)
      {
         SCIP_Bool iscompoftwocycles = FALSE;
         int ntwocyclesperm = 0;

         SCIP_CALL( getPermProperties(perms[components[i][j]], npermvars, &iscompoftwocycles, &ntwocyclesperm) );

         if ( ntwocyclescomp == -1 )
            ntwocyclescomp = ntwocyclesperm;

         /* no or different number of 2-cycles: permutations cannot generate orbitope */
         if ( ntwocyclescomp == 0 || ntwocyclescomp != ntwocyclesperm )
         {
            isorbitope = FALSE;
            break;
         }
      }

      if ( ! isorbitope )
         continue;

      /* iterate over permutations and check whether for each permutation there exists
       * another permutation whose 2-cycles intersect pairwise in exactly one element */

      /* whether a permutation was considered to contribute to orbitope */
      SCIP_CALL( SCIPallocBufferArray(scip, &usedperm, npermsincomponent[i]) );
      for (j = 0; j < npermsincomponent[i]; ++j)
         usedperm[j] = FALSE;
      nusedperms = 0;

      /* orbitope matrix for indices of variables in permvars array */
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx[j], npermsincomponent[i] + 1) );
      }

      /* order of columns of orbitopevaridx */
      SCIP_CALL( SCIPallocBufferArray(scip, &columnorder, npermsincomponent[i] + 1) );
      for (j = 0; j < npermsincomponent[i] + 1; ++j)
         columnorder[j] = npermsincomponent[i] + 2;

      /* fill first two columns of orbitopevaridx matrix */
      row = 0;
      for (j = 0; j < npermvars; ++j)
      {
         int permidx = components[i][0];

         /* avoid adding the same 2-cycle twice */
         if ( perms[permidx][j] > j )
         {
            orbitopevaridx[row][0] = j;
            orbitopevaridx[row++][1] = perms[permidx][j];
         }

         if ( row == ntwocyclescomp )
            break;
      }
      assert( row == ntwocyclescomp );

      usedperm[0] = TRUE;
      ++nusedperms;
      columnorder[0] = 0;
      columnorder[1] = 1;
      nfilledcols = 2;

      /* extend orbitopevaridx matrix to the left, i.e., find iteratively new permutations that
       * intersect the last added left column in each row in exactly one entry, starting with
       * column 0 */
      coltoextend = 0;
      for (j = 0; j < npermsincomponent[i]; ++j)
      {
         SCIP_Bool success = FALSE;
         SCIP_Bool infeasible = FALSE;

         if ( nusedperms == npermsincomponent[i] )
            break;

         if ( usedperm[j] )
            continue;

         SCIP_CALL( extendSubOrbitope(&orbitopevaridx, ntwocyclescomp, nfilledcols, coltoextend,
               perms[components[i][j]], TRUE, &success, &infeasible) );

         if ( infeasible )
         {
            isorbitope = FALSE;
            break;
         }
         else if ( success )
         {
            usedperm[j] = TRUE;
            ++nusedperms;
            coltoextend = nfilledcols;
            columnorder[nfilledcols++] = -1; /* mark column to be filled from the left */
            j = 0; /* reset j since previous permutations can now intersect with the latest added column */
         }
      }

      if ( ! isorbitope )
      {
         /* free data structures */
         SCIPfreeBufferArray(scip, &columnorder);
         for (j = 0; j < ntwocyclescomp; ++j)
            SCIPfreeBufferArray(scip, &orbitopevaridx[j]);
         SCIPfreeBufferArray(scip, &orbitopevaridx);
         SCIPfreeBufferArray(scip, &usedperm);

         continue;
      }

      coltoextend = 1;
      for (j = 0; j < npermsincomponent[i]; ++j)
      {
         SCIP_Bool success = FALSE;
         SCIP_Bool infeasible = FALSE;

         if ( nusedperms == npermsincomponent[i] )
            break;

         if ( usedperm[j] )
            continue;

         SCIP_CALL( extendSubOrbitope(&orbitopevaridx, ntwocyclescomp, nfilledcols, coltoextend,
               perms[components[i][j]], FALSE, &success, &infeasible) );

         if ( infeasible )
         {
            isorbitope = FALSE;
            break;
         }
         else if ( success )
         {
            usedperm[j] = TRUE;
            ++nusedperms;
            coltoextend = nfilledcols;
            columnorder[nfilledcols] = 1; /* mark column to be filled from the right */
            ++nfilledcols;
            j = 0; /* reset j since previous permutations can now intersect with the latest added column */
         }
      }

      if ( nusedperms < npermsincomponent[i] )
         isorbitope = FALSE;

      if ( ! isorbitope )
      {
         /* free data structures */
         SCIPfreeBufferArray(scip, &columnorder);
         for (j = 0; j < ntwocyclescomp; ++j)
            SCIPfreeBufferArray(scip, &orbitopevaridx[j]);
         SCIPfreeBufferArray(scip, &orbitopevaridx);
         SCIPfreeBufferArray(scip, &usedperm);

         continue;
      }

      /* we have found an orbitope, prepare data for orbitope conshdlr */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &vars[j], npermsincomponent[i]) );
      }

      /* prepare variable matrix (reorder columns of orbitopevaridx) */
      SCIP_CALL( generateOrbitopeVarsMatrix(&vars, ntwocyclescomp, npermsincomponent[i] + 1, permvars, npermvars,
            orbitopevaridx, columnorder) );

      /* @todo: add orbitope if full orbitope conshdlr is implemented */
      SCIPinfoMessage(scip, NULL, "Component %d is an orbitope with %d rows and %d columns.\n", i, ntwocyclescomp, npermsincomponent[i] + 1);

      /* free data structures */
      for (j = 0; j < ntwocyclescomp; ++j)
         SCIPfreeBufferArray(scip, &vars[j]);
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &columnorder);
      for (j = 0; j < ntwocyclescomp; ++j)
         SCIPfreeBufferArray(scip, &orbitopevaridx[j]);
      SCIPfreeBufferArray(scip, &orbitopevaridx);
      SCIPfreeBufferArray(scip, &usedperm);
   }

   return SCIP_OKAY;
}


/** add symresack constraints */
static
SCIP_RETCODE addSymresackConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOL*          presol              /**< symmetry breaking presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;
   int nperms;
   int** perms;
   int nsymresackcons = 0;
   SCIP_VAR** permvars;
   int npermvars;
   SCIP_Bool conssaddlp;
   int p;

   assert( scip != 0 );
   assert( presol != 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   perms = presoldata->perms;
   nperms = presoldata->nperms;
   permvars = presoldata->permvars;
   npermvars = presoldata->npermvars;
   conssaddlp = presoldata->conssaddlp;

   assert( nperms <= 0 || (nperms > 0 && perms != NULL) );
   assert( permvars != NULL );
   assert( npermvars > 0 );

   /* loop through perms and add symresack constraints */
   for (p = 0; p < nperms; ++p)
   {
      SCIP_CONS* cons;
      SCIP_Bool success = FALSE;

      SCIP_CALL( SCIPcreateConsSymresack(scip, &cons, "symresack", perms[p], permvars, npermvars,
            conssaddlp, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

      /* add the constraint only if the constraint is not trivial */
      if ( success )
      {
         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         presoldata->genconss[presoldata->ngenconss] = cons;
         ++presoldata->ngenconss;
         ++presoldata->nsymresacks;
         ++nsymresackcons;
         SCIPdebugMsg(scip, "Added symresack constraint: %u.\n", nsymresackcons);
      }
      else
      {
         /* otherwise the constraint was not generated */
         assert( cons == NULL );
      }
   }

   return SCIP_OKAY;
}


/** analyze generators and add symmetry breaking constraints */
static
SCIP_RETCODE addSymmetryBreakingConstraints(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_Bool decision;

   assert( scip != 0 );
   assert( presol != 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   /* exit if no or only trivial symmetry group is available */
   if ( presoldata->nperms < 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetBoolParam(scip, "presolving/" PRESOL_NAME"/addsymresacks", &decision) );

   if ( decision )
   {
      SCIP_CALL( addSymresackConss(scip, presol) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   SCIPfreeMemory(scip, &presoldata);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitSymbreak)
{
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   SCIPdebugMsg(scip, "Initpre method of symmetry breaking presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   /* check whether we should run */
#if 0
   SCIP_CALL( SCIPgetBoolParam(scip, "usesymbreak", &presoldata->enabled) );
#else
   presoldata->enabled = TRUE;
#endif

   if ( presoldata->enabled )
   {
      /* allow all problem specifications, since we handle them in the code above */
      SYMsetSpecRequirement(presoldata->symmetrypresol, SYM_SPEC_BINARY);
      SYMsetSpecRequirement(presoldata->symmetrypresol, SYM_SPEC_INTEGER);
      SYMsetSpecRequirement(presoldata->symmetrypresol, SYM_SPEC_REAL);
   }

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitSymbreak)
{
   SCIP_PRESOLDATA* presoldata;
   int i;
   SCIP_CONS** genconss;
   int ngenconss;

   assert( scip != NULL );
   assert( presol != NULL );

   SCIPdebugMessage("Exit method of symmetry breaking presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   ngenconss = presoldata->ngenconss;

   /* release constraints */
   genconss = presoldata->genconss;
   for (i = 0; i < ngenconss; ++i)
   {
      SCIP_CONS* cons = genconss[i];
      assert( cons != NULL );

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* reset data (necessary if another problem is read in the SCIP shell) */

   /* free pointers to symmetry group and binary variables */
   if ( presoldata->ngenconss > 0 )
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->genconss, presoldata->ngenconss);

   presoldata->genconss = NULL;
   presoldata->ngenconss = 0;

   /* free orbit structures */
   for (i = 0; i < presoldata->norbits; ++i)
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->orbits[i], presoldata->nvarsinorbits[i]);
   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->orbits, presoldata->norbits);
   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->nvarsinorbits, presoldata->norbits);

   presoldata->orbits = NULL;
   presoldata->nvarsinorbits = NULL;
   presoldata->norbits = -1;

   /* free components */
   for (i = 0; i < presoldata->ncomponents; ++i)
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->components[i], presoldata->npermsincomponent[i]);
    SCIPfreeBlockMemoryArrayNull(scip, &presoldata->components, presoldata->ncomponents);
    SCIPfreeBlockMemoryArrayNull(scip, &presoldata->npermsincomponent, presoldata->ncomponents);

    presoldata->components = NULL;
    presoldata->npermsincomponent = NULL;
    presoldata->ncomponents = -1;

    /* reset basic data */
    presoldata->addedconss = FALSE;
    presoldata->computedsymmetry = FALSE;
    presoldata->enabled = TRUE;
    presoldata->early = FALSE;
    presoldata->nperms = -1;
    presoldata->nsymresacks = 0;

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* check whether we have to run the presolver at the beginning of presolving */
   presoldata->early = ! SYMcomputeSymmetryPresolved(presoldata->symmetrypresol);

   if ( presoldata->early && presoldata->enabled )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Executing presolver <%s> early, since symmetries are computed early.\n\n", SCIPpresolGetName(presol));

      SCIP_CALL( SCIPsetIntParam(scip, "presolving/" PRESOL_NAME"/priority", 90000000) );
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int i;
   int noldfixedvars = *nfixedvars;
   int noldaggrvars = *naggrvars;
   int noldbdchgs = *nchgbds;
   int noldaddconss = *naddconss;


   assert( scip != NULL );
   assert( presol != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   /* possibly skip presolver */
   if ( ! presoldata->enabled )
      return SCIP_OKAY;

   /* skip presolving if we are not at the end */
   if ( ! presoldata->early && ! SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Presolving method of symmetry breaking presolver ...\n");

   /* get symmetry information, if not already computed */
   if ( ! presoldata->computedsymmetry )
   {
      SCIPdebugMsg(scip, "Symmetry breaking presolver: computing symmetry ...\n");
      assert( presoldata->nperms < 0 );

      /* get symmetries */
      SCIP_CALL( SCIPgetSymmetryGenerators(scip, presoldata->symmetrypresol, &(presoldata->npermvars),
            &(presoldata->permvars), &(presoldata->nperms), &(presoldata->perms)) );

      presoldata->computedsymmetry = TRUE;

      if ( presoldata->nperms <= 0 )
      {
         SCIPdebugMessage("Symmetry breaking presolver: no symmetry has been found, turning presolver off.\n");
         presoldata->enabled = FALSE;
         return SCIP_OKAY;
      }
      else if ( presoldata->nperms > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->genconss), presoldata->nperms) );

         SCIP_CALL( computeGroupOrbits(scip, presoldata) );

         SCIP_CALL( computeComponents(scip, presoldata) );

         SCIP_CALL( detectOrbitopes(scip, presoldata) );

      }
   }

   /* at this point, the symmetry group should be computed and nontrivial */
   assert( presoldata->nperms > 0 );
   *result = SCIP_DIDNOTFIND;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* if not already done, add symmetry breaking constraints */
   if ( ! presoldata->addedconss )
   {
      /* add symmetry breaking constraints */
      int noldngenconns = presoldata->ngenconss;

      SCIP_CALL( addSymmetryBreakingConstraints(scip, presol) );

      presoldata->addedconss = TRUE;
      *naddconss += presoldata->ngenconss - noldngenconns;
      SCIPdebugMsg(scip, "Added symmetry breaking constraints: %u.\n", presoldata->ngenconss - noldngenconns);

      /* if constraints have been added, loop through generated constraints and presolve each */
      for (i = 0; i < presoldata->ngenconss; ++i)
      {
         SCIP_CALL( SCIPpresolCons(scip, presoldata->genconss[i], nrounds, SCIP_PRESOLTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
               nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
               nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff has been detected */
         if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMessage("Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(presoldata->genconss[i]));
            return SCIP_OKAY;
         }
      }
      SCIPdebugMsg(scip, "Presolved %u constraints generated by symbreak.\n", presoldata->ngenconss);

      /* possibly stop */
      if ( SCIPisStopped(scip) )
         return SCIP_OKAY;
   }

   /* determine success */
   if ( noldaddconss + noldbdchgs + noldaggrvars + noldfixedvars > 0 )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the symbreak presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolSymbreak(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata = 0;
   SCIP_PRESOL* presol = 0;

   /* create presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* we must call the constructor explictly, because memory was m'alloced and not new'ed */
   presoldata->addedconss = FALSE;
   presoldata->computedsymmetry = FALSE;
   presoldata->enabled = TRUE;
   presoldata->early = FALSE;
   presoldata->nsymresacks = 0;
   presoldata->ngenconss = 0;
   presoldata->genconss = NULL;
   presoldata->nperms = -1;
   presoldata->norbits = -1;
   presoldata->nvarsinorbits = NULL;
   presoldata->orbits = NULL;
   presoldata->ncomponents = -1;
   presoldata->npermsincomponent = NULL;
   presoldata->components = NULL;

   /* determine cons_symmetries constraint handler (preuse presol) */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == 0 )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   presoldata->symmetrypresol = presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecSymbreak, presoldata) );

   /* set additional callbacks */
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSymbreak) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitSymbreak) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitSymbreak) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreSymbreak) );

   /* add symmetry breaking presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME"/conssaddlp",
         "Should the symmetry breaking constraints be added to the LP?",
         &presoldata->conssaddlp, TRUE, DEFAULT_CONSSADDLP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME"/addsymresacks",
         "Add inequalities for symresacks for each generator?",
         &presoldata->addsymresacks, TRUE, DEFAULT_ADDSYMRESACKS, NULL, NULL) );

   return SCIP_OKAY;
}
