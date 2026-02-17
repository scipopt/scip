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

/**@file   sym_sst.c
 * @ingroup DEFPLUGINS_SYM
 * @brief  symmetry handler for SST cuts
 * @author Christopher Hojny
 * @author Jasper van Doornmalen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_linear.h"
#include "scip/pub_cons.h"
#include "scip/pub_implics.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sym.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sym.h"
#include "scip/scip_var.h"
#include "scip/sym_sst.h"
#include "scip/symmetry.h"
#include "scip/type_implics.h"

/* symmetry handler properties */
#define SYM_NAME            "sym_sst"
#define SYM_DESC            "symmetry handler for SST cuts"
#define SYM_PRIORITY           -200000       /**< priority of try-add function */
#define SYM_PRESOLPRIORITY    -1000000       /**< priority of presolving method */
#define SYM_MAXPRESOLROUNDS          1       /**< maximum number of presolving rounds */

/* default value of parameters */
#define DEFAULT_ORBITRULE            1       /**< index of tie break rule for selecting orbit for SST constraints? */
#define DEFAULT_LEADERRULE           0       /**< index of rule for selecting leader variables for SST constraints? */
#define DEFAULT_LEADERVARTYPE        6       /**< bitset encoding allowed vars types for leaders (1 bin; 2 int; 4 cont);
                                              *   if multiple types are allowed, take the one with most affected vars */
#define DEFAULT_COMPUTENEWPERMS   TRUE       /**< Shall additional permutations of symmetry component be computed? */
#define DEFAULT_MAXNNEWPERMS       100       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer);
                                              *   -1: unbounded */
#define DEFAULT_ADDCONFLICTCUTS       TRUE   /**< Should SST constraints be added if we use a conflict based rule? */
#define DEFAULT_MIXEDCOMPONENTS       TRUE   /**< Should SST constraints be added if a symmetry component contains
                                              *   variables of different types? */

#define ISSSTBINACTIVE(x)          (((unsigned) x & SCIP_SSTTYPE_BINARY) != 0)
#define ISSSTINTACTIVE(x)          (((unsigned) x & SCIP_SSTTYPE_INTEGER) != 0)
#define ISSSTCONTACTIVE(x)         (((unsigned) x & SCIP_SSTTYPE_CONTINUOUS) != 0)

/** selection rules for leaders in SST cuts */
enum SST_LeaderRule
{
   SST_LEADERRULE_FIRSTINORBIT        = 0,       /**< first var in orbit */
   SST_LEADERRULE_LASTINORBIT         = 1,       /**< last var in orbit */
   SST_LEADERRULE_MAXCONFSINORBIT     = 2        /**< var with most conflicting vars in its orbit */
};
typedef enum SST_LeaderRule SST_LEADERRULE;

/** tie breaks for leader rule based on the leader's orbit */
enum SST_OrbitRule
{
   SST_ORBITRULE_MINORBIT            = 0,    /**< orbit of minimum size */
   SST_ORBITRULE_MAXORBIT            = 1,    /**< orbit of maximum size */
   SST_ORBITRULE_MAXCONFSINORBIT     = 2     /**< orbit with maximum number of vars in conflict with leader */
};
typedef enum SST_OrbitRule SST_ORBITRULE;

/** variable types for leader in Schreier-Sims cuts */
enum SST_Vartype
{
   SST_SSTTYPE_BINARY                 = 1,    /**< binary variables */
   SST_SSTTYPE_INTEGER                = 2,    /**< integer variables */
   SST_SSTTYPE_CONTINUOUS             = 4     /**< continuous variables */
};
typedef enum SST_Vartype SST_VARTYPE;

/** conflict data structure for SST cuts */
struct SST_ConflictData
{
   SCIP_VAR*             var;                /**< variable belonging to node */
   int                   orbitidx;           /**< orbit of variable w.r.t. current stabilizer subgroup
                                              *   or -1 if not affected by symmetry */
   int                   nconflictinorbit;   /**< number of variables the node's var is in conflict with */
   int                   orbitsize;          /**< size of the variable's orbit */
   int                   posinorbit;         /**< position of variable in its orbit */
   SCIP_Bool             active;             /**< whether variable has not been fixed by Schreier-Sims code */
   SCIP_CLIQUE**         cliques;            /**< list of setppc constraints */
   int                   ncliques;           /**< number of setppc constraints */
};
typedef struct SST_ConflictData SST_CONFLICTDATA;

/** compare function for sorting an array by the addresses of its members */
static
SCIP_DECL_SORTPTRCOMP(sortByPointerValue)
{
   /* @todo move to misc.c? Or even better: use a different sorting function, e.g., according to clique id */
   if ( elem1 < elem2 )
      return -1;
   else if ( elem1 > elem2 )
      return +1;
   return 0;
}

/*
 * Data structures
 */

/** symmetry component data */
struct SCIP_SymCompData
{
   SCIP_CONS**           sstconss;           /**< list of generated SST constraints */
   int                   nsstconss;          /**< number of generated SST constraints */
   int                   maxnsstconss;       /**< maximum number of conss in sstconss */
   SCIP_VAR**            orbitvars;          /**< variables used in SST cuts, orbits are stored consecutively,
                                              *   first element per orbit is the leader */
   int*                  orbitbegins;        /**< array indicating begin positions of orbits */
   int                   norbits;            /**< number of orbits */
   int                   lenorbitvars;       /**< length of orbitvars array */
   int                   lenorbitbegins;     /**< length of orbitbegins array */
};

/** symmetry handler data */
struct SCIP_SymhdlrData
{
   int                   leaderrule;         /**< rule to select leader */
   int                   orbitrule;          /**< rule to select orbit */
   int                   leadervartype;      /**< bitset encoding which variable types can be leaders;
                                              *   if multiple types are allowed, take the one with most affected vars */
   SCIP_Bool             computenewperms;    /**< Shall additional permutations of symmetry component be computed? */
   int                   maxnnewperms;       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer);
                                              *   -1: unbounded */
   SCIP_Bool             addconflictcuts;    /**< Should SST constraints be added if we use a conflict based rule? */
   SCIP_Bool             mixedcomponents;    /**< Should SST constraints be added if a symmetry component contains
                                              *   variables of different types? */
};

/*
 * Local methods
 */

/** creates transposed permutation matrix */
static
SCIP_RETCODE createPermsTranspose(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SYMTYPE           symtype,            /**< type of symmetries considered */
   int**                 perms,              /**< (signed) permutations for which SST constraints shall be added */
   int                   nperms,             /**< number of signed permutations */
   int                   npermvars,          /**< number of variables */
   int**                 newperms,           /**< additional permutations (or NULL) */
   int                   nnewperms,          /**< number of additional permutations */
   int***                permstrans          /**< pointer to hold transposed permutation matrix */
   )
{
   int p;
   int v;
   int i;
   int len;
   int ntotalperms;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(npermvars > 0);
   assert(newperms != NULL || nnewperms == 0);
   assert(permstrans != NULL);

   len = symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;
   ntotalperms = nperms + nnewperms;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, permstrans, len) );
   for( v = 0; v < len; ++v )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*permstrans)[v], ntotalperms) );
      for( p = 0; p < nperms; ++p )
         (*permstrans)[v][p] = perms[p][v];
      for( i = 0; i < nnewperms; ++i, ++p )
         (*permstrans)[v][p] = newperms[i][v];
   }

   return SCIP_OKAY;
}

/** computes number of variables affected by symmetries by type */
static
SCIP_RETCODE computeAffectedVartypesCount(
   int**                 perms,              /**< (signed) permutations */
   int                   nperms,             /**< number of signed permutations */
   SCIP_VAR**            permvars,           /**< variables on which the (signed) permutations act */
   int                   npermvars,          /**< number of variables */
   int*                  nmovedbinvars,      /**< pointer to store number of affected binary variables */
   int*                  nmovedintvars,      /**< pointer to store number of affected integer variables */
   int*                  nmovedcontvars      /**< pointer to store number of affected continuous variables */
   )
{
   int v;
   int p;

   assert(nperms >= 0);

   *nmovedbinvars = 0;
   *nmovedintvars = 0;
   *nmovedcontvars = 0;

   for( v = 0; v < npermvars; ++v )
   {
      for( p = 0; p < nperms; ++p )
      {
         if( perms[p][v] != v )
         {
            switch( SCIPgetSymInferredVarType(permvars[v]) )
            {
            case SCIP_VARTYPE_BINARY:
               ++(*nmovedbinvars);
               break;
            case SCIP_VARTYPE_INTEGER:
               ++(*nmovedintvars);
               break;
            case SCIP_VARTYPE_CONTINUOUS:
               ++(*nmovedcontvars);
               break;
            default:
               SCIPerrorMessage("unknown variable type\n");
               return SCIP_INVALIDDATA;
            } /*lint !e788*/
            break;
         }
      }
   }

   return SCIP_OKAY;
}

/** create conflict graph either for symmetric or for all variables
 *
 *  This routine just creates the graph, but does not add (symmetry) information to its nodes.
 *  This has to be done separately by the routine updateSymInfoConflictGraphSST().
 *
 *  The function returns with varconflicts as NULL when we do not create it.
 */
static
SCIP_RETCODE createConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SST_CONFLICTDATA**    varconflicts,       /**< pointer to store the variable conflict data */
   SCIP_VAR**            conflictvars,       /**< array of variables to encode in conflict graph */
   int                   nconflictvars,      /**< number of vars to encode in conflict graph */
   SCIP_HASHMAP*         conflictvarmap      /**< map of variables to indices in conflictvars array */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_VAR** cliquevars;
   SCIP_CLIQUE* clique;
   int* tmpncliques;
   int ncliques;
   int ncliquevars;
   int node;
   int c;
   int i;

#ifdef SCIP_DEBUG
   int varncliques = 0;
#endif

   assert(scip != NULL);
   assert(varconflicts != NULL);
   assert(conflictvars != NULL);
   assert(nconflictvars > 0);

   /* we set the pointer of varconflicts to NULL to illustrate that we didn't generate it */
   *varconflicts = NULL;

   /* get cliques for creating conflict structure */
   cliques = SCIPgetCliques(scip);
   ncliques = SCIPgetNCliques(scip);
   if( ncliques == 0 )
   {
      SCIPdebugMsg(scip, "No cliques present --> construction of conflict structure aborted.\n");
      return SCIP_OKAY;
   }

   /* construct variable conflicts */
   SCIPdebugMsg(scip, "Construction of conflict structure:\n");
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, varconflicts, nconflictvars) );
   for( i = 0; i < nconflictvars; ++i )
   {
      (*varconflicts)[i].ncliques = 0;
      (*varconflicts)[i].active = TRUE;
      (*varconflicts)[i].var = conflictvars[i];
      /* set remaining variable conflictdata at neutral entries */
      (*varconflicts)[i].cliques = NULL;
      (*varconflicts)[i].orbitidx = -1;
      (*varconflicts)[i].nconflictinorbit = 0;
      (*varconflicts)[i].orbitsize = -1;
      (*varconflicts)[i].posinorbit = -1;
   }

   /* Store, for each variable, the conflict cliques it is contained in.
    * In three steps:
    * (1.) Count the number of cliques it's contained in, per var, then
    * (2.) Create the array of this size, and
    * (3.) Fill the array with the cliques.
    * Starting with (1.):
    */
   for( c = 0; c < ncliques; ++c )
   {
      clique = cliques[c];
      assert(clique != NULL);

      cliquevars = SCIPcliqueGetVars(clique);
      ncliquevars = SCIPcliqueGetNVars(clique);
      assert(cliquevars != NULL);
      assert(ncliquevars > 0);

      SCIPdebugMsg(scip, "\tIdentify edges for clique ID: %u; Index: %d).\n", SCIPcliqueGetId(clique),
         SCIPcliqueGetIndex(clique));

      /* for all variables, list which cliques it is part of */
      for( i = 0; i < ncliquevars; ++i )
      {
         node = SCIPhashmapGetImageInt(conflictvarmap, cliquevars[i]);

         /* skip variables not in the conflictvars array (so not in hashmap, too) */
         if( node == INT_MAX )
            continue;
         assert(node >= 0);
         assert(node < nconflictvars);

         assert((*varconflicts)[node].var == cliquevars[i]);
         (*varconflicts)[node].active = TRUE;
         (*varconflicts)[node].ncliques++;
      }
   }

   /* (2.) allocate the arrays */
   for( i = 0; i < nconflictvars; ++i )
   {
      assert((*varconflicts)[i].ncliques >= 0);
      assert((*varconflicts)[i].cliques == NULL);
      if( (*varconflicts)[i].ncliques > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*varconflicts)[i].cliques, (*varconflicts)[i].ncliques) );
      }
   }

   /* (3.) fill the clique constraints */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &tmpncliques, nconflictvars) );
   for( c = 0; c < ncliques; ++c )
   {
      clique = cliques[c];
      assert(clique != NULL);

      cliquevars = SCIPcliqueGetVars(clique);
      ncliquevars = SCIPcliqueGetNVars(clique);
      assert(cliquevars != NULL);
      assert(ncliquevars > 0);

      SCIPdebugMsg(scip, "\tAdd edges for clique ID: %u; Index: %d).\n", SCIPcliqueGetId(clique),
         SCIPcliqueGetIndex(clique));

      /* for all variables, list which cliques it is part of */
      for( i = 0; i < ncliquevars; ++i )
      {
         node = SCIPhashmapGetImageInt(conflictvarmap, cliquevars[i]);

         /* skip variables not in the conflictvars array (so not in hashmap, too) */
         if( node == INT_MAX )
            continue;

         assert(node >= 0);
         assert(node < nconflictvars);
         assert((*varconflicts)[node].var == cliquevars[i]);

         /* add clique to the cliques */
         assert(tmpncliques[node] < (*varconflicts)[node].ncliques);
         assert((*varconflicts)[node].cliques != NULL);
         (*varconflicts)[node].cliques[tmpncliques[node]++] = clique;

#ifdef SCIP_DEBUG
         varncliques++;
#endif
      }
   }

   /* sort the variable cliques by the address, so checkSortedArraysHaveOverlappingEntry can detect intersections */
   for( i = 0; i < nconflictvars; ++i )
   {
      SCIPsortPtr((void**)(*varconflicts)[i].cliques, sortByPointerValue, (*varconflicts)[i].ncliques);
   }

#ifndef NDEBUG
   for( i = 0; i < nconflictvars; ++i )
   {
      assert(tmpncliques[i] == (*varconflicts)[i].ncliques);
   }
#endif

   SCIPfreeBufferArray(scip, &tmpncliques);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Construction of conflict graph terminated; %d variable-clique combinations detected.\n",
      varncliques);
#endif

   return SCIP_OKAY;
}

/** checks whether two arrays that are sorted with the same comparator have a common element */
static
SCIP_Bool checkSortedArraysHaveOverlappingEntry(
   void**                arr1,               /**< first array */
   int                   narr1,              /**< number of elements in first array */
   void**                arr2,               /**< second array */
   int                   narr2,              /**< number of elements in second array */
   SCIP_DECL_SORTPTRCOMP((*compfunc))        /**< comparator function that was used to sort arr1 and arr2; must define a total ordering */
   )
{
   /* @todo move to misc.c? */
   int it1;
   int it2;
   int cmp;

   assert(arr1 != NULL || narr1 == 0);
   assert(narr1 >= 0);
   assert(arr2 != NULL || narr2 == 0);
   assert(narr2 >= 0);
   assert(compfunc != NULL);

   /* there is no overlap if one of the two arrays is empty */
   if( narr1 <= 0 )
      return FALSE;
   if( narr2 <= 0 )
      return FALSE;

   it1 = 0;
   it2 = 0;

   while( TRUE )  /*lint !e716*/
   {
      cmp = compfunc(arr1[it1], arr2[it2]);
      if( cmp < 0 )
      {
         /* comparison function determines arr1[it1] < arr2[it2]
          * increase iterator for arr1
          */
         if( ++it1 >= narr1 )
            break;
         continue;
      }
      else if( cmp > 0 )
      {
         /* comparison function determines arr1[it1] > arr2[it2]
          * increase iterator for arr2
          */
         if( ++it2 >= narr2 )
            break;
         continue;
      }
      else
      {
         /* the entries arr1[it1] and arr2[it2] are the same with respect to the comparison function */
         assert(cmp == 0);
         return TRUE;
      }
   }

   /* no overlap detected */
   assert(it1 >= narr1 || it2 >= narr2);

   return FALSE;
}

/** add Schreier-Sims constraints for a specific orbit and update Schreier-Sims table */
static
SCIP_RETCODE addSSTConssOrbitAndUpdateSST(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< type of symmetries */
   int                   id,                 /**< identifier of symmetry component for that SST conss are added */
   SST_CONFLICTDATA*     varconflicts,       /**< conflict graph or NULL if useconflictgraph == FALSE */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            permvars,           /**< permvars array */
   int                   npermvars,          /**< number of variables in permvars */
   SCIP_Real*            permvardomaincenter,/**< domain center of permvars */
   int                   leaderrule,         /**< rule to select leader of SST conss */
   int                   orbitrule,          /**< rule to select orbit of SST conss */
   SCIP_Bool             addconflictcuts,    /**< Shall SST conss be added when a conflict-based rule is used? */
   SCIP_CONS***          sstconss,           /**< pointer to store SST conss */
   int*                  nsstconss,          /**< pointer to store number of SST conss */
   int*                  maxnsstconss,       /**< pointer to store maximum number of conss sstconss can hold */
   int*                  orbits,             /**< symmetry orbits */
   int*                  orbitbegins,        /**< array storing begin position for each orbit */
   int                   orbitidx,           /**< index of orbit for Schreier-Sims constraints */
   int                   orbitleaderidx,     /**< index of leader variable for Schreier-Sims constraints */
   SCIP_Shortbool*       orbitvarinconflict, /**< indicator whether orbitvar is in conflict with orbit leader */
   int                   norbitvarinconflict,/**< number of variables in conflict with orbit leader */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{ /*lint --e{613,641}*/
   SCIP_CONS* cons;
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_Real rhs = 0.0;
   int orbitsize;
   int posleader;
   int poscur;
   int ncuts = 0;
   SCIP_Bool addcuts = TRUE;
   int i;

   /* @todo Replace SST_CONFLICTDATA* by SST_CONFLICTDATA** */

   assert(scip != NULL);
   assert(permvars != NULL);
   assert(orbits != NULL);
   assert(orbitbegins != NULL);
   assert(orbitidx >= 0);
   assert(orbitleaderidx >= 0);
   assert(orbitvarinconflict != NULL || varconflicts == NULL);
   assert(norbitvarinconflict >= 0);
   assert(nchgbds != NULL);

   orbitsize = orbitbegins[orbitidx + 1] - orbitbegins[orbitidx];

   /* variables in conflict with leader are fixed and not treated by a cut; trailing -1 to not count the leader */
   if( leaderrule == SST_LEADERRULE_MAXCONFSINORBIT || orbitrule == SST_ORBITRULE_MAXCONFSINORBIT )
      addcuts = addconflictcuts;

   if( addcuts )
      ncuts = orbitsize - norbitvarinconflict - 1;

   /* (re-)allocate memory for SST constraints */
   if( ncuts > 0 )
   {
      if( *nsstconss == 0 )
      {
         assert(sstconss != NULL);
         assert(*maxnsstconss == 0);
         *maxnsstconss = 2 * ncuts;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, sstconss, *maxnsstconss) );
      }
      else if( *nsstconss + ncuts > *maxnsstconss )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, *maxnsstconss + 2 * ncuts);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, sstconss, *maxnsstconss, newsize) );
         *maxnsstconss = newsize;
      }
   }

   /* add SST constraints vars[0] >= vars[1], where vars[0] is always the leader */
   posleader = orbitbegins[orbitidx] + orbitleaderidx;
   vars[0] = permvars[orbits[posleader]];
   vals[0] = -1.0;
   vals[1] = 1.0;
   *nchgbds = 0;
   for( i = 0, poscur = orbitbegins[orbitidx]; i < orbitsize; ++i, ++poscur )
   {
      SCIP_CALL addcut;

      if( i == orbitleaderidx )
      {
         assert(orbitvarinconflict == NULL || ! orbitvarinconflict[i]);
         continue;
      }
      addcut = addcuts;

      vars[1] = permvars[orbits[poscur]];

      /* for signed permutations, we need to adapt the rhs of SST cuts (artificially center variables at 0) */
      if( symtype == SYM_SYMTYPE_SIGNPERM )
         rhs = permvardomaincenter[orbits[poscur]] - permvardomaincenter[orbits[posleader]];

      /* if the i-th variable in the orbit is in a conflict with the leader, fix it to 0 */
      if( varconflicts != NULL )
      {
         if( orbitvarinconflict[i] )
         {
            assert(SCIPvarIsBinary(vars[1]));
            assert(SCIPvarGetLbLocal(vars[1]) < 0.5);
            assert(varconflicts != NULL);

            /* if variable is fixed */
            if( SCIPvarGetUbLocal(vars[1]) > 0.5 )
            {
               SCIP_CALL( SCIPchgVarUb(scip, vars[1], 0.0) );
               ++(*nchgbds);

               /* deactivate the fixed variable (cannot contribute to a conflict anymore) */
               assert(varconflicts[orbits[poscur]].active);
               varconflicts[orbits[poscur]].active = FALSE;
               addcut = FALSE;
            }

            /* reset value */
            orbitvarinconflict[i] = FALSE;
         }
      }

      if( addcut )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SSTcut_symcomp_%d_%d_%d", id, orbits[posleader], orbits[poscur]);
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, -SCIPinfinity(scip), rhs,
               FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         (*sstconss)[(*nsstconss)++] = cons;
      }
   }

   return SCIP_OKAY;
}

/** update symmetry information of conflict graph */
static
SCIP_RETCODE updateSymInfoConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SST_CONFLICTDATA*     varconflicts,       /**< conflict structure */
   SCIP_VAR**            conflictvars,       /**< variables encoded in conflict structure */
   int                   nconflictvars,      /**< number of nodes/vars in conflict structure */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int                   norbits             /**< number of non-trivial orbits */
   )
{
   int i;
   int j;
   int ii;
   int jj;
   int r; /* r from orbit, the orbit index. */

   assert(scip != NULL);
   assert(varconflicts != NULL);
   assert(conflictvars != NULL);
   assert(nconflictvars > 0);
   assert(orbits != NULL);
   assert(orbitbegins != NULL);
   assert(norbits >= 0);

   /* initialize/reset variable information of nodes in conflict graph */
   for( i = 0; i < nconflictvars; ++i )
   {
      /* (re-)set node data */
      varconflicts[i].orbitidx = -1;
      varconflicts[i].nconflictinorbit = 0;
      varconflicts[i].orbitsize = -1;
      varconflicts[i].posinorbit = -1;
   }

   /* add orbit information to nodes of conflict graph */
   for( r = 0; r < norbits; ++r )
   {
      int posinorbit = 0;
      int orbitsize;

      orbitsize = orbitbegins[r + 1] - orbitbegins[r];
      assert(orbitsize >= 0);

      for( i = orbitbegins[r]; i < orbitbegins[r + 1]; ++i )
      {
         int pos;

         /* get variable and position in conflict graph */
         pos = orbits[i];
         assert(pos < nconflictvars);
         assert(varconflicts[pos].var == conflictvars[pos]);

         varconflicts[pos].orbitidx = r;
         varconflicts[pos].nconflictinorbit = 0;
         varconflicts[pos].orbitsize = orbitsize;
         varconflicts[pos].posinorbit = posinorbit++;
      }

      /* determine nconflictsinorbit
       *
       * For each pair of active variables in this orbit, check if it is part of a conflict clique.
       * Use that we store the cliques of this type in varconflicts[pos].cliques.
       * These lists are sorted (by the address of the constraint), so we only need to check for each i, j in the orbit
       * whether they are contained in the same clique.
       */
      for( i = orbitbegins[r]; i < orbitbegins[r + 1]; ++i )
      {
         ii = orbits[i];
         assert(varconflicts[ii].orbitidx == r);

         /* skip inactive variables */
         if( ! varconflicts[ii].active )
            continue;

         for( j = i + 1; j < orbitbegins[r + 1]; ++j )
         {
            jj = orbits[j];
            assert(varconflicts[jj].orbitidx == r);

            /* skip inactive variables */
            if( !varconflicts[jj].active )
               continue;

            /* Check if i and j are overlapping in some clique, where only one of the two could have value 1.
             * Use that cliques are sorted by the constraint address.
             *
             * @todo A better sorted order would be: First constraints with large variables (higher hitting probability)
             *  and then by a unique constraint identifier (address, or conspos).
             */
            if( checkSortedArraysHaveOverlappingEntry((void**)varconflicts[ii].cliques,
               varconflicts[ii].ncliques, (void**)varconflicts[jj].cliques, varconflicts[jj].ncliques,
               sortByPointerValue) )
            {
               /* there is overlap! */
               ++varconflicts[ii].nconflictinorbit;
               ++varconflicts[jj].nconflictinorbit;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** frees conflict graph */
static
SCIP_RETCODE freeConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SST_CONFLICTDATA**    varconflicts,       /**< conflict graph */
   int                   nvars               /**< number of nodes in conflict graph */
   )
{
   int i;
   int n;

   assert(scip != NULL);
   assert(varconflicts != NULL);
   assert(*varconflicts != NULL);
   assert(nvars >= 0);

   for( i = nvars - 1; i >= 0; --i )
   {
      n = (*varconflicts)[i].ncliques;
      SCIPfreeBlockMemoryArray(scip, &(*varconflicts)[i].cliques, n);
   }
   SCIPfreeBlockMemoryArray(scip, varconflicts, nvars);

   return SCIP_OKAY;
}

/** returns whether a permutation is an involution */
static
SCIP_Bool isInvolution(
   int*                  perm,               /**< permutation */
   int                   lenperm,            /**< length of permutation */
   SCIP_Bool*            istransposition     /**< pointer to store whether permutation is a transposition */
   )
{
   int lensupport = 0;
   int i;

   assert(perm != NULL);
   assert(lenperm > 0);
   assert(istransposition != NULL);
   *istransposition = FALSE;

   for( i = 0; i < lenperm; ++i )
   {
      if( perm[perm[i]] != i )
         return FALSE;
      if( perm[i] != i )
         ++lensupport;
   }

   if( lensupport == 2 )
      *istransposition = TRUE;
   return TRUE;
}

/** tries to generate involutions based on permutations in symmetry component
 *
 *  An involution is a permutation whose 2-fold application is the identity. Involutions
 *  are of particular interest, because their support is (in practice) often very small.
 *  Propagations based on involutions thus can be executed rather quickly. Moreover,
 *  when computing stabilizer subgroups by a filtering mechanism, it is more likely that
 *  a (sparse) involution is not filtered in contrast to dense non-involutions.
 *
 *  To create new involutions, we are given a list of involutions that have been found
 *  during symmetry detection. We then iterate through all pairs (p,q) of involutions in
 *  this list and generate the new involutions p * q (if p and q commute) and p * q * p
 *  (if p and q do not commute).
 */
static
SCIP_RETCODE tryGenerateInvolutions(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_SYMTYPE           symtype,            /**< symmetry type */
   int**                 perms,              /**< (signed) permutations matrix */
   int                   nperms,             /**< number of permutations */
   int                   npermvars,          /**< number of variables the permutations act on */
   int                   maxnnewinvolus,     /**< maximum number of involutions to be computed */
   int***                newperms,           /**< pointer to store new (signed) permutations */
   int*                  nnewperms,          /**< pointer to store number of new permutations */
   int*                  lennewperms         /**< pointer to store length of newperms */
   )
{
   int* tmpperm;
   int* perm1;
   int* perm2;
   int permlen;
   int p;
   int q;
   int i;
   SCIP_Bool commute;
   SCIP_Bool istransposition;
   SCIP_Bool dynamicmemsize;

   /* check whether we shall run */
   if( maxnnewinvolus <= 0 )
      return SCIP_OKAY;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms > 0);
   assert(npermvars > 0);
   assert(newperms != NULL);
   assert(nnewperms != NULL);
   assert(lennewperms != NULL);

   dynamicmemsize = maxnnewinvolus == -1;
   *lennewperms = dynamicmemsize ? 100 : maxnnewinvolus;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, newperms, *lennewperms) );
   *nnewperms = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpperm, npermvars) );

   permlen = symtype == (int)SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;

   /* try to generate new involutions by combining two involutions p and q
    *
    * If p and q commute, create involution p*q. Otherwise, create involutions p*q*p and q*p*q.
    */
   for( p = 0; p < nperms; ++p )
   {
      perm1 = perms[p];
      if( !isInvolution(perm1, npermvars, &istransposition) )
         continue;

      /* it seems promising to only combine involutions with at least two cycles each */
      if ( istransposition )
         continue;

      for( q = p + 1; q < nperms; ++q )
      {
         perm2 = perms[q];
         if( !isInvolution(perm2, npermvars, &istransposition) )
            continue;
         if( istransposition )
            continue;

         /* check whether perm1 and perm2 commute */
         commute = TRUE;
         for( i = 0; i < npermvars && commute; ++i )
         {
            if( perm1[perm2[i]] != perm2[perm1[i]] )
               commute = FALSE;
         }

         /* only consider involutions that have non-disjoint support */
         if( commute )
         {
            /* permutations commute, store perm1 * perm2 if we do not know it yet */
            for( i = 0; i < npermvars; ++i )
               tmpperm[i] = perm1[perm2[i]];

            if( isPermKnown(tmpperm, npermvars, perms, nperms) )
               continue;

            if ( isPermKnown(tmpperm, npermvars, *newperms, *nnewperms) )
               continue;

            if( dynamicmemsize && *nnewperms == *lennewperms )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, newperms, lennewperms, *lennewperms + 1) );
            }
            assert( *nnewperms < *lennewperms );

            /* recompute permutation, because we possibly also need the entries for negated variables */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*newperms)[*nnewperms], permlen) );
            for( i = 0; i < permlen; ++i )
               (*newperms)[*nnewperms][i] = perm1[perm2[i]];
            ++(*nnewperms);
         }
         else
         {
            /* permutations do not commute, compute perm1 * perm2 * perm1 */
            for( i = 0; i < npermvars; ++i )
               tmpperm[i] = perm1[perm2[perm1[i]]];

            /* do not store the permutation if it is already known
             * (also for signed permutations, it is sufficient to iterate over the first npermvars
             *  entries, because the permutation on the negated variables can be derived from these entries)
             */
            if( isPermKnown(tmpperm, npermvars, perms, nperms) )
               continue;

            if( isPermKnown(tmpperm, npermvars, *newperms, *nnewperms) )
               continue;

            /* we do not know the permutation yet, store it */
            if( dynamicmemsize && *nnewperms == *lennewperms )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, newperms, lennewperms, *lennewperms + 1) );
            }
            assert(*nnewperms < *lennewperms);

            /* recompute permutation, because we possibly also need the entries for negated variables */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*newperms)[*nnewperms], permlen) );
            for( i = 0; i < permlen; ++i )
               (*newperms)[*nnewperms][i] = perm1[perm2[perm1[i]]];
            ++(*nnewperms);

            if( !dynamicmemsize && *nnewperms == *lennewperms )
               break;

            /* compute perm2 * perm1 * perm2 */
            for( i = 0; i < npermvars; ++i )
               tmpperm[i] = perm2[perm1[perm2[i]]];

            /* do not store the permutation if it is already known */
            if( isPermKnown(tmpperm, npermvars, perms, nperms) )
               continue;

            if ( isPermKnown(tmpperm, npermvars, *newperms, *nnewperms) )
               continue;

            /* recompute permutation, because we possibly also need the entries for negated variables */
            if( dynamicmemsize && *nnewperms == *lennewperms )
            {
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, newperms, lennewperms, *lennewperms + 1) );
            }

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*newperms)[*nnewperms], permlen) );
            for( i = 0; i < permlen; ++i )
               (*newperms)[*nnewperms][i] = perm2[perm1[perm2[i]]];
            ++(*nnewperms);
         }

         if( !dynamicmemsize && *nnewperms >= *lennewperms )
            break;
      }
      if( !dynamicmemsize && *nnewperms >= *lennewperms )
         break;
   }

   SCIPfreeBufferArray(scip, &tmpperm);

   if( *nnewperms == 0 )
   {
      SCIPfreeBlockMemoryArray(scip, newperms, *lennewperms);
   }

   return SCIP_OKAY;
}

/** selection rule of next orbit/leader in orbit for Schreier-Sims constraints */
static
SCIP_RETCODE selectOrbitLeaderSSTConss(
   SCIP*                 scip,               /**< SCIP instance */
   SST_CONFLICTDATA*     varconflicts,       /**< variable conflicts structure, or NULL if we do not use it */
   SCIP_VAR**            conflictvars,       /**< variables encoded in conflict graph */
   int                   nconflictvars,      /**< number of variables encoded in conflict graph */
   int*                  orbits,             /**< orbits of stabilizer subgroup, expressed in terms of conflictvars */
   int*                  orbitbegins,        /**< array storing the begin position of each orbit in orbits */
   int                   norbits,            /**< number of orbits */
   int                   leaderrule,         /**< rule to select leader */
   int                   orbitrule,          /**< rule to select orbit */
   SCIP_VARTYPE          leadervartype,      /**< variable type of leader */
   int*                  orbitidx,           /**< pointer to index of selected orbit */
   int*                  leaderidx,          /**< pointer to leader in orbit */
   SCIP_Shortbool*       orbitvarinconflict, /**< array to store whether a var in the orbit is conflicting with leader */
   int*                  norbitvarinconflict,/**< pointer to store number of vars in the orbit in conflict with leader */
   SCIP_Bool*            success             /**< pointer to store whether orbit cut be selected successfully */
   )
{
   int varidx;
   int orbitcriterion;
   int curcriterion = INT_MIN;
   int orbitsize;
   int i;
   int leader = -1;

   assert(scip != NULL);
   assert(conflictvars != NULL);
   assert(nconflictvars > 0);
   assert(orbits != NULL);
   assert(orbitbegins != NULL);
   assert(norbits > 0);
   assert(orbitidx != NULL);
   assert(leaderidx != NULL);
   assert(orbitvarinconflict != NULL || varconflicts == NULL);
   assert(norbitvarinconflict != NULL);
   assert(success != NULL);

   *orbitidx = 0;
   *leaderidx = 0;
   *norbitvarinconflict = 0;
   *success = FALSE;

   /* terminate if leader or orbit rule cannot be checked */
   if( varconflicts == NULL && (leaderrule == (int) SCIP_LEADERRULE_MAXCONFLICTSINORBIT
         || orbitrule == (int) SST_ORBITRULE_MAXCONFSINORBIT) )
      return SCIP_OKAY;

   /* select the leader and its orbit */
   if( leaderrule == (int) SST_LEADERRULE_FIRSTINORBIT || leaderrule == (int) SST_LEADERRULE_LASTINORBIT )
   {
      orbitcriterion = INT_MIN;

      /* iterate over orbits and select the first one that meets the orbit rule */
      for( i = 0; i < norbits; ++i )
      {
         /* skip orbits containing vars different to the leader's vartype */
         /* Conflictvars is permvars! */
         if( SCIPgetSymInferredVarType(conflictvars[orbits[orbitbegins[i]]]) != leadervartype )
            continue;

         if( orbitrule == (int) SST_ORBITRULE_MINORBIT )
            curcriterion = orbitbegins[i] - orbitbegins[i + 1];
         else if( orbitrule == (int) SST_ORBITRULE_MAXORBIT )
            curcriterion = orbitbegins[i + 1] - orbitbegins[i];
         else
         {
            assert(orbitrule == (int) SST_ORBITRULE_MAXCONFSINORBIT);

            /* get first or last active variable in orbit */
            if( leaderrule == (int) SST_LEADERRULE_FIRSTINORBIT )
            {
               int cnt;

               cnt = orbitbegins[i];

               do
               {
                  varidx = orbits[cnt++];
               }
               while( SCIPvarGetProbindex(conflictvars[varidx]) == -1 && cnt < orbitbegins[i + 1]);
            }
            else
            {
               int cnt;

               cnt = orbitbegins[i + 1] - 1;

               do
               {
                  varidx = orbits[cnt--];
               }
               while( SCIPvarGetProbindex(conflictvars[varidx]) == -1 && cnt >= orbitbegins[i]);
            }

            /* skip inactive variables */
            if( SCIPvarGetProbindex(conflictvars[varidx]) == -1 )
               continue;

            assert(varconflicts[varidx].orbitidx == i);
            /* coverity[var_deref_op] */
            curcriterion = varconflicts[varidx].nconflictinorbit;
         }

         /* update selected orbit */
         if( curcriterion > orbitcriterion )
         {
            orbitcriterion = curcriterion;
            *orbitidx = i;
            *success = TRUE;

            if( leaderrule == (int) SST_LEADERRULE_FIRSTINORBIT )
               *leaderidx = 0;
            else
               *leaderidx = orbitbegins[i + 1] - orbitbegins[i] - 1;
         }
      }

      /* store variables in conflict with leader */
      if( *success && varconflicts != NULL )
      {
         leader = orbits[orbitbegins[*orbitidx] + *leaderidx];
         assert(leader < nconflictvars);

         if( orbitrule == (int) SST_ORBITRULE_MAXCONFSINORBIT && varconflicts[leader].ncliques > 0 )
         {
            /* Count how many active variables in the orbit conflict with "leader".
             * This is only needed if there are possible conflicts. */
            int varmapid;

            orbitsize = orbitbegins[*orbitidx + 1] - orbitbegins[*orbitidx];
            assert(varconflicts != NULL);
            assert(leader >= 0 && leader < nconflictvars);

            assert(orbitvarinconflict != NULL);

            for( i = 0; i < orbitsize; ++i )
            {
               /* skip the leader */
               if( i == *leaderidx )
                  continue;

               /* get variable index in conflict graph */
               varmapid = orbits[orbitbegins[*orbitidx] + i];

               /* only active variables */
               if( ! varconflicts[varmapid].active )
                  continue;

               /* check if leader and var have overlap */
               if( checkSortedArraysHaveOverlappingEntry((void**)varconflicts[leader].cliques,
                  varconflicts[leader].ncliques, (void**)varconflicts[varmapid].cliques,
                  varconflicts[varmapid].ncliques, sortByPointerValue) )
               {
                  /* there is overlap! */
                  orbitvarinconflict[i] = TRUE;
                  ++(*norbitvarinconflict);
               }
            }
         }
      }
   }
   else
   {
      /* Only three possible values for leaderrules, so it must be MAXCONFLICTSINORBIT.
       * In this case, the code must have computed the conflict graph. */
      assert(leaderrule == (int) SST_LEADERRULE_MAXCONFSINORBIT);
      assert(varconflicts != NULL);

      orbitcriterion = 0;

      /* iterate over variables and select the first one that meets the orbit rule */
      for( i = 0; i < nconflictvars; ++i )
      {
         /* skip vars different to the leader's vartype */
         if( SCIPgetSymInferredVarType(conflictvars[i]) != leadervartype )
            continue;

         /* skip variables not affected by symmetry */
         /* coverity[var_deref_op] */
         if( varconflicts[i].orbitidx == -1 )
            continue;

         curcriterion = varconflicts[i].nconflictinorbit;

         if( curcriterion > orbitcriterion )
         {
            orbitcriterion = curcriterion;
            *orbitidx = varconflicts[i].orbitidx;
            *leaderidx = varconflicts[i].posinorbit;
            *success = TRUE;
         }
      }

      /* store variables in conflict with leader */
      leader = orbits[orbitbegins[*orbitidx] + *leaderidx];
      assert(leader < nconflictvars);
      assert(norbitvarinconflict != NULL);

      if( *success && varconflicts[leader].ncliques > 0 )
      {
         /* count how many active variables in the orbit conflict with leader */
         int varmapid;

         orbitsize = orbitbegins[*orbitidx + 1] - orbitbegins[*orbitidx];
         assert(varconflicts != NULL);
         assert(leader >= 0 && leader < nconflictvars);

         assert(orbitvarinconflict != NULL);

         for( i = 0; i < orbitsize; ++i )
         {
            /* skip the leader */
            if( i == *leaderidx )
               continue;

            /* get variable index in conflict graph */
            varmapid = orbits[orbitbegins[*orbitidx] + i];

            /* only active variables */
            if( ! varconflicts[varmapid].active )
               continue;

            /* check if leader and var have overlap */
            if( checkSortedArraysHaveOverlappingEntry((void**)varconflicts[leader].cliques,
               varconflicts[leader].ncliques, (void**)varconflicts[varmapid].cliques,
               varconflicts[varmapid].ncliques, sortByPointerValue) )
            {
               /* there is overlap! */
               orbitvarinconflict[i] = TRUE;
               ++(*norbitvarinconflict);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** stores information about variables contained in the orbits of added SST cuts */
static
SCIP_RETCODE storeOrbitInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            permvars,           /**< variables permutations act on */
   SCIP_VAR***           orbitvars,          /**< pointer to array for storing variables contained in SST cuts */
   int**                 orbitbegins,        /**< pointer to array for storing begin positions of orbits */
   int*                  norbits,            /**< pointer to store number of orbits */
   int*                  lenorbitvars,       /**< pointer to store length of orbitvars array */
   int*                  lenorbitbegins,     /**< pointer to store length of orbitbeginstostore array */
   int*                  orbit,              /**< orbit of current round of SST cuts */
   int                   lenorbit,           /**< length of orbit */
   int                   orbitleaderidx      /**< position of leader within the orbit */
   )
{
   int nvars;
   int pos;
   int i;

   assert(scip != NULL);
   assert(permvars != NULL);
   assert(orbitvars != NULL);
   assert(orbitbegins != NULL);
   assert(norbits != NULL);
   assert(lenorbitbegins != NULL);
   assert(orbit != NULL);
   assert(lenorbit > 0);
   assert(0 <= orbitleaderidx && orbitleaderidx < lenorbit);

   nvars = *lenorbitbegins == 0 ? 0 : (*orbitbegins)[*norbits];

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitvars, lenorbitvars, nvars + lenorbit) );
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, orbitbegins, lenorbitbegins, *norbits + 2) );

   /* update orbitbegins array */
   if( *norbits == 0 )
      (*orbitbegins)[0] = 0;
   (*orbitbegins)[*norbits + 1] = (*orbitbegins)[*norbits] + lenorbit;

   /* store variables of orbit, the leader is the first variable */
   pos = (*orbitbegins)[*norbits];
   (*orbitvars)[pos++] = permvars[orbit[orbitleaderidx]];
   for( i = 0; i < lenorbit; ++i )
   {
      if( i == orbitleaderidx )
         continue;
      (*orbitvars)[pos++] = permvars[orbit[i]];
   }
   ++(*norbits);

   return SCIP_OKAY;
}

/** tries to add SST constraints */
static
SCIP_RETCODE tryAddSSTConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SYMTYPE           symtype,            /**< type of symmetries considered */
   int**                 perms,              /**< (signed) permutations for which SST constraints shall be added */
   int                   nperms,             /**< number of signed permutations */
   SCIP_VAR**            permvars,           /**< variables on which the (signed) permutations act */
   int                   npermvars,          /**< number of variables */
   SCIP_Real*            permvardomaincenter,/**< array of centers of variable domains */
   SCIP_HASHMAP*         permvarmap,         /**< map of variables to indices in permvars array */
   int                   id,                 /**< identifier of symmetry component */
   SCIP_CONS***          sstconss,           /**< pointer to hold SST constraints */
   int*                  nsstconss,          /**< pointer to store number of SST constraints */
   int*                  maxnsstconss,       /**< pointer to store length of sstconss array */
   SCIP_VAR***           orbitvars,          /**< pointer to array storing variables of orbits used for SST cuts */
   int**                 orbitbeginstostore, /**< pointer to array storing begin positions of orbits */
   int*                  norbitstostore,     /**< pointer to store number of orbits that are stored*/
   int*                  lenorbitvars,       /**< pointer to store length of orbitvars array */
   int*                  lenorbitbegins,     /**< pointer to store length of orbitbegins array */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SST_LEADERRULE        leaderrule,         /**< rule for selecting leaders */
   SST_ORBITRULE         orbitrule,          /**< rule for selecting orbits */
   SST_VARTYPE           leadervartype,      /**< rule for selecting variable type of leaders */
   SCIP_Bool             computenewperms,    /**< Shall additional permutations of symmetry component be computed? */
   int                   maxnnewperms,       /**< maximum number of additional permutations of symmetry component
                                              *   that is computed (used to have better approximation of stabilizer) */
   SCIP_Bool             addconflictcuts,    /**< Shall SST conss be added when a conflict-based rule is used? */
   SCIP_Bool             mixedcomponents,    /**< Shall SST conss be added if symmetry component contains
                                              *   different variable types? */
   SCIP_Bool*            success             /**< pointer to store whether SST constraints are created */
   )
{
   SST_CONFLICTDATA* varconflicts = NULL;
   int** permstrans;
   int nmovedpermvars;
   int nmovedbinpermvars;
   int nmovedintpermvars;
   int nmovedcontpermvars;
   int ntotalperms;
   int* orbits;
   int* orbitbegins;
   int norbits;
   int orbitidx;
   int orbitleaderidx;
   SCIP_Shortbool* orbitvarinconflict = NULL;
   int norbitvarinconflict;
   SCIP_Shortbool* inactiveperms;
   int ninactiveperms;
   int posleader;
   SCIP_VARTYPE selectedtype = SCIP_VARTYPE_CONTINUOUS;
   int nvarsselectedtype;
   SCIP_Bool conflictgraphcreated = FALSE;
   int norbitleadercomponent;
   SCIP_Shortbool* isaffected;
   SCIP_Shortbool* isproperperm;
   int i;
   int p;

   assert(scip != NULL);
   assert(permvars != NULL);
   assert(npermvars > 0);
   assert(nperms > 0);
   assert(orbitvars != NULL);
   assert(orbitbeginstostore != NULL);
   assert(norbitstostore != NULL);
   assert(lenorbitbegins != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( nchgbds != NULL )
      *nchgbds = 0;

   /* get number of affected vars */
   SCIP_CALL( computeAffectedVartypesCount(perms, nperms, permvars, npermvars,
         &nmovedbinpermvars, &nmovedintpermvars, &nmovedcontpermvars) );
   nmovedpermvars = nmovedbinpermvars + nmovedintpermvars + nmovedcontpermvars;
   assert(nmovedpermvars > 0);  /* nperms > 0 implies this */

   /* determine the leader's vartype */
   nvarsselectedtype = 0;
   if( ISSSTBINACTIVE(leadervartype) && nmovedbinpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_BINARY;
      nvarsselectedtype = nmovedbinpermvars;
   }

   if( ISSSTINTACTIVE(leadervartype) && nmovedintpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_INTEGER;
      nvarsselectedtype = nmovedintpermvars;
   }

   if( ISSSTCONTACTIVE(leadervartype) && nmovedcontpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_CONTINUOUS;
      nvarsselectedtype = nmovedcontpermvars;
   }

   /* terminate if no variables of a possible leader type is affected */
   if( nvarsselectedtype == 0 )
      return SCIP_OKAY;

   ntotalperms = nperms;
   if( computenewperms )
   {
      int** newperms = NULL;
      int nnewperms = 0;
      int lennewperms = 0;
      int permlen = 0;

      SCIP_CALL( tryGenerateInvolutions(scip, symtype, perms, nperms, npermvars, maxnnewperms,
            &newperms, &nnewperms, &lennewperms) );
      ntotalperms += nnewperms;

      SCIP_CALL( createPermsTranspose(scip, symtype, perms, nperms, npermvars, newperms, nnewperms, &permstrans) );

      /* compute whether permutation is a proper permutation, i.e., not signed */
      SCIP_CALL( SCIPallocBufferArray(scip, &isproperperm, ntotalperms) );
      for( p = 0; p < nperms; ++p )
      {
         isproperperm[p] = TRUE;
         if( symtype == SYM_SYMTYPE_PERM )
            continue;

         for( i = 0; i < npermvars; ++i )
         {
            if( perms[p][i] >= npermvars )
            {
               isproperperm[p] = FALSE;
               break;
            }
         }
      }

      for( p = 0; p < nnewperms; ++p )
      {
         isproperperm[nperms + p] = TRUE;
         if( symtype == SYM_SYMTYPE_PERM )
            continue;

         for( i = 0; i < npermvars; ++i )
         {
            if( newperms[p][i] >= npermvars )
            {
               isproperperm[nperms + p] = FALSE;
               break;
            }
         }
      }

      /* free additional permutations again */
      permlen = symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;
      for( p = nnewperms - 1; p >= 0; --p )
      {
         SCIPfreeBlockMemoryArray(scip, &newperms[p], permlen);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &newperms, lennewperms);
   }
   else
   {
      SCIP_CALL( createPermsTranspose(scip, symtype, perms, nperms, npermvars, NULL, 0, &permstrans) );

      /* compute whether permutation is a proper permutation, i.e., not signed */
      SCIP_CALL( SCIPallocBufferArray(scip, &isproperperm, nperms) );
      for( p = 0; p < nperms; ++p )
      {
         isproperperm[p] = TRUE;
         if( symtype == SYM_SYMTYPE_PERM )
            continue;

         for( i = 0; i < npermvars; ++i )
         {
            if( perms[p][i] >= npermvars )
            {
               isproperperm[p] = FALSE;
               break;
            }
         }
      }
   }
   assert(permstrans != NULL);

   /* compute variables that are affected by symmetry in component */
   SCIP_CALL( SCIPallocBufferArray(scip, &isaffected, npermvars) );
   for( i = 0; i < npermvars; ++i )
   {
      isaffected[i] = FALSE;
      for( p = 0; p < ntotalperms; ++p )
      {
         if( permstrans[i][p] != i )
         {
            isaffected[i] = TRUE;
            break;
         }
      }
   }

   /* @todo only create the conflict graph for the variables in the current component */
   /* possibly create conflict graph; graph is not created if no cliques are present */
   if( selectedtype == SCIP_VARTYPE_BINARY && (leaderrule == SST_LEADERRULE_MAXCONFSINORBIT
         || orbitrule == SST_ORBITRULE_MAXCONFSINORBIT) )
   {
      SCIP_CALL( createConflictGraphSST(scip, &varconflicts, permvars, npermvars, permvarmap) );
      conflictgraphcreated = varconflicts != NULL;
   }

   /* allocate data structures necessary for orbit computations and conflict graph */
   SCIP_CALL( SCIPallocBufferArray(scip, &inactiveperms, ntotalperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );

   if( conflictgraphcreated )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &orbitvarinconflict, npermvars) );
   }

   SCIPdebugMsg(scip, "Start selection of orbits and leaders for Schreier-Sims constraints.\n");
   SCIPdebugMsg(scip, "orbitidx\tleaderidx\torbitsize\n");

   /* initialize array indicating whether permutations shall not be considered for orbit permutations */
   ninactiveperms = 0;
   for( p = 0; p < ntotalperms; ++p )
   {
      /* possibly filter signed permutations */
      inactiveperms[p] = !isproperperm[p];
      if( inactiveperms[p] )
         ++ninactiveperms;
   }

   /* as long as the stabilizer is non-trivial, add Schreier-Sims constraints */
   norbitleadercomponent = 0;
   while( ninactiveperms < ntotalperms )
   {
      int nchanges = 0;

      /* compute orbits w.r.t. active perms */
      SCIP_CALL( SCIPcomputeOrbitsFilterSymNoComp(scip, npermvars, permstrans, ntotalperms, inactiveperms,
            isaffected, orbits, orbitbegins, &norbits, nmovedpermvars) );

      /* stop if we require pure components and a component contains variables of different types */
      if( !mixedcomponents )
      {
         for( p = 0; p < norbits; ++p )
         {
            /* stop if the first element of an orbits has the wrong vartype */
            if( SCIPgetSymInferredVarType(permvars[orbits[orbitbegins[p]]]) != selectedtype )
            {
               success = FALSE;
               break;
            }
         }
      }

      if( !success )
         break;

      /* update symmetry information of conflict graph */
      if( conflictgraphcreated )
      {
         SCIP_CALL( updateSymInfoConflictGraphSST(scip, varconflicts, permvars, npermvars, orbits, orbitbegins,
               norbits) );
      }

      /* possibly adapt the leader and orbit rule */
      if( leaderrule == SST_LEADERRULE_MAXCONFSINORBIT && !conflictgraphcreated )
         leaderrule = SST_LEADERRULE_FIRSTINORBIT;
      if( leaderrule == SST_LEADERRULE_MAXCONFSINORBIT && selectedtype != SCIP_VARTYPE_BINARY )
         leaderrule = SST_LEADERRULE_FIRSTINORBIT;
      if( orbitrule == SST_ORBITRULE_MAXCONFSINORBIT && !conflictgraphcreated )
         orbitrule = SST_ORBITRULE_MAXORBIT;
      if( orbitrule == SST_ORBITRULE_MAXCONFSINORBIT && selectedtype != SCIP_VARTYPE_BINARY )
         orbitrule = SST_ORBITRULE_MAXORBIT;

      /* select orbit and leader */
      SCIP_CALL( selectOrbitLeaderSSTConss(scip, varconflicts, permvars, npermvars, orbits, orbitbegins, norbits,
            leaderrule, orbitrule, selectedtype, &orbitidx, &orbitleaderidx, orbitvarinconflict, &norbitvarinconflict,
            success) );

      if( !success )
         break;

      assert(0 <= orbitidx && orbitidx < norbits);
      assert(0 <= orbitleaderidx && orbitleaderidx < orbitbegins[orbitidx + 1] - orbitbegins[orbitidx]);
      SCIPdebugMsg(scip, "%d\t\t%d\t\t%d\n", orbitidx, orbitleaderidx, orbitbegins[orbitidx + 1] - orbitbegins[orbitidx]);

      /* add Schreier-Sims constraints for the selected orbit and update Schreier-Sims table */
      SCIP_CALL( addSSTConssOrbitAndUpdateSST(scip, symtype, id, varconflicts, ntotalperms, permvars, npermvars,
            permvardomaincenter, leaderrule, orbitrule, addconflictcuts, sstconss, nsstconss, maxnsstconss,
            orbits, orbitbegins, orbitidx, orbitleaderidx, orbitvarinconflict, norbitvarinconflict, &nchanges) );

      SCIP_CALL( storeOrbitInformation(scip, permvars, orbitvars, orbitbeginstostore, norbitstostore, lenorbitvars,
            lenorbitbegins, &orbits[orbitidx], orbitbegins[orbitidx + 1] - orbitbegins[orbitidx], orbitleaderidx) );

      ++norbitleadercomponent;

      if( nchgbds != NULL )
         *nchgbds += nchanges;

      /* deactivate permutations that move the orbit leader */
      posleader = orbits[orbitbegins[orbitidx] + orbitleaderidx];
      for( p = 0; p < ntotalperms; ++p )
      {
         if( inactiveperms[p] )
            continue;

         if( permstrans[posleader][p] != posleader )
         {
            inactiveperms[p] = TRUE;
            ++ninactiveperms;
         }
      }
   }

   /* if Schreier-Sims constraints have been added, store that Schreier-Sims has been used for this component */
   if( norbitleadercomponent > 0 )
      *success = TRUE;

   if( conflictgraphcreated )
   {
      SCIPfreeBufferArray(scip, &orbitvarinconflict);
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   if( varconflicts != NULL )
   {
      /* nconflictvars at construction is npermvars */
      SCIP_CALL( freeConflictGraphSST(scip, &varconflicts, npermvars) );
   }
   SCIPfreeBufferArray(scip, &inactiveperms);

   SCIPfreeBufferArray(scip, &isproperperm);
   SCIPfreeBufferArray(scip, &isaffected);

   for( i = 0; i < (symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars); ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &permstrans[i], ntotalperms);
   }
   SCIPfreeBlockMemoryArray(scip, &permstrans, symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars);

   return SCIP_OKAY;
}

/*
 * Callback methods of symmetry handler
 */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryAddSST)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_CONS** sstconss = NULL;
   int nsstconss = 0;
   int maxnsstconss = 0;
   SCIP_VAR** orbitvars = NULL;
   int* orbitbegins = NULL;
   int lenorbitvars = 0;
   int lenorbitbegins = 0;
   int norbits = 0;
   int c;

   assert(success != NULL);
   assert(symhdlr != NULL);
   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms >= 0);
   assert(permvars != NULL || npermvars == 0);
   assert(nchgbds != NULL);

   *success = FALSE;
   *nchgbds = 0;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   /* try to add SST constraints */
   SCIP_CALL( tryAddSSTConss(scip, symtype, perms, nperms, permvars, npermvars, permvardomcenter, permvarmap, id,
         &sstconss, &nsstconss, &maxnsstconss, &orbitvars, &orbitbegins, &norbits, &lenorbitvars, &lenorbitbegins,
         nchgbds, (SST_LEADERRULE)symhdlrdata->leaderrule, (SST_ORBITRULE)symhdlrdata->orbitrule,
         (SST_VARTYPE)symhdlrdata->leadervartype, symhdlrdata->computenewperms, symhdlrdata->maxnnewperms,
         symhdlrdata->addconflictcuts, symhdlrdata->mixedcomponents, success) );

   /* in case of success, store information in symmetry component's data */
   if( !(*success) )
   {
      assert(sstconss == NULL);
      assert(nsstconss == 0);
      assert(maxnsstconss == 0);
      assert(orbitvars == NULL);
      assert(orbitbegins == NULL);
      assert(norbits == 0);
      assert(lenorbitvars == 0);
      assert(lenorbitbegins == 0);
      assert(*nchgbds == 0);

      return SCIP_OKAY;
   }

   assert(symcompdata != NULL);
   assert(naddedconss != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symcompdata)->sstconss, nsstconss) );
   for( c = 0; c < nsstconss; ++c )
      (*symcompdata)->sstconss[c] = sstconss[c];
   (*symcompdata)->nsstconss = nsstconss;
   (*symcompdata)->maxnsstconss = nsstconss;
   (*symcompdata)->orbitvars = orbitvars;
   (*symcompdata)->orbitbegins = orbitbegins;
   (*symcompdata)->norbits = norbits;
   (*symcompdata)->lenorbitvars = lenorbitvars;
   (*symcompdata)->lenorbitbegins = lenorbitbegins;
   *naddedconss = nsstconss;

   SCIPfreeBlockMemoryArray(scip, &sstconss, maxnsstconss);

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symhdlrExitSST)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata->sstconss == NULL || symdata->nsstconss > 0);

      for( c = 0; c < symdata->nsstconss; ++c )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &symdata->sstconss[c]) );
      }

      SCIPfreeBlockMemoryArray(scip, &symdata->orbitbegins, symdata->lenorbitbegins);
      SCIPfreeBlockMemoryArray(scip, &symdata->orbitvars, symdata->lenorbitvars);
      SCIPfreeBlockMemoryArrayNull(scip, &symdata->sstconss, symdata->maxnsstconss);
      SCIPfreeBlockMemory(scip, &symdata);
   }

   return SCIP_OKAY;
}

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
static
SCIP_DECL_SYMHDLRFREE(symhdlrFreeSST)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   assert(scip != NULL);
   assert(symhdlr != NULL);

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &symhdlrdata);

   return SCIP_OKAY;
}

/** presolving method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRESOL(symhdlrPresolSST)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   assert(result != NULL);

   *result = nsymcomps > 0 ? SCIP_DIDNOTFIND : SCIP_DIDNOTRUN;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = SCIPsymcompGetData(symcomps[s]);
      assert(symdata->nsstconss > 0);

      for( c = 0; c < symdata->nsstconss; ++c )
      {
         SCIP_CALL( SCIPpresolCons(scip, symdata->sstconss[c], nrounds, presoltiming, nnewfixedvars, nnewaggrvars,
               nnewchgvartypes, nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs,
               nnewchgsides, nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes, ndelconss, naddconss,
               nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff or unboundedness has been detected */
         if( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n",
               SCIPconsGetName(symdata->sstconss[c]));
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** symmetry component display method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRINT(symhdlrPrintSST)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symcompdata;
   int i;
   int v;

   assert(symcomp != NULL);

   symcompdata = SCIPsymcompGetData(symcomp);
   assert(symcompdata != NULL);

   assert(symcompdata->orbitvars != NULL);
   assert(symcompdata->orbitbegins != NULL);

   SCIPinfoMessage(scip, file, "handled by %d SST cuts\n", symcompdata->nsstconss);
   for( i = 0; i < symcompdata->norbits; ++i )
   {
      SCIPinfoMessage(scip, file, "\torbit %5d: leader <%s>; followers", i,
         SCIPvarGetName(symcompdata->orbitvars[symcompdata->orbitbegins[i]]));

      for( v = symcompdata->orbitbegins[i] + 1; v < symcompdata->orbitbegins[i + 1]; ++v )
         SCIPinfoMessage(scip, file, " <%s>", SCIPvarGetName(symcompdata->orbitvars[v]));
      SCIPinfoMessage(scip, file, "\n");
   }

   return SCIP_OKAY;
}

/** include symmetry handler for symresack constraints */
SCIP_RETCODE SCIPincludeSymhdlrSST(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &symhdlrdata) );

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, SYM_NAME, SYM_DESC, SYM_PRIORITY, 0, 0, SYM_PRESOLPRIORITY,
         -1, -1, FALSE, FALSE, 1.0, SYM_MAXPRESOLROUNDS, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST,
         symhdlrTryAddSST, NULL, symhdlrFreeSST, NULL, symhdlrExitSST, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, symhdlrPresolSST, symhdlrPrintSST, symhdlrdata) );

   /* add parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "symmetries/" SYM_NAME "/leaderrule",
         "rule to select the leader within the orbit (0: first, 1: last, 2: max. number of conflicts in orbit)",
         &symhdlrdata->leaderrule, TRUE, DEFAULT_LEADERRULE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "symmetries/" SYM_NAME "/orbitrule",
         "rule to select the orbit for SST cuts (0: min. size, 1: max. size, 2: max. number of conflicts)",
         &symhdlrdata->orbitrule, TRUE, DEFAULT_ORBITRULE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "symmetries/" SYM_NAME "/leadervartype",
         "bitset encoding allowed variable types for leader  (1: bin, 2: int, 4: cont)",
         &symhdlrdata->leadervartype, TRUE, DEFAULT_LEADERVARTYPE, 0, 7, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/computenewperms",
         "Shall additional permutations of symmetry component be computed?",
         &symhdlrdata->computenewperms, TRUE, DEFAULT_COMPUTENEWPERMS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "symmetries/" SYM_NAME "/maxnnewperms",
         "maximum number of additional permutations of symmetry component that is computed " \
         "(used to have better approximation of stabilizer); -1: unbounded",
         &symhdlrdata->maxnnewperms, TRUE, DEFAULT_MAXNNEWPERMS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/addconflictcuts",
         "Should SST constraints be added is a conflict-based rule is used?",
         &symhdlrdata->addconflictcuts, TRUE, DEFAULT_ADDCONFLICTCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "symmetries/" SYM_NAME "/mixedcomponents",
         "Should SST constraints be added is a symmetry component contains variables of different types?",
         &symhdlrdata->mixedcomponents, TRUE, DEFAULT_MIXEDCOMPONENTS, NULL, NULL) );

   return SCIP_OKAY;
}
