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

/**@file   prop_symmetry.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  propagator for handling symmetries
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 * @author Fabian Wegscheider
 * @author Jasper van Doornmalen
 *
 * This propagator combines the following symmetry handling functionalities:
 * - It allows to compute symmetries of the problem and to store this information in adequate form. The symmetry
 *   information can be accessed through external functions.
 * - It implements various methods to handle the symmetries:
 *    - orbital reduction, which generalizes orbital fixing. See symmetry_orbital.c
 *    - (dynamic) orbitopal reduction, which generalizes (dynamic) orbital fixing. See symmetry_orbitopal.c
 *    - static orbitopal fixing (for binary variable domains) for full orbitopes. See cons_orbitope.c
 *    - static orbitopal fixing (for binary variable domains) for packing-partitioning orbitopes. See cons_orbitope.c
 *    - (dynamic) lexicographic reduction. See symmetry_lexred.c
 *    - static lexicographic fixing for binary variable domains (i.e., symresack propagation). See cons_symresack.c
 *    - static lexicographic fixing for binary variable domains on involutions (i.e., orbisacks). See cons_orbisack.c
 *    - Symmetry breaking inequalities based on the Schreier-Sims Table (i.e., SST cuts).
 *    - Strong and weak symmetry breaking inequalities.
 *
 *
 * @section SYMCOMP Symmetry Computation
 *
 * The generic functionality of the compute_symmetry.h interface is used.
 * We do not copy symmetry information, since it is not clear how this information transfers. Moreover, copying
 * symmetry might inhibit heuristics. But note that solving a sub-SCIP might then happen without symmetry information!
 *
 *
 * @section SYMBREAK Symmetry handling by the (unified) symmetry handling constraints
 *
 * Many common methods are captured by a framework that dynamifies symmetry handling constraints. The ideas are
 * described in@n
 * J. van Doornmalen, C. Hojny, "A Unified Framework for Symmetry Handling", preprint, 2023,
 * https://doi.org/10.48550/arXiv.2211.01295.
 *
 * This paper shows that various symmetry handling methods are compatible under certain conditions, and provides
 * generalizations to common symmetry handling constraints from binary variable domains to arbitrary variable domains.
 * This includes symresack propagation, orbitopal fixing, and orbital fixing, that are generalized to
 * lexicographic reduction, orbitopal reduction and orbital reduction, respectively. For a description and
 * implementation, see symmetry_lexred.c, symmetry_orbitopal.c and symmetry_orbital.c, respectively.
 * The static counterparts on binary variable domains are cons_symresack.c and cons_orbisack.c for lexicographic
 * reduction (cf. symresack propagation), and cons_orbitope.c and cons_orbisack.c for orbitopal reduction
 * (cf. orbitopal fixing). We refer to the description of tryAddSymmetryHandlingMethods for the order in which these
 * methods are applied.
 *
 * @section SST Cuts derived from the Schreier Sims table
 *
 * SST cuts have been introduced by@n
 * D. Salvagnin: Symmetry Breaking Inequalities from the Schreier-Sims table. CPAIOR 2018 Proceedings, 521-529, 2018.
 *
 * These inequalities are computed as follows. Throughout these procedure a set of so-called leaders is maintained.
 * Initially the set of leaders is empty. In a first step, select a variable \f$x_i\f$ and compute its orbit w.r.t.
 * the symmetry group of the mixed-integer program. For each variable \f$x_j\f$ in the orbit of \f$x_i\f$, the
 * inequality \f$x_i \geq x_j\f$ is a valid symmetry handling inequality, which can be added to the mixed-integer
 * program. We call \f$x_i\f$ the leader of this inequality. Add the leader \f$x_i\f$ to the set of leaders and
 * compute the pointwise stabilizer of the leader set. In the next step, select a new variable, compute its orbit
 * w.r.t. the stabilizer group of the leaders, add the inequalities based on this orbit, and add the new leader
 * to the set of leaders. This procedure is iterated until the pointwise stabilizer group of the leaders has become
 * trivial.
 *
 * @todo Possibly turn off propagator in subtrees.
 * @todo Check application of conflict resolution.
 * @todo Check whether one should switch the role of 0 and 1
 * @todo Implement stabilizer computation?
 * @todo Implement isomorphism pruning?
 * @todo Implement particular preprocessing rules
 * @todo Separate permuted cuts (first experiments not successful)
 * @todo Allow the computation of local symmetries
 * @todo Order rows of orbitopes (in particular packing/partitioning) w.r.t. cliques in conflict graph.
 * @todo A dynamic variant for packing-partitioning orbitopal structures
 * @todo A dynamic variant for suborbitopes
 */
/* #define SCIP_OUTPUT */
/* #define SCIP_OUTPUT_COMPONENT */
/* #define SCIP_DISPLAY_SYM_CHECK */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
#include <scip/cons_knapsack.h>
#include <scip/cons_varbound.h>
#include <scip/cons_setppc.h>
#include <scip/cons_and.h>
#include <scip/cons_logicor.h>
#include <scip/cons_or.h>
#include <scip/cons_orbitope.h>
#include <scip/cons_symresack.h>
#include <scip/cons_xor.h>
#include <scip/cons_linking.h>
#include <scip/cons_bounddisjunction.h>
#include <scip/cons_indicator.h>
#include <scip/cons_nonlinear.h>
#include <scip/cons_sos1.h>
#include <scip/cons_sos2.h>
#include <scip/expr_pow.h>
#include <scip/expr_product.h>
#include <scip/pub_expr.h>
#include <scip/misc.h>
#include <scip/scip_datastructures.h>

#include <scip/prop_symmetry.h>
#include <symmetry/compute_symmetry.h>
#include <scip/event_shadowtree.h>
#include <scip/symmetry.h>
#include <scip/symmetry_graph.h>
#include <scip/symmetry_orbitopal.h>
#include <scip/symmetry_orbital.h>
#include <scip/symmetry_lexred.h>

#include <math.h>
#include <string.h>

/* propagator properties */
#define PROP_NAME            "symmetry"
#define PROP_DESC            "propagator for handling symmetry"
#define PROP_TIMING    SCIP_PROPTIMING_BEFORELP   /**< propagation timing mask */
#define PROP_PRIORITY          -1000000           /**< propagator priority */
#define PROP_FREQ                     1           /**< propagator frequency */
#define PROP_DELAY                FALSE           /**< should propagation method be delayed, if other propagators found reductions? */

#define PROP_PRESOL_PRIORITY  -10000000           /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers) */
#define PROP_PRESOLTIMING   SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1           /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */


/* default parameter values for symmetry computation */
#define DEFAULT_MAXGENERATORS        1500    /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
#define DEFAULT_CHECKSYMMETRIES     FALSE    /**< Should all symmetries be checked after computation? */
#define DEFAULT_DISPLAYNORBITVARS   FALSE    /**< Should the number of variables affected by some symmetry be displayed? */
#define DEFAULT_USECOLUMNSPARSITY   FALSE    /**< Should the number of conss a variable is contained in be exploited in symmetry detection? */
#define DEFAULT_DOUBLEEQUATIONS     FALSE    /**< Double equations to positive/negative version? */
#define DEFAULT_COMPRESSSYMMETRIES   TRUE    /**< Should non-affected variables be removed from permutation to save memory? */
#define DEFAULT_COMPRESSTHRESHOLD     0.5    /**< Compression is used if percentage of moved vars is at most the threshold. */
#define DEFAULT_SYMFIXNONBINARYVARS FALSE    /**< Disabled parameter */
#define DEFAULT_ENFORCECOMPUTESYMMETRY FALSE /**< always compute symmetries, even if they cannot be handled */
#define DEFAULT_SYMTYPE (int) SYM_SYMTYPE_PERM /**< type of symmetries to be computed */

/* default parameters for linear symmetry constraints */
#define DEFAULT_CONSSADDLP           TRUE    /**< Should the symmetry breaking constraints be added to the LP? */
#define DEFAULT_ADDSYMRESACKS        TRUE    /**< Add inequalities for symresacks for each generator? */
#define DEFAULT_DETECTDOUBLELEX      TRUE    /**< Should we check whether the components of the symmetry group can be handled by double lex matrices? */
#define DEFAULT_DETECTORBITOPES      TRUE    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
#define DEFAULT_DETECTSUBGROUPS      TRUE    /**< Should we try to detect orbitopes in subgroups of the symmetry group? */
#define DEFAULT_ADDWEAKSBCS          TRUE    /**< Should we add weak SBCs for enclosing orbit of symmetric subgroups? */
#define DEFAULT_ADDSTRONGSBCS       FALSE    /**< Should we add strong SBCs for enclosing orbit of symmetric subgroups if orbitopes are not used? */
#define DEFAULT_ADDCONSSTIMING          2    /**< timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving) */
#define DEFAULT_MAXNCONSSSUBGROUP  500000    /**< Maximum number of constraints up to which subgroup structures are detected */
#define DEFAULT_USEDYNAMICPROP       TRUE    /**< whether dynamic propagation should be used for full orbitopes */
#define DEFAULT_PREFERLESSROWS       TRUE    /**< Shall orbitopes with less rows be preferred in detection? */

/* default parameters for symmetry computation */
#define DEFAULT_SYMCOMPTIMING           2    /**< timing of symmetry computation (0 = before presolving, 1 = during presolving, 2 = at first call) */
#define DEFAULT_PERFORMPRESOLVING       0    /**< Run orbital fixing during presolving? (disabled parameter) */
#define DEFAULT_RECOMPUTERESTART        0    /**< Recompute symmetries after a restart has occurred? (0 = never) */

/* default parameters for Schreier Sims constraints */
#define DEFAULT_SSTTIEBREAKRULE   1          /**< index of tie break rule for selecting orbit for Schreier Sims constraints? */
#define DEFAULT_SSTLEADERRULE     0          /**< index of rule for selecting leader variables for Schreier Sims constraints? */
#define DEFAULT_SSTLEADERVARTYPE 14          /**< bitset encoding which variable types can be leaders (1: binary; 2: integer; 4: impl. int; 8: continuous);
                                              *   if multiple types are allowed, take the one with most affected vars */
#define DEFAULT_ADDCONFLICTCUTS       TRUE   /**< Should Schreier Sims constraints be added if we use a conflict based rule? */
#define DEFAULT_SSTADDCUTS            TRUE   /**< Should Schreier Sims constraints be added? */
#define DEFAULT_SSTMIXEDCOMPONENTS    TRUE   /**< Should Schreier Sims constraints be added if a symmetry component contains variables of different types? */

/* output table properties */
#define TABLE_NAME_SYMMETRY     "symmetry"
#define TABLE_DESC_SYMMETRY     "symmetry handling statistics"
#define TABLE_POSITION_SYMMETRY 7001                    /**< the position of the statistics table */
#define TABLE_EARLIEST_SYMMETRY SCIP_STAGE_SOLVING      /**< output of the statistics table is only printed from this stage onwards */


/* other defines */
#define MAXGENNUMERATOR          64000000    /**< determine maximal number of generators by dividing this number by the number of variables */
#define COMPRESSNVARSLB             25000    /**< lower bound on the number of variables above which compression could be performed */

/* macros for getting activeness of symmetry handling methods */
#define ISSYMRETOPESACTIVE(x)      (((unsigned) x & SYM_HANDLETYPE_SYMBREAK) != 0)
#define ISORBITALREDUCTIONACTIVE(x) (((unsigned) x & SYM_HANDLETYPE_ORBITALREDUCTION) != 0)
#define ISSSTACTIVE(x)             (((unsigned) x & SYM_HANDLETYPE_SST) != 0)

#define ISSSTBINACTIVE(x)          (((unsigned) x & SCIP_SSTTYPE_BINARY) != 0)
#define ISSSTINTACTIVE(x)          (((unsigned) x & SCIP_SSTTYPE_INTEGER) != 0)
#define ISSSTIMPLINTACTIVE(x)      (((unsigned) x & SCIP_SSTTYPE_IMPLINT) != 0)
#define ISSSTCONTACTIVE(x)         (((unsigned) x & SCIP_SSTTYPE_CONTINUOUS) != 0)

/* enable symmetry statistics */
#define SYMMETRY_STATISTICS 1

/** propagator data */
struct SCIP_PropData
{
   /* symmetry group information */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nbinpermvars;       /**< number of binary variables for permutations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int                   nperms;             /**< number of permutations */
   int                   nmaxperms;          /**< maximal number of permutations (needed for freeing storage) */
   int**                 perms;              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   int**                 permstrans;         /**< pointer to store transposed permutation generators as (npermvars x nperms) matrix */
   SCIP_HASHMAP*         permvarmap;         /**< map of variables to indices in permvars array */
   int                   nmovedpermvars;     /**< number of variables moved by any permutation */
   int                   nmovedbinpermvars;  /**< number of binary variables moved by any permutation */
   int                   nmovedintpermvars;  /**< number of integer variables moved by any permutation */
   int                   nmovedimplintpermvars; /**< number of implicitly integer variables moved by any permutation */
   int                   nmovedcontpermvars; /**< number of continuous variables moved by any permutation */
   SCIP_HASHMAP*         customsymopnodetypes; /**< types of operator nodes introduced
                                                *   by a user for symmetry detection */
   int                   nopnodetypes;       /**< current number of operator node types used for symmetry detection */
   SCIP_Real*            permvardomaincenter; /**< center of variable domains (needed for signed permutations) */
   int                   symtype;             /**< type of symmetries to be computed */

   /* components of symmetry group */
   int                   ncomponents;        /**< number of components of symmetry group */
   int                   ncompblocked;       /**< number of components that have been blocked */
   int*                  components;         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins;    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent;     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   unsigned*             componentblocked;   /**< array to store which symmetry methods have been applied to a component using
                                              *   the same bitset as for misc/usesymmetry */
   SCIP_Bool*            componenthassignedperm; /**< array to indicate whether a component has a signed permutation */

   /* further symmetry information */
   int                   nmovedvars;         /**< number of variables moved by some permutation */
   SCIP_Real             log10groupsize;     /**< log10 of size of symmetry group */
   SCIP_Bool             binvaraffected;     /**< whether binary variables are affected by some symmetry */

   /* for symmetry computation */
   int                   maxgenerators;      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
   SCIP_Bool             checksymmetries;    /**< Should all symmetries be checked after computation? */
   SCIP_Bool             displaynorbitvars;  /**< Whether the number of variables in non-trivial orbits shall be computed */
   SCIP_Bool             compresssymmetries; /**< Should non-affected variables be removed from permutation to save memory? */
   SCIP_Real             compressthreshold;  /**< Compression is used if percentage of moved vars is at most the threshold. */
   SCIP_Bool             compressed;         /**< Whether symmetry data has been compressed */
   SCIP_Bool             computedsymmetry;   /**< Have we already tried to compute symmetries? */
   int                   usesymmetry;        /**< encoding of active symmetry handling methods (for debugging) */
   SCIP_Bool             usecolumnsparsity;  /**< Should the number of conss a variable is contained in be exploited in symmetry detection? */
   SCIP_Bool             doubleequations;    /**< Double equations to positive/negative version? */
   SCIP_Bool             enforcecomputesymmetry; /**< always compute symmetries, even if they cannot be handled */

   /* for symmetry constraints */
   SCIP_Bool             triedaddsymmethods; /**< whether we already tried to add symmetry handling methods */
   SCIP_Bool             conssaddlp;         /**< Should the symmetry breaking constraints be added to the LP? */
   SCIP_Bool             addsymresacks;      /**< Add symresack constraints for each generator? */
   int                   addconsstiming;     /**< timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving) */
   SCIP_CONS**           genorbconss;        /**< list of generated orbitope/orbisack/symresack constraints */
   SCIP_CONS**           genlinconss;        /**< list of generated linear constraints */
   int                   ngenorbconss;       /**< number of generated orbitope/orbisack/symresack constraints */
   int                   genorbconsssize;    /**< size of generated orbitope/orbisack/symresack constraints array */
   int                   ngenlinconss;       /**< number of generated linear constraints */
   int                   genlinconsssize;    /**< size of linear constraints array */
   int                   nsymresacks;        /**< number of symresack constraints */
   SCIP_Bool             detectdoublelex;    /**< Should we check whether the components of the symmetry group can be handled by double lex matrices? */
   SCIP_Bool             detectorbitopes;    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
   SCIP_Bool             detectsubgroups;    /**< Should we try to detect orbitopes in subgroups of the symmetry group? */
   SCIP_Bool             addweaksbcs;        /**< Should we add weak SBCs for enclosing orbit of symmetric subgroups? */
   SCIP_Bool             addstrongsbcs;      /**< Should we add strong SBCs for enclosing orbit of symmetric subgroups if orbitopes are not used? */
   int                   norbitopes;         /**< number of orbitope constraints */
   SCIP_Bool*            isnonlinvar;        /**< array indicating whether variables appear non-linearly */
   SCIP_CONSHDLR*        conshdlr_nonlinear; /**< nonlinear constraint handler */
   int                   maxnconsssubgroup;  /**< maximum number of constraints up to which subgroup structures are detected */
   SCIP_Bool             usedynamicprop;     /**< whether dynamic propagation should be used for full orbitopes */
   SCIP_Bool             preferlessrows;     /**< Shall orbitopes with less rows be preferred in detection? */

   /* data necessary for symmetry computation order */
   int                   recomputerestart;   /**< Recompute symmetries after a restart has occurred? (0 = never, 1 = always, 2 = if symmetry reduction found) */
   int                   symcomptiming;      /**< timing for computation symmetries (0 = before presolving, 1 = during presolving, 2 = at first call) */
   int                   lastrestart;        /**< last restart for which symmetries have been computed */
   SCIP_Bool             symfoundreduction;  /**< whether symmetry handling propagation has found a reduction since the last time computing symmetries */

   /* data necessary for Schreier Sims constraints */
   SCIP_CONS**           sstconss;           /**< list of generated schreier sims conss */
   int                   nsstconss;          /**< number of generated schreier sims conss */
   int                   maxnsstconss;       /**< maximum number of conss in sstconss */
   int                   sstleaderrule;      /**< rule to select leader  */
   int                   ssttiebreakrule;    /**< tie break rule for leader selection */
   int                   sstleadervartype;   /**< bitset encoding which variable types can be leaders;
                                              *   if multiple types are allowed, take the one with most affected vars */
   int*                  leaders;            /**< index of orbit leaders in permvars */
   int                   nleaders;           /**< number of orbit leaders in leaders array */
   int                   maxnleaders;        /**< maximum number of leaders in leaders array */
   SCIP_Bool             addconflictcuts;    /**< Should Schreier Sims constraints be added if we use a conflict based rule? */
   SCIP_Bool             sstaddcuts;         /**< Should Schreier Sims constraints be added? */
   SCIP_Bool             sstmixedcomponents; /**< Should Schreier Sims constraints be added if a symmetry component contains variables of different types? */

   SCIP_EVENTHDLR*       shadowtreeeventhdlr;/**< pointer to event handler for shadow tree */
   SCIP_ORBITOPALREDDATA* orbitopalreddata;  /**< container for the orbitopal reduction data */
   SCIP_ORBITALREDDATA*  orbitalreddata;     /**< container for orbital reduction data */
   SCIP_LEXREDDATA*      lexreddata;         /**< container for lexicographic reduction propagation */
};

/** conflict data structure for SST cuts */
struct SCIP_ConflictData
{
   SCIP_VAR*             var;                /**< variable belonging to node */
   int                   orbitidx;           /**< orbit of variable w.r.t. current stabilizer subgroup
                                              *   or -1 if not affected by symmetry */
   int                   nconflictinorbit;   /**< number of variables the node's var is in conflict with */
   int                   orbitsize;          /**< size of the variable's orbit */
   int                   posinorbit;         /**< position of variable in its orbit */
   SCIP_Bool             active;             /**< whether variable has not been fixed by Schreier Sims code */
   SCIP_CLIQUE**         cliques;            /**< List of setppc constraints. */
   int                   ncliques;           /**< Number of setppc constraints. */
};
typedef struct SCIP_ConflictData SCIP_CONFLICTDATA;


/** compare function for sorting an array by the addresses of its members  */
static
SCIP_DECL_SORTPTRCOMP(sortByPointerValue)
{
   /* @todo move to misc.c? */
   if ( elem1 < elem2 )
      return -1;
   else if ( elem1 > elem2 )
      return +1;
   return 0;
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

   assert( arr1 != NULL || narr1 == 0 );
   assert( narr1 >= 0 );
   assert( arr2 != NULL || narr2 == 0 );
   assert( narr2 >= 0 );
   assert( compfunc != NULL );

   /* there is no overlap if one of the two arrays is empty */
   if ( narr1 <= 0 )
      return FALSE;
   if ( narr2 <= 0 )
      return FALSE;

   it1 = 0;
   it2 = 0;

   while ( TRUE )  /*lint !e716*/
   {
      cmp = compfunc(arr1[it1], arr2[it2]);
      if ( cmp < 0 )
      {
         /* comparison function determines arr1[it1] < arr2[it2]
          * increase iterator for arr1
          */
         if ( ++it1 >= narr1 )
            break;
         continue;
      }
      else if ( cmp > 0 )
      {
         /* comparison function determines arr1[it1] > arr2[it2]
          * increase iterator for arr2
          */
         if ( ++it2 >= narr2 )
            break;
         continue;
      }
      else
      {
         /* the entries arr1[it1] and arr2[it2] are the same with respect to the comparison function */
         assert( cmp == 0 );
         return TRUE;
      }
   }

   /* no overlap detected */
   assert( it1 >= narr1 || it2 >= narr2 );
   return FALSE;
}


/*
 * Display dialog callback methods
 */

/** displays the cycle of a symmetry */
static
SCIP_RETCODE displayCycleOfSymmetry(
   SCIP*                 scip,               /**< SCIP pointer */
   int*                  perm,               /**< symmetry */
   SYM_SYMTYPE           symtype,            /**< type of symmetry */
   int                   baseidx,            /**< variable index for which cycle is computed */
   SCIP_Bool*            covered,            /**< allocated array to store covered variables */
   int                   nvars,              /**< number of (non-negated) variables in symmetry */
   SCIP_VAR**            vars                /**< variables on which symmetry acts */
   )
{
   char* string;
   int varidx;
   int j;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( 0 <= baseidx );
   assert( (symtype == SYM_SYMTYPE_PERM && baseidx < nvars) ||
      (symtype == SYM_SYMTYPE_SIGNPERM && baseidx < 2 * nvars) );
   assert( covered != NULL );

   /* skip fixed points or elements already covered in previous cycle */
   if ( perm[baseidx] == baseidx || covered[baseidx] )
      return SCIP_OKAY;

   varidx = baseidx >= nvars ? baseidx - nvars : baseidx;
   string = (char*) SCIPvarGetName(vars[varidx]);
   SCIPinfoMessage(scip, NULL, "  (%s<%s>", baseidx >= nvars ? "negated " : "", string);
   j = perm[baseidx];
   covered[baseidx] = TRUE;
   while ( j != baseidx )
   {
      covered[j] = TRUE;
      varidx = j >= nvars ? j - nvars : j;
      string = (char*) SCIPvarGetName(vars[varidx]);
      SCIPinfoMessage(scip, NULL, ",%s<%s>", j >= nvars ? "negated " : "", string);
      j = perm[j];
   }
   SCIPinfoMessage(scip, NULL, ")\n");

   return SCIP_OKAY;
}

/** displays symmetry information without taking components into account */
static
SCIP_RETCODE displaySymmetriesWithoutComponents(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_Bool* covered;
   SYM_SYMTYPE symtype;
   int* perm;
   int permlen;
   int npermvars;
   int i;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->nperms > 0 );
   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );

   symtype = (SYM_SYMTYPE) propdata->symtype;
   npermvars = propdata->npermvars;
   permlen = symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;

   if ( symtype == SYM_SYMTYPE_SIGNPERM )
      SCIPinfoMessage(scip, NULL, "Display permutations as signed permutations (allowing translations)\n");

   SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, permlen) );

   for (p = 0; p < propdata->nperms; ++p)
   {
      SCIPinfoMessage(scip, NULL, "Permutation %d:\n", p);
      perm = propdata->perms[p];

      for (i = 0; i < permlen; ++i)
      {
         SCIP_CALL( displayCycleOfSymmetry(scip, perm, symtype, i, covered, npermvars, propdata->permvars) );
      }

      for (i = 0; i < permlen; ++i)
         covered[i] = FALSE;
   }

   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}

/** displays symmetry information taking components into account */
static
SCIP_RETCODE displaySymmetriesWithComponents(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_Bool* covered;
   SYM_SYMTYPE symtype;
   int* perm;
   int comppermlen;
   int permlen;
   int npermvars;
   int i;
   int p;
   int c;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->nperms > 0 );
   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->ncomponents > 0 );

   symtype = (SYM_SYMTYPE) propdata->symtype;
   npermvars = propdata->npermvars;
   permlen = symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, permlen) );

   for (c = 0; c < propdata->ncomponents; ++c)
   {
      int cnt;

      SCIPinfoMessage(scip, NULL, "Display symmetries of component %d.\n", c);
      if ( propdata->componenthassignedperm[c] )
         SCIPinfoMessage(scip, NULL, "   Symmetries are displayed as signed permutations (allowing translations).\n");
      else
         SCIPinfoMessage(scip, NULL, "   Symmetries are displayed as permutations.\n");

      comppermlen = propdata->componenthassignedperm[c] ? 2 * npermvars : npermvars;

      for (p = propdata->componentbegins[c], cnt = 0; p < propdata->componentbegins[c + 1]; ++p, ++cnt)
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d:\n", cnt);
         perm = propdata->perms[propdata->components[p]];

         for (i = 0; i < comppermlen; ++i)
         {
            SCIP_CALL( displayCycleOfSymmetry(scip, perm, symtype, i, covered, npermvars, propdata->permvars) );
         }

         for (i = 0; i < comppermlen; ++i)
            covered[i] = FALSE;
      }
   }

   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}

/** dialog execution method for the display symmetry information command */
static
SCIP_DECL_DIALOGEXEC(dialogExecDisplaySymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* add your dialog to history of dialogs that have been executed */
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   propdata = (SCIP_PROPDATA*)SCIPdialogGetData(dialog);
   assert( propdata != NULL );

   if ( propdata->nperms == -1 )
   {
      SCIPinfoMessage(scip, NULL, "Cannot display symmetries. Symmetries have not been computed yet.\n");
   }
   else if ( propdata->nperms == 0 )
   {
      SCIPinfoMessage(scip, NULL, "Cannot display symmetries. No symmetries detected.\n");
   }
   else if ( propdata->ncomponents < 0 )
   {
      SCIP_CALL( displaySymmetriesWithoutComponents(scip, propdata) );
   }
   else
   {
      SCIP_CALL( displaySymmetriesWithComponents(scip, propdata) );
   }

   /* next dialog will be root dialog again */
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/*
 * Table callback methods
 */

/** table data */
struct SCIP_TableData
{
   SCIP_PROPDATA*        propdata;           /** pass data of propagator for table output function */
};


/** output method of symmetry propagator statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputSymmetry)
{
   SCIP_TABLEDATA* tabledata;
   int nred;
   int ncutoff;
   SCIP_Real time;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );
   assert( tabledata->propdata != NULL );

   if ( tabledata->propdata->orbitopalreddata || tabledata->propdata->orbitalreddata
      || tabledata->propdata->lexreddata )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "Symmetry           :\n");
      if ( tabledata->propdata->orbitopalreddata )
      {
         SCIP_CALL( SCIPorbitopalReductionGetStatistics(scip, tabledata->propdata->orbitopalreddata, &nred, &ncutoff) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  orbitopal reduction: %10d reductions applied,"
            " %10d cutoffs\n", nred, ncutoff);
      }
      if ( tabledata->propdata->orbitalreddata )
      {
         SCIP_CALL( SCIPorbitalReductionGetStatistics(scip, tabledata->propdata->orbitalreddata, &nred, &ncutoff) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  orbital reduction: %10d reductions applied,"
            " %10d cutoffs\n", nred, ncutoff);
      }
      if ( tabledata->propdata->lexreddata )
      {
         SCIP_CALL( SCIPlexicographicReductionGetStatistics(scip, tabledata->propdata->lexreddata, &nred, &ncutoff) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  lexicographic red: %10d reductions applied,"
            " %10d cutoffs\n", nred, ncutoff);
      }
      if ( tabledata->propdata->shadowtreeeventhdlr )
      {
         time = SCIPgetShadowTreeEventHandlerExecutionTime(scip, tabledata->propdata->shadowtreeeventhdlr);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  shadow tree time : %10.2f s\n", time);
      }
   }

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeSymmetry)
{
   SCIP_TABLEDATA* tabledata;
   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   SCIPfreeBlockMemory(scip, &tabledata);

   return SCIP_OKAY;
}



/*
 * local data structures
 */

/** data structure to store arrays used for sorting colored component types */
struct SYM_Sortgraphcompvars
{
   int*                  components;         /**< array of components */
   int*                  colors;             /**< array of colors */
};
typedef struct SYM_Sortgraphcompvars SYM_SORTGRAPHCOMPVARS;

/** sorts variable indices according to their corresponding component in the graph
 *
 *  Variables are sorted first by the color of their component and then by the component index.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortGraphCompVars)
{
   SYM_SORTGRAPHCOMPVARS* data;

   data = (SYM_SORTGRAPHCOMPVARS*) dataptr;

   if ( data->colors[ind1] < data->colors[ind2] )
      return -1;
   else if ( data->colors[ind1] > data->colors[ind2] )
      return 1;

   if ( data->components[ind1] < data->components[ind2] )
      return -1;
   if ( data->components[ind1] > data->components[ind2] )
      return 1;

   return 0;
}

/** compares two symmetry detection graphs
 *
 *  Graphs are compared by their number of consnodes, then their constypes, then by their lhs,
 *  then by their rhs, then by their total number of nodes, then by the number of operator nodes,
 *  then by their number of value nodes, and then by their number of edges.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
int compareSymgraphs(
   SCIP*                 scip,               /**< SCIP pointer (or NULL) */
   SYM_GRAPH*            G1,                 /**< first graph in comparison */
   SYM_GRAPH*            G2                  /**< second graph in comparison */
   )
{
   int c;
   int perm1;
   int perm2;

   /* compare the number of constraint nodes */
   if ( G1->nconsnodes < G2->nconsnodes )
      return -1;
   if ( G1->nconsnodes > G2->nconsnodes )
      return 1;

   /* compare the constraint nodes of the two graphs */
   for (c = 0; c < G1->nconsnodes; ++c)
   {
      perm1 = G1->consnodeperm[c];
      perm2 = G2->consnodeperm[c];

      if ( SCIPconsGetHdlr(G1->conss[perm1]) < SCIPconsGetHdlr(G2->conss[perm2]) )
         return -1;
      if ( SCIPconsGetHdlr(G1->conss[perm1]) > SCIPconsGetHdlr(G2->conss[perm2]) )
         return 1;

      /* compare using SCIP functions when SCIP is available */
      if ( scip != NULL )
      {
         if ( SCIPisLT(scip, G1->lhs[perm1], G2->lhs[perm2]) )
            return -1;
         if ( SCIPisGT(scip, G1->lhs[perm1], G2->lhs[perm2]) )
            return 1;

         if ( SCIPisLT(scip, G1->rhs[perm1], G2->rhs[perm2]) )
            return -1;
         if ( SCIPisGT(scip, G1->rhs[perm1], G2->rhs[perm2]) )
            return 1;
      }
      else
      {
         if ( G1->lhs[perm1] < G2->lhs[perm2] )
            return -1;
         if ( G1->lhs[perm1] > G2->lhs[perm2] )
            return 1;

         if ( G1->rhs[perm1] < G2->rhs[perm2] )
            return -1;
         if ( G1->rhs[perm1] > G2->rhs[perm2] )
            return 1;
      }
   }

   /* compare number of remaining node types */
   if ( G1->nnodes < G2->nnodes )
      return -1;
   if ( G1->nnodes > G2->nnodes )
      return 1;

   if ( G1->nopnodes < G2->nopnodes )
      return -1;
   if ( G1->nopnodes > G2->nopnodes )
      return 1;

   if ( G1->nvalnodes < G2->nvalnodes )
      return -1;
   if ( G1->nvalnodes > G2->nvalnodes )
      return 1;

   if ( G1->nedges < G2->nedges )
      return -1;
   if ( G1->nedges > G2->nedges )
      return 1;

   return 0;
}

/** sorts symmetry detection graphs
 *
 *  Graphs are sorted by their number of consnodes, then their constypes, then by their lhs,
 *  then by their rhs, then by their total number of nodes, then by the number of operator nodes,
 *  then by their number of value nodes, and then by their number of edges.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortSymgraphs)
{
   SYM_GRAPH** data;
   SYM_GRAPH* G1;
   SYM_GRAPH* G2;

   data = (SYM_GRAPH**) dataptr;
   G1 = data[ind1];
   G2 = data[ind2];

   return compareSymgraphs(NULL, G1, G2);
}

/*
 * Local methods
 */

#ifndef NDEBUG
/** checks that symmetry data is all freed */
static
SCIP_Bool checkSymmetryDataFree(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert( propdata->permvarmap == NULL );
   assert( propdata->genorbconss == NULL );
   assert( propdata->genlinconss == NULL );
   assert( propdata->ngenlinconss == 0 );
   assert( propdata->ngenorbconss == 0 );
   assert( propdata->genorbconsssize == 0 );
   assert( propdata->genlinconsssize == 0 );
   assert( propdata->sstconss == NULL );
   assert( propdata->leaders == NULL );

   assert( propdata->permvardomaincenter == NULL );
   assert( propdata->permvars == NULL );
   assert( propdata->perms == NULL );
   assert( propdata->permstrans == NULL );
   assert( propdata->npermvars == 0 );
   assert( propdata->nbinpermvars == 0 );
   assert( propdata->nperms == -1 || propdata->nperms == 0 );
   assert( propdata->nmaxperms == 0 );
   assert( propdata->nmovedpermvars == -1 );
   assert( propdata->nmovedbinpermvars == 0 );
   assert( propdata->nmovedintpermvars == 0 );
   assert( propdata->nmovedimplintpermvars == 0 );
   assert( propdata->nmovedcontpermvars == 0 );
   assert( propdata->nmovedvars == -1 );
   assert( propdata->binvaraffected == FALSE );
   assert( propdata->isnonlinvar == NULL );

   assert( propdata->componenthassignedperm == NULL );
   assert( propdata->componentblocked == NULL );
   assert( propdata->componentbegins == NULL );
   assert( propdata->components == NULL );
   assert( propdata->ncomponents == -1 );
   assert( propdata->ncompblocked == 0 );

   return TRUE;
}
#endif


/** resets symmetry handling propagators that depend on the branch-and-bound tree structure */
static
SCIP_RETCODE resetDynamicSymmetryHandling(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   /* propagators managed by a different file */
   if ( propdata->orbitalreddata != NULL )
   {
      SCIP_CALL( SCIPorbitalReductionReset(scip, propdata->orbitalreddata) );
   }
   if ( propdata->orbitopalreddata != NULL )
   {
      SCIP_CALL( SCIPorbitopalReductionReset(scip, propdata->orbitopalreddata) );
   }
   if ( propdata->lexreddata != NULL )
   {
      SCIP_CALL( SCIPlexicographicReductionReset(scip, propdata->lexreddata) );
   }

   return SCIP_OKAY;
}


/** frees symmetry data */
static
SCIP_RETCODE freeSymmetryData(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );

   SCIP_CALL( resetDynamicSymmetryHandling(scip, propdata) );

   if ( propdata->permvarmap != NULL )
   {
      SCIPhashmapFree(&propdata->permvarmap);
   }

   /* release all variables contained in permvars array */
   for (i = 0; i < propdata->npermvars; ++i)
   {
      assert( propdata->permvars[i] != NULL );
      SCIP_CALL( SCIPreleaseVar(scip, &propdata->permvars[i]) );
   }

   /* free permstrans matrix*/
   if ( propdata->permstrans != NULL )
   {
      assert( propdata->nperms > 0 );
      assert( propdata->permvars != NULL );
      assert( propdata->npermvars > 0 );
      assert( propdata->nmaxperms > 0 );

      for (i = 0; i < propdata->npermvars; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &propdata->permstrans[i], propdata->nmaxperms);
      }
      SCIPfreeBlockMemoryArray(scip, &propdata->permstrans, propdata->npermvars);
   }

   /* free data of added orbitope/orbisack/symresack constraints */
   if ( propdata->genorbconss != NULL )
   {
      assert( propdata->ngenorbconss > 0 );

      /* release constraints */
      while ( propdata->ngenorbconss > 0 )
      {
         assert( propdata->genorbconss[propdata->ngenorbconss - 1] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genorbconss[--propdata->ngenorbconss]) );
      }
      assert( propdata->ngenorbconss == 0 );

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->genorbconss, propdata->genorbconsssize);
      propdata->genorbconsssize = 0;
   }

   /* free data of added constraints */
   if ( propdata->genlinconss != NULL )
   {
      /* release constraints */
      for (i = 0; i < propdata->ngenlinconss; ++i)
      {
         assert( propdata->genlinconss[i] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genlinconss[i]) );
      }

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->genlinconss, propdata->genlinconsssize);
      propdata->ngenlinconss = 0;
      propdata->genlinconsssize = 0;
   }

   if ( propdata->sstconss != NULL )
   {
      assert( propdata->nsstconss > 0 );

      /* release constraints */
      for (i = 0; i < propdata->nsstconss; ++i)
      {
         assert( propdata->sstconss[i] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->sstconss[i]) );
      }

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->sstconss, propdata->maxnsstconss);
      propdata->sstconss = NULL;
      propdata->nsstconss = 0;
      propdata->maxnsstconss = 0;
   }

   if ( propdata->leaders != NULL )
   {
      assert( propdata->maxnleaders > 0 );

      SCIPfreeBlockMemoryArray(scip, &propdata->leaders, propdata->maxnleaders);
      propdata->maxnleaders = 0;
      propdata->leaders = NULL;
      propdata->nleaders = 0;
   }

   /* free components */
   if ( propdata->ncomponents > 0 )
   {
      assert( propdata->componentblocked != NULL );
      assert( propdata->vartocomponent != NULL );
      assert( propdata->componentbegins != NULL );
      assert( propdata->components != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->componenthassignedperm, propdata->ncomponents);
      SCIPfreeBlockMemoryArray(scip, &propdata->componentblocked, propdata->ncomponents);
      SCIPfreeBlockMemoryArray(scip, &propdata->vartocomponent, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->componentbegins, propdata->ncomponents + 1);
      SCIPfreeBlockMemoryArray(scip, &propdata->components, propdata->nperms);

      propdata->ncomponents = -1;
      propdata->ncompblocked = 0;
   }

   /* free main symmetry data */
   if ( propdata->nperms > 0 )
   {
      int permlen;

      assert( propdata->permvars != NULL );

      if ( (SYM_SYMTYPE) propdata->symtype == SYM_SYMTYPE_SIGNPERM )
         permlen = 2 * propdata->npermvars;
      else
         permlen = propdata->npermvars;

      SCIPfreeBlockMemoryArray(scip, &propdata->permvars, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->permvardomaincenter, propdata->npermvars);

      if ( propdata->perms != NULL )
      {
         for (i = 0; i < propdata->nperms; ++i)
         {
            SCIPfreeBlockMemoryArray(scip, &propdata->perms[i], permlen);
         }
         SCIPfreeBlockMemoryArray(scip, &propdata->perms, propdata->nmaxperms);
      }

      SCIPfreeBlockMemoryArrayNull(scip, &propdata->isnonlinvar, propdata->npermvars);

      propdata->npermvars = 0;
      propdata->nbinpermvars = 0;
      propdata->nperms = -1;
      propdata->nmaxperms = 0;
      propdata->nmovedpermvars = -1;
      propdata->nmovedbinpermvars = 0;
      propdata->nmovedintpermvars = 0;
      propdata->nmovedimplintpermvars = 0;
      propdata->nmovedcontpermvars = 0;
      propdata->nmovedvars = -1;
      propdata->log10groupsize = -1.0;
      propdata->binvaraffected = FALSE;
      propdata->isnonlinvar = NULL;
   }
   propdata->nperms = -1;

   assert( checkSymmetryDataFree(propdata) );

   propdata->computedsymmetry = FALSE;
   propdata->compressed = FALSE;

   return SCIP_OKAY;
}


/** makes sure that the constraint array (potentially NULL) of given array size is sufficiently large */
static
SCIP_RETCODE ensureDynamicConsArrayAllocatedAndSufficientlyLarge(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS***          consarrptr,         /**< constraint array pointer */
   int*                  consarrsizeptr,     /**< constraint array size pointer */
   int                   consarrsizereq      /**< constraint array size required */
   )
{
   int newsize;

   assert( scip != NULL );
   assert( consarrptr != NULL );
   assert( consarrsizeptr != NULL );
   assert( consarrsizereq > 0 );
   assert( *consarrsizeptr >= 0 );
   assert( (*consarrsizeptr == 0) == (*consarrptr == NULL) );

   /* array is already sufficiently large */
   if ( consarrsizereq <= *consarrsizeptr )
      return SCIP_OKAY;

   /* compute new size */
   newsize = SCIPcalcMemGrowSize(scip, consarrsizereq);
   assert( newsize > *consarrsizeptr );

   /* allocate or reallocate */
   if ( *consarrptr == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, consarrptr, newsize) );
   }
   else
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, consarrptr, *consarrsizeptr, newsize) );
   }

   *consarrsizeptr = newsize;

   return SCIP_OKAY;
}

/** set symmetry data */
static
SCIP_RETCODE setSymmetryData(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries in perms */
   SCIP_VAR**            vars,               /**< vars present at time of symmetry computation */
   int                   nvars,              /**< number of vars present at time of symmetry computation */
   int                   nbinvars,           /**< number of binary vars present at time of symmetry computation */
   SCIP_VAR***           permvars,           /**< pointer to permvars array */
   int*                  npermvars,          /**< pointer to store number of permvars */
   int*                  nbinpermvars,       /**< pointer to store number of binary permvars */
   SCIP_Real**           permvardomaincenter, /**< pointer to store center points of variable domains */
   int**                 perms,              /**< permutations matrix (nperms x nvars) */
   int                   nperms,             /**< number of permutations */
   int*                  nmovedvars,         /**< pointer to store number of vars affected by symmetry (if usecompression) or NULL */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether a binary variable is affected by symmetry */
   SCIP_Bool             usecompression,     /**< whether symmetry data shall be compressed */
   SCIP_Real             compressthreshold,  /**< if percentage of moved vars is at most threshold, compression is done */
   SCIP_Bool*            compressed          /**< pointer to store whether compression has been performed */
   )
{
   SCIP_Real ub;
   SCIP_Real lb;
   int i;
   int p;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( permvars != NULL );
   assert( npermvars != NULL );
   assert( nbinpermvars != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( binvaraffected != NULL );
   assert( SCIPisGE(scip, compressthreshold, 0.0) );
   assert( SCIPisLE(scip, compressthreshold, 1.0) );
   assert( compressed != NULL );

   /* set default return values */
   *permvars = vars;
   *npermvars = nvars;
   *nbinpermvars = nbinvars;
   *binvaraffected = FALSE;
   *compressed = FALSE;

   /* if we possibly perform compression */
   if ( usecompression && SCIPgetNVars(scip) >= COMPRESSNVARSLB )
   {
      SCIP_Real percentagemovedvars;
      int* labelmovedvars;
      int* labeltopermvaridx;
      int nbinvarsaffected = 0;

      assert( nmovedvars != NULL );

      *nmovedvars = 0;

      /* detect number of moved vars and label moved vars */
      SCIP_CALL( SCIPallocBufferArray(scip, &labelmovedvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &labeltopermvaridx, nvars) );
      for (i = 0; i < nvars; ++i)
      {
         labelmovedvars[i] = -1;

         for (p = 0; p < nperms; ++p)
         {
            if ( perms[p][i] != i )
            {
               labeltopermvaridx[*nmovedvars] = i;
               labelmovedvars[i] = (*nmovedvars)++;

               if ( SCIPvarIsBinary(vars[i]) )
                  ++nbinvarsaffected;
               break;
            }
         }
      }

      if ( nbinvarsaffected > 0 )
         *binvaraffected = TRUE;

      /* check whether compression should be performed */
      percentagemovedvars = (SCIP_Real) *nmovedvars / (SCIP_Real) nvars;
      if ( *nmovedvars > 0 && SCIPisLE(scip, percentagemovedvars, compressthreshold) )
      {
         /* remove variables from permutations that are not affected by any permutation */
         for (p = 0; p < nperms; ++p)
         {
            /* iterate over labels and adapt permutation (possibly taking signed permutations into account) */
            for (i = 0; i < *nmovedvars; ++i)
            {
               assert( i <= labeltopermvaridx[i] );
               if ( perms[p][labeltopermvaridx[i]] >= nvars )
               {
                  int imgvaridx;

                  assert( symtype == SYM_SYMTYPE_SIGNPERM );

                  imgvaridx = perms[p][labeltopermvaridx[i]] - nvars;
                  perms[p][i] = labelmovedvars[imgvaridx] + *nmovedvars;

                  assert( 0 <= perms[p][i] && perms[p][i] < 2 * (*nmovedvars) );
               }
               else
                  perms[p][i] = labelmovedvars[perms[p][labeltopermvaridx[i]]];
            }

            if ( symtype == SYM_SYMTYPE_SIGNPERM )
            {
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &perms[p], 2 * nvars, 2 * (*nmovedvars)) );
            }
            else
            {
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &perms[p], nvars, *nmovedvars) );
            }
         }

         /* remove variables from permvars array that are not affected by any symmetry */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, permvars, *nmovedvars) );
         for (i = 0; i < *nmovedvars; ++i)
         {
            (*permvars)[i] = vars[labeltopermvaridx[i]];
         }
         *npermvars = *nmovedvars;
         *nbinpermvars = nbinvarsaffected;
         *compressed = TRUE;

         SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      }
      SCIPfreeBufferArray(scip, &labeltopermvaridx);
      SCIPfreeBufferArray(scip, &labelmovedvars);
   }
   else
   {
      /* detect whether binary variable is affected by symmetry and count number of binary permvars */
      for (i = 0; i < nbinvars; ++i)
      {
         for (p = 0; p < nperms && ! *binvaraffected; ++p)
         {
            if ( perms[p][i] != i )
            {
               if ( SCIPvarIsBinary(vars[i]) )
                  *binvaraffected = TRUE;
               break;
            }
         }
      }
   }

   /* store center points of variable domains */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, permvardomaincenter, *npermvars) );
   for (i = 0; i < *npermvars; ++i)
   {
      ub = SCIPvarGetUbGlobal((*permvars)[i]);
      lb = SCIPvarGetLbGlobal((*permvars)[i]);

      (*permvardomaincenter)[i] = 0.5 * (ub + lb);
   }

   return SCIP_OKAY;
}

/** returns whether a constraint handler can provide required symmetry information */
static
SCIP_Bool conshdlrCanProvideSymInformation(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SYM_SYMTYPE           symtype             /**< type of symmetries for which information are needed */
   )
{
   assert( conshdlr != NULL );

   switch ( symtype )
   {
   case SYM_SYMTYPE_PERM:
      return SCIPconshdlrSupportsPermsymDetection(conshdlr);
   default:
      assert( symtype == SYM_SYMTYPE_SIGNPERM );
      return SCIPconshdlrSupportsSignedPermsymDetection(conshdlr);
   } /*lint !e788*/
}

/** returns whether all constraint handlers with constraints can provide symmetry information */
static
SCIP_Bool conshdlrsCanProvideSymInformation(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype             /**< type of symmetries for which information are needed */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   SCIP_CONSHDLR* conshdlr;
   int nconshdlrs;
   int c;

   conshdlrs = SCIPgetConshdlrs(scip);
   assert( conshdlrs != NULL );

   nconshdlrs = SCIPgetNConshdlrs(scip);
   for (c = 0; c < nconshdlrs; ++c)
   {
      conshdlr = conshdlrs[c];
      assert( conshdlr != NULL );

      if ( ! conshdlrCanProvideSymInformation(conshdlr, symtype) && SCIPconshdlrGetNConss(conshdlr) > 0 )
      {
         char name[SCIP_MAXSTRLEN];

         if ( symtype == SYM_SYMTYPE_PERM )
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "CONSGETPERMSYMGRAPH");
         else
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "CONSGETSIGNEDPERMSYMGRAPH");

         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   Symmetry detection interrupted: constraints of type %s do not provide symmetry information.\n"
            "   If symmetries shall be detected, implement the %s callback.\n",
            SCIPconshdlrGetName(conshdlr), name);

         return FALSE;
      }
   }

   /* check whether all expressions provide sufficient symmetry information */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   if ( conshdlr != NULL && SCIPconshdlrGetNConss(conshdlr) > 0 )
   {
      SCIP_EXPRHDLR* exprhdlr;

      for (c = 0; c < SCIPgetNExprhdlrs(scip); ++c)
      {
         SCIP_Bool found = FALSE;
         exprhdlr = SCIPgetExprhdlrs(scip)[c];

         if ( SCIPexprhdlrHasGetSymData(exprhdlr) )
            continue;

         /* check whether exprhdlr is known by SCIP (and handles symmetries correctly) */
         if ( strcmp(SCIPexprhdlrGetName(exprhdlr), "var") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "sum") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "product") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "val") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "pow") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "signpow") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "exp") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "log") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "abs") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "sin") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "cos") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "entropy") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "erf") == 0
            || strcmp(SCIPexprhdlrGetName(exprhdlr), "varidx") == 0 )
            found = TRUE;

         /* there exists an unknown expression handler that does not provide symmetry information */
         if ( ! found )
         {
            SCIPwarningMessage(scip, "Expression handler %s does not implement the EXPRGETSYMDATA callback.\n"
               "Computed symmetries might be incorrect if the expression uses different constants or assigns\n"
               "different coefficients to its children.\n", SCIPexprhdlrGetName(SCIPgetExprhdlrs(scip)[c]));
         }
      }
   }

   return TRUE;
}

/** provides estimates for the number of nodes and edges in a symmetry detection graph */
static
SCIP_RETCODE estimateSymgraphSize(
   SCIP*                 scip,               /**< SCIP pointer */
   int*                  nopnodes,           /**< pointer to store estimate for number of operator nodes */
   int*                  nvalnodes,          /**< pointer to store estimate for number of value nodes */
   int*                  nconsnodes,         /**< pointer to store estimate for number of constraint nodes */
   int*                  nedges              /**< pointer to store estimate for number of edges */
   )
{
   SCIP_CONS** conss;
   SCIP_Bool success;
   int nvars;
   int nconss;
   int num;
   int c;

   assert( scip != NULL );
   assert( nopnodes != NULL );
   assert( nvalnodes != NULL );
   assert( nconsnodes != NULL );
   assert( nedges != NULL );

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   assert( conss != NULL || nconss == 0 );

   *nconsnodes = nconss;

   /* get estimate from different types of constraints */
   *nopnodes = 0;
   *nvalnodes = 0;
   for (c = 0; c < nconss; ++c)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "bounddisjunction") == 0 )
      {
         SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &num, &success) );

         if ( success )
         {
            *nopnodes += num;
            *nvalnodes += num;
         }
      }
      else if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "indicator") == 0 )
      {
         *nopnodes += 3;
         *nvalnodes += 1;
      }
      else if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "nonlinear") == 0 )
      {
         SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &num, &success) );

         /* use binary trees as a proxy for an expression tree */
         if ( success )
         {
            int depth;
            int numnodes;
            int expval;

            depth = (int) log2((double) num);
            expval = (int) exp2((double) (depth + 1));
            numnodes = MIN(expval, 100);

            *nopnodes += numnodes;
            *nvalnodes += MAX((int) 0.1 * numnodes, 1);
         }
      }
      else if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SOS1") == 0 )
      {
         SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &num, &success) );

         if ( success )
            *nopnodes += num;
      }
      else if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SOS2") == 0 )
      {
         SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &num, &success) );

         if ( success )
            *nopnodes += num - 1;
      }
   }

   /* use a staggered scheme for the number of edges since this can become large
    *
    * In most cases, edges represent variable coefficients from linear constraints.
    * For this reason, use number of variables as proxy.
    */
   if ( nvars <= 100000 )
      *nedges = 100 * nvars;
   else if ( nvars <= 1000000 )
      *nedges = 32 * nvars;
   else if ( nvars <= 16700000 )
      *nedges = 16 * nvars;
   else
      *nedges = INT_MAX / 10;

   return SCIP_OKAY;
}

/** checks whether computed symmetries are indeed symmetries */
static
SCIP_RETCODE checkSymmetriesAreSymmetries(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries to be checked */
   int**                 perms,              /**< array of permutations */
   int                   nperms,             /**< number of permutations */
   int                   npermvars,          /**< number of variables permutations act on */
   SYM_SPEC              fixedtype           /**< variable types that must be fixed by symmetries */
   )
{
   SYM_GRAPH** graphs;
   SCIP_CONS** conss;
   SCIP_VAR** symvars;
   SCIP_Bool success;
   int* graphperm;
   int* groupbegins;
   int ngroups = 1;
   int nsymvars;
   int nconss;
   int p;
   int c;
   int g;
#ifdef SCIP_DISPLAY_SYM_CHECK
   int permlen;
   SCIP_Bool* covered;
#endif

   assert( scip != NULL );
   assert( perms != NULL );
   assert( nperms > 0 );
   assert( npermvars > 0 );

   /* get symmetry detection graphs for all constraints */
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   assert( conss != NULL );

   symvars = SCIPgetVars(scip);
   nsymvars = SCIPgetNVars(scip);
   assert( nsymvars == npermvars );

   SCIP_CALL( SCIPallocBufferArray(scip, &graphs, nconss) );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CALL( SCIPcreateSymgraph(scip, symtype, &graphs[c], symvars, nsymvars, 10, 10, 1, 100) );

      switch ( symtype )
      {
      case SYM_SYMTYPE_PERM:
         SCIP_CALL( SCIPgetConsPermsymGraph(scip, conss[c], graphs[c], &success) );
         break;
      default:
         assert( symtype == SYM_SYMTYPE_SIGNPERM );
         SCIP_CALL( SCIPgetConsSignedPermsymGraph(scip, conss[c], graphs[c], &success) );
      } /*lint !e788*/

      SCIP_CALL( SCIPcomputeSymgraphColors(scip, graphs[c], fixedtype) );

      assert( success );
   }

   /* sort graphs for quicker comparisons */
   SCIP_CALL( SCIPallocBufferArray(scip, &graphperm, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &groupbegins, nconss + 1) );
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CALL( SCIPcreateSymgraphConsnodeperm(scip, graphs[c]) );
   }

   SCIPsort(graphperm, SYMsortSymgraphs, graphs, nconss);

   groupbegins[0] = 0;
   for (c = 1; c < nconss; ++c)
   {
      if ( compareSymgraphs(scip, graphs[graphperm[c]], graphs[graphperm[c-1]]) != 0 )
         groupbegins[ngroups++] = c;
   }
   groupbegins[ngroups] = nconss;

   /* remove information from symmetry detection graph that is not needed anymore */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CALL( SCIPfreeSymgraphConsnodeperm(scip, graphs[c]) );
   }

#ifdef SCIP_DISPLAY_SYM_CHECK
   permlen = symtype == SYM_SYMTYPE_SIGNPERM ? 2 * npermvars : npermvars;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, permlen) );
#endif

   /* iterate over all permutations and check whether they define symmetries */
   for (p = 0; p < nperms; ++p)
   {
      SYM_GRAPH* graph;
      SCIP_Bool found = TRUE;
      int d;
#ifdef SCIP_DISPLAY_SYM_CHECK
      int i;

      SCIPinfoMessage(scip, NULL, "Check whether permutation %d is a symmetry:\n", p);
      for (i = 0; i < permlen; ++i)
      {
         SCIP_CALL( displayCycleOfSymmetry(scip, perms[p], symtype, i, covered, npermvars, SCIPgetVars(scip)) );
      }

      for (i = 0; i < permlen; ++i)
         covered[i] = FALSE;
      SCIPinfoMessage(scip, NULL, "Check whether every constraint has a symmetric counterpart.\n");
#endif

      /* for every constraint, create permuted graph by copying nodes and edges */
      for (g = 0; g < ngroups; ++g)
      {
         for (c = groupbegins[g]; c < groupbegins[g+1]; ++c)
         {
#ifdef SCIP_DISPLAY_SYM_CHECK
            SCIPinfoMessage(scip, NULL, "Check whether constraint %d has a symmetric counterpart:\n",
               graphperm[c]);
            SCIP_CALL( SCIPprintCons(scip, conss[graphperm[c]], NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPcopySymgraph(scip, &graph, graphs[graphperm[c]], perms[p], fixedtype) );

            /* if adapted graph is equivalent to original graph, we don't need to check further graphs */
            if ( SYMcheckGraphsAreIdentical(scip, symtype, graph, graphs[graphperm[c]]) )
            {
#ifdef SCIP_DISPLAY_SYM_CHECK
               SCIPinfoMessage(scip, NULL, "\tconstraint is symmetric to itself\n");
#endif
               SCIP_CALL( SCIPfreeSymgraph(scip, &graph) );
               continue;
            }

            /* check whether graph has an isomorphic counterpart */
            found = FALSE;
            for (d = groupbegins[g]; d < groupbegins[g+1] && ! found; ++d)
            {
               found = SYMcheckGraphsAreIdentical(scip, symtype, graph, graphs[graphperm[d]]);

#ifdef SCIP_DISPLAY_SYM_CHECK
               SCIPinfoMessage(scip, NULL, "\tconstraint is %ssymmetric to constraint %d\n\t", !found ? "not " : "", d);
               SCIP_CALL( SCIPprintCons(scip, conss[graphperm[d]], NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
#endif
            }

            SCIP_CALL( SCIPfreeSymgraph(scip, &graph) );

            if ( ! found )
            {
#ifdef SCIP_DISPLAY_SYM_CHECK
               SCIPfreeBufferArray(scip, &covered);
#endif
               SCIPerrorMessage("permutation %d is not a symmetry\n", p);
               return SCIP_ERROR;
            }
         }
      }
   }

#ifdef SCIP_DISPLAY_SYM_CHECK
   SCIPfreeBufferArray(scip, &covered);
#endif

   SCIPfreeBufferArray(scip, &groupbegins);
   SCIPfreeBufferArray(scip, &graphperm);

   for (c = nconss - 1; c >= 0; --c)
   {
      SCIP_CALL( SCIPfreeSymgraph(scip, &graphs[c]) );
   }
   SCIPfreeBufferArray(scip, &graphs);

   return SCIP_OKAY;
}

/** computes symmetry group of a CIP */
static
SCIP_RETCODE computeSymmetryGroup(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries to be computed */
   SCIP_Bool             compresssymmetries, /**< Should non-affected variables be removed from permutation to save memory? */
   SCIP_Real             compressthreshold,  /**< if percentage of moved vars is at most threshold, compression is done */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_SPEC              fixedtype,          /**< variable types that must be fixed by symmetries */
   SCIP_Bool             checksymmetries,    /**< Should all symmetries be checked after computation? */
   SCIP_VAR***           permvars,           /**< pointer to permvars array */
   int*                  npermvars,          /**< pointer to store number of permvars */
   int*                  nbinpermvars,       /**< pointer to store number of binary permvars */
   SCIP_Real**           permvardomaincenter, /**< pointer to store center points of variable domains */
   int***                perms,              /**< pointer to store permutation matrix (nperms x nvars) */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations
                                              *   (needed for freeing storage) */
   int*                  nmovedvars,         /**< pointer to store number of vars affected
                                              *   by symmetry (if usecompression) or NULL */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether a binary variable is affected by symmetry */
   SCIP_Bool*            compressed,         /**< pointer to store whether compression has been performed */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   SCIP_Real*            symcodetime,        /**< pointer to store the time for symmetry code */
   SCIP_Bool*            success             /**< pointer to store whether symmetry computation was successful */
   )
{
   SCIP_CONS** conss;
   SYM_GRAPH* graph;
   int nconsnodes = 0;
   int nvalnodes = 0;
   int nopnodes = 0;
   int nedges = 0;
   int nconss;
   int c;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( npermvars != NULL );
   assert( nbinpermvars != NULL );
   assert( perms != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( nmovedvars != NULL );
   assert( binvaraffected != NULL );
   assert( compressed != NULL );
   assert( log10groupsize != NULL );
   assert( symcodetime != NULL );
   assert( success != NULL );

   /* init pointers */
   *permvars = NULL;
   *npermvars = 0;
   *nbinpermvars = 0;
   *perms = NULL;
   *nperms = 0;
   *nmaxperms = 0;
   *nmovedvars = -1;
   *binvaraffected = FALSE;
   *compressed = FALSE;
   *log10groupsize = 0;
   *success = FALSE;
   *symcodetime = 0.0;

   /* check whether all constraints can provide symmetry information */
   if ( ! conshdlrsCanProvideSymInformation(scip, symtype) )
      return SCIP_OKAY;

   /* get symmetry detection graphs from constraints */
   conss = SCIPgetConss(scip);
   assert( conss != NULL );

   nconss = SCIPgetNConss(scip);

   /* exit if no constraints or no variables are available */
   if ( nconss == 0 || SCIPgetNVars(scip) == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* get an estimate for the number of nodes and edges */
   SCIP_CALL( estimateSymgraphSize(scip, &nopnodes, &nvalnodes, &nconsnodes, &nedges) );

   /* create graph */
   SCIP_CALL( SCIPcreateSymgraph(scip, symtype, &graph, SCIPgetVars(scip), SCIPgetNVars(scip),
         nopnodes, nvalnodes, nconsnodes, nedges) );

   *success = TRUE;
   for (c = 0; c < nconss && *success; ++c)
   {
      if ( symtype == SYM_SYMTYPE_PERM )
      {
         SCIP_CALL( SCIPgetConsPermsymGraph(scip, conss[c], graph, success) );
      }
      else
      {
         assert( symtype == SYM_SYMTYPE_SIGNPERM );
         SCIP_CALL( SCIPgetConsSignedPermsymGraph(scip, conss[c], graph, success) );
      }

      /* terminate early if graph could not be returned */
      if ( ! *success )
      {
         SCIP_CALL( SCIPfreeSymgraph(scip, &graph) );

         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPcomputeSymgraphColors(scip, graph, fixedtype) );

   /* terminate early in case all variables are different */
   if ( (symtype == SYM_SYMTYPE_PERM && SCIPgetSymgraphNVarcolors(graph) == SCIPgetNVars(scip))
      || (symtype == SYM_SYMTYPE_SIGNPERM && SCIPgetSymgraphNVarcolors(graph) == 2 * SCIPgetNVars(scip)) )
   {
      SCIP_CALL( SCIPfreeSymgraph(scip, &graph) );
      return SCIP_OKAY;
   }

   /*
    * actually compute symmetries
    */
   SCIP_CALL( SYMcomputeSymmetryGenerators(scip, maxgenerators, graph, nperms, nmaxperms,
         perms, log10groupsize, symcodetime) );

   if ( checksymmetries && *nperms > 0 )
   {
      SCIP_CALL( checkSymmetriesAreSymmetries(scip, symtype, *perms, *nperms, SCIPgetNVars(scip), fixedtype) );
   }

   /* potentially store symmetries */
   if ( *nperms > 0 )
   {
      SCIP_VAR** vars;
      int nvars;

      nvars = SCIPgetNVars(scip);
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, SCIPgetVars(scip), nvars) ); /*lint !e666*/

      SCIP_CALL( setSymmetryData(scip, symtype, vars, nvars, SCIPgetNBinVars(scip), permvars, npermvars, nbinpermvars,
            permvardomaincenter, *perms, *nperms, nmovedvars, binvaraffected,
            compresssymmetries, compressthreshold, compressed) );
   }

   /* free symmetry graph */
   SCIP_CALL( SCIPfreeSymgraph(scip, &graph) );

   return SCIP_OKAY;
}

/** returns whether a symmetry is a non-standard permutation */
static
SCIP_Bool isNonstandardPerm(
   SCIP*                 scip,               /**< SCIP instance */
   int*                  symmetry,           /**< a symmetry encoded as a signed permutation */
   SCIP_VAR**            vars,               /**< array of variables the symmetry acts on */
   int                   nvars               /**< number of variables in vars */
   )
{
   int v;

   assert( symmetry != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );

   for (v = 0; v < nvars; ++v)
   {
      /* the symmetry is signed */
      if ( symmetry[v] >= nvars )
         return TRUE;

      /* the domain of symmetric variables is different */
      if ( !SCIPisEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetLbLocal(vars[symmetry[v]]))
         || !SCIPisEQ(scip, SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbLocal(vars[symmetry[v]])) )
      {
         assert( SCIPisEQ(scip, SCIPvarGetUbLocal(vars[v]) - SCIPvarGetLbLocal(vars[v]),
               SCIPvarGetUbLocal(vars[symmetry[v]]) - SCIPvarGetLbLocal(vars[symmetry[v]])) );
         return TRUE;
      }
   }

   return FALSE;
}

/** checks whether component contains non-standard permutations
 *
 *  If all symmetries are standard permutations, stores them as such.
 */
static
SCIP_RETCODE checkComponentsForNonstandardPerms(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int* components;
   int* componentbegins;
   int ncomponents;
   int i;
   int c;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->ncomponents > 0 );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );

   components = propdata->components;
   componentbegins = propdata->componentbegins;
   ncomponents = propdata->ncomponents;

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(propdata->componenthassignedperm), ncomponents) );

   /* stop if no non-standard permutations can exist */
   if ( (SYM_SYMTYPE) propdata->symtype == SYM_SYMTYPE_PERM )
      return SCIP_OKAY;

   /* for each component, check whether it has a non-standard permutation  */
   for (c = 0; c < ncomponents; ++c)
   {
      for (i = componentbegins[c]; i < componentbegins[c + 1]; ++i)
      {
         if ( isNonstandardPerm(scip, propdata->perms[components[i]], propdata->permvars, propdata->npermvars) )
         {
            propdata->componenthassignedperm[c] = TRUE;
            break;
         }
      }
   }

   return SCIP_OKAY;
}

/** ensures that the symmetry components are already computed */
static
SCIP_RETCODE ensureSymmetryComponentsComputed(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   /* symmetries must have been determined */
   assert( propdata->nperms >= 0 );

   /* stop if already computed */
   if ( propdata->ncomponents >= 0 )
      return SCIP_OKAY;

   /* compute components */
   assert( propdata->ncomponents == -1 );
   assert( propdata->components == NULL );
   assert( propdata->componentbegins == NULL );
   assert( propdata->vartocomponent == NULL );

#ifdef SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation started\n", SCIPgetSolvingTime(scip));
#endif

   SCIP_CALL( SCIPcomputeComponentsSym(scip, (SYM_SYMTYPE) propdata->symtype, propdata->perms, propdata->nperms,
         propdata->permvars, propdata->npermvars, FALSE, &propdata->components, &propdata->componentbegins,
         &propdata->vartocomponent, &propdata->componentblocked, &propdata->ncomponents) );

#ifdef SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation finished\n", SCIPgetSolvingTime(scip));
#endif

   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );
   assert( propdata->ncomponents > 0 );

   /* structure of symmetries can be simplified if they are standard permutations */
   SCIP_CALL( checkComponentsForNonstandardPerms(scip, propdata) );
   assert( propdata->componenthassignedperm != NULL );

   return SCIP_OKAY;
}


/** ensures that permvarmap is initialized */
static
SCIP_RETCODE ensureSymmetryPermvarmapComputed(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;

   assert( scip != NULL );
   assert( propdata != NULL );

   /* symmetries must have been determined */
   assert( propdata->nperms >= 0 );

   /* stop if already computed */
   if ( propdata->permvarmap != NULL )
      return SCIP_OKAY;

   /* create hashmap for storing the indices of variables */
   SCIP_CALL( SCIPhashmapCreate(&propdata->permvarmap, SCIPblkmem(scip), propdata->npermvars) );

   /* insert variables into hashmap  */
   for (v = 0; v < propdata->npermvars; ++v)
   {
      SCIP_CALL( SCIPhashmapInsertInt(propdata->permvarmap, propdata->permvars[v], v) );
   }

   return SCIP_OKAY;
}


/** ensures that permstrans is initialized */
static
SCIP_RETCODE ensureSymmetryPermstransComputed(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );

   /* symmetries must have been determined */
   assert( propdata->nperms >= 0 );

   /* stop if already computed */
   if ( propdata->permstrans != NULL )
      return SCIP_OKAY;

   /* transpose symmetries matrix here */
   assert( propdata->permstrans == NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permstrans, propdata->npermvars) );
   for (v = 0; v < propdata->npermvars; ++v)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->permstrans[v]), propdata->nmaxperms) );
      for (p = 0; p < propdata->nperms; ++p)
         propdata->permstrans[v][p] = propdata->perms[p][v];
   }

   return SCIP_OKAY;
}


/** ensures that movedpermvarscounts is initialized */
static
SCIP_RETCODE ensureSymmetryMovedpermvarscountsComputed(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );

   /* symmetries must have been determined */
   assert( propdata->nperms >= 0 );

   /* stop if already computed */
   if ( propdata->nmovedpermvars >= 0 )
      return SCIP_OKAY;
   assert( propdata->nmovedpermvars == -1 );

   propdata->nmovedpermvars = 0;
   propdata->nmovedbinpermvars = 0;
   propdata->nmovedintpermvars = 0;
   propdata->nmovedimplintpermvars = 0;
   propdata->nmovedcontpermvars = 0;

   for (p = 0; p < propdata->nperms; ++p)
   {
      for (v = 0; v < propdata->npermvars; ++v)
      {
         if ( propdata->perms[p][v] != v )
         {
            ++propdata->nmovedpermvars;

            switch ( SCIPvarGetType(propdata->permvars[v]) )
            {
            case SCIP_VARTYPE_BINARY:
               ++propdata->nmovedbinpermvars;
               break;
            case SCIP_VARTYPE_INTEGER:
               ++propdata->nmovedintpermvars;
               break;
            case SCIP_VARTYPE_IMPLINT:
               ++propdata->nmovedimplintpermvars;
               break;
            case SCIP_VARTYPE_CONTINUOUS:
               ++propdata->nmovedcontpermvars;
               break;
            default:
               SCIPerrorMessage("Variable provided with unknown vartype\n");
               return SCIP_ERROR;
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** returns whether any allowed symmetry handling method is effective for the problem instance */
static
SCIP_Bool testSymmetryComputationRequired(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   /* must always compute symmetry if it is enforced */
   if ( propdata->enforcecomputesymmetry )
      return TRUE;

   /* for dynamic symmetry handling or orbital reduction, branching must be possible */
   if ( propdata->usedynamicprop || ISORBITALREDUCTIONACTIVE(propdata->usesymmetry) )
   {
      /* @todo a proper test whether variables can be branched on or not */
      if ( SCIPgetNBinVars(scip) > 0 )
         return TRUE;
      if ( SCIPgetNIntVars(scip) > 0 )
         return TRUE;
      /* continuous variables can be branched on if nonlinear constraints exist */
      if ( ( SCIPgetNContVars(scip) > 0 || SCIPgetNImplVars(scip) > 0 )
         && SCIPconshdlrGetNActiveConss(propdata->conshdlr_nonlinear) > 0 )
         return TRUE;
   }

   /* for SST, matching leadervartypes */
   if ( ISSSTACTIVE(propdata->usesymmetry) )
   {
      if ( ISSSTBINACTIVE(propdata->sstleadervartype) && SCIPgetNBinVars(scip) > 0 ) /*lint !e641*/
         return TRUE;
      if ( ISSSTINTACTIVE(propdata->sstleadervartype) && SCIPgetNIntVars(scip) > 0 ) /*lint !e641*/
         return TRUE;
      if ( ISSSTIMPLINTACTIVE(propdata->sstleadervartype) && SCIPgetNImplVars(scip) > 0 ) /*lint !e641*/
         return TRUE;
      if ( ISSSTCONTACTIVE(propdata->sstleadervartype) && SCIPgetNContVars(scip) > 0 ) /*lint !e641*/
         return TRUE;
   }

   /* for static symmetry handling constraints, binary variables must be present */
   if ( ISSYMRETOPESACTIVE(propdata->usesymmetry) )
   {
      if ( SCIPgetNBinVars(scip) > 0 )
         return TRUE;
   }

   /* if all tests above fail, then the symmetry handling methods cannot achieve anything */
   return FALSE;
}

/** determines symmetry */
static
SCIP_RETCODE determineSymmetry(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SYM_SPEC              symspecrequire,     /**< symmetry specification for which we need to compute symmetries */
   SYM_SPEC              symspecrequirefixed /**< symmetry specification of variables which must be fixed by symmetries */
   )
{ /*lint --e{641}*/
   SCIP_Bool successful;
   SCIP_Real symcodetime = 0.0;
   int maxgenerators;
   unsigned int type = 0;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* do not compute symmetry if reoptimization is enabled */
   if ( SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   /* do not compute symmetry if Benders decomposition enabled */
   if ( SCIPgetNActiveBenders(scip) > 0 )
      return SCIP_OKAY;

   /* skip symmetry computation if no graph automorphism code was linked */
   if ( ! SYMcanComputeSymmetry() )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   Deactivated symmetry handling methods, since SCIP was built without symmetry detector (SYM=none).\n");

      return SCIP_OKAY;
   }

   /* do not compute symmetry if there are active pricers */
   if ( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* avoid trivial cases */
   nvars = SCIPgetNVars(scip);
   if ( nvars <= 0 )
      return SCIP_OKAY;

   /* do not compute symmetry if we cannot handle it */
   if ( !testSymmetryComputationRequired(scip, propdata) )
      return SCIP_OKAY;

   /* determine symmetry specification */
   if ( SCIPgetNBinVars(scip) > 0 )
      type |= (int) SYM_SPEC_BINARY;
   if ( SCIPgetNIntVars(scip) > 0 )
      type |= (int) SYM_SPEC_INTEGER;
   /* count implicit integer variables as real variables, since we cannot currently handle integral variables well */
   if ( SCIPgetNContVars(scip) > 0 || SCIPgetNImplVars(scip) > 0 )
      type |= (int) SYM_SPEC_REAL;

   /* skip symmetry computation if required variables are not present */
   if ( ! (type & symspecrequire) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation skipped: type (bin %c, int %c, cont %c) does not match requirements (bin %c, int %c, cont %c).\n",
         SCIPgetSolvingTime(scip),
         SCIPgetNBinVars(scip) > 0 ? '+' : '-',
         SCIPgetNIntVars(scip) > 0  ? '+' : '-',
         SCIPgetNContVars(scip) + SCIPgetNImplVars(scip) > 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_REAL) != 0 ? '+' : '-');

      return SCIP_OKAY;
   }

   /* skip computation if symmetry has already been computed */
   if ( propdata->computedsymmetry )
      return SCIP_OKAY;

   assert( propdata->npermvars == 0 );
   assert( propdata->permvars == NULL );
   assert( propdata->nperms < 0 );
   assert( propdata->nmaxperms == 0 );
   assert( propdata->perms == NULL );

   /* output message */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   (%.1fs) symmetry computation started: requiring (bin %c, int %c, cont %c), (fixed: bin %c, int %c, cont %c)\n",
      SCIPgetSolvingTime(scip),
      (symspecrequire & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
      (symspecrequire & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
      (symspecrequire & (int) SYM_SPEC_REAL) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_REAL) != 0 ? '+' : '-');

   /* output warning if we want to fix certain symmetry parts that we also want to compute */
   if ( symspecrequire & symspecrequirefixed )
      SCIPwarningMessage(scip, "Warning: some required symmetries must be fixed.\n");

   /* determine maximal number of generators depending on the number of variables */
   maxgenerators = propdata->maxgenerators;
   maxgenerators = MIN(maxgenerators, MAXGENNUMERATOR / nvars);

   /* actually compute (global) symmetry */
   SCIP_CALL( computeSymmetryGroup(scip, (SYM_SYMTYPE) propdata->symtype,
         propdata->compresssymmetries, propdata->compressthreshold,
         maxgenerators, symspecrequirefixed, propdata->checksymmetries, &propdata->permvars,
         &propdata->npermvars, &propdata->nbinpermvars, &propdata->permvardomaincenter,
         &propdata->perms, &propdata->nperms, &propdata->nmaxperms,
         &propdata->nmovedvars, &propdata->binvaraffected, &propdata->compressed,
         &propdata->log10groupsize, &symcodetime, &successful) );

   /* mark that we have computed the symmetry group */
   propdata->computedsymmetry = TRUE;

   /* store restart level */
   propdata->lastrestart = SCIPgetNRuns(scip);

   /* return if not successful */
   if ( ! successful )
   {
      assert( checkSymmetryDataFree(propdata) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) could not compute symmetry\n", SCIPgetSolvingTime(scip));

      return SCIP_OKAY;
   }

   /* return if no symmetries found */
   if ( propdata->nperms == 0 )
   {
      assert( checkSymmetryDataFree(propdata) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry present (symcode time: %.2f)\n", SCIPgetSolvingTime(scip), symcodetime);

      return SCIP_OKAY;
   }
   assert( propdata->nperms > 0 );
   assert( propdata->npermvars > 0 );

   /* display statistics */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) symmetry computation finished: %d generators found (max: ",
      SCIPgetSolvingTime(scip), propdata->nperms);

   /* display statistics: maximum number of generators */
   if ( maxgenerators == 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "-");
   else
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%d", maxgenerators);

   /* display statistics: log10 group size, number of affected vars*/
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", log10 of symmetry group size: %.1f", propdata->log10groupsize);

   if ( propdata->displaynorbitvars )
   {
      if ( propdata->nmovedvars == -1 )
      {
         SCIP_CALL( SCIPdetermineNVarsAffectedSym(scip, propdata->perms, propdata->nperms, propdata->permvars,
               propdata->npermvars, &(propdata->nmovedvars)) );
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", number of affected variables: %d)\n", propdata->nmovedvars);
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ") (symcode time: %.2f)\n", symcodetime);

   /* capture all variables while they are in the permvars array */
   for (i = 0; i < propdata->npermvars; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, propdata->permvars[i]) );
   }

   return SCIP_OKAY;
}


/*
 * Functions for symmetry constraints
 */


/** Checks whether given set of 2-cycle permutations forms an orbitope and if so, builds the variable index matrix.
 *
 *  If @p activevars == NULL, then the function assumes all permutations of the component are active and therefore all
 *  moved vars are considered.
 *
 *  We need to keep track of the number of generated columns, because we might not be able to detect all orbitopes.
 *  For example (1,2), (2,3), (3,4), (3,5) defines the symmetric group on {1,2,3,4,5}, but the generators we expect
 *  in our construction need shape (1,2), (2,3), (3,4), (4,5).
 *
 *  @pre @p orbitopevaridx has to be an initialized 2D array of size @p ntwocycles x @p nperms
 *  @pre @p columnorder has to be an initialized array of size nperms
 *  @pre @p nusedelems has to be an initialized array of size npermvars
 */
static
SCIP_RETCODE checkTwoCyclePermsAreOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            permvars,           /**< array of all permutation variables */
   int                   npermvars,          /**< number of permutation variables */
   int**                 perms,              /**< array of all permutations of the symmetry group */
   int*                  activeperms,        /**< indices of the relevant permutations in perms */
   int                   ntwocycles,         /**< number of 2-cycles in the permutations */
   int                   nactiveperms,       /**< number of active permutations */
   int**                 orbitopevaridx,     /**< pointer to store variable index matrix */
   int*                  columnorder,        /**< pointer to store column order */
   int*                  nusedelems,         /**< pointer to store how often each element was used */
   int*                  nusedcols,          /**< pointer to store number of columns used in orbitope (or NULL) */
   SCIP_Shortbool*       rowisbinary,        /**< pointer to store which rows are binary (or NULL) */
   SCIP_Bool*            isorbitope,         /**< buffer to store result */
   SCIP_Shortbool*       activevars          /**< bitset to store whether a variable is active (or NULL) */
   )
{  /*lint --e{571}*/
   SCIP_Bool* usedperm;
   SCIP_Bool foundperm = FALSE;
   int nusedperms = 0;
   int nfilledcols;
   int coltoextend;
   int ntestedperms = 0;
   int row = 0;
   int j;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( perms != NULL );
   assert( activeperms != NULL );
   assert( orbitopevaridx != NULL );
   assert( columnorder != NULL );
   assert( nusedelems != NULL );
   assert( isorbitope != NULL );
   assert( nactiveperms > 0 );
   assert( ntwocycles > 0 );
   assert( npermvars > 0 );
   assert( activevars == NULL || (0 <= nactiveperms && nactiveperms < npermvars) );

   *isorbitope = TRUE;
   if ( nusedcols != NULL )
      *nusedcols = 0;

   /* whether a permutation was considered to contribute to orbitope */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &usedperm, nactiveperms) );

   /* fill first two columns of orbitopevaridx matrix */

   /* look for the first active permutation which moves an active variable */
   while ( ! foundperm )
   {
      int permidx;

      assert( ntestedperms < nactiveperms );

      permidx = activeperms[ntestedperms];

      for (j = 0; j < npermvars; ++j)
      {
         if ( activevars != NULL && ! activevars[j] )
            continue;

         assert( activevars == NULL || activevars[perms[permidx][j]] );

         /* avoid adding the same 2-cycle twice */
         if ( perms[permidx][j] > j )
         {
            assert( SCIPvarIsBinary(permvars[j]) == SCIPvarIsBinary(permvars[perms[permidx][j]]) );

            if ( rowisbinary != NULL && SCIPvarIsBinary(permvars[j]) )
               rowisbinary[row] = TRUE;

            orbitopevaridx[row][0] = j;
            orbitopevaridx[row++][1] = perms[permidx][j];
            ++(nusedelems[j]);
            ++(nusedelems[perms[permidx][j]]);

            foundperm = TRUE;
         }

         if ( row == ntwocycles )
            break;
      }

      ++ntestedperms;
   }

   /* in the subgroup case it might happen that a generator has less than ntwocycles many 2-cyles */
   if ( row != ntwocycles )
   {
      *isorbitope = FALSE;
      SCIPfreeBufferArray(scip, &usedperm);
      return SCIP_OKAY;
   }

   usedperm[ntestedperms - 1] = TRUE;
   ++nusedperms;
   columnorder[0] = 0;
   columnorder[1] = 1;
   nfilledcols = 2;

   /* extend orbitopevaridx matrix to the left, i.e., iteratively find new permutations that
    * intersect the last added left column in each row in exactly one entry, starting with
    * column 0 */
   coltoextend = 0;
   for (j = ntestedperms; j < nactiveperms; ++j)
   {  /* lint --e{850} */
      SCIP_Bool success = FALSE;
      SCIP_Bool infeasible = FALSE;

      if ( nusedperms == nactiveperms )
         break;

      if ( usedperm[j] )
         continue;

      SCIP_CALL( SCIPextendSubOrbitope(orbitopevaridx, ntwocycles, nfilledcols, coltoextend,
            perms[activeperms[j]], TRUE, &nusedelems, permvars, NULL, &success, &infeasible) );

      if ( infeasible )
      {
         *isorbitope = FALSE;
         break;
      }
      else if ( success )
      {
         usedperm[j] = TRUE;
         ++nusedperms;
         coltoextend = nfilledcols;
         columnorder[nfilledcols++] = -1; /* mark column to be filled from the left */
         j = 0; /*lint !e850*/ /* reset j since previous permutations can now intersect with the latest added column */
      }
   }

   if ( ! *isorbitope ) /*lint !e850*/
   {
      SCIPfreeBufferArray(scip, &usedperm);
      return SCIP_OKAY;
   }

   coltoextend = 1;
   for (j = ntestedperms; j < nactiveperms; ++j)
   {  /*lint --e(850)*/
      SCIP_Bool success = FALSE;
      SCIP_Bool infeasible = FALSE;

      if ( nusedperms == nactiveperms )
         break;

      if ( usedperm[j] )
         continue;

      SCIP_CALL( SCIPextendSubOrbitope(orbitopevaridx, ntwocycles, nfilledcols, coltoextend,
            perms[activeperms[j]], FALSE, &nusedelems, permvars, NULL, &success, &infeasible) );

      if ( infeasible )
      {
         *isorbitope = FALSE;
         break;
      }
      else if ( success )
      {
         usedperm[j] = TRUE;
         ++nusedperms;
         coltoextend = nfilledcols;
         columnorder[nfilledcols] = 1; /* mark column to be filled from the right */
         ++nfilledcols;
         j = 0; /*lint !e850*/ /* reset j since previous permutations can now intersect with the latest added column */
      }
   }

   if ( activevars == NULL && nusedperms < nactiveperms ) /*lint !e850*/
      *isorbitope = FALSE;

   if ( nusedcols != NULL )
      *nusedcols = nfilledcols;
   assert( ! *isorbitope || activevars == NULL || nusedperms < nfilledcols );

   SCIPfreeBufferArray(scip, &usedperm);

   return SCIP_OKAY;
}

/** choose an order in which the generators should be added for subgroup detection */
static
SCIP_RETCODE chooseOrderOfGenerators(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int                   compidx,            /**< index of component */
   int**                 genorder,           /**< (initialized) buffer to store the resulting order of generator */
   int*                  ntwocycleperms      /**< pointer to store the number of 2-cycle permutations in component compidx */
   )
{
   int** perms;
   int* components;
   int* componentbegins;
   int* ntwocycles;
   int npermvars;
   int npermsincomp;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( compidx >= 0 );
   assert( compidx < propdata->ncomponents );
   assert( genorder != NULL );
   assert( *genorder != NULL );
   assert( ntwocycleperms != NULL );
   assert( propdata->computedsymmetry );
   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->ncomponents > 0 );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );

   perms = propdata->perms;
   npermvars = propdata->npermvars;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   npermsincomp = componentbegins[compidx + 1] - componentbegins[compidx];
   *ntwocycleperms = npermsincomp;

   SCIP_CALL( SCIPallocBufferArray(scip, &ntwocycles, npermsincomp) );

   for (i = 0; i < npermsincomp; ++i)
   {
      int* perm;
      int nbincycles;

      perm = perms[components[componentbegins[compidx] + i]];

      SCIP_CALL( SCIPisInvolutionPerm(perm, propdata->permvars, npermvars, &(ntwocycles[i]), &nbincycles, FALSE) );

      /* we skip permutations which do not purely consist of 2-cycles */
      if ( ntwocycles[i] == 0 )
      {
         /* we change the number of two cycles for this perm so that it will be sorted to the end */
         if ( propdata->preferlessrows )
            ntwocycles[i] = npermvars;
         else
            ntwocycles[i] = 0;
         --(*ntwocycleperms);
      }
      else if ( ! propdata->preferlessrows )
         ntwocycles[i] = - ntwocycles[i];
   }

   SCIPsortIntInt(ntwocycles, *genorder, npermsincomp);

   SCIPfreeBufferArray(scip, &ntwocycles);

   return SCIP_OKAY;
}


/** builds the graph for symmetric subgroup detection from the given permutation of generators
 *
 *  After execution, @p graphcomponents contains all permvars sorted by their color and component,
 *  @p graphcompbegins points to the indices where new components in @p graphcomponents start and
 *  @p compcolorbegins points to the indices where new colors in @p graphcompbegins start.
*/
static
SCIP_RETCODE buildSubgroupGraph(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  genorder,           /**< order in which the generators should be considered */
   int                   ntwocycleperms,     /**< number of 2-cycle permutations in this component */
   int                   compidx,            /**< index of the component */
   int**                 graphcomponents,    /**< buffer to store the components of the graph (ordered var indices) */
   int**                 graphcompbegins,    /**< buffer to store the indices of each new graph component */
   int**                 compcolorbegins,    /**< buffer to store at which indices a new color begins */
   int*                  ngraphcomponents,   /**< pointer to store the number of graph components */
   int*                  ncompcolors,        /**< pointer to store the number of different colors */
   int**                 usedperms,          /**< buffer to store the indices of permutations that were used */
   int*                  nusedperms,         /**< pointer to store the number of used permutations in the graph */
   int                   usedpermssize,      /**< initial size of usedperms */
   SCIP_Shortbool*       permused            /**< initialized buffer to store which permutations have been used
                                              *   (identified by index in component) */
   )
{
   SCIP_DISJOINTSET* vartocomponent;
   SCIP_DISJOINTSET* comptocolor;
   int** perms;
   int* components;
   int* componentbegins;
   int* componentslastperm;
   SYM_SORTGRAPHCOMPVARS graphcompvartype;
   int npermvars;
   int nextcolor;
   int nextcomp;
   int j;
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( graphcomponents != NULL );
   assert( graphcompbegins != NULL );
   assert( compcolorbegins != NULL );
   assert( ngraphcomponents != NULL );
   assert( ncompcolors != NULL );
   assert( genorder != NULL );
   assert( usedperms != NULL );
   assert( nusedperms != NULL );
   assert( usedpermssize > 0 );
   assert( permused != NULL );
   assert( ntwocycleperms >= 0 );
   assert( compidx >= 0 );
   assert( compidx < propdata->ncomponents );
   assert( propdata->computedsymmetry );
   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->ncomponents > 0 );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );
   assert( ! propdata->componentblocked[compidx] );

   perms = propdata->perms;
   npermvars = propdata->npermvars;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   *nusedperms = 0;

   assert( ntwocycleperms <= componentbegins[compidx + 1] - componentbegins[compidx] );

   SCIP_CALL( SCIPcreateDisjointset(scip, &vartocomponent, npermvars) );
   SCIP_CALL( SCIPcreateDisjointset(scip, &comptocolor, npermvars) );
   SCIP_CALL( SCIPallocBufferArray( scip, &componentslastperm, npermvars) );

   for (k = 0; k < npermvars; ++k)
      componentslastperm[k] = -1;

   for (j = 0; j < ntwocycleperms; ++j)
   {
      int* perm;
      int firstcolor = -1;

      /* use given order of generators */
      perm = perms[components[componentbegins[compidx] + genorder[j]]];
      assert( perm != NULL );

      /* iteratively handle each swap of perm until an invalid one is found or all edges have been added */
      for (k = 0; k < npermvars; ++k)
      {
         int comp1;
         int comp2;
         int color1;
         int color2;
         int img;

         img = perm[k];
         assert( perm[img] == k );

         if ( img <= k )
            continue;

         comp1 = SCIPdisjointsetFind(vartocomponent, k);
         comp2 = SCIPdisjointsetFind(vartocomponent, img);

         if ( comp1 == comp2 )
         {
            /* another permutation has already merged these variables into one component; store its color */
            if ( firstcolor < 0 )
            {
               assert( SCIPdisjointsetFind(comptocolor, comp1) == SCIPdisjointsetFind(comptocolor, comp2) );
               firstcolor = SCIPdisjointsetFind(comptocolor, comp1);
            }
            componentslastperm[comp1] = j;
            continue;
         }

         /* if it is the second time that the component is used for this generator,
          * it is not guaranteed that the group acts like the symmetric group, so skip it
          */
         if ( componentslastperm[comp1] == j || componentslastperm[comp2] == j )
            break;

         color1 = SCIPdisjointsetFind(comptocolor, comp1);
         color2 = SCIPdisjointsetFind(comptocolor, comp2);

         /* a generator is not allowed to connect two components of the same color, since they depend on each other */
         if ( color1 == color2 )
            break;

         componentslastperm[comp1] = j;
         componentslastperm[comp2] = j;

         if ( firstcolor < 0 )
            firstcolor = color1;
      }

      /* if the generator is invalid, delete the newly added edges, go to next generator */
      if ( k < npermvars )
         continue;

      /* if the generator only acts on already existing components, we don't have to store it */
      if ( firstcolor == -1 )
         continue;

      /* check whether we need to resize */
      if ( *nusedperms >= usedpermssize )
      {
         int newsize = SCIPcalcMemGrowSize(scip, (*nusedperms) + 1);
         assert( newsize > usedpermssize );

         SCIP_CALL( SCIPreallocBufferArray(scip, usedperms, newsize) );

         usedpermssize = newsize;
      }

      (*usedperms)[*nusedperms] = components[componentbegins[compidx] + genorder[j]];
      ++(*nusedperms);
      permused[genorder[j]] = TRUE;

      /* if the generator can be added, update the datastructures for graph components and colors */
      for (k = 0; k < npermvars; ++k)
      {
         int comp1;
         int comp2;
         int color1;
         int color2;
         int img;

         img = perm[k];
         assert( perm[img] == k );

         if ( img <= k )
            continue;

         comp1 = SCIPdisjointsetFind(vartocomponent, k);
         comp2 = SCIPdisjointsetFind(vartocomponent, img);

         /* components and colors don't have to be updated if the components are the same */
         if ( comp1 == comp2 )
            continue;

         color1 = SCIPdisjointsetFind(comptocolor, comp1);
         color2 = SCIPdisjointsetFind(comptocolor, comp2);

         if ( color1 != color2 )
         {
            SCIPdisjointsetUnion(comptocolor, firstcolor, color1, TRUE);
            SCIPdisjointsetUnion(comptocolor, firstcolor, color2, TRUE);
         }

         SCIPdisjointsetUnion(vartocomponent, comp1, comp2, FALSE);

         assert( SCIPdisjointsetFind(vartocomponent, k) == SCIPdisjointsetFind(vartocomponent, img) );
         assert( SCIPdisjointsetFind(comptocolor, SCIPdisjointsetFind(vartocomponent, k)) == firstcolor );
         assert( SCIPdisjointsetFind(comptocolor, SCIPdisjointsetFind(vartocomponent, img)) == firstcolor );
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, graphcomponents, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(graphcompvartype.components), npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(graphcompvartype.colors), npermvars) );

   /*
    * At this point, we have built the colored graph. Now we transform the information in the
    * disjoint sets to the arrays graphcomponents, graphcompbegins, and compcolorbegins (see above).
    */

   /* build the struct graphcompvartype which is used to sort the graphcomponents array */
   for (j = 0; j < npermvars; ++j)
   {
      int comp;

      comp = SCIPdisjointsetFind(vartocomponent, j);

      graphcompvartype.components[j] = comp;
      graphcompvartype.colors[j] = SCIPdisjointsetFind(comptocolor, comp);

      (*graphcomponents)[j] = j;
   }

   /* sort graphcomponents first by color, then by component */
   SCIPsort(*graphcomponents, SYMsortGraphCompVars, (void*) &graphcompvartype, npermvars);

   *ngraphcomponents = SCIPdisjointsetGetComponentCount(vartocomponent);
   *ncompcolors = SCIPdisjointsetGetComponentCount(comptocolor);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, graphcompbegins, (*ngraphcomponents) + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, compcolorbegins, (*ncompcolors) + 1) );

   nextcolor = 1;
   nextcomp = 1;
   (*graphcompbegins)[0] = 0;
   (*compcolorbegins)[0] = 0;

   /* find the starting indices of new components and new colors */
   for (j = 1; j < npermvars; ++j)
   {
      int idx1;
      int idx2;

      idx1 = (*graphcomponents)[j];
      idx2 = (*graphcomponents)[j-1];

      assert( graphcompvartype.colors[idx1] >= graphcompvartype.colors[idx2] );

      if ( graphcompvartype.components[idx1] != graphcompvartype.components[idx2] )
      {
         (*graphcompbegins)[nextcomp] = j;

         if ( graphcompvartype.colors[idx1] > graphcompvartype.colors[idx2] )
         {
            (*compcolorbegins)[nextcolor] = nextcomp;
            ++nextcolor;
         }

         ++nextcomp;
      }
   }
   assert( nextcomp == *ngraphcomponents );
   assert( nextcolor == *ncompcolors );

   (*compcolorbegins)[nextcolor] = *ngraphcomponents;
   (*graphcompbegins)[nextcomp] = npermvars;

   SCIPfreeBufferArray(scip, &(graphcompvartype.colors));
   SCIPfreeBufferArray(scip, &(graphcompvartype.components));
   SCIPfreeBufferArray(scip, &componentslastperm);
   SCIPfreeDisjointset(scip, &comptocolor);
   SCIPfreeDisjointset(scip, &vartocomponent);

   return SCIP_OKAY;
}

/** adds an orbitope constraint for a suitable color of the subgroup graph */
static
SCIP_RETCODE addOrbitopeSubgroup(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  usedperms,          /**< array of the permutations that build the orbitope */
   int                   nusedperms,         /**< number of permutations in usedperms */
   int*                  compcolorbegins,    /**< array indicating where a new graphcolor begins */
   int*                  graphcompbegins,    /**< array indicating where a new graphcomponent begins */
   int*                  graphcomponents,    /**< array of all variable indices sorted by color and comp */
   int                   graphcoloridx,      /**< index of the graph color */
   int                   nrows,              /**< number of rows in the orbitope  */
   int                   ncols,              /**< number of columns in the orbitope  */
   int*                  firstvaridx,        /**< buffer to store the index of the largest variable (or NULL) */
   int*                  compidxfirstrow,    /**< buffer to store the comp index for the first row (or NULL) */
   int**                 lexorder,           /**< pointer to array storing lexicographic order defined by sub orbitopes */
   int*                  nvarslexorder,      /**< number of variables in lexicographic order */
   int*                  maxnvarslexorder,   /**< maximum number of variables in lexicographic order */
   SCIP_Bool             mayinteract,        /**< whether orbitope's symmetries might interact with other symmetries */
   SCIP_Bool*            success             /**< whether the orbitope could be added */
   )
{  /*lint --e{571}*/
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR*** orbitopevarmatrix;
   SCIP_Shortbool* activevars;
   int** orbitopevaridx;
   int* columnorder;
   int* nusedelems;
   SCIP_CONS* cons;
   SCIP_Bool isorbitope;
   SCIP_Bool infeasible = FALSE;
#ifndef NDEBUG
   int nactivevars = 0;
#endif
   int ngencols = 0;
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( usedperms != NULL );
   assert( compcolorbegins != NULL );
   assert( graphcompbegins != NULL );
   assert( graphcomponents != NULL );
   assert( nusedperms > 0 );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( lexorder != NULL );
   assert( nvarslexorder != NULL );
   assert( maxnvarslexorder != NULL );

   *success = FALSE;

   /* create hashset to mark variables */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &activevars, propdata->npermvars) );

   /* orbitope matrix for indices of variables in permvars array */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx, nrows) );
   for (k = 0; k < nrows; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx[k], ncols) ); /*lint !e866*/
   }

   /* order of columns of orbitopevaridx */
   SCIP_CALL( SCIPallocBufferArray(scip, &columnorder, ncols) );
   for (k = 0; k < ncols; ++k)
      columnorder[k] = ncols + 1;

   /* count how often an element was used in the potential orbitope */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nusedelems, propdata->npermvars) );

   /* mark variables in this subgroup orbitope */
   for (k = compcolorbegins[graphcoloridx]; k < compcolorbegins[graphcoloridx+1]; ++k)
   {
      SCIP_VAR* firstvar;
      int compstart;
      int l;

      compstart = graphcompbegins[k];
      firstvar = propdata->permvars[graphcomponents[compstart]];

      if ( ! SCIPvarIsBinary(firstvar) )
         continue;

      for (l = 0; l < ncols; ++l)
      {
         int varidx;

         varidx = graphcomponents[compstart + l];
         assert( ! activevars[varidx] );

         activevars[varidx] = TRUE;
#ifndef NDEBUG
         ++nactivevars;
#endif
      }
   }
   assert( nactivevars == nrows * ncols );

   /* build the variable index matrix for the orbitope
    *
    * It is possible that we find an orbitope, but not using all possible columns. For example
    * (1,2), (2,3), (3,4), (3,5) defines the symmetric group on {1,2,3,4,5}, but the generators
    * we expect in our construction need shape (1,2), (2,3), (3,4), (4,5). For this reason,
    * we need to store how many columns have been generated.
    *
    * @todo ensure compatibility with more general generators
    */
   SCIP_CALL( checkTwoCyclePermsAreOrbitope(scip, propdata->permvars, propdata->npermvars,
         propdata->perms, usedperms, nrows, nusedperms, orbitopevaridx, columnorder,
         nusedelems, &ngencols, NULL, &isorbitope, activevars) );

   /* it might happen that we cannot detect the orbitope if it is generated by permutations with different
    *  number of 2-cycles.
    */
   if ( ! isorbitope )
   {
      SCIPfreeBufferArray(scip, &nusedelems);
      SCIPfreeBufferArray(scip, &columnorder);
      for (k = nrows - 1; k >= 0; --k)
      {
         SCIPfreeBufferArray(scip, &orbitopevaridx[k]);
      }
      SCIPfreeBufferArray(scip, &orbitopevaridx);
      SCIPfreeBufferArray(scip, &activevars);

      return SCIP_OKAY;
   }

   /* There are three possibilities for the structure of columnorder:
    * 1)  [0, 1, -1, -1, ..., -1]
    * 2)  [0, 1, 1, 1, ..., 1]
    * 3)  [0, 1, -1, -1, ...., -1, 1, 1, ..., 1]
    *
    * The '1'-columns will be added to the matrix first and in the last 2
    * cases the method starts from the right. So to store the variable index
    * that will be in the upper-left corner, we need either the entryin the
    * second column (case 1) or the entry in the last column (cases 2 and 3).
    */
   if ( firstvaridx != NULL )
   {
      if ( columnorder[ngencols-1] > -1 )
         *firstvaridx = orbitopevaridx[0][ngencols-1];
      else
         *firstvaridx = orbitopevaridx[0][1];
   }

   /* find corresponding graphcomponent of first variable (needed for weak sbcs) */
   if ( compidxfirstrow != NULL && firstvaridx != NULL )
   {
      *compidxfirstrow = -1;

      for (k = compcolorbegins[graphcoloridx]; k < compcolorbegins[graphcoloridx+1] && (*compidxfirstrow) < 0; ++k)
      {
         SCIP_VAR* firstvar;
         int compstart;
         int l;

         compstart = graphcompbegins[k];
         firstvar = propdata->permvars[graphcomponents[compstart]];

         if ( ! SCIPvarIsBinary(firstvar) )
            continue;

         /* iterate over all columns (elements in orbit), because we cannot see from ngencols which columns
          * have been left out
          */
         for (l = 0; l < ncols; ++l)
         {
            if ( graphcomponents[compstart + l] == *firstvaridx )
            {
               *compidxfirstrow = k;
               break;
            }
         }
      }
      assert( *compidxfirstrow > -1 );
   }

   /* prepare orbitope variable matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix, nrows) );
   for (k = 0; k < nrows; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevarmatrix[k], ngencols) );
   }

   /* build the matrix containing the actual variables of the orbitope */
   SCIP_CALL( SCIPgenerateOrbitopeVarsMatrix(scip, &orbitopevarmatrix, nrows, ngencols,
         propdata->permvars, propdata->npermvars, orbitopevaridx, columnorder,
         nusedelems, NULL, &infeasible, TRUE, lexorder, nvarslexorder, maxnvarslexorder) );

   assert( ! infeasible );
   assert( firstvaridx == NULL || propdata->permvars[*firstvaridx] == orbitopevarmatrix[0][0] );

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "suborbitope_%d_%d", graphcoloridx, propdata->norbitopes);

   SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, orbitopevarmatrix,
         SCIP_ORBITOPETYPE_FULL, nrows, ngencols, FALSE, mayinteract, FALSE, FALSE, propdata->conssaddlp,
         TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   *success = TRUE;

   /* do not release constraint here - will be done later */
   SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genorbconss,
      &propdata->genorbconsssize, propdata->ngenorbconss + 1) );
   propdata->genorbconss[propdata->ngenorbconss++] = cons;
   ++propdata->norbitopes;

   for (k = nrows - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &orbitopevarmatrix[k]);
   SCIPfreeBufferArray(scip, &orbitopevarmatrix);
   SCIPfreeBufferArray(scip, &nusedelems);
   SCIPfreeBufferArray(scip, &columnorder);
   for (k = nrows - 1; k >= 0; --k)
      SCIPfreeBufferArray(scip, &orbitopevaridx[k]);
   SCIPfreeBufferArray(scip, &orbitopevaridx);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

/** adds strong SBCs for a suitable color of the subgroup graph */
static
SCIP_RETCODE addStrongSBCsSubgroup(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  graphcompbegins,    /**< array indicating where a new graphcomponent begins */
   int*                  graphcomponents,    /**< array of all variable indices sorted by color and comp */
   int                   graphcompidx,       /**< index of the graph component */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing lexicographic order defined by sub orbitopes */
   int*                  nvarsorder,         /**< number of variables in lexicographic order */
   int*                  maxnvarsorder       /**< maximum number of variables in lexicographic order */
   )
{
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( graphcompbegins != NULL );
   assert( graphcomponents != NULL );
   assert( graphcompidx >= 0 );
   assert( ! storelexorder || lexorder != NULL );
   assert( ! storelexorder || nvarsorder != NULL );
   assert( ! storelexorder || maxnvarsorder != NULL );

   /* possibly store lexicographic order defined by strong SBCs */
   if ( storelexorder )
   {
      if ( *maxnvarsorder == 0 )
      {
         *maxnvarsorder = graphcompbegins[graphcompidx + 1] - graphcompbegins[graphcompidx + 1];
         *nvarsorder = 0;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, lexorder, *maxnvarsorder) );
      }
      else
      {
         assert( *nvarsorder == *maxnvarsorder );

         *maxnvarsorder += graphcompbegins[graphcompidx + 1] - graphcompbegins[graphcompidx + 1];

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, lexorder, *nvarsorder, *maxnvarsorder) );
      }

      (*lexorder)[*nvarsorder++] = graphcomponents[graphcompbegins[graphcompidx]];
   }

   /* add strong SBCs (lex-max order) for chosen graph component */
   for (k = graphcompbegins[graphcompidx]+1; k < graphcompbegins[graphcompidx+1]; ++k)
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_CONS* cons;
      SCIP_VAR* vars[2];
      SCIP_Real vals[2] = {1, -1};

      vars[0] = propdata->permvars[graphcomponents[k-1]];
      vars[1] = propdata->permvars[graphcomponents[k]];

      if ( storelexorder )
         (*lexorder)[*nvarsorder++] = graphcomponents[k];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "strong_sbcs_%s_%s", SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]));

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, 0.0,
            SCIPinfinity(scip), propdata->conssaddlp, propdata->conssaddlp, TRUE,
            TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* check whether we need to resize */
      SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genlinconss,
         &propdata->genlinconsssize, propdata->ngenlinconss + 1) );
      propdata->genlinconss[propdata->ngenlinconss] = cons;
      ++propdata->ngenlinconss;
   }

   return SCIP_OKAY;
}

/** adds weak SBCs for a suitable color of the subgroup graph */
static
SCIP_RETCODE addWeakSBCsSubgroup(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  compcolorbegins,    /**< array indicating where a new graphcolor begins */
   int*                  graphcompbegins,    /**< array indicating where a new graphcomponent begins */
   int*                  graphcomponents,    /**< array of all variable indices sorted by color and comp */
   int                   ncompcolors,        /**< number of colors in the graph  */
   int*                  chosencomppercolor, /**< array indicating which comp was handled per color */
   int*                  firstvaridxpercolor,/**< array indicating the largest variable per color */
   int                   symgrpcompidx,      /**< index of the component of the symmetry group */
   int*                  naddedconss,        /**< buffer to store the number of added constraints */
   SCIP_Bool             storelexorder,      /**< whether the lexicographic order induced by the orbitope shall be stored */
   int**                 lexorder,           /**< pointer to array storing lexicographic order defined by sub orbitopes */
   int*                  nvarsorder,         /**< number of variables in lexicographic order */
   int*                  maxnvarsorder       /**< maximum number of variables in lexicographic order */
   )
{  /*lint --e{571}*/
   SCIP_HASHMAP* varsinlexorder;
   SCIP_Shortbool* usedvars;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2] = {1, -1};
   SCIP_Shortbool* varfound;
   int* orbit[2];
   int orbitsize[2] = {1, 1};
   int activeorb = 0;
   int chosencolor = -1;
   int j;
   int k;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( compcolorbegins != NULL );
   assert( graphcompbegins != NULL );
   assert( graphcomponents != NULL );
   assert( firstvaridxpercolor != NULL );
   assert( chosencomppercolor != NULL );
   assert( naddedconss != NULL );
   assert( symgrpcompidx >= 0 );
   assert( symgrpcompidx < propdata->ncomponents );
   assert( ! storelexorder || lexorder != NULL );
   assert( ! storelexorder || nvarsorder != NULL );
   assert( ! storelexorder || maxnvarsorder != NULL );

   *naddedconss = 0;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &usedvars, propdata->npermvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &varfound, propdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbit[0], propdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbit[1], propdata->npermvars) );

   /* Store the entries in lexorder in a hashmap, for fast lookups. */
   if ( lexorder == NULL || *lexorder == NULL )
   {
      /* Lexorder does not exist, so do not create hashmap. */
      varsinlexorder = NULL;
   }
   else
   {
      assert( *maxnvarsorder >= 0 );
      assert( *nvarsorder >= 0 );

      SCIP_CALL( SCIPhashmapCreate(&varsinlexorder, SCIPblkmem(scip), *maxnvarsorder) );

      for (k = 0; k < *nvarsorder; ++k)
      {
         /* add element from lexorder to hashmap.
          * Use insert, as duplicate entries in lexorder is not permitted. */
         assert((*lexorder)[k] >= 0);
         assert( ! SCIPhashmapExists(varsinlexorder, (void*) (size_t) (*lexorder)[k]) ); /* Use int as pointer */
         SCIP_CALL( SCIPhashmapInsertInt(varsinlexorder, (void*) (size_t) (*lexorder)[k], k) );
      }
   }

   /* We will store the newest and the largest orbit and activeorb will be used to mark at which entry of the array
    * orbit the newly computed one will be stored. */
   if ( ncompcolors > 0 )
   {
      SCIP_CALL( ensureSymmetryPermstransComputed(scip, propdata) );
   }
   for (j = 0; j < ncompcolors; ++j)
   {
      int graphcomp;
      int graphcompsize;
      int varidx;

      /* skip color for which we did not add anything */
      if ( chosencomppercolor[j] < 0 )
         continue;

      assert( firstvaridxpercolor[j] >= 0 );

      graphcomp = chosencomppercolor[j];
      graphcompsize = graphcompbegins[graphcomp+1] - graphcompbegins[graphcomp];
      varidx = firstvaridxpercolor[j];
      assert(varidx >= 0);

      /* if the first variable was already contained in another orbit or if there are no variables left anyway, skip the
       * component */
      if ( varfound[varidx] || graphcompsize == propdata->npermvars )
         continue;

      /* If varidx is in lexorder, then it must be the first entry of lexorder. */
      if ( varsinlexorder != NULL
         && SCIPhashmapExists(varsinlexorder, (void*) (size_t) varidx)
         && lexorder != NULL && *lexorder != NULL && *maxnvarsorder > 0 && *nvarsorder > 0
         && (*lexorder)[0] != varidx )
         continue;

      /* mark all variables that have been used in strong SBCs */
      for (k = graphcompbegins[graphcomp]; k < graphcompbegins[graphcomp+1]; ++k)
      {
         assert( 0 <= graphcomponents[k] && graphcomponents[k] < propdata->npermvars );

         usedvars[graphcomponents[k]] = TRUE;
      }

      SCIP_CALL( SCIPcomputeOrbitVar(scip, propdata->npermvars, propdata->perms,
            propdata->permstrans, propdata->components, propdata->componentbegins,
            usedvars, varfound, varidx, symgrpcompidx,
            orbit[activeorb], &orbitsize[activeorb]) );

      assert( orbit[activeorb][0] ==  varidx );

      if ( orbitsize[activeorb] > orbitsize[1 - activeorb] ) /*lint !e514*/
      {
         /* if the new orbit is larger then the old largest one, flip activeorb */
         activeorb = 1 - activeorb;
         chosencolor = j;
      }

      /* reset array */
      for (k = graphcompbegins[graphcomp]; k < graphcompbegins[graphcomp+1]; ++k)
         usedvars[graphcomponents[k]] = FALSE;
   }

   /* check if we have found at least one non-empty orbit */
   if ( chosencolor > -1 )
   {
      /* flip activeorb again to avoid confusion, it is then at the largest orbit */
      activeorb = 1 - activeorb;

      assert( orbit[activeorb][0] == firstvaridxpercolor[chosencolor] );
      vars[0] = propdata->permvars[orbit[activeorb][0]];

      assert( chosencolor > -1 );
      SCIPdebugMsg(scip, "    adding %d weak sbcs for enclosing orbit of color %d.\n", orbitsize[activeorb]-1, chosencolor);

      *naddedconss = orbitsize[activeorb] - 1;

      /* add weak SBCs for rest of enclosing orbit */
      for (j = 1; j < orbitsize[activeorb]; ++j)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];

         vars[1] = propdata->permvars[orbit[activeorb][j]];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "weak_sbcs_%d_%s_%s", symgrpcompidx, SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]));

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, 0.0,
               SCIPinfinity(scip), propdata->conssaddlp, propdata->conssaddlp, TRUE,
               TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         /* check whether we need to resize */
         SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genlinconss,
            &propdata->genlinconsssize, propdata->ngenlinconss + 1) );
         propdata->genlinconss[propdata->ngenlinconss] = cons;
         ++propdata->ngenlinconss;
      }

      /* possibly store lexicographic order defined by weak SBCs */
      if ( storelexorder )
      {
         int varidx;

         varidx = orbit[activeorb][0];
         assert(varidx >= 0);

         if ( *maxnvarsorder == 0 )
         {
            *maxnvarsorder = 1;
            *nvarsorder = 0;

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, lexorder, *maxnvarsorder) );
            (*lexorder)[(*nvarsorder)++] = varidx;
         }
         else
         {
            assert( *nvarsorder == *maxnvarsorder );
            assert( varsinlexorder != NULL );
            assert( lexorder != NULL );
            assert( *lexorder != NULL );

            /* the leader of the weak inequalities has to be the first element in the lexicographic order */
            if ( varidx == (*lexorder)[0] )
            {
               /* lexorder is already ok!! */
               assert( SCIPhashmapExists(varsinlexorder, (void*) (size_t) varidx) );
            }
            else
            {
               /* Then varidx must not be in the lexorder,
                * We must add it at the front of the array, and maintain the current order. */
               assert( ! SCIPhashmapExists(varsinlexorder, (void*) (size_t) varidx) );

               ++(*maxnvarsorder);
               ++(*nvarsorder);

               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, lexorder, *nvarsorder, *maxnvarsorder) );

               /* Shift array by one position to the right */
               for (k = *maxnvarsorder - 1; k >= 1; --k)
                  (*lexorder)[k] = (*lexorder)[k - 1];

               (*lexorder)[0] = varidx;
            }
         }
      }
   }
   else
      SCIPdebugMsg(scip, "  no further weak sbcs are valid\n");

   SCIPfreeBufferArray(scip, &orbit[1]);
   SCIPfreeBufferArray(scip, &orbit[0]);
   if ( varsinlexorder != NULL )
      SCIPhashmapFree(&varsinlexorder);
   SCIPfreeBufferArray(scip, &varfound);
   SCIPfreeCleanBufferArray(scip, &usedvars);

   return SCIP_OKAY;
}


/** temporarily adapt symmetry data to new variable order given by Schreier Sims */
static
SCIP_RETCODE adaptSymmetryDataSST(
   SCIP*                 scip,               /**< SCIP instance */
   int**                 origperms,          /**< permutation matrix w.r.t. original variable ordering */
   int**                 modifiedperms,      /**< memory for permutation matrix w.r.t. new variable ordering */
   int                   nperms,             /**< number of permutations */
   SCIP_VAR**            origpermvars,       /**< array of permutation vars w.r.t. original variable ordering */
   SCIP_VAR**            modifiedpermvars,   /**< memory for array of permutation vars w.r.t. new variable ordering */
   int                   npermvars,          /**< length or modifiedpermvars array */
   int*                  leaders,            /**< leaders of Schreier Sims constraints */
   int                   nleaders            /**< number of leaders */
   )
{
   int* permvaridx;
   int* posinpermvar;
   int leader;
   int curposleader;
   int varidx;
   int lidx;
   int i;
   int l;
   int p;

   assert( scip != NULL );
   assert( origperms != NULL );
   assert( modifiedperms != NULL );
   assert( nperms > 0 );
   assert( origpermvars != NULL );
   assert( modifiedpermvars != NULL );
   assert( npermvars > 0 );
   assert( leaders != NULL );
   assert( nleaders > 0 );

   /* initialize map from position in lexicographic order to index of original permvar */
   SCIP_CALL( SCIPallocBufferArray(scip, &permvaridx, npermvars) );
   for (i = 0; i < npermvars; ++i)
      permvaridx[i] = i;

   /* initialize map from permvaridx to its current position in the reordered permvars array */
   SCIP_CALL( SCIPallocBufferArray(scip, &posinpermvar, npermvars) );
   for (i = 0; i < npermvars; ++i)
      posinpermvar[i] = i;

   /* Iterate over leaders and put the l-th leader to the l-th position of the lexicographic order.
    * We do this by swapping the l-th leader with the element at position l of the current permvars array. */
   for (l = 0; l < nleaders; ++l)
   {
      leader = leaders[l];
      curposleader = posinpermvar[leader];
      varidx = permvaridx[curposleader];
      lidx = permvaridx[l];

      /* swap the permvar at position l with the l-th leader */
      permvaridx[curposleader] = lidx;
      permvaridx[l] = varidx;

      /* update the position map */
      posinpermvar[lidx] = curposleader;
      posinpermvar[leader] = l;
   }

   /* update the permvars array to new variable order */
   for (i = 0; i < npermvars; ++i)
      modifiedpermvars[i] = origpermvars[permvaridx[i]];

   /* update the permutation to the new variable order */
   for (p = 0; p < nperms; ++p)
   {
      for (i = 0; i < npermvars; ++i)
         modifiedperms[p][i] = posinpermvar[origperms[p][permvaridx[i]]];
   }

   SCIPfreeBufferArray(scip, &permvaridx);
   SCIPfreeBufferArray(scip, &posinpermvar);

   return SCIP_OKAY;
}


/* returns the number of found orbitopes with at least three columns per graph component or 0
 * if the found orbitopes do not satisfy certain criteria for being used
 */
static
int getNOrbitopesInComp(
   SCIP_VAR**            permvars,           /**< array of variables affected by symmetry */
   int*                  graphcomponents,    /**< array of graph components */
   int*                  graphcompbegins,    /**< array indicating starting position of graph components */
   int*                  compcolorbegins,    /**< array indicating starting positions of potential orbitopes */
   int                   ncompcolors,        /**< number of components encoded in compcolorbegins */
   int                   symcompsize         /**< size of symmetry component for that we detect suborbitopes */
   )
{
   SCIP_Bool oneorbitopecriterion = FALSE;
   SCIP_Bool multorbitopecriterion = FALSE;
   int norbitopes = 0;
   int j;

   assert( graphcompbegins != NULL );
   assert( compcolorbegins != NULL );
   assert( ncompcolors >= 0 );
   assert( symcompsize > 0 );

   for (j = 0; j < ncompcolors; ++j)
   {
      SCIP_VAR* firstvar;
      int largestcompsize = 0;
      int nbinrows= 0;
      int k;

      /* skip trivial components */
      if ( graphcompbegins[compcolorbegins[j+1]] - graphcompbegins[compcolorbegins[j]] < 2 )
         continue;

      /* check whether components of this color build an orbitope (with > 2 columns) */
      for (k = compcolorbegins[j]; k < compcolorbegins[j+1]; ++k)
      {
         int compsize;

         compsize = graphcompbegins[k+1] - graphcompbegins[k];

         /* the first component that we are looking at for this color */
         if ( largestcompsize < 1 )
         {
            if ( compsize < 3 )
               break;

            largestcompsize = compsize;
         }
         else if ( compsize != largestcompsize )
            break;

         firstvar = permvars[graphcomponents[graphcompbegins[k]]];

         /* count number of binary orbits (comps) */
         if ( SCIPvarIsBinary(firstvar) )
            ++nbinrows;
      }

      /* we have found an orbitope */
      if ( k == compcolorbegins[j+1] )
      {
         SCIP_Real threshold;
         int ncols;

         ++norbitopes;
         ncols = graphcompbegins[compcolorbegins[j] + 1] - graphcompbegins[compcolorbegins[j]];

         threshold = 0.7 * (SCIP_Real) symcompsize;

         /* check whether criteria for adding orbitopes are satisfied */
         if ( nbinrows <= 2 * ncols || (nbinrows <= 8 * ncols && nbinrows < 100) )
            multorbitopecriterion = TRUE;
         else if ( nbinrows <= 3 * ncols || (SCIP_Real) nbinrows * ncols >= threshold )
            oneorbitopecriterion = TRUE;
      }
   }

   if ( (norbitopes == 1 && oneorbitopecriterion) || (norbitopes >= 2 && multorbitopecriterion) )
      return norbitopes;

   return 0;
}


/** checks whether subgroups of the components are symmetric groups and adds SBCs for them */
static
SCIP_RETCODE detectAndHandleSubgroups(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int                   cidx                /**< index of component which shall be handled */
   )
{
   int* genorder;
   int p;
#ifdef SCIP_DEBUG
   int norbitopes = 0;
   int nstrongsbcs = 0;
   int nweaksbcs = 0;
#endif
   int** modifiedperms;
   SCIP_VAR** modifiedpermvars;
   int* nvarsincomponent;

   int* graphcomponents;
   int* graphcompbegins;
   int* compcolorbegins;
   int* chosencomppercolor = NULL;
   int* firstvaridxpercolor = NULL;
   int* usedperms;
   int usedpermssize;
   int ngraphcomponents;
   int ncompcolors;
   int ntwocycleperms;
   int npermsincomp;
   int nusedperms;
   int ntrivialcolors = 0;
   int j;
   int* lexorder = NULL;
   int nvarslexorder = 0;
   int maxnvarslexorder = 0;
   SCIP_Shortbool* permused;
   SCIP_Bool allpermsused = FALSE;
   SCIP_Bool handlednonbinarysymmetry = FALSE;
   int norbitopesincomp;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->computedsymmetry );
   assert( propdata->nperms >= 0 );
   assert( 0 <= cidx && cidx < propdata->ncomponents );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );

   /* exit if no symmetry is present or component is blocked */
   if ( propdata->nperms == 0 || propdata->componentblocked[cidx] )
      return SCIP_OKAY;

   /* exit if instance is too large */
   if ( SCIPgetNConss(scip) > propdata->maxnconsssubgroup )
      return SCIP_OKAY;

   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->permvars != NULL );

   /* create array for permutation order */
   SCIP_CALL( SCIPallocBufferArray(scip, &genorder, propdata->nperms) );

   /* create arrays for modified permutations in case we adapt the lexicographic order because of suborbitopes */
   SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms, propdata->nperms) );
   for (p = 0; p < propdata->nperms; ++p)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms[p], propdata->npermvars) );
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &modifiedpermvars, propdata->npermvars) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &nvarsincomponent, propdata->npermvars) );
   for (p = 0; p < propdata->npermvars; ++p)
   {
      if ( propdata->vartocomponent[p] >= 0 )
         ++nvarsincomponent[propdata->vartocomponent[p]];
   }

   SCIPdebugMsg(scip, "starting subgroup detection routine for component %d\n", cidx);

   npermsincomp = propdata->componentbegins[cidx + 1] - propdata->componentbegins[cidx];

   /* set the first npermsincomp entries of genorder; the others are not used for this component */
   for (j = 0; j < npermsincomp; ++j)
      genorder[j] = j;

   SCIP_CALL( chooseOrderOfGenerators(scip, propdata, cidx, &genorder, &ntwocycleperms) );

   assert( ntwocycleperms >= 0 );
   assert( ntwocycleperms <= npermsincomp );

   SCIPdebugMsg(scip, "component %d has %d permutations consisting of 2-cycles\n", cidx, ntwocycleperms);

#ifdef SCIP_MORE_DEBUG
   SCIP_Bool* used;
   int perm;
   int p;
   int k;

   SCIP_CALL( SCIPallocBufferArray(scip, &used, propdata->npermvars) );
   for (p = propdata->componentbegins[cidx]; p < propdata->componentbegins[cidx+1]; ++p)
   {
      perm = propdata->components[p];

      for (k = 0; k < propdata->npermvars; ++k)
         used[k] = FALSE;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "permutation %d\n", perm);

      for (k = 0; k < propdata->npermvars; ++k)
      {
         if ( used[k] )
            continue;

         j = propdata->perms[perm][k];

         if ( k == j )
            continue;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "(%s,", SCIPvarGetName(propdata->permvars[k]));
         used[k] = TRUE;
         while (j != k)
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s,", SCIPvarGetName(propdata->permvars[j]));
            used[j] = TRUE;

            j = propdata->perms[perm][j];
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ")");
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
   }

   SCIPfreeBufferArray(scip, &used);
#endif

   if ( ntwocycleperms < 2 )
   {
      SCIPdebugMsg(scip, "  --> skip\n");
      goto FREEBASICMEM;
   }

   usedpermssize = ntwocycleperms / 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &usedperms, usedpermssize) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &permused, npermsincomp) );

   SCIP_CALL( buildSubgroupGraph(scip, propdata, genorder, ntwocycleperms, cidx,
         &graphcomponents, &graphcompbegins, &compcolorbegins, &ngraphcomponents,
         &ncompcolors, &usedperms, &nusedperms, usedpermssize, permused) );

   SCIPdebugMsg(scip, "  created subgroup detection graph using %d of the permutations\n", nusedperms);

   if ( nusedperms == npermsincomp )
      allpermsused = TRUE;

   assert( graphcomponents != NULL );
   assert( graphcompbegins != NULL );
   assert( compcolorbegins != NULL );
   assert( ngraphcomponents > 0 );
   assert( ncompcolors > 0 );
   assert( nusedperms <= ntwocycleperms );
   assert( ncompcolors < propdata->npermvars );

   if ( nusedperms == 0 )
   {
      SCIPdebugMsg(scip, "  -> skipping component, since less no permutation was used\n");

      SCIPfreeBufferArray(scip, &permused);
      SCIPfreeBufferArray(scip, &usedperms);

      goto FREEBASICMEM;
   }

   SCIPdebugMsg(scip, "  number of different colors in the graph: %d\n", ncompcolors);

   if ( propdata->addstrongsbcs || propdata->addweaksbcs )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &chosencomppercolor, ncompcolors) );
      SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercolor, ncompcolors) );

      /* Initialize the arrays with -1 to encode that we have not added orbitopes/strong SBCs
       * yet. In case we do not modify this entry, no weak inequalities are added based on
       * this component.
       */
      for (j = 0; j < ncompcolors; ++j)
      {
         chosencomppercolor[j] = -1;
         firstvaridxpercolor[j] = -1;
      }
   }

   norbitopesincomp = getNOrbitopesInComp(propdata->permvars, graphcomponents, graphcompbegins, compcolorbegins,
      ncompcolors, nvarsincomponent[cidx]);

   /* if there is just one orbitope satisfying the requirements, handle the full component by symresacks */
   if ( norbitopesincomp == 1 )
   {
      int k;

      for (k = 0; k < npermsincomp; ++k)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_comp%d_perm%d", cidx, k);

         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name,
               propdata->perms[propdata->components[propdata->componentbegins[cidx] + k]],
               propdata->permvars, propdata->npermvars, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons));

         /* do not release constraint here - will be done later */
         SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genorbconss,
               &propdata->genorbconsssize, propdata->ngenorbconss + 1) );
         propdata->genorbconss[propdata->ngenorbconss++] = cons;
         ++propdata->nsymresacks;

         if ( ! propdata->componentblocked[cidx] )
         {
            propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
            ++propdata->ncompblocked;
         }

         SCIPdebugMsg(scip, "  add symresack for permutation %d of component %d\n", k, cidx);
      }

      goto FREEMEMORY;
   }

   for (j = 0; j < ncompcolors; ++j)
   {
      int nbinarycomps = 0;
      int largestcolorcomp = -1;
      int largestcompsize = 0;
      int k;
      SCIP_Bool isorbitope = TRUE;
      SCIP_Bool orbitopeadded = FALSE;
      SCIP_Bool useorbitope;
#ifdef SCIP_DEBUG
      SCIP_Bool binaffected = FALSE;
      SCIP_Bool intaffected = FALSE;
      SCIP_Bool contaffected = FALSE;
#endif

      /* skip trivial components */
      if ( graphcompbegins[compcolorbegins[j+1]] - graphcompbegins[compcolorbegins[j]] < 2 )
      {
         if( chosencomppercolor != NULL )
            chosencomppercolor[j] = -1;

         ++ntrivialcolors;
         continue;
      }

      SCIPdebugMsg(scip, "    color %d has %d components with overall %d variables\n", j, compcolorbegins[j+1] - compcolorbegins[j],
         graphcompbegins[compcolorbegins[j+1]] - graphcompbegins[compcolorbegins[j]]);

      /* check whether components of this color might build an orbitope (with > 2 columns) */
      for (k = compcolorbegins[j]; k < compcolorbegins[j+1]; ++k)
      {
         SCIP_VAR* firstvar;
         int compsize;

         compsize = graphcompbegins[k+1] - graphcompbegins[k];

         /* the first component that we are looking at for this color */
         if ( largestcompsize < 1 )
         {
            if ( compsize < 3 )
            {
               isorbitope = FALSE;
               break;
            }

            largestcompsize = compsize;
            largestcolorcomp = k;
         }
         else if ( compsize != largestcompsize )
         {
            /* variable orbits (compsize) have not the same size, cannot define orbitope */
            isorbitope = FALSE;
            break;
         }

         firstvar = propdata->permvars[graphcomponents[graphcompbegins[k]]];

         /* count number of binary orbits (comps) */
         if ( SCIPvarIsBinary(firstvar) )
            ++nbinarycomps;
      }

#ifdef SCIP_DEBUG
      for (k = compcolorbegins[j]; k < compcolorbegins[j+1]; ++k)
      {
         SCIP_VAR* firstvar;

         firstvar = propdata->permvars[graphcomponents[graphcompbegins[k]]];

         if ( SCIPvarIsBinary(firstvar) )
            binaffected = TRUE;
         else if (SCIPvarIsIntegral(firstvar) )
            intaffected = TRUE;
         else
            contaffected = TRUE;
      }

      SCIPdebugMsg(scip, "      affected types (bin,int,cont): (%d,%d,%d)\n", binaffected, intaffected, contaffected);
#endif

      /* only use the orbitope if there are binary rows */
      useorbitope = FALSE;
      if ( norbitopesincomp > 0 && nbinarycomps > 0 )
         useorbitope = TRUE;

      if ( isorbitope && useorbitope )
      {
         int firstvaridx;
         int chosencomp;

         SCIPdebugMsg(scip, "      detected an orbitope with %d rows and %d columns\n", nbinarycomps, largestcompsize);

         assert( nbinarycomps > 0 );
         assert( largestcompsize > 2 );

         /* add the orbitope constraint for this color
          *
          * It might happen that we cannot generate the orbitope matrix if the orbitope is not generated by permutations
          * all having the same number of 2-cycles, e.g., the orbitope generated by (1,2)(4,5), (2,3), (5,6).
          */
         SCIP_CALL( addOrbitopeSubgroup(scip, propdata, usedperms, nusedperms, compcolorbegins,
               graphcompbegins, graphcomponents, j, nbinarycomps, largestcompsize, &firstvaridx, &chosencomp,
               &lexorder, &nvarslexorder, &maxnvarslexorder, allpermsused, &orbitopeadded) );

         if ( orbitopeadded )
         {
            if ( propdata->addstrongsbcs || propdata->addweaksbcs )
            {
               assert( chosencomppercolor != NULL );
               assert( firstvaridxpercolor != NULL );

               /* adapt the first variable per color to be compatible with the created orbiope (upper left variable) */
               assert( compcolorbegins[j] <= chosencomp && chosencomp < compcolorbegins[j+1] );
               assert( 0 <= firstvaridx && firstvaridx < propdata->npermvars );

               chosencomppercolor[j] = chosencomp;
               firstvaridxpercolor[j] = firstvaridx;
            }

            if ( ! propdata->componentblocked[cidx] )
            {
               propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
               ++propdata->ncompblocked;
            }

#ifdef SCIP_DEBUG
            ++norbitopes;
#endif
         }
      }

      /* if no (useable) orbitope was found, possibly add strong SBCs */
      if ( propdata->addstrongsbcs && ! orbitopeadded )
      {
         assert( largestcolorcomp >= 0 );
         assert( largestcolorcomp < ngraphcomponents );
         assert( largestcompsize > 0 );

         if( propdata->addweaksbcs )
         {
            assert( chosencomppercolor != NULL );
            assert( firstvaridxpercolor != NULL );

            chosencomppercolor[j] = largestcolorcomp;
            firstvaridxpercolor[j] = graphcomponents[graphcompbegins[largestcolorcomp]];
         }

         SCIPdebugMsg(scip, "      choosing component %d with %d variables and adding strong SBCs\n",
            largestcolorcomp, graphcompbegins[largestcolorcomp+1] - graphcompbegins[largestcolorcomp]);

         /* add the strong SBCs for the corresponding component */
         SCIP_CALL( addStrongSBCsSubgroup(scip, propdata, graphcompbegins, graphcomponents, largestcolorcomp,
               propdata->addsymresacks, &lexorder, &nvarslexorder, &maxnvarslexorder) );

         /* store whether symmetries on non-binary symmetries have been handled */
         if ( ! SCIPvarIsBinary(propdata->permvars[graphcomponents[graphcompbegins[largestcolorcomp]]]) )
            handlednonbinarysymmetry = TRUE;

         if ( ! propdata->componentblocked[cidx] )
         {
            propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
            ++propdata->ncompblocked;
         }

#ifdef SCIP_DEBUG
         nstrongsbcs += graphcompbegins[largestcolorcomp+1] - graphcompbegins[largestcolorcomp] - 1;
#endif
      }
      else if ( ! orbitopeadded )
      {
         SCIPdebugMsg(scip, "      no useable orbitope found and no SBCs added\n");

         /* mark the color as not handled */
         if ( propdata->addweaksbcs )
         {
            assert( chosencomppercolor != NULL );
            chosencomppercolor[j] = -1; /*lint !e613*/
         }
      }
   }

   SCIPdebugMsg(scip, "    skipped %d trivial colors\n", ntrivialcolors);

   /* possibly add weak SBCs for enclosing orbit of first component */
   if ( propdata->addweaksbcs && propdata->componentblocked[cidx] && nusedperms < npermsincomp )
   {
      int naddedconss;

      assert( firstvaridxpercolor != NULL );
      assert( chosencomppercolor != NULL );

      SCIP_CALL( addWeakSBCsSubgroup(scip, propdata, compcolorbegins, graphcompbegins,
            graphcomponents, ncompcolors, chosencomppercolor, firstvaridxpercolor,
            cidx, &naddedconss, propdata->addsymresacks, &lexorder, &nvarslexorder, &maxnvarslexorder) );

      assert( naddedconss < propdata->npermvars );

#ifdef SCIP_DEBUG
      nweaksbcs += naddedconss;
#endif
   }
   else
      SCIPdebugMsg(scip, "  don't add weak sbcs because all generators were used or the settings forbid it\n");

   /* if suborbitopes or strong group actions have been found, potentially add symresacks adapted to
    * variable order given by lexorder if no symmetries on non-binary variables have been handled
    */
   if ( nvarslexorder > 0 && propdata->addsymresacks && ! handlednonbinarysymmetry )
   {
      int k;

      SCIP_CALL( adaptSymmetryDataSST(scip, propdata->perms, modifiedperms, propdata->nperms,
            propdata->permvars, modifiedpermvars, propdata->npermvars, lexorder, nvarslexorder) );

      for (k = 0; k < npermsincomp; ++k)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];
         int* symresackperm;
         SCIP_Bool actsonbinary = FALSE;

         /* skip permutations that have been used to build an orbitope */
         if ( permused[k] )
            continue;

         /* skip permutations that do not act on binary variables */
         symresackperm = modifiedperms[propdata->components[propdata->componentbegins[cidx] + k]];
         for (j = 0; j < propdata->nperms && !actsonbinary; ++j)
         {
            if ( symresackperm[j] != j && SCIPvarIsBinary(modifiedpermvars[j]) )
               actsonbinary = TRUE;
         }

         if ( ! actsonbinary )
            continue;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_comp%d_perm%d", cidx, k);

         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name,
               symresackperm, modifiedpermvars, propdata->npermvars, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons));

         /* do not release constraint here - will be done later */
         SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genorbconss,
               &propdata->genorbconsssize, propdata->ngenorbconss + 1) );
         propdata->genorbconss[propdata->ngenorbconss++] = cons;
         ++propdata->nsymresacks;

         if ( ! propdata->componentblocked[cidx] )
         {
            propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
            ++propdata->ncompblocked;
         }

         SCIPdebugMsg(scip, "  add symresack for permutation %d of component %d adapted to suborbitope lexorder\n", k, cidx);
      }
   }

 FREEMEMORY:
   SCIPfreeBlockMemoryArrayNull(scip, &lexorder, maxnvarslexorder);

   SCIPfreeBufferArrayNull(scip, &firstvaridxpercolor);
   SCIPfreeBufferArrayNull(scip, &chosencomppercolor);
   SCIPfreeBlockMemoryArrayNull(scip, &compcolorbegins, ncompcolors + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graphcompbegins, ngraphcomponents + 1);
   SCIPfreeBlockMemoryArrayNull(scip, &graphcomponents, propdata->npermvars);
   SCIPfreeBufferArrayNull(scip, &permused);
   SCIPfreeBufferArrayNull(scip, &usedperms);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "total number of added (sub-)orbitopes: %d\n", norbitopes);
   SCIPdebugMsg(scip, "total number of added strong sbcs: %d\n", nstrongsbcs);
   SCIPdebugMsg(scip, "total number of added weak sbcs: %d\n", nweaksbcs);
#endif

 FREEBASICMEM:
   SCIPfreeBufferArray(scip, &nvarsincomponent);

   SCIPfreeBufferArray(scip, &modifiedpermvars);
   for (p = propdata->nperms - 1; p >= 0; --p)
   {
      SCIPfreeBufferArray(scip, &modifiedperms[p]);
   }
   SCIPfreeBufferArray(scip, &modifiedperms);
   SCIPfreeBufferArray(scip, &genorder);

   return SCIP_OKAY;
}


/*
 * Functions for symmetry constraints
 */


/** update symmetry information of conflict graph */
static
SCIP_RETCODE updateSymInfoConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONFLICTDATA*    varconflicts,       /**< conflict structure */
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

   assert( scip != NULL );
   assert( varconflicts != NULL );
   assert( conflictvars != NULL );
   assert( nconflictvars > 0 );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits >= 0 );

   /* initialize/reset variable information of nodes in conflict graph */
   for (i = 0; i < nconflictvars; ++i)
   {
      /* (re-)set node data */
      varconflicts[i].orbitidx = -1;
      varconflicts[i].nconflictinorbit = 0;
      varconflicts[i].orbitsize = -1;
      varconflicts[i].posinorbit = -1;
   }

   /* add orbit information to nodes of conflict graph */
   for (r = 0; r < norbits; ++r)
   {
      int posinorbit = 0;
      int orbitsize;

      orbitsize = orbitbegins[r + 1] - orbitbegins[r];
      assert( orbitsize >= 0 );

      for (i = orbitbegins[r]; i < orbitbegins[r + 1]; ++i)
      {
         int pos;

         /* get variable and position in conflict graph */
         pos = orbits[i];
         assert( pos < nconflictvars );
         assert( varconflicts[pos].var == conflictvars[pos] );

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
      for (i = orbitbegins[r]; i < orbitbegins[r + 1]; ++i)
      {
         ii = orbits[i];
         assert( varconflicts[ii].orbitidx == r );

         /* skip inactive variables */
         if ( ! varconflicts[ii].active )
            continue;

         for (j = i + 1; j < orbitbegins[r + 1]; ++j)
         {
            jj = orbits[j];
            assert( varconflicts[jj].orbitidx == r );

            /* skip inactive variables */
            if ( ! varconflicts[jj].active )
               continue;

            /* Check if i and j are overlapping in some clique, where only one of the two could have value 1.
             * Use that cliques are sorted by the constraint address.
             *
             * @todo A better sorted order would be: First constraints with large variables (higher hitting probability)
             *  and then by a unique constraint identifier (address, or conspos).
             */
            if ( checkSortedArraysHaveOverlappingEntry((void**)varconflicts[ii].cliques,
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
   SCIP_CONFLICTDATA**   varconflicts,       /**< pointer to store the variable conflict data */
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

   assert( scip != NULL );
   assert( varconflicts != NULL );
   assert( conflictvars != NULL );
   assert( nconflictvars > 0 );

   /* we set the pointer of varconflicts to NULL to illustrate that we didn't generate it */
   *varconflicts = NULL;

   /* get cliques for creating conflict structure */

   cliques = SCIPgetCliques(scip);
   ncliques = SCIPgetNCliques(scip);
   if ( ncliques == 0 )
   {
      SCIPdebugMsg(scip, "No cliques present --> construction of conflict structure aborted.\n");
      return SCIP_OKAY;
   }

   /* construct variable conflicts */
   SCIPdebugMsg(scip, "Construction of conflict structure:\n");
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, varconflicts, nconflictvars) );
   for (i = 0; i < nconflictvars; ++i)
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
   for (c = 0; c < ncliques; ++c)
   {
      clique = cliques[c];
      assert( clique != NULL );

      cliquevars = SCIPcliqueGetVars(clique);
      ncliquevars = SCIPcliqueGetNVars(clique);
      assert( cliquevars != NULL );
      assert( ncliquevars > 0 );

      SCIPdebugMsg(scip, "\tIdentify edges for clique ID: %d; Index: %d).\n", SCIPcliqueGetId(clique),
         SCIPcliqueGetIndex(clique));

      /* for all variables, list which cliques it is part of */
      for (i = 0; i < ncliquevars; ++i)
      {
         node = SCIPhashmapGetImageInt(conflictvarmap, cliquevars[i]);

         /* skip variables not in the conflictvars array (so not in hashmap, too) */
         if ( node == INT_MAX )
            continue;
         assert( node >= 0 );
         assert( node < nconflictvars );

         assert( (*varconflicts)[node].var == cliquevars[i] );
         (*varconflicts)[node].active = TRUE;
         (*varconflicts)[node].ncliques++;
      }
   }

   /* (2.) allocate the arrays */
   for (i = 0; i < nconflictvars; ++i)
   {
      assert( (*varconflicts)[i].ncliques >= 0 );
      assert( (*varconflicts)[i].cliques == NULL );
      if ( (*varconflicts)[i].ncliques > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*varconflicts)[i].cliques, (*varconflicts)[i].ncliques) );
      }
   }

   /* (3.) fill the clique constraints */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &tmpncliques, nconflictvars) );
   for (c = 0; c < ncliques; ++c)
   {
      clique = cliques[c];
      assert( clique != NULL );

      cliquevars = SCIPcliqueGetVars(clique);
      ncliquevars = SCIPcliqueGetNVars(clique);
      assert( cliquevars != NULL );
      assert( ncliquevars > 0 );

      SCIPdebugMsg(scip, "\tAdd edges for clique ID: %d; Index: %d).\n", SCIPcliqueGetId(clique),
         SCIPcliqueGetIndex(clique));

      /* for all variables, list which cliques it is part of */
      for (i = 0; i < ncliquevars; ++i)
      {
         node = SCIPhashmapGetImageInt(conflictvarmap, cliquevars[i]);

         /* skip variables not in the conflictvars array (so not in hashmap, too) */
         if ( node == INT_MAX )
            continue;

         assert( node >= 0 );
         assert( node < nconflictvars );
         assert( (*varconflicts)[node].var == cliquevars[i] );

         /* add clique to the cliques */
         assert( tmpncliques[node] < (*varconflicts)[node].ncliques );
         assert( (*varconflicts)[node].cliques != NULL );
         (*varconflicts)[node].cliques[tmpncliques[node]++] = clique;

#ifdef SCIP_DEBUG
         varncliques++;
#endif
      }
   }

   /* sort the variable cliques by the address, so checkSortedArraysHaveOverlappingEntry can detect intersections */
   for (i = 0; i < nconflictvars; ++i)
   {
      SCIPsortPtr((void**)(*varconflicts)[i].cliques, sortByPointerValue, (*varconflicts)[i].ncliques);
   }

#ifndef NDEBUG
   for (i = 0; i < nconflictvars; ++i)
   {
      assert( tmpncliques[i] == (*varconflicts)[i].ncliques );
   }
#endif

   SCIPfreeBufferArray(scip, &tmpncliques);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Construction of conflict graph terminated; %d variable-clique combinations detected.\n",
      varncliques);
#endif

   return SCIP_OKAY;
}

/** frees conflict graph */
static
SCIP_RETCODE freeConflictGraphSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONFLICTDATA**   varconflicts,       /**< conflict graph */
   int                   nvars               /**< number of nodes in conflict graph */
   )
{
   int i;
   int n;

   assert( scip != NULL );
   assert( varconflicts != NULL );
   assert( *varconflicts != NULL );
   assert( nvars >= 0 );

   for (i = nvars - 1; i >= 0; --i)
   {
      n = (*varconflicts)[i].ncliques;
      SCIPfreeBlockMemoryArray(scip, &(*varconflicts)[i].cliques, n);
   }
   SCIPfreeBlockMemoryArray(scip, varconflicts, nvars);

   return SCIP_OKAY;
}


/** adds symresack constraints */
static
SCIP_RETCODE addSymresackConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   int                   cidx                /**< index of component to be handled */
   )
{ /*lint --e{641}*/
   int* components;
   int* componentbegins;
   SCIP_VAR** permvars;
   SCIP_Bool conssaddlp;
   int** modifiedperms = NULL;
   SCIP_VAR** modifiedpermvars = NULL;
   int** perms;
   int nsymresackcons = 0;
   int npermvars;
   int nperms;
   int i;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->npermvars >= 0 );
   assert( propdata->nbinpermvars >= 0 );

   /* if no symmetries on binary variables are present */
   if ( propdata->nbinpermvars == 0 )
   {
      assert( propdata->binvaraffected == 0 );
      return SCIP_OKAY;
   }

   perms = propdata->perms;
   nperms = propdata->nperms;
   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   conssaddlp = propdata->conssaddlp;
   components = propdata->components;
   componentbegins = propdata->componentbegins;

   assert( nperms <= 0 || perms != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( 0 <= cidx && cidx < propdata->ncomponents );

   /* exit if component is already blocked by incompatible methods */
   if ( propdata->componentblocked[cidx] & (~SYM_HANDLETYPE_SST) )
      return SCIP_OKAY;
   if ( (propdata->componentblocked[cidx] & SYM_HANDLETYPE_SST) )
   {
      /* the leader must be binary for compatability */
      if ( (ISSSTINTACTIVE(propdata->sstleadervartype)
            || ISSSTIMPLINTACTIVE(propdata->sstleadervartype)
            || ISSSTCONTACTIVE(propdata->sstleadervartype)) )
         return SCIP_OKAY;
   }

   /* skip component if it has signed permutations */
   if ( propdata->componenthassignedperm[cidx] )
      return SCIP_OKAY;

   /* adapt natural variable order to a variable order that is compatible with Schreier Sims constraints */
   if ( propdata->nleaders > 0 && ISSSTBINACTIVE(propdata->sstleadervartype) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms, nperms) );
      for (p = 0; p < nperms; ++p)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &modifiedperms[p], npermvars) );
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &modifiedpermvars, npermvars) );

      for (i = 0; i < npermvars; ++i)
         modifiedpermvars[i] = permvars[i];

      SCIP_CALL( adaptSymmetryDataSST(scip, perms, modifiedperms, nperms, permvars, modifiedpermvars, npermvars,
            propdata->leaders, propdata->nleaders) );
   }

   /* loop through perms in component cidx and add symresack constraints */
   for (p = componentbegins[cidx]; p < componentbegins[cidx + 1]; ++p)
   {
      SCIP_CONS* cons;
      int permidx;
      char name[SCIP_MAXSTRLEN];

      permidx = components[p];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_component%d_perm%d", cidx, permidx);

      /* adapt permutation to leader */
      if ( propdata->nleaders > 0 && ISSSTBINACTIVE(propdata->sstleadervartype) )
      {
         assert( (propdata->componentblocked[cidx] & SYM_HANDLETYPE_SST) != 0 );
         assert( modifiedperms != NULL );
         assert( modifiedpermvars != NULL );

         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, modifiedperms[permidx], modifiedpermvars, npermvars, FALSE,
               conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[permidx], permvars, npermvars, FALSE,
               conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
      propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
      SCIP_CALL( SCIPaddCons(scip, cons) );

      /* do not release constraint here - will be done later */
      SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genorbconss,
            &propdata->genorbconsssize, propdata->ngenorbconss + 1) );
      propdata->genorbconss[propdata->ngenorbconss++] = cons;
      ++propdata->nsymresacks;
      ++nsymresackcons;
   }

   if ( propdata->nleaders > 0 && ISSSTBINACTIVE(propdata->sstleadervartype) )
   {
      assert( modifiedperms != NULL );
      assert( modifiedpermvars != NULL );

      SCIPfreeBufferArray(scip, &modifiedpermvars);
      for (p = nperms - 1; p >= 0; --p)
      {
         SCIPfreeBufferArray(scip, &modifiedperms[p]);
      }
      SCIPfreeBufferArray(scip, &modifiedperms);
   }

   SCIPdebugMsg(scip, "Added %d symresack constraints.\n", nsymresackcons);

   return SCIP_OKAY;
}


/** add Schreier Sims constraints for a specific orbit and update Schreier Sims table */
static
SCIP_RETCODE addSSTConssOrbitAndUpdateSST(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONFLICTDATA*    varconflicts,       /**< conflict graph or NULL if useconflictgraph == FALSE */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   SCIP_VAR**            permvars,           /**< permvars array */
   int*                  orbits,             /**< symmetry orbits */
   int*                  orbitbegins,        /**< array storing begin position for each orbit */
   int                   orbitidx,           /**< index of orbit for Schreier Sims constraints */
   int                   orbitleaderidx,     /**< index of leader variable for Schreier Sims constraints */
   SCIP_Shortbool*       orbitvarinconflict, /**< indicator whether orbitvar is in conflict with orbit leader */
   int                   norbitvarinconflict,/**< number of variables in conflict with orbit leader */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL) */
   )
{ /*lint --e{613,641}*/
   SCIP_CONS* cons;
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   int orbitsize;
   int posleader;
   int poscur;
   int ncuts = 0;
   SCIP_Bool addcuts = FALSE;
   int i;
#ifndef NDEBUG
   int j;
#endif

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( permvars != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( orbitidx >= 0 );
   assert( orbitleaderidx >= 0 );
   assert( orbitvarinconflict != NULL || varconflicts == NULL );
   assert( norbitvarinconflict >= 0 );
   assert( nchgbds != NULL );

   orbitsize = orbitbegins[orbitidx + 1] - orbitbegins[orbitidx];

   /* variables in conflict with leader are fixed and not treated by a cut; trailing -1 to not count the leader */
   if ( propdata->sstaddcuts )
      addcuts = TRUE;
   else if ( propdata->sstleaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT
      || propdata->ssttiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT )
      addcuts = propdata->addconflictcuts;

   if ( addcuts )
      ncuts = orbitsize - norbitvarinconflict - 1;

   /* (re-)allocate memory for Schreier Sims constraints and leaders */
   if ( ncuts > 0 )
   {
      if ( propdata->nsstconss == 0 )
      {
         assert( propdata->sstconss == NULL );
         assert( propdata->maxnsstconss == 0 );
         propdata->maxnsstconss = 2 * ncuts;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->sstconss), propdata->maxnsstconss) );
      }
      else if ( propdata->nsstconss + ncuts > propdata->maxnsstconss )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, propdata->maxnsstconss + 2 * ncuts);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(propdata->sstconss),
               propdata->maxnsstconss, newsize) );
         propdata->maxnsstconss = newsize;
      }
   }

   if ( propdata->nleaders == 0 )
   {
      propdata->maxnleaders = MIN(propdata->nperms, propdata->npermvars);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->leaders), propdata->maxnleaders) );
   }
   assert( propdata->nleaders < propdata->maxnleaders );

   /* add Schreier Sims constraints vars[0] >= vars[1], where vars[0] is always the leader */
   posleader = orbitbegins[orbitidx] + orbitleaderidx;
   vars[0] = permvars[orbits[posleader]];
   vals[0] = -1.0;
   vals[1] = 1.0;
   propdata->leaders[propdata->nleaders++] = orbits[posleader];
   *nchgbds = 0;
   for (i = 0, poscur = orbitbegins[orbitidx]; i < orbitsize; ++i, ++poscur)
   {
      if ( i == orbitleaderidx )
      {
         assert( orbitvarinconflict == NULL || ! orbitvarinconflict[i] );
         continue;
      }

      vars[1] = permvars[orbits[poscur]];
#ifndef NDEBUG
      for (j = 0; j < propdata->nleaders - 1; ++j)
      {
         assert( propdata->leaders[j] != orbits[poscur] );
      }
#endif

      /* if the i-th variable in the orbit is in a conflict with the leader, fix it to 0 */
      if ( varconflicts != NULL )
      {
         if ( orbitvarinconflict[i] )
         {
            assert( SCIPvarIsBinary(vars[1]) );
            assert( SCIPvarGetLbLocal(vars[1]) < 0.5 );
            assert( varconflicts != NULL );

            /* if variable is fixed */
            if ( SCIPvarGetUbLocal(vars[1]) > 0.5 )
            {
               SCIP_CALL( SCIPchgVarUb(scip, vars[1], 0.0) );
               ++(*nchgbds);

               /* deactivate the fixed variable (cannot contribute to a conflict anymore) */
               assert( varconflicts[orbits[poscur]].active );
               varconflicts[orbits[poscur]].active = FALSE;
            }

            /* reset value */
            orbitvarinconflict[i] = FALSE;
         }
         else if ( addcuts )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SSTcut_%d_%d", orbits[posleader], orbits[poscur]);
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, - SCIPinfinity(scip), 0.0,
                  FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            propdata->sstconss[propdata->nsstconss++] = cons;
         }
      }
      else if ( addcuts )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SSTcut_%d_%d", orbits[posleader], orbits[poscur]);
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, - SCIPinfinity(scip), 0.0,
               FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         propdata->sstconss[propdata->nsstconss++] = cons;
      }
   }

   return SCIP_OKAY;
}


/** selection rule of next orbit/leader in orbit for Schreier Sims constraints */
static
SCIP_RETCODE selectOrbitLeaderSSTConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONFLICTDATA*    varconflicts,       /**< variable conflicts structure, or NULL if we do not use it */
   SCIP_VAR**            conflictvars,       /**< variables encoded in conflict graph */
   int                   nconflictvars,      /**< number of variables encoded in conflict graph */
   int*                  orbits,             /**< orbits of stabilizer subgroup, expressed in terms of conflictvars */
   int*                  orbitbegins,        /**< array storing the begin position of each orbit in orbits */
   int                   norbits,            /**< number of orbits */
   int                   leaderrule,         /**< rule to select leader */
   int                   tiebreakrule,       /**< tie break rule to select leader */
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

   assert( scip != NULL );
   assert( conflictvars != NULL );
   assert( nconflictvars > 0 );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( norbits > 0 );
   assert( orbitidx != NULL );
   assert( leaderidx != NULL );
   assert( orbitvarinconflict != NULL || varconflicts == NULL );
   assert( norbitvarinconflict != NULL );
   assert( success != NULL );

   *orbitidx = 0;
   *leaderidx = 0;
   *norbitvarinconflict = 0;
   *success = FALSE;

   /* terminate if leader or tiebreak rule cannot be checked */
   if ( varconflicts == NULL && (leaderrule == (int) SCIP_LEADERRULE_MAXCONFLICTSINORBIT
         || tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT) )
      return SCIP_OKAY;

   /* select the leader and its orbit */
   if ( leaderrule == (int) SCIP_LEADERRULE_FIRSTINORBIT || leaderrule == (int) SCIP_LEADERRULE_LASTINORBIT )
   {
      orbitcriterion = INT_MIN;

      /* iterate over orbits and select the first one that meets the tiebreak rule */
      for (i = 0; i < norbits; ++i)
      {
         /* skip orbits containing vars different to the leader's vartype */
         /* Conflictvars is permvars! */
         if ( SCIPvarGetType(conflictvars[orbits[orbitbegins[i]]]) != leadervartype )
            continue;

         if ( tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MINORBIT )
            curcriterion = orbitbegins[i] - orbitbegins[i + 1];
         else if ( tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXORBIT )
            curcriterion = orbitbegins[i + 1] - orbitbegins[i];
         else
         {
            assert( tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT );

            /* get first or last active variable in orbit */
            if ( leaderrule == (int) SCIP_LEADERRULE_FIRSTINORBIT )
            {
               int cnt;

               cnt = orbitbegins[i];

               do
               {
                  varidx = orbits[cnt++];
               }
               while ( SCIPvarGetProbindex(conflictvars[varidx]) == -1 && cnt < orbitbegins[i + 1]);
            }
            else
            {
               int cnt;

               cnt = orbitbegins[i + 1] - 1;

               do
               {
                  varidx = orbits[cnt--];
               }
               while ( SCIPvarGetProbindex(conflictvars[varidx]) == -1 && cnt >= orbitbegins[i]);
            }

            /* skip inactive variables */
            if ( SCIPvarGetProbindex(conflictvars[varidx]) == -1 )
               continue;

            assert( varconflicts[varidx].orbitidx == i );
            /* coverity[var_deref_op] */
            curcriterion = varconflicts[varidx].nconflictinorbit;
         }

         /* update selected orbit */
         if ( curcriterion > orbitcriterion )
         {
            orbitcriterion = curcriterion;
            *orbitidx = i;
            *success = TRUE;

            if ( leaderrule == (int) SCIP_LEADERRULE_FIRSTINORBIT )
               *leaderidx = 0;
            else
               *leaderidx = orbitbegins[i + 1] - orbitbegins[i] - 1;
         }
      }

      /* store variables in conflict with leader */
      if ( *success && varconflicts != NULL )
      {
         leader = orbits[orbitbegins[*orbitidx] + *leaderidx];
         assert( leader < nconflictvars );

         if ( tiebreakrule == (int) SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT
            && varconflicts[leader].ncliques > 0 )
         {
            /* count how many active variables in the orbit conflict with "leader"
             * This is only needed if there are possible conflicts.
             */
            int varmapid;

            orbitsize = orbitbegins[*orbitidx + 1] - orbitbegins[*orbitidx];
            assert( varconflicts != NULL );
            assert( leader >= 0 && leader < nconflictvars );

            assert( orbitvarinconflict != NULL );

            for (i = 0; i < orbitsize; ++i)
            {
               /* skip the leader */
               if ( i == *leaderidx )
                  continue;

               /* get variable index in conflict graph */
               varmapid = orbits[orbitbegins[*orbitidx] + i];

               /* only active variables */
               if ( ! varconflicts[varmapid].active )
                  continue;

               /* check if leader and var have overlap */
               if ( checkSortedArraysHaveOverlappingEntry((void**)varconflicts[leader].cliques,
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
      /* only three possible values for leaderrules, so it must be MAXCONFLICTSINORBIT
       * In this case, the code must have computed the conflict graph.
       */
      assert( leaderrule == (int) SCIP_LEADERRULE_MAXCONFLICTSINORBIT );
      assert( varconflicts != NULL );

      orbitcriterion = 0;

      /* iterate over variables and select the first one that meets the tiebreak rule */
      for (i = 0; i < nconflictvars; ++i)
      {
         /* skip vars different to the leader's vartype */
         if ( SCIPvarGetType(conflictvars[i]) != leadervartype )
            continue;

         /* skip variables not affected by symmetry */
         /* coverity[var_deref_op] */
         if ( varconflicts[i].orbitidx == -1 )
            continue;

         curcriterion = varconflicts[i].nconflictinorbit;

         if ( curcriterion > orbitcriterion )
         {
            orbitcriterion = curcriterion;
            *orbitidx = varconflicts[i].orbitidx;
            *leaderidx = varconflicts[i].posinorbit;
            *success = TRUE;
         }
      }

      /* store variables in conflict with leader */
      leader = orbits[orbitbegins[*orbitidx] + *leaderidx];
      assert( leader < nconflictvars );
      assert( norbitvarinconflict != NULL );

      if ( *success && varconflicts[leader].ncliques > 0 )
      {
         /* count how many active variables in the orbit conflict with leader */
         int varmapid;

         orbitsize = orbitbegins[*orbitidx + 1] - orbitbegins[*orbitidx];
         assert( varconflicts != NULL );
         assert( leader >= 0 && leader < nconflictvars );

         assert( orbitvarinconflict != NULL );

         for (i = 0; i < orbitsize; ++i)
         {
            /* skip the leader */
            if ( i == *leaderidx )
               continue;

            /* get variable index in conflict graph */
            varmapid = orbits[orbitbegins[*orbitidx] + i];
            /* only active variables */
            if ( ! varconflicts[varmapid].active )
               continue;

            /* check if leader and var have overlap */
            if ( checkSortedArraysHaveOverlappingEntry((void**)varconflicts[leader].cliques,
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


/** add Schreier Sims constraints to the problem */
static
SCIP_RETCODE addSSTConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   SCIP_Bool             onlywithcontvars,   /**< only handle components that contain continuous variables with SST */
   int*                  nchgbds,            /**< pointer to store number of bound changes (or NULL) */
   int                   cidx                /**< index of component which shall be handled */
   )
{ /*lint --e{641}*/
   SCIP_CONFLICTDATA* varconflicts = NULL;
   SCIP_HASHMAP* permvarmap;
   SCIP_VAR** permvars;
   int** permstrans;
   int npermvars;
   int nmovedpermvars;
   int nmovedbinpermvars;
   int nmovedintpermvars;
   int nmovedimplintpermvars;
   int nmovedcontpermvars;
   int nperms;

   int* orbits;
   int* orbitbegins;
   int norbits;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   unsigned* componentblocked;

   int orbitidx;
   int orbitleaderidx;
   SCIP_Shortbool* orbitvarinconflict = NULL;
   int norbitvarinconflict;
   SCIP_Shortbool* inactiveperms;
   int ninactiveperms;
   int posleader;
   int leaderrule;
   int tiebreakrule;
   int leadervartype;
   SCIP_VARTYPE selectedtype = SCIP_VARTYPE_CONTINUOUS;
   int nvarsselectedtype;
   SCIP_Bool conflictgraphcreated = FALSE;
   SCIP_Bool mixedcomponents;
   int norbitleadercomponent;
   int* perm;
   SCIP_VARTYPE vartype;

   int i;
   int c;
   int p;
   SCIP_Bool success = TRUE;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->computedsymmetry );

   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   nperms = propdata->nperms;
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( nperms > 0 );

   SCIP_CALL( ensureSymmetryPermvarmapComputed(scip, propdata) );
   permvarmap = propdata->permvarmap;
   assert( permvarmap != NULL );

   SCIP_CALL( ensureSymmetryPermstransComputed(scip, propdata) );
   permstrans = propdata->permstrans;
   assert( permstrans != NULL );

   components = propdata->components;
   componentbegins = propdata->componentbegins;
   componentblocked = propdata->componentblocked;
   vartocomponent = propdata->vartocomponent;
   ncomponents = propdata->ncomponents;

   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( vartocomponent != NULL );
   assert( componentblocked != NULL );
   assert( ncomponents > 0 );
   assert( 0 <= cidx && cidx < ncomponents );

   /* exit if component is blocked */
   if ( componentblocked[cidx] )
      return SCIP_OKAY;

   /* skip component if it has signed permutations */
   if ( propdata->componenthassignedperm[cidx] )
      return SCIP_OKAY;

   leaderrule = propdata->sstleaderrule;
   tiebreakrule = propdata->ssttiebreakrule;
   leadervartype = propdata->sstleadervartype;
   mixedcomponents = propdata->sstmixedcomponents;

   /* if not already computed, get number of affected vars */
   SCIP_CALL( ensureSymmetryMovedpermvarscountsComputed(scip, propdata) );
   nmovedpermvars = propdata->nmovedpermvars;
   nmovedbinpermvars = propdata->nmovedbinpermvars;
   nmovedintpermvars = propdata->nmovedintpermvars;
   nmovedimplintpermvars = propdata->nmovedimplintpermvars;
   nmovedcontpermvars = propdata->nmovedcontpermvars;
   assert( nmovedpermvars > 0 );  /* nperms > 0 implies this */

   /* determine the leader's vartype */
   nvarsselectedtype = 0;
   if ( ISSSTBINACTIVE(leadervartype) && nmovedbinpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_BINARY;
      nvarsselectedtype = nmovedbinpermvars;
   }

   if ( ISSSTINTACTIVE(leadervartype) && nmovedintpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_INTEGER;
      nvarsselectedtype = nmovedintpermvars;
   }

   if ( ISSSTIMPLINTACTIVE(leadervartype) && nmovedimplintpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_IMPLINT;
      nvarsselectedtype = nmovedimplintpermvars;
   }

   if ( ISSSTCONTACTIVE(leadervartype) && nmovedcontpermvars > nvarsselectedtype )
   {
      selectedtype = SCIP_VARTYPE_CONTINUOUS;
      nvarsselectedtype = nmovedcontpermvars;
   }

   /* terminate if no variables of a possible leader type is affected */
   if ( nvarsselectedtype == 0 )
      return SCIP_OKAY;

   /* ignore this component if no continuous variables are contained */
   if ( onlywithcontvars )
   {
      for (p = componentbegins[cidx]; p < componentbegins[cidx + 1]; ++p)
      {
         perm = propdata->perms[p];
         for (i = 0; i < propdata->npermvars; ++i)
         {
            if ( perm[i] == i )
               continue;
            vartype = SCIPvarGetType(propdata->permvars[i]);
            if ( vartype == SCIP_VARTYPE_CONTINUOUS || vartype == SCIP_VARTYPE_IMPLINT )
               goto COMPONENTOK;
         }
      }
      /* loop terminated naturally, so component does not have continuous or implicitly integer variables. */
      return SCIP_OKAY;

   COMPONENTOK:
      ;
   }

   /* @todo online create the conflict graph for the variable in the current component */
   /* possibly create conflict graph; graph is not created if no cliques are present */
   if ( selectedtype == SCIP_VARTYPE_BINARY && (leaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT
         || tiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT) )
   {
      SCIP_CALL( createConflictGraphSST(scip, &varconflicts, permvars, npermvars, permvarmap) );
      conflictgraphcreated = varconflicts != NULL;
   }

   /* allocate data structures necessary for orbit computations and conflict graph */
   SCIP_CALL( SCIPallocBufferArray(scip, &inactiveperms, nperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );

   if ( conflictgraphcreated )
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &orbitvarinconflict, npermvars) );
   }

   SCIPdebugMsg(scip, "Start selection of orbits and leaders for Schreier Sims constraints.\n");
   SCIPdebugMsg(scip, "orbitidx\tleaderidx\torbitsize\n");

   if ( nchgbds != NULL )
      *nchgbds = 0;

   /* initialize array indicating whether permutations shall not be considered for orbit permutations */
   for (c = 0; c < ncomponents; ++c)
   {
      for (p = componentbegins[c]; p < componentbegins[c + 1]; ++p)
      {
         if ( c == cidx )
            inactiveperms[components[p]] = FALSE;
         else
            inactiveperms[components[p]] = TRUE;
      }
   }
   ninactiveperms = nperms - componentbegins[cidx + 1] + componentbegins[cidx];

   /* as long as the stabilizer is non-trivial, add Schreier Sims constraints */
   norbitleadercomponent = 0;
   while ( ninactiveperms < nperms )
   {
      int nchanges = 0;

      /* compute orbits w.r.t. active perms */
      SCIP_CALL( SCIPcomputeOrbitsFilterSym(scip, npermvars, permstrans, nperms, inactiveperms,
            orbits, orbitbegins, &norbits, components, componentbegins, vartocomponent,
            componentblocked, ncomponents, nmovedpermvars) );

      /* stop if we require pure components and a component contains variables of different types */
      if ( ! mixedcomponents )
      {
         for (p = 0; p < norbits; ++p)
         {
            /* stop if the first element of an orbits has the wrong vartype */
            if ( SCIPvarGetType(permvars[orbits[orbitbegins[p]]]) != selectedtype )
            {
               success = FALSE;
               break;
            }
         }
      }

      if ( ! success )
         break;

      /* update symmetry information of conflict graph */
      if ( conflictgraphcreated )
      {
         SCIP_CALL( updateSymInfoConflictGraphSST(scip, varconflicts, permvars, npermvars, orbits, orbitbegins,
               norbits) );
      }

      /* possibly adapt the leader and tie-break rule */
      if ( leaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT && ! conflictgraphcreated )
         leaderrule = SCIP_LEADERRULE_FIRSTINORBIT;
      if ( leaderrule == SCIP_LEADERRULE_MAXCONFLICTSINORBIT && selectedtype != SCIP_VARTYPE_BINARY )
         leaderrule = SCIP_LEADERRULE_FIRSTINORBIT;
      if ( tiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT && ! conflictgraphcreated )
         tiebreakrule = SCIP_LEADERTIEBREAKRULE_MAXORBIT;
      if ( tiebreakrule == SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT && selectedtype != SCIP_VARTYPE_BINARY )
         tiebreakrule = SCIP_LEADERTIEBREAKRULE_MAXORBIT;

      /* select orbit and leader */
      SCIP_CALL( selectOrbitLeaderSSTConss(scip, varconflicts, permvars, npermvars, orbits, orbitbegins,
            norbits, propdata->sstleaderrule, propdata->ssttiebreakrule, selectedtype, &orbitidx, &orbitleaderidx,
            orbitvarinconflict, &norbitvarinconflict, &success) );

      if ( ! success )
         break;

      assert( 0 <= orbitidx && orbitidx < norbits );
      assert( 0 <= orbitleaderidx && orbitleaderidx < orbitbegins[orbitidx + 1] - orbitbegins[orbitidx] );
      SCIPdebugMsg(scip, "%d\t\t%d\t\t%d\n", orbitidx, orbitleaderidx, orbitbegins[orbitidx + 1] - orbitbegins[orbitidx]);

      /* add Schreier Sims constraints for the selected orbit and update Schreier Sims table */
      SCIP_CALL( addSSTConssOrbitAndUpdateSST(scip, varconflicts, propdata, permvars,
            orbits, orbitbegins, orbitidx, orbitleaderidx, orbitvarinconflict, norbitvarinconflict, &nchanges) );

      ++norbitleadercomponent;

      if ( nchgbds != NULL )
         *nchgbds += nchanges;

      /* deactivate permutations that move the orbit leader */
      posleader = orbits[orbitbegins[orbitidx] + orbitleaderidx];
      for (p = 0; p < nperms; ++p)
      {
         if ( inactiveperms[p] )
            continue;

         if ( permstrans[posleader][p] != posleader )
         {
            inactiveperms[p] = TRUE;
            ++ninactiveperms;
         }
      }
   }

   /* if Schreier Sims constraints have been added, store that Schreier Sims has been used for this component */
   if ( norbitleadercomponent > 0 )
      componentblocked[cidx] |= SYM_HANDLETYPE_SST;

   if ( conflictgraphcreated )
   {
      SCIPfreeBufferArray(scip, &orbitvarinconflict);
   }
   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);
   if ( varconflicts != NULL )
   {
      /* nconflictvars at construction is npermvars */
      SCIP_CALL( freeConflictGraphSST(scip, &varconflicts, npermvars) );
   }
   SCIPfreeBufferArray(scip, &inactiveperms);

   return SCIP_OKAY;
}


/** orbitopal reduction */
static
SCIP_RETCODE addOrbitopesDynamic(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< propdata */
   int                   id,                 /**< ID for orbitope constraint (needed for name) */
   int**                 varidxmatrix,       /**< matrix containing variable indices in orbitope matrix */
   int                   nrows,              /**< number of rows of orbitope */
   int                   ncols,              /**< number of columns of orbitope */
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

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->usedynamicprop );
   assert( varidxmatrix != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( success != NULL );

   *success = FALSE;

   /* add linear constraints x_1 >= x_2 >= ... >= x_ncols for single-row orbitopes */
   if ( nrows == 1 )
   {
      /* restrict to the packing and partitioning rows */
      SCIP_CONS* cons;
      SCIP_VAR* consvars[2];
      SCIP_Real conscoefs[2] = { -1.0, 1.0 };

      /* for all adjacent column pairs, add linear constraint */
      SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genlinconss,
            &propdata->genlinconsssize, propdata->ngenlinconss + ncols - 1) );
      for (i = 0; i < ncols - 1; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_1row_comp_%d_col%d", id, i);

         consvars[0] = propdata->permvars[varidxmatrix[0][i]];
         consvars[1] = propdata->permvars[varidxmatrix[0][i + 1]];

         /* enforce, but do not check */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, consvars, conscoefs, -SCIPinfinity(scip), 0.0,
               propdata->conssaddlp, propdata->conssaddlp, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         propdata->genlinconss[propdata->ngenlinconss++] = cons;
      }

      *success = TRUE;
      return SCIP_OKAY;
   }

   /* for only 2 columns, the the component can be completely handled by lexicographic reduction */
   if ( ncols == 2 && propdata->lexreddata != NULL )
   {
      int* orbisackperm;

      /* If the component is an orbitope with 2 columns, then there is 1 generator of order 2. */
      orbisackperm = propdata->perms[propdata->components[propdata->componentbegins[id]]];

      SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, propdata->lexreddata,
            propdata->permvars, propdata->npermvars, orbisackperm, (SYM_SYMTYPE) propdata->symtype,
            propdata->permvardomaincenter, TRUE, success) );
      if ( *success )
         return SCIP_OKAY;
   }

   /* create orbitope variable matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &varmatrix, nrows) );
   for (i = 0; i < nrows; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &varmatrix[i], ncols) );
      for (j = 0; j < ncols; ++j)
         varmatrix[i][j] = propdata->permvars[varidxmatrix[i][j]];
   }

   pprows = NULL;
   SCIP_CALL( SCIPisPackingPartitioningOrbitope(scip, varmatrix, nrows, ncols, &pprows, &npprows, &type) );

   /* does it have at least 3 packing-partitioning rows? */
   ispporbitope = npprows >= 3;  /* (use same magic number as cons_orbitope.c) */

   if ( ispporbitope ) /* @todo if it's a pporbitope, we do it statically right now. */
   {
      /* restrict to the packing and partitioning rows */
      SCIP_CONS* cons;
      SCIP_VAR*** ppvarsarrayonlypprows;
      int r;

      assert( pprows != NULL );

      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsarrayonlypprows, npprows) );

      r = 0;
      for (i = 0; i < nrows; ++i)
      {
         if ( pprows[i] )
         {
            assert( r < npprows );
            ppvarsarrayonlypprows[r++] = varmatrix[i];
         }
      }
      assert( r == npprows );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_pp_comp_%d", id);
      SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, ppvarsarrayonlypprows, SCIP_ORBITOPETYPE_PACKING,
            npprows, ncols, FALSE, FALSE, FALSE, FALSE, propdata->conssaddlp,
            TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      /* check whether we need to resize */
      SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genlinconss,
            &propdata->genlinconsssize, propdata->ngenlinconss + 1) );
      /* @todo we add orbitopes to the dynamically sized array `genlinconss` instead of `genorbconss` to ensure
       * compatability with the static orbitope function, which allocates this array statically
       */
      propdata->genlinconss[propdata->ngenlinconss++] = cons;
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
      for (i = 0; i < nrows; ++i)
      {
         for (j = 0; j < ncols; ++j)
            orbitopevarmatrix[pos++] = varmatrix[i][j];
      }

      /* get column ordering */
      columnordering = SCIPorbitopalReductionGetDefaultColumnOrdering(propdata->orbitopalreddata);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_full_comp_%d", id);
      SCIP_CALL( SCIPorbitopalReductionAddOrbitope(scip, propdata->orbitopalreddata,
            SCIP_ROWORDERING_BRANCHING, columnordering,
            orbitopevarmatrix, nrows, ncols, success) );
      *success = TRUE;

      SCIPfreeBufferArray(scip, &orbitopevarmatrix);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &pprows, nrows);

   for (i = nrows - 1; i >= 0; --i)
   {
      SCIPfreeBufferArray(scip, &varmatrix[i]);
   }
   SCIPfreeBufferArray(scip, &varmatrix);

   return SCIP_OKAY;
}


/** applies pp-orbitope upgrade if at least 50% of the permutations in a component correspond to pp-orbisacks */
static
SCIP_RETCODE componentPackingPartitioningOrbisackUpgrade(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< propdata */
   int**                 componentperms,     /**< permutations in the component */
   int                   componentsize,      /**< number of permutations in the component */
   SCIP_Bool             hassignedperm,      /**< whether the component has a signed permutation */
   SCIP_Bool*            success             /**< whether the packing partitioning upgrade succeeded */
   )
{
   int c;
   int i;
   int j;
   int p;
   int* perm;
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_CONS* cons;
   SCIP_CONS** setppconsssort;
   int nsetppconss;
   int nsetppcvars;
   SCIP_VAR** setppcvars;
   int nsetppcconss;
   int** pporbisackperms;
   int npporbisackperms;
   SCIP_VAR* var;
   int varid;
   SCIP_CONS*** permvarssetppcconss;
   int* npermvarssetppcconss;
   int* maxnpermvarssetppcconss;
   int maxntwocycles;
   int ntwocycles;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( componentperms != NULL );
   assert( componentsize > 0 );
   assert( success != NULL );

   /* we did not upgrade yet */
   *success = FALSE;

   /* currently, we cannot handle signed permutations */
   if ( hassignedperm )
      return SCIP_OKAY;

   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
      return SCIP_OKAY;

   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);
   if ( nsetppcconss == 0 )
      return SCIP_OKAY;

   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   assert( setppcconss != NULL );

   SCIP_CALL( ensureSymmetryPermvarmapComputed(scip, propdata) );

   /* collect non-covering constraints and sort by pointer for easy intersection finding */
   SCIP_CALL( SCIPallocBufferArray(scip, &setppconsssort, nsetppcconss) );
   nsetppconss = 0;
   for (c = 0; c < nsetppcconss; ++c)
   {
      cons = setppcconss[c];

      /* only packing or partitioning constraints, no covering types */
      if ( SCIPgetTypeSetppc(scip, cons) == SCIP_SETPPCTYPE_COVERING )
         continue;

      setppconsssort[nsetppconss++] = cons;
   }
   SCIPsortPtr((void**) setppconsssort, sortByPointerValue, nsetppcconss);

   /* For each permvar, introduce an array of setppc constraints (initially NULL) for each variable,
    * and populate it with the setppc constraints that it contains. This array follows the ordering by cons ptr address.
    */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &permvarssetppcconss, propdata->npermvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &npermvarssetppcconss, propdata->npermvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &maxnpermvarssetppcconss, propdata->npermvars) );
   for (c = 0; c < nsetppconss; ++c)
   {
      assert( c >= 0 );
      assert( c < nsetppconss );
      cons = setppconsssort[c];
      assert( cons != NULL );

      setppcvars = SCIPgetVarsSetppc(scip, cons);
      nsetppcvars = SCIPgetNVarsSetppc(scip, cons);

      for (i = 0; i < nsetppcvars; ++i)
      {
         var = setppcvars[i];
         assert( var != NULL );
         varid = SCIPhashmapGetImageInt(propdata->permvarmap, (void*) var);
         assert( varid == INT_MAX || varid < propdata->npermvars );
         assert( varid >= 0 );
         if ( varid < propdata->npermvars )
         {
            SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip,
               &(permvarssetppcconss[varid]), &maxnpermvarssetppcconss[varid], npermvarssetppcconss[varid] + 1) );
            assert( npermvarssetppcconss[varid] < maxnpermvarssetppcconss[varid] );
            permvarssetppcconss[varid][npermvarssetppcconss[varid]++] = cons;
         }
      }
   }

   /* for all permutations, test involutions on binary variables and test if they are captured by setppc conss */
   SCIP_CALL( SCIPallocBufferArray(scip, &pporbisackperms, componentsize) );
   maxntwocycles = 0;
   npporbisackperms = 0;
   for (p = 0; p < componentsize; ++p)
   {
      perm = componentperms[p];
      ntwocycles = 0;

      /* check if the binary orbits are involutions */
      for (i = 0; i < propdata->npermvars; ++i)
      {
         j = perm[i];

         /* ignore fixed points in permutation */
         if ( i == j )
            continue;
         /* only check for situations where i and j are binary variables */
         assert( SCIPvarGetType(propdata->permvars[i]) == SCIPvarGetType(propdata->permvars[j]) );
         if ( SCIPvarGetType(propdata->permvars[i]) != SCIP_VARTYPE_BINARY )
            continue;
         /* the permutation must be an involution on binary variables */
         if ( perm[j] != i )
            goto NEXTPERMITER;
         /* i and j are a two-cycle, so we find this once for i and once for j. Only handle this once for i < j. */
         if ( i > j )
            continue;
         /* disqualify permutation if i and j are not in a common set packing constraint */
         if ( !checkSortedArraysHaveOverlappingEntry((void**) permvarssetppcconss[i], npermvarssetppcconss[i],
            (void**) permvarssetppcconss[j], npermvarssetppcconss[j], sortByPointerValue) )
            goto NEXTPERMITER;
         ++ntwocycles;
      }

      /* The permutation qualifies if all binary variables are either a reflection or in a 2-cycle. There must be at
       * least one binary 2-cycle, because otherwise the permutation is the identity, or it permutes
       * nonbinary variables.
       */
      if ( ntwocycles > 0 )
      {
         pporbisackperms[npporbisackperms++] = perm;
         if ( ntwocycles > maxntwocycles )
            maxntwocycles = ntwocycles;
      }

   NEXTPERMITER:
      ;
   }

   /* if at least 50% of such permutations are packing-partitioning type, apply packing upgrade */
   if ( npporbisackperms * 2 >= componentsize )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR** ppvarsblock;
      SCIP_VAR*** ppvarsmatrix;
      SCIP_VAR** row;
      int nrows;

      assert( npporbisackperms > 0 );
      assert( maxntwocycles > 0 );

      /* instead of allocating and re-allocating multiple times, recycle the ppvars array */
      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsblock, 2 * maxntwocycles) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ppvarsmatrix, maxntwocycles) );
      for (i = 0; i < maxntwocycles; ++i)
         ppvarsmatrix[i] = &(ppvarsblock[2 * i]);

      /* for each of these perms, create the packing orbitope matrix and add constraint*/
      for (p = 0; p < npporbisackperms; ++p)
      {
         perm = pporbisackperms[p];

         /* populate ppvarsmatrix */
         nrows = 0;
         for (i = 0; i < propdata->npermvars; ++i)
         {
            j = perm[i];

            /* ignore fixed points in permutation, and only consider rows with i < j */
            if ( i >= j )
               continue;
            /* only for situations where i and j are binary variables */
            assert( SCIPvarGetType(propdata->permvars[i]) == SCIPvarGetType(propdata->permvars[j]) );
            if ( SCIPvarGetType(propdata->permvars[i]) != SCIP_VARTYPE_BINARY )
               continue;
            assert( perm[j] == i );
            assert( checkSortedArraysHaveOverlappingEntry((void**) permvarssetppcconss[i], npermvarssetppcconss[i],
               (void**) permvarssetppcconss[j], npermvarssetppcconss[j], sortByPointerValue) );

            assert( nrows < maxntwocycles );
            row = ppvarsmatrix[nrows++];
            row[0] = propdata->permvars[i];
            row[1] = propdata->permvars[j];
            assert( row[0] != row[1] );
         }
         assert( nrows > 0 );

         /* create constraint, use same parameterization as in orbitope packing partitioning checker */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_pp_upgrade_lexred%d", p);
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, ppvarsmatrix, SCIP_ORBITOPETYPE_PACKING, nrows, 2,
            FALSE, FALSE, FALSE, FALSE,
            propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genlinconss,
            &propdata->genlinconsssize, propdata->ngenlinconss + 1) );
         /* @todo we add orbitopes to the dynamically sized array `genlinconss` instead of `genorbconss` to ensure
          * compatability with the static orbitope function, which allocates this array statically
          */
         propdata->genlinconss[propdata->ngenlinconss++] = cons;
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      SCIPfreeBufferArray(scip, &ppvarsmatrix);
      SCIPfreeBufferArray(scip, &ppvarsblock);

      *success = TRUE;
   }

   /* free pp orbisack array */
   SCIPfreeBufferArray(scip, &pporbisackperms);

   /* clean the non-clean arrays */
   for (varid = 0; varid < propdata->npermvars; ++varid)
   {
      assert( (permvarssetppcconss[varid] == NULL) == (maxnpermvarssetppcconss[varid] == 0) );
      assert( npermvarssetppcconss[varid] >= 0 );
      assert( maxnpermvarssetppcconss[varid] >= 0 );
      assert( npermvarssetppcconss[varid] <= maxnpermvarssetppcconss[varid] );
      if ( npermvarssetppcconss[varid] == 0 )
         continue;
      SCIPfreeBlockMemoryArray(scip, &permvarssetppcconss[varid], maxnpermvarssetppcconss[varid]);
      permvarssetppcconss[varid] = NULL;
      npermvarssetppcconss[varid] = 0;
      maxnpermvarssetppcconss[varid] = 0;
   }
   SCIPfreeCleanBufferArray(scip, &maxnpermvarssetppcconss);
   SCIPfreeCleanBufferArray(scip, &npermvarssetppcconss);
   SCIPfreeCleanBufferArray(scip, &permvarssetppcconss);
   SCIPfreeBufferArray(scip, &setppconsssort);

   return SCIP_OKAY;
}


/** dynamic permutation lexicographic reduction */
static
SCIP_RETCODE tryAddOrbitalRedLexRed(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< propdata */
   int                   cidx                /**< index of component */
   )
{
   int componentsize;
   int** componentperms;
   int p;

   SCIP_Bool checkorbired;
   SCIP_Bool checklexred;
   SCIP_Bool success;
   SCIP_PARAM* checkpporbisack;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( ISORBITALREDUCTIONACTIVE(propdata->usesymmetry)
      || (
         ISSYMRETOPESACTIVE(propdata->usesymmetry)
         && propdata->usedynamicprop
         && propdata->addsymresacks
      ) );
   assert( propdata->nperms > 0 );
   assert( 0 <= cidx && cidx < propdata->ncomponents );
   assert( propdata->componentblocked != NULL );

   /* exit if component is already blocked */
   if ( propdata->componentblocked[cidx] )
      return SCIP_OKAY;

   /* in this function orbital reduction or dynamic lexicographic reduction propagation must be enabled */
   checkorbired = ISORBITALREDUCTIONACTIVE(propdata->usesymmetry);
   checklexred = ISSYMRETOPESACTIVE(propdata->usesymmetry) && propdata->usedynamicprop && propdata->addsymresacks;
   assert( checkorbired || checklexred );

   SCIP_CALL( ensureSymmetryMovedpermvarscountsComputed(scip, propdata) );
   assert( propdata->nmovedpermvars );

   /* collect the permutations of this component */
   componentsize = propdata->componentbegins[cidx + 1] - propdata->componentbegins[cidx];
   SCIP_CALL( SCIPallocBufferArray(scip, &componentperms, componentsize) );
   for (p = 0; p < componentsize; ++p)
      componentperms[p] = propdata->perms[propdata->components[propdata->componentbegins[cidx] + p]];

   /* check if many component permutations contain many packing partitioning orbisacks
    *
    * 1. Get the checkpporbisack param from the parameter hashset. This returns NULL if it is not initialized,
    *    likely because the orbisack constraint handler is not loaded.
    * 2. If the param is not NULL, then we only do the packing-partitioning upgrade step if its value is TRUE.
    * Packing-partitioning orbitopes are only implemented for binary orbitopes, so binary variables must be moved.
    */
   checkpporbisack = SCIPgetParam(scip, "constraints/orbisack/checkpporbisack");
   if ( ( checkpporbisack == NULL || SCIPparamGetBool(checkpporbisack) == TRUE ) && propdata->nmovedbinpermvars > 0 )
   {
      SCIP_CALL( componentPackingPartitioningOrbisackUpgrade(scip, propdata,
            componentperms, componentsize, propdata->componenthassignedperm[cidx], &success) );

      if ( success )
      {
         propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
         goto FINISHCOMPONENT;
      }
   }

   /* handle component permutations with orbital reduction */
   if ( checkorbired && !propdata->componenthassignedperm[cidx] )
   {
      SCIP_CALL( SCIPorbitalReductionAddComponent(scip, propdata->orbitalreddata,
            propdata->permvars, propdata->npermvars, componentperms, componentsize, &success) );
      if ( success )
         propdata->componentblocked[cidx] |= SYM_HANDLETYPE_ORBITALREDUCTION;
   }

   /* handle component permutations with the dynamic lexicographic reduction propagator */
   if ( checklexred )
   {
      /* handle every permutation in the component with the dynamic lexicographic reduction propagator */
      for (p = 0; p < componentsize; ++p)
      {
         assert( componentperms[p] != NULL );
         SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, propdata->lexreddata,
               propdata->permvars, propdata->npermvars, componentperms[p],
               (SYM_SYMTYPE) propdata->symtype, propdata->permvardomaincenter, TRUE, &success) );
         if ( success )
            propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
      }
   }

 FINISHCOMPONENT:
   /* if it got blocked here */
   if ( propdata->componentblocked[cidx] )
      ++propdata->ncompblocked;

   SCIPfreeBufferArray(scip, &componentperms);

   return SCIP_OKAY;
}


/** displays statistics on the used symmetry handling methods */
static
SCIP_RETCODE SCIPdisplaySymmetryStatistics(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata            /**< data of symmetry propagator */
   )
{
   int ncomponentshandled;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "dynamic symmetry handling statistics:\n");
   if ( propdata->orbitopalreddata )
   {
      SCIP_CALL( SCIPorbitopalReductionPrintStatistics(scip, propdata->orbitopalreddata) );
   }
   if ( propdata->orbitalreddata )
   {
      SCIP_CALL( SCIPorbitalReductionPrintStatistics(scip, propdata->orbitalreddata) );
   }
   if ( propdata->lexreddata )
   {
      SCIP_CALL( SCIPlexicographicReductionPrintStatistics(scip, propdata->lexreddata) );
   }
   if ( propdata->ncomponents >= 0 )
   {
      /* report the number of handled components
       *
       * Since SST is compatible with static symresacks, the propdata->ncompblocked counter is not the number of
       * handled components. Compute this statistic based on the componentblocked array.
       */
      ncomponentshandled = 0;
      for (i = 0; i < propdata->ncomponents; ++i)
      {
         if ( propdata->componentblocked[i] )
            ++ncomponentshandled;
      }
      assert( propdata->ncompblocked <= ncomponentshandled );
      assert( ncomponentshandled <= propdata->ncomponents );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "handled %d out of %d symmetry components\n",
         ncomponentshandled, propdata->ncomponents);
   }

   return SCIP_OKAY;
}

/** handles orbitope action by static or dynamic symmetry handling methods */
static
SCIP_RETCODE handleOrbitope(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   int                   id,                 /**< ID of orbitope (used for constraint name) */
   int**                 varidxmatrix,       /**< matrix containing variable indices of orbitope */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   SCIP_Bool*            success             /**< pointer to store whether orbitope could be added successfully */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );
   assert( varidxmatrix != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( success != NULL );

   *success = FALSE;

   /* dynamic propagation */
   if ( propdata->usedynamicprop )
   {
      SCIP_CALL( addOrbitopesDynamic(scip, propdata, id, varidxmatrix, nrows, ncols, success) );
   }
   /* static variant only for binary variables */
   else if ( propdata->binvaraffected )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR*** orbitopematrix;
      SCIP_CONS* cons;
      int i;
      int j;
      int nbinrows = 0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "orbitope_component_%d", id);

      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, nrows) );
      for (i = 0; i < nrows; ++i)
      {
         /* skip rows without binary variables */
         if ( ! SCIPvarIsBinary(propdata->permvars[varidxmatrix[i][0]]) )
            continue;

         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[nbinrows], ncols) );
         for (j = 0; j < ncols; ++j)
         {
            assert( SCIPvarIsBinary(propdata->permvars[varidxmatrix[i][j]]) );
            orbitopematrix[nbinrows][j] = propdata->permvars[varidxmatrix[i][j]];
         }
         ++nbinrows;
      }

      if ( nbinrows > 0 )
      {
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, orbitopematrix, SCIP_ORBITOPETYPE_FULL,
               nbinrows, ncols, propdata->usedynamicprop /* @todo disable */, FALSE, FALSE, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genorbconss,
               &propdata->genorbconsssize, propdata->ngenorbconss + 1) );
         propdata->genorbconss[propdata->ngenorbconss++] = cons;
         ++propdata->norbitopes;

         *success = TRUE;
      }

      for (i = nbinrows - 1; i >= 0; --i)
      {
         SCIPfreeBufferArray(scip, &orbitopematrix[i]);
      }
      SCIPfreeBufferArray(scip, &orbitopematrix);
   }

   return SCIP_OKAY;
}

/** handles binary double lex matrix by adding static orbitope constraints
 *
 * @todo Extend method to general variable types and dynamic variable orders.
 */
static
SCIP_RETCODE handleDoublelLexMatrix(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   int                   id,                 /**< ID of double lex matrix (used for constraint names) */
   int**                 varidxmatrix,       /**< matrix containing variable indices of double lex matrix */
   int                   nrows,              /**< number of rows of matrix */
   int                   ncols,              /**< number of columns of matrix */
   int*                  rowsbegin,          /**< array indicating where a new row block begins */
   int*                  colsbegin,          /**< array indicating where a new column block begins */
   int                   nrowblocks,         /**< number of row blocks */
   int                   ncolblocks,         /**< number of column blocks */
   SCIP_Bool*            success             /**< pointer to store whether orbitope could be added successfully */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR*** orbitopematrix;
   SCIP_CONS* cons;
   int maxdim;
   int i;
   int p;
   int j;
   int col;
   int nbinrows;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( varidxmatrix != NULL );
   assert( nrows > 0 );
   assert( ncols > 0 );
   assert( rowsbegin != NULL );
   assert( colsbegin != NULL );
   assert( nrowblocks > 0 );
   assert( ncolblocks > 0 );
   assert( success != NULL );

   /* ensure that we can store orbitope constraints in probdata */
   SCIP_CALL( ensureDynamicConsArrayAllocatedAndSufficientlyLarge(scip, &propdata->genorbconss,
         &propdata->genorbconsssize, propdata->ngenorbconss + nrowblocks + ncolblocks) );

   maxdim = MAX(nrows, ncols);
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix, maxdim) );
   for (i = 0; i < maxdim; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopematrix[i], maxdim) );
   }

   /* add orbitopes corresponding to column blocks of doublelexmatrix */
   for (p = 0; p < ncolblocks; ++p)
   {
      nbinrows = 0;
      for (i = 0; i < nrows; ++i)
      {
         /* skip rows that do not contain binary variables */
         if ( ! SCIPvarIsBinary(propdata->permvars[varidxmatrix[i][colsbegin[p]]]) )
            continue;

         for (col = 0, j = colsbegin[p]; j < colsbegin[p + 1]; ++j, ++col)
         {
            assert( SCIPvarIsBinary(propdata->permvars[varidxmatrix[i][j]]) );
            orbitopematrix[nbinrows][col] = propdata->permvars[varidxmatrix[i][j]];
         }
         ++nbinrows;
      }

      if ( nbinrows > 0 )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "doublelex_cols_%d_%d", id, p);
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, orbitopematrix, SCIP_ORBITOPETYPE_FULL,
               nrows, colsbegin[p + 1] - colsbegin[p], FALSE, FALSE, TRUE, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         propdata->genorbconss[(propdata->ngenorbconss)++] = cons;
         /* do not release constraint here - will be done later */
      }
   }

   /* add orbitopes corresponding to row blocks of doublelexmatrix */
   for (p = 0; p < nrowblocks; ++p)
   {
      nbinrows = 0;
      for (i = 0; i < ncols; ++i)
      {
         /* skip rows that do not contain binary variables */
         if ( ! SCIPvarIsBinary(propdata->permvars[varidxmatrix[rowsbegin[p]][i]]) )
            continue;

         for (col = 0, j = rowsbegin[p]; j < rowsbegin[p + 1]; ++j, ++col)
         {
            assert( SCIPvarIsBinary(propdata->permvars[varidxmatrix[j][i]]) );
            orbitopematrix[nbinrows][col] = propdata->permvars[varidxmatrix[j][i]];
         }
         ++nbinrows;
      }

      if ( nbinrows > 0 )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "doublelex_rows_%d_%d", id, p);
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, name, orbitopematrix, SCIP_ORBITOPETYPE_FULL,
               ncols, rowsbegin[p + 1] - rowsbegin[p], FALSE, FALSE, TRUE, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         propdata->genorbconss[(propdata->ngenorbconss)++] = cons;
         /* do not release constraint here - will be done later */
      }
   }

   for (i = maxdim - 1; i >= 0; --i)
   {
      SCIPfreeBufferArray(scip, &orbitopematrix[i]);
   }
   SCIPfreeBufferArray(scip, &orbitopematrix);

   return SCIP_OKAY;
}

/** tries to handle symmetries of single lex matrices (orbitopes) or double lex matrices */
static
SCIP_RETCODE tryHandleSingleOrDoubleLexMatricesComponent(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   SCIP_Bool             detectsinglelex,    /**< whether single lex matrices shall be detected */
   int                   cidx                /**< index of component */
   )
{
   int** lexmatrix = NULL;
   int* lexrowsbegin = NULL;
   int* lexcolsbegin = NULL;
   int nrows;
   int ncols;
   int nrowmatrices;
   int ncolmatrices;
   int** perms;
   int compsize;
   int i;
   int p;
   SCIP_Bool isorbitope;
   SCIP_Bool success = FALSE;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( 0 <= cidx && cidx < propdata->ncomponents );

   /* exit if component is already blocked */
   if ( propdata->componentblocked[cidx] )
      return SCIP_OKAY;

   /* exit if component has non-standard permutations */
   if ( propdata->componenthassignedperm[cidx] )
      return SCIP_OKAY;

   /* exit if polyhedral methods are disabled when looking for double lex matrices */
   if ( !ISSYMRETOPESACTIVE(propdata->usesymmetry) && !detectsinglelex )
      return SCIP_OKAY;

   /* get permutations of component */
   compsize = propdata->componentbegins[cidx + 1] - propdata->componentbegins[cidx];
   SCIP_CALL( SCIPallocBufferArray(scip, &perms, compsize) );
   for (p = 0, i = propdata->componentbegins[cidx]; i < propdata->componentbegins[cidx + 1]; ++i)
      perms[p++] = propdata->perms[propdata->components[i]];

   SCIP_CALL( SCIPdetectSingleOrDoubleLexMatrices(scip, detectsinglelex, perms, compsize, propdata->npermvars,
         &success, &isorbitope, &lexmatrix, &nrows, &ncols,
         &lexrowsbegin, &lexcolsbegin, &nrowmatrices, &ncolmatrices) );

   SCIPfreeBufferArray(scip, &perms);

   /* possibly handle double lex matrix or orbitope */
   if ( success )
   {
      assert( lexmatrix != NULL );
      assert( nrows > 0 );
      assert( ncols > 0 );

      if ( isorbitope )
      {
         SCIP_CALL( handleOrbitope(scip, propdata, cidx, lexmatrix, nrows, ncols, &success) );
      }
      else
      {
         SCIP_Bool hasbinaryvar = FALSE;

         /* check whether a binary variable is contained in the matrix */
         for (i = 0; i < nrows && !hasbinaryvar; ++i)
         {
            for (p = 0; p < ncols; ++p)
            {
               if ( SCIPvarIsBinary(propdata->permvars[lexmatrix[i][p]]) )
               {
                  hasbinaryvar = TRUE;
                  break;
               }
            }
         }

         if ( hasbinaryvar )
         {
            SCIP_CALL( handleDoublelLexMatrix(scip, propdata, cidx, lexmatrix, nrows, ncols,
                  lexrowsbegin, lexcolsbegin, nrowmatrices, ncolmatrices, &success) );
         }
         else
            success = FALSE;
      }

      /* free memory not needed anymore */
      for (i = nrows - 1; i >= 0; --i)
      {
         SCIPfreeBlockMemoryArray(scip, &lexmatrix[i], ncols);
      }
      SCIPfreeBlockMemoryArray(scip, &lexmatrix, nrows);
      if ( ncolmatrices > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &lexcolsbegin, ncolmatrices);
      }
      if ( nrowmatrices > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &lexrowsbegin, nrowmatrices);
      }
   }

   if ( success )
   {
      propdata->componentblocked[cidx] |= SYM_HANDLETYPE_SYMBREAK;
      ++(propdata->ncompblocked);
   }

   return SCIP_OKAY;
}

/** tries to handle subgroups of component */
static
SCIP_RETCODE tryHandleSubgroups(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   int                   cidx                /**< index of component */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );
   assert( 0 <= cidx && cidx < propdata->ncomponents );

   /* exit if component is already blocked */
   if ( propdata->componentblocked[cidx] )
      return SCIP_OKAY;

   /* skip component if it has signed permutations */
   if ( propdata->componenthassignedperm[cidx] )
      return SCIP_OKAY;

   /* only run if subgroups shall be detected and we can handle them */
   if ( !propdata->usedynamicprop && ISSYMRETOPESACTIVE(propdata->usesymmetry) && propdata->detectsubgroups
      && propdata->binvaraffected && propdata->ncompblocked < propdata->ncomponents )
   {
      /* @todo also implement a dynamic variant */
      SCIP_CALL( detectAndHandleSubgroups(scip, propdata, cidx) );
   }

   return SCIP_OKAY;
}


/** tries to add symmetry handling methods to component of symmetry group
 *
 * For a component, we handle the symmetries as follows:
 * 1. If orbitope detection is enabled and the component is an orbitope: Apply one of the following:
 *   1.1. If dynamic symmetry handling methods are used:
 *     1.1.1. If the orbitope has a single row, add linear constraints x_1 >= x_2 ... >= x_n.
 *     1.1.2. If it has only two columns only, use lexicographic reduction; cf. symmetry_lexred.c
 *     1.1.3. If there are at least 3 binary rows with packing-partitioning constraints,
 *       use a static packing-partitioning orbitopal fixing; cf. cons_orbitope.c
 *       @todo make a dynamic adaptation for packing-partitioning orbitopes.
 *     1.1.4. If none of these standard cases apply, use dynamic orbitopal reduction; cf. symmetry_orbitopal.c
 *   1.2. If static symmetry handling methods are used: Use static orbitopal fixing (binary variables only);
 *     cf. cons_orbitope.c
 * 2. If no dynamic symmetry handling methods are used, and if (orbitopal) subgroup detection is enabled,
 *      detect those and add static orbitopes if necessary.
 * 3. Otherwise, if orbital reduction is enabled, or if dynamic methods are enabled and lexicographic reduction
 *     propagations can be applied:
 *   3.1. If orbital reduction is enabled: Use orbital reduction.
 *   3.2. And, if dynamic methods and lexicographic for single permutations reduction are enabled, use that.
 * 4. Otherwise, if possible, use SST cuts.
 * 5. Otherwise, if possible, add symresacks (lexicographic reduction on binary variables using a static ordering).
 */
static
SCIP_RETCODE tryAddSymmetryHandlingMethodsComponent(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry propagator */
   int                   cidx,               /**< index of component */
   int*                  nchgbds             /**< pointer to store number of bound changes (or NULL)*/
   )
{
   SCIP_Bool useorbitalredorlexred;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->ncomponents >= 0 );
   assert( 0 <= cidx && cidx < propdata->ncomponents );

   /* ignore blocked components */
   if ( propdata->componentblocked[cidx] )
      return SCIP_OKAY;

   /* detect if orbital reduction or lexicographic reduction shall be applied */
   useorbitalredorlexred = ISORBITALREDUCTIONACTIVE(propdata->usesymmetry)
      || ( ISSYMRETOPESACTIVE(propdata->usesymmetry) && propdata->usedynamicprop && propdata->addsymresacks );

   /* try to apply symmetry handling methods */
   if ( propdata->detectdoublelex || propdata->detectorbitopes )
   {
      SCIP_Bool detectsinglelex;

      detectsinglelex = propdata->detectdoublelex ? FALSE : TRUE;

      SCIP_CALL( tryHandleSingleOrDoubleLexMatricesComponent(scip, propdata, detectsinglelex, cidx) );
   }
   SCIP_CALL( tryHandleSubgroups(scip, propdata, cidx) );
   if ( ISSSTACTIVE(propdata->usesymmetry) )
   {
      SCIP_CALL( addSSTConss(scip, propdata, useorbitalredorlexred, nchgbds, cidx) );
   }
   if ( useorbitalredorlexred )
   {
      SCIP_CALL( tryAddOrbitalRedLexRed(scip, propdata, cidx) );
   }
   SCIP_CALL( addSymresackConss(scip, propdata, cidx) );

   return SCIP_OKAY;
}


/** determines problem symmetries and activates symmetry handling methods */
static
SCIP_RETCODE tryAddSymmetryHandlingMethods(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< symmetry breaking propagator */
   int*                  nchgbds,            /**< pointer to store number of bound changes (or NULL)*/
   SCIP_Bool*            earlyterm           /**< pointer to store whether we terminated early (or NULL) */
   )
{
   SCIP_PROPDATA* propdata;
   int c;

   assert( prop != NULL );
   assert( scip != NULL );

   if ( nchgbds != NULL )
      *nchgbds = 0;
   if ( earlyterm != NULL )
      *earlyterm = FALSE;

   /* only allow symmetry handling methods if strong and weak dual reductions are permitted */
   if ( !SCIPallowStrongDualReds(scip) || !SCIPallowWeakDualReds(scip) )
   {
      if ( earlyterm != NULL )
         *earlyterm = TRUE;
      return SCIP_OKAY;
   }

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* if no symmetries may be handled, stop here */
   if ( propdata->usesymmetry == 0 )
   {
      if ( earlyterm != NULL )
         *earlyterm = TRUE;
      return SCIP_OKAY;
   }

   /* if symmetry handling methods have already been added */
   if ( propdata->triedaddsymmethods )
   {
      assert( propdata->nperms >= 0 );

      if ( earlyterm != NULL )
         *earlyterm = TRUE;

      return SCIP_OKAY;
   }
   assert( !propdata->triedaddsymmethods );

   /* compute symmetries, if it is not computed before */
   if ( !propdata->computedsymmetry )
   {
      /* verify that no symmetry information is present */
      assert( checkSymmetryDataFree(propdata) );
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }

   /* stop if symmetry computation failed, the reason should be given inside determineSymmetry */
   if ( !propdata->computedsymmetry )
      return SCIP_OKAY;

   /* mark that symmetry handling methods are now tried to be added */
   propdata->triedaddsymmethods = TRUE;
   assert( propdata->nperms >= 0 );

   /* no symmetries present, so nothing to be handled */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   /* compute components of symmetry group */
   SCIP_CALL( ensureSymmetryComponentsComputed(scip, propdata) );
   assert( propdata->ncomponents > 0 );

   /* iterate over components and handle each by suitable symmetry handling methods */
   for (c = 0; c < propdata->ncomponents; ++c)
   {
      SCIP_CALL( tryAddSymmetryHandlingMethodsComponent(scip, propdata, c, nchgbds) );

      if ( SCIPisStopped(scip) || propdata->ncompblocked >= propdata->ncomponents )
         break;
   }

#ifdef SYMMETRY_STATISTICS
   SCIP_CALL( SCIPdisplaySymmetryStatistics(scip, propdata) );
#endif

   return SCIP_OKAY;
}


/** apply propagation methods for various symmetry handling constraints */
static
SCIP_RETCODE propagateSymmetry(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Bool*            infeasible,         /**< pointer for storing feasibility state */
   int*                  nred,               /**< pointer for number of reductions */
   SCIP_Bool*            didrun              /**< pointer for storing whether a propagator actually ran */
   )
{
   int nredlocal;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( infeasible != NULL );
   assert( nred != NULL );
   assert( didrun != NULL );

   *nred = 0;
   *infeasible = FALSE;
   *didrun = FALSE;

   /* apply orbitopal reduction */
   SCIP_CALL( SCIPorbitopalReductionPropagate(scip, propdata->orbitopalreddata, infeasible, &nredlocal, didrun) );
   *nred += nredlocal;
   if ( *infeasible )
      return SCIP_OKAY;

   /* apply orbital reduction */
   SCIP_CALL( SCIPorbitalReductionPropagate(scip, propdata->orbitalreddata, infeasible, &nredlocal, didrun) );
   *nred += nredlocal;
   if ( *infeasible )
      return SCIP_OKAY;

   /* apply dynamic lexicographic reduction */
   SCIP_CALL( SCIPlexicographicReductionPropagate(scip, propdata->lexreddata, infeasible, &nredlocal, didrun) );
   *nred += nredlocal;
   if ( *infeasible )
      return SCIP_OKAY;

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* get nonlinear conshdlr for future checks on whether there are nonlinear constraints */
   propdata->conshdlr_nonlinear = SCIPfindConshdlr(scip, "nonlinear");

   /* check whether we should run */
   if ( propdata->usesymmetry < 0 )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &propdata->usesymmetry) );
   }
   assert( propdata->usesymmetry >= 0 );

   /* terminate early if no symmetries will be handled */
   if ( propdata->usesymmetry == 0 )
      return SCIP_OKAY;

   /* add symmetry handling constraints if required  */
   if ( propdata->addconsstiming == 0 )
   {
      SCIPdebugMsg(scip, "Try to add symmetry handling constraints before presolving.\n");

      SCIP_CALL( tryAddSymmetryHandlingMethods(scip, prop, NULL, NULL) );
   }
   else if ( propdata->symcomptiming == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Symmetry computation before presolving:\n");

      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of propagator (called after presolving has been finished) */
static
SCIP_DECL_PROPEXITPRE(propExitpreSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Exitpre method of propagator <%s> ...\n", PROP_NAME);

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* terminate early if no symmetries will be handled */
   if ( propdata->usesymmetry == 0 )
      return SCIP_OKAY;

   /* guarantee that symmetries are computed (and handled) if the solving process has not been interrupted
    * and even if presolving has been disabled */
   if ( SCIPgetStatus(scip) == SCIP_STATUS_UNKNOWN )
   {
      SCIP_CALL( tryAddSymmetryHandlingMethods(scip, prop, NULL, NULL) );
   }

   /* if timing requests it, guarantee that symmetries are computed even if presolving is disabled */
   if ( propdata->symcomptiming <= 1 && SCIPgetStatus(scip) == SCIP_STATUS_UNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Symmetry computation at end of presolving:\n");

      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolSymmetry)
{
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Exitpre method of propagator <%s> ...\n", PROP_NAME);

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* reset symmetry handling propagators that depend on the branch-and-bound tree structure */
   SCIP_CALL( resetDynamicSymmetryHandling(scip, propdata) );

   return SCIP_OKAY;
} /*lint !e715*/


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;
   int noldngenconns;
   int nchanges;
   SCIP_Bool earlyterm;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* terminate early if no symmetries will be handled */
   if ( propdata->usesymmetry == 0 )
      return SCIP_OKAY;

   /* possibly create symmetry handling constraints */

   /* skip presolving if we are not at the end if addconsstiming == 2 */
   assert( 0 <= propdata->addconsstiming && propdata->addconsstiming <= SYM_COMPUTETIMING_AFTERPRESOL );
   if ( propdata->addconsstiming > SYM_COMPUTETIMING_DURINGPRESOL && ! SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   noldngenconns = propdata->ngenorbconss + propdata->nsstconss + propdata->ngenlinconss;

   SCIP_CALL( tryAddSymmetryHandlingMethods(scip, prop, &nchanges, &earlyterm) );

   /* if we actually tried to add symmetry handling constraints */
   if ( ! earlyterm ) /*lint !e774*/
   {
      *result = SCIP_DIDNOTFIND;

      if ( nchanges > 0 )
      {
         *result = SCIP_SUCCESS;
         *nchgbds += nchanges;
      }

      /* if symmetry handling constraints have been added, presolve each */
      if ( propdata->ngenorbconss > 0 || propdata->ngenlinconss > 0 || propdata->nsstconss > 0 )
      {
         /* at this point, the symmetry group should be computed and nontrivial */
         assert( propdata->nperms > 0 );
         assert( propdata->triedaddsymmethods );

         /* we have added at least one symmetry handling constraints, i.e., we were successful */
         *result = SCIP_SUCCESS;

         *naddconss += propdata->ngenorbconss + propdata->ngenlinconss + propdata->nsstconss - noldngenconns;
         SCIPdebugMsg(scip, "Added symmetry breaking constraints: %d.\n", *naddconss);

         /* if constraints have been added, loop through generated constraints and presolve each */
         for (i = 0; i < propdata->ngenorbconss; ++i)
         {
            SCIP_CALL( SCIPpresolCons(scip, propdata->genorbconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                  nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                  nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

            /* exit if cutoff or unboundedness has been detected */
            if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
            {
               SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->genorbconss[i]));
               return SCIP_OKAY;
            }
         }

         for (i = 0; i < propdata->ngenlinconss; ++i)
         {
            SCIP_CALL( SCIPpresolCons(scip, propdata->genlinconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                  nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                  nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

            /* exit if cutoff or unboundedness has been detected */
            if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
            {
               SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->genlinconss[i]));
               return SCIP_OKAY;
            }
         }
         SCIPdebugMsg(scip, "Presolved %d generated constraints.\n",
            propdata->ngenorbconss + propdata->ngenlinconss);

         for (i = 0; i < propdata->nsstconss; ++i)
         {
            SCIP_CALL( SCIPpresolCons(scip, propdata->sstconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                  nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                  nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

            /* exit if cutoff or unboundedness has been detected */
            if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
            {
               SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->sstconss[i]));
               return SCIP_OKAY;
            }
         }
         SCIPdebugMsg(scip, "Presolved %d generated Schreier Sims constraints.\n", propdata->nsstconss);
      }
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible;
   SCIP_Bool didrun;
   int nred;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* usesymmetry must be read in order for propdata to have initialized symmetry handling propagators */
   if ( propdata->usesymmetry < 0 )
      return SCIP_OKAY;

   SCIP_CALL( propagateSymmetry(scip, propdata, &infeasible, &nred, &didrun) );

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      propdata->symfoundreduction = TRUE;
      return SCIP_OKAY;
   }
   if ( nred > 0 )
   {
      assert( didrun );
      *result = SCIP_REDUCEDDOM;
      propdata->symfoundreduction = TRUE;
   }
   else if ( didrun )
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitSymmetry)
{
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Exiting propagator <%s>.\n", PROP_NAME);

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( freeSymmetryData(scip, propdata) );

   /* reset basic data */
   propdata->usesymmetry = -1;
   propdata->triedaddsymmethods = FALSE;
   propdata->nsymresacks = 0;
   propdata->norbitopes = 0;
   propdata->lastrestart = 0;
   propdata->symfoundreduction = FALSE;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator
 *
 *  @todo Implement reverse propagation.
 *
 *  Note that this is relatively difficult to obtain: One needs to include all bounds of variables that are responsible
 *  for creating the orbit in which the variables that was propagated lies. This includes all variables that are moved
 *  by the permutations which are involved in creating the orbit.
 */
static
SCIP_DECL_PROPRESPROP(propRespropSymmetry)
{  /*lint --e{715,818}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIPdebugMsg(scip, "Freeing symmetry propagator.\n");

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->customsymopnodetypes != NULL );

   SCIPhashmapFree(&propdata->customsymopnodetypes);

   assert( propdata->lexreddata != NULL );
   SCIP_CALL( SCIPlexicographicReductionFree(scip, &propdata->lexreddata) );

   assert( propdata->orbitalreddata != NULL );
   SCIP_CALL( SCIPorbitalReductionFree(scip, &propdata->orbitalreddata) );

   assert( propdata->orbitopalreddata != NULL );
   SCIP_CALL( SCIPorbitopalReductionFree(scip, &propdata->orbitopalreddata) );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


/*
 * External methods
 */

/** include symmetry propagator */
SCIP_RETCODE SCIPincludePropSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;
   SCIP_PROPDATA* propdata = NULL;
   SCIP_PROP* prop = NULL;
   SCIP_DIALOG* rootdialog;
   SCIP_DIALOG* displaymenu;
   SCIP_DIALOG* dialog;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   assert( propdata != NULL );

   propdata->npermvars = 0;
   propdata->nbinpermvars = 0;
   propdata->permvars = NULL;
   propdata->nperms = -1;
   propdata->nmaxperms = 0;
   propdata->perms = NULL;
   propdata->permstrans = NULL;
   propdata->permvarmap = NULL;
   propdata->permvardomaincenter = NULL;

   propdata->ncomponents = -1;
   propdata->ncompblocked = 0;
   propdata->components = NULL;
   propdata->componentbegins = NULL;
   propdata->vartocomponent = NULL;
   propdata->componentblocked = NULL;
   propdata->componenthassignedperm = NULL;

   propdata->log10groupsize = -1.0;
   propdata->nmovedvars = -1;
   propdata->binvaraffected = FALSE;
   propdata->computedsymmetry = FALSE;
   propdata->conshdlr_nonlinear = NULL;

   propdata->usesymmetry = -1;
   propdata->triedaddsymmethods = FALSE;
   propdata->genorbconss = NULL;
   propdata->genlinconss = NULL;
   propdata->ngenorbconss = 0;
   propdata->ngenlinconss = 0;
   propdata->genorbconsssize = 0;
   propdata->genlinconsssize = 0;
   propdata->nsymresacks = 0;
   propdata->norbitopes = 0;
   propdata->isnonlinvar = NULL;

   propdata->nmovedpermvars = -1;
   propdata->nmovedbinpermvars = 0;
   propdata->nmovedintpermvars = 0;
   propdata->nmovedimplintpermvars = 0;
   propdata->nmovedcontpermvars = 0;
   propdata->lastrestart = 0;
   propdata->symfoundreduction = FALSE;

   propdata->sstconss = NULL;
   propdata->nsstconss = 0;
   propdata->maxnsstconss = 0;
   propdata->leaders = NULL;
   propdata->nleaders = 0;
   propdata->maxnleaders = 0;

   SCIP_CALL( SCIPhashmapCreate(&propdata->customsymopnodetypes, SCIPblkmem(scip), 10) );
   propdata->nopnodetypes = (int) SYM_CONSOPTYPE_LAST;

   /* include constraint handler */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC,
         PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, propExecSymmetry, propdata) );
   assert( prop != NULL );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSymmetry) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSymmetry) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreSymmetry) );
   SCIP_CALL( SCIPsetPropExitpre(scip, prop, propExitpreSymmetry) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolSymmetry) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropSymmetry) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolSymmetry, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );

   /* include table */
   SCIP_CALL( SCIPallocBlockMemory(scip, &tabledata) );
   tabledata->propdata = propdata;
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_SYMMETRY, TABLE_DESC_SYMMETRY, TRUE,
         NULL, tableFreeSymmetry, NULL, NULL, NULL, NULL, tableOutputSymmetry,
         tabledata, TABLE_POSITION_SYMMETRY, TABLE_EARLIEST_SYMMETRY) );

   /* include display dialog */
   rootdialog = SCIPgetRootDialog(scip);
   assert(rootdialog != NULL);
   if( SCIPdialogFindEntry(rootdialog, "display", &displaymenu) != 1 )
   {
      SCIPerrorMessage("display sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( ! SCIPdialogHasEntry(displaymenu, "symmetries") );
   SCIP_CALL( SCIPincludeDialog(scip, &dialog,
      NULL, dialogExecDisplaySymmetry, NULL, NULL,
      "symmetry", "display generators of symmetry group in cycle notation, if available",
         FALSE, (SCIP_DIALOGDATA*)propdata) );
   SCIP_CALL( SCIPaddDialogEntry(scip, displaymenu, dialog) );
   SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );

   /* add parameters for computing symmetry */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxgenerators",
         "limit on the number of generators that should be produced within symmetry detection (0 = no limit)",
         &propdata->maxgenerators, TRUE, DEFAULT_MAXGENERATORS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/checksymmetries",
         "Should all symmetries be checked after computation?",
         &propdata->checksymmetries, TRUE, DEFAULT_CHECKSYMMETRIES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/displaynorbitvars",
         "Should the number of variables affected by some symmetry be displayed?",
         &propdata->displaynorbitvars, TRUE, DEFAULT_DISPLAYNORBITVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/doubleequations",
         "Double equations to positive/negative version?",
         &propdata->doubleequations, TRUE, DEFAULT_DOUBLEEQUATIONS, NULL, NULL) );

   /* add parameters for adding symmetry handling constraints */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/conssaddlp",
         "Should the symmetry breaking constraints be added to the LP?",
         &propdata->conssaddlp, TRUE, DEFAULT_CONSSADDLP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addsymresacks",
         "Add inequalities for symresacks for each generator?",
         &propdata->addsymresacks, TRUE, DEFAULT_ADDSYMRESACKS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/detectdoublelex",
         "Should we check whether the components of the symmetry group can be handled by double lex matrices?",
         &propdata->detectdoublelex, TRUE, DEFAULT_DETECTDOUBLELEX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/detectorbitopes",
         "Should we check whether the components of the symmetry group can be handled by orbitopes?",
         &propdata->detectorbitopes, TRUE, DEFAULT_DETECTORBITOPES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/detectsubgroups",
         "Should we try to detect symmetric subgroups of the symmetry group on binary variables?",
         &propdata->detectsubgroups, TRUE, DEFAULT_DETECTSUBGROUPS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addweaksbcs",
         "Should we add weak SBCs for enclosing orbit of symmetric subgroups?",
         &propdata->addweaksbcs, TRUE, DEFAULT_ADDWEAKSBCS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/addconsstiming",
         "timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving)",
         &propdata->addconsstiming, TRUE, DEFAULT_ADDCONSSTIMING, 0, 2, NULL, NULL) );

   /* add parameters for orbital reduction */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/ofsymcomptiming",
         "timing of symmetry computation (0 = before presolving, 1 = during presolving, 2 = at first call)",
         &propdata->symcomptiming, TRUE, DEFAULT_SYMCOMPTIMING, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/performpresolving",
         "run orbital fixing during presolving? (disabled)",
         NULL, TRUE, DEFAULT_PERFORMPRESOLVING, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/recomputerestart",
         "recompute symmetries after a restart has occurred? (0 = never)",
         &propdata->recomputerestart, TRUE, DEFAULT_RECOMPUTERESTART, 0, 0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/compresssymmetries",
         "Should non-affected variables be removed from permutation to save memory?",
         &propdata->compresssymmetries, TRUE, DEFAULT_COMPRESSSYMMETRIES, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "propagating/" PROP_NAME "/compressthreshold",
         "Compression is used if percentage of moved vars is at most the threshold.",
         &propdata->compressthreshold, TRUE, DEFAULT_COMPRESSTHRESHOLD, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usecolumnsparsity",
         "Should the number of conss a variable is contained in be exploited in symmetry detection?",
         &propdata->usecolumnsparsity, TRUE, DEFAULT_USECOLUMNSPARSITY, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxnconsssubgroup",
         "maximum number of constraints up to which subgroup structures are detected",
         &propdata->maxnconsssubgroup, TRUE, DEFAULT_MAXNCONSSSUBGROUP, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/usedynamicprop",
         "whether dynamified symmetry handling constraint methods should be used",
         &propdata->usedynamicprop, TRUE, DEFAULT_USEDYNAMICPROP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addstrongsbcs",
         "Should strong SBCs for enclosing orbit of symmetric subgroups be added if orbitopes are not used?",
         &propdata->addstrongsbcs, TRUE, DEFAULT_ADDSTRONGSBCS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/ssttiebreakrule",
         "rule to select the orbit in Schreier Sims inequalities (variable in 0: minimum size orbit; 1: maximum size orbit; 2: orbit with most variables in conflict with leader)",
         &propdata->ssttiebreakrule, TRUE, DEFAULT_SSTTIEBREAKRULE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/sstleaderrule",
         "rule to select the leader in an orbit (0: first var; 1: last var; 2: var having most conflicting vars in orbit)",
         &propdata->sstleaderrule, TRUE, DEFAULT_SSTLEADERRULE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/sstleadervartype",
         "bitset encoding which variable types can be leaders (1: binary; 2: integer; 4: impl. int; 8: continuous);" \
         "if multiple types are allowed, take the one with most affected vars",
         &propdata->sstleadervartype, TRUE, DEFAULT_SSTLEADERVARTYPE, 1, 15, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addconflictcuts",
         "Should Schreier Sims constraints be added if we use a conflict based rule?",
         &propdata->addconflictcuts, TRUE, DEFAULT_ADDCONFLICTCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/sstaddcuts",
         "Should Schreier Sims constraints be added?",
         &propdata->sstaddcuts, TRUE, DEFAULT_SSTADDCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/sstmixedcomponents",
         "Should Schreier Sims constraints be added if a symmetry component contains variables of different types?",
         &propdata->sstmixedcomponents, TRUE, DEFAULT_SSTMIXEDCOMPONENTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/symfixnonbinaryvars",
         "Whether all non-binary variables shall be not affected by symmetries if OF is active? (disabled)",
         NULL, TRUE, DEFAULT_SYMFIXNONBINARYVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/enforcecomputesymmetry",
         "Is only symmetry on binary variables used?",
         &propdata->enforcecomputesymmetry, TRUE, DEFAULT_ENFORCECOMPUTESYMMETRY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/preferlessrows",
         "Shall orbitopes with less rows be preferred in detection?",
         &propdata->preferlessrows, TRUE, DEFAULT_PREFERLESSROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/symtype",
         "Type of symmetries that shall be computed?",
         &propdata->symtype, TRUE, DEFAULT_SYMTYPE, 0, 1, NULL, NULL) );

   /* possibly add description */
   if ( SYMcanComputeSymmetry() )
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetName(), SYMsymmetryGetDesc()) );
      if ( SYMsymmetryGetAddName() != NULL )
      {
         SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetAddName(), SYMsymmetryGetAddDesc()) );
      }
   }

   /* depending functionality */
   SCIP_CALL( SCIPincludeEventHdlrShadowTree(scip, &propdata->shadowtreeeventhdlr) );
   assert( propdata->shadowtreeeventhdlr != NULL );

   SCIP_CALL( SCIPincludeOrbitopalReduction(scip, &propdata->orbitopalreddata) );
   assert( propdata->orbitopalreddata != NULL );

   SCIP_CALL( SCIPincludeOrbitalReduction(scip, &propdata->orbitalreddata, propdata->shadowtreeeventhdlr) );
   assert( propdata->orbitalreddata != NULL );

   SCIP_CALL( SCIPincludeLexicographicReduction(scip, &propdata->lexreddata, propdata->shadowtreeeventhdlr) );
   assert( propdata->lexreddata != NULL );

   return SCIP_OKAY;
}


/** return currently available symmetry group information */
SCIP_RETCODE SCIPgetSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP**        permvarmap,         /**< pointer to store hash map of permvars (or NULL) */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix (or NULL)*/
   int***                permstrans,         /**< pointer to store permutation generators as (npermvars x nperms) matrix (or NULL)*/
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of group size (or NULL) */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether binary variables are affected (or NULL) */
   int**                 components,         /**< pointer to store components of symmetry group (or NULL) */
   int**                 componentbegins,    /**< pointer to store begin positions of components in components array (or NULL) */
   int**                 vartocomponent,     /**< pointer to store assignment from variable to its component (or NULL) */
   int*                  ncomponents         /**< pointer to store number of components (or NULL) */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( perms != NULL || permstrans != NULL );
   assert( ncomponents != NULL || (components == NULL && componentbegins == NULL && vartocomponent == NULL) );

   /* find symmetry propagator */
   prop = SCIPfindProp(scip, "symmetry");
   if ( prop == NULL )
   {
      SCIPerrorMessage("Could not find symmetry propagator.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *npermvars = propdata->npermvars;
   *permvars = propdata->permvars;

   if ( permvarmap != NULL )
   {
      if ( propdata->nperms > 0 )
      {
         SCIP_CALL( ensureSymmetryPermvarmapComputed(scip, propdata) );
      }
      *permvarmap = propdata->permvarmap;
   }

   *nperms = propdata->nperms;
   if ( perms != NULL )
   {
      *perms = propdata->perms;
      assert( *perms != NULL || *nperms <= 0 );
   }

   if ( permstrans != NULL )
   {
      if ( propdata->nperms > 0 )
      {
         SCIP_CALL( ensureSymmetryPermstransComputed(scip, propdata) );
      }
      *permstrans = propdata->permstrans;
      assert( *permstrans != NULL || *nperms <= 0 );
   }

   if ( log10groupsize != NULL )
      *log10groupsize = propdata->log10groupsize;

   if ( binvaraffected != NULL )
      *binvaraffected = propdata->binvaraffected;

   if ( components != NULL || componentbegins != NULL || vartocomponent != NULL || ncomponents != NULL )
   {
      if ( propdata->nperms > 0 )
      {
         SCIP_CALL( ensureSymmetryComponentsComputed(scip, propdata) );
      }
   }

   if ( components != NULL )
      *components = propdata->components;

   if ( componentbegins != NULL )
      *componentbegins = propdata->componentbegins;

   if ( vartocomponent )
      *vartocomponent = propdata->vartocomponent;

   if ( ncomponents )
      *ncomponents = propdata->ncomponents;

   return SCIP_OKAY;
}


/** return number of the symmetry group's generators */
int SCIPgetSymmetryNGenerators(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );

   prop = SCIPfindProp(scip, PROP_NAME);
   if ( prop == NULL )
      return 0;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   if ( propdata->nperms < 0 )
      return 0;
   else
      return propdata->nperms;
}

/** creates new operator node type (used for symmetry detection) and returns its representation
 *
 *  If the operator node already exists, the function terminates with SCIP_INVALIDDATA.
 */
SCIP_RETCODE SCIPcreateSymOpNodeType(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           opnodename,         /**< name of new operator node type */
   int*                  nodetype            /**< pointer to store the node type */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( nodetype != NULL );

   prop = SCIPfindProp(scip, PROP_NAME);
   if ( prop == NULL )
   {
      SCIPerrorMessage("Cannot create operator node type, symmetry propagator has not been included.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->customsymopnodetypes != NULL );

   if ( SCIPhashmapExists(propdata->customsymopnodetypes, (void*) opnodename) )
   {
      SCIPerrorMessage("Cannot create operator node type %s, it already exists.\n", opnodename);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPhashmapInsertInt(propdata->customsymopnodetypes, (void*) opnodename, propdata->nopnodetypes) );
   *nodetype = propdata->nopnodetypes++;

   return SCIP_OKAY;
}

/** returns representation of an operator node type.
 *
 *  If the node type does not already exist, a new node type will be created.
 */
SCIP_RETCODE SCIPgetSymOpNodeType(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           opnodename,         /**< name of new operator node type */
   int*                  nodetype            /**< pointer to store the node type */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );

   prop = SCIPfindProp(scip, PROP_NAME);
   if ( prop == NULL )
   {
      SCIPerrorMessage("Cannot return operator node type, symmetry propagator has not been included.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->customsymopnodetypes != NULL );

   if ( ! SCIPhashmapExists(propdata->customsymopnodetypes, (void*) opnodename) )
   {
      SCIP_CALL( SCIPcreateSymOpNodeType(scip, opnodename, nodetype) );
   }
   else
      *nodetype = SCIPhashmapGetImageInt(propdata->customsymopnodetypes, (void*) opnodename);

   return SCIP_OKAY;
}
