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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_symmetry.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  propagator for handling symmetries
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 *
 * This propagator combines the following symmetry handling functionalities:
 * - It allows to compute symmetries of the problem and to store this information in adequate form. The symmetry
 *   information can be accessed through external functions.
 * - It allows to add the following symmetry breaking constraints:
 *    - symresack constraints, which separate minimal cover inequalities
 *    - orbitope constraints, if special symmetry group structures are detected
 * - It allows to apply orbital fixing.
 *
 *
 * @section SYMCOMP Symmetry Computation
 *
 * The following comments apply to symmetry computation.
 *
 * - The generic functionality of the compute_symmetry.h interface is used.
 * - We treat implicit integer variables as if they were continuous/real variables. The reason is that there is currently
 *   no distinction between implicit integer and implicit binary. Moreover, currently implicit integer variables hurt
 *   our code more than continuous/real variables (we basically do not handle integral variables at all).
 * - We do not copy symmetry information, since it is not clear how this information transfers. Moreover, copying
 *   symmetry might inhibit heuristics. But note that solving a sub-SCIP might then happen without symmetry information!
 *
 *
 * @section SYMBREAK Symmetry Handling Constraints
 *
 * The following comments apply to adding symmetry handling constraints.
 *
 * - The code automatically detects whether symmetry substructures like symresacks or orbitopes are present and possibly
 *   adds the corresponding constraints.
 * - If orbital fixing is active, only orbitopes are added (if present) and no symresacks.
 * - We try to compute symmetry as late as possible and then add constraints based on this information.
 * - Currently, we only allocate memory for pointers to symresack constraints for group generators. If further
 *   constraints are considered, we have to reallocate memory.
 *
 *
 * @section OF Orbital Fixing
 *
 * Orbital fixing is implemented as introduced by@n
 * F. Margot: Exploiting orbits in symmetric ILP. Math. Program., 98(1-3):3–21, 2003.
 *
 * The method computes orbits of variables with respect to the subgroup of the symmetry group that stabilizes the
 * variables globally fixed or branched to 1. Then one can fix all variables in an orbit to 0 or 1 if one of the other
 * variables in the orbit is fixed to 0 or 1, respectively. Different from Margot, the subgroup is obtained by filtering
 * out generators that do not individually stabilize the variables branched to 1.
 *
 * @pre All variable fixings applied by other components are required to be strict, i.e., if one variable is fixed to
 *      a certain value v, all other variables in the same variable orbit can be fixed to v as well, c.f.@n
 *      F. Margot: Symmetry in integer linear programming. 50 Years of Integer Programming, 647-686, Springer 2010.
 *
 * To illustrate this, consider the example \f$\max\{x_1 + x_2 : x_1 + x_2 \leq 1, Ay \leq b,
 * (x,y) \in \{0,1\}^{2 + n}\} \f$. Since \f$x_1\f$ and \f$x_2\f$ are independent from the remaining problem, the
 * setppc constraint handler may fix \f$(x_1,x_2) = (1,0)\f$. However, since both variables are symmetric, this setting
 * is not strict (if it was strict, both variables would have been set to the same value) and orbital fixing would
 * declare this subsolution as infeasible (there exists an orbit of non-branching variables that are fixed to different
 * values). To avoid this situation, we have to assume that all non-strict settings fix variables globally, i.e., we
 * can take care of it by taking variables into account that have been globally fixed to 1. In fact, it suffices to
 * consider one kind of global fixings since stabilizing one kind prevents an orbit to contain variables that have
 * been fixed globally to different values.
 *
 * @pre All non-strict settings are global settings, since otherwise, we cannot (efficiently) take care of them.
 *
 * @pre No non-strict setting algorithm is interrupted early (e.g., by a time or iteration limit), since this may lead to
 * wrong decisions by orbital fixing as well. For example, if cons_setppc in the above toy example starts by fixing
 * \f$x_2 = 0\f$ and is interrupted afterwards, orbital fixing detects that the orbit \f$\{x_1, x_2\}\f$ contains
 * one variable that is fixed to 0, and thus, it fixes \f$x_1\f$ to 0 as well. Thus, after these reductions, every
 * feasible solution has objective 0 which is not optimal. This situation would not occur if the non-strict setting is
 * complete, because then \f$x_1\f$ is globally fixed to 1, and thus, is stabilized in orbital fixing.
 *
 * Note that orbital fixing might lead to wrong results if it is called in repropagation of a node, because the path
 * from the node to the root might have been changed. Thus, the stabilizers of global 1-fixing and 1-branchings of the
 * initial propagation and repropagation might differ, which may cause conflicts. For this reason, orbital fixing cannot
 * be called in repropagation.
 *
 * @note If, besides orbital fixing, also symmetry handling constraints shall be added, orbital fixing is only applied
 *       to symmetry components that are not handled by orbitope constraints.
 *
 * @todo Possibly turn off propagator in subtrees.
 * @todo Check application of conflict resolution.
 * @todo Check whether one should switch the role of 0 and 1
 * @todo Implement stablizer computation?
 * @todo Implement isomorphism pruning?
 * @todo Implement particular preprocessing rules
 * @todo Separate permuted cuts (first experiments not successful)
 * @todo Allow the computation of local symmetries
 * @todo Order rows of orbitopes (in particular packing/partitioning) w.r.t. cliques in conflict graph.
 */
/* #define SCIP_OUTPUT */
/* #define SCIP_OUTPUT_COMPONENT */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
#include <scip/cons_knapsack.h>
#include <scip/cons_varbound.h>
#include <scip/cons_setppc.h>
#include <scip/cons_and.h>
#include <scip/cons_logicor.h>
#include <scip/cons_or.h>
#include "scip/cons_orbitope.h"
#include "scip/cons_symresack.h"
#include <scip/cons_xor.h>
#include <scip/cons_linking.h>
#include <scip/cons_bounddisjunction.h>
#include <scip/misc.h>

#include <scip/prop_symmetry.h>
#include <symmetry/compute_symmetry.h>
#include <scip/symmetry.h>

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
#define DEFAULT_DOUBLEEQUATIONS      FALSE   /**< Double equations to positive/negative version? */
#define DEFAULT_COMPRESSSYMMETRIES   TRUE    /**< Should non-affected variables be removed from permutation to save memory? */
#define DEFAULT_COMPRESSTHRESHOLD     0.5    /**< Compression is used if percentage of moved vars is at most the threshold. */
#define DEFAULT_SYMFIXNONBINARYVARS FALSE    /**< Whether all non-binary variables shall be not affected by symmetries if OF is active? */

/* default parameters for symmetry constraints */
#define DEFAULT_CONSSADDLP           TRUE    /**< Should the symmetry breaking constraints be added to the LP? */
#define DEFAULT_ADDSYMRESACKS        TRUE    /**< Add inequalities for symresacks for each generator? */
#define DEFAULT_DETECTORBITOPES      TRUE    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
#define DEFAULT_ADDCONSSTIMING          2    /**< timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving) */

/* default parameters for orbital fixing */
#define DEFAULT_OFSYMCOMPTIMING         2    /**< timing of symmetry computation for orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call) */
#define DEFAULT_PERFORMPRESOLVING   FALSE    /**< Run orbital fixing during presolving? */
#define DEFAULT_RECOMPUTERESTART    FALSE    /**< Recompute symmetries after a restart has occurred? */
#define DEFAULT_DISABLEOFRESTART    FALSE    /**< whether OF shall be disabled if OF has found a reduction and a restart occurs */


/* event handler properties */
#define EVENTHDLR_SYMMETRY_NAME    "symmetry"
#define EVENTHDLR_SYMMETRY_DESC    "filter global variable fixing event handler for orbital fixing"

/* output table properties */
#define TABLE_NAME_ORBITALFIXING        "orbitalfixing"
#define TABLE_DESC_ORBITALFIXING        "orbital fixing statistics"
#define TABLE_POSITION_ORBITALFIXING    7001                    /**< the position of the statistics table */
#define TABLE_EARLIEST_ORBITALFIXING    SCIP_STAGE_SOLVING      /**< output of the statistics table is only printed from this stage onwards */


/* other defines */
#define MAXGENNUMERATOR          64000000    /**< determine maximal number of generators by dividing this number by the number of variables */
#define SCIP_SPECIALVAL 1.12345678912345e+19 /**< special floating point value for handling zeros in bound disjunctions */
#define COMPRESSNVARSLB             25000    /**< lower bound on the number of variables above which compression could be performed */

/* macros for getting activeness of symmetry handling methods */
#define ISSYMRETOPESACTIVE(x)      (((unsigned) x & SYM_HANDLETYPE_SYMBREAK) != 0)
#define ISORBITALFIXINGACTIVE(x)   (((unsigned) x & SYM_HANDLETYPE_ORBITALFIXING) != 0)



/** propagator data */
struct SCIP_PropData
{
   /* symmetry group information */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nbinpermvars;       /**< number of binary variables for permuations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
#ifndef NDEBUG
   SCIP_Real*            permvarsobj;        /**< objective values of permuted variables (for debugging) */
#endif
   int                   nperms;             /**< number of permutations */
   int                   nmaxperms;          /**< maximal number of permutations (needed for freeing storage) */
   int**                 perms;              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   int**                 permstrans;         /**< pointer to store transposed permutation generators as (npermvars x nperms) matrix */
   SCIP_HASHMAP*         permvarmap;         /**< map of variables to indices in permvars array */

   /* components of symmetry group */
   int                   ncomponents;        /**< number of components of symmetry group */
   int*                  components;         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins;    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent;     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   SCIP_Shortbool*       componentblocked;   /**< array to store whether a component is blocked to be considered by
                                              *   further symmetry handling techniques */

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
   SCIP_Bool             symfixnonbinaryvars; /**< Whether all non-binary variables shall be not affected by symmetries if OF is active? */

   /* for symmetry constraints */
   SCIP_Bool             symconsenabled;     /**< Should symmetry constraints be added? */
   SCIP_Bool             triedaddconss;      /**< whether we already tried to add symmetry breaking constraints */
   SCIP_Bool             conssaddlp;         /**< Should the symmetry breaking constraints be added to the LP? */
   SCIP_Bool             addsymresacks;      /**< Add symresack constraints for each generator? */
   int                   addconsstiming;     /**< timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving) */
   SCIP_CONS**           genconss;           /**< list of generated constraints */
   int                   ngenconss;          /**< number of generated constraints */
   int                   nsymresacks;        /**< number of symresack constraints */
   SCIP_Bool             detectorbitopes;    /**< Should we check whether the components of the symmetry group can be handled by orbitopes? */
   int                   norbitopes;         /**< number of orbitope constraints */

   /* data necessary for orbital fixing */
   SCIP_Bool             ofenabled;          /**< Run orbital fixing? */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for handling global variable fixings */
   SCIP_Shortbool*       bg0;                /**< bitset to store variables globally fixed to 0 */
   int*                  bg0list;            /**< list of variables globally fixed to 0 */
   int                   nbg0;               /**< number of variables in bg0 and bg0list */
   SCIP_Shortbool*       bg1;                /**< bitset to store variables globally fixed or branched to 1 */
   int*                  bg1list;            /**< list of variables globally fixed or branched to 1 */
   int                   nbg1;               /**< number of variables in bg1 and bg1list */
   int*                  permvarsevents;     /**< stores events caught for permvars */
   SCIP_Shortbool*       inactiveperms;      /**< array to store whether permutations are inactive */
   int                   nmovedpermvars;     /**< number of variables moved by any permutation in a symmetry component that is handled by OF */
   SCIP_Bool             performpresolving;  /**< Run orbital fixing during presolving? */
   SCIP_Bool             recomputerestart;   /**< Recompute symmetries after a restart has occured? */
   int                   ofsymcomptiming;    /**< timing of orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call) */
   int                   lastrestart;        /**< last restart for which symmetries have been computed */
   int                   nfixedzero;         /**< number of variables fixed to 0 */
   int                   nfixedone;          /**< number of variables fixed to 1 */
   SCIP_Longint          nodenumber;         /**< number of node where propagation has been last applied */
   SCIP_Bool             offoundreduction;   /**< whether orbital fixing has found a reduction since the last time computing symmetries */
   SCIP_Bool             disableofrestart;   /**< whether OF shall be disabled if OF has found a reduction and a restart occurs */
};



/*
 * Event handler callback methods
 */

/** exec the event handler for handling global variable bound changes (necessary for orbital fixing)
 *
 *  Global variable fixings during the solving process might arise because parts of the tree are pruned or if certain
 *  preprocessing steps are performed that do not correspond to strict setting algorithms. Since these fixings might be
 *  caused by or be in conflict with orbital fixing, they can be in conflict with the symmetry handling decisions of
 *  orbital fixing in the part of the tree that is not pruned. Thus, we have to take global fixings into account when
 *  filtering out symmetries.
 */
static
SCIP_DECL_EVENTEXEC(eventExecSymmetry)
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR* var;
   int varidx;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_SYMMETRY_NAME) == 0 );
   assert( event != NULL );

   propdata = (SCIP_PROPDATA*) eventdata;
   assert( propdata != NULL );
   assert( propdata->permvarmap != NULL );
   assert( propdata->permstrans != NULL );
   assert( propdata->nperms > 0 );
   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );

   /* get fixed variable */
   var = SCIPeventGetVar(event);
   assert( var != NULL );
   assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );

   if ( ! SCIPhashmapExists(propdata->permvarmap, (void*) var) )
   {
      SCIPerrorMessage("Invalid variable.\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }
   varidx = SCIPhashmapGetImageInt(propdata->permvarmap, (void*) var);
   assert( 0 <= varidx && varidx < propdata->npermvars );

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_GUBCHANGED )
   {
      assert( SCIPisEQ(scip, SCIPeventGetNewbound(event), 0.0) );
      assert( SCIPisEQ(scip, SCIPeventGetOldbound(event), 1.0) );

      SCIPdebugMsg(scip, "Mark variable <%s> as globally fixed to 0.\n", SCIPvarGetName(var));
      assert( ! propdata->bg0[varidx] );
      propdata->bg0[varidx] = TRUE;
      propdata->bg0list[propdata->nbg0++] = varidx;
      assert( propdata->nbg0 <= propdata->npermvars );
   }

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_GLBCHANGED )
   {
      assert( SCIPisEQ(scip, SCIPeventGetNewbound(event), 1.0) );
      assert( SCIPisEQ(scip, SCIPeventGetOldbound(event), 0.0) );

      SCIPdebugMsg(scip, "Mark variable <%s> as globally fixed to 1.\n", SCIPvarGetName(var));
      assert( ! propdata->bg1[varidx] );
      propdata->bg1[varidx] = TRUE;
      propdata->bg1list[propdata->nbg1++] = varidx;
      assert( propdata->nbg1 <= propdata->npermvars );
   }

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


/** output method of orbital fixing propagator statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputOrbitalfixing)
{
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );
   assert( tabledata->propdata != NULL );

   if ( tabledata->propdata->nperms > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "Orbital fixing     :\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  vars fixed to 0  :%11d\n", tabledata->propdata->nfixedzero);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  vars fixed to 1  :%11d\n", tabledata->propdata->nfixedone);
   }

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeOrbitalfixing)
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

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyVartype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare the types of two variables according to objective, lower and upper bound, variable type, and column sparsity.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQVartype)
{
   SCIP* scip;
   SYM_VARTYPE* k1;
   SYM_VARTYPE* k2;

   scip = (SCIP*) userptr;
   k1 = (SYM_VARTYPE*) key1;
   k2 = (SYM_VARTYPE*) key2;

   /* first check objective coefficients */
   if ( ! SCIPisEQ(scip, k1->obj, k2->obj) )
      return FALSE;

   /* if still undecided, take lower bound */
   if ( ! SCIPisEQ(scip, k1->lb, k2->lb) )
      return FALSE;

   /* if still undecided, take upper bound */
   if ( ! SCIPisEQ(scip, k1->ub, k2->ub) )
      return FALSE;

   /* if still undecided, take variable type */
   if ( k1->type != k2->type )
      return FALSE;

   /* if still undecided, take number of conss var is contained in */
   if ( k1->nconss != k2->nconss )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValVartype)
{  /*lint --e{715}*/
   SYM_VARTYPE* k;

   k = (SYM_VARTYPE*) key;

   return SCIPhashFour(SCIPrealHashCode(k->obj), SCIPrealHashCode(k->lb), SCIPrealHashCode((double) k->nconss), SCIPrealHashCode(k->ub));
}

/** data struct to store arrays used for sorting rhs types */
struct SYM_Sortrhstype
{
   SCIP_Real*            vals;               /**< array of values */
   SYM_RHSSENSE*         senses;             /**< array of senses of rhs */
   int                   nrhscoef;           /**< size of arrays (for debugging) */
};
typedef struct SYM_Sortrhstype SYM_SORTRHSTYPE;

/** sorts rhs types - first by sense, then by value
 *
 *  Due to numerical issues, we first sort by sense, then by value.
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortRhsTypes)
{
   SYM_SORTRHSTYPE* data;
   SCIP_Real diffvals;

   data = (SYM_SORTRHSTYPE*) dataptr;
   assert( 0 <= ind1 && ind1 < data->nrhscoef );
   assert( 0 <= ind2 && ind2 < data->nrhscoef );

   /* first sort by senses */
   if ( data->senses[ind1] < data->senses[ind2] )
      return -1;
   else if ( data->senses[ind1] > data->senses[ind2] )
      return 1;

   /* senses are equal, use values */
   diffvals = data->vals[ind1] - data->vals[ind2];

   if ( diffvals < 0.0 )
      return -1;
   else if ( diffvals > 0.0 )
      return 1;

   return 0;
}

/** sorts matrix coefficients
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortMatCoef)
{
   SCIP_Real diffvals;
   SCIP_Real* vals;

   vals = (SCIP_Real*) dataptr;
   diffvals = vals[ind1] - vals[ind2];

   if ( diffvals < 0.0 )
      return -1;
   else if ( diffvals > 0.0 )
      return 1;

   return 0;
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
   assert( propdata->permvarsevents == NULL );
   assert( propdata->bg0list == NULL );
   assert( propdata->bg0 == NULL );
   assert( propdata->bg1list == NULL );
   assert( propdata->bg1 == NULL );
   assert( propdata->nbg0 == 0 );
   assert( propdata->nbg1 == 0 );
   assert( propdata->genconss == NULL );

   assert( propdata->permvars == NULL );
   assert( propdata->permvarsobj == NULL );
   assert( propdata->inactiveperms == NULL );
   assert( propdata->perms == NULL );
   assert( propdata->permstrans == NULL );
   assert( propdata->npermvars == 0 );
   assert( propdata->nbinpermvars == 0 );
   assert( propdata->nperms == -1 || propdata->nperms == 0 );
   assert( propdata->nmaxperms == 0 );
   assert( propdata->nmovedpermvars == 0 );
   assert( propdata->nmovedvars == -1 );
   assert( propdata->binvaraffected == FALSE );

   assert( propdata->componentblocked == NULL );
   assert( propdata->componentbegins == NULL );
   assert( propdata->components == NULL );
   assert( propdata->ncomponents == -1 );

   return TRUE;
}
#endif


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

   if ( propdata->permvarmap != NULL )
   {
      SCIPhashmapFree(&propdata->permvarmap);
   }

   /* drop events */
   if ( propdata->permvarsevents != NULL )
   {
      assert( propdata->permvars != NULL );
      assert( propdata->npermvars > 0 );

      for (i = 0; i < propdata->npermvars; ++i)
      {
         if ( SCIPvarGetType(propdata->permvars[i]) == SCIP_VARTYPE_BINARY )
         {
            /* If symmetry is computed before presolving, it might happen that some variables are turned into binary
             * variables, for which no event has been catched. Since there currently is no way of checking whether a var
             * event has been caught for a particular variable, we use the stored eventfilter positions. */
            if ( propdata->permvarsevents[i] >= 0 )
            {
               SCIP_CALL( SCIPdropVarEvent(scip, propdata->permvars[i], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                     propdata->eventhdlr, (SCIP_EVENTDATA*) propdata, propdata->permvarsevents[i]) );
            }
         }
      }
      SCIPfreeBlockMemoryArray(scip, &propdata->permvarsevents, propdata->npermvars);
   }

   /*  release variables */
   if ( propdata->binvaraffected )
   {
      for (i = 0; i < propdata->nbinpermvars; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &propdata->permvars[i]) );
      }
   }

   /* free lists for orbitopal fixing */
   if ( propdata->bg0list != NULL )
   {
      assert( propdata->bg0 != NULL );
      assert( propdata->bg1list != NULL );
      assert( propdata->bg1 != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->bg0list, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->bg0, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->bg1list, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->bg1, propdata->npermvars);

      propdata->nbg0 = 0;
      propdata->nbg1 = 0;
   }

   /* other data */
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->inactiveperms, propdata->nperms);

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

   /* free data of added constraints */
   if ( propdata->genconss != NULL )
   {
      assert( propdata->ngenconss > 0 || (ISORBITALFIXINGACTIVE(propdata->usesymmetry) && propdata->norbitopes == 0) );

      /* release constraints */
      for (i = 0; i < propdata->ngenconss; ++i)
      {
         assert( propdata->genconss[i] != NULL );
         SCIP_CALL( SCIPreleaseCons(scip, &propdata->genconss[i]) );
      }

      /* free pointers to symmetry group and binary variables */
      SCIPfreeBlockMemoryArray(scip, &propdata->genconss, propdata->nperms);
      propdata->ngenconss = 0;
   }

   /* free components */
   if ( propdata->ncomponents > 0 )
   {
      assert( propdata->componentblocked != NULL );
      assert( propdata->vartocomponent != NULL );
      assert( propdata->componentbegins != NULL );
      assert( propdata->components != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->componentblocked, propdata->ncomponents);
      SCIPfreeBlockMemoryArray(scip, &propdata->vartocomponent, propdata->npermvars);
      SCIPfreeBlockMemoryArray(scip, &propdata->componentbegins, propdata->ncomponents + 1);
      SCIPfreeBlockMemoryArray(scip, &propdata->components, propdata->nperms);

      propdata->ncomponents = -1;
   }

   /* free main symmetry data */
   if ( propdata->nperms > 0 )
   {
      assert( propdata->permvars != NULL );

      SCIPfreeBlockMemoryArray(scip, &propdata->permvars, propdata->npermvars);

      /* if orbital fixing runs exclusively, propdata->perms was already freed in determineSymmetry() */
      if ( propdata->perms != NULL )
      {
         for (i = 0; i < propdata->nperms; ++i)
         {
            SCIPfreeBlockMemoryArray(scip, &propdata->perms[i], propdata->npermvars);
         }
         SCIPfreeBlockMemoryArray(scip, &propdata->perms, propdata->nmaxperms);
      }

#ifndef NDEBUG
      SCIPfreeBlockMemoryArrayNull(scip, &propdata->permvarsobj, propdata->npermvars);
#endif

      propdata->npermvars = 0;
      propdata->nbinpermvars = 0;
      propdata->nperms = -1;
      propdata->nmaxperms = 0;
      propdata->nmovedpermvars = 0;
      propdata->nmovedvars = -1;
      propdata->log10groupsize = -1.0;
      propdata->binvaraffected = FALSE;
   }
   propdata->nperms = -1;

   assert( checkSymmetryDataFree(propdata) );

   propdata->computedsymmetry = FALSE;
   propdata->compressed = FALSE;

   return SCIP_OKAY;
}


/** deletes symmetry handling constraints */
static
SCIP_RETCODE delSymConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );

   if ( propdata->ngenconss == 0 )
   {
      if ( propdata->genconss != NULL )
         SCIPfreeBlockMemoryArray(scip, &propdata->genconss, propdata->nperms);
      propdata->triedaddconss = FALSE;

      return SCIP_OKAY;
   }
   assert( propdata->genconss != NULL );
   assert( propdata->nperms > 0 );
   assert( propdata->nperms >= propdata->ngenconss );

   for (i = 0; i < propdata->ngenconss; ++i)
   {
      assert( propdata->genconss[i] != NULL );

      SCIP_CALL( SCIPdelCons(scip, propdata->genconss[i]) );
      SCIP_CALL( SCIPreleaseCons(scip, &propdata->genconss[i]) );
   }

   /* free pointers to symmetry group and binary variables */
   SCIPfreeBlockMemoryArray(scip, &propdata->genconss, propdata->nperms);
   propdata->ngenconss = 0;
   propdata->triedaddconss = FALSE;

   return SCIP_OKAY;
}


/** determines whether variable should be fixed by permutations */
static
SCIP_Bool SymmetryFixVar(
   SYM_SPEC              fixedtype,          /**< bitset of variable types that should be fixed */
   SCIP_VAR*             var                 /**< variable to be considered */
   )
{
   if ( (fixedtype & SYM_SPEC_INTEGER) && SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
      return TRUE;
   if ( (fixedtype & SYM_SPEC_BINARY) && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      return TRUE;
   if ( (fixedtype & SYM_SPEC_REAL) &&
      (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT) )
      return TRUE;
   return FALSE;
}


/** Transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant.
 *
 *  @note @p constant needs to be initialized!
 */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( *vars != NULL );
   assert( *scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if ( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if ( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for (v = 0; v < *nvars; ++v)
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&(*vars)[v], &(*scalars)[v], constant) );
      }
   }
   return SCIP_OKAY;
}


/** fills in matrix elements into coefficient arrays */
static
SCIP_RETCODE collectCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             doubleequations,    /**< Double equations to positive/negative version? */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficients values (or NULL if all linear coefficient values are 1) */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             istransformed,      /**< whether the constraint is transformed */
   SYM_RHSSENSE          rhssense,           /**< identifier of constraint type */
   SYM_MATRIXDATA*       matrixdata,         /**< matrix data to be filled in */
   int*                  nconssforvar        /**< pointer to array to store for each var the number of conss */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real constant = 0.0;
   int nrhscoef;
   int nmatcoef;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( nlinvars == 0 || linvars != NULL );
   assert( lhs <= rhs );

   /* do nothing if constraint is empty */
   if ( nlinvars == 0 )
      return SCIP_OKAY;

   /* ignore redundant constraints */
   if ( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* duplicate variable and value array */
   nvars = nlinvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, linvars, nvars) );
   if ( linvals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &vals, linvals, nvars) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      for (j = 0; j < nvars; ++j)
         vals[j] = 1.0;
   }
   assert( vars != NULL );
   assert( vals != NULL );

   /* get active variables */
   SCIP_CALL( getActiveVariables(scip, &vars, &vals, &nvars, &constant, istransformed) );

   /* check whether constraint is empty after transformation to active variables */
   if ( nvars <= 0 )
   {
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
      return SCIP_OKAY;
   }

   /* handle constant */
   if ( ! SCIPisInfinity(scip, -lhs) )
      lhs -= constant;
   if ( ! SCIPisInfinity(scip, rhs) )
      rhs -= constant;

   /* check whether we have to resize; note that we have to add 2 * nvars since two inequalities may be added */
   if ( matrixdata->nmatcoef + 2 * nvars > matrixdata->nmaxmatcoef )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, matrixdata->nmatcoef + 2 * nvars);
      assert( newsize >= 0 );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matidx), matrixdata->nmaxmatcoef, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matrhsidx), matrixdata->nmaxmatcoef, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matvaridx), matrixdata->nmaxmatcoef, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(matrixdata->matcoef), matrixdata->nmaxmatcoef, newsize) );
      SCIPdebugMsg(scip, "Resized matrix coefficients from %u to %d.\n", matrixdata->nmaxmatcoef, newsize);
      matrixdata->nmaxmatcoef = newsize;
   }

   nrhscoef = matrixdata->nrhscoef;
   nmatcoef = matrixdata->nmatcoef;

   /* check lhs/rhs */
   if ( SCIPisEQ(scip, lhs, rhs) )
   {
      SCIP_Bool poscoef = FALSE;
      SCIP_Bool negcoef = FALSE;

      assert( ! SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      matrixdata->rhscoef[nrhscoef] = rhs;

      /* if we deal with special constraints */
      if ( rhssense >= SYM_SENSE_XOR )
         matrixdata->rhssense[nrhscoef] = rhssense;
      else
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_EQUATION;
      matrixdata->rhsidx[nrhscoef] = nrhscoef;

      for (j = 0; j < nvars; ++j)
      {
         assert( nmatcoef < matrixdata->nmaxmatcoef );

         matrixdata->matidx[nmatcoef] = nmatcoef;
         matrixdata->matrhsidx[nmatcoef] = nrhscoef;

         assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

         if ( nconssforvar != NULL )
            nconssforvar[SCIPvarGetProbindex(vars[j])] += 1;
         matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
         matrixdata->matcoef[nmatcoef++] = vals[j];
         if ( SCIPisPositive(scip, vals[j]) )
            poscoef = TRUE;
         else
            negcoef = TRUE;
      }
      nrhscoef++;

      /* add negative of equation; increases chance to detect symmetry, but might increase time to compute symmetry. */
      if ( doubleequations && poscoef && negcoef )
      {
         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;
            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
            matrixdata->matcoef[nmatcoef++] = -vals[j];
         }
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_EQUATION;
         matrixdata->rhsidx[nrhscoef] = nrhscoef;
         matrixdata->rhscoef[nrhscoef++] = -rhs;
      }
   }
   else
   {
#ifndef NDEBUG
      if ( rhssense == SYM_SENSE_BOUNDIS_TYPE_2 )
      {
         assert( ! SCIPisInfinity(scip, -lhs) );
         assert( ! SCIPisInfinity(scip, rhs) );
      }
#endif

      if ( ! SCIPisInfinity(scip, -lhs) )
      {
         matrixdata->rhscoef[nrhscoef] = -lhs;
         if ( rhssense >= SYM_SENSE_XOR )
         {
            assert( rhssense == SYM_SENSE_BOUNDIS_TYPE_2 );
            matrixdata->rhssense[nrhscoef] = rhssense;
         }
         else
            matrixdata->rhssense[nrhscoef] = SYM_SENSE_INEQUALITY;

         matrixdata->rhsidx[nrhscoef] = nrhscoef;

         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;
            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);

            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            if ( nconssforvar != NULL )
               nconssforvar[SCIPvarGetProbindex(vars[j])] += 1;

            matrixdata->matcoef[nmatcoef++] = -vals[j];
         }
         nrhscoef++;
      }

      if ( ! SCIPisInfinity(scip, rhs) )
      {
         matrixdata->rhscoef[nrhscoef] = rhs;
         if ( rhssense >= SYM_SENSE_XOR )
         {
            assert( rhssense == SYM_SENSE_BOUNDIS_TYPE_2 );
            matrixdata->rhssense[nrhscoef] = rhssense;
         }
         else
            matrixdata->rhssense[nrhscoef] = SYM_SENSE_INEQUALITY;

         matrixdata->rhsidx[nrhscoef] = nrhscoef;

         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;

            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            if ( nconssforvar != NULL )
               nconssforvar[SCIPvarGetProbindex(vars[j])] += 1;

            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
            matrixdata->matcoef[nmatcoef++] = vals[j];
         }
         nrhscoef++;
      }
   }
   matrixdata->nrhscoef = nrhscoef;
   matrixdata->nmatcoef = nmatcoef;

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** checks whether given permutations form a symmetry of a MIP
 *
 *  We need the matrix and rhs in the original order in order to speed up the comparison process. The matrix is needed
 *  in the right order to easily check rows. The rhs is used because of cache effects.
 */
static
SCIP_RETCODE checkSymmetriesAreSymmetries(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SPEC              fixedtype,          /**< variable types that must be fixed by symmetries */
   SYM_MATRIXDATA*       matrixdata,         /**< matrix data */
   int                   nperms,             /**< number of permutations */
   int**                 perms               /**< permutations */
   )
{
   SCIP_Real* permrow = 0;
   int* rhsmatbeg = 0;
   int oldrhs;
   int j;
   int p;

   SCIPdebugMsg(scip, "Checking whether symmetries are symmetries (generators: %u).\n", nperms);

   /* set up dense row for permuted row */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &permrow, matrixdata->npermvars) );

   /* set up map between rows and first entry in matcoef array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rhsmatbeg, matrixdata->nrhscoef) );
   for (j = 0; j < matrixdata->nrhscoef; ++j)
      rhsmatbeg[j] = -1;

   /* build map from rhs into matrix */
   oldrhs = -1;
   for (j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int rhs;

      rhs = matrixdata->matrhsidx[j];
      if ( rhs != oldrhs )
      {
         assert( 0 <= rhs && rhs < matrixdata->nrhscoef );
         rhsmatbeg[rhs] = j;
         oldrhs = rhs;
      }
   }

   /* create row */
   for (j = 0; j < matrixdata->npermvars; ++j)
      permrow[j] = 0.0;

   /* check all generators */
   for (p = 0; p < nperms; ++p)
   {
      int* P;
      int r1;
      int r2;

      SCIPdebugMsg(scip, "Verifying automorphism group generator #%d ...\n", p);
      P = perms[p];
      assert( P != NULL );

      for (j = 0; j < matrixdata->npermvars; ++j)
      {
         if ( SymmetryFixVar(fixedtype, matrixdata->permvars[j]) && P[j] != j )
         {
            SCIPdebugMsg(scip, "Permutation does not fix types %u, moving variable %d.\n", fixedtype, j);
            return SCIP_ERROR;
         }
      }

      /* check all constraints == rhs */
      for (r1 = 0; r1 < matrixdata->nrhscoef; ++r1)
      {
         int npermuted = 0;

         /* fill row into permrow (dense) */
         j = rhsmatbeg[r1];
         assert( 0 <= j && j < matrixdata->nmatcoef );
         assert( matrixdata->matrhsidx[j] == r1 ); /* note: row cannot be empty by construction */

         /* loop through row */
         while ( j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r1 )
         {
            int varidx;

            assert( matrixdata->matvaridx[j] < matrixdata->npermvars );
            varidx = P[matrixdata->matvaridx[j]];
            assert( 0 <= varidx && varidx < matrixdata->npermvars );
            if ( varidx != matrixdata->matvaridx[j] )
               ++npermuted;
            assert( SCIPisZero(scip, permrow[varidx]) );
            permrow[varidx] = matrixdata->matcoef[j];
            ++j;
         }

         /* if row is not affected by permutation, we do not have to check it */
         if ( npermuted > 0 )
         {
            /* check other rows (sparse) */
            SCIP_Bool found = FALSE;
            for (r2 = 0; r2 < matrixdata->nrhscoef; ++r2)
            {
               /* a permutation must map constraints of the same type and respect rhs coefficients */
               if ( matrixdata->rhssense[r1] == matrixdata->rhssense[r2] && SCIPisEQ(scip, matrixdata->rhscoef[r1], matrixdata->rhscoef[r2]) )
               {
                  j = rhsmatbeg[r2];
                  assert( 0 <= j && j < matrixdata->nmatcoef );
                  assert( matrixdata->matrhsidx[j] == r2 );
                  assert( matrixdata->matvaridx[j] < matrixdata->npermvars );

                  /* loop through row r2 and check whether it is equal to permuted row r */
                  while (j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r2 && SCIPisEQ(scip, permrow[matrixdata->matvaridx[j]], matrixdata->matcoef[j]) )
                     ++j;

                  /* check whether rows are completely equal */
                  if ( j >= matrixdata->nmatcoef || matrixdata->matrhsidx[j] != r2 )
                  {
                     /* perm[p] is indeed a symmetry */
                     found = TRUE;
                     break;
                  }
               }
            }

            assert( found );
            if ( ! found ) /*lint !e774*/
            {
               SCIPerrorMessage("Found permutation that is not a symmetry.\n");
               return SCIP_ERROR;
            }
         }

         /* reset permrow */
         j = rhsmatbeg[r1];
         while ( j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r1 )
         {
            int varidx;
            varidx = P[matrixdata->matvaridx[j]];
            permrow[varidx] = 0.0;
            ++j;
         }
      }
   }

   SCIPfreeBlockMemoryArray(scip, &rhsmatbeg, matrixdata->nrhscoef);
   SCIPfreeBlockMemoryArray(scip, &permrow, matrixdata->npermvars);

   return SCIP_OKAY;
}


/** returns the number of active constraints that can be handled by symmetry */
static
int getNSymhandableConss(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONSHDLR* conshdlr;
   int nhandleconss = 0;

   assert( scip != NULL );

   conshdlr = SCIPfindConshdlr(scip, "linear");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "linking");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "setppc");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "xor");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "and");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "or");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "logicor");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "knapsack");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "bounddisjunction");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);

   return nhandleconss;
}


/** set symmetry data */
static
SCIP_RETCODE setSymmetryData(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR**            vars,               /**< vars present at time of symmetry computation */
   int                   nvars,              /**< number of vars present at time of symmetry computation */
   int                   nbinvars,           /**< number of binary vars present at time of symmetry computation */
   SCIP_VAR***           permvars,           /**< pointer to permvars array */
   int*                  npermvars,          /**< pointer to store number of permvars */
   int*                  nbinpermvars,       /**< pointer to store number of binary permvars */
   int**                 perms,              /**< permutations matrix (nperms x nvars) */
   int                   nperms,             /**< number of permutations */
   int*                  nmovedvars,         /**< pointer to store number of vars affected by symmetry (if usecompression) or NULL */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether a binary variable is affected by symmetry */
   SCIP_Bool             usecompression,     /**< whether symmetry data shall be compressed */
   SCIP_Real             compressthreshold,  /**< if percentage of moved vars is at most threshold, compression is done */
   SCIP_Bool*            compressed          /**< pointer to store whether compression has been performed */
   )
{
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
            /* iterate over labels and adapt permutation */
            for (i = 0; i < *nmovedvars; ++i)
            {
               assert( i <= labeltopermvaridx[i] );
               perms[p][i] = labelmovedvars[perms[p][labeltopermvaridx[i]]];
            }

            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &perms[p], nvars, *nmovedvars) );
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

   return SCIP_OKAY;
}


/** computes symmetry group of a MIP */
static
SCIP_RETCODE computeSymmetryGroup(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_Bool             doubleequations,    /**< Double equations to positive/negative version? */
   SCIP_Bool             compresssymmetries, /**< Should non-affected variables be removed from permutation to save memory? */
   SCIP_Real             compressthreshold,  /**< Compression is used if percentage of moved vars is at most the threshold. */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_SPEC              fixedtype,          /**< variable types that must be fixed by symmetries */
   SCIP_Bool             local,              /**< Use local variable bounds? */
   SCIP_Bool             checksymmetries,    /**< Should all symmetries be checked after computation? */
   SCIP_Bool             usecolumnsparsity,  /**< Should the number of conss a variable is contained in be exploited in symmetry detection? */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   int*                  nbinpermvars,       /**< pointer to store number of binary variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   int*                  nmovedvars,         /**< pointer to store number of moved vars */
   SCIP_Bool*            binvaraffected,     /**< pointer to store wether a binary variable is affected by symmetry */
   SCIP_Bool*            compressed,         /**< pointer to store whether compression has been performed */
   SCIP_Bool*            success             /**< pointer to store whether symmetry computation was successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SYM_MATRIXDATA matrixdata;
   SCIP_HASHTABLE* vartypemap;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SYM_VARTYPE* uniquevararray;
   SYM_RHSSENSE oldsense = SYM_SENSE_UNKOWN;
   SYM_SORTRHSTYPE sortrhstype;
   SCIP_Real oldcoef = SCIP_INVALID;
   SCIP_Real val;
   int* nconssforvar = NULL;
   int nuniquevararray = 0;
   int nhandleconss;
   int nactiveconss;
   int nconss;
   int nvars;
   int nbinvars;
   int nvarsorig;
   int nallvars;
   int c;
   int j;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( nbinpermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( binvaraffected != NULL );
   assert( compressed != NULL );
   assert( success != NULL );
   assert( SYMcanComputeSymmetry() );

   /* init */
   *npermvars = 0;
   *nbinpermvars = 0;
   *permvars = NULL;
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *nmovedvars = -1;
   *binvaraffected = FALSE;
   *compressed = FALSE;
   *success = FALSE;

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nvarsorig = nvars;

   /* exit if no constraints or no variables are available */
   if ( nconss == 0 || nvars == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   conss = SCIPgetConss(scip);
   assert( conss != NULL );

   /* compute the number of active constraints */
   nactiveconss = SCIPgetNActiveConss(scip);

   /* exit if no active constraints are available */
   if ( nactiveconss == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* before we set up the matrix, check whether we can handle all constraints */
   nhandleconss = getNSymhandableConss(scip);
   assert( nhandleconss <= nactiveconss );
   if ( nhandleconss < nactiveconss )
   {
      /* In this case we found unkown constraints and we exit, since we cannot handle them. */
      *success = FALSE;
      *nperms = -1;
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "Detecting %ssymmetry on %d variables and %d constraints.\n", local ? "local " : "", nvars, nactiveconss);

   /* copy variables */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, SCIPgetVars(scip), nvars) ); /*lint !e666*/
   assert( vars != NULL );

   /* fill matrixdata */

   /* use a staggered scheme for allocating space for non-zeros of constraint matrix since it can become large */
   if ( nvars <= 100000 )
      matrixdata.nmaxmatcoef = 100 * nvars;
   else if ( nvars <= 1000000 )
      matrixdata.nmaxmatcoef = 32 * nvars;
   else if ( nvars <= 16700000 )
      matrixdata.nmaxmatcoef = 16 * nvars;
   else
      matrixdata.nmaxmatcoef = INT_MAX / 10;

   matrixdata.nmatcoef = 0;
   matrixdata.nrhscoef = 0;
   matrixdata.nuniquemat = 0;
   matrixdata.nuniquevars = 0;
   matrixdata.nuniquerhs = 0;
   matrixdata.npermvars = nvars;
   matrixdata.permvars = vars;
   matrixdata.permvarcolors = NULL;
   matrixdata.matcoefcolors = NULL;
   matrixdata.rhscoefcolors = NULL;

   /* prepare matrix data (use block memory, since this can become large) */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhscoef, 2 * nactiveconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhssense, 2 * nactiveconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhsidx, 2 * nactiveconss) );

   /* prepare temporary constraint data (use block memory, since this can become large);
    * also allocate memory for fixed vars since some vars might have been deactivated meanwhile */
   nallvars = nvars + SCIPgetNFixedVars(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvars, nallvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvals, nallvars) );

   /* allocate memory for getting the number of constraints that contain a variable */
   if ( usecolumnsparsity )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &nconssforvar, nvars) );
   }

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      const char* conshdlrname;
      SCIP_CONS* cons;
      SCIP_VAR** linvars;
      int nconsvars;

      /* get constraint */
      cons = conss[c];
      assert( cons != NULL );

      /* skip non-active constraints */
      if ( ! SCIPconsIsActive(cons) )
         continue;

      /* Skip conflict constraints if we are late in the solving process */
      if ( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsConflict(cons) )
         continue;

      /* get constraint handler */
      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( conshdlrname != NULL );

      /* check type of constraint */
      if ( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( collectCoefficients(scip, doubleequations, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
               SCIPgetNVarsLinear(scip, cons), SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons),
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "linking") == 0 )
      {
         SCIP_VAR** curconsvars;
         SCIP_Real* curconsvals;
         int i;

         /* get constraint variables and their coefficients */
         curconsvals = SCIPgetValsLinking(scip, cons);
         SCIP_CALL( SCIPgetBinvarsLinking(scip, cons, &curconsvars, &nconsvars) );
         /* SCIPgetBinVarsLinking returns the number of binary variables, but we also need the integer variable */
         nconsvars++;

         /* copy vars and vals for binary variables */
         for( i = 0; i < nconsvars - 1; i++ )
         {
            consvars[i] = curconsvars[i];
            consvals[i] = (SCIP_Real) curconsvals[i];
         }

         /* set final entry of vars and vals to the linking variable and its coefficient, respectively */
         consvars[nconsvars - 1] = SCIPgetLinkvarLinking(scip, cons);
         consvals[nconsvars - 1] = -1.0;

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata, nconssforvar) );
         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, NULL, nconsvars - 1, 1.0, 1.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "setppc") == 0 )
      {
         linvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( collectCoefficients(scip, doubleequations, linvars, 0, nconsvars, 1.0, 1.0, SCIPconsIsTransformed(cons), SYM_SENSE_EQUATION, &matrixdata, nconssforvar) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( collectCoefficients(scip, doubleequations, linvars, 0, nconsvars, -SCIPinfinity(scip), 1.0, SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( collectCoefficients(scip, doubleequations, linvars, 0, nconsvars, 1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
            break;
         default:
            SCIPerrorMessage("Unknown setppc type %d.\n", SCIPgetTypeSetppc(scip, cons));
            return SCIP_ERROR;
         }
      }
      else if ( strcmp(conshdlrname, "xor") == 0 )
      {
         SCIP_VAR** curconsvars;
         SCIP_VAR* var;

         /* get number of variables of XOR constraint (without integer variable) */
         nconsvars = SCIPgetNVarsXor(scip, cons);

         /* get variables of XOR constraint */
         curconsvars = SCIPgetVarsXor(scip, cons);
         for (j = 0; j < nconsvars; ++j)
         {
            assert( curconsvars[j] != NULL );
            consvars[j] = curconsvars[j];
            consvals[j] = 1.0;
         }

         /* intVar of xor constraint might have been removed */
         var = SCIPgetIntVarXor(scip, cons);
         if ( var != NULL )
         {
            consvars[nconsvars] = var;
            consvals[nconsvars++] = 2.0;
         }
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, (SCIP_Real) SCIPgetRhsXor(scip, cons),
               (SCIP_Real) SCIPgetRhsXor(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_XOR, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "and") == 0 )
      {
         SCIP_VAR** curconsvars;

         /* get number of variables of AND constraint (without resultant) */
         nconsvars = SCIPgetNVarsAnd(scip, cons);

         /* get variables of AND constraint */
         curconsvars = SCIPgetVarsAnd(scip, cons);

         for (j = 0; j < nconsvars; ++j)
         {
            assert( curconsvars[j] != NULL );
            consvars[j] = curconsvars[j];
            consvals[j] = 1.0;
         }

         assert( SCIPgetResultantAnd(scip, cons) != NULL );
         consvars[nconsvars] = SCIPgetResultantAnd(scip, cons);
         consvals[nconsvars++] = 2.0;
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_AND, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "or") == 0 )
      {
         SCIP_VAR** curconsvars;

         /* get number of variables of OR constraint (without resultant) */
         nconsvars = SCIPgetNVarsOr(scip, cons);

         /* get variables of OR constraint */
         curconsvars = SCIPgetVarsOr(scip, cons);

         for (j = 0; j < nconsvars; ++j)
         {
            assert( curconsvars[j] != NULL );
            consvars[j] = curconsvars[j];
            consvals[j] = 1.0;
         }

         assert( SCIPgetResultantOr(scip, cons) != NULL );
         consvars[nconsvars] = SCIPgetResultantOr(scip, cons);
         consvals[nconsvars++] = 2.0;
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_OR, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( collectCoefficients(scip, doubleequations, SCIPgetVarsLogicor(scip, cons), 0, SCIPgetNVarsLogicor(scip, cons),
               1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array and get active variables of constraint */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         for (j = 0; j < nconsvars; ++j)
            consvals[j] = (SCIP_Real) weights[j];
         assert( nconsvars <= nallvars );

         SCIP_CALL( collectCoefficients(scip, doubleequations, SCIPgetVarsKnapsack(scip, cons), consvals, nconsvars, -SCIPinfinity(scip),
               (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "varbound") == 0 )
      {
         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvals[0] = 1.0;

         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
               SCIPgetRhsVarbound(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata, nconssforvar) );
      }
      else if ( strcmp(conshdlrname, "bounddisjunction") == 0 )
      {
         /* To model bound disjunctions, we normalize each constraint
          * \f[
          *   (x_1 \{\leq,\geq\} b_1) \vee \ldots \vee (x_n \{\leq,\geq\} b_n)
          * \f]
          * to a constraint of type
          * \f[
          *   (x_1 \leq b'_1 \vee \ldots \vee (x_n \leq b'_n).
          * \f]
          *
          * If no variable appears twice in such a normalized constraint, we say this bound disjunction
          * is of type 1. If the bound disjunction has length two and both disjunctions contain the same variable,
          * we say the bound disjunction is of type 2. Further bound disjunctions are possible, but can currently
          * not be handled.
          *
          * Bound disjunctions of type 1 are modeled as the linear constraint
          * \f[
          *    b'_1 \cdot x_1 + \ldots +  b'_n \cdot x_n = 0
          * \f]
          * and bound disjunctions of type 2 are modeled as the linear constraint
          * \f[
          *    \min\{b'_1, b'_2\} \leq x_1 \leq \max\{b'_1, b'_2\}.
          * \f]
          * Note that problems arise if \fb'_i = 0\f for some variable \fx_i\f, because its coefficient in the
          * linear constraint is 0. To avoid this, we replace 0 by a special number.
          */
         SCIP_VAR** bounddisjvars;
         SCIP_BOUNDTYPE* boundtypes;
         SCIP_Real* bounds;
         SCIP_Bool repetition = FALSE;
         int nbounddisjvars;
         int k;

         /* collect coefficients for normalized constraint */
         nbounddisjvars = SCIPgetNVarsBounddisjunction(scip, cons);
         bounddisjvars = SCIPgetVarsBounddisjunction(scip, cons);
         boundtypes = SCIPgetBoundtypesBounddisjunction(scip, cons);
         bounds = SCIPgetBoundsBounddisjunction(scip, cons);

         /* copy data */
         for (j = 0; j < nbounddisjvars; ++j)
         {
            consvars[j] = bounddisjvars[j];

            /* normalize bounddisjunctions to SCIP_BOUNDTYPE_LOWER */
            if ( boundtypes[j] == SCIP_BOUNDTYPE_LOWER )
               consvals[j] = - bounds[j];
            else
               consvals[j] = bounds[j];

            /* special treatment of 0 values */
            if ( SCIPisZero(scip, consvals[j]) )
               consvals[j] = SCIP_SPECIALVAL;

            /* detect whether a variable appears in two literals */
            for (k = 0; k < j && ! repetition; ++k)
            {
               if ( consvars[j] == consvars[k] )
                  repetition = TRUE;
            }

            /* stop, we cannot handle bounddisjunctions with more than two variables that contain a variable twice */
            if ( repetition && nbounddisjvars > 2 )
            {
               *success = FALSE;

               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                  "   Deactivated symmetry handling methods, there exist constraints that cannot be handled by symmetry methods.\n");

               if ( usecolumnsparsity )
                  SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvars);

               SCIPfreeBlockMemoryArrayNull(scip, &consvals, nallvars);
               SCIPfreeBlockMemoryArrayNull(scip, &consvars, nallvars);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);
               SCIPfreeBlockMemoryArrayNull(scip, &vars, nvars);

               return SCIP_OKAY;
            }
         }
         assert( ! repetition || nbounddisjvars == 2 );

         /* if no variable appears twice */
         if ( ! repetition )
         {
            /* add information for bounddisjunction of type 1 */
            SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, nbounddisjvars, 0.0, 0.0,
                  SCIPconsIsTransformed(cons), SYM_SENSE_BOUNDIS_TYPE_1, &matrixdata, nconssforvar) );
         }
         else
         {
            /* add information for bounddisjunction of type 2 */
            SCIP_Real lhs;
            SCIP_Real rhs;

            lhs = MIN(consvals[0], consvals[1]);
            rhs = MAX(consvals[0], consvals[1]);

            consvals[0] = 1.0;

            SCIP_CALL( collectCoefficients(scip, doubleequations, consvars, consvals, 1, lhs, rhs,
                  SCIPconsIsTransformed(cons), SYM_SENSE_BOUNDIS_TYPE_2, &matrixdata, nconssforvar) );
         }
      }
      else
      {
         SCIPerrorMessage("Cannot determine symmetries for constraint <%s> of constraint handler <%s>.\n",
            SCIPconsGetName(cons), SCIPconshdlrGetName(conshdlr) );
         return SCIP_ERROR;
      }
   }
   assert( matrixdata.nrhscoef <= 2 * nactiveconss );
   assert( matrixdata.nrhscoef >= 0 );

   SCIPfreeBlockMemoryArray(scip, &consvals, nallvars);
   SCIPfreeBlockMemoryArray(scip, &consvars, nallvars);

   /* if no active constraint contains active variables */
   if ( matrixdata.nrhscoef == 0 )
   {
      *success = TRUE;

      if ( usecolumnsparsity )
         SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvars);

      /* free matrix data */
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

      SCIPfreeBlockMemoryArray(scip, &vars, nvars);

      return SCIP_OKAY;
   }

   /* sort matrix coefficients (leave matrix array intact) */
   SCIPsort(matrixdata.matidx, SYMsortMatCoef, (void*) matrixdata.matcoef, matrixdata.nmatcoef);

   /* sort rhs types (first by sense, then by value, leave rhscoef intact) */
   sortrhstype.vals = matrixdata.rhscoef;
   sortrhstype.senses = matrixdata.rhssense;
   sortrhstype.nrhscoef = matrixdata.nrhscoef;
   SCIPsort(matrixdata.rhsidx, SYMsortRhsTypes, (void*) &sortrhstype, matrixdata.nrhscoef);

   /* create map for variables to indices */
   SCIP_CALL( SCIPhashtableCreate(&vartypemap, SCIPblkmem(scip), 5 * nvars, SYMhashGetKeyVartype, SYMhashKeyEQVartype, SYMhashKeyValVartype, (void*) scip) );
   assert( vartypemap != NULL );

   /* allocate space for mappings to colors */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.permvarcolors, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniquevararray, nvars) );

   /* determine number of different coefficents */

   /* find non-equivalent variables: same objective, lower and upper bounds, and variable type */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;

      var = vars[j];
      assert( var != NULL );

      /* if the variable type should be fixed, just increase the color */
      if ( SymmetryFixVar(fixedtype, var) )
      {
         matrixdata.permvarcolors[j] = matrixdata.nuniquevars++;
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "Detected variable <%s> of fixed type %d - color %d.\n", SCIPvarGetName(var), SCIPvarGetType(var), matrixdata.nuniquevars - 1);
#endif
      }
      else
      {
         SYM_VARTYPE* vt;

         vt = &uniquevararray[nuniquevararray];
         assert( nuniquevararray <= matrixdata.nuniquevars );

         vt->obj = SCIPvarGetObj(var);
         if ( local )
         {
            vt->lb = SCIPvarGetLbLocal(var);
            vt->ub = SCIPvarGetUbLocal(var);
         }
         else
         {
            vt->lb = SCIPvarGetLbGlobal(var);
            vt->ub = SCIPvarGetUbGlobal(var);
         }
         vt->type = SCIPvarGetType(var);
         vt->nconss = usecolumnsparsity ? nconssforvar[j] : 0; /*lint !e613*/

         if ( ! SCIPhashtableExists(vartypemap, (void*) vt) )
         {
            SCIP_CALL( SCIPhashtableInsert(vartypemap, (void*) vt) );
            vt->color = matrixdata.nuniquevars;
            matrixdata.permvarcolors[j] = matrixdata.nuniquevars++;
            ++nuniquevararray;
#ifdef SCIP_OUTPUT
            SCIPdebugMsg(scip, "Detected variable <%s> of new type (probindex: %d, obj: %g, lb: %g, ub: %g, type: %d) - color %d.\n",
               SCIPvarGetName(var), SCIPvarGetProbindex(var), vt->obj, vt->lb, vt->ub, vt->type, matrixdata.nuniquevars - 1);
#endif
         }
         else
         {
            SYM_VARTYPE* vtr;

            vtr = (SYM_VARTYPE*) SCIPhashtableRetrieve(vartypemap, (void*) vt);
            matrixdata.permvarcolors[j] = vtr->color;
         }
      }
   }

   /* If every variable is unique, terminate. -> no symmetries can be present */
   if ( matrixdata.nuniquevars == nvars )
   {
      *success = TRUE;

      /* free matrix data */
      SCIPfreeBlockMemoryArray(scip, &uniquevararray, nvars);

      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.permvarcolors, nvars);
      SCIPhashtableFree(&vartypemap);

      if ( usecolumnsparsity )
         SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvars);

      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
      SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

      SCIPfreeBlockMemoryArray(scip, &vars, nvars);

      return SCIP_OKAY;
   }

   /* find non-equivalent matrix entries (use sorting to avoid too many map calls) */
   for (j = 0; j < matrixdata.nmatcoef; ++j)
   {
      int idx;

      idx = matrixdata.matidx[j];
      assert( 0 <= idx && idx < matrixdata.nmatcoef );

      val = matrixdata.matcoef[idx];
      assert( oldcoef == SCIP_INVALID || oldcoef <= val ); /*lint !e777*/

      if ( ! SCIPisEQ(scip, val, oldcoef) )
      {
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "Detected new matrix entry type %f - color: %d\n.", val, matrixdata.nuniquemat);
#endif
         matrixdata.matcoefcolors[idx] = matrixdata.nuniquemat++;
         oldcoef = val;
      }
      else
      {
         assert( matrixdata.nuniquemat > 0 );
         matrixdata.matcoefcolors[idx] = matrixdata.nuniquemat - 1;
      }
   }

   /* find non-equivalent rhs */
   oldcoef = SCIP_INVALID;
   for (j = 0; j < matrixdata.nrhscoef; ++j)
   {
      SYM_RHSSENSE sense;
      int idx;

      idx = matrixdata.rhsidx[j];
      assert( 0 <= idx && idx < matrixdata.nrhscoef );
      sense = matrixdata.rhssense[idx];
      val = matrixdata.rhscoef[idx];

      /* make sure that new senses are treated with new color */
      if ( sense != oldsense )
         oldcoef = SCIP_INVALID;
      oldsense = sense;
      assert( oldcoef == SCIP_INVALID || oldcoef <= val ); /*lint !e777*/

      /* assign new color to new type */
      if ( ! SCIPisEQ(scip, val, oldcoef) )
      {
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "Detected new rhs type %f, type: %u - color: %d\n", val, sense, matrixdata.nuniquerhs);
#endif
         matrixdata.rhscoefcolors[idx] = matrixdata.nuniquerhs++;
         oldcoef = val;
      }
      else
      {
         assert( matrixdata.nuniquerhs > 0 );
         matrixdata.rhscoefcolors[idx] = matrixdata.nuniquerhs - 1;
      }
   }
   assert( 0 < matrixdata.nuniquevars && matrixdata.nuniquevars <= nvars );
   assert( 0 < matrixdata.nuniquerhs && matrixdata.nuniquerhs <= matrixdata.nrhscoef );
   assert( 0 < matrixdata.nuniquemat && matrixdata.nuniquemat <= matrixdata.nmatcoef );

   SCIPdebugMsg(scip, "Number of detected different variables: %d (total: %d).\n", matrixdata.nuniquevars, nvars);
   SCIPdebugMsg(scip, "Number of detected different rhs types: %d (total: %d).\n", matrixdata.nuniquerhs, matrixdata.nrhscoef);
   SCIPdebugMsg(scip, "Number of detected different matrix coefficients: %d (total: %d).\n", matrixdata.nuniquemat, matrixdata.nmatcoef);

   /* do not compute symmetry if all variables are non-equivalent (unique) or if all matrix coefficients are different */
   if ( matrixdata.nuniquevars < nvars && matrixdata.nuniquemat < matrixdata.nmatcoef )
   {
      /* determine generators */
      SCIP_CALL( SYMcomputeSymmetryGenerators(scip, maxgenerators, &matrixdata, nperms, nmaxperms, perms, log10groupsize) );
      assert( *nperms <= *nmaxperms );

      /* SCIPisStopped() might call SCIPgetGap() which is only available after initpresolve */
      if ( checksymmetries && SCIPgetStage(scip) > SCIP_STAGE_INITPRESOLVE && ! SCIPisStopped(scip) )
      {
         SCIP_CALL( checkSymmetriesAreSymmetries(scip, fixedtype, &matrixdata, *nperms, *perms) );
      }

      if ( *nperms > 0 )
      {
         SCIP_CALL( setSymmetryData(scip, vars, nvars, nbinvars, permvars, npermvars, nbinpermvars, *perms, *nperms,
               nmovedvars, binvaraffected, compresssymmetries, compressthreshold, compressed) );
      }
      else
      {
         SCIPfreeBlockMemoryArray(scip, &vars, nvars);
      }
   }
   else
   {
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
   }
   *success = TRUE;

   /* free matrix data */
   SCIPfreeBlockMemoryArray(scip, &uniquevararray, nvarsorig);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.permvarcolors, nvarsorig);
   SCIPhashtableFree(&vartypemap);

   if ( usecolumnsparsity )
      SCIPfreeBlockMemoryArrayNull(scip, &nconssforvar, nvarsorig);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

   return SCIP_OKAY;
}


/** determines symmetry */
static
SCIP_RETCODE determineSymmetry(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SYM_SPEC              symspecrequire,     /**< symmetry specification for which we need to compute symmetries */
   SYM_SPEC              symspecrequirefixed /**< symmetry specification of variables which must be fixed by symmetries */
   )
{
   SCIP_Bool successful;
   int maxgenerators;
   int nhandleconss;
   int nconss;
   unsigned int type = 0;
   int nvars;
   int j;
   int p;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );
   assert( propdata->ofenabled || propdata->symconsenabled );

   /* do not compute symmetry if reoptimization is enabled */
   if ( SCIPisReoptEnabled(scip) )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      return SCIP_OKAY;
   }

   /* do not compute symmetry if Benders decomposition enabled */
   if ( SCIPgetNActiveBenders(scip) > 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      return SCIP_OKAY;
   }

   /* skip symmetry computation if no graph automorphism code was linked */
   if ( ! SYMcanComputeSymmetry() )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      nconss = SCIPgetNActiveConss(scip);
      nhandleconss = getNSymhandableConss(scip);

      /* print verbMessage only if problem consists of symmetry handable constraints */
      assert( nhandleconss <=  nconss );
      if ( nhandleconss < nconss )
         return SCIP_OKAY;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   Deactivated symmetry handling methods, since SCIP was built without symmetry detector (SYM=none).\n");

      return SCIP_OKAY;
   }

   /* do not compute symmetry if there are active pricers */
   if ( SCIPgetNActivePricers(scip) > 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;
      return SCIP_OKAY;
   }

   /* avoid trivial cases */
   nvars = SCIPgetNVars(scip);
   if ( nvars <= 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

   /* do not compute symmetry if there are no binary variables */
   if ( SCIPgetNBinVars(scip) == 0 )
   {
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

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

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

   /* if a restart occured, either disable orbital fixing... */
   if ( propdata->offoundreduction && propdata->disableofrestart  && SCIPgetNRuns(scip) > propdata->lastrestart )
      propdata->ofenabled = FALSE;
   /* ... or free symmetries after a restart to recompute them later */
   else if ( (propdata->offoundreduction || propdata->recomputerestart)
      && propdata->nperms > 0 && SCIPgetNRuns(scip) > propdata->lastrestart )
   {
      assert( propdata->npermvars > 0 );
      assert( propdata->permvars != NULL );
      assert( ! propdata->ofenabled || propdata->permvarmap != NULL );
      assert( ! propdata->ofenabled || propdata->bg0list != NULL );

      /* reset symmetry information */
      SCIP_CALL( delSymConss(scip, propdata) );
      SCIP_CALL( freeSymmetryData(scip, propdata) );
   }

   /* skip computation if symmetry has already been computed */
   if ( propdata->computedsymmetry )
      return SCIP_OKAY;

   /* skip symmetry computation if there are constraints that cannot be handled by symmetry */
   nconss = SCIPgetNActiveConss(scip);
   nhandleconss = getNSymhandableConss(scip);
   if ( nhandleconss < nconss )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation skipped: there exist constraints that cannot be handled by symmetry methods.\n",
         SCIPgetSolvingTime(scip));

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

   assert( propdata->npermvars == 0 );
   assert( propdata->permvars == NULL );
#ifndef NDEBUG
   assert( propdata->permvarsobj == NULL );
#endif
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
   SCIP_CALL( computeSymmetryGroup(scip, propdata->doubleequations, propdata->compresssymmetries, propdata->compressthreshold,
	 maxgenerators, symspecrequirefixed, FALSE, propdata->checksymmetries, propdata->usecolumnsparsity,
         &propdata->npermvars, &propdata->nbinpermvars, &propdata->permvars, &propdata->nperms, &propdata->nmaxperms,
         &propdata->perms, &propdata->log10groupsize, &propdata->nmovedvars, &propdata->binvaraffected,
         &propdata->compressed, &successful) );

   /* mark that we have computed the symmetry group */
   propdata->computedsymmetry = TRUE;

   /* store restart level */
   propdata->lastrestart = SCIPgetNRuns(scip);

   /* return if not successful */
   if ( ! successful )
   {
      assert( checkSymmetryDataFree(propdata) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) could not compute symmetry\n", SCIPgetSolvingTime(scip));

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

   /* return if no symmetries found */
   if ( propdata->nperms == 0 )
   {
      assert( checkSymmetryDataFree(propdata) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry present\n", SCIPgetSolvingTime(scip));

      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

   /* display statistics */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) symmetry computation finished: %d generators found (max: ",
      SCIPgetSolvingTime(scip), propdata->nperms);

   /* display statistics: maximum number of generators */
   if ( maxgenerators == 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "-");
   else
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%u", maxgenerators);

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
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ")\n");

   /* exit if no binary variables are affected by symmetry */
   if ( ! propdata->binvaraffected )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry on binary variables present.\n", SCIPgetSolvingTime(scip));

      /* free data and exit */
      SCIP_CALL( freeSymmetryData(scip, propdata) );

      /* disable OF and symmetry handling constraints */
      propdata->ofenabled = FALSE;
      propdata->symconsenabled = FALSE;

      return SCIP_OKAY;
   }

   assert( propdata->nperms > 0 );
   assert( propdata->npermvars > 0 );

   /* compute components */
   assert( propdata->components == NULL );
   assert( propdata->componentbegins == NULL );
   assert( propdata->vartocomponent == NULL );

#ifdef SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation started\n", SCIPgetSolvingTime(scip));
#endif

   /* we only need the components for orbital fixing and orbitope detection */
   if ( propdata->ofenabled || ( propdata->symconsenabled && propdata->detectorbitopes ) )
   {
      SCIP_CALL( SCIPcomputeComponentsSym(scip, propdata->perms, propdata->nperms, propdata->permvars,
            propdata->npermvars, FALSE, &propdata->components, &propdata->componentbegins,
            &propdata->vartocomponent, &propdata->componentblocked, &propdata->ncomponents) );
   }

#ifdef SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation finished\n", SCIPgetSolvingTime(scip));
#endif

   /* set up data for OF */
   if ( propdata->ofenabled )
   {
      int componentidx;
      int v;

      /* transpose symmetries matrix here if necessary */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permstrans, propdata->npermvars) );
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->permstrans[v]), propdata->nmaxperms) );
         for (p = 0; p < propdata->nperms; ++p)
         {
            if ( SCIPvarIsBinary(propdata->permvars[v]) )
               propdata->permstrans[v][p] = propdata->perms[p][v];
            else
               propdata->permstrans[v][p] = v; /* ignore symmetry information on non-binary variables */
         }
      }

      /* prepare array for active permutations */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->inactiveperms, propdata->nperms) );
      for (v = 0; v < propdata->nperms; ++v)
         propdata->inactiveperms[v] = FALSE;

      /* create hashmap for storing the indices of variables */
      assert( propdata->permvarmap == NULL );
      SCIP_CALL( SCIPhashmapCreate(&propdata->permvarmap, SCIPblkmem(scip), propdata->npermvars) );

      /* prepare data structures */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permvarsevents, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg0, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg0list, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg1, propdata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->bg1list, propdata->npermvars) );

      /* insert variables into hashmap  */
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPhashmapInsertInt(propdata->permvarmap, propdata->permvars[v], v) );

         propdata->bg0[v] = FALSE;
         propdata->bg1[v] = FALSE;
         propdata->permvarsevents[v] = -1;

         /* collect number of moved permvars that are handled by OF */
         componentidx = propdata->vartocomponent[v];
         if ( componentidx > -1 && ! propdata->componentblocked[componentidx] )
            propdata->nmovedpermvars += 1;

         /* Only catch binary variables, since integer variables should be fixed pointwise; implicit integer variables
          * are not branched on. */
         if ( SCIPvarGetType(propdata->permvars[v]) == SCIP_VARTYPE_BINARY )
         {
            /* catch whether binary variables are globally fixed; also store filter position */
            SCIP_CALL( SCIPcatchVarEvent(scip, propdata->permvars[v], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                  propdata->eventhdlr, (SCIP_EVENTDATA*) propdata, &propdata->permvarsevents[v]) );
         }
      }
      assert( propdata->nbg1 == 0 );
   }

#ifndef NDEBUG
   /* store objective coefficients for debug purposes */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->permvarsobj, propdata->npermvars) );
   for (j = 0; j < propdata->npermvars; ++j)
      propdata->permvarsobj[j] = SCIPvarGetObj(propdata->permvars[j]);
#endif

   /* capture binary variables and forbid multi-aggregation of symmetric variables
    *
    * note: binary variables are in the beginning of pervars
    */
   for (j = 0; j < propdata->nbinpermvars; ++j)
   {
      SCIP_CALL( SCIPcaptureVar(scip, propdata->permvars[j]) );

      if ( propdata->compressed )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->permvars[j]) );
      }
      else
      {
         for (p = 0; p < propdata->nperms; ++p)
         {
            if ( propdata->perms[p][j] != j )
            {
               SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, propdata->permvars[j]) );
               break;
            }
         }
      }
   }

   /* free original perms matrix if no symmetry constraints are added */
   if ( ! propdata->symconsenabled )
   {
      for (p = 0; p < propdata->nperms; ++p)
      {
         SCIPfreeBlockMemoryArray(scip, &(propdata->perms)[p], propdata->npermvars);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &propdata->perms, propdata->nmaxperms);
   }

   return SCIP_OKAY;
}



/*
 * Functions for symmetry constraints
 */


/** checks whether components of the symmetry group can be completely handled by orbitopes */
static
SCIP_RETCODE detectOrbitopes(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROPDATA*        propdata,           /**< pointer to data of symmetry propagator */
   int*                  components,         /**< array containing components of symmetry group */
   int*                  componentbegins,    /**< array containing begin positions of components in components array */
   int                   ncomponents         /**< number of components */
   )
{
   SCIP_VAR** permvars;
   int** perms;
   int npermvars;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( components != NULL );
   assert( componentbegins != NULL );
   assert( ncomponents > 0 );
   assert( propdata->nperms >= 0 );

   /* exit if no symmetry is present */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   assert( propdata->nperms > 0 );
   assert( propdata->perms != NULL );
   assert( propdata->nbinpermvars >= 0 );
   assert( propdata->npermvars >= 0 );
   assert( propdata->permvars != NULL );

   /* exit if no symmetry on binary variables is present */
   if ( propdata->nbinpermvars == 0 )
   {
      assert( ! propdata->binvaraffected );
      return SCIP_OKAY;
   }

   perms = propdata->perms;
   npermvars = propdata->npermvars;
   permvars = propdata->permvars;

   /* iterate over components */
   for (i = 0; i < ncomponents; ++i)
   {
      SCIP_VAR*** vars;
      SCIP_VAR*** varsallocorder;
      SCIP_CONS* cons;
      SCIP_Bool* usedperm;
      SCIP_Bool isorbitope = TRUE;
      SCIP_Bool infeasibleorbitope;
      int** orbitopevaridx;
      int* columnorder;
      int npermsincomponent;
      int ntwocyclescomp = INT_MAX;
      int nfilledcols;
      int nusedperms;
      int* nusedelems;
      int coltoextend;
      int j;
      int row;

      /* get properties of permutations */
      npermsincomponent = componentbegins[i + 1] - componentbegins[i];
      assert( npermsincomponent > 0 );
      for (j = componentbegins[i]; j < componentbegins[i + 1]; ++j)
      {
         SCIP_Bool iscompoftwocycles = FALSE;
         SCIP_Bool allvarsbinary = TRUE;
         int ntwocyclesperm = 0;

         SCIP_CALL( SCIPgetPropertiesPerm(perms[components[j]], permvars, npermvars, &iscompoftwocycles, &ntwocyclesperm, &allvarsbinary) );

         /* if we are checking the first permutation */
         if ( ntwocyclescomp == INT_MAX )
            ntwocyclescomp = ntwocyclesperm;

         /* no or different number of 2-cycles or not all vars binary: permutations cannot generate orbitope */
         if ( ntwocyclescomp == 0 || ntwocyclescomp != ntwocyclesperm || ! allvarsbinary )
         {
            isorbitope = FALSE;
            break;
         }
      }

      /* if no orbitope was detected */
      if ( ! isorbitope )
         continue;
      assert( ntwocyclescomp > 0 );
      assert( ntwocyclescomp < INT_MAX );

      /* iterate over permutations and check whether for each permutation there exists
       * another permutation whose 2-cycles intersect pairwise in exactly one element */

      /* whether a permutation was considered to contribute to orbitope */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &usedperm, npermsincomponent) );
      nusedperms = 0;

      /* orbitope matrix for indices of variables in permvars array */
      SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &orbitopevaridx[j], npermsincomponent + 1) ); /*lint !e866*/
      }

      /* order of columns of orbitopevaridx */
      SCIP_CALL( SCIPallocBufferArray(scip, &columnorder, npermsincomponent + 1) );
      for (j = 0; j < npermsincomponent + 1; ++j)
         columnorder[j] = npermsincomponent + 2;

      /* count how often an element was used in the potential orbitope */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &nusedelems, npermvars) );

      /* fill first two columns of orbitopevaridx matrix */
      row = 0;
      for (j = 0; j < npermvars; ++j)
      {
         int permidx;

         permidx = components[componentbegins[i]];

         /* avoid adding the same 2-cycle twice */
         if ( perms[permidx][j] > j )
         {
            orbitopevaridx[row][0] = j;
            orbitopevaridx[row++][1] = perms[permidx][j];
            nusedelems[j] += 1;
            nusedelems[perms[permidx][j]] += 1;
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

      /* extend orbitopevaridx matrix to the left, i.e., iteratively find new permutations that
       * intersect the last added left column in each row in exactly one entry, starting with
       * column 0 */
      coltoextend = 0;
      for (j = 0; j < npermsincomponent; ++j)
      {  /* lint --e{850} */
         SCIP_Bool success = FALSE;
         SCIP_Bool infeasible = FALSE;

         if ( nusedperms == npermsincomponent )
            break;

         if ( usedperm[j] )
            continue;

         SCIP_CALL( SCIPextendSubOrbitope(orbitopevaridx, ntwocyclescomp, nfilledcols, coltoextend,
               perms[components[componentbegins[i] + j]], TRUE, &nusedelems, &success, &infeasible) );

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
            j = 0; /*lint !e850*/ /* reset j since previous permutations can now intersect with the latest added column */
         }
      }

      if ( ! isorbitope ) /*lint !e850*/
         goto FREEDATASTRUCTURES;

      coltoextend = 1;
      for (j = 0; j < npermsincomponent; ++j)
      {  /*lint --e(850)*/
         SCIP_Bool success = FALSE;
         SCIP_Bool infeasible = FALSE;

         if ( nusedperms == npermsincomponent )
            break;

         if ( usedperm[j] )
            continue;

         SCIP_CALL( SCIPextendSubOrbitope(orbitopevaridx, ntwocyclescomp, nfilledcols, coltoextend,
               perms[components[componentbegins[i] + j]], FALSE, &nusedelems, &success, &infeasible) );

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
            j = 0; /*lint !e850*/ /* reset j since previous permutations can now intersect with the latest added column */
         }
      }

      if ( nusedperms < npermsincomponent ) /*lint !e850*/
         isorbitope = FALSE;

      if ( ! isorbitope )
         goto FREEDATASTRUCTURES;

      /* we have found a potential orbitope, prepare data for orbitope conshdlr */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, ntwocyclescomp) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varsallocorder, ntwocyclescomp) );
      for (j = 0; j < ntwocyclescomp; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &vars[j], npermsincomponent + 1) ); /*lint !e866*/
         varsallocorder[j] = vars[j]; /* to ensure that we can free the buffer in reverse order */
      }

      /* prepare variable matrix (reorder columns of orbitopevaridx) */
      infeasibleorbitope = FALSE;
      SCIP_CALL( SCIPgenerateOrbitopeVarsMatrix(&vars, ntwocyclescomp, npermsincomponent + 1, permvars, npermvars,
            orbitopevaridx, columnorder, nusedelems, &infeasibleorbitope) );

      if ( ! infeasibleorbitope )
      {
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &cons, "orbitope", vars, SCIP_ORBITOPETYPE_FULL,
               ntwocyclescomp, npermsincomponent + 1, TRUE, FALSE,
               propdata->conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         propdata->genconss[propdata->ngenconss++] = cons;
         ++propdata->norbitopes;

         propdata->componentblocked[i] = TRUE;
      }

      /* free data structures */
      for (j = ntwocyclescomp - 1; j >= 0; --j)
      {
         SCIPfreeBufferArray(scip, &varsallocorder[j]);
      }
      SCIPfreeBufferArray(scip, &varsallocorder);
      SCIPfreeBufferArray(scip, &vars);

   FREEDATASTRUCTURES:
      SCIPfreeBufferArray(scip, &nusedelems);
      SCIPfreeBufferArray(scip, &columnorder);
      for (j = ntwocyclescomp - 1; j >= 0; --j)
      {
         SCIPfreeBufferArray(scip, &orbitopevaridx[j]);
      }
      SCIPfreeBufferArray(scip, &orbitopevaridx);
      SCIPfreeBufferArray(scip, &usedperm);
   }

   return SCIP_OKAY;
}


/** adds symresack constraints */
static
SCIP_RETCODE addSymresackConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< symmetry breaking propagator */
   int*                  components,         /**< array containing components of symmetry group */
   int*                  componentbegins,    /**< array containing begin positions of components in components array */
   int                   ncomponents         /**< number of components */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** permvars;
   SCIP_Bool conssaddlp;
   int** perms;
   int nsymresackcons = 0;
   int npermvars;
   int i;
   int p;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
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
   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   conssaddlp = propdata->conssaddlp;

   assert( propdata->nperms <= 0 || perms != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );

   /* if components have not been computed */
   if ( ncomponents == -1 )
   {
      assert( ! propdata->ofenabled );
      assert( ! propdata->detectorbitopes );

      /* loop through perms and add symresack constraints */
      for (p = 0; p < propdata->nperms; ++p)
      {
         SCIP_CONS* cons;
         char name[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_perm%d", p);
         SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[p], permvars, npermvars, FALSE,
               conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         propdata->genconss[propdata->ngenconss++] = cons;
         ++propdata->nsymresacks;
         ++nsymresackcons;
      }
   }
   else
   {
      /* loop through components */
      for (i = 0; i < ncomponents; ++i)
      {
         /* skip components that were treated by different symmetry handling techniques */
         if ( propdata->componentblocked[i] )
            continue;

         /* loop through perms in component i and add symresack constraints */
         for (p = componentbegins[i]; p < componentbegins[i + 1]; ++p)
         {
            SCIP_CONS* cons;
            int permidx;
            char name[SCIP_MAXSTRLEN];

            permidx = components[p];

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symbreakcons_component%d_perm%d", i, permidx);
            SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, name, perms[permidx], permvars, npermvars, FALSE,
                  conssaddlp, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );

            /* do not release constraint here - will be done later */
            propdata->genconss[propdata->ngenconss++] = cons;
            ++propdata->nsymresacks;
            ++nsymresackcons;
         }
      }
   }

   SCIPdebugMsg(scip, "Added %d symresack constraints.\n", nsymresackcons);

   return SCIP_OKAY;
}


/** finds problem symmetries */
static
SCIP_RETCODE tryAddSymmetryHandlingConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< symmetry breaking propagator */
   SCIP_Bool*            earlyterm           /**< pointer to store whether we terminated early  (or NULL) */
   )
{
   SCIP_PROPDATA* propdata;

   assert( prop != NULL );
   assert( scip != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->symconsenabled );

   /* possibly compute symmetry */
   if ( propdata->ofenabled )
   {
      if ( propdata->symfixnonbinaryvars )
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
      }
      else
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
      }
   }
   else
   {
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_INTEGER | SYM_SPEC_REAL, 0) );
   }
   assert( propdata->binvaraffected || ! propdata->symconsenabled );

   /* if constraints have already been added */
   if ( propdata->triedaddconss )
   {
      assert( propdata->nperms > 0 );

      if ( earlyterm != NULL )
         *earlyterm = TRUE;

      return SCIP_OKAY;
   }

   if ( propdata->nperms <= 0 || ! propdata->symconsenabled )
      return SCIP_OKAY;

   assert( propdata->nperms > 0 );
   assert( propdata->binvaraffected );
   propdata->triedaddconss = TRUE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->genconss, propdata->nperms) );

   if ( propdata->detectorbitopes )
   {
      SCIP_CALL( detectOrbitopes(scip, propdata, propdata->components, propdata->componentbegins, propdata->ncomponents) );
   }

   /* disable orbital fixing if all components are handled by orbitopes */
   if ( propdata->ncomponents == propdata->norbitopes )
      propdata->ofenabled = FALSE;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
   {
      if ( propdata->ngenconss == 0 )
      {
         SCIPfreeBlockMemoryArrayNull(scip, &propdata->genconss, propdata->nperms);
      }
      return SCIP_OKAY;
   }

   /* add symmetry breaking constraints if orbital fixing is not used outside orbitopes */
   if ( ! propdata->ofenabled )
   {
      /* exit if no or only trivial symmetry group is available */
      if ( propdata->nperms <= 0 || ! propdata->binvaraffected )
         return SCIP_OKAY;

      if ( propdata->addsymresacks )
      {
         SCIP_CALL( addSymresackConss(scip, prop, propdata->components, propdata->componentbegins, propdata->ncomponents) );
      }
   }

   return SCIP_OKAY;
}



/*
 * Local methods for orbital fixing
 */


/** performs orbital fixing
 *
 *  Note that we do not have to distinguish between variables that have been fixed or branched to 1, since the
 *  stabilizer is with respect to the variables that have been branched to 1. Thus, if an orbit contains a variable that
 *  has been branched to 1, the whole orbit only contains variables that have been branched to 1 - and nothing can be
 *  fixed.
 */
static
SCIP_RETCODE performOrbitalFixing(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR**            permvars,           /**< variables */
   int                   npermvars,          /**< number of variables */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int                   norbits,            /**< number of orbits */
   SCIP_Bool*            infeasible,         /**< pointer to store whether problem is infeasible */
   int*                  nfixedzero,         /**< pointer to store number of variables fixed to 0 */
   int*                  nfixedone           /**< pointer to store number of variables fixed to 1 */
   )
{
   SCIP_Bool tightened;
   int i;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( infeasible != NULL );
   assert( nfixedzero != NULL );
   assert( nfixedone != NULL );
   assert( norbits > 0 );
   assert( orbitbegins[0] == 0 );

   *infeasible = FALSE;
   *nfixedzero = 0;
   *nfixedone = 0;

   /* check all orbits */
   for (i = 0; i < norbits; ++i)
   {
      SCIP_Bool havefixedone = FALSE;
      SCIP_Bool havefixedzero = FALSE;
      SCIP_VAR* var;
      int j;

      /* we only have nontrivial orbits */
      assert( orbitbegins[i+1] - orbitbegins[i] >= 2 );

      /* check all variables in the orbit */
      for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
      {
         assert( 0 <= orbits[j] && orbits[j] < npermvars );
         var = permvars[orbits[j]];
         assert( var != NULL );

         /* check whether variable is not binary (and not implicit integer!) */
         if ( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
         {
            /* skip orbit if there are non-binary variables */
            havefixedone = FALSE;
            havefixedzero = FALSE;
            break;
         }

         /* if variable is fixed to 1 -> can fix all variables in orbit to 1 */
         if ( SCIPvarGetLbLocal(var) > 0.5 )
            havefixedone = TRUE;

         /* check for zero-fixed variables */
         if ( SCIPvarGetUbLocal(var) < 0.5 )
            havefixedzero = TRUE;
      }

      /* check consistency */
      if ( havefixedone && havefixedzero )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }

      /* fix all variables to 0 if there is one variable fixed to 0 */
      if ( havefixedzero )
      {
         assert( ! havefixedone );

         for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         {
            assert( 0 <= orbits[j] && orbits[j] < npermvars );
            var = permvars[orbits[j]];
            assert( var != NULL );

            /* only variables that are not yet fixed to 0 */
            if ( SCIPvarGetUbLocal(var) > 0.5 )
            {
               SCIPdebugMsg(scip, "can fix <%s> (index %d) to 0.\n", SCIPvarGetName(var), orbits[j]);
               assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
               /* due to aggregation, var might already be fixed to 1, so do not put assert here */

               /* do not use SCIPinferBinvarProp(), since conflict analysis is not valid */
               SCIP_CALL( SCIPtightenVarUb(scip, var, 0.0, FALSE, infeasible, &tightened) );
               if ( *infeasible )
                  return SCIP_OKAY;
               if ( tightened )
                  ++(*nfixedzero);
            }
         }
      }

      /* fix all variables to 1 if there is one variable fixed to 1 */
      if ( havefixedone )
      {
         assert( ! havefixedzero );

         for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         {
            assert( 0 <= orbits[j] && orbits[j] < npermvars );
            var = permvars[orbits[j]];
            assert( var != NULL );

            /* only variables that are not yet fixed to 1 */
            if ( SCIPvarGetLbLocal(var) < 0.5)
            {
               SCIPdebugMsg(scip, "can fix <%s> (index %d) to 1.\n", SCIPvarGetName(var), orbits[j]);
               assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
               /* due to aggregation, var might already be fixed to 0, so do not put assert here */

               /* do not use SCIPinferBinvarProp(), since conflict analysis is not valid */
               SCIP_CALL( SCIPtightenVarLb(scip, var, 1.0, FALSE, infeasible, &tightened) );
               if ( *infeasible )
                  return SCIP_OKAY;
               if ( tightened )
                  ++(*nfixedone);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** Gets branching variables on the path to root
 *
 *  The variables are added to bg1 and bg1list, which are prefilled with the variables globally fixed to 1.
 */
static
SCIP_RETCODE computeBranchingVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varmap,             /**< map of variables to indices in vars array */
   SCIP_Shortbool*       bg1,                /**< bitset marking the variables globally fixed or branched to 1 */
   int*                  bg1list,            /**< array to store the variable indices globally fixed or branched to 1 */
   int*                  nbg1                /**< pointer to store the number of variables in bg1 and bg1list */
   )
{
   SCIP_NODE* node;

   assert( scip != NULL );
   assert( varmap != NULL );
   assert( bg1 != NULL );
   assert( bg1list != NULL );
   assert( nbg1 != NULL );
   assert( *nbg1 >= 0 );

   /* get current node */
   node = SCIPgetCurrentNode(scip);

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPprintNodeRootPath(scip, node, NULL) );
#endif

   /* follow path to the root (in the root no domains were changed due to branching) */
   while ( SCIPnodeGetDepth(node) != 0 )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_DOMCHG* domchg;
      SCIP_VAR* branchvar;
      int nboundchgs;
      int i;

      /* get domain changes of current node */
      domchg = SCIPnodeGetDomchg(node);

      /* If we stopped due to a solving limit, it might happen that a non-root node has no domain changes, in all other
       * cases domchg should not be NULL. */
      if ( domchg != NULL )
      {
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
            branchvar = SCIPboundchgGetVar(boundchg);

            /* we only consider binary variables */
            if ( SCIPvarGetType(branchvar) == SCIP_VARTYPE_BINARY )
            {
               /* if branching variable is not known (may have been created meanwhile,
                * e.g., by prop_inttobinary; may have been removed from symmetry data
                * due to compression), continue with parent node */
               if ( ! SCIPhashmapExists(varmap, (void*) branchvar) )
                  break;

               if ( SCIPvarGetLbLocal(branchvar) > 0.5 )
               {
                  int branchvaridx;

                  branchvaridx = SCIPhashmapGetImageInt(varmap, (void*) branchvar);
                  assert( branchvaridx < nvars );

                  /* the variable might already be fixed to 1 */
                  if ( ! bg1[branchvaridx] )
                  {
                     bg1[branchvaridx] = TRUE;
                     bg1list[(*nbg1)++] = branchvaridx;
                  }
               }
            }
         }
      }

      node = SCIPnodeGetParent(node);
   }

   return SCIP_OKAY;
}


/** propagates orbital fixing */
static
SCIP_RETCODE propagateOrbitalFixing(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata,           /**< data of symmetry breaking propagator */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the node is detected to be infeasible */
   int*                  nprop               /**< pointer to store the number of propagations */
   )
{
   SCIP_Shortbool* inactiveperms;
   SCIP_Shortbool* bg0;
   SCIP_Shortbool* bg1;
   SCIP_VAR** permvars;
   int* orbitbegins;
   int* orbits;
   int* components;
   int* componentbegins;
   int* vartocomponent;
   int ncomponents;
   int* bg0list;
   int nbg0;
   int* bg1list;
   int nbg1;
   int nactiveperms;
   int norbits;
   int npermvars;
   int nbinpermvars;
   int** permstrans;
   int nperms;
   int p;
   int v;
   int j;
   int componentidx;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->ofenabled );
   assert( infeasible != NULL );
   assert( nprop != NULL );

   *infeasible = FALSE;
   *nprop = 0;

   /* possibly compute symmetry */
   if ( propdata->symfixnonbinaryvars )
   {
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
   }
   else
   {
      SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
   }
   assert( propdata->binvaraffected || ! propdata->ofenabled );

   /* return if there is no symmetry available */
   nperms = propdata->nperms;
   if ( nperms <= 0 || ! propdata->ofenabled )
      return SCIP_OKAY;

   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->permvarmap != NULL );
   assert( propdata->permstrans != NULL );
   assert( propdata->inactiveperms != NULL );
   assert( propdata->components != NULL );
   assert( propdata->componentbegins != NULL );
   assert( propdata->vartocomponent != NULL );
   assert( propdata->ncomponents > 0 );

   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   nbinpermvars = propdata->nbinpermvars;
   permstrans = propdata->permstrans;
   inactiveperms = propdata->inactiveperms;
   components = propdata->components;
   componentbegins = propdata->componentbegins;
   vartocomponent = propdata->vartocomponent;
   ncomponents = propdata->ncomponents;

   /* init bitset for marking variables (globally fixed or) branched to 1 */
   assert( propdata->bg1 != NULL );
   assert( propdata->bg1list != NULL );
   assert( propdata->nbg1 >= 0 );
   assert( propdata->nbg1 <= npermvars );

   bg1 = propdata->bg1;
   bg1list = propdata->bg1list;
   nbg1 = propdata->nbg1;

   /* get branching variables */
   SCIP_CALL( computeBranchingVariables(scip, npermvars, propdata->permvarmap, bg1, bg1list, &nbg1) );
   assert( nbg1 >= propdata->nbg1 );

   /* reset inactive permutations */
   nactiveperms = nperms;
   for (p = 0; p < nperms; ++p)
      propdata->inactiveperms[p] = FALSE;

   /* get pointers for bg0 */
   assert( propdata->bg0 != NULL );
   assert( propdata->bg0list != NULL );
   assert( propdata->nbg0 >= 0 );
   assert( propdata->nbg0 <= npermvars );

   bg0 = propdata->bg0;
   bg0list = propdata->bg0list;
   nbg0 = propdata->nbg0;

   /* filter out permutations that move variables that are fixed to 0 */
   for (j = 0; j < nbg0 && nactiveperms > 0; ++j)
   {
      int* pt;

      v = bg0list[j];
      assert( 0 <= v && v < npermvars );
      assert( bg0[v] );

      componentidx = vartocomponent[v];

      /* skip unaffected variables and blocked components */
      if ( componentidx < 0 || propdata->componentblocked[componentidx] )
         continue;

      pt = permstrans[v];
      assert( pt != NULL );

      for (p = componentbegins[componentidx]; p < componentbegins[componentidx + 1]; ++p)
      {
         int img;
         int perm;

         perm = components[p];

         /* skip inactive permutations */
         if ( inactiveperms[perm] )
            continue;

         img = pt[perm];

         if ( img != v )
         {
#ifndef NDEBUG
            SCIP_VAR* varv = permvars[v];
            SCIP_VAR* varimg = permvars[img];

            /* check whether moved variables have the same type (might have been aggregated in the meanwhile) */
            assert( SCIPvarGetType(varv) == SCIPvarGetType(varimg) ||
               (SCIPvarIsBinary(varv) && SCIPvarIsBinary(varimg)) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_IMPLINT && SCIPvarGetType(varimg) == SCIP_VARTYPE_CONTINUOUS &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(varimg) == SCIP_VARTYPE_IMPLINT &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) );
            assert( SCIPisEQ(scip, propdata->permvarsobj[v], propdata->permvarsobj[img]) );
#endif

            /* we are moving a variable globally fixed to 0 to a variable not of this type */
            if ( ! bg0[img] )
            {
               inactiveperms[perm] = TRUE; /* mark as inactive */
               --nactiveperms;
            }
         }
      }
   }

   /* filter out permutations that move variables that are fixed to different values */
   for (j = 0; j < nbg1 && nactiveperms > 0; ++j)
   {
      int* pt;

      v = bg1list[j];
      assert( 0 <= v && v < npermvars );
      assert( bg1[v] );

      componentidx = vartocomponent[v];

      /* skip unaffected variables and blocked components */
      if ( componentidx < 0 || propdata->componentblocked[componentidx] )
         continue;

      pt = permstrans[v];
      assert( pt != NULL );

      for (p = componentbegins[componentidx]; p < componentbegins[componentidx + 1]; ++p)
      {
         int img;
         int perm;

         perm = components[p];

         /* skip inactive permutations */
         if ( inactiveperms[perm] )
            continue;

         img = pt[perm];

         if ( img != v )
         {
#ifndef NDEBUG
            SCIP_VAR* varv = permvars[v];
            SCIP_VAR* varimg = permvars[img];

            /* check whether moved variables have the same type (might have been aggregated in the meanwhile) */
            assert( SCIPvarGetType(varv) == SCIPvarGetType(varimg) ||
               (SCIPvarIsBinary(varv) && SCIPvarIsBinary(varimg)) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_IMPLINT && SCIPvarGetType(varimg) == SCIP_VARTYPE_CONTINUOUS &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(varimg) == SCIP_VARTYPE_IMPLINT &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) );
            assert( SCIPisEQ(scip, propdata->permvarsobj[v], propdata->permvarsobj[img]) );
#endif

            /* we are moving a variable globally fixed or branched to 1 to a variable not of this type */
            if ( ! bg1[img] )
            {
               inactiveperms[perm] = TRUE; /* mark as inactive */
               --nactiveperms;
            }
         }
      }
   }

   /* Clean bg1 list - need to do this after the main loop! (Not needed for bg0.)
    * Note that variables globally fixed to 1 are not resetted, since the loop starts at propdata->nbg1. */
   for (j = propdata->nbg1; j < nbg1; ++j)
      bg1[bg1list[j]] = FALSE;

   /* exit if no active permuations left */
   if ( nactiveperms == 0 )
      return SCIP_OKAY;

   /* compute orbits of binary variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, nbinpermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, nbinpermvars) );
   SCIP_CALL( SCIPcomputeOrbitsFilterSym(scip, nbinpermvars, permstrans, nperms, inactiveperms,
         orbits, orbitbegins, &norbits, components, componentbegins, vartocomponent, propdata->componentblocked, ncomponents, propdata->nmovedpermvars) );

   if ( norbits > 0 )
   {
      int nfixedzero = 0;
      int nfixedone = 0;

      SCIPdebugMsg(scip, "Perform orbital fixing on %d orbits (%d active perms).\n", norbits, nactiveperms);
      SCIP_CALL( performOrbitalFixing(scip, permvars, nbinpermvars, orbits, orbitbegins, norbits, infeasible, &nfixedzero, &nfixedone) );

      propdata->nfixedzero += nfixedzero;
      propdata->nfixedone += nfixedone;
      *nprop = nfixedzero + nfixedone;

      SCIPdebugMsg(scip, "Orbital fixings: %d 0s, %d 1s.\n", nfixedzero, nfixedone);
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

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

   /* check whether we should run */
   if ( propdata->usesymmetry < 0 )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &propdata->usesymmetry) );

      if ( ISSYMRETOPESACTIVE(propdata->usesymmetry) )
         propdata->symconsenabled = TRUE;
      else
         propdata->symconsenabled = FALSE;

      if ( ISORBITALFIXINGACTIVE(propdata->usesymmetry) )
         propdata->ofenabled = TRUE;
      else
         propdata->ofenabled = FALSE;
   }

   /* add symmetry handling constraints if required  */
   if ( propdata->symconsenabled && propdata->addconsstiming == 0 )
   {
      SCIPdebugMsg(scip, "Try to add symmetry handling constraints before presolving.");

      SCIP_CALL( tryAddSymmetryHandlingConss(scip, prop, NULL) );
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

   /* guarantee that symmetries are computed (and handled) if the solving process has not been interrupted
    * and even if presolving has been disabled */
   if ( propdata->symconsenabled && SCIPgetStatus(scip) == SCIP_STATUS_UNKNOWN )
   {
      SCIP_CALL( tryAddSymmetryHandlingConss(scip, prop, NULL) );
   }

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );
   assert( propdata->usesymmetry >= 0 );

   /* possibly create symmetry handling constraints */
   if ( propdata->symconsenabled )
   {
      int noldngenconns;
      SCIP_Bool earlyterm = FALSE;

      /* skip presolving if we are not at the end if addconsstiming == 2 */
      assert( 0 <= propdata->addconsstiming && propdata->addconsstiming <= 2 );
      if ( propdata->addconsstiming > 1 && ! SCIPisPresolveFinished(scip) )
         return SCIP_OKAY;

      /* possibly stop */
      if ( SCIPisStopped(scip) )
         return SCIP_OKAY;

      noldngenconns = propdata->ngenconss;

      SCIP_CALL( tryAddSymmetryHandlingConss(scip, prop, &earlyterm) );

      /* if we actually tried to add symmetry handling constraints */
      if ( ! earlyterm )
      {
         *result = SCIP_DIDNOTFIND;

         /* if symmetry handling constraints have been added, presolve each */
         if ( propdata->ngenconss > 0 )
         {
            /* at this point, the symmetry group should be computed and nontrivial */
            assert( propdata->nperms > 0 );
            assert( propdata->triedaddconss );

            /* we have added at least one symmetry handling constraints, i.e., we were successful */
            *result = SCIP_SUCCESS;

            *naddconss += propdata->ngenconss - noldngenconns;
            SCIPdebugMsg(scip, "Added symmetry breaking constraints: %d.\n", propdata->ngenconss - noldngenconns);

            /* if constraints have been added, loop through generated constraints and presolve each */
            for (i = 0; i < propdata->ngenconss; ++i)
            {
               SCIP_CALL( SCIPpresolCons(scip, propdata->genconss[i], nrounds, SCIP_PROPTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
                     nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
                     nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

               /* exit if cutoff or unboundedness has been detected */
               if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
               {
                  SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(propdata->genconss[i]));
                  return SCIP_OKAY;
               }
            }
            SCIPdebugMsg(scip, "Presolved %d generated constraints.\n", propdata->ngenconss);
         }
      }
   }

   /* run OF presolving */
   assert( 0 <= propdata->ofsymcomptiming && propdata->ofsymcomptiming <= 2 );
   if ( propdata->ofenabled && propdata->performpresolving && propdata->ofsymcomptiming <= 1 )
   {
      SCIP_Bool infeasible;
      int nprop;

      /* if we did not tried to add symmetry handling constraints */
      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;

      SCIPdebugMsg(scip, "Presolving <%s>.\n", PROP_NAME);

      SCIP_CALL( propagateOrbitalFixing(scip, propdata, &infeasible, &nprop) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         propdata->offoundreduction = TRUE;
      }
      else if ( nprop > 0 )
      {
         *result = SCIP_SUCCESS;
         *nfixedvars += nprop;
         propdata->offoundreduction = TRUE;
      }
   }
   else if ( propdata->ofenabled && propdata->ofsymcomptiming == 1 )
   {
      /* otherwise compute symmetry if timing requests it */
      if ( propdata->symfixnonbinaryvars )
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY, SYM_SPEC_INTEGER | SYM_SPEC_REAL) );
      }
      else
      {
         SCIP_CALL( determineSymmetry(scip, propdata, SYM_SPEC_BINARY | SYM_SPEC_REAL, SYM_SPEC_INTEGER) );
      }
      assert( propdata->binvaraffected || ! propdata->ofenabled );
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Longint nodenumber;
   int nprop = 0;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* do nothing if we are in a probing node */
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run again in repropagation, since the path to the root might have changed */
   if ( SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* if usesymmetry has not been read so far */
   if ( propdata->usesymmetry < 0 )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &propdata->usesymmetry) );
      if ( ISSYMRETOPESACTIVE(propdata->usesymmetry) )
         propdata->symconsenabled = TRUE;
      else
         propdata->symconsenabled = FALSE;

      if ( ISORBITALFIXINGACTIVE(propdata->usesymmetry) )
         propdata->ofenabled = TRUE;
      else
         propdata->ofenabled = FALSE;
   }

   /* do not propagate if orbital fixing is not enabled */
   if ( ! propdata->ofenabled )
      return SCIP_OKAY;

   /* return if there is no symmetry available */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   /* return if we already ran in this node */
   nodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   if ( nodenumber == propdata->nodenumber )
      return SCIP_OKAY;
   propdata->nodenumber = nodenumber;

   /* propagate */
   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "Propagating <%s>.\n", SCIPpropGetName(prop));

   SCIP_CALL( propagateOrbitalFixing(scip, propdata, &infeasible, &nprop) );

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      propdata->offoundreduction = TRUE;
   }
   else if ( nprop > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      propdata->offoundreduction = TRUE;
   }

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
   propdata->symconsenabled = FALSE;
   propdata->triedaddconss = FALSE;
   propdata->nsymresacks = 0;
   propdata->norbitopes = 0;
   propdata->ofenabled = FALSE;
   propdata->lastrestart = 0;
   propdata->nfixedzero = 0;
   propdata->nfixedone = 0;
   propdata->nodenumber = -1;
   propdata->offoundreduction = FALSE;

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

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   assert( propdata != NULL );

   propdata->npermvars = 0;
   propdata->nbinpermvars = 0;
   propdata->permvars = NULL;
#ifndef NDEBUG
   propdata->permvarsobj = NULL;
#endif
   propdata->nperms = -1;
   propdata->nmaxperms = 0;
   propdata->perms = NULL;
   propdata->permstrans = NULL;
   propdata->permvarmap = NULL;

   propdata->ncomponents = -1;
   propdata->components = NULL;
   propdata->componentbegins = NULL;
   propdata->vartocomponent = NULL;
   propdata->componentblocked = NULL;

   propdata->log10groupsize = -1.0;
   propdata->nmovedvars = -1;
   propdata->binvaraffected = FALSE;
   propdata->computedsymmetry = FALSE;

   propdata->usesymmetry = -1;
   propdata->symconsenabled = FALSE;
   propdata->triedaddconss = FALSE;
   propdata->genconss = NULL;
   propdata->ngenconss = 0;
   propdata->nsymresacks = 0;
   propdata->norbitopes = 0;

   propdata->ofenabled = FALSE;
   propdata->bg0 = NULL;
   propdata->bg0list = NULL;
   propdata->nbg0 = 0;
   propdata->bg1 = NULL;
   propdata->bg1list = NULL;
   propdata->nbg1 = 0;
   propdata->permvarsevents = NULL;
   propdata->inactiveperms = NULL;
   propdata->nmovedpermvars = 0;
   propdata->lastrestart = 0;
   propdata->nfixedzero = 0;
   propdata->nfixedone = 0;
   propdata->nodenumber = -1;
   propdata->offoundreduction = FALSE;

   /* create event handler */
   propdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(propdata->eventhdlr), EVENTHDLR_SYMMETRY_NAME, EVENTHDLR_SYMMETRY_DESC,
         eventExecSymmetry, NULL) );
   assert( propdata->eventhdlr != NULL );

   /* include constraint handler */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC,
         PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, propExecSymmetry, propdata) );
   assert( prop != NULL );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSymmetry) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSymmetry) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreSymmetry) );
   SCIP_CALL( SCIPsetPropExitpre(scip, prop, propExitpreSymmetry) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropSymmetry) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolSymmetry, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );

   /* include table */
   SCIP_CALL( SCIPallocBlockMemory(scip, &tabledata) );
   tabledata->propdata = propdata;
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_ORBITALFIXING, TABLE_DESC_ORBITALFIXING, TRUE,
         NULL, tableFreeOrbitalfixing, NULL, NULL, NULL, NULL, tableOutputOrbitalfixing,
         tabledata, TABLE_POSITION_ORBITALFIXING, TABLE_EARLIEST_ORBITALFIXING) );

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
         "propagating/" PROP_NAME "/detectorbitopes",
         "Should we check whether the components of the symmetry group can be handled by orbitopes?",
         &propdata->detectorbitopes, TRUE, DEFAULT_DETECTORBITOPES, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/addconsstiming",
         "timing of adding constraints (0 = before presolving, 1 = during presolving, 2 = after presolving)",
         &propdata->addconsstiming, TRUE, DEFAULT_ADDCONSSTIMING, 0, 2, NULL, NULL) );

   /* add parameters for orbital fixing */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/ofsymcomptiming",
         "timing of symmetry computation for orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call)",
         &propdata->ofsymcomptiming, TRUE, DEFAULT_OFSYMCOMPTIMING, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/performpresolving",
         "run orbital fixing during presolving?",
         &propdata->performpresolving, TRUE, DEFAULT_PERFORMPRESOLVING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/recomputerestart",
         "recompute symmetries after a restart has occured?",
         &propdata->recomputerestart, TRUE, DEFAULT_RECOMPUTERESTART, NULL, NULL) );

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

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/disableofrestart",
         "Shall orbital fixing be disabled if orbital fixing has found a reduction and a restart occurs?",
         &propdata->disableofrestart, TRUE, DEFAULT_DISABLEOFRESTART, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/symfixnonbinaryvars",
         "Whether all non-binary variables shall be not affected by symmetries if OF is active?",
         &propdata->symfixnonbinaryvars, TRUE, DEFAULT_SYMFIXNONBINARYVARS, NULL, NULL) );

   /* possibly add description */
   if ( SYMcanComputeSymmetry() )
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetName(), SYMsymmetryGetDesc()) );
   }

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
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether binary variables are affected */
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
      *permvarmap = propdata->permvarmap;

   *nperms = propdata->nperms;
   if ( perms != NULL )
   {
      *perms = propdata->perms;
      assert( *perms != NULL || *nperms <= 0 );
   }

   if ( permstrans != NULL )
   {
      *permstrans = propdata->permstrans;
      assert( *permstrans != NULL || *nperms <= 0 );
   }

   if ( log10groupsize != NULL )
      *log10groupsize = propdata->log10groupsize;

   if ( binvaraffected != NULL )
      *binvaraffected = propdata->binvaraffected;

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

/** return whether orbital fixing is enabled */
SCIP_Bool SCIPisOrbitalfixingEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );

   prop = SCIPfindProp(scip, PROP_NAME);
   if ( prop == NULL )
      return FALSE;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   return propdata->ofenabled;
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
