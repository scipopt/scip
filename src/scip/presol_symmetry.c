/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_symmetry.c
 * @brief  presolver for storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 *
 * This presolver computes symmetries of the problem and stores this information in adequate form. It does not
 * perform additional actions. The symmetry information can be accessed through external functions. However, the user
 * has to declare the type of symmetry that is needed before execution, see SYMsetSpecRequirement().
 *
 * @note We treat implict integer variables as if they were continuous/real variables. The reason is that there is
 * currently no distinction between implicit integer and implicit binary. Moreover, currently implicit integer variables
 * hurt our code more than continuous/real variables (we basically do not handle integral variables at all).
 *
 * @note We do not copy symmetry information, since it is not clear how this information transfers. Moreover, copying
 * symmetry might inhibit heuristics. But note that solving the a sub-SCIP might then happen without symmetry
 * information!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
#include <scip/cons_knapsack.h>
#include <scip/cons_varbound.h>
#include <scip/cons_setppc.h>
#include <scip/cons_and.h>
#include <scip/cons_logicor.h>
#include <scip/cons_or.h>
#include <scip/cons_xor.h>
#include <scip/cons_linking.h>
#include <scip/misc.h>

#include <scip/presol_symmetry.h>
#include <symmetry/compute_symmetry.h>

#include <string.h>

/* presolver properties */
#define PRESOL_NAME            "symmetry"
#define PRESOL_DESC            "presolver for computing and storing symmetry information about current problem"
#define PRESOL_PRIORITY               0      /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1      /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING            SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

/* default parameter values */
#define DEFAULT_MAXGENERATORS      1500      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
#define DEFAULT_CHECKSYMMETRIES   FALSE      /**< Should all symmetries be checked after computation? */
#define DEFAULT_DISPLAYNORBITVARS FALSE      /**< Should the number of variables affected by some symmetry be displayed? */

/* event handler properties */
#define EVENTHDLR_SYMMETRY_NAME    "symmetry"
#define EVENTHDLR_SYMMETRY_DESC    "filter global variable fixing event handler for orbital fixing"

/* other defines */
#define MAXGENNUMERATOR        64000000      /**< determine maximal number of generators by dividing this number by the number of variables */

#define SCIP_OUTPUT               FALSE
#define SCIP_OUTPUT_COMPONENT     FALSE

/* macros for getting activeness of symmetry handling methods */
#define ISSYMRETOPESACTIVE(x)      ((x & SYM_HANDLETYPE_SYMBREAK) != 0)
#define ISORBITALFIXINGACTIVE(x)   ((x & SYM_HANDLETYPE_ORBITALFIXING) != 0)


/** presolver data */
struct SCIP_PresolData
{
   int                   maxgenerators;      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
   SCIP_Bool             checksymmetries;    /**< Should all symmetries be checked after computation? */
   SCIP_Bool             displaynorbitvars;  /**< Whether the number of variables in non-trivial orbits shall be computed */
   int                   npermvars;          /**< number of variables for permutations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   SCIP_Real*            permvarsobj;        /**< objective values of permuted variables (for debugging) */
   int                   nperms;             /**< number of permutations */
   int                   nmaxperms;          /**< maximal number of permutations (needed for freeing storage) */
   int**                 perms;              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   int**                 permstrans;         /**< pointer to store transposed permutation generators as (npermvars x nperms) matrix */
   SCIP_Real             log10groupsize;     /**< log10 of size of symmetry group */
   int                   norbitvars;         /**< number of vars that are contained in a non-trivial orbit */
   SCIP_Bool             binvaraffected;     /**< whether binary variables are affected by some symmetry */
   SCIP_Bool             computedsym;        /**< Have we already tried to compute symmetries? */
   SCIP_Bool             successful;         /**< Was the computation of symmetries successful? */
   int                   usesymmetry;        /**< encoding of active symmetry handling methods  */

   /* components of symmetry group */
   int                   ncomponents;        /**< number of components of symmetry group */
   int*                  components;         /**< array containing the indices of permutations sorted by components */
   int*                  componentbegins;    /**< array containing in i-th position the first position of
                                              *   component i in components array */
   int*                  vartocomponent;     /**< array containing for each permvar the index of the component it is
                                              *   contained in (-1 if not affected) */
   SCIP_Shortbool*       componentblocked;   /**< array to store whether a component is blocked to be considered by
                                              *   further symmetry handling techniques */

   /* data necessary for orbital fixing */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for handling global variable fixings */
   SCIP_Shortbool*       bg0;                /**< bitset to store variables globally fixed to 0 */
   int*                  bg0list;            /**< list of variables globally fixed to 0 */
   int                   nbg0;               /**< number of variables in bg0 and bg0list */
   SCIP_Shortbool*       bg1;                /**< bitset to store variables globally fixed or branched to 1 */
   int*                  bg1list;            /**< list of variables globally fixed or branched to 1 */
   int                   nbg1;               /**< number of variables in bg1 and bg1list */
   SCIP_HASHMAP*         permvarmap;         /**< map of variables to indices in permvars array */
   int*                  permvarsevents;     /**< stores events caught for permvars */
};


/*
 * local data structures
 */

/* ------------------- map for variable types ------------------- */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyVartype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare the types of two variables according to objective, lower and upper bound, and variable type.
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

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValVartype)
{  /*lint --e{715}*/
   SYM_VARTYPE* k;

   k = (SYM_VARTYPE*) key;

   return SCIPhashTwo(SCIPcombineTwoInt(SCIPrealHashCode(k->obj), SCIPrealHashCode(k->lb)), SCIPrealHashCode(k->ub));
}


/* ------------------- sorting function for rhs types ------------------- */

/** data struct to store arrays used for sorting rhs types */
struct SYM_Sortrhstype
{
   SCIP_Real*            vals;               /**< array of values */
   SYM_RHSSENSE*         senses;             /**< array of senses of rhs */
   int                   nrhscoef;           /**< size of arrays (for debugging) */
};
typedef struct SYM_Sortrhstype SYM_SORTRHSTYPE;

/** sort rhs types - first by sense, then by value
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

/** sort matrix coefficients
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


/** fill in matrix elements into coefficient arrays */
static
SCIP_RETCODE collectCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficients values (or NULL if all linear coefficient values are 1) */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             istransformed,      /**< whether the constraint is transformed */
   SYM_RHSSENSE          rhssense,           /**< identifier of constraint type */
   SYM_MATRIXDATA*       matrixdata          /**< matrix data to be filled in */
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
      assert( ! SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      matrixdata->rhscoef[nrhscoef] = rhs;
      /* if we deal with special constraints */
      if ( (int) rhssense >= 3 )
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

         matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);
         matrixdata->matcoef[nmatcoef++] = vals[j];
      }
      nrhscoef++;
   }
   else
   {
      if ( ! SCIPisInfinity(scip, -lhs) )
      {
         matrixdata->rhscoef[nrhscoef] = -lhs;
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_INEQUALITY;
         matrixdata->rhsidx[nrhscoef] = nrhscoef;

         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;
            matrixdata->matvaridx[nmatcoef] = SCIPvarGetProbindex(vars[j]);

            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

            matrixdata->matcoef[nmatcoef++] = -vals[j];
         }
         nrhscoef++;
      }

      if ( ! SCIPisInfinity(scip, rhs) )
      {
         matrixdata->rhscoef[nrhscoef] = rhs;
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_INEQUALITY;
         matrixdata->rhsidx[nrhscoef] = nrhscoef;

         for (j = 0; j < nvars; ++j)
         {
            assert( nmatcoef < matrixdata->nmaxmatcoef );
            matrixdata->matidx[nmatcoef] = nmatcoef;
            matrixdata->matrhsidx[nmatcoef] = nrhscoef;

            assert( 0 <= SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) < SCIPgetNVars(scip) );

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

   /* set up dense arrow for permuted row */
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

/** compute symmetry group of MIP */
static
SCIP_RETCODE computeSymmetryGroup(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_SPEC              fixedtype,          /**< variable types that must be fixed by symmetries */
   SCIP_Bool             local,              /**< Use local variable bounds? */
   SCIP_Bool             checksymmetries,    /**< Should all symmetries be checked after computation? */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_Real**           permvarsobj,        /**< objective values of permuted variables */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   int***                permstrans,         /**< pointer to store permutation generators as (npermvars x nperms) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   int                   usesymmetry,        /**< identifier of active symmetry handling routines */
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
   int nuniquevararray = 0;
   int nhandleconss;
   int nactiveconss;
   int nconss;
   int nvars;
   int nallvars;
   int c;
   int j;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( permvarsobj != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( success != NULL );

   /* init */
   *npermvars = 0;
   *permvars = NULL;
   *permvarsobj = NULL;
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *permstrans = NULL;
   *log10groupsize = 0;
   *success = FALSE;

   /* skip if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
      return SCIP_OKAY;

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

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
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "Detecting %ssymmetry on %d variables and %d constraints.\n", local ? "local " : "", nvars, nactiveconss);

   /* copy variables */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, SCIPgetVars(scip), nvars) ); /*lint !e666*/
   assert( vars != NULL );

   /* fill matrixdata */
   matrixdata.nmaxmatcoef = 100 * nvars;
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
         SCIP_CALL( collectCoefficients(scip, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
               SCIPgetNVarsLinear(scip, cons), SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons),
               SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "linking") == 0 )
      {
         SCIP_VAR** curconsvars;
         int* curconsvals;
         int i;

         /* get constraint variables and their amount */
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
         consvars[nconsvars - 1] = SCIPgetIntvarLinking(scip, cons);
         consvals[nconsvars - 1] = -1;

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, nconsvars, 0.0, 0.0,
                        SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata) );
         SCIP_CALL( collectCoefficients(scip, consvars, NULL, nconsvars - 1, 1.0, 1.0,
                        SCIPconsIsTransformed(cons), SYM_SENSE_UNKOWN, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "setppc") == 0 )
      {
         linvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( collectCoefficients(scip, linvars, 0, nconsvars, 1.0, 1.0, SCIPconsIsTransformed(cons), SYM_SENSE_EQUATION, &matrixdata) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( collectCoefficients(scip, linvars, 0, nconsvars, -SCIPinfinity(scip), 1.0, SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( collectCoefficients(scip, linvars, 0, nconsvars, 1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata) );
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

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, nconsvars, (SCIP_Real) SCIPgetRhsXor(scip, cons),
               (SCIP_Real) SCIPgetRhsXor(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_XOR, &matrixdata) );
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

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_AND, &matrixdata) );
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

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, nconsvars, 0.0, 0.0,
               SCIPconsIsTransformed(cons), SYM_SENSE_OR, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( collectCoefficients(scip, SCIPgetVarsLogicor(scip, cons), 0, SCIPgetNVarsLogicor(scip, cons),
               1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata) );
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

         SCIP_CALL( collectCoefficients(scip, SCIPgetVarsKnapsack(scip, cons), consvals, nconsvars, -SCIPinfinity(scip),
               (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "varbound") == 0 )
      {
         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvals[0] = 1.0;

         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
               SCIPgetRhsVarbound(scip, cons), SCIPconsIsTransformed(cons), SYM_SENSE_INEQUALITY, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "bounddisjunction") == 0 )
      {
         /* currently assume bound disjunctions are o.k. for non local symmetry groups */
         if ( local )
         {
            /* @todo we need to handle bounddisjunctions if local symmetry groups are considered */
            SCIPerrorMessage("Cannot determine symmetries for constraint <%s> of constraint handler <%s>.\n",
               SCIPconsGetName(cons), SCIPconshdlrGetName(conshdlr) );
            return SCIP_ERROR;
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
   assert( matrixdata.nrhscoef > 0 ); /* cannot have empty rows! */

   SCIPfreeBlockMemoryArray(scip, &consvals, nallvars);
   SCIPfreeBlockMemoryArray(scip, &consvars, nallvars);

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

      /* if the variable type should be fixed just increase the color */
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

      /* updata data if nontrivial symmetry */
      if ( *nperms > 0 )
      {
         /* transpose symmetries matrix here if necessary */
         if ( ISORBITALFIXINGACTIVE(usesymmetry) )
         {
            int p;

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, permstrans, nvars) );
            for (j = 0; j < nvars; ++j)
            {
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*permstrans)[j], *nmaxperms) );
               for (p = 0; p < *nperms; ++p)
                  (*permstrans)[j][p] = (*perms)[p][j];
            }

            if ( ! ISSYMRETOPESACTIVE(usesymmetry) )
            {
               /* free original perms matrix */
               for (p = 0; p < *nperms; ++p)
               {
                  SCIPfreeBlockMemoryArray(scip, &(*perms)[p], nvars);
               }
               SCIPfreeBlockMemoryArrayNull(scip, perms, *nmaxperms);
               *perms = NULL;
            }
         }

         /* symmetric variables are not allowed to be multi-aggregated */
         for (j = 0; j < nvars; ++j)
         {
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, vars[j]) );
         }

#ifndef NDEBUG
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, permvarsobj, nvars) );
         for (j = 0; j < nvars; ++j)
            (*permvarsobj)[j] = SCIPvarGetObj(vars[j]);
#endif
      }
   }
   *success = TRUE;

   if ( *nperms > 0 )
   {
      /* copy variables */
      *permvars = vars;
      *npermvars = nvars;
   }
   else
   {
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
   }

   /* free matrix data */
   SCIPfreeBlockMemoryArray(scip, &uniquevararray, nvars);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoefcolors, matrixdata.nrhscoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoefcolors, matrixdata.nmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.permvarcolors, nvars);
   SCIPhashtableFree(&vartypemap);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nactiveconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

   return SCIP_OKAY;
}


/** compute components of symmetry group */
static
SCIP_RETCODE computeComponents(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   SCIP_DISJOINTSET* componentstovar = NULL;
   int** perms;
   int** permstrans;
   int* permtovarcomp;
   int* permtocomponent;
   int nperms;
   int npermvars;
   int ncomponents;
   int p;
   int i;
   int idx;

   assert( scip != NULL );
   assert( presoldata != NULL );

   assert( presoldata->ncomponents == -1 );
   assert( presoldata->components == NULL );
   assert( presoldata->componentbegins == NULL );
   assert( presoldata->vartocomponent == NULL );
   assert( presoldata->componentblocked == NULL );

#if SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation started\n", SCIPgetSolvingTime(scip));
#endif

   /* get data */
   nperms = presoldata->nperms;

   if ( nperms <= 0 )
      return SCIP_OKAY;

   npermvars = presoldata->npermvars;
   perms = presoldata->perms;
   permstrans = presoldata->permstrans;
   assert( npermvars > 0 );
   assert( (! ISORBITALFIXINGACTIVE(presoldata->usesymmetry) && perms != NULL)
      || (ISORBITALFIXINGACTIVE(presoldata->usesymmetry) && permstrans != NULL) );

   SCIP_CALL( SCIPdisjointsetCreate(&componentstovar, SCIPblkmem(scip), npermvars) );
   ncomponents = npermvars;

   /* init array that stores for each permutation the representative of its affected variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &permtovarcomp, nperms) );
   for (p = 0; p < nperms; ++p)
      permtovarcomp[p] = -1;

   /* find permutation components and store for each variable an affecting permutation (or -1)  */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->vartocomponent, npermvars) );
   for (i = 0; i < npermvars; ++i)
   {
      presoldata->vartocomponent[i] = -1;

      for (p = 0; p < nperms; ++p)
      {
         int img;

         img = ISORBITALFIXINGACTIVE(presoldata->usesymmetry) ? permstrans[i][p] : perms[p][i]; /*lint !e613*/

         /* perm p affects i -> possibly merge var components */
         if ( img != i )
         {
            int component1;
            int component2;
            int representative;

            component1 = SCIPdisjointsetFind(componentstovar, i);
            component2 = SCIPdisjointsetFind(componentstovar, img);
            presoldata->vartocomponent[i] = p;
            presoldata->vartocomponent[img] = p;

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
               --ncomponents;
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
               --ncomponents;
            }
            else if ( representative > component1 )
            {
               assert( representative == component2 );
               permtovarcomp[p] = component1;
            }
         }
      }

      /* reduce number of components by singletons */
      if ( presoldata->vartocomponent[i] == -1 )
         --ncomponents;
      else if ( SCIPvarIsBinary(presoldata->permvars[i]) )
         presoldata->binvaraffected = TRUE;
   }
   assert( ncomponents > 0 );
   presoldata->ncomponents = ncomponents;

   /* update permvartocomp array to final variable representatives */
   for (p = 0; p < nperms; ++p)
      permtovarcomp[p] = SCIPdisjointsetFind(componentstovar, permtovarcomp[p]);

   /* init components array by trivial natural order of permutations */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->components, nperms) );
   for (p = 0; p < nperms; ++p)
      presoldata->components[p] = p;

   /* get correct order of components array */
   SCIPsortIntInt(permtovarcomp, presoldata->components, nperms);

   /* determine componentbegins and store components for each permutation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->componentbegins, ncomponents + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permtocomponent, nperms) );

   presoldata->componentbegins[0] = 0;
   permtocomponent[presoldata->components[0]] = 0;
   idx = 0;

   for (p = 1; p < nperms; ++p)
   {
      if ( permtovarcomp[p] > permtovarcomp[p - 1] )
         presoldata->componentbegins[++idx] = p;

      assert( presoldata->components[p] >= 0 );
      assert( presoldata->components[p] < nperms );
      permtocomponent[presoldata->components[p]] = idx;
   }
   assert( ncomponents == idx + 1 );
   presoldata->componentbegins[++idx] = nperms;

   /* determine vartocomponent */
   for (i = 0; i < npermvars; ++i)
   {
      int permidx;
      permidx = presoldata->vartocomponent[i];
      assert( -1 <= permidx && permidx < nperms );

      if ( permidx != -1 )
      {
         assert( 0 <= permtocomponent[permidx] );
         assert( permtocomponent[permidx] < ncomponents );

         presoldata->vartocomponent[i] = permtocomponent[permidx];
      }
   }

   /* init componentblocked */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->componentblocked, ncomponents) );
   for (i = 0; i < ncomponents; ++i)
      presoldata->componentblocked[i] = FALSE;

   SCIPfreeBufferArray(scip, &permtocomponent);
   SCIPfreeBufferArray(scip, &permtovarcomp);
   SCIPdisjointsetFree(&componentstovar, SCIPblkmem(scip));

#if SCIP_OUTPUT_COMPONENT
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) component computation finished\n", SCIPgetSolvingTime(scip));
#endif

#if SCIP_OUTPUT
   printf("number of components: %d\n", presoldata->ncomponents);
   for (i = 0; i < ncomponents; ++i)
   {
      printf("Component %d contains the following permutations:\n\t", i);
      for (p = presoldata->componentbegins[i]; p < presoldata->componentbegins[i + 1]; ++p)
      {
         printf("%d, ", presoldata->components[p]);
      }
      printf("\n");
   }
#endif

   return SCIP_OKAY;
}


/** determine whether binary variable is effected (and potentially compute number of affected variables) */
static
SCIP_RETCODE determineBinvarAffected(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_Bool             completestatistic   /**< whether number of affected vars should be computed */
   )
{
   int** perms;
   int nperms;
   int nvars;
   SCIP_Shortbool* affected;
   int i;
   int p;
   int naffected = 0;

   assert( scip != NULL );
   assert( presoldata != NULL );

   if ( presoldata->binvaraffected && !completestatistic )
      return SCIP_OKAY;

   assert( presoldata->perms != NULL );
   assert( presoldata->nperms > 0 );
   assert( presoldata->npermvars > 0 );

   perms = presoldata->perms;
   nperms = presoldata->nperms;
   nvars = presoldata->npermvars;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &affected, nvars) );

   /* iterate over permutations and check which variables are affected by some symmetry */
   for (p = 0; p < nperms && (completestatistic || ! presoldata->binvaraffected); ++p)
   {
      for (i = 0; i < nvars; ++i)
      {
         if ( affected[i] )
            continue;

         if ( perms[p][i] != i )
         {
            if ( SCIPvarIsBinary(presoldata->permvars[i]) )
            {
               presoldata->binvaraffected = TRUE;

               if ( ! completestatistic )
                  break;
            }

            affected[i] = TRUE;
            ++naffected;
         }
      }
   }

   if ( completestatistic )
      presoldata->norbitvars = naffected;

   SCIPfreeBufferArray(scip, &affected);

   return SCIP_OKAY;
}


/** determine symmetry */
static
SCIP_RETCODE determineSymmetry(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SYM_SPEC              symspecrequire,     /**< symmetry specification for which we need to compute symmetries */
   SYM_SPEC              symspecrequirefixed /**< symmetry specification of variables which must be fixed by symmetries */
   )
{
   int maxgenerators;
   int type = 0;
   int nvars;

   assert( scip != NULL );
   assert( presoldata != NULL );

   assert( ! presoldata->computedsym );
   assert( presoldata->npermvars == 0 );
   assert( presoldata->permvars == NULL );
   assert( presoldata->permvarsobj == NULL );
   assert( presoldata->nperms == 0 );
   assert( presoldata->nmaxperms == 0 );
   assert( presoldata->perms == NULL );

   presoldata->computedsym = TRUE;

#ifndef NDEBUG
   {
      int usesymmetry;
      SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &usesymmetry) );
      assert( usesymmetry );
   }
#endif

   /* do not compute symmetry if there are active pricers */
   if ( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* avoid trivial cases */
   nvars = SCIPgetNVars(scip);
   if ( nvars <= 0 )
      return SCIP_OKAY;

   /* determine symmetry specification */
   if ( SCIPgetNBinVars(scip) > 0 )
      type |= (int) SYM_SPEC_BINARY;
   if ( SCIPgetNIntVars(scip) > 0 )
      type |= (int) SYM_SPEC_INTEGER;
   /* count implicit integer variables as real variables, since we cannot currently handle integral variables well */
   if ( SCIPgetNContVars(scip) > 0 || SCIPgetNImplVars(scip) > 0 )
      type |= (int) SYM_SPEC_REAL;

   /* skip symmetry computation if no graph automorphism code was linked */
   if ( ! SYMcanComputeSymmetry() )
   {
      int nconss = SCIPgetNActiveConss(scip);
      int nhandleconss = getNSymhandableConss(scip);

      /* print verbMessage only if problem consists of symmetry handable constraints */
      assert( nhandleconss <=  nconss );
      if ( nhandleconss < nconss )
         return SCIP_OKAY;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   Deactivated symmetry handling methods, since SCIP was built without symmetry detector (SYM=none).\n");
      return SCIP_OKAY;
   }
   /* skip symmetry computation if required variables are not present */
   else if ( ! (type & symspecrequire) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation skipped: type (bin %c, int %c, cont %c) does not match requirements (bin %c, int %c, cont %c)\n",
         SCIPgetSolvingTime(scip),
         SCIPgetNBinVars(scip) > 0 ? '+' : '-',
         SCIPgetNIntVars(scip) > 0  ? '+' : '-',
         SCIPgetNContVars(scip) + SCIPgetNImplVars(scip) > 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
         (symspecrequire & (int) SYM_SPEC_REAL) != 0 ? '+' : '-');
      return SCIP_OKAY;
   }
   /* skip symmetry computation if there are constraints that cannot be handled by symmetry */
   else if ( getNSymhandableConss(scip) < SCIPgetNActiveConss(scip) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation skipped: there exist constraints that cannot be handled by symmetry methods\n",
         SCIPgetSolvingTime(scip));
      return SCIP_OKAY;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   (%.1fs) symmetry computation started: requiring (bin %c, int %c, cont %c), (fixed: bin %c, int %c, cont %c)\n",
      SCIPgetSolvingTime(scip),
      (symspecrequire & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
      (symspecrequire & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
      (symspecrequire & (int) SYM_SPEC_REAL) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_BINARY) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_INTEGER) != 0 ? '+' : '-',
      (symspecrequirefixed & (int) SYM_SPEC_REAL) != 0 ? '+' : '-');

   if ( symspecrequire & symspecrequirefixed )
      SCIPwarningMessage(scip, "Warning: some required symmetries must be fixed.\n");

   /* actually compute (global) symmetry */
   /* determine maximal number of generators depending on the number of variables */
   maxgenerators = presoldata->maxgenerators;
   maxgenerators = MIN(maxgenerators, MAXGENNUMERATOR / nvars);

   SCIP_CALL( computeSymmetryGroup(scip, maxgenerators, symspecrequirefixed, FALSE, presoldata->checksymmetries,
         &presoldata->npermvars, &presoldata->permvars, &presoldata->permvarsobj, &presoldata->nperms,
         &presoldata->nmaxperms, &presoldata->perms, &presoldata->permstrans,
         &presoldata->log10groupsize, presoldata->usesymmetry, &presoldata->successful) );

   /* output statistics */
   if ( ! presoldata->successful )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) could not compute symmetry\n", SCIPgetSolvingTime(scip));
   else if ( presoldata->nperms == 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry present\n", SCIPgetSolvingTime(scip));
   else
   {
      assert( presoldata->nperms > 0 );

      if ( presoldata->displaynorbitvars )
      {
         SCIP_CALL( determineBinvarAffected(scip, presoldata, TRUE) );
      }
      else if ( ISSYMRETOPESACTIVE(presoldata->usesymmetry) )
      {
         SCIP_CALL( determineBinvarAffected(scip, presoldata, FALSE) );
      }

      /* display statistics: number of generators */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "   (%.1fs) symmetry computation finished: %d generators found (max: ",
         SCIPgetSolvingTime(scip), presoldata->nperms);

      /* display statistics: maximum number of generators*/
      if ( maxgenerators == 0 )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "-");
      else
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%u", maxgenerators);

      /* display statistics: log10 group size, number of affected vars*/
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", log10 of symmetry group size: %.1f", presoldata->log10groupsize);

      /* display statistics: number of affected vars*/
      if ( presoldata->displaynorbitvars )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ", number of affected variables: %d)\n", presoldata->norbitvars);
      else
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, ")\n");

      /* do not deactivate components if no binary variables are affected in the polyhedral setting */
      if ( ! presoldata->binvaraffected && presoldata->usesymmetry == 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) no symmetry on binary variables present\n", SCIPgetSolvingTime(scip));

         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}


/*
 * Event handler callback methods
 */

/** exec the event handler for handling global variable lower bound changes (necessary for orbital fixing)
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
   SCIP_PRESOLDATA* presoldata;
   SCIP_VAR* var;
   int varidx;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_SYMMETRY_NAME) == 0 );
   assert( event != NULL );

   presoldata = (SCIP_PRESOLDATA*) eventdata;
   assert( presoldata != NULL );
   assert( presoldata->permvarmap != NULL );
   assert( presoldata->permstrans != NULL );
   assert( presoldata->nperms > 0 );
   assert( presoldata->permvars != NULL );
   assert( presoldata->npermvars > 0 );

   /* get fixed variable */
   var = SCIPeventGetVar(event);
   assert( var != NULL );
   assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );

   if ( ! SCIPhashmapExists(presoldata->permvarmap, (void*) var) )
   {
      SCIPerrorMessage("Invalid variable.\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }
   varidx = SCIPhashmapGetImageInt(presoldata->permvarmap, (void*) var);
   assert( 0 <= varidx && varidx < presoldata->npermvars );

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_GUBCHANGED )
   {
      assert( SCIPisEQ(scip, SCIPeventGetNewbound(event), 0.0) );
      assert( SCIPisEQ(scip, SCIPeventGetOldbound(event), 1.0) );

      SCIPdebugMsg(scip, "Mark variable <%s> as globally fixed to 0.\n", SCIPvarGetName(var));
      assert( ! presoldata->bg0[varidx] );
      presoldata->bg0[varidx] = TRUE;
      presoldata->bg0list[presoldata->nbg0++] = varidx;
      assert( presoldata->nbg0 <= presoldata->npermvars );
   }

   if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_GLBCHANGED )
   {
      assert( SCIPisEQ(scip, SCIPeventGetNewbound(event), 1.0) );
      assert( SCIPisEQ(scip, SCIPeventGetOldbound(event), 0.0) );

      SCIPdebugMsg(scip, "Mark variable <%s> as globally fixed to 1.\n", SCIPvarGetName(var));
      assert( ! presoldata->bg1[varidx] );
      presoldata->bg1[varidx] = TRUE;
      presoldata->bg1list[presoldata->nbg1++] = varidx;
      assert( presoldata->nbg1 <= presoldata->npermvars );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitSymmetry)
{
   SCIP_PRESOLDATA* presoldata;
   int i;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   SCIPdebugMsg(scip, "Exiting symmetry presolver.\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   if ( presoldata->ncomponents > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->componentblocked, presoldata->ncomponents);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->vartocomponent, presoldata->npermvars);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->componentbegins, presoldata->ncomponents + 1);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->components, presoldata->nperms);
   }

   if ( ISORBITALFIXINGACTIVE(presoldata->usesymmetry) )
   {
      int v;

      if ( presoldata->permvarmap != NULL )
      {
         SCIPhashmapFree(&presoldata->permvarmap);
      }

      /* free variables */
      for (v = 0; v < presoldata->npermvars; ++v)
      {
         if ( SCIPvarGetType(presoldata->permvars[v]) == SCIP_VARTYPE_BINARY && presoldata->permvarsevents[v] >= 0 )
         {
            /* If symmetry is computed before presolving, it might happen that some variables are turned into binary
             * variables, for which no event has been catched. Since there currently is no way of checking whether a var
             * event has been caught for a particular variable, we use the stored eventfilter positions. */
            SCIP_CALL( SCIPdropVarEvent(scip, presoldata->permvars[v], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                  presoldata->eventhdlr, (SCIP_EVENTDATA*) presoldata, presoldata->permvarsevents[v]) );
         }
         SCIP_CALL( SCIPreleaseVar(scip, &presoldata->permvars[v]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg0list, presoldata->npermvars);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg0, presoldata->npermvars);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg1list, presoldata->npermvars);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg1, presoldata->npermvars);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvarsevents, presoldata->npermvars);

      /* free permstrans matrix*/
      assert( presoldata->permstrans != NULL || presoldata->nperms == 0 );
      for (i = 0; i < presoldata->npermvars; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &presoldata->permstrans[i], presoldata->nmaxperms);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permstrans, presoldata->npermvars);
   }

   if ( ISSYMRETOPESACTIVE(presoldata->usesymmetry) )
   {
      assert( presoldata->perms != NULL || presoldata->nperms == 0 );
      for (i = 0; i < presoldata->nperms; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &presoldata->perms[i], presoldata->npermvars);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->perms, presoldata->nmaxperms);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvars, presoldata->npermvars);
   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvarsobj, presoldata->npermvars);

   /* reset settings */
   presoldata->npermvars = 0;
   presoldata->nperms = 0;
   presoldata->nmaxperms = 0;
   presoldata->norbitvars = 0;
   presoldata->binvaraffected = FALSE;
   presoldata->computedsym = FALSE;
   presoldata->successful = FALSE;
   presoldata->ncomponents = -1;
   presoldata->nbg0 = 0;
   presoldata->nbg1 = 0;
   presoldata->permvarmap = NULL;
   presoldata->permvarsevents = NULL;

   return SCIP_OKAY;
}

/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreSymmetry)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* check whether we should run */
   SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &presoldata->usesymmetry) );

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeSymmetry)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   SCIPdebugMsg(scip, "Freeing symmetry presolver.\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   SCIPfreeBlockMemory(scip, &presoldata);

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecSymmetry)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}


/*
 * External methods
 */

/** include symmetry constraint handler */
SCIP_RETCODE SCIPincludePresolSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol = NULL;
   SCIP_PRESOLDATA* presoldata = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );
   assert( presoldata != NULL );

   presoldata->npermvars = 0;
   presoldata->permvars = NULL;
   presoldata->permvarsobj = NULL;
   presoldata->perms = NULL;
   presoldata->permstrans = NULL;
   presoldata->nperms = 0;
   presoldata->nmaxperms = 0;
   presoldata->norbitvars = 0;
   presoldata->binvaraffected = FALSE;
   presoldata->computedsym = FALSE;
   presoldata->successful = FALSE;
   presoldata->ncomponents = -1;
   presoldata->components = NULL;
   presoldata->componentbegins = NULL;
   presoldata->vartocomponent = NULL;
   presoldata->componentblocked = NULL;
   presoldata->bg0 = NULL;
   presoldata->bg0list = NULL;
   presoldata->nbg0 = 0;
   presoldata->bg1 = NULL;
   presoldata->bg1list = NULL;
   presoldata->nbg1 = 0;
   presoldata->permvarmap = NULL;
   presoldata->permvarsevents = NULL;

   /* create event handler if orbital fixing is active */
   presoldata->eventhdlr = NULL;
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(presoldata->eventhdlr), EVENTHDLR_SYMMETRY_NAME, EVENTHDLR_SYMMETRY_DESC,
            eventExecSymmetry, NULL) );
      assert( presoldata->eventhdlr != NULL );

   /* include constraint handler */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC,
         PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING, presolExecSymmetry, presoldata) );
   assert( presol != NULL );

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSymmetry) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitSymmetry) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreSymmetry) );

   /* add parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/maxgenerators",
         "limit on the number of generators that should be produced within symmetry detection (0 = no limit)",
         &presoldata->maxgenerators, TRUE, DEFAULT_MAXGENERATORS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/checksymmetries",
         "Should all symmetries be checked after computation?",
         &presoldata->checksymmetries, TRUE, DEFAULT_CHECKSYMMETRIES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/displaynorbitvars",
         "Should the number of variables affected by some symmetry be displayed?",
         &presoldata->displaynorbitvars, TRUE, DEFAULT_DISPLAYNORBITVARS, NULL, NULL) );

   /* possibly add description */
   if ( SYMcanComputeSymmetry() )
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetName(), SYMsymmetryGetDesc()) );
   }

   return SCIP_OKAY;
}


/** return symmetry group generators */
SCIP_RETCODE SCIPgetGeneratorsSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SPEC              symspecrequire,     /**< symmetry specification for which we need to compute symmetries */
   SYM_SPEC              symspecrequirefixed,/**< symmetry specification of variables which must be fixed by symmetries */
   SCIP_Bool             recompute,          /**< Have symmetries already been computed? */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
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
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;
   SCIP_Bool computedsym;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( perms != NULL || permstrans != NULL );
   assert( ncomponents != NULL || (components == NULL && componentbegins == NULL && vartocomponent == NULL) );

   /* find symmetry presolver */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == NULL )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* free symmetry information if we recompute symmetries */
   if ( recompute )
   {
      int i;

      if ( presoldata->ncomponents > 0 )
      {
         assert( presoldata->componentblocked != NULL );
         assert( presoldata->vartocomponent != NULL );
         assert( presoldata->componentbegins != NULL );
         assert( presoldata->components != NULL );
         SCIPfreeBlockMemoryArray(scip, &presoldata->componentblocked, presoldata->ncomponents);
         SCIPfreeBlockMemoryArray(scip, &presoldata->vartocomponent, presoldata->npermvars);
         SCIPfreeBlockMemoryArray(scip, &presoldata->componentbegins, presoldata->ncomponents + 1);
         SCIPfreeBlockMemoryArray(scip, &presoldata->components, presoldata->nperms);

         presoldata->ncomponents = -1;
      }

      /* free data needed for orbital fixing */
      if ( ISORBITALFIXINGACTIVE(presoldata->usesymmetry) )
      {
         int v;

         SCIPhashmapFree(&presoldata->permvarmap);

         /* free variables */
         for (v = 0; v < presoldata->npermvars; ++v)
         {
            if ( SCIPvarGetType(presoldata->permvars[v]) == SCIP_VARTYPE_BINARY && presoldata->permvarsevents[v] >= 0 )
            {
               /* If symmetry is computed before presolving, it might happen that some variables are turned into binary
                * variables, for which no event has been catched. Since there currently is no way of checking whether a var
                * event has been caught for a particular variable, we use the stored eventfilter positions. */
               SCIP_CALL( SCIPdropVarEvent(scip, presoldata->permvars[v], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                     presoldata->eventhdlr, (SCIP_EVENTDATA*) presoldata, presoldata->permvarsevents[v]) );
            }
            SCIP_CALL( SCIPreleaseVar(scip, &presoldata->permvars[v]) );
         }
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg0list, presoldata->npermvars);
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg0, presoldata->npermvars);
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg1list, presoldata->npermvars);
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->bg1, presoldata->npermvars);
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvarsevents, presoldata->npermvars);

         assert( presoldata->permstrans != NULL );
         for (i = 0; i < presoldata->npermvars; ++i)
         {
            SCIPfreeBlockMemoryArray(scip, &presoldata->permstrans[i], presoldata->nmaxperms);
         }
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permstrans, presoldata->npermvars);
      }

      /* free data needed for symretopes */
      if ( ISSYMRETOPESACTIVE(presoldata->usesymmetry) )
      {
         assert( presoldata->perms != NULL );
         for (i = 0; i < presoldata->nperms; ++i)
         {
            SCIPfreeBlockMemoryArray(scip, &presoldata->perms[i], presoldata->npermvars);
         }
         SCIPfreeBlockMemoryArrayNull(scip, &presoldata->perms, presoldata->nmaxperms);
      }

      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvars, presoldata->npermvars);
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvarsobj, presoldata->npermvars);

      /* reset settings */
      presoldata->npermvars = 0;
      presoldata->nperms = 0;
      presoldata->nmaxperms = 0;
      presoldata->norbitvars = 0;
      presoldata->binvaraffected = FALSE;
      presoldata->computedsym = FALSE;
      presoldata->successful = FALSE;
      presoldata->ncomponents = -1;
      presoldata->nbg0 = 0;
      presoldata->nbg1 = 0;
      presoldata->permvarmap = NULL;
      presoldata->permvarsevents = NULL;
   }

   /* if not already done before, compute symmetries; store old value (might get manipulated by determineSymmetry()) */
   computedsym = presoldata->computedsym;
   if ( ! computedsym )
   {
      if ( SCIPgetStage(scip) != SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING &&
           SCIPgetStage(scip) != SCIP_STAGE_EXITPRESOLVE && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVED &&
           SCIPgetStage(scip) != SCIP_STAGE_INITSOLVE && SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      {
         SCIPerrorMessage("Cannot call symmetry detection outside of presolving.\n");
         return SCIP_INVALIDCALL;
      }

      /* determine symmetry here */
      SCIP_CALL( determineSymmetry(scip, presoldata, symspecrequire, symspecrequirefixed) );
   }

   *npermvars = presoldata->npermvars;
   *permvars = presoldata->permvars;
   *nperms = presoldata->nperms;
   if ( perms != NULL )
   {
      *perms = presoldata->perms;
      assert( *perms != NULL || *nperms == 0 );
   }
   if ( permstrans != NULL )
   {
      *permstrans = presoldata->permstrans;
      assert( *permstrans != NULL || *nperms == 0 );
   }

   if ( log10groupsize != NULL )
      *log10groupsize = presoldata->log10groupsize;
   if ( binvaraffected != NULL )
      *binvaraffected = presoldata->binvaraffected;

   if ( ncomponents != NULL || components != NULL || componentbegins != NULL || vartocomponent != NULL )
   {
      /* components might have been already computed if orbitopes and orbital fixing are both used */
      if ( presoldata->ncomponents == -1 )
      {
         SCIP_CALL( computeComponents(scip, presoldata) );
      }

      if ( components != NULL )
         *components = presoldata->components;

      if ( componentbegins != NULL )
         *componentbegins = presoldata->componentbegins;

      if ( vartocomponent )
         *vartocomponent = presoldata->vartocomponent;

      if ( ncomponents )
         *ncomponents = presoldata->ncomponents;
   }

   /* if not already done before, set data for event handler if orbital fixing is active */
   if ( ! computedsym && ISORBITALFIXINGACTIVE(presoldata->usesymmetry) )
   {
      int v;

      /* create hashmap for storing the indices of variables */
      assert( presoldata->permvarmap == NULL );
      SCIP_CALL( SCIPhashmapCreate(&presoldata->permvarmap, SCIPblkmem(scip), presoldata->npermvars) );

      /* insert variables into hashmap and capture variables */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->permvarsevents, presoldata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->bg0, presoldata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->bg0list, presoldata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->bg1, presoldata->npermvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &presoldata->bg1list, presoldata->npermvars) );

      for (v = 0; v < presoldata->npermvars; ++v)
      {
         SCIP_CALL( SCIPhashmapInsertInt(presoldata->permvarmap, (*permvars)[v], v) );
         SCIP_CALL( SCIPcaptureVar(scip, (*permvars)[v]) );

         presoldata->bg0[v] = FALSE;
         presoldata->bg1[v] = FALSE;
         presoldata->permvarsevents[v] = -1;

         /* only catch binary variables, since integer variables should be fixed pointwise;
          * implicit integer variables are not branched on
          */
         if ( SCIPvarGetType((*permvars)[v]) == SCIP_VARTYPE_BINARY )
         {
            /* catch whether lower bounds are changed, i.e., binary variables are fixed to 1;
             * also store filter position
             */
            SCIP_CALL( SCIPcatchVarEvent(scip, (*permvars)[v], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
                  presoldata->eventhdlr, (SCIP_EVENTDATA*) presoldata, &presoldata->permvarsevents[v]) );
         }
      }
      assert( presoldata->nbg1 == 0 );
   }

   return SCIP_OKAY;
}


/** return objective coefficients of permuted variables at time of symmetry computation */
SCIP_RETCODE SCIPgetPermvarsObjSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           permvarsobj         /**< pointer to store objective coefficients of permuted variables (NULL if not available) */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   assert( scip != NULL );
   assert( permvarsobj != NULL );

   /* find symmetry presolver */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == NULL )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   *permvarsobj = presoldata->permvarsobj;

   return SCIP_OKAY;
}


/** block component of symmetry group to be considered by symmetry handling routines */
SCIP_RETCODE SCIPsetSymmetryComponentblocked(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   i                   /**< index of component to block */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* find symmetry presolver */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == NULL )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );
   assert( 0 <= i && i < presoldata->ncomponents );
   assert( presoldata->componentblocked != NULL );

   presoldata->componentblocked[i] = TRUE;

   return SCIP_OKAY;
}


/** get blocked status component of symmetry group */
SCIP_Shortbool SCIPgetSymmetryComponentblocked(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   i                   /**< index of component to check blocked status */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* find symmetry presolver */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == NULL )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return FALSE;
   }
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );
   assert( 0 <= 0 && i < presoldata->ncomponents );
   assert( presoldata->componentblocked != NULL );

   return presoldata->componentblocked[i];
}


/** return symmetry information on globally fixed variables */
SCIP_RETCODE SCIPgetSyminfoGloballyFixedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Shortbool**      bg0,                /**< pointer to store array indicating whether var is globally fixed to 0 */
   int**                 bg0list,            /**< pointer to store list of vars globally fixed to 0 */
   int**                 nbg0,               /**< pointer to store memory position of number of vars globally fixed to 0 */
   SCIP_Shortbool**      bg1,                /**< pointer to store array indicating whether var is globally fixed to 1 */
   int**                 bg1list,            /**< pointer to store list of vars globally fixed to 1 */
   int**                 nbg1,               /**< pointer to store memory position of number of vars globally fixed to 1 */
   SCIP_HASHMAP**        permvarmap          /**< pointer to store hash map of permvars */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   assert( scip != NULL );
   assert( permvarmap != NULL );
   assert( bg0 != NULL );
   assert( bg0list != NULL );
   assert( bg1 != NULL );
   assert( bg1list != NULL );

   /* find symmetry presolver */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == NULL )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   *permvarmap = presoldata->permvarmap;
   *bg0 = presoldata->bg0;
   *bg0list = presoldata->bg0list;
   *nbg0 = &(presoldata->nbg0);
   *bg1 = presoldata->bg1;
   *bg1list = presoldata->bg1list;
   *nbg1 = &(presoldata->nbg1);

   return SCIP_OKAY;
}
