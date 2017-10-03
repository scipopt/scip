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

/**@file   presol_symmetry.cpp
 * @brief  presovler for storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 *
 * This presolver computes symmetries of the problem and stores this information in adequate form. It does not
 * perform additional actions. The symmetry information can be accessed through external functions. However, the user
 * has to declare the type of symmetry that is needed before execution, see SYMsetSpecRequirement().
 *
 * @note We treat implict integer variables as if they were continuous/real variables. The reason is that there is
 * currently no distinction between implicit integer and implicit binary. Moreover, currently implicit integer variables
 * hurt our code more than continuous/real variables (we basically do not handle integral variables at all).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
#include <scip/cons_knapsack.h>
#include <scip/cons_varbound.h>
#include <scip/cons_setppc.h>
#include <scip/cons_logicor.h>
#include <scip/cons_xor.h>

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
#define DEFAULT_DETECTSYMPRESOL    TRUE      /**< Should the symmetry be detected within presolving (otherwise before presol)? */

/* other defines */
#define MAXGENNUMERATOR        64000000      /**< determine maximal number of generators by dividing this number by the number of variables */


/** presolver data */
struct SCIP_PresolData
{
   SCIP_Bool             detectsympresol;    /**< Should the symmetry be detected within presolving (otherwise before presol)? */
   int                   maxgenerators;      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
   int                   symspecrequire;     /**< symmetry specification for which we need to compute symmetries */
   int                   symspecrequirefixed;/**< symmetry specification of variables which must be fixed by symmetries */
   int                   npermvars;          /**< number of variables for permutations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int                   nperms;             /**< number of permutations */
   int                   nmaxperms;          /**< maximal number of permutations (needed for freeing storage) */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   SCIP_Bool             computedsym;        /**< Have we already tried to compute symmetries? */
   SCIP_Bool             successful;         /**< Was the computation of symmetries successful? */
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
   return (uint64_t) (SCIPrealHashCode(k->obj) + SCIPrealHashCode(k->lb) + SCIPrealHashCode(k->ub));  /*lint !e776*/
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
   else
   {
      /* senses are equal, use values */
      diffvals = data->vals[ind1] - data->vals[ind2];

      if ( diffvals < 0.0 )
         return -1;
      else if ( diffvals > 0.0 )
         return 1;
   }
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
   int                   fixedtype,          /**< bitset of variable types that should be fixed */
   SCIP_VAR*             var                 /**< variable to be considered */
   )
{
   if ( (fixedtype & (int) SYM_SPEC_INTEGER) && SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
      return TRUE;
   if ( (fixedtype & (int) SYM_SPEC_BINARY) && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      return TRUE;
   if ( (fixedtype & (int) SYM_SPEC_REAL) &&
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

   assert( scip != 0 );
   assert( vars != 0 );
   assert( scalars != 0 );
   assert( *vars != 0 );
   assert( *scalars != 0 );
   assert( nvars != 0 );
   assert( constant != 0 );

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
   SCIP_Bool             isxor,              /**< whether the constraint is an XOR constraint */
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
   if ( linvals != 0 )
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

   /* check whether we have to resize */
   if ( matrixdata->nmatcoef + nvars > matrixdata->nmaxmatcoef )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, matrixdata->nmatcoef + nvars);
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
      if ( isxor )
         matrixdata->rhssense[nrhscoef] = SYM_SENSE_XOR;
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


#ifdef SCIP_DEBUG
/** checks whether given permutations form a symmetry of a MIP
 *
 *  We need the matrix and rhs in the original order in order to speed up the comparison process. The matrix is needed
 *  in the right order to easily check rows. The rhs is used because of cache effects.
 */
static
SCIP_RETCODE checkSymmetriesAreSymmetries(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   fixedtype,          /**< variable types that must be fixed by symmetries */
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
      assert( P != 0 );

      for (j = 0; j < matrixdata->npermvars; ++j)
      {
         if ( SymmetryFixVar(fixedtype, matrixdata->permvars[j]) && P[j] != j )
         {
            SCIPdebugMsg(scip, "Permutation does not fix types %d, moving variable %d.\n", fixedtype, j);
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
#endif


/** compute symmetry group of MIP */
static
SCIP_RETCODE computeSymmetryGroup(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   int                   fixedtype,          /**< variable types that must be fixed by symmetries */
   SCIP_Bool             local,              /**< Use local variable bounds? */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
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
   int nhandleconss = 0;
   int nactiveconss = 0;
   int nconss;
   int nvars;
   int c;
   int j;

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( success != NULL );

   /* init */
   *npermvars = 0;
   *permvars = NULL;
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *success = FALSE;

   /* skip if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
   {
      SCIPwarningMessage(scip, "Cannot compute symmetry, since no third party software has been linked in.\n");
      return SCIP_OKAY;
   }

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
   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );
      if ( SCIPconsIsActive(conss[c]) )
         ++nactiveconss;
   }

   /* exit if no active constraints are available */
   if ( nactiveconss == 0 )
   {
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* before we set up the matrix, check whether we can handle all constraints */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "setppc");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "xor");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "logicor");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "knapsack");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "bounddisjunction");
   nhandleconss += SCIPconshdlrGetNActiveConss(conshdlr);
   if ( nhandleconss < nactiveconss )
   {
      SCIPwarningMessage(scip, "Cannot compute symmetry, since unkown constraints are present.\n");
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

   /* prepare temporary constraint data (use block memory, since this can become large) */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvars, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consvals, nvars) );

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
               SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "setppc") == 0 )
      {
         linvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( collectCoefficients(scip, linvars, 0, nconsvars, 1.0, 1.0, SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( collectCoefficients(scip, linvars, 0, nconsvars, -SCIPinfinity(scip), 1.0, SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( collectCoefficients(scip, linvars, 0, nconsvars, 1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
            break;
         default:
            SCIPerrorMessage("Unknown setppc type %d.\n", SCIPgetTypeSetppc(scip, cons));
            return SCIP_ERROR;
         }
      }
      else if ( strcmp(conshdlrname, "xor") == 0 )
      {
         SCIP_Bool consvarssuccess;

         /* get number of variables of XOR constraint (should include integer variable) */
         SCIP_CALL( SCIPgetConsNVars(scip, cons, &nconsvars, &consvarssuccess) );
         assert( consvarssuccess );
         assert( nconsvars <= nvars );
         assert( nconsvars == SCIPgetNVarsXor(scip, cons) + 1 );

         /* get variables of XOR constraint */
         SCIP_CALL( SCIPgetConsVars(scip, cons, consvars, nconsvars, &consvarssuccess) );
         assert( consvarssuccess );

         for (j = 0; j < nconsvars; ++j)
         {
            /* mark integer variable with a coefficient of 2 to distinguish it from the other variables */
            if ( SCIPvarGetType(consvars[j]) == SCIP_VARTYPE_INTEGER )
               consvals[j] = 2.0;
            else
               consvals[j] = 1.0;
         }

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, nconsvars, (SCIP_Real) SCIPgetRhsXor(scip, cons),
               (SCIP_Real) SCIPgetRhsXor(scip, cons), SCIPconsIsTransformed(cons), TRUE, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( collectCoefficients(scip, SCIPgetVarsLogicor(scip, cons), 0, SCIPgetNVarsLogicor(scip, cons),
               1.0, SCIPinfinity(scip), SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         linvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);
         assert( nconsvars <= nvars );
         assert( consvals != NULL );

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         for (j = 0; j < nconsvars; ++j)
            consvals[j] = (SCIP_Real) weights[j];

         SCIP_CALL( collectCoefficients(scip, linvars, consvals, nconsvars, -SCIPinfinity(scip),
               (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "varbound") == 0 )
      {
         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( collectCoefficients(scip, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
               SCIPgetRhsVarbound(scip, cons), SCIPconsIsTransformed(cons), FALSE, &matrixdata) );
      }
      else if ( strcmp(conshdlrname, "bounddisjunction") == 0 )
      {
         /* currently assume bound disjunctions are o.k. for non local symmetry groups */
         if ( ! local )
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

   SCIPfreeBlockMemoryArray(scip, &consvals, nvars);
   SCIPfreeBlockMemoryArray(scip, &consvars, nvars);

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

   for (j = 0; j < matrixdata.nrhscoef; ++j)
      printf("%f ", matrixdata.rhscoef[matrixdata.rhsidx[j]]);
   printf("\n");

   /* determine number of different coefficents */

   /* find non-equivalent variables: same objective, lower and upper bounds, and variable type */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;

      var = vars[j];
      assert( var != 0 );

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
         SCIPdebugMsg(scip, "Detected new new rhs type %f, type: %u - color: %d\n", val, sense, matrixdata.nuniquerhs);
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
   assert( matrixdata.nuniquevars > 0 );
   assert( matrixdata.nuniquerhs > 0 );
   assert( matrixdata.nuniquemat > 0 );

   SCIPdebugMsg(scip, "Number of detected different variables: %d (total: %d).\n", matrixdata.nuniquevars, nvars);
   SCIPdebugMsg(scip, "Number of detected different rhs types: %d (total: %d).\n", matrixdata.nuniquerhs, matrixdata.nrhscoef);
   SCIPdebugMsg(scip, "Number of detected different matrix coefficients: %d (total: %d).\n", matrixdata.nuniquemat, matrixdata.nmatcoef);

   /* determine generators */
   SCIP_CALL( SYMcomputeSymmetryGenerators(scip, maxgenerators, &matrixdata, nperms, nmaxperms, perms) );

#ifdef SCIP_DEBUG
   if ( ! SCIPisStopped(scip) )
   {
      SCIP_CALL( checkSymmetriesAreSymmetries(scip, fixedtype, &matrixdata, *nperms, *perms) );
   }
#endif

   /* output time */
   if ( ! local )
   {
      if ( maxgenerators == 0 )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Number of generators:\t\t\t%d \t(max: -)\n", *nperms);
      else
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Number of generators:\t\t\t%u \t(max: %u)\n", *nperms, maxgenerators);
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


/** determine symmetry */
static
SCIP_RETCODE determineSymmetry(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
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
   assert( presoldata->nperms == 0 );
   assert( presoldata->nmaxperms == 0 );
   assert( presoldata->perms == NULL );

   presoldata->computedsym = TRUE;

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

   /* skip symmetry computation if required variables are not present */
   if ( ! (type & presoldata->symspecrequire) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Skip symmetry computation, since type does not match requirements (%d bin, %d int, %d cont); required: (%d bin, %d int, %d cont).\n",
         SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNContVars(scip) + SCIPgetNImplVars(scip),
         (presoldata->symspecrequire & (int) SYM_SPEC_BINARY) != 0,
         (presoldata->symspecrequire & (int) SYM_SPEC_INTEGER) != 0,
         (presoldata->symspecrequire & (int) SYM_SPEC_REAL) != 0);
      return SCIP_OKAY;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Required symmetry:\t\t\t(%d bin, %d int, %d cont); (fixed: %d bin, %d int, %d cont)\n",
      (presoldata->symspecrequire & (int) SYM_SPEC_BINARY) != 0,
      (presoldata->symspecrequire & (int) SYM_SPEC_INTEGER) != 0,
      (presoldata->symspecrequire & (int) SYM_SPEC_REAL) != 0,
      (presoldata->symspecrequirefixed & (int) SYM_SPEC_BINARY) != 0,
      (presoldata->symspecrequirefixed & (int) SYM_SPEC_INTEGER) != 0,
      (presoldata->symspecrequirefixed & (int) SYM_SPEC_REAL) != 0);

   if ( presoldata->symspecrequire & presoldata->symspecrequirefixed )
      SCIPwarningMessage(scip, "Warning: some required symmetries must be fixed.\n");

   /* actually compute (global) symmetry */
   /* determine maximal number of generators depending on the number of variables */
   maxgenerators = presoldata->maxgenerators;
   maxgenerators = MIN(maxgenerators, MAXGENNUMERATOR / nvars);

   SCIP_CALL( computeSymmetryGroup(scip, maxgenerators, presoldata->symspecrequirefixed, FALSE,
         &presoldata->npermvars, &presoldata->permvars, &presoldata->nperms, &presoldata->nmaxperms, &presoldata->perms, &presoldata->successful) );

   if ( ! presoldata->successful )
      return SCIP_OKAY;

   if ( presoldata->nperms == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "No symmetry found.\n");
   }
   else
   {
      /* turn off some other presolving methods in order to be sure that they do not destroy symmetry afterwards */
      assert( presoldata->nperms > 0 );

      /* domcol avoids S_2-symmetries and may not be compatible with other symmetry handling methods. */
      SCIP_CALL( SCIPsetIntParam(scip, "presolving/domcol/maxrounds", 0) );
      SCIPinfoMessage(scip, NULL, "Turned off presolver <domcol>.\n");

      /* components creates sub-SCIPs on which no symmetry handling is installed, thus turn this off. */
      SCIP_CALL( SCIPsetIntParam(scip, "constraints/components/maxprerounds", 0) );
      SCIPinfoMessage(scip, NULL, "Turned off presolver <components>.\n\n");
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of presolver
 */

/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreSymmetry)
{
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   SCIPdebugMsg(scip, "Initialization of symmetry presolver.\n");

   /* compute symmetries if not requested during presolving */
   if ( ! presoldata->detectsympresol && ! presoldata->computedsym )
   {
      /* determine symmetry here in initpre, since other plugins specify their problem type in init() */
      SCIP_CALL( determineSymmetry(scip, presoldata) );
   }

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeSymmetry)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int i;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );
   assert( presoldata != NULL );

   SCIPdebugMsg(scip, "Freeing symmetry presolver.\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->permvars, presoldata->npermvars);
   for (i = 0; i < presoldata->nperms; ++i)
   {
      SCIPfreeBlockMemoryArray(scip, &presoldata->perms[i], presoldata->npermvars);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &presoldata->perms, presoldata->nmaxperms);

   SCIPfreeBlockMemory(scip, &presoldata);

   return SCIP_OKAY;
}

#if 0
/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSymmetry)
{
   SCIP_PRESOLDATA* sourcepresoldata;
   SCIP_PRESOLDATA* targetpresoldata;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );
   assert( sourcecons != NULL );

   SCIPdebugMsg(scip, "Transforming symmetry constraint <%s> ...\n", SCIPconsGetName(sourcecons));

   sourcepresoldata = SCIPconsGetData(sourcecons);
   assert( sourcepresoldata != NULL );
   SCIP_CALL( SCIPallocBlockMemory(scip, &targetpresoldata) );

   /* copy pointers to group */
   targetpresoldata->npermvars = sourcepresoldata->npermvars;
   targetpresoldata->permvars = NULL;
   targetpresoldata->computedsym = sourcepresoldata->computedsym;
   targetpresoldata->nperms = 0;
   targetpresoldata->nmaxperms = 0;
   targetpresoldata->perms = NULL;

   /* copy variables and set up variable map */
   if ( sourcepresoldata->npermvars > 0 )
   {
      assert( sourcepresoldata->permvars != NULL );

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetpresoldata->permvars), sourcepresoldata->permvars, sourcepresoldata->npermvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, sourcepresoldata->npermvars, sourcepresoldata->permvars, targetpresoldata->permvars) );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), presol, targetpresoldata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}
#endif

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


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreSymmetry)
{
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   /* skip if we are in a restart */
   if ( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   /* skip if we are exiting */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Exitpre method of symmetry presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* compute symmetries if requested during presolving */
   if ( presoldata->detectsympresol && ! presoldata->computedsym )
   {
      SCIP_CALL( determineSymmetry(scip, presoldata) );
   }

   return SCIP_OKAY;
}


/** copy method for presolver plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopySymmetry)
{
   assert( scip != NULL );
   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   /* Do not copy, since we do not know how to efficiently copy symmetry information, but declare copy to be valid in
    * order to not disturb sub-SCIP heuristics. */

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
   assert( presoldata != 0 );

   presoldata->symspecrequire = 0;
   presoldata->symspecrequirefixed = 0;
   presoldata->npermvars = 0;
   presoldata->permvars = NULL;
   presoldata->perms = NULL;
   presoldata->nperms = 0;
   presoldata->nmaxperms = 0;
   presoldata->computedsym = FALSE;
   presoldata->successful = FALSE;

   /* include constraint handler */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC,
         PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING, presolExecSymmetry, presoldata) );
   assert( presol != NULL );

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSymmetry) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreSymmetry) );
   SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreSymmetry) );

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolvers/" PRESOL_NAME"/detectsympresol",
         "Should the symmetry be detected after presolving (otherwise before presol)?",
         &presoldata->detectsympresol, TRUE, DEFAULT_DETECTSYMPRESOL, 0, 0) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolvers/" PRESOL_NAME"/maxgenerators", "limit on the number of generators that should be produced within symmetry detection (0 = no limit)",
         &presoldata->maxgenerators, TRUE, DEFAULT_MAXGENERATORS, 0, INT_MAX, 0, 0) );

   /* possibly add description */
   if ( SYMcanComputeSymmetry() )
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetName(), SYMsymmetryGetDesc()) );
   }

   return SCIP_OKAY;
}


/** return symmetry group generators */
SCIP_RETCODE SCIPgetSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< symmetry presolver */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert( presol != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( perms != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   if ( ! presoldata->computedsym )
   {
      if ( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVED )
      {
         SCIPerrorMessage("Cannot call symmetry detection outside of presolving.\n");
         return SCIP_INVALIDCALL;
      }

      /* determine symmetry here */
      SCIP_CALL( determineSymmetry(scip, presoldata) );
   }

   *npermvars = presoldata->npermvars;
   *permvars = presoldata->permvars;
   *nperms = presoldata->nperms;
   *perms = presoldata->perms;

   return SCIP_OKAY;
}


/** specify symmetry type for which we need symmetries */
void SYMsetSpecRequirement(
   SCIP_PRESOL*          presol,             /**< symmetry presolver */
   SYM_SPEC              type                /**< variable types the callee is interested in */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   presoldata->symspecrequire |= (int) type;
}


/** specify symmetry type which symmetry group must fix */
void SYMsetSpecRequirementFixed(
   SCIP_PRESOL*          presol,             /**< symmetry presolver */
   SYM_SPEC              fixedtype           /**< variable types that callee wants to have fixed */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   presoldata->symspecrequirefixed |= (int) fixedtype;
}


/** whether symmetry should be computed for presolved system */
SCIP_Bool SYMdetectSymmetryPresolved(
   SCIP_PRESOL*          presol              /**< symmetry presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert( presol != NULL );
   assert( strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   return presoldata->detectsympresol;
}
