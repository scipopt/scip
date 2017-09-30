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

/**@file   cons_symmetry.cpp
 * @brief  constraint handler for computing and storing symmetry information about current problem
 * @author Marc Pfetsch
 * @author Thomas Rehn
 *
 * This constraint handler computes symmetries of the problem and stores this information in adequate form. It does not
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

#include <scip/cons_symmetry.h>
#include <symmetry/compute_symmetry.h>

#include <string.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "symmetry"
#define CONSHDLR_DESC          "constraint handler for computing and storing symmetry information about current problem"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -9000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE


/* default parameter values */
#define DEFAULT_MAXGENERATORS      1500      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
#define DEFAULT_DETECTSYMPRESOL    TRUE      /**< Should the symmetry be detected within presolving (otherwise before presol)? */

/* other defines */
#define MAXGENNUMERATOR        64000000      /**< determine maximal number of generators by dividing this number by the number of variables */


/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             detectsympresol;    /**< Should the symmetry be detected within presolving (otherwise before presol)? */
   int                   maxgenerators;      /**< limit on the number of generators that should be produced within symmetry detection (0 = no limit) */
   int                   symspecrequire;     /**< symmetry specification for which we need to compute symmetries */
   int                   symspecrequirefixed;/**< symmetry specification of variables which must be fixed by symmetries */
};


/** constraint data */
struct SCIP_ConsData
{
   int                   npermvars;          /**< number of variables for permutations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int                   npermbinvars;       /**< number of binary variables for permutations */
   int*                  permbinvars;        /**< indices of binary variables for permutations */
   SCIP_HASHMAP*         varmap;             /**< map of variables to indices in permvars array */
   SCIP_Bool             computedsym;        /**< Have we already tried to compute symmetries? */
   int                   nperms;             /**< number of permutations */
   int                   nmaxperms;          /**< maximal number of permutations (needed for freeing storage) */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
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


/* ------------------- map for rhs types ------------------- */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyRhstype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare the types of two rhs according to value and sense.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQRhstype)
{
   SCIP* scip;
   SYM_RHSTYPE* k1;
   SYM_RHSTYPE* k2;

   scip = (SCIP*) userptr;
   k1 = (SYM_RHSTYPE*) key1;
   k2 = (SYM_RHSTYPE*) key2;

   /* first check value */
   if ( ! SCIPisEQ(scip, k1->val, k2->val) )
      return FALSE;

   /* if still undecided, take sense */
   if ( k1->sense != k2->sense )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValRhstype)
{  /*lint --e{715}*/
   SYM_RHSTYPE* k;

   k = (SYM_RHSTYPE*) key;
   return SCIPrealHashCode(k->val) + (uint64_t) k->sense;
}


/* ------------------- map for matrix coefficient types ------------------- */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyMattype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare the types of two rhs according to value and sense.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQMattype)
{
   SCIP* scip;
   SYM_MATTYPE* k1;
   SYM_MATTYPE* k2;

   scip = (SCIP*) userptr;
   k1 = (SYM_MATTYPE*) key1;
   k2 = (SYM_MATTYPE*) key2;

   /* first check value */
   if ( ! SCIPisEQ(scip, k1->val, k2->val) )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValMattype)
{  /*lint --e{715}*/
   SYM_MATTYPE* k;

   k = (SYM_MATTYPE*) key;
   return SCIPrealHashCode(k->val);
}



/*
 * Local methods
 */





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
   int                   npermvars,          /**< number of variables in permutations */
   SCIP_VAR**            permvars,           /**< variables in permutations */
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
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &permrow, npermvars) );

   /* set up map between rows and first entry in matcoef array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rhsmatbeg, matrixdata->nrhscoef) );
   for (j = 0; j < matrixdata->nrhscoef; ++j)
      rhsmatbeg[j] = matrixdata->nmatcoef + 1;

   /* build map from rhs into matrix */
   oldrhs = matrixdata->nrhscoef + 1;
   for (j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int rhs;

      rhs = matrixdata->matrhsidx[j];
      if ( rhs != oldrhs )
      {
         assert( rhs < matrixdata->nrhscoef );
         rhsmatbeg[rhs] = j;
         oldrhs = rhs;
      }
   }

   /* create row */
   for (j = 0; j < npermvars; ++j)
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

      for (j = 0; j < npermvars; ++j)
      {
         if ( SymmetryFixVar(fixedtype, permvars[j]) && P[j] != j )
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
         assert( j < matrixdata->nmatcoef );
         assert( matrixdata->matrhsidx[j] == r1 ); /* note: row cannot be empty by construction */

         /* loop through row */
         while ( j < matrixdata->nmatcoef && matrixdata->matrhsidx[j] == r1 )
         {
            int varidx;

            assert( matrixdata->matvaridx[j] < npermvars );
            varidx = P[matrixdata->matvaridx[j]];
            assert( varidx < npermvars );
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
                  assert( j < matrixdata->nmatcoef );
                  assert( matrixdata->matrhsidx[j] == r2 );
                  assert( matrixdata->matvaridx[j] < npermvars );

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
   SCIPfreeBlockMemoryArray(scip, &permrow, npermvars);

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
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SYM_MATRIXDATA matrixdata;
   SCIP_HASHTABLE* vartypemap;
   SCIP_HASHTABLE* rhstypemap;
   SCIP_HASHTABLE* mattypemap;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SYM_MATTYPE* uniquematarray;
   SYM_VARTYPE* uniquevararray;
   SYM_RHSTYPE* uniquerhsarray;
   SYM_RHSSENSE oldsense = SYM_SENSE_UNKOWN;
   SCIP_Real oldcoef = SCIP_INVALID;
   SCIP_Real val;
   int nuniquevararray = 0;
   int nhandleconss = 0;
   int nconss;
   int nvars;
   int c;
   int j;

#ifdef SCIP_DEBUG
   /* for debugging: store original coefficients */
   SCIP_Real* matcoeforig;
   SCIP_Real* rhscoeforig;
#endif

   assert( scip != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );

   /* init */
   *npermvars = 0;
   *permvars = NULL;
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;

   /* skip if no symmetry can be computed */
   if ( ! SYMcanComputeSymmetry() )
   {
      SCIPwarningMessage(scip, "Cannot compute symmetry, since no third party software has been linked in.\n");
      return SCIP_OKAY;
   }

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   /* exit if no constraints or no variables  are available */
   if ( nconss == 0 || nvars == 0 )
      return SCIP_OKAY;

   /* before we set up the matrix, check whether we can handle all constraints */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "setppc");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "xor");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "logicor");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "knapsack");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "bounddisjunction");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   conshdlr = SCIPfindConshdlr(scip, "symmetry");
   nhandleconss += SCIPconshdlrGetNConss(conshdlr);
   if ( nhandleconss < nconss )
   {
      SCIPwarningMessage(scip, "Cannot compute symmetry, since unkown constraints are present.\n");
      return SCIP_OKAY;
   }

   conss = SCIPgetConss(scip);
   assert( conss != NULL );

   SCIPdebugMsg(scip, "Detecting %ssymmetry on %d variables and %d constraints.\n", local ? "local " : "", nvars, nconss);

   /* copy variables */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, SCIPgetVars(scip), nvars) );
   assert( vars != NULL );

   /* fill matrixdata */
   matrixdata.nmaxmatcoef = 100 * nvars;
   matrixdata.nmatcoef = 0;
   matrixdata.nrhscoef = 0;
   matrixdata.rhstypemap = NULL;
   matrixdata.mattypemap = NULL;
   matrixdata.nuniquemat = 0;
   matrixdata.nuniquevars = 0;
   matrixdata.nuniquerhs = 0;
   matrixdata.npermvars = nvars;
   matrixdata.permvars = vars;
   matrixdata.permvarcolors = NULL;

   /* prepare matrix data (use block memory, since this can become large) */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhscoef, 2 * nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhssense, 2 * nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.rhsidx, 2 * nconss) );

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
         /* get number of variables of XOR constraint (should include integer variable) */
         SCIP_Bool success;

         SCIP_CALL( SCIPgetConsNVars(scip, cons, &nconsvars, &success) );
         assert( success );
         assert( nconsvars <= nvars );
         assert( nconsvars == SCIPgetNVarsXor(scip, cons) + 1 );

         /* get variables of XOR constraint */
         SCIP_CALL( SCIPgetConsVars(scip, cons, consvars, nconsvars, &success) );
         assert( success );

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
      else if ( strcmp(conshdlrname, "symmetry") == 0 )
      {
         /* skip */
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
   assert( matrixdata.nrhscoef <= 2 * nconss );
   assert( matrixdata.nrhscoef > 0 ); /* cannot have empty rows! */

   SCIPfreeBlockMemoryArray(scip, &consvals, nvars);
   SCIPfreeBlockMemoryArray(scip, &consvars, nvars);

#ifdef SCIP_DEBUG
   /* if symmetries have to be checked store original matrix and rhs coefficient arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &matcoeforig, matrixdata.matcoef, matrixdata.nmatcoef) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &rhscoeforig, matrixdata.rhscoef, matrixdata.nrhscoef) );
#endif

   /* determine nonredundant list of coefficients: first sort */
   SCIPsortRealInt(matrixdata.matcoef, matrixdata.matidx, matrixdata.nmatcoef);
   SCIPsortRealInt(matrixdata.rhscoef, matrixdata.rhsidx, matrixdata.nrhscoef);

#ifdef SCIP_DEBUG
   assert( matcoeforig != NULL && rhscoeforig !=  NULL );
   for (j = 0; j < matrixdata.nrhscoef; ++j)
      assert( SCIPisEQ(scip, matrixdata.rhscoef[j], rhscoeforig[matrixdata.rhsidx[j]]) );
#endif

   /* create maps for coefficients to indices */
   SCIP_CALL( SCIPhashtableCreate(&vartypemap, SCIPblkmem(scip), 5 * nvars, SYMhashGetKeyVartype, SYMhashKeyEQVartype, SYMhashKeyValVartype, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&rhstypemap, SCIPblkmem(scip), 5 * matrixdata.nrhscoef, SYMhashGetKeyRhstype, SYMhashKeyEQRhstype, SYMhashKeyValRhstype, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&mattypemap, SCIPblkmem(scip), 5 * matrixdata.nmatcoef, SYMhashGetKeyMattype, SYMhashKeyEQMattype, SYMhashKeyValMattype, (void*) scip) );
   assert( vartypemap != NULL );
   assert( rhstypemap != NULL );
   assert( mattypemap != NULL );
   matrixdata.rhstypemap = rhstypemap;
   matrixdata.mattypemap = mattypemap;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &matrixdata.permvarcolors, nvars) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniquematarray, matrixdata.nmatcoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniquerhsarray, matrixdata.nrhscoef) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniquevararray, nvars) );

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
      assert( 0 <= matrixdata.matidx[j] && matrixdata.matidx[j] < matrixdata.nmatcoef );

      val = matrixdata.matcoef[j];
      if ( ! SCIPisEQ(scip, val, oldcoef) )
      {
         SYM_MATTYPE* mt;

         mt = &uniquematarray[matrixdata.nuniquemat];
         mt->val = val;
         mt->color = matrixdata.nuniquemat++;

         assert( ! SCIPhashtableExists(mattypemap, (void*) mt) );
         SCIP_CALL( SCIPhashtableInsert(mattypemap, (void*) mt) );

         oldcoef = val;
#ifdef SCIP_OUTPUT
         SCIPdebugMsg(scip, "detected new matrix entry type %f - color: %d\n", val, matrixdata.nuniquemat - 1);
#endif
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
      val = matrixdata.rhscoef[j];
      if ( ! SCIPisEQ(scip, val, oldcoef) || oldsense != sense )
      {
         SYM_RHSTYPE* rt;

         rt = &uniquerhsarray[matrixdata.nuniquerhs];
         rt->val = val;
         rt->sense = sense;

         if ( ! SCIPhashtableExists(rhstypemap, (void*) rt) )
         {
            rt->color = matrixdata.nuniquerhs++;
            SCIP_CALL( SCIPhashtableInsert(rhstypemap, (void*) rt) );

#ifdef SCIP_OUTPUT
            SCIPdebugMsg(scip, "have new rhs type %f, type: %u - color: %u\n", val, sense, matrixdata.nuniquerhs - 1);
#endif
         }
         oldcoef = val;
         oldsense = sense;
      }
   }
   assert( matrixdata.nuniquevars > 0 );
   assert( matrixdata.nuniquerhs > 0 );
   assert( matrixdata.nuniquemat > 0 );

   SCIPdebugMsg(scip, "Number of detected different variables: %d (total: %d).\n", matrixdata.nuniquevars, nvars);
   SCIPdebugMsg(scip, "Number of detected different nonzero coefficients: rhs: %d (total %d), matrix: %d (total %d).\n",
      matrixdata.nuniquerhs, matrixdata.nrhscoef, matrixdata.nuniquemat, matrixdata.nmatcoef);

   /* determine generators */
   SCIP_CALL( SYMcomputeSymmetryGenerators(scip, maxgenerators, &matrixdata, nperms, nmaxperms, perms) );

#ifdef SCIP_DEBUG
   if ( ! SCIPisStopped(scip) )
   {
      SCIP_CALL( checkSymmetriesAreSymmetries(scip, fixedtype, nvars, vars, &matrixdata, *nperms, *perms) );
   }
   SCIPfreeBlockMemoryArray(scip, &rhscoeforig, matrixdata.nrhscoef);
   SCIPfreeBlockMemoryArray(scip, &matcoeforig, matrixdata.nmatcoef);
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
   SCIPfreeBlockMemoryArray(scip, &uniquerhsarray, matrixdata.nrhscoef);
   SCIPfreeBlockMemoryArray(scip, &uniquematarray, matrixdata.nmatcoef);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.permvarcolors, nvars);
   SCIPhashtableFree(&matrixdata.mattypemap);
   SCIPhashtableFree(&matrixdata.rhstypemap);
   SCIPhashtableFree(&vartypemap);

   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhsidx, 2 * nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhssense, 2 * nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.rhscoef, 2 * nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matvaridx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matrhsidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matidx, matrixdata.nmaxmatcoef);
   SCIPfreeBlockMemoryArrayNull(scip, &matrixdata.matcoef, matrixdata.nmaxmatcoef);

   /* copy variables */
   *permvars = vars;
   *npermvars = nvars;

   return SCIP_OKAY;
}


/** determine symmetry */
static
SCIP_RETCODE determineSymmetry(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< symmetries constraint handler data */
   SCIP_CONSDATA*        consdata            /**< symmetries constraint data */
   )
{
   int maxgenerators;
   int type = 0;
   int nvars;

   assert( scip != NULL );
   assert( consdata != NULL );

   assert( ! consdata->computedsym );
   assert( consdata->npermvars == 0 );
   assert( consdata->permvars == NULL );
   assert( consdata->npermbinvars == 0 );
   assert( consdata->permbinvars == NULL );
   assert( consdata->nperms == 0 );
   assert( consdata->nmaxperms == 0 );
   assert( consdata->perms == NULL );

   consdata->computedsym = TRUE;

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
   if ( ! (type & conshdlrdata->symspecrequire) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Skip symmetry computation, since type does not match requirements (%d bin, %d int, %d cont); required: (%d bin, %d int, %d cont).\n",
         SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNContVars(scip) + SCIPgetNImplVars(scip),
         (conshdlrdata->symspecrequire & (int) SYM_SPEC_BINARY) != 0,
         (conshdlrdata->symspecrequire & (int) SYM_SPEC_INTEGER) != 0,
         (conshdlrdata->symspecrequire & (int) SYM_SPEC_REAL) != 0);
      return SCIP_OKAY;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Required symmetry:\t\t\t(%d bin, %d int, %d cont); (fixed: %d bin, %d int, %d cont)\n",
      (conshdlrdata->symspecrequire & (int) SYM_SPEC_BINARY) != 0,
      (conshdlrdata->symspecrequire & (int) SYM_SPEC_INTEGER) != 0,
      (conshdlrdata->symspecrequire & (int) SYM_SPEC_REAL) != 0,
      (conshdlrdata->symspecrequirefixed & (int) SYM_SPEC_BINARY) != 0,
      (conshdlrdata->symspecrequirefixed & (int) SYM_SPEC_INTEGER) != 0,
      (conshdlrdata->symspecrequirefixed & (int) SYM_SPEC_REAL) != 0);

   if ( conshdlrdata->symspecrequire & conshdlrdata->symspecrequirefixed )
      SCIPwarningMessage(scip, "Warning: some required symmetries must be fixed.\n");

   /* actually compute (global) symmetry */
   /* determine maximal number of generators depending on the number of variables */
   maxgenerators = conshdlrdata->maxgenerators;
   maxgenerators = MIN(maxgenerators, MAXGENNUMERATOR / nvars);

   SCIP_CALL( computeSymmetryGroup(scip, maxgenerators, conshdlrdata->symspecrequirefixed, FALSE,
         &consdata->npermvars, &consdata->permvars, &consdata->nperms, &consdata->nmaxperms, &consdata->perms) );

   if ( consdata->nperms == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "No symmetry found.\n");
   }
   else
   {
      /* turn off some other presolving methods in order to be sure that they do not destroy symmetry afterwards */
      assert( consdata->nperms > 0 );

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
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSymmetry)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMsg(scip, "Freeing symmetry constraint handler.\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreSymmetry)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != 0 );

   /* compute symmetries if not requested during presolving */
   if ( ! conshdlrdata->detectsympresol )
   {
      SCIP_CONSDATA* consdata;
      for (c = 0; c < nconss; ++c)
      {
         assert( conss[c] != 0 );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != 0 );

         if ( consdata->computedsym )
            continue;

         SCIPdebugMsg(scip, "Initialization of symmetry constraint <%s>.\n", SCIPconsGetName(conss[c]));

         /* determine symmetry here in initpre, since other plugins specify their problem type in init() */
         SCIP_CALL( determineSymmetry(scip, conshdlrdata, consdata) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSymmetry)
{  /*lint --e{715}*/
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( consdata != NULL );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Deleting symmetry constraint <%s>.\n", SCIPconsGetName(cons));

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->permvars, (*consdata)->npermvars);
   if ( (*consdata)->varmap != 0 )
   {
      SCIPhashmapFree(&(*consdata)->varmap);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->permbinvars, (*consdata)->npermbinvars);
   for (i = 0; i < (*consdata)->nperms; ++i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->perms[i], (*consdata)->npermvars);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->perms, (*consdata)->nmaxperms);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSymmetry)
{
   SCIP_CONSDATA* sourceconsdata;
   SCIP_CONSDATA* targetconsdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );

   SCIPdebugMsg(scip, "Transforming symmetry constraint <%s> ...\n", SCIPconsGetName(sourcecons));

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert( sourceconsdata != NULL );
   SCIP_CALL( SCIPallocBlockMemory(scip, &targetconsdata) );

   /* copy pointers to group */
   targetconsdata->npermvars = sourceconsdata->npermvars;
   targetconsdata->permvars = NULL;
   targetconsdata->npermbinvars = sourceconsdata->npermbinvars;
   targetconsdata->permbinvars = NULL;
   targetconsdata->varmap = 0;
   targetconsdata->computedsym = sourceconsdata->computedsym;
   targetconsdata->nperms = 0;
   targetconsdata->nmaxperms = 0;
   targetconsdata->perms = NULL;

   /* copy variables and set up variable map */
   if ( sourceconsdata->npermvars > 0 )
   {
      int cnt = 0;
      int j;

      assert( sourceconsdata->permvars != NULL );

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetconsdata->permvars), sourceconsdata->permvars, sourceconsdata->npermvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, sourceconsdata->npermvars, sourceconsdata->permvars, targetconsdata->permvars) );

      if ( sourceconsdata->npermbinvars > 0 )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetconsdata->permbinvars), sourceconsdata->permbinvars, sourceconsdata->npermbinvars) );
      }

      /* create hashmap for storing the indices of variables */
      SCIP_CALL( SCIPhashmapCreate(&targetconsdata->varmap, SCIPblkmem(scip), 5 * sourceconsdata->npermvars) );

      /* insert variables and determine binary variables */
      for (j = 0; j < targetconsdata->npermvars; ++j)
      {
         SCIP_CALL( SCIPhashmapInsert(targetconsdata->varmap, targetconsdata->permvars[j], (void*)(size_t) j) );
         if ( SCIPvarGetType(targetconsdata->permvars[j]) == SCIP_VARTYPE_BINARY )
            targetconsdata->permbinvars[cnt++] = j;
      }
   }
   else
   {
      assert( targetconsdata->npermbinvars == 0 );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetconsdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSymmetry)
{  /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSymmetry)
{  /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for primal solutions */
static
SCIP_DECL_CONSCHECK(consCheckSymmetry)
{  /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSymmetry)
{  /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   /* do nothing */
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreSymmetry)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* skip if we are in a restart */
   if ( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   /* skip if we are exiting */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Exitpre method of symmetry constraint handler ...\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* compute symmetries if requested during presolving */
   if ( conshdlrdata->detectsympresol )
   {
      SCIP_CONSDATA* consdata;
      int c;

      for (c = 0; c < nconss; ++c)
      {
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );

         if ( consdata->computedsym )
            continue;

         SCIPdebugMsg(scip, "Exitpre method of symmetry constraint <%s>.\n", SCIPconsGetName(conss[c]));

         SCIP_CALL( determineSymmetry(scip, conshdlrdata, consdata) );
      }
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSymmetry)
{  /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* do nothing */

   return SCIP_OKAY;
}


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySymmetry)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   SCIP_CALL( SCIPincludeConshdlrSymmetry(scip) );
   *valid = TRUE;

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySymmetry)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourceconshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   /* Do not copy, since we do not know how to efficiently copy symmetry information, but declare copy to be valid in
    * order to not disturb sub-SCIP heuristics. */
   *valid = TRUE;

   return SCIP_OKAY;
}





/*
 * External methods
 */

/** include symmetry constraint handler */
SCIP_RETCODE SCIPincludeConshdlrSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert( conshdlrdata != 0 );
   conshdlrdata->symspecrequire = 0;
   conshdlrdata->symspecrequirefixed = 0;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSymmetry, consEnfopsSymmetry, consCheckSymmetry, consLockSymmetry, conshdlrdata) );
   assert( conshdlr != NULL );

   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSymmetry) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreSymmetry) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSymmetry) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySymmetry, consCopySymmetry) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSymmetry) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymmetry, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSymmetry) );

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME"/detectsympresol",
         "Should the symmetry be detected within presolving (otherwise before presol)?",
         &conshdlrdata->detectsympresol, TRUE, DEFAULT_DETECTSYMPRESOL, 0, 0) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME"/maxgenerators", "limit on the number of generators that should be produced within symmetry detection (0 = no limit)",
         &conshdlrdata->maxgenerators, TRUE, DEFAULT_MAXGENERATORS, 0, INT_MAX, 0, 0) );

   /* possibly add description */
   if ( SYMcanComputeSymmetry() )
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SYMsymmetryGetName(), SYMsymmetryGetDesc()) );
   }

   return SCIP_OKAY;
}


/** create symmetry constraint */
SCIP_RETCODE SCIPcreateConsSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = NULL;

   assert( scip != NULL );
   assert( cons != NULL );

   /* find the symmetry constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("symmetry constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMsg(scip, "Creating symmetry constraint <%s>.\n", name);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->npermvars = 0;
   consdata->permvars = NULL;
   consdata->npermbinvars = 0;
   consdata->permbinvars = NULL;
   consdata->varmap = NULL;
   consdata->computedsym = FALSE;
   consdata->perms = NULL;
   consdata->nperms = 0;
   consdata->nmaxperms = 0;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE));

   return SCIP_OKAY;
}


/** return symmetry group generators */
SCIP_RETCODE SCIPgetSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< symmetry constraint handler */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP**        permvarmap,         /**< map of variables to indices in permvars array */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms               /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;

   assert( conshdlr != NULL );
   assert( npermvars != NULL );
   assert( permvars != NULL );
   assert( nperms != NULL );
   assert( perms != NULL );

   if ( SCIPconshdlrGetNConss(conshdlr) <= 0 )
   {
      SCIPerrorMessage("Could not find symmetry handling constraint.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* arbitrarly get first constraint */
   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );
   assert( conss[0] != NULL );
   consdata = SCIPconsGetData(conss[0]);

   if ( ! consdata->computedsym )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      if ( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVED )
      {
         SCIPerrorMessage("Cannot call symmetry detection outside of presolving.\n");
         return SCIP_INVALIDCALL;
      }

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != 0 );

      /* determine symmetry here */
      SCIP_CALL( determineSymmetry(scip, conshdlrdata, consdata) );
   }

   *npermvars = consdata->npermvars;
   *permvars = consdata->permvars;
   *nperms = consdata->nperms;
   *perms = consdata->perms;

   if ( permvarmap != NULL )
      *permvarmap = consdata->varmap;

   return SCIP_OKAY;
}


/** specify symmetry type for which we need symmetries */
void SYMsetSpecRequirement(
   SCIP_CONSHDLR*        conshdlr,           /**< symmetry constraint handler */
   SYM_SPEC              type                /**< variable types the callee is interested in */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != 0 );

   conshdlrdata->symspecrequire |= (int) type;
}


/** specify symmetry type which symmetry group must fix */
void SYMsetSpecRequirementFixed(
   SCIP_CONSHDLR*        conshdlr,           /**< symmetry constraint handler */
   SYM_SPEC              fixedtype           /**< variable types that callee wants to have fixed */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != 0 );

   conshdlrdata->symspecrequirefixed |= (int) fixedtype;
}


/** whether symmetry should be computed for presolved system */
SCIP_Bool SYMdetectSymmetryPresolved(
   SCIP_CONSHDLR*        conshdlr            /**< symmetry constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != 0 );

   return conshdlrdata->detectsympresol;
}
