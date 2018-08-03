/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/*  This file was written by Giacomo Nannicini,                              */
/*    Copyright (C) 2012 Singapore University of Technology and Design       */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_gmi.c
 * @brief  Gomory Mixed-Integer Cuts
 * @author Giacomo Nannicini
 * @author Marc Pfetsch
 *
 * This file implements a Gomory Mixed-Integer (GMI) cuts generator that reads cuts from the simplex tableau, applying
 * the textbook formula:
 * \f[
 *    \sum_{j \in J_I : f_j \leq f_0} f_j x_j + \sum_{j \in J_I : f_j > f_0} f_0 \frac{1-f_j}{1 - f_0} x_j +
 *    \sum_{j \in J_C : a_j \geq 0}   a_j x_j - \sum_{j \in J_C : a_j < 0} f_0  \frac{a_j}{1-f_0} x_j \geq f_0.
 * \f]
 * Here, \f$J_I\f$ and \f$J_C \subseteq \{1, \ldots, n\}\f$ are the indices of integer and continuous non basic
 * variables, respectively. The tableaux row is given by \f$a_j\f$ and its right hand side is \f$a_0\f$. The values
 * \f$f_j\f$ for \f$j = 0, \ldots, n\f$ denote the fractional values of the tableaux row and rhs, i.e., \f$f_j = a_j -
 * \lfloor a_j \rfloor\f$.
 *
 * Here is a brief description of the simplex tableau that we can expect from the SCIP LP interfaces:
 *
 * - Nonbasic columns can be at lower or upper bound, or they can be nonbasic at zero if they are free. Nonbasic columns
 *   at the upper bound must be flipped. Nonbasic free variables at zero are currently untested in the cut generator,
 *   but they should be handled properly anyway.
 *
 * - Nonbasic rows can be at lower or upper bound, depending on whether the lower or upper bound of the row is
 *   attained. SCIP always adds slack/surplus variables with a coefficient of +1: the slack variable is nonnegative in
 *   case of a <= constraint, it is nonpositive in case of a >= or ranged constraint. Therefore, slack variables
 *   corresponding to >= or ranged constraints must be flipped if the row is at its lower bound. (Ranged constraints at
 *   the upper bound do not have to - * be flipped because the variable is nonpositive.)
 *
 * Generated cuts are modified and their numerical properties are checked before being added to the LP relaxation.
 * Default parameters for cut modification and checking procedures are taken from the paper
 *
 * G. Cornuejols, F. Margot, and G. Nannicini:@n
 * On the safety of Gomory cut generators.@n
 * Mathematical Programming Computation 5, No. 4 (2013), pp. 345-395.
 *
 * In addition to the routines described in the paper above, here we additionally check the support of the cutting
 * plane.
 *
 * @todo Check whether it is worth rescaling the cut to have integral coefficients on integer variables. This may lead
 * to an integral slack variable, that has stronger cut coefficients in subsequent rounds.
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/pub_misc.h"
#include "sepa_gmi.h"

#define SEPA_NAME              "gmi"
#define SEPA_DESC              "Gomory Mixed-Integer cuts separator"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of Gomory separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        30 /**< maximal number of Gomory separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          -1 /**< maximal number of Gomory cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT      -1 /**< maximal number of Gomory cuts separated per separation round in root node */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_SEPARATEROWS       TRUE /**< separate rows with integral slack? */

#define DEFAULT_AWAY              0.005 /**< minimal fractionality of a basic variable in order to try GMI cut - default */
#define DEFAULT_MIN_VIOLATION      0.00 /**< minimal violation to accept cut - default */
#define DEFAULT_EPS_COEFF         1e-11 /**< tolerance for zeroing out small coefficients - default */
#define DEFAULT_EPS_RELAX_ABS     1e-11 /**< absolute cut rhs relaxation - default */
#define DEFAULT_EPS_RELAX_REL     1e-13 /**< relative cut rhs relaxation - default */
#define DEFAULT_MAX_DYN          1.0e+6 /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients - default */
#define DEFAULT_MAX_SUPP_ABS       1000 /**< maximum cut support - absolute value in the formula - default */
#define DEFAULT_MAX_SUPP_REL        0.1 /**< maximum cut support - relative value in the formula - default */


/** separator data */
struct SCIP_SepaData
{
   int                   maxrounds;          /**< maximal number of Gomory separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of Gomory separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of Gomory cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of Gomory cuts separated per separation round in root node */
   int                   lastncutsfound;     /**< total number of cuts found after last call of separator */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Bool             separaterows;       /**< separate rows with integral slack? */
   SCIP_Real             away;               /**< minimal fractionality of a basis variable in order to try GMI cut */
   SCIP_Real             minviolation;       /**< minimal violation to accept cut */
   SCIP_Real             epscoeff;           /**< tolerance for zeroing out small coefficients */
   SCIP_Real             epsrelaxabs;        /**< absolute cut rhs relaxation */
   SCIP_Real             epsrelaxrel;        /**< relative cut rhs relaxation */
   SCIP_Real             maxdynamism;        /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients */
   int                   maxsuppabs;         /**< maximum cut support - absolute value in the formula */
   SCIP_Real             maxsupprel;         /**< maximum cut support - relative value in the formula */
};


/*
 * local methods
 */

/** Modify the cut to make it numerically safer, and packs it from dense format to sparse format.
 *
 *  See paper "On the safety of Gomory cut generators" by Cornuejols, Margot, Nannicini for more information. Returns
 *  TRUE if cut is accepted, FALSE if it is discarded.
 */
static
SCIP_Bool modifyAndPackCut(
   SCIP*                 scip,               /**< pointer to the SCIP environment */
   SCIP_SEPADATA*        sepadata,           /**< pointer to separator data */
   int                   ncols,              /**< number of columns in the LP */
   SCIP_COL**            cols,               /**< columns of the LP */
   SCIP_Real*            densecoefs,         /**< cut in dense format on input */
   SCIP_Real*            sparsecoefs,        /**< cut coefficients in sparse format on output */
   int*                  cutind,             /**< cut indices in sparse format on output */
   int*                  cutnz,              /**< pointer to store the number of nonzero elements in the cut in sparse format on output */
   SCIP_Real*            cutrhs              /**< pointer to store the rhs of the cut, initialized to original value, modified */
   )
{
   SCIP_COL* col;
   int i;
   int c;

   assert(scip != NULL);
   assert(cols != NULL);
   assert(densecoefs != NULL);
   assert(sparsecoefs != NULL);
   assert(cutind != NULL);
   assert(cutnz != NULL);
   assert(cutrhs != NULL);

   *cutnz = 0; /* this is the current position in the cut array */

   /* Check each cut coefficient. If it is small, try set it to zero. */
   for( c = 0; c < ncols; ++c )
   {
      col = cols[c];
      assert(col != NULL);
      i = SCIPcolGetLPPos(col);
      assert( 0 <= i );

      /* Cycle over small elements that are not zero. If the element is zero, it will be discarded anyway. */
      if( EPSZ(densecoefs[i], sepadata->epscoeff) && ! SCIPisZero(scip, densecoefs[i]) )
      {
         if( densecoefs[i] > 0.0 )
         {
            /* If we would have to modify the rhs by a multiple of infinity, discard the cut altogether. */
            if( SCIPisInfinity(scip, -SCIPcolGetLb(col)) )
               return FALSE;

            /* Zero out coefficient and modify rhs to preserve validity and possibly strengthen the cut. */
            *cutrhs -= densecoefs[i] * SCIPcolGetLb(cols[c]);
         }
         else if( densecoefs[i] < 0.0 )
         {
            /* If we would have to modify the rhs by a multiple of infinity, discard the cut altogether. */
            if( SCIPisInfinity(scip, SCIPcolGetUb(col)) )
               return FALSE;

            /* Zero out coefficient and modify rhs to preserve validity and possibly strengthen the cut. */
            *cutrhs -= densecoefs[i] * SCIPcolGetUb(cols[c]);
         }
      } /* if( EPSZ(densecoefs[i], sepadata->epscoeff) && ! SCIPisZero(densecoefs[i]) ) */
      else if( ! EPSZ(densecoefs[i], sepadata->epscoeff) )
      {
         /* cut coefficient is large enough - keep it and write in sparse form */
	 sparsecoefs[*cutnz] = densecoefs[i];
         cutind[*cutnz] = c;
         (*cutnz)++;
      }
   } /* for( c = 0; c < ncols; ++c ) */

   /* Relax rhs of the cut */
   *cutrhs += REALABS(*cutrhs) * sepadata->epsrelaxrel + sepadata->epsrelaxabs;

   return (*cutnz > 0) ? TRUE : FALSE;
}

/** Check the numerical properties of the cut.
 *
 *  See paper "On the safety of Gomory cut generators" by Cornuejols, Margot, Nannicini for more information. Returns
 *  TRUE if cut is accepted, FALSE if it is discarded.
 */
static
SCIP_Bool checkNumerics(
   SCIP*                 scip,               /**< pointer to the SCIP environment */
   SCIP_SEPADATA*        sepadata,           /**< pointer to separator data */
   int                   ncols,              /**< number of columns in the LP */
   SCIP_COL**            cols,               /**< columns of the LP */
   SCIP_Real*            cutcoefs,           /**< cut in sparse format */
   int*                  cutind,             /**< cut indices in sparse format */
   int                   cutnz,              /**< number of nonzero elements in the cut in sparse format */
   SCIP_Real             cutrhs,             /**< rhs of the cut */
   SCIP_Real*            cutact              /**< pointer to store activity of the cut at the current LP optimum will go here on output */
   )
{
   SCIP_Real violation;
   SCIP_Real mincoef;
   SCIP_Real maxcoef;
   int i;

   assert(scip != NULL);
   assert(cols != NULL);
   assert(cutcoefs != NULL);
   assert(cutind != NULL);
   assert(cutact != NULL);
   assert(cutnz > 0);

   /* Check maximum support */
   if( cutnz > ncols * sepadata->maxsupprel + sepadata->maxsuppabs )
   {
      SCIPdebugMsg(scip, "Cut too dense (%d > %d).\n", cutnz, (int) (ncols * sepadata->maxsupprel + sepadata->maxsuppabs));
      return FALSE;
   }

   /* Compute cut violation and dynamism */
   mincoef = SCIP_REAL_MAX;
   maxcoef = 0.0;
   *cutact = 0.0;

   for( i = 0; i < cutnz; ++i )
   {
      mincoef = MIN(mincoef, REALABS(cutcoefs[i])); /*lint !e666*/
      maxcoef = MAX(maxcoef, REALABS(cutcoefs[i])); /*lint !e666*/
      *cutact += cutcoefs[i] * SCIPcolGetPrimsol(cols[cutind[i]]);
   }

   /* Check dynamism */
   if( maxcoef > mincoef * sepadata->maxdynamism )
   {
      SCIPdebugMsg(scip, "Cut too dynamic (%g > %g).\n", maxcoef, mincoef * sepadata->maxdynamism);
      return FALSE;
   }

   /* Check minimum violation */
   violation = *cutact - cutrhs;
   if( REALABS(cutrhs) > 1.0 )
      violation /= REALABS(cutrhs);

   return (violation >= sepadata->minviolation) ? TRUE : FALSE;
}

/** Method to obtain a GMI in the space of the original variables from a row of the simplex tableau.
 *
 *  Returns TRUE if cut is successfully created, FALSE if no cut was generated or if it should be discarded. If the
 *  function returns FALSE, the contents of cutcoefs, cutind, cutnz, cutrhs, cutact may be garbage.
 */
static
SCIP_Bool getGMIFromRow(
   SCIP*                 scip,               /**< pointer to the SCIP environment */
   SCIP_SEPADATA*        sepadata,           /**< pointer to separator data */
   int                   ncols,              /**< number of columns in the LP */
   int                   nrows,              /**< number of rows in the LP */
   SCIP_COL**            cols,               /**< columns of the LP */
   SCIP_ROW**            rows,               /**< rows of the LP */
   SCIP_Real*            binvrow,            /**< row of the basis inverse */
   SCIP_Real*            binvarow,           /**< row of the simplex tableau */
   SCIP_Real             rowrhs,             /**< rhs of the tableau row, i.e., corresponding element in the LP solution */
   SCIP_Real*            cutcoefs,           /**< array for cut elements in sparse format - must be of size ncols */
   int*                  cutind,             /**< array for indices of nonzero cut coefficients - must be of size ncols */
   int*                  cutnz,              /**< pointer to store number of nonzero elements in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store cut rhs */
   SCIP_Real*            cutact,             /**< pointer to store cut activity at the current LP optimum - only meaningful if returns TRUE */
   SCIP_Real*            workcoefs           /**< working array of size ncols, allocated by caller for efficiency */
   )
{
   SCIP_COL* col;
   SCIP_ROW* row;
   SCIP_Real rowelem;
   SCIP_Real cutelem;
   SCIP_Real f0;
   SCIP_Real ratiof0compl;
   SCIP_Bool success;
   int i;
   int c;

   assert(scip != NULL);
   assert(cols != NULL);
   assert(rows != NULL);
   assert(binvrow != NULL);
   assert(binvarow != NULL);
   assert(cutcoefs != NULL);
   assert(cutind != NULL);
   assert(cutnz != NULL);
   assert(cutrhs != NULL);
   assert(cutact != NULL);
   assert(workcoefs != NULL);

   /* Compute cut fractionality f0 and f0/(1-f0). */
   f0 = SCIPfeasFrac(scip, rowrhs);
   ratiof0compl = f0/(1-f0);

   /* rhs of the cut is the fractional part of the LP solution for the basic variable */
   *cutrhs = -f0;

   /* clear cutcoefs */
   BMSclearMemoryArray(workcoefs, ncols);

   /* Generate cut coefficients for the original variables. We first use workcoefs to store the cut in dense form, then
    * we clean and pack the cut to sparse form in cutcoefs. */
   for( c = 0; c < ncols; ++ c)
   {
      col = cols[c];
      assert( col != NULL );

      /* Get simplex tableau element. */
      switch ( SCIPcolGetBasisStatus(col) )
      {
      case SCIP_BASESTAT_LOWER:
         /* Take element if nonbasic at lower bound. */
         rowelem = binvarow[c];
         break;
      case SCIP_BASESTAT_UPPER:
         /* Flip element if nonbasic at upper bound. */
         rowelem = -binvarow[c];
         break;
      case SCIP_BASESTAT_ZERO:
         /* Nonbasic free variable at zero: cut coefficient is zero, skip */
         continue;
      case SCIP_BASESTAT_BASIC:
      default:
         /* Basic variable: skip */
         continue;
      }

      /* Integer variables */
      if( SCIPcolIsIntegral(col) )
      {
         /* If cutelem < 0, then we know SCIPisZero(scip, cutelem) is true and hope it doesn't do much damage. */
         cutelem = SCIPfrac(scip, rowelem);

         if( cutelem > f0 )
         {
            /* cut element if f > f0 */
            cutelem = -((1.0 - cutelem) * ratiof0compl);
         }
         else
         {
            /* cut element if f <= f0 */
            cutelem = -cutelem;
         }
      }
      /* Continuous variables */
      else
      {
         if( rowelem < 0.0 )
         {
            /* cut element if f < 0 */
            cutelem = rowelem * ratiof0compl;
         }
         else
         {
            /* cut element if f >= 0 */
            cutelem = -rowelem;
         }
      }

      if( ! SCIPisZero(scip, cutelem) )
      {
         /* Unflip if necessary, and adjust rhs if at lower or upper bound. */
         if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
         {
            cutelem = -cutelem;
            *cutrhs += cutelem * SCIPcolGetUb(col);
         }
         else
            *cutrhs += cutelem * SCIPcolGetLb(col);

         /* Add coefficient to cut in dense form. */
         workcoefs[SCIPcolGetLPPos(col)] = cutelem;
      }
   } /* for( c = 0; c < ncols; ++c) */

   /* Generate cut coefficients for the slack variables. */
   for( c = 0; c < nrows; ++c )
   {
      row = rows[c];
      assert( row != NULL );

      /* Get simplex tableau element. */
      switch ( SCIProwGetBasisStatus(row) )
      {
      case SCIP_BASESTAT_LOWER:
         /* Take element if nonbasic at lower bound. */
         rowelem = binvrow[SCIProwGetLPPos(row)];
         /* But if this is a >= or ranged constraint at the lower bound, we have to flip the row element. */
         if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
            rowelem = -rowelem;
         break;
      case SCIP_BASESTAT_UPPER:
         /* Take element if nonbasic at upper bound - see notes at beginning of file: only nonpositive slack variables
          * can be nonbasic at upper, therefore they should be flipped twice and we can take the element directly. */
         rowelem = binvrow[SCIProwGetLPPos(row)];
         break;
      case SCIP_BASESTAT_ZERO:
         /* Nonbasic free variable at zero: cut coefficient is zero, skip */
         SCIPdebugMsg(scip, "Free nonbasic slack variable, this should not happen!\n");
         continue;
      case SCIP_BASESTAT_BASIC:
      default:
         /* Basic variable: skip */
         continue;
      }

      /* Check if row is integral and will stay integral through the Branch-and-Cut tree; if so, strengthen
       * coefficient */
      if( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) )
      {
         /* If cutelem < 0, then we know SCIPisZero(scip, cutelem) is true and hope it doesn't do much damage. */
         cutelem = SCIPfrac(scip, rowelem);

         if( cutelem > f0 )
         {
            /* cut element if f > f0 */
            cutelem = -((1.0 - cutelem) * ratiof0compl);
         }
         else
         {
            /* cut element if f <= f0 */
            cutelem = -cutelem;
         }
      }
      else
      {
	 if( rowelem < 0.0 )
	 {
	    /* cut element if f < 0 */
	    cutelem = rowelem * ratiof0compl;
	 }
	 else
	 {
	    /* cut element if f >= 0 */
	    cutelem = -rowelem;
	 }
      }

      if( ! SCIPisZero(scip, cutelem) )
      {
         /* Coefficient is large enough, we can keep it. */
         SCIP_COL** rowcols;
         SCIP_Real* rowvals;

         SCIP_Real act;
         SCIP_Real rlhs;
         SCIP_Real rrhs;
         SCIP_Real rhsslack;

         /* get lhs/rhs */
         rlhs = SCIProwGetLhs(row);
         rrhs = SCIProwGetRhs(row);
         assert( SCIPisLE(scip, rlhs, rrhs) );
         assert( ! SCIPisInfinity(scip, rlhs) || ! SCIPisInfinity(scip, rrhs) );

         /* If the slack variable is fixed, we can ignore this cut coefficient. */
         if( SCIPisFeasZero(scip, rrhs - rlhs) )
            continue;

         act = SCIPgetRowLPActivity(scip, row);
         rhsslack = rrhs - act;

         /* Unflip slack variable and adjust rhs if necessary. */
         if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
         {
            /* If >= or ranged constraint, flip element back to original */
            assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );
            cutelem = -cutelem;
         }

         rowcols = SCIProwGetCols(row);
         rowvals = SCIProwGetVals(row);

         /* Eliminate slack variable. */
         for( i = 0; i < SCIProwGetNLPNonz(row); ++i )
            workcoefs[SCIPcolGetLPPos(rowcols[i])] -= cutelem * rowvals[i];

         if ( SCIPisFeasZero(scip, rhsslack) )
            *cutrhs -= cutelem * (rrhs - SCIProwGetConstant(row));
         else
         {
            assert( SCIPisFeasZero(scip, act - rlhs) );
            *cutrhs -= cutelem * (rlhs - SCIProwGetConstant(row));
         }
      }
   } /* for( c = 0; c < nrows; ++ c) */

   /* Initialize cut activity. */
   *cutact = 0.0;

   /* Modify cut to make it numerically safer, and check that it is numerically safe. */
   success = modifyAndPackCut(scip, sepadata, ncols, cols, workcoefs, cutcoefs, cutind, cutnz, cutrhs);
   if ( success )
   {
      success = checkNumerics(scip, sepadata, ncols, cols, cutcoefs, cutind, *cutnz, *cutrhs, cutact);
      SCIPdebugMsg(scip, "checkNumerics returned: %u.\n", success);
      return success;
   }
   SCIPdebugMsg(scip, "modifyAndPackCut was not successful.\n");

   return FALSE;
}


/*
 * Callback methods
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyGMI)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaGMI(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeGMI)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpGMI)
{  /*lint --e{715}*/
   char cutname[SCIP_MAXSTRLEN];
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* binvrow;
   SCIP_Real* binvarow;
   SCIP_Real* cutcoefs;
   SCIP_Real* workcoefs;
   SCIP_Real cutrhs;
   int* cutind;
   int* basisind;
   int nvars;
   int ncols;
   int nrows;
   int ncalls;
   int depth;
   int maxsepacuts;
   int ncuts;
   int cutnz;
   int c;
   int i;
   int j;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* Only call separator, if we are not close to terminating. */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* Only call separator, if an optimal LP solution is at hand. */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* Only call separator, if the LP solution is basic. */
   if( ! SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* Only call separator, if there are fractional variables. */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* Only call the Gomory cut separator a given number of times at each node. */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot) || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* exit if LP is trivial */
   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &workcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutind, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvarow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );

   /* get basis indices */
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
      maxsepacuts = sepadata->maxsepacutsroot;
   else
      maxsepacuts = sepadata->maxsepacuts;

   if( maxsepacuts == -1 )
      maxsepacuts = INT_MAX;

   /* For all basic columns belonging to integer variables, try to generate a Gomory cut. */
   ncuts = 0;
   for( i = 0; i < nrows && ncuts < maxsepacuts && ! SCIPisStopped(scip) && *result != SCIP_CUTOFF; ++i )
   {
      SCIP_Bool tryrow;
      SCIP_Real primsol;

      tryrow = FALSE;
      c = basisind[i];
      primsol = SCIP_INVALID;

      SCIPdebugMsg(scip, "Row %d basic variable %d with value %f\n", i, basisind[i], (c >= 0) ? SCIPcolGetPrimsol(cols[c]) : SCIPgetRowActivity(scip, rows[-c-1]));
      if( c >= 0 )
      {
         SCIP_VAR* var;
         assert(c < ncols);
         assert(cols[c] != NULL);
         var = SCIPcolGetVar(cols[c]);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            primsol = SCIPcolGetPrimsol(cols[c]);
            assert(SCIPgetVarSol(scip, var) == primsol); /*lint !e777*/

            if( (SCIPfeasFrac(scip, primsol) >= sepadata->away) && (SCIPfeasFrac(scip, primsol) <= 1.0 - sepadata->away) )
            {
               SCIPdebugMsg(scip, "trying Gomory cut for col <%s> [%g] row %i\n", SCIPvarGetName(var), primsol, i);
               tryrow = TRUE;
            }
         }
      }
      else if( sepadata->separaterows )
      {
         SCIP_ROW* row;
         assert(0 <= -c-1 && -c-1 < nrows);
         row = rows[-c-1];
         if( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) )
         {
            /* Compute value of the slack variable (we only care about the correct fractionality) */
            if ( SCIPisInfinity(scip, SCIProwGetRhs(row)) )
               primsol = SCIProwGetLhs(row) - SCIPgetRowLPActivity(scip, row);
            else
               primsol = SCIProwGetRhs(row) - SCIPgetRowLPActivity(scip, row);

            if( (SCIPfeasFrac(scip, primsol) >= sepadata->away) && (SCIPfeasFrac(scip, primsol) <= 1.0 - sepadata->away) )
            {
               SCIPdebugMsg(scip, "trying Gomory cut for row <%s> [%g]\n", SCIProwGetName(row), primsol);
               SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
               tryrow = TRUE;
            }
         }
      }

      if( tryrow )
      {
         SCIP_Real cutact;
         SCIP_Bool success;
         SCIP_Bool cutislocal;

         /* get the row of B^-1 for this basic integer variable with fractional solution value */
         SCIP_CALL( SCIPgetLPBInvRow(scip, i, binvrow, NULL, NULL) );

         /* get the tableau row for this basic integer variable with fractional solution value */
         SCIP_CALL( SCIPgetLPBInvARow(scip, i, binvrow, binvarow, NULL, NULL) );

         /* this is an approximation (one could also pass over coefficients and check whether local rows have been used): */
         cutislocal = (depth != 0) ? TRUE : FALSE;

         /* create a GMI cut out of the simplex tableau row */
         success = getGMIFromRow(scip, sepadata, ncols, nrows, cols, rows, binvrow, binvarow, primsol, cutcoefs, cutind, &cutnz, &cutrhs, &cutact, workcoefs);

         SCIPdebugMsg(scip, " -> success = %u: %g <= %g\n", success, cutact, cutrhs);

         /* if successful, add the row as a cut */
         if( success )
         {
            SCIP_ROW* cut;

            /* construct cut name */
            if( c >= 0 )
               (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "gmi%d_x%d", SCIPgetNLPs(scip), c);
            else
               (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "gmi%d_s%d", SCIPgetNLPs(scip), -c-1);

            /* create empty cut */
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );

            /* cache the row extension and only flush them if the cut gets added */
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

            /* collect all non-zero coefficients */
            for( j = 0; j < cutnz; ++j )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, SCIPcolGetVar(cols[cutind[j]]), cutcoefs[j]) );
            }

            if( SCIProwGetNNonz(cut) == 0 )
            {
               assert(SCIPisFeasNegative(scip, cutrhs));
               SCIPdebugMsg(scip, " -> Gomory cut detected infeasibility with cut 0 <= %f\n", cutrhs);
               *result = SCIP_CUTOFF;
               break;
            }

            /* Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound
               changes); the latter cuts will be handeled internally in sepastore. */
            if( SCIProwGetNNonz(cut) == 1 || SCIPisCutEfficacious(scip, NULL, cut) )
            {
               SCIP_Bool infeasible;

               SCIPdebugMsg(scip, " -> Gomory cut for <%s>: act=%f, rhs=%f, eff=%f\n",
                  c >= 0 ? SCIPvarGetName(SCIPcolGetVar(cols[c])) : SCIProwGetName(rows[-c-1]),
                  cutact, cutrhs, SCIPgetCutEfficacy(scip, NULL, cut));

               SCIPdebugMsg(scip, " -> found Gomory cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                  cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                  SCIPgetCutEfficacy(scip, NULL, cut),
                  SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                  SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));

               /* flush all changes before adding the cut */
               SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

               SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );

               /* add global cuts that are not implicit bound changes to the cut pool */
               if( ! cutislocal && SCIProwGetNNonz(cut) > 1 )
               {
                  SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               }

               if ( infeasible )
                  *result = SCIP_CUTOFF;
               else
                  *result = SCIP_SEPARATED;
               ncuts++;
            }

            /* release the row */
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &binvarow);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &workcoefs);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &cutind);

   SCIPdebugMsg(scip, "end searching Gomory cuts: found %d cuts.\n", ncuts);

   sepadata->lastncutsfound = SCIPgetNCutsFound(scip);

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the GMI MIR cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaGMI(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->lastncutsfound = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpGMI, NULL, sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyGMI) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeGMI) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gmi/maxrounds",
         "maximal number of gmi separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gmi/maxroundsroot",
         "maximal number of gmi separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gmi/maxsepacuts",
         "maximal number of gmi cuts separated per separation round (-1: unlimited)",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gmi/maxsepacutsroot",
         "maximal number of gmi cuts separated per separation round in the root node (-1: unlimited)",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gmi/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gmi/separaterows",
         "separate rows with integral slack",
         &sepadata->separaterows, FALSE, DEFAULT_SEPARATEROWS, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/away",
         "minimal fractionality of a basic variable in order to try GMI cut",
         &sepadata->away, FALSE, DEFAULT_AWAY, 0.0, 0.5, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/minviolation",
         "minimal violation to accept cut",
         &sepadata->minviolation, FALSE, DEFAULT_MIN_VIOLATION, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/epscoeff",
         "tolerance for zeroing out small coefficients",
         &sepadata->epscoeff, FALSE, DEFAULT_EPS_COEFF, 0.0, 0.01, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/epsrelaxabs",
         "absolute cut rhs relaxation",
         &sepadata->epsrelaxabs, FALSE, DEFAULT_EPS_RELAX_ABS, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/epsrelaxrel",
         "relative cut rhs relaxation",
         &sepadata->epsrelaxrel, FALSE, DEFAULT_EPS_RELAX_REL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/maxdynamism",
         "maximal valid range max(|weights|)/min(|weights|) of cut coefficients",
         &sepadata->maxdynamism, FALSE, DEFAULT_MAX_DYN, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gmi/maxsuppabs",
         "maximum cut support - absolute value in the formula",
         &sepadata->maxsuppabs, FALSE, DEFAULT_MAX_SUPP_ABS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gmi/maxsupprel",
         "maximum cut support - relative value in the formula",
         &sepadata->maxsupprel, FALSE, DEFAULT_MAX_SUPP_REL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
