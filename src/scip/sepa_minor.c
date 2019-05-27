#define SCIP_STATISTIC
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

/**@file   sepa_minor.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  minor separator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_minor.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_iterator.h"
#include "nlpi/nlpi_ipopt.h"

#define SEPA_NAME              "minor"
#define SEPA_DESC              "separator template"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXMINORS           100 /**< default maximum number for minors (0: no limit) */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_VAR**            minors;             /**< variables of 2x2 minors; each minor is stored like (auxvar_x^2,auxvar_y^2,auxvar_xy) */
   int                   nminors;            /**< total number of minors */
   int                   minorssize;         /**< size of minors array */
   int                   maxminors;          /**< maximum number for minors (0: no limit) */
   SCIP_Bool             detectedminors;     /**< has the minor detection beeing called? */
};

/*
 * Local methods
 */

/** helper method to store a 2x2 minor in the separation data */
static
SCIP_RETCODE sepadataAddMinor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR*             auxvarxx,           /**< auxiliary variable for x*x */
   SCIP_VAR*             auxvaryy,           /**< auxiliary variable for y*y */
   SCIP_VAR*             auxvarxy            /**< auxiliary variable for x*y */
   )
{
   assert(sepadata != NULL);
   assert(auxvarxx != NULL);
   assert(auxvaryy != NULL);
   assert(auxvarxy != NULL);
   assert(auxvarxx != auxvaryy);
   assert(auxvarxx != auxvarxy);
   assert(auxvaryy != auxvarxy);

   SCIPdebugMsg(scip, "store 2x2 minor: %s %s %s\n", SCIPvarGetName(auxvarxx), SCIPvarGetName(auxvaryy), SCIPvarGetName(auxvarxy));

   /* reallocate if necessary */
   if( sepadata->minorssize < 3 * (sepadata->nminors + 1) )
   {
      int newsize = SCIPcalcMemGrowSize(scip, 3 * (sepadata->nminors + 1));
      assert(newsize > 3 * (sepadata->nminors + 1));

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->minors), sepadata->minorssize, newsize) );
      sepadata->minorssize = newsize;
   }

   /* store minor */
   sepadata->minors[3 * sepadata->nminors] = auxvarxx;
   sepadata->minors[3 * sepadata->nminors + 1] = auxvaryy;
   sepadata->minors[3 * sepadata->nminors + 2] = auxvarxy;
   ++(sepadata->nminors);

   /* capture variables */
   SCIPcaptureVar(scip, auxvarxx);
   SCIPcaptureVar(scip, auxvaryy);
   SCIPcaptureVar(scip, auxvarxy);

   return SCIP_OKAY;
}

/** helper method to clear separation data */
static
SCIP_RETCODE sepadataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   int i;

   assert(sepadata != NULL);

   SCIPdebugMsg(scip, "clear separation data\n");

   /* release captured variables */
   for( i = 0; i < 3 * sepadata->nminors; ++i )
   {
      assert(sepadata->minors[i] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &sepadata->minors[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &sepadata->minors, sepadata->minorssize);

   /* reset counters */
   sepadata->nminors = 0;
   sepadata->minorssize = 0;

   return SCIP_OKAY;
}

/** method to detect and store minors; the cut generated in separatePoint() looks like
 *
 *    v_1^2 X_ii + 2 v_1 v_2 X_ij + v_2^2 X_jj >= 0
 *
 *  The cut is only useful if X_ii and X_jj are locked downwards. Depending on the sign of
 *  v_1 v_2, X_ij needs to be up- or down-locked.
 *
 *  @todo It is checked whether bilinear terms are locked in both directions. We can detect the
 *        expressions here and check their locks depending on sign(v_1 v_2) in separatePoints().
 */
static
SCIP_RETCODE detectMinors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_HASHMAP* exprmap;
   SCIP_HASHMAP* quadmap;
   SCIP_VAR** xs;
   SCIP_VAR** ys;
   SCIP_VAR** auxvars;
   int* perm = NULL;
   int nbilinterms = 0;
   int nquadterms = 0;
   int c;
   int i;

#ifdef SCIP_STATISTIC
   SCIP_Real totaltime = -SCIPgetTotalTime(scip);
#endif

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* check whether minor detection has been called already */
   if( sepadata->detectedminors )
      return SCIP_OKAY;

   assert(sepadata->minors == NULL);
   assert(sepadata->nminors == 0);

   /* we assume that the auxiliary variables in the expression constraint handler have been already generated */
   sepadata->detectedminors = TRUE;

   /* check whether there are expression constraints available */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   if( conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "call detectMinors()\n");

   /* allocate memory */
   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&exprmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&quadmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxvars, SCIPgetNVars(scip)) );

   for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSEXPR_EXPR* expr;

      cons = SCIPconshdlrGetConss(conshdlr)[c];
      assert(cons != NULL);
      assert(SCIPgetExprConsExpr(scip, cons));

      SCIP_CALL( SCIPexpriteratorInit(it, SCIPgetExprConsExpr(scip, cons), SCIP_CONSEXPRITERATOR_DFS, FALSE) );
      SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_ENTEREXPR);

      for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
      {
         SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
         SCIP_CONSEXPR_EXPR** children;
         SCIP_VAR* auxvar;

         SCIPdebugMsg(scip, "visit expression %p in constraint %s\n", (void*)expr, SCIPconsGetName(cons));

         /* check whether the expression as an auxiliary variable */
         auxvar = SCIPgetConsExprExprAuxVar(expr);
         if( auxvar == NULL )
            continue;

         /* check whether the expression has been considered in another constraint */
         if( SCIPhashmapExists(exprmap, (void*)expr) )
            continue;

         exprhdlr = SCIPgetConsExprExprHdlr(expr);
         assert(exprhdlr != NULL);
         children = SCIPgetConsExprExprChildren(expr);

         /* check for expr = (x)^2 */
         if( SCIPgetConsExprExprNChildren(expr) == 1 && exprhdlr == SCIPgetConsExprExprHdlrPower(conshdlr)
            && SCIPgetConsExprExprPowExponent(expr) == 2.0
            && SCIPisConsExprExprVar(children[0])
            && SCIPgetConsExprExprNLocksNeg(expr) > 0 )
         {
            SCIP_VAR* quadvar;

            assert(children[0] != NULL);

            quadvar = SCIPgetConsExprExprVarVar(children[0]);
            assert(quadvar != NULL);
            assert(!SCIPhashmapExists(quadmap, (void*)quadvar));
            SCIPdebugMsg(scip, "found %s = (%s)^2\n", SCIPvarGetName(auxvar), SCIPvarGetName(quadvar));

            /* hash the quadratic variable to its corresponding auxiliary variable */
            SCIP_CALL( SCIPhashmapInsert(quadmap, (void*)quadvar, auxvar) );
            ++nquadterms;

            /* add expression to the map to not reconsider it */
            SCIP_CALL( SCIPhashmapInsert(exprmap, (void*)expr, NULL) );
         }
         /* check for expr = x * y */
         else if( SCIPgetConsExprExprNChildren(expr) == 2 && exprhdlr == SCIPgetConsExprExprHdlrProduct(conshdlr)
            && SCIPisConsExprExprVar(children[0]) && SCIPisConsExprExprVar(children[1])
            && SCIPgetConsExprExprNLocksNeg(expr) > 0 && SCIPgetConsExprExprNLocksPos(expr) > 0 )
         {
            assert(children[0] != NULL);
            assert(children[1] != NULL);

            xs[nbilinterms] = SCIPgetConsExprExprVarVar(children[0]);
            ys[nbilinterms] = SCIPgetConsExprExprVarVar(children[1]);
            auxvars[nbilinterms] = auxvar;
            SCIPdebugMsg(scip, "found %s = %s * %s\n", SCIPvarGetName(auxvar), SCIPvarGetName(xs[nbilinterms]), SCIPvarGetName(ys[nbilinterms]));
            ++nbilinterms;

            /* add expression to the map to not reconsider it */
            SCIP_CALL( SCIPhashmapInsert(exprmap, (void*)expr, NULL) );
         }
      }
   }
   assert(nbilinterms < SCIPgetNVars(scip));
   SCIPdebugMsg(scip, "stored %d bilinear terms in total\n", nbilinterms);

   /* permute bilinear terms if there are be too many of them; the motivation for this is that we don't want to
    * prioritize variables because of the order in the bilinear terms where they appear; however, variables that
    * appear more often in bilinear terms might be more important than others so the corresponding bilinear terms
    * are more likely to be chosen
    */
   if( sepadata->maxminors > 0 && sepadata->maxminors < nbilinterms && sepadata->maxminors < SQR(nquadterms) )
   {
      SCIP_RANDNUMGEN* randnumgen;

      /* TODO use global seed */
      SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 0, 0) );
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nbilinterms) );

      for( i = 0; i < nbilinterms; ++i )
         perm[i] = i;

      /* permute array */
      SCIPrandomPermuteIntArray(randnumgen, perm, 0, nbilinterms);

      SCIPfreeRandom(scip, &randnumgen);
   }

   /* store 2x2 principal minors */
   for( i = 0; i < nbilinterms && (sepadata->maxminors == 0 || sepadata->nminors < sepadata->maxminors); ++i )
   {
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_VAR* auxvar;

      if( perm == NULL )
      {
         x = xs[i];
         y = ys[i];
         auxvar = auxvars[i];
      }
      else
      {
         x = xs[perm[i]];
         y = ys[perm[i]];
         auxvar = auxvars[perm[i]];
      }

      assert(x != NULL);
      assert(y != NULL);
      assert(auxvar != NULL);
      assert(x != y);

      if( SCIPhashmapExists(quadmap, (void*)x) && SCIPhashmapExists(quadmap, (void*)y) )
      {
         SCIP_VAR* auxvarxx;
         SCIP_VAR* auxvaryy;

         auxvarxx = (SCIP_VAR*)SCIPhashmapGetImage(quadmap, (void*)x);
         assert(auxvarxx != NULL);
         auxvaryy = (SCIP_VAR*)SCIPhashmapGetImage(quadmap, (void*)y);
         assert(auxvaryy != NULL);

         /* store minor into te separation data */
         SCIP_CALL( sepadataAddMinor(scip, sepadata, auxvarxx, auxvaryy, auxvar) );
      }
   }
   SCIPdebugMsg(scip, "found %d 2x2 minors in total\n", sepadata->nminors);

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &perm);
   SCIPfreeBufferArray(scip, &auxvars);
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);
   SCIPhashmapFree(&quadmap);
   SCIPhashmapFree(&exprmap);
   SCIPexpriteratorFree(&it);

#ifdef SCIP_STATISTIC
   totaltime += SCIPgetTotalTime(scip);
   SCIPstatisticMessage("MINOR DETECT %s %f %d\n", SCIPgetProbName(scip), totaltime, sepadata->nminors);
#endif

   return SCIP_OKAY;
}

/** separates cuts for stored minors */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_Real minefficacy;
   int i;

   assert(sepa != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* TODO add a parameter for this */
   minefficacy = 1e-4;

   /* check whether there are some minors available */
   if( sepadata->nminors == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   for( i = 0; i < sepadata->nminors && (*result != SCIP_CUTOFF); ++i )
   {
      SCIP_VAR* xx;
      SCIP_VAR* yy;
      SCIP_VAR* xy;
      SCIP_Real solxx;
      SCIP_Real solyy;
      SCIP_Real solxy;
      SCIP_Real determinant;

      xx = sepadata->minors[3*i];
      assert(xx != NULL);
      solxx = SCIPgetSolVal(scip, sol, xx);

      yy = sepadata->minors[3*i + 1];
      assert(yy != NULL);
      solyy = SCIPgetSolVal(scip, sol, yy);

      xy = sepadata->minors[3*i + 2];
      assert(xy != NULL);
      solxy = SCIPgetSolVal(scip, sol, xy);

      determinant = solxx * solyy - SQR(solxy);
      SCIPdebugMsg(scip, "(%s,%s,%s) = (%g,%g,%g) implies determinant = %g\n", SCIPvarGetName(xx), SCIPvarGetName(yy),
         SCIPvarGetName(xy), solxx, solyy, solxy, determinant);

      /* check whether minor is positive semi-definit in the current solution */
      if( SCIPisFeasLT(scip, determinant, 0.0) )
      {
         SCIP_Real eigenvals[2];
         SCIP_Real matrix[4];
         int k;

         matrix[0] = solxx;
         matrix[1] = solxy;
         matrix[2] = solxy;
         matrix[3] = solyy;

         /* use LAPACK to compute the eigenvalues and eigenvectors */
         if( LapackDsyev(TRUE, 2, matrix, eigenvals) != SCIP_OKAY )
         {
            SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors of augmented quadratic form matrix.\n");
            continue;
         }

         /* at least one eigenvalue needs to be negative */
         assert(SCIPisLT(scip, eigenvals[0], 0.0) || SCIPisLT(scip, eigenvals[1], 0.0));

         SCIPdebugMsg(scip, "eigenvec[0] = (%g,%g) with eigenvalue[0] = %g\n", matrix[0], matrix[1], eigenvals[0]);
         assert(SCIPisRelEQ(scip, solxx * matrix[0] + solxy * matrix[1], eigenvals[0] * matrix[0]));
         assert(SCIPisRelEQ(scip, solxy * matrix[0] + solyy * matrix[1], eigenvals[0] * matrix[1]));

         SCIPdebugMsg(scip, "eigenvec[1] = (%g,%g) with eigenvalue[1] = %g\n", matrix[2], matrix[3], eigenvals[1]);
         assert(SCIPisRelEQ(scip, solxx * matrix[2] + solxy * matrix[3], eigenvals[1] * matrix[2]));
         assert(SCIPisRelEQ(scip, solxy * matrix[2] + solyy * matrix[3], eigenvals[1] * matrix[3]));

         /* try to generate a cut for every negative eigenvalue */
         for( k = 0; k < 2 && (*result != SCIP_CUTOFF); ++k )
         {
            SCIP_ROWPREP* rowprep;
            SCIP_VAR* vars[3];
            SCIP_Real coefs[3];

            if( !SCIPisFeasLT(scip, eigenvals[k], 0.0) )
               continue;

            if( matrix[2*k] * matrix[2*k + 1] > 0.0 )
               printf("XXXX %g\n", matrix[2*k] * matrix[2*k + 1]);

            /* create rowprep */
            SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, FALSE) );

            /* cuts reads as v_1^2 X_ii + 2 v_1 v_2 X_ij + v_2^2 X_jj >= 0 */
            vars[0] = xx;
            coefs[0] = SQR(matrix[2*k]);
            vars[1] = xy;
            coefs[1] = 2.0 * matrix[2*k] * matrix[2*k + 1];
            vars[2] = yy;
            coefs[2] = SQR(matrix[2*k+1]);

            SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, 3, vars, coefs) );
            SCIPdebug( SCIPprintRowprep(scip, rowprep, NULL) );
            SCIPdebugMsg(scip, "cut violation %g minefficacy = %g\n", SCIPgetRowprepViolation(scip, rowprep, sol, NULL), minefficacy);

            /* check cut violation */
            if( SCIPgetRowprepViolation(scip, rowprep, sol, NULL) > minefficacy )
            {
               SCIP_ROW* row;
               SCIP_Bool infeasible;

               /* create, add, and release row */
               SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &row) );

               /* update result pointer */
               *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
            }

            /* free rowprep */
            SCIPfreeRowprep(scip, &rowprep);
         }
      }
   }

   /* TODO generate a/several cut/cuts */

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMinor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaMinor(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->minors == NULL);
   assert(sepadata->nminors == 0);
   assert(sepadata->minorssize == 0);

   /* free separator data */
   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitMinor NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitMinor NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMinor)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* clear separation data */
   SCIP_CALL( sepadataClear(scip, sepadata) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMinor)
{  /*lint --e{715}*/

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, sepa) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, NULL, allowlocal, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMinor)
{  /*lint --e{715}*/

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, sepa) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, sol, allowlocal, result) );

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the minor separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMinor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata = NULL;
   SCIP_SEPA* sepa = NULL;

   /* create minor separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpMinor, sepaExecsolMinor,
         sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMinor) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMinor) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitMinor) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitMinor) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolMinor) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolMinor) );

   /* add minor separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxminors",
         "maximum number for minors (0: no limit)",
         &sepadata->maxminors, FALSE, DEFAULT_MAXMINORS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
