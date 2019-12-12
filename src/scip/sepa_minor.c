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
 * @brief  principal minor separator
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
#define DEFAULT_MINCUTVIOL         1e-4 /**< default minimum required violation of a cut */
#define DEFAULT_RANDSEED            157 /**< default random seed */

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
   SCIP_Real             mincutviol;         /**< minimum required violation of a cut */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
};

/*
 * Local methods
 */

/** helper method to store a 2x2 minor in the separation data */
static
SCIP_RETCODE sepadataAddMinor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR*             x,                  /**< x variable */
   SCIP_VAR*             y,                  /**< y variable */
   SCIP_VAR*             auxvarxx,           /**< auxiliary variable for x*x */
   SCIP_VAR*             auxvaryy,           /**< auxiliary variable for y*y */
   SCIP_VAR*             auxvarxy            /**< auxiliary variable for x*y */
   )
{
   assert(sepadata != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);
   assert(auxvarxx != NULL);
   assert(auxvaryy != NULL);
   assert(auxvarxy != NULL);
   assert(auxvarxx != auxvaryy);
   assert(auxvarxx != auxvarxy);
   assert(auxvaryy != auxvarxy);

   SCIPdebugMsg(scip, "store 2x2 minor: %s %s %s for x=%s y=%y\n", SCIPvarGetName(auxvarxx), SCIPvarGetName(auxvaryy),
      SCIPvarGetName(auxvarxy), SCIPvarGetName(x), SCIPvarGetName(y));

   /* reallocate if necessary */
   if( sepadata->minorssize < 5 * (sepadata->nminors + 1) )
   {
      int newsize = SCIPcalcMemGrowSize(scip, 5 * (sepadata->nminors + 1));
      assert(newsize > 5 * (sepadata->nminors + 1));

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->minors), sepadata->minorssize, newsize) );
      sepadata->minorssize = newsize;
   }

   /* store minor */
   sepadata->minors[5 * sepadata->nminors] = x;
   sepadata->minors[5 * sepadata->nminors + 1] = y;
   sepadata->minors[5 * sepadata->nminors + 2] = auxvarxx;
   sepadata->minors[5 * sepadata->nminors + 3] = auxvaryy;
   sepadata->minors[5 * sepadata->nminors + 4] = auxvarxy;
   ++(sepadata->nminors);

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, x) );
   SCIP_CALL( SCIPcaptureVar(scip, y) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxx) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvaryy) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxy) );

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
   for( i = 0; i < 5 * sepadata->nminors; ++i )
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

/** helper method to get the variables associated to a minor */
static
SCIP_RETCODE getMinorVars(
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   int                   idx,                /**< index of the stored minor */
   SCIP_VAR**            x,                  /**< pointer to store x variable */
   SCIP_VAR**            y,                  /**< pointer to store x variable */
   SCIP_VAR**            auxvarxx,           /**< pointer to store auxiliary variable for x*x */
   SCIP_VAR**            auxvaryy,           /**< pointer to store auxiliary variable for y*y */
   SCIP_VAR**            auxvarxy            /**< pointer to store auxiliary variable for x*y */
   )
{
   assert(sepadata != NULL);
   assert(idx >= 0 && idx < sepadata->nminors);
   assert(auxvarxx != NULL);
   assert(auxvaryy != NULL);
   assert(auxvarxy != NULL);

   *x = sepadata->minors[5 * idx];
   *y = sepadata->minors[5 * idx + 1];
   *auxvarxx = sepadata->minors[5 * idx + 2];
   *auxvaryy = sepadata->minors[5 * idx + 3];
   *auxvarxy = sepadata->minors[5 * idx + 4];

   return SCIP_OKAY;
}

/** method to detect and store principal minors */
static
SCIP_RETCODE detectMinors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
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
            && SCIPisConsExprExprVar(children[0]) )
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
            && SCIPisConsExprExprVar(children[0]) && SCIPisConsExprExprVar(children[1]) )
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
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nbilinterms) );

      for( i = 0; i < nbilinterms; ++i )
         perm[i] = i;

      /* permute array */
      SCIPrandomPermuteIntArray(sepadata->randnumgen, perm, 0, nbilinterms);
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
         SCIP_CALL( sepadataAddMinor(scip, sepadata, x, y, auxvarxx, auxvaryy, auxvar) );
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

/** helper method to compute eigenvectors and eigenvalues */
static
SCIP_RETCODE getEigenValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x,                  /**< solution value of x */
   SCIP_Real             y,                  /**< solution value of y */
   SCIP_Real             xx,                 /**< solution value of x*x */
   SCIP_Real             yy,                 /**< solution value of y*y */
   SCIP_Real             xy,                 /**< solution value of x*y */
   SCIP_Real*            eigenvals,          /**< array to store eigenvalues (at least of size 3) */
   SCIP_Real*            eigenvecs,          /**< array to store eigenvalues (at least of size 9) */
   SCIP_Bool*            success             /**< pointer to store whether eigenvalue computation was successful */
   )
{
   assert(eigenvals != NULL);
   assert(eigenvecs != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* construct matrix */
   eigenvecs[0] = 1.0;
   eigenvecs[1] = x;
   eigenvecs[2] = y;
   eigenvecs[3] = x;
   eigenvecs[4] = xx;
   eigenvecs[5] = xy;
   eigenvecs[6] = y;
   eigenvecs[7] = xy;
   eigenvecs[8] = yy;

   /* use LAPACK to compute the eigenvalues and eigenvectors */
   if( LapackDsyev(TRUE, 3, eigenvecs, eigenvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors of augmented quadratic form matrix.\n");
      *success = FALSE;
   }

   return SCIP_OKAY;
}

/** helper generate and add a cut */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< solution to separate (might be NULL) */
   SCIP_VAR*             x,                  /**< x variable */
   SCIP_VAR*             y,                  /**< y variable */
   SCIP_VAR*             xx,                 /**< auxiliary variable for x*x */
   SCIP_VAR*             yy,                 /**< auxiliary variable for y*y */
   SCIP_VAR*             xy,                 /**< auxiliary variable for x*y */
   SCIP_Real*            eigenvec,           /**< array containing an eigenvector */
   SCIP_Real             eigenval,           /**< eigenvalue */
   SCIP_Real             mincutviol,         /**< minimal required violation */
   SCIP_RESULT*          result              /**< pointer to update the result */
   )
{
   SCIP_VAR* vars[5] = {x, y, xx, yy, xy};
   SCIP_Real coefs[5];
   SCIP_Real constant;
   SCIP_ROWPREP* rowprep;
   SCIP_Bool success;

   assert(x != NULL);
   assert(y != NULL);
   assert(xx != NULL);
   assert(yy != NULL);
   assert(xy != NULL);
   assert(eigenvec != NULL);
   assert(mincutviol >= 0.0);
   assert(result != NULL);

   /* check whether the resulting cut is violated enough */
   if( !SCIPisFeasLT(scip, eigenval, -mincutviol) )
      return SCIP_OKAY;

   /* the resulting cut reads as v_0^2 + 2v_0v_1 * x + 2v_0v_2 * y + v_1^2 * xx + v_2^2 * yy + 2v_1v_2 * xy */
   constant = SQR(eigenvec[0]);
   coefs[0] = 2.0 * eigenvec[0] * eigenvec[1];
   coefs[1] = 2.0 * eigenvec[0] * eigenvec[2];
   coefs[2] = SQR(eigenvec[1]);
   coefs[3] = SQR(eigenvec[2]);
   coefs[4] = 2.0 * eigenvec[1] * eigenvec[2];

   /* create rowprep */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, FALSE) );
   SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, 5, vars, coefs) );
   SCIPaddRowprepConstant(rowprep, constant);
   SCIPdebug( SCIPprintRowprep(scip, rowprep, NULL) );
   SCIPdebugMsg(scip, "cut violation %g mincutviol = %g\n", SCIPgetRowprepViolation(scip, rowprep, sol, NULL), mincutviol);

   /* cleanup coefficient and side, esp treat epsilon to integral values; don't consider scaling up here */
   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, SCIP_CONSEXPR_CUTMAXRANGE, 0.0, NULL, &success) );

   /* check cut violation */
   if( success && SCIPgetRowprepViolation(scip, rowprep, sol, NULL) > mincutviol )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;
      char name[SCIP_MAXSTRLEN];

      /* set name of rowprep */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "minor_%s_%s_%s", SCIPvarGetName(xx), SCIPvarGetName(yy),
         SCIPvarGetName(xy));
      memcpy(rowprep->name, name, (unsigned long)SCIP_MAXSTRLEN);

      /* create, add, and release row */
      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );

      /* update result pointer */
      *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   }

   /* free rowprep */
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** separates cuts for stored principal minors */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
   int i;

   assert(sepa != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* check whether there are some minors available */
   if( sepadata->nminors == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   for( i = 0; i < sepadata->nminors && (*result != SCIP_CUTOFF); ++i )
   {
      SCIP_Real eigenvals[3];
      SCIP_Real eigenvecs[9];
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_VAR* xx;
      SCIP_VAR* yy;
      SCIP_VAR* xy;
      SCIP_Real solx;
      SCIP_Real soly;
      SCIP_Real solxx;
      SCIP_Real solyy;
      SCIP_Real solxy;
      SCIP_Bool success;
      int k;

      /* get variables of the i-th minor */
      SCIP_CALL( getMinorVars(sepadata, i, &x, &y, &xx, &yy, &xy) );
      assert(x != NULL);
      assert(y != NULL);
      assert(xx != NULL);
      assert(yy != NULL);
      assert(xy != NULL);

      /* get current solution values */
      solx = SCIPgetSolVal(scip, sol, x);
      soly = SCIPgetSolVal(scip, sol, y);
      solxx = SCIPgetSolVal(scip, sol, xx);
      solyy = SCIPgetSolVal(scip, sol, yy);
      solxy = SCIPgetSolVal(scip, sol, xy);
      SCIPdebugMsg(scip, "solution values (x,y,xx,yy,xy)=(%g,%g,%g,%g,%g)\n", solx, soly, solxx, solyy, solxy);

      /* compute eigenvalues and eigenvectors */
      SCIP_CALL( getEigenValues(scip, solx, soly, solxx, solyy, solxy, eigenvals, eigenvecs, &success) );
      if( !success )
         continue;

      /* try to generate a cut for each negative eigenvalue */
      for( k = 0; k < 3 && (*result != SCIP_CUTOFF); ++k )
      {
         SCIPdebugMsg(scip, "eigenvalue = %g  eigenvector = (%g,%g,%g)\n", eigenvals[k], eigenvecs[3*k], eigenvecs[3*k + 1], eigenvecs[3*k + 2]);
         SCIP_CALL( addCut(scip, sepa, sol, x, y, xx, yy, xy, &eigenvecs[3*k], eigenvals[k], sepadata->mincutviol, result) );
         SCIPdebugMsg(scip, "result: %u\n", *result);
      }
   }

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
static
SCIP_DECL_SEPAINIT(sepaInitMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->randnumgen == NULL);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}


/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->randnumgen != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}


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
   SCIP_CALL( detectMinors(scip, SCIPsepaGetData(sepa)) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, NULL, result) );

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
   SCIP_CALL( detectMinors(scip, SCIPsepaGetData(sepa)) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, sol, result) );

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

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/mincutviol",
         "minimum required violation of a cut",
         &sepadata->mincutviol, FALSE, DEFAULT_MINCUTVIOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
