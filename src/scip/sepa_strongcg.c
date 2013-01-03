/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_strongcg.c
 * @brief  Strong CG Cuts (Letchford & Lodi)
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_strongcg.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "strongcg"
#define SEPA_DESC              "Strong CG cuts separator (Letchford and Lodi)"
#define SEPA_PRIORITY             -2000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of strong CG separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        20 /**< maximal number of strong CG separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of strong CG cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of strong CG cuts separated per separation round in root node */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAXWEIGHTRANGE    1e+04 /**< maximal valid range max(|weights|)/min(|weights|) of row weights */

#define MAKECUTINTEGRAL        /* try to scale all cuts to integral coefficients */
/*#define MAKEINTCUTINTEGRAL*/     /* try to scale cuts without continuous variables to integral coefficients */
#define FORCECUTINTEGRAL       /* discard cut if conversion to integral coefficients failed */
#define SEPARATEROWS           /* separate rows with integral slack */

#define BOUNDSWITCH              0.9999
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE
#define MAKECONTINTEGRAL          FALSE
#define MINFRAC                    0.05
#define MAXFRAC                    0.95

#define MAXAGGRLEN(nvars)          (0.1*(nvars)+1000) /**< maximal length of base inequality */


/** separator data */
struct SCIP_SepaData
{
   SCIP_Real             maxweightrange;     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   int                   maxrounds;          /**< maximal number of strong CG separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of strong CG separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of strong CG cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of strong CG cuts separated per separation round in root node */
   int                   lastncutsfound;     /**< total number of cuts found after last call of separator */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};


/*
 * local methods
 */

/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   SCIP_VAR**            cutvars,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real val;
   SCIP_Real absval;
   SCIP_Real cutsqrnorm;
   SCIP_Real act;
   SCIP_Real norm;
   int len;
   int v;

   assert(nvars == 0 || cutcoefs != NULL);
   assert(nvars == 0 || varsolvals != NULL);
   assert(cutvars != NULL);
   assert(cutvals != NULL);
   assert(cutlen != NULL);
   assert(cutact != NULL);
   assert(cutnorm != NULL);

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch( normtype )
   {
   case 'e':
      cutsqrnorm = 0.0;
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            absval = ABS(val);
            norm = MAX(norm, absval);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 's':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += ABS(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 'd':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm = 1.0;
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyStrongcg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaStrongcg(scip) );
 
   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeStrongcg)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpStrongcg)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* varsolvals;
   SCIP_Real* binvrow;
   SCIP_Real* cutcoefs;
   SCIP_Real cutrhs;
   SCIP_Real cutact;
   SCIP_Real maxscale;
   SCIP_Longint maxdnom;
   int* basisind;
   int nvars;
   int ncols;
   int nrows;
   int ncalls;
   int depth;
   int maxdepth;
   int maxsepacuts;
   int ncuts;
   int c;
   int i;
   SCIP_Bool success;
   SCIP_Bool cutislocal;
   char normtype;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call the strong CG cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if the LP solution is basic */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

#if 0 /* if too many columns, separator is usually very slow: delay it until no other cuts have been found */
   if( ncols >= 50*nrows )
      return SCIP_OKAY;
   if( ncols >= 5*nrows )
   {
      int ncutsfound;

      ncutsfound = SCIPgetNCutsFound(scip);
      if( ncutsfound > sepadata->lastncutsfound || !SCIPsepaWasLPDelayed(sepa) )
      {
         sepadata->lastncutsfound = ncutsfound;
         *result = SCIP_DELAYED;
         return SCIP_OKAY;
      }
   }
#endif

   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* set the maximal denominator in rational representation of strong CG cut and the maximal scale factor to
    * scale resulting cut to integral values to avoid numerical instabilities
    */
   /**@todo find better but still stable strong CG cut settings: look at dcmulti, gesa3, khb0525, misc06, p2756 */
   maxdepth = SCIPgetMaxDepth(scip);
   if( depth == 0 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/4 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/2 )
   {
      maxdnom = 100;
      maxscale = 100.0;
   }
   else
   {
      maxdnom = 10;
      maxscale = 10.0;
   }

   *result = SCIP_DIDNOTFIND;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   varsolvals = NULL; /* allocate this later, if needed */

   /* get basis indices */
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
      maxsepacuts = sepadata->maxsepacutsroot;
   else
      maxsepacuts = sepadata->maxsepacuts;

   SCIPdebugMessage("searching strong CG cuts: %d cols, %d rows, maxdnom=%"SCIP_LONGINT_FORMAT", maxscale=%g, maxcuts=%d\n",
      ncols, nrows, maxdnom, maxscale, maxsepacuts);

   /* for all basic columns belonging to integer variables, try to generate a strong CG cut */
   ncuts = 0;
   for( i = 0; i < nrows && ncuts < maxsepacuts && !SCIPisStopped(scip); ++i )
   {
      SCIP_Bool tryrow;

      tryrow = FALSE;
      c = basisind[i];
      if( c >= 0 )
      {
         SCIP_VAR* var;

         assert(c < ncols);
         var = SCIPcolGetVar(cols[c]);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_Real primsol;

            primsol = SCIPcolGetPrimsol(cols[c]);
            assert(SCIPgetVarSol(scip, var) == primsol); /*lint !e777*/

            if( SCIPfeasFrac(scip, primsol) >= MINFRAC )
            {
               SCIPdebugMessage("trying strong CG cut for col <%s> [%g]\n", SCIPvarGetName(var), primsol);
               tryrow = TRUE;
            }
         }
      }
#ifdef SEPARATEROWS
      else
      {
         SCIP_ROW* row;

         assert(0 <= -c-1 && -c-1 < nrows); 
         row = rows[-c-1];
         if( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) )
         {
            SCIP_Real primsol;

            primsol = SCIPgetRowActivity(scip, row);
            if( SCIPfeasFrac(scip, primsol) >= MINFRAC )
            {
               SCIPdebugMessage("trying strong CG cut for row <%s> [%g]\n", SCIProwGetName(row), primsol);
               tryrow = TRUE;
            }
         }
      }
#endif

      if( tryrow )
      {
         /* get the row of B^-1 for this basic integer variable with fractional solution value */
         SCIP_CALL( SCIPgetLPBInvRow(scip, i, binvrow) );

#ifdef SCIP_DEBUG
         /* initialize variables, that might not have been initialized in SCIPcalcMIR if success == FALSE */
         cutact = 0.0;
         cutrhs = SCIPinfinity(scip);
#endif
         /* create a strong CG cut out of the weighted LP rows using the B^-1 row as weights */
         SCIP_CALL( SCIPcalcStrongCG(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, (int) MAXAGGRLEN(nvars), sepadata->maxweightrange, MINFRAC, MAXFRAC,
               binvrow, 1.0, cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
         assert(ALLOWLOCAL || !cutislocal);
         SCIPdebugMessage(" -> success=%u: %g <= %g\n", success, cutact, cutrhs);
               
         /* if successful, convert dense cut into sparse row, and add the row as a cut */
         if( success && SCIPisFeasGT(scip, cutact, cutrhs) )
         {
            SCIP_VAR** cutvars;
            SCIP_Real* cutvals;
            SCIP_Real cutnorm;
            int cutlen;

            /* if this is the first successful cut, get the LP solution for all COLUMN variables */
            if( varsolvals == NULL )
            {
               int v;

               SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
               for( v = 0; v < nvars; ++v )
               {
                  if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
                     varsolvals[v] = SCIPvarGetLPSol(vars[v]);
               }
            }
            assert(varsolvals != NULL);

            /* get temporary memory for storing the cut as sparse row */
            SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, nvars) );

            /* store the cut as sparse row, calculate activity and norm of cut */
            SCIP_CALL( storeCutInArrays(scip, nvars, vars, cutcoefs, varsolvals, normtype,
                  cutvars, cutvals, &cutlen, &cutact, &cutnorm) );

            SCIPdebugMessage(" -> strong CG cut for <%s>: act=%f, rhs=%f, norm=%f, eff=%f\n",
               c >= 0 ? SCIPvarGetName(SCIPcolGetVar(cols[c])) : SCIProwGetName(rows[-c-1]),
               cutact, cutrhs, cutnorm, (cutact - cutrhs)/cutnorm);

            if( SCIPisPositive(scip, cutnorm) && SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm) )
            {
               SCIP_ROW* cut;
               char cutname[SCIP_MAXSTRLEN];

               /* create the cut */
               if( c >= 0 )
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "scg%d_x%d", SCIPgetNLPs(scip), c);
               else
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "scg%d_s%d", SCIPgetNLPs(scip), -c-1);
               SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );
               SCIP_CALL( SCIPaddVarsToRow(scip, cut, cutlen, cutvars, cutvals) );
               /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/

               assert(success);
#ifdef MAKECUTINTEGRAL
               /* try to scale the cut to integral values */
               SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                     maxdnom, maxscale, MAKECONTINTEGRAL, &success) );
#else
#ifdef MAKEINTCUTINTEGRAL
               /* try to scale the cut to integral values if there are no continuous variables
                *  -> leads to an integral slack variable that can later be used for other cuts
                */
               {
                  int k;
                  for( k = 0; k < cutlen && SCIPvarIsIntegral(cutvars[k]); ++k )
                  {}
                  if( k == cutlen )
                  {
                     SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                           maxdnom, maxscale, MAKECONTINTEGRAL, &success) );
                  }
               }
#endif
#endif

#ifndef FORCECUTINTEGRAL
               success = TRUE;
#endif

               if( success )
               {
                  if( !SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIPdebugMessage(" -> strong CG cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                        cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                        SCIPgetCutEfficacy(scip, NULL, cut));
                     /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/
                     success = FALSE;
                  }
                  else
                  {
                     SCIPdebugMessage(" -> found strong CG cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                        cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                        SCIPgetCutEfficacy(scip, NULL, cut),
                        SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                        SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
                     /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE) );
                     if( !cutislocal )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                     }
                     *result = SCIP_SEPARATED;
                     ncuts++;
                  }
               }
               else
               {
                  SCIPdebugMessage(" -> strong CG cut <%s> couldn't be scaled to integral coefficients: act=%f, rhs=%f, norm=%f, eff=%f\n",
                     cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, NULL, cut));
               }

               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &cutvals);
            SCIPfreeBufferArray(scip, &cutvars);
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArrayNull(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutcoefs);

   SCIPdebugMessage("end searching strong CG cuts: found %d cuts\n", ncuts);

   sepadata->lastncutsfound = SCIPgetNCutsFound(scip);

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the Strong CG cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaStrongcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
   sepadata->lastncutsfound = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpStrongcg, NULL,
         sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyStrongcg) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeStrongcg) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxrounds",
         "maximal number of strong CG separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxroundsroot",
         "maximal number of strong CG separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxsepacuts",
         "maximal number of strong CG cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxsepacutsroot",
         "maximal number of strong CG cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/strongcg/maxweightrange",
         "maximal valid range max(|weights|)/min(|weights|) of row weights",
         &sepadata->maxweightrange, TRUE, DEFAULT_MAXWEIGHTRANGE, 1.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/strongcg/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

