/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_gomory.c,v 1.43 2005/02/08 14:22:30 bzfpfend Exp $"

/**@file   sepa_gomory.c
 * @brief  Gomory MIR Cuts
 * @author Tobias Achterberg
 */

/**@todo try k-Gomory-cuts (s. Cornuejols: K-Cuts: A Variation of Gomory Mixed Integer Cuts from the LP Tableau) */
/**@todo Gomory cuts don't seem to work with SOPLEX, but c-MIR cuts do! */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_gomory.h"


#define SEPA_NAME              "gomory"
#define SEPA_DESC              "Gomory MIR cuts separator"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                    10
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of gomory separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of gomory separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of gomory cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of gomory cuts separated per separation round in root node */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAXWEIGHTRANGE    1e+04 /**< maximal valid range max(|weights|)/min(|weights|) of row weights */

#define BOUNDSWITCH              0.9999
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE
#define MAKECONTINTEGRAL          FALSE
#define MINFRAC                    0.05


/** separator data */
struct SepaData
{
   Real             maxweightrange;     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   int              maxrounds;          /**< maximal number of gomory separation rounds per node (-1: unlimited) */
   int              maxroundsroot;      /**< maximal number of gomory separation rounds in the root node (-1: unlimited) */
   int              maxsepacuts;        /**< maximal number of gomory cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of gomory cuts separated per separation round in root node */
   Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};




/*
 * local methods
 */

/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm */
static
RETCODE storeCutInArrays(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of problem variables */
   VAR**            vars,               /**< problem variables */
   Real*            cutcoefs,           /**< dense coefficient vector */
   Real*            varsolvals,         /**< dense variable LP solution vector */
   char             normtype,           /**< type of norm to use for efficacy norm calculation */
   COL**            cutcols,            /**< array to store columns of sparse cut vector */
   Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*             cutlen,             /**< pointer to store number of nonzero entries in cut */
   Real*            cutact,             /**< pointer to store activity of cut */
   Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   Real val;
   Real absval;
   Real cutsqrnorm;
   Real act;
   Real norm;
   int len;
   int v;

   assert(nvars == 0 || cutcoefs != NULL);
   assert(nvars == 0 || varsolvals != NULL);
   assert(cutcols != NULL);
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
            assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN);
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutcols[len] = SCIPvarGetCol(vars[v]);
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
            assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN);
            act += val * varsolvals[v];
            absval = REALABS(val);
            norm = MAX(norm, absval);
            cutcols[len] = SCIPvarGetCol(vars[v]);
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
            assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN);
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutcols[len] = SCIPvarGetCol(vars[v]);
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
            assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN);
            act += val * varsolvals[v];
            norm = 1.0;
            cutcols[len] = SCIPvarGetCol(vars[v]);
            cutvals[len] = val;
            len++;
         }
      }
      break;
   default:
      errorMessage("invalid efficacy norm parameter '%c'\n", normtype);
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

/** destructor of separator to free user data (called when SCIP is exiting) */
static
DECL_SEPAFREE(sepaFreeGomory)
{  /*lint --e{715}*/
   SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** initialization method of separator (called when problem solving starts) */
#define sepaInitGomory NULL


/** deinitialization method of separator (called when problem solving exits) */
#define sepaExitGomory NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolGomory NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolGomory NULL


/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecGomory)
{  /*lint --e{715}*/
   SEPADATA* sepadata;
   VAR** vars;
   COL** cols;
   ROW** rows;
   Real* varsolvals;
   Real* binvrow;
   Real* cutcoefs;
   Real cutrhs;
   Real cutact;
   Real maxscale;
   Longint maxdnom;
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
   Bool success;
   Bool cutislocal;
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

   /* only call the gomory cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get variables data */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get LP data */
   CHECK_OKAY( SCIPgetLPColsData(scip, &cols, &ncols) );
   CHECK_OKAY( SCIPgetLPRowsData(scip, &rows, &nrows) );
   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

   /* get the type of norm to use for efficacy calculations */
   CHECK_OKAY( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* set the maximal denominator in rational representation of gomory cut and the maximal scale factor to
    * scale resulting cut to integral values to avoid numerical instabilities
    */
   /**@todo find better but still stable gomory cut settings: look at dcmulti, gesa3, khb0525, misc06, p2756 */
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
   CHECK_OKAY( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &basisind, nrows) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &binvrow, nrows) );
   varsolvals = NULL; /* allocate this later, if needed */

   /* get basis indices */
   CHECK_OKAY( SCIPgetLPBasisInd(scip, basisind) );

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
      maxsepacuts = sepadata->maxsepacutsroot;
   else
      maxsepacuts = sepadata->maxsepacuts;

   debugMessage("searching gomory cuts: %d cols, %d rows, maxdnom=%lld, maxscale=%g, maxcuts=%d\n",
      ncols, nrows, maxdnom, maxscale, maxsepacuts);

   /* for all basic columns belonging to integer variables, try to generate a gomory cut */
   ncuts = 0;
   for( i = 0; i < nrows && ncuts < maxsepacuts; ++i )
   {
      c = basisind[i];
      if( c >= 0 )
      {
         VAR* var;

         assert(c < ncols);
         var = SCIPcolGetVar(cols[c]);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            Real primsol;

            primsol = SCIPcolGetPrimsol(cols[c]);
            assert(SCIPgetVarSol(scip, var) == primsol); /*lint !e777*/

            if( SCIPfeasFrac(scip, primsol) >= MINFRAC )
            {
               debugMessage("trying gomory cut for <%s> [%g]\n", SCIPvarGetName(var), primsol);

               /* get the row of B^-1 for this basic integer variable with fractional solution value */
               CHECK_OKAY( SCIPgetLPBInvRow(scip, i, binvrow) );

               /* create a MIR cut out of the weighted LP rows using the B^-1 row as weights */
               CHECK_OKAY( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, sepadata->maxweightrange, MINFRAC,
                     binvrow, 1.0, cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
               assert(ALLOWLOCAL || !cutislocal);
               debugMessage(" -> success=%d: %g <= %g\n", success, cutact, cutrhs);

               /* if successful, convert dense cut into sparse row, and add the row as a cut */
               if( success && SCIPisFeasGT(scip, cutact, cutrhs) )
               {
                  COL** cutcols;
                  Real* cutvals;
                  Real cutnorm;
                  int cutlen;

                  /* if this is the first successful cut, get the LP solution for all COLUMN variables */
                  if( varsolvals == NULL )
                  {
                     int v;

                     CHECK_OKAY( SCIPallocBufferArray(scip, &varsolvals, nvars) );
                     for( v = 0; v < nvars; ++v )
                     {
                        if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
                           varsolvals[v] = SCIPvarGetLPSol(vars[v]);
                     }
                  }
                  assert(varsolvals != NULL);

                  /* get temporary memory for storing the cut as sparse row */
                  CHECK_OKAY( SCIPallocBufferArray(scip, &cutcols, nvars) );
                  CHECK_OKAY( SCIPallocBufferArray(scip, &cutvals, nvars) );

                  /* store the cut as sparse row, calculate activity and norm of cut */
                  CHECK_OKAY( storeCutInArrays(scip, nvars, vars, cutcoefs, varsolvals, normtype,
                        cutcols, cutvals, &cutlen, &cutact, &cutnorm) );

                  debugMessage(" -> gomory cut for <%s>: act=%f, rhs=%f, norm=%f, eff=%f\n",
                     SCIPvarGetName(var), cutact, cutrhs, cutnorm, (cutact - cutrhs)/cutnorm);

                  if( SCIPisPositive(scip, cutnorm) && SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm) )
                  {
                     ROW* cut;
                     char cutname[MAXSTRLEN];

                     /* create the cut */
                     sprintf(cutname, "gom%d_%d", SCIPgetNLPs(scip), c);
                     CHECK_OKAY( SCIPcreateRow(scip, &cut, cutname, cutlen, cutcols, cutvals, -SCIPinfinity(scip), cutrhs, 
                           cutislocal, FALSE, sepadata->dynamiccuts) );
                     debug(SCIPprintRow(scip, cut, NULL));

                     /* try to scale the cut to integral values */
                     CHECK_OKAY( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                           maxdnom, maxscale, MAKECONTINTEGRAL, &success) );

                     if( success )
                     {
                        if( !SCIPisCutEfficacious(scip, cut) )
                        {
                           debugMessage(" -> gomory cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                              cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                              SCIPgetCutEfficacy(scip, cut));
                           debug(SCIPprintRow(scip, cut, NULL));
                           success = FALSE;
                        }
                        else
                        {
                           debugMessage(" -> found gomory cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                              cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                              SCIPgetCutEfficacy(scip, cut), SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                              SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
                           debug(SCIPprintRow(scip, cut, NULL));
                           CHECK_OKAY( SCIPaddCut(scip, cut, FALSE) );
                           if( !cutislocal )
                           {
                              CHECK_OKAY( SCIPaddPoolCut(scip, cut) );
                           }
                           *result = SCIP_SEPARATED;
                           ncuts++;
                        }
                     }
                     else
                     {
                        debugMessage(" -> gomory cut <%s> couldn't be scaled to integral coefficients: act=%f, rhs=%f, norm=%f, eff=%f\n",
                           cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, cut));
                     }

                     /* release the row */
                     CHECK_OKAY( SCIPreleaseRow(scip, &cut) );
                  }

                  /* free temporary memory */
                  SCIPfreeBufferArray(scip, &cutvals);
                  SCIPfreeBufferArray(scip, &cutcols);
               }
            }
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArrayNull(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutcoefs);

   debugMessage("end searching gomory cuts: found %d cuts\n", ncuts);

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the Gomory MIR cut separator and includes it in SCIP */
RETCODE SCIPincludeSepaGomory(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SEPADATA* sepadata;

   /* create separator data */
   CHECK_OKAY( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_DELAY,
         sepaFreeGomory, sepaInitGomory, sepaExitGomory,
         sepaInitsolGomory, sepaExitsolGomory, sepaExecGomory,
         sepadata) );

   /* add separator parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/gomory/maxrounds",
         "maximal number of gomory separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/gomory/maxroundsroot",
         "maximal number of gomory separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/gomory/maxsepacuts",
         "maximal number of gomory cuts separated per separation round",
         &sepadata->maxsepacuts, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddIntParam(scip,
         "separating/gomory/maxsepacutsroot",
         "maximal number of gomory cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "separating/gomory/maxweightrange",
         "maximal valid range max(|weights|)/min(|weights|) of row weights",
         &sepadata->maxweightrange, DEFAULT_MAXWEIGHTRANGE, 1.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddBoolParam(scip,
         "separating/gomory/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

