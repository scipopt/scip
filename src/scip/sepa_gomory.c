/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_gomory.c
 * @brief  Gomory Cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_gomory.h"


#define SEPA_NAME         "gomory"
#define SEPA_DESC         "gomory cuts separator"
#define SEPA_PRIORITY     0
#define SEPA_FREQ         8



/*
 * Callback methods
 */

static
DECL_SEPAEXEC(SCIPsepaExecGomory)
{
   VAR** vars;
   COL** cols;
   ROW** rows;
   Real* varsol;
   Real* binvrow;
   Real* cutcoef;
   Real cutrhs;
   Bool success;
   Longint maxdnom;
   int* basisind;
   int nvars;
   int ncols;
   int nrows;
   int actdepth;
   int maxdepth;
   int c;
   int i;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

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

   /* set the maximal denominator in rational representation of gomory cut to avoid numerical instabilities */
   actdepth = SCIPgetActDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   if( actdepth == 0 )
      maxdnom = 1000000;
   else if( actdepth <= maxdepth/4 )
      maxdnom = 100;
   else if( actdepth <= maxdepth/2 )
      maxdnom = 10;
   else
      maxdnom = 1;

   debugMessage("searching gomory cuts: %d cols, %d rows, maxdnom=%lld\n", ncols, nrows, maxdnom);

   *result = SCIP_DIDNOTFIND;

   /* allocate temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &cutcoef, nvars) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &basisind, nrows) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &binvrow, nrows) );
   varsol = NULL; /* allocate this later, if needed */

   /* get basis indices */
   CHECK_OKAY( SCIPgetLPBasisInd(scip, basisind) );

   /* for all basic columns belonging to integer variables, try to generate a gomory cut */
   for( i = 0; i < nrows; ++i )
   {
      c = basisind[i];
      if( c >= 0 )
      {
         VAR* var;

         assert(c < ncols);
         var = SCIPcolGetVar(cols[c]);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINOUS )
         {
            Real primsol;

            primsol = SCIPcolGetPrimsol(cols[c]);
            assert(SCIPgetVarSol(scip, var) == primsol);

            if( !SCIPisIntegral(scip, primsol) )
            {
               /* get the row of B^-1 for this basic integer variable with fractional solution value */
               CHECK_OKAY( SCIPgetLPBInvRow(scip, i, binvrow) );

               /* create a MIR cut out of the weighted LP rows using the B^-1 row as weights */
               CHECK_OKAY( SCIPcalcMIR(scip, 0.05, binvrow, cutcoef, &cutrhs, &success) );

               /* if successful, convert dense cut into sparse row, and add the row as a cut */
               if( success )
               {
                  COL** cutcols;
                  Real* cutvals;
                  Real cutact;
                  Real cutsqrnorm;
                  Real cutnorm;
                  Real val;
                  int cutlen;
                  int v;

                  /* if this is the first successful cut, get the LP solution for all COLUMN variables */
                  if( varsol == NULL )
                  {
                     CHECK_OKAY( SCIPcaptureBufferArray(scip, &varsol, nvars) );
                     for( v = 0; v < nvars; ++v )
                     {
                        if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
                           varsol[v] = SCIPvarGetLPSol(vars[v]);
                     }
                  }
                  assert(varsol != NULL);

                  /* get temporary memory for storing the cut as sparse row */
                  CHECK_OKAY( SCIPcaptureBufferArray(scip, &cutcols, nvars) );
                  CHECK_OKAY( SCIPcaptureBufferArray(scip, &cutvals, nvars) );

                  /* store the cut as sparse row, calculate activity of cut */
                  cutlen = 0;
                  cutact = 0.0;
                  cutsqrnorm = 0.0;
                  for( v = 0; v < nvars; ++v )
                  {
                     val = cutcoef[v];
                     if( !SCIPisZero(scip, val) )
                     {
                        assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN);
                        cutact += val * varsol[v];
                        cutsqrnorm += SQR(val);
                        cutcols[cutlen] = SCIPvarGetCol(vars[v]);
                        cutvals[cutlen] = val;
                        cutlen++;
                     }
                  }
                  cutnorm = SQRT(cutsqrnorm);

                  if( SCIPisPositive(scip, cutnorm)
                     && SCIPisFeasGT(scip, cutact, cutrhs)
                     && SCIPisCutViolated(scip, cutact/cutnorm, cutrhs/cutnorm) )
                  {
                     ROW* cut;
                     char cutname[255];
                     Bool success;

                     /* create the cut */
                     sprintf(cutname, "gom%d_%d", SCIPgetNLPs(scip), c);
                     CHECK_OKAY( SCIPcreateRow(scip, &cut, cutname, cutlen, cutcols, cutvals, -SCIPinfinity(scip), cutrhs, 
                                 TRUE, FALSE, TRUE) );
                     /*debugMessage(" -> found potential gomory cut <%s>: activity=%f, rhs=%f:", cutname, cutact, cutrhs);
                       debug(SCIProwPrint(cut, NULL));*/

                     /* try to scale the cut to integral values */
                     CHECK_OKAY( SCIPmakeRowRational(scip, cut, maxdnom, &success) );
                     
                     /* if scaling was successful, add the cut */
                     if( success )
                     {
                        cutact = SCIPgetRowLPActivity(scip, cut);
                        cutrhs = SCIProwGetRhs(cut);
                        cutnorm = SCIProwGetNorm(cut);
                        if( SCIPisPositive(scip, cutnorm)
                           && SCIPisFeasGT(scip, cutact, cutrhs)
                           && SCIPisCutViolated(scip, cutact/cutnorm, cutrhs/cutnorm) )
                        {
                           debugMessage(" -> found integral gomory cut <%s>: activity=%f, rhs=%f:",
                              cutname, cutact, cutrhs);
                           debug(SCIProwPrint(cut, NULL));
                           CHECK_OKAY( SCIPaddCut(scip, cut, (cutact-cutrhs)/cutnorm/(cutlen+1)) );
                           *result = SCIP_SEPARATED;
                        }
                     }

                     /* release the row */
                     CHECK_OKAY( SCIPreleaseRow(scip, &cut) );
                  }

                  /* free temporary memory */
                  CHECK_OKAY( SCIPreleaseBufferArray(scip, &cutvals) );
                  CHECK_OKAY( SCIPreleaseBufferArray(scip, &cutcols) );
               }
            }
         }
      }
   }

   /* free temporary memory */
   if( varsol != NULL )
   {
      CHECK_OKAY( SCIPreleaseBufferArray(scip, &varsol) );
   }
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &binvrow) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &basisind) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &cutcoef) );

   debugMessage("end searching gomory cuts\n");

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the gomory separator and includes it in SCIP */
RETCODE SCIPincludeSepaGomory(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ,
                  NULL, NULL, NULL, SCIPsepaExecGomory,
                  NULL) );

   return SCIP_OKAY;
}

