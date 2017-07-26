
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_model.c
 * @brief  model-separator. Separates triangle-inequalities in SparseApprox Problem
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "../../SparseApprox/src/sepa_model.h"

#include <assert.h>
#include <string.h>

#include "probdata_spa.h"
#include "scip/cons_logicor.h"

#define SEPA_NAME              "model"
#define SEPA_DESC              "separator to separate triangle-inequalities in cycle-clustering application"
#define SEPA_PRIORITY          53687091
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

/** copy method for separator plugins (called when SCIP copies plugins) */

static
SCIP_DECL_SEPACOPY(sepaCopyModel)
{   /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaModel(scip) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpModel)
{
   SCIP_VAR**** edgevars;
   SCIP_VAR*** binvars;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   int nbins;
   int ncluster;
   int i;
   int j;
   int c1;
   int c2;
   int l;
   char consname[SCIP_MAXSTRLEN];
   int** sign;
   SCIP_Real violation;

   edgevars = SCIPspaGetEdgevars(scip);
   binvars = SCIPspaGetBinvars(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);

   assert(nbins > 0);
   assert(ncluster > 0 && ncluster < nbins);
   assert(NULL != edgevars);

   *result = SCIP_DIDNOTFIND;
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign, 3) );
   for( i = 0; i < 3; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign[i], 3) );
      for( j = 0; j < 3; ++j )
      {
         sign[i][j] = 1;
      }
      sign[i][i] = -1;
   }

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == edgevars[i][j] )
            continue;
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            for( l = 0; l < 3; ++l )
            {
               violation = sign[l][0] * SCIPvarGetLPSol(edgevars[i][j][0]) + sign[l][1] * SCIPvarGetLPSol(binvars[i][c1]) + sign[l][2] * SCIPvarGetLPSol(binvars[j][c1]) - 1;
               if(SCIPisPositive(scip, violation) )
               {
                  (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_var%d", i, j, c1, l);
                  SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, consname, 0, NULL, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                  if( sign[l][0] > 0 )
                     var = SCIPvarGetNegatedVar(edgevars[i][j][0] );
                  else
                     var = edgevars[i][j][0];
                  SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );

                  if( sign[l][1] > 0 )
                     var = SCIPvarGetNegatedVar(binvars[i][c1]);
                  else
                     var = binvars[i][c1];
                  SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );

                  if( sign[l][2] > 0 )
                     var = SCIPvarGetNegatedVar(binvars[j][c1]);
                  else
                     var = binvars[j][c1];
                  SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );

                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               for( c2 = 0; c2 < ncluster; ++c2 )
               {
                  SCIP_VAR* edge;
                  if( c2 == c1 + 1 || ( c2 == 0 && c1 == ncluster -1) )
                     edge = edgevars[i][j][1];
                  else if( c2 == c1 - 1 || ( c1 == 0 && c2 == ncluster -1) )
                     edge = edgevars[j][i][1];
                  else
                     continue;

                  violation = sign[l][0] * SCIPvarGetLPSol(edge) + sign[l][1] * SCIPvarGetLPSol(binvars[i][c1]) + sign[l][2] * SCIPvarGetLPSol(binvars[j][c2]) - 1;
                  if( SCIPisPositive(scip, violation) )
                  {
                     (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_%d_var%d", i, j, c1, c2, l);
                     SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, consname, 0, NULL, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                     if( sign[l][0] > 0 )
                        var = SCIPvarGetNegatedVar(edge);
                     else
                        var = edge;
                     SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );

                     if( sign[l][1] > 0 )
                        var = SCIPvarGetNegatedVar(binvars[i][c1]);
                     else
                        var = binvars[i][c1];
                     SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );

                     if( sign[l][2] > 0 )
                        var = SCIPvarGetNegatedVar(binvars[j][c2]);
                     else
                        var = binvars[j][c2];
                     SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );

                     SCIP_CALL( SCIPaddCons(scip, cons) );
                     SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                  }
               }
            }
         }
      }
   }

   for( i = 0; i < 3; ++i )
   {
      SCIPfreeMemoryArray(scip, &sign[i]);
   }
   SCIPfreeMemoryArray(scip, &sign);
   return SCIP_OKAY;
}


/** creates the Model separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaModel(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPA* sepa;


   /* include separator */

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpModel, NULL,
      NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyModel) );


   return SCIP_OKAY;
}
