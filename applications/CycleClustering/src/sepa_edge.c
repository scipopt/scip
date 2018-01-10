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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_edge.c
 * @brief  edge-separator. Separates triangle-inequalities in cycle clustering problem
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "sepa_edge.h"

#include "probdata_cyc.h"

#define SEPA_NAME                "edge"
#define SEPA_DESC              "separator to separate triangle-inequalities in cycle-clustering application"
#define SEPA_PRIORITY              2000
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                     500
#define MAXROUNDS                    15

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyEdge)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaEdge(scip) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpEdge)
{  /*lint --e{715}*/
   SCIP_VAR**** edgevars;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Real** sign;
   SCIP_Real* violation;
   SCIP_ROW** cuts;
   int nbins;
   int ncuts;
   int i;
   int j;
   int k;
   int l;
   int rounds;
   int size;

   rounds = SCIPsepaGetNCallsAtNode(sepa);
   if( rounds >= MAXROUNDS )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   edgevars = SCIPcycGetEdgevars(scip);
   nbins = SCIPcycGetNBins(scip);

   assert(nbins > 0);
   assert(NULL != edgevars);

   ncuts = 0;
   size = MAXCUTS;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocClearMemoryArray(scip, &violation, size) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cuts, size) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign, 3) );
   for( i = 0; i < 3; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign[i], 3) ); /*lint !e866*/
      for( j = 0; j < 3; ++j )
      {
         sign[i][j] = 1.0;
      }
      sign[i][i] = -1.0;
   }

   if( SCIPcycGetNCluster(scip) > 3 )
   {
      /* separate edges by the valid inequality y_ij1 + y_ik1 - y_jk0 <= 1 and y_ij0+y_ik1-y_jk1 <= 1 */
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < nbins; ++j )
         {
            for( k = 0; k < nbins; ++k )
            {
               if( NULL == edgevars[i][j] || NULL == edgevars[i][k] || NULL == edgevars[j][k] )
                  continue;
               if( (i != j && i != k && j > k) )
               {
                  violation[ncuts] = SCIPvarGetLPSol(edgevars[i][j][1]) + SCIPvarGetLPSol(edgevars[i][k][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) - 1;

                  if( violation[ncuts] > 0)
                  {
                     (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                     SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                     SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][0], -1.0) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][j][1], 1.0) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][k][1], 1.0) );
                     SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncuts]) );

                     if( ncuts >= size - 1 )
                     {
                        SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                        SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                        size += MAXCUTS;
                     }

                     ncuts++;
                  }

                  violation[ncuts] = SCIPvarGetLPSol(edgevars[j][i][1]) + SCIPvarGetLPSol(edgevars[k][i][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) - 1;

                  if( violation[ncuts] > 0)
                  {
                     (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                     SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                     SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][0], -1.0) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][i][1], 1.0) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][i][1], 1.0) );
                     SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncuts]) );

                     if( ncuts >= size - 1 )
                     {
                        SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                        SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                        size += MAXCUTS;
                     }

                     ncuts++;
                  }
               }
            }
         }
      }

      /* separate edges by the valid inequality y_ij0 + y_ik0 - y_jk0 <= 1 */
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < i; ++j )
         {
            for( k = 0; k < j; ++k )
            {
               if( NULL == edgevars[i][j] || NULL == edgevars[j][k] || NULL == edgevars[i][k] || NULL == edgevars[k][i] || NULL == edgevars[k][j] || NULL == edgevars[j][i] )
                  continue;
               for( l = 0; l < 3; ++l )
               {
                  violation[ncuts] = sign[l][0] * SCIPvarGetLPSol(edgevars[i][j][0]) + sign[l][1] * SCIPvarGetLPSol(edgevars[i][k][0]) + sign[l][2] * SCIPvarGetLPSol(edgevars[j][k][0]) - 1;
                  violation[ncuts] += 0.5 * sign[l][0] * (SCIPvarGetLPSol(edgevars[i][j][1]) + SCIPvarGetLPSol(edgevars[j][i][1]));
                  violation[ncuts] += 0.5 * sign[l][1] * (SCIPvarGetLPSol(edgevars[i][k][1]) + SCIPvarGetLPSol(edgevars[k][i][1]));
                  violation[ncuts] += 0.5 * sign[l][2] * (SCIPvarGetLPSol(edgevars[k][j][1]) + SCIPvarGetLPSol(edgevars[j][k][1]));

                  if( violation[ncuts] > 0 )
                  {
                     (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d_%f%f%f", i, j, k, sign[l][0], sign[l][1], sign[l][2] );
                     SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][j][0], sign[l][0]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][j][1], 0.5 * sign[l][0]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][i][1], 0.5 * sign[l][0]) );

                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][k][0], sign[l][1]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][i][1], 0.5 * sign[l][1]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][k][1], 0.5 * sign[l][1]) );

                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][0], sign[l][2]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][j][1], 0.5 * sign[l][2]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][1], 0.5 * sign[l][2]) );

                     if( ncuts >= size - 1 )
                     {
                        SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                        SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                        size += MAXCUTS;
                     }

                     ncuts++;
                  }
               }
            }
         }
      }
   }

   if( SCIPcycGetNCluster(scip) == 3 )
   {
      for( i = 0; i < nbins; ++i )
      {
         for (j = 0; j < nbins; ++j )
         {
            for (k = 0; k < nbins; ++k )
            {
               if (NULL == edgevars[i][j] || NULL == edgevars[j][k] || NULL == edgevars[i][k])
                  continue;

               violation[ncuts] = SCIPvarGetLPSol(edgevars[i][j][1]) + SCIPvarGetLPSol(edgevars[j][k][1]) - SCIPvarGetLPSol(edgevars[k][i][1]) - 1;

               if( violation[ncuts] > 0 )
               {
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i, j, k);
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][j][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][i][1], -1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncuts]) );
                  if( ncuts >= size - 1 )
                  {
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                     size += MAXCUTS;
                  }

                  ncuts++;
               }

               violation[ncuts] = SCIPvarGetLPSol(edgevars[i][j][1]) + SCIPvarGetLPSol(edgevars[i][k][1]) +  SCIPvarGetLPSol(edgevars[j][k][1]) +  SCIPvarGetLPSol(edgevars[k][j][1]) - 2;

               if( violation[ncuts] > 0 )
               {
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i, j, k);
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 2.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][j][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][j][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[i][k][1], 1.0) );

                  SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncuts]) );

                  if( ncuts >= size - 1 )
                  {
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                     size += MAXCUTS;
                  }

                  ncuts++;
               }

               violation[ncuts] = SCIPvarGetLPSol(edgevars[j][i][1]) + SCIPvarGetLPSol(edgevars[k][i][1]) +  SCIPvarGetLPSol(edgevars[j][k][1]) +  SCIPvarGetLPSol(edgevars[k][j][1]) - 2;

               if( violation[ncuts] > 0 )
               {
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i, j, k);
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 2.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][i][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[j][k][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][j][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[k][i][1], 1.0) );

                  SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncuts]) );
                  if( ncuts >= size - 1 )
                  {
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                     size += MAXCUTS;
                  }

                  ncuts++;
               }
            }
         }
      }
   }

   /* apply the cuts with the highest violation */
   SCIPsortDownRealPtr(violation, ((void **) cuts), ncuts);
   for( i = 0; i < MIN(MAXCUTS, ncuts); ++i )
   {
      if( SCIPisCutEfficacious(scip, NULL, cuts[i]) )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cuts[i]) );
         *result = SCIP_SEPARATED;
      }
   }

   /* free memory */
   for( i = 0; i < ncuts; ++i )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(cuts[i])) );
   }
   SCIPfreeMemoryArray(scip, &cuts);
   SCIPfreeMemoryArray(scip, &violation);
   for( i = 0; i < 3; ++i )
   {
      SCIPfreeMemoryArray(scip, &sign[i]);
   }
   SCIPfreeMemoryArray(scip, &sign);

   return SCIP_OKAY;
}

/** creates the Edge separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaEdge(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPA* sepa;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpEdge, NULL,
      NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyEdge) );

   return SCIP_OKAY;
}
