
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

/**@file   sepa_partition.c
 * @brief  partition-separator. Separates triangle-inequalities in SparseApprox Problem
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "sepa_partition.h"
#include "probdata_spa.h"
#include "scip/cons_linear.h"

#define SEPA_NAME              "partition"
#define SEPA_DESC              "separator to separate triangle-inequalities in cycle-clustering application"
#define SEPA_PRIORITY              1500
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                     500

static
SCIP_RETCODE createPartitionCut(
   SCIP*                 scip,
   SCIP_SEPA*            sepa,
   SCIP_VAR****          edgevars,
   int*                  S,
   int*                  T,
   int                   nS,
   int                   nT,
   SCIP_RESULT*          result,
   int*                  ncuts
)
{
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Real violation;
   int i;
   int j;

   /* sort the pointers */
   SCIPsortInt(S, nS);
   SCIPsortInt(T, nT);
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "PartitionCut_%d_%d", nS, nT);
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), MIN(nS,nT), FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
   violation = MIN(nS,nT);
   for( i = 0; i < nS; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[S[i]][S[j]][0], -1.0) );
         violation += SCIPvarGetLPSol(edgevars[S[i]][S[j]][0]);
      }
   }
   for( i = 0; i < nT; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[T[i]][T[j]][0], -1.0) );
         violation += SCIPvarGetLPSol(edgevars[T[i]][T[j]][0]);
      }
   }
   for( i = 0; i < nS; ++i )
   {
      for( j = 0; j < nT; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[S[i]][T[j]][1], 1.0) );
         violation -= SCIPvarGetLPSol(edgevars[S[i]][T[j]][1]);
      }
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
   if( SCIPisCutEfficacious(scip, NULL, cut) )
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      *result = SCIP_SEPARATED;
      *ncuts += 1;
      SCIPdebug(SCIPprintRow(scip, cut, NULL));
      SCIPdebug(SCIPdebugMessagePrint(scip, "violation: %f \n", violation));

   }
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   return SCIP_OKAY;
}

/** copy method for separator plugins (called when SCIP copies plugins) */

static
SCIP_DECL_SEPACOPY(sepaCopyPartition)
{   /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaPartition(scip) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpPartition)
{
   SCIP_VAR**** edgevars;

   int nbins;
   int i;
   int j;
   int k;
   int l,l2;
   int rounds;
   int ncuts;
   SCIP_Real* fractionality;
   int*       idx;
   SCIP_Real triangleval;
   SCIP_Real violation;

   rounds = SCIPsepaGetNCallsAtNode(sepa);
   ncuts = 0;
   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPspaGetNrBins(scip);

   assert(nbins > 0);
   assert(NULL != edgevars);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocMemoryArray(scip, &fractionality, nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &idx, nbins) );
   for( i = 0; i < nbins; ++i )
   {
      idx[i] = i;
      fractionality[i] = 0;
      for( j = 0; j < nbins; ++j )
      {
         if( NULL == edgevars[i][j] )
            continue;
         fractionality[i] += MIN(SCIPvarGetLPSol(edgevars[i][j][1]), 1 - SCIPvarGetLPSol(edgevars[i][j][1])) + MIN(1 - SCIPvarGetLPSol(edgevars[MAX(i,j)][MIN(i,j)][0]), SCIPvarGetLPSol(edgevars[MAX(i,j)][MIN(i,j)][0]));
      }
   }
   SCIPsortDownRealInt(fractionality, idx, nbins);
   /* we try to construct partition inequalities from triangle-inequalities that are tight*/
   for( i = 0; i < nbins && ncuts < MAXCUTS; ++i )
   {
      int candidate = idx[i];
      for( j = 0; j < nbins && ncuts < MAXCUTS; ++j )
      {
         for( k = 0; k < nbins && ncuts < MAXCUTS; ++k )
         {
            int S[nbins];
            int T[nbins];
            int nS = 0;
            int nT = 0;
            SCIP_Real addval;
            SCIP_Real bestval;
            if( NULL == edgevars[candidate][j] || NULL == edgevars[candidate][k] || NULL == edgevars[j][k] || candidate == j || i ==k || j == k )
               continue;
            if( (candidate != j && candidate != k && j > k) )
            {
               nS = 1;
               S[0] = candidate;
               nT = 2;
               T[0] = j;
               T[1] = k;
               triangleval = SCIPvarGetLPSol(edgevars[candidate][j][1]) + SCIPvarGetLPSol(edgevars[candidate][k][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) - 1;
               if( SCIPisLE(scip, triangleval, 0.0) )
               {
                  /*                   add a bin to T*/
                  bestval = -SCIPinfinity(scip);
                  T[2] = -1;
                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[candidate][l] || NULL == edgevars[MAX(l,j)][MIN(l,j)] || NULL == edgevars[MAX(l,k)][MIN(l,k)] )
                        continue;

                     addval = SCIPvarGetLPSol(edgevars[candidate][l][1]) - SCIPvarGetLPSol(edgevars[MAX(l,j)][MIN(l,j)][0]) - SCIPvarGetLPSol(edgevars[MAX(l,k)][MIN(l,k)][0]);
                     if( addval > bestval )
                     {
                        bestval = addval;
                        T[2] = l;
                     }
                  }
                  nT++;
                  if( bestval > 0 )
                     SCIP_CALL( createPartitionCut(scip, sepa, edgevars, S, T, nS, nT, result, &ncuts) );
                  bestval = -SCIPinfinity(scip);
                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[MAX(candidate,l)][MIN(candidate,l)] || NULL == edgevars[l][T[0]] || NULL == edgevars[l][T[1]] || NULL == edgevars[l][T[2]] )
                        continue;
                     addval = -SCIPvarGetLPSol(edgevars[MAX(candidate,l)][MIN(candidate,l)][0]) + SCIPvarGetLPSol(edgevars[l][T[0]][1]) + SCIPvarGetLPSol(edgevars[l][T[1]][1]) + SCIPvarGetLPSol(edgevars[l][T[1]][1]);
                     if( addval > bestval)
                     {
                        bestval = addval;
                        S[1] = l;
                     }
                  }
               }
               nS = 2;
               S[0] = k;
               S[1] = j;
               nT = 1;
               T[0] = candidate;
               triangleval = SCIPvarGetLPSol(edgevars[j][candidate][1]) + SCIPvarGetLPSol(edgevars[k][candidate][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) - 1;
               if( SCIPisLE(scip, triangleval, 0.0) )
               {
                  /*                   add a bin to S*/
                  bestval = -SCIPinfinity(scip);
                  S[2] = -1;
                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[l][candidate] || NULL == edgevars[MAX(l,j)][MIN(l,j)] || NULL == edgevars[MAX(l,k)][MIN(l,k)] )
                        continue;

                     addval = SCIPvarGetLPSol(edgevars[l][candidate][1]) - SCIPvarGetLPSol(edgevars[MAX(l,j)][MIN(l,j)][0]) - SCIPvarGetLPSol(edgevars[MAX(l,k)][MIN(l,k)][0]);
                     if( addval > bestval )
                     {
                        bestval = addval;
                        S[2] = l;
                     }
                  }
                  nS++;
                  if( bestval > 0 )
                     SCIP_CALL( createPartitionCut(scip, sepa, edgevars, S, T, nS, nT, result, &ncuts) );
                  bestval = -SCIPinfinity(scip);
                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[MAX(l,candidate)][MIN(l,candidate)] || NULL == edgevars[l][S[0]] || NULL == edgevars[l][S[1]] || NULL == edgevars[l][S[2]] )
                        continue;
                     addval = -SCIPvarGetLPSol(edgevars[MAX(l,candidate)][MIN(l,candidate)][0]) + SCIPvarGetLPSol(edgevars[S[0]][l][1]) + SCIPvarGetLPSol(edgevars[S[1]][l][1]) + SCIPvarGetLPSol(edgevars[S[2]][l][1]);
                     if( addval > bestval)
                     {
                        bestval = addval;
                        T[1] = l;
                     }
                  }
                  if( bestval > 0 )
                     SCIP_CALL( createPartitionCut(scip, sepa, edgevars, S, T, nS, nT, result, &ncuts) );
               }
            }
         }
      }
   }

   SCIPfreeMemoryArray(scip, &fractionality);
   SCIPfreeMemoryArray(scip, &idx);

   return SCIP_OKAY;
}


/** creates the Partition separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaPartition(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPA* sepa;


   /* include separator */

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpPartition, NULL,
      NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyPartition) );


   return SCIP_OKAY;
}
