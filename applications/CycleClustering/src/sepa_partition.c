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

/**@file   sepa_partition.c
 * @brief  partition-separator. Searches for two partitions of size 2 and 3 (extension of triangle-inequalities).
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "sepa_partition.h"

#include "probdata_cyc.h"
#include "scip/cons_linear.h"

#define SEPA_NAME           "partition"
#define SEPA_DESC              "separator to separate triangle-inequalities in cycle-clustering application"
#define SEPA_PRIORITY              1500
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                     500
#define MAXROUNDS                    15


/** Given two partitions S, T creates the corresponding cut and adds it do SCIP */
static
SCIP_RETCODE createPartitionCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< Separator */
   SCIP_VAR****          edgevars,           /**< The edge-variables */
   int*                  firstpart,          /**< The first partition */
   int*                  secondpart,         /**< The second partition */
   int                   nfirst,             /**< Number of states in first partition */
   int                   nsecond,            /**< Number of states in second partition */
   SCIP_RESULT*          result,             /**< Result pointer */
   int*                  ncuts               /**< Number of generated cuts */
   )
{
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];
   int i;
   int j;

   /* sort the pointers */
   SCIPsortInt(firstpart, nfirst);
   SCIPsortInt(secondpart, nsecond);

   /* create cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "PartitionCut_%d_%d", nfirst, nsecond);
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), (SCIP_Real) MIN(nfirst, nsecond), FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

   for( i = 0; i < nfirst; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[firstpart[i]][firstpart[j]][0], -1.0) );
      }
   }
   for( i = 0; i < nsecond; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[secondpart[i]][secondpart[j]][0], -1.0) );
      }
   }
   for( i = 0; i < nfirst; ++i )
   {
      for( j = 0; j < nsecond; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[firstpart[i]][secondpart[j]][1], 1.0) );
      }
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

   if( SCIPisCutEfficacious(scip, NULL, cut) )
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      *result = SCIP_SEPARATED;
      *ncuts += 1;
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyPartition)
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
   SCIP_VAR**** edgevars;
   SCIP_Real* fractionality;
   int*       idx;
   SCIP_Real triangleval;
   SCIP_Real addval;
   SCIP_Real bestval;
   SCIP_Real fval;
   SCIP_Real incval;
   int* firstpart;
   int* secondpart;
   int nfirst = 0;
   int nsecond = 0;
   int nbins;
   int i;
   int j;
   int k;
   int l;
   int ncuts;
   int rounds;

   ncuts = 0;
   edgevars = SCIPcycGetEdgevars(scip);
   nbins = SCIPcycGetNBins(scip);
   rounds = SCIPsepaGetNCallsAtNode(sepa);

   assert(nbins > 0);
   assert(NULL != edgevars);

   *result = SCIP_DIDNOTFIND;

   if( rounds >= MAXROUNDS )
      {
         *result =  SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

   SCIP_CALL( SCIPallocMemoryArray(scip, &fractionality, nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &idx, nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &firstpart, nbins) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &secondpart, nbins) );

   for( i = 0; i < nbins; ++i )
   {
      idx[i] = i;
      fractionality[i] = 0;

      for( j = 0; j < nbins; ++j )
      {
         if( NULL == edgevars[i][j] )
            continue;
         fval = SCIPvarGetLPSol(edgevars[i][j][1]);
         incval = SCIPvarGetLPSol(edgevars[MAX(i,j)][MIN(i,j)][0]);
         fractionality[i] += MIN(fval, 1 - fval) + MIN(1 - incval, incval);
      }
   }

   /* sort by fractionality of edgevars */
   SCIPsortDownRealInt(fractionality, idx, nbins);

   /* we try to construct partition inequalities from triangle-inequalities that are satisfied at equality*/
   for( i = 0; i < nbins && ncuts < MAXCUTS; ++i )
   {
      int candidate = idx[i];

      for( j = 0; j < nbins && ncuts < MAXCUTS; ++j )
      {
         for( k = 0; k < nbins && ncuts < MAXCUTS; ++k )
         {
            if( NULL == edgevars[candidate][j] || NULL == edgevars[candidate][k] || NULL == edgevars[j][k] || candidate == j || i ==k || j == k )
               continue;

            if( (candidate != j && candidate != k && j > k) )
            {
               nfirst = 1;
               firstpart[0] = candidate;
               nsecond = 2;
               secondpart[0] = j;
               secondpart[1] = k;

               triangleval = SCIPvarGetLPSol(edgevars[candidate][j][1]) + SCIPvarGetLPSol(edgevars[candidate][k][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) - 1;

               if( SCIPisLE(scip, triangleval, 0.0) )
               {
                  /* add a state to second partition*/
                  bestval = -SCIPinfinity(scip);
                  secondpart[2] = -1;
                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[candidate][l] || NULL == edgevars[MAX(l,j)][MIN(l,j)] || NULL == edgevars[MAX(l,k)][MIN(l,k)] )
                        continue;

                     addval = SCIPvarGetLPSol(edgevars[candidate][l][1]) - SCIPvarGetLPSol(edgevars[MAX(l,j)][MIN(l,j)][0]) - SCIPvarGetLPSol(edgevars[MAX(l,k)][MIN(l,k)][0]);

                     if( addval > bestval )
                     {
                        bestval = addval;
                        secondpart[2] = l;
                     }
                  }
                  nsecond++;

                  /*if partition with 3 in second and one in first violates inequality then add cut */
                  if( bestval > 0 )
                     SCIP_CALL( createPartitionCut(scip, sepa, edgevars, firstpart, secondpart, nfirst, nsecond, result, &ncuts) );

                  /* add one more state to first partition */
                  bestval = -SCIPinfinity(scip);
                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[MAX(candidate,l)][MIN(candidate,l)] || NULL == edgevars[l][secondpart[0]] || NULL == edgevars[l][secondpart[1]] || NULL == edgevars[l][secondpart[2]] )
                        continue;
                     addval = -SCIPvarGetLPSol(edgevars[MAX(candidate,l)][MIN(candidate,l)][0]) + SCIPvarGetLPSol(edgevars[l][secondpart[0]][1]) + SCIPvarGetLPSol(edgevars[l][secondpart[1]][1]) + SCIPvarGetLPSol(edgevars[l][secondpart[1]][1]);
                     if( addval > bestval)
                     {
                        bestval = addval;
                        firstpart[1] = l;
                     }
                  }

               }

               /* now try to find partition with 3 in first and 2 in second set */
               nfirst = 2;
               firstpart[0] = k;
               firstpart[1] = j;
               nsecond = 1;
               secondpart[0] = candidate;

               triangleval = SCIPvarGetLPSol(edgevars[j][candidate][1]) + SCIPvarGetLPSol(edgevars[k][candidate][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) - 1;

               if( SCIPisLE(scip, triangleval, 0.0) )
               {
                  /* add a state to first set*/
                  bestval = -SCIPinfinity(scip);
                  firstpart[2] = -1;

                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[l][candidate] || NULL == edgevars[MAX(l,j)][MIN(l,j)] || NULL == edgevars[MAX(l,k)][MIN(l,k)] )
                        continue;

                     addval = SCIPvarGetLPSol(edgevars[l][candidate][1]) - SCIPvarGetLPSol(edgevars[MAX(l,j)][MIN(l,j)][0]) - SCIPvarGetLPSol(edgevars[MAX(l,k)][MIN(l,k)][0]);
                     if( addval > bestval )
                     {
                        bestval = addval;
                        firstpart[2] = l;
                     }
                  }
                  nfirst++;

                  if( bestval > 0 )
                     SCIP_CALL( createPartitionCut(scip, sepa, edgevars, firstpart, secondpart, nfirst, nsecond, result, &ncuts) );

                  bestval = -SCIPinfinity(scip);

                  for( l = 0; l < nbins; ++l )
                  {
                     if( NULL == edgevars[MAX(l,candidate)][MIN(l,candidate)] || NULL == edgevars[l][firstpart[0]] || NULL == edgevars[l][firstpart[1]] || NULL == edgevars[l][firstpart[2]] )
                        continue;
                     addval = -SCIPvarGetLPSol(edgevars[MAX(l,candidate)][MIN(l,candidate)][0]) + SCIPvarGetLPSol(edgevars[firstpart[0]][l][1]) + SCIPvarGetLPSol(edgevars[firstpart[1]][l][1]) + SCIPvarGetLPSol(edgevars[firstpart[2]][l][1]);
                     if( addval > bestval)
                     {
                        bestval = addval;
                        secondpart[1] = l;
                     }
                  }

                  if( bestval > 0 )
                     SCIP_CALL( createPartitionCut(scip, sepa, edgevars, firstpart, secondpart, nfirst, nsecond, result, &ncuts) );
               }
            }
         }
      }
   }

   SCIPfreeMemoryArray(scip, &fractionality);
   SCIPfreeMemoryArray(scip, &idx);
   SCIPfreeMemoryArray(scip, &firstpart);
   SCIPfreeMemoryArray(scip, &secondpart);

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
