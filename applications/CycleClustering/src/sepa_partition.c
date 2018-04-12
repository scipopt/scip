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

#define SEPA_NAME              "partition"
#define SEPA_DESC              "separator to separate partition-inequalities in cycle-clustering application"
#define SEPA_PRIORITY              1500
#define SEPA_FREQ                     5
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                    5000
#define MAXROUNDS                    20


/** Given two partitions S, T creates the corresponding cut and adds it do SCIP */
static
SCIP_RETCODE createPartitionCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< Separator */
   SCIP_ROW***           cuts,               /**< Array to store generated cut */
   int*                  cutsize,            /**< Size of the cut array */
   int*                  ncutscreated,       /**< Number of created cuts */
   int*                  firstpart,          /**< The first partition */
   int*                  secondpart,         /**< The second partition */
   int                   nfirst,             /**< Number of states in first partition */
   int                   nsecond,            /**< Number of states in second partition */
   SCIP_Real**           violations,         /**< Array to stor the violation of each cut */
   SCIP_Real             violation           /**< Violation of the cut that should be created */
   )
{
   SCIP_VAR**** edgevars;
   char cutname[SCIP_MAXSTRLEN];
   int i;
   int j;
   int inda;
   int indb;

   edgevars = SCIPcycGetEdgevars(scip);

   assert(NULL != edgevars);

   if( *cutsize - 1 <= *ncutscreated )
   {
      *cutsize = *cutsize * 2;
      SCIP_CALL( SCIPreallocMemoryArray(scip, cuts, (int) *cutsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, violations, (int) *cutsize) );
   }

   (*violations)[*ncutscreated] = violation;

   /* create cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "PartitionCut_%d_%d", nfirst, nsecond);
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &((*cuts)[*ncutscreated]), sepa, cutname, -SCIPinfinity(scip), (SCIP_Real) MIN(nfirst, nsecond), FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, (*cuts)[*ncutscreated]) );

   for( i = 0; i < nfirst; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         inda = MAX(firstpart[i], firstpart[j]);
         indb = MIN(firstpart[i], firstpart[j]);
         SCIP_CALL( SCIPaddVarToRow(scip, (*cuts)[*ncutscreated], edgevars[inda][indb][0], -1.0) );
      }
   }
   for( i = 0; i < nsecond; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         inda = MAX(secondpart[i], secondpart[j]);
         indb = MIN(secondpart[i], secondpart[j]);
         SCIP_CALL( SCIPaddVarToRow(scip, (*cuts)[*ncutscreated], edgevars[inda][indb][0], -1.0) );
      }
   }
   for( i = 0; i < nfirst; ++i )
   {
      for( j = 0; j < nsecond; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, (*cuts)[*ncutscreated], edgevars[firstpart[i]][secondpart[j]][1], 1.0) );
      }
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, (*cuts)[*ncutscreated]) );

   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, (*cuts)[*ncutscreated], NULL) ) );
   (*ncutscreated)++;

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
   SCIP_DIGRAPH* edgegraph;
   int* idx;
   int states[5];
   SCIP_Real violation;
   SCIP_Real violationchg;
   SCIP_Real bestvalue;
   SCIP_Real lpvalforward;
   SCIP_Real lpvalincluster;
   SCIP_Real* violations;
   SCIP_ROW** cuts;
   int cutsize;
   int ncutscreated;
   int ncutsapplied;
   int maxcutscreated;
   int* firstpart;
   int* secondpart;
   int nfirst;
   int nsecond;
   int nstates;
   int rounds;
   int i;
   int j;
   int k;
   int l;
   SCIP_Bool usecutselection;


   /* get necessary probdata */
   edgevars = SCIPcycGetEdgevars(scip);
   edgegraph = SCIPcycGetEdgeGraph(scip);
   nstates = SCIPcycGetNBins(scip);
   rounds = SCIPsepaGetNCallsAtNode(sepa);
   cutsize = MAXCUTS;
   ncutscreated = 0;
   maxcutscreated = 10 * MAXCUTS;
   SCIP_CALL( SCIPgetBoolParam(scip, "cycleclustering/usecutselection", &usecutselection) );


   assert(nstates > 0);
   assert(NULL != edgevars);
   assert(NULL != edgegraph);

   *result = SCIP_DIDNOTFIND;

   if( rounds >= MAXROUNDS )
   {
      *result =  SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &fractionality, nstates) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &idx, nstates) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &firstpart, nstates) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &secondpart, nstates) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cuts, cutsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &violations, cutsize) );

   /* sort edges by decreasing fractionality of lp-solution */
   for( i = 0; i < nstates; ++i )
   {
      idx[i] = i;
      fractionality[i] = 0;

      for( j = 0; j < SCIPdigraphGetNSuccessors(edgegraph, i); ++j )
      {
         states[0] = i;
         states[1] = SCIPdigraphGetSuccessors(edgegraph, i)[j];

         assert(NULL != edgevars[states[0]][states[1]]);

         lpvalforward = SCIPvarGetLPSol(edgevars[states[0]][states[1]][1]);
         lpvalincluster = SCIPvarGetLPSol(edgevars[MAX(states[0],states[1])][MIN(states[0],states[1])][0]);
         fractionality[states[0]] += MIN(lpvalforward, 1 - lpvalforward) + MIN(1 - lpvalincluster, lpvalincluster);
      }
   }

   /* sort by fractionality of edgevars */
   SCIPsortDownRealInt(fractionality, idx, nstates);

   /* we try to construct partition inequalities from triangle-inequalities that are almost satisfied at equality */
   for( i = 0; i < nstates && ncutscreated < maxcutscreated; ++i )
   {
      states[0] = idx[i];

      for( j = 0; j < SCIPdigraphGetNSuccessors(edgegraph, states[0]) && ncutscreated < maxcutscreated; ++j )
      {
         states[1] = SCIPdigraphGetSuccessors(edgegraph, states[0])[j];

         for( k = 0; k < SCIPdigraphGetNSuccessors(edgegraph, states[1]) && ncutscreated < maxcutscreated; ++k )
         {
            states[2] = SCIPdigraphGetSuccessors(edgegraph, states[1])[k];

            if( NULL == edgevars[states[0]][states[1]] || NULL == edgevars[states[0]][states[2]] || NULL == edgevars[states[1]][states[2]] )
               continue;

            if( states[1] > states[2] )
            {
               /* first case, construct partition with 2 predecessors and 3 successors */
               nfirst = 1;
               firstpart[0] = states[0];
               firstpart[1] = -1;
               nsecond = 2;
               secondpart[0] = states[1];
               secondpart[1] = states[2];
               secondpart[2] = -1;

               /* Get violation of trianlge inequality for these three states */
               violation = SCIPvarGetLPSol(edgevars[states[0]][states[1]][1]);
               violation += SCIPvarGetLPSol(edgevars[states[0]][states[2]][1]);
               violation -= SCIPvarGetLPSol(edgevars[states[1]][states[2]][0]);
               violation -= 1;

               if( SCIPisGE(scip, violation, -0.2) )
               {
                  /* add a state to second partition*/
                  bestvalue = -SCIPinfinity(scip);
                  secondpart[2] = -1;
                  for( l = 0; l < SCIPdigraphGetNSuccessors(edgegraph, states[2]); ++l )
                  {
                     states[3] = SCIPdigraphGetSuccessors(edgegraph, states[2])[l];
                     if( NULL == edgevars[states[0]][states[3]] || NULL == edgevars[states[1]][states[3]] || NULL == edgevars[states[2]][states[3]] )
                        continue;

                     violationchg = SCIPvarGetLPSol(edgevars[states[0]][states[3]][1]);
                     violationchg -= SCIPvarGetLPSol(edgevars[MAX(states[1],states[3])][MIN(states[1],states[3])][0]);
                     violationchg -= SCIPvarGetLPSol(edgevars[MAX(states[2],states[3])][MIN(states[2],states[3])][0]);

                     if( violationchg > bestvalue )
                     {
                        bestvalue = violationchg;
                        secondpart[2] = states[3];
                     }
                  }

                  states[3] = secondpart[2];

                  /* if we did not find a state that we can add we can stop */
                  if( states[3] == -1 )
                     continue;
                  nsecond++;

                  violation += bestvalue;

                  /* add one more state to first partition */
                  bestvalue = -SCIPinfinity(scip);
                  for( l = 0; l < SCIPdigraphGetNSuccessors(edgegraph, states[3]); ++l )
                  {
                     states[4] = SCIPdigraphGetSuccessors(edgegraph, states[3])[l];

                     if( NULL == edgevars[states[0]][states[4]] || NULL == edgevars[states[1]][states[4]] || NULL == edgevars[states[2]][states[4]] || NULL == edgevars[states[3]][states[4]] )
                        continue;

                     /* compute what has changed from the violation of the 1-4 inequality */
                     violationchg = -SCIPvarGetLPSol(edgevars[MAX(states[0], states[4])][MIN(states[0],states[4])][0]) - 1.0;
                     violationchg += SCIPvarGetLPSol(edgevars[states[4]][secondpart[0]][1]);
                     violationchg += SCIPvarGetLPSol(edgevars[states[4]][secondpart[1]][1]);
                     violationchg += SCIPvarGetLPSol(edgevars[states[4]][secondpart[2]][1]);

                     /* create cut if inequality is violated by lp-solution */
                     if( SCIPisPositive(scip, violation + violationchg) )
                     {
                        firstpart[1] = states[4];
                        nfirst = 2;
                        SCIP_CALL( createPartitionCut(scip, sepa, &cuts, &cutsize, &ncutscreated, firstpart, secondpart, nfirst, nsecond, &violations, violation + violationchg) );
                        break;
                     }
                  }
               }

               /* now try to find partition with 3 in first and 2 in second set */
               nfirst = 2;
               firstpart[0] = states[1];
               firstpart[1] = states[2];
               firstpart[2] = -1;
               nsecond = 1;
               secondpart[0] = states[0];
               secondpart[1] = -1;

               violation = SCIPvarGetLPSol(edgevars[states[1]][states[0]][1]);
               violation += SCIPvarGetLPSol(edgevars[states[2]][states[0]][1]);
               violation -= SCIPvarGetLPSol(edgevars[states[1]][states[2]][0]);
               violation -= 1;

               if( SCIPisGE(scip, violation, -0.2) )
               {
                  /* add a state to second partition*/
                  bestvalue = -SCIPinfinity(scip);
                  firstpart[2] = -1;
                  for( l = 0; l < SCIPdigraphGetNSuccessors(edgegraph, states[2]); ++l )
                  {
                     states[3] = SCIPdigraphGetSuccessors(edgegraph, states[2])[l];
                     if( NULL == edgevars[states[0]][states[3]] || NULL == edgevars[states[1]][states[3]] || NULL == edgevars[states[2]][states[3]] )
                        continue;

                     violationchg = SCIPvarGetLPSol(edgevars[states[3]][states[0]][1]);
                     violationchg -= SCIPvarGetLPSol(edgevars[MAX(states[1],states[3])][MIN(states[1],states[3])][0]);
                     violationchg -= SCIPvarGetLPSol(edgevars[MAX(states[2],states[3])][MIN(states[2],states[3])][0]);

                     if( violationchg > bestvalue )
                     {
                        bestvalue = violationchg;
                        firstpart[2] = states[3];
                     }
                  }

                  states[3] = firstpart[2];

                  if( states[3] == -1 )
                     continue;
                  nfirst++;

                  violation += bestvalue;

                  /* add one more state to second partition */
                  bestvalue = -SCIPinfinity(scip);
                  for( l = 0; l < SCIPdigraphGetNSuccessors(edgegraph, states[3]); ++l )
                  {
                     states[4] = SCIPdigraphGetSuccessors(edgegraph, states[3])[l];

                     if( NULL == edgevars[states[0]][states[4]] || NULL == edgevars[states[1]][states[4]] || NULL == edgevars[states[2]][states[4]] || NULL == edgevars[states[3]][states[4]] )
                        continue;

                     violationchg = -SCIPvarGetLPSol(edgevars[MAX(states[0], states[4])][MIN(states[0],states[4])][0]) - 1.0;
                     violationchg += SCIPvarGetLPSol(edgevars[firstpart[0]][states[4]][1]);
                     violationchg += SCIPvarGetLPSol(edgevars[firstpart[1]][states[4]][1]);
                     violationchg += SCIPvarGetLPSol(edgevars[firstpart[2]][states[4]][1]);

                     if( SCIPisPositive(scip, violation + violationchg) )
                     {
                        secondpart[1] = states[4];
                        nsecond = 2;
                        SCIP_CALL( createPartitionCut(scip, sepa, &cuts, &cutsize, &ncutscreated, firstpart, secondpart, nfirst, nsecond, &violations, violation + violationchg) );
                        break;
                     }
                  }
               }
            }
         }
      }
   }

   /* apply the cuts with the highest violation or use cut-selection */
   if( usecutselection )
      SCIP_CALL( SCIPselectCuts(scip, cuts, NULL, 0.8, 0.0, 0.1, 0.5, 0.5, 1.0, 0.1, 0.1, ncutscreated, 0, MAXCUTS, &ncutsapplied ) );
   else
   {
      SCIPsortDownRealPtr(violations, ((void **) cuts), ncutscreated);
      ncutsapplied = MIN(ncutscreated, MAXCUTS);
   }
   for( j = 0; j < ncutsapplied; ++j )
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cuts[j]) );
      *result = SCIP_SEPARATED;
   }

   SCIPfreeMemoryArray(scip, &fractionality);
   SCIPfreeMemoryArray(scip, &idx);
   SCIPfreeMemoryArray(scip, &firstpart);
   SCIPfreeMemoryArray(scip, &secondpart);
   for( i = 0; i < ncutscreated; ++i )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(cuts[i])) );
   }
   SCIPfreeMemoryArray(scip, &cuts);
   SCIPfreeMemoryArray(scip, &violations);

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
