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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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
#define MAXCUTS                    2000 /**< maximal number of cuts that can be added to cut pool */
#define MAXCUTSCREATED            10000 /**< maximal number of cuts to select from */
#define MAXROUNDS                    20 /**< maximal number of separation rounds per node */
#define MAXTRIANGLEDISTANCE        -0.2 /**< maximal negative violation of triangle-inequality to construct cut from */


/** Given two partitions S, T creates the corresponding cut and adds it do SCIP */
static
SCIP_RETCODE createPartitionCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_ROW***           cuts,               /**< array to store generated cut */
   int*                  cutsize,            /**< size of the cut array */
   int*                  ncutscreated,       /**< number of created cuts */
   int*                  firstpart,          /**< the first partition */
   int*                  secondpart,         /**< the second partition */
   int                   nfirst,             /**< number of states in first partition */
   int                   nsecond,            /**< number of states in second partition */
   SCIP_Real**           violations,         /**< array to stor the violation of each cut */
   SCIP_Real             violation           /**< violation of the cut that should be created */
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
      SCIP_CALL( SCIPreallocBufferArray(scip, cuts, (int) *cutsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, violations, (int) *cutsize) );
   }

   (*violations)[*ncutscreated] = violation;

   /* create cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "PartitionCut_%d_%d", nfirst, nsecond);
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &((*cuts)[*ncutscreated]), sepa, cutname, -SCIPinfinity(scip),
      (SCIP_Real) MIN(nfirst, nsecond), FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, (*cuts)[*ncutscreated]) );

   for( i = 0; i < nfirst; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         inda = MAX(firstpart[i], firstpart[j]);
         indb = MIN(firstpart[i], firstpart[j]);
         SCIP_CALL( SCIPaddVarToRow(scip, (*cuts)[*ncutscreated], getEdgevar(edgevars, inda, indb, 0), -1.0) );
      }
   }

   for( i = 0; i < nsecond; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         inda = MAX(secondpart[i], secondpart[j]);
         indb = MIN(secondpart[i], secondpart[j]);
         SCIP_CALL( SCIPaddVarToRow(scip, (*cuts)[*ncutscreated], getEdgevar(edgevars, inda, indb, 0), -1.0) );
      }
   }

   for( i = 0; i < nfirst; ++i )
   {
      for( j = 0; j < nsecond; ++j )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, (*cuts)[*ncutscreated],
            getEdgevar(edgevars, firstpart[i], secondpart[j], 1), 1.0) );
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
   SCIP_Real goodscorefac;
   SCIP_Real badscorefac;
   SCIP_Real goodmaxparall;
   SCIP_Real maxparall;
   SCIP_Real dircutoffdist;
   SCIP_Real efficacyweight;
   SCIP_Real objparalweight;
   SCIP_Real intsuppweight;
   SCIP_Real* violations;
   SCIP_ROW** cuts;
   int cutsize;
   int ncutscreated;
   int ncutsapplied;
   int* firstpart;
   int* secondpart;
   int** successors;
   int* nsuccessors;
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

   SCIP_CALL( SCIPgetBoolParam(scip, "cycleclustering/usecutselection", &usecutselection) );

   assert(nstates > 0);
   assert(NULL != edgevars);
   assert(NULL != edgegraph);

   *result = SCIP_DIDNOTFIND;

   if( SCIPcycGetNCluster(scip) == 3 || rounds >= MAXROUNDS )
   {
      *result =  SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &successors, 5) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsuccessors, 5) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &fractionality, nstates) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idx, nstates) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &firstpart, nstates) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &secondpart, nstates) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cuts, cutsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violations, cutsize) );

   /* sort edges by decreasing fractionality of lp-solution */
   for( i = 0; i < nstates; ++i )
   {
      idx[i] = i;
      fractionality[i] = 0;
      successors[0] = SCIPdigraphGetSuccessors(edgegraph, i);
      nsuccessors[0] = SCIPdigraphGetNSuccessors(edgegraph, i);

      for( j = 0; j < nsuccessors[0]; ++j )
      {
         states[0] = i;
         states[1] = successors[0][j];

         lpvalforward = SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[1], 1));
         lpvalincluster = SCIPvarGetLPSol(getEdgevar(edgevars, MAX(states[0],states[1]), MIN(states[0],states[1]), 0));
         fractionality[states[0]] += MIN(lpvalforward, 1 - lpvalforward) + MIN(1 - lpvalincluster, lpvalincluster);
      }
   }

   /* sort by fractionality of edgevars */
   SCIPsortDownRealInt(fractionality, idx, nstates);

   /* we try to construct partition inequalities from triangle-inequalities that are almost satisfied at equality */
   for( i = 0; i < nstates && ncutscreated < MAXCUTSCREATED; ++i )
   {
      states[0] = idx[i];
      successors[0] = SCIPdigraphGetSuccessors(edgegraph, states[0]);
      nsuccessors[0] = SCIPdigraphGetNSuccessors(edgegraph, states[0]);

      for( j = 0; j < nsuccessors[0] && ncutscreated < MAXCUTSCREATED; ++j )
      {
         states[1] = successors[0][j];
         successors[1] = SCIPdigraphGetSuccessors(edgegraph, states[1]);
         nsuccessors[1] = SCIPdigraphGetNSuccessors(edgegraph, states[1]);

         for( k = 0; k < nsuccessors[1] && ncutscreated < MAXCUTSCREATED; ++k )
         {
            states[2] = successors[1][k];
            successors[2] = SCIPdigraphGetSuccessors(edgegraph, states[2]);
            nsuccessors[2] = SCIPdigraphGetNSuccessors(edgegraph, states[2]);

            /* check if all edges in triangle exist */
            if( !edgesExist(edgevars, states, 3) )
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

               /* get violation of trianlge inequality for these three states */
               violation = SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[1], 1));
               violation += SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[2], 1));
               violation -= SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[2], 0));
               violation -= 1;

               if( SCIPisGE(scip, violation, MAXTRIANGLEDISTANCE) )
               {
                  /* add a state to second partition*/
                  bestvalue = -SCIPinfinity(scip);
                  secondpart[2] = -1;
                  for( l = 0; l < nsuccessors[2]; ++l )
                  {
                     states[3] = successors[2][l];
                     if( !edgesExist(edgevars, states, 4) )
                        continue;

                     violationchg = SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[3], 1));
                     violationchg -= SCIPvarGetLPSol(getEdgevar(edgevars,
                        MAX(states[1],states[3]), MIN(states[1],states[3]), 0));
                     violationchg -= SCIPvarGetLPSol(getEdgevar(edgevars,
                        MAX(states[2],states[3]), MIN(states[2],states[3]), 0));

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

                  successors[3] = SCIPdigraphGetSuccessors(edgegraph, states[3]);
                  nsuccessors[3] = SCIPdigraphGetNSuccessors(edgegraph, states[3]);

                  nsecond++;
                  violation += bestvalue;

                  /* add one more state to first partition */
                  bestvalue = -SCIPinfinity(scip);
                  for( l = 0; l < nsuccessors[3]; ++l )
                  {
                     states[4] = successors[3][l];

                     if( !edgesExist(edgevars, states, 5) )
                        continue;

                     /* compute what has changed from the violation of the 1-4 inequality */
                     violationchg = -SCIPvarGetLPSol(getEdgevar(edgevars,
                        MAX(states[0], states[4]), MIN(states[0],states[4]), 0)) - 1.0;
                     violationchg += SCIPvarGetLPSol(getEdgevar(edgevars, states[4], secondpart[0], 1));
                     violationchg += SCIPvarGetLPSol(getEdgevar(edgevars, states[4], secondpart[1], 1));
                     violationchg += SCIPvarGetLPSol(getEdgevar(edgevars, states[4], secondpart[2], 1));

                     /* create cut if inequality is violated by lp-solution */
                     if( SCIPisPositive(scip, violation + violationchg) )
                     {
                        firstpart[1] = states[4];
                        nfirst = 2;
                        SCIP_CALL( createPartitionCut(scip, sepa, &cuts, &cutsize, &ncutscreated, firstpart, secondpart,
                           nfirst, nsecond, &violations, violation + violationchg) );

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

               violation = SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[0], 1));
               violation += SCIPvarGetLPSol(getEdgevar(edgevars, states[2], states[0], 1));
               violation -= SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[2], 0));
               violation -= 1;

               if( SCIPisGE(scip, violation, MAXTRIANGLEDISTANCE) )
               {
                  /* add a state to second partition*/
                  bestvalue = -SCIPinfinity(scip);
                  firstpart[2] = -1;
                  for( l = 0; l < nsuccessors[2]; ++l )
                  {
                     states[3] = successors[2][l];
                     if( !edgesExist(edgevars, states, 4) )
                        continue;

                     violationchg = SCIPvarGetLPSol(getEdgevar(edgevars, states[3], states[0], 1));
                     violationchg -= SCIPvarGetLPSol(getEdgevar(edgevars,
                        MAX(states[1],states[3]), MIN(states[1],states[3]), 0));
                     violationchg -= SCIPvarGetLPSol(getEdgevar(edgevars,
                        MAX(states[2],states[3]), MIN(states[2],states[3]), 0));

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
                  for( l = 0; l < nsuccessors[3]; ++l )
                  {
                     states[4] = successors[3][l];

                     if( !edgesExist(edgevars, states, 5) )
                        continue;

                     violationchg = -SCIPvarGetLPSol(getEdgevar(edgevars,
                        MAX(states[0], states[4]), MIN(states[0],states[4]), 0)) - 1.0;
                     violationchg += SCIPvarGetLPSol(getEdgevar(edgevars, firstpart[0], states[4], 1));
                     violationchg += SCIPvarGetLPSol(getEdgevar(edgevars, firstpart[1], states[4], 1));
                     violationchg += SCIPvarGetLPSol(getEdgevar(edgevars, firstpart[2], states[4], 1));
                     violationchg += SCIPvarGetLPSol(edgevars[firstpart[2]][states[4]][1]);

                     if( SCIPisPositive(scip, violation + violationchg) )
                     {
                        secondpart[1] = states[4];
                        nsecond = 2;
                        SCIP_CALL( createPartitionCut(scip, sepa, &cuts, &cutsize, &ncutscreated, firstpart,
                           secondpart, nfirst, nsecond, &violations, violation + violationchg) );

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
   {
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/goodscorefac", &goodscorefac) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/badscorefac", &badscorefac) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/goodmaxparall", &goodmaxparall) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/maxparall", &maxparall) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/dircutoffdist", &dircutoffdist) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/efficacyweight", &efficacyweight) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/objparalweight", &objparalweight) );
      SCIP_CALL( SCIPgetRealParam(scip, "cycleclustering/intsuppweight", &intsuppweight) );

      SCIP_CALL( SCIPselectCuts(scip, cuts, NULL, goodscorefac, badscorefac, goodmaxparall, maxparall, dircutoffdist,
         efficacyweight, objparalweight, intsuppweight, ncutscreated, 0, MAXCUTS, &ncutsapplied ) );
   }
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

   SCIPfreeBlockMemoryArray(scip, &fractionality, nstates);
   SCIPfreeBlockMemoryArray(scip, &idx, nstates);
   SCIPfreeBlockMemoryArray(scip, &firstpart, nstates);
   SCIPfreeBlockMemoryArray(scip, &secondpart, nstates);

   for( i = 0; i < ncutscreated; ++i )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(cuts[i])) );
   }

   SCIPfreeBufferArray(scip, &cuts);
   SCIPfreeBufferArray(scip, &violations);
   SCIPfreeBufferArray(scip, &nsuccessors);
   SCIPfreeBufferArray(scip, &successors);

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
      SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpPartition, NULL, NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyPartition) );

   return SCIP_OKAY;
}
