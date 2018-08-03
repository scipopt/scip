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
#define SEPA_PRIORITY              5000
#define SEPA_FREQ                     5
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                    2000 /**< maximal number of cuts that can be added to cut pool */
#define MAXCUTSCREATED            10000 /**< maximal number of cuts to select from */
#define MAXROUNDS                    20

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
   SCIP_VAR**** edgevars;                    /* edgevars from probdata */
   SCIP_DIGRAPH* edgegraph;                  /* edgegraph from probdata */
   char cutname[SCIP_MAXSTRLEN];             /* name of the cut */
   SCIP_Real** sign;                         /* matrix of sign-permuations */
   SCIP_Real* violation;                     /* array of violations */
   SCIP_ROW** cuts;                          /* array of generated cuts */
   SCIP_Real goodscorefac;                   /* parameters for cut-selection */
   SCIP_Real badscorefac;
   SCIP_Real goodmaxparall;
   SCIP_Real maxparall;
   SCIP_Real dircutoffdist;
   SCIP_Real efficacyweight;
   SCIP_Real objparalweight;
   SCIP_Real intsuppweight;
   int* succs1;                              /* successors of first state */
   int* succs2;                              /* successors of second state */
   int nstates;                              /* number of states */
   int ncutscreated;                         /* number of generated cuts */
   int ncutsapplied;                         /* number of cuts that were added to the pool */
   int size;                                 /* size of the cuts-array */
   int rounds;                               /* number of separation rounds */
   int states[3];                            /* states in a triangle */
   int nsuccs1;                              /* number of successors */
   int nsuccs2;
   int j;                                    /* running indices */
   int k;
   int l;
   SCIP_Bool usecutselection;                /* should cut selection be uses */

   rounds = SCIPsepaGetNCallsAtNode(sepa);
   if( rounds >= MAXROUNDS )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   edgegraph = SCIPcycGetEdgeGraph(scip);
   edgevars = SCIPcycGetEdgevars(scip);
   nstates = SCIPcycGetNBins(scip);

   assert(nstates > 0);
   assert(NULL != edgevars);
   assert(NULL != edgegraph);

   ncutscreated = 0;
   ncutsapplied = 0;
   size = MAXCUTS;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetBoolParam(scip, "cycleclustering/usecutselection", &usecutselection) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &violation, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cuts, size) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign, 3) );

   /* for some inequalities, you can swap the minus sign to any variable of the triangle */
   for( j = 0; j < 3; ++j )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &sign[j], 3) ); /*lint !e866*/

      for( k = 0; k < 3; ++k )
      {
         sign[j][k] = 1.0;
      }
      sign[j][j] = -1.0;
   }

   if( SCIPcycGetNCluster(scip) > 3 )
   {
      /* separate edges by the valid inequality z_ij1 + z_ik1 - y_jk0 <= 1 and z_ji + z_ki - y_jk1 <= 1
       * (minus sign can be anywhere)
       */
      for( states[0] = 0; states[0] < nstates && ncutscreated < MAXCUTSCREATED; ++states[0] )
      {
         succs1 = SCIPdigraphGetSuccessors(edgegraph, states[0]);
         nsuccs1 = SCIPdigraphGetNSuccessors(edgegraph, states[0]);

         for( j = 0; j < nsuccs1 && ncutscreated < MAXCUTSCREATED; ++j )
         {
            states[1] = succs1[j];
            succs2 = SCIPdigraphGetSuccessors(edgegraph, states[1]);
            nsuccs2 = SCIPdigraphGetNSuccessors(edgegraph, states[1]);

            for( k = 0; k < nsuccs2 && ncutscreated < MAXCUTSCREATED; ++k )
            {
               states[2] = succs2[k];

               if( !edgesExist(edgevars, states, 3) )
                  continue;

               if( (states[0] != states[1] && states[0] != states[2] && states[1] > states[2]) )
               {
                  /* permute the minus sign */
                  for( l = 0; l < 3 ; ++l )
                  {
                     violation[ncutscreated] = sign[l][0]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[1], 1));
                     violation[ncutscreated] += sign[l][1]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[2], 1));
                     violation[ncutscreated] += sign[l][2]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[2], 0)) - 1;

                     if( violation[ncutscreated] > 0 )
                     {
                        (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "trianglefw_%d_%d_%d_%d", states[0], states[1], states[2], l );
                        SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncutscreated]), sepa, cutname,
                           -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );

                        SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncutscreated]) );

                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[1], states[2], 0), sign[l][2]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[0], states[1], 1), sign[l][0]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[0], states[2], 1), sign[l][1]) );

                        SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncutscreated]) );

                        if( ncutscreated >= size - 1 )
                        {
                           SCIP_CALL( SCIPreallocBufferArray(scip, &violation, (int) (size + MAXCUTS)) );
                           SCIP_CALL( SCIPreallocBufferArray(scip, &cuts, (int) (size + MAXCUTS)) );
                           size += MAXCUTS;
                        }

                        ncutscreated++;
                     }

                     violation[ncutscreated] = sign[l][0]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[0], 1));
                     violation[ncutscreated] += sign[l][1]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[2], states[0], 1));
                     violation[ncutscreated] += sign[l][2]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[2], 0)) - 1;

                     if( violation[ncutscreated] > 0)
                     {
                        (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "trianglebw_%d_%d_%d_%d", states[0], states[1], states[2], l );
                        SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncutscreated]), sepa, cutname,
                           -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );

                        SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncutscreated]) );

                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[1], states[2], 0), sign[l][2]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[1], states[0], 1), sign[l][0]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[2], states[0], 1), sign[l][1]) );

                        SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncutscreated]) );

                        if( ncutscreated >= size - 1 )
                        {
                           SCIP_CALL( SCIPreallocBufferArray(scip, &violation, (int) (size + MAXCUTS)) );
                           SCIP_CALL( SCIPreallocBufferArray(scip, &cuts, (int) (size + MAXCUTS)) );
                           size += MAXCUTS;
                        }

                        ncutscreated++;
                     }
                  }
               }

               if( states[0] > states[1] && states[1] > states[2] )
               {
                  for( l = 0; l < 3; ++l )
                  {
                     violation[ncutscreated] = sign[l][0]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[1], 0));
                     violation[ncutscreated] += sign[l][1]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[2], 0));
                     violation[ncutscreated] += sign[l][2]
                        * SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[2], 0)) - 1;

                     if( violation[ncutscreated] > 0 )
                     {
                        (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", states[0], states[1], states[2]);
                        SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncutscreated]), sepa, cutname,
                           -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[0], states[1], 0), sign[l][0]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[0], states[2], 0), sign[l][1]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                           getEdgevar(edgevars, states[1], states[2], 0), sign[l][2]) );

                        if( ncutscreated >= size - 1 )
                        {
                           SCIP_CALL( SCIPreallocBufferArray(scip, &violation, (int) (size + MAXCUTS)) );
                           SCIP_CALL( SCIPreallocBufferArray(scip, &cuts, (int) (size + MAXCUTS)) );
                           size += MAXCUTS;
                        }

                        ncutscreated++;
                     }
                  }
               }
            }
         }
      }
   }

   if( SCIPcycGetNCluster(scip) == 3 )
   {
      for( states[0] = 0; states[0] < nstates; ++states[0] )
      {
         succs1 = SCIPdigraphGetSuccessors(edgegraph, states[0]);
         nsuccs1 = SCIPdigraphGetNSuccessors(edgegraph, states[0]);

         for( j = 0; j < nsuccs1 && ncutscreated < MAXCUTSCREATED; ++j )
         {
            states[1] = succs1[j];
            succs2 = SCIPdigraphGetSuccessors(edgegraph, states[1]);
            nsuccs2 = SCIPdigraphGetNSuccessors(edgegraph, states[1]);

            for( k = 0; k < nsuccs2 && ncutscreated < MAXCUTSCREATED; ++k )
            {
               states[2] = succs2[k];

               if( !edgesExist(edgevars, states, 3) )
                  continue;

               violation[ncutscreated] = SCIPvarGetLPSol(getEdgevar(edgevars, states[0], states[1], 1));
               violation[ncutscreated] += SCIPvarGetLPSol(getEdgevar(edgevars, states[1], states[2], 1));
               violation[ncutscreated] -= SCIPvarGetLPSol(getEdgevar(edgevars, states[2], states[0], 1)) - 1;

               if( violation[ncutscreated] > 0 )
               {
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", states[0], states[1], states[2]);
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncutscreated]), sepa, cutname,
                     -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );

                  SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncutscreated]) );

                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                     getEdgevar(edgevars, states[0], states[1], 1), 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                     getEdgevar(edgevars, states[1], states[2], 1), 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncutscreated],
                     getEdgevar(edgevars, states[2], states[0], 1), -1.0) );

                  SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncutscreated]) );

                  if( ncutscreated >= size - 1 )
                  {
                     SCIP_CALL( SCIPreallocBufferArray(scip, &violation, (int) (size + MAXCUTS)) );
                     SCIP_CALL( SCIPreallocBufferArray(scip, &cuts, (int) (size + MAXCUTS)) );
                     size += MAXCUTS;
                  }

                  ncutscreated++;
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
      SCIPsortDownRealPtr(violation, ((void **) cuts), ncutscreated);
      ncutsapplied = MIN(ncutscreated, MAXCUTS);
   }

   for( j = 0; j < ncutsapplied; ++j )
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cuts[j]) );
      *result = SCIP_SEPARATED;
   }

   /* free memory */
   for( j = 0; j < ncutscreated; ++j )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(cuts[j])) );
   }
   SCIPfreeBufferArray(scip, &cuts);
   SCIPfreeBufferArray(scip, &violation);
   for( j = 0; j < 3; ++j )
   {
      SCIPfreeMemoryArray(scip, &sign[j]);
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
