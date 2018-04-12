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
#define SEPA_PRIORITY              5000
#define SEPA_FREQ                     5
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXCUTS                    2000
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
   int nstates;                              /* number of states */
   int ncuts;                                /* number of generated cuts */
   int ncutsapplied;                         /* number of cuts that were added to the pool */
   int size;                                 /* size of the cuts-array */
   int rounds;                               /* number of separation rounds */
   int state1;                               /* states in triangle */
   int state2;
   int state3;
   int j;                                    /* running indices */
   int k;
   int l;
   SCIP_Bool usecutselection;

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

   ncuts = 0;
   ncutsapplied = 0;
   size = MAXCUTS;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetBoolParam(scip, "cycleclustering/usecutselection", &usecutselection) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &violation, size) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cuts, size) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign, 3) );

   /** for some inequalities, you can swap the minus sign to any variable of the triangle */
   for( j = 0; j < 3; ++j )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign[j], 3) ); /*lint !e866*/
      for( k = 0; k < 3; ++k )
      {
         sign[j][k] = 1.0;
      }
      sign[j][j] = -1.0;
   }

   if( SCIPcycGetNCluster(scip) > 3 )
   {
      /* separate edges by the valid inequality z_ij1 + z_ik1 - y_jk0 <= 1 and z_ji + z_ki - y_jk1 <= 1 (minus sign can be anywhere) */
      for( state1 = 0; state1 < nstates; ++state1 )
      {
         for( j = 0; j < SCIPdigraphGetNSuccessors(edgegraph, state1) ; ++j )
         {
            state2 = SCIPdigraphGetSuccessors(edgegraph, state1)[j];

            for( k = 0; k < SCIPdigraphGetNSuccessors(edgegraph, state2); ++k )
            {
               state3 = SCIPdigraphGetSuccessors(edgegraph, state2)[k];

               if( NULL == edgevars[state1][state2] || NULL == edgevars[state1][state3] || NULL == edgevars[state2][state3] )
                  continue;
               if( (state1 != state2 && state1 != state3 && state2 > state3) )
               {
                  /* permute the minus sign */
                  for( l = 0; l < 3 ; ++l )
                  {
                     violation[ncuts] = sign[l][0] * SCIPvarGetLPSol(edgevars[state1][state2][1]) + sign[l][1] * SCIPvarGetLPSol(edgevars[state1][state3][1]) + sign[l][2] * SCIPvarGetLPSol(edgevars[state2][state3][0]) - 1;

                     if( violation[ncuts] > 0)
                     {
                        (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "trianglefw_%d_%d_%d_%d", state1, state2, state3, l );
                        SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                        SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state2][state3][0], sign[l][2]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state1][state2][1], sign[l][0]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state1][state3][1], sign[l][1]) );
                        SCIP_CALL( SCIPflushRowExtensions(scip, cuts[ncuts]) );

                        if( ncuts >= size - 1 )
                        {
                           SCIP_CALL( SCIPreallocMemoryArray(scip, &violation, (int) (size + MAXCUTS)) );
                           SCIP_CALL( SCIPreallocMemoryArray(scip, &cuts, (int) (size + MAXCUTS)) );
                           size += MAXCUTS;
                        }

                        ncuts++;
                     }

                     violation[ncuts] = sign[l][0] * SCIPvarGetLPSol(edgevars[state2][state1][1]) + sign[l][1] * SCIPvarGetLPSol(edgevars[state3][state1][1]) + sign[l][2] * SCIPvarGetLPSol(edgevars[state2][state3][0]) - 1;

                     if( violation[ncuts] > 0)
                     {
                        (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "trianglebw_%d_%d_%d_%d", state1, state2, state3, l );
                        SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                        SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state2][state3][0], sign[l][2]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state2][state1][1], sign[l][0]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state3][state1][1], sign[l][1]) );
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

               if( state1 > state2 && state2 > state3 )
               {
                  for( l = 0; l < 3; ++l )
                  {
                     violation[ncuts] = sign[l][0] * SCIPvarGetLPSol(edgevars[state1][state2][0]) + sign[l][1] * SCIPvarGetLPSol(edgevars[state1][state3][0]) + sign[l][2] * SCIPvarGetLPSol(edgevars[state2][state3][0]) - 1;

                     if( violation[ncuts] > 0 )
                     {
                        (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d_%f%f%f", state1, state2, state3, sign[l][0], sign[l][1], sign[l][2] );
                        SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state1][state2][0], sign[l][0]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state1][state3][0], sign[l][1]) );
                        SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state2][state3][0], sign[l][2]) );

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
   }

   if( SCIPcycGetNCluster(scip) == 3 )
   {
     for( state1 = 0; state1 < nstates; ++state1 )
      {
         for( j = 0; j < SCIPdigraphGetNSuccessors(edgegraph, state1) ; ++j )
         {
            state2 = SCIPdigraphGetSuccessors(edgegraph, state1)[j];

            for( k = 0; k < SCIPdigraphGetNSuccessors(edgegraph, state2); ++k )
            {
               state3 = SCIPdigraphGetSuccessors(edgegraph, state2)[k];

               if (NULL == edgevars[state1][state2] || NULL == edgevars[state2][state3] || NULL == edgevars[state1][state3])
                  continue;

               violation[ncuts] = SCIPvarGetLPSol(edgevars[state1][state2][1]) + SCIPvarGetLPSol(edgevars[state2][state3][1]) - SCIPvarGetLPSol(edgevars[state3][state1][1]) - 1;

               if( violation[ncuts] > 0 )
               {
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", state1, state2, state3);
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &(cuts[ncuts]), sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cuts[ncuts]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state1][state2][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state2][state3][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cuts[ncuts], edgevars[state3][state1][1], -1.0) );
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

   /* apply the cuts with the highest violation or use cut-selection */
   if( usecutselection )
      SCIP_CALL( SCIPselectCuts(scip, cuts, NULL, 0.8, 0.0, 0.1, 0.5, 0.5, 1.0, 0.1, 0.1, ncuts, 0, MAXCUTS, &ncutsapplied ) );
   else
   {
      SCIPsortDownRealPtr(violation, ((void **) cuts), ncuts);
      ncutsapplied = MIN(ncuts, MAXCUTS);
   }
   for( j = 0; j < ncutsapplied; ++j )
   {
      SCIP_CALL( SCIPaddPoolCut(scip, cuts[j]) );
      *result = SCIP_SEPARATED;
   }

   /* free memory */
   for( j = 0; j < ncuts; ++j )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(cuts[j])) );
   }
   SCIPfreeMemoryArray(scip, &cuts);
   SCIPfreeMemoryArray(scip, &violation);
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
