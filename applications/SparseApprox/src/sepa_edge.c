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

/**@file   sepa_edge.c
 * @brief  edge-separator
 *  @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_edge.h"
#include "probdata_spa.h"
#include "scip/cons_linear.h"

#define SEPA_NAME              "edge"
#define SEPA_DESC              "separator for graph partitioning in sparse approximation project"
#define SEPA_PRIORITY               5000
#define SEPA_FREQ                   1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define MAXROUNDS                 5
#define MAXCUTS                   500

/** copy method for separator plugins (called when SCIP copies plugins) */

static
SCIP_DECL_SEPACOPY(sepaCopyEdge)
{   /*lint --e{715}*/
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
{
   SCIP_VAR**** edgevars;

   int nbins;
   int ncluster;
   int ncuts;
   int i;
   int j;
   int k;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_ROW* cut;
   SCIP_Bool infeasible;
   int rounds;

   rounds = SCIPsepaGetNCallsAtNode(sepa);
   if( rounds >= MAXROUNDS && SCIPgetDepth(scip) == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);

   assert(nbins > 0 && ncluster > 0);
   assert(NULL != edgevars);

   ncuts = 0;

   *result = SCIP_DIDNOTFIND;

   /* separate edges by the valied inequality y_ij1 + y_ik1 - y_jk0 <= 1 */
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < nbins; ++j )
      {
         for( k = 0; k < j; ++k )
         {
            if( (i > j) || (i < k) )
            {
               if( SCIPvarGetLPSol(edgevars[i][j][1]) + SCIPvarGetLPSol(edgevars[i][k][1]) - SCIPvarGetLPSol(edgevars[j][k][0]) > 1 )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][1], 1) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][k][1], 1) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], -1) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
                     if( !infeasible )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                        *result = SCIP_SEPARATED;
                     }
                     else
                        *result = SCIP_CUTOFF;
                     ncuts++;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
            }
         }
      }
   }

   /* separate edges by the valied inequality y_ij1 + y_jk1 - y_ki0 <= 1 */
   if ( ncluster == 3 )
   {
      for( i = 0; i < nbins; ++i )
      {
         for( j = 0; j < i; ++j )
         {
            for( k = 0; k < j; ++k )
            {
               if( SCIPvarGetLPSol(edgevars[i][j][1]) + SCIPvarGetLPSol(edgevars[j][k][1]) - SCIPvarGetLPSol(edgevars[k][i][1]) > 1 )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][1], 1) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][1], 1) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[k][i][1], -1) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
                     if( !infeasible )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                        *result = SCIP_SEPARATED;
                     }
                     else
                        *result = SCIP_CUTOFF;
                     ncuts++;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
               if( SCIPvarGetLPSol(edgevars[k][j][1]) + SCIPvarGetLPSol(edgevars[j][i][1]) - SCIPvarGetLPSol(edgevars[i][k][1]) > 1 )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[k][j][1], 1) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][i][1], 1) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][k][1], -1) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
                     if( !infeasible )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                        *result = SCIP_SEPARATED;
                     }
                     else
                        *result = SCIP_CUTOFF;
                     ncuts++;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
            }
         }
      }
   }
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
