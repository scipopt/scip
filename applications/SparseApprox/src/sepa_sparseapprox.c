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

/**@file   sepa_SparseApprox.c
 * @brief  SparseApprox separator
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_sparseapprox.h"
#include "probdata_spa.h"
#include "scip/cons_linear.h"

#define SEPA_NAME              "sparsapprox"
#define SEPA_DESC              "separator for graph partitioning in sparse approximation project"
#define SEPA_PRIORITY               5000
#define SEPA_FREQ                   10
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


/** copy method for separator plugins (called when SCIP copies plugins) */

static
SCIP_DECL_SEPACOPY(sepaCopySparseApprox)
{   /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaSparseApprox(scip) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpSparseApprox)
{
   SCIP_VAR***** edgevars;
   SCIP_VAR*** binvars;
   int i;               /* the bin in the first partition */
   int j;
   int w;
   int k;
   int nbins;
   int ncluster;
   int ncuts = 0;
   SCIP_ROW* cut;
   SCIP_Bool infeasible;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Real lpval;
   SCIP_Real lpsum;
   int* candidate;         /* at each iteration this array holds the candidates that could got into the second partition */
   int* partition;         /* at each iteration, holds the bins in the second partition */
   int counter;
   int partitionsize;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   assert(nbins > 0);
   assert(ncluster > 0 && ncluster <= nbins);

   /* construct the 2 partition inequalities heuristically */

   /* get the variables from scip */
   edgevars = SCIPspaGetEdgevars(scip);


   SCIP_CALL( SCIPallocClearMemoryArray(scip, &candidate, nbins) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &partition, nbins) );


   for( i = 0; i < nbins; ++i )
   {
      /* construct the set candidate */
      counter = 0;
      for( j = 0; j < i; ++j )
      {
         partition[j] = -1;
         candidate[j] = -1;
         if( j == i )
            continue;
         /* only add bins where there is a fractional uncut value on the edge */
         if( NULL == edgevars[i][j][0][0])
            continue;
         lpval = SCIPvarGetLPSol(edgevars[i][j][0][0]);
         if( lpval > 0 && lpval < 1 && candidate[counter - 1] != j )
         {
            candidate[counter] = j;
            ++counter;
         }
      }
      /* if there are no fractional edges between i and any other bin continue */
      if( 0 == counter )
         continue;
      /* construcht the partition set */
      partition[0] = candidate[0];
      partitionsize = 1;

      for( j = 1; j < counter; ++j )
      {
         SCIP_Bool addtopartition = TRUE;
         /* at first, only add edges with uncut value 0 inside the partition */
         for( w = 1; w < partitionsize; ++w )
         {
            if( candidate[j] < partition[w] && NULL != edgevars[partition[w]][candidate[j]][0][0] && SCIPvarGetLPSol(edgevars[partition[w]][candidate[j]][0][0]) != 0 )
               addtopartition = FALSE;
            if( candidate[j] > partition[w] && NULL != edgevars[candidate[j]][partition[w]][0][0] && SCIPvarGetLPSol(edgevars[candidate[j]][partition[w]][0][0]) != 0 )
               addtopartition = FALSE;

         }
         if( addtopartition )
         {
            partition[partitionsize] = candidate[j];
            partitionsize++;
         }
      }
      assert( nbins > partitionsize - 1 );
      assert( nbins > counter - 1 );
      /* check if the current lpsol violates this constraint. This is the case if the sum over the uncut edges between i and the partition is greater than 1 */
      lpsum = 0;
      for( j = 0; j < partitionsize; ++j )
      {
         if( i > partition[j] && NULL != edgevars[i][partition[j]][0][0] )
            lpsum += SCIPvarGetLPSol(edgevars[i][partition[j]][0][0]);
         else if( NULL != edgevars[partition[j]][i][0][0] )
            lpsum += SCIPvarGetLPSol(edgevars[partition[j]][i][0][0]);
      }
      if( lpsum > 1 )
      {
         (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "2partsize1_%d_%d", i+1, ncuts );
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
         for( j = 0; j < partitionsize; ++j )
         {
            if( i < partition[j] && NULL != edgevars[i][partition[j]][0][0] )
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][partition[j]][0][0], 1.0) );
            else if( NULL != edgevars[partition[j]][i][0][0] )
               SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[partition[j]][i][0][0], 1.0) );
         }
         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
         if( SCIPisCutEfficacious(scip, NULL, cut) )
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
            if( !infeasible )
            {
               SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               *result = SCIP_SEPARATED;
            }
         }
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         ncuts++;
      }
   }
   /* if we did not find a violated constraint yet, redo the above with the modification that the uncut edges in the partition can be non-zero */
   if( *result != SCIP_SEPARATED )
   {
      for( i = 0; i < nbins; ++i )
      {
         /* construct the set w as candidate */
         counter = 0;
         for( j = 0; j < i; ++j )
         {
            partition[j] = -1;
            candidate[j] = -1;
            if( j == i )
               continue;

            lpval = 0;
            if( NULL != edgevars[i][j][0][0] )
               lpval = SCIPvarGetLPSol(edgevars[i][j][0][0]);
            if( lpval > 0 && lpval < 1 && candidate[counter - 1] != j )
            {
               candidate[counter] = j;
               ++counter;
            }

         }
         if( 0 == counter )
            continue;
         else
         {
            partition[0] = candidate[0];
            partitionsize = 1;
         }
         for( j = 1; j < counter; ++j )
         {
            /* we now need that the uncut value within the partition is smaller than that between i and the partition instead of 0 */
            SCIP_Real sum = 0.0;  /* this is the uncut value within the partition */
            for( w = 0; w < partitionsize; ++w )
            {
               if( candidate[j] < partition[w] && NULL != edgevars[partition[w]][candidate[j]][0][0] )
                  sum += SCIPvarGetLPSol(edgevars[partition[w]][candidate[j]][0][0]);
               else if( NULL != edgevars[candidate[j]][partition[w]][0][0] )
                  sum += SCIPvarGetLPSol(edgevars[candidate[j]][partition[w]][0][0]);
            }
            /* lpval is the uncut value between i and the partition */
            lpval = 0;
            if ( candidate[j] < 1 && NULL != edgevars[i][candidate[j]][0][0] )
               lpval += SCIPvarGetLPSol(edgevars[i][candidate[j]][0][0]);
            else if( NULL != edgevars[candidate[j]][i][0][0] )
               lpval += SCIPvarGetLPSol(edgevars[candidate[j]][i][0][0]);
            if( lpval - sum > 0 )
            {
               partition[partitionsize] = candidate[j];
               partitionsize++;
            }
         }
         assert( nbins > partitionsize - 1 );
         assert( nbins > counter - 1 );
         /* check if the current lpsol violates this constraint */
         lpsum = 0;
         for( j = 0; j < partitionsize; ++j )
         {
            if( i > partition[j] && NULL != edgevars[i][partition[j]][0][0] )
               lpsum += SCIPvarGetLPSol(edgevars[i][partition[j]][0][0]);
            else if( NULL != edgevars[partition[j]][i][0][0] )
               lpsum += SCIPvarGetLPSol(edgevars[partition[j]][i][0][0]);
         }
         if( lpsum > 1 )
         {
            (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "2partsize2_%d_&d", i+1, ncuts );
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
            for( j = 0; j < partitionsize; ++j )
            {
               if( i < partition[j] && NULL != edgevars[i][partition[j]][0][0] )
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][partition[j]][0][0], 1.0));
               else if( NULL != edgevars[partition[j]][i][0][0] )
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[partition[j]][i][0][0], 1.0) );
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
            if( SCIPisCutEfficacious(scip, NULL, cut) )
            {
               SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
               if( !infeasible )
               {
                  SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                  *result = SCIP_SEPARATED;
               }
            }
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            ncuts++;
         }
      }
   }

   if( *result != SCIP_SEPARATED )
      *result = SCIP_DIDNOTFIND;

   SCIPfreeMemoryArray(scip, &candidate);
   SCIPfreeMemoryArray(scip, &partition);

   return SCIP_OKAY;
}


/** creates the SparseApprox separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaSparseApprox(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPA* sepa;


   /* include separator */

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpSparseApprox, NULL,
      NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopySparseApprox) );


   return SCIP_OKAY;
}
