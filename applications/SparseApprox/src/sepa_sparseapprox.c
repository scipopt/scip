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
#define SEPA_PRIORITY               1000000
#define SEPA_FREQ                   1.0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */



/*
 * Callback methods of separator
 */


/** copy method for separator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_SEPACOPY(sepaCopySparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaCopySparseApprox NULL
#endif

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeSparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeSparseApprox NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitSparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitSparseApprox NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitSparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitSparseApprox NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolSparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolSparseApprox NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolSparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolSparseApprox NULL
#endif


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpSparseApprox)
{
   SCIP_VAR***** edgevars;
   int i;               /* the bin in the first partition */
   int j;
   int w;
   int k;
   int nbins;
   int ncluster;
   SCIP_CONS* cons;
   char consname[SCIP_MAXSTRLEN];
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


   /* construct the 2 partition inequalities heuristically */

   /* get the variables from scip */
   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);

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
         lpval = 0;

         /* only add bins where there is a fractional uncut value on the edge */
         for( k = 0; k < ncluster; ++k )
         {
            lpval += SCIPvarGetLPSol(edgevars[i][j][k][k]);
         }
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
      else
      {
         partition[0] = candidate[0];
         partitionsize = 1;
      }
      for( j = 1; j < counter; ++j )
      {
         SCIP_Bool addtopartition = TRUE;
         /* at first, only add edges with uncut value 0 inside the partition */
         for( w = 1; w < partitionsize; ++w )
         {
            for( k = 0; k < ncluster; ++k )
            {
               if( candidate[j] < partition[w] && SCIPvarGetLPSol(edgevars[partition[w]][candidate[j]][k][k]) != 0 )
                  addtopartition = FALSE;
               if( candidate[j] > partition[w] && SCIPvarGetLPSol(edgevars[candidate[j]][partition[w]][k][k]) != 0 )
                  addtopartition = FALSE;
            }

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
         for( k = 0; k < ncluster; ++k )
         {
            if( i > partition[j])
               lpsum += SCIPvarGetLPSol(edgevars[i][partition[j]][k][k]);
            else
               lpsum += SCIPvarGetLPSol(edgevars[partition[j]][i][k][k]);
         }
      }
      if( lpsum > 1 )
      {
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "2partsize1_%d", i+1 );
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(scip), 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         for( j = 0; j < partitionsize; ++j )
         {
            for( k = 0; k < ncluster; ++k )
            {
               if( i < partition[j] )
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons, edgevars[i][partition[j]][k][k], 1.0) );
               else
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons, edgevars[partition[j]][i][k][k], 1.0) );
            }
         }
         SCIP_CALL( SCIPaddCons(scip, cons ) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         *result = SCIP_CONSADDED;
      }
   }
   /* if we did not find a violated constraint yet, redo the above with the modification that the uncut edges in the partition can be non-zero */
   if( *result != SCIP_CONSADDED )
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
            for( k = 0; k < ncluster; ++k )
            {
               lpval += SCIPvarGetLPSol(edgevars[i][j][k][k]);
            }
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
               for( k = 0; k < ncluster; ++k )
               {
                  if( candidate[j] < partition[w] )
                     sum += SCIPvarGetLPSol(edgevars[partition[w]][candidate[j]][k][k]);
                  else
                     sum += SCIPvarGetLPSol(edgevars[candidate[j]][partition[w]][k][k]);
               }
            }
            /* lpval is the uncut value between i and the partition */
            lpval = 0;
            for( k = 0; k < ncluster; ++k )
            {
               if ( candidate[j] < 1 )
                  lpval += SCIPvarGetLPSol(edgevars[i][candidate[j]][k][k]);
               else
                  lpval += SCIPvarGetLPSol(edgevars[candidate[j]][i][k][k]);
            }
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
            for( k = 0; k < ncluster; ++k )
            {
               if( i > partition[j])
                  lpsum += SCIPvarGetLPSol(edgevars[i][partition[j]][k][k]);
               else
                  lpsum += SCIPvarGetLPSol(edgevars[partition[j]][i][k][k]);
            }
         }
         if( lpsum > 1 )
         {
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "2partsize1_%d", i+1 );
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, 0.0, 1.0, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            for( j = 0; j < partitionsize; ++j )
            {
               for( k = 0; k < ncluster; ++k )
               {
                  if( i < partition[j] )
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, edgevars[i][partition[j]][k][k], 1.0) );
                  else
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, edgevars[partition[j]][i][k][k], 1.0) );
               }
            }
            SCIP_CALL( SCIPaddCons(scip, cons ) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            *result = SCIP_CONSADDED;
         }
      }
   }
   if( *result != SCIP_CONSADDED )
      *result = SCIP_DIDNOTFIND;

   SCIPfreeMemoryArray(scip, &candidate);
   SCIPfreeMemoryArray(scip, &partition);

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolSparseApprox)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of SparseApprox separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExecsolSparseApprox NULL
#endif


/*
 * separator specific interface methods
 */

/** creates the SparseApprox separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaSparseApprox(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create SparseApprox separator data */
   sepadata = NULL;
   /* TODO: (optional) create separator specific data here */

   sepa = NULL;

   /* include separator */
#if 0
   /* use SCIPincludeSepa() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaCopySparseApprox, sepaFreeSparseApprox, sepaInitSparseApprox, sepaExitSparseApprox, sepaInitsolSparseApprox, sepaExitsolSparseApprox, sepaExeclpSparseApprox, sepaExecsolSparseApprox,
      sepadata) );
#else
   /* use SCIPincludeSepaBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpSparseApprox, sepaExecsolSparseApprox,
      sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopySparseApprox) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeSparseApprox) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitSparseApprox) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitSparseApprox) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolSparseApprox) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolSparseApprox) );
#endif

   /* add SparseApprox separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
