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

/**@file   branch_multinode.c
 * @brief  mutlinode branching rule for the set-partitioning part in cycle clustering application.
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "branch_multinode.h"

#include "probdata_cyc.h"
#include "scip/branch_relpscost.h"


#define BRANCHRULE_NAME            "multinode"
#define BRANCHRULE_DESC            "multinode branching creates a child for every variable of a setpartitioning constraint"
#define BRANCHRULE_PRIORITY        10000000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** get the branching candidates viable for multinode branching */
static
SCIP_RETCODE getBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< the address of the branching candidates */
   SCIP_Real*            branchcandssol,     /**< pointer to solution values of the candidates */
   SCIP_Real*            branchcandsfrac,    /**< pointer to fractionalities of the candidates */
   int*                  ncands              /**< number of branching candidates */
   )
{
   SCIP_VAR*** binvars;
   int nbins;
   int ncluster;
   int i;
   int k;

   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   binvars = SCIPcycGetBinvars(scip);

   /* all binvars that are in the lp, and have fractional values are viable candidates */
   for( i = 0; i < nbins; ++i )
   {
      for( k = 0; k < ncluster; ++k )
      {
         if( SCIPvarGetStatus(binvars[i][k]) ==  SCIP_VARSTATUS_COLUMN
            && !SCIPisFeasIntegral(scip, SCIPvarGetLPSol(binvars[i][k])) )
         {
            (branchcands)[*ncands] = binvars[i][k];
            (branchcandssol)[*ncands] = SCIPvarGetLPSol(binvars[i][k]);
            (branchcandsfrac)[*ncands] = MAX(1-(branchcandssol)[*ncands], (branchcandssol)[*ncands]);
            (*ncands)++;
         }
      }
   }

   return SCIP_OKAY;
}

/** branch on a selected bin -> Create at most |Cluster| children */
static
SCIP_RETCODE branchOnBin(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   row,                /**< the row in the binvar-matrix (not lp-row) to be branched on */
   SCIP_RESULT*          result              /**< pointer to store result of branching */
   )
{
   SCIP_VAR*** binvars;
   SCIP_Bool* branched;
   SCIP_Real priority;
   SCIP_Real estimate;
   SCIP_Real minestzero = SCIP_INVALID;
   SCIP_Real tmp;
   SCIP_NODE* node;
   int ncands = 0;
   int k;
   int ncluster;

   binvars = SCIPcycGetBinvars(scip);
   ncluster = SCIPcycGetNCluster(scip);
   assert(NULL != binvars[row]);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &branched, ncluster) );

   /* check all candidates */
   for( k = 0; k < ncluster; ++k )
   {
      if( (SCIPvarGetStatus(binvars[row][k]) == SCIP_VARSTATUS_LOOSE ||
         SCIPvarGetStatus(binvars[row][k]) == SCIP_VARSTATUS_COLUMN) &&
         !SCIPisZero(scip, SCIPvarGetLPSol(binvars[row][k])) && !SCIPisEQ(scip, SCIPvarGetLPSol(binvars[row][k]), 1.0) )
      {
         ncands++;
         priority = SCIPcalcNodeselPriority(scip, binvars[row][k], SCIP_BRANCHDIR_UPWARDS, 1.0);
         estimate = SCIPcalcChildEstimate(scip, binvars[row][k], 1.0);
         tmp = SCIPcalcChildEstimate(scip, binvars[row][k], 0.0);
         minestzero = MIN(tmp, minestzero);

         /* branch all viable candidates upwards */
         SCIP_CALL( SCIPcreateChild(scip, &node, priority, estimate) );
         SCIP_CALL( SCIPchgVarLbNode(scip, node, binvars[row][k], 1.0) );

         branched[k] = TRUE;

         *result = SCIP_BRANCHED;
         if( ncands > 2 )
            break;
      }
   }

   /* create one child, were all the before upwards branched variables are now set to 0. Only do so if at least one
    * upwards branching was done and if not all the variables were branched upwards
    */
   if( ncands > 0 && ncands < ncluster )
   {
      SCIP_CALL( SCIPcreateChild(scip, &node, (SCIP_Real)ncands, minestzero) );
      for( k = 0; k < ncluster; ++k )
      {
         if( branched[k] )
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, node, binvars[row][k], 0.0) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &branched);

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMultinode)
{  /*lint --e{715}*/
   int i;
   int k;
   int nbins;
   int ncluster;
   int ncands;

   SCIP_Real* score;
   SCIP_Real max;
   int maxrow;
   SCIP_VAR*** binvars;
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandssol;
   SCIP_Real* branchcandsfrac;

   binvars = SCIPcycGetBinvars(scip);
   nbins = SCIPcycGetNBins(scip);
   ncluster = SCIPcycGetNCluster(scip);
   *result = SCIP_DIDNOTRUN;

   assert(nbins > 0);
   assert(ncluster > 0 && ncluster <= nbins);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &score, nbins) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &branchcands, (SCIP_Longint) nbins * ncluster) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &branchcandssol, (SCIP_Longint) nbins * ncluster) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &branchcandsfrac, (SCIP_Longint) nbins * ncluster) );

   ncands = 0;
   /* get the candidates */
   SCIP_CALL( getBranchCands( scip, branchcands, branchcandssol, branchcandsfrac, &ncands) );
   if( ncands != 0 )
   {
      /* compute the relpcost for the candidates */
      SCIP_CALL( SCIPexecRelpscostBranching(scip, branchcands, branchcandssol, branchcandsfrac, ncands, FALSE,  result) );

      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM);

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM )
      {
         maxrow = -1;
         max = -SCIPinfinity(scip);

         /* compute the best candidate from pseudocostscore */
         for( i = 0; i < nbins; ++i )
         {
            score[i] = 0;

            for( k = 0; k < ncluster; ++k )
            {
               /* Do not branch on variables that are already fixed locally */
               if( SCIPisEQ(scip, SCIPvarGetUbLocal(binvars[i][k]), SCIPvarGetLbLocal(binvars[i][k])) )
               {
                  score[i] = -SCIPinfinity(scip);
                  break;
               }

               if( SCIPvarGetStatus(binvars[i][k]) == SCIP_VARSTATUS_COLUMN &&
                  !SCIPisZero(scip, SCIPvarGetLPSol(binvars[i][k])) &&
                  !SCIPisEQ(scip, SCIPvarGetLPSol(binvars[i][k]), 1.0) )
               {
                  score[i] += SCIPgetVarPseudocostScore(scip, binvars[i][k], SCIPvarGetLPSol(binvars[i][k]));
               }
            }

            if( SCIPisLT(scip, max, score[i]) && !SCIPisInfinity(scip, -score[i]) )
            {
               max = score[i];
               maxrow = i;
            }
         }

         if( -1 != maxrow )
         {
            SCIP_CALL( branchOnBin(scip, maxrow, result) );
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &score);
   SCIPfreeBufferArray(scip, &branchcands);
   SCIPfreeBufferArray(scip, &branchcandssol);
   SCIPfreeBufferArray(scip, &branchcandsfrac);

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the mutlinode branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMultinode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   branchruledata = NULL;

   /* include branching rule */
   /* use SCIPincludeBranchruleBasic() plus setter functions if you want to set callbacks one-by-one
    * and your code should compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
      BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMultinode) );

   return SCIP_OKAY;
}
