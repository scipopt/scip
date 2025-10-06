/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_coloring.c
 * @brief  default branching rule for the vertex coloring problem
 * @author Gerald Gamrath
 * @author Julian Meffert
 *
 * This file implements the standard branching rule for the coloring algorithm.
 *
 * As we use column generation, we may not branch on the variables themselves,
 * but on some sort of constraints that we introduce in the pricing problem.
 *
 * In our case, we choose two nodes v and w, which are not adjacent in the current graph, and
 * consider the following two constraints: SAME(v,w) and DIFFER(v,w).  SAME(v,w) requires that both
 * nodes v and w get the same color, whereas DIFFER(v,w) forbids this. For each pair of nodes, each
 * feasible solution fulfills exactly one of these constraints. Hence, splitting the solution space
 * into two parts, one fulfilling SAME(v,w) and the other DIFFER(v,w), does not cut off any feasible
 * solution and can therefore be used as the branching rule.
 *
 * The branching is done as follows: Given the optimal (fractional) solution of the current
 * branch-and-bound node, choose the least/most fractional variable and the corresponding stable set
 * s1. Now choose two nodes v, w and another stable set s2, such that v is part of both stable sets,
 * whereas w is part of exactly one of the stable sets.  Create two children of the current node,
 * one with the restriction SAME(v,w), the other one with restriction DIFFER(v,w). Therefore, each
 * node gets a constraint of type @c cons_storeGraph, which enforces the branching decision and
 * assures that each coloring of the nodes in the respective subgraph assigns to both nodes the same
 * color/different colors by fixing stable sets to 0 that violate this constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "branch_coloring.h"

#define BRANCHRULE_NAME            "coloring"
#define BRANCHRULE_DESC            "branching rule template"
#define BRANCHRULE_PRIORITY        50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define BRANCHRULE_STRATEGIES          "ml" /**< possible variable selection strategies m=most fractional, l=least fractional */
#define BRANCHRULE_STRATEGY_DEFAULT     'l' /**< default variable selection strategy */


/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   char strategy;                /* determines the variable selection,
                                    l: for least fractional variable,
                                    m: for most fractional variable */
};


/*
 * Callback methods of branching rule
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpColoring)
{
   /* array of candidates for branching + fractionalities of candidates + length of array */
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   /* variables for finding the most fractional column */
   SCIP_Real fractionality;
   SCIP_Real bestfractionality;
   int bestcand;
   /* array of variables in a constraint + length of array */
   SCIP_VAR** vars;
   int nvars;
   /* the variables for 2 stable sets, needed to find the two nodes for branching */
   SCIP_VAR* s1;
   SCIP_VAR* s2;
   /* the 2 stable sets: array with all nodes and arraylength for each of them */
   int* set1;
   int setlength1;
   int* set2;
   int setlength2;
   /* the 2 nodes, for which the branching is done by DIFFER and SAME */
   int node1;
   int node2;
   /* the constraint belonging to node1 */
   SCIP_CONS* cons1;
   /* the nodes in the branch&bound-tree which are created */
   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   /* the constraints for the created b&b-nodes */
   SCIP_CONS* conssame;
   SCIP_CONS* consdiffer;
   /* the constraint of the processed b&b-node */
   SCIP_CONS* currentcons;
   /* the variable selection strategy */
   char strategy;
   /* the branching rule data */
   SCIP_BRANCHRULEDATA*  branchruledata;

   int i;
   int j;
   int k;
   int l;
   int setindex;

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   branchruledata = SCIPbranchruleGetData(branchrule);
   strategy = branchruledata->strategy;

   bestcand = -1;

   switch( strategy )
   {
   case 'l':
      /* search the least fractional candidate */
      bestfractionality = 1;
      for( i = 0; i < nlpcands; ++i )
      {
         assert(lpcands[i] != NULL);
         fractionality = lpcandsfrac[i];
         fractionality = MIN( fractionality, 1.0-fractionality );
         if ( fractionality < bestfractionality )
         {
            bestfractionality = fractionality;
            bestcand = i;
         }
      }
      break;

   case 'm':
      /* search the most fractional candidate */
      bestfractionality = 0;
      for( i = 0; i < nlpcands; ++i ) {
         assert(lpcands[i] != NULL);
         fractionality = lpcandsfrac[i];
         fractionality = MIN( fractionality, 1.0-fractionality );
         if ( fractionality > bestfractionality )
         {
            bestfractionality = fractionality;
            bestcand = i;
         }
      }
      break;
   default:
      SCIPABORT();
      return SCIP_PARAMETERWRONGVAL;
   }

   assert(bestcand >= 0);
   assert(SCIPisFeasPositive(scip, bestfractionality));

   /* s1 = column belonging to bestcand */
   s1 = lpcands[bestcand];
   setindex = (int)(size_t) SCIPvarGetData(s1);

   /* get stable set corresponding to variable s1 */
   COLORprobGetStableSet(scip, setindex, &set1, &setlength1);

   node1 = -1;
   node2 = -1;
   s2 = NULL;
   /* search for two nodes node1, node2 and column s2 (s2 != s1) such that:
      the node1-constraint is covered by s1 and s2
      the node2-constraint is covered by exactly one of the columns s1,s2 */
   for ( i = 0; ((i < setlength1) && (node2 == -1)); i++ )
   {
      node1 = COLORconsGetRepresentative(scip, set1[i]);
      /* search for other set containing the node */
      cons1 = COLORprobGetConstraint(scip, node1);
      vars = SCIPgetVarsSetppc(scip, cons1);
      nvars = SCIPgetNVarsSetppc(scip, cons1);
      for ( j  = 0; j < nvars; j++ )
      {
         if ( vars[j] != s1 && !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[j])) )
         {
            s2 = vars[j];
            setindex = (int)(size_t) SCIPvarGetData(s2);
            /* get Stable Set corresponding to Variable s2 */
            COLORprobGetStableSet(scip, setindex, &set2, &setlength2);
            /* for all nodes in set1 */
            for ( k = 0; k < setlength1; k++ )
            {
               /* set node2 = current node in set1 */
               node2 = COLORconsGetRepresentative(scip, set1[k]);
               if ( node2 == node1)
               {
                  node2 = -1;
               }
               else
               {
                  /* check whether node2 is in set2 */
                  for ( l = 0; l < setlength2; l++ )
                  {
                     if ( COLORconsGetRepresentative(scip, set2[l]) == node2 )
                     {
                        /* node2 is in both sets -> no branching-candidate */
                        node2 = -1;
                        break;  /* for l */
                     }
                  }
                  /* if node2 found, get out of for-loops */
                  if ( node2 != -1 )
                  {
                     break; /* for k */
                  }
               }
            }
            if ( node2 != -1 )
            {
               break;  /* for j */
            }
            for ( k = 0; k < setlength2; k++ )
            {
               /* set node2 = current node in set1 */
               node2 = COLORconsGetRepresentative(scip, set2[k]);
               if ( node2 == node1)
               {
                  node2 = -1;
               }
               else
               {
                  /* check whether node2 is in set2 */
                  for ( l = 0; l < setlength1; l++ )
                  {
                     if ( COLORconsGetRepresentative(scip, set1[l]) == node2 )
                     {
                        /* node2 is in both sets -> no branching-candidate */
                        node2 = -1;
                        break;  /* for l */
                     }
                  }
                  /* if node2 found, get out of for-loops */
                  if ( node2 != -1 )
                  {
                     break; /* for k */
                  }
               }
            }
            if ( node2 != -1 )
            {
               break;  /* for j */
            }
         }
      }
   }

   assert(node2 != -1);
   assert(node1 != -1);
   assert(node1 == COLORconsGetRepresentative(scip, node1));
   assert(node2 == COLORconsGetRepresentative(scip, node2));
   assert(SCIPconsIsEnabled(COLORprobGetConstraint(scip,node1)));
   assert(SCIPconsIsEnabled(COLORprobGetConstraint(scip,node2)));
   assert(!tcliqueIsEdge(COLORconsGetCurrentGraph(scip), node1, node2));

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* create corresponding constraints */
   currentcons = COLORconsGetActiveStoreGraphCons(scip);
   SCIP_CALL( COLORcreateConsStoreGraph(scip, &conssame,   "same",   currentcons, COLOR_CONSTYPE_SAME,   node1, node2, childsame) );
   SCIP_CALL( COLORcreateConsStoreGraph(scip, &consdiffer, "differ", currentcons, COLOR_CONSTYPE_DIFFER, node1, node2, childdiffer) );

   /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childsame, conssame, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdiffer, consdiffer, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &conssame) );
   SCIP_CALL( SCIPreleaseCons(scip, &consdiffer) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}/*lint !e715*/


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsColoring)
{
   /* the 2 nodes, for which the branching is done by DIFFER and SAME */
   int node1;
   int node2;
   /* the nodes in the branch&bound-tree which are created */
   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   /* the constraints for the created b&b-nodes */
   SCIP_CONS* conssame;
   SCIP_CONS* consdiffer;
   /* the constraint of the processed b&b-node */
   SCIP_CONS* currentcons;

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* search for two nodes node1, node2 such that:
      node1 and node2 are neither in the same union nor adjacent */
   for ( node1 = 0; node1 < COLORprobGetNNodes(scip); ++node1 )
   {
      if ( node1 != COLORconsGetRepresentative(scip, node1) )
      {
         continue;
      }
      for ( node2 = node1+1; node2 < COLORprobGetNNodes(scip); ++node2 )
      {
         if ( node2 != COLORconsGetRepresentative(scip, node2) )
         {
            continue;
         }
         if ( (node2 != node1) && !tcliqueIsEdge(COLORconsGetCurrentGraph(scip), node1, node2))
         {
            /* create the b&b-tree child-nodes of the current node */
            SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
            SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

            /* create corresponding constraints */
            currentcons = COLORconsGetActiveStoreGraphCons(scip);
            SCIP_CALL( COLORcreateConsStoreGraph(scip, &conssame,   "same",   currentcons, COLOR_CONSTYPE_SAME,   node1, node2, childsame) );
            SCIP_CALL( COLORcreateConsStoreGraph(scip, &consdiffer, "differ", currentcons, COLOR_CONSTYPE_DIFFER, node1, node2, childdiffer) );

            /* add constraints to nodes */
            SCIP_CALL( SCIPaddConsNode(scip, childsame, conssame, NULL) );
            SCIP_CALL( SCIPaddConsNode(scip, childdiffer, consdiffer, NULL) );

            /* release constraints */
            SCIP_CALL( SCIPreleaseCons(scip, &conssame) );
            SCIP_CALL( SCIPreleaseCons(scip, &consdiffer) );

            *result = SCIP_BRANCHED;

            return SCIP_OKAY;
         }
      }
   }

   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyColoring)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   return SCIP_OKAY;
}


/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeColoring)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeBlockMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}/*lint !e715*/



/*
 * branching rule specific interface methods
 */

/** creates the coloring branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleColoring(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   assert(scip != NULL);

   /* create branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   branchrule = NULL;
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
	 BRANCHRULE_MAXBOUNDDIST, branchruledata) );
   assert(branchrule != NULL);

   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpColoring) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsColoring) );
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyColoring) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeColoring) );

   SCIP_CALL( SCIPaddCharParam(scip, "branching/" BRANCHRULE_NAME "/strategy",
         "variable selection strategy, 'l'east fractional or 'm'ost fractional variable",
         &branchruledata->strategy, FALSE, BRANCHRULE_STRATEGY_DEFAULT, BRANCHRULE_STRATEGIES, NULL, NULL) );

   return SCIP_OKAY;
}
