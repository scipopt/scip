/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_strongcoloring.c
 * @brief  branching rule performing strong branching for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements an additional branching rule for the coloring algorithm.
 *
 * We are looking for two nodes v and w, which are not adjacent in the current graph, and consider
 * the following two constraints: SAME(v,w) and DIFFER(v,w). More information about the meaning of
 * these constraints can be found in the documentation of the branching rule in branch_coloring.c.
 *
 * This branching rule puts some more effort into the choice of the two nodes and performs a
 * strongbranching. This means that for every possible choice of two nodes, it solves the LPs of the
 * created children and computes a score with respect to the increase of the lower bound in both
 * nodes. After that, it takes the combination of nodes yielding the best score. The interesting
 * point is that the strongbranching is not performed for each variable, as it is done in some
 * default branching rules of SCIP and supported by the LP-solver, but is done for a constraint,
 * since we are branching on constraints. Look at executeStrongBranching() to see how it is
 * done. There are also some improvements, since testing all possible combination of nodes is very
 * expensive.  The first possibility to avoid this is to stop the computation of scores once a
 * possible branching is found that has only one feasible child. This results in more restrictions
 * in this child without increasing the number of unprocessed nodes.
 *
 * The second improvement is to compute a priority for all possible combinations, w.r.t. the
 * fractional values of the variables. Then, only the first best k combinations are investigated by
 * strongbranching.
 *
 * This code is not optimized and in most cases inferior to the standard branching rule. It is only
 * a demonstration of how to perform strongbranching on constraints!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "branch_strongcoloring.h"
#include "pricer_coloring.h"

#define BRANCHRULE_NAME            "strongcoloring"
#define BRANCHRULE_DESC            "branching rule template"
#define BRANCHRULE_PRIORITY        15000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define EPS 0.0000001

/* default values for parameters */
#define DEFAULT_BRANCHINGMODE 2 
#define DEFAULT_FIXINGSSCOREMODE 3
#define DEFAULT_MAXPRICINGROUNDS -1
#define DEFAULT_USETCLIQUE TRUE
#define DEFAULT_LOOKAHEAD 10



/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   int branchingmode;            /* determines the branchingmode, 0: for fullstrong branching, 
                                    1: strong branching, take first possible branching with only one child-node 
                                    2: strong branching with prior sorting of candidates w.r.t. the fractional value of concerned sets */
   int length;                   /* length of the arrays samevalue and differvalue, length = n*(n-1)/2 with n = NNodes*/
   SCIP_Real* samevalue;         /* value of variables, that would be fixed to 0 for same(i,j), index = nodes2index(i,j) */
   SCIP_Real* differvalue;       /* value of variables, that would be fixed to 0 for differ(i,j), index = nodes2index(i,j) */
   SCIP_Real* combinedvalue;     /* combination of samevalue and differvalue, computed by computeScore*/
   int* permutation;             /* permutation of the indexes of the array combinedvalue, s.t. it is sorted */
   SCIP_Bool usetclique;         /* should the exact pricing with the tclique-algorithm be used for the strongbranchings? */
   int maxpricingrounds;         /* maximal number of pricing rounds used for each probing node in the strongbranching */
   int lookahead;                /* number of candidates to be considered in branchingmode 2 */
   int fixingsscoremode;         /* determines the weightings of the two factors for prior sorting by fractional LP value */
   
};




/*
 * Local methods
 */

/** computes a score for the two improvements that are achieved in the two sons for a branching decision */
static
double computeScore(
   SCIP_Real             val1,               /**< the first value */ 
   SCIP_Real             val2,               /**< the second value */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
      return 0.2 * MAX( val1, val2 ) + 0.8 * MIN( val1, val2 );
}

/** computes a score for the fractional values of the variables that would be fixed to zero for a same- or differ-branching */
static 
SCIP_Real computeFixingsScore(
   SCIP_Real             samevalue,          /**< value of the fractional variables fixed to 0 for a same-branching*/
   SCIP_Real             differvalue,        /**< value of the fractional variables fixed to 0 for a differ-branching*/
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   if ( branchruledata->fixingsscoremode == 1 )
   {
      return 3*samevalue+differvalue; 
   }
   if ( branchruledata->fixingsscoremode == 2 )
   {
      return 2*samevalue+differvalue; 
   }
   if ( branchruledata->fixingsscoremode == 3 )
   {
      return samevalue+10*differvalue; 
   }
   if ( branchruledata->fixingsscoremode == 4 )
   {
      if ( samevalue == -1 && differvalue == -1 )
         return -1;
      return samevalue*differvalue; 
   }
   return samevalue*differvalue; 
}

/** for given nodes node1, node2, compute the corresponding index in the arrays branchruledata->same-/differvalue */
static 
int nodes2index(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node1,              /**< the first node */
   int                   node2               /**< the second node */
   )
{
   int ind;
   int nnodes;
   int i;

   assert(scip != NULL);
   assert(node1 >= 0 && node2 >= 0);

   /* node 1 has to be smaller than node 2 */
   if ( node1 > node2 )
   {
      ind = node1;
      node1 = node2;
      node2 = ind;
   }
   nnodes = COLORprobGetNNodes(scip);
   assert(node1 < nnodes && node2 < nnodes);
   ind = 0;
   for ( i = 0; i < node1; i++ )
      ind += (nnodes - i - 1);
   ind += ( node2-node1-1); 
   return ind;
}

/** for given index of the arrays branchruledata->same-/differvalue, compute the two nodes, the index represents */
static 
void index2nodes(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   ind,                /**< the given index in the arrays */
   int*                  node1,              /**< return value: the first node */
   int*                  node2               /**< return value: the second node */
   )
{
   int nnodes;
   int value;

   assert(scip != NULL);
   assert(node2 != NULL && node1 != NULL);

   nnodes = COLORprobGetNNodes(scip);
   *node1 = 0;
   value = 0;
   while ( value + nnodes - 1 - *node1 <= ind )
   {
      value += (nnodes - 1 - *node1);
      *node1 = *node1 + 1;
   }
   *node2 = *node1 + 1 + (ind - value);
}

/** computes for each pair of nodes (i,j) two values, one for same (i,j), the other for differ(i,j) which are the sum of 
    the values of variables with fractional parts, that would be fixed for this decision 
    asd */
static
SCIP_RETCODE computeBranchingPriorities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< the data of the branching rule */
   )
{
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   TCLIQUE_GRAPH* graph;
   int nlpcands;
   int i;
   int j;
   int k;
   int node1;
   int node2;
   SCIP_VAR* var;
   int setindex;
   int* set;
   int setlength;
   int nnodes;

   assert(scip != NULL);
   assert(branchruledata != NULL);

   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, &nlpcands, NULL) );
   nnodes = COLORprobGetNNodes(scip);
   graph = COLORconsGetCurrentGraph(scip);

   assert(graph != NULL);
   assert(nnodes >= 0);

   /* fill array samevalue, differvalue with zeroes, or -1 for impossible branchings */
   for ( i = 0; i < branchruledata->length; i++ )
   {
      index2nodes(scip, i, &node1, &node2);
      /* there is an edge between node1 and node2 --> no branching possible --> set value to -1 */
      if ( tcliqueIsEdge(graph, node1, node2) )
      {
         branchruledata->samevalue[i] = -1;
         branchruledata->differvalue[i] = -1;
         continue;
      }
      branchruledata->samevalue[i] = 0;
      branchruledata->differvalue[i] = 0;
   }

   /* for all branching candidates (variables with fractional value) check for which branching decisions they would be 
      fixed to 0 and add the fractional part to the related entry in the array samevalue or differvalue */
   for ( i = 0; i < nlpcands; i++ )
   {
      assert(SCIPisFeasPositive(scip, lpcandsfrac[i]));
      var = lpcands[i];
      setindex = (int)(size_t) SCIPvarGetData(var);
      COLORprobGetStableSet(scip, setindex, &set, &setlength);
      for ( j = 0; j < setlength; j++ )
      {
         node1 = set[j];
         /* if node1 is part of a union and not its representant, continue */
         if ( COLORconsGetRepresentative(scip, node1) != node1 )
         {
            continue;
         }
         k = 0;
         for ( node2 = nnodes-1; node2 >= 0; node2-- )
         {
            /* if k is a node, which is part of, but not representant of a union, increment k */
            while ( k < setlength && COLORconsGetRepresentative(scip, set[k]) != set[k] )
            {
               k++;
            }
            /* node1 is equal to node2 -> increment k and continue */
            if ( node2 == node1 )
            {
               assert(k == j);
               k++;
               continue;
            }
            /* if node2 is part of a union and not its representant, continue */
            if ( COLORconsGetRepresentative(scip, node2) != node2 )
               continue;
            /* if there is an edge between node1 and node2 in the current graph, continue */
            if ( branchruledata->differvalue[nodes2index(scip, node1, node2)] == -1 )
            {
               continue;
            }
            /* node2 is also in the set --> the variable would be fixed to 0 for differ(node1, node2) */
            if ( k < setlength && node2 == set[k] )
            {
               branchruledata->differvalue[nodes2index(scip, node1, node2)] += lpcandsfrac[i];
               assert(COLORprobIsNodeInStableSet(scip, setindex, node1) && COLORprobIsNodeInStableSet(scip, setindex, node2));
               k++;
            }
            /* node2 is not in the set --> the variable would be fixed to 0 for same(node1, node2) */
            else
            {
               branchruledata->samevalue[nodes2index(scip, node1, node2)] += lpcandsfrac[i];
               assert(COLORprobIsNodeInStableSet(scip, setindex, node1) && !COLORprobIsNodeInStableSet(scip, setindex, node2));
            }
         }
         assert(k == setlength);
      }
   }

   return SCIP_OKAY;

}



/** computes the lower bound that would a child node with the given branching decision would have */
static 
SCIP_Real executeStrongBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   COLOR_CONSTYPE        constype,           /**< the type of the contraint: SAME or DIFFER */
   int                   node1,              /**< the first node for the branching constraint */
   int                   node2,              /**< the second node for the branching constraint */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< the data of the branching rule */
   )
{

   SCIP_NODE* newnode;
   SCIP_CONS* currentcons;
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   
   SCIP_Real newLb;

   assert(scip != NULL);

   /* get the constraint of the current Node in the B&B-Tree */
   currentcons = COLORconsGetActiveStoreGraphCons(scip);

   /* start Probing */
   SCIPstartProbing(scip);
 
  /* create new probing node and add store graph cons to it with same(node1, node2) */
   SCIPnewProbingNode(scip);
   newnode = SCIPgetCurrentNode(scip);
   SCIP_CALL( COLORcreateConsStoreGraph(scip, &cons, "probingcons", currentcons, constype, node1, node2, newnode) );
   SCIP_CALL( SCIPaddConsNode(scip, newnode, cons, NULL) ); 
   /* propagate the new b&b-node, i.e. fix vars to 0 that don't contain both node1 and node2 */
   SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );
   /* solve the LP using pricing */
   SCIP_CALL( SCIPsolveProbingLPWithPricing(scip, FALSE, FALSE, branchruledata->maxpricingrounds, &cutoff) );
   assert(!cutoff);
   /* get the changed objective value */
   newLb = SCIPgetLPObjval(scip);

   SCIP_CALL( SCIPdelCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPendProbing(scip);

   return newLb;
}


/** index comparison method two values in a real array */
static
SCIP_DECL_SORTINDCOMP(consdataCompValues)
{  
   SCIP_Real* values;
   
   values = (SCIP_Real*)dataptr;

   assert(values != NULL);

   if ( values[ind1] > values[ind2] )
   {
      return -1;
   }
   if ( values[ind1] < values[ind2] )
   {
      return 1;
   }
   return 0;
}


/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyStrongcoloring)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpStrongcoloring)
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

   int i;
   int j;
   int nnodes;

   SCIP_Bool* wasnode1;
   SCIP_Bool* wasnode2;
   SCIP_Bool start;
   TCLIQUE_GRAPH* graph;
   SCIP_Real currLb;
   SCIP_Real sameLb;
   SCIP_Real differLb;

   SCIP_Real bestscore;
   SCIP_Real bestdiffer;
   SCIP_Real bestsame;
   SCIP_Real score;
   int bestnode2;
   int bestnode1;

   SCIP_BRANCHRULEDATA* branchruledata;

#ifndef NDEBUG
   SCIP_NODE* node;
#endif

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   nnodes = COLORprobGetNNodes(scip);
   graph = COLORconsGetCurrentGraph(scip);

   if ( branchruledata->branchingmode == 2 )
   {
      SCIP_CALL( computeBranchingPriorities(scip, branchruledata) );

      for ( i = 0; i < branchruledata->length; i++ )
      {
         branchruledata->combinedvalue[i] = computeFixingsScore(branchruledata->samevalue[i], branchruledata->differvalue[i], branchruledata);
      }
      /* get permutation of indexes, so that the array is sorted */
      /** @todo could be improved by only getting the k best indexes */
      SCIPsort(branchruledata->permutation, consdataCompValues, branchruledata->combinedvalue, branchruledata->length);

      bestscore = -1;
      bestnode1 = -1;
      bestnode2 = -1;
      bestdiffer = -1;
      bestsame = -1;

      for ( i = 0; i < branchruledata->lookahead && i < branchruledata->length; i++ )
      {
         index2nodes(scip, branchruledata->permutation[i], &node1, &node2);
         currLb = SCIPgetLPObjval(scip);

         /* SAME */
         sameLb = executeStrongBranching(scip, COLOR_CONSTYPE_SAME, node1, node2, branchruledata);
         if ( sameLb-currLb > 1000 )
         {
            sameLb = currLb + 1000;
         }
         
         /* DIFFER */
         differLb = executeStrongBranching(scip, COLOR_CONSTYPE_DIFFER, node1, node2, branchruledata);
         if ( differLb-currLb > 1000 )
         {
            differLb = currLb + 1000;
         }

         score = computeScore( sameLb-currLb, differLb-currLb, branchruledata );
         assert( !SCIPisFeasZero(scip, score) || (SCIPisFeasZero(scip, 0.2 * (sameLb-currLb)) && SCIPisFeasZero(scip, 0.2 * (differLb-currLb))
               && (SCIPisFeasZero(scip, sameLb-currLb) || SCIPisFeasZero(scip, differLb-currLb))) );

         if ( score > bestscore )
         {
            bestscore = score;
            bestnode1 = node1;
            bestnode2 = node2;
            bestdiffer = differLb-currLb;
            bestsame = sameLb-currLb;
         }
         if ( bestdiffer > 999 || bestsame > 999 )
         {
            break;
         }
      }

   }
   else
   {
      assert(branchruledata->branchingmode == 0 || branchruledata->branchingmode == 1);
      /* create array wasnode1 and wasnode2 and fill them with FALSE */
      SCIP_CALL( SCIPallocBufferArray(scip, &wasnode1, nnodes) );
      BMSclearMemoryArray(wasnode1, nnodes);
      SCIP_CALL( SCIPallocBufferArray(scip, &wasnode2, nnodes) );

      bestscore = -1;
      bestnode1 = -1;
      bestnode2 = -1;
      bestdiffer = -1;
      bestsame = -1;

      SCIPsetBoolParam(scip, "pricers/coloring/usetclique", branchruledata->usetclique);
#ifndef NDEBUG
      node = SCIPgetCurrentNode(scip);
#endif
      currentcons = COLORconsGetActiveStoreGraphCons(scip);
       
      start = TRUE;
      for ( i = SCIPgetDepth(scip)%nnodes; (start || (i != SCIPgetDepth(scip)%nnodes)); i=((i+1)%nnodes) )
      {
         start = FALSE;
         node1 = COLORconsGetRepresentative(scip, i);
         /* check whether node1 was already tested */
         if ( wasnode1[node1] == TRUE )
         {
            continue;
         }
         else
         {
            wasnode1[node1] = TRUE;
         }
         BMSclearMemoryArray(wasnode2, nnodes);

         for ( j = i+1; j < nnodes; j++ )
         {
            node2 = COLORconsGetRepresentative(scip, j);
            if ( node2 == node1 || tcliqueIsEdge(graph, node1, node2) || node2 < i )
            {
               continue;
            }
            else
            {
               /* check whether node2 was already tested */
               if ( wasnode2[node2] == TRUE ) continue;
               else wasnode2[node2] = TRUE;
            
               currLb = SCIPgetLPObjval(scip);

               assert(currentcons == COLORconsGetActiveStoreGraphCons(scip));
               assert(node == SCIPgetCurrentNode(scip));

               /* compute lower bounds for possible branchings */

               /* SAME */
               sameLb = executeStrongBranching(scip, COLOR_CONSTYPE_SAME, node1, node2, branchruledata);
               if ( sameLb-currLb > 1000 )
               {
                  sameLb = currLb + 1000;
               }

               /* DIFFER */
               differLb = executeStrongBranching(scip, COLOR_CONSTYPE_DIFFER, node1, node2, branchruledata);
               if ( differLb-currLb > 1000 )
               {
                  differLb = currLb + 1000;
               }

               score = computeScore( sameLb-currLb, differLb-currLb, branchruledata );
               if ( score > bestscore )
               {
                  bestscore = score;
                  bestnode1 = node1;
                  bestnode2 = node2;
                  bestdiffer = differLb-currLb;
                  bestsame = sameLb-currLb;
               }
               if ( (branchruledata->branchingmode == 1) && (bestdiffer > 999 || bestsame > 999) )
               {
                  break;
               }

            }
         }
         if ( (branchruledata->branchingmode == 1) && (bestdiffer > 999 || bestsame > 999) )
         {
            break;
         }
      }
   
      SCIPsetBoolParam(scip, "pricers/coloring/usetclique", TRUE);
      assert(node == SCIPgetCurrentNode(scip));
      assert(currentcons == COLORconsGetActiveStoreGraphCons(scip));

      SCIPfreeBufferArray(scip, &wasnode2);
      SCIPfreeBufferArray(scip, &wasnode1);

   }
      
   assert(!SCIPisSumNegative(scip, bestscore));
  
   node1 = bestnode1;
   node2 = bestnode2;

   /* branchingmode >= 1 --> only create nodes, that do not have a LP solution that is much bigger than the lower bound */
   if ( branchruledata->branchingmode >= 1 && branchruledata->usetclique == TRUE )
   {
      *result = SCIP_CUTOFF;
      currentcons = COLORconsGetActiveStoreGraphCons(scip);

      if ( bestdiffer <= 999 )
      {
         /* create the b&b-tree child-nodes of the current node */
         SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );
      
         /* create corresponding constraints */
         SCIP_CALL( COLORcreateConsStoreGraph(scip, &consdiffer, "differ", currentcons, COLOR_CONSTYPE_DIFFER, node1, node2, childdiffer) );
      
         /* add constraints to nodes */
         SCIP_CALL( SCIPaddConsNode(scip, childdiffer, consdiffer, NULL) );
      
         /* release constraints */
         SCIP_CALL( SCIPreleaseCons(scip, &consdiffer) );
      
         *result = SCIP_BRANCHED;
      }
   
      if ( bestsame <= 999 )
      {
         /* create the b&b-tree child-nodes of the current node */
         SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
      
         /* create corresponding constraints */
         SCIP_CALL( COLORcreateConsStoreGraph(scip, &conssame,   "same",   currentcons, COLOR_CONSTYPE_SAME,   node1, node2, childsame) );
      
         /* add constraints to nodes */
         SCIP_CALL( SCIPaddConsNode(scip, childsame, conssame, NULL) );
      
         /* release constraints */
         SCIP_CALL( SCIPreleaseCons(scip, &conssame) );
      
         *result = SCIP_BRANCHED;
      }
   }
   /* create both children */
   else
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
   }

   return SCIP_OKAY;
}


/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeStrongcoloring)
{  
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitStrongcoloring)
{  
   SCIP_BRANCHRULEDATA* branchruledata;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get memory for the arrays */
   branchruledata->length = (COLORprobGetNNodes(scip)*(COLORprobGetNNodes(scip)-1))/2;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(branchruledata->samevalue), branchruledata->length) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(branchruledata->differvalue), branchruledata->length) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(branchruledata->combinedvalue), branchruledata->length) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(branchruledata->permutation), branchruledata->length) );

   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitStrongcoloring)
{  
   SCIP_BRANCHRULEDATA* branchruledata;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* free arrays */
   SCIPfreeMemoryArray(scip, &(branchruledata->samevalue));
   SCIPfreeMemoryArray(scip, &(branchruledata->differvalue));
   SCIPfreeMemoryArray(scip, &(branchruledata->combinedvalue));
   SCIPfreeMemoryArray(scip, &(branchruledata->permutation));
   
   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the coloring branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleStrongcoloring(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   assert(scip != NULL);

   /* create branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   branchrule = NULL;
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
	 BRANCHRULE_MAXBOUNDDIST, branchruledata) );
   assert(branchrule != NULL);
   
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyStrongcoloring) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeStrongcoloring) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpStrongcoloring) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitStrongcoloring) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitStrongcoloring) );


   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/strongcoloring/lookahead",
         "number of candidates to be considered in branchingmode 2",
         &branchruledata->lookahead, TRUE, DEFAULT_LOOKAHEAD, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/strongcoloring/usetclique",
         "should the exact pricing with the tclique-algorithm be used for the strongbranchings?",
         &branchruledata->usetclique, FALSE, DEFAULT_USETCLIQUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/strongcoloring/maxpricingrounds",
         "maximal number of pricing rounds used for each probing node in the strongbranching",
         &branchruledata->maxpricingrounds, TRUE, DEFAULT_MAXPRICINGROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/strongcoloring/branchingmode",
         "determines the branchingmode, 0: fullstrong branching, 1: strong branching, take first possible branching with only one child-node, 2: strong branching with prior sorting of candidates w.r.t. the fractional value of concerned sets */",
         &branchruledata->branchingmode, FALSE, DEFAULT_BRANCHINGMODE, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/strongcoloring/fixingsscoremode",
         "determines the weightings of the two factors for prior sorting by fractional LP value",
         &branchruledata->fixingsscoremode, TRUE, DEFAULT_FIXINGSSCOREMODE, 0, 4, NULL, NULL) );

   return SCIP_OKAY;
}
