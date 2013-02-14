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

/**@file   pricer_coloring.c
 * @brief  variable pricer for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements the pricer for the coloring algorithm.
 *
 * It computes maximal stable sets in the current graph whose corresponding variables can improve
 * the current LP solution.  This is done by computing a maximum weighted stable set in the current
 * graph with dual-variables of the node constraints as weights. A variable can improve the
 * solution, if the weight of the corresponding stable set is larger than 1, since it then has
 * negative reduced costs, which are given by (1 - weight of the set).
 *
 * The pricer first tries to compute such a stable set using a a greedy-method. If it fails, the tclique-algorithm is
 * used on the complementary graph. This is a branch-and-bound based algorithm for maximal cliques,
 * included in SCIP.  In this case, not only the best solution is added to the LP, but also all other
 * stable sets found during the branch-and-bound process that could improve the current LP solution
 * are added, limited to a maximal number that can be changed by a parameter.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "pricer_coloring.h"
#include "reader_col.h"
#include "cons_storeGraph.h"
#include <stdio.h>
#include <stdlib.h>


#define PRICER_NAME            "coloring"
#define PRICER_DESC            "pricer for coloring"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

/* defines for rounding for tclique */
#define MAXDNOM                1000LL 
#define MINDELTA               1e-03 
#define MAXDELTA               1e-09 
#define MAXSCALE               1000.0 


/* default values for parameters */
#define DEFAULT_MAXVARSROUND     0
#define DEFAULT_USETCLIQUE      TRUE
#define DEFAULT_USEGREEDY       TRUE
#define DEFAULT_ONLYBEST        TRUE
#define DEFAULT_MAXROUNDSROOT   -1
#define DEFAULT_MAXROUNDSNODE   -1
#define DEFAULT_MAXTCLIQUENODES INT_MAX



/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   SCIP*            scip;                    /* SCIP data structure */
   int              maxvarsround;            /* maximal number of variables created each round */
   int              oldmaxvarsround;         /* maximal number of variables created each round, old value before parameter was changed */
   int              nstablesetsfound;        /* number of improving stable sets found in the current round so far */
   SCIP_CONS**      constraints;             /* array containing all node constraints */
   SCIP_Real        scalefactor;             /* the factor used for scaling the rational values to integers for the tclique-weights */
   SCIP_Real*       pi;                      /* array of the dual solutions */
   SCIP_Bool        onlybest;                /* determines whether the maxvarsround variables with the best reduced costs should be added 
                                                (onlybest = true) or the first maxvarsround variables which are found are added (false) */
   SCIP_Bool        usegreedy;               /* determines whether a greedy method is used for finding variables with neg. reduced costs */
   SCIP_Bool        usetclique;              /* determines whether the tclique method is used for finding improving variables */
   int**            improvingstablesets;     /* array to store the maxvarsround stable sets with the most negative reduced costs */
   int*             nstablesetnodes;         /* array which stores the lengths of the stable sets in improvingstablesets */
   int              actindex;                /* the index at which the current stable set was inserted into improvingstablesets */
   SCIP_NODE*       bbnode;                  /* the current B&B-tree node, used for limiting the number of pricing rounds */
   int              noderounds;              /* the number of remaining pricing rounds at the current node */
   int              maxroundsroot;           /* maximum number of pricing rounds in the root, -1 for infinity, attention: positive value may lead to a non-optimal solution */
   int              maxroundsnode;           /* maximum number of pricing rounds in the B&B-nodes, -1 for infinity, attention: positive value may lead to a non-optimal solution */
   int              maxtcliquenodes;         /* maximum number of nodes used in the tclique algorithm for solving the stable set problem */
   SCIP_Real        lowerbound;              /* lower bound computed by the pricer */
};


/*
 * Local methods
 */

/** returns whether the graph has an uncolored node 
 */
static
SCIP_Bool hasUncoloredNode(
   TCLIQUE_GRAPH*        graph,              /**< the graph that should be colored */
   SCIP_Bool*            colored             /**< array of booleans, colored[i] == TRUE iff node i is colored */
   )
{
   int i;

   assert(graph != NULL);
   assert(colored != NULL);

   for ( i = 0; i < tcliqueGetNNodes(graph); i++)
   {
      /* node not yet colored */
      if (!colored[i])
      {
	return TRUE;
      }
   }
   return FALSE;
}

/** sorts the nodes from 0 to nnodes-1 w.r.t. the given weights */
static
SCIP_RETCODE sortNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< the weights for sorting */
   int                   nnodes,             /**< the number of nodes */
   int*                  sortednodes         /**< the array that will be overwritten with the sorted node numbers */
   )
{
   int i;
   SCIP_Real* values;

   assert(scip != NULL);
   assert(weights != NULL);

   /* create array with indices and copy the weights-array */
   SCIP_CALL( SCIPallocBufferArray(scip, &values, nnodes) ); 
   for ( i = 0; i < nnodes; i++)
   {
      sortednodes[i] = i;
      values[i] = weights[i];
   }

   /* sort the nodes w.r.t. the computed values */
   SCIPsortDownRealInt(values, sortednodes, nnodes);
   SCIPfreeBufferArray(scip, &values);

   return SCIP_OKAY;
}

/** computes a stable set with a greedy-method.  attention: the weight of the maximum stable set is not computed! */
static 
SCIP_RETCODE greedyStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< pointer to graph data structure */
   SCIP_Bool*            colored,            /**< array for marking yet colored nodes */
   int*                  maxstablesetnodes,  /**< pointer to store nodes of the maximum weight stableset */
   int*                  nmaxstablesetnodes  /**< pointer to store number of nodes in the maximum weight stableset */
   )
{
   SCIP_Bool indnode;
   int nnodes;
   int i;
   int j;
   int* degrees;
   int* sortednodes;
   SCIP_Real* values;    /* values for sorting the nodes: deg(v)+w(v)*nnodes  */

   assert(scip != NULL);
   assert(graph != NULL);
   assert(maxstablesetnodes != NULL);
   assert(nmaxstablesetnodes != NULL);

   /* get number of nodes */
   nnodes = tcliqueGetNNodes(graph);
   *nmaxstablesetnodes = 0;

   /* get the  degrees for the nodes in the graph */
   degrees = tcliqueGetDegrees(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &values, nnodes) );   
   SCIP_CALL( SCIPallocBufferArray(scip, &sortednodes, nnodes) );

   /* set values to the nodes which are used for sorting them */
   /* value = degree of the node + weight of the node * number of nodes, therefore the yet colored nodes
      (which have weight 0) have lower values than the not yet colored nodes which have weight 1 */
   for ( i = 0; i < nnodes; i++ )
   {
      sortednodes[i] = i;
      values[i] = ( colored[i] == TRUE ? degrees[i] : degrees[i]+nnodes );
   }

   /* sort the nodes w.r.t. the computed values */
   SCIPsortDownRealInt(values, sortednodes, nnodes);

   /* insert first node */
   maxstablesetnodes[0] = sortednodes[0];
   (*nmaxstablesetnodes) = 1;
   for ( i = 1; i < nnodes; i++)
   {
      /* check whether node is independent to nodes in the set */
      indnode = TRUE;
      for ( j = 0; j < (*nmaxstablesetnodes); j++ )
      {
         if ( tcliqueIsEdge(graph, sortednodes[i], maxstablesetnodes[j]) )
         {
            indnode = FALSE;
            break;
         }
      }
      if ( indnode == TRUE )
      {
         /* node is independent, thus add it to the set */
         maxstablesetnodes[*nmaxstablesetnodes] = sortednodes[i];
         (*nmaxstablesetnodes) = (*nmaxstablesetnodes)+1;
      }

   }
   SCIPfreeBufferArray(scip, &sortednodes);
   SCIPfreeBufferArray(scip, &values);   
   
   return SCIP_OKAY;
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = EPSFLOOR(sval, 0.0);
   upval = EPSCEIL(sval, 0.0);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** get integral number with error in the bounds which corresponds to given value scaled by a given scalar;
 *  should be used in connection with isIntegralScalar()  
 */
static 
SCIP_Longint getIntegralVal(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;
   SCIP_Longint intval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = EPSFLOOR(sval, 0.0);
   upval = EPSCEIL(sval, 0.0);

   if( SCIPrelDiff(sval, upval) >= mindelta )
      intval = (SCIP_Longint) upval;
   else
      intval = (SCIP_Longint) downval;
   
   return intval;
}



/** generates improving variables using a stable set found by the algorithm for maximum weight clique,
 *  decides whether to stop generating cliques with the algorithm for maximum weight clique
 */
static
TCLIQUE_NEWSOL(tcliqueNewsolPricer)
{
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(acceptsol != NULL);
   assert(stopsolving != NULL);

   pricerdata = (SCIP_PRICERDATA*)tcliquedata;

   assert(pricerdata != NULL);
   assert(pricerdata->scip != NULL);
   assert(pricerdata->nstablesetsfound >= 0);
   assert(pricerdata->scalefactor > 0);

   *acceptsol = FALSE;
   *stopsolving = FALSE;

   /* if the stable set was already created in a former pricing round, we don't have to add it a second time */
   if ( !COLORprobStableSetIsNew(pricerdata->scip, cliquenodes, ncliquenodes) )
      return;

   /* compute the index, at which the new stable set will be stored in the improvingstablesets-array */
   pricerdata->actindex = (pricerdata->actindex+1)%(pricerdata->maxvarsround);
   
   /* found maxvarsround variables */
   if ( pricerdata->nstablesetnodes[pricerdata->actindex] == -1 )
   {
      /* if we are looking for the best stable sets, continue at the beginning 
         and overwrite the stable set with least improvement */
      if ( pricerdata->onlybest )
      {
         pricerdata->actindex = 0;
      }
      /* if we only want to find maxvarsround variables (onlybest = false), we can stop now */
      else
      {
         *stopsolving = TRUE;
         return;         
      }
   }

   /* write the new improving stable set into the improvingstablesets-array */
   pricerdata->nstablesetnodes[pricerdata->actindex] = ncliquenodes;
   for ( i = 0; i < ncliquenodes; i++ )
   {
      pricerdata->improvingstablesets[pricerdata->actindex][i] = cliquenodes[i];
   }

   /* accept the solution as new incumbent */
   *acceptsol = TRUE;
   
}


/*
 * Callback methods of variable pricer
 */

/** copy method for pricer plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRICERCOPY(pricerCopyColoring)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(pricer != NULL);
   assert(strcmp(SCIPpricerGetName(pricer), PRICER_NAME) == 0);

   return SCIP_OKAY;
}


/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeColoring)
{ 
   SCIP_PRICERDATA* pricerdata;  

   assert(scip != NULL);
  
   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   /* free memory for pricerdata*/
   if ( pricerdata != NULL )
   {
      SCIPfreeMemory(scip, &pricerdata);
   }
   
   SCIPpricerSetData(pricer, NULL);
   return SCIP_OKAY;
}



/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolColoring)
{  
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   
   /* set maximal number of variables to be priced in each round */
   SCIPsetIntParam(scip, "pricers/coloring/maxvarsround", MAX(5,COLORprobGetNStableSets(scip))*MAX(50,COLORprobGetNNodes(scip))/50);

   pricerdata->bbnode = NULL;

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pi), COLORprobGetNNodes(scip)) );

   return SCIP_OKAY;
}



/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolColoring)
{  
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* free memory */
   for ( i = 0; i < pricerdata->maxvarsround; i++ )
   {
      SCIPfreeMemoryArray(scip, &(pricerdata->improvingstablesets[i]));
   }
   SCIPfreeMemoryArray(scip, &(pricerdata->improvingstablesets));
   SCIPfreeMemoryArray(scip, &(pricerdata->nstablesetnodes));
   SCIPfreeMemoryArray(scip, &(pricerdata->pi));

   return SCIP_OKAY;
}




/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostColoring)
{  
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */

   TCLIQUE_GRAPH*   graph;                 /* the current graph */
   TCLIQUE_GRAPH*   cgraph;                /* the complementary graph, used for tclique-algorithm */
   int              nnodes;                /* number of nodes in the graph */

   int*             sortednodes;           /* array of the nodes, sorted in specific way, atm by decreasing dual-solution*/
   SCIP_Real        maxstablesetweightreal;/* weigth of the maximal stable set computed by the greedy */
   SCIP_Bool        indnode;               /* boolean for greedy: is node independant? */

   int*             maxstablesetnodes;     /* pointer to store nodes of the maximum weight clique */
   int              nmaxstablesetnodes;    /* number of nodes in the maximum weight clique */
   TCLIQUE_WEIGHT   maxstablesetweight;    /* weight of the maximum weight clique */
   TCLIQUE_STATUS   status;                /* status of clique-computation */
   SCIP_Real        maxredcost;

   SCIP_VAR*        var;                   /* pointer to the new created variable */
   int              setnumber;             /* index of the new created variable */

   /* variables used for scaling the rational dual-solutions to integers */
   SCIP_Bool        scalesuccess;
   SCIP_Bool        weightsIntegral;

   int              i;
   int              j;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* count down number of remaining pricing rounds at the current node */
   if ( pricerdata->bbnode == SCIPgetCurrentNode(scip) )
   {
      if ( pricerdata->noderounds > 0 )
         pricerdata->noderounds--;
   }
   else
   {
      if ( pricerdata->bbnode == NULL )
      {
         pricerdata->noderounds = pricerdata->maxroundsroot;
         pricerdata->lowerbound = - SCIPinfinity(scip);
      }
      else
      {
         pricerdata->noderounds = pricerdata->maxroundsnode;
         pricerdata->lowerbound = - SCIPinfinity(scip);
      }
      pricerdata->bbnode = SCIPgetCurrentNode(scip);
   }
   /* stop pricing if limit for pricing rounds reached */
   if ( pricerdata->noderounds == 0 )
   {
      SCIPdebugMessage("maxrounds reached, pricing interrupted\n");

      /* set result and lowerbound pointer */
      *result = SCIP_DIDNOTRUN;
      *lowerbound = pricerdata->lowerbound;

      return SCIP_OKAY;
   }

   /* set result pointer */
   *result = SCIP_SUCCESS;

   /* get graph and number of nodes */
   graph = COLORconsGetCurrentGraph(scip);
   assert(graph != NULL);
   nnodes = tcliqueGetNNodes(graph);

   assert(SCIPgetNVars(scip) == COLORprobGetNStableSets(scip));

   /* get constraints */
   pricerdata->constraints = COLORprobGetConstraints(scip);

   /* get dual solutions and save them in pi */
   for ( i = 0; i < nnodes; i++)
   {
      pricerdata->pi[i] = SCIPgetDualsolSetppc(scip, pricerdata->constraints[i]);
   }
   pricerdata->nstablesetsfound = 0;
   /* ......greedy-heuristic........ */
   if ( pricerdata->usegreedy )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &sortednodes, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &maxstablesetnodes, nnodes) );
      SCIP_CALL( sortNodes(scip, pricerdata->pi, nnodes, sortednodes) );

      SCIPdebugMessage("starting greedy...\n");

      /* insert first node */
      maxstablesetnodes[0] = sortednodes[0];
      nmaxstablesetnodes = 1;
      maxstablesetweightreal = pricerdata->pi[sortednodes[0]];
      
      for ( i = 1; i < nnodes; i++ )
      {
         /* test if node is independant to nodes in stable set */
         indnode = TRUE;
         for ( j = 0; j < nmaxstablesetnodes; j++ )
         {
            if ( tcliqueIsEdge(graph, sortednodes[i], maxstablesetnodes[j]) )
            {
               indnode = FALSE;
               break;
            }
         }
         /* if node is independant to nodes in stable set, insert it into stable set*/
         if ( indnode )
         {
            maxstablesetnodes[nmaxstablesetnodes] = sortednodes[i];
            nmaxstablesetnodes = nmaxstablesetnodes+1;
            maxstablesetweightreal = maxstablesetweightreal + pricerdata->pi[sortednodes[i]];
         }
      }
      
      
      SCIPdebugMessage("value of the greedy-heuristik: %f \n", maxstablesetweightreal);
      setnumber = -1;
      if ( SCIPisFeasGT(scip, maxstablesetweightreal, 1.0) && COLORprobStableSetIsNew(scip, maxstablesetnodes, nmaxstablesetnodes) )
      {
         SCIP_CALL( COLORprobAddNewStableSet(scip, maxstablesetnodes, nmaxstablesetnodes, &setnumber) );
         
         assert(setnumber >= 0);
         pricerdata->nstablesetnodes[pricerdata->nstablesetsfound] = nmaxstablesetnodes;
         for ( i = 0; i < nmaxstablesetnodes; i++ )
         {
            pricerdata->improvingstablesets[pricerdata->nstablesetsfound][i] = maxstablesetnodes[i];
         }
         pricerdata->nstablesetsfound += 1;
         
         /* create variable for the stable set and add it to SCIP */
         SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0, 1, 1, SCIP_VARTYPE_BINARY, 
               TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setnumber) );

         COLORprobAddVarForStableSet(scip, setnumber, var);
         SCIPvarMarkDeletable(var);
         SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
         SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );
         
         /* add variable to the constraints in which it appears */
         for ( i = 0; i < nmaxstablesetnodes; i++ )
         {
            /* add variable to node constraints of nodes in the set */
            SCIP_CALL( SCIPaddCoefSetppc(scip, pricerdata->constraints[maxstablesetnodes[i]], var) );
         }
      }
      
      SCIPfreeBufferArray(scip, &sortednodes);
      SCIPfreeBufferArray(scip, &maxstablesetnodes);

      SCIPdebugMessage("%d vars created via greedy\n", pricerdata->nstablesetsfound);
   } 
   

   /* solve with tclique-algorithm */   
   /* only use tclique if the greedy found no improving stable set */
   if ( pricerdata->nstablesetsfound == 0 && pricerdata->usetclique )
   {
      SCIPdebugMessage("starting tclique algorithm...\n");
      maxredcost = 0;
      /* get the complementary graph from the current cons */
      cgraph = COLORconsGetComplementaryGraph(scip);
      SCIP_CALL( SCIPallocBufferArray(scip, &maxstablesetnodes, nnodes) );
      /* get dual solutions and set weight of nodes */
      weightsIntegral = TRUE;
      for ( i = 0; i < nnodes; i++ )
      {
         pricerdata->pi[i] = SCIPgetDualsolSetppc(scip, pricerdata->constraints[i]);
         
         if( !isIntegralScalar(pricerdata->pi[i], 1.0, -MINDELTA, MAXDELTA) )
         {
            weightsIntegral = FALSE;
         }
      }
      /* are weigths integral? */
      if( weightsIntegral )
      {
         pricerdata->scalefactor = 1.0;
         scalesuccess = TRUE;
      }
      else
      {
         /* compute factor, which makes the weights integral */
         scalesuccess = FALSE;
         SCIP_CALL( SCIPcalcIntegralScalar(pricerdata->pi, nnodes, -MINDELTA, MAXDELTA, MAXDNOM, MAXSCALE, 
               &(pricerdata->scalefactor), &scalesuccess) );
      }
      assert(scalesuccess);
      /* change the weights for the nodes in the graph to the dual solution value * scalefactor */
      for ( i = 0; i < nnodes; i++ )
      {
         tcliqueChangeWeight(cgraph, i, getIntegralVal(pricerdata->pi[i], pricerdata->scalefactor, -MINDELTA, MAXDELTA));
      }
      /* clear the improvingstablesets array */
      pricerdata->actindex = -1;
      for ( i = 0; i < pricerdata->maxvarsround-pricerdata->nstablesetsfound; i++ )
      {
         pricerdata->nstablesetnodes[i] = 0;
      }
      for ( ; i < pricerdata->maxvarsround; i++ )
      {
         pricerdata->nstablesetnodes[i] = -1;
      }

      /* compute maximal clique */
      tcliqueMaxClique(NULL, NULL, NULL, NULL, cgraph, tcliqueNewsolPricer, (TCLIQUE_DATA*)pricerdata, maxstablesetnodes,
         &(nmaxstablesetnodes), &maxstablesetweight, 0,
         (int)getIntegralVal(pricerdata->scalefactor, 1.0, -MINDELTA, MAXDELTA), pricerdata->maxtcliquenodes, 0, INT_MAX, -1,
         NULL, &status);
      assert(status == TCLIQUE_OPTIMAL);

      /* if only the best vaiable should be priced per round, take the one which is given as return value from
         tcliqueMaxClique and put it into improvingstablesets array so that it will be inserted into the LP */
      if ( pricerdata->onlybest && pricerdata->maxvarsround == 1 )
      {
         pricerdata->nstablesetnodes[0] = nmaxstablesetnodes;
         for ( i = 0; i < nmaxstablesetnodes; i++ )
         {
            pricerdata->improvingstablesets[0][i] = maxstablesetnodes[i];
         }
      }
      
      SCIPfreeBufferArray(scip, &maxstablesetnodes);

      /* insert all variables in the array improvingstablesets into the LP */
      for ( i = 0; i < pricerdata->maxvarsround; i++ )
      { 
         if ( pricerdata->nstablesetnodes[i] > 0 )
         {
            maxstablesetweightreal = 0;
            for ( j = 0; j < pricerdata->nstablesetnodes[i]; j++ )
            {
               maxstablesetweightreal += pricerdata->pi[pricerdata->improvingstablesets[i][j]];
            }
            if ( maxredcost < maxstablesetweightreal )
            {
               maxredcost = maxstablesetweightreal;
            }
            if ( SCIPisFeasGT(scip, maxstablesetweightreal, 1.0) )
            {
               setnumber = -1;
               /* insert new variable */
               SCIP_CALL( COLORprobAddNewStableSet(pricerdata->scip, pricerdata->improvingstablesets[i], 
                     pricerdata->nstablesetnodes[i], &setnumber) );
               /* only insert, if there yet is no variable for this stable set */
               if ( setnumber >= 0  )
               {
                  /* create variable for the stable set and add it to SCIP */
                  SCIP_CALL( SCIPcreateVar(pricerdata->scip, &var, NULL, 0, 1, 1, SCIP_VARTYPE_BINARY, 
                        TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setnumber) );

                  COLORprobAddVarForStableSet(pricerdata->scip, setnumber, var);
                  SCIPvarMarkDeletable(var);
                  SCIP_CALL( SCIPaddPricedVar(pricerdata->scip, var, 1.0) );
                  SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

                  pricerdata->nstablesetsfound += 1;
                  /* add variable to the constraints in which it appears */
                  for ( j = 0; j < pricerdata->nstablesetnodes[i]; j++ )
                  {
                     /* add variable to node constraints of nodes in the set */
                     SCIP_CALL( SCIPaddCoefSetppc(pricerdata->scip, 
                           pricerdata->constraints[pricerdata->improvingstablesets[i][j]], var) );
                  }
               }
            }
         }
      }
      if ( SCIPisFeasGT(scip, maxredcost, 1.0) )
      {
         if ( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            pricerdata->lowerbound = MAX( pricerdata->lowerbound, (SCIPgetLPObjval(scip) + ((1.0 - maxredcost) * SCIPgetPrimalbound(scip))) );
         }
      }
   }

   
   return SCIP_OKAY;
}


/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasColoring)
{  
   TCLIQUE_GRAPH* graph;
   int nnodes;                  /* number of nodes */
   int* maxstablesetnodes;      /* array containig the nodes of the max stable set */
   int nmaxstablesetnodes;      /* number of nodes in stable set */
   int setnumber;               /* number of already found stable sets */
   SCIP_VAR* var;               /* var for the actual stable set */
   SCIP_CONS** constraints;     /* array of added constraints */
   SCIP_Bool* colored;          /* array for marking of yet colored nodes
                                   colored_i = true iff node i is already colored */
   int**              stablesets;
   int*               nstablesetelements;
   int                nstablesets;
   int i;
   int j;
   assert(scip != NULL);

   graph = COLORconsGetCurrentGraph(scip);
   assert(graph != NULL);

   nnodes = COLORprobGetNNodes(scip);
   assert(nnodes > 0);

   /* get the node-constraits */
   constraints = COLORprobGetConstraints(scip);
   assert(constraints != NULL);

   /* get all yet computed stable sets */
   COLORprobGetStableSets(scip, &stablesets, &nstablesetelements, &nstablesets);
   assert(stablesets != NULL && nstablesetelements != NULL);
   assert(nstablesets >= 0);
   assert(nnodes == tcliqueGetNNodes(graph));

   /* allocate memory for arrays */
   SCIP_CALL( SCIPallocBufferArray( scip, &colored, nnodes) );
   SCIP_CALL( SCIPallocBufferArray( scip, &maxstablesetnodes, nnodes) );
   nmaxstablesetnodes = 0;

   /* fill colored-array with FALSE */
   BMSclearMemoryArray(colored, nnodes);

   /* go through all stable sets and set colored to true for nodes in them */
   for ( i = 0; i < nstablesets; i++ )
   {
      if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal( COLORprobGetVarForStableSet(scip, i))) 
         && (SCIPgetNNodes(scip) == 0 || SCIPvarIsInLP(COLORprobGetVarForStableSet(scip, i))
            || SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ) )
      {
         for ( j = 0; j < nstablesetelements[i]; j++ )
         {
            colored[stablesets[i][j]] = TRUE;
         }
      }
   }

   /* create maximal Stable Sets until all Nodes are covered */
   while ( hasUncoloredNode(graph, colored) )
   {
      greedyStableSet(scip, graph, colored, maxstablesetnodes, &nmaxstablesetnodes);
      SCIPsortDownInt(maxstablesetnodes, nmaxstablesetnodes);
      SCIP_CALL( COLORprobAddNewStableSet(scip, maxstablesetnodes, nmaxstablesetnodes, &setnumber) );
      assert(setnumber != -1);
      
      /* create variable for the stable set and add it to SCIP*/
      SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0, 1, 1, SCIP_VARTYPE_BINARY, 
            TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*) (size_t) setnumber) );
      COLORprobAddVarForStableSet(scip, setnumber, var);
      SCIPvarMarkDeletable(var);
      SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
      SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

      for ( i = 0; i < nmaxstablesetnodes; i++ )
      {
         /* add variable to node constraints of nodes in the set */
         SCIP_CALL( SCIPaddCoefSetppc(scip, constraints[maxstablesetnodes[i]], var) );
         /* mark node as colored */
         colored[maxstablesetnodes[i]] = TRUE;
      }
   }
   /* free memory */
   SCIPfreeBufferArray(scip, &maxstablesetnodes);
   SCIPfreeBufferArray(scip, &colored);
   return SCIP_OKAY;

}

/** method to call, when the maximal number of variables priced in each round is changed */
static
SCIP_DECL_PARAMCHGD(paramChgdMaxvarsround)
{  
   SCIP_PARAMDATA* paramdata;
   SCIP_PRICERDATA* pricerdata;
   int i;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);
   pricerdata = (SCIP_PRICERDATA*) paramdata;

   if( pricerdata->maxvarsround == pricerdata->oldmaxvarsround )
      return SCIP_OKAY;
   
   if ( pricerdata->maxvarsround <= 1 ) 
      pricerdata->maxvarsround = 2;

   if ( pricerdata->maxvarsround == pricerdata->oldmaxvarsround && pricerdata->nstablesetnodes != NULL )
      return SCIP_OKAY;

   /* free old memory */
   if ( pricerdata -> oldmaxvarsround > 0 )
   {
      /* free memory */
      for ( i = 0; i < pricerdata->oldmaxvarsround; i++ )
      {
         SCIPfreeMemoryArray(scip, &(pricerdata->improvingstablesets[i]));
      }
      SCIPfreeMemoryArray(scip, &(pricerdata->improvingstablesets));
      SCIPfreeMemoryArray(scip, &(pricerdata->nstablesetnodes));
   }
   
   /* allocate memory of the new size */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nstablesetnodes), pricerdata->maxvarsround) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->improvingstablesets), pricerdata->maxvarsround) );
   for ( i = 0; i < pricerdata->maxvarsround; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->improvingstablesets[i]), COLORprobGetNNodes(scip)) );
   }

   SCIPdebugMessage("maxvarsround changed from %d to %d\n", pricerdata->oldmaxvarsround, pricerdata->maxvarsround);

   pricerdata->oldmaxvarsround = pricerdata->maxvarsround;

   return SCIP_OKAY;
}


/*
 * variable pricer specific interface methods
 */

/** creates the coloring variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerColoring(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
   pricerdata->scip = scip;

   pricerdata->maxvarsround = 0;
   pricerdata->oldmaxvarsround = 0;


   pricer = NULL;
   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostColoring, pricerFarkasColoring, pricerdata) );
   assert(pricer != NULL);
   
   /* include non-fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPricerCopy(scip, pricer, pricerCopyColoring) );
   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeColoring) );
   SCIP_CALL( SCIPsetPricerInitsol(scip, pricer, pricerInitsolColoring) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolColoring) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxvarsround",
         "maximum number of variables that the coloring variable pricer creates each round",
         &pricerdata->maxvarsround, TRUE, DEFAULT_MAXVARSROUND, -1, INT_MAX, paramChgdMaxvarsround, (SCIP_PARAMDATA*)pricerdata) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "pricers/coloring/usetclique",
         "should the tclique-algorithm be used to solve the pricing-problem to optimality?\n \
             WARNING: computed (optimal) solutions are not necessarily optimal if this is set to FALSE",
         &pricerdata->usetclique, TRUE, DEFAULT_USETCLIQUE, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip,
         "pricers/coloring/usegreedy",
         "should a greedy method be used to compute improving stable sets before potential use of tclique",
         &pricerdata->usegreedy, FALSE, DEFAULT_USEGREEDY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "pricers/coloring/onlybest",
         "should the best variables be addded to the problem instead of adding the first found variables?",
         &pricerdata->onlybest, FALSE, DEFAULT_ONLYBEST, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxroundsroot",
         "maximum number of pricing rounds in the root node (-1: no limit)",
         &pricerdata->maxroundsroot, TRUE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxroundsnode",
         "maximum number of pricing rounds in each node (except root node)(-1: no limit)",
         &pricerdata->maxroundsnode, TRUE, DEFAULT_MAXROUNDSNODE, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxtcliquenodes",
         "maximum number of B&B-nodes used in the tclique-algorithm",
         &pricerdata->maxtcliquenodes, TRUE, DEFAULT_MAXTCLIQUENODES, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
