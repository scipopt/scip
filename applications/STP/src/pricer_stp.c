
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

/**@file   pricer_stp.c
 * @brief  stp variable pricer
 * @author Daniel Rehfeldt
 */

#include <assert.h>
#include <string.h>
#include "scip/cons_linear.h"
#include <stdio.h>
#include <stdlib.h>
#include "pricer_stp.h"
#include "grph.h"


#define PRICER_NAME            "stp"
#define PRICER_DESC            "pricer for stp"
#define PRICER_PRIORITY        1
#define PRICER_DELAY           TRUE    /* only call pricer if all problem variables have non-negative reduced costs */

/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_CONS**           edgecons;           /**< array containing all edge constraints */
   SCIP_CONS**           pathcons;           /**< array containing all path constraints */
   SCIP_Real*            pi;                 /**< array of the dual solutions to the edge constraints*/
   SCIP_Real*            mi;	             /**< array of the dual solutions to the path constraint*/
   SCIP_Real             lowerbound;         /**< lower bound computed by the pricer */
   int* 	         realterms;          /**< array of all terminals except the root */
   int*                  ncreatedvars;       /**< array (1,...,realnterms); in ncreatedvars[t] the number of
                                                created variables (during pricing) pertaining to terminal t is stored */
   int                   root;               /**< root node */
   int                   nedges;             /**< number of edges */
   int                   realnterms;         /**< number of terminals except the root */
   SCIP_Bool             bigt;               /**< stores whether the 'T' model is being used */
};


/*
 * Local methods
 */


/*
 * Callback methods of variable pricer
 */

/** copy method for pricer plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRICERCOPY(pricerCopyStp)
{ /*lint --e{715}*/
   SCIPdebugPrintf("pricerCopy \n");
   assert(pricer != NULL);
   assert(strcmp(SCIPpricerGetName(pricer), PRICER_NAME) == 0);

   return SCIP_OKAY;
}

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
/**! [SnippetPricerFreeSTP] */
static
SCIP_DECL_PRICERFREE(pricerFreeStp)
{
   SCIP_PRICERDATA* pricerdata;
   SCIPdebugPrintf("pricerfree \n");
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
/**! [SnippetPricerFreeSTP] */

/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitStp)
{
   SCIP_PRICERDATA* pricerdata;
   SCIPdebugPrintf("pricer Init \n");
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   pricerdata->edgecons = SCIPprobdataGetEdgeConstraints(scip);
   pricerdata->pathcons = SCIPprobdataGetPathConstraints(scip);
   pricerdata->root = SCIPprobdataGetRoot(scip);
   pricerdata->nedges = SCIPprobdataGetNEdges(scip);
   pricerdata->realnterms = SCIPprobdataGetRNTerms(scip);
   pricerdata->realterms = SCIPprobdataGetRTerms(scip);
   pricerdata->bigt = SCIPprobdataIsBigt(scip);
   pricerdata->lowerbound = 0;
   return SCIP_OKAY;
}

/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolStp)
{
   SCIP_PRICERDATA* pricerdata;
   SCIPdebugPrintf("pricerinitsol \n");
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* allocate memory */
   if( !pricerdata->bigt )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pi), SCIPprobdataGetNEdges(scip) * SCIPprobdataGetRNTerms(scip)) );
   }
   else
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pi), SCIPprobdataGetNEdges(scip)) );
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->mi), SCIPprobdataGetRNTerms(scip)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pricerdata->ncreatedvars, SCIPprobdataGetRNTerms(scip)) );
   BMSclearMemoryArray(pricerdata->ncreatedvars, SCIPprobdataGetRNTerms(scip));

   return SCIP_OKAY;
}

/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolStp)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* free memory */
   SCIPfreeMemoryArray(scip, &(pricerdata->mi));
   SCIPfreeMemoryArray(scip, &(pricerdata->pi));
   SCIPfreeMemoryArray(scip, &pricerdata->ncreatedvars);

   return SCIP_OKAY;
}

/** method for either Farkas or Redcost pricing */
static
SCIP_RETCODE pricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_Real*            lowerbound,         /**< lowerbound pointer */
   SCIP_Bool             farkas              /**< TRUE: Farkas pricing; FALSE: Redcost pricing */
   )
{
   SCIP_PRICERDATA* pricerdata; /* the data of the pricer */
   SCIP_PROBDATA* probdata;
   GRAPH* graph;
   SCIP_VAR* var;
   PATH* path;
   SCIP_Real* edgecosts;  /* edgecosts of the current subproblem */
   char varname[SCIP_MAXSTRLEN];
   SCIP_Real newlowerbound = -SCIPinfinity(scip);
   SCIP_Real redcost;   /* reduced cost */
   int tail;
   int e;
   int t;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIPdebugMessage("solstat=%d\n", SCIPgetLPSolstat(scip));

   if( !farkas && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      newlowerbound = SCIPgetSolTransObj(scip, NULL);

   SCIPdebug( SCIP_CALL( SCIPprintSol(scip, NULL, NULL, FALSE) ) );

# if 0
   if ( pricerdata->lowerbound <= 4 )
   {
      char label[SCIP_MAXSTRLEN];
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "X%g.gml", pricerdata->lowerbound);
      SCIP_CALL( SCIPprobdataPrintGraph(scip, label , NULL, TRUE) );
      pricerdata->lowerbound++;
   }
#endif
   /* get the graph*/
   graph = SCIPprobdataGetGraph(probdata);

   /* get dual solutions and save them in mi and pi */
   for( t = 0; t < pricerdata->realnterms; ++t )
   {
      if( farkas )
      {
	 pricerdata->mi[t] = SCIPgetDualfarkasLinear(scip, pricerdata->pathcons[t]);
      }
      else
      {
         pricerdata->mi[t] = SCIPgetDualsolLinear(scip, pricerdata->pathcons[t]);
         assert(!SCIPisNegative(scip, pricerdata->mi[t]));
      }
   }

   for( e = 0; e < pricerdata->nedges; ++e )
   {
      if( !pricerdata->bigt )
      {
         for( t = 0; t < pricerdata->realnterms; ++t )
         {
            if( farkas )
	    {
               pricerdata->pi[t * pricerdata->nedges + e] = SCIPgetDualfarkasLinear(
                  scip, pricerdata->edgecons[t * pricerdata->nedges + e]);
	    }
            else
	    {
               pricerdata->pi[t * pricerdata->nedges + e] = SCIPgetDualsolLinear(
                  scip, pricerdata->edgecons[t * pricerdata->nedges + e]);
	    }
         }
      }
      else
      {
         if( farkas )
	 {
	    pricerdata->pi[e] = SCIPgetDualfarkasLinear(
               scip, pricerdata->edgecons[e]);
	 }
	 else
	 {
	    pricerdata->pi[e] = SCIPgetDualsolLinear(
               scip, pricerdata->edgecons[e]);
	 }
      }
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &path, graph->knots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgecosts, pricerdata->nedges) );

   if( pricerdata->bigt )
   {
      for( e = 0; e < pricerdata->nedges; ++e )
      {
         edgecosts[e] = (-pricerdata->pi[e]);
      }
   }
   /* find shortest r-t (r root, t terminal) paths and create corresponding variables iff reduced cost < 0 */
   for( t = 0; t < pricerdata->realnterms; ++t )
   {
      for( e = 0; e < pricerdata->nedges; ++e )
      {
	 if( !pricerdata->bigt )
	 {
            edgecosts[e] = (-pricerdata->pi[t * pricerdata->nedges + e]);
	 }

         assert(!SCIPisNegative(scip, edgecosts[e]));
      }

      for( i = 0; i < graph->knots; i++ )
         graph->mark[i] = 1;

      graph_path_exec(scip, graph, FSP_MODE, pricerdata->root, edgecosts, path);

      /* compute reduced cost of shortest path to terminal t */
      redcost = 0.0;
      tail = pricerdata->realterms[t];
      while( tail != pricerdata->root )
      {
         redcost += edgecosts[path[tail].edge];
	 tail = graph->tail[path[tail].edge];
      }
      redcost -= pricerdata->mi[t];

      if( !farkas && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         newlowerbound += redcost;
      }
      /* check if reduced cost < 0 */
      if( SCIPisNegative(scip, redcost) )
      {
	 /* create variable to the shortest path (having reduced cost < 0) */
         var = NULL;
	 sprintf(varname, "PathVar%d_%d", t, pricerdata->ncreatedvars[t]);
         ++(pricerdata->ncreatedvars[t]);

         SCIP_CALL( SCIPcreateVarBasic(scip, &var, varname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddPricedVar(scip, var, -redcost) );
         tail = pricerdata->realterms[t];
         while( tail != pricerdata->root )
         {
            /* add variable to constraints */
	    if( !pricerdata->bigt )
	    {
	       SCIP_CALL( SCIPaddCoefLinear(scip, pricerdata->edgecons[t * pricerdata->nedges + path[tail].edge], var, 1.0) );
	    }
	    else
	    {
	       SCIP_CALL( SCIPaddCoefLinear(scip, pricerdata->edgecons[path[tail].edge], var, 1.0) );
	    }

	    tail = graph->tail[path[tail].edge];
         }
         SCIP_CALL( SCIPaddCoefLinear(scip, pricerdata->pathcons[t], var, 1.0) );
      }
   }

   if( !farkas && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      *lowerbound = newlowerbound;

   SCIPfreeMemoryArray(scip, &edgecosts);
   SCIPfreeMemoryArray(scip, &path);

   return SCIP_OKAY;
}

/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostStp)
{ /*lint --e{715}*/
   SCIP_CALL( pricing(scip, pricer, lowerbound, FALSE) );

   /* set result pointer */
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** Farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasStp)
{ /*lint --e{715}*/
   SCIP_CALL( pricing(scip, pricer,  NULL, TRUE) );

   return SCIP_OKAY;
}

/*
 * variable pricer specific interface methods
 */

/** creates the stp variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{

   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;
   SCIPdebugPrintf("include Pricer \n");

   /* create stp variable pricer data */
   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
   pricerdata->scip = scip;

   /* include variable pricer */
   pricer = NULL;
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostStp, pricerFarkasStp, pricerdata) );
   assert(pricer != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPricerCopy(scip, pricer, pricerCopyStp) );
   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeStp) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitStp) );
   SCIP_CALL( SCIPsetPricerInitsol(scip, pricer, pricerInitsolStp) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolStp) );

   return SCIP_OKAY;
}
