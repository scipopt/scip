/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_mcf.c,v 1.25 2008/06/18 18:35:13 bzfpfend Exp $"

/*#define SCIP_DEBUG*/
/**@file   sepa_mcf.c
 * @brief  multi-commodity-flow network cut separator
 * @author Tobias Achterberg
 *
 * We try to identify a multi-commodity flow structure in the LP relaxation of the
 * following type:
 *
 *  (1)  \sum_{a \in \delta^+(v)} f_a^k  - \sum_{a \in \delta^-(v)} f_a^k  <=  -d_v^k   for all v \in V and k \in K
 *  (2)  \sum_{k \in K} f_a^k - c_a x_a                                    <=  0        for all a \in A
 *
 * Constraints (1) are flow conservation constraints, which say that for each commodity k and node v the
 * outflow (\delta^+(v)) minus the inflow (\delta^-(v)) of a node v must not exceed the negative of the demand of
 * node v in commodity k. To say it the other way around, inflow minus outflow must be at least equal to the demand.
 * Constraints (2) are the arc capacity constraints, which say that the sum of all flow over an arc a must not
 * exceed its capacity c_a x_a, with x being a binary or integer variable.
 * c_a x_a does not need to be a single product of a capacity and an integer variable; we also accept general scalar
 * products.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_mcf.h"


#define SEPA_NAME              "mcf"
#define SEPA_DESC              "multi-commodity-flow network cut separator"
#define SEPA_PRIORITY            -10000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_NCLUSTERS             8 /**< number of clusters to generate in the shrunken network */
#define DEFAULT_MAXROWFAC          1e+4 /**< maximal row aggregation factor */
#define DEFAULT_MAXTESTDELTA         -1 /**< maximal number of different deltas to try (-1: unlimited) */
#define DEFAULT_FIXINTEGRALRHS     TRUE /**< should an additional variable be complemented if f0 = 0? */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MODELTYPE             0 /**< model type of network (0: auto, 1:directed, 2:undirected) */
#define DEFAULT_FLOWTYPE              0 /**< type of flow variables in the network (0: auto, 1:continuous/integer, 2:binary) */

#define BOUNDSWITCH                 0.5
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE
#define MINFRAC                    0.05
#define MAXFRAC                    0.95

#define MINNODES                      2 /**< minimal number of nodes in network to keep it for separation */
#define MINARCS                       1 /**< minimal number of arcs in network to keep it for separation */




/*
 * Data structures
 */

/** model type of the network */
enum SCIP_McfModeltype
{
   SCIP_MCFMODELTYPE_AUTO       = 0,         /**< model type should be detected automatically */
   SCIP_MCFMODELTYPE_DIRECTED   = 1,         /**< directed network where each arc has its own capacity */
   SCIP_MCFMODELTYPE_UNDIRECTED = 2          /**< directed network where anti-parallel arcs share the capacity */
};
typedef enum SCIP_McfModeltype SCIP_MCFMODELTYPE;

/** type of flow variables in the model */
enum SCIP_McfFlowtype
{
   SCIP_MCFFLOWTYPE_AUTO       = 0,          /**< flow variable type should be detected automatically */
   SCIP_MCFFLOWTYPE_CONTINUOUS = 1,          /**< continuous or integer flow variables: coefficients in capacity constraints are +1 */
   SCIP_MCFFLOWTYPE_BINARY     = 2           /**< binary flow variables: coefficients in capacity are equal to the demand */
};
typedef enum SCIP_McfFlowtype SCIP_MCFFLOWTYPE;

/** extracted multi-commodity-flow network */
struct SCIP_McfNetwork
{
   SCIP_ROW***           nodeflowrows;       /**< nodeflowrows[v][k]: flow conservation constraint for node v and
                                              *   commodity k; NULL if this node does not exist in the commodity */
   SCIP_Real**           nodeflowscales;     /**< scaling factors to convert nodeflowrows[v][k] into a +/-1 >= row */
   SCIP_ROW**            arccapacityrows;    /**< arccapacity[a]: capacity constraint on arc a;
                                              *   NULL if uncapacitated */
   SCIP_Real*            arccapacityscales;  /**< scaling factors to convert arccapacity[a] into a <= row with
                                              *   positive entries for the flow variables */
   int*                  arcsources;         /**< source node ids of arcs */
   int*                  arctargets;         /**< target node ids of arcs */
   int*                  colcommodity;       /**< commodity number of each column, or -1 */
   int                   nnodes;             /**< number of nodes in the graph */
   int                   narcs;              /**< number of arcs in the graph */
   int                   ncommodities;       /**< number of commodities */
   SCIP_MCFMODELTYPE     modeltype;          /**< detected model type of the network */
   SCIP_MCFFLOWTYPE      flowtype;           /**< detected flow variable type of the network */
};
typedef struct SCIP_McfNetwork SCIP_MCFNETWORK;

/** separator data */
struct SCIP_SepaData
{
   SCIP_MCFNETWORK**     mcfnetworks;        /**< array of multi-commodity-flow network structures */
   int                   nmcfnetworks;       /**< number of multi-commodity-flow networks (-1: extraction not yet done) */
   SCIP_Real             maxrowfac;          /**< maximal row aggregation factor */
   int                   nclusters;          /**< number of clusters to generate in the shrunken network */
   int                   maxtestdelta;	     /**< maximal number of different deltas to try (-1: unlimited) */
   int                   modeltype;          /**< model type of the network */
   int                   flowtype;           /**< type of flow variables in the network */
   SCIP_Bool             fixintegralrhs;     /**< should an additional variable be complemented if f0 = 0? */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};

/** internal MCF extraction data to pass to subroutines */
struct mcfdata
{
   unsigned char*        flowrowsigns;       /**< potential or actual sides of rows to be used as flow conservation constraint */
   SCIP_Real*            flowrowscalars;     /**< scalar of rows to transform into +/-1 coefficients */
   SCIP_Real*            flowrowscores;      /**< score value indicating how sure we are that this is indeed a flow conservation constraint */
   unsigned char*        capacityrowsigns;   /**< potential or actual sides of rows to be used as capacity constraint */
   SCIP_Real*            capacityrowscores;  /**< score value indicating how sure we are that this is indeed a capacity constraint */
   int*                  flowcands;          /**< list of row indices that are candidates for flow conservation constraints */
   int                   nflowcands;         /**< number of elements in flow candidate list */
   int*                  capacitycands;      /**< list of row indices that are candidates for capacity constraints */
   int                   ncapacitycands;     /**< number of elements in capacity candidate list */
   SCIP_Bool*            plusflow;           /**< is column c member of a flow row with coefficient +1? */
   SCIP_Bool*            minusflow;          /**< is column c member of a flow row with coefficient -1? */
   int                   ncommodities;       /**< number of commodities */
   int*                  commoditysigns;     /**< +1: regular, -1: all arcs have opposite direction; 0: undecided */
   int                   commoditysignssize; /**< size of commoditysigns array */
   int*                  colcommodity;       /**< commodity number of each column, or -1 */
   int*                  rowcommodity;       /**< commodity number of each row, or -1 */
   int*                  colarcid;           /**< arc id of each flow column, or -1 */
   int*                  rowarcid;           /**< arc id of each capacity row, or -1 */
   int*                  rownodeid;          /**< node id of each flow conservation row, or -1 */
   int*                  arcsources;         /**< source node ids of arcs */
   int*                  arctargets;         /**< target node ids of arcs */
   int*                  firstoutarcs;       /**< for each node the first arc id for which the node is the source node */
   int*                  firstinarcs;        /**< for each node the first arc id for which the node is the target node */
   int*                  nextoutarcs;        /**< for each arc the next outgoing arc in the adjacenty list */
   int*                  nextinarcs;         /**< for each arc the next outgoing arc in the adjacenty list */
   int*                  newcols;            /**< columns of current commodity that have to be inspected for incident flow conservation rows */
   int                   nnewcols;           /**< number of newcols */
   int                   narcs;              /**< number of arcs in the extracted graph */
   int                   nnodes;             /**< number of nodes in the extracted graph */
   int                   ninconsistencies;   /**< number of inconsistencies between the commodity graphs */
   SCIP_ROW**            capacityrows;       /**< capacity row for each arc */
   int                   capacityrowssize;   /**< size of array */
   SCIP_Bool*            colisincident;      /**< temporary memory for column collection */
   int*                  zeroarcarray;       /**< temporary array of zeros */
   SCIP_MCFMODELTYPE     modeltype;          /**< model type that is used for this network extraction */
   SCIP_MCFFLOWTYPE      flowtype;           /**< flow variable type that is used for this network extraction */
};
typedef struct mcfdata MCFDATA;              /**< internal MCF extraction data to pass to subroutines */

/** data structure to put on the arc heap */
struct arcentry
{
   int                   arcid;              /**< index of the arc */
   SCIP_Real             weight;             /**< weight of the arc in the separation problem */
};
typedef struct arcentry ARCENTRY;

/** arc priority queue */
struct arcqueue
{
   SCIP_PQUEUE*          pqueue;             /**< priority queue of elements */
   ARCENTRY*             arcentries;         /**< elements on the heap */
};
typedef struct arcqueue ARCQUEUE;

/** partitioning of the nodes into clusters */
struct nodepartition
{
   int*                  representatives;    /**< mapping of node ids to their representatives within their cluster */
   int*                  nodeclusters;       /**< cluster for each node id */
   int*                  clusternodes;       /**< node ids sorted by cluster */
   int*                  clusterbegin;       /**< first entry in clusternodes for each cluster (size: nclusters+1) */
   int                   nclusters;          /**< number of clusters */
};
typedef struct nodepartition NODEPARTITION;



/*
 * Local methods
 */

#define LHSPOSSIBLE     1                    /**< we may use the constraint as lhs <= a*x */
#define RHSPOSSIBLE     2                    /**< we may use the constraint as a*x <= rhs */
#define LHSASSIGNED     4                    /**< we have chosen to use the constraint as lhs <= a*x */
#define RHSASSIGNED     8                    /**< we have chosen to use the constraint as a*x <= rhs */
#define DISCARDED      16                    /**< we have chosen to not use the constraint */
#define UNDIRECTED     32                    /**< the capacity candidate has two flow variables for a commodity */


/** creates an empty MCF network data structure */
static
SCIP_RETCODE mcfnetworkCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK**     mcfnetwork          /**< MCF network structure */
   )
{
   assert(mcfnetwork != NULL);

   SCIP_CALL( SCIPallocMemory(scip, mcfnetwork) );
   (*mcfnetwork)->nodeflowrows = NULL;
   (*mcfnetwork)->nodeflowscales = NULL;
   (*mcfnetwork)->arccapacityrows = NULL;
   (*mcfnetwork)->arccapacityscales = NULL;
   (*mcfnetwork)->arcsources = NULL;
   (*mcfnetwork)->arctargets = NULL;
   (*mcfnetwork)->colcommodity = NULL;
   (*mcfnetwork)->nnodes = 0;
   (*mcfnetwork)->narcs = 0;
   (*mcfnetwork)->ncommodities = 0;

   return SCIP_OKAY;
}

/** frees MCF network data structure */
static
SCIP_RETCODE mcfnetworkFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK**     mcfnetwork          /**< MCF network structure */
   )
{
   assert(mcfnetwork != NULL);

   if( *mcfnetwork != NULL )
   {
      int v;
      int a;

      for( v = 0; v < (*mcfnetwork)->nnodes; v++ )
      {
         int k;

         for( k = 0; k < (*mcfnetwork)->ncommodities; k++ )
         {
            if( (*mcfnetwork)->nodeflowrows[v][k] != NULL )
            {
               SCIP_CALL( SCIPreleaseRow(scip, &(*mcfnetwork)->nodeflowrows[v][k]) );
            }
         }
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowrows[v]);
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowscales[v]);
      }
      for( a = 0; a < (*mcfnetwork)->narcs; a++ )
      {
         if( (*mcfnetwork)->arccapacityrows[a] != NULL )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &(*mcfnetwork)->arccapacityrows[a]) );
         }
      }
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowrows);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowscales);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arccapacityrows);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arccapacityscales);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arcsources);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arctargets);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->colcommodity);

      SCIPfreeMemory(scip, mcfnetwork);
   }

   return SCIP_OKAY;
}

/** fills the MCF network structure with the MCF data */
static
SCIP_RETCODE mcfnetworkFill(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK*      mcfnetwork,         /**< MCF network structure */
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   int*                  compnodeid,         /**< temporary storage for v -> compv mapping; must be set to -1 for all v */
   int*                  compnodes,          /**< array of node ids of the component */
   int                   ncompnodes,         /**< number of nodes in the component */
   int*                  comparcs,           /**< array of arc ids of the component */
   int                   ncomparcs           /**< number of arcs in the component */
   )
{
   unsigned char*    flowrowsigns     = mcfdata->flowrowsigns;
   SCIP_Real*        flowrowscalars   = mcfdata->flowrowscalars;
   unsigned char*    capacityrowsigns = mcfdata->capacityrowsigns;
   int*              flowcands        = mcfdata->flowcands;
   int               nflowcands       = mcfdata->nflowcands;
   int               ncommodities     = mcfdata->ncommodities;
   int*              commoditysigns   = mcfdata->commoditysigns;
   int*              colcommodity     = mcfdata->colcommodity;
   int*              rowcommodity     = mcfdata->rowcommodity;
   int*              rownodeid        = mcfdata->rownodeid;
   SCIP_ROW**        capacityrows     = mcfdata->capacityrows;
   SCIP_MCFMODELTYPE modeltype        = mcfdata->modeltype;
   SCIP_MCFFLOWTYPE  flowtype         = mcfdata->flowtype;

   SCIP_Real* comdemands;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;
   int* compcommodity;
   int ncompcommodities;
   int v;
   int a;
   int k;
   int i;
   int c;

   assert(mcfnetwork != NULL);
   assert(modeltype != SCIP_MCFMODELTYPE_AUTO);
   assert(flowtype != SCIP_MCFFLOWTYPE_AUTO);
   assert(2 <= ncompnodes && ncompnodes <= mcfdata->nnodes);
   assert(1 <= ncomparcs && ncomparcs <= mcfdata->narcs);
   assert(ncommodities > 0);

#ifndef NDEBUG
   /* v -> compv mapping must be all -1 */
   for( v = 0; v < mcfdata->nnodes; v++ )
      assert(compnodeid[v] == -1);
#endif

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &comdemands, ncommodities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compcommodity, ncommodities) );

   /* initialize demand array */
   BMSclearMemoryArray(comdemands, ncommodities);

   /* initialize k -> compk mapping */
   for( k = 0; k < ncommodities; k++ )
      compcommodity[k] = -1;

   /* get LP rows and cols data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* generate v -> compv mapping */
   for( i = 0; i < ncompnodes; i++ )
   {
      v = compnodes[i];
      assert(0 <= v && v < mcfdata->nnodes);
      compnodeid[v] = i;
   }

   /* generate k -> compk mapping */
   ncompcommodities = 0;
   for( i = 0; i < nflowcands; i++ )
   {
      int r;
      int rv;

      r = flowcands[i];
      assert(0 <= r && r < nrows);

      rv = rownodeid[r];
      if( rv >= 0 && compnodeid[rv] >= 0 )
      {
         k = rowcommodity[r];
         assert(0 <= k && k < ncommodities);
         if( compcommodity[k] == -1 )
         {
            compcommodity[k] = ncompcommodities;
            ncompcommodities++;
         }
      }
   }

   /**@todo model type and flow type may be different for each component */
   /* record model and flow type */
   mcfnetwork->modeltype = modeltype;
   mcfnetwork->flowtype = flowtype;

   /* record network size */
   mcfnetwork->nnodes = ncompnodes;
   mcfnetwork->narcs = ncomparcs;
   mcfnetwork->ncommodities = ncompcommodities;

   /* allocate memory for arrays and initialize with default values */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->nodeflowrows, mcfnetwork->nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->nodeflowscales, mcfnetwork->nnodes) );
   for( v = 0; v < mcfnetwork->nnodes; v++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->nodeflowrows[v], mcfnetwork->ncommodities) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->nodeflowscales[v], mcfnetwork->ncommodities) );
      for( k = 0; k < mcfnetwork->ncommodities; k++ )
      {
         mcfnetwork->nodeflowrows[v][k] = NULL;
         mcfnetwork->nodeflowscales[v][k] = 0.0;
      }
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->arccapacityrows, mcfnetwork->narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->arccapacityscales, mcfnetwork->narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->arcsources, mcfnetwork->narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->arctargets, mcfnetwork->narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfnetwork->colcommodity, ncols) );
   for( a = 0; a < mcfnetwork->narcs; a++ )
   {
      mcfnetwork->arccapacityrows[a] = NULL;
      mcfnetwork->arccapacityscales[a] = 0.0;
      mcfnetwork->arcsources[a] = -1;
      mcfnetwork->arctargets[a] = -1;
   }
   BMSclearMemoryArray(mcfnetwork->colcommodity, mcfnetwork->ncommodities);

   /* fill in existing node data */
   for( i = 0; i < nflowcands; i++ )
   {
      int r;
      int rv;

      r = flowcands[i];
      assert(0 <= r && r < nrows);

      rv = rownodeid[r];
      if( rv >= 0 && compnodeid[rv] >= 0 )
      {
         SCIP_Real scale;
         int rk;

         v = compnodeid[rv];
         rk = rowcommodity[r];
         assert(v < mcfnetwork->nnodes);
         assert(0 <= rk && rk < ncommodities);
         assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);

         k = compcommodity[rk];
         assert(0 <= k && k < mcfnetwork->ncommodities);

         /* fill in node -> row assignment */
         SCIP_CALL( SCIPcaptureRow(scip, rows[r]) );
         mcfnetwork->nodeflowrows[v][k] = rows[r];
         scale = flowrowscalars[r];
         if( (flowrowsigns[r] & LHSASSIGNED) != 0 )
            scale *= -1.0;
         if( commoditysigns[rk] == -1 )
            scale *= -1.0;
         mcfnetwork->nodeflowscales[v][k] = scale;
      }
   }

   /* fill in existing arc data */
   for( a = 0; a < mcfnetwork->narcs; a++ )
   {
      SCIP_ROW* capacityrow;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      int rowlen;
      int globala;
      int r;
      int j;

      globala = comparcs[a];
      capacityrow = capacityrows[globala];
      assert(capacityrow != NULL);

      r = SCIProwGetLPPos(capacityrow);
      assert(0 <= r && r < nrows);
      assert((capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);
      assert(mcfdata->rowarcid[r] == globala);

      SCIP_CALL( SCIPcaptureRow(scip, capacityrow) );
      mcfnetwork->arccapacityrows[a] = capacityrow;
      mcfnetwork->arccapacityscales[a] = 1.0;

      /* Scale constraint such that the coefficients for the flow variables are normalized in such a way that coefficients in
       * multiple capacity constraints that belong to the same commodity are (hopefully) equal.
       * This is needed for binary flow variable models in which the demand of each commodity is stored as the coefficient in
       * the capacity constraints. Since it may happen (e.g., due to presolve) that different capacity constraints are scaled
       * differently, we need to find scaling factors to make the demand coefficients of each commodity equal.
       * To do this, we scale the first capacity constraint with +1 and then store the coefficients of the flow variables
       * as target demands for the commodities. Then, we scale the other constraints in such a way that these demands are hit, if possible.
       * Note that for continuous flow variable models, the coefficients in the capacity constraints are usually +1.0.
       * This is achieved automatically by our scaling routine.
       */
      rowcols = SCIProwGetCols(capacityrow);
      rowvals = SCIProwGetVals(capacityrow);
      rowlen = SCIProwGetNLPNonz(capacityrow);
      for( j = 0; j < rowlen; j++ )
      {
         c = SCIPcolGetLPPos(rowcols[j]);
         assert(0 <= c && c < SCIPgetNLPCols(scip));
         k = colcommodity[c];
         if( k >= 0 )
         {
            if( comdemands[k] != 0.0 )
            {
               /* update the scaling factor */
               mcfnetwork->arccapacityscales[a] = comdemands[k]/rowvals[j];
               break;
            }
         }
      }

      /* use negative scaling if we use the left hand side, use positive scaling if we use the right hand side */
      mcfnetwork->arccapacityscales[a] = ABS(mcfnetwork->arccapacityscales[a]);
      if( (capacityrowsigns[r] & LHSASSIGNED) != 0 )
         mcfnetwork->arccapacityscales[a] *= -1.0;

      /* record the commodity demands */
      for( j = 0; j < rowlen; j++ )
      {
         c = SCIPcolGetLPPos(rowcols[j]);
         assert(0 <= c && c < SCIPgetNLPCols(scip));
         k = colcommodity[c];
         if( k >= 0 && comdemands[k] == 0.0 )
            comdemands[k] = mcfnetwork->arccapacityscales[a] * rowvals[j];
      }

      /* copy the source/target node assignment */
      if( mcfdata->arcsources[globala] >= 0 )
      {
         assert(mcfdata->arcsources[globala] < mcfdata->nnodes);
         assert(0 <= compnodeid[mcfdata->arcsources[globala]] && compnodeid[mcfdata->arcsources[globala]] < mcfnetwork->nnodes);
         mcfnetwork->arcsources[a] = compnodeid[mcfdata->arcsources[globala]];
      }
      if( mcfdata->arctargets[globala] >= 0 )
      {
         assert(mcfdata->arctargets[globala] < mcfdata->nnodes);
         assert(0 <= compnodeid[mcfdata->arctargets[globala]] && compnodeid[mcfdata->arctargets[globala]] < mcfnetwork->nnodes);
         mcfnetwork->arctargets[a] = compnodeid[mcfdata->arctargets[globala]];
      }
   }

   /* translate colcommodity array */
   for( c = 0; c < ncols; c++ )
   {
      k = colcommodity[c];
      if( k >= 0 )
         mcfnetwork->colcommodity[c] = compcommodity[k];
      else
         mcfnetwork->colcommodity[c] = -1;
   }

   /* reset v -> compv mapping */
   for( i = 0; i < ncompnodes; i++ )
   {
      assert(0 <= compnodes[i] && compnodes[i] < mcfdata->nnodes);
      assert(compnodeid[compnodes[i]] == i);
      compnodeid[compnodes[i]] = -1;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &compcommodity);
   SCIPfreeBufferArray(scip, &comdemands);

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** displays the MCF network */
static
void mcfnetworkPrint(
   SCIP_MCFNETWORK*      mcfnetwork          /**< MCF network structure */
   )
{
   if( mcfnetwork == NULL )
      printf("MCF network is empty\n");
   else
   {
      int v;
      int a;

      for( v = 0; v < mcfnetwork->nnodes; v++ )
      {
         int k;

         printf("node %2d:\n", v);
         for( k = 0; k < mcfnetwork->ncommodities; k++ )
         {
            printf("  commodity %2d: ", k);
            if( mcfnetwork->nodeflowrows[v][k] != NULL )
            {
               printf("<%s> [%+g]\n", SCIProwGetName(mcfnetwork->nodeflowrows[v][k]), mcfnetwork->nodeflowscales[v][k]);
               /*SCIProwPrint(mcfnetwork->nodeflowrows[v][k], NULL);*/
            }
            else
               printf("-\n");
         }
      }

      for( a = 0; a < mcfnetwork->narcs; a++ )
      {
         printf("arc %2d [%2d -> %2d]: ", a, mcfnetwork->arcsources[a], mcfnetwork->arctargets[a]);
         if( mcfnetwork->arccapacityrows[a] != NULL )
         {
            printf("<%s> [%+g]\n", SCIProwGetName(mcfnetwork->arccapacityrows[a]), mcfnetwork->arccapacityscales[a]);
            /*SCIProwPrint(mcfnetwork->arccapacityrows[a], NULL);*/
         }
         else
            printf("-\n");
      }
   }
}

/** displays commodities and its members */
static
void printCommodities(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char* flowrowsigns     = mcfdata->flowrowsigns;
   unsigned char* capacityrowsigns = mcfdata->capacityrowsigns;
   int            ncommodities     = mcfdata->ncommodities;
   int*           commoditysigns   = mcfdata->commoditysigns;
   int*           colcommodity     = mcfdata->colcommodity;
   int*           rowcommodity     = mcfdata->rowcommodity;
   int*           colarcid         = mcfdata->colarcid;
   int*           rownodeid        = mcfdata->rownodeid;
   int            nnodes           = mcfdata->nnodes;
   SCIP_ROW**     capacityrows     = mcfdata->capacityrows;

   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int k;
   int c;
   int r;
   int a;

   cols = SCIPgetLPCols(scip);
   ncols = SCIPgetNLPCols(scip);
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);

   for( k = 0; k < ncommodities; k++ )
   {
      printf("commodity %d (sign: %+d):\n", k, commoditysigns[k]);

      for( c = 0; c < ncols; c++ )
      {
         if( colcommodity[c] == k )
            printf(" col <%s>: arc %d\n", SCIPvarGetName(SCIPcolGetVar(cols[c])), colarcid != NULL ? colarcid[c] : -1);
      }
      for( r = 0; r < nrows; r++ )
      {
         if( rowcommodity[r] == k )
            printf(" row <%s>: node %d [sign:%+d]\n", SCIProwGetName(rows[r]), rownodeid != NULL ? rownodeid[r] : -1,
                   (flowrowsigns[r] & RHSASSIGNED) != 0 ? +1 : -1);
      }
      printf("\n");
   }

   if( rownodeid != NULL )
   {
      int v;

      for( v = 0; v < nnodes; v++ )
      {
         printf("node %d:\n", v);
         for( r = 0; r < nrows; r++ )
         {
            if( rownodeid[r] == v )
               printf(" row <%s> [sign:%+d]\n", SCIProwGetName(rows[r]), (flowrowsigns[r] & RHSASSIGNED) != 0 ? +1 : -1);
         }
         printf("\n");
      }
   }

   printf("capacities:\n");
   for( a = 0; a < mcfdata->narcs; a++ )
   {
      printf("  arc %d: ", a);
      if( capacityrows[a] != NULL )
      {
         r = SCIProwGetLPPos(capacityrows[a]);
         assert(0 <= r && r < nrows);
         if( (capacityrowsigns[r] & LHSASSIGNED) != 0 )
            printf(" row <%s> [sign:-1]\n", SCIProwGetName(rows[r]));
         else if( (capacityrowsigns[r] & RHSASSIGNED) != 0 )
            printf(" row <%s> [sign:+1]\n", SCIProwGetName(rows[r]));
      }
      else
         printf(" -\n");
   }
   printf("\n");

   printf("unused columns:\n");
   for( c = 0; c < ncols; c++ )
   {
      if( colcommodity[c] == -1 )
      {
         SCIP_VAR* var = SCIPcolGetVar(cols[c]);
         printf(" col <%s> [%g,%g]\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      }
   }
   printf("\n");

   printf("unused rows:\n");
   for( r = 0; r < nrows; r++ )
   {
      if( rowcommodity[r] == -1 && (capacityrowsigns == NULL || (capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) == 0) )
      {
         printf(" row <%s>\n", SCIProwGetName(rows[r]));
         SCIPdebug( SCIPprintRow(scip, rows[r], NULL) );
      }
   }
   printf("\n");
}
#endif

/** comparator method for flow and capacity row candidates */
static
SCIP_DECL_SORTINDCOMP(compCands)
{
   SCIP_Real* rowscores = (SCIP_Real*)dataptr;

   if( rowscores[ind2] < rowscores[ind1] )
      return -1;
   else if( rowscores[ind2] > rowscores[ind1] )
      return +1;
   else
      return 0;
}

/** extracts flow conservation from the LP */
static
SCIP_RETCODE extractFlowRows(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char* flowrowsigns;
   SCIP_Real*     flowrowscalars;
   SCIP_Real*     flowrowscores;
   int*           flowcands;

   SCIP_ROW** rows;
   int nrows;
   int r;

   SCIP_Real maxdualflow;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* allocate temporary memory for extraction data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowrowsigns, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowrowscalars, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowrowscores, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowcands, nrows) );
   flowrowsigns      = mcfdata->flowrowsigns;
   flowrowscalars    = mcfdata->flowrowscalars;
   flowrowscores     = mcfdata->flowrowscores;
   flowcands         = mcfdata->flowcands;

   maxdualflow = 0.0;
   for( r = 0; r < nrows; r++ )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real rowlhs;
      SCIP_Real rowrhs;
      int rowlen;
      int nbinvars;
      int nintvars;
      int nimplintvars;
      int ncontvars;
      SCIP_Real coef;
      SCIP_Bool hasposcoef;
      SCIP_Bool hasnegcoef;
      SCIP_Real absdualsol;
      int i;

      row = rows[r];
      assert(SCIProwGetLPPos(row) == r);

      /* get dual solution, if available */
      absdualsol = SCIProwGetDualsol(row);
      if( absdualsol == SCIP_INVALID )
         absdualsol = 0.0;
      absdualsol = ABS(absdualsol);

      flowrowsigns[r] = 0;
      flowrowscalars[r] = 0.0;
      flowrowscores[r] = 0.0;

      rowlen = SCIProwGetNNonz(row);
      if( rowlen == 0 )
         continue;
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
      rowlhs = SCIProwGetLhs(row);
      rowrhs = SCIProwGetRhs(row);

      /* identify flow conservation constraints */
      coef = ABS(rowvals[0]);
      hasposcoef = FALSE;
      hasnegcoef = FALSE;
      nbinvars = 0;
      nintvars = 0;
      nimplintvars = 0;
      ncontvars = 0;
      for( i = 0; i < rowlen; i++ )
      {
         SCIP_Real absval = ABS(rowvals[i]);
         if( !SCIPisEQ(scip, absval, coef) )
            break;

#if 0
         if( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) != SCIP_VARTYPE_CONTINUOUS )
            break; /*????????????????????*/
#endif
         hasposcoef = hasposcoef || (rowvals[i] > 0.0);
         hasnegcoef = hasnegcoef || (rowvals[i] < 0.0);
         switch( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) )
         {
         case SCIP_VARTYPE_BINARY:
            nbinvars++;
            break;
         case SCIP_VARTYPE_INTEGER:
            nintvars++;
            break;
         case SCIP_VARTYPE_IMPLINT:
            nimplintvars++;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            ncontvars++;
            break;
         default:
            SCIPerrorMessage("unknown variable type\n");
            SCIPABORT();
         }
      }
      if( i == rowlen )
      {
         /* Flow conservation constraints should always be a*x <= -d.
          * If lhs and rhs are finite, both sides are still valid candidates.
          */
         if( !SCIPisInfinity(scip, -rowlhs) )
            flowrowsigns[r] |= LHSPOSSIBLE;
         if( !SCIPisInfinity(scip, rowrhs) )
            flowrowsigns[r] |= RHSPOSSIBLE;
         flowrowscalars[r] = 1.0/coef;
         flowcands[mcfdata->nflowcands] = r;
         mcfdata->nflowcands++;
      }

      /* calculate flow row score */
      if( (flowrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) != 0 )
      {
         flowrowscalars[r] = 1.0;

         /* row is an equation: score +1000 */
         if( (flowrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) == (LHSPOSSIBLE | RHSPOSSIBLE) )
            flowrowscores[r] += 1000.0;

         /* row does not need to be scaled: score +1000 */
         if( SCIPisEQ(scip, flowrowscalars[r], 1.0) )
            flowrowscores[r] += 1000.0;

         /* row has positive and negative coefficients: score +500 */
         if( hasposcoef && hasnegcoef )
            flowrowscores[r] += 500.0;

         /* all variables are of the same type:
          *    continuous: score +1000
          *    integer:    score  +500
          *    binary:     score  +100
          */
         if( ncontvars == rowlen )
            flowrowscores[r] += 1000.0;
         else if( nintvars + nimplintvars == rowlen )
            flowrowscores[r] += 500.0;
         else if( nbinvars == rowlen )
            flowrowscores[r] += 100.0;

         /* the longer the row, the earlier we want to process it: score +10*len/(len+10) */
         flowrowscores[r] += 10.0*rowlen/(rowlen+10.0);

         assert(flowrowscores[r] > 0.0);

         /* update maximum dual solution value for additional score tie breaking */
         maxdualflow = MAX(maxdualflow, absdualsol);

         /**@todo go through list of several model types, depending on the current model type throw away invalid constraints
          *       instead of assigning a low score
          */
      }
   }

   /* apply additional score tie breaking using the dual solutions */
   if( SCIPisPositive(scip, maxdualflow) )
   {
      int i;

      for( i = 0; i < mcfdata->nflowcands; i++ )
      {
         SCIP_Real dualsol;

         r = flowcands[i];
         assert(0 <= r && r < nrows);
         dualsol = SCIProwGetDualsol(rows[r]);
         if( dualsol == SCIP_INVALID )
            dualsol = 0.0;
         else if( flowrowsigns[r] == (LHSPOSSIBLE | RHSPOSSIBLE) )
            dualsol = ABS(dualsol);
         else if( flowrowsigns[r] == RHSPOSSIBLE )
            dualsol = -dualsol;
         flowrowscores[r] += dualsol/maxdualflow + 1.0;
         assert(flowrowscores[r] > 0.0);
      }
   }

   /* sort candidates by score */
   SCIPsortInd(mcfdata->flowcands, compCands, (void*)flowrowscores, mcfdata->nflowcands);

   SCIPdebugMessage("flow conservation candidates:\n");
   for( r = 0; r < mcfdata->nflowcands; r++ )
   {
      //SCIPdebug(SCIPprintRow(scip, rows[mcfdata->flowcands[r]], NULL));
      SCIPdebugMessage("%4d [score: %2g]: %s\n", r, flowrowscores[mcfdata->flowcands[r]], SCIProwGetName(rows[mcfdata->flowcands[r]]));
   }

   return SCIP_OKAY;
}

/** extracts flow conservation and capacity rows from the LP */
static
SCIP_RETCODE extractCapacityRows(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int*              colcommodity = mcfdata->colcommodity;
   int               ncommodities = mcfdata->ncommodities;
   SCIP_MCFMODELTYPE modeltype    = mcfdata->modeltype;

   unsigned char* capacityrowsigns;
   SCIP_Real*     capacityrowscores;
   int*           capacitycands;

   SCIP_ROW** rows;
   int nrows;
   int r;

   SCIP_Real maxdualcapacity;
   int maxcolspercommodity;
   int *ncolspercommodity;
   int ndirectedcands;
   int nundirectedcands;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* allocate temporary memory for extraction data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->capacityrowsigns, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->capacityrowscores, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->capacitycands, nrows) );
   capacityrowsigns  = mcfdata->capacityrowsigns;
   capacityrowscores = mcfdata->capacityrowscores;
   capacitycands     = mcfdata->capacitycands;

   /* allocate temporary memory for model type identification */
   SCIP_CALL( SCIPallocBufferArray(scip, &ncolspercommodity, ncommodities) );

   /* identify model type and set the maximal number of flow variables per capacity constraint and commodity */
   switch( modeltype )
   {
   case SCIP_MCFMODELTYPE_AUTO:
      maxcolspercommodity = 2; /* will be set to 1 later if we detect that the network is directed */
      break;
   case SCIP_MCFMODELTYPE_DIRECTED:
      maxcolspercommodity = 1;
      break;
   case SCIP_MCFMODELTYPE_UNDIRECTED:
      maxcolspercommodity = 2;
      break;
   default:
      SCIPerrorMessage("invalid parameter value %d for model type\n", modeltype);
      return SCIP_INVALIDDATA;
   }

   /**@todo instead of discarding rows with too many entries per commodity, decrease its score */
   maxdualcapacity = 0.0;
   ndirectedcands = 0;
   nundirectedcands = 0;
   for( r = 0; r < nrows; r++ )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real rowlhs;
      SCIP_Real rowrhs;
      int rowlen;
      int nposflowcoefs;
      int nnegflowcoefs;
      int nposcapacitycoefs;
      int nnegcapacitycoefs;
      int nbadcoefs;
      SCIP_Real sameflowcoef;
      SCIP_Real sameabsflowcoef;
      SCIP_Real absdualsol;
      unsigned int rowsign;
      int i;

      row = rows[r];
      assert(SCIProwGetLPPos(row) == r);

      /* get dual solution, if available */
      absdualsol = SCIProwGetDualsol(row);
      if( absdualsol == SCIP_INVALID )
         absdualsol = 0.0;
      absdualsol = ABS(absdualsol);

      capacityrowsigns[r] = 0;
      capacityrowscores[r] = 0.0;

      rowlen = SCIProwGetNNonz(row);
      if( rowlen == 0 )
         continue;
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
      rowlhs = SCIProwGetLhs(row);
      rowrhs = SCIProwGetRhs(row);

      /* reset commodity counting array */
      BMSclearMemoryArray(ncolspercommodity, ncommodities);

      /* identify capacity constraints */
      nposflowcoefs = 0;
      nnegflowcoefs = 0;
      nposcapacitycoefs = 0;
      nnegcapacitycoefs = 0;
      nbadcoefs = 0;
      sameflowcoef = 0.0;
      sameabsflowcoef = 0.0;
      rowsign = 0;
      if( !SCIPisInfinity(scip, -rowlhs) )
         rowsign |= LHSPOSSIBLE;
      if( !SCIPisInfinity(scip, rowrhs) )
         rowsign |= RHSPOSSIBLE;
      for( i = 0; i < rowlen; i++ )
      {
         int c;
         int k;

         c = SCIPcolGetLPPos(rowcols[i]);
         assert(0 <= c && c < SCIPgetNLPCols(scip));

         /* check if this is a flow variable */
         k = colcommodity[c];
         assert(-1 <= k && k < ncommodities);
         if( k >= 0 )
         {
            SCIP_Real abscoef;

            abscoef = ABS(rowvals[i]);
            if( sameflowcoef == 0.0 )
               sameflowcoef = rowvals[i];
            else if( !SCIPisEQ(scip, sameflowcoef, rowvals[i]) )
               sameflowcoef = SCIP_REAL_MAX;
            if( sameabsflowcoef == 0.0 )
               sameabsflowcoef = abscoef;
            else if( !SCIPisEQ(scip, sameabsflowcoef, abscoef) )
               sameabsflowcoef = SCIP_REAL_MAX;

            if( rowvals[i] > 0.0 )
               nposflowcoefs++;
            else
               nnegflowcoefs++;

            ncolspercommodity[k]++;
            if( ncolspercommodity[k] > maxcolspercommodity )
            {
               /* the row has too many entries in this commodities: discard capacity row */
               break;
            }
            else if( ncolspercommodity[k] == 2 )
               capacityrowsigns[r] |= UNDIRECTED;
         }
         else if( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) != SCIP_VARTYPE_CONTINUOUS )
         {
            /* an integer variable which is not a flow variable can be used as capacity variable */
            if( rowvals[i] > 0.0 )
               nposcapacitycoefs++;
            else
               nnegcapacitycoefs++;
         }
         else
         {
            /* a continuous variable which is not a flow variable cannot be used for anything: this is bad! */
            nbadcoefs++;
         }
      }

      /* check if this is a valid capacity constraint */
      if( i == rowlen && rowsign != 0 && nposflowcoefs + nnegflowcoefs > 0 )
      {
         capacityrowsigns[r] |= rowsign;
         capacitycands[mcfdata->ncapacitycands] = r;
         mcfdata->ncapacitycands++;

         /* calculate capacity row score */
         capacityrowscores[r] = 1.0;

         /* row is of type f - c*x <= b: score +1000 */
         if( (capacityrowsigns[r] & RHSPOSSIBLE) != 0 && nnegflowcoefs == 0 && nposcapacitycoefs == 0 && nnegcapacitycoefs > 0 )
            capacityrowscores[r] += 1000.0;
         if( (capacityrowsigns[r] & LHSPOSSIBLE) != 0 && nposflowcoefs == 0 && nposcapacitycoefs > 0 && nnegcapacitycoefs == 0 )
            capacityrowscores[r] += 1000.0;

         /* row has no continuous variables that are not flow variables: score +1000 */
         if( nbadcoefs == 0 )
            capacityrowscores[r] += 1000.0;

         /* almost all commodities are covered: score +500*(nposflow+nnegflow)/ncommodities */
         capacityrowscores[r] += 500.0 * (nposflowcoefs+nnegflowcoefs)/(SCIP_Real)ncommodities;

         /* all coefficients of flow variables are +1 or all are -1: score +500 */
         if( SCIPisEQ(scip, ABS(sameflowcoef), 1.0) )
            capacityrowscores[r] += 500.0;

         /* all coefficients of flow variables are equal: score +200 */
         if( sameflowcoef != 0.0 && sameflowcoef != SCIP_REAL_MAX )
            capacityrowscores[r] += 200.0;

         /* all coefficients of flow variables are +1 or -1: score +100 */
         if( SCIPisEQ(scip, sameabsflowcoef, 1.0) )
            capacityrowscores[r] += 100.0;

         /* all coefficients of flow variables are equal in their absolute values: score +50 */
         if( sameflowcoef != 0.0 && sameflowcoef != SCIP_REAL_MAX )
            capacityrowscores[r] += 50.0;

         /* flow coefficients are mostly of the same sign: score +20*max(npos,nneg)/(npos+nneg) */
         capacityrowscores[r] += 20.0 * MAX(nposflowcoefs, nnegflowcoefs)/(SCIP_Real)(nposflowcoefs+nnegflowcoefs);

         /* capacity coefficients are mostly of the same sign: score +10*max(npos,nneg)/(npos+nneg+1) */
         capacityrowscores[r] += 10.0 * MAX(nposcapacitycoefs, nnegcapacitycoefs)/(SCIP_Real)(nposcapacitycoefs+nnegcapacitycoefs+1.0);

         /* row is a <= row with non-negative right hand side: score +10 */
         if( (capacityrowsigns[r] & RHSPOSSIBLE) != 0 && !SCIPisNegative(scip, rowrhs)  )
            capacityrowscores[r] += 10.0;

         /* row is an inequality: score +10 */
         if( SCIPisInfinity(scip, -rowlhs) != SCIPisInfinity(scip, rowrhs) )
            capacityrowscores[r] += 10.0;

         assert(capacityrowscores[r] > 0.0);

         /* update maximum dual solution value for additional score tie breaking */
         maxdualcapacity = MAX(maxdualcapacity, absdualsol);

         /* if the model type should be detected automatically, count the number of directed and undirected capacity candidates */
         if( modeltype == SCIP_MCFMODELTYPE_AUTO )
         {
            assert(maxcolspercommodity == 2);
            if( (capacityrowsigns[r] & UNDIRECTED) != 0 )
               nundirectedcands++;
            else
               ndirectedcands++;
         }
      }
   }

   /* if the model type should be detected automatically, decide it by a majority vote */
   if( modeltype == SCIP_MCFMODELTYPE_AUTO )
   {
      if( ndirectedcands > nundirectedcands )
         modeltype = SCIP_MCFMODELTYPE_DIRECTED;
      else
         modeltype = SCIP_MCFMODELTYPE_UNDIRECTED;

      SCIPdebugMessage("detected model type: %s (%d directed candidates, %d undirected candidates)\n",
                       modeltype == SCIP_MCFMODELTYPE_DIRECTED ? "directed" : "undirected", ndirectedcands, nundirectedcands);

      if( modeltype == SCIP_MCFMODELTYPE_DIRECTED )
      {
         int i;

         /* discard all undirected arcs */
         for( i = 0; i < mcfdata->ncapacitycands; i++ )
         {
            r = capacitycands[i];
            assert(0 <= r && r < nrows);
            if( (capacityrowsigns[r] & UNDIRECTED) != 0 )
            {
               SCIPdebugMessage(" -> discarded capacity row %d <%s>: undirected flow\n", r, SCIProwGetName(rows[r]));
               capacityrowsigns[r] |= DISCARDED;
            }
         }
      }

      /* record the detected model type */
      mcfdata->modeltype = modeltype;
   }

   /* apply additional score tie breaking using the dual solutions */
   if( SCIPisPositive(scip, maxdualcapacity) )
   {
      int i;

      for( i = 0; i < mcfdata->ncapacitycands; i++ )
      {
         SCIP_Real dualsol;

         r = capacitycands[i];
         assert(0 <= r && r < nrows);
         dualsol = SCIProwGetDualsol(rows[r]);
         if( dualsol == SCIP_INVALID )
            dualsol = 0.0;
         else if( capacityrowsigns[r] == (LHSPOSSIBLE | RHSPOSSIBLE) )
            dualsol = ABS(dualsol);
         else if( capacityrowsigns[r] == RHSPOSSIBLE )
            dualsol = -dualsol;
         capacityrowscores[r] += dualsol/maxdualcapacity;
         assert(capacityrowscores[r] > 0.0);
      }
   }

   /* sort candidates by score */
   SCIPsortInd(mcfdata->capacitycands, compCands, (void*)capacityrowscores, mcfdata->ncapacitycands);

   SCIPdebugMessage("capacity candidates:\n");
   for( r = 0; r < mcfdata->ncapacitycands; r++ )
   {
      //SCIPdebug(SCIPprintRow(scip, rows[mcfdata->capacitycands[r]], NULL));
      SCIPdebugMessage("%4d [score: %2g]: %s\n", r, capacityrowscores[mcfdata->capacitycands[r]], SCIProwGetName(rows[mcfdata->capacitycands[r]]));
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &ncolspercommodity);

   return SCIP_OKAY;
}

/** creates a new commodity */
static
SCIP_RETCODE createNewCommodity(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   /* get memory for commoditysigns array */
   assert(mcfdata->ncommodities <= mcfdata->commoditysignssize);
   if( mcfdata->ncommodities == mcfdata->commoditysignssize )
   {
      mcfdata->commoditysignssize = MAX(2*mcfdata->commoditysignssize, mcfdata->ncommodities+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &mcfdata->commoditysigns, mcfdata->commoditysignssize) );
   }
   assert(mcfdata->ncommodities < mcfdata->commoditysignssize);

   /* create commodity */
   SCIPdebugMessage("**** creating new commodity %d ****\n", mcfdata->ncommodities);
   mcfdata->commoditysigns[mcfdata->ncommodities] = 0;
   mcfdata->ncommodities++;

   return SCIP_OKAY;
}

/** adds the given flow row and all involved columns to the current commodity */
static
void addFlowrowToCommodity(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             row,                /**< flow row to add to current commodity */
   unsigned char         rowsign             /**< possible flow row signs to use */
   )
{
   unsigned char* flowrowsigns   = mcfdata->flowrowsigns;
   SCIP_Bool*     plusflow       = mcfdata->plusflow;
   SCIP_Bool*     minusflow      = mcfdata->minusflow;
   int            ncommodities   = mcfdata->ncommodities;
   int*           commoditysigns = mcfdata->commoditysigns;
   int*           colcommodity   = mcfdata->colcommodity;
   int*           rowcommodity   = mcfdata->rowcommodity;
   int*           newcols        = mcfdata->newcols;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int rowscale;
   int r;
   int k;
   int i;

   k = ncommodities-1;
   assert(k >= 0);

   r = SCIProwGetLPPos(row);
   assert(r >= 0);
   assert(rowcommodity[r] == -1);
   assert((flowrowsigns[r] | rowsign) == flowrowsigns[r]);
   assert((rowsign & (LHSPOSSIBLE | RHSPOSSIBLE)) == rowsign);
   assert(rowsign != 0);

   /* if the row is only useable as flow row in one direction, we cannot change the sign
    * of the whole commodity anymore
    */
   if( (flowrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) != (LHSPOSSIBLE | RHSPOSSIBLE) )
      commoditysigns[k] = +1; /* we cannot switch directions */

   /* decide the sign (direction) of the row */
   if( rowsign == LHSPOSSIBLE )
      rowsign = LHSASSIGNED;
   else if( rowsign == RHSPOSSIBLE )
      rowsign = RHSASSIGNED;
   else
   {
      SCIP_Real dualsol = SCIProwGetDualsol(row);

      assert(rowsign == (LHSPOSSIBLE | RHSPOSSIBLE));

      /* if we have a valid non-zero dual solution, choose the side which is tight */
      if( !SCIPisZero(scip, dualsol) && dualsol != SCIP_INVALID )
      {
         if( dualsol > 0.0 )
            rowsign = LHSASSIGNED;
         else
            rowsign = RHSASSIGNED;
      }
      else
      {
         SCIP_Real rowlhs = SCIProwGetLhs(row);
         SCIP_Real rowrhs = SCIProwGetRhs(row);

         /* choose row sign such that we get a*x <= -d with d non-negative */
         if( rowrhs < 0.0 )
            rowsign = RHSASSIGNED;
         else if( rowlhs > 0.0 )
            rowsign = LHSASSIGNED;
         else
            rowsign = RHSASSIGNED; /* if we are still undecided, choose rhs */
      }
   }
   if( rowsign == RHSASSIGNED )
      rowscale = +1;
   else
      rowscale = -1;
   flowrowsigns[r] |= rowsign;

   SCIPdebugMessage("adding flow row %d <%s> with sign %+d to commodity %d [score:%g]\n",
                    r, SCIProwGetName(row), rowscale, k, mcfdata->flowrowscores[r]);
   SCIPdebug( SCIPprintRow(scip, row, NULL); ); /*???????????????????*/

   /* add row to commodity */
   rowcommodity[r] = k;
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   for( i = 0; i < rowlen; i++ )
   {
      SCIP_Real val;
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* assign column to commodity */
      if( colcommodity[c] == -1 )
      {
         assert(!plusflow[c]);
         assert(!minusflow[c]);
         assert(mcfdata->nnewcols < SCIPgetNLPCols(scip));
         colcommodity[c] = k;
         newcols[mcfdata->nnewcols] = c;
         mcfdata->nnewcols++;
      }
      assert(colcommodity[c] == k);

      /* update plusflow/minusflow */
      val = rowscale * rowvals[i];
      if( val > 0.0 )
      {
         assert(!plusflow[c]);
         plusflow[c] = TRUE;
      }
      else
      {
         assert(!minusflow[c]);
         minusflow[c] = TRUE;
      }
   }
}

/** deletes the last commodity, which must only have a single node, and removes the flow row again from the system */
static
void deleteSingleNodeCommodity(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             row                 /**< single flow row of the commodity */
   )
{
   unsigned char* flowrowsigns = mcfdata->flowrowsigns;
   SCIP_Bool*     plusflow     = mcfdata->plusflow;
   SCIP_Bool*     minusflow    = mcfdata->minusflow;
   int            ncommodities = mcfdata->ncommodities;
   int*           colcommodity = mcfdata->colcommodity;
   int*           rowcommodity = mcfdata->rowcommodity;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int r;
   int k;
   int i;

   k = ncommodities-1;
   assert(k >= 0);

   r = SCIProwGetLPPos(row);
   assert(r >= 0);
   assert(rowcommodity[r] == k);
   assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);

   SCIPdebugMessage("deleting commodity %d with single flow row %d <%s>\n", k, r, SCIProwGetName(row));

   /* remove the lhs/rhs assignment */
   flowrowsigns[r] &= ~(LHSASSIGNED | RHSASSIGNED);

   /* remove row from commodity */
   rowcommodity[r] = -1;
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   for( i = 0; i < rowlen; i++ )
   {
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* remove column from commodity */
      assert(colcommodity[c] == k);
      colcommodity[c] = -1;

      /* reset plusflow/minusflow */
      assert(plusflow[c] != minusflow[c]);
      plusflow[c] = FALSE;
      minusflow[c] = FALSE;
   }

   /* get rid of commodity */
   mcfdata->ncommodities--;
}

/** checks whether the given row fits into the given commodity and returns the possible flow row signs */
static
void getFlowrowFit(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             row,                /**< flow row to check */
   int                   k,                  /**< commodity that the flow row should enter */
   unsigned char*        rowsign             /**< pointer to store the possible flow row signs */
   )
{
   unsigned char* flowrowsigns = mcfdata->flowrowsigns;
   SCIP_Bool*     plusflow     = mcfdata->plusflow;
   SCIP_Bool*     minusflow    = mcfdata->minusflow;
   int*           colcommodity = mcfdata->colcommodity;
   int*           rowcommodity = mcfdata->rowcommodity;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int flowrowsign;
   int r;
   int j;

   *rowsign = 0;

   r = SCIProwGetLPPos(row);
   assert(0 <= r && r < SCIPgetNLPRows(scip));

   /* ignore rows that are already used */
   if( rowcommodity[r] != -1 )
      return;

   /* check if row is an available flow row */
   flowrowsign = flowrowsigns[r];
   assert((flowrowsign & (LHSPOSSIBLE | RHSPOSSIBLE | DISCARDED)) == flowrowsign);
   if( (flowrowsign & DISCARDED) != 0 )
      return;
   if( (flowrowsign & (LHSPOSSIBLE | RHSPOSSIBLE)) == 0 )
      return;

   /* check whether the row fits w.r.t. the signs of the coefficients */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   for( j = 0; j < rowlen && flowrowsign != 0; j++ )
   {
      int rowc;

      rowc = SCIPcolGetLPPos(rowcols[j]);
      assert(0 <= rowc && rowc < SCIPgetNLPCols(scip));

      /* check if column already belongs to the same commodity */
      if( colcommodity[rowc] == k )
      {
         /* column only fits if it is not yet present with the same sign */
         if( plusflow[rowc] )
         {
            /* column must not be included with positive sign */
            if( rowvals[j] > 0.0 )
               flowrowsign &= ~RHSPOSSIBLE;
            else
               flowrowsign &= ~LHSPOSSIBLE;
         }
         if( minusflow[rowc] )
         {
            /* column must not be included with negative sign */
            if( rowvals[j] > 0.0 )
               flowrowsign &= ~LHSPOSSIBLE;
            else
               flowrowsign &= ~RHSPOSSIBLE;
         }
      }
      else if( colcommodity[rowc] != -1 )
      {
         /* column does not fit if it already belongs to a different commodity */
         flowrowsign = 0;
      }
   }

   if( flowrowsign == 0 )
   {
      /* we can discard the row, since it can also not be member of a different commodity */
      SCIPdebugMessage(" -> discard flow row %d <%s>\n", r, SCIProwGetName(row));
      flowrowsigns[r] |= DISCARDED;
   }

   *rowsign = flowrowsign;
}

/** returns a flow conservation row that fits into the current commodity, or NULL */
static
void getNextFlowrow(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW**            nextrow,            /**< pointer to store next row */
   unsigned char*        nextrowsign         /**< pointer to store possible signs of next row */
   )
{
   SCIP_Real* flowrowscores = mcfdata->flowrowscores;
   SCIP_Bool* plusflow      = mcfdata->plusflow;
   SCIP_Bool* minusflow     = mcfdata->minusflow;
   int*       newcols       = mcfdata->newcols;
   int        ncommodities  = mcfdata->ncommodities;

   SCIP_COL** cols;
   int k;

   assert(nextrow != NULL);
   assert(nextrowsign != NULL);

   *nextrow = NULL;
   *nextrowsign = 0;

   k = ncommodities-1;

   cols = SCIPgetLPCols(scip);
   assert(cols != NULL);

   /* check if there are any columns left in the commodity that have not yet been inspected for incident flow rows */
   while( mcfdata->nnewcols > 0 )
   {
      SCIP_COL* col;
      SCIP_ROW** colrows;
      int collen;
      SCIP_ROW* bestrow;
      unsigned char bestrowsign;
      SCIP_Real bestscore;
      int c;
      int i;

      /* pop next new column from stack */
      c = newcols[mcfdata->nnewcols-1];
      mcfdata->nnewcols--;
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* check if this columns already as both signs */
      assert(plusflow[c] || minusflow[c]);
      if( plusflow[c] && minusflow[c] )
         continue;

      /* check whether column is incident to a valid flow row that fits into the current commodity */
      bestrow = NULL;
      bestrowsign = 0;
      bestscore = 0.0;
      col = cols[c];
      colrows = SCIPcolGetRows(col);
      collen = SCIPcolGetNLPNonz(col);
      for( i = 0; i < collen; i++ )
      {
         SCIP_ROW* row;
         unsigned char flowrowsign;

         row = colrows[i];

         /* check if row fits into the current commodity */
         getFlowrowFit(scip, mcfdata, row, k, &flowrowsign);

         /* do we have a winner? */
         if( flowrowsign != 0 )
         {
            int r = SCIProwGetLPPos(row);
            assert(0 <= r && r < SCIPgetNLPRows(scip));
            assert(flowrowscores[r] > 0.0);

            if( flowrowscores[r] > bestscore )
            {
               bestrow = row;
               bestrowsign = flowrowsign;
               bestscore = flowrowscores[r];
            }
         }
      }

      /* if there was a valid row for this column, pick the best one
       * Note: This is not the overall best row, only the one for the first column that has a valid row.
       *       However, picking the overall best row seems to be too expensive
       */
      if( bestrow != NULL )
      {
         assert(bestscore > 0.0);
         assert(bestrowsign != 0);
         *nextrow = bestrow;
         *nextrowsign = bestrowsign;
         break;
      }
   }
}

/** extracts flow conservation rows and puts them into commodities */
static
SCIP_RETCODE extractFlow(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int* flowcands = mcfdata->flowcands;

   SCIP_Bool* plusflow;
   SCIP_Bool* minusflow;
   int* colcommodity;
   int* rowcommodity;

   SCIP_ROW** rows;
   int nrows;
   int ncols;
   int i;
   int c;
   int r;

   /* get LP data */
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->plusflow, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->minusflow, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->colcommodity, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->rowcommodity, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->newcols, ncols) );
   plusflow = mcfdata->plusflow;
   minusflow = mcfdata->minusflow;
   colcommodity = mcfdata->colcommodity;
   rowcommodity = mcfdata->rowcommodity;

   /* 3. Extract network structure of flow conservation constraints:
    *    (a) Initialize plusflow[c] = minusflow[c] = FALSE for all columns c and other local data.
    */
   BMSclearMemoryArray(plusflow, ncols);
   BMSclearMemoryArray(minusflow, ncols);
   for( c = 0; c < ncols; c++ )
      colcommodity[c] = -1;
   for( r = 0; r < nrows; r++ )
      rowcommodity[r] = -1;

   /*    (b) As long as there are flow conservation candidates left:
    *        (i) Create new commodity and use first flow conservation constraint as newrow.
    *       (ii) Add newrow to commodity, update pluscom/minuscom accordingly.
    *      (iii) For the newly added columns search for an incident flow conservation constraint. Pick the one of highest ranking.
    *       (iv) If found, set newrow to this row and goto (ii).
    */
   for( i = 0; i < mcfdata->nflowcands; i++ )
   {
      SCIP_ROW* newrow;
      SCIP_ROW* lastrow;
      unsigned char newrowsign;
      int nnodes;

      r = flowcands[i];
      assert(0 <= r && r < nrows);
      newrow = rows[r];

      /* check if row fits into a new commodity */
      getFlowrowFit(scip, mcfdata, newrow, mcfdata->ncommodities, &newrowsign);
      if( newrowsign == 0 )
         continue;

      /* start new commodity */
      SCIP_CALL( createNewCommodity(scip, mcfdata) );
      nnodes = 0;
      lastrow = NULL;

      /* fill commodity with flow conservation constraints */
      do
      {
         /* add new row to commodity */
         addFlowrowToCommodity(scip, mcfdata, newrow, newrowsign);
         lastrow = newrow;
         nnodes++;

         /* get next row to add */
         getNextFlowrow(scip, mcfdata, &newrow, &newrowsign);
      }
      while( newrow != NULL );

      SCIPdebugMessage(" -> finished commodity %d: identified %d nodes\n", mcfdata->ncommodities-1, nnodes);

      /* if the commodity has only one node, delete it again */
      if( nnodes == 1 )
         deleteSingleNodeCommodity(scip, mcfdata, lastrow);
   }

   return SCIP_OKAY;
}

/** identifies capacity constraints for the arcs and assigns arc ids to columns and capacity constraints */
static
SCIP_RETCODE extractCapacities(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char*    capacityrowsigns  = mcfdata->capacityrowsigns;
   SCIP_Real*        capacityrowscores = mcfdata->capacityrowscores;
   int*              colcommodity      = mcfdata->colcommodity;
   int*              rowcommodity      = mcfdata->rowcommodity;

   int* colarcid;
   int* rowarcid;

   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;

   int r;
   int c;
   int i;

   assert(mcfdata->narcs == 0);

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* allocate temporary memory for extraction data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->colarcid, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->rowarcid, nrows) );
   colarcid = mcfdata->colarcid;
   rowarcid = mcfdata->rowarcid;

   /* initialize arcid arrays */
   for( c = 0; c < ncols; c++ )
      colarcid[c] = -1;
   for( r = 0; r < nrows; r++ )
      rowarcid[r] = -1;

   /**@todo use capacity candidates in their score order
    *       instead of looping through the used flow columns and their column vectors
    */
   /* for each column, search for a capacity constraint */
   for( c = 0; c < ncols; c++ )
   {
      SCIP_ROW* bestcapacityrow;
      SCIP_Real bestscore;
      SCIP_ROW** colrows;
      int collen;

      /* ignore columns that are not flow variables */
      if( colcommodity[c] == -1 )
         continue;

      /* ignore columns that are already assigned to an arc */
      if( colarcid[c] >= 0 )
         continue;

      /* scan the column to search for valid capacity constraints */
      bestcapacityrow = NULL;
      bestscore = 0.0;
      colrows = SCIPcolGetRows(cols[c]);
      collen = SCIPcolGetNLPNonz(cols[c]);
      for( i = 0; i < collen; i++ )
      {
         r = SCIProwGetLPPos(colrows[i]);
         assert(0 <= r && r < nrows);

         /* row must not be already assigned */
         assert((capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) == 0);

         /* ignore rows that are not capacity candidates */
         if( (capacityrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) == 0 )
            continue;

         /* ignore discarded rows */
         if( (capacityrowsigns[r] & DISCARDED) != 0 )
            continue;

         /* ignore rows that are already used as flow conservation constraints */
         if( rowcommodity[r] != -1 )
            continue;

         /* check if this capacity candidate has better score */
         assert(capacityrowscores[r] > 0.0);
         if( capacityrowscores[r] > bestscore )
         {
            bestcapacityrow = colrows[i];
            bestscore = capacityrowscores[r];
         }
      }

      /* if no capacity row has been found, leave the column unassigned */
      if( bestcapacityrow != NULL )
      {
         SCIP_COL** rowcols;
         int rowlen;

         /* store the row */
         assert(mcfdata->narcs <= mcfdata->capacityrowssize);
         if( mcfdata->narcs == mcfdata->capacityrowssize )
         {
            mcfdata->capacityrowssize = MAX(2*mcfdata->capacityrowssize, mcfdata->narcs+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, &mcfdata->capacityrows, mcfdata->capacityrowssize) );
         }
         assert(mcfdata->narcs < mcfdata->capacityrowssize);
         mcfdata->capacityrows[mcfdata->narcs] = bestcapacityrow;

         /* assign the capacity row to a new arc id */
         r = SCIProwGetLPPos(bestcapacityrow);
         assert(0 <= r && r < nrows);
         rowarcid[r] = mcfdata->narcs;

         /* decide which sign to use */
         if( (capacityrowsigns[r] & RHSPOSSIBLE) != 0 )
            capacityrowsigns[r] |= RHSASSIGNED;
         else
         {
            assert((capacityrowsigns[r] & LHSPOSSIBLE) != 0);
            capacityrowsigns[r] |= LHSASSIGNED;
         }

         SCIPdebugMessage("assigning capacity row %d <%s> with sign %+d to arc %d [score:%g]\n",
                          r, SCIProwGetName(bestcapacityrow), (capacityrowsigns[r] & RHSASSIGNED) != 0 ? +1 : -1, mcfdata->narcs,
                          mcfdata->capacityrowscores[r]);

         /* assign all involved flow variables to the new arc id */
         /*SCIPdebugMessage(" -> flow:");*/
         rowcols = SCIProwGetCols(bestcapacityrow);
         rowlen = SCIProwGetNLPNonz(bestcapacityrow);
         for( i = 0; i < rowlen; i++ )
         {
            int rowc;

            rowc = SCIPcolGetLPPos(rowcols[i]);
            assert(0 <= rowc && rowc < ncols);

            /* due to aggregations in preprocessing it may happen that a flow variable appears in multiple capacity constraints;
             * in this case, assign it to the first that has been found, and count the number of occurrences
             */
            if( colcommodity[rowc] >= 0 )
            {
               if( colarcid[rowc] == -1 )
               {
                  /*SCIPdebug( printf(" x%d<%s>[%d]", rowc, SCIPvarGetName(SCIPcolGetVar(rowcols[i])), colcommodity[rowc]) );*/
                  colarcid[rowc] = mcfdata->narcs;
               }
            }
         }
         /*SCIPdebug( printf("\n") );*/

         /* increase number of arcs */
         mcfdata->narcs++;
      }
      else
      {
         SCIPdebugMessage("no capacity row found for column x%d <%s> in commodity %d\n", c, SCIPvarGetName(SCIPcolGetVar(cols[c])), colcommodity[c]);
      }
   }

   return SCIP_OKAY;
}

/** collects all flow columns of all commodities (except the one of the base row) that are incident to the node described by the given flow row */
static
void collectIncidentFlowCols(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             flowrow,            /**< flow conservation constraint that defines the node */
   int                   basecommodity       /**< commodity of the base row */
   )
{
   int*           colcommodity  = mcfdata->colcommodity;
   int*           colarcid      = mcfdata->colarcid;
   int*           newcols       = mcfdata->newcols;
   SCIP_ROW**     capacityrows  = mcfdata->capacityrows;
   SCIP_Bool*     colisincident = mcfdata->colisincident;

   SCIP_COL** rowcols;
   int rowlen;
   int i;

#ifndef NDEBUG
   /* check that the marker array is correctly initialized */
   for( i = 0; i < SCIPgetNLPCols(scip); i++ )
      assert(!colisincident[i]);
#endif

   /* loop through all flow columns in the flow conservation constraint */
   rowcols = SCIProwGetCols(flowrow);
   rowlen = SCIProwGetNLPNonz(flowrow);
   mcfdata->nnewcols = 0;
   for( i = 0; i < rowlen; i++ )
   {
      SCIP_COL** capacityrowcols;
      int capacityrowlen;
      int arcid;
      int c;
      int j;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* get arc id of the column in the flow conservation constraint */
      arcid = colarcid[c];
      if( arcid == -1 )
         continue;
      assert(arcid < mcfdata->narcs);

      /* collect flow variables in the capacity constraint of this arc */
      assert(capacityrows[arcid] != NULL);
      capacityrowcols = SCIProwGetCols(capacityrows[arcid]);
      capacityrowlen = SCIProwGetNLPNonz(capacityrows[arcid]);
      for( j = 0; j < capacityrowlen; j++ )
      {
         int caprowc;

         caprowc = SCIPcolGetLPPos(capacityrowcols[j]);
         assert(0 <= caprowc && caprowc < SCIPgetNLPCols(scip));

         /* ignore columns that do not belong to a commodity, i.e., are not flow variables */
         if( colcommodity[caprowc] == -1 )
         {
            assert(colarcid[caprowc] == -1);
            continue;
         }
         assert(colarcid[caprowc] <= arcid); /* colarcid < arcid if column belongs to multiple arcs, for example, due to an aggregation in presolving */

         /* ignore columns in the same commodity as the base row */
         if( colcommodity[caprowc] == basecommodity )
            continue;

         /* if not already done, collect the column */
         if( !colisincident[caprowc] )
         {
            assert(mcfdata->nnewcols < SCIPgetNLPCols(scip));
            colisincident[caprowc] = TRUE;
            newcols[mcfdata->nnewcols] = caprowc;
            mcfdata->nnewcols++;
         }
      }
   }
}

/** compares given row against a base node flow row and calculates a similarity score;
 *  score is 0.0 if the rows are incompatible
 */
static
SCIP_RETCODE getNodeSilimarityScore(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   int                   baserowlen,         /**< length of base node flow row */
   int*                  basearcpattern,     /**< arc patern of base node flow row */
   SCIP_ROW*             row,                /**< row to compare against base node flow row */
   SCIP_Real*            score,              /**< pointer to store the similarity score */
   SCIP_Bool*            invertcommodity     /**< pointer to store whether the arcs in the commodity of the row have
                                              *   to be inverted for the row to be compatible to the base row */
   )
{
   unsigned char*    flowrowsigns   = mcfdata->flowrowsigns;
   int*              commoditysigns = mcfdata->commoditysigns;
   int               narcs          = mcfdata->narcs;
   int*              rowcommodity   = mcfdata->rowcommodity;
   int*              colarcid       = mcfdata->colarcid;
   int*              arcpattern     = mcfdata->zeroarcarray;
   SCIP_MCFMODELTYPE modeltype      = mcfdata->modeltype;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int rowcom;
   int rowcomsign;
   SCIP_Bool incompatible;
   SCIP_Real overlap;
   int* overlappingarcs;
   int noverlappingarcs;
   int r;
   int i;

   *score = 0.0;
   *invertcommodity = FALSE;

#ifndef NDEBUG
   for( i = 0; i < narcs; i++ )
      assert(arcpattern[i] == 0);
#endif

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &overlappingarcs, narcs) );

   r = SCIProwGetLPPos(row);
   assert(0 <= r && r < SCIPgetNLPRows(scip));
   assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);
   rowcom = rowcommodity[r];
   assert(0 <= rowcom && rowcom < mcfdata->ncommodities);
   rowcomsign = commoditysigns[rowcom];
   assert(-1 <= rowcomsign && rowcomsign <= +1);

   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   incompatible = FALSE;
   noverlappingarcs = 0;
   for( i = 0; i < rowlen; i++ )
   {
      int c;
      int arcid;
      int valsign;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      arcid = colarcid[c];
      if( arcid == -1 )
         continue;
      assert(arcid < narcs);

      /* get the sign of the coefficient in the flow conservation constraint */
      valsign = (rowvals[i] > 0.0 ? +1 : -1);
      if( (flowrowsigns[r] & LHSASSIGNED) != 0 )
         valsign *= -1;

      /* check if this arc is also member of the base row */
      if( basearcpattern[arcid] != 0 )
      {
         /* check if the sign of the arc matches in the directed case */
         if( modeltype == SCIP_MCFMODELTYPE_DIRECTED )
         {
            int validcomsign;

            if( valsign * basearcpattern[arcid] > 0 )
               validcomsign = +1;
            else
               validcomsign = -1;

            if( rowcomsign == 0 )
            {
               /* the first entry decices whether we have to invert the commodity */
               rowcomsign = validcomsign;
            }
            else if( rowcomsign != validcomsign )
            {
               /* the signs do not fit: this is incompatible */
               incompatible = TRUE;
               break;
            }
         }
         else
         {
            /* in the undirected case, we ignore the sign of the coefficient */
            valsign = +1;
         }

         /* store overlapping arc pattern */
         if( arcpattern[arcid] == 0 )
         {
            overlappingarcs[noverlappingarcs] = arcid;
            noverlappingarcs++;
         }
         arcpattern[arcid] += valsign;
      }
   }

   /* calculate the weighted overlap and reset the zeroarcarray */
   overlap = 0.0;
   for( i = 0; i < noverlappingarcs; i++ )
   {
      SCIP_Real basenum;
      SCIP_Real arcnum;
      int arcid;

      arcid = overlappingarcs[i];
      assert(0 <= arcid && arcid < narcs);
      assert(modeltype == SCIP_MCFMODELTYPE_UNDIRECTED || rowcomsign * basearcpattern[arcid] * arcpattern[arcid] > 0);

      basenum = ABS(basearcpattern[arcid]);
      arcnum = ABS(arcpattern[arcid]);
      assert(basenum != 0.0);
      assert(arcnum != 0.0);

      if( basenum > arcnum )
         overlap += arcnum/basenum;
      else
         overlap += basenum/arcnum;

      arcpattern[arcid] = 0;
   }

   /* calculate the score: maximize overlap and use minimal number of non-overlapping entries as tie breaker */
   if( !incompatible && overlap > 0.0 )
   {
      assert(overlap <= rowlen);
      assert(overlap <= baserowlen);
      assert(noverlappingarcs >= 1);
      *score = overlap - (rowlen + baserowlen - 2.0*overlap)/(narcs+1.0);
      if( noverlappingarcs >= 2 ) /* only one overlapping arc is very dangerous, since this can also be the other end node of the arc */
         *score += 1000.0;
      *score = MAX(*score, 1e-6); /* score may get negative due to many columns in row without an arcid */
      *invertcommodity = (rowcomsign == -1);
   }

   SCIPdebugMessage(" -> node similarity: row <%s>: incompatible=%d overlap=%g rowlen=%d baserowlen=%d score=%g\n",
                    SCIProwGetName(row), incompatible, overlap, rowlen, baserowlen, *score);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &overlappingarcs);
   
   return SCIP_OKAY;
}


/** assigns node ids to flow conservation constraints */
static
SCIP_RETCODE extractNodes(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char*    flowrowsigns   = mcfdata->flowrowsigns;
   int               ncommodities   = mcfdata->ncommodities;
   int*              commoditysigns = mcfdata->commoditysigns;
   int               narcs          = mcfdata->narcs;
   int*              flowcands      = mcfdata->flowcands;
   int               nflowcands     = mcfdata->nflowcands;
   int*              rowcommodity   = mcfdata->rowcommodity;
   int*              colarcid       = mcfdata->colarcid;
   int*              newcols        = mcfdata->newcols;
   SCIP_MCFMODELTYPE modeltype      = mcfdata->modeltype;
   int*              rownodeid;
   SCIP_Bool*        colisincident;
   SCIP_Bool*        rowprocessed;

   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;

   int* arcpattern;
   SCIP_ROW** bestflowrows;
   SCIP_Real* bestscores;
   SCIP_Bool* bestinverted;
   int r;
   int c;
   int n;

   assert(mcfdata->nnodes == 0);
   assert(modeltype != SCIP_MCFMODELTYPE_AUTO);

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->rownodeid, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->colisincident, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->zeroarcarray, narcs) );
   BMSclearMemoryArray(mcfdata->zeroarcarray, narcs);
   rownodeid = mcfdata->rownodeid;
   colisincident = mcfdata->colisincident;

   /* allocate temporary local memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &arcpattern, narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bestflowrows, ncommodities) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bestscores, ncommodities) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bestinverted, ncommodities) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &rowprocessed, nrows) );

   /* initialize temporary memory */
   for( r = 0; r < nrows; r++ )
      rownodeid[r] = -1;
   for( c = 0; c < ncols; c++ )
      colisincident[c] = FALSE;

   /* process all flow conservation constraints that have been used */
   for( n = 0; n < nflowcands; n++ )
   {
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      int rowlen;
      int rowscale;
      int basecommodity;
      int i;

      r = flowcands[n];
      assert(0 <= r && r < nrows);

      /* ignore rows that are not used as flow conservation constraint */
      basecommodity = rowcommodity[r];
      if( basecommodity == -1 )
         continue;
      assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);
      assert(mcfdata->rowarcid[r] == -1);

      /* skip rows that are already assigned to a node */
      if( rownodeid[r] >= 0 )
         continue;

      /* assign row to new node id */
      SCIPdebugMessage("assigning row %d <%s> of commodity %d to node %d [score: %g]\n",
                       r, SCIProwGetName(rows[r]), basecommodity, mcfdata->nnodes, mcfdata->flowrowscores[r]);
      rownodeid[r] = mcfdata->nnodes;

      /* get the arc pattern of the flow row */
      BMSclearMemoryArray(arcpattern, narcs);
      rowcols = SCIProwGetCols(rows[r]);
      rowvals = SCIProwGetVals(rows[r]);
      rowlen = SCIProwGetNLPNonz(rows[r]);
      if( (flowrowsigns[r] & RHSASSIGNED) != 0 )
         rowscale = +1;
      else
         rowscale = -1;
      if( commoditysigns[basecommodity] == -1 )
         rowscale *= -1;

      for( i = 0; i < rowlen; i++ )
      {
         int arcid;

         c = SCIPcolGetLPPos(rowcols[i]);
         assert(0 <= c && c < ncols);
         arcid = colarcid[c];
         if( arcid >= 0 )
         {
            /* due to presolving we may have multiple flow variables of the same arc in the row */
            if( modeltype == SCIP_MCFMODELTYPE_UNDIRECTED || rowscale * rowvals[i] > 0.0 )
            {
               assert(arcpattern[arcid] >= 0);
               arcpattern[arcid]++;
            }
            else
            {
               assert(arcpattern[arcid] <= 0);
               arcpattern[arcid]--;
            }
         }
      }

      /* initialize arrays to store best flow rows */
      for( i = 0; i < ncommodities; i++ )
      {
         bestflowrows[i] = NULL;
         bestscores[i] = 0.0;
         bestinverted[i] = FALSE;
      }

      /* collect columns that are member of incident arc capacity constraints */
      collectIncidentFlowCols(scip, mcfdata, rows[r], basecommodity);

      /* initialize rowprocessed array */
      BMSclearMemoryArray(rowprocessed, nrows);

      /* identify flow conservation constraints in other commodities that match this node;
       * search for flow rows in the column vectors of the indicent columns
       */
      for( i = 0; i < mcfdata->nnewcols; i++ )
      {
         SCIP_ROW** colrows;
         int collen;
         int j;

         c = newcols[i];
         assert(0 <= c && c < ncols);
         assert(mcfdata->colcommodity[c] >= 0);
         assert(mcfdata->colcommodity[c] != basecommodity);

         /* clean up the marker array */
         assert(colisincident[c]);
         colisincident[c] = FALSE;

         /* scan column vector for flow conservation constraints */
         colrows = SCIPcolGetRows(cols[c]);
         collen = SCIPcolGetNLPNonz(cols[c]);
         for( j = 0; j < collen; j++ )
         {
            int colr;
            int rowcom;
            SCIP_Real score;
            SCIP_Bool invertcommodity;

            colr = SCIProwGetLPPos(colrows[j]);
            assert(0 <= colr && colr < nrows);

            /* ignore rows that have already been processed */
            if( rowprocessed[colr] )
               continue;
            rowprocessed[colr] = TRUE;

            /* ignore rows that are not flow conservation constraints in the network */
            rowcom = rowcommodity[colr];
            assert(rowcom != basecommodity);
            if( rowcom == -1 )
               continue;
            assert(rowcom == mcfdata->colcommodity[c]);
            assert((flowrowsigns[colr] & (LHSASSIGNED | RHSASSIGNED)) != 0);
            assert(mcfdata->rowarcid[colr] == -1);

            /* ignore rows that are already assigned to a node */
            if( rownodeid[colr] >= 0 )
               continue;

            /* compare row against arc pattern and calculate score */
            SCIP_CALL( getNodeSilimarityScore(scip, mcfdata, rowlen, arcpattern, colrows[j], &score, &invertcommodity) );
            if( score > bestscores[rowcom] )
            {
               bestflowrows[rowcom] = colrows[j];
               bestscores[rowcom] = score;
               bestinverted[rowcom] = invertcommodity;
            }
         }
      }

      /* for each commodity, pick the best flow conservation constraint to define this node */
      for( i = 0; i < ncommodities; i++ )
      {
         int comr;

         if( bestflowrows[i] == NULL )
            continue;

         comr = SCIProwGetLPPos(bestflowrows[i]);
         assert(0 <= comr && comr < nrows);
         assert(rowcommodity[comr] == i);
         assert((flowrowsigns[comr] & (LHSASSIGNED | RHSASSIGNED)) != 0);
         assert(rownodeid[comr] == -1);

         /* assign flow row to current node */
         SCIPdebugMessage(" -> assigning row %d <%s> of commodity %d to node %d [invert:%d]\n",
                          comr, SCIProwGetName(rows[comr]), i, mcfdata->nnodes, bestinverted[i]);
         rownodeid[comr] = mcfdata->nnodes;

         /* fix the direction of the arcs of the commodity */
         if( bestinverted[i] )
         {
            assert(commoditysigns[i] != +1);
            commoditysigns[i] = -1;
         }
         else
         {
            assert(commoditysigns[i] != -1);
            commoditysigns[i] = +1;
         }
      }
      
      /* increase number of nodes */
      mcfdata->nnodes++;
   }

   /* free local temporary memory */
   SCIPfreeMemoryArray(scip, &rowprocessed);
   SCIPfreeMemoryArray(scip, &bestinverted);
   SCIPfreeMemoryArray(scip, &bestscores);
   SCIPfreeMemoryArray(scip, &bestflowrows);
   SCIPfreeMemoryArray(scip, &arcpattern);

   return SCIP_OKAY;
}

/* if there are still undecided commodity signs, fix them to +1 */
static
void fixCommoditySigns(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int* commoditysigns = mcfdata->commoditysigns;
   int k;

   for( k = 0; k < mcfdata->ncommodities; k++ )
   {
      if( commoditysigns[k] == 0 )
         commoditysigns[k] = +1;
   }
}

/** cleans up the network: gets rid of commodities without arcs or with at most one node */
static
SCIP_RETCODE cleanupNetwork(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int*       flowcands      = mcfdata->flowcands;
   int        nflowcands     = mcfdata->nflowcands;
   int*       colcommodity   = mcfdata->colcommodity;
   int*       rowcommodity   = mcfdata->rowcommodity;
   int*       colarcid       = mcfdata->colarcid;
   int*       rowarcid       = mcfdata->rowarcid;
   int*       rownodeid      = mcfdata->rownodeid;
   int        ncommodities   = mcfdata->ncommodities;
   int*       commoditysigns = mcfdata->commoditysigns;
   int        narcs          = mcfdata->narcs;
   int        nnodes         = mcfdata->nnodes;
   SCIP_ROW** capacityrows   = mcfdata->capacityrows;

   SCIP_ROW** rows;
   int nrows;
   int ncols;

   int* nnodespercom;
   int* narcspercom;
   SCIP_Bool* arcisincom;
   int* perm;
   int permsize;
   int maxnnodes;
   int nnodesthreshold;
   int newncommodities;

   int i;
   int a;
   int k;

   printf/*???????????????SCIPdebugMessage*/("network before cleanup has %d nodes, %d arcs, and %d commodities\n", nnodes, narcs, ncommodities);

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   ncols = SCIPgetNLPCols(scip);

   /* allocate temporary memory */
   permsize = ncommodities;
   permsize = MAX(permsize, narcs);
   permsize = MAX(permsize, nnodes);
   SCIP_CALL( SCIPallocBufferArray(scip, &nnodespercom, ncommodities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &narcspercom, ncommodities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &arcisincom, ncommodities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, permsize) );
   BMSclearMemoryArray(nnodespercom, ncommodities);
   BMSclearMemoryArray(narcspercom, ncommodities);

   /**@todo remove nodes without any incoming and outgoing arcs */

   /* count the number of nodes in each commodity */
   for( i = 0; i < nflowcands; i++ )
   {
      int r;

      r = flowcands[i];
      assert(0 <= r && r < nrows);
      assert((rownodeid[r] >= 0) == (rowcommodity[r] >= 0));
      if( rowcommodity[r] >= 0 )
      {
         assert(rowcommodity[r] < ncommodities);
         nnodespercom[rowcommodity[r]]++;
      }
   }

   /* count the number of arcs in each commodity */
   for( a = 0; a < narcs; a++ )
   {
      SCIP_COL** rowcols;
      int rowlen;
      int r;
      int j;

      r = SCIProwGetLPPos(capacityrows[a]);
      assert(0 <= r && r < nrows);
      assert(rowarcid[r] == a);

      /* identify commodities which are touched by this arc capacity constraint */
      BMSclearMemoryArray(arcisincom, ncommodities);
      rowcols = SCIProwGetCols(rows[r]);
      rowlen = SCIProwGetNLPNonz(rows[r]);
      for( j = 0; j < rowlen; j++ )
      {
         int c;

         c = SCIPcolGetLPPos(rowcols[j]);
         assert(0 <= c && c < ncols);
         if( colcommodity[c] >= 0 && colarcid[c] == a )
         {
            assert(colcommodity[c] < ncommodities);
            arcisincom[colcommodity[c]] = TRUE;
         }
      }

      /* increase arc counters of touched commodities */
      for( k = 0; k < ncommodities; k++ )
      {
         if( arcisincom[k] )
            narcspercom[k]++;
      }
   }

   /* calculate maximal number of nodes per commodity */
   maxnnodes = 0;
   for( k = 0; k < ncommodities; k++ )
      maxnnodes = MAX(maxnnodes, nnodespercom[k]);

   /* we want to keep only commodities that have at least a certain size relative
    * to the largest commodity
    */
   nnodesthreshold = 0.2 * maxnnodes;
   nnodesthreshold = MAX(nnodesthreshold, 1);
   nnodesthreshold = MIN(nnodesthreshold, 10);
   SCIPdebugMessage(" -> node threshold: %d\n", nnodesthreshold);

   /* discard trivial commodities */
   newncommodities = 0;
   for( k = 0; k < ncommodities; k++ )
   {
      printf/*?????????????SCIPdebugMessage*/(" -> commodity %d: %d nodes, %d arcs\n", k, nnodespercom[k], narcspercom[k]);

      /* only keep commodities of a certain size that have at least one arc */
      if( nnodespercom[k] >= nnodesthreshold && narcspercom[k] >= 1 )
      {
         assert(newncommodities <= k);
         perm[k] = newncommodities;
         commoditysigns[newncommodities] = commoditysigns[k];
         newncommodities++;
      }
      else
         perm[k] = -1;
   }

   if( newncommodities < ncommodities )
   {
      SCIP_Bool* arcisused;
      SCIP_Bool* nodeisused;
      int newnarcs;
      int newnnodes;
      int c;
      int v;

      printf/*?????????????SCIPdebugMessage*/(" -> discarding %d of %d commodities\n", ncommodities - newncommodities, ncommodities);

      SCIP_CALL( SCIPallocBufferArray(scip, &arcisused, narcs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nodeisused, nnodes) );

      /* update data structures to new commodity ids */
      BMSclearMemoryArray(arcisused, narcs);
      BMSclearMemoryArray(nodeisused, nnodes);
      for( c = 0; c < ncols; c++ )
      {
         if( colcommodity[c] >= 0 )
         {
            assert(-1 <= colarcid[c] && colarcid[c] < narcs);
            assert(colcommodity[c] < mcfdata->ncommodities);
            colcommodity[c] = perm[colcommodity[c]];
            assert(colcommodity[c] < newncommodities);
            if( colcommodity[c] == -1 )
            {
               /* we are lazy and do not update plusflow and minusflow */
               colarcid[c] = -1;
            }
            else if( colarcid[c] >= 0 )
               arcisused[colarcid[c]] = TRUE;
         }
      }
      for( i = 0; i < nflowcands; i++ )
      {
         int r;

         r = flowcands[i];
         assert(0 <= r && r < nrows);
         assert((rownodeid[r] >= 0) == (rowcommodity[r] >= 0));
         if( rowcommodity[r] >= 0 )
         {
            assert(0 <= rownodeid[r] && rownodeid[r] < nnodes);
            assert(rowcommodity[r] < mcfdata->ncommodities);
            rowcommodity[r] = perm[rowcommodity[r]];
            assert(rowcommodity[r] < newncommodities);
            if( rowcommodity[r] == -1 )
            {
               /* we are lazy and do not update flowrowsigns */
               rownodeid[r] = -1;
            }
            else
               nodeisused[rownodeid[r]] = TRUE;
         }
      }
      mcfdata->ncommodities = newncommodities;
      ncommodities = newncommodities;

      /* discard unused arcs */
      newnarcs = 0;
      for( a = 0; a < narcs; a++ )
      {
         int r;

         if( arcisused[a] )
         {
            assert(newnarcs <= a);
            perm[a] = newnarcs;
            capacityrows[newnarcs] = capacityrows[a];
            newnarcs++;
         }
         else
         {
            /* we are lazy and do not update capacityrowsigns */
            perm[a] = -1;
         }
         r = SCIProwGetLPPos(capacityrows[a]);
         assert(0 <= r && r < nrows);
         assert(rowarcid[r] == a);
         rowarcid[r] = perm[a];
      }

      /* update remaining data structures to new arc ids */
      if( newnarcs < narcs )
      {
         printf/*?????????????SCIPdebugMessage*/(" -> discarding %d of %d arcs\n", narcs - newnarcs, narcs);

         for( c = 0; c < ncols; c++ )
         {
            if( colarcid[c] >= 0 )
            {
               colarcid[c] = perm[colarcid[c]];
               assert(colarcid[c] >= 0); /* otherwise colarcid[c] was set to -1 in the colcommodity update */
            }
         }
         mcfdata->narcs = newnarcs;
         narcs = newnarcs;
      }
#ifndef NDEBUG
      for( a = 0; a < narcs; a++ )
      {
         int r;
         r = SCIProwGetLPPos(capacityrows[a]);
         assert(0 <= r && r < nrows);
         assert(rowarcid[r] == a);
      }
#endif

      /* discard unused nodes */
      newnnodes = 0;
      for( v = 0; v < nnodes; v++ )
      {
         if( nodeisused[v] )
         {
            assert(newnnodes <= v);
            perm[v] = newnnodes;
            newnnodes++;
         }
         else
            perm[v] = -1;
      }

      /* update data structures to new node ids */
      if( newnnodes < nnodes )
      {
         printf/*?????????????SCIPdebugMessage*/(" -> discarding %d of %d nodes\n", nnodes - newnnodes, nnodes);

         for( i = 0; i < nflowcands; i++ )
         {
            int r;

            r = flowcands[i];
            assert(0 <= r && r < nrows);
            assert((rownodeid[r] >= 0) == (rowcommodity[r] >= 0));
            if( rowcommodity[r] >= 0 )
            {
               assert(rowcommodity[r] < ncommodities);
               rownodeid[r] = perm[rownodeid[r]];
               assert(rownodeid[r] >= 0); /* otherwise we would have deleted the commodity in the rowcommodity update above */
            }
         }
         mcfdata->nnodes = newnnodes;
         nnodes = newnnodes;
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &nodeisused);
      SCIPfreeBufferArray(scip, &arcisused);
   }
   
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &arcisincom);
   SCIPfreeBufferArray(scip, &narcspercom);
   SCIPfreeBufferArray(scip, &nnodespercom);

   printf/*???????????????SCIPdebugMessage*/("network after cleanup has %d nodes, %d arcs, and %d commodities\n", nnodes, narcs, ncommodities);

   return SCIP_OKAY;
}

/** identifies the (at most) two nodes which contain the given flow variable */
static
void getIncidentNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_COL*             col,                /**< flow column */
   int*                  sourcenode,         /**< pointer to store the source node of the flow column */
   int*                  targetnode          /**< pointer to store the target node of the flow column */
   )
{
   unsigned char*    flowrowsigns     = mcfdata->flowrowsigns;
   int*              commoditysigns   = mcfdata->commoditysigns;
   int*              rowcommodity     = mcfdata->rowcommodity;
   int*              rownodeid        = mcfdata->rownodeid;

   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   int collen;
   int i;

   assert(sourcenode != NULL);
   assert(targetnode != NULL);

   *sourcenode = -1;
   *targetnode = -1;

   /* search for flow conservation rows in the column vector */
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   collen = SCIPcolGetNLPNonz(col);
   for( i = 0; i < collen; i++ )
   {
      int r;

      r = SCIProwGetLPPos(colrows[i]);
      assert(0 <= r && r < SCIPgetNLPRows(scip));

      if( rownodeid[r] >= 0 )
      {
         int v;
         int k;
         int scale;

         v = rownodeid[r];
         k = rowcommodity[r];
         assert(0 <= v && v < mcfdata->nnodes);
         assert(0 <= k && k < mcfdata->ncommodities);
         assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);

         /* check whether the flow row is inverted */
         scale = +1;
         if( (flowrowsigns[r] & LHSASSIGNED) != 0 )
            scale *= -1;
         if( commoditysigns[k] == -1 )
            scale *= -1;

         /* decide whether this node is source or target */
         if( scale * colvals[i] > 0.0 )
         {
            assert(*sourcenode == -1);
            *sourcenode = v;
            if( *targetnode >= 0 )
               break;
         }
         else
         {
            assert(*targetnode == -1);
            *targetnode = v;
            if( *sourcenode >= 0 )
               break;
         }
      }
   }
}

/** for each arc identifies a source and target node */
static
SCIP_RETCODE identifySourcesTargets(
   SCIP*                 scip,               /**< SCIP data structure */
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int*              colarcid         = mcfdata->colarcid;
   int               narcs            = mcfdata->narcs;
   int               nnodes           = mcfdata->nnodes;
   SCIP_ROW**        capacityrows     = mcfdata->capacityrows;
   SCIP_MCFMODELTYPE modeltype        = mcfdata->modeltype;
   int*              arcsources;
   int*              arctargets;
   int*              firstoutarcs;
   int*              firstinarcs;
   int*              nextoutarcs;
   int*              nextinarcs;

   int *sourcenodecnt;
   int *targetnodecnt;
   int *touchednodes;
   int ntouchednodes;

   int v;
   int a;

   /* allocate memory in mcfdata */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->arcsources, narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->arctargets, narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->firstoutarcs, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->firstinarcs, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->nextoutarcs, narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->nextinarcs, narcs) );
   arcsources   = mcfdata->arcsources;
   arctargets   = mcfdata->arctargets;
   firstoutarcs = mcfdata->firstoutarcs;
   firstinarcs  = mcfdata->firstinarcs;
   nextoutarcs  = mcfdata->nextoutarcs;
   nextinarcs   = mcfdata->nextinarcs;

   /* initialize adjacency lists */
   for( v = 0; v < nnodes; v++ )
   {
      firstoutarcs[v] = -1;
      firstinarcs[v] = -1;
   }
   for( a = 0; a < narcs; a++ )
   {
      nextoutarcs[a] = -1;
      nextinarcs[a] = -1;
   }

   /* allocate temporary memory for source and target node identification */
   SCIP_CALL( SCIPallocBufferArray(scip, &sourcenodecnt, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &targetnodecnt, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &touchednodes, nnodes) );
   BMSclearMemoryArray(sourcenodecnt, nnodes);
   BMSclearMemoryArray(targetnodecnt, nnodes);

   mcfdata->ninconsistencies = 0;

   /* search for source and target nodes */
   for( a = 0; a < narcs; a++ )
   {
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      int rowlen;
      int bestsourcev;
      int besttargetv;
      int bestsourcecnt;
      int besttargetcnt;
      int totalnodecnt;
      int r;
      int i;

      r = SCIProwGetLPPos(capacityrows[a]);
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert((mcfdata->capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);
      assert(mcfdata->rowarcid[r] == a);

#ifndef NDEBUG
      for( i = 0; i < nnodes; i++ )
      {
         assert(sourcenodecnt[i] == 0);
         assert(targetnodecnt[i] == 0);
      }
#endif

      /* check the flow variables of the capacity row for flow conservation constraints */
      rowcols = SCIProwGetCols(capacityrows[a]);
      rowvals = SCIProwGetVals(capacityrows[a]);
      rowlen = SCIProwGetNLPNonz(capacityrows[a]);
      ntouchednodes = 0;
      totalnodecnt = 0;
      for( i = 0; i < rowlen; i++ )
      {
         int c;

         c = SCIPcolGetLPPos(rowcols[i]);
         assert(0 <= c && c < SCIPgetNLPCols(scip));
         if( colarcid[c] >= 0 )
         {
            int sourcev;
            int targetv;

            /* identify the (at most) two nodes which contain this flow variable */
            getIncidentNodes(scip, mcfdata, rowcols[i], &sourcev, &targetv);

            /* count the nodes */
            if( sourcev >= 0 )
            {
               if( sourcenodecnt[sourcev] == 0 && targetnodecnt[sourcev] == 0 )
               {
                  touchednodes[ntouchednodes] = sourcev;
                  ntouchednodes++;
               }
               sourcenodecnt[sourcev]++;
               totalnodecnt++;
            }
            if( targetv >= 0 )
            {
               if( sourcenodecnt[targetv] == 0 && targetnodecnt[targetv] == 0 )
               {
                  touchednodes[ntouchednodes] = targetv;
                  ntouchednodes++;
               }
               targetnodecnt[targetv]++;
               totalnodecnt++;
            }
         }
      }

      /* perform a majority vote on source and target node */
      bestsourcev = -1;
      besttargetv = -1;
      bestsourcecnt = 0;
      besttargetcnt = 0;
      for( i = 0; i < ntouchednodes; i++ )
      {
         v = touchednodes[i];
         assert(0 <= v && v < nnodes);
            
         if( modeltype == SCIP_MCFMODELTYPE_DIRECTED )
         {
            /* in the directed model, we distinguish between source and target */
            if( sourcenodecnt[v] >= targetnodecnt[v] )
            {
               if( sourcenodecnt[v] > bestsourcecnt )
               {
                  bestsourcev = v;
                  bestsourcecnt = sourcenodecnt[v];
               }
            }
            else
            {
               if( targetnodecnt[v] > besttargetcnt )
               {
                  besttargetv = v;
                  besttargetcnt = targetnodecnt[v];
               }
            }
         }
         else
         {
            int nodecnt = sourcenodecnt[v] + targetnodecnt[v];

            /* in the undirected model, we use source for the maximum and target for the second largest number of total hits */
            if( nodecnt > bestsourcecnt )
            {
               besttargetv = bestsourcev;
               besttargetcnt = bestsourcecnt;
               bestsourcev = v;
               bestsourcecnt = nodecnt;
            }
            else if( nodecnt > besttargetcnt )
            {
               besttargetv = v;
               besttargetcnt = nodecnt;
            }
         }

         /* clear the nodecnt arrays */
         sourcenodecnt[v] = 0;
         targetnodecnt[v] = 0;
      }

      /* assign the incident nodes */
      assert(bestsourcev == -1 || bestsourcev != besttargetv);
      arcsources[a] = bestsourcev;
      arctargets[a] = besttargetv;
      SCIPdebugMessage("arc %d: %d -> %d (sourcecnt=%d, targetcnt=%d, %d inconsistencies)\n",
                       a, bestsourcev, besttargetv, bestsourcecnt, besttargetcnt, totalnodecnt - bestsourcecnt - besttargetcnt);

      /* update adjacency lists */
      if( bestsourcev != -1 )
      {
         nextoutarcs[a] = firstoutarcs[bestsourcev];
         firstoutarcs[bestsourcev] = a;
      }
      if( besttargetv != -1 )
      {
         nextinarcs[a] = firstinarcs[besttargetv];
         firstinarcs[besttargetv] = a;
      }

      /* update the number of inconsistencies */
      assert(totalnodecnt >= bestsourcecnt + besttargetcnt);
      mcfdata->ninconsistencies += totalnodecnt - bestsourcecnt - besttargetcnt;
   }

   printf/*???????????????SCIPdebugMessage*/("extracted network has %d inconsistencies\n", mcfdata->ninconsistencies);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &touchednodes);
   SCIPfreeBufferArray(scip, &targetnodecnt);
   SCIPfreeBufferArray(scip, &sourcenodecnt);

   return SCIP_OKAY;
}

#define UNKNOWN 0  /**< node has not yet been seen */
#define ONSTACK 1  /**< node is currently on the processing stack */
#define VISITED 2  /**< node has been visited and assigned to some component */

/** returns lists of nodes and arcs in the connected component of the given startv */
static
SCIP_RETCODE identifyComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   int*                  nodevisited,        /**< array to mark visited nodes */
   int                   startv,             /**< node for which the connected component should be generated */
   int*                  compnodes,          /**< array to store node ids of the component */
   int*                  ncompnodes,         /**< pointer to store the number of nodes in the component */
   int*                  comparcs,           /**< array to store arc ids of the component */
   int*                  ncomparcs           /**< pointer to store the number of arcs in the component */
   )
{
   int* arcsources   = mcfdata->arcsources;
   int* arctargets   = mcfdata->arctargets;
   int* firstoutarcs = mcfdata->firstoutarcs;
   int* firstinarcs  = mcfdata->firstinarcs;
   int* nextoutarcs  = mcfdata->nextoutarcs;
   int* nextinarcs   = mcfdata->nextinarcs;
   int  nnodes       = mcfdata->nnodes;

   int* stacknodes;
   int nstacknodes;

   assert(nodevisited != NULL);
   assert(0 <= startv && startv < nnodes);
   assert(nodevisited[startv] == UNKNOWN);
   assert(compnodes != NULL);
   assert(ncompnodes != NULL);
   assert(comparcs != NULL);
   assert(ncomparcs != NULL);

   *ncompnodes = 0;
   *ncomparcs = 0;

   /* allocate temporary memory for node stack */
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknodes, nnodes) );

   /* put startv on stack */
   stacknodes[0] = startv;
   nstacknodes = 1;
   nodevisited[startv] = ONSTACK;

   /* perform depth-first search */
   while( nstacknodes > 0 )
   {
      int v;
      int a;

      /* pop first element from stack */
      v = stacknodes[nstacknodes-1];
      nstacknodes--;
      assert(0 <= v && v < nnodes);
      assert(nodevisited[v] == ONSTACK);
      nodevisited[v] = VISITED;

      /* put node into component */
      assert(*ncompnodes < nnodes);
      compnodes[*ncompnodes] = v;
      (*ncompnodes)++;

      /* go through the list of outgoing arcs */
      for( a = firstoutarcs[v]; a != -1; a = nextoutarcs[a] )
      {
         int targetv;

         assert(0 <= a && a < mcfdata->narcs);

         targetv = arctargets[a];

         /* check if we have already visited the target node */
         if( targetv != -1 && nodevisited[targetv] == VISITED )
            continue;

         /* put arc to component */
         assert(*ncomparcs < mcfdata->narcs);
         comparcs[*ncomparcs] = a;
         (*ncomparcs)++;

         /* push target node to stack */
         if( targetv != -1 && nodevisited[targetv] == UNKNOWN )
         {
            assert(nstacknodes < nnodes);
            stacknodes[nstacknodes] = targetv;
            nstacknodes++;
            nodevisited[targetv] = ONSTACK;
         }
      }

      /* go through the list of ingoing arcs */
      for( a = firstinarcs[v]; a != -1; a = nextinarcs[a] )
      {
         int sourcev;

         assert(0 <= a && a < mcfdata->narcs);

         sourcev = arcsources[a];

         /* check if we have already seen the source node */
         if( sourcev != -1 && nodevisited[sourcev] == VISITED )
            continue;

         /* put arc to component */
         assert(*ncomparcs < mcfdata->narcs);
         comparcs[*ncomparcs] = a;
         (*ncomparcs)++;

         /* push source node to stack */
         if( sourcev != -1 && nodevisited[sourcev] == UNKNOWN )
         {
            assert(nstacknodes < nnodes);
            stacknodes[nstacknodes] = sourcev;
            nstacknodes++;
            nodevisited[sourcev] = ONSTACK;
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &stacknodes);

   return SCIP_OKAY;
}

/** extracts MCF network structures from the current LP */
static
SCIP_RETCODE mcfnetworkExtract(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_MCFNETWORK***    mcfnetworks,        /**< pointer to store array of MCF network structures */
   int*                  nmcfnetworks        /**< pointer to store number of MCF networks */
   )
{
   MCFDATA mcfdata;

   SCIP_MCFMODELTYPE modeltype = sepadata->modeltype;
   SCIP_MCFFLOWTYPE  flowtype  = sepadata->flowtype;

   int* nodevisited;
   int* compnodeid;
   int* compnodes;
   int* comparcs;
   int ncompnodes;
   int ncomparcs;
   int mcfnetworkssize;
   int v;

   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;

   assert(mcfnetworks != NULL);
   assert(nmcfnetworks != NULL);

   *mcfnetworks = NULL;
   *nmcfnetworks = 0;
   mcfnetworkssize = 0;

   /* Algorithm to identify multi-commodity-flow network with capacity constraints
    *
    * 1. Identify candidate rows for flow conservation constraints in the LP.
    * 2. Sort flow conservation candidates by a ranking on how sure we are that it is indeed a constraint of the desired type.
    * 3. Extract network structure of flow conservation constraints:
    *    (a) Initialize plusflow[c] = minusflow[c] = FALSE for all columns c and other local data.
    *    (b) As long as there are flow conservation candidates left:
    *        (i) Create new commodity and use first flow conservation constraint as newrow.
    *       (ii) Add newrow to commodity, update pluscom/minuscom accordingly.
    *      (iii) For the newly added columns search for an incident flow conservation constraint. Pick the one of highest ranking.
    *       (iv) If found, set newrow to this row and goto (ii).
    *        (v) If only very few flow rows have been used, discard the commodity immediately.
    * 4. Identify candidate rows for capacity constraints in the LP.
    * 5. Sort capacity constraint candidates by a ranking on how sure we are that it is indeed a constraint of the desired type.
    * 6. Identify capacity constraints for the arcs and assign arc ids to columns and capacity constraints.
    * 7. Assign node ids to flow conservation constraints.
    */

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* initialize local extraction data */
   mcfdata.flowrowsigns = NULL;
   mcfdata.flowrowscalars = NULL;
   mcfdata.flowrowscores = NULL;
   mcfdata.capacityrowsigns = NULL;
   mcfdata.capacityrowscores = NULL;
   mcfdata.flowcands = NULL;
   mcfdata.nflowcands = 0;
   mcfdata.capacitycands = NULL;
   mcfdata.ncapacitycands = 0;
   mcfdata.plusflow = NULL;
   mcfdata.minusflow = NULL;
   mcfdata.ncommodities = 0;
   mcfdata.commoditysigns = NULL;
   mcfdata.commoditysignssize = 0;
   mcfdata.colcommodity = NULL;
   mcfdata.rowcommodity = NULL;
   mcfdata.colarcid = NULL;
   mcfdata.rowarcid = NULL;
   mcfdata.rownodeid = NULL;
   mcfdata.arcsources = NULL;
   mcfdata.arctargets = NULL;
   mcfdata.firstoutarcs = NULL;
   mcfdata.firstinarcs = NULL;
   mcfdata.nextoutarcs = NULL;
   mcfdata.nextinarcs = NULL;
   mcfdata.newcols = NULL;
   mcfdata.nnewcols = 0;
   mcfdata.narcs = 0;
   mcfdata.nnodes = 0;
   mcfdata.ninconsistencies = 0;
   mcfdata.capacityrows = NULL;
   mcfdata.capacityrowssize = 0;
   mcfdata.colisincident = NULL;
   mcfdata.zeroarcarray = NULL;
   mcfdata.modeltype = modeltype;
   mcfdata.flowtype = flowtype;

   /* 1. Identify candidate rows for flow conservation constraints in the LP.
    * 2. Sort flow conservation candidates by a ranking on how sure we are that it is indeed a constraint of the desired type.
    */
   SCIP_CALL( extractFlowRows(scip, &mcfdata) );

   if( mcfdata.nflowcands > 0 )
   {
      /* 3. Extract network structure of flow conservation constraints. */
      SCIPdebugMessage("****** extracting flow ******\n");
      SCIP_CALL( extractFlow(scip, &mcfdata) );
#ifdef SCIP_DEBUG
      printCommodities(scip, &mcfdata);
#endif
      
      /* 4. Identify candidate rows for capacity constraints in the LP.
       * 5. Sort capacity constraint candidates by a ranking on how sure we are that it is indeed a constraint of the desired type.
       */
      SCIP_CALL( extractCapacityRows(scip, &mcfdata) );

      if( mcfdata.ncapacitycands > 0 )
      {
         /* 6. Identify capacity constraints for the arcs and assign arc ids to columns and capacity constraints. */
         SCIPdebugMessage("****** extracting capacities ******\n");
         SCIP_CALL( extractCapacities(scip, &mcfdata) );

         /* 7. Assign node ids to flow conservation constraints. */
         SCIPdebugMessage("****** extracting nodes ******\n");
         SCIP_CALL( extractNodes(scip, &mcfdata) );

         /* if there are still undecided commodity signs, fix them to +1 */
         fixCommoditySigns(scip, &mcfdata);

         /**@todo Now we need to assign arcids to flow variables that have not yet been assigned to an arc because
          *       there was no valid capacity constraint.
          *       Go through the list of nodes and generate uncapacitated arcs in the network for the flow variables
          *       that do not yet have an arc assigned, such that the commodities still match.
          */

         /* clean up the network: get rid of commodities without arcs or with at most one node */
         SCIP_CALL( cleanupNetwork(scip, &mcfdata) );

         /* assign source and target nodes to arc */
         SCIP_CALL( identifySourcesTargets(scip, &mcfdata) );

         /**@todo auto-detect type of flow variables */
         if( mcfdata.flowtype == SCIP_MCFFLOWTYPE_AUTO )
            mcfdata.flowtype = SCIP_MCFFLOWTYPE_CONTINUOUS;
      }
   }

#ifdef SCIP_DEBUG
   printCommodities(scip, &mcfdata);
#endif

   /* allocate temporary memory for component finding */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, mcfdata.nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compnodes, mcfdata.nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comparcs, mcfdata.narcs) );
   BMSclearMemoryArray(nodevisited, mcfdata.nnodes);

   /* allocate temporary memory for v -> compv mapping */
   SCIP_CALL( SCIPallocBufferArray(scip, &compnodeid, mcfdata.nnodes) );
   for( v = 0; v < mcfdata.nnodes; v++ )
      compnodeid[v] = -1;

   /* search components and create a network structure for each of them */
   for( v = 0; v < mcfdata.nnodes; v++ )
   {
      /* ignore nodes that have been already assigned to a component */
      assert(nodevisited[v] == UNKNOWN || nodevisited[v] == VISITED);
      if( nodevisited[v] == VISITED )
         continue;

      /* identify nodes and arcs of this component */
      SCIP_CALL( identifyComponent(scip, &mcfdata, nodevisited, v, compnodes, &ncompnodes, comparcs, &ncomparcs) );
      assert(ncompnodes >= 1);
      assert(compnodes[0] == v);
      assert(nodevisited[v] == VISITED);

      /* ignore network component if it is trivial */
      if( ncompnodes >= MINNODES && ncomparcs >= MINARCS )
      {
         /* make sure that we have enough memory for the new network pointer */
         assert(*nmcfnetworks <= mcfnetworkssize);
         if( *nmcfnetworks == mcfnetworkssize )
         {
            mcfnetworkssize = MAX(2*mcfnetworkssize, *nmcfnetworks+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, mcfnetworks, mcfnetworkssize) );
         }
         assert(*nmcfnetworks < mcfnetworkssize);

         /* create network data structure */
         SCIP_CALL( mcfnetworkCreate(scip, &(*mcfnetworks)[*nmcfnetworks]) );
         assert((*mcfnetworks)[*nmcfnetworks] != NULL);
         
         /* fill sparse network structure */
         SCIP_CALL( mcfnetworkFill(scip, (*mcfnetworks)[*nmcfnetworks], &mcfdata, compnodeid, compnodes, ncompnodes, comparcs, ncomparcs) );

         (*nmcfnetworks)++;
      }
      else
      {
         SCIPdebugMessage(" -> discarded component with %d nodes and %d arcs\n", ncompnodes, ncomparcs);
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &compnodeid);
   SCIPfreeBufferArray(scip, &comparcs);
   SCIPfreeBufferArray(scip, &compnodes);
   SCIPfreeBufferArray(scip, &nodevisited);

   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &mcfdata.arcsources);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.arctargets);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.firstoutarcs);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.firstinarcs);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.nextoutarcs);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.nextinarcs);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.zeroarcarray);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.colisincident);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacityrows);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.rownodeid);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.rowarcid);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.colarcid);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.newcols);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.rowcommodity);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.colcommodity);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.commoditysigns);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.minusflow);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.plusflow);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacitycands);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowcands);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacityrowscores);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacityrowsigns);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowrowscores);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowrowscalars);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowrowsigns);

   return SCIP_OKAY;
}


/** comparison method for weighted arcs */
static
SCIP_DECL_SORTPTRCOMP(compArcs)
{
   ARCENTRY* arc1 = (ARCENTRY*)elem1;
   ARCENTRY* arc2 = (ARCENTRY*)elem2;

   if( arc1->weight > arc2->weight )
      return -1;
   else if( arc1->weight < arc2->weight )
      return +1;
   else
      return 0;
}

/** creates a priority queue and fills it with the given arc entries */
static
SCIP_RETCODE arcqueueCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK*      mcfnetwork,         /**< MCF network structure */
   ARCQUEUE**            arcqueue            /**< pointer to arc priority queue */
   )
{
   int a;

   assert(mcfnetwork != NULL);
   assert(arcqueue != NULL);

   SCIP_CALL( SCIPallocMemory(scip, arcqueue) );

   /* calculate weights for arcs */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*arcqueue)->arcentries, mcfnetwork->narcs) );
   for( a = 0; a < mcfnetwork->narcs; a++ )
   {
      SCIP_Real slack;
      SCIP_Real dualsol;
      SCIP_Real scale;

      slack = SCIPgetRowFeasibility(scip, mcfnetwork->arccapacityrows[a]);
      slack = MAX(slack, 0.0); /* can only be negative due to numerics */
      dualsol = SCIProwGetDualsol(mcfnetwork->arccapacityrows[a]);
      scale = ABS(mcfnetwork->arccapacityscales[a]);

      /* shrink rows with large slack first, and from the tight rows, shrink the ones with small dual solution value first */
      (*arcqueue)->arcentries[a].arcid = a;
      (*arcqueue)->arcentries[a].weight = scale * slack - ABS(dualsol)/scale;
      SCIPdebugMessage("arc %2d <%s>: slack=%g, dualsol=%g -> weight=%g\n", a, SCIProwGetName(mcfnetwork->arccapacityrows[a]),
                       slack, dualsol, (*arcqueue)->arcentries[a].weight);
   }

   /* create priority queue */
   SCIP_CALL( SCIPpqueueCreate(&(*arcqueue)->pqueue, mcfnetwork->narcs, 2.0, compArcs) );

   /* fill priority queue with arc data */
   for( a = 0; a < mcfnetwork->narcs; a++ )
   {
      SCIP_CALL( SCIPpqueueInsert((*arcqueue)->pqueue, (void*)&(*arcqueue)->arcentries[a]) );
   }

   return SCIP_OKAY;
}

/** frees memory of an arc queue */
static
void arcqueueFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ARCQUEUE**            arcqueue            /**< pointer to arc priority queue */
   )
{
   assert(arcqueue != NULL);
   assert(*arcqueue != NULL);

   SCIPpqueueFree(&(*arcqueue)->pqueue);
   SCIPfreeMemoryArray(scip, &(*arcqueue)->arcentries);
   SCIPfreeMemory(scip, arcqueue);
}

/** returns whether there are any arcs left on the queue */
static
SCIP_Bool arcqueueIsEmpty(
   ARCQUEUE*             arcqueue            /**< arc priority queue */
   )
{
   assert(arcqueue != NULL);

   return (SCIPpqueueFirst(arcqueue->pqueue) == NULL);
}

/** removes the top element from the arc priority queue and returns the arcid */
static
int arcqueueRemove(
   ARCQUEUE*             arcqueue            /**< arc priority queue */
   )
{
   ARCENTRY* arcentry;

   assert(arcqueue != NULL);

   arcentry = (ARCENTRY*)SCIPpqueueRemove(arcqueue->pqueue);
   if( arcentry != NULL )
      return arcentry->arcid;
   else
      return -1;
}

/** returns the representative node in the cluster of the given node */
static
int nodepartitionGetRepresentative(
   NODEPARTITION*        nodepartition,      /**< node partition data structure */
   int                   v                   /**< node id to get representative for */
   )
{
   int* representatives;

   assert(nodepartition != NULL);

   /* we apply a union find algorithm */
   representatives = nodepartition->representatives;
   while( v != representatives[v] )
   {
      representatives[v] = representatives[representatives[v]];
      v = representatives[v];
   }

   return v;
}

/** joins two clusters given by their representative nodes */
static
void nodepartitionJoin(
   NODEPARTITION*        nodepartition,      /**< node partition data structure */
   int                   rep1,               /**< representative of first cluster */
   int                   rep2                /**< representative of second cluster */
   )
{
   assert(nodepartition != NULL);
   assert(rep1 != rep2);
   assert(nodepartition->representatives[rep1] == rep1);
   assert(nodepartition->representatives[rep2] == rep2);

   /* make sure that the smaller representative survives
    *  -> node 0 is always a representative
    */
   if( rep1 < rep2 )
      nodepartition->representatives[rep2] = rep1;
   else
      nodepartition->representatives[rep1] = rep2;
}

/** partitions nodes into a small number of clusters */
static
SCIP_RETCODE nodepartitionCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK*      mcfnetwork,         /**< MCF network structure */
   NODEPARTITION**       nodepartition,      /**< pointer to node partition data structure */
   int                   nclusters           /**< number of clusters to generate */
   )
{
   ARCQUEUE* arcqueue;
   int* clustersize;
   int nclustersleft;
   int v;
   int c;
   int pos;

   assert(mcfnetwork != NULL);
   assert(nodepartition != NULL);

   /* allocate and initialize memory */
   SCIP_CALL( SCIPallocMemory(scip, nodepartition) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*nodepartition)->representatives, mcfnetwork->nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*nodepartition)->nodeclusters, mcfnetwork->nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*nodepartition)->clusternodes, mcfnetwork->nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*nodepartition)->clusterbegin, nclusters+1) );
   (*nodepartition)->nclusters = 0;

   /* we start with each node being in its own cluster */
   for( v = 0; v < mcfnetwork->nnodes; v++ )
      (*nodepartition)->representatives[v] = v;

   /* create priority queue for arcs */
   SCIP_CALL( arcqueueCreate(scip, mcfnetwork, &arcqueue) );

   /* loop over arcs in order of their weights */
   nclustersleft = mcfnetwork->nnodes;
   while( !arcqueueIsEmpty(arcqueue) && nclustersleft > nclusters )
   {
      int a;
      int sourcev;
      int targetv;
      int sourcerep;
      int targetrep;
      
      /* get the next arc */
      a = arcqueueRemove(arcqueue);
      assert(0 <= a && a < mcfnetwork->narcs);
      assert(arcqueue->arcentries[a].arcid == a);

      /* get the source and target node of the arc */
      sourcev = mcfnetwork->arcsources[a];
      targetv = mcfnetwork->arctargets[a];
      assert(sourcev < mcfnetwork->nnodes);
      assert(targetv < mcfnetwork->nnodes);
      if( sourcev == -1 || targetv == -1 )
         continue;

      /* identify the representatives of the two nodes */
      sourcerep = nodepartitionGetRepresentative(*nodepartition, sourcev);
      targetrep = nodepartitionGetRepresentative(*nodepartition, targetv);
      assert(0 <= sourcerep && sourcerep < mcfnetwork->nnodes);
      assert(0 <= targetrep && targetrep < mcfnetwork->nnodes);
      
      /* there is nothing to do if the two nodes are already in the same cluster */
      if( sourcerep == targetrep )
         continue;

      /* shrink the arc by joining the two clusters */
      SCIPdebugMessage("shrinking arc %d <%s> with weight %g (s=%g y=%g): join representatives %d and %d\n",
                       a, SCIProwGetName(mcfnetwork->arccapacityrows[a]), arcqueue->arcentries[a].weight,
                       SCIPgetRowFeasibility(scip, mcfnetwork->arccapacityrows[a]),
                       SCIProwGetDualsol(mcfnetwork->arccapacityrows[a]), sourcerep, targetrep);
      nodepartitionJoin(*nodepartition, sourcerep, targetrep);
      nclustersleft--;
   }

   /* node 0 must be a representative due to our join procedure */
   assert((*nodepartition)->representatives[0] == 0);

   /* if there have been too few arcs to shrink the graph to the required number of clusters, join clusters with first cluster
    * to create a larger disconnected cluster
    */
   if( nclustersleft > nclusters )
   {
      for( v = 1; v < mcfnetwork->nnodes && nclustersleft > nclusters; v++ )
      {
         int rep;

         rep = nodepartitionGetRepresentative(*nodepartition, v);
         if( rep != 0 )
         {
            nodepartitionJoin(*nodepartition, 0, rep);
            nclustersleft--;
         }
      }
   }
   assert(nclustersleft <= nclusters);

   /* extract the clusters */
   SCIP_CALL( SCIPallocBufferArray(scip, &clustersize, nclusters) );
   BMSclearMemoryArray(clustersize, nclusters);
   for( v = 0; v < mcfnetwork->nnodes; v++ )
   {
      int rep;

      /* get cluster of node */
      rep = nodepartitionGetRepresentative(*nodepartition, v);
      assert(rep <= v); /* due to our joining procedure */
      if( rep == v )
      {
         /* node is its own representative: this is a new cluster */
         c = (*nodepartition)->nclusters;
         (*nodepartition)->nclusters++;
      }
      else
         c = (*nodepartition)->nodeclusters[rep];
      assert(0 <= c && c < nclusters);

      /* assign node to cluster */
      (*nodepartition)->nodeclusters[v] = c;
      clustersize[c]++;
   }

   /* fill the clusterbegin array */
   pos = 0;
   for( c = 0; c < (*nodepartition)->nclusters; c++ )
   {
      (*nodepartition)->clusterbegin[c] = pos;
      pos += clustersize[c];
   }
   assert(pos == mcfnetwork->nnodes);
   (*nodepartition)->clusterbegin[(*nodepartition)->nclusters] = mcfnetwork->nnodes;

   /* fill the clusternodes array */
   BMSclearMemoryArray(clustersize, (*nodepartition)->nclusters);
   for( v = 0; v < mcfnetwork->nnodes; v++ )
   {
      c = (*nodepartition)->nodeclusters[v];
      assert(0 <= c && c < (*nodepartition)->nclusters);
      pos = (*nodepartition)->clusterbegin[c] + clustersize[c];
      assert(pos < (*nodepartition)->clusterbegin[c+1]);
      (*nodepartition)->clusternodes[pos] = v;
      clustersize[c]++;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &clustersize);

   /* free arc queue */
   arcqueueFree(scip, &arcqueue);

   return SCIP_OKAY;
}

/** frees node partition data */
static
void nodepartitionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   NODEPARTITION**       nodepartition       /**< pointer to node partition data structure */
   )
{
   assert(nodepartition != NULL);
   assert(*nodepartition != NULL);

   SCIPfreeMemoryArray(scip, &(*nodepartition)->representatives);
   SCIPfreeMemoryArray(scip, &(*nodepartition)->nodeclusters);
   SCIPfreeMemoryArray(scip, &(*nodepartition)->clusternodes);
   SCIPfreeMemoryArray(scip, &(*nodepartition)->clusterbegin);
   SCIPfreeMemory(scip, nodepartition);
}

#ifdef SCIP_DEBUG
static
void nodepartitionPrint(
   NODEPARTITION*        nodepartition       /**< node partition data structure */
   )
{
   int c;

   for( c = 0; c < nodepartition->nclusters; c++ )
   {
      int i;

      printf("cluster %d:", c);
      for( i = nodepartition->clusterbegin[c]; i < nodepartition->clusterbegin[c+1]; i++ )
         printf(" %d", nodepartition->clusternodes[i]);
      printf("\n");
   }
}
#endif

/** returns whether given node v is in a cluster that belongs to the partition S */
static
SCIP_Bool nodeInPartition(
   NODEPARTITION*        nodepartition,      /**< node partition data structure */
   unsigned int          partition,          /**< bit field to mark clusters in partition S */
   int                   v                   /**< node to check */
   )
{
   int cluster;
   unsigned int clusterbit;

   assert(nodepartition != NULL);

   /* if the node does not exist, it is not in the partition */
   if( v < 0 )
      return FALSE;

   cluster = nodepartition->nodeclusters[v];
   assert(0 <= cluster && cluster < nodepartition->nclusters);
   clusterbit = (1 << cluster);

   return ((partition & clusterbit) != 0);
}

/** adds given cut to LP if violated */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< coefficients of active variables in cut */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Bool             cutislocal,         /**< is the cut only locally valid? */
   int*                  ncuts               /**< pointer to count the number of added cuts */
   )
{
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_VAR** vars;
   int nvars;
   int v;
   
   assert(scip != NULL);
   assert(sepadata != NULL);      
   assert(cutcoefs != NULL);
   assert(ncuts != NULL);

   /* get active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars == 0 || vars != NULL);

   /* create the cut */
   sprintf(cutname, "mcf%d_%d", SCIPgetNLPs(scip), *ncuts);
   SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), cutrhs, 
                                 cutislocal, FALSE, sepadata->dynamiccuts) );

   /* add coefficients */
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
   for( v = 0; v < nvars; v++ )
   {
      if( SCIPisZero(scip, cutcoefs[v]) )
         continue;

      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[v], cutcoefs[v]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

   /* check efficacy */
   if( SCIPisCutEfficacious(scip, sol, cut) )
   {
      SCIPdebugMessage(" -> found MCF cut <%s>: rhs=%f, act=%f eff=%f\n",
                       cutname, cutrhs, SCIPgetRowSolActivity(scip, cut, sol), SCIPgetCutEfficacy(scip, sol, cut));
      SCIPdebug(SCIPprintRow(scip, cut, NULL));
      SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE) );
      if( !cutislocal )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      }
      (*ncuts)++;
   }

   /* release the row */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;   
}

/** enumerates cuts between subsets of the clusters */
static
SCIP_RETCODE generateClusterCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_MCFNETWORK*      mcfnetwork,         /**< MCF network structure */
   NODEPARTITION*        nodepartition,      /**< node partition data structure */
   int*                  ncuts               /**< pointer to count the number of added cuts */
   )
{
   SCIP_ROW***       nodeflowrows      = mcfnetwork->nodeflowrows;
   SCIP_Real**       nodeflowscales    = mcfnetwork->nodeflowscales;
   SCIP_ROW**        arccapacityrows   = mcfnetwork->arccapacityrows;
   SCIP_Real*        arccapacityscales = mcfnetwork->arccapacityscales;
   int*              colcommodity      = mcfnetwork->colcommodity;
   int*              arcsources        = mcfnetwork->arcsources;
   int*              arctargets        = mcfnetwork->arctargets;
   int               nnodes            = mcfnetwork->nnodes;
   int               narcs             = mcfnetwork->narcs;
   int               ncommodities      = mcfnetwork->ncommodities;
   SCIP_MCFMODELTYPE modeltype         = mcfnetwork->modeltype;

   int nclusters = nodepartition->nclusters;

   int maxtestdelta = sepadata->maxtestdelta;

   int nrows;
   int nvars;

   SCIP_Real* cutcoefs;
   SCIP_Real* deltas;
   int deltassize;
   int ndeltas;
   SCIP_Real* cmirweights;
   SCIP_Real* comcutdemands;
   SCIP_Real* comdemands;
   unsigned int partition;
   unsigned int allpartitions;

   assert((unsigned int)nclusters <= 8*sizeof(unsigned int));

   if( maxtestdelta == -1 )
      maxtestdelta = INT_MAX;

   nrows = SCIPgetNLPRows(scip);
   nvars = SCIPgetNVars(scip);

   /* Our system has the following form:
    *  (1)  \sum_{a \in \delta^+(v)} f_a^k  - \sum_{a \in \delta^-(v)} f_a^k  <=  -d_v^k   for all v \in V and k \in K
    *  (2)  \sum_{k \in K} f_a^k - c_a x_a                                    <=  0        for all a \in A
    *
    * The partitioning yields two clusters:
    *
    *                A^+
    *   cluster S  ------>  cluster T
    *              <------
    *                A^-
    *
    * Now we look at all commodities in which we have to route flow from T to S:
    *   K^+ = {k : d^k_S := sum_{v \in S} d_v^k > 0}
    *
    * Then, the base constraint of the c-MIR cut is the sum of those flow conservation constraints and the
    * capacity constraints for arcs A^-:
    *
    *   sum_{k \in K^+} sum_{v \in S} (sum_{a \in \delta^+(v)} f_a^k  - sum_{a \in \delta^-(v)} f_a^k)           <=  sum_{k \in K^+} sum_{v \in S} -d_v^k
    *                                                           + sum_{a \in A^-} sum_{k \in K} f_a^k - c_a x_a  <=  0
    */

   deltassize = 16;
   SCIP_CALL( SCIPallocMemoryArray(scip, &deltas, deltassize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cmirweights, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comcutdemands, ncommodities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comdemands, ncommodities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );

   /* loop over all possible partitions of the clusters;
    * cluster i is in S iff bit i of 'partition' is 1
    */
   if( mcfnetwork->modeltype == SCIP_MCFMODELTYPE_DIRECTED )
   {
      /* in the directed case, we try all (non-trivial) partitions */
      allpartitions = (1 << nclusters) - 1;
   }
   else
   {
      /* In the undirected case, we only need to try half of the partitions, since they are pairwise
       * equivalent (S-T is equivalent to T-S). Thus, we fix the last cluster to belong to partition T.
       */
      assert(mcfnetwork->modeltype == SCIP_MCFMODELTYPE_UNDIRECTED);
      allpartitions = (1 << (nclusters-1));
   }
   for( partition = 1; partition <= allpartitions-1; partition++ )
   {
      int v;
      int a;
      int d;

      SCIPdebugMessage("generating cuts for partition 0x%x\n", partition);

      /* clear memory */
      BMSclearMemoryArray(cmirweights, nrows);
      BMSclearMemoryArray(comcutdemands, ncommodities);
      BMSclearMemoryArray(comdemands, ncommodities);
      ndeltas = 0;

      /* Identify commodities with positive T -> S demand */
      for( v = 0; v < nnodes; v++ )
      {
         int k;

         /* check if node belongs to S */
         if( !nodeInPartition(nodepartition, partition, v) )
            continue;

         /* update commodity demand */
         for( k = 0; k < ncommodities; k++ )
         {
            SCIP_Real rhs;

            if( nodeflowrows[v][k] == NULL )
               continue;

            if( nodeflowscales[v][k] > 0.0 )
               rhs = SCIProwGetRhs(nodeflowrows[v][k]);
            else
               rhs = SCIProwGetLhs(nodeflowrows[v][k]);

            comcutdemands[k] += rhs * nodeflowscales[v][k];
         }
      }

#ifdef SCIP_DEBUG
      for( v = 0; v < ncommodities; v++ )
      {
         SCIPdebugMessage(" -> commodity %d: cutdemand=%g\n", v, comcutdemands[v]);
      }
#endif

      /* set weights of capacity rows that go from T to S, i.e., a \in A^- */
      for( a = 0; a < narcs; a++ )
      {
         SCIP_COL** rowcols;
         SCIP_Real* rowvals;
         int rowlen;
         int r;
         int j;

         assert(arccapacityrows[a] != NULL);
         assert(arcsources[a] < nnodes);
         assert(arctargets[a] < nnodes);

         /* check if this is an arc of our cut */
         if( modeltype == SCIP_MCFMODELTYPE_DIRECTED )
         {
            /* in the directed case, check if arc goes from T to S */
            if( nodeInPartition(nodepartition, partition, arcsources[a]) || !nodeInPartition(nodepartition, partition, arctargets[a]) )
               continue;
         }
         else
         {
            /* in the undirected case, check if the arc has endpoints in S and T */
            if( nodeInPartition(nodepartition, partition, arcsources[a]) && nodeInPartition(nodepartition, partition, arctargets[a]) )
            {
               /* both endpoints are in S */
               continue;
            }
            if( !nodeInPartition(nodepartition, partition, arcsources[a]) && !nodeInPartition(nodepartition, partition, arctargets[a]) )
            {
               /* both endpoints are in T */
               continue;
            }
         }

         /* use capacity row in c-MIR cut */
         r = SCIProwGetLPPos(arccapacityrows[a]);
         assert(r < nrows);
         if( r == -1 ) /* row might have been removed from LP in the meantime */
            continue;
         assert(cmirweights[r] == 0.0);

         /* if one of the arc nodes is unknown, we only use the capacity row if it does not have slack */
         if( arcsources[a] == -1 || arctargets[a] == -1 )
         {
            SCIP_Real feasibility;

            feasibility = SCIPgetRowSolFeasibility(scip, arccapacityrows[a], sol);
            if( SCIPisFeasPositive(scip, feasibility) )
               continue;
         }
         
         cmirweights[r] = arccapacityscales[a];
         SCIPdebugMessage(" -> arc %d, capacity row <%s>: weight=%g slack=%g dual=%g\n", a, SCIProwGetName(arccapacityrows[a]), cmirweights[r],
                          SCIPgetRowFeasibility(scip, arccapacityrows[a]), SCIProwGetDualsol(arccapacityrows[a]));
         SCIPdebug(SCIPprintRow(scip, arccapacityrows[a], NULL));

         /* extract useful deltas for c-MIR scaling and update the demand value for commodities (in binary flow model) */
         rowcols = SCIProwGetCols(arccapacityrows[a]);
         rowvals = SCIProwGetVals(arccapacityrows[a]);
         rowlen = SCIProwGetNLPNonz(arccapacityrows[a]);
         for( j = 0; j < rowlen; j++ )
         {
            SCIP_Real coef;
            SCIP_Bool exists;
            int left;
            int right;
            int k;
            int c;

            coef = rowvals[j] * arccapacityscales[a];
            coef = ABS(coef);
            
            /* binary search if we already know this coefficient */
            exists = FALSE;
            left = 0;
            right = ndeltas-1;
            while( left <= right )
            {
               int mid = (left+right)/2;
               if( SCIPisEQ(scip, deltas[mid], coef) )
               {
                  exists = TRUE;
                  break;
               }
               else if( coef < deltas[mid] )
                  right = mid-1;
               else
                  left = mid+1;
            }

            /* insert new candidate value */
            if( !exists )
            {
               assert(right == left-1);
               assert(ndeltas <= deltassize);
               if( ndeltas == deltassize )
               {
                  deltassize *= 2;
                  SCIP_CALL( SCIPreallocMemoryArray(scip, &deltas, deltassize) );
               }
               if( left < ndeltas )
                  BMScopyMemoryArray(&deltas[left+1], &deltas[left], ndeltas - left);
               deltas[left] = coef;
               ndeltas++;
            }

            /* update commodity demands */
            c = SCIPcolGetLPPos(rowcols[j]);
            assert(0 <= c && c < SCIPgetNLPCols(scip));
            k = colcommodity[c];
            if( k >= 0 )
               comdemands[k] = coef;
         }
      }

      /* set weights of node flow conservation constraints in c-MIR aggregation */
      for( v = 0; v < nnodes; v++ )
      {
         int k;

         /* check if node belongs to S */
         if( !nodeInPartition(nodepartition, partition, v) )
            continue;

         /* aggregate flow conservation constraints of the 'active' commodities */
         for( k = 0; k < ncommodities; k++ )
         {
            SCIP_Real scale;
            int r;

            /* if commodity was not hit by the capacity constraints of the cut in the graph, ignore the commodity */
            if ( comdemands[k] == 0.0 )
               continue;

            scale = comdemands[k];
            if( modeltype == SCIP_MCFMODELTYPE_DIRECTED )
            {
               /* in the directed case, use rows of commodities with positive demand (negative -d_k) */
               if( !SCIPisNegative(scip, comcutdemands[k]) )
                  continue;
            }
            else
            {
               /* in the undirected case, use rows of commodities with non-zero demand */
               if( SCIPisZero(scip, comcutdemands[k]) )
                  continue;

               /* if the demand (-d_k) is negative (i.e., points into the wrong direction), we use the flow constraint with opposite sign */
               if( comcutdemands[k] > 0.0 )
                  scale *= -1.0;
            }
            if( nodeflowrows[v][k] == NULL )
               continue;
            
            r = SCIProwGetLPPos(nodeflowrows[v][k]);
            assert(r < nrows);
            if( r >= 0 ) /* row might have been removed from LP in the meantime */
            {
               SCIP_Real feasibility;

               assert(cmirweights[r] == 0.0);

               /* ignore rows with slack */
               feasibility = SCIPgetRowSolFeasibility(scip, nodeflowrows[v][k], sol);
               if( !SCIPisFeasPositive(scip, feasibility) )
               {
                  cmirweights[r] = scale * nodeflowscales[v][k];
                  SCIPdebugMessage(" -> node %d, commodity %d, flow row <%s>: scale=%g weight=%g slack=%g dual=%g\n",
                                   v, k, SCIProwGetName(nodeflowrows[v][k]), scale, cmirweights[r], 
                                   SCIPgetRowFeasibility(scip, nodeflowrows[v][k]), SCIProwGetDualsol(nodeflowrows[v][k]));
                  SCIPdebug(SCIPprintRow(scip, nodeflowrows[v][k], NULL));
               }
            }
         }
      }

      /* try out deltas to generate c-MIR cuts */
      /**@todo use deltas as in the c-MIR separator from the mksetcoefs returned by SCIPcalcMIR() */
      for( d = 0; d < maxtestdelta && d < ndeltas; d++ )
      {
         SCIP_Real cutrhs;
         SCIP_Real cutact;
         SCIP_Bool success;
         SCIP_Bool cutislocal;

         /* do not use too small deltas */
         if( SCIPisFeasZero(scip, deltas[d]) )
            continue;

         SCIP_CALL( SCIPcalcMIR(scip, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, sepadata->fixintegralrhs, NULL, NULL, 
                                sepadata->maxrowfac, MINFRAC, MAXFRAC, cmirweights, 1.0/deltas[d], NULL, cutcoefs, &cutrhs, &cutact, 
                                &success, &cutislocal) );
         assert(ALLOWLOCAL || !cutislocal);

#if 0
#ifdef SCIP_DEBUG
         SCIPdebugMessage("    =>");
         for( a = 0; a < nvars; a++ )
         {
            if( !SCIPisZero(scip, cutcoefs[a]) )
               printf(" %+g<%s>", cutcoefs[a], SCIPvarGetName(SCIPgetVars(scip)[a]));
         }
         printf(" <= %g  (act: %g)\n", cutrhs, cutact);
#endif
#endif

         if( success && SCIPisFeasGT(scip, cutact, cutrhs) )
         {
            SCIPdebugMessage(" -> delta = %g  -> rhs: %g, act: %g\n", deltas[d], cutrhs, cutact);
            SCIP_CALL( addCut(scip, sepadata, sol, cutcoefs, cutrhs, cutislocal, ncuts) );
         }
      }
   }

   /* free local memory */
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &comdemands);
   SCIPfreeBufferArray(scip, &comcutdemands);
   SCIPfreeBufferArray(scip, &cmirweights);
   SCIPfreeMemoryArray(scip, &deltas);

   return SCIP_OKAY;
}

/** searches and adds MCF network cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_MCFNETWORK** mcfnetworks;
   int nmcfnetworks;
   int ncuts;
   int i;

   assert(result != NULL);

   printf("************** [%6.2f] MCF start\n", SCIPgetSolvingTime(scip)); /*????????????????????*/
   *result = SCIP_DIDNOTRUN;
   ncuts = 0;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* get or extract network flow structure */
   if( sepadata->nmcfnetworks == -1 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( mcfnetworkExtract(scip, sepadata, &sepadata->mcfnetworks, &sepadata->nmcfnetworks) );
      /*?????????????#ifdef SCIP_DEBUG*/
      printf/*??????????????SCIPdebugMessage*/("extracted %d networks\n", sepadata->nmcfnetworks);
      for( i = 0; i < sepadata->nmcfnetworks; i++ )
      {
         printf/*??????????????SCIPdebugMessage*/(" -> extracted network %d has %d nodes, %d arcs, and %d commodities (modeltype: %s, flowtype: %s)\n",
                                                  i, sepadata->mcfnetworks[i]->nnodes, sepadata->mcfnetworks[i]->narcs,
                                                  sepadata->mcfnetworks[i]->ncommodities,
                                                  sepadata->mcfnetworks[i]->modeltype == SCIP_MCFMODELTYPE_DIRECTED ? "directed" : "undirected",
                                                  sepadata->mcfnetworks[i]->flowtype == SCIP_MCFFLOWTYPE_CONTINUOUS ? "continuous" : "binary");
         SCIPdebug( mcfnetworkPrint(sepadata->mcfnetworks[i]) );
      }
      /*?????????????#endif*/
   }
   assert(sepadata->nmcfnetworks >= 0);
   assert(sepadata->mcfnetworks != NULL || sepadata->nmcfnetworks == 0);
   mcfnetworks = sepadata->mcfnetworks;
   nmcfnetworks = sepadata->nmcfnetworks;

   if( nmcfnetworks > 0 )
   {
      /* separate cuts */
      *result = SCIP_DIDNOTFIND;

      for( i = 0; i < nmcfnetworks; i++ )
      {
         SCIP_MCFNETWORK* mcfnetwork;
         NODEPARTITION* nodepartition;

         mcfnetwork = mcfnetworks[i];

         /* if the network does not have at least 2 nodes and 1 arc, we did not create it */
         assert(mcfnetwork->nnodes >= 2);
         assert(mcfnetwork->narcs >= 1);

         printf("************** [%6.2f] MCF create node partition\n", SCIPgetSolvingTime(scip)); /*????????????????????*/
         /* partition nodes into a small number of clusters */
         SCIP_CALL( nodepartitionCreate(scip, mcfnetwork, &nodepartition, sepadata->nclusters) );
         SCIPdebug( nodepartitionPrint(nodepartition) );

         printf("************** [%6.2f] MCF generate cuts\n", SCIPgetSolvingTime(scip)); /*????????????????????*/
         /* enumerate cuts between subsets of the clusters */
         SCIP_CALL( generateClusterCuts(scip, sepadata, sol, mcfnetwork, nodepartition, &ncuts) );

         /* free node partition */
         nodepartitionFree(scip, &nodepartition);

         printf/*?????????????SCIPdebugMessage*/("MCF network has %d nodes, %d arcs, %d commodities. Found %d MCF network cuts.\n",
                                                 mcfnetwork->nnodes, mcfnetwork->narcs, mcfnetwork->ncommodities, ncuts);

         /* adjust result code */
         if( ncuts > 0 )
            *result = SCIP_SEPARATED;
      }
   }
   printf("************** [%6.2f] MCF end\n", SCIPgetSolvingTime(scip)); /*????????????????????*/

   return SCIP_OKAY;
}




/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMcf)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->mcfnetworks == NULL);
   assert(sepadata->nmcfnetworks == -1);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitMcf NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitMcf NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolMcf NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMcf)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   int i;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* free MCF networks */
   for( i = 0; i < sepadata->nmcfnetworks; i++ )
   {
      SCIP_CALL( mcfnetworkFree(scip, &sepadata->mcfnetworks[i]) );
   }
   SCIPfreeMemoryArrayNull(scip, &sepadata->mcfnetworks);
   sepadata->nmcfnetworks = -1;

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMcf)
{  /*lint --e{715}*/
   /* separate cuts on the LP solution */
   SCIP_CALL( separateCuts(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMcf)
{  /*lint --e{715}*/
   /* separate cuts on the given primal solution */
   SCIP_CALL( separateCuts(scip, sepa, sol, result) );

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the mcf separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMcf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

#if 1
   /* disabled, because separator is not yet finished */
   return SCIP_OKAY;
#endif

   /* create cmir separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
   sepadata->mcfnetworks = NULL;
   sepadata->nmcfnetworks = -1;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeMcf, sepaInitMcf, sepaExitMcf, 
         sepaInitsolMcf, sepaExitsolMcf,
         sepaExeclpMcf, sepaExecsolMcf,
         sepadata) );

   /* add mcf separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mcf/nclusters",
         "number of clusters to generate in the shrunken network",
         &sepadata->nclusters, TRUE, DEFAULT_NCLUSTERS, 2, 8*sizeof(unsigned int), NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/mcf/maxrowfac",
         "maximal row aggregation factor",
         &sepadata->maxrowfac, TRUE, DEFAULT_MAXROWFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mcf/maxtestdelta",
         "maximal number of different deltas to try (-1: unlimited)",
         &sepadata->maxtestdelta, TRUE, DEFAULT_MAXTESTDELTA, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mcf/modeltype",
         "model type of network (0: auto, 1:directed, 2:undirected)",
         &sepadata->modeltype, TRUE, DEFAULT_MODELTYPE, 0, 2, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mcf/flowtype",
         "type of flow variables in the network (0: auto, 1:continuous/integer, 2:binary)",
         &sepadata->flowtype, TRUE, DEFAULT_FLOWTYPE, 0, 2, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/mcf/fixintegralrhs",
         "should an additional variable be complemented if f0 = 0?",
         &sepadata->fixintegralrhs, TRUE, DEFAULT_FIXINTEGRALRHS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/mcf/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
