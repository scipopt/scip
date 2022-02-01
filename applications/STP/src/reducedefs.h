/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reducedefs.h
 * @brief  includes reductions definitions and inline methods used for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef APPLICATIONS_STP_SRC_REDUCEDEFS_H_
#define APPLICATIONS_STP_SRC_REDUCEDEFS_H_

#include "scip/scip.h"
#include "portab.h"


#ifdef __cplusplus
extern "C" {
#endif


#define STP_REDUCTION_NONE      0
#define STP_REDUCTION_BASIC     1
#define STP_REDUCTION_ADVANCED  2

#define STP_DAMODE_HOPS         -9991
#define STP_DAMODE_FAST          0
#define STP_DAMODE_MEDIUM        1
#define STP_DAMODE_EXTENSIVE     2


/** lightweight minimum spanning tree structure that allows to add vertices to given MST on complete graph (in CSR format) */
typedef struct dynamic_complete_minimum_spanning_tree DCMST;

/** auxiliary data structure for ruling out all 1-hop stars of a given node */
typedef struct node_one_hop_star STAR;

/** SD distance graph data */
typedef struct special_distance_graph SDGRAPH;

/** SD neighbors */
typedef struct special_distance_neighbors SDN;

/** link-cut tree for bottleneck operations */
typedef struct bottleneck_link_cut_tree BLCTREE;

/** primal solution data retained during reduction process */
typedef struct reduction_solution_storage REDSOL;

/** INTERNAL primal solution data retained during reduction loop */
typedef struct reduction_local_solution_storage REDSOLLOCAL;


enum EXTRED_MODE { extred_none = 0, extred_fast = 1, extred_full = 2 };


/** reduction parameters */
typedef struct reduction_parameters
{
   SCIP_Bool             dualascent;         /**< do dual-ascent reduction? */
   SCIP_Bool             boundreduce;        /**< do bound-based reduction? */
   SCIP_Bool             nodereplacing;      /**< should node replacement (by edges) be performed? */
   int                   reductbound;        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   int                   reductbound_min;    /**< absolute minimum */
   SCIP_Bool             userec;             /**< use recombination heuristic? */
   SCIP_Bool             fullreduce;         /**< use full reductions? (including extended techniques) */
   SCIP_Bool             usestrongreds;      /**< allow strong reductions? */
} RPARAMS;


// todo: parameter so that after restart only standard SD is performed...
// todo: need to adapt reductbound?? make smaller... */
/** bi-decomposition reduction parameters */
typedef struct bidecomposition_reduction_parameters
{
   int                   depth;              /**< current depth */
   int                   maxdepth;           /**< maximum recursive depth of decomposition */
   SCIP_Bool             newLevelStarted;    /**< no level? */
} BIDECPARAMS;


/** reduction information and some buffers */
typedef struct reduction_base
{
   RPARAMS*              redparameters;      /**< parameters */
   BIDECPARAMS*          bidecompparams;     /**< bidecomposition parameters or NULL */
   int*                  solnode;            /**< solution nodes array (or NULL) */
   REDSOL*               redsol;             /**< primal solution container */
 /* buffer: */
   PATH* vnoi;
   PATH* path;
   SCIP_Real* nodearrreal;
   int* heap;
   int* state;
   int* vbase;
   int* nodearrint;
   int* edgearrint;
   int* nodearrint2;
   STP_Bool* nodearrchar;
} REDBASE;



/** Stores data for computation of special distance/bottleneck distance computations  */
typedef struct special_distance_storage
{
   SDPROFIT*             sdprofit;           /**< SD bias for nodes  (or NULL) */
   SDGRAPH*              sdgraph;            /**< special distance graph on terminals      */
   TPATHS*               terminalpaths;      /**< terminal paths                 */
   SDN*                  sdneighbors;        /**< neighbors */
   BLCTREE*              blctree;            /**< bottleneck tree (or NULL) */
   SCIP_Bool             isBiased;           /**< are the SDs biased? */
   SCIP_Bool             hasNeigborUpdate;   /**< with neighbor update? NOTE: does not allow certain methods */
} SD;


/** lightweight store for implied profit */
struct special_distance_implied_profit
{
   SCIP_Real* RESTRICT   nodes_bias;         /**< bias per node */
   SCIP_Real* RESTRICT   nodes_bias2;        /**< second best bias per node */
   int* RESTRICT         nodes_biassource;   /**< source terminal per node */
   int* RESTRICT         nodes_biassource2;  /**< second source terminal per node */
};


/** reduced cost reduction parameters */
typedef struct reduce_costs_reduction_parameters
{
   int                   damode;             /**< mode */
   enum EXTRED_MODE      extredMode;         /**< mode of extended reductions */
   SCIP_Bool             useRec;             /**< use recombination heuristic? */
   SCIP_Bool             useSlackPrune;      /**< use slack-prune heuristic? */
   SCIP_Bool             nodereplacing;      /**< should node replacement (by edges) be performed? */
   /* PC/MW only values: */
   SCIP_Bool             pcmw_solbasedda;    /**< rerun Da based on best primal solution */
   SCIP_Bool             pcmw_useMultRoots;  /**< vary root for DA? (if possible) */
   SCIP_Bool             pcmw_markroots;     /**< should terminals proven to be part of an opt. sol. be marked as such? */
   SCIP_Bool             pcmw_fastDa;        /**< run dual ascent heuristic in fast mode? */
} RPDA;


/** single special distance for PC */
typedef struct single_special_distance_pc
{
   PATH*  pathtail;
   PATH*  pathhead;
   int*    heap;
   int*    statetail;
   int*    statehead;
   int*    memlbltail;
   int*    memlblhead;
   int*    pathmaxnodetail;
   int*    pathmaxnodehead;
} SD1PC;



/** gets profit for given node */
inline static
SCIP_Real reduce_sdprofitGetProfit(
   const SDPROFIT*      sdprofit,           /**< the SD profit */
   int                  node,               /**< node to get profit for */
   int                  nonsource1,         /**< node that should not be a source */
   int                  nonsource2          /**< node that should not be a source */
)
{
   const int source1 = sdprofit->nodes_biassource[node];

   assert(nonsource1 != nonsource2 || nonsource1 == -1);
   assert(GE(sdprofit->nodes_bias[node], 0.0));
   assert(LE(sdprofit->nodes_bias[node], FARAWAY));
   assert(GE(sdprofit->nodes_bias2[node], 0.0));
   assert(LE(sdprofit->nodes_bias2[node], FARAWAY));
   assert(GE(sdprofit->nodes_bias[node], sdprofit->nodes_bias2[node]));

   if( source1 != nonsource1 && source1 != nonsource2 )
   {
      return sdprofit->nodes_bias[node];
   }
   else
   {
      const int source2 = sdprofit->nodes_biassource2[node];

      if( source2 != nonsource1 && source2 != nonsource2 )
      {
         return sdprofit->nodes_bias2[node];
      }
   }

   return 0.0;
}


/** gets biased distance */
inline static
SCIP_Real reduce_sdprofitGetBiasedDist(
   const SDPROFIT*      sdprofit,           /**< the SD profit */
   int                  node,               /**< node along which to get biased distance */
   SCIP_Real            edgecost,           /**< edge cost */
   SCIP_Real            nodedist,           /**< node distance */
   int                  nonsource1,         /**< node that should not be a source */
   int                  nonsource2          /**< node that should not be a source */
)
{
   SCIP_Real distnew = nodedist + edgecost;
   const SCIP_Real profit = reduce_sdprofitGetProfit(sdprofit, node, nonsource1, nonsource2);
   SCIP_Real bias = MIN(edgecost, profit);
   if( nodedist < bias )
      bias = nodedist;

   distnew -= bias;

   assert(GE(profit, 0.0));
   assert(GE(bias, 0.0));
   assert(GE(edgecost, 0.0));
   assert(GE(nodedist, 0.0));
   assert(GE(distnew, 0.0));

   return distnew;
}



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_REDUCEDEFS_H_ */
