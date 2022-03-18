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

/**@file   extreducedefs.h
 * @brief  includes extended reductions definitions and inline methods used for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_EXTREDUCEDEFS_H_
#define APPLICATIONS_STP_SRC_EXTREDUCEDEFS_H_


#ifdef __cplusplus
extern "C" {
#endif


#define STP_EXT_CLOSENODES_MAXN 50
#define STP_EXT_MAXSTACKSIZE  5000
#define STP_EXT_MAXNCOMPS     1000
#define STP_EXT_MAXDFSDEPTH 6
#define STP_EXT_MAXDFSDEPTH_GUARD (STP_EXT_MAXDFSDEPTH + 2)
#define STP_EXT_MIDDFSDEPTH 5
#define STP_EXT_MINDFSDEPTH 4
#define STP_EXT_MAXGRAD 9
#define STP_EXT_MAXEDGES 500
#define STP_EXT_DEPTH_MIDNEDGES 50000
#define STP_EXT_DEPTH_MAXNEDGES 100000
#define STP_EXTTREE_MAXNEDGES 25
#define STP_EXTTREE_MAXNLEAVES 20
#define STP_EXTTREE_MAXNLEAVES_GUARD (STP_EXTTREE_MAXNLEAVES + STP_EXT_MAXGRAD)
#define EXT_EDGE_WRAPPED -10


#define EXT_STATE_NONE     0
#define EXT_STATE_EXPANDED 1
#define EXT_STATE_MARKED   2

#define STP_MLDISTS_ID_EMPTY      -2
#define STP_MLDISTS_ID_UNSET      -1
#define STP_MLDISTS_DIST_UNSET    -1.0

#include "scip/scip.h"
#include "reduce.h"
#include "graph.h"

/** structure for storing distances in the extension tree */
typedef struct multi_level_distances_storage MLDISTS;

/** structure for storing data for algorithms based on contraction of extension tree */
typedef struct extension_tree_contraction CONTRACT;


/** path root state */
typedef struct pathroot_state
{
   int pathroot_id;
   int pathroot_nrecomps;
} PRSTATE;

/** distance data */
typedef struct distance_data
{
   //    SCIP_Bool*  pathrootsd_isdirty;
   SD* sdistdata;                    /* can be NULL */
   DHEAP* dheap;
   RANGE* closenodes_range;          /* of size nnodes */
   int* closenodes_indices;          /* of size closenodes_totalsize */
   int* closenodes_prededges;        /* of size closenodes_totalsize; NULL in opt mode */
   SCIP_Real* closenodes_distances;  /* of size closenodes_totalsize */
   int* pathroot_blocksizes;         /* of size nedges / 2 */
   int* pathroot_blocksizesmax;      /* of size nedges / 2 */
   PRSTATE** pathroot_blocks;        /* of size nedges / 2 */
   SCIP_Bool* pathroot_isdirty;      /* of size nnodes */
   int* pathroot_nrecomps;           /* of size nnodes */
   DIJK* dijkdata;
   int closenodes_totalsize;
   int closenodes_maxnumber;
   SCIP_Bool hasPathReplacement;
} DISTDATA;


/** permanent extension data */
typedef struct extension_data_permanent
{
/* non-owned data: */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator (initialized to NULL)  */
   REDCOST*              redcostdata;        /**< reduced cost data (initialized to NULL) */
   DISTDATA*             distdata_default;   /**< default distance data (initialized to NULL) */
   DISTDATA*             distdata_biased;    /**< biased distance data (initialized to NULL) */
/* owned data: */
   CONTRACT*             contration;         /**< contraction data */
   DCMST*                dcmst;              /**< dynamic MST */
   CSRDEPO*              msts_comp;          /**< storage for MSTs with extending node of the given level */
   CSRDEPO*              msts_levelbase;     /**< storage for MSTs without extending node of the level below */
   MLDISTS*              sds_horizontal;     /**< SDs from deepest leaves to remaining ones */
   MLDISTS*              sds_vertical;       /**< SDs between leaves of same depth */
   MLDISTS*              sdsbias_horizontal; /**< as above, but for biased SDs; can be NULL */
   MLDISTS*              sdsbias_vertical;   /**< as above, but for biased SDs; can be NULL */
   STP_Bool*             edgedeleted;        /**< (non-owned!) edge array to mark which directed edge can be removed */
   SCIP_Bool*            isterm;             /**< marks whether node is a terminal (or proper terminal for PC) */
   SCIP_Real*            bottleneckDistNode; /**< needs to be set to -1.0 (size nnodes) */
   SCIP_Real*            pcSdToNode;         /**< needs to be set to -1.0 for all nodes, or NULL if not Pc */
   int*                  tree_deg;           /**< -1 for forbidden nodes (e.g. PC terminals), nnodes for tail, 0 otherwise; in method: position ( > 0) for nodes in tree */
   const int*            result;             /**< solution array or NULL... NON-OWNED! */
   STP_Vectype(int)*     nodes_implications; /**< implied nodes for each node */
   int                   nnodes;             /**< number of nodes */
   int                   tree_maxnleaves;
   int                   tree_maxdepth;
   int                   tree_maxnedges;
   SCIP_Bool             redcostEqualAllow;  /**< delete also for equality of reduced costs? */
   SCIP_Bool             useSdBias;          /**< use biased bottleneck Steiner distance? (only for pseudo-elimination) */
   SCIP_Bool             solIsValid;         /**< use primal solution? */
   enum EXTRED_MODE      mode;               /**< mode */
} EXTPERMA;


/** Reduction data; just used internally.
 *  Stores information for SDs between tree vertices,
 *  reduced costs, and pseudo-ancestor conflicts. */
typedef struct reduction_data
{
   CONTRACT*             contration;         /**< contraction data */
   DCMST* const dcmst;
   CSRDEPO* const msts_comp;
   CSRDEPO* const msts_levelbase;
   MLDISTS* const sds_horizontal;
   MLDISTS* const sds_vertical;
   MLDISTS* const sdsbias_horizontal;      /**< can be NULL */
   MLDISTS* const sdsbias_vertical;        /**< can be NULL */
   const STP_Bool* const edgedeleted;
   int* const pseudoancestor_mark;
   STP_Vectype(int)* nodes_implications;
   SCIP_Real* const redcost_treecosts; /* for each level */
   SCIP_Bool* const redcost_noReversedTree; /* for each level */
   SCIP_Real* const redcost_treenodeswaps; /* for each level and each node */
   const int redcost_nlevels;
   const SCIP_Bool redcost_allowEquality;
   const SCIP_Bool sdsbias_use;
} REDDATA;


/** PC/MW data; just used internally. */
typedef struct pcmw_specific_data
{
  SCIP_Real* const pcSdToNode;                 /**< all entries need to be set to -1.0 */
  int* const pcSdCands;
  SCIP_Real tree_innerPrize;
  int nPcSdCands;
  int pcSdStart;
} PCDATA;


/** initial extension component
 *  NOTE: it is vital the the first edge of the star component comes from the root!
 *  (will thus be constantly asserted) */
typedef struct initial_extension_component
{
   int*                  compedges;          /**< edges of the component */
   int*                  extleaves;          /**< leaves to extend from */
   int                   nextleaves;         /**< number of extension nodes */
   int                   ncompedges;         /**< number of edges of the component */
   int                   comproot;           /**< component root */
   int                   genstar_centeredge; /**< center-edge or -1 */
   SCIP_Bool             allowReversion;     /**< allow change of comproot? (with extleaves = \{comproot\}) */
} EXTCOMP;


/** extension data; just used internally */
typedef struct extension_data
{
   int* const extstack_data;
   int* const extstack_start;
   int* const extstack_state;
   int* const tree_leaves;
   int* const tree_innerNodes;
   int* const tree_edges;
   int* const tree_deg;                         /**< -1 for forbidden nodes (e.g. PC terminals), nnodes for current tail,
                                                      0 otherwise; in method: position ( > 0) for nodes in tree */
   SCIP_Real* const tree_bottleneckDistNode;    /**< needs to be set to -1.0 (for all nodes) */
   int* const tree_parentNode;
   SCIP_Real* const tree_parentEdgeCost;        /**< of size nnodes */
   const SCIP_Bool* const node_isterm;          /**< marks whether node is a terminal (or proper terminal for PC) */
   REDDATA* const reddata;
   DISTDATA* const distdata;
   DISTDATA* const distdata_biased;       /**< can be NULL */
   const REDCOST* const redcostdata;
   PCDATA* const pcdata;
   SCIP_Bool* const sdeq_edgesIsForbidden;
   STP_Vectype(int) sdeq_resetStack;
   SCIP_Bool sdeq_hasForbiddenEdges;
   SCIP_Real tree_cost;
   int tree_nDelUpArcs;
   int tree_root;
   int tree_starcenter;
   int tree_nedges;
   int tree_depth;
   int tree_nleaves;
   int tree_ninnerNodes;
   int extstack_ncomponents;
   int ncostupdatestalls;           /**< cost update stalls counter */
   int genstar_centeredge;           /* center edge or -1  */
   const int extstack_maxncomponents;
   const int extstack_maxsize;
   const int tree_maxnleaves;
   const int tree_maxdepth;
   const int tree_maxnedges;
   enum EXTRED_MODE      mode;               /**< mode */
   const EXTCOMP* const extcomp;
} EXTDATA;



/* inline methods
 */


/** prize-collecting problem? */
static inline
SCIP_Bool extProbIsPc(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(graph && extdata);
   assert(graph_pc_isPc(graph) == (extdata->pcdata->pcSdToNode != NULL));

   return (extdata->pcdata->pcSdToNode != NULL);
}


/** currently at initial component? */
static inline
SCIP_Bool extIsAtInitialComp(
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(extdata);
   assert(extdata->extstack_ncomponents > 0);

   return ((extdata->extstack_ncomponents - 1) == 0);
}


/** is the initial component a single edge? */
static inline
SCIP_Bool extInitialCompIsEdge(
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(extdata);
   assert(extdata->extstack_start[1] >= 1);
   assert(extdata->extstack_start[1] == 1 || extdata->extstack_start[1] >= 3);
   assert(extdata->extstack_start[1] != 1 || extdata->genstar_centeredge == -1);

   return (extdata->extstack_start[1] == 1);
}


/** is the initial component a star? */
static inline
SCIP_Bool extInitialCompIsStar(
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(extdata);
   assert(extdata->extstack_start[1] >= 1);
   assert(extdata->extstack_start[1] == 1 || extdata->extstack_start[1] >= 3);

   return (extdata->extstack_start[1] >= 3 && extdata->genstar_centeredge == -1);
}


/** is the initial component a star? */
static inline
SCIP_Bool extInitialCompIsGenStar(
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(extdata);

   return (extdata->genstar_centeredge != -1);
}


/** currently at initial star? */
static inline
SCIP_Bool extIsAtInitialStar(
   const EXTDATA*        extdata             /**< extension data */
)
{
   return (extInitialCompIsStar(extdata) && extIsAtInitialComp(extdata));
}


/** currently at initial general star? */
static inline
SCIP_Bool extIsAtInitialGenStar(
   const EXTDATA*        extdata             /**< extension data */
)
{
   return (extInitialCompIsGenStar(extdata) && extIsAtInitialComp(extdata));
}


/** returns current position in the stack */
static inline
int extStackGetPosition(
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(extdata);
   assert(extdata->extstack_ncomponents > 0);

   return (extdata->extstack_ncomponents - 1);
}


/** returns root of top component on the stack */
static inline
int extStackGetTopRoot(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int stackpos = extStackGetPosition(extdata);
   const int comproot = graph->tail[extdata->extstack_data[extdata->extstack_start[stackpos]]];

   assert(comproot >= 0);
   assert(extdata->tree_deg[comproot] >= 1 || comproot == extdata->tree_root);

   return comproot;
}


/** returns start of outgoing edges */
static inline
int extStackGetOutEdgesStart(
   const EXTDATA*        extdata,            /**< extension data */
   int                   stackpos            /**< current position */
)
{
   int start;

   assert(extdata);
   assert(stackpos >= 0 && stackpos <= extStackGetPosition(extdata));

   if( stackpos == 0 && extInitialCompIsStar(extdata) )
   {
      start = extdata->extstack_start[stackpos] + 1;
   }
   else if( stackpos == 0 && extInitialCompIsGenStar(extdata) )
   {
      start = extdata->extstack_start[stackpos] + 2;
   }
   else
   {
      start = extdata->extstack_start[stackpos];
   }

   return start;
}


/** returns end of outgoing edges */
static inline
int extStackGetOutEdgesEnd(
   const EXTDATA*        extdata,            /**< extension data */
   int                   stackpos            /**< current position */
)
{
   assert(extdata);
   assert(stackpos >= 0 && stackpos <= extStackGetPosition(extdata));

   return (extdata->extstack_start[stackpos + 1]);
}


/** returns start of outgoing edges of top */
static inline
int extStackGetTopOutEdgesStart(
   const EXTDATA*        extdata,            /**< extension data */
   int                   stackpos            /**< current position */
)
{
   assert(extdata);
   assert(stackpos == extStackGetPosition(extdata));

   return extStackGetOutEdgesStart(extdata, stackpos);
}


/** returns end of outgoing edges of top */
static inline
int extStackGetTopOutEdgesEnd(
   const EXTDATA*        extdata,            /**< extension data */
   int                   stackpos            /**< current position */
)
{
   assert(extdata);
   assert(stackpos == extStackGetPosition(extdata));

   return extStackGetOutEdgesEnd(extdata, stackpos);
}


/** is component at top position wrapped? */
static inline
SCIP_Bool extStackTopIsWrapped(
   const EXTDATA*        extdata             /**< extension data */
   )
{
   const int pos = extStackGetPosition(extdata);

   assert(extdata->extstack_start[pos + 1] - extdata->extstack_start[pos] >= 1);

   if( extdata->extstack_start[pos + 1] - extdata->extstack_start[pos] == 1
    && extdata->extstack_data[extdata->extstack_start[pos]] == EXT_EDGE_WRAPPED  )
   {
      return TRUE;
   }

   assert(extdata->extstack_data[extdata->extstack_start[pos]] != EXT_EDGE_WRAPPED);

   return FALSE;
}


/** is component at top position a singleton edge? */
static inline
SCIP_Bool extStackTopIsSingleton(
   const EXTDATA*        extdata             /**< extension data */
   )
{
   const int pos = extStackGetPosition(extdata);
   assert(extdata->extstack_start[pos + 1] - extdata->extstack_start[pos] >= 1);

   return (extdata->extstack_start[pos + 1] - extdata->extstack_start[pos] == 1);
}


/** Finds position of given leaf in leaves data.
 *  Returns -1 if leaf could not be found. */
static inline
int extLeafFindPos(
   const EXTDATA*        extdata,            /**< extension data */
   int                   leaf                /**< leaf to find */
)
{
   int i;
   const int* const tree_leaves = extdata->tree_leaves;

   assert(tree_leaves);
   assert(extdata->tree_nleaves > 1);
   assert(leaf >= 0 && extdata->tree_deg[leaf] >= 1);

   for( i = extdata->tree_nleaves - 1; i >= 0; i-- )
   {
      const int currleaf = tree_leaves[i];

      if( currleaf == leaf )
         break;
   }

   return i;
}


/** can we extend the tree from given leaf? */
inline static
SCIP_Bool extLeafIsExtendable(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   int                   leaf                /**< the leaf */
)
{
   assert(graph && isterm);
   assert(leaf >= 0 && leaf < graph->knots);

   return (!isterm[leaf] && graph->grad[leaf] <= STP_EXT_MAXGRAD);
}


/** Gets proper SD value. I.e., non-negative, but possible FARAWAY */
static inline
SCIP_Real extSdGetProper(
   SCIP_Real             sd_in               /**< the special distance */
)
{
   assert(LE(sd_in, FARAWAY));
   assert(EQ(sd_in, -1.0) || GE(sd_in, 0.0));

   return (sd_in >= -0.5) ? sd_in : FARAWAY;
}


/** is given SD non-trivial? */
static inline
SCIP_Real extSdIsNonTrivial(
   SCIP_Real             specialDist         /**< SD */
  )
{
   assert(specialDist >= 0 || EQ(specialDist, -1.0));
   assert(LT(specialDist, FARAWAY) && "todo: might fail for edge replacement");

   return (specialDist >= -0.5);
}


/** does reduction data have biased SD information? */
static inline
SCIP_Bool extReddataHasBiasedSds(
   const REDDATA*        reddata             /**< reduction data */
)
{
   assert(reddata);
   assert((reddata->sdsbias_vertical != NULL) == (reddata->sdsbias_horizontal != NULL));
   assert(!reddata->sdsbias_use || (reddata->sdsbias_vertical != NULL));

   return (reddata->sdsbias_use);
}



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_EXTREDUCEDEFS_H_ */
