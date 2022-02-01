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

/**@file   graphdefs.h
 * @brief  includes graph definitions used for Steiner tree problems
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef APPLICATIONS_STP_SRC_GRAPHDEFS_H_
#define APPLICATIONS_STP_SRC_GRAPHDEFS_H_

#include "scip/scip.h"
#include "misc_stp.h"
#include "portab.h"


#ifdef __cplusplus
extern "C" {
#endif



#define STP_SPG                      0
#define STP_SAP                      1
#define STP_PCSPG                    2
#define STP_RPCSPG                   3
#define STP_NWSPG                    4
#define STP_DCSTP                    5
#define STP_NWPTSPG                  6
#define STP_RSMT                     7
#define STP_OARSMT                   8
#define STP_MWCSP                    9
#define STP_DHCSTP                   10
#define STP_GSTP                     11
#define STP_RMWCSP                   12
#define STP_BRMWCSP                  13

#define EAT_FREE     -1
#define EAT_LAST     -2
#define EAT_HIDE     -3

#define STP_TERM           0        /**< terminal */
#define STP_TERM_NONE     -1        /**< non-terminal */
#define STP_TERM_PSEUDO   -2        /**< pseudo-terminal (for PC/MW variants) */
#define STP_TERM_NONLEAF  -3        /**< non-leaf (pseudo-) terminal (for PC/MW variants) */

#define STP_CENTER_OK    0           /**< do nothing */
#define STP_CENTER_DEG   1           /**< find maximum degree */
#define STP_CENTER_SUM   2           /**< find the minimum distance sum */
#define STP_CENTER_MIN   3           /**< find the minimum largest distance */
#define STP_CENTER_ALL   4           /**< find the minimum distance sum to all knots */

#define TERM2EDGE_NOTERM      -1    /**< for PC/MW: vertex is no terminal */
#define TERM2EDGE_FIXEDTERM   -2    /**< for PC/MW: vertex is fixed terminal; artificial root is also considered a fixed terminal */
#define TERM2EDGE_NONLEAFTERM -3    /**< for PC/MW: vertex is non-leaf terminal */

#define SDSTAR_BASE_UNSET  -1
#define SDSTAR_BASE_KILLED -2

#define STP_DELPSEUDO_MAXGRAD   7
#define STP_DELPSEUDO_MAXNEDGES 21


/* ((((edge) % 2) == 0) ? ((edge) + 1) : ((edge) - 1)) without branch */
#define flipedge(edge) ( ((edge) + 1) - 2 * ((edge) % 2) )
#define flipedge_Uint(edge) ( (((unsigned int) edge) + 1) - 2 * (((unsigned int) edge) % 2) )

#define CONNECT      0
#define UNKNOWN    (-1)
#define FARAWAY            1e15
#define BLOCKED            1e10              /**< used for temporarily blocking an edge */
#define BLOCKED_MINOR      (BLOCKED - 1.0)   /**< used for permanently blocking an edge;
                                              * different from BLOCKED because of weird prize sum in reduce_base.c */

#define EDGE_BLOCKED       0
#define EDGE_MODIFIABLE    1

#define MST_MODE   0
#define FSP_MODE   1
#define BSP_MODE   2

#define Is_term(a)         ((a) >= 0)
#define Is_pseudoTerm(a)   ((a) == STP_TERM_PSEUDO)
#define Is_nonleafTerm(a)  ((a) == STP_TERM_NONLEAF)
#define Is_anyTerm(a)      ((a) >= 0 || (a) == STP_TERM_PSEUDO || (a) == STP_TERM_NONLEAF )
#define Edge_anti(a) ((((a) % 2) > 0) ? (a) - 1 : (a) + 1)
#define Edge_even(a) ((((a) % 2) == 0) ? (a) : (a) - 1)


/* stp file format */
#define STP_FILE_MAGIC       0x33d32945
#define STP_FILE_VERSION_MAJOR   1
#define STP_FILE_VERSION_MINOR   0

typedef enum { FF_BEA, FF_STP, FF_PRB, FF_GRD } FILETYPE;


/** fixed graph components */
typedef struct fixed_graph_components FIXED;


/** node ancestors resulting from pseudo-elimination */
typedef struct pseudo_ancestors PSEUDOANS;


/** depository for several CSR storages */
typedef struct compressed_sparse_storages_depository CSRDEPO;

/** bottleneck Steiner distance (implied) profit */
typedef struct special_distance_implied_profit SDPROFIT;


/** helper needed for extracting and reinserting subgraph */
typedef struct subgraph_extraction_insertion SUBINOUT;


/** CSR like graph storage */
typedef struct csr_storage
{
   int*                  start;              /**< start position for each node */
   int*                  head;               /**< edge head array */
   SCIP_Real*            cost;               /**< edge cost array */
   int*                  edge_id;            /**< edge ids */
   int                   nedges_max;         /**< maximum number of edges (real number given by start[nnodes]) */
   int                   nnodes;             /**< number of nodes */
} CSR;


/** for dynamic CSR */
typedef struct csr_range
{
   int start;
   int end;
} RANGE;


/** dynamic CSR */
typedef struct dynamic_csr_storage
{
   RANGE*                range;              /**< CSR range */
   int*                  head;               /**< edge head array */
   int*                  edgeid;             /**< gets id from CSR edge */
   int*                  id2csredge;         /**< gets CRS edge from id */
   SCIP_Real*            cost;               /**< edge cost array */
   SCIP_Real*            cost2;              /**< second edge cost array, initialized with NULL and never freed! */
   SCIP_Real*            cost3;              /**< third edge cost array, initialized with NULL and never freed! */
   int                   nedges;             /**< number of edges */
   int                   nnodes;             /**< number of nodes */
} DCSR;


/** singleton ancestors for given undirected edge */
typedef struct singleton_ancestors_edge
{
   IDX*                  ancestors;          /**< ancestors */
   IDX*                  revancestors;       /**< reverse ancestors */
   int*                  pseudoancestors;    /**< pseudo ancestors */
   int                   npseudoancestors;   /**< number of pseudo ancestors */
   int                   edge;               /**< edge index */
} SINGLETONANS;


/** Steiner graph data structure */
typedef struct
{
   /* Nodes */
   int                   norgmodelknots;     /**< Number of nodes in the original model               */
   int                   norgmodelterms;     /**< Number of terminals in the original model           */
   int                   ksize;              /**< Count of allocated knot slots                       */
   int                   knots;              /**< Count of nodes in graph                             */
   int                   orgknots;           /**< Count of nodes prior to graph reduction             */
   int                   terms;              /**< Count of terminals                                  */
   int                   layers;             /**< Count of different networks                         */
   int                   orgsource;          /**< root of unreduced graph                             */
   int                   source;             /**< The root                                            */
   int* RESTRICT         term;               /**< Array [0..nodes-1] of networknumber for             */
                                             /**< knot [i], -1 if [i] is never a terminal             */
   int* RESTRICT         mark;               /**< Array [0..nodes-1], normally TRUE or FALSE          */
                                             /**< to mark nodes for inclusion in the shortest         */
                                             /**< path / minimum spanning tree routine                */
   int* RESTRICT         grad;               /**< Array [0..nodes-1] with degree of knot [i]          */
   int* RESTRICT         inpbeg;             /**< Array [0..nodes-1] with starting slot index         */
                                             /**< for the ieat array, -1 if not used                  */
   int* RESTRICT         outbeg;             /**< Array [0..nodes-1] with starting slot index         */
                                             /**< for the oeat array, -1 if not used                  */
   int* RESTRICT         maxdeg;             /**< For HCDSTP: Array [0..nodes-1] containing the maximal degrees
                                                   of all nodes               */
   int*                  term2edge;          /**< (R)PCSTP and (R)MWCSP: Array [0..nodes-1] of edge to twin terminal or -1 */

   SCIP_Real*            prize;              /**< For NWSTP, (R)PCSTP and (R)MWCSP: Array [0..nodes-1] of node costs       */
   SCIP_Real*            costbudget;         /**< budget cost value for (R)BMWCSP:  Array [0..nodes-1] */
   SCIP_Real             budget;             /**< budget value for (R)BMWCSP */

   /* Edges */
   int                   norgmodeledges;     /**< Count of edges prior to graph transformation                               */
   int                   hoplimit;           /**< maximal number of edges allowed for a solution to be feasible
                                               (only used for HCDSTPs)                                                      */
   int                   esize;              /**< Count of allocated edge slots                                             */
   int                   edges;              /**< Count of edges in the graph                                               */
   int                   orgedges;
   SCIP_Real*            cost;               /**< Array [0..edges-1] of positive edge costs                                 */
   SCIP_Real*            cost_org_pc;        /**< Array [0..edges-1] of positive edge costs for non-transformed PC/MW variants   */
   int* RESTRICT         tail;               /**< Array [0..edges-1] of node-number of tail of edge [i]                     */
   int* RESTRICT         head;               /**< Array [0..edges-1] of node-number of head of edge [i]                     */
   int* RESTRICT         orgtail;            /**< Array [0..edges-1] of node-number of tail of edge [i] prior to reduction  */
   int* RESTRICT         orghead;            /**< Array [0..edges-1] of node-number of head of edge [i] prior to reduction  */
   int* RESTRICT         rootedgeprevs;      /**< Array [0..edges-1] for PC and MW problems */

   /* Nodes/Edges */
   int* RESTRICT         ieat;               /**< Array [0..edges-1], incoming edge allocation table          */
   int* RESTRICT         oeat;               /**< Array [0..edges-1], outgoing edge allocation table          */

   /* History */
   IDX**                 ancestors;          /**< list of ancestor edges to each edge (to keep track of reductions)         */
   IDX**                 pcancestors;        /**< list of ancestor edges to each node (to keep track of PC/MW reductions )  */
   PSEUDOANS*            pseudoancestors;    /**< pseudo ancestors */
   FIXED*                fixedcomponents;    /**< fixed components */
   int*                  contracttrace;      /**< used to trace node contractions */

   /* Data for min cut computation todo extra struct */
   int                   mincut_nnodes;
   int                   mincut_nedges;
   int* RESTRICT         mincut_dist;        /**< dist[i] : Distance-label of node i          */
   int* RESTRICT         mincut_head;        /**< head[i] : Head of active queue with label i */
   int* RESTRICT         mincut_head_inact;  /**< head[i] : Head of inactive queue with label i */
   int* RESTRICT         mincut_numb;        /**< numb[i] : numb[i] nodes with label i        */
   int* RESTRICT         mincut_prev;
   int* RESTRICT         mincut_next;
   int* RESTRICT         mincut_temp;
   int* RESTRICT         mincut_e;           /**< e[i] : Excess of node i                     */
   int* RESTRICT         mincut_x;           /**< x[k] : Actual flow on arc k                 */
   int* RESTRICT         mincut_r;           /**< r[k] : residual capacity of arc k                    */

   /* Data for sp and mst computation */
   int* RESTRICT         path_heap;
   int* RESTRICT         path_state;

   /* Data for grid problems */
   int                   grid_dim;
   int*                  grid_ncoords;
   int**                 grid_coordinates;

   /* Global information */
   int                   stp_type;           /**< Steiner problem variant  */
   SCIP_Bool             is_packed;          /**< graph already packed?    */
   SCIP_Bool             withInexactReductions;
   SCIP_Bool             extended;           /**< For (R)PCSTP and (R)MWCSP: signifies whether problem is in extended
                                                 form (TRUE) or not (FALSE) */
   /* other adjacency storages */
   CSR*                  csr_storage;        /**< CSR structure or NULL */
   DCSR*                 dcsr_storage;       /**< Dynamic CSR structure or NULL */
} GRAPH;

typedef struct presolve_info
{
   SCIP_Real fixed;
   SCIP_Real upper;
   SCIP_Real lower;
   int    time;
} PRESOL;


/** ONE segment of a path */
typedef struct shortest_path
{
   SCIP_Real             dist;               /* Distance to the end of the path             */
   signed int            edge;               /* Incoming edge to go along                   */
} PATH;

/** Steiner nodes to terminal paths */
typedef struct nodes_to_terminal_paths TPATHS;

/** heap entry */
typedef struct dijkstra_heap_entry
{
   SCIP_Real             key;
   int                   node;
} DENTRY;

/** Dijkstra heap */
typedef struct dijkstra_heap
{
   int                   capacity;           /**< maximum size */
   int                   size;               /**< size */
   int*                  position;           /**< position of an index in range 0 to capacity */
   DENTRY*               entries;            /**< number of components  */
} DHEAP;


/** Dijkstra data */
typedef struct dijkstra_data
{
   /* temporary arrays: */
   int*                  visitlist;          /**< stores all visited nodes */
   DHEAP*                dheap;              /**< Dijkstra heap, initially cleaned */
   SCIP_Real*            node_distance;      /**< distances array for each node, initially set to FARAWAY */
   STP_Bool*             node_visited;       /**< stores whether a node has been visited, initially set to FALSE */
   /* long-term arrays: */
   SCIP_Real*            node_bias;          /**< bias of node for PC problem */
   int                   nvisits;            /**< number of visited nodes, initially set to -1 */
   int                   edgelimit;          /**< number of edges to consider */
} DIJK;


/** Voronoi data */
typedef struct voronoi_storage
{
   SCIP_Real*            nodes_dist;         /**< distance to base for each node */
   int*                  nodes_predEdge;     /**< predecessor edge (incoming) to each node */
   int*                  nodes_base;         /**< base of each node*/
   int                   nnodes;             /**< number of nodes */
   SCIP_Bool             usingBufferArrays;  /**< are buffer arrays being used? */
} VNOI;


/** Stores data for computation of special distance/bottleneck distance clique computations  */
typedef struct special_distance_clique
{
   DIJK*                 dijkdata;           /**< temporary data */
   SCIP_Real*            sds;                /**< array for SDs of clique */
   int*                  cliquenodes;        /**< nodes */
   const int*            cliqueToNodeMap;    /**< makes clique nodes to original nodes; NON-OWNED! */
   int                   centernode;         /**< center node or, if there is none, -1 */
   int                   ncliquenodes;       /**< number of nodes */
} SDCLIQUE;


#ifdef __cplusplus
}
#endif



#endif /* APPLICATIONS_STP_SRC_GRAPHDEFS_H_ */
