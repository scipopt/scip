/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grph.h
 * @brief  includes various files containing graph methods used for Steiner tree problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef _GRAPH_H_
#define _GRAPH_H_

#define EAT_FREE     -1
#define EAT_LAST     -2
#define EAT_HIDE     -3

#define GRAPH_HAS_COORDINATES     1
#define GRAPH_IS_GRIDGRAPH        2
#define GRAPH_IS_DIRECTED         4

#define STP_SPG                      0
#define STP_SAP                      1
#define STP_PCSPG                    2
#define STP_RPCSPG                   3
#define STP_NWSPG                    4
#define STP_DCSTP                    5
#define STP_REVENUES_BUDGET_HOPCONS  6
#define STP_RSMT                     7
#define STP_OARSMT                   8
#define STP_MWCSP                    9
#define STP_DHCSTP                   10
#define STP_GSTP                     11
#define STP_RMWCSP                   12

typedef unsigned char STP_Bool;

#include "scip/scip.h"
#include "misc_stp.h"

typedef struct
{
   /* Nodes */
   int                   norgmodelknots;     /**< Number of nodes in the original model               */
   int                   flags;              /**< To store attributes                                 */
   int                   ksize;              /**< Count of allocated knot slots                       */
   int                   knots;              /**< Count of nodes in graph                             */
   int                   orgknots;           /**< Count of nodes prior to graph reduction             */
   int                   terms;              /**< Count of terminals                                  */
   int                   layers;             /**< Count of different networks                         */
   int                   orgsource;          /**< root of unreduced graph                             */
   int*                  locals;             /**< Array [0..layers-1] of count of terminals           */
                                             /**< in network [i]                                      */
   int* RESTRICT         source;             /**< Array [0..layers-1] of knot number of the           */
                                             /**< root of network [i], -1 if unknown                  */
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

   SCIP_Real*            prize;              /**< For (R)PCSTP and (R)MWCSP: Array [0..nodes-1] of node costs       */

   /* Edges */
   IDX*                  fixedges;           /**< list of fixed edges*/
   IDX**                 ancestors;          /**< list of ancestor edges to each edge (to keep track of reductions)         */
   IDX**                 pcancestors;        /**< list of ancestor edges to each node (to keep track of PC/MW reductions )  */
   int                   norgmodeledges;     /**< Count of edges prior to graph transformation                               */
   int                   hoplimit;           /**< maximal number of edges allowed for a solution to be feasible
                                               (only used for HCDSTPs)                                                      */
   int                   esize;              /**< Count of allocated edge slots                                             */
   int                   edges;              /**< Count of edges in the graph                                               */
   int                   orgedges;
   SCIP_Real*            cost;               /**< Array [0..edges-1] of positive edge costs                                  */
   int* RESTRICT         tail;               /**< Array [0..edges-1] of node-number of tail of edge [i]                     */
   int* RESTRICT         head;               /**< Array [0..edges-1] of node-number of head of edge [i]                     */
   int* RESTRICT         orgtail;            /**< Array [0..edges-1] of node-number of tail of edge [i] prior to reduction  */
   int* RESTRICT         orghead;            /**< Array [0..edges-1] of node-number of head of edge [i] prior to reduction  */

   /* Nodes/Edges */
   int* RESTRICT         ieat;               /**< Array [0..edges-1], incoming edge allocation table          */
   int* RESTRICT         oeat;               /**< Array [0..edges-1], outgoing edge allocation table          */

   /* Data for min cut computation */
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
   int                   stp_type;          /**< Steiner problem variant                                                */
   SCIP_Bool             extended;          /**< For (R)PCSTP and (R)MWCSP: signifies whether problem is in extended
                                                 form (TRUE) or not (FALSE)                                             */

} GRAPH;

typedef struct presolve_info
{
   SCIP_Real fixed;
   SCIP_Real upper;
   SCIP_Real lower;
   int    time;
} PRESOL;

/* ONE segment of a path
 */
typedef struct shortest_path
{
   SCIP_Real             dist;               /* Distance to the end of the path             */
   signed int            edge;               /* Incoming edge to go along                   */
} PATH;

/* ((((edge) % 2) == 0) ? ((edge) + 1) : ((edge) - 1)) without branch */
#define flipedge(edge) ( ((edge) + 1) - 2 * ((edge) % 2) )

#define PATH_NIL    ((PATH*)0)

#define CONNECT      0
#define UNKNOWN    (-1)
#define FARAWAY      1e15
#define BLOCKED     1e10

#define EDGE_BLOCKED      0
#define EDGE_MODIFIABLE    1

#define MST_MODE   0
#define FSP_MODE   1
#define BSP_MODE   2

#define NO_CHANGE    -10

#define Is_term(a)   ((a) >= 0)
#define Is_pterm(a)  ((a) == -2)
#define Is_gterm(a)  ((a) == -2 || (a) >= 0 )
#define Edge_anti(a) ((((a) % 2) > 0) ? (a) - 1 : (a) + 1)

#define STP_MAGIC       0x33d32945
#define VERSION_MAJOR   1
#define VERSION_MINOR   0

typedef enum { FF_BEA, FF_STP, FF_PRB, FF_GRD } FILETYPE;

/* grphbase.c
 */
extern void   graph_flags(GRAPH*, int);
extern void   graph_show(const GRAPH*);
extern void   graph_ident(const GRAPH*);
extern void   graph_knot_add(GRAPH*, int);
extern void   graph_knot_chg(GRAPH*, int, int);
extern void   graph_knot_contract_dir(GRAPH*, int, int);
extern void   graph_get_csr(const GRAPH*, int* RESTRICT, int* RESTRICT, int* RESTRICT, int*);
extern void   graph_edge_add(SCIP*, GRAPH*, int, int, double, double);
extern void   graph_edge_del(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_edge_hide(GRAPH*, int);
extern void   graph_uncover(GRAPH*);
extern void   graph_subtractprize(SCIP*, GRAPH*, SCIP_Real, int);
extern void   graph_trail(const GRAPH*, int);
extern void   graph_free(SCIP*, GRAPH*, SCIP_Bool);
extern void   graph_getNVET(const GRAPH*, int*, int*, int*);
extern void   graph_listSetSolNode(const GRAPH*, STP_Bool*, IDX*);
extern SCIP_RETCODE   graph_trail_arr(SCIP*, const GRAPH*, int);
extern SCIP_RETCODE   graph_init(SCIP*, GRAPH**, int, int, int, int);
extern SCIP_RETCODE   graph_2org(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_2trans(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_2orgcheck(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_2transcheck(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_init_history(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_resize(SCIP*, GRAPH*, int, int, int);
extern SCIP_RETCODE   graph_knot_contract(SCIP*, GRAPH*, int*, int, int);
extern SCIP_RETCODE   graph_knot_contractpc(SCIP*, GRAPH*, int*, int, int, int);
extern SCIP_RETCODE   graph_RerootSol(SCIP*, GRAPH*, int*, int);
extern SCIP_RETCODE   graph_PcToSap(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_MwToRmw(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_PcSapCopy(SCIP*, GRAPH*, GRAPH**, SCIP_Real*);
extern SCIP_RETCODE   graph_PcSapCopyShift(SCIP*, GRAPH*, GRAPH**, SCIP_Real*);
extern SCIP_RETCODE   graph_PcRSapCopy(SCIP*, GRAPH*, GRAPH**, int*, int, int);
extern SCIP_RETCODE   graph_edge_reinsert(SCIP*, GRAPH*, int, int, int, SCIP_Real, IDX*, IDX*, IDX*, IDX*);
extern SCIP_RETCODE   graph_MwcsToSap(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE   graph_prize_transform(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_rootprize_transform(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_maxweight_transform(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE   graph_rootmaxweight_transform(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_grid_create(SCIP*, GRAPH**, int**, int, int, int);
extern SCIP_RETCODE   graph_obstgrid_create(SCIP*, GRAPH**, int**, int**, int, int, int, int);
extern SCIP_RETCODE   graph_grid_coordinates(SCIP*, int**, int**, int*, int, int);
extern SCIP_RETCODE   graph_copy(SCIP*, const GRAPH*, GRAPH**);
extern SCIP_RETCODE   graph_pack(SCIP*, GRAPH*, GRAPH**, SCIP_Bool);
extern int    graph_edge_redirect(SCIP*, GRAPH*, int, int, int, SCIP_Real);
extern int    graph_valid(const GRAPH*);
extern int    graph_pterm_delete(SCIP*, GRAPH*, int);
extern SCIP_Bool graph_sol_valid(SCIP*, const GRAPH*, const int*);
extern SCIP_Real graph_computeSolVal(const SCIP_Real*, const int*, SCIP_Real, int);

/* grphpath.c
 */
extern void   graph_path_exit(SCIP*, GRAPH*);
extern void   graph_path_exec(SCIP*, const GRAPH*, const int, int, const SCIP_Real*, PATH*);
extern void   graph_path_execX(SCIP*, const GRAPH*, int, const SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_invroot(SCIP*, const GRAPH*, int, const SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_st(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, int*, int, SCIP_RANDNUMGEN*, STP_Bool*);
extern void   graph_path_st_rmw(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, const STP_Bool*, STP_Bool*);
extern void   graph_path_st_pcmw(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, const STP_Bool*, STP_Bool*);
extern void   graph_path_st_pcmw_full(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_reduce(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_extend(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, STP_Bool*, SCIP_Bool*);
extern void   graph_path_st_rpc(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   voronoi(SCIP* scip, const GRAPH*, SCIP_Real*, SCIP_Real*, STP_Bool*, int*, PATH*);
extern void   get2next(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   get3next(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   get4next(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   getnext3terms(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   getnext4terms(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   getnext4tterms(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   voronoi_mw(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   voronoi_terms(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   voronoi_inout(const GRAPH*);
extern void   voronoi_term(const GRAPH*, double*, double*, double*, PATH*, int*, int*, int*, int*, int);
extern void   heap_add(int*, int*, int*, int, PATH*);
extern void   voronoi_repair(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, PATH*, int*, int, UF*);
extern void   voronoi_repair_mult(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, int*, int*, STP_Bool*, UF*, PATH*);
extern void   voronoiSteinerTreeExt(SCIP*, const GRAPH*, SCIP_Real*, int*, STP_Bool*, PATH*);
extern void   sdpaths(SCIP*, const GRAPH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int, int, int);
extern void   voronoi_mw_radius(SCIP* scip, const GRAPH*, PATH*, const SCIP_Real*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE   voronoi_extend(SCIP*, const GRAPH*, SCIP_Real*, PATH*, SCIP_Real**, int**, int**, STP_Bool*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE   graph_path_init(SCIP*, GRAPH*);
extern SCIP_RETCODE   voronoi_dist(SCIP*, const GRAPH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, int*, PATH*);
extern SCIP_RETCODE   voronoi_radius(SCIP* scip, const GRAPH*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*);



/* grphmcut.c
 */
#if 0
extern void   graph_mincut_exit(SCIP*, GRAPH*);
extern void   graph_mincut_exec(GRAPH*, int, int, int, const int*, int*, const int*, const int*, const int*, int);
extern SCIP_RETCODE   graph_mincut_init(SCIP*, GRAPH*);
#else
extern void   graph_mincut_exit(SCIP*, GRAPH*);
extern void   graph_mincut_exec(const GRAPH*, const int, const int, const int, const int, const int, const int*, const int*, int* RESTRICT, const int*, const int*, const int*, const SCIP_Bool);
extern SCIP_RETCODE   graph_mincut_init(SCIP*, GRAPH*);
#endif

/* grphload.c
 */
extern SCIP_RETCODE graph_load(SCIP*, GRAPH**, const char*, PRESOL*);

/* grphsave.c
 */
extern void graph_save(const GRAPH*, const char*, FILETYPE);
extern void SCIPwriteStp(SCIP*, const GRAPH*, FILE*, SCIP_Real);

/* grph2fig.c
 */
extern void graph_writefig(const GRAPH*, const char*, const double *, int);
extern void graph_bfscoord(GRAPH* g);
extern void graph_boxcoord(GRAPH* g);

/* reduce.c
 */
extern SCIP_RETCODE level0(SCIP*, GRAPH*);
extern SCIP_RETCODE reduceStp(SCIP*, GRAPH**, SCIP_Real*, int, SCIP_Bool, SCIP_Bool, int*, SCIP_Bool);
extern SCIP_RETCODE redLoopStp(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Real, SCIP_Bool, SCIP_Bool, SCIP_Bool, int, int*, SCIP_Bool);
extern SCIP_RETCODE redLoopPc(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Bool, SCIP_Bool, int, SCIP_Bool);
extern SCIP_RETCODE redLoopMw(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, STP_Bool, STP_Bool, STP_Bool, int, SCIP_Bool);
extern SCIP_RETCODE reduce(SCIP*, GRAPH**, SCIP_Real*, int, int, SCIP_Bool);

/* reduce_alt.c
 */

extern void    reduce_alt_ans(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    sdsp_reduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, int*);
extern SCIP_RETCODE    sdsp_sap_reduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    sd_red(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, SCIP_Bool, int*);
extern SCIP_RETCODE    sdpc_reduction(SCIP*, GRAPH*, PATH*, SCIP_Real*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    sd2_reduction(SCIP*, GRAPH*, SCIP_Real*, int*, int*);
extern SCIP_RETCODE    getSD(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int, int, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    bd3_reduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    bdr_reduction(SCIP*, GRAPH*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    nv_reduction(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    nv_reductionAdv(SCIP*, GRAPH*, PATH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    sl_reduction(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*, STP_Bool*, int*, int*);
extern SCIP_RETCODE    ledge_reduction(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    ansadvReduction(SCIP*, GRAPH*, SCIP_Real*, int*, int*);
extern SCIP_RETCODE    ansadv2Reduction(SCIP*, GRAPH*, SCIP_Real*, int*, int*);
extern SCIP_RETCODE    cnsAdvReduction(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    nnpReduction(SCIP*, GRAPH*, SCIP_Real*, int*, int*, int*, int*, int, STP_Bool*);
extern SCIP_RETCODE    npvReduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    chain2Reduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);

/* reduce_bnd.c
 */
extern SCIP_RETCODE    da_reduce(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*, int*, int, SCIP_RANDNUMGEN*, SCIP_Bool, int*);
extern SCIP_RETCODE    da_reduceSlackPrune(SCIP*, SCIP_VAR**, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*, STP_Bool*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    da_reduceSlackPruneMw(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, STP_Bool*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    da_reducePcMw(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, STP_Bool*, int*, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    bound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    bound_reduceMw(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
#if 0
extern SCIP_RETCODE    bound_reducePrune(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);
#endif
extern SCIP_RETCODE    bound_reducePrune(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    hopbound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    hcrbound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    hcrcbound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);

/* reduce_simple.c
 */
extern SCIP_RETCODE    contractZeroEdges(SCIP*, GRAPH*);
extern SCIP_RETCODE    degree_test(SCIP*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE    degree_test_hc(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    degree_test_mw(SCIP*, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    degree_test_pc(SCIP*, GRAPH*, SCIP_Real*, int*, int*, SCIP_Bool);
extern SCIP_RETCODE    degree_test_sap(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    rptReduction(SCIP*, GRAPH*, SCIP_Real*, int*);

/* validate.c
 */
extern SCIP_RETCODE    SCIPvalidateStpSol(SCIP*, const GRAPH*, const double*, SCIP_Bool*);

#endif /* !_GRAPH_H_ */
