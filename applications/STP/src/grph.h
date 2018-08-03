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
   int* RESTRICT         rootedgeprevs;      /**< Array [0..edges-1] for PC and MW problems */

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
#define flipedge_Uint(edge) ( (((unsigned int) edge) + 1) - 2 * (((unsigned int) edge) % 2) )

#define PATH_NIL    ((PATH*)0)
#define CONNECT      0
#define UNKNOWN    (-1)
#define FARAWAY      1e15
#define BLOCKED      1e10

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

#define VERSION_SCIPJACK "1.3"
#define STP_MAGIC       0x33d32945
#define VERSION_MAJOR   1
#define VERSION_MINOR   0

typedef enum { FF_BEA, FF_STP, FF_PRB, FF_GRD } FILETYPE;

/* grphbase.c
 */
extern void   graph_pc_knot2nonTerm(GRAPH*, int);
extern void   graph_pc_updateTerm2edge(GRAPH*, const GRAPH*, int, int, int, int);
extern void   graph_pc_subtractPrize(SCIP*, GRAPH*, SCIP_Real, int);
extern void   graph_pc_chgPrize(SCIP*, GRAPH*, SCIP_Real, int);
extern void   graph_show(const GRAPH*);
extern void   graph_knot_add(GRAPH*, int);
extern void   graph_knot_chg(GRAPH*, int, int);
extern void   graph_knot_del(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_knot_contract_dir(GRAPH*, int, int);
extern void   graph_get_csr(const GRAPH*, int* RESTRICT, int* RESTRICT, int* RESTRICT, int*);
extern void   graph_edge_add(SCIP*, GRAPH*, int, int, double, double);
extern void   graph_edge_del(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_edge_hide(GRAPH*, int);
extern void   graph_edge_printInfo(SCIP*, const GRAPH*, int);
extern void   graph_uncover(GRAPH*);
extern void   graph_trail(const GRAPH*, int);
extern void   graph_free(SCIP*, GRAPH**, SCIP_Bool);
extern void   graph_free_history(SCIP*, GRAPH*);
extern void   graph_free_historyDeep(SCIP*, GRAPH*);
extern void   graph_get_NVET(const GRAPH*, int*, int*, int*);
extern void   graph_sol_setNodeList(const GRAPH*, STP_Bool*, IDX*);
extern void   graph_pc_2org(GRAPH*);
extern void   graph_pc_2trans(GRAPH*);
extern void   graph_pc_2orgcheck(GRAPH*);
extern void   graph_pc_2transcheck(GRAPH*);
extern void   graph_pc_adaptSap(SCIP*, SCIP_Real, GRAPH*, SCIP_Real*);
extern void   graph_pc_presolExit(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_init(SCIP*, GRAPH*, int, int);
extern SCIP_RETCODE   graph_pc_2pc(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_2rpc(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_2mw(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE   graph_pc_2rmw(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_mw2rmw(SCIP*, GRAPH*, SCIP_Real);
extern SCIP_RETCODE   graph_pc_getSap(SCIP*, GRAPH*, GRAPH**, SCIP_Real*);
extern SCIP_RETCODE   graph_pc_getSapShift(SCIP*, GRAPH*, GRAPH**, SCIP_Real*);
extern SCIP_RETCODE   graph_pc_getRsap(SCIP*, GRAPH*, GRAPH**, int*, int, int);
extern SCIP_RETCODE   graph_pc_contractEdgeAncestors(SCIP*, GRAPH*, int, int, int);
extern SCIP_RETCODE   graph_pc_contractEdge(SCIP*, GRAPH*, int*, int, int, int);
extern SCIP_RETCODE   graph_resize(SCIP*, GRAPH*, int, int, int);
extern SCIP_RETCODE   graph_knot_contract(SCIP*, GRAPH*, int*, int, int);
extern SCIP_RETCODE   graph_knot_contractLowdeg2High(SCIP*, GRAPH*, int*, int, int);
extern SCIP_RETCODE   graph_sol_reroot(SCIP*, GRAPH*, int*, int);
extern SCIP_RETCODE   graph_sol_getOrg(SCIP*, const GRAPH*, const GRAPH*, const int*, int*);
extern SCIP_RETCODE   graph_sol_markPcancestors(SCIP*, IDX**, const int*, const int*, int, STP_Bool*, STP_Bool*, int*, int*, int*);
extern SCIP_RETCODE   graph_edge_reinsert(SCIP*, GRAPH*, int, int, int, SCIP_Real, IDX*, IDX*, IDX*, IDX*, SCIP_Bool);
extern SCIP_RETCODE   graph_knot_delPseudo(SCIP*, GRAPH*, const SCIP_Real*, const SCIP_Real*, const SCIP_Real*, int, SCIP_Bool*);
extern SCIP_RETCODE   graph_grid_create(SCIP*, GRAPH**, int**, int, int, int);
extern SCIP_RETCODE   graph_obstgrid_create(SCIP*, GRAPH**, int**, int**, int, int, int, int);
extern SCIP_RETCODE   graph_grid_coordinates(SCIP*, int**, int**, int*, int, int);
extern SCIP_RETCODE   graph_copy_data(SCIP*, const GRAPH*, GRAPH*);
extern SCIP_RETCODE   graph_copy(SCIP*, const GRAPH*, GRAPH**);
extern SCIP_RETCODE   graph_pack(SCIP*, GRAPH*, GRAPH**, SCIP_Bool);
extern SCIP_RETCODE   graph_trail_arr(SCIP*, const GRAPH*, int);
extern SCIP_RETCODE   graph_get_edgeConflicts(SCIP*, const GRAPH*);
extern SCIP_RETCODE   graph_init(SCIP*, GRAPH**, int, int, int);
extern SCIP_RETCODE   graph_init_history(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_termsReachable(SCIP*, const GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pc_presolInit(SCIP*, GRAPH*);
extern int    graph_edge_redirect(SCIP*, GRAPH*, int, int, int, SCIP_Real, SCIP_Bool);
extern int    graph_pc_deleteTerm(SCIP*, GRAPH*, int);
extern SCIP_Bool graph_valid(const GRAPH*);
extern SCIP_Bool graph_pc_term2edgeConsistent(const GRAPH*);
extern SCIP_Bool graph_sol_unreduced(SCIP*, const GRAPH*, const int*);
extern SCIP_Bool graph_sol_valid(SCIP*, const GRAPH*, const int*);
extern SCIP_Bool graph_pc_isPcMw(const GRAPH*);
extern SCIP_Real graph_sol_getObj(const SCIP_Real*, const int*, SCIP_Real, int);
extern SCIP_Real graph_pc_getPosPrizeSum(SCIP*, const GRAPH*);

/* grphpath.c
 */
extern void   graph_path_exit(SCIP*, GRAPH*);
extern void   graph_path_exec(SCIP*, const GRAPH*, const int, int, const SCIP_Real*, PATH*);
extern void   graph_path_execX(SCIP*, const GRAPH*, int, const SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_invroot(SCIP*, const GRAPH*, int, const SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_st(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, int*, int, SCIP_RANDNUMGEN*, STP_Bool*);
extern void   graph_path_st_rmw(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_full(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_reduce(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_extend(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, STP_Bool*, SCIP_Bool*);
extern void   graph_path_st_rpc(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_voronoi(SCIP* scip, const GRAPH*, SCIP_Real*, SCIP_Real*, STP_Bool*, int*, PATH*);
extern void   graph_get2next(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_get3next(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_get4next(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_get3nextTerms(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_get4nextTerms(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_voronoiMw(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_voronoiTerms(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   voronoi_inout(const GRAPH*);
extern void   voronoi_term(const GRAPH*, double*, double*, double*, PATH*, int*, int*, int*, int*, int);
extern void   heap_add(int*, int*, int*, int, PATH*);
extern void   graph_voronoiRepair(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, PATH*, int*, int, UF*);
extern void   graph_voronoiRepairMult(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, int*, int*, STP_Bool*, UF*, PATH*);
extern void   voronoiSteinerTreeExt(SCIP*, const GRAPH*, SCIP_Real*, int*, STP_Bool*, PATH*);
extern void   graph_sdPaths(SCIP*, const GRAPH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int, int, int);
extern void   graph_path_PcMwSd(SCIP*, const GRAPH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int*, int, int, int);
extern void   graph_voronoiWithRadiusMw(SCIP* scip, const GRAPH*, PATH*, const SCIP_Real*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE   graph_voronoiExtend(SCIP*, const GRAPH*, SCIP_Real*, PATH*, SCIP_Real**, int**, int**, STP_Bool*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE   graph_path_init(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_voronoiWithDist(SCIP*, const GRAPH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, int*, PATH*);
extern SCIP_RETCODE   graph_voronoiWithRadius(SCIP* scip, const GRAPH*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE   graph_get4nextTTerms(SCIP*, const GRAPH*, SCIP_Real*, PATH*, int*, int*, int*);

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

/* reduce.c
 */
extern SCIP_RETCODE level0(SCIP*, GRAPH*);
extern SCIP_RETCODE level0save(SCIP*, GRAPH*);
extern SCIP_RETCODE reduceStp(SCIP*, GRAPH**, SCIP_Real*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE redLoopStp(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Real, SCIP_Bool, SCIP_Bool, SCIP_Bool, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE redLoopPc(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Bool, SCIP_Bool, int, SCIP_Bool);
extern SCIP_RETCODE redLoopMw(SCIP*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, STP_Bool, STP_Bool, STP_Bool, int, SCIP_Bool);
extern SCIP_RETCODE reduce(SCIP*, GRAPH**, SCIP_Real*, int, int, SCIP_Bool);

/* reduce_alt.c
 */
extern void    reduce_ans(SCIP*, GRAPH*, int*, int*);
extern void    reduce_ansAdv(SCIP*, GRAPH*, int*, int*, SCIP_Bool);
extern void    reduce_ansAdv2(SCIP*, GRAPH*, int*, int*);
extern void    reduce_nnp(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    reduce_sdsp(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, int*);
extern SCIP_RETCODE    reduce_sdspSap(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_sd(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, SCIP_Bool, int*);
extern SCIP_RETCODE    reduce_sdPc(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_getSd(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int, int, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    reduce_getSdPcMw(SCIP*, const GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE    reduce_nts(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_bd34(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_bdr(SCIP*, GRAPH*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_nv(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_nvAdv(SCIP*, GRAPH*, PATH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_sl(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*, STP_Bool*, int*, int*);
extern SCIP_RETCODE    reduce_ledge(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_cnsAdv(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    reduce_npv(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_chain2(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);

/* reduce_bnd.c
 */
extern SCIP_RETCODE    reduce_deleteConflictEdges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_check3Tree(SCIP*, const GRAPH*, int, const SCIP_Real*, const SCIP_Real*, const PATH*, const int*, SCIP_Real, const int*, int, int*, SCIP_Real*, SCIP_Bool*, unsigned int*, int*, SCIP_Bool*);
extern SCIP_RETCODE    reduce_da(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*, int*, int, SCIP_RANDNUMGEN*, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    reduce_daSlackPrune(SCIP*, SCIP_VAR**, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*, STP_Bool*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    reduce_daSlackPruneMw(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, STP_Bool*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    reduce_daPcMw(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, STP_Bool*, int*, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_RANDNUMGEN*, SCIP_Real);
extern SCIP_RETCODE    reduce_bound(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundMw(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundPrune(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, const int*, const int*, int*, int);
extern SCIP_RETCODE    reduce_boundHop(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopR(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopRc(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);
extern int reduce_extendedEdge(SCIP*, GRAPH*, const PATH*, const SCIP_Real*, const SCIP_Real*, const int*, SCIP_Real, int, int*, STP_Bool*);

/* reduce_simple.c
 */
extern SCIP_RETCODE    reduce_contractZeroEdges(SCIP*, GRAPH*, SCIP_Bool);
extern SCIP_RETCODE    reduce_simple(SCIP*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE    reduce_simple_hc(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_mw(SCIP*, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_pc(SCIP*, GRAPH*, SCIP_Real*, int*, int*, SCIP_Bool);
extern SCIP_RETCODE    reduce_simple_sap(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_rpt(SCIP*, GRAPH*, SCIP_Real*, int*);

/* validate.c
 */
extern SCIP_RETCODE    SCIPStpValidateSol(SCIP*, const GRAPH*, const double*, SCIP_Bool*);

#endif /* !_GRAPH_H_ */
