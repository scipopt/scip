/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grph.h
 * @brief  includes various files containing graph methods used for Steiner problems
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

#define STP_UNDIRECTED               0
#define STP_DIRECTED                 1
#define STP_PRIZE_COLLECTING         2
#define STP_ROOTED_PRIZE_COLLECTING  3
#define STP_NODE_WEIGHTS             4
#define STP_DEG_CONS                 5
#define STP_REVENUES_BUDGET_HOPCONS  6
#define STP_GRID                     7
#define STP_OBSTACLES_GRID           8
#define STP_MAX_NODE_WEIGHT          9
#define STP_HOP_CONS                 10
#define GSTP                         11
#define PRIZE                        12

#include "scip/scip.h"
#include "misc_stp.h"

typedef struct
{
   /* Knots
    */
   int     norgmodelknots;
   int     flags;  /**< To store attributes                         */
   int     ksize;  /**< Count of allocated knot slots               */
   int     knots;  /**< Count of knots in graph                     */
   int     orgknots;
   int     terms;  /**< Count of terminals                          */
   int     layers; /**< Count of different networks                 */
   int*    locals; /**< Array [0..layers-1] of count of terminals   */
                   /**< in network [i]                              */
   int*    source; /**< Array [0..layers-1] of knot number of the   */
                   /**< root of network [i], -1 if unknown          */
   int*    term;   /**< Array [0..knots-1] of networknumber for     */
                   /**< knot [i], -1 if [i] is never a terminal     */
   int*    mark;   /**< Array [0..knots-1], normaly TRUE or FALSE   */
                   /**< to mark knots for inclusion in the shortest */
                   /**< path / minimum spanning tree routine        */
   int*    grad;   /**< Array [0..knots-1] with degree of knot [i]  */
   int*    inpbeg; /**< Array [0..knots-1] with starting slot index */
                   /**< for the ieat array, -1 if not used          */
   int*    outbeg; /**< Array [0..knots-1] with starting slot index */
                   /**< for the oeat array, -1 if not used          */
   int*    maxdeg; /**< Array [0..knots-1] containing the maximal
                      degrees of all nodes (only used for Degree-
                      Constraint STPs)                            */
   /* Edges
    */
   IDX*    fixedges;  /**< list of fixed edges*/
   IDX**   ancestors; /**< list of ancestor edges to each edge (required to keep track of reductions) */
   IDX**   pcancestors; /**< list of ancestor edges to each node (required to keep track of reductions in PC) */
   int     norgmodeledges;
   int     hoplimit;  /**< maximal number of edges allowed for a solution to be feasible
                         (only for problem type STP_HOP_CONS) */
   int     esize;  /**< Count of allocated edge slots               */
   int     edges;  /**< Count of edges in the graph                 */
   int     orgedges;
   SCIP_Real* cost;   /**< Array [0..edges-1] of positiv edge costs  */
   SCIP_Real* prize;   /**< Array [0..edges-1] of positiv node costs */
   int*    tail;   /**< Array [0..edges-1] of knot-number of tail   */
                   /**< of edge [i]                                 */
   int*    head;   /**< Array [0..edges-1] of knot-number of head   */
                   /**< of edge [i]                                 */
   int*    orgtail;
   int*    orghead;
   /* Knots/Edges
    */
   int*    ieat;   /**< Array [0..edges-1], incomming edge          */
                   /**< allocation table                            */
   int*    oeat;   /**< Array [0..edges-1], outgoing edge           */
                   /**< allocation table                            */

   /* data for min cut computation
    */
   int*    mincut_dist;    /**< dist[i] : Distance-label of Knot i          */
   int*    mincut_head;    /**< head[i] : Head of Active Queue with Label i */
   int*    mincut_numb;    /**< numb[i] : numb[i] Knots with Label i        */
   int*    mincut_prev;
   int*    mincut_next;
   int*    mincut_temp;
   int*    mincut_e;       /**< e[i] : Excess of Knot i                     */
   int*    mincut_x;       /**< x[k] : Actual Flow on Arc k                 */
   int*    mincut_r;       /**< r[k] : Capacity of Arc k                    */
   /* data for math and mst computation
    */
   int* path_heap;
   int* path_state;

   /* data for grid problems
    */
   int       grid_dim;
   int*      grid_ncoords;
   int**     grid_coordinates;

   /* type of Steiner problem  */
   int stp_type;

} GRAPH;

typedef struct presolve_info
{
   double fixed;
   double upper;
   double lower;
   int    time;
} PRESOL;

/* ONE segment of a path
 */
typedef struct shortest_path
{
   double       dist;         /* Distance to the end of the path             */
   signed int   edge;         /* First edge to go                            */
} PATH;


#define flipedge(edge) (((edge % 2) == 0) ? edge + 1 : edge - 1)

#define PATH_NIL    ((PATH*)0)

#define CONNECT      0
#define UNKNOWN    (-1)
#define FARAWAY      1e15
#define BLOCKED     1e10

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

extern SCIP_RETCODE graph_init(SCIP*, GRAPH**, int, int, int, int);
extern SCIP_RETCODE   pcgraphorg(SCIP*, GRAPH*);
extern SCIP_RETCODE   pcgraphtrans(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_init_history(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_resize(SCIP*, GRAPH*, int, int, int);
extern void   graph_free(SCIP*, GRAPH*, SCIP_Bool);
extern SCIP_RETCODE   graph_prize_transform(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_rootprize_transform(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_maxweight_transform(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE graph_grid_create(SCIP*, GRAPH**, int**, int, int, int);
extern SCIP_RETCODE graph_obstgrid_create(SCIP*, GRAPH**, int**, int**, int, int, int, int);
extern SCIP_RETCODE   graph_grid_coordinates(SCIP*, int**, int**, int*, int, int);
#if 0
extern GRAPH* graph_copy(const GRAPH*);
extern GRAPH* graph_init2(int, int, int, int);
extern GRAPH* graph_pack2(SCIP*, GRAPH*, SCIP_Bool);
#endif
extern SCIP_RETCODE graph_copy(SCIP*, GRAPH*, GRAPH**);
extern void   graph_flags(GRAPH*, int);
extern void   graph_show(const GRAPH*);
extern void   graph_ident(const GRAPH*);
extern void   graph_knot_add(GRAPH*, int);
extern void   graph_knot_chg(GRAPH*, int, int);
extern SCIP_RETCODE   graph_knot_contract(SCIP*, GRAPH*, int, int);
extern SCIP_RETCODE   graph_knot_contractpc(SCIP*, GRAPH*, int, int, int);
extern void   graph_knot_contract_dir(GRAPH*, int, int);
extern void   graph_edge_add(SCIP*, GRAPH*, int, int, double, double);
extern void   graph_edge_del(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_edge_hide(GRAPH*, int);
extern int    graph_edge_redirect(SCIP*, GRAPH*, int, int, int, SCIP_Real);
extern SCIP_RETCODE   graph_edge_reinsert(SCIP*, GRAPH*, int, int, int, SCIP_Real, IDX*, IDX*, IDX*, IDX*);
extern void   graph_uncover(GRAPH*);
extern void   prize_subtract(SCIP*, GRAPH*, SCIP_Real, int);
extern SCIP_RETCODE  graph_pack(SCIP*, GRAPH*, GRAPH**, SCIP_Bool);
extern void   graph_trail(const GRAPH*, int);
extern int    graph_valid(const GRAPH*);
extern char    graph_valid2(SCIP*, const GRAPH*, SCIP_Real*);
extern SCIP_Bool graph_sol_valid(SCIP*, const GRAPH*, int*);

/* grphpath.c
 */
extern SCIP_RETCODE   graph_path_init(SCIP*, GRAPH*);
extern void   graph_path_exit(SCIP*, GRAPH*);
extern void   graph_path_exec(SCIP*, const GRAPH*, int, int, SCIP_Real*, PATH*);
extern void   graph_path_execX(SCIP*, const GRAPH*, int, SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_st(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, int*, int, char*);
#if 0
extern void   calculate_distances(SCIP*, const GRAPH*, PATH**, double*, int);
#endif
extern void   voronoi(SCIP* scip, const GRAPH*, SCIP_Real*, SCIP_Real*, char*, int*, PATH*);
extern void   get2next(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   get3next(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   getnext3terms(SCIP*, const GRAPH*, SCIP_Real*, SCIP_Real*, PATH*, int*, int*, int*);
extern void   voronoi_terms(SCIP*, const GRAPH*, SCIP_Real*, PATH*, int*, int*, int*);
extern SCIP_RETCODE   voronoi_dist(SCIP*, const GRAPH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, PATH*);
extern SCIP_RETCODE   voronoi_radius(SCIP* scip, const GRAPH*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*);
extern void   voronoi_inout(const GRAPH*);
extern void   voronoi_term(const GRAPH*, double*, double*, double*, PATH*, int*, int*, int*, int*, int);
extern void   voronoi_hop(const GRAPH*, double*, double*, double*, PATH*, int*, int*, int*, int*, int*);
extern void   heap_add(int*, int*, int*, int, PATH*);
extern void   voronoi_repair(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, PATH*, int*, int, UF*);
extern void   voronoi_repair_mult(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, int*, int*, char*, UF*, PATH*);
extern void   voronoi_slrepair(SCIP*, const GRAPH*, SCIP_Real*, PATH*, int*, int*, int*, int, int);
extern SCIP_RETCODE voronoi_extend(SCIP*, const GRAPH*, SCIP_Real*, PATH*, VLIST**, char*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE voronoi_extend2(SCIP*, const GRAPH*, SCIP_Real*, PATH*, SCIP_Real**, int**, int**, char*, int*, int*, int*, int, int, int);
extern void   sdpaths(SCIP*, const GRAPH*, PATH*, SCIP_Real*, int*, int*, int*, int*, int, int, int);
extern void   graph_path_length(const GRAPH*, const PATH*);

/* grphmcut.c
 */
extern SCIP_RETCODE   graph_mincut_init(SCIP*, GRAPH*);
extern void   graph_mincut_exit(SCIP*, GRAPH*);
extern void   graph_mincut_exec(GRAPH*, int, int, const int*, int*, int);

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
extern void level0(SCIP*, GRAPH*);
extern SCIP_RETCODE reduce(SCIP*, GRAPH**, SCIP_Real*, int, int);
extern SCIP_RETCODE bound_reduce(SCIP*, GRAPH*, PATH*, double*, double*, double*, double*, int*, int*, int*, int*);
extern SCIP_RETCODE hopbound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE hcrbound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE hcrcbound_reduce(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);

/* sdtest.c
 */
extern SCIP_RETCODE    sd_reduction(SCIP*, GRAPH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int, unsigned int*);

extern SCIP_RETCODE    sdsp_reduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    sd_red(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    sdpc_reduction(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    sd2_reduction(SCIP*, GRAPH*, SCIP_Real*, int*, int*);
extern SCIP_RETCODE    getSD(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, int*, int*, int*, int*, int*, int, int, int, SCIP_Bool);
extern SCIP_RETCODE    bd3_reduction(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    nv_reduction(SCIP*, GRAPH*, PATH*,double*, int*, int*, int*, int*);
extern SCIP_RETCODE    sl_reduction(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*);
extern SCIP_RETCODE    ledge_reduction(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*);
#if 0
extern SCIP_RETCODE    sd_reduction_dir(SCIP*, GRAPH*, double**, double**, double**, double**, double*, int*, int*, int*, int*);
extern SCIP_RETCODE    nsv_reduction(SCIP*, GRAPH*, SCIP_Real*, SCIP_Real*, int*);
extern SCIP_RETCODE    nv_reduction_optimal(SCIP*, GRAPH*, double*, int*, int);
#endif

/* dirreduce.c
 */
extern SCIP_RETCODE degree_test_dir(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE degree_test_pc(SCIP*, GRAPH*, SCIP_Real*, int*);

/* validate.c
 */
extern SCIP_RETCODE    SCIPvalidateStpSol(SCIP*, const GRAPH*, const double*, SCIP_Bool*);

#endif /* !_GRAPH_H_ */
