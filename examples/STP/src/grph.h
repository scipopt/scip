/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Header                                                        */
/*   File....: grph.h                                                        */
/*   Name....: Graph Routines                                                */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#define EAT_FREE     -1
#define EAT_LAST     -2
#define EAT_HIDE     -3

#define GRAPH_HAS_COORDINATES     1
#define GRAPH_IS_GRIDGRAPH        2
#define GRAPH_IS_DIRECTED         4

#define STP_UNDIRECTED              0
#define STP_DIRECTED                1
#define STP_PRIZE_COLLECTING        2
#define STP_NODE_WEIGHTS            3
#define STP_DEG_CONS                4
#define STP_REVENUES_BUDGET_HOPCONS 5
#define STP_GRID                    6
#define STP_MAX_NODE_WEIGHT         7

#include "scip/scip.h"
#include "misc_stp.h"
typedef struct
{
   /* Knots
    */
   int     flags;  /* To store attributes                         */
   int     ksize;  /* Count of allocated knot slots               */
   int     knots;  /* Count of knots in graph                     */
   int     terms;  /* Count of terminals                          */
   int     layers; /* Count of different networks                 */
   int*    locals; /* Array [0..layers-1] of count of terminals   */
                   /* in network [i]                              */
   int*    source; /* Array [0..layers-1] of knot number of the   */
                   /* root of network [i], -1 if unknown          */
   int*    term;   /* Array [0..knots-1] of networknumber for     */
                   /* knot [i], -1 if [i] is never a terminal     */
   double* prize;  /* Array [0..terms-1] of terminal prizes	  */
   int*    mark;   /* Array [0..knots-1], normaly TRUE or FALSE   */
                   /* to mark knots for inclusion in the shortest */
                   /* path / minimum spanning tree routine        */
   int*    grad;   /* Array [0..knots-1] with degree of knot [i]  */
   int*    inpbeg; /* Array [0..knots-1] with starting slot index */
                   /* for the ieat array, -1 if not used          */
   int*    outbeg; /* Array [0..knots-1] with starting slot index */
                   /* for the oeat array, -1 if not used          */
   int*    maxdeg; /* Array [0..knots-1] containing the maximal
                      degrees of all nodes (only used for Degree-
                      Constraint STPs)                            */
   /* Edges
    */
   int     esize;  /* Count of allocated edge slots               */
   int     edges;  /* Count of edges in the graph                 */
   SCIP_Real* cost;   /* Array [0..edges-1] of positiv edge costs    */
   int*    tail;   /* Array [0..edges-1] of knot-number of tail   */
                   /* of edge [i]                                 */
   int*    head;   /* Array [0..edges-1] of knot-number of head   */
                   /* of edge [i]                                 */
   /* Knots/Edges
    */
   int*    ieat;   /* Array [0..edges-1], incomming edge          */
                   /* allocation table                            */
   int*    oeat;   /* Array [0..edges-1], outgoing edge           */
                   /* allocation table                            */

   int*    xpos;   /* Array [0..knots-1] with X Coordinates of    */
                   /* Knot [i]                                    */
   int*    ypos;   /* Array [0..knots-1] with Y Coordinates of    */
                   /* Knot [i]                                    */
   /* data for min cut computation
    */
   int*    mincut_dist;    /* dist[i] : Distance-label of Knot i          */
   int*    mincut_head;    /* head[i] : Head of Active Queue with Label i */
   int*    mincut_numb;    /* numb[i] : numb[i] Knots with Label i        */
   int*    mincut_prev;
   int*    mincut_next;
   int*    mincut_temp;
   int*    mincut_e;       /* e[i] : Excess of Knot i                     */
   int*    mincut_x;       /* x[k] : Actual Flow on Arc k                 */
   int*    mincut_r;       /* r[k] : Capacity of Arc k                    */
   /* data for math and mst computation
    */
   int* path_heap;
   int* path_state;

   /* data for grid problems
    */
   int       grid_dim;
   int*      grid_ncoords;
   int**     grid_coordinates;

   /* type of Steiner Tree Problem  */
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

/* ONE segment of a voronoi path
 */
typedef struct voronoi_path
{
   double       dist;         /* Distance to the end of the path             */
   signed int   edge;         /* First edge to go                            */
   signed int   base;         /* Voronoi base                            */
} VNOI;

#define flipedge(edge) (((edge % 2) == 0) ? edge + 1 : edge - 1)

#define PATH_NIL    ((PATH*)0)

#define CONNECT      0
#define UNKNOWN    (-1)
#define FARAWAY      1e15

#define MST_MODE   0
#define FSP_MODE   1

#define NO_CHANGE    -10

#define Is_term(a)   ((a) >= 0)
#define Edge_anti(a) ((((a) % 2) > 0) ? (a) - 1 : (a) + 1)

#define STP_MAGIC       0x33d32945
#define VERSION_MAJOR   1
#define VERSION_MINOR   0

typedef enum { FF_BEA, FF_STP, FF_PRB, FF_GRD } FILETYPE;

/* grphbase.c
 */
extern GRAPH *graph_init(int, int, int, int);
extern void   graph_resize(GRAPH*, int, int, int);
extern void   graph_free(GRAPH*);
extern void   graph_prize_transform(GRAPH*);
extern void   graph_maxweight_transform(GRAPH*, double*);
extern GRAPH* graph_grid_create(int**, int, int, int);
extern void   graph_grid_coordinates(int**, int**, int*, int, int);
extern GRAPH* graph_copy(const GRAPH*);
extern void   graph_flags(GRAPH*, int);
extern void   graph_show(const GRAPH*);
extern void   graph_ident(const GRAPH*);
extern void   graph_knot_add(GRAPH*, int, int, int);
extern void   graph_knot_chg(GRAPH*, int, int, int, int);
extern void   graph_knot_contract(GRAPH*, int, int);
extern void   graph_edge_add(GRAPH*, int, int, double, double);
extern void   graph_edge_del(GRAPH*, int);
extern void   graph_edge_hide(GRAPH*, int);
extern void   graph_uncover(GRAPH*);
extern GRAPH* graph_pack(GRAPH*);
extern void   graph_trail(const GRAPH*, int i);
extern int    graph_valid(const GRAPH*);
extern char    graph_sol_valid(const GRAPH*, int*);

/* grphpath.c
 */
extern void   graph_path_init(GRAPH*);
extern void   graph_path_exit(GRAPH*);
extern void   graph_path_exec(const GRAPH*, int, int, SCIP_Real*, PATH*);
extern void   graph_path_exec2(const GRAPH*, int, int, const double*, PATH*, char*, int*, int*);
extern void   voronoi(SCIP* scip, const GRAPH*, SCIP_Real*, SCIP_Real*, char*, int*, PATH*);
extern void   heap_add(int*, int*, int*, int, PATH*);
extern void   voronoi_repair(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, PATH*, int*, int, UF*);
extern void   voronoi_repair_mult(SCIP*, const GRAPH*, SCIP_Real*, int*, int*, int*, int*, char*, UF*, PATH*);
extern SCIP_RETCODE  voronoi_extend(SCIP*, const GRAPH*, SCIP_Real*, PATH*, VLIST**, char*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE  voronoi_extend2(SCIP*, const GRAPH*, SCIP_Real*, PATH*, SCIP_Real**, int**, int**, char*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE  voronoi_extend3(SCIP*, const GRAPH*, const double*, PATH*, GNODE***, int**, int**, char*, int*, int*, int*, int, int, int);



extern void   graph_path_length(const GRAPH*, const PATH*);

/* grphmcut.c
 */
extern void   graph_mincut_init(GRAPH*);
extern void   graph_mincut_exit(GRAPH*);
extern void   graph_mincut_exec(GRAPH*, int, int, const int*, int*, int);

/* grphload.c
 */
extern GRAPH* graph_load(const char*, PRESOL*);

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
extern double reduce(GRAPH*, int, SCIP*);

/* sdtest.c
 */
extern int    sd_reduction(GRAPH*);
extern int    bd3_reduction(GRAPH*);
extern int    nsv_reduction(GRAPH*, double*);

/* validate.c
 */
extern int    validate(const GRAPH*, const double*);

#endif /* !_GRAPH_H_ */
