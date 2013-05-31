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
   int*    mark;   /* Array [0..knots-1], normaly TRUE or FALSE   */
                   /* to mark knots for inclusion in the shortest */
                   /* path / minimum spanning tree routine        */
   int*    grad;   /* Array [0..knots-1] with degree of knot [i]  */
   int*    inpbeg; /* Array [0..knots-1] with starting slot index */
                   /* for the ieat array, -1 if not used          */
   int*    outbeg; /* Array [0..knots-1] with starting slot index */
                   /* for the oeat array, -1 if not used          */
   /* Edges
    */
   int     esize;  /* Count of allocated edge slots               */
   int     edges;  /* Count of edges in the graph                 */
   double* cost;   /* Array [0..edges-1] of positiv edge costs    */
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

/* grphpath.c
 */
extern void   graph_path_init(const GRAPH*);
extern void   graph_path_exit(void);
extern void   graph_path_exec(const GRAPH*, int, int, const double*, PATH*);
extern void   graph_path_length(const GRAPH*, const PATH*);

/* grphmcut.c
 */
extern void   graph_mincut_init(const GRAPH*);
extern void   graph_mincut_exit(void);
extern void   graph_mincut_exec(const GRAPH*, int, int, const int*, int*, int);

/* grphload.c
 */
extern GRAPH* graph_load(const char*, PRESOL*);

/* grphsave.c
 */
extern void graph_save(const GRAPH*, const char*, FILETYPE);

/* grph2fig.c
 */
extern void graph_writefig(const GRAPH*, const char*, const double *, int);
extern void graph_bfscoord(GRAPH* g);
extern void graph_boxcoord(GRAPH* g);

/* reduce.c
 */
extern double reduce(GRAPH*, int);

/* sdtest.c
 */
extern int    sd_reduction(GRAPH*);
extern int    bd3_reduction(GRAPH*);
extern int    nsv_reduction(GRAPH*, double*);

/* validate.c
 */
extern int    validate(const GRAPH*, const double*);

#endif /* !_GRAPH_H_ */
