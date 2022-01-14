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

/**@file   dptermsinterns.h
 * @brief  Dynamic programming internals for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 *
 * Internal methods and data structures for DP.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_DPTERMSINTERNS_H_
#define APPLICATIONS_STP_SRC_DPTERMSINTERNS_H_

#include "scip/scip.h"
#include "scip/rbtree.h"
#include "graph.h"
#include "stpvector.h"
#include "stpbitset.h"
#include "stpprioqueue.h"

//#define STP_DPTERM_USEDA

/** dynamic programming search tree */
typedef struct dynamic_programming_search_tree DPSTREE;



/*
 * Data structures
 */


/** trace for reconstructing a sub-solution */
typedef struct solution_trace
{
   int                   prevs[2];           /**< marker to get ancestor solutions (0,1,2 ancestors possible) */
   SCIP_Real             cost;               /**< solution cost */
#ifdef STP_DPTERM_USEDA
   SCIP_Real             redcost;            /**< reduced solution cost */
#endif
   int                   root;               /**< solution root */
} SOLTRACE;


/** sub-solution with extension */
typedef struct dynamic_programming_subsolution
{
   SCIP_RBTREE_HOOKS;                        /**< for red-black tree */
   STP_Bitset            bitkey;             /**< key marking the terminals in sub-solution */
   STP_Vectype(SOLTRACE) traces;             /**< traces of solution */
} DPSUBSOL;


/** compressed graph with less information */
typedef struct dynamic_programming_graph
{
   int*                  terminals;          /**< array of terminals; in {0,1,...,nnodes - 1} */
   int*                  nodes_termId;       /**< per node: terminal (0,1,..), or -1 if non-terminal */
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
   int                   nterms;             /**< number of terminals */
} DPGRAPH;


/** additional data */
typedef struct dynamic_programming_misc
{
   STP_Vectype(SOLTRACE)  global_traces;
   STP_Vectype(STP_Bitset) global_termbits;
   STP_Vectype(int)      global_termbitscount;
   STP_Vectype(int)      global_starts;
   int                   opt_prev[2];
   SCIP_Real             opt_obj;
   int                   opt_root;
   int                   global_size;
} DPMISC;



/** reduced cost data */
typedef struct dynamic_programming_reduced_costs
{
   SCIP_Real*            csr_redcosts;
   SCIP_Real*            nodes_rootdist;
   SCIP_Real             cutoffbound;
   SCIP_Real             upperbound;
} DPREDCOST;



/** solver */
typedef struct dynamic_programming_solver
{
   STP_Vectype(int)      solnodes;           /**< (final) solution nodes */
   DPGRAPH*              dpgraph;            /**< graph */
   DPSUBSOL*             soltree_root;       /**< root of solution tree */
   DPSTREE*              dpstree;            /**< tree for finding solution combinations */
   DPMISC*               dpmisc;             /**< this and that */
   STP_PQ*               solpqueue;          /**< sub-solutions */
   DHEAP*                dheap;              /**< heap of size nnodes */
#ifdef STP_DPTERM_USEDA
   DPREDCOST*            dpredcosts;
#endif
} DPSOLVER;


/*
 * Macro hacks
 */


/* NOTE: needed to find element in a red-black tree */
#define SUBSOL_LT(key,subsol)  stpbitset_GT(key, subsol->bitkey)
#define SUBSOL_GT(key,subsol)  stpbitset_LT(key, subsol->bitkey)

static inline
SCIP_DEF_RBTREE_FIND(findSubsol, STP_Bitset, DPSUBSOL, SUBSOL_LT, SUBSOL_GT) /*lint !e123*/



/*
 * Inline methods
 */


/** initializes */
static inline
SCIP_RETCODE dpterms_dpsubsolInit(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSUBSOL**            subsol              /**< solution */
)
{
   DPSUBSOL* sub;
   SCIP_CALL( SCIPallocBlockMemory(scip, subsol) );
   sub = *subsol;
   sub->bitkey = NULL;
   sub->traces = NULL;

   return SCIP_OKAY;
}


/** frees */
static inline
void dpterms_dpsubsolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DPSUBSOL**            subsol              /**< solution */
)
{
   DPSUBSOL* sub = *subsol;

   if( sub->bitkey )
      stpbitset_free(scip, &(sub->bitkey));

   if( sub->traces )
   {
      StpVecFree(scip, sub->traces);
   }

   SCIPfreeBlockMemory(scip, subsol);
}



/*
 *
 */


/* dpterms_util.c
 */
extern SCIP_Bool        dpterms_intersectsEqualNaive(SCIP*, STP_Bitset, STP_Bitset, STP_Vectype(int), DPMISC*);
extern STP_Vectype(int) dpterms_collectIntersectsNaive(SCIP*, STP_Bitset, STP_Bitset, DPMISC*);
extern SCIP_RETCODE     dpterms_streeInit(SCIP*, int, int, DPSTREE**);
extern void             dpterms_streeFree(SCIP*, DPSTREE**);
extern SCIP_RETCODE     dpterms_streeInsert(SCIP*, STP_Bitset, STP_Bitset, int64_t, DPSTREE*);
extern STP_Vectype(int) dpterms_streeCollectIntersects(SCIP*, STP_Bitset, STP_Bitset, DPSTREE*);


/* dpterms_core.c
 */
extern SCIP_RETCODE     dpterms_coreSolve(SCIP*, GRAPH*, DPSOLVER*, SCIP_Bool*);



#endif /* APPLICATIONS_STP_SRC_DPTERMSINTERNS_H_ */
