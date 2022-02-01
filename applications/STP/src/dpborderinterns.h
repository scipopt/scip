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

/**@file   dpborderinterns.h
 * @brief  Dynamic programming internals for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 *
 * Internal methods and data structures for DP.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_DPBORDERINTERNS_H_
#define APPLICATIONS_STP_SRC_DPBORDERINTERNS_H_


#include "scip/scip.h"
#include "graph.h"
#include "stpvector.h"
#include "dpborder_hashmap.h"

#define BPBORDER_MAXNPARTITIONS 50000000
#define BPBORDER_MAXBORDERSIZE  16
#define DPB_Ptype char
#define DPBORDER_GROWTH_FACTOR 4


/** nodes sequence structure */
typedef struct dynamic_programming_border_nodes_sequence
{
   int*                  nodessquence;       /**< ordering of the nodes */
   uint64_t              maxnpartitions;     /**< maximum number of partitions */
   int                   maxbordersize;      /**< maximum size of any border */
   int                   nnodes;             /**< number of nodes of underlying graph */
} DPBSEQUENCE;


/** nodes sequence structure */
typedef struct dynamic_programming_border_level
{
   STP_Vectype(int)      bordernodesMapToOrg;/**< maps border nodes to original nodes */
   int                   globalstartidx;     /**< start position of level in global data */
   int                   nbordernodes;       /**< size of border */
   int                   extnode;            /**< extension nodes */
   SCIP_Bool             exnodeIsTerm;       /**< is the extension node a terminal? */
} DPBLEVEL;


/** single partition */
typedef struct dynamic_programming_border_partition
{
   DPB_Ptype*            partchars;          /**< partition characters */
   int                   partsize;           /**< size of partition */
   DPB_Ptype             delimiter;          /**< delimiter */
} DPBPART;



/** DP border structure */
struct dynamic_programming_border
{
   DPBHASHMAP            hashmap;            /**< hash map */
   DPBSEQUENCE*          dpbsequence;        /**< ordering of nodes */
   STP_Vectype(DPBLEVEL*) borderlevels;      /**< data for each border */
   SCIP_Bool*            nodes_isBorder;     /**< marks whether node is in current border */
   int*                  nodes_outdeg;       /**< degree w.r.t. not yet visited nodes */
   int*                  bordercharmap;      /**< maps last border chars to current border chars */
   SCIP_Real*            borderchardists;    /**< distance for last border nodes (chars) to extension nodes */
   STP_Vectype(int)      bordernodes;        /**< current border nodes */
   STP_Vectype(int)      prevbordernodes;    /**< nodes that are in previous but not current border */
   DPB_Ptype*            global_partitions;  /**< partitions */
   STP_Vectype(int)      global_partstarts;  /**< CSR like starts of partitions in array "global_partitions" */
   STP_Vectype(SCIP_Real) global_partcosts;  /**< costs of each partition */
   STP_Vectype(int)      global_predparts;   /**< predecessor partitions; of size global_npartitions */
   STP_Vectype(SCIP_Bool) global_partsUseExt; /**< partition uses extension node? */
   SCIP_Real             global_obj;         /**< objective */
   int                   global_npartitions; /**< number of global partitions */
   int                   global_partcap;     /**< capacity of array global_partitions */
   int                   global_optposition; /**< index of best solution partition */
   int                   ntermsvisited;      /**< number of already visited nodes */
   int                   nterms;             /**< number of terminals */
   int                   nnodes;             /**< number of nodes of underlying graph */
   DPB_Ptype             extborderchar;      /**< -1 if extnode is not contained! */
};

/*
 * Inline methods
 */


/** gets border delimiter for given iteration */
static inline
int dpborder_getDelimiter(
   const DPBORDER*       dpborder,           /**< border */
   int                   iteration           /**< iteration number */
)
{
   assert(iteration >= 0);
   assert(dpborder->borderlevels[iteration]->nbordernodes > 0);
   assert(dpborder->borderlevels[iteration]->nbordernodes <= BPBORDER_MAXBORDERSIZE);

   return dpborder->borderlevels[iteration]->nbordernodes;
}


/** gets border delimiter */
static inline
int dpborder_getTopDelimiter(
   const DPBORDER*       dpborder            /**< border */
)
{
   const int pos = StpVecGetSize(dpborder->borderlevels) - 1;

   assert(pos >= 0);
   assert(dpborder->borderlevels[pos]->nbordernodes >= 0);
   assert(dpborder->borderlevels[pos]->nbordernodes <= BPBORDER_MAXBORDERSIZE);

   return dpborder->borderlevels[pos]->nbordernodes;
}


/** gets top level */
static inline
DPBLEVEL* dpborder_getTopLevel(
   const DPBORDER*       dpborder            /**< border */
)
{
   const int pos = StpVecGetSize(dpborder->borderlevels) - 1;
   assert(pos >= 0);
   assert(dpborder->borderlevels[pos]);

   return dpborder->borderlevels[pos];
}


/** gets previous level */
static inline
DPBLEVEL* dpborder_getPredLevel(
   const DPBORDER*       dpborder            /**< border */
)
{
   const int pos = StpVecGetSize(dpborder->borderlevels) - 1;
   assert(pos >= 1);
   assert(dpborder->borderlevels[pos - 1]);

   return dpborder->borderlevels[pos - 1];
}


/*
 * Internal interface methods
 */

extern SCIP_Real dpborder_partGetConnectionCost(const DPBORDER*, const DPBPART*, const int*, int);
extern int dpborder_partglobalGetCard(int, int, const DPBORDER*);
extern int dpborder_partGetIdxNew(SCIP*, const DPBPART*, const int*, int, DPBORDER*);
extern int dpborder_partGetIdxNewExclusive(SCIP*, const DPBPART*, DPBORDER*);
extern STP_Vectype(int)  dpborder_partGetCandstarts(SCIP*, const DPBPART*, const DPBORDER*);
extern SCIP_Bool  dpborder_partIsValid(const DPBPART*);
extern void  dpborder_partPrint(const DPBPART*);
extern void  dpborder_markSolNodes(const DPBORDER*, STP_Bool* RESTRICT);
extern SCIP_RETCODE  dpborder_dpbsequenceInit(SCIP*, const GRAPH*, DPBSEQUENCE**);
extern void          dpborder_dpbsequenceFree(SCIP*, DPBSEQUENCE**);
extern void          dpborder_dpbsequenceCopy(const DPBSEQUENCE*, DPBSEQUENCE*);
extern SCIP_RETCODE  dpborder_dpblevelInit(SCIP*, DPBLEVEL**);
extern void          dpborder_dpblevelFree(SCIP*, DPBLEVEL**);
extern SCIP_RETCODE  dpborder_coreComputeOrderingSimple(SCIP*, const GRAPH*, DPBORDER*);
extern SCIP_RETCODE  dpborder_coreUpdateOrdering(SCIP*, const GRAPH*, DPBORDER*);
extern SCIP_RETCODE  dpborder_coreSolve(SCIP*, const GRAPH*, DPBORDER*, SCIP_Bool*);

#endif /* APPLICATIONS_STP_SRC_DPBORDERINTERNS_H_ */
