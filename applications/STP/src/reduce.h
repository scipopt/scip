/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce.h
 * @brief  includes various reduction methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_REDUCE_H_
#define APPLICATIONS_STP_SRC_REDUCE_H_


#include "scip/scip.h"
#include "graph.h"

#define STP_REDUCTION_NONE      0
#define STP_REDUCTION_BASIC     1
#define STP_REDUCTION_ADVANCED  2


/** lightweight minimum spanning tree structure that allows to add vertices to given MST on complete graph (in CSR format) */
typedef struct dynamic_complete_minimum_spanning_tree DCMST;

/** auxiliary data structure for ruling out all 1-hop stars of a given node */
typedef struct node_one_hop_star STAR;

/** SD distance graph data */
typedef struct special_distance_graph
{
   GRAPH*                distgraph;          /**< (complete) distance graph */
   PATH*                 sdmst;              /**< MST on sdgraph */
   SCIP_Real*            mstcosts;           /**< maximum MST edge costs in descending order */
   int*                  nodes_id;           /**< number of each node in original graph */
   int*                  halfedge_isInMst;   /**< signifies whether edge of original graph is part of MST
                                                  NOTE: operates on edges / 2! */
   SCIP_Real             mstmaxcost;         /**< maximum edge cost */
   int                   nnodesorg;          /**< number of nodes of original graph */
   int                   nedgesorg;          /**< number of edges of original graph */
} SDGRAPH;


/** reduced cost result data */
typedef struct reduce_costs_data
{
   SCIP_Real*            redEdgeCost;        /**< reduced costs */
   SCIP_Real*            rootToNodeDist;     /**< shortest path distances from root  */
   PATH*                 nodeTo3TermsPaths;  /**< paths to three nearest terminals */
   int*                  nodeTo3TermsBases;  /**< three nearest terminals */
   SCIP_Real             cutoff;             /**< reduced cost cutoff value or -1.0 if not used */
   SCIP_Real             dualBound;          /**< dual bound or -1.0 if not used */
   int                   redCostRoot;        /**< graph root for reduced cost calculation */
#ifndef NDEBUG
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
#endif
} REDCOST;



/** reduction parameters */
typedef struct reduction_parameters
{
   SCIP_Bool             dualascent;         /**< do dual-ascent reduction? */
   SCIP_Bool             boundreduce;        /**< do bound-based reduction? */
   SCIP_Bool             nodereplacing;      /**< should node replacement (by edges) be performed? */
   int                   reductbound;        /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             userec;             /**< use recombination heuristic? */
   SCIP_Bool             fullreduce;         /**< use full reductions? (including extended techniques) */
} RPARAMS;




/** reduced cost reduction parameters */
typedef struct reduce_costs_parameters
{
   int                   prevrounds;         /**< number of reduction rounds that have been performed already */
   SCIP_Bool             useRec;             /**< use recombination heuristic? */
   SCIP_Bool             useExtRed;          /**< use extended tests? */
   SCIP_Bool             nodereplacing;      /**< should node replacement (by edges) be performed? */
   /* PC/MW only values: */
   SCIP_Bool             pcmw_solbasedda;    /**< rerun Da based on best primal solution */
   SCIP_Bool             pcmw_useMultRoots;  /**< vary root for DA? (if possible) */
   SCIP_Bool             pcmw_markroots;     /**< should terminals proven to be part of an opt. sol. be marked as such? */
   SCIP_Bool             pcmw_fastDa;        /**< run dual ascent heuristic in fast mode? */
} RPDA;


/* reduce.c
 */
extern SCIP_RETCODE reduceLevel0(SCIP*, GRAPH*);
extern SCIP_RETCODE reduceLevel0infeas(SCIP*, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE reduceLevel0RpcRmw(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE reduceLevel0RpcRmwInfeas(SCIP*, GRAPH*, SCIP_Real*, SCIP_Bool*);
extern SCIP_RETCODE reduceStp(SCIP*, GRAPH*, SCIP_Real*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reducePc(SCIP*, const int*, GRAPH*, SCIP_Real*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduceMw(SCIP*, GRAPH*, SCIP_Real*, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduceHc(SCIP*, GRAPH*, SCIP_Real*, int);
extern SCIP_RETCODE reduceSap(SCIP*, GRAPH*, SCIP_Real*, int);
extern SCIP_RETCODE reduceNw(SCIP*, GRAPH*, SCIP_Real*, int);
extern SCIP_RETCODE redLoopStp(SCIP*, const RPARAMS*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*);
extern SCIP_RETCODE redLoopPc(SCIP*, const int*, GRAPH*, PATH*, PATH*,  GNODE**, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, SCIP_Bool, SCIP_Bool, SCIP_Bool, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE redLoopMw(SCIP*, GRAPH*, PATH*, GNODE**, SCIP_Real*, int*, int*, int*, int*, STP_Bool*, SCIP_Real*, STP_Bool, STP_Bool, STP_Bool, int, SCIP_Bool);
extern SCIP_RETCODE reduce(SCIP*, GRAPH*, SCIP_Real*, int, int, SCIP_Bool);


/* reduce_alt.c
 */
extern SCIP_RETCODE    reduce_ans(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_ansAdv(SCIP*, GRAPH*, int*, SCIP_Bool);
extern SCIP_RETCODE    reduce_ansAdv2(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_nnp(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_nv(SCIP*, GRAPH*, PATH*, double*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_nvAdv(SCIP*, const int*, GRAPH*, PATH*, SCIP_Real*, double*, int*, int*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_sl(SCIP*, const int*, GRAPH*, PATH*, double*, int*, int*, int*, int*, STP_Bool*, int*, int*);
extern SCIP_RETCODE    reduce_cnsAdv(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    reduce_npv(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_chain2(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int);


/* reduce_sd.c
 */
extern SCIP_RETCODE    reduce_ledge(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_sdsp(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, int*);
extern SCIP_RETCODE    reduce_sdStar(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdStarBiased(SCIP*, int, const int*,  GRAPH*, int*);
extern SCIP_RETCODE    reduce_sdStarPc(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdStarPc2(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalk(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdWalk_csr(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalkTriangle(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalkExt(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdWalkExt2(SCIP*, int, const int*, GRAPH*, int*,  SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdspSap(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_sd(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, SCIP_Bool, int*);
extern SCIP_RETCODE    reduce_sdPc(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_getSd(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int, int, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    reduce_getSdPcMw(SCIP*, const GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE    reduce_nts(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_bd34(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, SCIP_Real*);
extern SCIP_RETCODE    reduce_bd34WithSd(SCIP*, GRAPH*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);


/* reduce_bnd.c
 */
extern SCIP_RETCODE    reduce_bound(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundMw(SCIP*, GRAPH*, PATH*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundPruneHeur(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, const int*, const int*, int*, int);
extern SCIP_RETCODE    reduce_boundHop(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopR(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopRc(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);


/* reduce_da.c
 */
extern SCIP_RETCODE    reduce_da(SCIP*, GRAPH*, const RPDA*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*, SCIP_RANDNUMGEN*);
extern SCIP_RETCODE    reduce_daSlackPrune(SCIP*, SCIP_VAR**, GRAPH*, PATH*, GNODE**, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*, STP_Bool*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    reduce_daPcMw(SCIP*, GRAPH*, const RPDA*, PATH*, GNODE**, SCIP_Real*, int*, int*, int*, STP_Bool*, int*, SCIP_RANDNUMGEN*, SCIP_Real);


/* reduce_ext.c
 */
extern SCIP_RETCODE    reduce_deleteConflictEdges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_extendedCheck3Tree(SCIP*, const GRAPH*, int, const SCIP_Real*, const SCIP_Real*, const PATH*, const int*, SCIP_Real, const int*, int, SCIP_Real*, SCIP_Bool*, unsigned int*, int*, SCIP_Bool*);
extern int             reduce_extendedEdge(SCIP*, GRAPH*, const PATH*, const SCIP_Real*, const SCIP_Real*, const int*, SCIP_Real, int, STP_Bool*, SCIP_Bool);


/* reduce_simple.c
 */
extern SCIP_RETCODE    reduce_simple(SCIP*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern SCIP_RETCODE    reduce_simple_hc(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_sap(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_deleteMultiedges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_contract0Edges(SCIP*, GRAPH*, SCIP_Bool);
extern SCIP_RETCODE    reduce_aritculations(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_fixedConflicts(SCIP*, const int*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_rpt(SCIP*, GRAPH*, SCIP_Real*, int*);
extern void            reduce_identifyNonLeafTerms(SCIP*, GRAPH*);


/* reduce_pcsimple.c
 */
extern SCIP_RETCODE    reduce_simple_mw(SCIP*, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_pc(SCIP*, const int*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern void            reduce_removeDeg0NonLeafTerms(SCIP*, GRAPH*, SCIP_Real*);


/* reduce_util.c
 */
extern SCIP_RETCODE    reduce_applyPseudoDeletions(SCIP*, const REDCOST*, const SCIP_Bool*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_dcmstInit(SCIP*, int, DCMST**);
extern void            reduce_dcmstFree(SCIP*, DCMST**);
SCIP_Bool              reduce_dcmstMstIsValid(SCIP*, const CSR*);
extern void            reduce_dcmstAddNode(SCIP*, const CSR*, const SCIP_Real*, DCMST*, CSR*);
extern void            reduce_dcmstAddNodeInplace(SCIP*, const SCIP_Real*, DCMST*, CSR*);
extern void            reduce_dcmstGet0NodeMst(SCIP*, CSR*);
extern void            reduce_dcmstGet1NodeMst(SCIP*, CSR*);
extern void            reduce_dcmstGet2NodeMst(SCIP*, SCIP_Real, CSR*);
extern void            reduce_dcmstGet3NodeMst(SCIP*, SCIP_Real, SCIP_Real, SCIP_Real, CSR*);
extern SCIP_Real       reduce_dcmstGetExtWeight(SCIP*, const CSR*, const SCIP_Real*, DCMST*);
extern SCIP_Real       reduce_dcmstGetWeight(SCIP*, const CSR*);
extern int             reduce_dcmstGetMaxnnodes(const DCMST*);
extern SCIP_Real*      reduce_dcmstGetAdjcostBuffer(const DCMST*);
extern SCIP_RETCODE    reduce_starInit(SCIP*, int, STAR**);
extern void            reduce_starFree(SCIP*, STAR**);
extern void            reduce_starReset(const GRAPH*, int, STAR*);
extern const int*      reduce_starGetNext(STAR*, int*);
extern const int*      reduce_starGetRuledOutEdges(STAR*, int*);
extern int             reduce_starGetCenter(const STAR*);
extern void            reduce_starSetRuledOut(STAR*);
extern void            reduce_starSetFailed(STAR*);
extern SCIP_Bool       reduce_starAllAreChecked(const STAR*);
extern SCIP_RETCODE    reduce_redcostdataInit(SCIP*, int, int, SCIP_Real, int, REDCOST*);
extern void            reduce_redcostdataFreeMembers(SCIP*, REDCOST*);

#endif /* APPLICATIONS_STP_SRC_REDUCE_H_ */
