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
#include "mincut.h"
#include "bidecomposition.h"
#include "stpvector.h"
#include "redcosts.h"
#include "reducedefs.h"
#include "extreducedefs.h"
#include "termsepadefs.h"


#ifdef __cplusplus
extern "C" {
#endif


/* reduce_base.c
 */
extern int reduce_getMinNreductions(const GRAPH*, int);
extern SCIP_RETCODE reduce_baseInit(SCIP*, const GRAPH*, REDBASE**);
extern void reduce_baseFree(SCIP*, REDBASE**);
extern SCIP_RETCODE reduce_stp(SCIP*, GRAPH*, REDSOL*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduce_pc(SCIP*, REDSOL*, GRAPH*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduce_mw(SCIP*, REDSOL*, GRAPH*, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduce_hc(SCIP*, GRAPH*, SCIP_Real*, int);
extern SCIP_RETCODE reduce_sap(SCIP*, GRAPH*, SCIP_Bool, SCIP_Real*, int);
extern SCIP_RETCODE reduce_nw(SCIP*, GRAPH*, SCIP_Real*, int);
extern SCIP_RETCODE reduce_dc(SCIP*, GRAPH*, SCIP_Real*, int);
extern SCIP_RETCODE reduce_redLoopStp(SCIP*, GRAPH*, REDBASE*);
extern SCIP_RETCODE reduce_redLoopPc(SCIP*, REDSOL*, GRAPH*, PATH*, PATH*, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*,
                              SCIP_Bool, SCIP_Bool, SCIP_Bool, int, SCIP_Bool, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduce_redLoopMw(SCIP*, REDSOL*, GRAPH*, PATH*, SCIP_Real*, int*, int*, int*, STP_Bool*,
                              STP_Bool, STP_Bool, STP_Bool, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE reduce_exec(SCIP*, GRAPH*, REDSOL*, int, int, SCIP_Bool);

/* reduce_sol.c
 */
extern SCIP_RETCODE reduce_sollocalInit(SCIP*, const GRAPH*, REDSOLLOCAL**);
extern void         reduce_sollocalFree(SCIP*, REDSOLLOCAL**);
extern int*         reduce_sollocalGetSolnode(REDSOLLOCAL*);
extern void reduce_sollocalSetOffset(SCIP_Real, REDSOLLOCAL*);
extern void reduce_sollocalUpdateUpperBound(SCIP_Real, REDSOLLOCAL*);
extern SCIP_RETCODE reduce_sollocalUpdateNodesol(SCIP*, const int*, GRAPH*, REDSOLLOCAL*);
extern SCIP_Real reduce_sollocalGetUpperBound(const REDSOLLOCAL*);
extern SCIP_Real reduce_sollocalGetUpperBoundWithOffset(const REDSOLLOCAL*);
extern SCIP_Bool reduce_sollocalHasUpperBound(const REDSOLLOCAL*);
extern SCIP_Bool reduce_sollocalUsesNodesol(const REDSOLLOCAL*);
extern SCIP_RETCODE reduce_sollocalRebuildTry(SCIP*, GRAPH*, REDSOLLOCAL*);
extern SCIP_RETCODE reduce_solInit(SCIP*, const GRAPH*, SCIP_Bool, REDSOL**);
extern void         reduce_solFree(SCIP*, REDSOL**);
extern SCIP_RETCODE reduce_solInitLocal(SCIP*, const GRAPH*, REDSOL*, REDSOLLOCAL**);
extern void         reduce_solFinalizeLocal(SCIP*, const GRAPH*, REDSOL*);
extern void         reduce_solReInitLocal(const GRAPH*, REDSOL*);
extern void         reduce_solPack(const GRAPH*, const int*, int, REDSOL*);
extern SCIP_RETCODE reduce_solLevelAdd(SCIP*, const GRAPH*, REDSOL*);
extern SCIP_RETCODE reduce_solLevelTopUpdate(SCIP*, const GRAPH*, REDSOL*);
extern void reduce_solLevelTopFinalize(SCIP*, GRAPH*, REDSOL*);
extern void reduce_solLevelTopRemove(SCIP*, REDSOL*);
extern void reduce_solLevelTopClean(SCIP*, REDSOL*);
extern void reduce_solLevelTopTransferSolBack(const int*, REDSOL*);
extern SCIP_RETCODE reduce_solLevelTopTransferSolTo(const int*, REDSOL*);
extern void reduce_solSetOffset(SCIP_Real, REDSOL*);
extern SCIP_Real reduce_solGetOffset(const REDSOL*);
extern SCIP_Real reduce_solGetUpperBoundWithOffset(const REDSOL*);
extern const int* reduce_solGetNodesolPointer(const REDSOL*);
extern SCIP_Bool reduce_solUsesNodesol(const REDSOL*);
extern SCIP_Real* reduce_solGetOffsetPointer(REDSOL*);
extern SCIP_RETCODE         reduce_solGetEdgesol(SCIP*, GRAPH*, REDSOL*, SCIP_Real*, int*);
extern SCIP_RETCODE         reduce_solAddNodesol(const GRAPH*, const int*, REDSOL*);
extern void         reduce_solGetNodesol(const GRAPH*, REDSOL*, int*);


/* reduce_alt.c
 */
extern SCIP_RETCODE    reduce_impliedProfitBased(SCIP*, int, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_impliedProfitBasedRpc(SCIP*, GRAPH*, REDSOLLOCAL*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_ans(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_ansAdv(SCIP*, GRAPH*, int*, SCIP_Bool);
extern SCIP_RETCODE    reduce_ansAdv2(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_nnp(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_nv(SCIP*, GRAPH*, PATH*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_nvAdv(SCIP*, const int*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_sl(SCIP*, const int*, GRAPH*, PATH*, SCIP_Real*, int*, int*, STP_Bool*, int*, int*);
extern SCIP_RETCODE    reduce_nsvImplied(SCIP*, const SD*, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_nsvImpliedRecord(SCIP*, const SD*, GRAPH*, STP_Vectype(int)*);
extern SCIP_RETCODE    reduce_cnsAdv(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    reduce_npv(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_chain2(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int);


/*
 * reduce_sd.c
 */
extern SCIP_RETCODE    reduce_sdEdgeCliqueStar(SCIP*, int, GRAPH*, int*);
extern SCIP_RETCODE    reduce_sdImpLongEdge(SCIP*, const int*, GRAPH*, SD*, int*);
extern SCIP_RETCODE    reduce_sdsp(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*, int*, int, SCIP_Bool);
extern SCIP_RETCODE    reduce_sdStar(SCIP*, int, SCIP_Bool, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdStarBiased(SCIP*, int, SCIP_Bool, GRAPH*, int*);
extern SCIP_RETCODE    reduce_sdStarBiasedWithProfit(SCIP*, int, const SDPROFIT*, SCIP_Bool, GRAPH*, int*);
extern SCIP_RETCODE    reduce_sdStarPc(SCIP*, int, const int*, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdStarPc2(SCIP*, int, SCIP_Bool, GRAPH*, SCIP_Real*, int*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalk(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdWalk_csr(SCIP*, int, const int*, GRAPH*, int*, SCIP_Real*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalkTriangle(SCIP*, int, SCIP_Bool, GRAPH*, int*, SCIP_Real*, int*, STP_Bool*, DHEAP*, int*);
extern SCIP_RETCODE    reduce_sdWalkExt(SCIP*, int, SCIP_Bool, GRAPH*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdWalkExt2(SCIP*, int, const int*, GRAPH*, int*,  SCIP_Real*, int*, int*, int*, STP_Bool*, int*);
extern SCIP_RETCODE    reduce_sdspSap(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int);
extern SCIP_RETCODE    reduce_sd(SCIP*, GRAPH*, REDBASE*, int*);
extern SCIP_RETCODE    reduce_sdBiased(SCIP*, SD*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_sdBiasedNeighbor(SCIP*, SD*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_sdPc(SCIP*, GRAPH*, PATH*, int*, int*, int*, int*, int*, int*);
extern SCIP_Real       reduce_sdGetSd(const GRAPH*, int, int, SCIP_Real, SCIP_Real, SD*);
extern SCIP_Real       reduce_sdGetSdIntermedTerms(const GRAPH*, int, int, SCIP_Real, SCIP_Real, SD*);
extern SCIP_RETCODE    reduce_getSdByPaths(SCIP*, GRAPH*, PATH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int, int, int, SCIP_Bool, SCIP_Bool);
extern SCIP_RETCODE    reduce_bdk(SCIP*, int, GRAPH*, int*);
extern SCIP_RETCODE    reduce_bdkBiased(SCIP*, int, GRAPH*, int*);
extern SCIP_RETCODE    reduce_bdkWithSd(SCIP*, int, SD*, GRAPH*, int*);


/*
 * reduce_path.c
 */
extern SCIP_RETCODE    reduce_pathreplace(SCIP*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_pathreplaceExt(SCIP*, GRAPH*, EXTPERMA*, int*);


/* reduce_sdcomp.c
 */
extern SCIP_RETCODE    reduce_bd34(SCIP*, GRAPH*, PATH*, PATH*, int*, int*, int*, int*, int*, int*, int, SCIP_Real*);
extern SCIP_RETCODE    reduce_bd34WithSd(SCIP*, GRAPH*, SDGRAPH*, PATH*, int*, int*);


/* reduce_bnd.c
 */
extern SCIP_RETCODE    reduce_bound(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundMw(SCIP*, GRAPH*, PATH*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundPruneHeur(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, const int*, const int*, int*, int);
extern SCIP_RETCODE    reduce_boundHop(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopR(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*, int*, int*);
extern SCIP_RETCODE    reduce_boundHopRc(SCIP*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, SCIP_Bool);
extern SCIP_RETCODE    reduce_boundHopDa(SCIP*, GRAPH*, int*, SCIP_RANDNUMGEN*);

/* reduce_da.c
 */
extern SCIP_RETCODE    reduce_da(SCIP*, GRAPH*, const RPDA*, REDSOLLOCAL*, SCIP_Real*, int*, SCIP_RANDNUMGEN*);
extern SCIP_RETCODE    reduce_dapaths(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_daSlackPrune(SCIP*, GRAPH*, int, SCIP_Bool, int*, int*, int*, SCIP_Real*, SCIP_Bool*, SCIP_Bool*);
extern SCIP_RETCODE    reduce_daPcMw(SCIP*, GRAPH*, const RPDA*, REDSOLLOCAL*, PATH*, SCIP_Real*, int*, int*, int*, STP_Bool*, int*, SCIP_RANDNUMGEN*, SCIP_Real);

/* reduce_ext.c
 */
extern SCIP_RETCODE    reduce_deleteConflictEdges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_extendedCheck3Tree(SCIP*, const GRAPH*, int, const SCIP_Real*, const SCIP_Real*, const PATH*, const int*, SCIP_Real, const int*, int, SCIP_Real*, SCIP_Bool*, unsigned int*, int*, SCIP_Bool*);
extern int             reduce_extendedEdge(SCIP*, GRAPH*, const PATH*, const SCIP_Real*, const SCIP_Real*, const int*, SCIP_Real, int, STP_Bool*, SCIP_Bool);


/* reduce_simple.c
 */
extern void            reduce_nodesDeg1(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_simple(SCIP*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern void    reduce_simple_dc(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_simple_hc(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_sap(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_deleteMultiedges(SCIP*, GRAPH*);
extern SCIP_RETCODE    reduce_contract0Edges(SCIP*, GRAPH*, int*, SCIP_Bool);
extern SCIP_RETCODE    reduce_fixedConflicts(SCIP*, const int*, GRAPH*, int*);
extern SCIP_RETCODE    reduce_cutEdgeTryPrune(SCIP*, int, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE    reduce_rpt(SCIP*, GRAPH*, SCIP_Real*, int*);
extern void            reduce_identifyNonLeafTerms(SCIP*, GRAPH*);
extern SCIP_RETCODE reduce_unconnected(SCIP*, GRAPH*);
extern SCIP_RETCODE reduce_unconnectedForDirected(SCIP*, GRAPH*);
extern SCIP_RETCODE reduce_unconnectedInfeas(SCIP*, SCIP_Bool, GRAPH*, SCIP_Bool*);


/* reduce_sepa.c
 */
extern SCIP_RETCODE    reduce_articulations(SCIP*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_bidecomposition(SCIP*, GRAPH*, REDBASE*, int*, SCIP_Bool*);
extern SCIP_RETCODE    reduce_bidecompositionExact(SCIP*, GRAPH*, REDBASE*, int*, int*);
extern SCIP_RETCODE    reduce_nonTerminalComponents(SCIP*, const CUTNODES*, GRAPH*, SCIP_Real*, int*);



/* reduce_pcsimple.c
 */
extern SCIP_RETCODE    reduce_simple_mw(SCIP*, GRAPH*, int*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_simple_pc(SCIP*, const int*, GRAPH*, SCIP_Real*, int*, int*, int*);
extern void            reduce_removeDeg0NonLeafTerms(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE reduce_unconnectedRpcRmw(SCIP*, GRAPH*, SCIP_Real*);
extern SCIP_RETCODE reduce_unconnectedRpcRmwInfeas(SCIP*, GRAPH*, SCIP_Real*, SCIP_Bool*);


/* reduce_util.c
 */
extern void            reduce_impliedNodesGet(SCIP*, const GRAPH*, STP_Vectype(int)*);
extern void            reduce_impliedNodesRepair(SCIP*, const GRAPH*, int, int, STP_Vectype(int)*);
extern SCIP_Bool       reduce_impliedNodesIsValid(const GRAPH*, const STP_Vectype(int)*);
extern SCIP_RETCODE    reduce_applyPseudoDeletions(SCIP*, const SCIP_Bool*, REDCOST*, GRAPH*, SCIP_Real*, int*);
extern SCIP_RETCODE    reduce_blctreeInit(SCIP*, GRAPH*, BLCTREE**);
extern void            reduce_blctreeFree(SCIP*, BLCTREE**);
extern int             reduce_blctreeGetMstNedges(const BLCTREE*);
extern void            reduce_blctreeGetMstEdges(const GRAPH*, const BLCTREE*, int*);
extern void            reduce_blctreeGetMstEdgesToCutDist(const GRAPH*, const BLCTREE*, SCIP_Real* RESTRICT, SCIP_Real* RESTRICT);
extern void            reduce_blctreeGetMstBottlenecks(const GRAPH*, const BLCTREE*, SCIP_Real*);
extern void            reduce_blctreeGetMstEdgesState(const GRAPH*, const BLCTREE*, SCIP_Bool*);
extern SCIP_RETCODE    reduce_blctreeRebuild(SCIP*, GRAPH*, BLCTREE*);
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
extern void            reduce_starResetWithEdges(const GRAPH*, const STP_Vectype(int), STAR*);
extern const int*      reduce_starGetNext(STAR*, int*);
extern const int*      reduce_starGetNextAndPosition(STAR*, int*, int*);
extern const int*      reduce_starGetRuledOutEdges(STAR*, int*);
extern int             reduce_starGetCenter(const STAR*);
extern void            reduce_starCurrentSetRuledOut(STAR*);
extern void            reduce_starCurrentSetFailed(STAR*);
extern SCIP_Bool       reduce_starHasPromisingEdges(const STAR*);
extern SCIP_Bool       reduce_starAllAreChecked(const STAR*);


/* reduce_sdutil.c
 */
extern SCIP_RETCODE    reduce_sdneighborInit(SCIP*, const GRAPH*, SDN**);

extern void            reduce_sdneighborGetCloseTerms(const GRAPH*, const SDN*, int, SCIP_Real, int* RESTRICT, SCIP_Real* RESTRICT, int* RESTRICT);
extern void            reduce_sdneighborFree(SCIP*, SDN**);
extern const SCIP_Bool*    reduce_sdneighborGetBlocked(const SDN*);
extern SCIP_RETCODE    reduce_sdUpdateWithSdNeighbors(SCIP*, GRAPH*, SD*, int*);
extern SCIP_RETCODE    reduce_sdprofitInit(SCIP*, const GRAPH*, SDPROFIT**);
extern SCIP_RETCODE    reduce_sdprofitInit1stOnly(SCIP*, const GRAPH*, const SCIP_Real*, SDPROFIT**);
extern void            reduce_sdprofitFree(SCIP*, SDPROFIT**);
extern SCIP_RETCODE    reduce_sdprofitUpdateFromBLC(SCIP*, const GRAPH*, const BLCTREE*, SCIP_Bool, SDPROFIT*);
extern SCIP_RETCODE    reduce_sdprofitBuildFromBLC(SCIP*, const GRAPH*, const BLCTREE*, SCIP_Bool, SDPROFIT*);
extern void            reduce_sdprofitPrintStats(const GRAPH*, const SDPROFIT*);
extern SCIP_RETCODE     reduce_sdInit(SCIP*, GRAPH*, SD**);
extern SCIP_RETCODE     reduce_sdInitBiased(SCIP*, GRAPH*, SD**);
extern SCIP_RETCODE     reduce_sdInitBiasedBottleneck(SCIP*, GRAPH*, SD**);
extern SCIP_RETCODE     reduce_sdRepair(SCIP*, int, SCIP_Bool, GRAPH*, SD*);
extern SCIP_RETCODE     reduce_sdRepairSetUp(SCIP*, const GRAPH*, SD*);
extern SCIP_RETCODE     reduce_sdAddNeighborSd(SCIP*, const GRAPH*, SD*);
extern void             reduce_sdFree(SCIP*, SD**);
extern SCIP_RETCODE     reduce_sdGetSdsCliquegraph(SCIP*, const GRAPH*, int, const int*, DIJK*, SD*, GRAPH*);


/* reduce_sdgraph.c
 */
extern SCIP_RETCODE     reduce_sdgraphInit(SCIP*, const GRAPH*, SDGRAPH**);
extern SCIP_RETCODE     reduce_sdgraphInitFromDistGraph(SCIP*, const GRAPH*, GRAPH*, int*, SDGRAPH**);
extern SCIP_RETCODE     reduce_sdgraphInitBiased(SCIP*, const GRAPH*, const SDPROFIT*, SDGRAPH**);
extern SCIP_RETCODE     reduce_sdgraphInitBiasedFromTpaths(SCIP*, GRAPH*, const SDPROFIT*, const TPATHS*, SDGRAPH**);
extern SCIP_Real        reduce_sdgraphGetMaxCost(const SDGRAPH*);
extern const SCIP_Real* reduce_sdgraphGetOrderedMstCosts(const SDGRAPH*);
extern const int*       reduce_sdgraphGetOrgnodesToSdMap(const SDGRAPH*);
extern void             reduce_sdgraphInitOrderedMstCosts(SDGRAPH*);
extern const STP_Bool*  reduce_sdgraphGetMstHalfMark(const SDGRAPH*);
extern SCIP_Bool        reduce_sdgraphHasMstHalfMark(const SDGRAPH*);
extern SCIP_Bool        reduce_sdgraphEdgeIsInMst(const SDGRAPH*, int);
extern SCIP_Bool        reduce_sdgraphHasOrderedMstCosts(const SDGRAPH*);
extern SCIP_Real        reduce_sdgraphGetSd(int, int, SDGRAPH*);
extern void             reduce_sdgraphFree(SCIP*, SDGRAPH**);
extern void             reduce_sdgraphFreeFromDistGraph(SCIP*, SDGRAPH**);
extern void             reduce_sdgraphInsertEdge(SCIP*, int, int, SCIP_Real, int, int* RESTRICT, SDGRAPH*, SCIP_Bool*);
extern SCIP_RETCODE     reduce_sdgraphMstBuild(SCIP*, const GRAPH*, SDGRAPH*);
extern void             reduce_sdgraphMstSortCosts(SDGRAPH*);


/* reduce_termsepa.c
 */
extern SCIP_RETCODE     reduce_termcompInit(SCIP*, const GRAPH*, COMPBUILDER*, TERMCOMP**);
void                    reduce_termcompFree(SCIP*, TERMCOMP**);
extern SCIP_RETCODE     reduce_termcompInitTbottleneck(SCIP*, const int*, TERMCOMP*);
extern SCIP_RETCODE     reduce_termcompBuildSubgraph(SCIP*, GRAPH*, TERMCOMP*);
extern SCIP_RETCODE     reduce_termcompBuildSubgraphWithSds(SCIP*, GRAPH*, EXTPERMA*, TERMCOMP*);
extern SCIP_RETCODE     reduce_termcompChangeSubgraphToBottleneck(SCIP*, GRAPH*, TERMCOMP*, SCIP_Bool*);
extern void             reduce_termsepaGetNextComp(SCIP*, const GRAPH*, TERMSEPAS*, COMPBUILDER*, SCIP_Bool*);
extern SCIP_RETCODE     reduce_compbuilderInit(SCIP*, const GRAPH*, COMPBUILDER**);
void                    reduce_compbuilderFree(SCIP*, COMPBUILDER**);
SCIP_Real               reduce_compbuilderGetSubNodesRatio(const COMPBUILDER*);
void                    reduce_compbuilderPrintSeparators(const GRAPH*, const COMPBUILDER*);
void                    reduce_termcompChangeSubgraphToOrgCosts(const GRAPH*, TERMCOMP*);


/* reduce_termsepada.c
 */
extern SCIP_RETCODE     reduce_termsepaDaWithExperma(SCIP*, GRAPH*, EXTPERMA*, SCIP_Bool*, int*);
extern SCIP_RETCODE     reduce_termsepaDa(SCIP*, GRAPH*, int*);


/* reduce_termsepafull.c
 */
extern SCIP_RETCODE     reduce_termsepaFull(SCIP*, GRAPH*, int*, REDBASE*, int*);

#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_REDUCE_H_ */
