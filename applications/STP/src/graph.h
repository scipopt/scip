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

/**@file   graph.h
 * @brief  includes various files containing graph methods used for Steiner tree problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_GRAPH_H_
#define APPLICATIONS_STP_GRAPH_H_

#include "scip/scip.h"
#include "misc_stp.h"
#include "graphdefs.h"
#include "redcosts.h"
#include "reducedefs.h"
#include "portab.h"
#include "stpvector.h"


#ifdef __cplusplus
extern "C" {
#endif


//#define STP_RPC_FIXEDPROPER // only for testing



/* graph_history.c
 */
extern SCIP_RETCODE   graph_initContractTracing(SCIP*, GRAPH*);
extern int            graph_contractTrace(int, const GRAPH*);
extern SCIP_Bool      graph_knot_hasContractTrace(int, const GRAPH*);
extern void          graph_fixed_resetMoved(GRAPH*);
extern SCIP_RETCODE  graph_singletonAncestors_init(SCIP*, const GRAPH*, int, SINGLETONANS*);
extern void          graph_singletonAncestors_freeMembers(SCIP*, SINGLETONANS*);
extern SCIP_Bool     graph_valid_ancestors(SCIP*, const GRAPH*);
extern SCIP_RETCODE   graph_initAncestors(SCIP*, GRAPH*);
/* Pseudo ancestors */
extern SCIP_RETCODE   graph_initPseudoAncestors(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_initPseudoAncestorsSized(SCIP*, int, GRAPH*);
extern void           graph_freePseudoAncestors(SCIP*, GRAPH*);
extern void           graph_edge_delPseudoAncestors(SCIP*, int, GRAPH*);
extern void           graph_knot_delPseudoAncestors(SCIP*, int, GRAPH*);
extern void           graph_edge_printPseudoAncestors(const GRAPH*, int);
extern void           graph_knot_printPseudoAncestors(const GRAPH*, int);
extern int            graph_edge_nPseudoAncestors(const GRAPH*, int);
extern int            graph_knot_nPseudoAncestors(const GRAPH*, int);
extern const int*     graph_edge_getPseudoAncestors(const GRAPH*, int);
extern IDX*           graph_edge_getAncestors(const GRAPH*, int);
extern const int*     graph_knot_getPseudoAncestors(const GRAPH*, int);
extern int            graph_getNpseudoAncestors(const GRAPH*);
extern void      graph_addPseudoAncestor(GRAPH*, int*);
extern void      graph_addPseudoAncestors(int, GRAPH*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendCopyEdge(SCIP*, int, int, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendCopyNode(SCIP*, int, int, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendCopyNodeToEdge(SCIP*, int, int, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendCopyEdgeToNode(SCIP*, int, int, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendCopySingToEdge(SCIP*, int, const SINGLETONANS*, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendCopyArrayToEdge(SCIP*, int, const int*, int, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendMoveEdge(SCIP*, int, int, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_appendMoveNode(SCIP*, int, int, SCIP_Bool, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_pseudoAncestors_addToEdge(SCIP*, int, int, GRAPH*);
extern SCIP_RETCODE   graph_pseudoAncestors_addToNode(SCIP*, int, int, GRAPH*);
extern SCIP_RETCODE   graph_free_pseudoAncestorsBlock(SCIP*, int, GRAPH*);
//extern SCIP_RETCODE   graph_checkConflict1_pseudoAncestors(SCIP*, const GRAPH*, int, SCIP_Bool*);
extern int            graph_pseudoAncestorsGetHashArraySize(const PSEUDOANS*);
extern void           graph_pseudoAncestors_hashNode(const PSEUDOANS*, int, int*);
extern void           graph_pseudoAncestors_hashEdge(const PSEUDOANS*, int, int*);
extern void           graph_pseudoAncestors_unhashNode(const PSEUDOANS*, int, int*);
extern void           graph_pseudoAncestors_unhashEdge(const PSEUDOANS*, int, int*);
extern void           graph_pseudoAncestors_hashNodeDirty(const PSEUDOANS*, int, SCIP_Bool, SCIP_Bool*, int*);
extern void           graph_pseudoAncestors_hashEdgeDirty(const PSEUDOANS*, int, SCIP_Bool, SCIP_Bool*, int*);
extern void           graph_pseudoAncestors_unhashNodeDirty(const PSEUDOANS*, int, int*);
extern void           graph_pseudoAncestors_unhashEdgeDirty(const PSEUDOANS*, int, int*);
extern SCIP_Bool      graph_pseudoAncestors_edgeIsHashed(const PSEUDOANS*, int, const int*);
extern SCIP_Bool      graph_pseudoAncestors_nodeIsHashed(const PSEUDOANS*, int, const int*);
extern SCIP_Bool      graph_pseudoAncestors_edgesInConflict(SCIP*, const GRAPH*, const int*, int);
extern SCIP_Bool      graph_valid_pseudoAncestors(SCIP*, const GRAPH*);
/* Fixed components */
extern SCIP_RETCODE   graph_init_fixed(SCIP*, GRAPH*);
extern void           graph_free_fixed(SCIP*, GRAPH*);
extern void           graph_free_fixedEdgesOnly(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_fixed_add(SCIP*, IDX*, const int*, int, GRAPH*);
extern SCIP_RETCODE   graph_fixed_addEdge(SCIP*, int, GRAPH*);
extern SCIP_RETCODE   graph_fixed_addNodePc(SCIP*, int, GRAPH*);
extern SCIP_RETCODE   graph_fixed_moveNodePc(SCIP*, int, GRAPH*);
extern SCIP_RETCODE   graph_copyFixed(SCIP*, GRAPH*, SCIP_Bool, GRAPH*);
extern IDX*           graph_getFixedges(const GRAPH*);
extern const int*     graph_getFixpseudonodes(SCIP*, const GRAPH*);
extern int            graph_getNfixpseudonodes(const GRAPH*);

/* graph_util.c
 */
extern SCIP_RETCODE   graph_termsReachable(SCIP*, const GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_findCentralTerminal(SCIP*, const GRAPH*, int, int*);
extern SCIP_RETCODE   graph_trail_arr(SCIP*, const GRAPH*, int, SCIP_Bool*);
extern SCIP_RETCODE   graph_trail_costAware(SCIP*, const GRAPH*, int, SCIP_Bool*);
/* Dijkstra heap: */
extern SCIP_RETCODE graph_heap_create(SCIP*, int, int*, DENTRY*, DHEAP**);
extern void   graph_heap_free(SCIP*, SCIP_Bool, SCIP_Bool, DHEAP**);
extern void   graph_heap_deleteMin(int*, SCIP_Real*, DHEAP*);
extern void   graph_heap_deleteMinGetNode(int*, DHEAP*);
extern int    graph_heap_deleteMinReturnNode(DHEAP*);
extern void   graph_heap_clean(SCIP_Bool, DHEAP*);
extern void   graph_heap_correct(int, SCIP_Real, DHEAP*);
extern SCIP_Bool graph_heap_isClean(const DHEAP*);
/* CSR storage: */
extern SCIP_RETCODE   graph_init_csr(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_init_csrWithEdgeId(SCIP*, GRAPH*);
extern void           graph_free_csr(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_csr_alloc(SCIP*, int, int, CSR**);
extern SCIP_RETCODE   graph_csr_allocWithEdgeId(SCIP*, int, int, CSR**);
extern void           graph_csr_copy(const CSR*, CSR*);
extern void           graph_csr_build(const GRAPH*, const SCIP_Real*, CSR*);
extern void           graph_csr_buildCosts(const GRAPH*, const CSR*, const SCIP_Real*, SCIP_Real* RESTRICT);
extern void           graph_csr_chgCosts(const GRAPH*, const SCIP_Real*, CSR*);
extern void           graph_csr_print(const CSR*);
extern int           graph_csr_getNedges(const CSR*);
extern void           graph_csr_free(SCIP*, CSR**);
extern SCIP_Bool      graph_csr_costsAreInSync(const GRAPH*, const CSR*, const SCIP_Real*);
extern SCIP_Bool      graph_csr_isValid(const CSR*, SCIP_Bool verbose);
extern SCIP_Bool      graph_valid_csr(const GRAPH*, SCIP_Bool verbose);
extern void           graph_buildOrgNodesToReducedMap(const GRAPH*, int*);
/* CSR depository: */
extern SCIP_RETCODE   graph_csrdepo_init(SCIP*, int, int, CSRDEPO**);
extern void           graph_csrdepo_free(SCIP*, CSRDEPO**);
extern void           graph_csrdepo_print(const CSRDEPO*);
extern void           graph_csrdepo_getCSR(const CSRDEPO*, int, CSR*);
extern void           graph_csrdepo_getTopCSR(const CSRDEPO*, CSR*);
extern int            graph_csrdepo_getNcsrs(const CSRDEPO*);
extern int            graph_csrdepo_getDataSize(const CSRDEPO*);
extern void           graph_csrdepo_clean(CSRDEPO*);
extern void           graph_csrdepo_removeTop(CSRDEPO*);
extern void           graph_csrdepo_addEmptyTop(CSRDEPO*, int, int);
extern void           graph_csrdepo_addEmptyTopTree(CSRDEPO*, int);
extern SCIP_Bool      graph_csrdepo_isEmpty(const CSRDEPO*);
extern SCIP_Bool      graph_csrdepo_hasEmptyTop(const CSRDEPO*);
extern void           graph_csrdepo_getEmptyTop(const CSRDEPO*, CSR*);
extern void           graph_csrdepo_emptyTopSetMarked(CSRDEPO*);

/* Dynamic CSR storage: */
extern SCIP_RETCODE   graph_init_dcsr(SCIP*, GRAPH*);
extern void           graph_update_dcsr(SCIP*, GRAPH*);
extern void           graph_dcsr_deleteEdge(DCSR*, int, int);
extern void           graph_dcsr_deleteEdgeBi(SCIP*, DCSR*, int);
extern void           graph_free_dcsr(SCIP*, GRAPH*);
extern SCIP_Bool      graph_valid_dcsr(const GRAPH*, SCIP_Bool verbose);
/* misc: */
extern SCIP_RETCODE  graph_dijkLimited_init(SCIP*, const GRAPH*, DIJK**);
extern SCIP_RETCODE  graph_dijkLimited_initPcShifts(SCIP*, const GRAPH*, DIJK*);
extern void          graph_dijkLimited_reset(const GRAPH*, DIJK*);
extern void          graph_dijkLimited_clean(const GRAPH*, DIJK*);
extern void          graph_dijkLimited_free(SCIP*, DIJK**);

/* graph_base.c
 */
extern void   graph_mark(GRAPH*);
extern void   graph_show(const GRAPH*);
extern void   graph_uncover(GRAPH*);
extern void   graph_free(SCIP*, GRAPH**, SCIP_Bool);
extern void   graph_freeHistory(SCIP*, GRAPH*);
extern void   graph_freeHistoryDeep(SCIP*, GRAPH*);
extern void   graph_getIsTermArray(const GRAPH*, SCIP_Bool*);
extern void   graph_getTerms(const GRAPH*, int*);
extern SCIP_RETCODE   graph_getTermsRandom(SCIP*, const GRAPH*, int*);
extern void   graph_getEdgeCosts(const GRAPH*, SCIP_Real* RESTRICT, SCIP_Real* RESTRICT);
extern void   graph_getEdgeRevCosts(const GRAPH*, SCIP_Real* RESTRICT);
extern void   graph_getCsr(const GRAPH*, int* RESTRICT, int* RESTRICT, int* RESTRICT, int*);
extern SCIP_RETCODE   graph_resize(SCIP*, GRAPH*, int, int, int);
extern SCIP_RETCODE   graph_copy(SCIP*, const GRAPH*, GRAPH**);
extern SCIP_RETCODE   graph_copyPseudoAncestors(SCIP*, const GRAPH*, GRAPH*);
extern SCIP_RETCODE   graph_copyData(SCIP*, const GRAPH*, GRAPH*);
extern SCIP_RETCODE   graph_pack(SCIP*, GRAPH*, GRAPH**, REDSOL*, SCIP_Bool);
extern SCIP_RETCODE   graph_init(SCIP*, GRAPH**, int, int, int);
extern SCIP_RETCODE   graph_initHistory(SCIP*, GRAPH*);
extern SCIP_Bool      graph_isMarked(const GRAPH*);
extern SCIP_Bool      graph_isSetUp(const GRAPH*);
extern SCIP_RETCODE   graph_buildCompleteGraph(SCIP*, GRAPH**, int);
extern SCIP_Bool graph_valid(SCIP*, const GRAPH*);
extern SCIP_Bool graph_validInput(SCIP*, const GRAPH*);
extern SCIP_Bool graph_knotIsNWLeaf(const GRAPH*, int);

/* graph_sub.c
 */
extern SCIP_RETCODE   graph_subgraphExtract(SCIP*, GRAPH*, SUBINOUT*, GRAPH**);
extern SCIP_RETCODE   graph_subinoutInit(SCIP*, const GRAPH*, SUBINOUT**);
extern void           graph_subinoutFree(SCIP*, SUBINOUT**);
extern void           graph_subinoutClean(SCIP*, SUBINOUT*);
extern SCIP_Bool      graph_subinoutUsesNewHistory(const SUBINOUT*);
extern SCIP_RETCODE   graph_subinoutActivateEdgeMap(const GRAPH*, SUBINOUT*);
extern void           graph_subinoutActivateNewHistory(SUBINOUT*);
const int*            graph_subinoutGetSubToOrgEdgeMap(const SUBINOUT*);
const int*            graph_subinoutGetSubToOrgNodeMap(const SUBINOUT*);
const int*            graph_subinoutGetOrgToSubNodeMap(const SUBINOUT*);
const int*            graph_subinoutGetContractionRecord(const SUBINOUT*);
int                   graph_knot_getContractionRecordAncestor(int, const SUBINOUT*);
extern SCIP_RETCODE   graph_subgraphCompleteNewHistory(SCIP*, const int*, GRAPH*, GRAPH*);
extern SCIP_RETCODE   graph_subgraphReinsert(SCIP*, SUBINOUT*, GRAPH*, GRAPH**);
extern void           graph_subgraphFree(SCIP*, GRAPH**);

/* graph_grid.c
 */
extern SCIP_RETCODE   graph_grid_create(SCIP*, GRAPH**, int**, int, int, int);
extern SCIP_RETCODE   graph_obstgrid_create(SCIP*, GRAPH**, int**, int**, int, int, int, int);
extern SCIP_RETCODE   graph_grid_coordinates(SCIP*, int**, int**, int*, int, int);

/* graph_edge.c
 */
extern void   graph_edge_add(SCIP*, GRAPH*, int, int, double, double);
extern void   graph_edge_addBi(SCIP*, GRAPH*, int, int, double);
extern void   graph_edge_addSubgraph(SCIP*, const GRAPH*, const int*, int, GRAPH*);
extern void   graph_edge_del(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_edge_delFull(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_edge_delBlocked(SCIP*, GRAPH*, const SCIP_Bool*, SCIP_Bool);
extern void   graph_edge_delHistory(SCIP*, GRAPH*, int);
extern void   graph_edge_hide(GRAPH*, int);
extern int    graph_edge_redirect(SCIP*, GRAPH*, int, int, int, SCIP_Real, SCIP_Bool, SCIP_Bool);

/* graph_node.c
 */
extern SCIP_RETCODE   graph_knot_contract(SCIP*, GRAPH*, int*, int, int);
extern SCIP_RETCODE   graph_knot_contractFixed(SCIP*, GRAPH*, int*, int, int, int);
extern SCIP_RETCODE   graph_knot_contractLowdeg2High(SCIP*, GRAPH*, int*, int, int);
extern SCIP_RETCODE   graph_knot_replaceDeg2(SCIP*, int, SCIP_Real, int, GRAPH*, SCIP_Bool*);
extern SCIP_RETCODE   graph_edge_reinsert(SCIP*, GRAPH*, int, int, int, SCIP_Real, int, SINGLETONANS*, SINGLETONANS*, int*, SCIP_Bool*);
extern SCIP_Bool graph_knot_isInRange(const GRAPH*, int);
extern void   graph_knot_add(GRAPH*, int);
extern void   graph_knot_chg(GRAPH*, int, int);
extern void   graph_knot_del(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_knot_delFull(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_knot_contract_dir(GRAPH*, int, int);

/* graph_delpseudo.c
 */
extern SCIP_RETCODE   graph_edge_delPseudo(SCIP*, GRAPH*, const SCIP_Real*, const SCIP_Real*, const SCIP_Real*, int, SCIP_Real*, SCIP_Bool*);
extern SCIP_RETCODE   graph_edge_delPseudoPath(SCIP*, GRAPH*, int, int, int, SCIP_Real*);
extern SCIP_RETCODE   graph_knot_delPseudo(SCIP*, GRAPH*, const SCIP_Real*, const SCIP_Real*, const SCIP_Real*, int, REDCOST*, SCIP_Bool*);
extern SCIP_RETCODE   graph_knot_delPseudoCheckIfPossible(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, const SCIP_Real*, int, SCIP_Bool*);

/* graph_stats.c
 */
extern SCIP_RETCODE   graph_printEdgeConflicts(SCIP*, const GRAPH*);
extern void   graph_edge_printInfo(const GRAPH*, int);
extern SCIP_Bool graph_edge_isBlocked(const GRAPH*, int);
extern SCIP_Bool graph_edge_isDeleted(const GRAPH*, int);
extern SCIP_Bool graph_edge_isInRange(const GRAPH*, int);
extern SCIP_Bool graph_hasMultiEdges(SCIP*, const GRAPH*, SCIP_Bool);
extern SCIP_Bool graph_isAlmostUniform(const GRAPH*);
extern void   graph_knot_printInfo(const GRAPH*, int);
extern void   graph_printInfo(const GRAPH*);
extern void   graph_printInfoReduced(const GRAPH*);
extern int    graph_get_nNodes(const GRAPH*);
extern int    graph_get_nEdges(const GRAPH*);
extern int    graph_get_nTerms(const GRAPH*);
extern void   graph_get_nVET(const GRAPH*, int*, int*, int*);
extern SCIP_Real graph_get_avgDeg(const GRAPH*);
extern SCIP_Bool graph_typeIsSpgLike(const GRAPH*);
extern SCIP_Bool graph_typeIsUndirected(const GRAPH*);
extern SCIP_Bool graph_typeIsDirected(const GRAPH*);


/* graph_trans.c
 */
extern SCIP_RETCODE   graph_transNw(SCIP*, PRESOL*, GRAPH*);
extern SCIP_RETCODE   graph_transNw2sap(SCIP*, PRESOL*, GRAPH*);
extern SCIP_RETCODE   graph_transNw2pc(SCIP*, PRESOL*, GRAPH*);
extern void   graph_transGstpClean(PRESOL*, GRAPH*);
extern SCIP_RETCODE   graph_transPc(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_transPc2Spg(SCIP*, PRESOL*, GRAPH*);
extern SCIP_RETCODE   graph_transRpc(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_transRpc2SpgTrivial(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_transRpc2FixedProper(SCIP*, PRESOL*, GRAPH*);
extern SCIP_RETCODE   graph_transMw(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_transRmw(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_transPcmw2rooted(SCIP*, GRAPH*, SCIP_Real, SCIP_Bool);
extern SCIP_RETCODE   graph_transPcGetSap(SCIP*, GRAPH*, GRAPH**, SCIP_Real*);
extern SCIP_RETCODE   graph_transPcGetRsap(SCIP*, GRAPH*, GRAPH**, const int*, int, int);
extern SCIP_RETCODE   graph_transRpcGetSpg(SCIP*, const GRAPH*, SCIP_Real, SCIP_Real*, int**, GRAPH**);
extern SCIP_Bool      graph_transRpcToSpgIsStable(const GRAPH*, SCIP_Real);


/* graph_pcbase.c
 */
extern void   graph_pc_knotToNonTermProperty(GRAPH*, int);
extern void   graph_pc_knotToFixedTermProperty(GRAPH*, int);
extern void   graph_pc_termToNonLeafTerm(SCIP*, GRAPH*, int, SCIP_Bool);
extern void   graph_pc_termToNonTerm(SCIP*, GRAPH*, int);
extern void   graph_pc_fixedTermToNonTerm(SCIP*, GRAPH*, int);
extern void   graph_pc_knotToFixedTerm(SCIP*, GRAPH*, int, SCIP_Real*);
extern void   graph_pc_nonTermToFixedTerm(GRAPH*, int, SCIP_Real*);
extern void   graph_pc_updateSubgraphEdge(const GRAPH*, const int*, int, GRAPH*);
extern void   graph_pc_enforcePseudoTerm(SCIP*, GRAPH*, int);
extern void   graph_pc_enforceNonLeafTerm(SCIP*, GRAPH*, int);
extern SCIP_Bool graph_pc_nonLeafTermIsEnforced(SCIP*, const GRAPH*, int);
extern void   graph_pc_enforceNode(SCIP*, GRAPH*, int, SCIP_Real*);
extern void   graph_pc_subtractPrize(SCIP*, GRAPH*, SCIP_Real, int);
extern void   graph_pc_knotChgPrize(GRAPH*, SCIP_Real, int);
extern void   graph_pc_2org(SCIP*, GRAPH*);
extern void   graph_pc_2trans(SCIP*, GRAPH*);
extern void   graph_pc_2orgcheck(SCIP*, GRAPH*);
extern void   graph_pc_2transcheck(SCIP*, GRAPH*);
extern int    graph_pc_getNorgEdges(const GRAPH*);
extern void   graph_pc_getReductionRatios(const GRAPH*, SCIP_Real*, SCIP_Real*);
extern void   graph_pc_getOrgCosts(SCIP*, const GRAPH*, SCIP_Real*);
extern void   graph_pc_getOrgCostsCsr(SCIP*, const GRAPH*, CSR*);
extern void   graph_pc_markOrgGraph(GRAPH*);
extern void   graph_pc_adaptSap(SCIP_Real, GRAPH*, SCIP_Real*);
extern void   graph_pc_presolExit(SCIP*, GRAPH*);
extern void   graph_pc_getBiased(const GRAPH*, SCIP_Real* RESTRICT, SCIP_Real* RESTRICT);
extern void   graph_pc_termMarkProper(const GRAPH*, int*);
extern void   graph_pc_shiftNonLeafCosts(SCIP*, GRAPH*);
extern SCIP_Real graph_pc_getNonLeafTermOffset(SCIP*, const GRAPH*);
extern SCIP_RETCODE   graph_pc_presolInit(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_initSubgraph(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_finalizeSubgraph(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_pc_initPrizes(SCIP*, GRAPH*, int);
extern SCIP_RETCODE   graph_pc_initTerm2Edge(SCIP*, GRAPH*, int);
extern SCIP_RETCODE   graph_pc_contractNodeAncestors(SCIP*, GRAPH*, int, int, int);
extern SCIP_RETCODE   graph_pc_contractEdge(SCIP*, GRAPH*, int*, int, int, int);
extern SCIP_RETCODE   graph_pc_contractEdgeUnordered(SCIP*, GRAPH*, int*, int, int);
extern int    graph_pc_deleteTerm(SCIP*, GRAPH*, int, SCIP_Real*);
extern int    graph_pc_realDegree(const GRAPH*, int, SCIP_Bool);
extern int    graph_pc_getRoot2PtermEdge(const GRAPH*, int);
extern int    graph_pc_nFixedTerms(const GRAPH*);
extern int    graph_pc_nNonFixedTerms(const GRAPH*);
extern int    graph_pc_nProperPotentialTerms(const GRAPH*);
extern int    graph_pc_nNonLeafTerms(const GRAPH*);
extern int    graph_pc_getTwinTerm(const GRAPH*, int);
extern SCIP_Bool graph_pc_costsEqualOrgCosts(SCIP*, const GRAPH*, const SCIP_Real*);
extern SCIP_Bool graph_pc_edgeIsExtended(const GRAPH*, int);
extern SCIP_Bool graph_pc_knotIsFixedTerm(const GRAPH*, int);
extern SCIP_Bool graph_pc_knotIsPropPotTerm(const GRAPH*, int);
extern SCIP_Bool graph_pc_knotHasMaxPrize(const GRAPH*, int);
extern SCIP_Bool graph_pc_knotIsDummyTerm(const GRAPH*, int);
extern SCIP_Bool graph_pc_termIsNonLeafTerm(const GRAPH*, int);
extern SCIP_Bool graph_pc_knotIsNonLeafTerm(const GRAPH*, int);
extern SCIP_Bool graph_pc_evalTermIsNonLeaf(SCIP*, const GRAPH*, int);
extern SCIP_Bool graph_pc_term2edgeIsConsistent(SCIP*, const GRAPH*);
extern SCIP_Bool graph_pc_transOrgAreConistent(SCIP*, const GRAPH*, SCIP_Bool);
extern SCIP_Real graph_pc_getPosPrizeSum(SCIP*, const GRAPH*);
extern SCIP_Bool graph_pc_isPc(const GRAPH*);
extern SCIP_Bool graph_pc_isMw(const GRAPH*);
extern SCIP_Bool graph_pc_isPcMw(const GRAPH*);
extern SCIP_Bool graph_pc_isRootedPcMw(const GRAPH*);
extern SCIP_Bool graph_pc_isUnrootedPcMw(const GRAPH*);
extern SCIP_Real graph_pc_solGetObj(SCIP*, const GRAPH*, const int*, SCIP_Real);


/* graph_path.c
 */
extern void   graph_path_exit(SCIP*, GRAPH*);
extern void   graph_path_exec(SCIP*, const GRAPH*, int, int, const SCIP_Real*, PATH*);
extern void   graph_pathInLimitedExec(const GRAPH*, const SCIP_Real*, const SCIP_Bool*, int, DIJK*, SCIP_Real*);
extern void   graph_path_execX(SCIP*, const GRAPH*, int, const SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_invroot(SCIP*, const GRAPH*, int, const SCIP_Real*, SCIP_Real*, int*);
extern void   graph_path_st(const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern SCIP_RETCODE graph_path_st_dc(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*, int*, STP_Bool*);
extern void   graph_path_st_rpcmw(GRAPH*, SCIP_Real*, int*, const SCIP_Real*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern SCIP_RETCODE   graph_path_st_brmwcs(SCIP*, GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*, SCIP_Bool*);
extern void   graph_path_st_pcmw(GRAPH*, SCIP_Real*, int*, const SCIP_Real*, const SCIP_Real*, SCIP_Bool, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_full(GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_reduce(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real*, int*, int, STP_Bool*);
extern void   graph_path_st_pcmw_extend(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Bool, PATH*, STP_Bool*, SCIP_Bool*);
extern void   graph_path_st_pcmw_extendBiased(SCIP*, GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, STP_Bool*, SCIP_Bool*);
extern void   graph_path_st_pcmw_extendOut(SCIP*, const GRAPH*, int, STP_Bool*, SCIP_Real*, int*, STP_Bool*, DHEAP*, SCIP_Bool*);
extern void   graph_pathHeapAdd(const PATH*, int, int* RESTRICT, int* RESTRICT, int*);
extern void   graph_path_PcMwSd(SCIP*, const GRAPH*, PATH*, SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE   graph_path_init(SCIP*, GRAPH*);
extern SCIP_Bool      graph_path_exists(const GRAPH*);
extern void   graph_sdPaths(const GRAPH*, PATH*, const SCIP_Real*, SCIP_Real, int*, int*, int*, int*, int, int, int);

/* graph_tpath.c
 */
extern void   graph_add2ndTermPaths(const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*);
extern void   graph_add3rdTermPaths(const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*);
extern void   graph_add4thTermPaths(const GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*);
extern void   graph_get2nextTermPaths(GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*);
extern void   graph_get3nextTermPaths(GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*);
extern void   graph_get4nextTermPaths(GRAPH*, const SCIP_Real*, const SCIP_Real*, PATH*, int*, int*);
extern SCIP_RETCODE   graph_get4nextTTerms(SCIP*, GRAPH*, const SCIP_Real*, PATH*, int*, int*, int*);
extern SCIP_RETCODE     graph_tpathsInit(SCIP*, GRAPH*, TPATHS**);
extern SCIP_RETCODE     graph_tpathsInitBiased(SCIP*, const SDPROFIT*, GRAPH*, TPATHS**);
extern SCIP_RETCODE     graph_tpathsRecomputeBiased(const SDPROFIT*, GRAPH*, TPATHS*);
extern SCIP_RETCODE     graph_tpathsRepair(SCIP*, int, SCIP_Bool, const GRAPH*, TPATHS*);
extern SCIP_RETCODE     graph_tpathsRepairSetUp(const GRAPH*, TPATHS*);
extern void             graph_tpathsFree(SCIP*, TPATHS**);
extern void             graph_tpathsPrintForNode(const GRAPH*, const SDPROFIT*, const TPATHS*, int);
extern void graph_tpathsGetProfitNodes(SCIP*, const GRAPH*, const TPATHS*, const SDPROFIT*, int, int, STP_Vectype(int));
extern void             graph_tpathsAdd1st(const GRAPH*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsAdd2nd(const GRAPH*, const SCIP_Real*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsAdd3rd(const GRAPH*, const SCIP_Real*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsAdd4th(const GRAPH*, const SCIP_Real*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsSetAll2(GRAPH*, const SCIP_Real*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsSetAll3(GRAPH*, const SCIP_Real*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsSetAll4(GRAPH*, const SCIP_Real*, const SCIP_Real*, const SDPROFIT*, TPATHS*);
extern void             graph_tpathsGetClosestTerm(const GRAPH*, const TPATHS*, int,
                           int* RESTRICT, int* RESTRICT, SCIP_Real* RESTRICT);
extern void             graph_tpathsGet3CloseTerms(const GRAPH*, const TPATHS*, int, SCIP_Real,
                           int* RESTRICT, int* RESTRICT, SCIP_Real* RESTRICT, int* RESTRICT);
extern void             graph_tpathsGet4CloseTerms(const GRAPH*, const TPATHS*, int, SCIP_Real,
                           int* RESTRICT, int* RESTRICT, SCIP_Real* RESTRICT, int* RESTRICT);
extern void             graph_tpathsGet4CloseTermsLE(const GRAPH*, const TPATHS*, int, SCIP_Real,
                           int* RESTRICT, int* RESTRICT, SCIP_Real* RESTRICT, int* RESTRICT);

/* graph_sdpath.c
 */
extern void   graph_sdStar(SCIP*, const GRAPH*, SCIP_Bool, int, int, int*, SCIP_Real*, int*, int*, DHEAP*, STP_Bool*, SCIP_Bool*);
extern SCIP_RETCODE   graph_sdStarBiased(SCIP*, const GRAPH*, const SDPROFIT*, int, int*, DIJK*, SCIP_Bool*);
extern SCIP_RETCODE   graph_sdCloseNodesBiased(SCIP*, const GRAPH*, const SDPROFIT*, int, DIJK*);
extern SCIP_Bool   graph_sdWalksConnected(SCIP*, const GRAPH*, const int*, const SCIP_Real*, const STP_Bool*, int, int, SCIP_Real*, int*, int*, STP_Bool*, SCIP_Bool);
extern SCIP_Bool graph_sdWalks(SCIP*, const GRAPH*, const SCIP_Real*, const int*, SCIP_Real, int, int, int, SCIP_Real*, int*, int*, int*, int*, STP_Bool*);
extern SCIP_Bool graph_sdWalks_csr(SCIP*, const GRAPH*, const int*, SCIP_Real, int, int, int, SCIP_Real*, int*, int*, DHEAP*, STP_Bool*);
extern SCIP_Bool graph_sdWalksTriangle(SCIP*, const GRAPH*, const int*, const int*, SCIP_Real, int, int, int, SCIP_Real*, SCIP_Real*, int*, int*, DHEAP*, STP_Bool*);
extern SCIP_Bool graph_sdWalksExt(SCIP*, const GRAPH*, const SCIP_Real*, SCIP_Real, int, int, int, int, SCIP_Real*, int*, int*, int*, int*, int*, int*, STP_Bool*);
extern SCIP_Bool graph_sdWalksExt2(SCIP*, const GRAPH*, const SCIP_Real*, const int*, SCIP_Real, int, int, int, int, SCIP_Real*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, STP_Bool*);
extern SCIP_RETCODE   graph_sdComputeCliqueStar(SCIP*, const GRAPH*, const SDPROFIT*, SDCLIQUE*);


/* graph_vnoi.c
 */
extern SCIP_RETCODE graph_vnoiInit(SCIP*, const GRAPH*, SCIP_Bool, VNOI**);
extern void graph_vnoiFree(SCIP*, VNOI**);
extern SCIP_RETCODE graph_vnoiCompute(SCIP*, const GRAPH*, VNOI*);
extern SCIP_RETCODE graph_vnoiComputeImplied(SCIP*, const GRAPH*, const SDPROFIT*, VNOI*);
extern void   graph_voronoi(const GRAPH*, const SCIP_Real*, const SCIP_Real*, const STP_Bool*, int*, PATH*);
extern void   graph_voronoiTerms(const GRAPH*, const SCIP_Bool*, int* RESTRICT, PATH* RESTRICT);
extern void   graph_voronoiRepair(SCIP*, const GRAPH*, const SCIP_Real*, const SCIP_Real*, int*, int*, PATH*, int*, int, UF*);
extern void   graph_voronoiRepairMult(SCIP*, const GRAPH*, const SCIP_Real*, const STP_Bool*, int* RESTRICT, int* RESTRICT, int* RESTRICT, int* RESTRICT, UF* RESTRICT, PATH* RESTRICT);
extern void   graph_voronoiWithRadiusMw(SCIP*, const GRAPH*, PATH*, const SCIP_Real*, SCIP_Real*, int*, int*, int*);
extern void   graph_voronoiMw(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, int*, int*, int*);
extern void   graph_add1stTermPaths(const GRAPH*, const SCIP_Real*, PATH*, int*, int*);
extern void   voronoiSteinerTreeExt(SCIP*, const GRAPH*, SCIP_Real*, int*, STP_Bool*, PATH*);
extern SCIP_RETCODE   graph_voronoiExtend(SCIP*, const GRAPH*, const SCIP_Real*, PATH*, SCIP_Real**, int**, int**, STP_Bool*, int*, int*, int*, int, int, int);
extern SCIP_RETCODE   graph_voronoiWithDist(SCIP*, const GRAPH*, const SCIP_Bool*, const SCIP_Real*, double*, int*, int*, int*, int*, PATH*);
extern SCIP_RETCODE   graph_voronoiWithRadius(SCIP*, const GRAPH*, GRAPH*, PATH*, SCIP_Real*, SCIP_Real*, SCIP_Real*, int*, int*, int*);


/* graph_mcut.c
 */
#if 0
extern void   graph_mincut_exit(SCIP*, GRAPH*);
extern void   graph_mincut_exec(GRAPH*, int, int, int, const int*, int*, const int*, const int*, const int*, int);
extern SCIP_RETCODE   graph_mincut_init(SCIP*, GRAPH*);
#else
extern void   graph_mincut_exit(SCIP*, GRAPH*);
extern void   graph_mincut_exec(const GRAPH*, const int, const int, const int, const int, const int, const int*, const int*, int* RESTRICT, const int*, const int*, const int*, const SCIP_Bool);
extern SCIP_RETCODE   graph_mincut_init(SCIP*, GRAPH*);
extern SCIP_RETCODE   graph_mincut_reInit(SCIP*, int, int, GRAPH*);
extern SCIP_Bool  graph_mincut_isInitialized(const GRAPH*);
extern void    graph_mincut_setDefaultVals(GRAPH*);
#endif

/* graph_load.c
 */
extern SCIP_RETCODE graph_load(SCIP*, GRAPH**, const char*, PRESOL*);

/* graph_save.c
 */
extern void graph_save(SCIP*, const GRAPH*, const char*, FILETYPE);
extern void graph_writeReductionStats(const GRAPH*, const char*, const char*);
extern void graph_writeReductionRatioStats(const GRAPH*, const char*, const char*);
extern SCIP_RETCODE graph_writeReductionRatioStatsLive(SCIP*, GRAPH*, const char*);
extern SCIP_RETCODE graph_writeGml(const GRAPH*, const char*, const SCIP_Bool*);
extern SCIP_RETCODE graph_writeGmlSub(const GRAPH*, const char*, const SCIP_Bool*);
extern void graph_writeStp(SCIP*, const GRAPH*, FILE*, SCIP_Real);
extern void graph_writeStpByName(SCIP*, const GRAPH*, const char*, SCIP_Real);
extern void graph_writeStpOrg(SCIP*, const GRAPH*, const char*);


/* validate.c
 */
extern SCIP_RETCODE    SCIPStpValidateSol(SCIP*, const GRAPH*, const double*, SCIP_Bool, SCIP_Bool*);



#ifdef __cplusplus
}
#endif


#endif /* !APPLICATIONS_STP_GRAPH_H_ */
