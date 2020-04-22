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

/**@file   solstp.h
 * @brief  includes methods for Steiner tree problem solutions
 * @author Daniel Rehfeldt
 *
 * Methods for manipulating solutions (i.e. trees) to Steiner tree problems, such as pruning.
 * Also includes methods for obtaining information about solutions.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_SOLSTP_H_
#define APPLICATIONS_STP_SRC_SOLSTP_H_

#include "scip/scip.h"
#include "graph.h"
#include "shortestpath.h"

extern void solstp_pcConnectDummies(const GRAPH*, int, int* RESTRICT, STP_Bool* RESTRICT);
extern SCIP_Real solstp_pcGetObjCsr(const GRAPH*, const CSR*, const int*, const STP_Bool*);
extern int solstp_pcGetSolRoot(SCIP*, const GRAPH*, const STP_Bool*);

extern void   solstp_setNodeList(const GRAPH*, STP_Bool*, IDX*);
extern void   solstp_setVertexFromEdge(const GRAPH*, const int*, STP_Bool*);
extern void   solstp_print(const GRAPH*, const int*);
extern SCIP_RETCODE   solstp_markPcancestors(SCIP*, IDX**, const int*, const int*, int, STP_Bool*, STP_Bool*, int*, int*, int*);
extern SCIP_Bool solstp_isUnreduced(SCIP*, const GRAPH*, const int*);
extern SCIP_Bool solstp_isValid(SCIP*, const GRAPH*, const int*);
extern SCIP_Bool solstp_containsNode(const GRAPH*, const int*, int);
extern SCIP_Bool stpsol_pruningIsPossible(const GRAPH*, const int*, const STP_Bool*);
extern SCIP_Real solstp_getObj(const GRAPH*, const int*, SCIP_Real, int);
extern SCIP_Real solstp_getObjCsr(const GRAPH*, const CSR*, const int*, const STP_Bool*);
extern int       solstp_getNedges(const GRAPH*, const int*);
extern int       solstp_getNedgesBounded(const GRAPH*, const int*, int);
extern void      solstp_getTrivialSol(const GRAPH*, int*);
extern void      solstp_convertCsrToGraph(SCIP*, const GRAPH*, const CSR*, const int*, STP_Bool* RESTRICT, int* RESTRICT);
extern SCIP_RETCODE   solstp_getOrg(SCIP*, const GRAPH*, const GRAPH*, const int*, int*);
extern SCIP_RETCODE   solstp_reroot(SCIP*, GRAPH*, int*, int);
SCIP_RETCODE       solstp_prune(SCIP*, const GRAPH*, int*, STP_Bool*);
SCIP_RETCODE       solstp_pruneFromTmHeur(SCIP*, const GRAPH*, const SCIP_Real*, int* RESTRICT, STP_Bool* RESTRICT);
SCIP_RETCODE       solstp_pruneFromTmHeur_csr(SCIP*, const GRAPH*, SPATHS*, int* RESTRICT);
SCIP_RETCODE       solstp_pruneFromNodes(SCIP*, const GRAPH*, int*, STP_Bool*);
SCIP_RETCODE       solstp_pruneFromEdges(SCIP*, const GRAPH*, int*);


#endif /* APPLICATIONS_STP_SRC_SOLSTP_H_ */
