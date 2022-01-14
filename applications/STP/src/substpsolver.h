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

/**@file   substpsolver.h
 * @brief  Solver for Steiner tree (sub-) problems
 * @author Daniel Rehfeldt
 *
 * This file implements methods to solve Steiner tree subproblems to optimality, and to obtain an optimal Steiner tree.
 * Branch-and-cut or dynamic programming is used for solving the subproblems.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_SUBSTPSOLVER_H_
#define APPLICATIONS_STP_SRC_SUBSTPSOLVER_H_


#include "scip/scip.h"
#include "graph.h"
#include "probdata_stp.h"


/** Steiner tree sub-problem */
typedef struct sub_steiner_tree_problem SUBSTP;

#ifdef __cplusplus
extern "C" {
#endif


extern SCIP_RETCODE substpsolver_init(SCIP*, GRAPH*, SUBSTP**);
extern SCIP_RETCODE substpsolver_initDP(SCIP*, GRAPH*, SUBSTP**);
extern SCIP_RETCODE substpsolver_initBC(SCIP*, GRAPH*, SUBSTP**);
extern void         substpsolver_free(SCIP*, SUBSTP**);
extern SCIP_RETCODE substpsolver_transferHistory(const int*, GRAPH*, SUBSTP*);
extern SCIP_RETCODE substpsolver_initHistory(SUBSTP*);
extern SCIP_RETCODE substpsolver_solve(SCIP*, SUBSTP*, SCIP_Bool*);
extern int          substpsolver_getNsubedges(const SUBSTP*);
extern SCIP_RETCODE substpsolver_getSolution(SUBSTP*, int*);
extern SCIP_RETCODE substpsolver_setMute(SUBSTP*);
extern SCIP_RETCODE substpsolver_setProbIsIndependent(SUBSTP*);
extern SCIP_RETCODE substpsolver_setProbNoSubDP(SUBSTP*);
extern SCIP_RETCODE substpsolver_setProbFullPresolve(SUBSTP*);
extern SCIP_RETCODE substpsolver_getObjFromGraph(SCIP*, const GRAPH*, SCIP_Real*);


#ifdef __cplusplus
}
#endif

#endif /* APPLICATIONS_STP_SRC_SUBSTPSOLVER_H_ */
