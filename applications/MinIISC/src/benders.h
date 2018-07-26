/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benders.h
 * @brief  run Benders algorithm
 * @author Marc Pfetsch
 *
 * Run Benders algorithm using an oracle for solving the subproblems and solving the master problem to optimality.
 */

#ifndef __BENDERS_H__
#define __BENDERS_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Benders subproblem oracle solving status */
enum BENDERS_Status
{
   BENDERS_STATUS_UNKNOWN          =  0,     /**< the solving status is not yet known */
   BENDERS_STATUS_ADDEDCUT         =  1,     /**< a Benders cut has been added */
   BENDERS_STATUS_SUCESS           =  2,     /**< the solution is optimal, no further Benders cut has to be generated */
   BENDERS_STATUS_TIMELIMIT        =  3,     /**< the time limit has been reached */
   BENDERS_STATUS_USERINTERRUPT    =  4,     /**< the user has interrupted the solution of the subproblem */
   BENDERS_STATUS_ERROR            =  5      /**< an error occured during the solution of the subproblem */
};
typedef enum BENDERS_Status BENDERS_STATUS;

typedef struct BENDERS_Data BENDERS_DATA;    /**< user defined data to pass to the oracle */

/** user callback method for a Benders subproblem oracle
 *  input:
 *   - masterscip:       SCIP pointer of Benders master problem
 *   - nmastervars:      number of variables in master problem
 *   - mastervars:       variables in master problem
 *   - mastersolution:   solution of Benders master problem
 *   - data:             user data for oracle
 *   - timelimit:        time limit for subproblem
 *   - ntotalcuts:       total number of cuts
 *  output:
 *   - ncuts:            number of cuts added
 *   - status:           status
 *
 *  The oracle should take the given solution and possibly add a Benders Cut to the master problem.
 */
#define BENDERS_CUTORACLE(x) SCIP_RETCODE x (SCIP* masterscip, int nmastervars, SCIP_VAR** mastervars, SCIP_Real* mastersolution, BENDERS_DATA* data, SCIP_Real timelimit, SCIP_Longint ntotalcuts, int* ncuts, BENDERS_STATUS* status)



/** run Benders algorithm using an oracle for the subproblems */
extern
SCIP_RETCODE runBenders(
   SCIP*                 masterscip,         /**< master SCIP instance */
   BENDERS_CUTORACLE((*Oracle)),             /**< oracle for the Benders subproblem */
   BENDERS_DATA*         data,               /**< user data for oracle */
   SCIP_Real             timelimit,          /**< time limit read from arguments */
   SCIP_Real             memlimit,           /**< memory limit read from arguments */
   int                   dispfreq,           /**< display frequency */
   SCIP_Bool             usereopt,           /**< use reoptimization */
   SCIP_Bool             solvemasterapprox,  /**< Solve master problem approximately? */
   SCIP_Longint          masterstallnodes,   /**< stall nodes for master problem if solvemasterapprox is true */
   SCIP_Real             mastergaplimit,     /**< gap limit for master problem if solvemasterapprox is true */
   SCIP_VERBLEVEL        verblevel,          /**< verbosity level for output */
   SCIP_STATUS*          status              /**< status of optimization */
   );

#ifdef __cplusplus
}
#endif

#endif
