/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_nlp.h,v 1.10 2010/04/21 18:23:18 bzfviger Exp $"

/**@file   heur_nlp.h
 * @brief  NLP local search primal heuristic
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef HEUR_NLP_H_
#define HEUR_NLP_H_

#include "scip/scip.h"
#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** NLPI initialization method of NLP heuristic
 * 
 *  input:
 *  - scip            : SCIP main data structure
 *  - nlpi            : NLP solver interface
 *  - problem         : NLP solver problem
 *  - varmap          : variable mapping SCIP to NLPI
 */
#define SCIP_DECL_HEURNLPNLPIINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_NLPI* nlpi, SCIP_NLPIPROBLEM* problem, SCIP_HASHMAP* varmap)

/** method to check whether a constraint handler has nonlinear constraints
 * 
 *  input:
 *  - scip            : SCIP main data structure
 *  - fixedint        : should the case be considered where all discrete variable are fixed?
 *  - result          : buffer to store whether the NLPIINIT routine would add constraints to the NLP
 */
#define SCIP_DECL_HEURNLPHAVECONS(x) SCIP_RETCODE x (SCIP* scip, SCIP_Bool fixedint, SCIP_Bool* result)

/** creates the NLP local search primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurNlp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes an NLPI initialization method into the NLP heuristic
 * can be used by constraint handlers to register a function that inserts their constraints into an NLPI */
extern
SCIP_RETCODE SCIPincludeHeurNlpNlpiInit(
   SCIP*                   scip,               /**< SCIP data structure */
   SCIP_DECL_HEURNLPHAVECONS((*havecons)),     /**< method to call for checking if potential constraints for the NLP are present */
   SCIP_DECL_HEURNLPNLPIINIT((*nlpiinit)),     /**< method to call for initializing NLP */
   const char*             conshdlrname        /**< name of the constraint handler */
   );

/** updates the starting point for the NLP heuristic
 * 
 * Is called by a constraint handler that handles nonlinear constraints when a check on feasibility of a solution fails.
 */
extern
SCIP_RETCODE SCIPheurNlpUpdateStartpoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< NLP heuristic */
   SCIP_SOL*             solcand,            /**< solution candidate */
   SCIP_Real             violation           /**< constraint violation of solution candidate */
   );

/** main procedure of the NLP heuristic */
extern
SCIP_RETCODE SCIPapplyNlpHeur(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< result data structure                                          */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver                                 */
   SCIP_Real             timelimit,          /**< time limit for NLP solver                                      */
   SCIP_Longint*         iterused            /**< buffer to store number of iterations used by NLP solver, or NULL if not of interest */
   );

#ifdef __cplusplus
}
#endif

#endif /*HEUR_NLP_H_*/
