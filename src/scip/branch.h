/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch.h,v 1.21 2003/12/01 16:14:27 bzfpfend Exp $"

/**@file   branch.h
 * @brief  internal methods for branching rules and branching candidate storage
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BRANCH_H__
#define __BRANCH_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_misc.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_scip.h"
#include "type_branch.h"
#include "pub_branch.h"



/*
 * branching candidate storage methods
 */

/** creates a branching candidate storage */
extern
RETCODE SCIPbranchcandCreate(
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

/** frees branching candidate storage */
extern
RETCODE SCIPbranchcandFree(
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

/** gets branching candidates for LP solution branching (fractional variables) */
extern
RETCODE SCIPbranchcandGetLPCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   );

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
extern
RETCODE SCIPbranchcandGetPseudoCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   );

/** updates branching candidate list for a given variable */
extern
RETCODE SCIPbranchcandUpdateVar(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that changed its bounds */
   );




/*
 * branching rules
 */

/** creates a branching rule */
extern
RETCODE SCIPbranchruleCreate(
   BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

/** frees memory of branching rule */
extern
RETCODE SCIPbranchruleFree(
   BRANCHRULE**     branchrule,         /**< pointer to branching rule data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes branching rule */
extern
RETCODE SCIPbranchruleInit(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** deinitializes branching rule */
extern
RETCODE SCIPbranchruleExit(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** executes branching rule for fractional LP solution */
extern
RETCODE SCIPbranchruleExecLPSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** executes branching rule for not completely fixed pseudo solution */
extern
RETCODE SCIPbranchruleExecPseudoSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of branching rule */
extern
void SCIPbranchruleSetPriority(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the branching rule */
   );




/*
 * branching methods
 */

/** calculates the branching score out of the downward and upward gain prediction */
extern
Real SCIPbranchGetScore(
   const SET*       set,                /**< global SCIP settings */
   Real             downgain,           /**< prediction of objective gain for branching downwards */
   Real             upgain              /**< prediction of objective gain for branching upwards */
   );


#endif
