/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch.h
 * @brief  datastructures and methods for branching methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BRANCH_H__
#define __BRANCH_H__


typedef struct BranchCand BRANCHCAND;   /**< branching candidate storage */
typedef struct BranchRule BRANCHRULE;   /**< branching method data structure */
typedef struct BranchRuleData BRANCHRULEDATA; /**< branching method specific data */


/** destructor of branching method to free user data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    branchrule      : the branching rule itself
 */
#define DECL_BRANCHFREE(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** initialization method of branching rule (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    branchrule      : the branching rule itself
 */
#define DECL_BRANCHINIT(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** deinitialization method of branching rule (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    branchrule      : the branching rule itself
 */
#define DECL_BRANCHEXIT(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule)

/** branching execution method for fractional LP solutions
 *
 *  input:
 *    scip            : SCIP main data structure
 *    branchrule      : the branching rule itself
 *    result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : the current node was detected to be infeasible
 *    SCIP_BRANCHED   : branching was applied
 *    SCIP_REDUCEDDOM : a domain was reduced that rendered the actual LP solution infeasible
 *    SCIP_SEPARATED  : a cutting plane was generated
 *    SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXECLP(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule, RESULT* result)

/** branching execution method for not completely fixed pseudo solutions
 *
 *  input:
 *    scip            : SCIP main data structure
 *    branchrule      : the branching rule itself
 *    result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : the current node was detected to be infeasible
 *    SCIP_BRANCHED   : branching was applied
 *    SCIP_REDUCEDDOM : a domain was reduced that rendered the actual pseudo solution infeasible
 *    SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXECPS(x) RETCODE x (SCIP* scip, BRANCHRULE* branchrule, RESULT* result)



#include "scip.h"
#include "retcode.h"
#include "set.h"
#include "var.h"



/* branching candidate storage methods */


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




/* branching rules */

/** creates a branching rule */
extern
RETCODE SCIPbranchruleCreate(
   BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialise branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialise branching rule */
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
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** executes branching rule for not completely fixed pseudo solution */
extern
RETCODE SCIPbranchruleExecPseudoSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** gets name of branching rule */
extern
const char* SCIPbranchruleGetName(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets priority of branching rule */
extern
int SCIPbranchruleGetPriority(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets user data of branching rule */
extern
BRANCHRULEDATA* SCIPbranchruleGetData(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** sets user data of branching rule; user has to free old data in advance! */
extern
void SCIPbranchruleSetData(
   BRANCHRULE*      branchrule,         /**< branching rule */
   BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   );

/** is branching rule initialized? */
extern
Bool SCIPbranchruleIsInitialized(
   BRANCHRULE*      branchrule          /**< branching rule */
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
