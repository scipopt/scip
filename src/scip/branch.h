/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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
 *    branchrule      : the branching rule itself
 *    scip            : SCIP main data structure
 */
#define DECL_BRANCHFREE(x) RETCODE x (BRANCHRULE* branchrule, SCIP* scip)

/** initialization method of branching rule (called at problem creation)
 *
 *  input:
 *    branchrule      : the branching rule itself
 *    scip            : SCIP main data structure
 */
#define DECL_BRANCHINIT(x) RETCODE x (BRANCHRULE* branchrule, SCIP* scip)

/** deinitialization method of branching rule (called at problem destruction)
 *
 *  input:
 *    branchrule      : the branching rule itself
 *    scip            : SCIP main data structure
 */
#define DECL_BRANCHEXIT(x) RETCODE x (BRANCHRULE* branchrule, SCIP* scip)

/** branching execution method for fractional LP solutions
 *
 *  input:
 *    branchrule      : the branching rule itself
 *    scip            : SCIP main data structure
 *    result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : the current node was detected to be infeasible
 *    SCIP_BRANCHED   : branching was applied
 *    SCIP_REDUCEDDOM : a domain was reduced that rendered the actual LP solution infeasible
 *    SCIP_SEPARATED  : a cutting plane was generated (only if "lpvalid" is TRUE)
 *    SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXLP(x) RETCODE x (BRANCHRULE* branchrule, SCIP* scip, RESULT* result)

/** branching execution method for not completely fixed pseudo solutions
 *
 *  input:
 *    branchrule      : the branching rule itself
 *    scip            : SCIP main data structure
 *    result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : the current node was detected to be infeasible
 *    SCIP_BRANCHED   : branching was applied
 *    SCIP_REDUCEDDOM : a domain was reduced that rendered the actual pseudo solution infeasible
 *    SCIP_DIDNOTRUN  : the branching rule was skipped
 */
#define DECL_BRANCHEXPS(x) RETCODE x (BRANCHRULE* branchrule, SCIP* scip, RESULT* result)



#include "scip.h"
#include "retcode.h"
#include "set.h"
#include "var.h"



/* branching candidate storage methods */

extern
RETCODE SCIPbranchcandCreate(           /**< creates a branching candidate storage */
   BRANCHCAND**     branchcand,         /**< pointer to store branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   );

extern
RETCODE SCIPbranchcandFree(             /**< frees branching candidate storage */
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

extern
RETCODE SCIPbranchcandGetLPCands(       /**< gets branching candidates for LP solution branching (fractional variables) */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   );

extern
RETCODE SCIPbranchcandGetPseudoCands(   /**< gets branching candidates for pseudo solution branching (nonfixed variables) */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   );

extern
RETCODE SCIPbranchcandUpdateVar(        /**< updates branching candidate list for a given variable */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that changed its bounds */
   );




/* branching rules */

extern
RETCODE SCIPbranchruleCreate(           /**< creates a branching rule */
   BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE((*branchfree)),      /**< destructor of branching rule */
   DECL_BRANCHINIT((*branchinit)),      /**< initialise branching rule */
   DECL_BRANCHEXIT((*branchexit)),      /**< deinitialise branching rule */
   DECL_BRANCHEXLP((*branchexlp)),      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps)),      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

extern
RETCODE SCIPbranchruleFree(             /**< frees memory of branching rule */
   BRANCHRULE**     branchrule,         /**< pointer to branching rule data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPbranchruleInit(             /**< initializes branching rule */
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPbranchruleExit(             /**< deinitializes branching rule */
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPbranchruleExecLPSol(        /**< executes branching rule for fractional LP solution */
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
RETCODE SCIPbranchruleExecPseudoSol(    /**< executes branching rule for not completely fixed pseudo solution */
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
const char* SCIPbranchruleGetName(      /**< gets name of branching rule */
   BRANCHRULE*      branchrule          /**< branching rule */
   );

extern
int SCIPbranchruleGetPriority(          /**< gets priority of branching rule */
   BRANCHRULE*      branchrule          /**< branching rule */
   );

extern
BRANCHRULEDATA* SCIPbranchruleGetData(  /**< gets user data of branching rule */
   BRANCHRULE*      branchrule          /**< branching rule */
   );

extern
void SCIPbranchruleSetData(             /**< sets user data of branching rule; user has to free old data in advance! */
   BRANCHRULE*      branchrule,         /**< branching rule */
   BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   );

extern
Bool SCIPbranchruleIsInitialized(       /**< is branching rule initialized? */
   BRANCHRULE*      branchrule          /**< branching rule */
   );


#endif
