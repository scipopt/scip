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
typedef struct Branch BRANCH;           /**< branching method data structure */
typedef struct BranchData BRANCHDATA;   /**< branching method specific data */


/** destructor of branching method to free user data (called when SCIP is exiting)
 *
 *  input:
 *    branch          : the branching method itself
 *    scip            : SCIP main data structure
 */
#define DECL_BRANCHFREE(x) RETCODE x (BRANCH* branch, SCIP* scip)

/** initialization method of branching method (called at problem creation)
 *
 *  input:
 *    branch          : the branching method itself
 *    scip            : SCIP main data structure
 */
#define DECL_BRANCHINIT(x) RETCODE x (BRANCH* branch, SCIP* scip)

/** deinitialization method of branching method (called at problem destruction)
 *
 *  input:
 *    branch          : the branching method itself
 *    scip            : SCIP main data structure
 */
#define DECL_BRANCHEXIT(x) RETCODE x (BRANCH* branch, SCIP* scip)

/** branching execution method for fractional LP solutions
 *
 *  input:
 *    branch          : the branching method itself
 *    scip            : SCIP main data structure
 *    result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *    SCIP_BRANCHED   : branching was applied
 *    SCIP_REDUCEDDOM : a domain was reduced that rendered the actual LP solution infeasible
 *    SCIP_SEPARATED  : a cutting plane was generated (only if "lpvalid" is TRUE)
 *    SCIP_DIDNOTRUN  : the branching method was skipped
 */
#define DECL_BRANCHEXLP(x) RETCODE x (BRANCH* branch, SCIP* scip, RESULT* result)

/** branching execution method for not completely fixed pseudo solutions
 *
 *  input:
 *    branch          : the branching method itself
 *    scip            : SCIP main data structure
 *    result          : pointer to store the result of the branching call
 *
 *  possible return values for *result:
 *    SCIP_BRANCHED   : branching was applied
 *    SCIP_REDUCEDDOM : a domain was reduced that rendered the actual pseudo solution infeasible
 *    SCIP_DIDNOTRUN  : the branching method was skipped
 */
#define DECL_BRANCHEXPS(x) RETCODE x (BRANCH* branch, SCIP* scip, RESULT* result)



#include "scip.h"
#include "retcode.h"
#include "set.h"
#include "var.h"



/* branching candidate storage methods */

extern
RETCODE SCIPbranchcandCreate(           /**< creates a branching candidate storage */
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
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
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   );




/* branching methods */

extern
RETCODE SCIPbranchCreate(               /**< creates a branching method */
   BRANCH**         branch,             /**< pointer to store branching method */
   const char*      name,               /**< name of branching method */
   const char*      desc,               /**< description of branching method */
   DECL_BRANCHFREE((*branchfree)),      /**< destructor of branching method */
   DECL_BRANCHINIT((*branchinit)),      /**< initialise branching method */
   DECL_BRANCHEXIT((*branchexit)),      /**< deinitialise branching method */
   DECL_BRANCHEXLP((*branchexlp)),      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps)),      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHDATA*      branchdata          /**< branching method data */
   );

extern
RETCODE SCIPbranchFree(                 /**< frees memory of branching method */
   BRANCH**         branch,             /**< pointer to branching method data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPbranchInit(                 /**< initializes branching method */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPbranchExit(                 /**< deinitializes branching method */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPbranchExecLPSol(            /**< executes branching method for fractional LP solution */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
RETCODE SCIPbranchExecPseudoSol(        /**< executes branching method for not completely fixed pseudo solution */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
const char* SCIPbranchGetName(          /**< gets name of branching method */
   BRANCH*          branch              /**< branching method */
   );

extern
BRANCHDATA* SCIPbranchGetData(          /**< gets user data of branching method */
   BRANCH*          branch              /**< branching method */
   );

extern
void SCIPbranchSetData(                 /**< sets user data of branching method; user has to free old data in advance! */
   BRANCH*          branch,             /**< branching method */
   BRANCHDATA*      branchdata          /**< new branching method user data */
   );

extern
Bool SCIPbranchIsInitialized(           /**< is branching method initialized? */
   BRANCH*          branch              /**< branching method */
   );


#endif
