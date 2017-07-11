/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_benders.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for benders decomposition methods
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_BENDERS_H__
#define __SCIP_TYPE_BENDERS_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

enum SCIP_BendersEnfoType
{
    LP      = 1,
    RELAX   = 2,
    PSEUDO  = 3,
    CHECK   = 4
};
typedef enum SCIP_BendersEnfoType SCIP_BENDERSENFOTYPE;

typedef struct SCIP_Benders SCIP_BENDERS;           /**< variable benders data */
typedef struct SCIP_BendersData SCIP_BENDERSDATA;   /**< locally defined variable benders data */


/** copy method for benders plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 *  - valid           : was the copying process valid? 
 */
#define SCIP_DECL_BENDERSCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_Bool* valid)

/** destructor of variable benders to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** initialization method of variable benders (called after problem was transformed and benders is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** deinitialization method of variable benders (called before transformed problem is freed and benders is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSINITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** presolving deinitialization method of constraint handler (called after presolving has been finished)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSEXITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** solving process initialization method of variable benders (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The variable benders may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** solving process deinitialization method of variable benders (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The variable benders should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the variable benders itself
 */
#define SCIP_DECL_BENDERSEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** the execution method for Benders' decomposition. The execution method is called by the constraint handler and then
 *  solves each of the subproblems.
 *
 *  This method is called to verify a given solution with the Benders' decomposition subproblems.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - sol             : the solution that will be checked in the subproblem. Can be NULL.
 *  TODO: Need to update the parameters for this callback
 */
#define SCIP_DECL_BENDERSEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubXyz callback.
 *
 *  When creating the subproblems, they must be registered with the Benders' decomposition structure. This is done by
 *  calling SCIPaddBendersSubproblem.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - probnumber      : the subproblem problem number
 */
#define SCIP_DECL_BENDERSCREATESUB(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, int probnumber)

/** the solving method for a single Benders' decomposition subproblem. The solving methods are separated so that they
 *  can be called in parallel.
 *
 *  NOTE: The solving methods must be thread safe.
 *
 *  This method is called from within the execution method.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - sol             : the solution that will be checked in the subproblem. Can be NULL.
 *  - probnumber      : the subproblem problem number
 *  - infeasible      : pointer to store whether the problem is infeasible
 *  TODO: Need to update the parameters for this callback
 */
#define SCIP_DECL_BENDERSSOLVESUB(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_SOL* sol, int probnumber,\
  SCIP_Bool* infeasible)

/** the post-solve method for Benders' decomposition. The post-solve method is called after the subproblems have
 * been solved but before they are freed.
 *
 * This provides the opportunity for the user to clean up and data structures
 * that should not exist beyond the current iteration. Also, the user has full access to the master and subproblems. So
 * it is possible to construct solution for the master problem in the method.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - infeasible      : indicates whether at least one subproblem is infeasible
 *  TODO: Need to update the parameters for this callback
 */
#define SCIP_DECL_BENDERSPOSTSOLVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_Bool infeasible)

/** frees the subproblem so that it can be resolved in the next iteration. In the SCIP case, this involves freeing the
 *  transformed problem using SCIPfreeTransform()
 *
 *  NOTE: The freeing methods must be thread safe.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - probnumber      : the subproblem problem number
 *  TODO: Need to update the parameters for this callback
 */
#define SCIP_DECL_BENDERSFREESUB(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, int probnumber)

/** the variable mapping from the subproblem to the master problem.
 *
 *  This method is called when generating the cuts. The cuts are generated by using the solution to the subproblem to
 *  eliminate a solution to the master problem.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition structure
 *  - var             : the variable for which the corresponding variable in the master or subproblem is required
 *  - probnumber      : the number of the subproblem that the desired variable belongs to, -1 for the master problem
 */
#define SCIP_DECL_BENDERSGETVAR(x) SCIP_VAR* x (SCIP* scip, SCIP_BENDERS* benders, SCIP_VAR* var, int probnumber)

#ifdef __cplusplus
}
#endif

#endif
