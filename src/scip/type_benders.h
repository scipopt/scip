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
    LP      = 1,     /**< the Benders' subproblems are solved during the enforcement of an LP solution */
    RELAX   = 2,     /**< the Benders' subproblems are solved during the enforcement of a relaxation solution */
    PSEUDO  = 3,     /**< the Benders' subproblems are solved during the enforcement of a pseudo solution */
    CHECK   = 4      /**< the Benders' subproblems are solved during the checking of a solution for feasibility */
};
typedef enum SCIP_BendersEnfoType SCIP_BENDERSENFOTYPE;  /**< indicates the callback in cons_benders and cons_benderslp that triggered the subproblem solve */

typedef struct SCIP_Benders SCIP_BENDERS;           /**< Benders' decomposition data */
typedef struct SCIP_BendersData SCIP_BENDERSDATA;   /**< locally defined Benders' decomposition data */


/** copy method for Benders' decomposition plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders          : the Benders' decomposition itself
 *  - valid           : was the copying process valid?
 */
#define SCIP_DECL_BENDERSCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_Bool* valid)

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** initialization method of Benders' decomposition (called after problem was transformed and benders is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** deinitialization method of Benders' decomposition (called before transformed problem is freed and benders is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** presolving initialization method of the Benders' decomposition (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSINITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** presolving deinitialization method of the Benders' decomposition (called after presolving has been finished)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSEXITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The Benders' decomposition may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The Benders' decomposition should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition itself
 */
#define SCIP_DECL_BENDERSEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must register the SCIP instance for the subproblem with the Benders' decomposition core by calling
 *  SCIPaddBendersSubproblem. Typically, the user must create the SCIP instances for the subproblems. These can be
 *  created within a reader or probdata and then registered with the Benders' decomposition core during the call of this
 *  callback. If there are any settings required for solving the subproblems, then they should be set here. However,
 *  some settings will be overridden by the standard solving method included in the Benders' decomposition framework.
 *  If a special solving method is desired, the user can implement the bendersSolvesubXyz callback.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - probnumber      : the subproblem problem number
 */
#define SCIP_DECL_BENDERSCREATESUB(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, int probnumber)

/** called before the subproblem solving loop for Benders' decomposition. The pre subproblem solve function gives the
 *  user an oppportunity to perform any global set up for the Benders' decomposition.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 */
#define SCIP_DECL_BENDERSPRESUBSOLVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders)

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
 */
#define SCIP_DECL_BENDERSSOLVESUB(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_SOL* sol, int probnumber,\
  SCIP_Bool* infeasible)

/** the post-solve method for Benders' decomposition. The post-solve method is called after the subproblems have
 * been solved but before they have been freed. After the solving of the Benders' decomposition subproblems, the
 * subproblem solving data is freed in the SCIP_DECL_BENDERSFREESUB callback. However, it is not necessary to implement
 * SCIP_DECL_BENDERSFREESUB.
 *
 * If SCIP_DECL_BENDERSFREESUB is not implemented, then the Benders' decomposition framework will perform a default
 * freeing of the subproblems. If a subproblem is an LP, then they will be in probing mode for the subproblem
 * solve. So the freeing process involves ending the probing mode. If the subproblem is a MIP, then the subproblem is
 * solved by calling SCIPsolve. As such, the transformed problem must be freed after each subproblem solve.
 *
 * This callback provides the opportunity for the user to clean up any data structures that should not exist beyond the current
 * iteration.
 * The user has full access to the master and subproblems in this callback. So it is possible to construct solution for
 * the master problem in the method.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - sol             : the solution that was checked by solving the subproblems. Can be NULL.
 *  - infeasible      : indicates whether at least one subproblem is infeasible
 */
#define SCIP_DECL_BENDERSPOSTSOLVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_SOL* sol, SCIP_Bool infeasible)

/** frees the subproblem so that it can be resolved in the next iteration. As stated above, it is not necessary to
 *  implement this callback. If the callback is implemented, the subproblems should be freed by calling
 *  SCIPfreeTransform(). However, if the subproblems are LPs, then it could be more efficient to put the subproblem
 *  into probing mode prior to solving and then exiting the probing mode during the callback. To put the subproblem into
 *  probing mode, the subproblem must be in SCIP_STAGE_SOLVING. This can be achieved by using eventhandlers.
 *
 *  If SCIP_DECL_BENDERSFREESUB is not implemented, then the Benders' decomposition framework will perform a default
 *  freeing of the subproblems. If a subproblem is an LP, then they will be in probing mode for the subproblem
 *  solve. So the freeing process involves ending the probing mode. If the subproblem is a MIP, then the subproblem is
 *  solved by calling SCIPsolve. As such, the transformed problem must be freed after each subproblem solve.
 *
 *  NOTE: The freeing methods must be thread safe.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition data structure
 *  - probnumber      : the subproblem problem number
 */
#define SCIP_DECL_BENDERSFREESUB(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, int probnumber)

/** the variable mapping from the subproblem to the master problem. It is neccessary to have a mapping between every
 *  master problem variable and its counterpart in the subproblem. This mapping must go both ways: from master to sub
 *  and sub to master.
 *
 *  This method is called when generating the cuts. The cuts are generated by using the solution to the subproblem to
 *  eliminate a solution to the master problem.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition structure
 *  - var             : the variable for which the corresponding variable in the master or subproblem is required
 *  - mappedvar       : pointer to store the variable that is mapped to var
 *  - probnumber      : the number of the subproblem that the desired variable belongs to, -1 for the master problem
 */
#define SCIP_DECL_BENDERSGETVAR(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_VAR* var,\
   SCIP_VAR** mappedvar, int probnumber)

#ifdef __cplusplus
}
#endif

#endif
