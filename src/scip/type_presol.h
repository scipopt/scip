/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_presol.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PRESOL_H__
#define __SCIP_TYPE_PRESOL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Presol SCIP_PRESOL;           /**< presolver data structure */
typedef struct SCIP_PresolData SCIP_PRESOLDATA;   /**< presolver specific data */


/** copy method for presolver plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** destructor of presolver to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** initialization method of presolver (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** deinitialization method of presolver (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** presolving initialization method of presolver (called when presolving is about to begin)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 *  - isunbounded     : was the problem already declared to be unbounded
 *  - isinfeasible    : was the problem already declared to be infeasible
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_FEASIBLE   : no infeasibility or unboundedness could be found
 */
#define SCIP_DECL_PRESOLINITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol, SCIP_Bool isunbounded, \
      SCIP_Bool isinfeasible, SCIP_RESULT* result)

/** presolving deinitialization method of presolver (called after presolving has been finished)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 *  - isunbounded     : was the problem already declared to be unbounded
 *  - isinfeasible    : was the problem already declared to be infeasible
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_FEASIBLE   : no infeasibility or unboundedness could be found
 */
#define SCIP_DECL_PRESOLEXITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol, SCIP_Bool isunbounded, \
      SCIP_Bool isinfeasible, SCIP_RESULT* result)

/** execution method of presolver
 *
 *  The presolver should go through the variables and constraints and tighten the domains or
 *  constraints. Each tightening should increase the given total numbers of changes.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 *  - nrounds         : number of presolving rounds already done
 *  - nnewfixedvars   : number of variables fixed since the last call to the presolver
 *  - nnewaggrvars    : number of variables aggregated since the last call to the presolver
 *  - nnewchgvartypes : number of variable type changes since the last call to the presolver
 *  - nnewchgbds      : number of variable bounds tightened since the last call to the presolver
 *  - nnewholes       : number of domain holes added since the last call to the presolver
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolver
 *  - nnewaddconss    : number of added constraints since the last call to the presolver
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolver
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolver
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolver
 *
 *  @note the counters state the changes since the last call including the changes of this presolver during its last
 *        last call
 *
 *  input/output:
 *  - nfixedvars      : pointer to total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to total number of variable bounds tightened of all presolvers
 *  - naddholes       : pointer to total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to total number of deleted constraints of all presolvers
 *  - naddconss       : pointer to total number of added constraints of all presolvers
 *  - nupgdconss      : pointer to total number of upgraded constraints of all presolvers
 *  - nchgcoefs       : pointer to total number of changed coefficients of all presolvers
 *  - nchgsides       : pointer to total number of changed left/right hand sides of all presolvers
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_SUCCESS    : the presolver found a reduction
 *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolver was skipped
 *  - SCIP_DELAYED    : the presolver was skipped, but should be called again
 */
#define SCIP_DECL_PRESOLEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol, int nrounds, \
      int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
      int nnewdelconss, int nnewaddconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides, \
      int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes, \
      int* ndelconss, int* naddconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
