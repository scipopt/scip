/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_presol.h,v 1.8 2005/02/07 18:12:02 bzfpfend Exp $"

/**@file   type_presol.h
 * @brief  type definitions for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_PRESOL_H__
#define __TYPE_PRESOL_H__


typedef struct Presol PRESOL;           /**< presolver data structure */
typedef struct PresolData PRESOLDATA;   /**< presolver specific data */


/** destructor of presolver to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define DECL_PRESOLFREE(x) RETCODE x (SCIP* scip, PRESOL* presol)

/** initialization method of presolver (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define DECL_PRESOLINIT(x) RETCODE x (SCIP* scip, PRESOL* presol)

/** deinitialization method of presolver (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define DECL_PRESOLEXIT(x) RETCODE x (SCIP* scip, PRESOL* presol)

/** presolving initialization method of presolver (called when presolving is about to begin)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
 */
#define DECL_PRESOLINITPRE(x) RETCODE x (SCIP* scip, PRESOL* presol, RESULT* result)

/** presolving deinitialization method of presolver (called after presolving has been finished)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
 */
#define DECL_PRESOLEXITPRE(x) RETCODE x (SCIP* scip, PRESOL* presol, RESULT* result)

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
 *  - nnewchgbds      : number of variable bounds tightend since the last call to the presolver
 *  - nnewholes       : number of domain holes added since the last call to the presolver
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolver
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolver
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolver
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolver
 *
 *  input/output:
 *  - nfixedvars      : pointer to total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to total number of variable bounds tightend of all presolvers
 *  - naddholes       : pointer to total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to total number of deleted constraints of all presolvers
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
#define DECL_PRESOLEXEC(x) RETCODE x (SCIP* scip, PRESOL* presol, int nrounds,              \
   int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
   int nnewdelconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides,                 \
   int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes,        \
   int* ndelconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, RESULT* result)



#include "def.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_scip.h"


#endif
