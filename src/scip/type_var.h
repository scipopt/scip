/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_var.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for problem variables
 * @author Tobias Achterberg
 *
 *  This file defines the interface for user variable data implemented in C. Each variable can be equipped with a
 *  variable data struct. This data can be accessed via the function SCIPgetVardata() at any time after it is created
 *  and before it is deleted.
 *
 *  - \ref scip::ObjVardata "Corresponding C interface"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_VAR_H__
#define __SCIP_TYPE_VAR_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** status of problem variables */
enum SCIP_Varstatus
{
   SCIP_VARSTATUS_ORIGINAL   = 0,       /**< variable belongs to original problem */
   SCIP_VARSTATUS_LOOSE      = 1,       /**< variable is a loose variable of the transformed problem */
   SCIP_VARSTATUS_COLUMN     = 2,       /**< variable is a column of the transformed problem */
   SCIP_VARSTATUS_FIXED      = 3,       /**< variable is fixed to specific value in the transformed problem */
   SCIP_VARSTATUS_AGGREGATED = 4,       /**< variable is aggregated to x = a*y + c in the transformed problem */
   SCIP_VARSTATUS_MULTAGGR   = 5,       /**< variable is aggregated to x = a_1*y_1 + ... + a_k*y_k + c */
   SCIP_VARSTATUS_NEGATED    = 6        /**< variable is the negation of an original or transformed variable */
};
typedef enum SCIP_Varstatus SCIP_VARSTATUS;

/** variable type */
enum SCIP_Vartype
{
   SCIP_VARTYPE_BINARY     = 0,              /**< binary variable: \f$ x \in \{0,1\} \f$ */
   SCIP_VARTYPE_INTEGER    = 1,              /**< integer variable: \f$ x \in \{lb, \dots, ub\} \f$ */
#if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 6) /* _attribute__ ((deprecated)) within enums only supported by selected compilers */
   SCIP_VARTYPE_IMPLINT SCIP_DEPRECATED = 2, /**< @deprecated use `SCIPcreateVarImpl()` or `SCIPvarChgImplType()` to control implied integrality */
#else
   SCIP_VARTYPE_IMPLINT    = 2,              /**< @deprecated use `SCIPcreateVarImpl()` or `SCIPvarChgImplType()` to control implied integrality */
#endif
   SCIP_VARTYPE_CONTINUOUS = 3               /**< continuous variable: \f$ lb \leq x \leq ub \f$ */
};
typedef enum SCIP_Vartype SCIP_VARTYPE;

/** alternative to `SCIP_VARTYPE_IMPLINT` that comes without SCIP_DEPRECATED attribute
 * @deprecated still, the use of an implicit integral variable type is deprecated
 */
/** @todo remove / refactor usages of this macro */
#define SCIP_DEPRECATED_VARTYPE_IMPLINT ((SCIP_VARTYPE) 2)

/* CIP format variable characters */
#define SCIP_VARTYPE_BINARY_CHAR 'B'
#define SCIP_VARTYPE_INTEGER_CHAR 'I'
#define SCIP_VARTYPE_CONTINUOUS_CHAR 'C'
#define SCIP_DEPRECATED_VARTYPE_IMPLINT_CHAR 'M'

/** implied integral type */
enum SCIP_ImplintType
{
  SCIP_IMPLINTTYPE_NONE   = 0,               /**< The variable is not implied integral by other variables */
  SCIP_IMPLINTTYPE_WEAK   = 1,               /**< The constraint handlers enforce that if the problem is relaxed
                                              *   to have integrality constraints for the non-implied integral variables
                                              *   only, there exists an optimal solution where all weakly and strongly
                                              *   implied integral variables have integer solution values.
                                              *   For infeasible problems, when relaxing integrality of all implied
                                              *   integer variables, the problem remains infeasible.
                                              *   For unbounded problems, when enforcing integrality of all implied
                                              *   integer variables, the problem remains unbounded.
                                              *
                                              *   @note This notion of implied integrality is fragile and may break
                                              *   if extra constraints are added.
                                              *
                                              *   Example: The variable z is a weakly implied integral if it only occurs
                                              *   in the constraint 4x + 3y + z <= 10, where x and y are integer and
                                              *   z has objective 0. */
  SCIP_IMPLINTTYPE_STRONG = 2                /**< The constraint handlers enforce that if the problem is relaxed
                                              *   to have integrality constraints for the non-implied integral variables
                                              *   only, in every feasible solution all strongly implied integral
                                              *   variables have integer solution values.
                                              *
                                              *   @note This notion of implied integrality remains intact under the
                                              *   addition of additional constraints to the problem.
                                              *
                                              *   Example: The variable z is strongly implied integral if we have the
                                              *   constraint: 4x + 3y + z = 10, where x and y are integer variables. */
};
typedef enum SCIP_ImplintType SCIP_IMPLINTTYPE;

/** domain change data type */
enum SCIP_DomchgType
{
   SCIP_DOMCHGTYPE_DYNAMIC = 0,         /**< dynamic bound changes with size information of arrays */
   SCIP_DOMCHGTYPE_BOTH    = 1,         /**< static domain changes: number of entries equals size of arrays */
   SCIP_DOMCHGTYPE_BOUND   = 2          /**< static domain changes without any hole changes */
};
typedef enum SCIP_DomchgType SCIP_DOMCHGTYPE;

/** bound change type */
enum SCIP_BoundchgType
{
   SCIP_BOUNDCHGTYPE_BRANCHING = 0,     /**< bound change was due to a branching decision */
   SCIP_BOUNDCHGTYPE_CONSINFER = 1,     /**< bound change was due to an inference of a constraint (domain propagation) */
   SCIP_BOUNDCHGTYPE_PROPINFER = 2      /**< bound change was due to an inference of a domain propagator */
};
typedef enum SCIP_BoundchgType SCIP_BOUNDCHGTYPE;

/** types of variable locks */
#define NLOCKTYPES 2                    /**< number of lock types */
enum SCIP_LockType
{
   SCIP_LOCKTYPE_MODEL    = 0,          /**< variable locks for model and check constraints */
   SCIP_LOCKTYPE_CONFLICT = 1           /**< variable locks for conflict constraints */
};
typedef enum SCIP_LockType SCIP_LOCKTYPE;

typedef struct SCIP_DomChgBound SCIP_DOMCHGBOUND; /**< static domain change data for bound changes */
typedef struct SCIP_DomChgBoth SCIP_DOMCHGBOTH;   /**< static domain change data for bound and hole changes */
typedef struct SCIP_DomChgDyn SCIP_DOMCHGDYN;     /**< dynamic domain change data for bound and hole changes */
typedef union SCIP_DomChg SCIP_DOMCHG;            /**< changes in domains of variables */
typedef struct SCIP_BoundChg SCIP_BOUNDCHG;       /**< changes in bounds of variables */
typedef struct SCIP_BdChgIdx SCIP_BDCHGIDX;       /**< bound change index in path from root to current node */
typedef struct SCIP_BdChgInfo SCIP_BDCHGINFO;     /**< bound change information to track bound changes from root to current node */
typedef struct SCIP_BranchingData SCIP_BRANCHINGDATA; /**< data for branching decision bound changes */
typedef struct SCIP_InferenceData SCIP_INFERENCEDATA; /**< data for inferred bound changes */
typedef struct SCIP_HoleChg SCIP_HOLECHG;         /**< changes in holelist of variables */
typedef struct SCIP_Hole SCIP_HOLE;               /**< hole in a domain of an integer variable */
typedef struct SCIP_Holelist SCIP_HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct SCIP_Dom SCIP_DOM;                 /**< datastructures for storing domains of variables */
typedef struct SCIP_Original SCIP_ORIGINAL;       /**< original variable information */
typedef struct SCIP_Loose SCIP_LOOSE;             /**< loose variable information */
typedef struct SCIP_Aggregate SCIP_AGGREGATE;     /**< aggregation information */
typedef struct SCIP_AggregateExact SCIP_AGGREGATEEXACT; /**< exact aggregation information */
typedef struct SCIP_Multaggr SCIP_MULTAGGR;       /**< multiple aggregation information */
typedef struct SCIP_MultaggrExact SCIP_MULTAGGREXACT; /**< exact multiple aggregation information */
typedef struct SCIP_Negate SCIP_NEGATE;           /**< negation information */
typedef struct SCIP_Var SCIP_VAR;                 /**< variable of the problem */
typedef struct SCIP_VarData SCIP_VARDATA;         /**< user variable data */
typedef struct SCIP_VarDataExact SCIP_VARDATAEXACT;/**< exact data (obj, bounds, ...) */
typedef struct SCIP_DomExact SCIP_DOMEXACT;       /**< exact domain (using rationals) */

/** frees user data of original variable (called when the original variable is freed)
 *
 *  This method should free the user data of the original variable.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - var             : original variable the data to free is belonging to
 *  - vardata         : pointer to the user variable data to free
 */
#define SCIP_DECL_VARDELORIG(x) SCIP_RETCODE x (SCIP* scip, SCIP_VAR* var, SCIP_VARDATA** vardata)

/** creates transformed variable for original user variable
 *
 *  Because the original variable and the user data of the original variable should not be
 *  modified during the solving process, a transformed variable is created as a copy of
 *  the original variable. If the user variable data is never modified during the solving
 *  process anyways, it is enough to simply copy the user data's pointer. This is the
 *  default implementation, which is used when NULL is given as the VARTRANS method.
 *  If the user data may be modified during the solving process (e.g. during preprocessing),
 *  the VARTRANS method must be given and has to copy the user variable data to a different
 *  memory location.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sourcevar       : original variable
 *  - sourcedata      : source variable data to transform
 *  - targetvar       : transformed variable
 *  - targetdata      : pointer to store created transformed variable data
 */
#define SCIP_DECL_VARTRANS(x) SCIP_RETCODE x (SCIP* scip, SCIP_VAR* sourcevar, SCIP_VARDATA* sourcedata, SCIP_VAR* targetvar, SCIP_VARDATA** targetdata)

/** frees user data of transformed variable (called when the transformed variable is freed)
 *
 *  This method has to be implemented, if the VARTRANS method is not a simple pointer
 *  copy operation like in the default VARTRANS implementation. It should free the
 *  user data of the transformed variable, that was created in the VARTRANS method.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - var             : transformed variable the data to free is belonging to
 *  - vardata         : pointer to the user variable data to free
 */
#define SCIP_DECL_VARDELTRANS(x) SCIP_RETCODE x (SCIP* scip, SCIP_VAR* var, SCIP_VARDATA** vardata)

/** copies variable data of source SCIP variable for the target SCIP variable
 *
 *  This method should copy the variable data of the source SCIP and create a target variable data for target
 *  variable. This callback is optional. If the copying process was successful, the target variable gets this variable
 *  data assigned. In case the result pointer is set to SCIP_DIDNOTRUN, the target variable will have no variable data at
 *  all.
 *
 *  The variable map and the constraint map can be used via the function SCIPgetVarCopy() and SCIPgetConsCopy(),
 *  respectively, to get for certain variables and constraints of the source SCIP the counter parts in the target
 *  SCIP. You should be very careful in using these two methods since they could lead to infinite loop.
 *
 *  input:
 *  - scip            : target SCIP data structure
 *  - sourcescip      : source SCIP main data structure
 *  - sourcevar       : variable of the source SCIP
 *  - sourcedata      : variable data of the source variable which should get copied
 *  - varmap,         : a hashmap which stores the mapping of source variables to corresponding target variables
 *  - consmap,        : a hashmap which stores the mapping of source constraints to corresponding target constraints
 *  - targetvar       : variable of the (target) SCIP (targetvar is the copy of sourcevar)
 *  - targetdata      : pointer to store created copy of the variable data for the (target) SCIP
 *
 *  output:
 *  - result          : pointer to store the result of the call
 *
 *  possible return values for *result:
 *  - SCIP_DIDNOTRUN  : the copying process was not performed 
 *  - SCIP_SUCCESS    : the copying process was successfully performed
 */
#define SCIP_DECL_VARCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP* sourcescip, SCIP_VAR* sourcevar, SCIP_VARDATA* sourcedata, \
      SCIP_HASHMAP* varmap, SCIP_HASHMAP* consmap, SCIP_VAR* targetvar, SCIP_VARDATA** targetdata, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
