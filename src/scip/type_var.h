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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_var.h,v 1.14 2005/05/31 17:20:25 bzfpfend Exp $"

/**@file   type_var.h
 * @brief  type definitions for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_VAR_H__
#define __TYPE_VAR_H__



/** status of problem variables */
enum Varstatus
{
   SCIP_VARSTATUS_ORIGINAL   = 0,       /**< variable belongs to original problem */
   SCIP_VARSTATUS_LOOSE      = 1,       /**< variable is a loose variable of the transformed problem */
   SCIP_VARSTATUS_COLUMN     = 2,       /**< variable is a column of the transformed problem */
   SCIP_VARSTATUS_FIXED      = 3,       /**< variable is fixed to specific value in the transformed problem */
   SCIP_VARSTATUS_AGGREGATED = 4,       /**< variable is aggregated to x = a*y + c in the transformed problem */
   SCIP_VARSTATUS_MULTAGGR   = 5,       /**< variable is aggregated to x = a_1*y_1 + ... + a_k*y_k + c */
   SCIP_VARSTATUS_NEGATED    = 6        /**< variable is the negation of an original or transformed variable */
};
typedef enum Varstatus VARSTATUS;

/** variable type */
enum Vartype
{
   SCIP_VARTYPE_BINARY     = 0,         /**< binary variable: x in {0,1} */
   SCIP_VARTYPE_INTEGER    = 1,         /**< integer variable: x in {lb, ..., ub} */
   SCIP_VARTYPE_IMPLINT    = 2,         /**< implicit integer variable: continuous variable, that is always integral */
   SCIP_VARTYPE_CONTINUOUS = 3          /**< continuous variable: x in [lb,ub] */
};
typedef enum Vartype VARTYPE;

/** domain change data type */
enum DomchgType
{
   SCIP_DOMCHGTYPE_DYNAMIC = 0,         /**< dynamic bound changes with size information of arrays */
   SCIP_DOMCHGTYPE_BOTH    = 1,         /**< static domain changes: number of entries equals size of arrays */
   SCIP_DOMCHGTYPE_BOUND   = 2          /**< static domain changes without any hole changes */
};
typedef enum DomchgType DOMCHGTYPE;

/** bound change type */
enum BoundchgType
{
   SCIP_BOUNDCHGTYPE_BRANCHING = 0,     /**< bound change was due to a branching decision */
   SCIP_BOUNDCHGTYPE_CONSINFER = 1,     /**< bound change was due to an inference of a constraint (domain propagation) */
   SCIP_BOUNDCHGTYPE_PROPINFER = 2      /**< bound change was due to an inference of a domain propagator */
};
typedef enum BoundchgType BOUNDCHGTYPE;

typedef struct DomChgBound DOMCHGBOUND; /**< static domain change data for bound changes */
typedef struct DomChgBoth DOMCHGBOTH;   /**< static domain change data for bound and hole changes */
typedef struct DomChgDyn DOMCHGDYN;     /**< dynamic domain change data for bound and hole changes */
typedef union DomChg DOMCHG;            /**< changes in domains of variables */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct BdChgIdx BDCHGIDX;       /**< bound change index in path from root to current node */
typedef struct BdChgInfo BDCHGINFO;     /**< bound change information to track bound changes from root to current node */
typedef struct BranchingData BRANCHINGDATA; /**< data for branching decision bound changes */
typedef struct InferenceData INFERENCEDATA; /**< data for inferred bound changes */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct Dom DOM;                 /**< datastructures for storing domains of variables */
typedef struct VBounds VBOUNDS;         /**< variable bounds of a variable x in the form x <= c*y or x >= c*y */
typedef struct Implics IMPLICS;         /**< implications in the form x <= 0 or x >= 1 ==> y <= b or y >= b for x binary, NULL if x nonbinary */
typedef struct Original ORIGINAL;       /**< original variable information */
typedef struct Aggregate AGGREGATE;     /**< aggregation information */
typedef struct Multaggr MULTAGGR;       /**< multiple aggregation information */
typedef struct Negate NEGATE;           /**< negation information */
typedef struct Var VAR;                 /**< variable of the problem */
typedef struct VarData VARDATA;         /**< user variable data */


/** frees user data of original variable (called when the original variable is freed)
 *
 *  This method should free the user data of the original variable.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    var             : original variable the data to free is belonging to
 *    vardata         : pointer to the user variable data to free
 */
#define DECL_VARDELORIG(x) RETCODE x (SCIP* scip, VAR* var, VARDATA** vardata)

/** creates transformed variable for original user variable
 *
 *  Because the original variable and the user data of the original variable should not be
 *  modified during the solving process, a transformed variable is created as a copy of
 *  the original variable. If the user variable data is never modified during the solving
 *  process anyways, it is enough to simple copy the user data's pointer. This is the
 *  default implementation, which is used when a NULL is given as VARTRANS method.
 *  If the user data may be modified during the solving process (e.g. during preprocessing),
 *  the VARTRANS method must be given and has to copy the user variable data to a different
 *  memory location.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sourcevar       : original variable
 *    sourcedata      : source variable data to transform
 *    targetvar       : transformed variable
 *    targetdata      : pointer to store created transformed variable data
 */
#define DECL_VARTRANS(x) RETCODE x (SCIP* scip, VAR* sourcevar, VARDATA* sourcedata, VAR* targetvar, VARDATA** targetdata)

/** frees user data of transformed variable (called when the transformed variable is freed)
 *
 *  This method has to be implemented, if the VARTRANS method is not a simple pointer
 *  copy operation like in the default VARTRANS implementation. It should free the
 *  user data of the transformed variable, that was created in the VARTRANS method.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    var             : transformed variable the data to free is belonging to
 *    vardata         : pointer to the user variable data to free
 */
#define DECL_VARDELTRANS(x) RETCODE x (SCIP* scip, VAR* var, VARDATA** vardata)

#endif
