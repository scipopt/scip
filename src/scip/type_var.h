/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_var.h,v 1.4 2004/02/04 17:27:51 bzfpfend Exp $"

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
   SCIP_BOUNDCHGTYPE_INFERENCE = 1      /**< bound change was due to an inference (e.g. domain propagation) */
};
typedef enum BoundchgType BOUNDCHGTYPE;

typedef struct DomChgBound DOMCHGBOUND; /**< static domain change for bound changes */
typedef struct DomChgBoth DOMCHGBOTH;   /**< static domain change for bound and hole changes */
typedef struct DomChgDyn DOMCHGDYN;     /**< dynamic domain change for bound and hole changes */
typedef union DomChg DOMCHG;            /**< changes in domains of variables */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct BranchingData BRANCHINGDATA; /**< data for branching decision bound changes */
typedef struct InferenceData INFERENCEDATA; /**< data for inferred bound changes */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct Dom DOM;                 /**< datastructures for storing domains of variables */
typedef struct Aggregate AGGREGATE;     /**< aggregation information */
typedef struct Multaggr MULTAGGR;       /**< multiple aggregation information */
typedef struct Negate NEGATE;           /**< negation information */
typedef struct Var VAR;                 /**< variable of the problem */


#endif
