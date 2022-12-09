/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry computations
 * @author Marc Pfetsch
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SYMMETRY_H_
#define __SCIP_TYPE_SYMMETRY_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** symmetry type specification */
#define SYM_SPEC_INTEGER                UINT32_C(0x00000001)  /**< need symmetries for integer variables only */
#define SYM_SPEC_BINARY                 UINT32_C(0x00000002)  /**< need symmetries for binary variables only */
#define SYM_SPEC_REAL                   UINT32_C(0x00000004)  /**< need symmetries also for continuous variables */

typedef uint32_t SYM_SPEC;              /**< types of variables handled by symmetry */

/** symmetry timings */
#define SYM_COMPUTETIMING_BEFOREPRESOL    0  /**< compute symmetries before presolving */
#define SYM_COMPUTETIMING_DURINGPRESOL    1  /**< compute symmetries during presolving */
#define SYM_COMPUTETIMING_AFTERPRESOL     2  /**< compute symmetries after presolving */

/** define sense of rhs */
enum SYM_Rhssense
{
   SYM_SENSE_UNKOWN     = 0,                 /**< unknown sense */
   SYM_SENSE_INEQUALITY = 1,                 /**< linear inequality */
   SYM_SENSE_EQUATION   = 2,                 /**< linear equation */
   SYM_SENSE_XOR        = 3,                 /**< XOR constraint */
   SYM_SENSE_AND        = 4,                 /**< AND constraint */
   SYM_SENSE_OR         = 5,                 /**< OR constrant */
   SYM_SENSE_BOUNDIS_TYPE_1 = 6,             /**< bounddisjunction type 1 */
   SYM_SENSE_BOUNDIS_TYPE_2 = 7              /**< bounddisjunction type 2 */
};
typedef enum SYM_Rhssense SYM_RHSSENSE;

/** define type of nodes in symmetry detection expression trees */
enum SYM_Nodetype
{
   SYM_NODETYPE_OPERATOR = 0,                /**< operator node */
   SYM_NODETYPE_VAR      = 1,                /**< variable node */
   SYM_NODETYPE_COEF     = 2,                /**< coefficient node */
   SYM_NODETYPE_VAL      = 3,                /**< numerical value node */
   SYM_NODETYPE_RHS      = 4                 /**< rhs node */
};
typedef enum SYM_Nodetype SYM_NODETYPE;

/** define the type of constraints used in symmetry detection */
enum SYM_Constype
{
   SYM_CONSTYPE_LINEAR   = 0,                /**< linear constraint */
   SYM_CONSTYPE_SIMPLE   = 1,                /**< simple constraint */
   SYM_CONSTYPE_EXPR     = 2,                /**< constraints given by expression tree */
   SYM_CONSTYPE_OBJ      = 3                 /**< objective */
};
typedef enum SYM_Constype SYM_CONSTYPE;

/** define type of simple constraints/operators in symmetry detection */
enum SYM_Consoptype
{
   SYM_CONSOPTYPE_UNKNOWN  = 0,              /**< unkown constraint type */
   SYM_CONSOPTYPE_AND      = 1,              /**< constraint of type and */
   SYM_CONSOPTYPE_BDDISJ   = 2,              /**< constraint of type bounddisjunction */
   SYM_CONSOPTYPE_CARD     = 3,              /**< constraint of type cardinality */
   SYM_CONSOPTYPE_INDICATOR = 4,             /**< constraint of type indicator */
   SYM_CONSOPTYPE_OR       = 5,              /**< constraint of type or */
   SYM_CONSOPTYPE_PSEUDOBOOL = 6,            /**< constraint of type pseudoboolean */
   SYM_CONSOPTYPE_SOS1     = 7,              /**< constraint of type SOS1 */
   SYM_CONSOPTYPE_SOS2     = 8,              /**< constraint of type SOS2 */
   SYM_CONSOPTYPE_XOR      = 9,              /**< constraint of type xor */
   SYM_CONSOPTYPE_GEQ      = 10,             /**< needed to encode >= in bounddisjunctions */
   SYM_CONSOPTYPE_EQ       = 11,             /**< needed to encode == in indicator constraints */
   SYM_CONSOPTYPE_TUPLE    = 12,             /**< needed to encode pairs in SOS2 constraints */
   SYM_CONSOPTYPE_OBJ      = 13,             /**< needed to model the objective */
   SYM_CONSOPTYPE_NONLINEAR = 14,            /**< constraint of type nonlinear */
   SYM_CONSOPTYPE_POWER    = 15,             /**< needed to distinguish power and signpower */
   SYM_CONSOPTYPE_SIGNPOWER = 16,            /**< needed to distinguish power and signpower */
   SYM_CONSOPTYPE_BIPROD   = 17,             /**< needed to indicate product of two variables */
   SYM_CONSOPTYPE_SUM      = 18              /**< needed to indicate sums if sum-expr undefined */
};
typedef enum SYM_Consoptype SYM_CONSOPTYPE;

/** define type of flips that are allowed for children of operators */
enum SYM_Fliptype
{
   SYM_FLIPTYPE_NONE     = 0,                /**< no flips are permitted */
   SYM_FLIPTYPE_ODD      = 1,                /**< flips for odd functions are permitted */
   SYM_FLIPTYPE_EVEN     = 2,                /**< flips for even functions are permitted */
   SYM_FLIPTYPE_SHIFT_ODD = 3,               /**< flips for shift-odd functions are permitted */
   SYM_FLIPTYPE_BIPROD   = 4                 /**< flips for products of two variables */
};
typedef enum SYM_Fliptype SYM_FLIPTYPE;

/* type of symmetry handling codes */
#define SYM_HANDLETYPE_NONE             UINT32_C(0x00000000)  /**< no symmetry handling */
#define SYM_HANDLETYPE_SYMBREAK         UINT32_C(0x00000001)  /**< symmetry breaking inequalities */
#define SYM_HANDLETYPE_ORBITALFIXING    UINT32_C(0x00000002)  /**< orbital fixing */
#define SYM_HANDLETYPE_SST              UINT32_C(0x00000004)  /**< Schreier Sims cuts */
#define SYM_HANDLETYPE_SYMCONS (SYM_HANDLETYPE_SYMBREAK | SYM_HANDLETYPE_SST)

typedef uint32_t SYM_HANDLETYPE;        /**< type of symmetry handling */

typedef struct SYM_Vartype SYM_VARTYPE;      /**< data of variables that are considered to be equivalent */
typedef struct SYM_Optype SYM_OPTYPE;        /**< data of operators that are considered to be equivalent */
typedef struct SYM_Consttype SYM_CONSTTYPE;  /**< data of constants that are considered to be equivalent */
typedef struct SYM_Rhstype SYM_RHSTYPE;      /**< data of constraint sides that are considered to be equivalent */
typedef struct SYM_Matrixdata SYM_MATRIXDATA;/**< data for symmetry group computation on linear constraints */
typedef struct SYM_Reflsymdata SYM_REFLSYMDATA ;/**< data for reflection symmetry group computation */
typedef struct SYM_Exprdata SYM_EXPRDATA;    /**< data for symmetry group computation on nonlinear constraints */
typedef struct SYM_Consinfo SYM_CONSINFO;    /**< information about a constraint used in symmetry computation */
typedef struct SYM_Node SYM_NODE;            /**< data to encode a node of a symmetry detection graph */
typedef struct SYM_Edge SYM_EDGE;            /**< data to encode an edge of a symmetry detection graph */

/** selection rules for leaders in SST cuts */
enum SCIP_LeaderRule
{
   SCIP_LEADERRULE_FIRSTINORBIT        = 0,       /**< first var in orbit */
   SCIP_LEADERRULE_LASTINORBIT         = 1,       /**< last var in orbit */
   SCIP_LEADERRULE_MAXCONFLICTSINORBIT = 2        /**< var with most conflicting vars in its orbit */
};
typedef enum SCIP_LeaderRule SCIP_LEADERRULE;

/** tie breaks for leader rule based on the leader's orbit */
enum SCIP_LeaderTiebreakRule
{
   SCIP_LEADERTIEBREAKRULE_MINORBIT            = 0,    /**< orbit of minimum size */
   SCIP_LEADERTIEBREAKRULE_MAXORBIT            = 1,    /**< orbit of maximum size */
   SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT = 2     /**< orbit with maximum number of vars in conflict with leader */
};

/** variable types for leader in Schreier Sims cuts */
enum SCIP_SSTType
{
   SCIP_SSTTYPE_BINARY                 = 1,    /**< binary variables */
   SCIP_SSTTYPE_INTEGER                = 2,    /**< integer variables */
   SCIP_SSTTYPE_IMPLINT                = 4,    /**< implicitly integer variables */
   SCIP_SSTTYPE_CONTINUOUS             = 8     /**< continuous variables */
};

typedef enum SCIP_SSTType SCIP_SSTTYPE;

/** type of orbitope constraint: full, packing, or partitioning orbitope */
enum SCIP_OrbitopeType
{
   SCIP_ORBITOPETYPE_FULL         = 0,       /**< constraint is a full orbitope constraint:         rowsum(x) unrestricted */
   SCIP_ORBITOPETYPE_PARTITIONING = 1,       /**< constraint is a partitioning orbitope constraint: rowsum(x) == 1 */
   SCIP_ORBITOPETYPE_PACKING      = 2        /**< constraint is a packing orbitope constraint:      rowsum(x) <= 1 */
};
typedef enum SCIP_OrbitopeType SCIP_ORBITOPETYPE;

/** conditions to recompute symmetries after a restart */
enum SCIP_RecomputesymType
{
   SCIP_RECOMPUTESYM_NEVER         = 0,       /**< never recompute symmetries */
   SCIP_RECOMPUTESYM_ALWAYS        = 1,       /**< always recompute symmetries */
   SCIP_RECOMPUTESYM_OFFOUNDRED    = 2        /**< only if orbital fixing found a reduction in previous run */
};
typedef enum SCIP_RecomputesymType SCIP_RECOMPUTESYMTYPE;

#ifdef __cplusplus
}
#endif

#endif
