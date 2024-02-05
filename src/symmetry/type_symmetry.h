/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/** define symmetry types detectable by SCIP */
enum SYM_Symtype
{
   SYM_SYMTYPE_PERM      = 0,                /**< permutation symmetries */
   SYM_SYMTYPE_SIGNPERM  = 1                 /**< signed permutation symmetries */
};
typedef enum SYM_Symtype SYM_SYMTYPE;

/** define type of nodes in symmetry detection expression trees */
enum SYM_Nodetype
{
   SYM_NODETYPE_OPERATOR = 0,                /**< operator node */
   SYM_NODETYPE_VAL      = 1,                /**< numerical value node */
   SYM_NODETYPE_CONS     = 2,                /**< constraint node */
   SYM_NODETYPE_VAR      = 3                 /**< variable node */
};
typedef enum SYM_Nodetype SYM_NODETYPE;

/** define type of simple constraints/operators in symmetry detection */
enum SYM_Consoptype
{
   SYM_CONSOPTYPE_UNKNOWN     = 0,           /**< unknown constraint type */
   SYM_CONSOPTYPE_BDDISJ      = 1,           /**< constraint of type bounddisjunction */
   SYM_CONSOPTYPE_EQ          = 2,           /**< encodes == in indicator constraints for activation variable */
   SYM_CONSOPTYPE_SOS2_TUPLE  = 3,           /**< encodes pairs in SOS2 constraints */
   SYM_CONSOPTYPE_SUM         = 4,           /**< indicates sums if sum-expr undefined */
   SYM_CONSOPTYPE_SLACK       = 5,           /**< indicates slack vars in indicator constraints */
   SYM_CONSOPTYPE_COEF        = 6,           /**< indicates coefficients from parent expressions */
   SYM_CONSOPTYPE_SQDIFF      = 7,           /**< indicates a squared difference */
   SYM_CONSOPTYPE_CARD_TUPLE  = 8,           /**< encodes pairs in cardinality constraints */
   SYM_CONSOPTYPE_LAST        = 9            /**< number of predefined enum types, needs to always
                                              *   hold the biggest value */
};
typedef enum SYM_Consoptype SYM_CONSOPTYPE;

/* type of symmetry handling codes */
#define SYM_HANDLETYPE_NONE             UINT32_C(0x00000000)  /**< no symmetry handling */
#define SYM_HANDLETYPE_SYMBREAK         UINT32_C(0x00000001)  /**< symmetry breaking inequalities */
#define SYM_HANDLETYPE_ORBITALREDUCTION UINT32_C(0x00000002)  /**< orbital reduction */
#define SYM_HANDLETYPE_SST              UINT32_C(0x00000004)  /**< Schreier Sims cuts */
#define SYM_HANDLETYPE_SYMCONS (SYM_HANDLETYPE_SYMBREAK | SYM_HANDLETYPE_SST)

typedef uint32_t SYM_HANDLETYPE;        /**< type of symmetry handling */

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

#ifdef __cplusplus
}
#endif

#endif
