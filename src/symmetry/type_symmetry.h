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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SYMMETRY_H_
#define __SCIP_TYPE_SYMMETRY_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** symmetry type specification */
enum SYM_Spec
{
   SYM_SPEC_INTEGER = 1,                     /**< need symmetries for integer variables only */
   SYM_SPEC_BINARY  = 2,                     /**< need symmetries for binary variables only */
   SYM_SPEC_REAL    = 4                      /**< need symmetries also for continuous variables */
};
typedef enum SYM_Spec SYM_SPEC;

/** define sense of rhs */
enum SYM_Rhssense
{
   SYM_SENSE_UNKOWN     = 0,                 /**< unknown sense */
   SYM_SENSE_INEQUALITY = 1,                 /**< linear inequality */
   SYM_SENSE_EQUATION   = 2,                 /**< linear equation */
   SYM_SENSE_XOR        = 3,                 /**< XOR constraint */
   SYM_SENSE_AND        = 4,                 /**< AND constraint */
   SYM_SENSE_OR         = 5                  /**< OR constrant */
};
typedef enum SYM_Rhssense SYM_RHSSENSE;

/* type of symmetry handling codes */
#define SYM_HANDLETYPE_NONE             UINT64_C(0x00000000)  /**< no symmetry handling */
#define SYM_HANDLETYPE_SYMBREAK         UINT64_C(0x00000001)  /**< symmetry breaking inequalities */
#define SYM_HANDLETYPE_ORBITALFIXING    UINT64_C(0x00000002)  /**< orbital fixing */

typedef unsigned int SYM_HANDLETYPE;         /**< type of symmetry handling */

typedef struct SYM_Vartype SYM_VARTYPE;      /**< data of variables that are considered to be equivalent */
typedef struct SYM_Matrixdata SYM_MATRIXDATA;/**< data for symmetry group computation */

#ifdef __cplusplus
}
#endif

#endif
