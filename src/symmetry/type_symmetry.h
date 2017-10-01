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
   SYM_SENSE_XOR        = 3                  /**< XOR constraint */
};
typedef enum SYM_Rhssense SYM_RHSSENSE;

/** data of variables that are considered to be equivalent */
struct SYM_Vartype
{
   SCIP_Real             obj;                /**< objective of variable */
   SCIP_Real             lb;                 /**< lower bound of variable */
   SCIP_Real             ub;                 /**< upper bound of variable */
   SCIP_VARTYPE          type;               /**< type of variable */
   int                   color;              /**< store color */
};
typedef struct SYM_Vartype SYM_VARTYPE;

/** data of rhs that are considered to be equivalent */
struct SYM_Rhstype
{
   SCIP_Real             val;                /**< value of rhs */
   SYM_RHSSENSE          sense;              /**< sense of rhs */
   int                   color;              /**< store color */
};
typedef struct SYM_Rhstype SYM_RHSTYPE;

/** data for symmetry group computation */
struct SYM_Matrixdata
{
   SCIP_Real*            matcoef;            /**< nonzero coefficients appearing in the matrix */
   SCIP_Real*            rhscoef;            /**< rhs coefficients */
   SYM_RHSSENSE*         rhssense;           /**< sense of rhs */
   int*                  matrhsidx;          /**< indices of rhs corresponding to matrix entries */
   int*                  matvaridx;          /**< indices of variables for matrix entries */
   int*                  matidx;             /**< indices in mat(rhs/var)idx array corresponding to matrix coefficients */
   int*                  rhsidx;             /**< indices in rhstype array corresponding to rhs coefficients */
   int                   nmatcoef;           /**< number of coefficients in matrix */
   int                   nrhscoef;           /**< number of coefficients in rhs */
   int                   nmaxmatcoef;        /**< maximal number of matrix coefficients (will be increase on demand) */
   SCIP_Real*            matcoeforig;        /**< original coefficients in matrix */
   SCIP_Real*            rhscoeforig;        /**< original coefficients in rhs */
   int                   nuniquevars;        /**< number of unique variable types */
   int                   nuniquerhs;         /**< number of unique rhs types */
   int                   nuniquemat;         /**< number of unique matrix coefficients */
   SCIP_HASHTABLE*       rhstypemap;         /**< hash table for colors attached to non-equivalent rhs */
   int                   npermvars;          /**< number of variables for permutations */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int*                  permvarcolors;      /**< array for storing the colors of the individual variables */
   int*                  matcoefcolors;      /**< array for storing the colors of all matrix coefficients */
};
typedef struct SYM_Matrixdata SYM_MATRIXDATA;

#ifdef __cplusplus
}
#endif

#endif
