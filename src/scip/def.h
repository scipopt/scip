/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   def.h
 * @brief  comon defines and data types used in all packages of SCIP
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DEF_H__
#define __DEF_H__


#include <limits.h>

#include "retcode.h"



#define SCIP_VERSION                 20 /**< SCIP version number (multiplied by 100 to get integer number) */



#ifndef NDEBUG
#define CHECK_OKAY(x) { RETCODE _restat_;                                                                       \
                        if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                        {                                                                                   \
                          printf("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_);    \
                          return _restat_;                                                                  \
                        }                                                                                   \
                      }
#define ALLOC_OKAY(x) { if( NULL == (x) )                                                                   \
                        {                                                                                   \
                          printf("[%s:%d] No memory in function call\n", __FILE__, __LINE__);               \
                          return SCIP_NOMEMORY;                                                             \
                        }                                                                                   \
                      }
#else
#define CHECK_OKAY(x) { int _restat_; if( (_restat_ = (x)) != SCIP_OKAY ) return _restat_; }
#define ALLOC_OKAY(x) { if( NULL == (x) ) return SCIP_NOMEMORY; }
#endif


#ifndef SQR
#define SQR(x)        ((x)*(x))
#define SQRT(x)       (sqrt(x))
#endif

#ifndef ABS
#define ABS(x)        ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MAX
#define MAX(x,y)      ((x) >= (y) ? (x) : (y))     /**< returns maximum of x and y */
#define MIN(x,y)      ((x) <= (y) ? (x) : (y))     /**< returns minimum of x and y */
#endif

#define EPSEQ(x,y,eps) (ABS((x)-(y)) <= (eps))
#define EPSLT(x,y,eps) ((x)-(y) < -(eps))
#define EPSLE(x,y,eps) ((x)-(y) <= (eps))
#define EPSGT(x,y,eps) ((x)-(y) > (eps))
#define EPSGE(x,y,eps) ((x)-(y) >= -(eps))
#define EPSZ(x,eps)    (ABS(x) <= (eps))
#define EPSP(x,eps)    ((x) > (eps))
#define EPSN(x,eps)    ((x) < -(eps))



/*
 * Boolean values
 */

typedef int Bool;                       /**< type used for boolean values */
#ifndef TRUE
#define TRUE (0==0)                     /**< boolean value TRUE */
#define FALSE (0==1)                    /**< boolean value FALSE */
#endif


/*
 * Long Integer values
 */

#ifndef LLONG_MAX
#define LLONG_MAX	9223372036854775807LL
#define LLONG_MIN	(-LLONG_MAX - 1LL)
#endif

typedef long long Longint;              /**< type used for long integer values */
#define LONGINT_MAX          LLONG_MAX
#define LONGINT_MIN          LLONG_MIN



/*
 * Floating point values
 */

typedef double Real;                    /**< type used for floating point values */
#define SCIP_DEFAULT_EPSILON     1e-09  /**< default upper bound for floating points to be considered zero */
#define SCIP_DEFAULT_SUMEPSILON  1e-07  /**< default upper bound for sums of floating points to be considered zero */
#define SCIP_DEFAULT_FEASTOL     1e-06  /**< default LP feasibility tolerance */
#define SCIP_DEFAULT_INFINITY  1.0E+20  /**< default value considered to be infinity */
#define SCIP_INVALID           1.0E+99  /**< floating point value is not valid */


/*
 * Pointers
 */

#ifndef NULL
#define NULL ((void*)0)                 /**< zero pointer */
#endif


/*
 * Memory settings
 */

#define SCIP_SAFEMEMORY                 /**< use memory leakage detection in debug mode */
#define SCIP_BLOCKMEMORY                /**< use block memory */
#define SCIP_HASHSIZE_NAMES      131101 /**< size of hash table in name tables */
#define SCIP_HASHSIZE_CUTPOOLS   131101 /**< size of hash table in cut pools */



/*
 * Global debugging settings
 */

/*#define DEBUG 1*/
/*#define TODOMESSAGE 1*/


#endif

