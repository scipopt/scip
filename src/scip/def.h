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
#pragma ident "@(#) $Id: def.h,v 1.72 2005/02/14 13:35:42 bzfpfend Exp $"

/**@file   def.h
 * @brief  common defines and data types used in all packages of SCIP
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DEF_H__
#define __DEF_H__


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "scip/type_retcode.h"



#define SCIP_VERSION                 78 /**< SCIP version number (multiplied by 100 to get integer number) */



#define CHECK_ABORT_QUIET(x) do { if( (x) != SCIP_OKAY ) abort(); } while( FALSE )
#define CHECK_OKAY_QUIET(x)  do { RETCODE _restat_; if( (_restat_ = (x)) != SCIP_OKAY ) return _restat_; } while( FALSE )
#define ALLOC_OKAY_QUIET(x)  do { if( NULL == (x) ) return SCIP_NOMEMORY; } while( FALSE )

#define CHECK_ABORT(x) do                                                                                     \
                       {                                                                                      \
                          RETCODE _restat_;                                                                   \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             printf("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_);   \
                             fflush(stdout);                                                                  \
                             abort();                                                                         \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )


#ifndef NDEBUG
#define CHECK_OKAY(x)  do                                                                                     \
                       {                                                                                      \
                          RETCODE _restat_;                                                                   \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             printf("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_);   \
                             fflush(stdout);                                                                  \
                             return _restat_;                                                                 \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

#define ALLOC_OKAY(x)  do                                                                                     \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             printf("[%s:%d] ERROR: No memory in function call\n", __FILE__, __LINE__);       \
                             fflush(stdout);                                                                  \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )
#else
#define CHECK_OKAY(x)  CHECK_OKAY_QUIET(x)
#define ALLOC_OKAY(x)  ALLOC_OKAY_QUIET(x)
#endif


/*
 * Boolean values
 */

#define Bool unsigned int               /**< type used for boolean values */
#ifndef TRUE
#define TRUE  1                         /**< boolean value TRUE */
#define FALSE 0                         /**< boolean value FALSE */
#endif


/*
 * Long Integer values
 */

#ifndef LLONG_MAX
#define LLONG_MAX	9223372036854775807LL
#define LLONG_MIN	(-LLONG_MAX - 1LL)
#endif

#define Longint long long                    /**< type used for long integer values */
#define LONGINT_MAX          LLONG_MAX
#define LONGINT_MIN          LLONG_MIN
#define LONGINT_FORMAT          "%lld"


/*
 * Floating point values
 */

#define Real double                          /**< type used for floating point values */
#define REAL_MAX         (Real)DBL_MAX
#define REAL_MIN        -(Real)DBL_MAX
#define REAL_FORMAT              "%lf"

#define SCIP_DEFAULT_INFINITY         1e+20  /**< default value considered to be infinity */
#define SCIP_DEFAULT_EPSILON          1e-09  /**< default upper bound for floating points to be considered zero */
#define SCIP_DEFAULT_SUMEPSILON       1e-06  /**< default upper bound for sums of floating points to be considered zero */
#define SCIP_DEFAULT_FEASTOL          1e-06  /**< default feasibility tolerance for constraints */
#define SCIP_DEFAULT_DUALFEASTOL      1e-09  /**< default feasibility tolerance for reduced costs */
#define SCIP_DEFAULT_BOUNDSTREPS      1e-04  /**< default minimal improve for strengthening bounds */
#define SCIP_DEFAULT_PSEUDOCOSTEPS    1e-01  /**< default minimal variable distance value to use for pseudo cost updates */
#define SCIP_DEFAULT_PSEUDOCOSTDELTA  1e-04  /**< default minimal objective distance value to use for pseudo cost updates */
#define SCIP_MAXEPSILON               1e-03  /**< maximum value for any numerical epsilon */
#define SCIP_MINEPSILON               1e-20  /**< minimum value for any numerical epsilon */
#define SCIP_INVALID                  1e+99  /**< floating point value is not valid */


#define REALABS(x)        (fabs(x))
#define EPSEQ(x,y,eps)    (REALABS((x)-(y)) <= (eps))
#define EPSLT(x,y,eps)    ((x)-(y) < -(eps))
#define EPSLE(x,y,eps)    ((x)-(y) <= (eps))
#define EPSGT(x,y,eps)    ((x)-(y) > (eps))
#define EPSGE(x,y,eps)    ((x)-(y) >= -(eps))
#define EPSZ(x,eps)       (REALABS(x) <= (eps))
#define EPSP(x,eps)       ((x) > (eps))
#define EPSN(x,eps)       ((x) < -(eps))
#define EPSFLOOR(x,eps)   (floor((x)+(eps)))
#define EPSCEIL(x,eps)    (ceil((x)-(eps)))
#define EPSFRAC(x,eps)    ((x)-EPSFLOOR(x,eps))
#define EPSISINT(x,eps)   (EPSCEIL(x,eps)-(x) <= (eps))


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

#ifndef MAX3
#define MAX3(x,y,z) ((x) >= (y) ? MAX(x,z) : MAX(y,z))
#define MIN3(x,y,z) ((x) <= (y) ? MIN(x,z) : MIN(y,z))
#endif


/*
 * Pointers
 */

#ifndef NULL
#define NULL ((void*)0)                 /**< zero pointer */
#endif


/*
 * Strings
 */

#define MAXSTRLEN                 1024  /**< maximum string length in SCIP */


/*
 * Memory settings
 */

#define SCIP_HASHSIZE_NAMES      131101 /**< size of hash table in name tables */
#define SCIP_HASHSIZE_CUTPOOLS   131101 /**< size of hash table in cut pools */
#define SCIP_HASHSIZE_PARAMS       4099 /**< size of hash table in cut pools */
#define SCIP_HASHSIZE_VBC        131101 /**< size of hash map for node -> nodenum mapping used for VBC output */

/*#define NOBLOCKMEM*/



/*
 * Global debugging settings
 */

/*#define DEBUG*/


#endif
