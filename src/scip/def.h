/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   def.h
 * @brief  common defines and data types used in all packages of SCIP
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEF_H__
#define __SCIP_DEF_H__


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

/*
 * Boolean values
 */

#ifndef SCIP_Bool
#define SCIP_Bool unsigned int                    /**< type used for boolean values */
#ifndef TRUE
#define TRUE  1                         /**< boolean value TRUE */
#define FALSE 0                         /**< boolean value FALSE */
#endif
#endif

/*
 * Define the marco EXTERN and some functions depending if the OS is Windows or not
 */
#if defined(_WIN32) || defined(_WIN64)

#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#define getcwd _getcwd

#define EXTERN __declspec(dllexport)

#else
#define EXTERN extern
#endif



#include "scip/type_retcode.h"
#include "scip/pub_message.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SCIP_VERSION                301 /**< SCIP version number (multiplied by 100 to get integer number) */
#define SCIP_SUBVERSION               0 /**< SCIP sub version number */
#define SCIP_COPYRIGHT   "Copyright (c) 2002-2013 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)"


/*
 * CIP format variable characters
 */

#define SCIP_VARTYPE_BINARY_CHAR 'B'
#define SCIP_VARTYPE_INTEGER_CHAR 'I'
#define SCIP_VARTYPE_IMPLINT_CHAR 'M'
#define SCIP_VARTYPE_CONTINUOUS_CHAR 'C'

/*
 * Long Integer values
 */

#ifndef LLONG_MAX
#define LLONG_MAX        9223372036854775807LL
#define LLONG_MIN        (-LLONG_MAX - 1LL)
#endif

#define SCIP_Longint long long                         /**< type used for long integer values */
#define SCIP_LONGINT_MAX          LLONG_MAX
#define SCIP_LONGINT_MIN          LLONG_MIN
#ifndef SCIP_LONGINT_FORMAT
#if defined(_WIN32) || defined(_WIN64)
#define SCIP_LONGINT_FORMAT           "I64d"
#else
#define SCIP_LONGINT_FORMAT           "lld"
#endif
#endif


/*
 * Floating point values
 */

#define SCIP_Real double                               /**< type used for floating point values */
#define SCIP_REAL_MAX         (SCIP_Real)DBL_MAX
#define SCIP_REAL_MIN        -(SCIP_Real)DBL_MAX
#define SCIP_REAL_FORMAT               "lf"

#define SCIP_DEFAULT_INFINITY         1e+20  /**< default value considered to be infinity */
#define SCIP_DEFAULT_EPSILON          1e-09  /**< default upper bound for floating points to be considered zero */
#define SCIP_DEFAULT_SUMEPSILON       1e-06  /**< default upper bound for sums of floating points to be considered zero */
#define SCIP_DEFAULT_FEASTOL          1e-06  /**< default feasibility tolerance for constraints */
#define SCIP_DEFAULT_LPFEASTOL        1e-06  /**< default primal feasibility tolerance of LP solver */
#define SCIP_DEFAULT_DUALFEASTOL      1e-06  /**< default feasibility tolerance for reduced costs */
#define SCIP_DEFAULT_BARRIERCONVTOL   1e-10  /**< default convergence tolerance used in barrier algorithm */
#define SCIP_DEFAULT_BOUNDSTREPS       0.05  /**< default minimal relative improve for strengthening bounds */
#define SCIP_DEFAULT_PSEUDOCOSTEPS    1e-01  /**< default minimal variable distance value to use for pseudo cost updates */
#define SCIP_DEFAULT_PSEUDOCOSTDELTA  1e-04  /**< default minimal objective distance value to use for pseudo cost updates */
#define SCIP_DEFAULT_RECOMPFAC        1e+07  /**< default minimal decrease factor that causes the recomputation of a value (e.g., pseudo objective) instead of an update */
#define SCIP_DEFAULT_HUGEVAL          1e+15  /**< values larger than this are considered huge and should be handled separately (e.g., in activity computation) */
#define SCIP_MAXEPSILON               1e-03  /**< maximum value for any numerical epsilon */
#define SCIP_MINEPSILON               1e-20  /**< minimum value for any numerical epsilon */
#define SCIP_INVALID                  1e+99  /**< floating point value is not valid */
#define SCIP_UNKNOWN                  1e+98  /**< floating point value is not known (in primal solution) */


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
#define EPSROUND(x,eps)   (ceil((x)-0.5+(eps)))
#define EPSFRAC(x,eps)    ((x)-EPSFLOOR(x,eps))
#define EPSISINT(x,eps)   (EPSFRAC(x,eps) <= (eps))


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
#define MAX3(x,y,z) ((x) >= (y) ? MAX(x,z) : MAX(y,z))  /**< returns maximum of x, y, and z */
#define MIN3(x,y,z) ((x) <= (y) ? MIN(x,z) : MIN(y,z))  /**< returns minimum of x, y, and z */
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

#define SCIP_MAXSTRLEN             1024 /**< maximum string length in SCIP */

/*
 * Memory settings
 */

#define SCIP_HASHSIZE_PARAMS         4099 /**< size of hash table in parameter name tables */
#define SCIP_HASHSIZE_NAMES          131101 /**< size of hash table in name tables */
#define SCIP_HASHSIZE_CUTPOOLS       131101 /**< size of hash table in cut pools */
#define SCIP_HASHSIZE_CLIQUES        131101 /**< size of hash table in clique tables */
#define SCIP_HASHSIZE_NAMES_SMALL    8011   /**< size of hash table in name tables for small problems */
#define SCIP_HASHSIZE_CUTPOOLS_SMALL 8011   /**< size of hash table in cut pools for small problems */
#define SCIP_HASHSIZE_CLIQUES_SMALL  8011   /**< size of hash table in clique tables for small problems */
#define SCIP_HASHSIZE_VBC            131101 /**< size of hash map for node -> nodenum mapping used for VBC output */

/*#define BMS_NOBLOCKMEM*/


/*
 * Global debugging settings
 */

/*#define DEBUG*/


/*
 * Defines for handling SCIP return codes
 */

/** this macro is used to stop SCIP in debug mode such that errors can be debugged;
 *
 *  @note In optimized mode this macro has no effect. That means, in case of an error it has to be ensured that code
 *        terminates with an error code or continues safely.
 */
#define SCIPABORT() assert(FALSE)

#define SCIP_CALL_ABORT_QUIET(x)  do { if( (x) != SCIP_OKAY ) SCIPABORT(); } while( FALSE )
#define SCIP_CALL_QUIET(x)        do { SCIP_RETCODE _restat_; if( (_restat_ = (x)) != SCIP_OKAY ) return _restat_; } while( FALSE )
#define SCIP_ALLOC_ABORT_QUIET(x) do { if( NULL == (x) ) SCIPABORT(); } while( FALSE )
#define SCIP_ALLOC_QUIET(x)       do { if( NULL == (x) ) return SCIP_NOMEMORY; } while( FALSE )

#define SCIP_CALL_ABORT(x) do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             SCIPABORT();                                                                     \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_ALLOC_ABORT(x) do                                                                                \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n", __FILE__, __LINE__);            \
                             SCIPABORT();                                                                     \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             return _restat_;                                                                 \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_ALLOC(x)  do                                                                                     \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

/*
 * Define to mark deprecated API functions
 */

#if defined(_WIN32) || defined(_WIN64)
#  define SCIP_DEPRECATED __declspec(deprecated)
#elif defined(__GNUC__) && defined(__linux__)
#  define SCIP_DEPRECATED __attribute__ ((deprecated))
#else
#  define SCIP_DEPRECATED
#endif

#ifdef __cplusplus
}
#endif

#endif
