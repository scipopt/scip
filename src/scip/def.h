/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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


#define SCIP_VERSION                 10 /**< SCIP version number (multiplied by 100 to get integer number) */



#define CHECK_OKAY(x) { int _restat_; if( (_restat_ = (x)) < SCIP_OKAY ) return _restat_; }
#define ALLOC_OKAY(x) { if( NULL == (x) ) return SCIP_NOMEMORY; }

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



/*
 * Boolean values
 */

typedef int Bool;                       /**< type used for boolean values */
#ifndef TRUE
#define TRUE (0==0)                     /**< boolean value TRUE */
#define FALSE (0==1)                    /**< boolean value FALSE */
#endif


/*
 * Floating point values
 */

typedef double Real;                    /**< type used for floating point values */
#define SCIP_DEFAULT_EPSILON     1e-09  /**< default upper bound for floating points to be considered zero */
#define SCIP_DEFAULT_INFINITY  1.0E+20  /**< default value considered to be infinity */
#define SCIP_INVALID           1.0E+99  /**< floating point value is not valid */


/*
 * Pointers
 */

#ifndef NULL
#define NULL ((void*)0)                 /**< zero pointer */
#endif


/*
 * Message Output
 */

#define SCIP_DEFAULT_VERBLEVEL    SCIP_VERBLEVEL_NORMAL


/*
 * Dynamic Memory
 */

#define SCIP_DEFAULT_MEMGROWFAC       1.2
#define SCIP_DEFAULT_MEMGROWINIT      4
#define SCIP_DEFAULT_BUFGROWFAC       2.0
#define SCIP_DEFAULT_BUFGROWINIT  65536
#define SCIP_DEFAULT_TREEGROWFAC      2.0
#define SCIP_DEFAULT_TREEGROWINIT 65536
#define SCIP_DEFAULT_PATHGROWFAC      2.0
#define SCIP_DEFAULT_PATHGROWINIT   256


/*
 * Pricing
 */

#define SCIP_DEFAULT_MAXPRICEVARS    16      /**< maximal number of variables priced in per pricing round */
#define SCIP_PRICE_SCALELOOSE         1.0E-4 /**< scaling for pricing score of loose variables */


/*
 * Debugging
 */

#define DEBUG 1



#endif

