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
#define SCIP_HASHSIZE_NAMES       23663


/*
 * Pricing
 */

#define SCIP_DEFAULT_USEPRICING   FALSE /**< activate pricing of variables */
#define SCIP_DEFAULT_MAXPRICEVARS    16 /**< maximal number of variables priced in per pricing round */



/*
 * Cut Separation
 */

#define SCIP_DEFAULT_MAXSEPACUTS     32 /**< maximal number of cuts separated per separation round */



/*
 * Primal Solutions
 */

#define SCIP_DEFAULT_MAXSOL         256 /**< maximal number of solutions to store in the solution storage */


/*
 * Tree
 */

#define SCIP_DEFAULT_NODELIMIT  1000000 /**< maximal number of nodes to create */


/*
 * Display
 */

#define SCIP_DEFAULT_DISPWIDTH      140 /**< maximal number of characters in a node information line */
#define SCIP_DEFAULT_DISPFREQ     10000 /**< frequency for displaying node information lines */
#define SCIP_DEFAULT_DISPHEADERFREQ  15 /**< frequency for displaying header lines (every n'th node information line) */



/*
 * Block Memory
 */

#define SCIP_SAFEMEMORY                 /**< use memory leakage detection in debug mode */
#define SCIP_BLOCKMEMORY                /**< use block memory */



/*
 * Debugging
 */

/*#define DEBUG 1*/
/*#define TODOMESSAGE 1*/


#endif

