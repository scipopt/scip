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

#define SCIP_VERSION  0.1               /**< SCIP version number */

#define SCIP_MAXNCOL       0x0000fffff  /**< maximal number of columns; 20 bits available */
#define SCIP_MAXNROW       0x0000fffff  /**< maximal number of rows; 20 bits available */
#define SCIP_MAXNCHILDREN  0x00000ffff  /**< maximal number of children per node; 16 bits available */
#define SCIP_MAXNADDEDROWS 0x0000fffff  /**< maximal number of added rows per node; 20 bits available */

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
#define SCIP_DEFAULT_EPSZERO     1e-09  /**< default upper bound for floating points to be considered zero */
#define SCIP_DEFAULT_INFINITY  1.0E+20  /**< default value considered to be infinity */

/*
 * Pointers
 */

#ifndef NULL
#define NULL 0                          /**< zero pointer */
#endif




#endif

