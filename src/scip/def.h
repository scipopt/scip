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
 * @brief  comon defines used in all packages of SCIP
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DEF_H__
#define __DEF_H__

#define SCIP_VERSION_MAJOR  0           /**< major part of SCIP version */
#define SCIP_VERSION_MINOR  1           /**< minor part of SCIP version */

typedef int Bool;                       /**< type used for boolean values */

#ifndef NULL
#define NULL 0                          /**< zero pointer */
#endif

#ifndef TRUE
#define TRUE (0==0)                     /**< boolean value TRUE */
#define FALSE (0==1)                    /**< boolean value FALSE */
#endif

#define SCIP_INFINITY    1.0E+20        /**< value considered to be infinity */

#define SCIP_MAXNCOL       0x0000fffff  /**< maximal number of columns; 20 bits available */
#define SCIP_MAXNROW       0x0000fffff  /**< maximal number of rows; 20 bits available */
#define SCIP_MAXNCHILDREN  0x00000ffff  /**< maximal number of children per node; 16 bits available */
#define SCIP_MAXNADDEDROWS 0x0000fffff  /**< maximal number of added rows per node; 20 bits available */


enum Setting                            /**< possible settings for enabling/disabling algorithms and other features */
{
   SCIP_UNDEFINED = 0,                  /**< undefined setting */
   SCIP_DISABLED  = 1,                  /**< feature is disabled */
   SCIP_AUTO      = 2,                  /**< feature is set to automatic mode */
   SCIP_ENABLED   = 3                   /**< feature is enabled */
};
typedef enum Setting SETTING;



#define CHECK_OKAY(x) { int _restat_; if( (_restat_ = (x)) < SCIP_OKAY ) return _restat_; }
#define ALLOC_OKAY(x) { if( NULL == (x) ) return SCIP_NOMEMORY; }
#define ALLOC_NULL(x) { if( NULL == (x) ) return NULL; }

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

#endif

