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

/**@file   set.h
 * @brief  global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SET_H__
#define __SET_H__


typedef struct Set SET;                 /**< global SCIP settings */


#include "def.h"
#include "sort.h"


/** global SCIP settings */
struct Set
{
   Real             epsZero;            /**< absolute values smaller than this are considered zero */
   Real             infinity;           /**< values larger than this are considered infinity */
   Real             memGrowFac;         /**< memory growing factor for dynamically allocated arrays */
   int              memGrowInit;        /**< initial size of dynamically allocated arrays */
   Real             bufGrowFac;         /**< memory growing factor for buffer arrays */
   int              bufGrowInit;        /**< initial size of buffer arrays */
   Real             treeGrowFac;        /**< memory growing factor for tree array */
   int              treeGrowInit;       /**< initial size of tree array */
   Real             pathGrowFac;        /**< memory growing factor for path array */
   int              pathGrowInit;       /**< initial size of path array */
   DECL_SORTPTRCOMP((*nodecmp));        /**< compares two nodes regarding the order in the leaf list */
};


extern
RETCODE SCIPsetCreate(                  /**< creates global SCIP settings */
   SET**            set                 /**< pointer to SCIP settings */
   );

extern
int SCIPcalcMemGrowSize(                /**< calculate memory size for dynamically allocated arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

extern
int SCIPcalcBufGrowSize(                /**< calculate memory size for buffer arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

extern
int SCIPcalcPathGrowSize(               /**< calculate memory size for path array */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

extern
Bool SCIPisEQ(                          /**< checks, if values are in range of epsZero */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisL(                           /**< checks, if val1 is (more than epsZero) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisLE(                          /**< checks, if val1 is not (more than epsZero) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisG(                           /**< checks, if val1 is (more than epsZero) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisGE(                          /**< checks, if val1 is not (more than epsZero) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisZero(                        /**< checks, if value is in range epsZero of 0.0 */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Real SCIPinfinity(                      /**< returns infinity value */
   const SET*       set                 /**< global SCIP settings */
   );

extern
Bool SCIPisInfinity(                    /**< checks, if value is infinite */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   );

extern
Bool SCIPisPos(                         /**< checks, if value is greater than epsZero */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisNeg(                         /**< checks, if value is lower than -epsZero */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

#endif
