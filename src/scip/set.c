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

/**@file   set.c
 * @brief  global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "memory.h"
#include "tree.h"
#include "set.h"


#define DEFAULT_MEMGROWFAC        1.2
#define DEFAULT_MEMGROWINIT       4
#define DEFAULT_BUFGROWFAC        2.0
#define DEFAULT_BUFGROWINIT   65536
#define DEFAULT_TREEGROWFAC       2.0
#define DEFAULT_TREEGROWINIT  65536
#define DEFAULT_PATHGROWFAC       2.0
#define DEFAULT_PATHGROWINIT    256
#define DEFAULT_NODECMP      &SCIPnodeCmpLowerbound


RETCODE SCIPsetCreate(                  /**< creates global SCIP settings */
   SET**            set                 /**< pointer to SCIP settings */
   )
{
   assert(set != NULL);

   ALLOC_OKAY( allocMemory(*set) );

   (*set)->epsZero = SCIP_DEFAULT_EPSZERO;
   (*set)->infinity = SCIP_DEFAULT_INFINITY;
   (*set)->memGrowFac = DEFAULT_MEMGROWFAC;
   (*set)->memGrowInit = DEFAULT_MEMGROWINIT;
   (*set)->bufGrowFac = DEFAULT_BUFGROWFAC;
   (*set)->bufGrowInit = DEFAULT_BUFGROWINIT;
   (*set)->treeGrowFac = DEFAULT_TREEGROWFAC;
   (*set)->treeGrowInit = DEFAULT_TREEGROWINIT;
   (*set)->pathGrowFac = DEFAULT_PATHGROWFAC;
   (*set)->pathGrowInit = DEFAULT_PATHGROWINIT;
   (*set)->nodecmp = DEFAULT_NODECMP;

   return SCIP_OKAY;
}

static
int calcGrowSize(                       /**< calculate memory size for dynamically allocated arrays */
   int              initsize,           /**< initial size of array */
   Real             growfac,            /**< growing factor of array */
   int              num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(initsize >= 0);
   assert(growfac >= 1.0);

   /* calculate the size with this loop, such that the resulting numbers are allways the same (-> block memory) */
   size = initsize;
   while( size < num )
      size = growfac * size + 1;

   return size;
}

int SCIPcalcMemGrowSize(                /**< calculate memory size for dynamically allocated arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->memGrowInit, set->memGrowFac, num);
}

int SCIPcalcBufGrowSize(                /**< calculate memory size for buffer arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->bufGrowInit, set->bufGrowFac, num);
}

int SCIPcalcPathGrowSize(               /**< calculate memory size for path array */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->pathGrowInit, set->pathGrowFac, num);
}


Bool SCIPisEQ(                          /**< checks, if values are in range of epsZero */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( ABS(val1-val2) < set->epsZero );
}

Bool SCIPisL(                           /**< checks, if val1 is (more than epsZero) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 < val2 - set->epsZero );
}

Bool SCIPisLE(                          /**< checks, if val1 is not (more than epsZero) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 <= val2 + set->epsZero );
}

Bool SCIPisG(                           /**< checks, if val1 is (more than epsZero) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 > val2 + set->epsZero );
}

Bool SCIPisGE(                          /**< checks, if val1 is not (more than epsZero) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   return( val1 >= val2 - set->epsZero );
}

Bool SCIPisZero(                        /**< checks, if value is in range epsZero of 0.0 */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( ABS(val) <= set->epsZero );
}

Real SCIPinfinity(                      /**< returns infinity value */
   const SET*       set                 /**< global SCIP settings */
   )
{
   return set->infinity;
}

Bool SCIPisInfinity(                    /**< checks, if value is infinite */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   )
{
   return( val >= set->infinity );
}

Bool SCIPisPos(                         /**< checks, if value is greater than epsZero */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( val > set->epsZero );
}

Bool SCIPisNeg(                         /**< checks, if value is lower than -epsZero */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   )
{
   return( val < -set->epsZero );
}

