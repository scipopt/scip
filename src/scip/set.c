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
#include "set.h"

#define DEFAULT_EPSZERO       1e-09
#define DEFAULT_MEMGROWFAC    1.2
#define DEFAULT_MEMGROWADD    4
#define DEFAULT_MEMGROWINIT   4

SET* SCIPsetCreate(                     /**< creates global SCIP settings */
   void
   )
{
   SET* set;
   
   ALLOC_NULL( allocMemory(set) );
   set->epsZero = DEFAULT_EPSZERO;
   set->memGrowFac = DEFAULT_MEMGROWFAC;
   set->memGrowAdd = DEFAULT_MEMGROWADD;
   set->memGrowInit = DEFAULT_MEMGROWINIT;

   return set;
}

int SCIPcalcMemGrowSize(                /**< calculate memory size for dynamically allocated arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(set->memGrowAdd >= 1);
   assert(set->memGrowInit >= 0);

   /* calculate the size with this loop, such that the resulting numbers are allways the same (-> block memory) */
   size = set->memGrowInit;
   while( size < num )
      size = set->memGrowFac * size + set->memGrowAdd;

   return size;
}
