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

/**@file   prob.c
 * @brief  Methods and datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "prob.h"


/** main problem to solve */
struct Prob
{
   COL**            cols;               /**< array with problem variables */
   CONSLIST*        conslist;           /**< list of constraints of the problem */
   int              colssize;           /**< available slots in cols vector */
   int              ncols;              /**< number of variables in the problem (number of used slots in cols vector) */
};



/*
 * dymanic memory arrays
 */

static
RETCODE probEnsureColsMem(              /**< resizes cols array to be able to store at least num entries */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->colssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(prob->cols, newsize) );
      prob->colssize = newsize;
   }
   assert(num <= prob->colssize);

   return SCIP_OKAY;
}



/*
 * problem modification
 */

RETCODE SCIPprobAddCol(                 /**< adds variable to the problem */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< variable to add */
   )
{
   assert(prob != NULL);
   assert(set != NULL);
   assert(col != NULL);

   CHECK_OKAY( probEnsureColsMem(prob, set, prob->ncols+1) );
   prob->cols[prob->ncols] = col;
   prob->ncols++;
   SCIPcolCapture(col);

   return SCIP_OKAY;
}

RETCODE SCIPprobAddConstraint(          /**< adds constraint to the problem */
   PROB*            prob,               /**< problem data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(mem != NULL);
   assert(set != NULL);
   assert(cons != NULL);

   CHECK_OKAY( SCIPconslistAdd(&(prob->conslist), mem, cons) );

   return SCIP_OKAY;
}



RETCODE SCIPprobCreate(                 /**< creates problem data structure */
   PROB**           prob                /**< pointer to problem data structure */
   )
{
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(*prob) );

   return SCIP_OKAY;
}
