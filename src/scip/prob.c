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
#include <string.h>

#include "def.h"
#include "prob.h"


/** main problem to solve */
struct Prob
{
   char*            name;               /**< problem name */
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
 * problem information
 */

const char* SCIPprobGetName(            /**< gets problem name */
   const PROB*      prob                /**< problem data */
   )
{
   assert(prob != NULL);
   return prob->name;
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
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);

   CHECK_OKAY( SCIPconslistAdd(&(prob->conslist), memhdr, cons) );

   return SCIP_OKAY;
}



RETCODE SCIPprobCreate(                 /**< creates problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name                /**< problem name */
   )
{
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(*prob) );
   ALLOC_OKAY( allocMemoryArray((*prob)->name, strlen(name)+1) );
   copyMemoryArray((*prob)->name, name, strlen(name)+1);

   return SCIP_OKAY;
}

RETCODE SCIPprobFree(                   /**< frees problem data structure */
   PROB**           prob                /**< pointer to problem data structure */
   )
{
   int c;

   assert(prob != NULL);
   assert(*prob != NULL);

   freeMemoryArray((*prob)->name);
   freeMemoryArrayNull((*prob)->cols);
   /* conslist doesn't need to be freed, because it's in block memory
      SCIPconslistFree(&(*prob)->conslist, memhdr);
   */
   freeMemory(*prob);
   
   return SCIP_OKAY;
}

RETCODE SCIPprobDuplicate(              /**< duplicates problem data */
   PROB**           prob,               /**< pointer to target problem data structure */
   MEMHDR*          memhdr,             /**< block memory of new problem data */
   PROB*            source              /**< problem to duplicate */
   )
{
   errorMessage("Method not correctly implemented");

   /* copy all columns to new block memory,
    * copy all constraints to new block memory,
    * initialize cols and conslist correclty
    */

   return SCIPprobCreate(prob, source->name);
}
