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

/**@file   prob.h
 * @brief  Methods and datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PROB_H__
#define __PROB_H__


typedef struct Prob PROB;               /**< main problem to solve */


#include "retcode.h"
#include "constraint.h"
#include "lp.h"


/*
 * problem modification
 */

extern
RETCODE SCIPprobAddCol(                 /**< adds variable to the problem */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   COL*             col                 /**< variable to add */
   );

extern
RETCODE SCIPprobAddConstraint(          /**< adds constraint to the problem */
   PROB*            prob,               /**< problem data */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   );



/*
 * problem creation
 */

extern
RETCODE SCIPprobCreate(                 /**< creates problem data structure */
   PROB**           prob                /**< pointer to problem data structure */
   );


#endif
