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


/**< objective sense: minimization or maximization */
enum Objsense
{
   SCIP_OBJSENSE_MAXIMIZE = -1,         /**< maximization of objective function */
   SCIP_OBJSENSE_MINIMIZE = +1          /**< minimization of objective function (the default) */
};
typedef enum Objsense OBJSENSE;

typedef struct Prob PROB;               /**< main problem to solve */


#include "memory.h"
#include "retcode.h"
#include "cons.h"
#include "lp.h"
#include "var.h"
#include "stat.h"


/** main problem to solve */
struct Prob
{
   char*            name;               /**< problem name */
   VAR**            fixedvars;          /**< array with fixed and aggregated variables */
   VAR**            vars;               /**< array with non-fixed variables ordered binary, integer, implicit, continous */
   CONSLIST*        conslist;           /**< list of constraints of the problem */
   OBJSENSE         objsense;           /**< objective sense */
   Real             objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   Real             objlim;             /**< objective limit for non-fixed variables */
   int              fixedvarssize;      /**< available slots in fixedvars array */
   int              nfixedvars;         /**< number of fixed and aggregated variables in the problem */
   int              varssize;           /**< available slots in vars array */
   int              nvars;              /**< number of variables in the problem (number of used slots in vars array) */
   int              nbin;               /**< number of binary variables */
   int              nint;               /**< number of general integer variables */
   int              nimpl;              /**< number of implicit integer variables */
   int              ncont;              /**< number of continous variables */
   int              ncons;              /**< number of constraints in the problem (number of used slots in cons array) */
};


/*
 * problem creation
 */

extern
RETCODE SCIPprobCreate(                 /**< creates problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name                /**< problem name */
   );

extern
RETCODE SCIPprobFree(                   /**< frees problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's the original problem) */
   );

extern
RETCODE SCIPprobTransform(              /**< transform problem data into normalized form */
   PROB*            source,             /**< problem to transform */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB**           target              /**< pointer to target problem data structure */
   );

extern
RETCODE SCIPprobActivate(               /**< activates constraints in the problem */
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPprobDeactivate(             /**< deactivates constraints in the problem */
   PROB*            prob                /**< problem data */
   );


/*
 * problem modification
 */

extern
RETCODE SCIPprobAddVar(                 /**< adds variable to the problem and captures it */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable to add */
   );

extern
RETCODE SCIPprobAddCons(                /**< adds constraint to the problem and captures it */
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   );

extern
void SCIPprobSetObjsense(               /**< sets objective sense: minimization or maximization */
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   );

extern
void SCIPprobSetObjlim(                 /**< sets limit on objective function, such that only solutions better than this
                                           limit are accepted */
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< objective limit */
   );

extern
Real SCIPprobExternObjval(              /**< returns the external value of the given internal objective value */
   PROB*            prob,               /**< problem data */
   Real             objval              /**< internal objective value */
   );


/*
 * problem information
 */

extern
const char* SCIPprobGetName(            /**< gets problem name */
   const PROB*      prob                /**< problem data */
   );



#endif
