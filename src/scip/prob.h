/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
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


/** objective sense: minimization or maximization */
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
   VAR**            vars;               /**< array with active variables ordered binary, integer, implicit, continous */
   HASHTABLE*       varnames;           /**< hash table storing variable's names */
   CONS**           conss;              /**< array with constraints of the problem */
   HASHTABLE*       consnames;          /**< hash table storing constraints' names */
   OBJSENSE         objsense;           /**< objective sense */
   Real             objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   Real             objlim;             /**< objective limit for non-fixed variables */
   int              fixedvarssize;      /**< available slots in fixedvars array */
   int              nfixedvars;         /**< number of fixed and aggregated variables in the problem */
   int              varssize;           /**< available slots in vars array */
   int              nvars;              /**< number of mutable variables in the problem (used slots in vars array) */
   int              nbin;               /**< number of binary variables */
   int              nint;               /**< number of general integer variables */
   int              nimpl;              /**< number of implicit integer variables */
   int              ncont;              /**< number of continous variables */
   int              consssize;          /**< available slots in conss array */
   int              nconss;             /**< number of constraints in the problem (number of used slots in conss array) */
};


/*
 * problem creation
 */

/** creates problem data structure */
extern
RETCODE SCIPprobCreate(
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name                /**< problem name */
   );

/** frees problem data structure */
extern
RETCODE SCIPprobFree(
   PROB**           prob,               /**< pointer to problem data structure */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's the original problem) */
   );

/** transform problem data into normalized form */
extern
RETCODE SCIPprobTransform(
   PROB*            source,             /**< problem to transform */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB**           target              /**< pointer to target problem data structure */
   );

/** activates constraints in the problem */
extern
RETCODE SCIPprobActivate(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** deactivates constraints in the problem */
extern
RETCODE SCIPprobDeactivate(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );


/*
 * problem modification
 */

/** adds variable to the problem and captures it */
extern
RETCODE SCIPprobAddVar(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable to add */
   );

/** changes the type of a variable in the problem */
extern
RETCODE SCIPprobChgVarType(
   PROB*            prob,               /**< problem data */
   VAR*             var,                /**< variable to add */
   VARTYPE          vartype             /**< new type of variable */
   );

/** adds constraint to the problem and captures it */
extern
RETCODE SCIPprobAddCons(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   );

/** releases and removes constraint from the problem; if the user has not captured the constraint for his own use, the
 *  constraint may be invalid after the call
 */
extern
RETCODE SCIPprobDelCons(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   );

/** sets objective sense: minimization or maximization */
extern
void SCIPprobSetObjsense(
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   );

/** sets limit on objective function, such that only solutions better than this limit are accepted */
extern
void SCIPprobSetObjlim(
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< objective limit */
   );

/** returns the external value of the given internal objective value */
extern
Real SCIPprobExternObjval(
   PROB*            prob,               /**< problem data */
   Real             objval              /**< internal objective value */
   );


/*
 * problem information
 */

/** gets problem name */
extern
const char* SCIPprobGetName(
   PROB*            prob                /**< problem data */
   );

/** returns variable of the problem with given name */
extern
VAR* SCIPprobFindVar(
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   );

/** returns constraint of the problem with given name */
extern
CONS* SCIPprobFindCons(
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   );

/** displays actual pseudo solution */
extern
void SCIPprobPrintPseudoSol(
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   );



#endif
