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
typedef struct ProbData PROBDATA;       /**< user problem data set by the reader */


/** frees user problem data
 *
 *  input:
 *    scip            : SCIP main data structure
 *    probdata        : pointer to the user problem data to free
 */
#define DECL_PROBDELETE(x) RETCODE x (SCIP* scip, PROBDATA** probdata)

/** transforms user problem data into data belonging to the transformed problem
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sourcedata      : source problem data to transform
 *    targetdata      : pointer to store created transformed problem data
 */
#define DECL_PROBTRANS(x) RETCODE x (SCIP* scip, PROBDATA* sourcedata, PROBDATA** targetdata)



#include "memory.h"
#include "retcode.h"
#include "cons.h"
#include "tree.h"
#include "lp.h"
#include "var.h"
#include "stat.h"
#include "branch.h"


/** main problem to solve */
struct Prob
{
   char*            name;               /**< problem name */
   DECL_PROBDELETE  ((*probdelete));    /**< frees user problem data */
   DECL_PROBTRANS   ((*probtrans));     /**< transforms user problem data into data belonging to the transformed problem */
   PROBDATA*        probdata;           /**< user problem data set by the reader */
   VAR**            fixedvars;          /**< array with fixed and aggregated variables */
   VAR**            vars;               /**< array with active variables ordered binary, integer, implicit, continous */
   HASHTABLE*       varnames;           /**< hash table storing variable's names */
   CONS**           conss;              /**< array with constraints of the problem */
   HASHTABLE*       consnames;          /**< hash table storing constraints' names */
   OBJSENSE         objsense;           /**< objective sense */
   Real             objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   Real             objlim;             /**< objective limit as external value */
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
   int              maxnconss;          /**< maximum number of constraints existing at the same time */
   unsigned int     transformed:1;      /**< TRUE iff problem is the transformed problem */
};


/*
 * problem creation
 */

/** creates problem data structure */
extern
RETCODE SCIPprobCreate(
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name,               /**< problem name */
   DECL_PROBDELETE  ((*probdelete)),    /**< frees user problem data */
   DECL_PROBTRANS   ((*probtrans)),     /**< transforms user problem data into data belonging to the transformed problem */
   PROBDATA*        probdata,           /**< user problem data set by the reader */
   Bool             transformed         /**< is this the transformed problem? */
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
   TREE*            tree,               /**< branch and bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PROB**           target              /**< pointer to target problem data structure */
   );




/*
 * problem modification
 */

/** sets user problem data */
extern
void SCIPprobSetData(
   PROB*            prob,               /**< problem */
   PROBDATA*        probdata            /**< user problem data to use */
   );

/** adds variable to the problem and captures it */
extern
RETCODE SCIPprobAddVar(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< variable to add */
   );

/** changes the type of a variable in the problem */
extern
RETCODE SCIPprobChgVarType(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var,                /**< variable to add */
   VARTYPE          vartype             /**< new type of variable */
   );

/** informs problem, that the given variable was fixed */
extern
RETCODE SCIPprobVarFixed(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< variable to add */
   );

/** adds constraint to the problem and captures it; a local constraint is automatically upgraded into a global constraint */
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

/** resets maximum number of constraints to current number of constraints */
extern
void SCIPprobResetMaxNConss(
   PROB*            prob                /**< problem data */
   );

/** sets objective sense: minimization or maximization */
extern
void SCIPprobSetObjsense(
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   );

/** increases objective offset */
extern
void SCIPprobIncObjoffset(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             incval              /**< value to add to objective offset */
   );

/** sets limit on objective function, such that only solutions better than this limit are accepted */
extern
void SCIPprobSetExternObjlim(
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< external objective limit */
   );

/** sets limit on objective function as transformed internal objective value */
extern
void SCIPprobSetInternObjlim(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objlim              /**< transformed internal objective limit */
   );




/*
 * problem information
 */

/** gets problem name */
extern
const char* SCIPprobGetName(
   PROB*            prob                /**< problem data */
   );

/** gets user problem data */
extern
PROBDATA* SCIPprobGetData(
   PROB*            prob                /**< problem */
   );

/** returns the external value of the given internal objective value */
extern
Real SCIPprobExternObjval(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objval              /**< internal objective value */
   );

/** returns the internal value of the given external objective value */
extern
Real SCIPprobInternObjval(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objval              /**< external objective value */
   );

/** gets limit on objective function in external space */
extern
Real SCIPprobGetExternObjlim(
   PROB*            prob                /**< problem data */
   );

/** gets limit on objective function as transformed internal objective value */
extern
Real SCIPprobGetInternObjlim(
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
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

/** outputs problem statistics */
extern
void SCIPprobPrintStatistics(
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   );



#endif
