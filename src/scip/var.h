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

/**@file   var.h
 * @brief  Methods and datastructures for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __VAR_H__
#define __VAR_H__



/** status of problem variables */
enum Varstatus
{
   SCIP_VARSTATUS_ORIGINAL   = 0,       /**< variable belongs to original problem */
   SCIP_VARSTATUS_LOOSE      = 1,       /**< variable is a loose variable of the transformed problem */
   SCIP_VARSTATUS_COLUMN     = 2,       /**< variable is a column of the transformed problem */
   SCIP_VARSTATUS_FIXED      = 3,       /**< variable is fixed to specific value in the transformed problem */
   SCIP_VARSTATUS_AGGREGATED = 4,       /**< variable is aggregated to $x = a*y + c$ in the transformed problem */
   SCIP_VARSTATUS_MULTAGGR   = 5        /**< variable is aggregated to $x = a_1*y_1 + ... + a_k*y_k + c$ */
};
typedef enum Varstatus VARSTATUS;

/** variable type */
enum Vartype
{
   SCIP_VARTYPE_BINARY    = 0,          /**< binary variable: $x \in \{0,1\}$ */
   SCIP_VARTYPE_INTEGER   = 1,          /**< integer variable: $x \in \{lb, \ldots, \ub\}$ */
   SCIP_VARTYPE_IMPLINT   = 2,          /**< implicit integer variable: continous variable, that is always integral */
   SCIP_VARTYPE_CONTINOUS = 3           /**< continous variable: $x \in [lb,ub] */
};
typedef enum Vartype VARTYPE;

typedef struct DomChg DOMCHG;           /**< changes in domains of variables (fixed sized arrays) */
typedef struct DomChgDyn DOMCHGDYN;     /**< changes in domains of variables (dynamically sized arrays) */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct Dom DOM;                 /**< datastructures for storing domains of variables */
typedef struct Aggregate AGGREGATE;     /**< aggregation information */
typedef struct Multaggr MULTAGGR;       /**< multiple aggregation information */
typedef struct Var VAR;                 /**< variable of the problem */



#include "def.h"
#include "memory.h"
#include "retcode.h"
#include "lp.h"
#include "stat.h"
#include "tree.h"
#include "prob.h"
#include "branch.h"
#include "event.h"


/** tracks changes of the variable's domains (fixed sized arrays) */
struct DomChg
{
   BOUNDCHG*        boundchg;           /**< array with changes in bounds of variables */
   HOLECHG*         holechg;            /**< array with changes in hole lists */
   int              nboundchg;          /**< number of bound changes */
   int              nholechg;           /**< number of hole list changes */
};

/** domain of a variable */
struct Dom
{
   HOLELIST*        holelist;           /**< list of holes */
   Real             lb;                 /**< lower bounds of variables */
   Real             ub;                 /**< upper bounds of variables */
};

/** aggregation information: $x = a*y + c$ */
struct Aggregate
{
   VAR*             var;                /**< variable $y$ in aggregation */
   Real             scalar;             /**< multiplier $a$ in aggregation */
   Real             constant;           /**< constant shift $c$ in aggregation */
};

/** multiple aggregation information: $x = a_1*y_1 + ... + a_k*y_k + c$ */
struct Multaggr
{
   VAR**            vars;               /**< variables $y$ in multiple aggregation */
   Real*            scalars;            /**< multipliers $a$ in multiple aggregation */
   Real             constant;           /**< constant shift $c$ in multiple aggregation */
   int              nvars;              /**< number of variables in aggregation */
};

/** variable of the problem */
struct Var
{
   union
   {
      VAR*          transvar;           /**< pointer to representing transformed variable (for original variables) */
      COL*          col;                /**< LP column (for column variables) */
      AGGREGATE     aggregate;          /**< aggregation information (for aggregated variables) */
      MULTAGGR      multaggr;           /**< multiple aggregation information (for multiple aggregated variables) */
   } data;
   char*            name;               /**< name of the variable */
   VAR**            parentvars;         /**< parent variables in the aggregation tree */
   EVENTFILTER*     eventfilter;        /**< event filter for events concerning this variable; not for ORIGINAL vars */
   DOM              dom;                /**< domain of variable */
   Real             obj;                /**< objective function value of variable */
   int              index;              /**< consecutively numbered variable identifier */
   int              probindex;          /**< array position in problems vars array, or -1 if not assigned to a problem */
   int              parentvarssize;     /**< available slots in parentvars array */
   int              nparentvars;        /**< number of parent variables in aggregation tree (used slots of parentvars) */
   int              pseudocandindex;    /**< array position in pseudo branching candidates array, or -1 */
   int              eventqueueindexlb;  /**< array position in event queue of lower bound change event, or -1 */
   int              eventqueueindexub;  /**< array position in event queue of upper bound change event, or -1 */
   int              nuses;              /**< number of times, this variable is referenced */
   int              nlocksdown;         /**< number of locks for rounding down; if zero, rounding down is always feasible */
   int              nlocksup;           /**< number of locks for rounding up; if zero, rounding up is always feasible */
   unsigned int     vartype:2;          /**< type of variable: binary, integer, implicit integer, continous */
   unsigned int     varstatus:3;        /**< status of variable: original, transformed, column, fixed, aggregated */
};



/*
 * domain change methods
 */

/** frees fixed size domain change data */
extern
RETCODE SCIPdomchgFree(
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   );

/** applies domain change */
extern
RETCODE SCIPdomchgApply(
   const DOMCHG*    domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** undoes domain change */
extern
RETCODE SCIPdomchgUndo(
   const DOMCHG*    domchg,             /**< domain change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );


/*
 * dynamic size attachment methods
 */

/** creates a dynamic size attachment for a domain change data structure */
extern
RETCODE SCIPdomchgdynCreate(
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   );

/** frees a dynamic size attachment for a domain change data structure */
extern
void SCIPdomchgdynFree(
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   );

/** attaches dynamic size information to domain change data */
extern
void SCIPdomchgdynAttach(
   DOMCHGDYN*       domchgdyn,          /**< dynamic size information */
   DOMCHG**         domchg              /**< pointer to static domain change */
   );

/** detaches dynamic size information and shrinks domain change data */
extern
RETCODE SCIPdomchgdynDetach(
   DOMCHGDYN*       domchgdyn,          /**< dynamic size information */
   MEMHDR*          memhdr              /**< block memory */
   );

/** frees attached domain change data and detaches dynamic size attachment */
extern
RETCODE SCIPdomchgdynDiscard(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr              /**< block memory */
   );

/** adds bound change to domain changes */
extern
RETCODE SCIPdomchgdynAddBoundchg(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** adds hole change to domain changes */
extern
RETCODE SCIPdomchgdynAddHolechg(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   );

/** gets pointer to domain change data the dynamic size information references */
extern
DOMCHG** SCIPdomchgdynGetDomchgPtr(
   DOMCHGDYN*       domchgdyn           /**< dynamically sized domain change data structure */
   );



/*
 * methods for variables 
 */

/** creates and captures an original problem variable */
extern
RETCODE SCIPvarCreate(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

/** creates and captures a loose variable belonging to the transformed problem */
extern
RETCODE SCIPvarCreateTransformed(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

/** increases usage counter of variable */
extern
void SCIPvarCapture(
   VAR*             var                 /**< variable */
   );

/** decreases usage counter of variable, and frees memory if necessary */
extern
RETCODE SCIPvarRelease(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (may be NULL, if it's not a column variable) */
   );

/** adds a hole to the variables domain */
extern
RETCODE SCIPvarAddHole(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   );

/** copies original variable into loose transformed variable, that is captured */
extern
RETCODE SCIPvarTransform(
   VAR**            transvar,           /**< pointer to store the transformed variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR*             origvar             /**< original problem variable */
   );

/** converts transformed variable into column variable and creates LP column */
extern
RETCODE SCIPvarColumn(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat                /**< problem statistics */
   );

/** converts variable into fixed variable, updates LP respectively */
extern
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   Real             fixedval            /**< value to fix variable at */
   );

/** converts variable into aggregated variable, updates LP respectively */
extern
RETCODE SCIPvarAggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   VAR*             aggvar,             /**< variable $y$ in aggregation $x = a*y + c$ */
   Real             scalar,             /**< multiplier $a$ in aggregation $x = a*y + c$ */
   Real             constant            /**< constant shift $c$ in aggregation $x = a*y + c$ */
   );

/** changes type of variable; cannot be called, if var belongs to a problem */
extern
RETCODE SCIPvarChgType(
   VAR*             var,                /**< problem variable to change */
   VARTYPE          vartype             /**< new type of variable */
   );

/**< increases lock number for rounding down; tells variable, that rounding its value down will make the solution
 *   infeasible
 */
extern
void SCIPvarForbidRoundDown(
   VAR*             var                 /**< problem variable */
   );

/**< increases lock number for rounding up; tells variable, that rounding its value up will make the solution infeasible */
extern
void SCIPvarForbidRoundUp(
   VAR*             var                 /**< problem variable */
   );

/**< increases lock number for rounding down and up; tells variable, that rounding value in either direction will make
 *   the solution infeasible
 */
extern
void SCIPvarForbidRound(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding down; cancels a prior forbidRoundDown() */
extern
void SCIPvarAllowRoundDown(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding up; cancels a prior forbidRoundUp() */
extern
void SCIPvarAllowRoundUp(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding down & up; cancels a prior forbidRound() */
extern
void SCIPvarAllowRound(
   VAR*             var                 /**< problem variable */
   );

/** gets number of locks for rounding down */
extern
int SCIPvarGetNLocksDown(
   VAR*             var                 /**< problem variable */
   );

/** gets number of locks for rounding up */
extern
int SCIPvarGetNLocksUp(
   VAR*             var                 /**< problem variable */
   );

/** is it possible, to round variable down and stay feasible? */
extern
Bool SCIPvarMayRoundDown(
   VAR*             var                 /**< problem variable */
   );

/** is it possible, to round variable up and stay feasible? */
extern
Bool SCIPvarMayRoundUp(
   VAR*             var                 /**< problem variable */
   );

/** changes lower bound of variable */
extern
RETCODE SCIPvarChgLb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   );

/** changes upper bound of variable */
extern
RETCODE SCIPvarChgUb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   );

/** changes bound of variable */
extern
RETCODE SCIPvarChgBd(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** adjust lower bound to integral value, if variable is integral */
extern
void SCIPvarAdjustLb(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            lb                  /**< pointer to lower bound to adjust */
   );

/** adjust upper bound to integral value, if variable is integral */
extern
void SCIPvarAdjustUb(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            ub                  /**< pointer to upper bound to adjust */
   );

/** changes objective value of variable */
extern
RETCODE SCIPvarChgObj(
   VAR*             var,                /**< variable to change, must not be member of the problem */
   Real             newobj              /**< new objective value for variable */
   );

/** get name of variable */
extern
const char* SCIPvarGetName(
   VAR*             var                 /**< problem variable */
   );

/** gets status of variable */
extern
VARSTATUS SCIPvarGetStatus(
   VAR*             var                 /**< problem variable */
   );

/** gets type of variable */
extern
VARTYPE SCIPvarGetType(
   VAR*             var                 /**< problem variable */
   );

/** gets unique index of variable */
extern
int SCIPvarGetIndex(
   VAR*             var                 /**< problem variable */
   );

/** gets position of variable in problem, or -1 if variable is not active */
extern
int SCIPvarGetProbIndex(
   VAR*             var                 /**< problem variable */
   );

/** gets corresponding transformed variable of an original variable */
extern
VAR* SCIPvarGetTransformed(
   VAR*             var                 /**< problem variable */
   );

/** gets column of COLUMN variable */
extern
COL* SCIPvarGetCol(
   VAR*             var                 /**< problem variable */
   );

/** gets objective function value of variable */
extern
Real SCIPvarGetObj(
   VAR*             var                 /**< problem variable */
   );

/** gets lower bound of variable */
extern
Real SCIPvarGetLb(
   VAR*             var                 /**< problem variable */
   );

/** gets upper bound of variable */
extern
Real SCIPvarGetUb(
   VAR*             var                 /**< problem variable */
   );

/** gets primal LP solution value of variable */
extern
Real SCIPvarGetLPSol(
   VAR*             var                 /**< problem variable */
   );

/** gets pseudo solution value of variable at actual node */
extern
Real SCIPvarGetPseudoSol(
   VAR*             var                 /**< problem variable */
   );

/** gets solution value of variable at actual node: if LP was solved at the node, the method returns the LP primal
 *  solution value, otherwise the pseudo solution
 */
extern
Real SCIPvarGetSol(
   VAR*             var,                /**< problem variable */
   TREE*            tree                /**< branch-and-bound tree */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
RETCODE SCIPvarAddToRow(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

/** includes event handler in variable's event filter */
extern
RETCODE SCIPvarCatchEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );


/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
extern
DECL_HASHGETKEY(SCIPhashGetKeyVar);



#endif
