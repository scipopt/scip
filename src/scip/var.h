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
   SCIP_VARTYPE_IMPLINT   = 2,          /**< implicit integer variable: continous variable, that is allways integral */
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

extern
void SCIPdomchgFree(                    /**< frees fixed size domain change data */
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPdomchgApply(                /**< applies domain change */
   const DOMCHG*    domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

extern
RETCODE SCIPdomchgUndo(                 /**< undoes domain change */
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

extern
RETCODE SCIPdomchgdynCreate(            /**< creates a dynamic size attachment for a domain change data structure */
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
void SCIPdomchgdynFree(                 /**< frees a dynamic size attachment for a domain change data structure */
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
void SCIPdomchgdynAttach(               /**< attaches dynamic size information to domain change data */
   DOMCHGDYN*       domchgdyn,          /**< dynamic size information */
   DOMCHG**         domchg              /**< pointer to static domain change */
   );

extern
RETCODE SCIPdomchgdynDetach(            /**< detaches dynamic size information and shrinks domain change data */
   DOMCHGDYN*       domchgdyn,          /**< dynamic size information */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
void SCIPdomchgdynDiscard(              /**< frees attached domain change data and detaches dynamic size attachment */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   );

extern
DOMCHG** SCIPdomchgdynGetDomchgPtr(     /**< gets pointer to domain change data the dynamic size information references */
   DOMCHGDYN*       domchgdyn           /**< dynamically sized domain change data structure */
   );



/*
 * methods for variables 
 */

extern
RETCODE SCIPvarCreate(                  /**< creates and captures an original problem variable */
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

extern
RETCODE SCIPvarCreateTransformed(       /**< creates and captures a loose variable belonging to the transformed problem */
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

extern
void SCIPvarCapture(                    /**< increases usage counter of variable */
   VAR*             var                 /**< variable */
   );

extern
RETCODE SCIPvarRelease(                 /**< decreases usage counter of variable, and frees memory if necessary */
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (may be NULL, if it's not a column variable) */
   );

extern
RETCODE SCIPvarAddHole(                 /**< adds a hole to the variables domain */
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   );

extern
RETCODE SCIPvarTransform(               /**< copies original variable into loose transformed variable, that is captured */
   VAR**            transvar,           /**< pointer to store the transformed variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR*             origvar             /**< original problem variable */
   );

extern
RETCODE SCIPvarColumn(                  /**< converts transformed variable into column variable and creates LP column */
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat                /**< problem statistics */
   );

extern
RETCODE SCIPvarFix(                     /**< converts variable into fixed variable, updates LP respectively */
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   Real             fixedval            /**< value to fix variable at */
   );

extern
RETCODE SCIPvarAggregate(               /**< converts variable into aggregated variable, updates LP respectively */
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   VAR*             aggvar,             /**< variable $y$ in aggregation $x = a*y + c$ */
   Real             scalar,             /**< multiplier $a$ in aggregation $x = a*y + c$ */
   Real             constant            /**< constant shift $c$ in aggregation $x = a*y + c$ */
   );

extern
RETCODE SCIPvarChgType(                 /**< changes type of variable; cannot be called, if var belongs to a problem */
   VAR*             var,                /**< problem variable to change */
   VARTYPE          vartype             /**< new type of variable */
   );

extern
void SCIPvarForbidRoundDown(            /**< increases lock number for rounding down; tells variable, that rounding its
                                         *   value down will make the solution infeasible */
   VAR*             var                 /**< problem variable */
   );

extern
void SCIPvarForbidRoundUp(              /**< increases lock number for rounding up; tells variable, that rounding its
                                         *   value up will make the solution infeasible */
   VAR*             var                 /**< problem variable */
   );

extern
void SCIPvarForbidRound(                /**< increases lock number for rounding down and up; tells variable, that rounding
                                         *   value in either direction will make the solution infeasible */
   VAR*             var                 /**< problem variable */
   );

extern
void SCIPvarAllowRoundDown(             /**< decreases lock number for rounding down; cancels a prior forbidRoundDown() */
   VAR*             var                 /**< problem variable */
   );

extern
void SCIPvarAllowRoundUp(               /**< decreases lock number for rounding up; cancels a prior forbidRoundUp() */
   VAR*             var                 /**< problem variable */
   );

extern
void SCIPvarAllowRound(                 /**< decreases lock number for rounding down & up; cancels a prior forbidRound() */
   VAR*             var                 /**< problem variable */
   );

extern
int SCIPvarGetNLocksDown(               /**< gets number of locks for rounding down */
   VAR*             var                 /**< problem variable */
   );

extern
int SCIPvarGetNLocksUp(                 /**< gets number of locks for rounding up */
   VAR*             var                 /**< problem variable */
   );

extern
Bool SCIPvarMayRoundDown(               /**< is it possible, to round variable down and stay feasible? */
   VAR*             var                 /**< problem variable */
   );

extern
Bool SCIPvarMayRoundUp(                 /**< is it possible, to round variable up and stay feasible? */
   VAR*             var                 /**< problem variable */
   );

extern
RETCODE SCIPvarChgLb(                   /**< changes lower bound of variable */
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

extern
RETCODE SCIPvarChgUb(                   /**< changes upper bound of variable */
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

extern
RETCODE SCIPvarChgBd(                   /**< changes bound of variable */
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

extern
void SCIPvarAdjustLb(                   /**< adjust lower bound to integral value, if variable is integral */
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            lb                  /**< pointer to lower bound to adjust */
   );

extern
void SCIPvarAdjustUb(                   /**< adjust upper bound to integral value, if variable is integral */
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            ub                  /**< pointer to upper bound to adjust */
   );

extern
RETCODE SCIPvarChgObj(                  /**< changes objective value of variable */
   VAR*             var,                /**< variable to change, must not be member of the problem */
   Real             newobj              /**< new objective value for variable */
   );

extern
const char* SCIPvarGetName(             /**< get name of variable */
   VAR*             var                 /**< problem variable */
   );

extern
VARSTATUS SCIPvarGetStatus(             /**< gets status of variable */
   VAR*             var                 /**< problem variable */
   );

extern
VARTYPE SCIPvarGetType(                 /**< gets type of variable */
   VAR*             var                 /**< problem variable */
   );

extern
int SCIPvarGetIndex(                    /**< gets unique index of variable */
   VAR*             var                 /**< problem variable */
   );

extern
int SCIPvarGetProbIndex(                /**< gets position of variable in problem, or -1 if variable is not active */
   VAR*             var                 /**< problem variable */
   );

extern
VAR* SCIPvarGetTransformed(             /**< gets corresponding transformed variable of an original variable */
   VAR*             var                 /**< problem variable */
   );

extern
COL* SCIPvarGetCol(                     /**< gets column of COLUMN variable */
   VAR*             var                 /**< problem variable */
   );

extern
Real SCIPvarGetObj(                     /**< gets objective function value of variable */
   VAR*             var                 /**< problem variable */
   );

extern
Real SCIPvarGetLb(                      /**< gets lower bound of variable */
   VAR*             var                 /**< problem variable */
   );

extern
Real SCIPvarGetUb(                      /**< gets upper bound of variable */
   VAR*             var                 /**< problem variable */
   );

extern
Real SCIPvarGetLPSol(                   /**< get primal LP solution value of variable */
   VAR*             var                 /**< problem variable */
   );

extern
Real SCIPvarGetPseudoSol(               /**< get pseudo solution value of variable at actual node */
   VAR*             var                 /**< problem variable */
   );

extern
Real SCIPvarGetSol(                     /**< get solution value of variable at actual node: if LP was solved at the node,
                                           the method returns the LP primal solution value, otherwise the pseudo solution */
   VAR*             var,                /**< problem variable */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPvarAddToRow(                /**< resolves variable to columns and adds them with the coefficient to the row */
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

extern
RETCODE SCIPvarCatchEvent(              /**< includes event handler in variable's event filter */
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

extern
DECL_HASHGETKEY(SCIPhashGetKeyVar);     /**< gets the key (i.e. the name) of the given variable */



#endif
