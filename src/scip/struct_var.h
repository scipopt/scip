/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_var.h,v 1.6 2004/02/05 14:12:44 bzfpfend Exp $"

/**@file   struct_var.h
 * @brief  datastructures for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_VAR_H__
#define __STRUCT_VAR_H__


#include "def.h"
#include "type_history.h"
#include "type_event.h"
#include "type_var.h"
#include "type_cons.h"


/** hole in a domain */
struct Hole
{
   Real             left;               /**< left bound of open interval defining the hole (left,right) */
   Real             right;              /**< right bound of open interval defining the hole (left,right) */
};

/** list of domain holes */
struct Holelist
{
   HOLE             hole;               /**< this hole */
   HOLELIST*        next;               /**< next hole in list */
};

/** change in a hole list */
struct HoleChg
{
   HOLELIST**       ptr;                /**< changed list pointer */
   HOLELIST*        newlist;            /**< new value of list pointer */
   HOLELIST*        oldlist;            /**< old value of list pointer */
};

/** data for branching decision bound changes */
struct BranchingData
{
   Real             lpsolval;           /**< sol val of var in last LP prior to bound change, or SCIP_INVALID if unknown */
};

/** data for inferred bound changes */
struct InferenceData
{
   VAR*             var;                /**< variable that was changed (parent of var, or var itself) */
   CONS*            cons;               /**< constraint that inferred this bound change, or NULL */
};

/** change in one bound of a variable */
struct BoundChg
{
   VAR*             var;                /**< active variable to change the bounds for */
   union
   {
      BRANCHINGDATA branchingdata;      /**< data for branching decisions */
      INFERENCEDATA inferencedata;      /**< data for inferred bound changes */
   } data;
   Real             newbound;           /**< new value for bound */
   Real             oldbound;           /**< old value for bound */
   unsigned int     depth:16;           /**< depth in the tree, where this bound change took place */
   unsigned int     index:14;           /**< bound change index for each node representing the order of changes */
   unsigned int     boundtype:1;        /**< type of bound for var: lower or upper bound */
   unsigned int     boundchgtype:1;     /**< bound change type: branching decision or inferred bound change */
};

/** tracks changes of the variable's domains (static arrays, bound changes only) */
struct DomChgBound
{
   unsigned int     domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   unsigned int     nboundchgs:30;      /**< number of bound changes */
   BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
};

/** tracks changes of the variable's domains (static arrays, bound and hole changes) */
struct DomChgBoth
{
   unsigned int     domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   unsigned int     nboundchgs:30;      /**< number of bound changes */
   BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
   HOLECHG*         holechgs;           /**< array with changes in hole lists */
   int              nholechgs;          /**< number of hole list changes */
};

/** tracks changes of the variable's domains (dynamic arrays) */
struct DomChgDyn
{
   unsigned int     domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   unsigned int     nboundchgs:30;      /**< number of bound changes */
   BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
   HOLECHG*         holechgs;           /**< array with changes in hole lists */
   int              nholechgs;          /**< number of hole list changes */
   int              boundchgssize;      /**< size of bound changes array */
   int              holechgssize;       /**< size of hole changes array */
};

/** tracks changes of the variable's domains */
union DomChg
{
   DOMCHGBOUND      domchgbound;        /**< bound changes */
   DOMCHGBOTH       domchgboth;         /**< bound and hole changes */
   DOMCHGDYN        domchgdyn;          /**< bound and hole changes with dynamic arrays */
};

/** domain of a variable */
struct Dom
{
   HOLELIST*        holelist;           /**< list of holes */
   Real             lb;                 /**< lower bounds of variables */
   Real             ub;                 /**< upper bounds of variables */
};

/** aggregation information: x = a*y + c */
struct Aggregate
{
   VAR*             var;                /**< variable y in aggregation */
   Real             scalar;             /**< multiplier a in aggregation */
   Real             constant;           /**< constant shift c in aggregation */
};

/** multiple aggregation information: x = a_1*y_1 + ... + a_k*y_k + c */
struct Multaggr
{
   VAR**            vars;               /**< variables y in multiple aggregation */
   Real*            scalars;            /**< multipliers a in multiple aggregation */
   Real             constant;           /**< constant shift c in multiple aggregation */
   int              nvars;              /**< number of variables in aggregation */
   int              varssize;           /**< size of vars and scalars arrays */
};

/** negation information: x' = c - x */
struct Negate
{
   Real             constant;           /**< constant shift c in negation */
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
      NEGATE        negate;             /**< negation information (for negated variables) */
   } data;
   char*            name;               /**< name of the variable */
   VAR**            parentvars;         /**< parent variables in the aggregation tree */
   VAR*             negatedvar;         /**< pointer to the variables negation: x' = lb + ub - x, or NULL if not created */
   EVENTFILTER*     eventfilter;        /**< event filter for events concerning this variable; not for ORIGINAL vars */
   VAR*             infervar;           /**< variable that was assigned (parent of var, or var itself) */
   CONS*            infercons;          /**< constraint that deduced the assignment (binary variables only), or NULL */
   HISTORY*         lphistory;          /**< branching history information for downwards and upwards branching on LP */
   DOM              glbdom;             /**< domain of variable in global problem */
   DOM              locdom;             /**< domain of variable in current subproblem */
   Real             obj;                /**< objective function value of variable */
   Real             branchingpriority;  /**< priority of the variable to choose as branching variable */
   int              index;              /**< consecutively numbered variable identifier */
   int              probindex;          /**< array position in problems vars array, or -1 if not assigned to a problem */
   int              pseudocandindex;    /**< array position in pseudo branching candidates array, or -1 */
   int              eventqueueindexobj; /**< array position in event queue of objective change event, or -1 */
   int              eventqueueindexlb;  /**< array position in event queue of lower bound change event, or -1 */
   int              eventqueueindexub;  /**< array position in event queue of upper bound change event, or -1 */
   int              parentvarssize;     /**< available slots in parentvars array */
   int              nparentvars;        /**< number of parent variables in aggregation tree (used slots of parentvars) */
   int              nuses;              /**< number of times, this variable is referenced */
   int              nlocksdown;         /**< number of locks for rounding down; if zero, rounding down is always feasible */
   int              nlocksup;           /**< number of locks for rounding up; if zero, rounding up is always feasible */
   unsigned int     inferdepth:16;      /**< depth in the tree, where the last bound change took place */
   unsigned int     inferindex:15;      /**< bound change index for each node representing the order of changes */
   unsigned int     initial:1;          /**< TRUE iff var's column should be present in the initial root LP */
   unsigned int     removeable:1;       /**< TRUE iff var's column is removeable from the LP (due to aging or cleanup) */
   unsigned int     vartype:2;          /**< type of variable: binary, integer, implicit integer, continuous */
   unsigned int     varstatus:3;        /**< status of variable: original, transformed, column, fixed, aggregated */
   unsigned int     historyflag:2;      /**< temporary flag used in branching history update */
};


#endif
