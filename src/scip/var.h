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
   SCIP_VARSTATUS_AGGREGATED = 4        /**< variable is aggregated to $x = a*y + c$ in the transformed problem */
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

/** type of bound: lower or upper bound */
enum BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum BoundType BOUNDTYPE;

typedef struct DomChg DOMCHG;           /**< changes in domains of variables (fixed sized arrays) */
typedef struct DomChgDyn DOMCHGDYN;     /**< changes in domains of variables (dynamically sized arrays) */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct Dom DOM;                 /**< datastructures for storing domains of variables */
typedef struct Aggregate AGGREGATE;     /**< aggregation information */
typedef struct Var VAR;                 /**< variable of the problem */



#include "def.h"
#include "memory.h"
#include "retcode.h"
#include "lp.h"
#include "stat.h"


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

/** variable of the problem */
struct Var
{
   union
   {
      VAR*          transvar;           /**< pointer to representing transformed variable (for original variables) */
      COL*          col;                /**< LP column (for column variables) */
      AGGREGATE     aggregate;          /**< aggregation information (for aggregated variables) */
   } data;
   VAR*             origvar;            /**< pointer to original problem variable this var represents, or NULL */
   char*            name;               /**< name of the column */
   DOM              dom;                /**< domain of variable */
   Real             obj;                /**< objective function value of variable */
   unsigned int     vartype:2;          /**< type of variable: binary, integer, implicit integer, continous */
   unsigned int     varstatus:3;        /**< status of variable: original, transformed, column, fixed, aggregated */
};




extern
RETCODE SCIPdomchgdynCreate(            /**< creates a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn           /**< pointer to dynamically sized domain change data structure */
   );

extern
void SCIPdomchgdynFree(                 /**< frees a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn           /**< pointer to dynamically sized domain change data structure */
   );

extern
RETCODE SCIPdomchgdynCopy(              /**< copies data from fixed size domain change into dynamically sized one */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   DOMCHG*          domchg              /**< static domain change */
   );

extern
RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   );

extern
RETCODE SCIPdomchgCreate(               /**< creates domain change data (fixed size) from dynamically sized data */
   DOMCHG**         domchg,             /**< pointer to fixed size domain change data */
   MEMHDR*          memhdr,             /**< block memory */
   const DOMCHGDYN* domchgdyn           /**< dynamically sized domain change data structure */
   );

extern
void SCIPdomchgFree(                    /**< frees fixed size domain change data */
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPdomchgApply(                /**< applies domain change */
   const DOMCHG*    domchg,             /**< domain change to apply */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPdomchgUndo(                 /**< undoes domain change */
   const DOMCHG*    domchg,             /**< domain change to remove */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );



extern
RETCODE SCIPvarCreate(                  /**< creates an original problem variable */
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

extern
RETCODE SCIPvarCreateTransformed(       /**< creates a variable belonging only to the transformed problem */
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

extern
void SCIPvarFree(                       /**< frees a variable */
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
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
RETCODE SCIPvarTransform(               /**< copies original variable into transformed variable */
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   VAR**            transvar            /**< pointer to transformed variable */
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
RETCODE SCIPvarChgLb(                   /**< changes lower bound of variable */
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   );

extern
RETCODE SCIPvarChgUb(                   /**< changes upper bound of variable */
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   );



#endif
