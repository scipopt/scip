/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: implics.h,v 1.4 2005/08/12 11:06:21 bzfpfend Exp $"

/**@file   implics.h
 * @brief  methods for implications, variable bounds, and cliques
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_IMPLICS_H__
#define __SCIP_IMPLICS_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_implics.h"
#include "scip/type_branch.h"
#include "scip/pub_implics.h"

#ifdef NDEBUG
#include "scip/struct_implics.h"
#endif


/*
 * Methods for Variable Bounds
 */

/** frees a variable bounds data structure */
extern
void SCIPvboundsFree(
   VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BLKMEM*          blkmem              /**< block memory */
   );

/** adds a variable bound to the variable bounds data structure */
extern
RETCODE SCIPvboundsAdd(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   BOUNDTYPE        vboundtype,         /**< type of variable bound (LOWER or UPPER) */
   VAR*             var,                /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   Real             coef,               /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   Real             constant            /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   );

/** removes from variable x a variable bound x >=/<= b*z + d with binary or integer z */
extern
RETCODE SCIPvboundsDel(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             vbdvar              /**< variable z    in x >=/<= b*z + d */
   );

/** reduces the number of variable bounds stored in the given variable bounds data structure */
void SCIPvboundsShrink(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   int              newnvbds            /**< new number of variable bounds */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets number of variable bounds contained in given variable bounds data structure */
extern
int SCIPvboundsGetNVbds(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

/** gets array of variables contained in given variable bounds data structure */
extern
VAR** SCIPvboundsGetVars(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

/** gets array of coefficients contained in given variable bounds data structure */
extern
Real* SCIPvboundsGetCoefs(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

/** gets array of constants contained in given variable bounds data structure */
extern
Real* SCIPvboundsGetConstants(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPvboundsGetNVbds(vbounds)     ((vbounds)->len)
#define SCIPvboundsGetVars(vbounds)      ((vbounds)->vars)
#define SCIPvboundsGetCoefs(vbounds)     ((vbounds)->coefs)
#define SCIPvboundsGetConstants(vbounds) ((vbounds)->constants)

#endif




/*
 * Methods for Implications
 */

/** frees an implications data structure */
extern
void SCIPimplicsFree(
   IMPLICS**        implics,            /**< pointer of implications data structure to free */
   BLKMEM*          blkmem              /**< block memory */
   );

/** adds an implication x == 0/1 -> y <= b or y >= b to the implications data structure;
 *  the implication must be non-redundant
 */
extern
RETCODE SCIPimplicsAdd(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             varfixing,          /**< FALSE if implication for x == 0 has to be added, TRUE for x == 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   Bool*            conflict            /**< pointer to store whether implication causes a conflict for variable x */
   );

/** removes the implication  x <= 0 or x >= 1  ==>  y <= b  or  y >= b  from the implications data structure */
extern
RETCODE SCIPimplicsDel(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             varfixing,          /**< FALSE if y should be removed from implications for x <= 0, TRUE for x >= 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype            /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets number of implications for a given binary variable fixing */
extern
int SCIPimplicsGetNImpls(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets number of implications on binary variables for a given binary variable fixing */
extern
int SCIPimplicsGetNBinImpls(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with implied variables for a given binary variable fixing */
extern
VAR** SCIPimplicsGetVars(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with implication types for a given binary variable fixing */
extern
BOUNDTYPE* SCIPimplicsGetTypes(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with implication bounds for a given binary variable fixing */
extern
Real* SCIPimplicsGetBounds(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

/** gets array with unique implication identifiers for a given binary variable fixing */
extern
int* SCIPimplicsGetIds(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPimplicsGetNImpls(implics, varfixing)       ((implics)->nimpls[varfixing])
#define SCIPimplicsGetNBinImpls(implics, varfixing)    ((implics)->nbinimpls[varfixing])
#define SCIPimplicsGetVars(implics, varfixing)         ((implics)->vars[varfixing])
#define SCIPimplicsGetTypes(implics, varfixing)        ((implics)->types[varfixing])
#define SCIPimplicsGetBounds(implics, varfixing)       ((implics)->bounds[varfixing])
#define SCIPimplicsGetIds(implics, varfixing)          ((implics)->ids[varfixing])

#endif




/*
 * methods for cliques
 */

/** adds a single variable to the given clique */
extern
RETCODE SCIPcliqueAddVar(
   CLIQUE*          clique,             /**< clique data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to the clique */
   Bool             value,              /**< value of the variable in the clique */
   Bool*            doubleentry,        /**< pointer to store whether the variable and value occurs twice in the clique */
   Bool*            oppositeentry       /**< pointer to store whether the variable with opposite value is in the clique */
   );

/** removes a single variable from the given clique */
extern
RETCODE SCIPcliqueDelVar(
   CLIQUE*          clique,             /**< clique data structure */
   VAR*             var,                /**< variable to remove from the clique */
   Bool             value               /**< value of the variable in the clique */
   );

/** frees a clique list data structure */
extern
void SCIPcliquelistFree(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem              /**< block memory */
   );

/** adds a clique to the clique list */
extern
RETCODE SCIPcliquelistAdd(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             value,              /**< value of the variable for which the clique list should be extended */
   CLIQUE*          clique              /**< clique that should be added to the clique list */
   );

/** removes a clique from the clique list */
extern
RETCODE SCIPcliquelistDel(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem,             /**< block memory */
   Bool             value,              /**< value of the variable for which the clique list should be reduced */
   CLIQUE*          clique              /**< clique that should be deleted from the clique list */
   );

/** removes all listed entries from the cliques */
extern
void SCIPcliquelistRemoveFromCliques(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   VAR*             var                 /**< active problem variable the clique list belongs to */
   );

/** creates a clique table data structure */
extern
RETCODE SCIPcliquetableCreate(
   CLIQUETABLE**    cliquetable         /**< pointer to store clique table data structure */
   );

/** frees a clique table data structure */
extern
RETCODE SCIPcliquetableFree(
   CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   BLKMEM*          blkmem              /**< block memory */
   );

/** adds a clique to the clique table; performs implications if the clique contains the same variable twice */
extern
RETCODE SCIPcliquetableAdd(
   CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int              nvars,              /**< number of variables in the clique */
   Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*             nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   );

/** removes all empty and single variable cliques from the clique table, and converts all two variable cliques
 *  into implications
 */
extern
RETCODE SCIPcliquetableCleanup(
   CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the number of cliques stored in the clique list */
extern
int SCIPcliquelistGetNCliques(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   Bool             value               /**< value of the variable for which the cliques should be returned */
   );

/** returns the cliques stored in the clique list, or NULL if the clique list is empty */
extern
CLIQUE** SCIPcliquelistGetCliques(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   Bool             value               /**< value of the variable for which the cliques should be returned */
   );

/** checks whether variable is contained in all cliques of the cliquelist */
extern
void SCIPcliquelistCheck(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   VAR*             var                 /**< variable, the clique list belongs to */
   );

/** gets the number of cliques stored in the clique table */
extern
int SCIPcliquetableGetNCliques(
   CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcliquelistGetNCliques(cliquelist, value) ((cliquelist) != NULL ? (cliquelist)->ncliques[value] : 0)
#define SCIPcliquelistGetCliques(cliquelist, value)  ((cliquelist) != NULL ? (cliquelist)->cliques[value] : NULL)
#define SCIPcliquelistCheck(cliquelist, var)         /**/
#define SCIPcliquetableGetNCliques(cliquetable)      ((cliquetable)->ncliques)

#endif


#endif
