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
#pragma ident "@(#) $Id: implics.h,v 1.1 2005/08/08 13:20:35 bzfpfend Exp $"

/**@file   implics.h
 * @brief  methods for implications, variable bounds, and clique tables
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
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_implics.h"

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


#endif
