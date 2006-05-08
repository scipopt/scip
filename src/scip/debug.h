/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: debug.h,v 1.19 2006/05/08 10:12:09 bzfpfend Exp $"

/**@file   debug.h
 * @brief  methods for debugging
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEBUG_H__
#define __SCIP_DEBUG_H__

/** uncomment this define to activate debugging on given solution */
/*#define SCIP_DEBUG_SOLUTION "modell.sol"*/


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_lp.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"


#ifdef SCIP_DEBUG_SOLUTION

/** checks whether given row is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< row to check for validity */
   );

/** checks whether given global lower bound is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckLbGlobal(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lb                  /**< lower bound */
   );

/** checks whether given global upper bound is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckUbGlobal(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             ub                  /**< upper bound */
   );

/** checks whether given local bound implication is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckInference(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< local node where this bound change was applied */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** informs solution debugger, that the given node will be freed */
extern
SCIP_RETCODE SCIPdebugRemoveNode(
   SCIP_NODE*            node                /**< node that will be freed */
   );

/** checks whether given variable bound is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckVbound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable x in x <= b*z + d  or  x >= b*z + d */
   SCIP_BOUNDTYPE        vbtype,             /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             vbvar,              /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbcoef,             /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbconstant          /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   );

/** checks whether given implication is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckImplic(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound           /**< bound b    in implication y <= b or y >= b */
   );

#if 0 /*?????????????????????? update this method to the general conflict constraints */
/** checks whether given conflict is valid for the debugging solution */
extern
SCIP_RETCODE SCIPdebugCheckConflict(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_VAR**            conflictset,        /**< variables in the conflict set */
   int                   nliterals           /**< number of literals in the conflict set */
   );
#endif

/** creates the debugging propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPdebugIncludeProp(
   SCIP*                 scip                /**< SCIP data structure */
   );

#else

#define SCIPdebugCheckRow(set,row) SCIP_OKAY
#define SCIPdebugCheckLbGlobal(set,var,lb) SCIP_OKAY
#define SCIPdebugCheckUbGlobal(set,var,ub) SCIP_OKAY
#define SCIPdebugCheckInference(blkmem,set,node,var,newbound,boundtype) SCIP_OKAY
#define SCIPdebugRemoveNode(node) SCIP_OKAY
#define SCIPdebugCheckVbound(set,var,vbtype,vbvar,vbcoef,vbconstant) SCIP_OKAY
#define SCIPdebugCheckImplic(set,var,varfixing,implvar,impltype,implbound) SCIP_OKAY
#define SCIPdebugCheckConflict(blkmem,set,node,conflictset,nliterals) SCIP_OKAY
#define SCIPdebugIncludeProp(scip) SCIP_OKAY

#endif


#endif
