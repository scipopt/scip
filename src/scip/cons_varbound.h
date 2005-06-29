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
#pragma ident "@(#) $Id: cons_varbound.h,v 1.9 2005/06/29 11:08:05 bzfpfend Exp $"

/**@file   cons_varbound.h
 * @brief  constraint handler for varbound constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_VARBOUND_H__
#define __CONS_VARBOUND_H__


#include "scip/scip.h"


/** creates the handler for varbound constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConshdlrVarbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a varbound constraint: lhs <= x + c*y <= rhs */
extern
RETCODE SCIPcreateConsVarbound(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   VAR*             var,                /**< variable x that has variable bound */
   VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   Real             vbdcoef,            /**< coefficient c of bounding variable y */
   Real             lhs,                /**< left hand side of variable bound inequality */
   Real             rhs,                /**< right hand side of variable bound inequality */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             dynamic,            /**< is constraint subject to aging? */
   Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   );

/** gets the dual solution of the varbound constraint in the current LP */
extern
Real SCIPgetDualsolVarbound(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets the dual farkas value of the varbound constraint in the current infeasible LP */
extern
Real SCIPgetDualfarkasVarbound(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

#endif
