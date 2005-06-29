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
#pragma ident "@(#) $Id: cons_logicor.h,v 1.18 2005/06/29 11:08:04 bzfpfend Exp $"

/**@file   cons_logicor.h
 * @brief  constraint handler for logicor constraints
 *         (equivalent to set covering, but algorithms are suited for depth first search)
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_LOGICOR_H__
#define __CONS_LOGICOR_H__


#include "scip/scip.h"


/** creates the handler for logic or constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConshdlrLogicor(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a logic or constraint */
extern
RETCODE SCIPcreateConsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable during node processing (subject to col generation)? */
   Bool             dynamic,            /**< is constraint subject to aging? */
   Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   );

/** gets the dual solution of the logic or constraint in the current LP */
extern
Real SCIPgetDualsolLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets the dual farkas value of the logic or constraint in the current infeasible LP */
extern
Real SCIPgetDualfarkasLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets array of variables in logic or constraint */
extern
VAR** SCIPgetVarsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets number of variables in logic or constraint */
extern
int SCIPgetNVarsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );


#endif
