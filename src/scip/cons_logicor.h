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
#pragma ident "@(#) $Id: cons_logicor.h,v 1.12 2004/11/02 11:27:35 bzfpfend Exp $"

/**@file   cons_logicor.h
 * @brief  constraint handler for logicor constraints
 *         (equivalent to set covering, but algorithms are suited for depth first search)
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_LOGICOR_H__
#define __CONS_LOGICOR_H__


#include "scip.h"


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
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** gets the dual solution of the logic or constraint in the current LP */
extern
Real SCIPgetDualsolLogicor(
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
