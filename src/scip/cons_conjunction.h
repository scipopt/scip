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
#pragma ident "@(#) $Id: cons_conjunction.h,v 1.6 2005/05/31 17:20:11 bzfpfend Exp $"

/**@file   cons_conjunction.h
 * @brief  constraint handler for conjunction constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_CONJUNCTION_H__
#define __CONS_CONJUNCTION_H__


#include "scip/scip.h"


/** creates the handler for conjunction constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConshdlrConjunction(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a conjunction constraint */
extern
RETCODE SCIPcreateConsConjunction(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nconss,             /**< number of initial constraints in conjunction */
   CONS**           conss,              /**< initial constraint in conjunction */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             dynamic             /**< is constraint subject to aging? */
   );

/** adds constraint to the conjunction of constraints */
extern
RETCODE SCIPaddConsElemConjunction(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< conjunction constraint */
   CONS*            addcons             /**< additional constraint in conjunction */
   );

#endif
