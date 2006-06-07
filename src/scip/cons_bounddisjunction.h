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
#pragma ident "@(#) $Id: cons_bounddisjunction.h,v 1.2 2006/06/07 08:21:00 bzfpfend Exp $"

/**@file   cons_bounddisjunction.h
 * @brief  constraint handler for bound disjunction constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_BOUNDDISJUNCTION_H__
#define __SCIP_CONS_BOUNDDISJUNCTION_H__


#include "scip/scip.h"


/** creates the handler for bound disjunction constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrBounddisjunction(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a bound disjunction constraint */
extern
SCIP_RETCODE SCIPcreateConsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< variables of the literals in the constraint */
   SCIP_BOUNDTYPE*       boundtypes,         /**< types of bounds of the literals (lower or upper bounds) */
   SCIP_Real*            bounds,             /**< bounds of the literals */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable during node processing (subject to col generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup? */
   );

/** gets number of variables in bound disjunction constraint */
extern
int SCIPgetNVarsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of variables in bound disjunction constraint */
extern
SCIP_VAR** SCIPgetVarsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of bound types in bound disjunction constraint */
extern
SCIP_BOUNDTYPE* SCIPgetBoundtypesBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of bounds in bound disjunction constraint */
extern
SCIP_Real* SCIPgetBoundsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );


#endif
