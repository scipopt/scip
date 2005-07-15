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
#pragma ident "@(#) $Id: cons_setppc.h,v 1.21 2005/07/15 17:20:07 bzfpfend Exp $"

/**@file   cons_setppc.h
 * @brief  constraint handler for the set partitioning / packing / covering constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SETPPC_H__
#define __SCIP_CONS_SETPPC_H__


#include "scip/scip.h"


/** type of setppc constraint: set partitioning, set packing, or set covering */
enum SetppcType
{
   SCIP_SETPPCTYPE_PARTITIONING = 0,     /**< constraint is a set partitioning constraint: sum(x) == 1 */
   SCIP_SETPPCTYPE_PACKING      = 1,     /**< constraint is a set packing constraint:      sum(x) <= 1 */
   SCIP_SETPPCTYPE_COVERING     = 2      /**< constraint is a set covering constraint:     sum(x) >= 1 */
};
typedef enum SetppcType SETPPCTYPE;

/** creates the handler for set partitioning / packing / covering constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConshdlrSetppc(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a set partitioning constraint */
extern
RETCODE SCIPcreateConsSetpart(
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

/** creates and captures a set packing constraint */
extern
RETCODE SCIPcreateConsSetpack(
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

/** creates and captures a set covering constraint */
extern
RETCODE SCIPcreateConsSetcover(
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

/** adds coefficient in set partitioning / packing / covering constraint */
extern
RETCODE SCIPaddCoefSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var                 /**< variable to add to the constraint */
   );

/** gets the dual solution of the set partitioning / packing / covering constraint in the current LP */
extern
Real SCIPgetDualsolSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets the dual farkas value of the set partitioning / packing / covering constraint in the current infeasible LP */
extern
Real SCIPgetDualfarkasSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets array of variables in set partitioning / packing / covering constraint */
extern
VAR** SCIPgetVarsSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets number of variables in set partitioning / packing / covering constraint */
extern
int SCIPgetNVarsSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets type of set partitioning / packing / covering constraint */
extern
SETPPCTYPE SCIPgetTypeSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );


#endif
