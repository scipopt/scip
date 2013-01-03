/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_samediff.h
 * @brief  Constraint handler stores the local branching decision data
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This constraint handler is used to store the branching decision of the \ref BRANCHING "Ryan/Foster branching rule"
 * which is implemented in \ref branch_ryanfoster.c.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SAMEDIFF_H__
#define __SCIP_CONS_SAMEDIFF_H__


#include "scip/scip.h"

/* type of constraint: differ or same */
enum ConsType
{
   DIFFER = 0,                               /**< constraint representing the branching decision differ(i,j) */
   SAME   = 1,                               /**< constraint representing the branching decision same(i,j) */
};
typedef enum ConsType CONSTYPE;

/** creates the handler for element constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrSamediff(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a samediff constraint */
extern
SCIP_RETCODE SCIPcreateConsSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   itemid1,            /**< item id one */
   int                   itemid2,            /**< item id two */
   CONSTYPE              type,               /**< stores whether the items have to be in the SAME or DIFFER packing */
   SCIP_NODE*            node,               /**< the node in the B&B-tree at which the cons is sticking */
   SCIP_Bool             local               /**< is constraint only valid locally? */
   );

/** returns item id one */
extern
int SCIPgetItemid1Samediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

/** returns item id two */
extern
int SCIPgetItemid2Samediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

/** return constraint type SAME or DIFFER */
extern
CONSTYPE SCIPgetTypeSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

#endif
