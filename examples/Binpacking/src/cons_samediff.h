/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_samediff.h
 * @brief  Constraint handler stores the local branching decision data
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This constraint handler is used to store the branching decision of the \ref BINPACKING_BRANCHING "Ryan/Foster branching rule"
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
   SAME   = 1                                /**< constraint representing the branching decision same(i,j) */
};
typedef enum ConsType CONSTYPE;

/** creates the handler for element constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSamediff(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a samediff constraint */
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
int SCIPgetItemid1Samediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

/** returns item id two */
int SCIPgetItemid2Samediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

/** return constraint type SAME or DIFFER */
CONSTYPE SCIPgetTypeSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

#endif
