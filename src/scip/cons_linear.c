/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_linear.c
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_linear.h"


DECL_CONSINIT(SCIPconsInit_Linear)
{
   assert(self != NULL);
   assert(strcmp(self->name, CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(memhdr != NULL);

   return SCIP_OKAY;
}

DECL_CONSEXIT(SCIPconsExit_Linear)
{
   assert(self != NULL);
   assert(strcmp(self->name, CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(memhdr != NULL);

   return SCIP_OKAY;
}

DECL_CONSFREE(SCIPconsFree_Linear)
{
   assert(self != NULL);
   assert(strcmp(self->name, CONSHDLR_NAME) == 0);
   assert(memhdr != NULL);

   return SCIP_OKAY;
}

DECL_CONSCHCK(SCIPconsChck_Linear)
{
   assert(self != NULL);
   assert(strcmp(self->name, CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);
   assert(psol != NULL);

   return SCIP_OKAY;
}

DECL_CONSPROP(SCIPconsProp_Linear)
{
   assert(self != NULL);
   assert(strcmp(self->name, CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);

   return SCIP_OKAY;
}


RETCODE SCIPconsCreate_Linear(          /**< creates a linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             lhs,                /**< left hand side of row (for ranged rows) */
   Real             epsilon,            /**< maximal normed violation of row */
   ROWTYPE          rowtype,            /**< type of row */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model               /**< is constraint necessary for feasibility? */
   )
{
   ROW* row;

   CHECK_OKAY( SCIPcreateRow(scip, &row, name, len, col, val, rhs, lhs, epsilon, rowtype) );

   CHECK_OKAY( SCIPcaptureRow(scip, row) );

   CHECK_OKAY( SCIPcreateConstraint(scip, cons, &ConsHdlr_Linear, (CONSDATA*)row, original, model) );

   return SCIP_OKAY;
}
