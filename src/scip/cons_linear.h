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

/**@file   cons_linear.h
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_LINEAR_H__
#define __CONS_LINEAR_H__


#include "scip.h"



extern
RETCODE SCIPincludeConsHdlrLinear(      /**< creates the handler for linear constraints and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPcreateConsLinear(           /**< creates and captures a linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             model,              /**< is constraint necessary for feasibility? */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   );

extern
RETCODE SCIPcreateConsLinearLPRow(      /**< creates and captures a linear constraint from an LP row, captures the row */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   ROW*             row,                /**< LP row */
   Bool             model               /**< is constraint necessary for feasibility? */
   );

extern
RETCODE SCIPaddCoefConsLinear(          /**< adds coefficient in linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   );

extern
RETCODE SCIPgetLhsConsLinear(           /**< gets left hand side of linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real*            lhs                 /**< pointer to store left hand side */
   );

extern
RETCODE SCIPgetRhsConsLinear(           /**< gets right hand side of linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real*            rhs                 /**< pointer to store right hand side */
   );

extern
RETCODE SCIPchgLhsConsLinear(           /**< changes left hand side of linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             lhs                 /**< new left hand side */
   );

extern
RETCODE SCIPchgRhsConsLinear(           /**< changes right hand side of linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             rhs                 /**< new right hand side */
   );


#endif
