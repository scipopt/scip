/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_linear.c
 * @brief  miscellaneous methods for linear constraints
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/scip.h"
#include "scip/pub_misc_nonlinear.h"
#include "scip/scipdefplugins.h"


/** returns the right-hand side of an arbitrary SCIP constraint that can be represented as a single nonlinear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_Real SCIPconsNonlinearGetRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which right-hand side is queried */
   SCIP_Bool*            success             /**< pointer to store whether a valid right-hand side was returned */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(success != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   *success = TRUE;
   rhs = SCIP_INVALID;

   if( strcmp(conshdlrname, "nonlinear") == 0 )
   {
      rhs = SCIPgetRhsNonlinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "quadratic") == 0 )
   {
      rhs = SCIPgetRhsQuadratic(scip, cons);
   }
   else if( strcmp(conshdlrname, "abspower") == 0 )
   {
      rhs = SCIPgetRhsAbspower(scip, cons);
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return rhs for constraint of type <%s>\n", conshdlrname);
      *success = FALSE;
   }

   return rhs;
}

/** returns the left-hand side of an arbitrary SCIP constraint that can be represented as a single nonlinear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_Real SCIPconsNonlinearGetLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get left-hand side for */
   SCIP_Bool*            success             /**< pointer to store whether a valid left-hand side was returned */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_Real lhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(success != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   *success = TRUE;
   lhs = SCIP_INVALID;

   if( strcmp(conshdlrname, "nonlinear") == 0 )
   {
      lhs = SCIPgetLhsNonlinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "quadratic") == 0 )
   {
      lhs = SCIPgetLhsQuadratic(scip, cons);
   }
   else if( strcmp(conshdlrname, "abspower") == 0 )
   {
      lhs = SCIPgetLhsAbspower(scip, cons);
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return lhs for constraint of type <%s>\n", conshdlrname);
      *success = FALSE;
   }

   return lhs;
}

/** adds the given variable to the input constraint. */
SCIP_RETCODE SCIPconsNonlinearAddLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which row is queried */
   SCIP_VAR*             var,                /**< variable of the constraint entry */
   SCIP_Real             val                 /**< the coefficient of the constraint entry */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "nonlinear") == 0 )
   {
      SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, val) );
   }
   else if( strcmp(conshdlrname, "quadratic") == 0 )
   {
      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, var, val) );
   }
   else if( strcmp(conshdlrname, "abspower") == 0 )
   {
      SCIPerrorMessage("Sorry, can't add coefficient for constraint of type <%s>\n", conshdlrname);
      return SCIP_ERROR;
   }
   else
   {
      SCIPerrorMessage("Sorry, can't add coefficient for constraint of type <%s>\n", conshdlrname);
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
