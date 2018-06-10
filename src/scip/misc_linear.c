/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc_linear.c
 * @brief  miscellaneous methods for linear constraints
 * @author Jakob Witzig
 * @author Ambros Gleixner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/scip.h"
#include "scip/pub_misc_linear.h"
#include "scip/cons_setppc.h"
#include "scip/scipdefplugins.h"


/** returns the right-hand side of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_Real SCIPconsGetRhs(
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

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      rhs = SCIPgetRhsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch( SCIPgetTypeSetppc(scip, cons) )
      {
      case SCIP_SETPPCTYPE_PARTITIONING: /* fall through intended */
      case SCIP_SETPPCTYPE_PACKING:
         rhs = 1.0;
         break;

      case SCIP_SETPPCTYPE_COVERING:
         rhs = SCIPinfinity(scip);
         break;
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      rhs = SCIPinfinity(scip);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      rhs = SCIPgetCapacityKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      rhs = SCIPgetRhsVarbound(scip, cons);
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return rhs for constraint of type <%s>\n", conshdlrname);
      *success = FALSE;
   }

   return rhs;
}

/** returns the left-hand side of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_Real SCIPconsGetLhs(
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

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      lhs = SCIPgetLhsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch( SCIPgetTypeSetppc(scip, cons) )
      {
      case SCIP_SETPPCTYPE_PARTITIONING: /* fall through intended */
      case SCIP_SETPPCTYPE_COVERING:
         lhs = 1.0;
         break;

      case SCIP_SETPPCTYPE_PACKING:
         lhs = -SCIPinfinity(scip);
         break;
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      lhs = 1.0;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      lhs = -SCIPinfinity(scip);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      lhs = SCIPgetLhsVarbound(scip, cons);
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return lhs for constraint of type <%s>\n", conshdlrname);
      *success = FALSE;
   }

   return lhs;
}

/** returns the value array of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the involved values
 */
SCIP_RETCODE SCIPgetConsVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the coefficients are wanted */
   SCIP_Real*            vals,               /**< array to store the coefficients of the constraint */
   int                   varssize,           /**< available slots in vals array needed to check if the array is large enough */
   SCIP_Bool*            success             /**< pointer to store whether the coefficients are successfully copied */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vals != NULL);
   assert(success != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   conshdlrname = SCIPconshdlrGetName(conshdlr);

   *success = TRUE;

   SCIP_CALL( SCIPgetConsNVars(scip, cons, &nvars, success) );

   if( !(*success) )
   {
      SCIPwarningMessage(scip, "Cannot return value array for constraint of type\n", conshdlrname);
      return SCIP_OKAY;
   }

   if( varssize < nvars )
   {
      SCIPwarningMessage(scip, "Cannot return value array for constraint of type (insufficient memory provided)\n", conshdlrname);
      *success = FALSE;
      return SCIP_OKAY;
   }

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      SCIP_Real* linvals;

      linvals = SCIPgetValsLinear(scip, cons);
      assert(linvals != NULL);

      for( i = 0; i < nvars; i++ )
      {
         vals[i] = linvals[i];
      }
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      for( i = 0; i < nvars; i++ )
      {
         vals[i] = 1.0;
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      for( i = 0; i < nvars; i++ )
      {
         vals[i] = 1.0;
      }
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      SCIP_Longint* weights;

      weights = SCIPgetWeightsKnapsack(scip, cons);
      assert(weights != NULL);

      for( i = 0; i < nvars; i++ )
      {
         vals[i] = (SCIP_Real)weights[i];
      }
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      assert(nvars == 2);

      vals[0] = 1.0;
      vals[1] = SCIPgetVbdcoefVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIP_Real* weights;

      weights = SCIPgetWeightsSOS1(scip, cons);
      assert(weights != NULL);

      for( i = 0; i < nvars; i++ )
      {
         vals[i] = weights[i];
      }
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIP_Real* weights;

      weights = SCIPgetWeightsSOS2(scip, cons);
      assert(weights != NULL);

      for( i = 0; i < nvars; i++ )
      {
         vals[i] = weights[i];
      }
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return value array for constraint of type <%s>\n", conshdlrname);
      *success = FALSE;
   }

   return SCIP_OKAY;
}

/** returns the dual farkas sol of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the dual farkas solution
 */
void SCIPconsGetDualfarkas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get the dual farkas solution for */
   SCIP_Real*            dualfarkas,         /**< pointer to store the dual farkas solution */
   SCIP_Bool*            success             /**< pointer to store whether the dual farkas solution is successfully returned */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   *success = TRUE;

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      *dualfarkas = SCIPgetDualfarkasLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      *dualfarkas = SCIPgetDualfarkasSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      *dualfarkas = SCIPgetDualfarkasLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      *dualfarkas = SCIPgetDualfarkasKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      *dualfarkas = SCIPgetDualfarkasVarbound(scip, cons);
   }
   /* these are Benders' specific constraint handlers */
   else if( strcmp(conshdlrname, "origbranch") == 0 || strcmp(conshdlrname, "masterbranch") == 0 )
   {
      *dualfarkas = 0.0;
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return dual farkas solution for constraint of type <%s>\n", conshdlrname);
      *dualfarkas = 0.0;
      *success = FALSE;
   }
}

/** returns the dual sol of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *
 *  @note The success pointer indicates if the individual contraint handler was able to return the dual solution
 */
void SCIPconsGetDualsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get the dual solution for */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution */
   SCIP_Bool*            success             /**< pointer to store whether the dual solution is successfully returned */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   *success = TRUE;

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      *dualsol = SCIPgetDualsolLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      *dualsol = SCIPgetDualsolSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      *dualsol = SCIPgetDualsolLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      *dualsol = SCIPgetDualsolKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      *dualsol = SCIPgetDualsolVarbound(scip, cons);
   }
   /* these are Benders' specific constraint handlers */
   else if( strcmp(conshdlrname, "origbranch") == 0 || strcmp(conshdlrname, "masterbranch") == 0 )
   {
      *dualsol = 0.0;
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return dual solution for constraint of type <%s>\n", conshdlrname);
      *dualsol = 0.0;
      *success = FALSE;
   }
}

/** returns the row of an arbitrary SCIP constraint that can be represented as a single linear constraint
 *  or NULL of no row is awailable
 */
SCIP_ROW* SCIPconsGetRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint for which row is queried */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetRowLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return SCIPgetRowSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPgetRowLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetRowKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return SCIPgetRowVarbound(scip, cons);
   }
   else
   {
      SCIPwarningMessage(scip, "Cannot return row for constraint of type <%s>\n", conshdlrname);
   }

   return NULL;
}
