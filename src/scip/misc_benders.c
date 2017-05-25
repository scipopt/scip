/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    benders_misc.c
 * @brief   various SCIP helper methods for Benders' decomposition
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/misc_benders.h"
#include "scip/scipdefplugins.h"
#include <string.h>

/** returns TRUE if variable is relevant, FALSE otherwise */
SCIP_Bool BDisVarRelevant(
   SCIP_VAR*             var                 /**< variable to test */
   )
{
   assert(var != NULL);
   return (SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED)
      && (SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
}

/** returns the type of an arbitrary SCIP constraint */
consType BDconsGetType(
   SCIP_CONS*            cons                /**< constraint to get type for */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return linear;
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch ( SCIPgetTypeSetppc(NULL, cons) ) {
      case SCIP_SETPPCTYPE_COVERING:
         return setcovering;
      case SCIP_SETPPCTYPE_PACKING:
         return setpacking;
      case SCIP_SETPPCTYPE_PARTITIONING:
         return setpartitioning;
      default:
         return unknown;
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return logicor;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return knapsack;
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return varbound;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      return sos1;
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      return sos2;
   }
   return unknown;
}

/** returns the rhs of an arbitrary SCIP constraint */
SCIP_Real BDconsGetRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   )
{
   int nconsvars;
   SCIP_VAR** consvars;
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   SCIP_Longint * w;
   SCIP_Real rhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetRhsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch ( SCIPgetTypeSetppc(scip, cons) ) {
      case SCIP_SETPPCTYPE_PARTITIONING: /* fall through desired */
      case SCIP_SETPPCTYPE_PACKING:
         rhs = 1.0;

         nconsvars = SCIPgetNVarsSetppc(scip, cons);
         consvars = SCIPgetVarsSetppc(scip, cons);

         for( i = 0; i < nconsvars; i++ )
         {
            if( SCIPvarIsNegated(consvars[i]) )
               rhs -= 1.0;
         }

         return rhs;

      case SCIP_SETPPCTYPE_COVERING:
         return SCIPinfinity(scip);
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPinfinity(scip);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      rhs = SCIPgetCapacityKnapsack(scip, cons);

      /* copy Longint array to SCIP_Real array */
      w = SCIPgetWeightsKnapsack(scip, cons);
      consvars = SCIPgetVarsKnapsack(scip, cons);
      nconsvars = SCIPgetNVarsKnapsack(scip, cons);

      for( i = 0; i < nconsvars; i++ )
      {
         if( SCIPvarIsNegated(consvars[i]) )
            rhs -= w[i];
      }

      return rhs;
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      rhs = SCIPgetRhsVarbound(scip, cons);

      if( SCIPvarIsNegated(SCIPgetVarVarbound(scip, cons)) && !SCIPisInfinity(scip, rhs) )
         rhs -= 1.0;

      if( SCIPvarIsNegated(SCIPgetVbdvarVarbound(scip, cons)) && !SCIPisInfinity(scip, rhs) )
         rhs -= SCIPgetVbdcoefVarbound(scip, cons);

      return rhs;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED");
   }
   return -SCIPinfinity(scip);
}

/** returns the lhs of an arbitrary SCIP constraint */
SCIP_Real BDconsGetLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   int nconsvars;
   SCIP_VAR** consvars;
   SCIP_Real lhs;
   int i;


   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetLhsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch ( SCIPgetTypeSetppc(scip, cons) ) {
      case SCIP_SETPPCTYPE_PARTITIONING: /* fall through desired */
      case SCIP_SETPPCTYPE_COVERING:
         lhs = 1.0;

         nconsvars = SCIPgetNVarsSetppc(scip, cons);
         consvars = SCIPgetVarsSetppc(scip, cons);

         for( i = 0; i < nconsvars; i++ )
         {
            if( SCIPvarIsNegated(consvars[i]) )
               lhs -= 1.0;
         }

         return lhs;

      case SCIP_SETPPCTYPE_PACKING:
         return -SCIPinfinity(scip);
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      lhs = 1.0;

      nconsvars = SCIPgetNVarsLogicor(scip, cons);
      consvars = SCIPgetVarsLogicor(scip, cons);

      for( i = 0; i < nconsvars; i++ )
      {
         if( SCIPvarIsNegated(consvars[i]) )
            lhs -= 1.0;
      }

      return lhs;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return -SCIPinfinity(scip);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      lhs = SCIPgetLhsVarbound(scip, cons);

      if( SCIPvarIsNegated(SCIPgetVarVarbound(scip, cons)) && !SCIPisInfinity(scip, -lhs) )
         lhs -= 1.0;

      if( SCIPvarIsNegated(SCIPgetVbdvarVarbound(scip, cons)) && !SCIPisInfinity(scip, -lhs) )
         lhs -= SCIPgetVbdcoefVarbound(scip, cons);

      return lhs;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED");
   }
   return SCIPinfinity(scip);
}

/** returns the dual farkas sol of an arbitrary SCIP constraint */
extern
SCIP_Real BDconsGetDualfarkas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetDualfarkasLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return SCIPgetDualfarkasSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPgetDualfarkasLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetDualfarkasKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return SCIPgetDualfarkasVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "origbranch") == 0 )
   {
      SCIPdebugMessage("origbranch: return dualfarkas 0\n");
      return 0.0;
   }
   else if( strcmp(conshdlrname, "masterbranch") == 0 )
   {
      SCIPdebugMessage("masterbranch: return dualsol 0\n");
      return 0.0;
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED");
   }
   return SCIPinfinity(scip);
}

/** returns the dual sol of an arbitrary SCIP constraint */
extern
SCIP_Real BDconsGetDualsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetDualsolLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return SCIPgetDualsolSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPgetDualsolLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetDualsolKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return SCIPgetDualsolVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "origbranch") == 0 )
   {
      SCIPdebugMessage("origbranch: return Dualsol 0\n");
      return 0.0;
   }
   else if( strcmp(conshdlrname, "masterbranch") == 0 )
   {
      SCIPdebugMessage("masterbranch: return dualsol 0\n");
      return 0.0;
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED");
   }
   return SCIPinfinity(scip);
}



/** returns the row of an arbitrary SCIP constraint */
extern
SCIP_ROW* BDconsGetRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get left hand side for */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;

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
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
      return NULL;
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
      return NULL;
   }
   else if( strcmp(conshdlrname, "origbranch") == 0 )
   {
      SCIPdebugMessage("origbranch: return Dualsol 0\n");
      return NULL;
   }
   else if( strcmp(conshdlrname, "masterbranch") == 0 )
   {
      SCIPdebugMessage("masterbranch: return dualsol 0\n");
      return NULL;
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED");
      return NULL;
   }
   return NULL;
}



/** adds a coefficient to an arbitrary SCIP constraint */
extern
SCIP_RETCODE BDconsAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);
   assert(var != NULL);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, val) );
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      SCIP_CALL( SCIPaddCoefSetppc(scip, cons, var) );
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      SCIP_CALL( SCIPaddCoefLogicor(scip, cons, var) );
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      SCIP_CALL( SCIPaddCoefKnapsack(scip, cons, var, val) );
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      SCIPdebugMessage("WARNING: varbound NOT IMPLEMENTED\n");
      return SCIP_OKAY;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
      return SCIP_OKAY;
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
      return SCIP_OKAY;
   }
   else if( strcmp(conshdlrname, "origbranch") == 0 )
   {
      SCIPdebugMessage("origbranch: return Dualsol 0\n");
      return SCIP_OKAY;
   }
   else if( strcmp(conshdlrname, "masterbranch") == 0 )
   {
      SCIPdebugMessage("masterbranch: return dualsol 0\n");
      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED");
   }
   return SCIP_OKAY;
}



/** returns the number of variables in an arbitrary SCIP constraint */
int BDconsGetNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get number of variables */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;

   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetNVarsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return SCIPgetNVarsSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPgetNVarsLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetNVarsKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return 2;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      return SCIPgetNVarsSOS1(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      return SCIPgetNVarsSOS2(scip, cons);
   }
   else
   {
      SCIPdebugMessage("WARNING: NOT IMPLEMENTED <%s>\n", conshdlrname);
      return 0;
   }
}

/** returns the variable array of an arbitrary SCIP constraint */
SCIP_RETCODE BDconsGetVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get variables from */
   SCIP_VAR**            vars,               /**< array where variables are stored */
   int                   nvars               /**< size of storage array */
   )
{

   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vars != NULL);
   assert(nvars > 0);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      if( nvars < SCIPgetNVarsLinear(scip, cons) )
         return SCIP_INVALIDDATA;

      BMScopyMemoryArray(vars, SCIPgetVarsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons));
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      SCIP_VAR** consvars;
      int nconsvars;

      consvars = SCIPgetVarsSetppc(scip, cons);
      nconsvars = SCIPgetNVarsSetppc(scip, cons);

      if( nvars < nconsvars )
         return SCIP_INVALIDDATA;

      for( i = 0; i < nconsvars; i++ )
         if( !SCIPvarIsNegated(consvars[i]) )
            vars[i] = consvars[i];
         else
            vars[i] = SCIPvarGetNegatedVar(consvars[i]);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      SCIP_VAR** consvars;
      int nconsvars;

      consvars = SCIPgetVarsLogicor(scip, cons);
      nconsvars = SCIPgetNVarsLogicor(scip, cons);

      if( nvars < nconsvars )
         return SCIP_INVALIDDATA;

      for( i = 0; i < nconsvars; i++ )
         if( !SCIPvarIsNegated(consvars[i]) )
            vars[i] = consvars[i];
         else
            vars[i] = SCIPvarGetNegatedVar(consvars[i]);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      SCIP_VAR** consvars;
      int nconsvars;

      consvars = SCIPgetVarsKnapsack(scip, cons);
      nconsvars = SCIPgetNVarsKnapsack(scip, cons);

      if( nvars < nconsvars )
         return SCIP_INVALIDDATA;

      for( i = 0; i < nconsvars; i++ )
         if( !SCIPvarIsNegated(consvars[i]) )
            vars[i] = consvars[i];
         else
            vars[i] = SCIPvarGetNegatedVar(consvars[i]);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      if( nvars < 2 )
         return SCIP_INVALIDDATA;

      if( !SCIPvarIsNegated(SCIPgetVarVarbound(scip, cons)) )
         vars[0] = SCIPgetVarVarbound(scip, cons);
      else
         vars[0] = SCIPvarGetNegatedVar(SCIPgetVarVarbound(scip, cons));

      if( !SCIPvarIsNegated(SCIPgetVbdvarVarbound(scip, cons)) )
         vars[1] = SCIPgetVbdvarVarbound(scip, cons);
      else
         vars[1] = SCIPvarGetNegatedVar(SCIPgetVbdvarVarbound(scip, cons));
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      if( nvars < SCIPgetNVarsSOS1(scip, cons) )
         return SCIP_INVALIDDATA;

      BMScopyMemoryArray(vars, SCIPgetVarsSOS1(scip, cons), nvars);
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      if( nvars < SCIPgetNVarsSOS2(scip, cons) )
         return SCIP_INVALIDDATA;

      BMScopyMemoryArray(vars, SCIPgetVarsSOS2(scip, cons), nvars);
   }
   else
   {
      SCIPwarningMessage(scip, "WARNING: NOT IMPLEMENTED <%s>\n", conshdlrname);
   }
   return SCIP_OKAY;
}


/**
 * Returns the value array of an arbitrary SCIP constraint
 * @todo SOS1 & SOS2 not implemented yet
 */
SCIP_RETCODE BDconsGetVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get values from */
   SCIP_Real*            vals,               /**< array where values are stored */
   int                   nvals               /**< size of storage array */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   int i;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(vals != NULL);
   assert(nvals > 0);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      nvars = SCIPgetNVarsLinear(scip, cons);
      if( nvals < nvars )
         return SCIP_INVALIDDATA;

      BMScopyMemoryArray(vals, SCIPgetValsLinear(scip, cons), nvars);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      SCIP_VAR** vars;
      vars = SCIPgetVarsSetppc(scip, cons);
      nvars = SCIPgetNVarsSetppc(scip, cons);
      if( nvals < nvars )
         return SCIP_INVALIDDATA;

      for( i = 0; i < nvals; i++ )
         if( !SCIPvarIsNegated(vars[i]) )
            vals[i] = 1.0;
         else
            vals[i] = -1.0;
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      SCIP_VAR** vars;
      vars = SCIPgetVarsLogicor(scip, cons);
      nvars = SCIPgetNVarsLogicor(scip, cons);
      if( nvals < nvars )
         return SCIP_INVALIDDATA;

      for( i = 0; i < nvals; i++ )
      {
         if( !SCIPvarIsNegated(vars[i]) )
            vals[i] = 1.0;
         else
            vals[i] = -1.0;
      }
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      SCIP_VAR** vars;

      /* copy Longint array to SCIP_Real array */
      SCIP_Longint * w = SCIPgetWeightsKnapsack(scip, cons);
      vars = SCIPgetVarsKnapsack(scip, cons);
      nvars = SCIPgetNVarsKnapsack(scip, cons);
      if( nvals < nvars )
         return SCIP_INVALIDDATA;

      for( i = 0; i < nvars; i++ )
      {
         if( !SCIPvarIsNegated(vars[i]) )
            vals[i] = w[i];
         else
            vals[i] = -w[i];
      }
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      nvars = 2;
      if( nvals < nvars )
         return SCIP_INVALIDDATA;

      if( !SCIPvarIsNegated(SCIPgetVarVarbound(scip, cons)) )
         vals[0] = 1.0;
      else
         vals[0] = -1.0;

      if( !SCIPvarIsNegated(SCIPgetVbdvarVarbound(scip, cons)) )
         vals[1] = SCIPgetVbdcoefVarbound(scip, cons);
      else
         vals[1] = -SCIPgetVbdcoefVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      /* store constraint */
      SCIPerrorMessage("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      /* store constraint */
      SCIPdebugMessage("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPdebugMessage("WARNING: UNKNOWN NOT IMPLEMENTED: %s\n", conshdlrname);
   }
   return SCIP_OKAY;
}


/** returns true if the constraint should be a master constraint and false otherwise */
SCIP_Bool BDgetConsIsSetppc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SETPPCTYPE*      setppctype          /**< returns the type of the constraints */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int i;

   int nvars;
   SCIP_Bool relevant = TRUE;
   assert(scip != NULL);
   assert(cons != NULL);
   assert(setppctype != NULL);

   *setppctype = SCIP_SETPPCTYPE_PACKING;
   SCIPdebugMessage("cons %s is ", SCIPconsGetName(cons));

   if( BDconsGetType(cons) == setcovering || BDconsGetType(cons) == setpartitioning || BDconsGetType(cons) == logicor )
   {
      SCIPdebugPrintf("setcov, part or logicor.\n");
      return TRUE;
   }
   nvars = BDconsGetNVars(scip, cons);
   vars = NULL;
   vals = NULL;
   if( nvars > 0 )
   {
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals, nvars) );
      SCIP_CALL_ABORT( BDconsGetVars(scip, cons, vars, nvars) );
      SCIP_CALL_ABORT( BDconsGetVals(scip, cons, vals, nvars) );
   }

   /* check vars and vals for integrality */
   for( i = 0; i < nvars && relevant; ++i )
   {
      assert(vars != NULL);
      assert(vals != NULL);

      if( !SCIPvarIsBinary(vars[i]) )
      {
         SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[i]) );
         relevant = FALSE;
      }
      if( !SCIPisEQ(scip, vals[i], 1.0) )
      {
         SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[i]), vals[i] );
         relevant = FALSE;
      }
   }

   if( relevant )
   {
      SCIP_Real rhs = BDconsGetRhs(scip, cons);
      SCIP_Real lhs = BDconsGetLhs(scip, cons);
      SCIPdebugPrintf("(lhs %.2f, rhs %.2f)", lhs, rhs);

      if( SCIPisEQ(scip, lhs, 1.0) && SCIPisEQ(scip, rhs, 1.0) )
      {
         *setppctype = SCIP_SETPPCTYPE_PARTITIONING;
      }
      else if( SCIPisEQ(scip, lhs, 1.0) && SCIPisGE(scip, rhs, nvars*1.0) )
      {
         *setppctype = SCIP_SETPPCTYPE_COVERING;
      }
      else if( SCIPisLE(scip, lhs, 0.0) && SCIPisEQ(scip, rhs, 1.0) )
      {
         *setppctype = SCIP_SETPPCTYPE_PACKING;
      }
      else
      {
         relevant = FALSE;
      }
   }

   /* free temporary data  */
   SCIPfreeBufferArrayNull(scip, &vals);
   SCIPfreeBufferArrayNull(scip, &vars);

   SCIPdebugPrintf("%s master\n", relevant ? "in" : "not in");
   return relevant;
}

/** returns TRUE or FALSE, depending whether we are in the root node or not */
SCIP_Bool BDisRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   return (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip));
}
