/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_trivial.c,v 1.13 2004/04/27 15:50:02 bzfpfend Exp $"

/**@file   presol_trivial.c
 * @brief  trivial presolver: round fractional bounds on integer variables, fix variables with equal bounds
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "presol_trivial.h"


#define PRESOL_NAME            "trivial"
#define PRESOL_DESC            "trivial presolver: round fractional bounds on integers, fix variables with equal bounds"
#define PRESOL_PRIORITY        +9000000




/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
#define presolFreeTrivial NULL


/** initialization method of presolver (called after problem was transformed) */
#define presolInitTrivial NULL


/** deinitialization method of presolver (called before transformed problem is freed) */
#define presolExitTrivial NULL


/** presolving execution method */
static
DECL_PRESOLEXEC(presolExecTrivial)
{  /*lint --e{715}*/
   VAR** vars;
   Real lb;
   Real ub;
   Bool infeasible;
   Bool fixed;
   int nvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* get the problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* scan the variables for trivial bound reductions
    * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
    */
   for( v = nvars-1; v >= 0; --v )
   {
      /* get variable's bounds */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);

      /* is variable integral? */
      if( SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS )
      {
         Real newlb;
         Real newub;
         
         /* round fractional bounds on integer variables */
         newlb = SCIPceil(scip, lb);
         newub = SCIPfloor(scip, ub);

         /* check bounds on variable for infeasibility */
         if( newlb > newub + 0.5 )
         {
            SCIPmessage(scip, SCIP_VERBLEVEL_NORMAL,
               "problem infeasible: integral variable <%s> has bounds [%g,%g] rounded to [%g,%g]\n",
               SCIPvarGetName(vars[v]), lb, ub, newlb, newub);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* fix variables with equal bounds */
         if( newlb > newub - 0.5 )
         {
            debugMessage("fixing integral variable <%s>: [%g,%g] -> [%g,%g]\n",
               SCIPvarGetName(vars[v]), lb, ub, newlb, newub);
            CHECK_OKAY( SCIPfixVar(scip, vars[v], newlb, &infeasible, &fixed) );
            if( infeasible )
            {
               debugMessage(" -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            assert(fixed);
            (*nfixedvars)++;
         }
         else
         {
            /* round fractional bounds */
            if( !SCIPisFeasEQ(scip, lb, newlb) )
            {
               debugMessage("rounding lower bound of integral variable <%s>: [%g,%g] -> [%g,%g]\n",
                  SCIPvarGetName(vars[v]), lb, ub, newlb, ub);
               CHECK_OKAY( SCIPchgVarLb(scip, vars[v], newlb) );
               (*nchgbds)++;
            }
            if( !SCIPisFeasEQ(scip, ub, newub) )
            {
               debugMessage("rounding upper bound of integral variable <%s>: [%g,%g] -> [%g,%g]\n",
                  SCIPvarGetName(vars[v]), newlb, ub, newlb, newub);
               CHECK_OKAY( SCIPchgVarUb(scip, vars[v], newub) );
               (*nchgbds)++;
            }
         }
      }
      else
      {
         /* check bounds on continuous variable for infeasibility */
         if( SCIPisFeasGT(scip, lb, ub) )
         {
            SCIPmessage(scip, SCIP_VERBLEVEL_NORMAL,
               "problem infeasible: continuous variable <%s> has bounds [%g,%g]\n",
               SCIPvarGetName(vars[v]), lb, ub);
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* fix variables with equal bounds */
         if( SCIPisFeasEQ(scip, lb, ub) )
         {
            debugMessage("fixing continuous variable <%s>: [%g,%g]\n", SCIPvarGetName(vars[v]), lb, ub);
            CHECK_OKAY( SCIPfixVar(scip, vars[v], lb, &infeasible, &fixed) );
            if( infeasible )
            {
               debugMessage(" -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            assert(fixed);
            (*nfixedvars)++;
         }
      }
   }

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the trivial presolver and includes it in SCIP */
RETCODE SCIPincludePresolTrivial(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PRESOLDATA* presoldata;

   /* create trivial presolver data */
   presoldata = NULL;

   /* include presolver */
   CHECK_OKAY( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY,
                  presolFreeTrivial, presolInitTrivial, presolExitTrivial, presolExecTrivial,
                  presoldata) );

   return SCIP_OKAY;
}
