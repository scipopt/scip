/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_integral.c
 * @brief  constraint handler for integral constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_integral.h"


#define CONSHDLR_NAME          "integral"
#define CONSHDLR_DESC          "Integrality constraint"
#define CONSHDLR_SEPAPRIORITY  -1000000
#define CONSHDLR_ENFOPRIORITY         0
#define CONSHDLR_CHECKPRIORITY        0
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS        FALSE /**< the constraint handler is called without constraints */



/*
 * Callback methods
 */

static
DECL_CONSENFOLP(consEnfolpIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   debugMessage("Enfolp method of integrality constraint\n");

   /* call branching methods */
   CHECK_OKAY( SCIPbranchLP(scip, result) );

   return SCIP_OKAY;
}

static
DECL_CONSCHECK(consCheckIntegral)
{
   VAR** vars;
   Real solval;
   int nbin;
   int nint;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   CHECK_OKAY( SCIPgetVars(scip, &vars, NULL, &nbin, &nint, NULL, NULL) );

   *result = SCIP_FEASIBLE;

   if( checkintegrality )
   {
      for( v = 0; v < nbin + nint && *result == SCIP_FEASIBLE; ++v )
      {
         CHECK_OKAY( SCIPgetSolVal(scip, sol, vars[v], &solval) );
         if( !SCIPisIntegral(scip, solval) )
            *result = SCIP_INFEASIBLE;
      }
   }

   return SCIP_OKAY;
}





/*
 * constraint specific interface methods
 */

/** creates the handler for integrality constraint and includes it in SCIP */
RETCODE SCIPincludeConsHdlrIntegral(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  NULL, NULL, NULL, 
                  NULL, NULL, 
                  NULL, consEnfolpIntegral, NULL, consCheckIntegral, NULL,
                  NULL, NULL,
                  NULL) );

   return SCIP_OKAY;
}
