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
#define CONSHDLR_CHCKPRIORITY         0
#define CONSHDLR_NEEDSCONS        FALSE /**< the constraint handler is called without constraints */



/*
 * Callback methods
 */

static
DECL_CONSENLP(consEnlpIntegral)
{
   VAR** lpcands;
   int nlpcands;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   debugMessage("Enlp method of integrality constraint\n");

   /* get branching candidates */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, &nlpcands) );
   if( nlpcands == 0 )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
   *result = SCIP_INFEASIBLE;

   /* call branching methods */
   todoMessage("variable selection rules");

   /* if no branching method succeeded in choosing a branching, just branch on the first fractional variable */
   if( *result == SCIP_INFEASIBLE )
   {
      NODE* node;
      VAR* var;
      Real primsol;

      var = lpcands[0];
      assert(var != NULL);

      CHECK_OKAY( SCIPgetActVarSol(scip, var, &primsol) );
      assert(!SCIPisIntegral(scip, primsol));

      debugMessage("stupid branching on variable <%s> with value %f\n", SCIPvarGetName(var), primsol);
      
      /* create left child */
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgNodeUb(scip, node, var, SCIPfloor(scip, primsol)) );
      
      /* create right child */
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      CHECK_OKAY( SCIPchgNodeLb(scip, node, var, SCIPceil(scip, primsol)) );
      
      *result = SCIP_BRANCHED;
   }

   assert(*result != SCIP_INFEASIBLE);

   return SCIP_OKAY;
}

static
DECL_CONSCHCK(consChckIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   todoMessage("Chck method of integrality constraint");

   return SCIP_OKAY;
}





/*
 * constraint specific interface methods
 */

RETCODE SCIPincludeConsHdlrIntegral(      /**< creates the handler for integrality constraint and includes it in SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHCKPRIORITY, CONSHDLR_NEEDSCONS,
                  NULL, NULL, NULL, 
                  NULL, NULL, 
                  NULL, consEnlpIntegral, NULL, consChckIntegral, NULL,
                  NULL) );

   return SCIP_OKAY;
}
