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


/** constraint data for integrality constraint */
struct ConsData
{
   int              dummy;              /**< dummy to not have an empty struct */
};


/*
 * Callback methods
 */

DECL_CONSINIT(SCIPconsInitIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "initialise integrality constraint handler");

   return SCIP_OKAY;
}

DECL_CONSEXIT(SCIPconsExitIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   infoMessage(SCIPverbLevel(scip), SCIP_VERBLEVEL_FULL, "exit integrality constraint handler");

   return SCIP_OKAY;
}

DECL_CONSFREE(SCIPconsFreeIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata == NULL);

   /* nothing to do here, because the integraliy constraint has no constraint data */

   return SCIP_OKAY;
}

DECL_CONSTRAN(SCIPconsTranIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(SCIPstage(scip) == SCIP_STAGE_SOLVING);
   assert(sourcedata == NULL);
   assert(targetdata != NULL);

   debugMessage("Tran method of integrality constraint\n");

   *targetdata = NULL;
   
   return SCIP_OKAY;
}

DECL_CONSSEPA(SCIPconsSepaIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);

   debugMessage("Sepa method of integrality constraint\n");

   /* nothing to do here, because the integrality constraint handler doesn't separate anything */

   return SCIP_DIDNOTRUN;
}

DECL_CONSENFO(SCIPconsEnfoIntegral)
{
   VAR** vars;
   VAR* var;
   int nbin;
   int nint;
   int v;
   Real primsol;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);

   debugMessage("Enfo method of integrality constraint\n");

   CHECK_OKAY( SCIPgetVars(scip, &vars, NULL, &nbin, &nint, NULL, NULL) );

   todoMessage("variable selection rules");

   for( v = 0; v < nbin+nint; ++v )
   {
      var = vars[v];
      primsol = SCIPvarGetPrimsol(var);
      if( !SCIPisIntegral(scip, primsol) )
      {
         NODE* node;

         debugMessage("branching on variable <%s> with value %f\n", SCIPvarGetName(var), primsol);

         /* create left child */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         CHECK_OKAY( SCIPchgNodeUb(scip, node, var, SCIPfloor(scip, primsol)) );

         /* create right child */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         CHECK_OKAY( SCIPchgNodeLb(scip, node, var, SCIPceil(scip, primsol)) );

         return SCIP_BRANCHED;
      }
   }

   return SCIP_FEASIBLE;
}

DECL_CONSCHCK(SCIPconsChckIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   todoMessage("Chck method of integrality constraint");

   return SCIP_OKAY;
}

DECL_CONSPROP(SCIPconsPropIntegral)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPgetConsHdlrName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   todoMessage("Prop method of integrality constraint");

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
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHCKPRIORITY,
                  SCIPconsInitIntegral, SCIPconsExitIntegral, SCIPconsFreeIntegral,
                  SCIPconsTranIntegral, SCIPconsSepaIntegral, SCIPconsEnfoIntegral,
                  SCIPconsChckIntegral, SCIPconsPropIntegral,
                  NULL) );

   return SCIP_OKAY;
}
