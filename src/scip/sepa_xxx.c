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
#pragma ident "@(#) $Id: sepa_xxx.c,v 1.6 2004/07/14 14:05:09 bzfpfend Exp $"

/**@file   sepa_xxx.c
 * @brief  xxx separator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_xxx.h"


#define SEPA_NAME              "xxx"
#define SEPA_DESC              "separator template"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10




/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SepaData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
DECL_SEPAFREE(sepaFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeXxx NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
DECL_SEPAINIT(sepaInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitXxx NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
DECL_SEPAEXIT(sepaExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitXxx NULL
#endif


/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx separator not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * separator specific interface methods
 */

/** creates the xxx separator and includes it in SCIP */
RETCODE SCIPincludeSepaXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SEPADATA* sepadata;

   /* create xxx separator data */
   sepadata = NULL;
   /* TODO: (optional) create separator specific data here */

   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ,
         sepaFreeXxx, sepaInitXxx, sepaExitXxx, sepaExecXxx,
         sepadata) );

   /* add xxx separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
