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
#pragma ident "@(#) $Id: presol_xxx.c,v 1.7 2004/02/04 17:27:33 bzfpfend Exp $"

/**@file   presol_xxx.c
 * @brief  xxx presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "presol_xxx.h"


#define PRESOL_NAME            "xxx"
#define PRESOL_DESC            "presolver template"
#define PRESOL_PRIORITY        0




/*
 * Data structures
 */

/* TODO: fill in the necessary presolver data */

/** presolver data */
struct PresolData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of presolver to free user data (called when SCIP is exiting) */
#if 0
static
DECL_PRESOLFREE(presolFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx presolver not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolFreeXxx NULL
#endif


/** initialization method of presolver (called when problem solving starts) */
#if 0
static
DECL_PRESOLINIT(presolInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx presolver not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitXxx NULL
#endif


/** deinitialization method of presolver (called when problem solving exits) */
#if 0
static
DECL_PRESOLEXIT(presolExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx presolver not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitXxx NULL
#endif


/** execution method of presolver */
static
DECL_PRESOLEXEC(presolExecXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx presolver not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the xxx presolver and includes it in SCIP */
RETCODE SCIPincludePresolXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PRESOLDATA* presoldata;

   /* create xxx presolver data */
   presoldata = NULL;
   /* TODO: (optional) create presolver specific data here */

   /* include presolver */
   CHECK_OKAY( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY,
                  presolFreeXxx, presolInitXxx, presolExitXxx, presolExecXxx,
                  presoldata) );

   /* add xxx presolver parameters */
   /* TODO: (optional) add presolver specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
