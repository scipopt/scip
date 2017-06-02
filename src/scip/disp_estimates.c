/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp_estimates.c
 * @brief  estimates display column
 * @author Pierre Le Bodic
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/disp_estimates.h"
#include "scip/event_treesizeprediction.h"


#define DISP_NAME               "estimates"                
#define DISP_DESC               "estimates display column"
#define DISP_HEADER             "est n" 
#define DISP_WIDTH              7       /**< the width of the display column */
#define DISP_PRIORITY           110000  /**< the priority of the display column */
#define DISP_POSITION           30100   /**< the relative position of the display column */
#define DISP_STRIPLINE          TRUE    /**< the default for whether the display column should be separated 
                                         *   with a line from its right neighbor */


/*
 * Data structures
 */

/* TODO: fill in the necessary display column data */

/** display column data */
struct SCIP_DispData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of display column
 */

/* TODO: Implement all necessary display column methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for dialog plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_DISPCOPY(dispCopyEstimates)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of estimates display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispCopyEstimates NULL
#endif

/** destructor of display column to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_DISPFREE(dispFreeEstimates)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of estimates display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispFreeEstimates NULL
#endif


/** initialization method of display column (called after problem was transformed) */
#if 0
static
SCIP_DECL_DISPINIT(dispInitEstimates)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of estimates display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitEstimates NULL
#endif


/** deinitialization method of display column (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_DISPEXIT(dispExitEstimates)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of estimates display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitEstimates NULL
#endif


/** solving process initialization method of display column (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_DISPINITSOL(dispInitsolEstimates)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of estimates display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitsolEstimates NULL
#endif


/** solving process deinitialization method of display column (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_DISPEXITSOL(dispExitsolEstimates)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of estimates display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitsolEstimates NULL
#endif




/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputEstimates)
{  /*lint --e{715}*/
   SCIP_Longint estimate;

   assert(scip != NULL);
   /* We call the estimation method */
   estimate = SCIPtreeSizeGetEstimateRemaining(scip) + SCIPgetNNodes(scip);
   
   if( estimate != -1)
      SCIPdispLongint(SCIPgetMessagehdlr(scip), file, estimate, DISP_WIDTH);
   else
      SCIPinfoMessage(scip, file, "   -   ");
   
   return SCIP_OKAY;
}





/*
 * display column specific interface methods
 */

/** creates the estimates display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDispEstimates(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata;

   /* create estimates display column data */
   dispdata = NULL;
   /* TODO: (optional) create display column specific data here */

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO, 
         dispCopyEstimates,
         dispFreeEstimates, dispInitEstimates, dispExitEstimates, 
         dispInitsolEstimates, dispExitsolEstimates, dispOutputEstimates, 
         dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   /* add estimates display column parameters */
   /* TODO: (optional) add display column specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
