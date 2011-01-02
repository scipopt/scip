/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: disp_xxx.c,v 1.8 2011/01/02 11:10:47 bzfheinz Exp $"

/**@file   disp_xxx.c
 * @ingroup DISPLAYS
 * @brief  xxx display column
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/disp_xxx.h"


#define DISP_NAME               "xxx"                
#define DISP_DESC               "xxx display column"
#define DISP_HEADER             "xxx" 
#define DISP_WIDTH              14      /**< the width of the display column */
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
SCIP_DECL_DISPCOPY(dispCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define dispCopyXxx NULL
#endif

/** destructor of display column to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_DISPFREE(dispFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispFreeXxx NULL
#endif


/** initialization method of display column (called after problem was transformed) */
#if 0
static
SCIP_DECL_DISPINIT(dispInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitXxx NULL
#endif


/** deinitialization method of display column (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_DISPEXIT(dispExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitXxx NULL
#endif


/** solving process initialization method of display column (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_DISPINITSOL(dispInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitsolXxx NULL
#endif


/** solving process deinitialization method of display column (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_DISPEXITSOL(dispExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitsolXxx NULL
#endif




/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * display column specific interface methods
 */

/** creates the xxx display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDispXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata;

   /* create xxx display column data */
   dispdata = NULL;
   /* TODO: (optional) create display column specific data here */

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO, 
         dispCopyXxx,
         dispFreeXxx, dispInitXxx, dispExitXxx, 
         dispInitsolXxx, dispExitsolXxx, dispOutputXxx, 
         dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   /* add xxx display column parameters */
   /* TODO: (optional) add display column specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
