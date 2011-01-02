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
#pragma ident "@(#) $Id: nodesel_xxx.c,v 1.23 2011/01/02 11:10:44 bzfheinz Exp $"

/**@file   nodesel_xxx.c
 * @ingroup NODESELECTORS
 * @brief  xxx node selector
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/nodesel_xxx.h"


#define NODESEL_NAME            "xxx"
#define NODESEL_DESC            "node selector template"
#define NODESEL_STDPRIORITY     0
#define NODESEL_MEMSAVEPRIORITY 0




/*
 * Data structures
 */

/* TODO: fill in the necessary node selector data */

/** node selector data */
struct SCIP_NodeselData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of node selector
 */

/* TODO: Implement all necessary node selector methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for node selector plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_NODESELCOPY(nodeselCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define nodeselCopyXxx NULL
#endif

/** destructor of node selector to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_NODESELFREE(nodeselFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselFreeXxx NULL
#endif


/** initialization method of node selector (called after problem was transformed) */
#if 0
static
SCIP_DECL_NODESELINIT(nodeselInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselInitXxx NULL
#endif


/** deinitialization method of node selector (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_NODESELEXIT(nodeselExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitXxx NULL
#endif


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselInitsolXxx NULL
#endif


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_NODESELEXITSOL(nodeselExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitsolXxx NULL
#endif


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return 0;
}




/*
 * node selector specific interface methods
 */

/** creates the xxx node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   /* create xxx node selector data */
   nodeseldata = NULL;
   /* TODO: (optional) create node selector specific data here */

   /* include node selector */
   SCIP_CALL( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
         nodeselCopyXxx,
         nodeselFreeXxx, nodeselInitXxx, nodeselExitXxx, 
         nodeselInitsolXxx, nodeselExitsolXxx, nodeselSelectXxx, nodeselCompXxx,
         nodeseldata) );

   /* add xxx node selector parameters */
   /* TODO: (optional) add node selector specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
