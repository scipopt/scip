/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cutsel_hello.c
 * @ingroup DEFPLUGINS_CUTSEL
 * @brief  hello cut selector
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip_cutsel.h"
#include "scip/cutsel_hello.h"


#define CUTSEL_NAME              "hello"
#define CUTSEL_DESC              "cut selector template"
#define CUTSEL_PRIORITY                 0


/*
 * Data structures
 */

/* TODO: fill in the necessary cut selector data */

/** cut selector data */
struct SCIP_CutselData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of cut selector
 */

static
SCIP_DECL_CUTSELCOPY(cutselCopyHello)
{  /*lint --e{715}*/
    assert(scip != NULL);
    assert(cutsel != NULL);
    assert(strcmp(SCIPcutselGetName(cutsel), CUTSEL_NAME) == 0);

   printf("We are being copied, Mark\n");

    /* call inclusion method of node selector */
    SCIP_CALL( SCIPincludeCutselHello(scip) );

    return SCIP_OKAY;
}

/* TODO: Implement all necessary cut selector methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for cut selector plugin (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CUTSELCOPY(cutselCopyHello)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of hello cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
//#define cutselCopyHello NULL
#endif

/** destructor of cut selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_CUTSELFREE(cutselFreeHello)
{  /*lint --e{715}*/
   printf("We are being freed, Mark\n");

   return SCIP_OKAY;
}


/** initialization method of cut selector (called after problem was transformed) */
static
SCIP_DECL_CUTSELINIT(cutselInitHello)
{  /*lint --e{715}*/
   printf("We are initializing, Mark\n");

   return SCIP_OKAY;
}


/** deinitialization method of cut selector (called before transformed problem is freed) */
static
SCIP_DECL_CUTSELEXIT(cutselExitHello)
{  /*lint --e{715}*/
  printf("We are exiting, Mark\n");

   return SCIP_OKAY;
}


/** solving process initialization method of cut selector (called when branch and bound process is about to begin) */
static
SCIP_DECL_CUTSELINITSOL(cutselInitsolHello)
{  /*lint --e{715}*/
   printf("We are initialising at the very beginning, Mark\n");

   return SCIP_OKAY;
}


/** solving process deinitialization method of cut selector (called before branch and bound process data is freed) */
static
SCIP_DECL_CUTSELEXITSOL(cutselExitsolHello)
{  /*lint --e{715}*/
   printf("We are ending the whole thing, Mark\n");

   return SCIP_OKAY;
}


/** cut selection method of cut selector */
static
SCIP_DECL_CUTSELSELECT(cutselSelectHello)
{  /*lint --e{715}*/
   //printf("Hello Mark, I am a cut selector plugin\n");

   *nselectedcuts = 0;

   return SCIP_OKAY;
}


/*
 * cut selector specific interface methods
 */

/** creates the hello cut selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeCutselHello(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CUTSELDATA* cutseldata;
   SCIP_CUTSEL* cutsel;

   /* create hello cut selector data */
   cutseldata = NULL;
   /* TODO: (optional) create cut selector specific data here */

   cutsel = NULL;

   /* include cut selector */
#if 1
   /* use SCIPincludeCutsel() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeCutsel(scip, CUTSEL_NAME, CUTSEL_DESC, CUTSEL_PRIORITY,
         cutselCopyHello, cutselFreeHello, cutselInitHello, cutselExitHello, cutselInitsolHello, cutselExitsolHello, cutselSelectHello,
         cutseldata) );
#else
   /* use SCIPincludeCutselBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeCutselBasic(scip, &cutsel, CUTSEL_NAME, CUTSEL_DESC, CUTSEL_PRIORITY, cutselSelectHello,
            cutseldata) );

   assert(cutsel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetCutselCopy(scip, cutsel, cutselCopyHello) );
   SCIP_CALL( SCIPsetCutselFree(scip, cutsel, cutselFreeHello) );
   SCIP_CALL( SCIPsetCutselInit(scip, cutsel, cutselInitHello) );
   SCIP_CALL( SCIPsetCutselExit(scip, cutsel, cutselExitHello) );
   SCIP_CALL( SCIPsetCutselInitsol(scip, cutsel, cutselInitsolHello) );
   SCIP_CALL( SCIPsetCutselExitsol(scip, cutsel, cutselExitsolHello) );
#endif

   /* add hello cut selector parameters */
   /* TODO: (optional) add cut selector specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
