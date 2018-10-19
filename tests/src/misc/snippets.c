/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   snippet.c
 * @brief  unit tests for code snippets
 * @author Ksenia Bestuzheva
 */

#include<stdio.h>

#include "scip/type_dialog.h"
#include "scip/pub_misc.h"
#include "scip/dialog_default.h"
#include "scip/scip.h"

#include "include/scip_test.h"

struct SCIP_DialogData
{
};

/** GLOBAL VARIABLES **/
static SCIP* scip;

/** methods **/

static
SCIP_DECL_DIALOGFREE(dialogFreeDrawgraph)
{
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != NULL);

   SCIPfreeMemory(scip, &dialogdata);

   SCIPdialogSetData(dialog, NULL);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPincludeDialogDrawgraph(
   SCIP*  scip
   )
{
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* root;

   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
   }
   assert( root != NULL );

   SCIP_DIALOGDATA* dialogdata;
   SCIP_CALL( SCIPallocMemory(scip, &dialogdata) );

   if( !SCIPdialogHasEntry(root, "drawgraph") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, NULL, NULL, dialogFreeDrawgraph,
                                                "drawgraph", "draws the graph for the current problem instance", FALSE, dialogdata) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   return SCIP_OKAY;
}

/* TEST SUITE */
static
void setup(void)
{
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   cr_assert_not_null(scip);
}

static
void teardown(void)
{

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(snippets, .init = setup, .fini = teardown);

Test(snippets, dialog, .description = "tests Dialog code")
{
   SCIP_CALL( SCIPincludeDialogDrawgraph(scip) );
}
