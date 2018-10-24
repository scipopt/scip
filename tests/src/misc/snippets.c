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


/**! [SnippetDialogFree] */
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
/**! [SnippetDialogFree] */

/**! [SnippetDialogInclude] */
SCIP_RETCODE SCIPincludeDialogDrawgraph(
   SCIP*        scip,
   SCIP_DIALOG* root
   )
{
   SCIP_DIALOG* dialog;

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
/**! [SnippetDialogInclude] */

static
SCIP_DECL_DISPFREE(dispFreeMydisplaycolumn)
{
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);

   SCIPfreeMemory(scip, &dispdata);

   SCIPdispSetData(disp, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_TABLEFREE(tableFreeMystatisticstable)
{
   SCIP_TABLEDATA* tabledata;

   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);

   SCIPfreeMemory(scip, &tabledata);

   SCIPtableSetData(table, NULL);

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

Test(snippets, dialog, .description = "tests example Dialog code")
{
   /**! [SnippetDialogCreate] */
   SCIP_DIALOG* root;

   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
   }
   assert( root != NULL );
   /**! [SnippetDialogCreate] */

   SCIP_CALL( SCIPincludeDialogDrawgraph(scip, root) );
}

Test(snippets, memory, .description = "tests an example of memory allocation and freeing")
{
   /**! [SnippetArrayAllocAndFree] */
   int nparams;
   int* array;

   nparams = SCIPgetNParams(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &array, nparams) );

   /* do something ... */

   SCIPfreeBlockMemoryArray(scip, &array, nparams);
   /**! [SnippetArrayAllocAndFree] */
}

Test(snippets, display, .description = "tests an example of display code")
{
   
}

Test(snippets, table, .description = "tests an example of table code")
{

}
