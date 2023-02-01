/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   snippet.c
 * @brief  unit tests for code snippets
 * @author Ksenia Bestuzheva
 */

#include<stdio.h>

#include "scip/type_dialog.h"
#include "scip/type_disp.h"
#include "scip/pub_misc.h"
#include "scip/dialog_default.h"
#include "scip/scip.h"

#include "include/scip_test.h"

#define DISP_WIDTH              14                /**< the width of the display column */
#define DISP_PRIORITY           110000            /**< the priority of the display column */
#define DISP_POSITION           30100             /**< the relative position of the display column */
#define DISP_STRIPLINE          TRUE              /**< the default for whether the display column should be separated */

#define TABLE_POSITION          -1000             /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_INIT   /**< output of the statistics table is only printed from this stage onwards */

/** dialog data */
struct SCIP_DialogData
{
   int some_dialog_data;
   /* some data */
};

/** display column data */
struct SCIP_DispData
{
   int some_disp_data;
   /* some data */
};

/** statistics table data */
struct SCIP_TableData
{
   int some_table_data;
   /* some data */
};

/** methods **/


/** Dialog methods **/

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
static
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


/** Display methods **/

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputMydisplaycolumn)
{
   SCIPinfoMessage(scip, file, "\nOutput method of mydisplaycolumn\n");

   return SCIP_OKAY;
}

/**! [SnippetDispFree] */
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
/**! [SnippetDispFree] */

/** creates my display column and includes it in SCIP */
static
SCIP_RETCODE SCIPincludeDispMydisplaycolumn(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_DISPDATA* dispdata;

   /* create mydisplaycolumn display column data */
   SCIP_CALL( SCIPallocMemory(scip, &dispdata) );

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, "my display column", "an example of display code", "mydisplaycolumn", SCIP_DISPSTATUS_AUTO,
            NULL,
            dispFreeMydisplaycolumn, NULL, NULL,
            NULL, NULL, dispOutputMydisplaycolumn,
            dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   return SCIP_OKAY;
}

/** Table methods **/

/**! [SnippetTableFree] */
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
/**! [SnippetTableFree] */

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputMystatisticstable)
{
   SCIPinfoMessage(scip, file, "output method of my statistics table\n");

   return SCIP_OKAY;
}

/** creates the statistics table and includes it in SCIP */
static
SCIP_RETCODE SCIPincludeMytable(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_TABLEDATA* tabledata;

   /* create statistics table data */
   SCIP_CALL( SCIPallocMemory(scip, &tabledata) );

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, "my statistics table", "an example of a statistics table", TRUE,
          NULL, tableFreeMystatisticstable, NULL, NULL,
          NULL, NULL, tableOutputMystatisticstable,
          tabledata, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   return SCIP_OKAY;
}

/** GLOBAL VARIABLES **/
static SCIP* scip;

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

Test(snippets, dialog, .description = "tests an example of dialog code")
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
   SCIPincludeDispMydisplaycolumn(scip);
}

Test(snippets, table, .description = "tests an example of table code")
{
   SCIPincludeMytable(scip);
}
