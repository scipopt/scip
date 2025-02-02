/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file    scip_datatree.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for managing data trees
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip_datatree.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/pub_message.h"
#include "scip/set.h"
#include "scip/datatree.h"
#include "scip/pub_datatree.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_datatree.h"


/** creates a new SCIP_DATATREE */
SCIP_RETCODE SCIPcreateDatatree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE**       datatree,           /**< buffer to store created data tree */
   int                   capacity            /**< desired capacity, or -1 for default */
)
{
   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPdatatreeCreate(datatree, scip->mem->setmem, capacity) );

   return SCIP_OKAY;
}

/** creates a new SCIP_DATATREE and inserts it into a SCIP_DATATREE object */
SCIP_RETCODE SCIPcreateDatatreeInTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree where to insert new tree */
   SCIP_DATATREE**       newtree,            /**< buffer to store pointer to created data tree */
   const char*           name,               /**< name of entry to add */
   int                   capacity            /**< capacity of new tree, or -1 for default */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(newtree != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPdatatreeCreate(newtree, scip->mem->setmem, capacity) );

   SCIP_CALL( SCIPdatatreeInsertTree(datatree, scip->set, scip->mem->setmem, name, *newtree) );

   return SCIP_OKAY;
}

/** inserts a bool value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeBool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_Bool             value               /**< value to add */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPdatatreeInsertBool(datatree, scip->set, scip->mem->setmem, name, value) );

   return SCIP_OKAY;
}

/** inserts an int value into a SCIP_DATATREE object
 *
 *  The value will be stored as SCIP_Longint.
 */
SCIP_RETCODE SCIPinsertDatatreeInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   int                   value               /**< value to add */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPdatatreeInsertLong(datatree, scip->set, scip->mem->setmem, name, (SCIP_Longint)value) );

   return SCIP_OKAY;
}

/** inserts a long value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeLong(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_Longint          value               /**< value to add */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPdatatreeInsertLong(datatree, scip->set, scip->mem->setmem, name, value) );

   return SCIP_OKAY;
}

/** inserts a SCIP_Real value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeReal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_Real             value               /**< value to add */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPdatatreeInsertReal(datatree, scip->set, scip->mem->setmem, name, value) );

   return SCIP_OKAY;
}

/** inserts a string value into a SCIP_DATATREE object
 *
 *  The string value will be copied.
 */
SCIP_RETCODE SCIPinsertDatatreeString(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const char*           value               /**< value to add */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(value != NULL);

   SCIP_CALL( SCIPdatatreeInsertString(datatree, scip->set, scip->mem->setmem, name, value) );

   return SCIP_OKAY;
}

/** inserts a SCIP_Bool array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeBoolArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const SCIP_Bool*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues >= 0);

   SCIP_CALL( SCIPdatatreeInsertBoolArray(datatree, scip->set, scip->mem->setmem, name, values, nvalues) );

   return SCIP_OKAY;
}

/** inserts an int array into a SCIP_DATATREE object
 *
 *  The value will be stored as array of SCIP_Longint.
 */
SCIP_RETCODE SCIPinsertDatatreeIntArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const int*            values,             /**< values of entry */
   int                   nvalues             /**< number of values */
)
{
   SCIP_Longint* longvalues = NULL;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues >= 0);

   /* turn int array into long array */
   SCIP_CALL( SCIPallocBufferArray(scip, &longvalues, nvalues) );
   for( int i = 0; i < nvalues; i++ )
      longvalues[i] = (SCIP_Longint)values[i];

   SCIP_CALL( SCIPdatatreeInsertLongArray(datatree, scip->set, scip->mem->setmem, name, longvalues, nvalues) );

   SCIPfreeBufferArray(scip, &longvalues);

   return SCIP_OKAY;
}

/** inserts a SCIP_Real array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeRealArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const SCIP_Real*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues >= 0);

   SCIP_CALL( SCIPdatatreeInsertRealArray(datatree, scip->set, scip->mem->setmem, name, values, nvalues) );

   return SCIP_OKAY;
}

/** inserts a SCIP_Longint array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeLongArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const SCIP_Longint*   values,             /**< values of entry */
   int                   nvalues             /**< number of values */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues >= 0);

   SCIP_CALL( SCIPdatatreeInsertLongArray(datatree, scip->set, scip->mem->setmem, name, values, nvalues) );

   return SCIP_OKAY;
}

/** inserts a string array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPinsertDatatreeStringArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const char* const*    values,             /**< values of entry */
   int                   nvalues             /**< number of values */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues >= 0);

   SCIP_CALL( SCIPdatatreeInsertStringArray(datatree, scip->set, scip->mem->setmem, name, values, nvalues) );

   return SCIP_OKAY;
}

/** inserts a data tree value into a SCIP_DATATREE object
 *
 *  The data tree assumes ownership of value.
 */
SCIP_RETCODE SCIPinsertDatatreeTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_DATATREE*        value               /**< value to add */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(name != NULL);
   assert(value != NULL);

   SCIP_CALL( SCIPdatatreeInsertTree( datatree, scip->set, scip->mem->setmem, name, value ) );

   return SCIP_OKAY;
}

/** frees a SCIP_DATATREE object */
void SCIPfreeDatatree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE**       datatree            /**< pointer to data tree to free */
)
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIPdatatreeFree(datatree, scip->mem->setmem);
}

/** writes a SCIP_DATATREE object as JSON to a file */
SCIP_RETCODE SCIPwriteDatatreeJson(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file to write to, or NULL for stdout */
   SCIP_DATATREE*        datatree            /**< data tree to write */
)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPdatatreeWriteJson(datatree, scip->messagehdlr, file) );

   return SCIP_OKAY;
}

/** prints a generic table from a data tree */
SCIP_RETCODE SCIPprintDatatreeAsTable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   FILE*                 file,               /**< output file */
   const char*           sectionname,        /**< section name to process, e.g., "plugins" */
   const char*           tablename           /**< table name to process, e.g., "heuristics" */
)
{
   SCIP_DATATREE* section;
   SCIP_DATATREE* firstsectionitem = NULL;

   assert(scip != NULL);
   assert(datatree != NULL);

   /* Get the section name */
   SCIP_CALL( SCIPdatatreeGetTree( datatree, sectionname, &section ) );

   /* Print the table header */
   SCIPinfoMessage(scip, file, "%-19s:", tablename);

   if( section->nitems == 0 )
   {
      SCIPinfoMessage(scip, file, "No data available in section '%s'\n", sectionname);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPdatatreeGetTree( section, section->items[0].name, &firstsectionitem ) );

   /* Dynamically determine all keys from the first item */
   for( int i = 0; i < firstsectionitem->nitems; ++i )
   {
      size_t len = strlen(firstsectionitem->items[i].name);
      if( len > 14 )
      {
         SCIPinfoMessage(scip, file, " %-15.14s", firstsectionitem->items[i].name);
      }
      else
      {
         SCIPinfoMessage(scip, file, " %-15s", firstsectionitem->items[i].name);
      }
   }
   SCIPinfoMessage(scip, file, "\n");

   /* Print each item's data */
   for( int i = 0; i < section->nitems; ++i )
   {
      SCIP_DATATREE* sectionitem = NULL;
      SCIPinfoMessage(scip, file, "%-19s:", section->items[i].name);

      SCIP_CALL( SCIPdatatreeGetTree(section, section->items[i].name, &sectionitem) );
      assert(sectionitem != NULL);
      assert(sectionitem->nitems == firstsectionitem->nitems);

      /* Print the values of all keys for this item */
      for( int j = 0; j < sectionitem->nitems; ++j )
      {
         SCIP_DATATREEITEM* item = &sectionitem->items[j];

         assert(strcmp(item->name, firstsectionitem->items[j].name) == 0);

         switch( item->value.type )
         {
            case SCIP_DATATREE_BOOL:
            {
               SCIP_Bool boolval;
               SCIP_CALL( SCIPdatatreeGetBool(sectionitem, item->name, &boolval) );
               SCIPinfoMessage(scip, file, " %-15s", boolval ? "true" : "false");
               break;
            }
            case SCIP_DATATREE_LONG:
            {
               SCIP_Longint longval;
               SCIP_CALL( SCIPdatatreeGetLong(sectionitem, item->name, &longval) );
               SCIPinfoMessage(scip, file, " %-15" SCIP_LONGINT_FORMAT, longval);
               break;
            }
            case SCIP_DATATREE_REAL:
            {
               SCIP_Real realval;
               SCIP_CALL( SCIPdatatreeGetReal(sectionitem, item->name, &realval) );
               SCIPinfoMessage(scip, file, " %-15.2f", realval);
               break;
            }
            case SCIP_DATATREE_STRING:
            {
               const char* stringval;
               SCIP_CALL( SCIPdatatreeGetString(sectionitem, item->name, &stringval) );
               size_t len = strlen(stringval);
               if( len > 14 )
               {
                  SCIPinfoMessage(scip, file, " %-15.14s", stringval);
               }
               else
               {
                  SCIPinfoMessage(scip, file, " %-15s", stringval);
               }
               break;
            }
            case SCIP_DATATREE_DATATREE:
               SCIPwarningMessage(scip, "Table value type not supported in SCIPprintDatatreeAsTable\n");
               break;
            case SCIP_DATATREE_BOOLARRAY:
            case SCIP_DATATREE_LONGARRAY:
            case SCIP_DATATREE_REALARRAY:
            case SCIP_DATATREE_STRINGARRAY:
               SCIPwarningMessage(scip, "Table array value type not supported in SCIPprintDatatreeAsTable\n");
               break;
            default:
               SCIPinfoMessage(scip, file, " %-15s", "N/A");
               break;
         }
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   return SCIP_OKAY;
}
