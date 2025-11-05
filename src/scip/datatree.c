/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file    datatree.c
 * @ingroup OTHER_CFILES
 * @brief  methods for managing data trees
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/datatree.h"
#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/pub_datatree.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/syncstore.h"
#include "scip/struct_datatree.h"


/** default capacity of a table store */
#define DATATREE_DEFAULT_CAPACITY 4

/** ensure a SCIP_DATATREE object can store an additional item
 *
 * This function ensures that a SCIP_DATATREE object can store an additional item
 * by increasing the itemssize, if necessary.
 */
static
SCIP_RETCODE datatreeExpand(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int newsize;

   assert(datatree != NULL);
   assert(datatree->nitems <= datatree->itemssize);

   if( datatree->nitems < datatree->itemssize )
      return SCIP_OKAY;

   newsize = SCIPsetCalcMemGrowSize(set, datatree->itemssize + 1);

   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &datatree->items, datatree->itemssize, newsize) );
   datatree->itemssize = newsize;

   return SCIP_OKAY;
}

/** gives name to table valuetype enum */
static
const char* datatreeValueTypeToString(
   SCIP_DATATREE_VALUETYPE type
   )
{
   switch( type )
   {
      case SCIP_DATATREE_BOOL:
         return "bool";
      case SCIP_DATATREE_LONG:
         return "long";
      case SCIP_DATATREE_REAL:
         return "real";
      case SCIP_DATATREE_STRING:
         return "string";
      case SCIP_DATATREE_BOOLARRAY:
         return "boolarray";
      case SCIP_DATATREE_LONGARRAY:
         return "longarray";
      case SCIP_DATATREE_REALARRAY:
         return "realarray";
      case SCIP_DATATREE_STRINGARRAY:
         return "stringarray";
      case SCIP_DATATREE_DATATREE:
         return "datatree";
      default:
         return "unknown";
   }
}

/** creates a new SCIP_DATATREE with a given capacity for items */
SCIP_RETCODE SCIPdatatreeCreate(
   SCIP_DATATREE**       datatree,           /**< buffer to store pointer to created datatree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   capacity            /**< initial capacity */
   )
{
   assert(datatree != NULL);
   assert(blkmem != NULL);
   assert(capacity == -1 || capacity > 0);

   if( capacity == -1 )
      capacity = DATATREE_DEFAULT_CAPACITY;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, datatree) );
   (*datatree)->nitems = 0;
   (*datatree)->itemssize = capacity;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*datatree)->items, capacity) );

   return SCIP_OKAY;
}

/** frees a SCIP_DATATREE object */
void SCIPdatatreeFree(
   SCIP_DATATREE**       datatree,           /**< pointer to datatree to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(datatree != NULL);
   assert(*datatree != NULL);

   for( i = 0; i < (*datatree)->nitems; ++i )
   {
      SCIP_DATATREEITEM* item = &(*datatree)->items[i];
      switch( item->value.type )
      {
         case SCIP_DATATREE_BOOL:
         case SCIP_DATATREE_LONG:
         case SCIP_DATATREE_REAL:
            break;
         case SCIP_DATATREE_STRING:
            BMSfreeBlockMemoryArray(blkmem, &item->value.data.as_string, strlen(item->value.data.as_string) + 1);
            break;
         case SCIP_DATATREE_BOOLARRAY:
            BMSfreeBlockMemoryArray(blkmem, &item->value.data.as_boolarray, item->value.nvalues);
            break;
         case SCIP_DATATREE_LONGARRAY:
            BMSfreeBlockMemoryArray(blkmem, &item->value.data.as_longarray, item->value.nvalues);
            break;
         case SCIP_DATATREE_REALARRAY:
            BMSfreeBlockMemoryArray(blkmem, &item->value.data.as_realarray, item->value.nvalues);
            break;
         case SCIP_DATATREE_STRINGARRAY:
         {
            int j;
            for( j = 0; j < item->value.nvalues; ++j )
            {
               assert(item->value.data.as_stringarray[j] != NULL);
               BMSfreeBlockMemoryArray(blkmem, &item->value.data.as_stringarray[j], strlen(item->value.data.as_stringarray[j]) + 1);
            }
            BMSfreeBlockMemoryArray(blkmem, &item->value.data.as_stringarray, item->value.nvalues);
            break;
         }
         case SCIP_DATATREE_DATATREE:
            SCIPdatatreeFree(&item->value.data.as_dtree, blkmem);
            break;
         default:
            SCIPerrorMessage("Unknown type\n");
            SCIPABORT();
      }
      BMSfreeBlockMemoryArray(blkmem, &item->name, strlen(item->name) + 1);
   }

   BMSfreeBlockMemoryArray(blkmem, &(*datatree)->items, (*datatree)->itemssize);
   BMSfreeBlockMemory(blkmem, datatree);
}

/** inserts a SCIP_Bool value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertBool(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_Bool             value               /**< value of entry */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   item->value.type = SCIP_DATATREE_BOOL;
   item->value.data.as_bool = value;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a long value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertLong(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_Longint          value               /**< value of entry */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   item->value.type = SCIP_DATATREE_LONG;
   item->value.data.as_long = value;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a SCIP_Real value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertReal(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_Real             value               /**< value of entry */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   item->value.type = SCIP_DATATREE_REAL;
   item->value.data.as_real = value;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a string value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertString(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const char*           value               /**< value of entry */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   item->value.type = SCIP_DATATREE_STRING;
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->value.data.as_string, value, strlen(value) + 1) );

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a SCIP_Bool array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertBoolArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const SCIP_Bool*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues > 0);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->value.data.as_boolarray, values, nvalues) );
   item->value.nvalues = nvalues;
   item->value.type = SCIP_DATATREE_BOOLARRAY;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a SCIP_Longint array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertLongArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const SCIP_Longint*   values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues > 0);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->value.data.as_longarray, values, nvalues) );
   item->value.nvalues = nvalues;
   item->value.type = SCIP_DATATREE_LONGARRAY;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a SCIP_Real array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertRealArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const SCIP_Real*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues > 0);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->value.data.as_realarray, values, nvalues) );
   item->value.nvalues = nvalues;
   item->value.type = SCIP_DATATREE_REALARRAY;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a string array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertStringArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const char* const*    values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   )
{
   SCIP_DATATREEITEM* item;
   int i;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues > 0);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &item->value.data.as_stringarray, nvalues) );
   for( i = 0; i < nvalues; ++i )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->value.data.as_stringarray[i], values[i], strlen(values[i]) + 1) );
   }
   item->value.nvalues = nvalues;
   item->value.type = SCIP_DATATREE_STRINGARRAY;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** inserts a datatree value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertTree(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_DATATREE*        value               /**< value of entry */
   )
{
   SCIP_DATATREEITEM* item;

   assert(datatree != NULL);
   assert(name != NULL);

   /* allocate memory for the new item */
   SCIP_CALL( datatreeExpand(datatree, set, blkmem) );

   item = &datatree->items[datatree->nitems];

   /* copy the name */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &item->name, name, strlen(name) + 1) );

   /* set the value */
   item->value.type = SCIP_DATATREE_DATATREE;
   item->value.data.as_dtree = value;

   datatree->nitems++;

   return SCIP_OKAY;
}

/** writes a SCIP_DATATREE object as JSON to file */
SCIP_RETCODE SCIPdatatreeWriteJson(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to write to, or NULL for stdout */
   )
{
   int i;

   assert(datatree != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file, "{");
   for( i = 0; i < datatree->nitems; ++i )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "\"%s\": ", datatree->items[i].name);

      switch( datatree->items[i].value.type )
      {
         case SCIP_DATATREE_BOOL:
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%s", datatree->items[i].value.data.as_bool ? "true" : "false");
            break;
         }
         case SCIP_DATATREE_LONG:
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%" SCIP_LONGINT_FORMAT, datatree->items[i].value.data.as_long);
            break;
         }
         case SCIP_DATATREE_REAL:
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%f", datatree->items[i].value.data.as_real);
            break;
         }
         case SCIP_DATATREE_STRING:
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "\"%s\"", datatree->items[i].value.data.as_string);
            break;
         }
         case SCIP_DATATREE_DATATREE:
         {
            SCIP_CALL( SCIPdatatreeWriteJson( datatree->items[i].value.data.as_dtree, messagehdlr, file) );
            break;
         }
         case SCIP_DATATREE_BOOLARRAY:
         case SCIP_DATATREE_LONGARRAY:
         case SCIP_DATATREE_REALARRAY:
         case SCIP_DATATREE_STRINGARRAY:
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "[");
            for( int j = 0; j < datatree->items[i].value.nvalues; ++j )
            {
               switch( datatree->items[i].value.type )
               {
                  case SCIP_DATATREE_BOOLARRAY:
                  {
                     SCIPmessageFPrintInfo(messagehdlr, file, "%s", datatree->items[i].value.data.as_boolarray[j] ? "true" : "false");
                     break;
                  }

                  case SCIP_DATATREE_LONGARRAY:
                  {
                     SCIPmessageFPrintInfo(messagehdlr, file, "%" SCIP_LONGINT_FORMAT, datatree->items[i].value.data.as_longarray[j]);
                     break;
                  }

                  case SCIP_DATATREE_REALARRAY:
                  {
                     SCIPmessageFPrintInfo(messagehdlr, file, "%f", datatree->items[i].value.data.as_realarray[j]);
                     break;
                  }

                  case SCIP_DATATREE_STRINGARRAY:
                  {
                     SCIPmessageFPrintInfo(messagehdlr, file, "\"%s\"", datatree->items[i].value.data.as_stringarray[j]);
                     break;
                  }

                  default:
                  {
                     SCIPABORT();
                     return SCIP_ERROR;
                  }
               }
               if( j < datatree->items[i].value.nvalues - 1 )
               {
                  SCIPmessageFPrintInfo(messagehdlr, file, ", ");
               }
            }
            SCIPmessageFPrintInfo(messagehdlr, file, "]");
            break;
         }
         default:
         {
            SCIPerrorMessage("Unknown SCIP_TABLEVALUETYPE\n");
            return SCIP_ERROR;
         }
      }

      if( i < datatree->nitems - 1 )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, ", ");
      }
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "}");

   return SCIP_OKAY;
}

/** gets a bool value from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetBool(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Bool*            value               /**< buffer to store value */
   )
{
   int i;

   assert(datatree != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_BOOL )
         {
            SCIPerrorMessage("Value for key %s is not of type bool, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }
         *value = datatree->items[i].value.data.as_bool;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a long value from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetLong(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Longint*         value               /**< buffer to store value */
   )
{
   int i;

   assert(datatree != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_LONG )
         {
            SCIPerrorMessage("Value for key %s is not of type long, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }
         *value = datatree->items[i].value.data.as_long;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a SCIP_Real value from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetReal(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Real*            value               /**< buffer to store value */
   )
{
   int i;

   assert(datatree != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_REAL )
         {
            SCIPerrorMessage("Value for key %s is not of type real, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }
         *value = datatree->items[i].value.data.as_real;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a string value from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetString(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   const char**          value               /**< buffer to store pointer to string */
   )
{
   int i;

   assert(datatree != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_STRING )
         {
            SCIPerrorMessage("Value for key %s is not of type string, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }
         *value = datatree->items[i].value.data.as_string;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a bool array from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetBoolArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Bool**           values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   )
{
   int i;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_BOOLARRAY )
         {
            SCIPerrorMessage("Value for key %s is not of type bool array, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }

         *values = datatree->items[i].value.data.as_boolarray;
         *nvalues = datatree->items[i].value.nvalues;

         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a SCIP_Longint array from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetLongArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Longint**        values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   )
{
   int i;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_LONGARRAY )
         {
            SCIPerrorMessage("Value for key %s is not of type long array, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }

         *values = datatree->items[i].value.data.as_longarray;
         *nvalues = datatree->items[i].value.nvalues;

         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a SCIP_Real array from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetRealArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Real**           values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   )
{
   int i;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_REALARRAY )
         {
            SCIPerrorMessage("Value for key %s is not of type real array, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }

         *values = datatree->items[i].value.data.as_realarray;
         *nvalues = datatree->items[i].value.nvalues;

         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}


/** gets a string array from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetStringArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   char***               values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   )
{
   int i;

   assert(datatree != NULL);
   assert(name != NULL);
   assert(values != NULL);
   assert(nvalues != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_STRINGARRAY )
         {
            SCIPerrorMessage("Value for key %s is not of type string array, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }

         *values = datatree->items[i].value.data.as_stringarray;
         *nvalues = datatree->items[i].value.nvalues;

         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}

/** gets a datatree value from a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeGetTree(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_DATATREE**       value               /**< buffer to store pointer to datatree */
   )
{
   int i;

   assert(datatree != NULL);

   for( i = 0; i < datatree->nitems; ++i )
   {
      if( strcmp(datatree->items[i].name, name) == 0 )
      {
         if( datatree->items[i].value.type != SCIP_DATATREE_DATATREE )
         {
            SCIPerrorMessage("Value for key %s is not of type datatree, but of type %s\n", name, datatreeValueTypeToString(datatree->items[i].value.type));
            return SCIP_ERROR;
         }
         *value = datatree->items[i].value.data.as_dtree;
         return SCIP_OKAY;
      }
   }

   SCIPerrorMessage("Key %s not found\n", name);
   return SCIP_ERROR;
}
