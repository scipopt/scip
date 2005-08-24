/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader.c,v 1.29 2005/08/24 17:26:55 bzfpfend Exp $"

/**@file   reader.c
 * @brief  interface for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <strings.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/reader.h"

#include "scip/struct_reader.h"



/** creates a reader */
SCIP_RETCODE SCIPreaderCreate(
   SCIP_READER**         reader,             /**< pointer to store reader */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   SCIP_DECL_READERREAD  ((*readerread)),    /**< read method */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   assert(reader != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(extension != NULL);
   assert(readerread != NULL);

   SCIP_ALLOC( BMSallocMemory(reader) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*reader)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*reader)->desc, desc, strlen(desc)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*reader)->extension, extension, strlen(extension)+1) );
   (*reader)->readerfree = readerfree;
   (*reader)->readerread = readerread;
   (*reader)->readerdata = readerdata;

   return SCIP_OKAY;
}

/** frees memory of reader */
SCIP_RETCODE SCIPreaderFree(
   SCIP_READER**         reader,             /**< pointer to reader data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(reader != NULL);
   assert(*reader != NULL);
   assert(set != NULL);

   /* call destructor of reader */
   if( (*reader)->readerfree != NULL )
   {
      SCIP_CALL( (*reader)->readerfree(set->scip, *reader) );
   }

   BMSfreeMemoryArray(&(*reader)->name);
   BMSfreeMemoryArray(&(*reader)->desc);
   BMSfreeMemoryArray(&(*reader)->extension);
   BMSfreeMemory(reader);

   return SCIP_OKAY;
}

/** returns TRUE, if reader is responsible for files with the given extension */
static
SCIP_Bool readerIsApplicable(
   SCIP_READER*          reader,             /**< reader */
   const char*           extension           /**< extension of the input file name */
   )
{
   assert(reader != NULL);
   assert(reader->extension != NULL);

   return (extension != NULL && strcasecmp(reader->extension, extension) == 0)
      || (extension == NULL && *(reader->extension) == '\0');
}

/** reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPreaderRead(
   SCIP_READER*          reader,             /**< reader */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           filename,           /**< name of the input file */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   char* tmpfilename;
   char* path;
   char* name;
   char* extension;
   char* compression;

   assert(reader != NULL);
   assert(reader->readerread != NULL);
   assert(set != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   /* get path, name and extension from filename */
   SCIP_ALLOC( BMSduplicateMemoryArray(&tmpfilename, filename, strlen(filename)+1) );
   SCIPsplitFilename(tmpfilename, &path, &name, &extension, &compression);

   /* check, if reader is applicable on the given file */
   if( readerIsApplicable(reader, extension) )
   {
      /* call reader to read problem */
      SCIP_CALL( reader->readerread(set->scip, reader, filename, result) );
   }
   else
      *result = SCIP_DIDNOTRUN;

   BMSfreeMemoryArray(&tmpfilename);

   return SCIP_OKAY;
}

/** gets user data of reader */
SCIP_READERDATA* SCIPreaderGetData(
   SCIP_READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->readerdata;
}

/** sets user data of reader; user has to free old data in advance! */
void SCIPreaderSetData(
   SCIP_READER*          reader,             /**< reader */
   SCIP_READERDATA*      readerdata          /**< new reader user data */
   )
{
   assert(reader != NULL);

   reader->readerdata = readerdata;
}

/** gets name of reader */
const char* SCIPreaderGetName(
   SCIP_READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->name;
}

/** gets description of reader */
const char* SCIPreaderGetDesc(
   SCIP_READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->desc;
}

/** gets file extension of reader */
const char* SCIPreaderGetExtension(
   SCIP_READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->extension;
}

