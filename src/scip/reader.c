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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader.c,v 1.24 2005/02/14 13:35:49 bzfpfend Exp $"

/**@file   reader.c
 * @brief  interface for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/reader.h"

#include "scip/struct_reader.h"



/** creates a reader */
RETCODE SCIPreaderCreate(
   READER**         reader,             /**< pointer to store reader */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   DECL_READERREAD  ((*readerread)),    /**< read method */
   READERDATA*      readerdata          /**< reader data */
   )
{
   assert(reader != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(extension != NULL);
   assert(readerread != NULL);

   ALLOC_OKAY( allocMemory(reader) );
   ALLOC_OKAY( duplicateMemoryArray(&(*reader)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*reader)->desc, desc, strlen(desc)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*reader)->extension, extension, strlen(extension)+1) );
   (*reader)->readerfree = readerfree;
   (*reader)->readerread = readerread;
   (*reader)->readerdata = readerdata;

   return SCIP_OKAY;
}

/** frees memory of reader */
RETCODE SCIPreaderFree(
   READER**         reader,             /**< pointer to reader data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(reader != NULL);
   assert(*reader != NULL);
   assert(set != NULL);

   /* call destructor of reader */
   if( (*reader)->readerfree != NULL )
   {
      CHECK_OKAY( (*reader)->readerfree(set->scip, *reader) );
   }

   freeMemoryArray(&(*reader)->name);
   freeMemoryArray(&(*reader)->desc);
   freeMemoryArray(&(*reader)->extension);
   freeMemory(reader);

   return SCIP_OKAY;
}

/** returns TRUE, if reader is responsible for files with the given extension */
static
Bool readerIsApplicable(
   READER*          reader,             /**< reader */
   const char*      extension           /**< extension of the input file name */
   )
{
   assert(reader != NULL);
   assert(reader->extension != NULL);

   return (extension != NULL && strcasecmp(reader->extension, extension) == 0)
      || (extension == NULL && *(reader->extension) == '\0');
}

/** reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
RETCODE SCIPreaderRead(
   READER*          reader,             /**< reader */
   SET*             set,                /**< global SCIP settings */
   const char*      filename,           /**< name of the input file */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   char* tmpfilename;
   char* path;
   char* name;
   char* extension;

   assert(reader != NULL);
   assert(reader->readerread != NULL);
   assert(set != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   /* get path, name and extension from filename */
   ALLOC_OKAY( duplicateMemoryArray(&tmpfilename, filename, strlen(filename)+1) );
   SCIPsplitFilename(tmpfilename, &path, &name, &extension);

   /* check, if reader is applicable on the given file */
   if( readerIsApplicable(reader, extension) )
   {
      /* call reader to read problem */
      CHECK_OKAY( reader->readerread(set->scip, reader, filename, result) );
   }
   else
      *result = SCIP_DIDNOTRUN;

   freeMemoryArray(&tmpfilename);

   return SCIP_OKAY;
}

/** gets user data of reader */
READERDATA* SCIPreaderGetData(
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->readerdata;
}

/** sets user data of reader; user has to free old data in advance! */
void SCIPreaderSetData(
   READER*          reader,             /**< reader */
   READERDATA*      readerdata          /**< new reader user data */
   )
{
   assert(reader != NULL);

   reader->readerdata = readerdata;
}

/** gets name of reader */
const char* SCIPreaderGetName(
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->name;
}

/** gets description of reader */
const char* SCIPreaderGetDesc(
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->desc;
}

/** gets file extension of reader */
const char* SCIPreaderGetExtension(
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->extension;
}

