/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader.c
 * @brief  interface for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader.h"



/** input file reader */
struct Reader
{
   const char*      name;               /**< name of reader */
   const char*      desc;               /**< description of reader */
   const char*      extension;          /**< file extension that reader processes */
   DECL_READERFREE  ((*readerfree));    /**< destructor of reader */
   DECL_READERREAD  ((*readerread));    /**< read method */
   READERDATA*      readerdata;         /**< reader data */
};



/* reader methods */

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
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(reader != NULL);
   assert(*reader != NULL);

   /* call destructor of reader */
   if( (*reader)->readerfree != NULL )
   {
      CHECK_OKAY( (*reader)->readerfree(scip, *reader) );
   }

   freeMemoryArray(&(*reader)->name);
   freeMemoryArray(&(*reader)->desc);
   freeMemoryArray(&(*reader)->extension);
   freeMemory(reader);

   return SCIP_OKAY;
}

/** splits filename into path, name, and extension */
static
void splitFilename(
   char*            filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**           path,               /**< pointer to store path */
   char**           name,               /**< pointer to store name */
   char**           extension           /**< pointer to store extension */
   )
{
   char* lastslash;
   char* lastdot;

   assert(filename != NULL);
   assert(path != NULL);
   assert(name != NULL);
   assert(extension != NULL);

   lastslash = strrchr(filename, '/');
   if( lastslash == NULL )
   {
      *path = NULL;
      *name = filename;
   }
   else
   {
      *path = filename;
      *name = lastslash+1;
      *lastslash = '\0';
   }

   lastdot = strrchr(*name, '.');
   if( lastdot == NULL )
      *extension = NULL;
   else
   {
      *extension = lastdot+1;
      *lastdot = '\0';
   }
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
   assert(extension != NULL);

   return (strcasecmp(reader->extension, extension) == 0);
}

/** reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
RETCODE SCIPreaderRead(
   READER*          reader,             /**< reader */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      filename,           /**< name of the input file */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   RETCODE retcode;
   char* tmpfilename;
   char* path;
   char* name;
   char* extension;

   assert(reader != NULL);
   assert(reader->readerread != NULL);
   assert(scip != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   /* get path, name and extension from filename */
   ALLOC_OKAY( duplicateMemoryArray(&tmpfilename, filename, strlen(filename)+1) );
   splitFilename(tmpfilename, &path, &name, &extension);

   /* check, if reader is applicable on the given file */
   if( readerIsApplicable(reader, extension) )
   {
      /* call reader to read problem */
      CHECK_OKAY( reader->readerread(scip, reader, filename, result) );
   }
   else
      *result = SCIP_DIDNOTRUN;

   freeMemoryArray(&tmpfilename);

   return SCIP_OKAY;
}

/** gets name of reader */
const char* SCIPreaderGetName(
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->name;
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

