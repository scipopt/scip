/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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
   DECL_READERFREE((*readerfree));      /**< destructor of reader */
   DECL_READERINIT((*readerinit));      /**< initialise reader */
   DECL_READEREXIT((*readerexit));      /**< deinitialise reader */
   DECL_READERREAD((*readerread));      /**< read method */
   READERDATA*      readerdata;         /**< reader data */
   unsigned int     initialized:1;      /**< is reader initialized? */
};



/* reader methods */

RETCODE SCIPreaderCreate(               /**< creates a reader */
   READER**         reader,             /**< pointer to store reader */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE((*readerfree)),      /**< destructor of reader */
   DECL_READERINIT((*readerinit)),      /**< initialise reader */
   DECL_READEREXIT((*readerexit)),      /**< deinitialise reader */
   DECL_READERREAD((*readerread)),      /**< read method */
   READERDATA*      readerdata          /**< reader data */
   )
{
   assert(reader != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(extension != NULL);
   assert(readerread != NULL);

   ALLOC_OKAY( allocMemory(*reader) );
   ALLOC_OKAY( duplicateMemoryArray((*reader)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*reader)->desc, desc, strlen(desc)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*reader)->extension, extension, strlen(extension)+1) );
   (*reader)->readerfree = readerfree;
   (*reader)->readerinit = readerinit;
   (*reader)->readerexit = readerexit;
   (*reader)->readerread = readerread;
   (*reader)->readerdata = readerdata;
   (*reader)->initialized = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPreaderFree(                 /**< frees memory of reader */
   READER**         reader,             /**< pointer to reader data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(reader != NULL);
   assert(*reader != NULL);
   assert(!(*reader)->initialized);

   /* call destructor of reader */
   if( (*reader)->readerfree != NULL )
   {
      CHECK_OKAY( (*reader)->readerfree(*reader, scip) );
   }

   freeMemoryArray((*reader)->name);
   freeMemoryArray((*reader)->desc);
   freeMemoryArray((*reader)->extension);
   freeMemory(*reader);

   return SCIP_OKAY;
}

RETCODE SCIPreaderInit(                 /**< initializes reader */
   READER*          reader,             /**< reader */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(reader != NULL);
   assert(scip != NULL);

   if( reader->initialized )
   {
      char s[255];
      sprintf(s, "reader <%s> already initialized", reader->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( reader->readerinit != NULL )
   {
      CHECK_OKAY( reader->readerinit(reader, scip) );
   }
   reader->initialized = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPreaderExit(                 /**< deinitializes reader */
   READER*          reader,             /**< reader */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(reader != NULL);
   assert(scip != NULL);

   if( !reader->initialized )
   {
      char s[255];
      sprintf(s, "reader <%s> not initialized", reader->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( reader->readerexit != NULL )
   {
      CHECK_OKAY( reader->readerexit(reader, scip) );
   }
   reader->initialized = FALSE;

   return SCIP_OKAY;
}

static
void splitFilename(                     /**< splits filename into path, name, and extension */
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

static
Bool readerIsApplicable(                /**< returns TRUE, if reader is responsible for files with the given extension */
   READER*          reader,             /**< reader */
   const char*      extension           /**< extension of the input file name */
   )
{
   assert(reader != NULL);
   assert(reader->extension != NULL);
   assert(extension != NULL);

   return (strcasecmp(reader->extension, extension) == 0);
}

RETCODE SCIPreaderRead(                 /**< reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
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
   ALLOC_OKAY( duplicateMemoryArray(tmpfilename, filename, strlen(filename)+1) );
   splitFilename(tmpfilename, &path, &name, &extension);

   /* check, if reader is applicable on the given file */
   if( readerIsApplicable(reader, extension) )
   {
      /* create new problem */
      CHECK_OKAY( SCIPcreateProb(scip, name) );

      /* call reader to read problem */
      CHECK_OKAY( reader->readerread(reader, scip, filename, result) );
   }
   else
      *result = SCIP_DIDNOTRUN;

   freeMemoryArray(tmpfilename);

   return SCIP_OKAY;
}

const char* SCIPreaderGetName(          /**< gets name of reader */
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->name;
}

READERDATA* SCIPreaderGetData(          /**< gets user data of reader */
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->readerdata;
}

void SCIPreaderSetData(                 /**< sets user data of reader; user has to free old data in advance! */
   READER*          reader,             /**< reader */
   READERDATA*      readerdata          /**< new reader user data */
   )
{
   assert(reader != NULL);

   reader->readerdata = readerdata;
}

Bool SCIPreaderIsInitialized(           /**< is reader initialized? */
   READER*          reader              /**< reader */
   )
{
   assert(reader != NULL);

   return reader->initialized;
}
