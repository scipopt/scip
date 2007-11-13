/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader.c,v 1.37 2007/11/13 17:21:48 bzfheinz Exp $"

/**@file   reader.c
 * @brief  interface for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h>
#endif

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/reader.h"
#include "scip/prob.h"
#include "scip/pub_var.h"
#include "scip/var.h"
#include "scip/pub_cons.h"
#include "scip/cons.h"

#include "scip/struct_reader.h"



/** creates a reader */
SCIP_RETCODE SCIPreaderCreate(
   SCIP_READER**         reader,             /**< pointer to store reader */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   SCIP_DECL_READERREAD  ((*readerread)),    /**< read method */
   SCIP_DECL_READERWRITE ((*readerwrite)),   /**< write method */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   assert(reader != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(extension != NULL);

   SCIP_ALLOC( BMSallocMemory(reader) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*reader)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*reader)->desc, desc, strlen(desc)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*reader)->extension, extension, strlen(extension)+1) );
   (*reader)->readerfree = readerfree;
   (*reader)->readerread = readerread;
   (*reader)->readerwrite = readerwrite;
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
   SCIP_RETCODE retcode;
   char* tmpfilename;
   char* path;
   char* name;
   char* extension;
   char* compression;

   assert(reader != NULL);
   assert(set != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   /* get path, name and extension from filename */
   SCIP_ALLOC( BMSduplicateMemoryArray(&tmpfilename, filename, strlen(filename)+1) );
   SCIPsplitFilename(tmpfilename, &path, &name, &extension, &compression);

   /* check, if reader is applicable on the given file */
   if( readerIsApplicable(reader, extension) && reader->readerread != NULL )
   {
      /* call reader to read problem */
      retcode = reader->readerread(set->scip, reader, filename, result);
   }
   else
   {
      *result = SCIP_DIDNOTRUN;
      retcode = SCIP_OKAY;
   }

   BMSfreeMemoryArray(&tmpfilename);

   /* check for reader errors */
   if( retcode == SCIP_READERROR || retcode == SCIP_NOFILE || retcode == SCIP_PARSEERROR )
      return retcode;

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}


/* reset the variable name to the given one */
static
void resetVarname(
   SCIP_VAR*             var,                /**< variable */
   const char*           name                /**< variable name */
   )
{
   const char * oldname;

   assert( var != NULL );
   assert( name != NULL );
   
   /* get pointer to temporary generic name and free the memory */
   oldname = SCIPvarGetName(var);
   BMSfreeMemory(&oldname);
   
   /* reset name */
   SCIPvarSetNamePointer(var, name);
}


/** writes problem data to file with given reader or returns SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPreaderWrite(
   SCIP_READER*          reader,             /**< reader */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format */
   SCIP_Bool             genericnames,       /**< using generic variable and constraint names? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_RETCODE retcode;
   int i;
   
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR** fixedvars;
   int nfixedvars;
   SCIP_CONS** conss;
   int nconss;

   SCIP_VAR* var;
   SCIP_CONS* cons;

   int size;
   char* name;
   const char* consname;
   const char** varnames;
   const char** fixedvarnames;
   const char** consnames;
   
   assert(reader != NULL);
   assert(set != NULL);
   assert(extension != NULL);
   assert(result != NULL);

   vars = prob->vars;
   nvars = prob->nvars;
   fixedvars = prob->fixedvars;
   nfixedvars = prob->nfixedvars;
   conss = prob->conss;
   nconss = prob->nconss;
   
   /* check, if reader is applicable on the given file */
   if( readerIsApplicable(reader, extension) && reader->readerwrite != NULL )
   {
      if( genericnames )
      {
         /* save variable and constraint names and replace these names by generic names */

         /* alloac memory for saving the original variable and constraint names */
         SCIP_ALLOC( BMSallocMemoryArray(&varnames, nvars ) );
         SCIP_ALLOC( BMSallocMemoryArray(&fixedvarnames, nfixedvars ) );
         SCIP_ALLOC( BMSallocMemoryArray(&consnames, nconss ) );

         /* compute length of the generic variable names */
         size = nvars % 10 + 2;
         
         for( i = 0; i < nvars; ++i )
         {
            var = vars[i];
            varnames[i] = SCIPvarGetName(var);
            
            SCIP_ALLOC( BMSallocMemoryArray(&name, size ) );
            sprintf(name, "x%d", i);
            SCIPvarSetNamePointer(var, name);
         }  
         
         /* compute length of the generic variable names */
         size = nfixedvars % 10 + 2;

         for( i = 0; i < nfixedvars; ++i )
         {
            var = fixedvars[i];
            fixedvarnames[i] = SCIPvarGetName(var);
            
            SCIP_ALLOC( BMSallocMemoryArray(&name, size ) );
            sprintf(name, "y%d", i);
            SCIPvarSetNamePointer(var, name);
         }
         
         /* compute length of the generic constraint names */
         size = nconss % 10 + 2;
         
         for( i = 0; i < nconss; ++i )
         {
            cons = conss[i];
            consnames[i] = SCIPconsGetName(cons);

            SCIP_ALLOC( BMSallocMemoryArray(&name, size ) );
            sprintf(name, "c%d", i);
            SCIPconsSetNamePointer(cons, name);
         }
      }

      /* call reader to write problem */
      retcode = reader->readerwrite(set->scip, reader, file, prob->name, prob->probdata, prob->transformed,
         prob->transformed ? SCIP_OBJSENSE_MINIMIZE : prob->objsense, prob->objscale, prob->objoffset,
         vars, nvars, prob->nbinvars, prob->nintvars, prob->nimplvars, prob->ncontvars, 
         fixedvars, nfixedvars, prob->startnvars, 
         conss, nconss, prob->maxnconss, prob->startnconss, genericnames, result);
         
      if( genericnames )
      {
         /* reset variable and constraint names to original names */
         
         for( i = 0; i < nvars; ++i )
            resetVarname(vars[i], varnames[i]);
               
         for( i = 0; i < nfixedvars; ++i )
            resetVarname(fixedvars[i], fixedvarnames[i]);

         for( i = 0; i < nconss; ++i )
         {
            cons = conss[i];

            /* get pointer to temporary generic name and free the memory */
            consname = SCIPconsGetName(cons);
            BMSfreeMemory(&consname);
            
            /* reset name */
            SCIPconsSetNamePointer(cons, consnames[i]);
         }
         
         /* free memory */
         BMSfreeMemoryArray(&varnames);
         BMSfreeMemoryArray(&fixedvarnames);
         BMSfreeMemoryArray(&consnames);
      }
   }
   else
   {
      *result = SCIP_DIDNOTRUN;
      retcode = SCIP_OKAY;
   }
   
   /* check for reader errors */
   if( retcode == SCIP_WRITEERROR )
      return retcode;
   
   SCIP_CALL( retcode );

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

