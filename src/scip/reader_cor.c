/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_cor.c
 * @brief  COR file reader (MPS format of the core problem for stochastic programs)
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/reader_mps.h"
#include "scip/reader_cor.h"

#define READER_NAME             "correader"
#define READER_DESC             "file reader for core problem of stochastic programs in the COR file format"
#define READER_EXTENSION        "cor"

#define DEFAULT_LINEARIZE_ANDS         TRUE  /**< should possible \"and\" constraint be linearized when writing the mps file? */
#define DEFAULT_AGGRLINEARIZATION_ANDS TRUE  /**< should an aggregated linearization for and constraints be used? */

#define SCIP_DEFAULT_ARRAYSIZE      100

/** COR reading data */
struct SCIP_ReaderData
{
   const char**          varnames;
   const char**          consnames;
   int                   varnamessize;
   int                   consnamessize;
   int                   nvarnames;
   int                   nconsnames;
   SCIP_Bool             read;

   /* Required for the MPS reading and writing functions */
   SCIP_Bool             linearizeands;
   SCIP_Bool             aggrlinearizationands;
};

/** creates the reader data  */
static
SCIP_RETCODE createReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data structure */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);

   readerdata->read = FALSE;
   readerdata->nvarnames = 0;
   readerdata->nconsnames = 0;
   readerdata->varnamessize = SCIP_DEFAULT_ARRAYSIZE;
   readerdata->consnamessize = SCIP_DEFAULT_ARRAYSIZE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->varnames, readerdata->varnamessize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->consnames, readerdata->consnamessize) );

   return SCIP_OKAY;
}

/** creates the reader data  */
static
SCIP_RETCODE freeReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(readerdata != NULL);

   for( i = readerdata->nvarnames - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &readerdata->varnames[i], strlen(readerdata->varnames[i]) + 1);

   for( i = readerdata->nconsnames - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &readerdata->consnames[i], strlen(readerdata->consnames[i]) + 1);

   SCIPfreeBlockMemoryArray(scip, &readerdata->consnames, readerdata->consnamessize);
   SCIPfreeBlockMemoryArray(scip, &readerdata->varnames, readerdata->varnamessize);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderCor(scip) );

   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
/**! [SnippetReaderFreeCor] */
static
SCIP_DECL_READERFREE(readerFreeCor)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( freeReaderdata(scip, readerdata) );

   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}
/**! [SnippetReaderFreeCor] */


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCor)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadCor(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCor)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPwriteCor(scip, reader, file, name, transformed, objsense, objscale, objoffset, vars, nvars, nbinvars,
         nintvars, nimplvars, ncontvars, fixedvars, nfixedvars, conss, nconss, genericnames, result) );

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the cor file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );
   SCIP_CALL( createReaderdata(scip, readerdata) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCor) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeCor) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCor) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCor) );

   /* add lp-reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/linearize-and-constraints",
         "should possible \"and\" constraint be linearized when writing the mps file?",
         &readerdata->linearizeands, TRUE, DEFAULT_LINEARIZE_ANDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/aggrlinearization-ands",
         "should an aggregated linearization for and constraints be used?",
         &readerdata->aggrlinearizationands, TRUE, DEFAULT_AGGRLINEARIZATION_ANDS, NULL, NULL) );

   return SCIP_OKAY;
}



/** reads problem from file */
SCIP_RETCODE SCIPreadCor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPreadMps(scip, reader, filename, result, &readerdata->varnames, &readerdata->consnames,
         &readerdata->varnamessize, &readerdata->consnamessize, &readerdata->nvarnames, &readerdata->nconsnames) );

   if( (*result) == SCIP_SUCCESS )
      readerdata->read = TRUE;

   return SCIP_OKAY;
}

/** writes problem to file */
SCIP_RETCODE SCIPwriteCor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
                                              * extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of general integer variables */
   int                   nimplvars,          /**< number of implicit integer variables */
   int                   ncontvars,          /**< number of continuous variables */
   SCIP_VAR**            fixedvars,          /**< array with fixed and aggregated variables */
   int                   nfixedvars,         /**< number of fixed and aggregated variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_Bool             genericnames,       /**< should generic names be used for the variables and constraints? */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   if( genericnames )
   {
      SCIP_CALL( SCIPwriteMps(scip, reader, file, name, transformed, objsense, objscale, objoffset, vars,
            nvars, nbinvars, nintvars, nimplvars, ncontvars, fixedvars, nfixedvars, conss, nconss, result) );
   }
   else
   {
      SCIPwarningMessage(scip, "COR format is an MPS format with generic variable and constraint names\n");

      if( transformed )
      {
         SCIPwarningMessage(scip, "write transformed problem with generic variable and constraint names\n");
         SCIP_CALL( SCIPprintTransProblem(scip, file, "cor", TRUE) );
      }
      else
      {
         SCIPwarningMessage(scip, "write original problem with generic variable and constraint names\n");
         SCIP_CALL( SCIPprintOrigProblem(scip, file, "cor", TRUE) );
      }
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/*
 * Interface method for the tim and sto readers
 */
SCIP_Bool SCIPcorHasRead(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->read;
}

int SCIPcorGetNVarNames(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->nvarnames;
}

int SCIPcorGetNConsNames(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->nconsnames;
}

const char* SCIPcorGetVarName(
   SCIP_READER*          reader,             /**< the file reader itself */
   int                   i                   /**< the index of the constraint that is requested */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(i >= 0 && i < readerdata->nvarnames);

   return readerdata->varnames[i];
}

const char* SCIPcorGetConsName(
   SCIP_READER*          reader,             /**< the file reader itself */
   int                   i                   /**< the index of the constraint that is requested */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(i >= 0 && i < readerdata->nconsnames);

   return readerdata->consnames[i];
}
