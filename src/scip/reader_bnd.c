/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_bnd.c
 * @brief  file reader for variable bounds
 * @author Ambros Gleixner
 *
 * This reader allows to read a file containing new bounds for variables of the current problem.  Each line of the file
 * should have format
 *
 *    \<variable name\> \<lower bound\> \<upper bound\>
 *
 * where infinite bounds can be written as inf, +inf or -inf.  Note that only a subset of the variables may appear in
 * the file.  Lines with unknown variable names are ignored.  The writing functionality is currently not supported.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/reader_bnd.h"

#define READER_NAME             "bndreader"
#define READER_DESC             "file reader for variable bounds"
#define READER_EXTENSION        "bnd"


/*
 * Local methods of reader
 */

/** reads a given bound file, problem has to be in problem stage */
static
SCIP_RETCODE readBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname               /**< name of the input file */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Bool unknownvariablemessage;
   SCIP_Bool usevartable;
   int lineno;

   assert(scip != NULL);
   assert(fname != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );

   if( !usevartable )
   {
      SCIPerrorMessage("Cannot read bounds file if vartable is disabled. Make sure parameter 'misc/usevartable' is set to TRUE.\n");
      return SCIP_READERROR;
   }

   /* open input file */
   file = SCIPfopen(fname, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", fname);
      SCIPprintSysError(fname);
      return SCIP_NOFILE;
   }

   /* read the file */
   error = FALSE;
   unknownvariablemessage = FALSE;
   lineno = 0;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char lbstring[SCIP_MAXSTRLEN];
      char ubstring[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* parse the line */
      nread = sscanf(buffer, "%s %s %s\n", varname, lbstring, ubstring);
      if( nread < 2 )
      {
         SCIPerrorMessage("invalid input line %d in bounds file <%s>: <%s>\n", lineno, fname, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPwarningMessage(scip, "unknown variable <%s> in line %d of bounds file <%s>\n", varname, lineno, fname);
            SCIPwarningMessage(scip, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the lower bound value */
      if( strncasecmp(lbstring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(lbstring, "+inf", 4) == 0 || strncasecmp(lbstring, "inf", 3) == 0 )
         lb = SCIPinfinity(scip);
      else if( strncasecmp(lbstring, "-inf", 4) == 0 )
         lb = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(lbstring, "%lf", &lb);
         if( nread != 1 )
         {
            SCIPerrorMessage("invalid lower bound value <%s> for variable <%s> in line %d of bounds file <%s>\n",
               lbstring, varname, lineno, fname);
            error = TRUE;
            break;
         }
      }

      /* cast the upper bound value */
      if( strncasecmp(ubstring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(ubstring, "+inf", 4) == 0 || strncasecmp(ubstring, "inf", 3) == 0 )
         ub = SCIPinfinity(scip);
      else if( strncasecmp(ubstring, "-inf", 4) == 0 )
         ub = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(ubstring, "%lf", &ub);
         if( nread != 1 )
         {
            SCIPerrorMessage("invalid lower bound value <%s> for variable <%s> in line %d of bounds file <%s>\n",
               ubstring, varname, lineno, fname);
            error = TRUE;
            break;
         }
      }

      /* set the bounds of the variable */
      if( lb <= SCIPvarGetUbGlobal(var) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
         SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
         SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      }
   }

   /* close input file */
   SCIPfclose(file);

   /* return error if necessary */
   if ( error )
      return SCIP_READERROR;

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyBnd)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderBnd(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader
 *
 *  In order to determine the type of the file, we have to open it. Thus, it has to be opened
 *  twice. This might be removed, but is likely to not hurt the performance too much.
 */
static
SCIP_DECL_READERREAD(readerReadBnd)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of bounds file is only possible after a problem was created\n");
      return SCIP_READERROR;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of bounds file is only possible during problem creation stage\n");
      return SCIP_READERROR;
   }

   /* read bounds file */
   SCIP_CALL( readBounds(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * bnd file reader specific interface methods
 */

/** includes the bnd file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderBnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyBnd) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadBnd) );

   return SCIP_OKAY;
}
