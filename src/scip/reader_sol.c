/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_sol.c,v 1.3 2006/04/10 16:15:27 bzfpfend Exp $"

/**@file   reader_sol.c
 * @brief  file reader for partial solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <strings.h>

#include "scip/reader_sol.h"


#define READER_NAME             "solreader"
#define READER_DESC             "file reader for (partial) solutions"
#define READER_EXTENSION        "sol"



/*
 * local methods
 */

/** reads the given solution file */
static
SCIP_RETCODE readSol(
   SCIP*                 scip,               /**< SCIP data structure */   
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Bool unknownvariablemessage;
   int lineno;
   int nfixed;

   assert(scip != NULL);
   assert(filename != NULL);

   /* open input file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }   

   /* read the file */
   error = FALSE;
   unknownvariablemessage = FALSE;
   lineno = 0;
   nfixed = 0;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* the lines "solution status: ..." and "objective value: ..." may preceed the solution information */
      if( strncasecmp(buffer, "solution", 8) == 0 || strncasecmp(buffer, "objective", 9) == 0 )
         continue;

      /* parse the line */
      nread = sscanf(buffer, "%s %s %s\n", varname, valuestring, objstring);
      if( nread < 2 )
      {
         SCIPwarningMessage("invalid input line %d in solution file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPwarningMessage("unknown variable <%s> in line %d of solution file <%s>\n", varname, lineno, filename);
            SCIPwarningMessage("  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the value */
      if( strncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(valuestring, "+inf", 4) == 0 || strncasecmp(valuestring, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(valuestring, "%lf", &value);
         if( nread != 1 )
         {
            SCIPwarningMessage("invalid solution value <%s> for variable <%s> in line %d of solution file <%s>\n",
               valuestring, varname, lineno, filename);
            error = TRUE;
            break;
         }
      }

      /* fix the variable */
      SCIP_CALL( SCIPfixVar(scip, var, value, &infeasible, &fixed) );
      if( infeasible )
      {
         SCIPwarningMessage("infeasible solution value of <%s>[%g,%g] to %g in line %d of solution file <%s>\n",
            varname, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), value, lineno, filename);
         error = TRUE;
         break;
      }
      if( fixed )
         nfixed++;
   }

   /* close input file */
   SCIPfclose(file);

   /* display result */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "fixed %d variables from solution file <%s>\n", nfixed, filename);

   if( error )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSol NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSol)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPwarningMessage("reading of solution file is only possible after a problem was created\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( readSol(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * sol file reader specific interface methods
 */

/** includes the sol file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create sol reader data */
   readerdata = NULL;

   /* include sol reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeSol, readerReadSol, readerdata) );

   return SCIP_OKAY;
}

