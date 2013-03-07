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

/**@file   reader_sol.c
 * @brief  file reader for primal solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/reader_sol.h"
#include "xml/xml.h"

#define READER_NAME             "solreader"
#define READER_DESC             "file reader for primal solutions"
#define READER_EXTENSION        "sol"


/*
 * Local methods of reader
 */

/** reads a given SCIP solution file, problem has to be transformed in advance */
static
SCIP_RETCODE readSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname               /**< name of the input file */
   )
{
   SCIP_SOL* sol;
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Bool unknownvariablemessage;
   SCIP_Bool stored;
   SCIP_Bool usevartable;
   int lineno;

   assert(scip != NULL);
   assert(fname != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );

   if( !usevartable )
   {
      SCIPerrorMessage("Cannot read solution file if vartable is disabled. Make sure parameter 'misc/usevartable' is set to TRUE.\n");
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

   /* create zero solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   
   /* read the file */
   error = FALSE;
   unknownvariablemessage = FALSE;
   lineno = 0;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real value;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* there are some lines which may preceed the solution information */
      if( strncasecmp(buffer, "solution status:", 16) == 0 || strncasecmp(buffer, "objective value:", 16) == 0 ||
         strncasecmp(buffer, "Log started", 11) == 0 || strncasecmp(buffer, "Variable Name", 13) == 0 ||
         strncasecmp(buffer, "All other variables", 19) == 0 || strncasecmp(buffer, "\n", 1) == 0 || 
         strncasecmp(buffer, "NAME", 4) == 0 || strncasecmp(buffer, "ENDATA", 6) == 0 )    /* allow parsing of SOL-format on the MIPLIB 2003 pages */
         continue;

      /* parse the line */
      nread = sscanf(buffer, "%s %s %s\n", varname, valuestring, objstring);
      if( nread < 2 )
      {
         SCIPerrorMessage("Invalid input line %d in solution file <%s>: <%s>.\n", lineno, fname, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> in line %d of solution file <%s>\n", 
               varname, lineno, fname);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
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
            SCIPerrorMessage("Invalid solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
               valuestring, varname, lineno, fname);
            error = TRUE;
            break;
         }
      }

      /* set the solution value of the variable, if not multiaggregated */
      if( SCIPisTransformed(scip) && SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n", SCIPvarGetName(var));
      }
      else
      {
         SCIP_RETCODE retcode;
         retcode =  SCIPsetSolVal(scip, sol, var, value);

         if( retcode == SCIP_INVALIDDATA )
         {
            if( SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_FIXED )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored conflicting solution value for fixed variable <%s>\n",
                  SCIPvarGetName(var));
            }
            else
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n",
                  SCIPvarGetName(var));
            }
         }
         else
         {
            SCIP_CALL( retcode );
         }
      }
   }

   /* close input file */
   SCIPfclose(file);

   if( !error )
   {
      /* add and free the solution */
      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, &stored) );

         /* display result */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "primal solution from solution file <%s> was %s\n",
            fname, stored ? "accepted" : "rejected - solution is infeasible or objective too poor");
      }
      else
      {
         /* add primal solution to solution candidate storage, frees the solution afterwards */
         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );

         /* display result */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "primal solution from solution file <%s> was %s\n",
            fname, stored ? "accepted as candidate, will be checked when solving starts" : "rejected - solution objective too poor");
      }

      return SCIP_OKAY;
   }
   else
   {
      /* free solution */
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      return SCIP_READERROR;
   }
}

/** reads a given xml solution file */
static
SCIP_RETCODE readXMLSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_Bool unknownvariablemessage;
   SCIP_SOL* sol;
   SCIP_Bool error;
   XML_NODE* start;
   const XML_NODE* varsnode;
   const XML_NODE* varnode;
   const char* tag;

   assert( scip != NULL );
   assert( filename != NULL );

   /* read xml file */
   start = xmlProcess(filename);

   if( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing the XML solution file.\n");
      return SCIP_READERROR;
   }
   
   /* create zero solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   
   error = FALSE;
   
   /* find variable sections */
   tag = "variables";
   varsnode = xmlFindNodeMaxdepth(start, tag, 0, 3);
   if( varsnode == NULL )
   {
      /* free xml data */
      xmlFreeNode(start);

      SCIPerrorMessage("Variable section not found.\n");
      return SCIP_READERROR;
   }

   /* loop through all variables */
   unknownvariablemessage = FALSE;
   for( varnode = xmlFirstChild(varsnode); varnode != NULL; varnode = xmlNextSibl(varnode) )
   {
      const char* varname;
      const char* varvalue;
      SCIP_VAR* var;
      SCIP_Real value;
      int nread;

      /* find variable name */
      varname = xmlGetAttrval(varnode, "name");
      if( varname == NULL )
      {
         SCIPerrorMessage("Attribute \"name\" of variable not found.\n");
         error = TRUE;
         break;
      }
      
      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> of solution file <%s>\n", 
               varname, filename);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* find value of variable */
      varvalue = xmlGetAttrval(varnode, "value");
      if( varvalue == NULL )
      {
         SCIPerrorMessage("Attribute \"value\" of variable not found.\n");
         error = TRUE;
         break;
      }

      /* cast the value */
      if( strncasecmp(varvalue, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(varvalue, "+inf", 4) == 0 || strncasecmp(varvalue, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(varvalue, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(varvalue, "%lf", &value);
         if( nread != 1 )
         {
            SCIPwarningMessage(scip, "invalid solution value <%s> for variable <%s> in XML solution file <%s>\n", varvalue, varname, filename);
            error = TRUE;
            break;
         }
      }

      /* set the solution value of the variable */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var, value) );
   }

   if( !error )
   {
      SCIP_Bool stored;

      /* add and free the solution */
      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, &stored) );

         /* display result */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "primal solution from solution file <%s> was %s\n",
            filename, stored ? "accepted" : "rejected - solution is infeasible or objective too poor");
      }
      else
      {
         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );

         /* display result */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "primal solution from solution file <%s> was %s\n",
            filename, stored ? "accepted as candidate, will be checked when solving starts" : "rejected - solution objective too poor");
      }
   }
   else
   {
      /* free solution */
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      /* free xml data */
      xmlFreeNode(start);

      return SCIP_READERROR;
   }

   /* free xml data */
   xmlFreeNode(start);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderSol(scip) );
 
   return SCIP_OKAY;
}


/** problem reading method of reader 
 *
 *  In order to determine the type of the file, we have to open it. Thus, it has to be opened
 *  twice. This might be removed, but is likely to not hurt the performance too much.
 */
static
SCIP_DECL_READERREAD(readerReadSol)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   char buffer[SCIP_MAXSTRLEN];
   char *s;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of solution file is only possible after a problem was created\n");
      return SCIP_READERROR;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
         "primal solution from solution file <%s> was ignored - problem is already solved to optimality\n",
         filename);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* open input file in order to determine type */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* get next line */
   if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
   {
      SCIPerrorMessage("cannot parse file.\n");
      return SCIP_READERROR;
   }
   /* close file */
   SCIPfclose(file);

   /* decide whether it is xml */
   s = buffer;
   
   /* skip spaces */
   while( isspace((unsigned char)*s) )
      ++s;
   if( s[0] == '<' && s[1] == '?' && s[2] == 'x' && s[3] == 'm' && s[4] == 'l' )
   {
      /* read XML solution and add it to the solution pool */
      SCIP_CALL( readXMLSol(scip, filename) );
   }
   else
   {
      /* read the solution and add it to the solution pool */
      SCIP_CALL( readSol(scip, filename) );
   }

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
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySol) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSol) );

   return SCIP_OKAY;
}

