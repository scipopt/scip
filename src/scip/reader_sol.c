/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_sol.c,v 1.22 2011/01/02 11:10:43 bzfheinz Exp $"

/**@file   reader_sol.c
 * @ingroup FILEREADERS 
 * @brief  file reader for primal solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
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

/** reads a given xml solution file */
static
SCIP_RETCODE readXMLSol(
   SCIP*                 scip,              /**< SCIP data structure */
   const char*           filename           /**< name of the input file */
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
   start = xml_process(filename);

   if ( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing the XML solution file.\n");
      return SCIP_READERROR;
   }

   /* create primal solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   error = FALSE;

   /* find variable sections */
   tag = "variables";
   varsnode = xml_find_node(start, tag);
   if ( varsnode == NULL )
   {
      /* free xml data */
      xml_free_node(start);

      SCIPerrorMessage("Variable section not found.\n");
      return SCIP_READERROR;
   }

   /* loop through all variables */
   unknownvariablemessage = FALSE;
   for (varnode = xml_first_child(varsnode); varnode != NULL; varnode = xml_next_sibl(varnode))
   {
      const char* varname;
      const char* varvalue;
      SCIP_VAR* var;
      SCIP_Real value;
      int nread;

      /* find variable name */
      varname = xml_get_attrval(varnode, "name");
      if ( varname == NULL )
      {
         SCIPerrorMessage("Attribute \"name\" of variable not found.\n");
         error = TRUE;
         break;
      }
      
      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if ( var == NULL )
      {
         if ( !unknownvariablemessage )
         {
            SCIPwarningMessage("unknown variable <%s> in XML solution file <%s>\n", varname, filename);
            SCIPwarningMessage("  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* find value of variable */
      varvalue = xml_get_attrval(varnode, "value");
      if ( varvalue == NULL )
      {
         SCIPerrorMessage("Attribute \"value\" of variable not found.\n");
         error = TRUE;
         break;
      }

      /* cast the value */
      if ( strncasecmp(varvalue, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(varvalue, "+inf", 4) == 0 || strncasecmp(varvalue, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(varvalue, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(varvalue, "%lf", &value);
         if ( nread != 1 )
         {
            SCIPwarningMessage("invalid solution value <%s> for variable <%s> in XML solution file <%s>\n", varvalue, varname, filename);
            error = TRUE;
            break;
         }
      }

      /* set the solution value of the variable */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var, value) );
   }

   if ( ! error )
   {
      SCIP_Bool stored;

      /* add and free the solution */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, &stored) );

      /* display result */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "primal solution from solution file <%s> was %s\n",
         filename, stored ? "accepted" : "rejected - solution is infeasible or objective too poor");
   }
   else
   {
      /* free solution */
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      /* free xml data */
      xml_free_node(start);

      return SCIP_READERROR;
   }

   /* free xml data */
   xml_free_node(start);

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


/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSol NULL


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

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPwarningMessage("reading of solution file is only possible after a problem was created\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
         "primal solution from solution file <%s> was ignored - problem is already solved to optimality\n",
         filename);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* transform the problem such that adding primal solutions is possible */
   SCIP_CALL( SCIPtransformProb(scip) );

   /* open input file in order to determine type */
   file = SCIPfopen(filename, "r");
   if ( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* get next line */
   if ( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
   {
      SCIPerrorMessage("cannot parse file.\n");
      return SCIP_READERROR;
   }
   /* close file */
   SCIPfclose(file);

   /* decide whether it is xml */
   s = buffer;
   
   /* skip spaces */
   while ( isspace(*s) )
      ++s;
   if ( s[0] == '<' && s[1] == '?' && s[2] == 'x' && s[3] == 'm' && s[4] == 'l' )
   {
      /* read XML solution and add it to the solution pool */
      SCIP_CALL( readXMLSol(scip, filename) );
   }
   else
   {
      /* read the solution and add it to the solution pool */
      SCIP_CALL( SCIPreadSol(scip, filename) );
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteSol NULL


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
         readerCopySol, readerFreeSol, readerReadSol, readerWriteSol, 
         readerdata) );

   return SCIP_OKAY;
}

