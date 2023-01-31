/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_lop.c
 * @brief  linear ordering file reader
 * @author Marc Pfetsch
 *
 * This file implements the reader/parser used to read linear ordering problems. For more details see \ref READER. The
 * data should be given in LOLIB format, see <a href="http://grafo.etsii.urjc.es/optsicom/lolib/">LOLIB</a>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "cons_lop.h"
#include "reader_lop.h"
#include <scip/pub_fileio.h>

#define READER_NAME             "lopreader"
#define READER_DESC             "file reader for linear ordering problems"
#define READER_EXTENSION        "lop"

#define READER_STRLEN           65536

/* ----------------- auxiliary functions ------------------------ */

/** Get next number from file */
static
SCIP_RETCODE getNextNumber(
   SCIP_FILE*            file,               /**< file to read from */
   char*                 buffer,             /**< buffer to store lines in */
   char**                s,                  /**< pointer to string pointer */
   SCIP_Real*            value               /**< pointer to store the read value */
   )
{
   assert( file != NULL );
   assert( buffer != NULL );
   assert( s != NULL );
   assert( value != NULL );

   *value = SCIP_INVALID; /* for debugging */

   /* skip whitespace */
   while ( isspace(**s) )
      ++(*s);

   /* if we reached the end of the line, read new line */
   if ( **s == '\n' || **s == '\0' )
   {
      /* read line into buffer */
      if ( SCIPfgets(buffer, READER_STRLEN, file) == NULL )
      {
         SCIPerrorMessage("Error reading from file.\n");
         return SCIP_READERROR;
      }
      *s = buffer;

      /* skip whitespace */
      while ( isspace(**s) )
         ++(*s);
   }

   /* check whether we found a number */
   if ( isdigit(**s) || **s == '-' || **s == '+' )
   {
      *value = atof(*s);

      /* skip number */
      while ( isdigit(**s) || **s == '-' || **s == '+' )
         ++(*s);
   }
   else
   {
      SCIPerrorMessage("Error reading from file.\n");
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}



/** read weight matrix from file (in LOLIB format)
 *
 *  Format:
 *  comment line
 *  n = # of elements
 *  weight matrix (n times n double matrix)
 */
static
SCIP_RETCODE LOPreadFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   int*                  n,                  /**< pointer to store the number of elements on exit */
   SCIP_Real***          W                   /**< pointer to store the weight matrix on exit */
   )
{
   char buffer[READER_STRLEN];
   SCIP_FILE* file;
   char* s;
   char* nstr;
   int i;
   int j;

   assert( n != NULL );
   assert( W != NULL );

   /* open file */
   file = SCIPfopen(filename, "r");
   if ( file == NULL )
   {
      SCIPerrorMessage("Could not open file <%s>.\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* skip lines as comments until we found the first line that only contains the number of elements */
   *n = -1;
   do
   {
      /* read line */
      if ( SCIPfgets(buffer, READER_STRLEN, file) == NULL )
      {
         SCIPerrorMessage("Error reading file <%s>.\n", filename);
         return SCIP_READERROR;
      }

      /* skip whitespace */
      s = buffer;
      while( isspace(*s) )
         ++s;

      /* check whether rest of line only contains whitespace or numbers */
      nstr = s;
      while ( *s != '\0' && (isspace(*s) || isdigit(*s)) )
         ++s;

      /* if the line only contains a number, use this as the number of elements */
      if ( *s == '\0' || *s == '\n' )
         *n = atoi(nstr);
   }
   while ( ! SCIPfeof(file) && *n < 0 );

   if ( *n <= 0 )
   {
      SCIPerrorMessage("Reading of number of elements failed.\n");
      return SCIP_READERROR;
   }
   assert( *n > 0 );

   /* set up matrix */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, W, *n) );
   for (i = 0; i < *n; ++i)
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*W)[i]), *n) ); /*lint !e866*/

   /* read weight matrix */
   for (i = 0; i < *n; ++i)
   {
      for (j = 0; j < *n; ++j)
      {
	 SCIP_Real val;

         SCIP_CALL( getNextNumber(file, buffer, &s, &val) );
         assert( val != SCIP_INVALID );
	 (*W)[i][j] = val;
      }
   }
   (void) SCIPfclose(file);

   return SCIP_OKAY;
}

/** get problem name
 *
 *  Returns NULL on error
 */
static
SCIP_RETCODE getProblemName(
   const char*           filename,           /**< input filename */
   char*                 probname,           /**< output problemname */
   int                   maxSize             /**< maximum size of probname */
   )
{
   int i = 0;
   int j = 0;
   int l;

   /* first find end of string */
   while ( filename[i] != 0)
      ++i;
   l = i;

   /* go back until '.' or '/' or '\' appears */
   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/') && (filename[i] != '\\'))
      --i;

   /* if we found '.', search for '/' or '\\' */
   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/') && (filename[i] != '\\'))
	 --i;
   }

   /* correct counter */
   if ((filename[i] == '/') || (filename[i] == '\\'))
      ++i;

   /* copy name */
   while ( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if (j > maxSize-1)
	 return SCIP_ERROR;
   }
   probname[j] = 0;

   return SCIP_OKAY;
}



/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(LOPreaderRead)
{  /*lint --e{715}*/
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR*** vars;
   SCIP_CONS* cons;
   SCIP_Real** W;
   int n;
   int i;
   int j;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPinfoMessage(scip, NULL, "File name:\t\t%s\n", filename);

   /* read problem */
   SCIP_CALL( LOPreadFile(scip, filename, &n, &W) );

   /* generate problem name from filename */
   SCIP_CALL( getProblemName(filename, name, SCIP_MAXSTRLEN) );

   SCIPinfoMessage(scip, NULL, "Problem name:\t\t%s\n", name);
   SCIPinfoMessage(scip, NULL, "Number of elements:\t%d\n\n", n);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, name) );

   /* set maximization */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* generate variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, n) );
   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(vars[i]), n) );
      for (j = 0; j < n; ++j)
      {
	 if (j != i)
	 {
	    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x#%d#%d", i, j);
	    SCIP_CALL( SCIPcreateVar(scip, &(vars[i][j]), name, 0.0, 1.0, W[i][j], SCIP_VARTYPE_BINARY,
		  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
	    SCIP_CALL( SCIPaddVar(scip, vars[i][j]) );
	 }
	 else
	    vars[i][j] = NULL;
      }
   }

   /* generate linear ordering constraint */
   SCIP_CALL( SCIPcreateConsLOP(scip, &cons, "LOP", n, vars, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
	 FALSE, FALSE, FALSE, FALSE));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* print small model */
   if ( n <= 10 )
   {
      SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
      SCIPinfoMessage(scip, NULL, "\n");
   }

   /* free memory */
   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
      {
	 if (j != i)
	 {
	    SCIP_CALL( SCIPreleaseVar(scip, &(vars[i][j])) );
	 }
      }
      SCIPfreeBlockMemoryArray(scip, &(vars[i]), n);
      SCIPfreeBlockMemoryArray(scip, &(W[i]), n);
   }
   SCIPfreeBlockMemoryArray(scip, &vars, n);
   SCIPfreeBlockMemoryArray(scip, &W, n);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the linear ordering file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderLOP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );
   assert( reader != NULL );

   SCIP_CALL( SCIPsetReaderRead(scip, reader, LOPreaderRead) );

   return SCIP_OKAY;
}

/**@} */
