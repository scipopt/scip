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

/**@file   reader_cnf.c
 * @brief  CNF file reader
 * @author Thorsten Koch
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_cnf.h"
#include "cons_logicor.h"


#define READER_NAME             "CNFreader"
#define READER_DESC             "CNF file reader"
#define READER_EXTENSION        "cnf"



/*
 * CNF reader internal methods
 */

static
void readError(
   int              linecount,          /**< line number of error */
   const char*      errormsg            /**< error message */
   )
{
   char s[MAXSTRLEN];

   sprintf(s, "read error in line <%d>: %s", linecount, errormsg);
   errorMessage(s);
}

static
void readWarning(
   int              linecount,          /**< line number of error */
   const char*      warningmsg          /**< warning message */
   )
{
   char s[MAXSTRLEN];

   sprintf(s, "warning in line <%d>: %s", linecount, warningmsg);
   warningMessage(s);
}

/** reads the next non-empty non-comment line of a CNF file */
static
RETCODE readCNFLine(
   FILE*            file,               /**< input file */
   char*            buffer,             /**< buffer for storing the input line */
   int              size,               /**< size of the buffer */
   int*             linecount           /**< pointer to the line number counter */
   )
{
   char* line;
   int linelen;

   assert(file != NULL);
   assert(buffer != NULL);
   assert(size >= 2);
   assert(linecount != NULL);

   do
   {
      (*linecount)++;
      line = fgets(buffer, size, file);
      if( line != NULL )
      {
         linelen = strlen(line);
         if( linelen == size-1 )
         {
            char s[MAXSTRLEN];
            sprintf(s, "line too long (exceeds %d characters)", size-2);
            readError(*linecount, s);
            return SCIP_PARSEERROR;
         }
      }
      else
         linelen = 0;
   }
   while( line != NULL && (*line == 'c' || *line == '\n') );

   if( linelen >= 2 && line[linelen-2] == '\n' )
      line[linelen-2] = '\0';
   else if( linelen == 0 )
      *buffer = '\0';

   assert((line != NULL) ^ (*buffer == '\0'));
 
   return SCIP_OKAY;
}

/* Read SAT formula in "CNF File Format".
 * 
 *  The specification is taken from the
 *
 *  Satisfiability Suggested Format
 *
 *  Online available at http://www.intellektik.informatik.tu-darmstadt.de/SATLIB/Benchmarks/SAT/satformat.ps
 *
 *  The method reads all files of CNF format. Other formats (SAT, SATX, SATE) are not supported.
 */  
static
RETCODE readCNF(
   SCIP*            scip,               /**< SCIP data structure */   
   FILE*            file                /**< input file */
   )
{
   RETCODE retcode;
   VAR** vars;
   VAR** clausevars;
   CONS* cons;
   int* varsign;
   char* tok;
   char* nexttok;
   char line[MAXSTRLEN];
   char format[MAXSTRLEN];
   char varname[MAXSTRLEN];
   char s[MAXSTRLEN];
   Bool dynamiccols;
   Bool dynamicrows;
   Bool negative;
   int linecount;
   int clauselen;
   int clausenum;
   int nvars;
   int nclauses;
   int varnum;
   int v;

   assert(scip != NULL);
   assert(file != NULL);

   retcode = SCIP_OKAY;

   linecount = 0;

   /* read header */
   CHECK_OKAY( readCNFLine(file, line, sizeof(line), &linecount) );
   if( *line != 'p' )
   {
      readError(linecount, "problem declaration line expected");
      return SCIP_PARSEERROR;
   }
   if( sscanf(line, "p %8s %d %d", format, &nvars, &nclauses) != 3 )
   {
      readError(linecount, "invalid problem declaration (must be 'p cnf <nvars> <nclauses>')");
      return SCIP_PARSEERROR;
   }
   if( strcmp(format, "cnf") != 0 )
   {
      sprintf(s, "invalid format tag <%s> (must be 'cnf')", format);
      readError(linecount, s);
      return SCIP_PARSEERROR;
   }
   if( nvars <= 0 )
   {
      sprintf(s, "invalid number of variables <%d> (must be positive)", nvars);
      readError(linecount, s);
      return SCIP_PARSEERROR;
   }
   if( nclauses <= 0 )
   {
      sprintf(s, "invalid number of clauses <%d> (must be positive)", nclauses);
      readError(linecount, s);
      return SCIP_PARSEERROR;
   }

   /* get parameter values */
   CHECK_OKAY( SCIPgetBoolParam(scip, "reader/cnf/dynamiccols", &dynamiccols) );
   CHECK_OKAY( SCIPgetBoolParam(scip, "reader/cnf/dynamicrows", &dynamicrows) );

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &vars, nvars) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &clausevars, nvars) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &varsign, nvars) );

   /* create the variables */
   for( v = 0; v < nvars; ++v )
   {
      sprintf(varname, "x%d", v+1);
      CHECK_OKAY( SCIPcreateVar(scip, &vars[v], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, dynamiccols) );
      CHECK_OKAY( SCIPaddVar(scip, vars[v]) );
      varsign[v] = 0;
   }

   /* read clauses */
   clausenum = 0;
   clauselen = 0;
   do
   {
      retcode = readCNFLine(file, line, sizeof(line), &linecount);
      if( retcode != SCIP_OKAY )
         goto TERMINATE;

      if( *line != '\0' && *line != '%' )
      {
         tok = strtok_r(line, " \f\n\r\t", &nexttok);
         while( tok != NULL )
         {
            /* parse literal and check for errors */
            if( sscanf(tok, "%d", &v) != 1 )
            {
               sprintf(s, "invalid literal <%s>", tok);
               readError(linecount, s);
               retcode = SCIP_PARSEERROR;
               goto TERMINATE;
            }

            /* interpret literal number: v == 0: end of clause, v < 0: negated literal, v > 0: positive literal */
            if( v == 0 )
            {
               /* end of clause: construct clause and add it to SCIP */
               if( clauselen == 0 )
                  readWarning(linecount, "empty clause detected in line -- problem infeasible");

               clausenum++;
               sprintf(s, "c%d", clausenum);
               CHECK_OKAY( SCIPcreateConsLogicOr(scip, &cons, s, clauselen, clausevars, 
                              !dynamicrows, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicrows) );
               CHECK_OKAY( SCIPaddCons(scip, cons) );
               CHECK_OKAY( SCIPreleaseCons(scip, &cons) );
               clauselen = 0;
            }
            else if( v >= -nvars && v <= nvars )
            {
               if( clauselen >= nvars )
               {
                  readError(linecount, "too many literals in clause");
                  retcode = SCIP_PARSEERROR;
                  goto TERMINATE;
               }
         
               /* add literal to clause */
               varnum = ABS(v)-1;
               if( v < 0 )
               {
                  CHECK_OKAY( SCIPgetNegatedVar(scip, vars[varnum], &clausevars[clauselen]) );
                  varsign[varnum]--;
               }
               else
               {
                  clausevars[clauselen] = vars[varnum];
                  varsign[varnum]++;
               }
               clauselen++;
            }
            else
            {
               sprintf(s, "invalid variable number <%d>", ABS(v));
               readError(linecount, s);
               retcode = SCIP_PARSEERROR;
               goto TERMINATE;
            }

            /* get next token */
            tok = strtok_r(NULL, " \f\n\r\t", &nexttok);
         }
      }
   }
   while( *line != '\0' && *line != '%' );

   /* check for additional literals */
   if( clauselen > 0 )
   {
      sprintf(s, "found %d additional literals after last clause", clauselen);
      warningMessage(s);
   }

   /* check number of clauses */
   if( clausenum != nclauses )
   {
      sprintf(s, "expected %d clauses, but found %d", nclauses, clausenum);
      warningMessage(s);
   }

 TERMINATE:
   /* change objective values and release variables */
   CHECK_OKAY( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   for( v = 0; v < nvars; ++v )
   {
      CHECK_OKAY( SCIPchgVarObj(scip, vars[v], varsign[v]) );
      CHECK_OKAY( SCIPreleaseVar(scip, &vars[v]) );
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &varsign) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &clausevars) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &vars) );

   return retcode;
}



/*
 * Callback methods
 */

static
DECL_READERREAD(SCIPreaderReadCNF)
{
   FILE* f;
   RETCODE retcode;
   char s[MAXSTRLEN];
   int linecount;
   int nvars;
   int nclauses;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(filename != NULL);
   assert(result != NULL);

   /* open file */
   f = fopen(filename, "r");
   if( f == NULL )
   {
      char s[1024];
      sprintf(s, "cannot open file <%s> for reading", filename);
      errorMessage(s);
      perror(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   CHECK_OKAY( SCIPcreateProb(scip, filename, NULL, NULL, NULL) );

   /* read CNF file */
   retcode = readCNF(scip, f);

   /* close file */
   fclose(f);

   *result = SCIP_SUCCESS;

   return retcode;
}




/*
 * CNF file reader specific interface methods
 */

/** includes the CNF file reader in SCIP */
RETCODE SCIPincludeReaderCNF(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* include CNF reader */
   CHECK_OKAY( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
                  NULL, SCIPreaderReadCNF, NULL) );

   /* add CNF reader parameters */
   CHECK_OKAY( SCIPaddBoolParam(scip,
                  "reader/cnf/dynamiccols", "should columns be added and removed dynamically to the LP?",
                  NULL, FALSE, NULL, NULL) );
   CHECK_OKAY( SCIPaddBoolParam(scip,
                  "reader/cnf/dynamicrows", "should rows be added and removed dynamically to the LP?",
                  NULL, FALSE, NULL, NULL) );
   
   return SCIP_OKAY;
}

