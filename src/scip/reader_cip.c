/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   reader_cip.c
 * @ingroup DEFPLUGINS_READER
 * @brief  CIP file reader
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/rational.h"
#include "scip/cons_linear.h"
#include "scip/cons_exactlinear.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/pub_var.h"
#include "scip/reader_cip.h"
#include "scip/scip_exact.h"
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include "scip/scip_var.h"


#define READER_NAME             "cipreader"
#define READER_DESC             "file reader for CIP (Constraint Integer Program) format"
#define READER_EXTENSION        "cip"

#define DEFAULT_CIP_WRITEFIXEDVARS  TRUE     /**< Should fixed and aggregated variables be written when writing? */


/** CIP reading data */
struct SCIP_ReaderData
{
   SCIP_Bool             writefixedvars;     /**< Should fixed and aggregated variables be written when writing? */
};


/** Section of the in CIP files */
enum CipSection
{
   CIP_START,            /**< start tag */
   CIP_STATISTIC,        /**< statistics section */
   CIP_OBJECTIVE,        /**< objective */
   CIP_VARS,             /**< list of (free) variables */
   CIP_FIXEDVARS,        /**< list of fixed variables */
   CIP_CONSTRAINTS,      /**< constraints */
   CIP_END               /**< end of file tag */
};
typedef enum CipSection CIPSECTION;          /**< Section of the in CIP files */


/*
 * Data structures
 */

/** CIP reading data */
struct CipInput
{
   SCIP_FILE*            file;               /**< input file */
   char*                 strbuf;             /**< string buffer for input lines */
   int                   len;                /**< length of strbuf */
   int                   readingsize;        /**< size of block in which len is increased if necessary */
   int                   linenumber;         /**< number of line in input file */
   CIPSECTION            section;            /**< current section */
   SCIP_Bool             haserror;           /**< some error occurred */
   SCIP_Bool             endfile;            /**< we have reached the end of the file */
   SCIP_Real             objoffset;          /**< real objective offset */
   SCIP_Real             objscale;           /**< real objective scale */
   SCIP_RATIONAL*        objoffsetexact;     /**< exact objective offset */
   SCIP_RATIONAL*        objscaleexact;      /**< exact objective scale */
};
typedef struct CipInput CIPINPUT;            /**< CIP reading data */


/*
 * Local methods for reading/parsing
 */

/** get next input line; this are all characters until the next semicolon */
static
SCIP_RETCODE getInputString(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput            /**< CIP parsing data */
   )
{
   char* endline;
   char* endcharacter;
   char* windowsendlinechar;

   assert(cipinput != NULL);

   /* read next line */
   cipinput->endfile = (SCIPfgets(cipinput->strbuf, cipinput->len, cipinput->file) == NULL);

   if( cipinput->endfile )
   {
      /* clear the line for safety reason */
      BMSclearMemoryArray(cipinput->strbuf, cipinput->len);
      return SCIP_OKAY;
   }

   cipinput->linenumber++;
   endline = strchr(cipinput->strbuf, '\n');
   endcharacter = strchr(cipinput->strbuf, ';');

   while( endline == NULL || (endcharacter == NULL && cipinput->section == CIP_CONSTRAINTS && strncmp(cipinput->strbuf, "END", 3) != 0 ) )
   {
      int pos;

      /* we refill the buffer from the '\n' character */
      if( endline == NULL )
         pos = cipinput->len - 1;
      else
         pos = (int) (endline - cipinput->strbuf);

      /* don't erase the '\n' from all buffers for constraints */
      if( endline != NULL && cipinput->section == CIP_CONSTRAINTS )
         pos++;

      /* if necessary reallocate memory */
      if( pos + cipinput->readingsize >= cipinput->len )
      {
         cipinput->len = SCIPcalcMemGrowSize(scip, pos + cipinput->readingsize);
         SCIP_CALL( SCIPreallocBufferArray(scip, &(cipinput->strbuf), cipinput->len) );
      }

      /* read next line */
      cipinput->endfile = (SCIPfgets(&(cipinput->strbuf[pos]), cipinput->len - pos, cipinput->file) == NULL);

      if( cipinput->endfile )
      {
         /* clear the line for safety reason */
         BMSclearMemoryArray(cipinput->strbuf, cipinput->len);
         return SCIP_OKAY;
      }

      cipinput->linenumber++;
      endline = strrchr(&cipinput->strbuf[pos], '\n');
      endcharacter = strchr(&cipinput->strbuf[pos], ';');
   }
   assert(endline != NULL);

   /*SCIPdebugMsg(scip, "read line: %s\n", cipinput->strbuf);*/

   /* check for windows "carriage return" endline character */
   windowsendlinechar = strrchr(cipinput->strbuf, '\r');
   if( windowsendlinechar != NULL && windowsendlinechar + 1 == endline )
      --endline;
   else
      /* if the assert should not hold we found a windows "carriage return" which was not at the end of the line */
      assert(windowsendlinechar == NULL);

   if( cipinput->section == CIP_CONSTRAINTS && endcharacter != NULL && endline - endcharacter != 1 )
   {
      SCIPerrorMessage("Constraint line has to end with ';\\n' (line: %d).\n", cipinput->linenumber);
      cipinput->haserror = TRUE;
      return SCIP_OKAY; /* return error at hightest level */
   }

   *endline = '\0';

   return SCIP_OKAY;
}

/** read the problem name out of the statistics */
static
void getStart(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput            /**< CIP parsing data */
   )
{
   char* buf;

   assert(scip != NULL);

   buf = cipinput->strbuf;

   if( strncmp(buf, "STATISTICS", 9) == 0 )
   {
      cipinput->section = CIP_STATISTIC;
      return;
   }

   if( strncmp(buf, "VARIABLES", 8) == 0 || strncmp(buf, "FIXED", 5) == 0 || strncmp(buf, "CONSTRAINTS", 11) == 0 || strncmp(buf, "OBJECTIVE", 9) == 0 )
   {
      SCIPerrorMessage("Syntax Error: File has to start with 'STATISTICS' section.\n");
      cipinput->haserror = TRUE;
   }
}


/** read the problem name out of the statistics */
static
SCIP_RETCODE getStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput            /**< CIP parsing data */
   )
{
   char* buf;

   buf = cipinput->strbuf;

   if( strncmp(buf, "OBJECTIVE", 9) == 0 )
   {
      cipinput->section = CIP_OBJECTIVE;
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "parse statistics\n");

   if( strncmp(buf, "  Problem name", 14) == 0 )
   {
      char* name;
      char* s;

      name = strchr(buf, ':');

      if( name == NULL )
      {
         SCIPwarningMessage(scip, "did not find problem name (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
         return SCIP_OKAY;  /* no error, might work with empty problem name */
      }

      /* skip ':' */
      ++name;

      /* make sure that we terminate the string at comments ('#') or newline ('\r', '\n')*/
      if( NULL != (s = strpbrk(name, "#\r\n")) )
         *s = '\0';

      /* remove white space (tabs, ' ') in front of the name */
      SCIP_CALL( SCIPskipSpace(&name) );

      /* set problem name */
      SCIP_CALL( SCIPsetProbName(scip, name) );

      SCIPdebugMsg(scip, "problem name <%s>\n", name);
   }

   return SCIP_OKAY;
}

/** read objective sense, offset, and scale */
static
SCIP_RETCODE getObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput            /**< CIP parsing data */
   )
{
   SCIP_Bool success;
   char* buf;
   char* name;

   buf = cipinput->strbuf;

   if( strncmp(buf, "VARIABLES", 8) == 0 )
      cipinput->section = CIP_VARS;
   else if( strncmp(buf, "FIXED", 5) == 0 )
      cipinput->section = CIP_FIXEDVARS;
   else if( strncmp(buf, "CONSTRAINTS", 11) == 0 )
      cipinput->section = CIP_CONSTRAINTS;
   else if( strncmp(buf, "END", 3) == 0 )
      cipinput->section = CIP_END;

   if( cipinput->section != CIP_OBJECTIVE )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "parse objective information\n");

   /* remove white space */
   SCIP_CALL( SCIPskipSpace(&buf) );

   if( SCIPstrncasecmp(buf, "Sense", 5) == 0 )
   {
      SCIP_OBJSENSE objsense;

      name = strchr(buf, ':');

      if( name == NULL )
      {
         SCIPwarningMessage(scip, "did not find objective sense (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
         return SCIP_OKAY; /* no error - might work with default */
      }

      /* skip ':' */
      ++name;

      /* remove white space in front of the name */
      SCIP_CALL( SCIPskipSpace(&name) );

      if( SCIPstrncasecmp(name, "min", 3) == 0 )
         objsense = SCIP_OBJSENSE_MINIMIZE;
      else if( SCIPstrncasecmp(name, "max", 3) == 0 )
         objsense = SCIP_OBJSENSE_MAXIMIZE;
      else
      {
         SCIPwarningMessage(scip, "unknown objective sense '%s' (line: %d):\n%s\n", name, cipinput->linenumber, cipinput->strbuf);
         return SCIP_OKAY; /* no error - might work with default */
      }

      /* set problem name */
      SCIP_CALL( SCIPsetObjsense(scip, objsense) );
      SCIPdebugMsg(scip, "objective sense <%s>\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");
   }
   else if( SCIPstrncasecmp(buf, "Offset", 6) == 0 )
   {
      char* endptr;

      name = strchr(buf, ':');

      if( name == NULL )
      {
         SCIPwarningMessage(scip, "did not find offset (line: %d)\n", cipinput->linenumber);
         return SCIP_OKAY;
      }

      /* skip ':' */
      ++name;

      /* read exact offset */
      if( cipinput->objoffsetexact != NULL )
      {
         success = SCIPparseRational(scip, name, cipinput->objoffsetexact, &endptr);

         if( success )
            SCIPrationalDebugMessage("read exact objoffset %q\n", cipinput->objoffsetexact);
         else
            SCIPrationalSetReal(cipinput->objoffsetexact, 0.0);
      }
      /* read real offset */
      else
      {
         success = SCIPparseReal(scip, name, &cipinput->objoffset, &endptr);

         if( success )
            SCIPdebugMsg(scip, "read real objoffset %g\n", cipinput->objoffset);
         else
            cipinput->objoffset = 0.0;
      }

      if( !success )
      {
         SCIPwarningMessage(scip, "could not parse offset (line: %d)\n%s\n", cipinput->linenumber, cipinput->strbuf);
         return SCIP_OKAY;
      }
   }
   else if( SCIPstrncasecmp(buf, "Scale", 5) == 0 )
   {
      char* endptr;

      name = strchr(buf, ':');

      if( name == NULL )
      {
         SCIPwarningMessage(scip, "did not find scale (line: %d)\n", cipinput->linenumber);
         return SCIP_OKAY;
      }

      /* skip ':' */
      ++name;

      /* read exact scale */
      if( cipinput->objscaleexact != NULL )
      {
         success = SCIPparseRational(scip, name, cipinput->objscaleexact, &endptr);

         if( success )
            SCIPrationalDebugMessage("read exact objscale %q\n", cipinput->objscaleexact);
         else
            SCIPrationalSetReal(cipinput->objscaleexact, 1.0);
      }
      /* read real scale */
      else
      {
         success = SCIPparseReal(scip, name, &cipinput->objscale, &endptr);

         if( success )
            SCIPdebugMsg(scip, "read real objscale %g\n", cipinput->objscale);
         else
            cipinput->objscale = 1.0;
      }

      if( !success )
      {
         SCIPwarningMessage(scip, "could not parse objective scale (line: %d)\n%s\n", cipinput->linenumber, cipinput->strbuf);
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** read variable */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput,           /**< CIP parsing data */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable           /**< is var's column removable from the LP (due to aging or cleanup)? */
   )
{
   SCIP_Bool success;
   SCIP_VAR* var;
   char* buf;
   char* endptr;

   buf = cipinput->strbuf;

   if( strncmp(buf, "FIXED", 5) == 0 )
      cipinput->section = CIP_FIXEDVARS;
   else if( strncmp(buf, "CONSTRAINTS", 4) == 0 )
      cipinput->section = CIP_CONSTRAINTS;
   else if( strncmp(buf, "END", 3) == 0 )
      cipinput->section = CIP_END;

   if( cipinput->section != CIP_VARS )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "parse variable\n");

   /* parse the variable */
   SCIP_CALL( SCIPparseVar(scip, &var, buf, initial, removable, NULL, NULL, NULL, NULL, NULL, &endptr, &success) );

   if( !success )
   {
      SCIPerrorMessage("syntax error in variable information (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }

   /* scale exact objective */
   if( cipinput->objscaleexact != NULL )
   {
      if( !SCIPrationalIsEQReal(cipinput->objscaleexact, 1.0) )
      {
         SCIP_RATIONAL* newobjval;

         SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &newobjval) );

         SCIPrationalMult(newobjval, SCIPvarGetObjExact(var), cipinput->objscaleexact);
         SCIP_CALL( SCIPchgVarObjExact(scip, var, newobjval) );

         SCIPrationalFreeBuffer(SCIPbuffer(scip), &newobjval);
      }
   }
   /* scale real objective */
   else
   {
      if( cipinput->objscale != 1.0 ) /*lint !e777*/
      {
         SCIP_CALL( SCIPchgVarObj(scip, var, SCIPvarGetObj(var) * cipinput->objscale) );
      }
   }

   SCIP_CALL( SCIPaddVar(scip, var) );

   SCIPdebug( SCIP_CALL( SCIPprintVar(scip, var, NULL) ) );

   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** read fixed variable */
static
SCIP_RETCODE getFixedVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput            /**< CIP parsing data */
   )
{
   SCIP_Bool success;
   SCIP_VAR* var;
   char* buf;
   char* endptr;
   char name[SCIP_MAXSTRLEN];

   buf = cipinput->strbuf;

   if( strncmp(buf, "CONSTRAINTS", 11) == 0 )
      cipinput->section = CIP_CONSTRAINTS;
   else if( strncmp(buf, "END", 3) == 0 )
      cipinput->section = CIP_END;

   if( cipinput->section != CIP_FIXEDVARS )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "parse fixed variable\n");

   /* parse the variable */
   SCIP_CALL( SCIPparseVar(scip, &var, buf, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL, &endptr, &success) );

   if( !success )
   {
      SCIPerrorMessage("syntax error in variable information (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }

   /* skip intermediate stuff */
   buf = endptr;

   while( *buf != '\0' && (*buf == ' ' || *buf == ',') )
      ++buf;

   /* check whether variable is fixed */
   if( strncmp(buf, "fixed:", 6) == 0 )
   {
      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIPdebug( SCIP_CALL( SCIPprintVar(scip, var, NULL) ) );
   }
   else if( strncmp(buf, "negated:", 8) == 0 )
   {
      SCIP_CONS* lincons = NULL;
      SCIP_VAR* negvar;
      SCIP_VAR* vars[2];

      buf += 8;

      /* we can just parse the next variable (ignoring all other information in between) */
      SCIP_CALL( SCIPparseVarName(scip, buf, &negvar, &endptr) );

      if( negvar == NULL )
      {
         SCIPerrorMessage("could not parse negated variable (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
         cipinput->haserror = TRUE;
         return SCIP_OKAY;
      }

      assert(SCIPvarIsBinary(var));
      assert(SCIPvarIsBinary(negvar));

      SCIP_CALL( SCIPaddVar(scip, var) );

      SCIPdebugMsg(scip, "creating negated variable <%s> (of <%s>) ...\n", SCIPvarGetName(var), SCIPvarGetName(negvar) );
      SCIPdebug( SCIP_CALL( SCIPprintVar(scip, var, NULL) ) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "neg_%s", SCIPvarGetName(var) );
      vars[0] = var;
      vars[1] = negvar;

      /* add exact linear constraint for negation */
      if( SCIPisExact(scip) )
      {
         SCIP_RATIONAL** vals;

         SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &vals, 2) );

         SCIPrationalSetReal(vals[0], 1.0);
         SCIPrationalSetReal(vals[1], 1.0);
         SCIP_CALL( SCIPcreateConsExactLinear(scip, &lincons, name, 2, vars, vals, vals[0], vals[0], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );

         SCIPrationalFreeBufferArray(SCIPbuffer(scip), &vals, 2);
      }
      /* add real linear constraint for negation */
      else
      {
         SCIP_Real vals[2];

         vals[0] = 1.0;
         vals[1] = 1.0;
         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, 2, vars, vals, 1.0, 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
      }

      SCIPdebugMsg(scip, "coupling constraint:\n");
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   }
   else if( strncmp(buf, "aggregated:", 11) == 0 )
   {
      /* handle (multi-)aggregated variables */
      SCIP_CONS* lincons = NULL;
      SCIP_VAR** vars;
      const char* str;
      int nvarssize = 20;
      int requsize;
      int nvars;

      buf += 11;

      /* special handling of variables that seem to be slack variables of indicator constraints */
      str = SCIPvarGetName(var);
      if( strncmp(str, "indslack", 8) == 0 )
      {
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "indlin");
         (void)strncat(name, str+8, SCIP_MAXSTRLEN-7);
      }
      else if( strncmp(str, "t_indslack", 10) == 0 )
      {
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "indlin");
         (void)strncat(name, str+10, SCIP_MAXSTRLEN-7);
      }
      else
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPvarGetName(var) );

      SCIPdebugMsg(scip, "parsing aggregated variable <%s> ...\n", SCIPvarGetName(var));

      /* add exact linear constraint for (multi-)aggregation */
      if( SCIPisExact(scip) )
      {
         SCIP_RATIONAL** vals;
         SCIP_RATIONAL* rhs;

         SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &rhs) );

         /* parse exact constant */
         if( !SCIPparseRational(scip, buf, rhs, &endptr) )
         {
            SCIPerrorMessage("expected constant when aggregated variable information (line: %d):\n%s\n", cipinput->linenumber, buf);
            cipinput->haserror = TRUE;
            SCIPrationalFreeBuffer(SCIPbuffer(scip), &rhs);
            return SCIP_OKAY;
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvarssize) );
         SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &vals, nvarssize) );

         /* check whether constant is 0.0 */
         str = endptr;
         SCIP_CALL( SCIPskipSpace((char**)&str) );

         /* if next char is '<' we found a variable -> constant is 0 */
         if( *str != '<' )
         {
            buf = endptr;

            SCIPrationalDebugMessage("constant: %q\n", rhs);

            SCIPrationalMultReal(rhs, rhs, -1.0);
         }
         /* otherwise keep buf */
         else
            SCIPrationalSetReal(rhs, 0.0);

         vars[0] = var;
         SCIPrationalSetReal(vals[0], -1.0);
         --nvarssize;

         /* parse exact linear sum to get variables and coefficients */
         SCIP_CALL( SCIPparseVarsLinearsumExact(scip, buf, vars + 1, vals + 1, &nvars, nvarssize, &requsize, &endptr, &success) );
         if( success && requsize > nvarssize )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requsize + 1) );
            SCIP_CALL( SCIPrationalReallocBufferArray(SCIPbuffer(scip), &vals, nvarssize + 1, requsize + 1) );
            nvarssize = requsize;
            SCIP_CALL( SCIPparseVarsLinearsumExact(scip, buf, vars + 1, vals + 1, &nvars, nvarssize, &requsize, &endptr, &success) );
            assert(!success || requsize <= nvarssize);
         }

         if( success )
         {
            /* add aggregation constraint */
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPcreateConsExactLinear(scip, &lincons, name, nvars + 1, vars, vals, rhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         }

         SCIPrationalFreeBufferArray(SCIPbuffer(scip), &vals, nvarssize + 1);
         SCIPfreeBufferArray(scip, &vars);
         SCIPrationalFreeBuffer(SCIPbuffer(scip), &rhs);
      }
      /* add real linear constraint for (multi-)aggregation */
      else
      {
         SCIP_Real* vals;
         SCIP_Real rhs;

         /* parse real constant */
         if( !SCIPparseReal(scip, buf, &rhs, &endptr) )
         {
            SCIPerrorMessage("expected constant when aggregated variable information (line: %d):\n%s\n", cipinput->linenumber, buf);
            cipinput->haserror = TRUE;
            return SCIP_OKAY;
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvarssize) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvarssize) );

         /* check whether constant is 0.0 */
         str = endptr;
         SCIP_CALL( SCIPskipSpace((char**)&str) );

         /* if next char is '<' we found a variable -> constant is 0 */
         if( *str != '<' )
         {
            buf = endptr;

            SCIPdebugMsg(scip, "constant: %f\n", rhs);

            rhs *= -1.0;
         }
         /* otherwise keep buf */
         else
            rhs = 0.0;

         vars[0] = var;
         vals[0] = -1.0;
         --nvarssize;

         /* parse linear sum to get variables and coefficients */
         SCIP_CALL( SCIPparseVarsLinearsum(scip, buf, vars + 1, vals + 1, &nvars, nvarssize, &requsize, &endptr, &success) );
         if( success && requsize > nvarssize )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requsize + 1) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &vals, requsize + 1) );
            nvarssize = requsize;
            SCIP_CALL( SCIPparseVarsLinearsum(scip, buf, vars + 1, vals + 1, &nvars, nvarssize, &requsize, &endptr, &success) );
            assert(!success || requsize <= nvarssize);
         }

         if( success )
         {
            /* add aggregation constraint */
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, nvars + 1, vars, vals, rhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         }

         SCIPfreeBufferArray(scip, &vals);
         SCIPfreeBufferArray(scip, &vars);
      }

      if( success )
      {
         SCIPdebugMsg(scip, "coupling constraint:\n");
         SCIPdebugPrintCons(scip, lincons, NULL);
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      }
      else
      {
         SCIPwarningMessage(scip, "Could not read (multi-)aggregated variable <%s>: dependent variables unkown - consider changing the order (line: %d):\n%s\n",
            SCIPvarGetName(var), cipinput->linenumber, buf);
      }
   }
   else
   {
      SCIPerrorMessage("unknown section when parsing variables (line: %d):\n%s\n", cipinput->linenumber, buf);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** read constraint */
static
SCIP_RETCODE getConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput,           /**< CIP parsing data */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             dynamic,            /**< Is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONS* cons;
   char* buf;
   char* copybuf;
   SCIP_RETCODE retcode;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool success;
   int len;

   buf = cipinput->strbuf;

   if( strncmp(buf, "END", 3) == 0 )
   {
      cipinput->section = CIP_END;
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "parse constraints in line %d\n", cipinput->linenumber);

   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;

   /* get length of line and check for correct ending of constraint line */
   len = (int)strlen(buf);
   if( len < 1 )
   {
      SCIPerrorMessage("syntax error: expected constraint in line %d.\n", cipinput->linenumber);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }
   if ( buf[len - 1] != ';' )
   {
      SCIPerrorMessage("syntax error: line has to end with ';' (line: %d)\n", cipinput->linenumber);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }

   /* copy buffer for working purpose */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &copybuf, buf, len) );
   copybuf[len - 1] = '\0';

   /* parse the constraint */
   retcode = SCIPparseCons(scip, &cons, copybuf,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE, &success);

   /* free temporary buffer */
   SCIPfreeBufferArray(scip, &copybuf);

   SCIP_CALL( retcode );

   if( !success )
   {
      SCIPerrorMessage("syntax error when reading constraint (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCip)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);

   SCIP_STRINGEQ( SCIPreaderGetName(reader), READER_NAME, SCIP_INVALIDCALL );

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderCip(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeCip)
{
   SCIP_READERDATA* readerdata;

   SCIP_STRINGEQ( SCIPreaderGetName(reader), READER_NAME, SCIP_INVALIDCALL );

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCip)
{  /*lint --e{715}*/
   CIPINPUT cipinput;
   SCIP_Bool initialconss;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamiccols;
   SCIP_Bool dynamicrows;
   SCIP_Bool initialvar;
   SCIP_Bool removablevar;
   SCIP_RETCODE retcode = SCIP_OKAY;

   if( NULL == (cipinput.file = SCIPfopen(filename, "r")) )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   cipinput.len = 131071;
   SCIP_CALL( SCIPallocBufferArray(scip, &(cipinput.strbuf), cipinput.len) );

   cipinput.linenumber = 0;
   cipinput.section = CIP_START;
   cipinput.haserror = FALSE;
   cipinput.endfile = FALSE;
   cipinput.readingsize = 65535;
   cipinput.objoffset = 0.0;
   cipinput.objscale = 1.0;
   if( SCIPisExact(scip) )
   {
      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &cipinput.objoffsetexact) );
      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &cipinput.objscaleexact) );

      SCIPrationalSetReal(cipinput.objoffsetexact, 0.0);
      SCIPrationalSetReal(cipinput.objscaleexact, 1.0);
   }
   else
   {
      cipinput.objoffsetexact = NULL;
      cipinput.objscaleexact = NULL;
   }

   SCIP_CALL_TERMINATE( retcode, SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL), TERMINATE );

   SCIP_CALL_TERMINATE( retcode, SCIPgetBoolParam(scip, "reading/initialconss", &initialconss), TERMINATE );
   SCIP_CALL_TERMINATE( retcode, SCIPgetBoolParam(scip, "reading/dynamiccols", &dynamiccols), TERMINATE );
   SCIP_CALL_TERMINATE( retcode, SCIPgetBoolParam(scip, "reading/dynamicconss", &dynamicconss), TERMINATE );
   SCIP_CALL_TERMINATE( retcode, SCIPgetBoolParam(scip, "reading/dynamicrows", &dynamicrows), TERMINATE );

   initialvar = !dynamiccols;
   removablevar = dynamiccols;

   while( cipinput.section != CIP_END && !cipinput.haserror )
   {
      /* get next input string */
      SCIP_CALL_TERMINATE( retcode, getInputString(scip, &cipinput), TERMINATE );

      if( cipinput.endfile )
         break;

      switch( cipinput.section )
      {
      case CIP_START:
         getStart(scip, &cipinput);
         break;
      case CIP_STATISTIC:
         SCIP_CALL_TERMINATE( retcode, getStatistics(scip, &cipinput), TERMINATE );
         break;
      case CIP_OBJECTIVE:
         SCIP_CALL_TERMINATE( retcode, getObjective(scip, &cipinput), TERMINATE );
         break;
      case CIP_VARS:
         SCIP_CALL_TERMINATE( retcode, getVariable(scip, &cipinput, initialvar, removablevar), TERMINATE );
         break;
      case CIP_FIXEDVARS:
         SCIP_CALL_TERMINATE( retcode, getFixedVariable(scip, &cipinput), TERMINATE );
         break;
      case CIP_CONSTRAINTS:
         SCIP_CALL_TERMINATE( retcode, getConstraint(scip, &cipinput, initialconss, dynamicconss, dynamicrows), TERMINATE );
         break;
      default:
         SCIPerrorMessage("invalid CIP state\n");
         SCIPABORT();
         retcode = SCIP_INVALIDDATA;  /*lint !e527*/
         goto TERMINATE;
      } /*lint !e788*/
   }

   if( cipinput.haserror )
      goto TERMINATE;

   /* offset exact objective */
   if( cipinput.objoffsetexact != NULL )
   {
      if( !SCIPrationalIsZero(cipinput.objoffsetexact) )
      {
         SCIPrationalMult(cipinput.objoffsetexact, cipinput.objoffsetexact, cipinput.objscaleexact);
         SCIP_CALL_TERMINATE( retcode, SCIPaddOrigObjoffsetExact(scip, cipinput.objoffsetexact), TERMINATE );
      }
   }
   /* offset real objective */
   else
   {
      if( cipinput.objoffset != 0.0 ) /*lint !e777*/
      {
         cipinput.objoffset *= cipinput.objscale;
         SCIP_CALL_TERMINATE( retcode, SCIPaddOrigObjoffset(scip, cipinput.objoffset), TERMINATE );
      }
   }

   if( cipinput.section != CIP_END )
   {
      SCIPerrorMessage("unexpected EOF\n");
      cipinput.haserror = TRUE;
   }

 TERMINATE:
   /* close file stream */
   SCIPfclose(cipinput.file);

   if( cipinput.objscaleexact != NULL )
      SCIPrationalFreeBuffer(SCIPbuffer(scip), &cipinput.objscaleexact);
   if( cipinput.objoffsetexact != NULL )
      SCIPrationalFreeBuffer(SCIPbuffer(scip), &cipinput.objoffsetexact);
   SCIPfreeBufferArray(scip, &cipinput.strbuf);

   if( cipinput.haserror || retcode == SCIP_INVALIDDATA )
      return SCIP_READERROR;

   /* successfully parsed cip format */
   if( retcode == SCIP_OKAY )
      *result = SCIP_SUCCESS;

   return retcode;
}

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{  /*lint --e{715}*/
   if( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{  /*lint --e{715}*/
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCip)
{  /*lint --e{715}*/
   SCIP_HASHTABLE* varhash = NULL;
   SCIP_READERDATA* readerdata;
   int i;

   assert(reader != NULL);

   SCIP_STRINGEQ( SCIPreaderGetName(reader), READER_NAME, SCIP_INVALIDCALL );

   SCIPinfoMessage(scip, file, "STATISTICS\n");
   SCIPinfoMessage(scip, file, "  Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "  Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "  Constraints      : %d initial, %d maximal\n", startnconss, maxnconss);

   SCIPinfoMessage(scip, file, "OBJECTIVE\n");
   SCIPinfoMessage(scip, file, "  Sense            : %s\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");
   /* write exact offset */
   if( objoffsetexact != NULL )
   {
      assert(SCIPisExact(scip));

      if( !SCIPrationalIsZero(objoffsetexact) )
      {
         SCIPinfoMessage(scip, file, "  Offset           : ");
         SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, objoffsetexact);
         SCIPinfoMessage(scip, file, "\n");
      }
   }
   /* write real offset */
   else
   {
      if( objoffset != 0.0 ) /*lint !e777*/
         SCIPinfoMessage(scip, file, "  Offset           : %+.15g\n", objoffset);
   }
   /* write exact scale */
   if( objscaleexact != NULL )
   {
      assert(SCIPisExact(scip));

      if( !SCIPrationalIsEQReal(objscaleexact, 1.0) )
      {
         SCIPinfoMessage(scip, file, "  Scale            : ");
         SCIPrationalMessage(SCIPgetMessagehdlr(scip), file, objscaleexact);
         SCIPinfoMessage(scip, file, "\n");
      }
   }
   /* write real scale */
   else
   {
      if( objscale != 1.0 ) /*lint !e777*/
         SCIPinfoMessage(scip, file, "  Scale            : %.15g\n", objscale);
   }

   if ( nfixedvars > 0 )
   {
      /* set up hash table for variables that have been written property (used for writing out fixed vars in the right order) */
      SCIP_CALL( SCIPhashtableCreate(&varhash, SCIPblkmem(scip), nvars + nfixedvars, hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   }

   if ( nvars + nfixedvars > 0 )
   {
      SCIPinfoMessage(scip, file, "VARIABLES\n");
   }

   if( nvars > 0 )
   {
      for( i = 0; i < nvars; ++i )
      {
         SCIP_VAR* var;

         var = vars[i];
         assert( var != NULL );
         SCIP_CALL( SCIPprintVar(scip, var, file) );
         if ( varhash != NULL )
         {
            /* add free variable to hashtable */
            if ( ! SCIPhashtableExists(varhash, (void*) var) )
            {
               SCIP_CALL( SCIPhashtableInsert(varhash, (void*) var) );
            }
         }
      }
   }

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   if( readerdata->writefixedvars && nfixedvars > 0 )
   {
      int nwritten = 0;

      SCIPinfoMessage(scip, file, "FIXED\n");

      /* loop through variables until each has been written after the variables that it depends on have been written; this
       * requires several runs over the variables, but the depth (= number of loops) is usually small. */
      while ( nwritten < nfixedvars )
      {
         SCIPdebugMsg(scip, "written %d of %d fixed variables.\n", nwritten, nfixedvars);
         for (i = 0; i < nfixedvars; ++i)
         {
            SCIP_VAR* var;
            SCIP_VAR* tmpvar;

            var = fixedvars[i];
            assert( var != NULL );

            /* skip variables already written */
            if ( SCIPhashtableExists(varhash, (void*) var) )
               continue;

            switch ( SCIPvarGetStatus(var) )
            {
            case SCIP_VARSTATUS_FIXED:

               /* fixed variables can simply be output and added to the hashtable */
               SCIP_CALL( SCIPprintVar(scip, var, file) );
               assert( ! SCIPhashtableExists(varhash, (void*) var) );
               SCIP_CALL( SCIPhashtableInsert(varhash, (void*) var) );
               ++nwritten;

               break;

            case SCIP_VARSTATUS_NEGATED:

               tmpvar = SCIPvarGetNegationVar(var);
               assert( tmpvar != NULL );
               assert( var == SCIPvarGetNegatedVar(tmpvar) );

               /* if the negated variable has been written, we can write the current variable */
               if ( SCIPhashtableExists(varhash, (void*) tmpvar) )
               {
                  SCIP_CALL( SCIPprintVar(scip, var, file) );
                  assert( ! SCIPhashtableExists(varhash, (void*) var) );
                  SCIP_CALL( SCIPhashtableInsert(varhash, (void*) var) );
                  ++nwritten;
               }
               break;

            case SCIP_VARSTATUS_AGGREGATED:

               tmpvar = SCIPvarGetAggrVar(var);
               assert( tmpvar != NULL );

               /* if the aggregating variable has been written, we can write the current variable */
               if ( SCIPhashtableExists(varhash, (void*) tmpvar) )
               {
                  SCIP_CALL( SCIPprintVar(scip, var, file) );
                  assert( ! SCIPhashtableExists(varhash, (void*) var) );
                  SCIP_CALL( SCIPhashtableInsert(varhash, (void*) var) );
                  ++nwritten;
               }
               break;

            case SCIP_VARSTATUS_MULTAGGR:
            {
               SCIP_VAR** aggrvars;
               int naggrvars;
               int j;

               /* get the active representation */
               SCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );

               naggrvars = SCIPvarGetMultaggrNVars(var);
               aggrvars = SCIPvarGetMultaggrVars(var);
               assert(aggrvars != NULL || naggrvars == 0);

               for (j = 0; j < naggrvars; ++j)
               {
                  if( !SCIPhashtableExists(varhash, (void*) aggrvars[j]) ) /*lint !e613*/
                     break;
               }

               /* if all multi-aggregating variables have been written, we can write the current variable */
               if ( j >= naggrvars )
               {
                  SCIP_CALL( SCIPprintVar(scip, var, file) );
                  assert( ! SCIPhashtableExists(varhash, (void*) var) );
                  SCIP_CALL( SCIPhashtableInsert(varhash, (void*) var) );
                  ++nwritten;
               }
               break;
            }

            case SCIP_VARSTATUS_ORIGINAL:
            case SCIP_VARSTATUS_LOOSE:
            case SCIP_VARSTATUS_COLUMN:
               SCIPerrorMessage("Only fixed variables are allowed to be present in fixedvars list.\n");
               SCIPABORT();
               return SCIP_ERROR; /*lint !e527*/
            }
         }
      }
   }

   if( nconss > 0 )
   {
      SCIPinfoMessage(scip, file, "CONSTRAINTS\n");

      for( i = 0; i < nconss; ++i )
      {
         SCIP_CALL( SCIPprintCons(scip, conss[i], file) );
         SCIPinfoMessage(scip, file, ";\n");
      }
   }
   SCIPinfoMessage(scip, file, "END\n");

   *result = SCIP_SUCCESS;

   if( nfixedvars > 0 )
      SCIPhashtableFree(&varhash);
   else
      assert(varhash == NULL);

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the cip file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create cip reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* reader is safe to use in exact solving mode */
   SCIPreaderMarkExact(reader);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCip) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeCip) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCip) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCip) );

   /* add cip reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/cipreader/writefixedvars", "should fixed and aggregated variables be printed (if not, re-parsing might fail)",
         &readerdata->writefixedvars, FALSE, DEFAULT_CIP_WRITEFIXEDVARS, NULL, NULL) );

   return SCIP_OKAY;
}
