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

/**@file   reader_cip.c
 * @brief  CIP file reader
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 *
 * @todo Ensure order of fixed/(multi-)aggregated/negated variables in output to ensure correct reading; all used
 * variables in aggregations have to be defined previously.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <ctype.h>

#include "scip/reader_cip.h"
#include "scip/cons_linear.h"

#define READER_NAME             "cipreader"
#define READER_DESC             "file reader for CIP (Constraint Integer Program) format"
#define READER_EXTENSION        "cip"

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
   SCIP_Bool             aggregatedvars;     /**< whether there are (multi-)aggregated variables in the input */
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
         pos = endline - cipinput->strbuf;

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
      endline = strrchr(cipinput->strbuf, '\n');
      endcharacter = strchr(cipinput->strbuf, ';');
   }
   assert(endline != NULL);

   /*SCIPdebugMessage("read line: %s\n", cipinput->strbuf);*/

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
SCIP_RETCODE getStart(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput            /**< CIP parsing data */
   )
{
   if( strncmp(cipinput->strbuf, "STATISTICS", 9) == 0 )
   {
      cipinput->section = CIP_STATISTIC;
      return SCIP_OKAY;
   }
   
   return SCIP_OKAY;
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

   SCIPdebugMessage("parse statistics\n");

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
      while( isspace((unsigned char)* name) )
         ++name;

      /* set problem name */
      SCIP_CALL( SCIPsetProbName(scip, name) );

      SCIPdebugMessage("problem name <%s>\n", name);
   }

   return SCIP_OKAY;
}

/** read objective sense, offset, and scale */
static
SCIP_RETCODE getObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput,           /**< CIP parsing data */
   SCIP_Real*            objscale,           /**< buffer where to multiply with objective scale */
   SCIP_Real*            objoffset           /**< buffer where to add with objective offset */
   )
{
   char* buf;
   char* name;

   assert(objscale != NULL);
   assert(objoffset != NULL);

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

   SCIPdebugMessage("parse objective information\n");

   if( strncmp(buf, "  Sense", 7) == 0 )
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
      while( isspace((unsigned char)* name) )
         ++name;

      if( strncmp(name, "minimize", 3) == 0 )
         objsense = SCIP_OBJSENSE_MINIMIZE;
      else if( strncmp(name, "maximize", 3) == 0 )
         objsense = SCIP_OBJSENSE_MAXIMIZE;
      else
      {
         SCIPwarningMessage(scip, "unknown objective sense '%s' (line: %d):\n%s\n", name, cipinput->linenumber, cipinput->strbuf);
         return SCIP_OKAY; /* no error - might work with default */
      }

      /* set problem name */
      SCIP_CALL( SCIPsetObjsense(scip, objsense) );
      SCIPdebugMessage("objective sense <%s>\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");
   }
   else if( strncmp(buf, "  Offset", 7) == 0 )
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

      /* remove white space in front of the name */
      while(isspace((unsigned char)*name))
         name++;

      *objoffset += strtod(name, &endptr);
   }
   else if( strncmp(buf, "  Scale", 7) == 0 )
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

      /* remove white space in front of the name */
      while(isspace((unsigned char)*name))
         name++;

      *objscale *= strtod(name, &endptr);
   }

   return SCIP_OKAY;
}

/** read variable */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   CIPINPUT*             cipinput,           /**< CIP parsing data */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_Real             objscale            /**< objective scale */
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
   
   SCIPdebugMessage("parse variable\n");

   /* parse the variable */
   SCIP_CALL( SCIPparseVar(scip, &var, buf, initial, removable, NULL, NULL, NULL, NULL, NULL, &endptr, &success) );

   if( !success )
   {
      SCIPerrorMessage("syntax error in variable information (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
      cipinput->haserror = TRUE;
      return SCIP_OKAY;
   }
   
   if( objscale != 1.0 )
   {
      SCIP_CALL( SCIPchgVarObj(scip, var, SCIPvarGetObj(var) * objscale) );
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

   SCIPdebugMessage("parse fixed variables\n");

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

   while ( *buf != '\0' && (*buf == ' ' || *buf == ',') )
      ++buf;

   /* check whether variable is fixed */
   if ( strncmp(buf, "fixed:", 6) == 0 )
   {
      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIPdebug( SCIP_CALL( SCIPprintVar(scip, var, NULL) ) );
   }
   else if ( strncmp(buf, "negated:", 8) == 0 )
   {
      SCIP_CONS* lincons;
      SCIP_VAR* negvar;
      SCIP_Real vals[2];
      SCIP_VAR* vars[2];

      buf += 8;

      /* we can just parse the next variable (ignoring all other information in between) */
      SCIP_CALL( SCIPparseVarName(scip, buf, &negvar, &endptr) );

      if ( negvar == NULL )
      {
         SCIPerrorMessage("could not parse negated variable (line: %d):\n%s\n", cipinput->linenumber, cipinput->strbuf);
         cipinput->haserror = TRUE;
         return SCIP_OKAY;
      }
      assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
      assert( SCIPvarGetType(negvar) == SCIP_VARTYPE_BINARY );

      SCIP_CALL( SCIPaddVar(scip, var) );

      SCIPdebugMessage("creating negated variable <%s> (of <%s>) ...\n", SCIPvarGetName(var), SCIPvarGetName(negvar) );
      SCIPdebug( SCIP_CALL( SCIPprintVar(scip, var, NULL) ) );

      /* add linear constraint for negation */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "neg_%s", SCIPvarGetName(var) );
      vars[0] = var;
      vars[1] = negvar;
      vals[0] = 1.0;
      vals[1] = 1.0;
      SCIPdebugMessage("coupling constraint:\n");
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, 2, vars, vals, 1.0, 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   }
   else
   {
      if ( ! cipinput->aggregatedvars )
      {
         cipinput->aggregatedvars = TRUE;
         SCIPwarningMessage(scip, "the CIP-input contains (multi-)aggregated variables - this is not supported yet. Note that this might lead to parsing errors later.\n");
      }
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
   int len;

   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;

   SCIP_Bool success;

   buf = cipinput->strbuf;

   if( strncmp(buf, "END", 3) == 0 )
   {
      cipinput->section = CIP_END;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("parse constraints in line %d\n", cipinput->linenumber);

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
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &copybuf, buf, len) );
   copybuf[len - 1] = '\0';

   /* parse the constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, copybuf,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE, &success) );

   /* free temporary buffer */
   SCIPfreeMemoryArray(scip, &copybuf);

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
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderCip(scip) );
 
   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCip)
{  /*lint --e{715}*/

   CIPINPUT cipinput;
   SCIP_Real objscale;
   SCIP_Real objoffset;
   SCIP_Bool initialconss;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamiccols;
   SCIP_Bool dynamicrows;

   SCIP_Bool initialvar;
   SCIP_Bool removablevar;

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
   cipinput.aggregatedvars = FALSE;

   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &initialconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &dynamiccols) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &dynamicrows) );

   initialvar = !dynamiccols;
   removablevar = dynamiccols;

   objscale = 1.0;
   objoffset = 0.0;

   while( cipinput.section != CIP_END && !cipinput.haserror )
   {
      /* get next input string */
      SCIP_CALL( getInputString(scip, &cipinput) );

      if( cipinput.endfile )
         break;
      
      switch( cipinput.section )
      {
      case CIP_START:
         SCIP_CALL( getStart(scip, &cipinput) );
         break;
      case CIP_STATISTIC:
         SCIP_CALL( getStatistics(scip, &cipinput) );
         break;
      case CIP_OBJECTIVE:
         SCIP_CALL( getObjective(scip, &cipinput, &objscale, &objoffset) );
         break;
      case CIP_VARS:
         SCIP_CALL( getVariable(scip, &cipinput, initialvar, removablevar, objscale) );
         break;
      case CIP_FIXEDVARS:
         SCIP_CALL( getFixedVariable(scip, &cipinput) );
         break;
      case CIP_CONSTRAINTS:
         SCIP_CALL( getConstraint(scip, &cipinput, initialconss, dynamicconss, dynamicrows) );
         break;
      default:
         SCIPerrorMessage("invalid CIP state\n");
         SCIPABORT();
      } /*lint !e788*/ 
   }

   if( !SCIPisZero(scip, objoffset) && !cipinput.haserror )
   {
      SCIP_VAR* objoffsetvar;

      objoffset *= objscale;
      SCIP_CALL( SCIPcreateVar(scip, &objoffsetvar, "objoffset", objoffset, objoffset, 1.0, SCIP_VARTYPE_CONTINUOUS,
         TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, objoffsetvar) );
      SCIP_CALL( SCIPreleaseVar(scip, &objoffsetvar) );
      SCIPdebugMessage("added variables <objoffset> for objective offset of <%g>\n", objoffset);
   }

   /* close file stream */
   SCIPfclose(cipinput.file);
   
   if( cipinput.section != CIP_END && !cipinput.haserror )
   {
      SCIPerrorMessage("unexpected EOF\n");
   }

   SCIPfreeBufferArray(scip, &cipinput.strbuf);

   if( cipinput.haserror )
      return SCIP_READERROR;

   /* successfully parsed cip format */
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCip)
{  /*lint --e{715}*/
   int i;

   SCIPinfoMessage(scip, file, "STATISTICS\n");
   SCIPinfoMessage(scip, file, "  Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "  Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "  Constraints      : %d initial, %d maximal\n", startnconss, maxnconss);

   SCIPinfoMessage(scip, file, "OBJECTIVE\n");
   SCIPinfoMessage(scip, file, "  Sense            : %s\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");
   if( !SCIPisZero(scip, objoffset) )
      SCIPinfoMessage(scip, file, "  Offset           : %+.15g\n", objoffset);
   if( !SCIPisEQ(scip, objscale, 1.0) )
      SCIPinfoMessage(scip, file, "  Scale            : %.15g\n", objscale);

   if( nvars > 0 )
   {
      SCIPinfoMessage(scip, file, "VARIABLES\n");
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPprintVar(scip, vars[i], file) );
      }
   }

   if( nfixedvars > 0 )
   {
      SCIPinfoMessage(scip, file, "FIXED\n");

      /* first print fixed variables to increase chance that variables are defined for aggregated variables */
      for( i = 0; i < nfixedvars; ++i )
      {
         assert( SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_AGGREGATED ||
            SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_NEGATED );

         if (SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_FIXED )
            SCIP_CALL( SCIPprintVar(scip, fixedvars[i], file) );
      }

      /* next print aggregated or negated variables */
      for( i = 0; i < nfixedvars; ++i )
      {
         if ( SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_AGGREGATED || SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_NEGATED )
            SCIP_CALL( SCIPprintVar(scip, fixedvars[i], file) );
      }

      /* finally print multi-aggregated variables */
      for( i = 0; i < nfixedvars; ++i )
      {
         if ( SCIPvarGetStatus(fixedvars[i]) == SCIP_VARSTATUS_MULTAGGR )
            SCIP_CALL( SCIPprintVar(scip, fixedvars[i], file) );
      }
   }

   if( nconss > 0 )
   {
      SCIPinfoMessage(scip, file, "CONSTRAINTS\n");

      for( i = 0; i < nconss; ++i )
      {
         /* in case the transformed is written only constraint are posted which are enabled in the current node */
         assert(!transformed || SCIPconsIsEnabled(conss[i]));

         SCIP_CALL( SCIPprintCons(scip, conss[i], file) );
         SCIPinfoMessage(scip, file, ";\n");
      }
   }

   *result = SCIP_SUCCESS;

   SCIPinfoMessage(scip, file, "END\n");
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

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCip) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCip) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCip) );

   return SCIP_OKAY;
}
