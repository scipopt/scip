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
#pragma ident "@(#) $Id: reader_pip.c,v 1.18 2011/02/24 19:57:43 bzfwinkm Exp $"

/**@file   reader_pip.c
 * @ingroup FILEREADERS 
 * @brief  file reader for polynomial mixed-integer programs in PIP format
 * @author Stefan Vigerske
 *
 * @todo Test for uniqueness of variable names (after cutting down).
 * 
 * reader_lp has been used as a starting point for this reader.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h>
#endif
#include <ctype.h>

#include "scip/reader_pip.h"
#include "scip/cons_linear.h"
#include "scip/cons_quadratic.h"
#include "scip/pub_misc.h"

#define READER_NAME             "pipreader"
#define READER_DESC             "file reader for polynomial mixed-integer programs in PIP format"
#define READER_EXTENSION        "pip"


/*
 * Data structures
 */
#define PIP_MAX_LINELEN        65536
#define PIP_MAX_PUSHEDTOKENS   2
#define PIP_INIT_VARSSIZE      256
#define PIP_INIT_MONOMIALSSIZE 128
#define PIP_INIT_FACTORSSIZE   16

/** Section in PIP File */
enum PipSection
{
   PIP_START, PIP_OBJECTIVE, PIP_CONSTRAINTS, PIP_BOUNDS, PIP_GENERALS, PIP_BINARIES, PIP_END
};
typedef enum PipSection PIPSECTION;

enum PipExpType
{
   PIP_EXP_NONE, PIP_EXP_UNSIGNED, PIP_EXP_SIGNED
};
typedef enum PipExpType PIPEXPTYPE;

enum PipSense
{
   PIP_SENSE_NOTHING, PIP_SENSE_LE, PIP_SENSE_GE, PIP_SENSE_EQ
};
typedef enum PipSense PIPSENSE;

/** PIP reading data */
struct PipInput
{
   SCIP_FILE*           file;
   char                 linebuf[PIP_MAX_LINELEN+1];
   char                 probname[PIP_MAX_LINELEN];
   char                 objname[PIP_MAX_LINELEN];
   char*                token;
   char*                tokenbuf;
   char*                pushedtokens[PIP_MAX_PUSHEDTOKENS];
   int                  npushedtokens;
   int                  linenumber;
   int                  linepos;
   PIPSECTION           section;
   SCIP_OBJSENSE        objsense;
   SCIP_Bool            haserror;
};
typedef struct PipInput PIPINPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=*^";
static const char commentchars[] = "\\";




/*
 * Local methods (for reading)
 */

/** issues an error message and marks the PIP data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   const char*           msg                 /**< error message */
   )
{
   char formatstr[256];

   assert(pipinput != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error in line %d: %s ('%s')\n",
      pipinput->linenumber, msg, pipinput->token);
   if( pipinput->linebuf[strlen(pipinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", pipinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", pipinput->linebuf);
   }
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", pipinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, formatstr, "^");
   pipinput->section  = PIP_END;
   pipinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   return pipinput->haserror;
}

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char                  c                   /**< input character */
   )
{
   return (c == '\0') || (strchr(delimchars, c) != NULL);
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char                  c                   /**< input character */
   )
{
   return (strchr(tokenchars, c) != NULL);
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char                  c,                  /**< input character */
   char                  nextc,              /**< next input character */
   SCIP_Bool             firstchar,          /**< is the given character the first char of the token? */
   SCIP_Bool*            hasdot,             /**< pointer to update the dot flag */
   PIPEXPTYPE*           exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit(c) )
      return TRUE;
   else if( (*exptype == PIP_EXP_NONE) && !(*hasdot) && (c == '.') && isdigit(nextc) )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == PIP_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = PIP_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit(nextc) )
      {
         *exptype = PIP_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == PIP_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = PIP_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   int i;

   assert(pipinput != NULL);

   /* clear the line */
   BMSclearMemoryArray(pipinput->linebuf, PIP_MAX_LINELEN);

   /* read next line */
   pipinput->linepos = 0;
   pipinput->linebuf[PIP_MAX_LINELEN-2] = '\0';
   if( SCIPfgets(pipinput->linebuf, sizeof(pipinput->linebuf), pipinput->file) == NULL )
      return FALSE;
   pipinput->linenumber++;
   if( pipinput->linebuf[PIP_MAX_LINELEN-2] != '\0' )
   {
      SCIPerrorMessage("Error: line %d exceeds %d characters\n", pipinput->linenumber, PIP_MAX_LINELEN-2);
      pipinput->haserror = TRUE;
      return FALSE;
   }
   pipinput->linebuf[PIP_MAX_LINELEN-1] = '\0';
   pipinput->linebuf[PIP_MAX_LINELEN-2] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(pipinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
      }
   }

   return TRUE;
}

/** swaps the addresses of two pointers */
static
void swapPointers(
   char**                pointer1,           /**< first pointer */
   char**                pointer2            /**< second pointer */
   )
{
   char* tmp;

   tmp = *pointer1;
   *pointer1 = *pointer2;
   *pointer2 = tmp;
}

/** reads the next token from the input file into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   SCIP_Bool hasdot;
   PIPEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(pipinput != NULL);
   assert(pipinput->linepos < PIP_MAX_LINELEN);

   /* check the token stack */
   if( pipinput->npushedtokens > 0 )
   {
      swapPointers(&pipinput->token, &pipinput->pushedtokens[pipinput->npushedtokens-1]);
      pipinput->npushedtokens--;
      SCIPdebugMessage("(line %d) read token again: '%s'\n", pipinput->linenumber, pipinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = pipinput->linebuf;
   while( isDelimChar(buf[pipinput->linepos]) )
   {
      if( buf[pipinput->linepos] == '\0' )
      {
         if( !getNextLine(pipinput) )
         {
            pipinput->section = PIP_END;
            SCIPdebugMessage("(line %d) end of file\n", pipinput->linenumber);
            return FALSE;
         }
         assert(pipinput->linepos == 0);
      }
      else
         pipinput->linepos++;
   }
   assert(pipinput->linepos < PIP_MAX_LINELEN);
   assert(!isDelimChar(buf[pipinput->linepos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = PIP_EXP_NONE;
   if( isValueChar(buf[pipinput->linepos], buf[pipinput->linepos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < PIP_MAX_LINELEN);
         assert(!isDelimChar(buf[pipinput->linepos]));
         pipinput->token[tokenlen] = buf[pipinput->linepos];
         tokenlen++;
         pipinput->linepos++;
      }
      while( isValueChar(buf[pipinput->linepos], buf[pipinput->linepos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < PIP_MAX_LINELEN);
         pipinput->token[tokenlen] = buf[pipinput->linepos];
         tokenlen++;
         pipinput->linepos++;
         if( tokenlen == 1 && isTokenChar(pipinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[pipinput->linepos]) && !isTokenChar(buf[pipinput->linepos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1
         && (pipinput->token[tokenlen-1] == '<' || pipinput->token[tokenlen-1] == '>' || pipinput->token[tokenlen-1] == '=')
         && buf[pipinput->linepos] == '=' )
      {
         pipinput->linepos++;
      }
      else if( pipinput->token[tokenlen-1] == '=' && (buf[pipinput->linepos] == '<' || buf[pipinput->linepos] == '>') )
      {
         pipinput->token[tokenlen-1] = buf[pipinput->linepos];
         pipinput->linepos++;
      }
   }
   assert(tokenlen < PIP_MAX_LINELEN);
   pipinput->token[tokenlen] = '\0';

   SCIPdebugMessage("(line %d) read token: '%s'\n", pipinput->linenumber, pipinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);
   assert(pipinput->npushedtokens < PIP_MAX_PUSHEDTOKENS);

   swapPointers(&pipinput->pushedtokens[pipinput->npushedtokens], &pipinput->token);
   pipinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);
   assert(pipinput->npushedtokens < PIP_MAX_PUSHEDTOKENS);

   swapPointers(&pipinput->pushedtokens[pipinput->npushedtokens], &pipinput->tokenbuf);
   pipinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   swapPointers(&pipinput->token, &pipinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   SCIP_Bool iscolon;

   assert(pipinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(pipinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(pipinput) )
   {
      iscolon = (strcmp(pipinput->token, ":") == 0);
      pushToken(pipinput);
   }

   /* reinstall the previous token by swapping back the token buffer */
   swapTokenBuffer(pipinput);

   /* check for ':' */
   if( iscolon )
      return FALSE;

   if( strcasecmp(pipinput->token, "MINIMIZE") == 0
      || strcasecmp(pipinput->token, "MINIMUM") == 0
      || strcasecmp(pipinput->token, "MIN") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: OBJECTIVE\n", pipinput->linenumber);
      pipinput->section = PIP_OBJECTIVE;
      pipinput->objsense = SCIP_OBJSENSE_MINIMIZE;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "MAXIMIZE") == 0
      || strcasecmp(pipinput->token, "MAXIMUM") == 0
      || strcasecmp(pipinput->token, "MAX") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: OBJECTIVE\n", pipinput->linenumber);
      pipinput->section = PIP_OBJECTIVE;
      pipinput->objsense = SCIP_OBJSENSE_MAXIMIZE;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "SUBJECT") == 0 )
   {
      /* check if the next token is 'TO' */
      swapTokenBuffer(pipinput);
      if( getNextToken(pipinput) )
      {
         if( strcasecmp(pipinput->token, "TO") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", pipinput->linenumber);
            pipinput->section = PIP_CONSTRAINTS;
            return TRUE;
         }
         else
            pushToken(pipinput);
      }
      swapTokenBuffer(pipinput);
   }

   if( strcasecmp(pipinput->token, "SUCH") == 0 )
   {
      /* check if the next token is 'THAT' */
      swapTokenBuffer(pipinput);
      if( getNextToken(pipinput) )
      {
         if( strcasecmp(pipinput->token, "THAT") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", pipinput->linenumber);
            pipinput->section = PIP_CONSTRAINTS;
            return TRUE;
         }
         else
            pushToken(pipinput);
      }
      swapTokenBuffer(pipinput);
   }

   if( strcasecmp(pipinput->token, "st") == 0
      || strcasecmp(pipinput->token, "S.T.") == 0
      || strcasecmp(pipinput->token, "ST.") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", pipinput->linenumber);
      pipinput->section = PIP_CONSTRAINTS;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "BOUNDS") == 0
      || strcasecmp(pipinput->token, "BOUND") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: BOUNDS\n", pipinput->linenumber);
      pipinput->section = PIP_BOUNDS;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "GENERAL") == 0
      || strcasecmp(pipinput->token, "GENERALS") == 0
      || strcasecmp(pipinput->token, "GEN") == 0
      || strcasecmp(pipinput->token, "INTEGER") == 0
      || strcasecmp(pipinput->token, "INTEGERS") == 0
      || strcasecmp(pipinput->token, "INT") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: GENERALS\n", pipinput->linenumber);
      pipinput->section = PIP_GENERALS;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "BINARY") == 0
      || strcasecmp(pipinput->token, "BINARIES") == 0
      || strcasecmp(pipinput->token, "BIN") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: BINARIES\n", pipinput->linenumber);
      pipinput->section = PIP_BINARIES;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "END") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: END\n", pipinput->linenumber);
      pipinput->section = PIP_END;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   PIPINPUT*             pipinput,           /**< PIP reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(pipinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( pipinput->token[1] == '\0' )
   {
      if( *pipinput->token == '+' )
         return TRUE;
      else if( *pipinput->token == '-' )
      {
         *sign *= -1;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isValue(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(pipinput != NULL);
   assert(value != NULL);

   if( strcasecmp(pipinput->token, "INFINITY") == 0 || strcasecmp(pipinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(pipinput->token, &endptr);
      if( endptr != pipinput->token && *endptr == '\0' )
      {
         *value = val;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is an equation sense */
static
SCIP_Bool isSense(
   PIPINPUT*              pipinput,           /**< PIP reading data */
   PIPSENSE*              sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(pipinput != NULL);

   if( strcmp(pipinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = PIP_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(pipinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = PIP_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(pipinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = PIP_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns the variable with the given name, or creates a new variable if it does not exist */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 name,               /**< name of the variable */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Bool*            created             /**< pointer to store whether a new variable was created, or NULL */
   )
{
   assert(name != NULL);
   assert(var != NULL);

   *var = SCIPfindVar(scip, name);
   if( *var == NULL )
   {
      SCIP_VAR* newvar;
      SCIP_Bool dynamiccols;
      SCIP_Bool initial;
      SCIP_Bool removable;

      SCIP_CALL( SCIPgetBoolParam(scip, "reading/pipreader/dynamiccols", &dynamiccols) );
      initial = !dynamiccols;
      removable = dynamiccols;

      /* create new variable of the given name */
      SCIPdebugMessage("creating new variable: <%s>\n", name);
      SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
            initial, removable, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, newvar) );
      *var = newvar;

      /* because the variable was added to the problem, it is captured by SCIP and we can safely release it right now
       * without making the returned *var invalid
       */
      SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

      if( created != NULL )
         *created = TRUE;
   }
   else if( created != NULL )
      *created = FALSE;

   return SCIP_OKAY;
}

/** reads the header of the file */
static
SCIP_RETCODE readStart(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( !getNextToken(pipinput) )
         return SCIP_OKAY;
   }
   while( !isNewSection(pipinput) );

   return SCIP_OKAY;
}

/** ensure that an array of monomials can hold a minimum number of entries */
static
SCIP_RETCODE ensureMonomialsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRDATA_MONOMIAL*** monomials,      /**< pointer to current array of monomials */
   int*                  monomialssize,      /**< current size of monomials array at input; new size at exit */
   int                   minnmonomials       /**< required minimal size of monomials array */
   )
{
   int newsize;
   
   assert(scip != NULL);
   assert(monomials != NULL);
   assert(monomialssize != NULL);
   assert(*monomials != NULL || *monomialssize == 0);
   
   if( minnmonomials <= *monomialssize )
      return SCIP_OKAY;
   
   newsize = SCIPcalcMemGrowSize(scip, minnmonomials);
   
   if( *monomials != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, monomials, newsize) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, monomials, newsize) );
   }
   *monomialssize = newsize;
   
   return SCIP_OKAY;
}

/** ensure that arrays of exponents and variable indices can hold a minimum number of entries */
static
SCIP_RETCODE ensureFactorsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           exponents,          /**< pointer to current array of exponents */
   int**                 varidxs,            /**< pointer to current array of variable indices */
   int*                  factorssize,        /**< current size of arrays at input; new size at exit */
   int                   minnfactors         /**< required minimal size of arrays */
   )
{
   int newsize;
   
   assert(scip != NULL);
   assert(exponents != NULL);
   assert(varidxs != NULL);
   assert(factorssize != NULL);
   assert(*exponents != NULL || *factorssize == 0);
   assert(*varidxs   != NULL || *factorssize == 0);
   assert((*exponents != NULL) == (*varidxs != NULL));
   
   if( minnfactors <= *factorssize )
      return SCIP_OKAY;
   
   newsize = SCIPcalcMemGrowSize(scip, minnfactors);
   
   if( *exponents != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, exponents, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, varidxs,   newsize) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, exponents, newsize) );
      SCIP_CALL( SCIPallocBufferArray(scip, varidxs,   newsize) );
   }
   *factorssize = newsize;
   
   return SCIP_OKAY;
}

/** gives index of variable in vars array, inserts it at the end if not existing yet */
static
SCIP_RETCODE getVariableIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to current array of variables */
   int*                  varssize,           /**< current size of variables array at input; new size at exit */
   int*                  nvars,              /**< number of variables stored in array */
   SCIP_HASHMAP*         varhash,            /**< hashmap variables -> indices */
   SCIP_VAR*             var,                /**< the variable which index we need */
   int*                  varidx              /**< pointer to store index of variable in *vars */
   )
{
   assert(scip != NULL);
   assert(varssize != NULL);
   assert(vars != NULL || *varssize == 0);
   assert(nvars != NULL);
   assert(*nvars <= *varssize);
   assert(varhash != NULL);
   assert(var != NULL);
   assert(varidx != NULL);
   
   /* check if we saw this variable before */
   if( SCIPhashmapExists(varhash, (void*)var) )
   {
      *varidx = (int)(size_t)SCIPhashmapGetImage(varhash, (void*)var);
      assert(*varidx >= 0);
      assert(*varidx < *nvars);
      
      return SCIP_OKAY;
   }
   
   /* since variable is new, add it to the end of vars array and into hashmap */
   
   /* ensure enough space in vars array */
   if( *nvars + 1 > *varssize )
   {
      *varssize = SCIPcalcMemGrowSize(scip, *nvars + 1);
      if( vars == NULL )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, vars, *varssize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, *varssize) );
      }
   }
   assert(*vars != NULL);  /*lint !e613*/
   
   (*vars)[*nvars] = var;  /*lint !e613*/
   SCIP_CALL( SCIPhashmapInsert(varhash, (void*)var, (void*)(size_t)*nvars) );
   *varidx = *nvars;
   
   ++*nvars;
   
   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   PIP_MAX_LINELEN */
   SCIP_EXPRTREE**       exprtree,           /**< pointer to store constraint function as polynomial expression */
   int*                  degree,             /**< pointer to store degree of polynomial */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_EXPR* expression;
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;
   int nextcoefsign;
   int monomialdegree;
   SCIP_EXPR** varexprs;
   int i;
   
   SCIP_VAR** vars;
   int varssize;
   int nvars;
   SCIP_HASHMAP* varhash;
   
   SCIP_EXPRDATA_MONOMIAL** monomials;
   int monomialssize;
   int nmonomials;
   
   int nfactors;
   int factorssize;
   SCIP_Real* exponents;
   int* varidxs;
   
   assert(scip != NULL);
   assert(pipinput != NULL);
   assert(name != NULL);
   assert(exprtree != NULL);
   assert(degree != NULL);
   assert(newsection != NULL);

   *name = '\0';
   *exprtree = NULL;
   *degree = 0;
   *newsection = FALSE;

   /* read the first token, which may be the name of the line */
   if( getNextToken(pipinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(pipinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* remember the token in the token buffer */
      swapTokenBuffer(pipinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(pipinput) )
      {
         if( strcmp(pipinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the line name */
            (void)strncpy(name, pipinput->tokenbuf, PIP_MAX_LINELEN);
            name[PIP_MAX_LINELEN - 1] = '\0';
            SCIPdebugMessage("(line %d) read constraint name: '%s'\n", pipinput->linenumber, name);
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            pushToken(pipinput);
            pushBufferToken(pipinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(pipinput);
      }
   }

   /* initialize buffer for storing the variables */
   varssize = PIP_INIT_VARSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), SCIPcalcHashtableSize(PIP_INIT_VARSSIZE)) );
   
   /* initialize buffer for storing the monomials */
   monomialssize = PIP_INIT_MONOMIALSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &monomials, monomialssize) );
   
   /* initialize buffer for storing the factors in a monomial */
   factorssize = PIP_INIT_FACTORSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, factorssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varidxs,   factorssize) );

   /* read the coefficients */
   coefsign = +1;
   nextcoefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   nmonomials = 0;
   nvars = 0;
   nfactors = 0;
   monomialdegree = 0;
   while( getNextToken(pipinput) )
   {
      SCIP_VAR* var;
      int varidx;
      SCIP_Bool issense;
      SCIP_Bool issign;
      SCIP_Bool isnewsection;
      SCIP_Real exponent;
      
      issign = FALSE;   /* fix compiler warning */
      issense = FALSE;  /* fix lint warning */
      if( (isnewsection = isNewSection(pipinput)) ||  /*lint !e820*/ 
          (issense = isSense(pipinput, NULL))     ||  /*lint !e820*/
          (nfactors > 0 && (issign = isSign(pipinput, &nextcoefsign))) )  /*lint !e820*/
      {
         /* finish the current monomial */
         SCIP_CALL( ensureMonomialsSize(scip, &monomials, &monomialssize, nmonomials + 1) );
         SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip),
            &monomials[nmonomials], coefsign * coef, nfactors, varidxs, exponents) );
         ++nmonomials;
         
         if( monomialdegree > *degree )
            *degree = monomialdegree;

         /* reset variables */
         nfactors = 0;
         coef = 1.0;
         coefsign = +1;
         havesign = FALSE;
         havevalue = FALSE;
         monomialdegree = 0;
         
         if( isnewsection )
         {
            *newsection = TRUE;
            break;
         }
         
         if( issense )
         {
            /* put the sense back onto the token stack */
            pushToken(pipinput);
            break;
         }
         
         if( issign )
         {
            coefsign = nextcoefsign;
            SCIPdebugMessage("(line %d) read coefficient sign: %+d\n", pipinput->linenumber, coefsign);
            havesign = TRUE;
            nextcoefsign = +1;
            continue;
         }
      }
      
      /* check if we read a sign */
      if( isSign(pipinput, &coefsign) )
      {
         SCIPdebugMessage("(line %d) read coefficient sign: %+d\n", pipinput->linenumber, coefsign);

         if( nfactors > 0 || havevalue )
         {
            syntaxError(scip, pipinput, "sign can only be at beginning of monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }
         
         havesign = TRUE;
         continue;
      }
      
      /* check if we are in between factors of a monomial */
      if( strcmp(pipinput->token, "*") == 0 )
      {
         if( nfactors == 0 )
         {
            syntaxError(scip, pipinput, "cannot have '*' before first variable in monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }
         
         continue;
      }
      
      /* all but the first monomial need a sign */
      if( nmonomials > 0 && !havesign )
      {
         syntaxError(scip, pipinput, "expected sign ('+' or '-') or sense ('<' or '>')");
         goto TERMINATE_READPOLYNOMIAL;
      }

      /* check if we are at an exponent for the last variable */
      if( strcmp(pipinput->token, "^") == 0 )
      {
         if( !getNextToken(pipinput) || !isValue(scip, pipinput, &exponent) )
         {
            syntaxError(scip, pipinput, "expected exponent value after '^'");
            goto TERMINATE_READPOLYNOMIAL;
         }
         if( nfactors == 0 )
         {
            syntaxError(scip, pipinput, "cannot have '^' before first variable in monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }
         exponents[nfactors-1] = exponent;
         if( SCIPisIntegral(scip, exponent) && exponent > 0.0 )
            monomialdegree += (int)exponent - 1; /* -1, because we added +1 when we put the variable into varidxs */
         else
            monomialdegree = SCIP_EXPR_DEGREEINFINITY;

         SCIPdebugMessage("(line %d) read exponent value %g for variable %s\n", pipinput->linenumber, exponent,
            SCIPvarGetName(vars[varidxs[nfactors-1]]));
         continue;
      }
      
      /* check if we read a value */
      if( isValue(scip, pipinput, &coef) )
      {
         SCIPdebugMessage("(line %d) read coefficient value: %g with sign %+d\n", pipinput->linenumber, coef, coefsign);
         
         if( havevalue )
         {
            syntaxError(scip, pipinput, "two consecutive values");
            goto TERMINATE_READPOLYNOMIAL;
         }
         
         if( nfactors > 0 )
         {
            syntaxError(scip, pipinput, "coefficients can only be at the beginning of a monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }
         
         havevalue = TRUE;
         continue;
      }

      /* the token is a variable name: get the corresponding variable (or create a new one) */
      SCIP_CALL( getVariable(scip, pipinput->token, &var, NULL) );
      
      /* get the index of the variable in the vars array, or add there if not in it yet */
      SCIP_CALL( getVariableIndex(scip, &vars, &varssize, &nvars, varhash, var, &varidx) );
      
      SCIP_CALL( ensureFactorsSize(scip, &exponents, &varidxs, &factorssize, nfactors + 1) );
      
      exponents[nfactors] = 1.0;
      varidxs[nfactors]   = varidx;
      ++nfactors;
      ++monomialdegree;
   }
   
   if( nfactors > 0 )
   {
      syntaxError(scip, pipinput, "string ended before monomial has finished");
      goto TERMINATE_READPOLYNOMIAL;
   }
   
   /* create variable expressions */
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[i], SCIP_EXPR_VARIDX, i) );
   }
   
   /* create polynomial expression, let polynomial take over ownership of monomials */
   SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &expression, nvars, varexprs,
      nmonomials, monomials, 0.0, FALSE) );
   
   SCIPfreeBufferArray(scip, &varexprs);
   
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), exprtree, expression, 0, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(*exprtree, nvars, vars) );
   
   SCIPdebugMessage("read polynomial of degree %d: ", *degree);
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(*exprtree, NULL) ) );
   SCIPdebugPrintf("\n");
   
TERMINATE_READPOLYNOMIAL:
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &monomials);
   SCIPfreeBufferArray(scip, &exponents);
   SCIPfreeBufferArray(scip, &varidxs);
   SCIPhashmapFree(&varhash);

   return SCIP_OKAY;
}

/** given an expression tree that holds a polynomial expression of degree at most two,
 * gives the coefficients of the constant, linear, and quadratic part of this expression
 */
static
void getLinearAndQuadraticCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree holding polynomial expression */
   SCIP_Real*            constant,           /**< buffer to store constant monomials */
   int*                  nlinvars,           /**< buffer to store number of linear coefficients */
   SCIP_VAR**            linvars,            /**< array to fill with linear variables */
   SCIP_Real*            lincoefs,           /**< array to fill with coefficients of linear variables */
   int*                  nquadterms,         /**< buffer to store number of quadratic terms */ 
   SCIP_VAR**            quadvars1,          /**< array to fill with quadratic variables */
   SCIP_VAR**            quadvars2,          /**< array to fill with quadratic variables */
   SCIP_Real*            quadcoefs           /**< array to fill with coefficients of quadratic terms */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPRDATA_MONOMIAL** monomials;
   int nmonomials;
   int varidx;
   int i;

   expr = SCIPexprtreeGetRoot(exprtree);
   assert(expr != NULL);
   assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);
   assert(SCIPexprGetNChildren(expr) == SCIPexprtreeGetNVars(exprtree));

   nmonomials = SCIPexprGetNMonomials(expr);
   monomials  = SCIPexprGetMonomials(expr);

   *constant = 0.0;
   *nlinvars = 0;
   *nquadterms = 0;
   for( i = 0; i < nmonomials; ++i )
   {
      assert(SCIPexprGetMonomialNFactors(monomials[i]) >= 0);
      assert(SCIPexprGetMonomialNFactors(monomials[i]) <= 2);
      assert(SCIPexprGetMonomialExponents(monomials[i]) != NULL    || SCIPexprGetMonomialNFactors(monomials[i]) == 0);
      assert(SCIPexprGetMonomialChildIndices(monomials[i]) != NULL || SCIPexprGetMonomialNFactors(monomials[i]) == 0);
      
      if( SCIPexprGetMonomialNFactors(monomials[i]) == 0 )
      {
         /* constant monomial */
         *constant += SCIPexprGetMonomialCoef(monomials[i]);
      }
      else if( SCIPexprGetMonomialNFactors(monomials[i]) == 1 && SCIPexprGetMonomialExponents(monomials[i])[0] == 1.0 )
      {
         /* linear monomial */
         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */

         lincoefs[*nlinvars] = SCIPexprGetMonomialCoef(monomials[i]);
         linvars[*nlinvars]  = SCIPexprtreeGetVars(exprtree)[varidx];
         ++*nlinvars;
      }
      else if( SCIPexprGetMonomialNFactors(monomials[i]) == 1 )
      {
         /* square monomial */
         assert(SCIPexprGetMonomialExponents(monomials[i])[0] == 2.0);

         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */

         quadcoefs[*nquadterms] = SCIPexprGetMonomialCoef(monomials[i]);
         quadvars1[*nquadterms] = SCIPexprtreeGetVars(exprtree)[varidx];
         quadvars2[*nquadterms] = quadvars1[*nquadterms];
         ++*nquadterms;
      }
      else
      {
         /* bilinear monomial */
         assert(SCIPexprGetMonomialExponents(monomials[i])[0] == 1.0);
         assert(SCIPexprGetMonomialExponents(monomials[i])[1] == 1.0);

         quadcoefs[*nquadterms] = SCIPexprGetMonomialCoef(monomials[i]);

         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */
         quadvars1[*nquadterms] = SCIPexprtreeGetVars(exprtree)[varidx];
         
         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[1];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */
         quadvars2[*nquadterms] = SCIPexprtreeGetVars(exprtree)[varidx];
         
         ++*nquadterms;
      }
   }
}

/** reads the objective section */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   char name[PIP_MAX_LINELEN];
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   int degree;
   SCIP_Bool newsection;
   int varidx;
   int nmonomials;

   assert(pipinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readPolynomial(scip, pipinput, name, &exprtree, &degree, &newsection) );
   if( !hasError(pipinput) )
   {
      int i;

      expr = SCIPexprtreeGetRoot(exprtree);
      assert(expr != NULL);
      assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);

      nmonomials = SCIPexprGetNMonomials(expr);

      assert(degree >= 0);
      if( degree == 1 )
      {
         SCIP_Real coef;
         SCIP_VAR* var;
         SCIP_EXPRDATA_MONOMIAL** monomials;
         
         assert(SCIPexprtreeGetVars(exprtree) != NULL);
         assert(SCIPexprGetNChildren(expr) == SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetPolynomialConstant(expr) == 0.0);
         
         monomials  = SCIPexprGetMonomials(expr);
         
         for( i = 0; i < nmonomials; ++i )
         {
            assert(SCIPexprGetMonomialNFactors(monomials[i]) == 1);
            assert(SCIPexprGetMonomialExponents(monomials[i]) != NULL);
            assert(SCIPexprGetMonomialExponents(monomials[i])[0] == 1.0);
            assert(SCIPexprGetMonomialChildIndices(monomials[i]) != NULL);
            
            varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
            assert(varidx >= 0);
            assert(varidx < SCIPexprGetNChildren(expr));
            assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
            assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */
            
            coef = SCIPexprGetMonomialCoef(monomials[i]);
            var = SCIPexprtreeGetVars(exprtree)[varidx];
            
            SCIP_CALL( SCIPchgVarObj(scip, var, SCIPvarGetObj(var) + coef) );
         }
      }
      else if( degree == 2 )
      {
         /* insert dummy variable and constraint to represent quadratic part of objective */
         
         SCIP_VAR*  quadobjvar;
         SCIP_CONS* quadobjcons;
         SCIP_Real  lhs;
         SCIP_Real  rhs;
         
         SCIP_Real constant;
         int nlinvars;
         SCIP_VAR** linvars;
         SCIP_Real* lincoefs;
         int nquadterms;
         SCIP_VAR** quadvars1;
         SCIP_VAR** quadvars2;
         SCIP_Real* quadcoefs;
         
         SCIP_CALL( SCIPallocBufferArray(scip, &linvars,   nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs,  nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, nmonomials) );

         getLinearAndQuadraticCoefs(scip, exprtree, &constant, &nlinvars, linvars, lincoefs, &nquadterms, quadvars1, quadvars2, quadcoefs);

         SCIP_CALL( SCIPcreateVar(scip, &quadobjvar, "quadobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, 
               SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, quadobjvar) );

         if ( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
         {
            lhs = -SCIPinfinity(scip);
            rhs = -constant;
         }
         else
         {
            lhs = -constant;
            rhs = SCIPinfinity(scip);
         }

         SCIP_CALL( SCIPcreateConsQuadratic(scip, &quadobjcons, "quadobj", nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, quadobjcons, quadobjvar, -1.0) );

         SCIP_CALL( SCIPaddCons(scip, quadobjcons) );
         SCIPdebugMessage("(line %d) added constraint <%s> to represent quadratic objective: ", pipinput->linenumber, SCIPconsGetName(quadobjcons));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, quadobjcons, NULL) ) );

         SCIP_CALL( SCIPreleaseCons(scip, &quadobjcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &quadobjvar) );
         
         /* free memory */
         SCIPfreeBufferArray(scip, &linvars);
         SCIPfreeBufferArray(scip, &lincoefs);
         SCIPfreeBufferArray(scip, &quadvars1);
         SCIPfreeBufferArray(scip, &quadvars2);
         SCIPfreeBufferArray(scip, &quadcoefs);
      }
      else if( degree > 2 )
      {
         SCIPerrorMessage("polynomial constraints of degree > 2 not yet supported by SCIP\n");
         SCIP_CALL( SCIPexprtreeFree(&exprtree) );
         return SCIP_ERROR;
      }
   }
   
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** reads the constraints section 
 */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   char name[PIP_MAX_LINELEN];
   SCIP_CONS* cons;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   int degree;

   SCIP_Real constant;
   
   int nlinvars;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   
   int nquadcoefs;
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
   
   SCIP_Bool newsection;
   PIPSENSE sense;
   SCIP_Real sidevalue;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamicrows;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   int sidesign;
   int nmonomials;

   assert(pipinput != NULL);

   /* read polynomial */
   SCIP_CALL( readPolynomial(scip, pipinput, name, &exprtree, &degree, &newsection) );
   if ( hasError(pipinput) )
      goto TERMINATE;
   if ( newsection )
   {
      if ( exprtree != NULL )
         syntaxError(scip, pipinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if ( !getNextToken(pipinput) || !isSense(pipinput, &sense) )
   {
      syntaxError(scip, pipinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if ( !getNextToken(pipinput) )
   {
      syntaxError(scip, pipinput, "missing right hand side");
      goto TERMINATE;
   }
   if ( isSign(pipinput, &sidesign) )
   {
      if( !getNextToken(pipinput) )
      {
         syntaxError(scip, pipinput, "missing value of right hand side");
         goto TERMINATE;
      }
   }
   if ( !isValue(scip, pipinput, &sidevalue) )
   {
      syntaxError(scip, pipinput, "expected value as right hand side");
      goto TERMINATE;
   }
   sidevalue *= sidesign;

   /* create and add the linear constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/pipreader/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/pipreader/dynamicrows", &dynamicrows) );
   initial = !dynamicrows;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows;
   
   if( degree > 2 )
   {
      SCIPerrorMessage("polynomial constraints of degree > 2 not supported by SCIP yet\n");
      return SCIP_ERROR;
   }
   else
   {
      expr = SCIPexprtreeGetRoot(exprtree);
      assert(expr != NULL);
      assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);
      nmonomials = SCIPexprGetNMonomials(expr);
      
      SCIP_CALL( SCIPallocBufferArray(scip, &linvars,   nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs,  nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, nmonomials) );

      getLinearAndQuadraticCoefs(scip, exprtree, &constant, &nlinvars, linvars, lincoefs, &nquadcoefs, quadvars1, quadvars2, quadcoefs);

      /* assign the left and right hand side, depending on the constraint sense */
      switch ( sense )
      {
      case PIP_SENSE_GE:
         lhs = sidevalue - constant;
         rhs = SCIPinfinity(scip);
         break;
      case PIP_SENSE_LE:
         lhs = -SCIPinfinity(scip);
         rhs = sidevalue - constant;
         break;
      case PIP_SENSE_EQ:
         lhs = sidevalue - constant;
         rhs = sidevalue - constant;
         break;
      case PIP_SENSE_NOTHING:
      default:
         SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
         return SCIP_INVALIDDATA;
      }
      
      if( nquadcoefs == 0 )
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nlinvars, linvars, lincoefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsQuadratic(scip, &cons, name, nlinvars, linvars, lincoefs,
            nquadcoefs, quadvars1, quadvars2, quadcoefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
      }
      
      /* free memory */
      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &lincoefs);
      SCIPfreeBufferArray(scip, &quadvars1);
      SCIPfreeBufferArray(scip, &quadvars2);
      SCIPfreeBufferArray(scip, &quadcoefs);
      
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("(line %d) created constraint: ", pipinput->linenumber);
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   
TERMINATE:
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** reads the bounds section */
static
SCIP_RETCODE readBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   while( getNextToken(pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Real lb;
      SCIP_Real ub;
      int sign;
      SCIP_Bool hassign;
      PIPSENSE leftsense;

      /* check if we reached a new section */
      if( isNewSection(pipinput) )
         return SCIP_OKAY;

      /* default bounds are [0,+inf] */
      lb = 0.0;
      ub = SCIPinfinity(scip);
      leftsense = PIP_SENSE_NOTHING;

      /* check if the first token is a sign */
      sign = +1;
      hassign = isSign(pipinput, &sign);
      if( hassign && !getNextToken(pipinput) )
      {
         syntaxError(scip, pipinput, "expected value");
         return SCIP_OKAY;
      }

      /* the first token must be either a value or a variable name */
      if( isValue(scip, pipinput, &value) )
      {
         /* first token is a value: the second token must be a sense */
         if( !getNextToken(pipinput) || !isSense(pipinput, &leftsense) )
         {
            syntaxError(scip, pipinput, "expected bound sense '<=', '=', or '>='");
            return SCIP_OKAY;
         }

         /* update the bound corresponding to the sense */
         switch( leftsense )
         {
         case PIP_SENSE_GE:
            ub = sign * value;
            break;
         case PIP_SENSE_LE:
            lb = sign * value;
            break;
         case PIP_SENSE_EQ:
            lb = sign * value;
            ub = sign * value;
            break;
         case PIP_SENSE_NOTHING:
         default:
            SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
            return SCIP_INVALIDDATA;
         }
      }
      else if( hassign )
      {
         syntaxError(scip, pipinput, "expected value");
         return SCIP_OKAY;
      }
      else
         pushToken(pipinput);

      /* the next token must be a variable name */
      if( !getNextToken(pipinput) )
      {
         syntaxError(scip, pipinput, "expected variable name");
         return SCIP_OKAY;
      }
      SCIP_CALL( getVariable(scip, pipinput->token, &var, NULL) );

      /* the next token might be another sense, or the word "free" */
      if( getNextToken(pipinput) )
      {
         PIPSENSE rightsense;

         if( isSense(pipinput, &rightsense) )
         {
            /* check, if the senses fit */
            if( leftsense == PIP_SENSE_NOTHING
               || (leftsense == PIP_SENSE_LE && rightsense == PIP_SENSE_LE)
               || (leftsense == PIP_SENSE_GE && rightsense == PIP_SENSE_GE) )
            {
               if( !getNextToken(pipinput) )
               {
                  syntaxError(scip, pipinput, "expected value or sign");
                  return SCIP_OKAY;
               }

               /* check if the next token is a sign */
               sign = +1;
               hassign = isSign(pipinput, &sign);
               if( hassign && !getNextToken(pipinput) )
               {
                  syntaxError(scip, pipinput, "expected value");
                  return SCIP_OKAY;
               }

               /* the next token must be a value */
               if( !isValue(scip, pipinput, &value) )
               {
                  syntaxError(scip, pipinput, "expected value");
                  return SCIP_OKAY;
               }

               /* update the bound corresponding to the sense */
               switch( rightsense )
               {
               case PIP_SENSE_GE:
                  lb = sign * value;
                  break;
               case PIP_SENSE_LE:
                  ub = sign * value;
                  break;
               case PIP_SENSE_EQ:
                  lb = sign * value;
                  ub = sign * value;
                  break;
               case PIP_SENSE_NOTHING:
               default:
                  SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
                  return SCIP_INVALIDDATA;
               }
            }
            else
            {
               syntaxError(scip, pipinput, "the two bound senses do not fit");
               return SCIP_OKAY;
            }
         }
         else if( strcasecmp(pipinput->token, "FREE") == 0 )
         {
            if( leftsense != PIP_SENSE_NOTHING )
            {
               syntaxError(scip, pipinput, "variable with bound is marked as 'free'");
               return SCIP_OKAY;
            }
            lb = -SCIPinfinity(scip);
            ub = SCIPinfinity(scip);
         }
         else
         {
            /* the token was no sense: push it back to the token stack */
            pushToken(pipinput);
         }
      }

      /* change the bounds of the variable if bounds have been given (do not destroy earlier specification of bounds) */
      if ( lb != 0.0 )
	 SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      /*lint --e{777}*/
      if ( ub != SCIPinfinity(scip) )
	 SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      SCIPdebugMessage("(line %d) new bounds: <%s>[%g,%g]\n", pipinput->linenumber, SCIPvarGetName(var),
	 SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

/** reads the generals section */
static
SCIP_RETCODE readGenerals(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   while( getNextToken(pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;
      SCIP_Bool infeasible;

      /* check if we reached a new section */
      if( isNewSection(pipinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, pipinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, pipinput, "unknown variable in generals section");
         return SCIP_OKAY;
      }

      /* mark the variable to be integral */
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible) );
      /* don't assert feasibility here because the presolver will and should detect a infeasibility */
   }

   return SCIP_OKAY;
}

/** reads the binaries section */
static
SCIP_RETCODE readBinaries(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   while( getNextToken(pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;
      SCIP_Bool infeasible;

      /* check if we reached a new section */
      if( isNewSection(pipinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, pipinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, pipinput, "unknown variable in binaries section");
         return SCIP_OKAY;
      }

      /* mark the variable to be binary and change its bounds appropriately */
      if( SCIPvarGetLbGlobal(var) < 0.0 )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );
      }
      if( SCIPvarGetUbGlobal(var) > 1.0 )
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 1.0) );
      }
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
      /* don't assert feasibility here because the presolver will and should detect a infeasibility */
   }

   return SCIP_OKAY;
}

/** reads a PIP file
 */
static
SCIP_RETCODE readPIPFile(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   const char*           filename            /**< name of the input file */
   )
{
   assert(pipinput != NULL);

   /* open file */
   pipinput->file = SCIPfopen(filename, "r");
   if( pipinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* parse the file */
   pipinput->section = PIP_START;
   while( pipinput->section != PIP_END && !hasError(pipinput) )
   {
      switch( pipinput->section )
      {
      case PIP_START:
         SCIP_CALL( readStart(scip, pipinput) );
         break;

      case PIP_OBJECTIVE:
         SCIP_CALL( readObjective(scip, pipinput) );
         break;

      case PIP_CONSTRAINTS:
         SCIP_CALL( readConstraints(scip, pipinput) );
         break;

      case PIP_BOUNDS:
         SCIP_CALL( readBounds(scip, pipinput) );
         break;

      case PIP_GENERALS:
         SCIP_CALL( readGenerals(scip, pipinput) );
         break;

      case PIP_BINARIES:
         SCIP_CALL( readBinaries(scip, pipinput) );
         break;

      case PIP_END: /* this is already handled in the while() loop */
      default:
         SCIPerrorMessage("invalid PIP file section <%d>\n", pipinput->section);
         return SCIP_INVALIDDATA;
      }
   }

   /* close file */
   SCIPfclose(pipinput->file);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyPip)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderPip(scip) );
 
   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreePip NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadPip)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadPip(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWritePip NULL


/*
 * reader specific interface methods
 */

/** includes the pip file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderPip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create lp reader data */
   readerdata = NULL;

   /* include lp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyPip,
         readerFreePip, readerReadPip, readerWritePip, 
         readerdata) );

   /* add lp reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dynamiccols", "should columns be added and removed dynamically to the PIP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dynamicrows", "should rows be added and removed dynamically to the PIP?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}


/* reads problem from file */
SCIP_RETCODE SCIPreadPip(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   PIPINPUT pipinput;
   int i;

   /* initialize PIP input data */
   pipinput.file = NULL;
   pipinput.linebuf[0] = '\0';
   pipinput.probname[0] = '\0';
   pipinput.objname[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &pipinput.token, PIP_MAX_LINELEN) );
   pipinput.token[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &pipinput.tokenbuf, PIP_MAX_LINELEN) );
   pipinput.tokenbuf[0] = '\0';
   for( i = 0; i < PIP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((pipinput.pushedtokens)[i]), PIP_MAX_LINELEN) );
   }

   pipinput.npushedtokens = 0;
   pipinput.linenumber = 0;
   pipinput.linepos = 0;
   pipinput.section = PIP_START;
   pipinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   pipinput.haserror = FALSE;

   /* read the file */
   SCIP_CALL( readPIPFile(scip, &pipinput, filename) );

   /* free dynamically allocated memory */
   SCIPfreeMemoryArray(scip, &pipinput.token);
   SCIPfreeMemoryArray(scip, &pipinput.tokenbuf);
   for( i = 0; i < PIP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeMemoryArray(scip, &pipinput.pushedtokens[i]);
   }

   /* evaluate the result */
   if( pipinput.haserror )
      return SCIP_READERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, pipinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}

