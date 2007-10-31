/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_lp.c,v 1.30 2007/10/31 09:26:31 bzfheinz Exp $"

/**@file   reader_lp.c
 * @brief  LP file reader
 * @author Tobias Achterberg
 * @author Marc Pfetsch
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

#include "scip/reader_lp.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"

#define READER_NAME             "lpreader"
#define READER_DESC             "file reader for MIPs in ILOG's LP file format"
#define READER_EXTENSION        "lp"


/*
 * Data structures
 */

#define LP_MAX_LINELEN       65536
#define LP_MAX_PUSHEDTOKENS  2
#define LP_INIT_COEFSSIZE    8192

/** Section in LP File */
enum LpSection
{
   LP_START, LP_OBJECTIVE, LP_CONSTRAINTS, LP_BOUNDS, LP_GENERALS, LP_BINARIES, LP_SEMICONTINUOUS, LP_SOS, LP_END
};
typedef enum LpSection LPSECTION;

enum LpExpType
{
   LP_EXP_NONE, LP_EXP_UNSIGNED, LP_EXP_SIGNED
};
typedef enum LpExpType LPEXPTYPE;

enum LpSense
{
   LP_SENSE_NOTHING, LP_SENSE_LE, LP_SENSE_GE, LP_SENSE_EQ
};
typedef enum LpSense LPSENSE;

/** LP reading data */
struct LpInput
{
   SCIP_FILE*           file;
   char                 linebuf[LP_MAX_LINELEN];
   char                 probname[LP_MAX_LINELEN];
   char                 objname[LP_MAX_LINELEN];
   char*                token;
   char*                tokenbuf;
   char*                pushedtokens[LP_MAX_PUSHEDTOKENS];
   int                  npushedtokens;
   int                  linenumber;
   int                  linepos;
   LPSECTION            section;
   SCIP_OBJSENSE        objsense;
   SCIP_Bool            inlazyconstraints;
   SCIP_Bool            inusercuts;
   SCIP_Bool            haserror;
};
typedef struct LpInput LPINPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=";
static const char commentchars[] = "\\";




/*
 * Local methods (for reading)
 */

/** issues an error message and marks the LP data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   const char*           msg                 /**< error message */
   )
{
   char formatstr[255];

   assert(lpinput != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error in line %d: %s ('%s')\n",
      lpinput->linenumber, msg, lpinput->token);
   if( lpinput->linebuf[strlen(lpinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", lpinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", lpinput->linebuf);
   }
   sprintf(formatstr, "         %%%ds\n", lpinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, formatstr, "^");
   lpinput->section  = LP_END;
   lpinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   return lpinput->haserror;
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
   LPEXPTYPE*            exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit(c) )
      return TRUE;
   else if( (*exptype == LP_EXP_NONE) && !(*hasdot) && (c == '.') )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == LP_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = LP_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit(nextc) )
      {
         *exptype = LP_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == LP_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = LP_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   int i;

   assert(lpinput != NULL);

   /* clear the line */
   BMSclearMemoryArray(lpinput->linebuf, LP_MAX_LINELEN);

   /* read next line */
   lpinput->linepos = 0;
   lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';
   if( SCIPfgets(lpinput->linebuf, sizeof(lpinput->linebuf), lpinput->file) == NULL )
      return FALSE;
   lpinput->linenumber++;
   if( lpinput->linebuf[LP_MAX_LINELEN-2] != '\0' )
   {
      SCIPerrorMessage("Error: line %d exceeds %d characters\n", lpinput->linenumber, LP_MAX_LINELEN-2);
      lpinput->haserror = TRUE;
      return FALSE;
   }
   lpinput->linebuf[LP_MAX_LINELEN-1] = '\0';
   lpinput->linebuf[LP_MAX_LINELEN-2] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(lpinput->linebuf, commentchars[i]);
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
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Bool hasdot;
   LPEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(lpinput != NULL);
   assert(lpinput->linepos < LP_MAX_LINELEN);

   /* check the token stack */
   if( lpinput->npushedtokens > 0 )
   {
      swapPointers(&lpinput->token, &lpinput->pushedtokens[lpinput->npushedtokens-1]);
      lpinput->npushedtokens--;
      SCIPdebugMessage("(line %d) read token again: '%s'\n", lpinput->linenumber, lpinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = lpinput->linebuf;
   while( isDelimChar(buf[lpinput->linepos]) )
   {
      if( buf[lpinput->linepos] == '\0' )
      {
         if( !getNextLine(lpinput) )
         {
            lpinput->section = LP_END;
            SCIPdebugMessage("(line %d) end of file\n", lpinput->linenumber);
            return FALSE;
         }
         assert(lpinput->linepos == 0);
      }
      else
         lpinput->linepos++;
   }
   assert(lpinput->linepos < LP_MAX_LINELEN);
   assert(!isDelimChar(buf[lpinput->linepos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = LP_EXP_NONE;
   if( isValueChar(buf[lpinput->linepos], buf[lpinput->linepos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < LP_MAX_LINELEN);
         assert(!isDelimChar(buf[lpinput->linepos]));
         lpinput->token[tokenlen] = buf[lpinput->linepos];
         tokenlen++;
         lpinput->linepos++;
      }
      while( isValueChar(buf[lpinput->linepos], buf[lpinput->linepos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < LP_MAX_LINELEN);
         lpinput->token[tokenlen] = buf[lpinput->linepos];
         tokenlen++;
         lpinput->linepos++;
         if( tokenlen == 1 && isTokenChar(lpinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[lpinput->linepos]) && !isTokenChar(buf[lpinput->linepos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1
         && (lpinput->token[tokenlen-1] == '<' || lpinput->token[tokenlen-1] == '>' || lpinput->token[tokenlen-1] == '=')
         && buf[lpinput->linepos] == '=' )
      {
         lpinput->linepos++;
      }
      else if( lpinput->token[tokenlen-1] == '=' && (buf[lpinput->linepos] == '<' || buf[lpinput->linepos] == '>') )
      {
         lpinput->token[tokenlen-1] = buf[lpinput->linepos];
         lpinput->linepos++;
      }
   }
   assert(tokenlen < LP_MAX_LINELEN);
   lpinput->token[tokenlen] = '\0';

   SCIPdebugMessage("(line %d) read token: '%s'\n", lpinput->linenumber, lpinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);
   assert(lpinput->npushedtokens < LP_MAX_PUSHEDTOKENS);

   swapPointers(&lpinput->pushedtokens[lpinput->npushedtokens], &lpinput->token);
   lpinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);
   assert(lpinput->npushedtokens < LP_MAX_PUSHEDTOKENS);

   swapPointers(&lpinput->pushedtokens[lpinput->npushedtokens], &lpinput->tokenbuf);
   lpinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   swapPointers(&lpinput->token, &lpinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Bool iscolon;

   assert(lpinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(lpinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(lpinput) )
   {
      iscolon = (strcmp(lpinput->token, ":") == 0);
      pushToken(lpinput);
   }

   /* reinstall the previous token by swapping back the token buffer */
   swapTokenBuffer(lpinput);

   /* check for ':' */
   if( iscolon )
      return FALSE;

   if( strcasecmp(lpinput->token, "MINIMIZE") == 0
      || strcasecmp(lpinput->token, "MINIMUM") == 0
      || strcasecmp(lpinput->token, "MIN") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: OBJECTIVE\n", lpinput->linenumber);
      lpinput->section = LP_OBJECTIVE;
      lpinput->objsense = SCIP_OBJSENSE_MINIMIZE;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "MAXIMIZE") == 0
      || strcasecmp(lpinput->token, "MAXIMUM") == 0
      || strcasecmp(lpinput->token, "MAX") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: OBJECTIVE\n", lpinput->linenumber);
      lpinput->section = LP_OBJECTIVE;
      lpinput->objsense = SCIP_OBJSENSE_MAXIMIZE;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "SUBJECT") == 0 )
   {
      /* check if the next token is 'TO' */
      swapTokenBuffer(lpinput);
      if( getNextToken(lpinput) )
      {
         if( strcasecmp(lpinput->token, "TO") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", lpinput->linenumber);
            lpinput->section = LP_CONSTRAINTS;
            lpinput->inlazyconstraints = FALSE;
            lpinput->inusercuts = FALSE;
            return TRUE;
         }
         else
            pushToken(lpinput);
      }
      swapTokenBuffer(lpinput);
   }

   if( strcasecmp(lpinput->token, "SUCH") == 0 )
   {
      /* check if the next token is 'THAT' */
      swapTokenBuffer(lpinput);
      if( getNextToken(lpinput) )
      {
         if( strcasecmp(lpinput->token, "THAT") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", lpinput->linenumber);
            lpinput->section = LP_CONSTRAINTS;
            lpinput->inlazyconstraints = FALSE;
            lpinput->inusercuts = FALSE;
            return TRUE;
         }
         else
            pushToken(lpinput);
      }
      swapTokenBuffer(lpinput);
   }

   if( strcasecmp(lpinput->token, "st") == 0
      || strcasecmp(lpinput->token, "S.T.") == 0
      || strcasecmp(lpinput->token, "ST.") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", lpinput->linenumber);
      lpinput->section = LP_CONSTRAINTS;
      lpinput->inlazyconstraints = FALSE;
      lpinput->inusercuts = FALSE;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "LAZY") == 0 )
   {
      /* check if the next token is 'CONSTRAINTS' */
      swapTokenBuffer(lpinput);
      if( getNextToken(lpinput) )
      {
         if( strcasecmp(lpinput->token, "CONSTRAINTS") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS (lazy)\n", lpinput->linenumber);
            lpinput->section = LP_CONSTRAINTS;
            lpinput->inlazyconstraints = TRUE;
            lpinput->inusercuts = FALSE;
            return TRUE;
         }
         else
            pushToken(lpinput);
      }
      swapTokenBuffer(lpinput);
   }

   if( strcasecmp(lpinput->token, "USER") == 0 )
   {
      /* check if the next token is 'CUTS' */
      swapTokenBuffer(lpinput);
      if( getNextToken(lpinput) )
      {
         if( strcasecmp(lpinput->token, "CUTS") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS (user cuts)\n", lpinput->linenumber);
            lpinput->section = LP_CONSTRAINTS;
            lpinput->inlazyconstraints = FALSE;
            lpinput->inusercuts = TRUE;
            return TRUE;
         }
         else
            pushToken(lpinput);
      }
      swapTokenBuffer(lpinput);
   }

   if( strcasecmp(lpinput->token, "BOUNDS") == 0
      || strcasecmp(lpinput->token, "BOUND") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: BOUNDS\n", lpinput->linenumber);
      lpinput->section = LP_BOUNDS;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "GENERAL") == 0
      || strcasecmp(lpinput->token, "GENERALS") == 0
      || strcasecmp(lpinput->token, "GEN") == 0
      || strcasecmp(lpinput->token, "INTEGER") == 0
      || strcasecmp(lpinput->token, "INTEGERS") == 0
      || strcasecmp(lpinput->token, "INT") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: GENERALS\n", lpinput->linenumber);
      lpinput->section = LP_GENERALS;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "BINARY") == 0
      || strcasecmp(lpinput->token, "BINARIES") == 0
      || strcasecmp(lpinput->token, "BIN") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: BINARIES\n", lpinput->linenumber);
      lpinput->section = LP_BINARIES;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "SEMI-CONTINUOUS") == 0
      || strcasecmp(lpinput->token, "SEMIS") == 0
      || strcasecmp(lpinput->token, "SEMI") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: SEMICONTINUOUS\n", lpinput->linenumber);
      lpinput->section = LP_SEMICONTINUOUS;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "SOS") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: SOS\n", lpinput->linenumber);
      lpinput->section = LP_SOS;
      return TRUE;
   }

   if( strcasecmp(lpinput->token, "END") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: END\n", lpinput->linenumber);
      lpinput->section = LP_END;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   LPINPUT*              lpinput,            /**< LP reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(lpinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( lpinput->token[1] == '\0' )
   {
      if( *lpinput->token == '+' )
         return TRUE;
      else if( *lpinput->token == '-' )
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
   LPINPUT*              lpinput,            /**< LP reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(lpinput != NULL);
   assert(value != NULL);

   if( strcasecmp(lpinput->token, "INFINITY") == 0 || strcasecmp(lpinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(lpinput->token, &endptr);
      if( endptr != lpinput->token && *endptr == '\0' )
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
   LPINPUT*              lpinput,            /**< LP reading data */
   LPSENSE*              sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(lpinput != NULL);

   if( strcmp(lpinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = LP_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(lpinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = LP_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(lpinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = LP_SENSE_EQ;
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

      SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamiccols", &dynamiccols) );
      initial = !dynamiccols;
      removable = dynamiccols;

      /* create new variable of the given name */
      SCIPdebugMessage("creating new variable: <%s>\n", name);
      SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
            initial, removable, NULL, NULL, NULL, NULL) );
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
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( !getNextToken(lpinput) )
         return SCIP_OKAY;
   }
   while( !isNewSection(lpinput) );

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   LP_MAX_LINELEN */
   SCIP_VAR***           vars,               /**< pointer to store the array with variables (must be freed by caller) */
   SCIP_Real**           coefs,              /**< pointer to store the array with coefficients (must be freed by caller) */
   int*                  ncoefs,             /**< pointer to store the number of coefficients */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;
   int coefssize;

   assert(lpinput != NULL);
   assert(name != NULL);
   assert(vars != NULL);
   assert(coefs != NULL);
   assert(ncoefs != NULL);
   assert(newsection != NULL);

   *vars = NULL;
   *coefs = NULL;
   *name = '\0';
   *ncoefs = 0;
   *newsection = FALSE;

   /* read the first token, which may be the name of the line */
   if( getNextToken(lpinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(lpinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* remember the token in the token buffer */
      swapTokenBuffer(lpinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(lpinput) )
      {
         if( strcmp(lpinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the line name */
            strncpy(name, lpinput->tokenbuf, SCIP_MAXSTRLEN);
            name[SCIP_MAXSTRLEN-1] = '\0';
            SCIPdebugMessage("(line %d) read constraint name: '%s'\n", lpinput->linenumber, name);
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            pushToken(lpinput);
            pushBufferToken(lpinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(lpinput);
      }
   }

   /* initialize buffers for storing the coefficients */
   coefssize = LP_INIT_COEFSSIZE;
   SCIP_CALL( SCIPallocMemoryArray(scip, vars, coefssize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, coefs, coefssize) );

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   *ncoefs = 0;
   while( getNextToken(lpinput) )
   {
      SCIP_VAR* var;

      /* check if we reached a new section */
      if( isNewSection(lpinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* check if we reached an equation sense */
      if( isSense(lpinput, NULL) )
      {
         /* put the sense back onto the token stack */
         pushToken(lpinput);
         break;
      }

      /* check if we read a sign */
      if( isSign(lpinput, &coefsign) )
      {
         SCIPdebugMessage("(line %d) read coefficient sign: %+d\n", lpinput->linenumber, coefsign);
         havesign = TRUE;
         continue;
      }

      /* all but the first coefficient need a sign */
      if( *ncoefs > 0 && !havesign )
      {
         syntaxError(scip, lpinput, "expected sign ('+' or '-') or sense ('<' or '>')");
         return SCIP_OKAY;
      }

      /* check if we read a value */
      if( isValue(scip, lpinput, &coef) )
      {
         SCIPdebugMessage("(line %d) read coefficient value: %g with sign %+d\n", lpinput->linenumber, coef, coefsign);
         if( havevalue )
         {
            syntaxError(scip, lpinput, "two consecutive values");
            return SCIP_OKAY;
         }
         havevalue = TRUE;
         continue;
      }

      /* the token is a variable name: get the corresponding variable (or create a new one) */
      SCIP_CALL( getVariable(scip, lpinput->token, &var, NULL) );

      /* insert the coefficient */
      SCIPdebugMessage("(line %d) read coefficient: %+g<%s>\n", lpinput->linenumber, coefsign * coef, SCIPvarGetName(var));
      if( !SCIPisZero(scip, coef) )
      {
         /* resize the vars and coefs array if needed */
         if( *ncoefs >= coefssize )
         {
            coefssize *= 2;
            coefssize = MAX(coefssize, (*ncoefs)+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, vars, coefssize) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, coefs, coefssize) );
         }
         assert(*ncoefs < coefssize);

         /* add coefficient */
         (*vars)[*ncoefs] = var;
         (*coefs)[*ncoefs] = coefsign * coef;
         (*ncoefs)++;
      }

      /* reset the flags and coefficient value for the next coefficient */
      coefsign = +1;
      coef = 1.0;
      havesign = FALSE;
      havevalue = FALSE;
   }

   return SCIP_OKAY;
}

/** reads the objective section */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   char name[LP_MAX_LINELEN];
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int ncoefs;
   SCIP_Bool newsection;

   assert(lpinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, lpinput, name, &vars, &coefs, &ncoefs, &newsection) );
   if( !hasError(lpinput) )
   {
      int i;

      /* set the objective values */
      for( i = 0; i < ncoefs; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(scip, vars[i], coefs[i]) );
      }
   }

   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &vars);
   SCIPfreeMemoryArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** reads the constraints section */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   char name[LP_MAX_LINELEN];
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Bool newsection;
   LPSENSE sense;
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
   int ncoefs;
   int sidesign;

   assert(lpinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, lpinput, name, &vars, &coefs, &ncoefs, &newsection) );
   if( hasError(lpinput) )
      goto TERMINATE;
   if( newsection )
   {
      if( ncoefs > 0 )
         syntaxError(scip, lpinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if( !getNextToken(lpinput) || !isSense(lpinput, &sense) )
   {
      syntaxError(scip, lpinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if( !getNextToken(lpinput) )
   {
      syntaxError(scip, lpinput, "missing right hand side");
      goto TERMINATE;
   }
   if( isSign(lpinput, &sidesign) )
   {
      if( !getNextToken(lpinput) )
      {
         syntaxError(scip, lpinput, "missing value of right hand side");
         goto TERMINATE;
      }
   }
   if( !isValue(scip, lpinput, &sidevalue) )
   {
      syntaxError(scip, lpinput, "expected value as right hand side");
      goto TERMINATE;
   }
   sidevalue *= sidesign;

   /* assign the left and right hand side, depending on the constraint sense */
   switch( sense )
   {
   case LP_SENSE_GE:
      lhs = sidevalue;
      rhs = SCIPinfinity(scip);
      break;
   case LP_SENSE_LE:
      lhs = -SCIPinfinity(scip);
      rhs = sidevalue;
      break;
   case LP_SENSE_EQ:
      lhs = sidevalue;
      rhs = sidevalue;
      break;
   case LP_SENSE_NOTHING:
   default:
      SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
      return SCIP_INVALIDDATA;
   }

   /* create and add the linear constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamicrows", &dynamicrows) );
   initial = !dynamicrows && !lpinput->inlazyconstraints && !lpinput->inusercuts;
   separate = TRUE;
   enforce = !lpinput->inusercuts;
   check = !lpinput->inusercuts;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows || lpinput->inusercuts;
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, lhs, rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugMessage("(line %d) created constraint%s: ", lpinput->linenumber,
      lpinput->inlazyconstraints ? " (lazy)" : (lpinput->inusercuts ? " (user cut)" : ""));
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

 TERMINATE:
   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &vars);
   SCIPfreeMemoryArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** reads the bounds section */
static
SCIP_RETCODE readBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   while( getNextToken(lpinput) )
   {
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Real lb;
      SCIP_Real ub;
      int sign;
      SCIP_Bool hassign;
      LPSENSE leftsense;

      /* check if we reached a new section */
      if( isNewSection(lpinput) )
         return SCIP_OKAY;

      /* default bounds are [0,+inf] */
      lb = 0.0;
      ub = SCIPinfinity(scip);
      leftsense = LP_SENSE_NOTHING;

      /* check if the first token is a sign */
      sign = +1;
      hassign = isSign(lpinput, &sign);
      if( hassign && !getNextToken(lpinput) )
      {
         syntaxError(scip, lpinput, "expected value");
         return SCIP_OKAY;
      }

      /* the first token must be either a value or a variable name */
      if( isValue(scip, lpinput, &value) )
      {
         /* first token is a value: the second token must be a sense */
         if( !getNextToken(lpinput) || !isSense(lpinput, &leftsense) )
         {
            syntaxError(scip, lpinput, "expected bound sense '<=', '=', or '>='");
            return SCIP_OKAY;
         }

         /* update the bound corresponding to the sense */
         switch( leftsense )
         {
         case LP_SENSE_GE:
            ub = sign * value;
            break;
         case LP_SENSE_LE:
            lb = sign * value;
            break;
         case LP_SENSE_EQ:
            lb = sign * value;
            ub = sign * value;
            break;
         case LP_SENSE_NOTHING:
         default:
            SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
            return SCIP_INVALIDDATA;
         }
      }
      else if( hassign )
      {
         syntaxError(scip, lpinput, "expected value");
         return SCIP_OKAY;
      }
      else
         pushToken(lpinput);

      /* the next token must be a variable name */
      if( !getNextToken(lpinput) )
      {
         syntaxError(scip, lpinput, "expected variable name");
         return SCIP_OKAY;
      }
      SCIP_CALL( getVariable(scip, lpinput->token, &var, NULL) );

      /* the next token might be another sense, or the word "free" */
      if( getNextToken(lpinput) )
      {
         LPSENSE rightsense;

         if( isSense(lpinput, &rightsense) )
         {
            /* check, if the senses fit */
            if( leftsense == LP_SENSE_NOTHING
               || (leftsense == LP_SENSE_LE && rightsense == LP_SENSE_LE)
               || (leftsense == LP_SENSE_GE && rightsense == LP_SENSE_GE) )
            {
               if( !getNextToken(lpinput) )
               {
                  syntaxError(scip, lpinput, "expected value or sign");
                  return SCIP_OKAY;
               }

               /* check if the next token is a sign */
               sign = +1;
               hassign = isSign(lpinput, &sign);
               if( hassign && !getNextToken(lpinput) )
               {
                  syntaxError(scip, lpinput, "expected value");
                  return SCIP_OKAY;
               }

               /* the next token must be a value */
               if( !isValue(scip, lpinput, &value) )
               {
                  syntaxError(scip, lpinput, "expected value");
                  return SCIP_OKAY;
               }

               /* update the bound corresponding to the sense */
               switch( rightsense )
               {
               case LP_SENSE_GE:
                  lb = sign * value;
                  break;
               case LP_SENSE_LE:
                  ub = sign * value;
                  break;
               case LP_SENSE_EQ:
                  lb = sign * value;
                  ub = sign * value;
                  break;
               case LP_SENSE_NOTHING:
               default:
                  SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
                  return SCIP_INVALIDDATA;
               }
            }
            else
            {
               syntaxError(scip, lpinput, "the two bound senses do not fit");
               return SCIP_OKAY;
            }
         }
         else if( strcasecmp(lpinput->token, "FREE") == 0 )
         {
            if( leftsense != LP_SENSE_NOTHING )
            {
               syntaxError(scip, lpinput, "variable with bound is marked as 'free'");
               return SCIP_OKAY;
            }
            lb = -SCIPinfinity(scip);
            ub = SCIPinfinity(scip);
         }
         else
         {
            /* the token was no sense: push it back to the token stack */
            pushToken(lpinput);
         }
      }

      /* change the bounds of the variable */
      SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      SCIPdebugMessage("(line %d) new bounds: <%s>[%g,%g]\n", lpinput->linenumber, SCIPvarGetName(var), lb, ub);
   }

   return SCIP_OKAY;
}

/** reads the generals section */
static
SCIP_RETCODE readGenerals(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   while( getNextToken(lpinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;

      /* check if we reached a new section */
      if( isNewSection(lpinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpinput, "unknown variable in generals section");
         return SCIP_OKAY;
      }

      /* mark the variable to be integral */
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER) );
   }

   return SCIP_OKAY;
}

/** reads the binaries section */
static
SCIP_RETCODE readBinaries(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   while( getNextToken(lpinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;

      /* check if we reached a new section */
      if( isNewSection(lpinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpinput, "unknown variable in binaries section");
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
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY) );
   }

   return SCIP_OKAY;
}

/** reads the semicontinuous section */
static
SCIP_RETCODE readSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   while( getNextToken(lpinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(lpinput) )
         return SCIP_OKAY;

      /* semi-continuous variables are not yet supported by SCIP */
      syntaxError(scip, lpinput, "semi-continuous variables not yet supported by SCIP");
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** reads the sos section
 *
 *  The format is as follows:
 *
 *  SOS
 *  <constraint name>: [S1|S2]:: {<variable name>:<weight>}
 *  ...
 *  <constraint name>: [S1|S2]:: {<variable name>:<weight>}
 * */
static
SCIP_RETCODE readSos(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Bool initial, separate, enforce, check, propagate;
   SCIP_Bool local, modifiable, dynamic, removable;
   char name[SCIP_MAXSTRLEN];
   int cnt = 0;

   assert(lpinput != NULL);

   /* standard settings for SOS constraints: */
   initial = TRUE;
   separate = FALSE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;
   removable = FALSE;

   while ( getNextToken(lpinput) )
   {
      int type = -1;
      SCIP_CONS* cons;
      SCIP_Bool nextSection = FALSE;

      /* check if we reached a new section */
      if ( isNewSection(lpinput) )
         return SCIP_OKAY;

      /* check for an SOS constraint name */
      *name = '\0';

      /* remember the token in the token buffer */
      swapTokenBuffer(lpinput);

      /* get the next token and check, whether it is a colon */
      if ( getNextToken(lpinput) )
      {
         if ( strcmp(lpinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the constraint name */
            strncpy(name, lpinput->tokenbuf, SCIP_MAXSTRLEN);
            name[SCIP_MAXSTRLEN-1] = '\0';
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse it next */
            pushToken(lpinput);
            pushBufferToken(lpinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it next */
         pushBufferToken(lpinput);
      }

      /* get type */
      if ( ! getNextToken(lpinput) )
      {
	 syntaxError(scip, lpinput, "expected SOS type: 'S1::' or 'S2::'");
	 return SCIP_OKAY;
      }
      /* check whether constraint name was left out */
      if ( strcmp(lpinput->token, ":") == 0 )
      {
	 /* we have to push twice ':' and once the type: */
	 pushToken(lpinput);
	 lpinput->token[0] = ':';
	 lpinput->token[1] = '\0';
	 pushToken(lpinput);
	 swapTokenBuffer(lpinput);

	 /* set artificial name */
	 snprintf(name, SCIP_MAXSTRLEN, "SOS%d", ++cnt);
      }

      /* check whether it is type 1 or type 2 */
      if ( strcmp(lpinput->token, "S1") == 0 )
      {
	 type = 1;
	 SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
				       local, modifiable, dynamic, removable) );
      }
      else if ( strcmp(lpinput->token, "S2") == 0 )
      {
	 type = 2;
	 SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
				       local, modifiable, dynamic, removable) );
      }
      else
      {
	 syntaxError(scip, lpinput, "SOS constraint type other than 1 or 2 appeared");
	 return SCIP_OKAY;
      }
      assert( type == 1 || type == 2 );

      SCIPdebugMessage("created SOS%d constraint <%s>\n", type, name);

      /* make sure that a colons follows */
      if ( ! getNextToken(lpinput) || strcmp(lpinput->token, ":") != 0 )
      {
	 syntaxError(scip, lpinput, "SOS constraint type has to be followed by two colons");
	 return SCIP_OKAY;
      }

      /* make sure that another colons follows */
      if ( ! getNextToken(lpinput) || strcmp(lpinput->token, ":") != 0 )
      {
	 syntaxError(scip, lpinput, "SOS constraint type has to be followed by two colons");
	 return SCIP_OKAY;
      }

      /* parse elements of SOS constraint */
      while ( getNextToken(lpinput) && ! nextSection )
      {
	 SCIP_VAR* var;
	 SCIP_Real weight;

	 /* check if we reached a new section */
	 if ( isNewSection(lpinput) )
	 {
	    nextSection = TRUE;
	    break;
	 }

	 /* remember the token in the token buffer */
	 swapTokenBuffer(lpinput);

	 /* get variable and colon */
	 var = SCIPfindVar(scip, lpinput->tokenbuf);

	 /* if token is a variable name */
	 if ( var == NULL )
	 {
	    pushBufferToken(lpinput);
	    break;
	 }
	 else
	 {
	    SCIPdebugMessage("found variable <%s>\n", SCIPvarGetName(var));
	    if ( ! getNextToken(lpinput) || strcmp(lpinput->token, ":") != 0 )
	    {
	       syntaxError(scip, lpinput, "expected colon and weight.");
	       return SCIP_OKAY;
	    }
	    /* check next token */
	    if ( ! getNextToken(lpinput) )
	    {
	       /* push back token, since it could be the name of a new constraint */
	       pushToken(lpinput);
	       pushBufferToken(lpinput);
	       break;
	    }
	    else
	    {
	       /* get weight */
	       if ( ! isValue(scip, lpinput, &weight) )
	       {
		  /* push back token, since it could be the name of a new constraint */
		  pushToken(lpinput);
		  pushBufferToken(lpinput);
		  break;
	       }
	       else
	       {
		  /* we now know that we have a variable/weight pair -> add variable*/
		  switch (type)
		  {
		  case 1: SCIP_CALL( SCIPaddVarSOS1(scip, cons, var, weight) ); break;
		  case 2: SCIP_CALL( SCIPaddVarSOS2(scip, cons, var, weight) ); break;
		  default: abort(); // should not happen
		  }
		  SCIPdebugMessage("added variable <%s> with weight %g.\n", SCIPvarGetName(var), weight);
	       }
	    }
	 }
      }

      /* add the SOS constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("(line %d) added constraint <%s>: ", lpinput->linenumber, SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** reads an LP file */
static
SCIP_RETCODE readLPFile(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   const char*           filename            /**< name of the input file */
   )
{
   assert(lpinput != NULL);

   /* open file */
   lpinput->file = SCIPfopen(filename, "r");
   if( lpinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* parse the file */
   lpinput->section = LP_START;
   while( lpinput->section != LP_END )
   {
      switch( lpinput->section )
      {
      case LP_START:
         SCIP_CALL( readStart(scip, lpinput) );
         break;

      case LP_OBJECTIVE:
         SCIP_CALL( readObjective(scip, lpinput) );
         break;

      case LP_CONSTRAINTS:
         SCIP_CALL( readConstraints(scip, lpinput) );
         break;

      case LP_BOUNDS:
         SCIP_CALL( readBounds(scip, lpinput) );
         break;

      case LP_GENERALS:
         SCIP_CALL( readGenerals(scip, lpinput) );
         break;

      case LP_BINARIES:
         SCIP_CALL( readBinaries(scip, lpinput) );
         break;

      case LP_SEMICONTINUOUS:
         SCIP_CALL( readSemicontinuous(scip, lpinput) );
         break;

      case LP_SOS:
         SCIP_CALL( readSos(scip, lpinput) );
         break;

      case LP_END: /* this is already handled in the while() loop */
      default:
         SCIPerrorMessage("invalid LP file section <%d>\n", lpinput->section);
         return SCIP_INVALIDDATA;
      }
   }

   /* close file */
   SCIPfclose(lpinput->file);

   return SCIP_OKAY;
}


/*
 * Local methods (for writing)
 */

/** transforms given variables, scalars and constant to the corresponding active variables, scalars and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */ 
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize) );
            
      if( requiredsize > *nvars )
      {
         *nvars = requiredsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, *nvars ) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, *nvars ) );
         
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
      }
   }
   return SCIP_OKAY;
}


/* print row in LP format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   int*                  consnumber,         /**< pointer to constraint counter (or NULL) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /** row type ("=", "<=", or ">=") */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvars,              /**< number of variables */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int v;

   assert( scip != NULL );
   assert( strcmp(type, "=") == 0 || strcmp(type, "<=") || strcmp(type, ">=") );
   
   /* print row name */
   if( consnumber == NULL )
      SCIPinfoMessage(scip, file, " %s%s: ", rowname, rownameextension );
   else
   {
      (*consnumber)++;
      SCIPinfoMessage(scip, file, " c%d: ", *consnumber);
   }
   
   /* print coefficients */
   for( v = 0; v < nvars; ++v )
   {
      SCIPinfoMessage(scip, file, "%+g", vals[v]);
      if( consnumber != NULL )
         SCIPinfoMessage(scip, file, "x%d ", SCIPvarGetProbindex(vars[v]) + 1 );
      else
         SCIPinfoMessage(scip, file, " %s ", SCIPvarGetName(vars[v]));
   }
   
   /* print left hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;
   
   SCIPinfoMessage(scip, file, "%s %g\n", type, rhs );
}



/** prints given linear constraint information in LP format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   int*                  consnumber,         /**< pointer to constraint counter (or NULL) */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */ 
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( rowname != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;
   
   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars );
   if( vals != NULL )
      SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );
   
      for( v = 0; v < nactivevars; ++v )
      {
         activevals[v] = 1.0;
      }
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );
   
   /* print row(s) in LP format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );
      
      /* equal constrain */
      printRow(scip, file, consnumber, rowname, "", "=", activevars, activevals, nactivevars, rhs - activeconstant);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRow(scip, file, consnumber, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", ">=", 
            activevars, activevals, nactivevars, lhs - activeconstant);
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         printRow(scip, file, consnumber, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "<=", 
            activevars, activevals, nactivevars, rhs - activeconstant);
      }
   }
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}



/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeLp NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadLp)
{  /*lint --e{715}*/
   
   SCIP_CALL( SCIPreadLp(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteLp)
{  /*lint --e{715}*/
   
   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "/ STATISTICS\n");
   SCIPinfoMessage(scip, file, "/   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "/   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   
   SCIPinfoMessage(scip, file, "/   Constraints      : %d initial, %d maximal\n", startnconss, maxnconss);
   
   SCIP_CALL( SCIPwriteLp(scip, file, FALSE, objsense, vars, nvars, conss, nconss, result) );
   
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the lp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderLp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create lp reader data */
   readerdata = NULL;

   /* include lp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeLp, readerReadLp, readerWriteLp, readerdata) );

   /* add lp reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/lpreader/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/lpreader/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/lpreader/dynamicrows", "should rows be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}


/* reads problem from file */
SCIP_RETCODE SCIPreadLp(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   )
{
   LPINPUT lpinput;
   int i;

   /* initialize LP input data */
   lpinput.file = NULL;
   lpinput.linebuf[0] = '\0';
   lpinput.probname[0] = '\0';
   lpinput.objname[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpinput.token, LP_MAX_LINELEN) );
   lpinput.token[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpinput.tokenbuf, LP_MAX_LINELEN) );
   lpinput.tokenbuf[0] = '\0';
   for( i = 0; i < LP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &lpinput.pushedtokens[i], LP_MAX_LINELEN) );
   }

   lpinput.npushedtokens = 0;
   lpinput.linenumber = 0;
   lpinput.linepos = 0;
   lpinput.section = LP_START;
   lpinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   lpinput.inlazyconstraints = FALSE;
   lpinput.inusercuts = FALSE;
   lpinput.haserror = FALSE;

   /* read the file */
   SCIP_CALL( readLPFile(scip, &lpinput, filename) );

   /* free dynamically allocated memory */
   SCIPfreeMemoryArray(scip, &lpinput.token);
   SCIPfreeMemoryArray(scip, &lpinput.tokenbuf);
   for( i = 0; i < LP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeMemoryArray(scip, &lpinput.pushedtokens[i]);
   }

   /* evaluate the result */
   if( lpinput.haserror )
      return SCIP_PARSEERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, lpinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/* writes problem to file */
SCIP_RETCODE SCIPwriteLp(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file, or NULL if standard output should be used */
   SCIP_Bool          genericnames,       /**< use generic variable and row names? */
   SCIP_OBJSENSE      objsense,           /**< objective sense */
   SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, 
                                           *   continuous */
   int                nvars,              /**< number of mutable variables in the problem */
   SCIP_CONS**        conss,              /**< array with constraints of the problem */
   int                nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*       result              /**< pointer to store the result of the file writing call */
   )
{
   int c,v;
   int counter = 0;
   
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;
   const char* consname;
   SCIP_Bool transformed;
   int* consnumber;
   
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   
   assert( scip != NULL );

   if( genericnames )
      consnumber = &counter;
   else
      consnumber = NULL;
   
   /* print objective sense */
   SCIPinfoMessage(scip, file, "%s\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "Minimize" : "Maximize");
   
   /* print objective function */
   SCIPinfoMessage(scip, file, " obj:");
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];

      if(SCIPisZero(scip, SCIPvarGetObj(var)) )
         continue;
      
      if( genericnames )
         SCIPinfoMessage(scip, file, " %+gx%d", SCIPvarGetObj(var), SCIPvarGetProbindex(var) + 1 ); 
      else
         SCIPinfoMessage(scip, file, " %+g %s", SCIPvarGetObj(var), SCIPvarGetName(var) ); 
      
      if( (v + 1)  % 10 == 0 )
         SCIPinfoMessage(scip, file, "\n     ");
   }
   
   /* print "Subsect to" section */
   SCIPinfoMessage(scip, file, "\n\nSubject to\n");

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL);

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );
      
      consname = SCIPconsGetName(cons);
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      transformed = SCIPconsIsTransformed(cons);
      
      if( strcasecmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file, consname,  consnumber,
               SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
               SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcasecmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);
         
         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printLinearCons(scip, file, consname,  consnumber,
                  consvars, NULL, nconsvars, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printLinearCons(scip, file, consname, consnumber,
                  consvars, NULL, nconsvars, -SCIP_DEFAULT_INFINITY, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printLinearCons(scip, file, consname, consnumber,
                  consvars, NULL, nconsvars, 1.0, SCIP_DEFAULT_INFINITY, transformed) );
            break;
         }
      }
      else if( strcasecmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file, consname, consnumber,
               SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons), 
               1.0, SCIP_DEFAULT_INFINITY, transformed) );
      }
      else if( strcasecmp(conshdlrname, "knapsack") == 0 )
      {
         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         SCIP_Longint* weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
         {
            consvals[v] = weights[v];
         }
         
         SCIP_CALL( printLinearCons(scip, file, consname, consnumber,
               consvars, consvals, nconsvars, 
               -SCIP_DEFAULT_INFINITY, SCIPgetCapacityKnapsack(scip, cons), transformed) );
         
         SCIPfreeBufferArray(scip, &consvals);

      }
      else if( strcasecmp(conshdlrname, "varbounds") == 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( printLinearCons(scip, file, consname, consnumber,
               consvars, consvals, 2, 
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );
         
         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else
      {
         SCIPwarningMessage("constraint handler <%s> can not print requested format\n", conshdlrname );
      }
   }

   /* print "Bounds" section */
   SCIPinfoMessage(scip, file, "\n\nBounds\n");
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );
      
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
      
      if( !SCIPisInfinity(scip, -lb) || !SCIPisInfinity(scip, ub) )
      {
         /* print lower bound as far this one is not infinity */
         if( !SCIPisInfinity(scip, lb) )
            SCIPinfoMessage(scip, file, " %g <=", lb);
         
         /* print variable name */
         if( genericnames )
            SCIPinfoMessage(scip, file, " x%d ", SCIPvarGetProbindex(var) + 1 ); 
         else
            SCIPinfoMessage(scip, file, " %s ", SCIPvarGetName(var) ); 
         
         /* print upper bound as far this one is not infinity */
         if( !SCIPisInfinity(scip, ub) )
            SCIPinfoMessage(scip, file, "<= %g", ub);
         
         SCIPinfoMessage(scip, file, "\n");
      }
   }
   
   /* print generals section */
   if( nvars > 0 )
   {
      SCIPinfoMessage(scip, file, "\n\nGenerals\n");
      
      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         assert( var != NULL );

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS)
            continue;

         /* print variable name */
         if( genericnames )
            SCIPinfoMessage(scip, file, " x%d", SCIPvarGetProbindex(var) + 1 ); 
         else
            SCIPinfoMessage(scip, file, " %s", SCIPvarGetName(var) ); 
         
         if( (v + 1)  % 10 == 0 )
            SCIPinfoMessage(scip, file, "\n");
      }
   }
   
   /* end of lp format */
   SCIPmessageFPrintInfo(file, "\n\nEnd");
   
   *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}


