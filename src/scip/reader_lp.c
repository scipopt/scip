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

/**@file   reader_lp.c
 * @brief  LP file reader
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Stefan Heinz
 * @author Stefan Vigerske
 * @author Lars Schewe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strncasecmp() */
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
#include "scip/cons_indicator.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_soc.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/pub_misc.h"

#define READER_NAME             "lpreader"
#define READER_DESC             "file reader for MIPs in IBM CPLEX's LP file format"
#define READER_EXTENSION        "lp"


/*
 * Data structures
 */
#define LP_MAX_LINELEN       65536
#define LP_MAX_PUSHEDTOKENS  2
#define LP_INIT_COEFSSIZE    8192
#define LP_INIT_QUADCOEFSSIZE 16
#define LP_MAX_PRINTLEN      561       /**< the maximum length of any line is 560 + '\\0' = 561*/
#define LP_MAX_NAMELEN       256       /**< the maximum length for any name is 255 + '\\0' = 256 */
#define LP_PRINTLEN          100

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
   char                 linebuf[LP_MAX_LINELEN+1];
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
   SCIP_Bool            comment;
   SCIP_Bool            endline;
};
typedef struct LpInput LPINPUT;

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
   char formatstr[256];

   assert(lpinput != NULL);

   SCIPerrorMessage("Syntax error in line %d ('%s'): %s \n", lpinput->linenumber, lpinput->token, msg);
   if( lpinput->linebuf[strlen(lpinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", lpinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", lpinput->linebuf);
   }
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", lpinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, (const char*)formatstr, "^");
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
   switch (c)
   {
   case ' ':
   case '\f':
   case '\n':
   case '\r':
   case '\t':
   case '\v':
   case '\0':
      return TRUE;
   default:
      return FALSE;
   }
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char                  c                   /**< input character */
   )
{
   switch (c)
   {
   case '-':
   case '+':
   case ':':
   case '<':
   case '>':
   case '=':
   case '[':
   case ']':
   case '*':
   case '^':
      return TRUE;
   default:
      return FALSE;
   }
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

   if( isdigit((unsigned char)c) )
      return TRUE;
   else if( (*exptype == LP_EXP_NONE) && !(*hasdot) && (c == '.') && isdigit((unsigned char)nextc) )
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
      else if( isdigit((unsigned char)nextc) )
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
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   int i;

   assert(lpinput != NULL);

   /* if we previously detected a comment we have to parse the remaining line away if there is something left */
   if( !lpinput->endline && lpinput->comment )
   {
      SCIPdebugMessage("Throwing rest of comment away.\n");

      do
      {
         lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';
         (void)SCIPfgets(lpinput->linebuf, sizeof(lpinput->linebuf), lpinput->file);
      }
      while( lpinput->linebuf[LP_MAX_LINELEN-2] != '\0' );

      lpinput->comment = FALSE;
      lpinput->endline = TRUE;
   }

   /* read next line */
   lpinput->linepos = 0;
   lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';

   if( SCIPfgets(lpinput->linebuf, sizeof(lpinput->linebuf), lpinput->file) == NULL )
   {
      /* clear the line, this is really necessary here! */
      BMSclearMemoryArray(lpinput->linebuf, LP_MAX_LINELEN);

      return FALSE;
   }

   lpinput->linenumber++;

   /* if line is too long for our buffer correct the buffer and correct position in file */
   if( lpinput->linebuf[LP_MAX_LINELEN-2] != '\0' )
   {
      char* last;

      /* buffer is full; erase last token since it might be incomplete */
      lpinput->endline = FALSE;
      last = strrchr(lpinput->linebuf, ' ');

      if( last == NULL )
      {
         SCIPwarningMessage(scip, "we read %d characters from the file; this might indicate a corrupted input file!",
            LP_MAX_LINELEN - 2);
         lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';
         SCIPdebugMessage("the buffer might be corrupted\n");
      }
      else
      {
         SCIPfseek(lpinput->file, -(long) strlen(last) - 1, SEEK_CUR);
         SCIPdebugMessage("correct buffer, reread the last %ld characters\n", (long) strlen(last) + 1);
         *last = '\0';
      }
   }
   else
   {
      /* found end of line */
      lpinput->endline = TRUE;
   }
   lpinput->linebuf[LP_MAX_LINELEN-1] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
   lpinput->comment = FALSE;

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(lpinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

         lpinput->comment = TRUE;
         break;
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
   SCIP*                 scip,               /**< SCIP data structure */
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
         if( !getNextLine(scip, lpinput) )
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

      /* if the token is a power sign '^', skip a following '2'
       * if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1 && lpinput->token[tokenlen-1] == '^' && buf[lpinput->linepos] == '2' )
      {
         lpinput->linepos++;
      }
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
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Bool iscolon;

   assert(lpinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(lpinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(scip, lpinput) )
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
      if( getNextToken(scip, lpinput) )
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
      if( getNextToken(scip, lpinput) )
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

   if( strcasecmp(lpinput->token, "ST") == 0
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
      if( getNextToken(scip, lpinput) )
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
      if( getNextToken(scip, lpinput) )
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
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( !getNextToken(scip, lpinput) )
         return SCIP_OKAY;
   }
   while( !isNewSection(scip, lpinput) );

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   SCIP_Bool             isobjective,        /**< indicates whether we are currently reading the coefficients of the objective */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   LP_MAX_LINELEN */
   SCIP_VAR***           vars,               /**< pointer to store the array with variables (must be freed by caller) */
   SCIP_Real**           coefs,              /**< pointer to store the array with coefficients (must be freed by caller) */
   int*                  ncoefs,             /**< pointer to store the number of coefficients */
   SCIP_VAR***           quadvars1,          /**< pointer to store the array with first variables in quadratic terms (must be freed by caller) */
   SCIP_VAR***           quadvars2,          /**< pointer to store the array with second variables in quadratic terms (must be freed by caller) */
   SCIP_Real**           quadcoefs,          /**< pointer to store the array with coefficients in quadratic terms (must be freed by caller) */
   int*                  nquadcoefs,         /**< pointer to store the number of quadratic coefficients */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;
   int coefssize;
   SCIP_Bool inquadpart;
   int quadcoefssize;
   SCIP_VAR* firstquadvar;

   assert(lpinput != NULL);
   assert(name != NULL);
   assert(vars != NULL);
   assert(coefs != NULL);
   assert(ncoefs != NULL);
   assert(quadvars1 != NULL);
   assert(quadvars2 != NULL);
   assert(quadcoefs != NULL);
   assert(nquadcoefs != NULL);
   assert(newsection != NULL);

   *vars = NULL;
   *coefs = NULL;
   *quadvars1 = NULL;
   *quadvars2 = NULL;
   *quadcoefs = NULL;
   *name = '\0';
   *ncoefs = 0;
   *nquadcoefs = 0;
   *newsection = FALSE;
   inquadpart = FALSE;

   /* read the first token, which may be the name of the line */
   if( getNextToken(scip, lpinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* remember the token in the token buffer */
      swapTokenBuffer(lpinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(scip, lpinput) )
      {
         if( strcmp(lpinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the line name */
            (void)SCIPmemccpy(name, lpinput->tokenbuf, '\0', LP_MAX_LINELEN);

            name[LP_MAX_LINELEN - 1] = '\0';
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
   
   quadcoefssize = LP_INIT_QUADCOEFSSIZE;
   SCIP_CALL( SCIPallocMemoryArray(scip, quadvars1, quadcoefssize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, quadvars2, quadcoefssize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, quadcoefs, quadcoefssize) );

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   firstquadvar = NULL;
   *ncoefs = 0;
   *nquadcoefs = 0;
   while( getNextToken(scip, lpinput) )
   {
      SCIP_VAR* var;
      
      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
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
      
      /* check if we start a quadratic part */
      if( strcmp(lpinput->token, "[") == 0 )
      {
         if( inquadpart )
         {
            syntaxError(scip, lpinput, "cannot start quadratic part while already in quadratic part.");
            return SCIP_OKAY;
         }
         if( havesign && coefsign != +1 )
         {
            syntaxError(scip, lpinput, "cannot have '-' in front of quadratic part.");
            return SCIP_OKAY;
         }
         if( havevalue )
         {
            syntaxError(scip, lpinput, "cannot have value in front of quadratic part.");
            return SCIP_OKAY;
         }

         SCIPdebugMessage("(line %d) start quadratic part\n", lpinput->linenumber);
         inquadpart = TRUE;
         continue;
      }
      
      /* check if we end a quadratic part */
      if( *lpinput->token == ']' )
      {
         if( !inquadpart )
         {
            syntaxError(scip, lpinput, "cannot end quadratic part before starting one.");
            return SCIP_OKAY;
         }
         if( havesign || havevalue || firstquadvar != NULL )
         {
            if( firstquadvar == NULL )
            {
               syntaxError(scip, lpinput, "expected value or first quadratic variable.");
            }
            else
            {
               syntaxError(scip, lpinput, "expected second quadratic variable.");
            }
            return SCIP_OKAY;
         }
         
         SCIPdebugMessage("(line %d) end quadratic part\n", lpinput->linenumber);
         inquadpart = FALSE;
         
         if( isobjective )
         {
            /* quadratic part in objective has to end with '/2' */
            if( !getNextToken(scip, lpinput) )
            {
               syntaxError(scip, lpinput, "expected '/2' after end of quadratic part in objective.");
               return SCIP_OKAY;
            }
            if( strcmp(lpinput->token, "/2") == 0 )
            {
               SCIPdebugMessage("(line %d) saw '/2' after quadratic part in objective\n", lpinput->linenumber);
            }
            else if( strcmp(lpinput->token, "/") == 0 )
            {
               /* maybe it says '/ 2' */
               if( !getNextToken(scip, lpinput) || strcmp(lpinput->token, "2") != 0 )
               {
                  syntaxError(scip, lpinput, "expected '/2' after end of quadratic part in objective.");
                  return SCIP_OKAY;
               }
               SCIPdebugMessage("(line %d) saw '/ 2' after quadratic part in objective\n", lpinput->linenumber);
            }
            else
            {
               syntaxError(scip, lpinput, "expected '/2' after end of quadratic part in objective.");
               return SCIP_OKAY;
            }
         }
         
         continue;
      }

      /* check if we are in between two quadratic variables */
      if( strcmp(lpinput->token, "*") == 0 )
      {
         if( !inquadpart )
         {
            syntaxError(scip, lpinput, "cannot have '*' outside of quadratic part.");
            return SCIP_OKAY;
         }
         if( firstquadvar == NULL )
         {
            syntaxError(scip, lpinput, "cannot have '*' before first variable in quadratic term.");
            return SCIP_OKAY;
         }
         
         continue;
      }
      
      /* all but the first coefficient need a sign */
      if( !inquadpart && *ncoefs > 0 && !havesign )
      {
         syntaxError(scip, lpinput, "expected sign ('+' or '-') or sense ('<' or '>').");
         return SCIP_OKAY;
      }
      if(  inquadpart && *nquadcoefs > 0 && !havesign )
      {
         syntaxError(scip, lpinput, "expected sign ('+' or '-').");
         return SCIP_OKAY;
      }

      /* check if the last variable should be squared */
      if( strcmp(lpinput->token, "^") == 0 )
      {
         if( !inquadpart )
         {
            syntaxError(scip, lpinput, "cannot have squares ('^2') outside of quadratic part.");
            return SCIP_OKAY;
         }
         if( firstquadvar == NULL )
         {
            syntaxError(scip, lpinput, "cannot have square '^2' before variable.");
            return SCIP_OKAY;
         }
         
         var = firstquadvar;
      }
      else
      {
         /* check if we read a value */
         if( isValue(scip, lpinput, &coef) )
         {
            SCIPdebugMessage("(line %d) read coefficient value: %g with sign %+d\n", lpinput->linenumber, coef, coefsign);
            if( havevalue )
            {
               syntaxError(scip, lpinput, "two consecutive values.");
               return SCIP_OKAY;
            }
            havevalue = TRUE;
            continue;
         }

         /* the token is a variable name: get the corresponding variable (or create a new one) */
         SCIP_CALL( getVariable(scip, lpinput->token, &var, NULL) );
      }

      if( !inquadpart )
      {
         /* insert the linear coefficient */
         SCIPdebugMessage("(line %d) read linear coefficient: %+g<%s>\n", lpinput->linenumber, coefsign * coef, SCIPvarGetName(var));
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
      }
      else
      {
         if( firstquadvar == NULL )
         {
            /* if first quadratic variable read, store it and continue; expect second one in next round */
            firstquadvar = var;
            continue;
         }
         
         /* insert the quadratic coefficient */
         SCIPdebugMessage("(line %d) read quadratic coefficient: %+g<%s><%s>\n", lpinput->linenumber, (isobjective ? 0.5 : 1) * coefsign * coef, SCIPvarGetName(firstquadvar), SCIPvarGetName(var));
         if( !SCIPisZero(scip, coef) )
         {
            /* resize the vars and coefs array if needed */
            if( *nquadcoefs >= quadcoefssize )
            {
               quadcoefssize *= 2;
               quadcoefssize = MAX(quadcoefssize, (*nquadcoefs)+1);
               SCIP_CALL( SCIPreallocMemoryArray(scip, quadvars1, quadcoefssize) );
               SCIP_CALL( SCIPreallocMemoryArray(scip, quadvars2, quadcoefssize) );
               SCIP_CALL( SCIPreallocMemoryArray(scip, quadcoefs, quadcoefssize) );
            }
            assert(*nquadcoefs < quadcoefssize);

            /* add coefficient */
            (*quadvars1)[*nquadcoefs] = firstquadvar;
            (*quadvars2)[*nquadcoefs] = var;
            (*quadcoefs)[*nquadcoefs] = coefsign * coef;
            if( isobjective )
               (*quadcoefs)[*nquadcoefs] /= 2.0;
            (*nquadcoefs)++;
         }
      }

      /* reset the flags and coefficient value for the next coefficient */
      coefsign = +1;
      coef = 1.0;
      havesign = FALSE;
      havevalue = FALSE;
      firstquadvar = NULL;
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
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
   int nquadcoefs;
   SCIP_Bool newsection;

   assert(lpinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, lpinput, TRUE, name, &vars, &coefs, &ncoefs, &quadvars1, &quadvars2, &quadcoefs, &nquadcoefs, &newsection) );
   if( !hasError(lpinput) )
   {
      int i;

      /* set the linear objective values */
      for( i = 0; i < ncoefs; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(scip, vars[i], SCIPvarGetObj(vars[i]) + coefs[i]) );
      }

      /* insert dummy variable and constraint to represent quadratic part of objective */
      if( nquadcoefs > 0 )
      {
         SCIP_VAR*  quadobjvar;
         SCIP_CONS* quadobjcons;
         SCIP_Real  lhs;
         SCIP_Real  rhs;
         SCIP_Real  minusone;

         SCIP_CALL( SCIPcreateVar(scip, &quadobjvar, "quadobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, quadobjvar) );

         if( lpinput->objsense == SCIP_OBJSENSE_MINIMIZE )
         {
            lhs = -SCIPinfinity(scip);
            rhs = 0.0;
         }
         else
         {
            lhs = 0.0;
            rhs = SCIPinfinity(scip);
         }

         minusone = -1.0;
         SCIP_CALL( SCIPcreateConsQuadratic(scip, &quadobjcons, "quadobj", 1, &quadobjvar, &minusone, nquadcoefs, quadvars1, quadvars2, quadcoefs, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPaddCons(scip, quadobjcons) );
         SCIPdebugMessage("(line %d) added constraint <%s> to represent quadratic objective: ", lpinput->linenumber, SCIPconsGetName(quadobjcons));
         SCIPdebugPrintCons(scip, quadobjcons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &quadobjcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &quadobjvar) );
      }
   }

   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &vars);
   SCIPfreeMemoryArrayNull(scip, &coefs);
   SCIPfreeMemoryArrayNull(scip, &quadvars1);
   SCIPfreeMemoryArrayNull(scip, &quadvars2);
   SCIPfreeMemoryArrayNull(scip, &quadcoefs);

   return SCIP_OKAY;
}

/** create indicator constraint */
static
SCIP_RETCODE createIndicatorConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   const char*           name,               /**< name of indicator constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable */
   SCIP_Real             binvalue            /**< value of indicator part (0/1) */
   )
{
   char name2[LP_MAX_LINELEN];
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlincoefs;
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
   int nquadcoefs;
   SCIP_Bool newsection;
   LPSENSE linsense;
   SCIP_Real linsidevalue;
   int linsidesign;
   SCIP_CONS* cons;
   SCIP_Real linrhs;
   int j;
   SCIP_Bool linConsEQ;

   SCIP_Bool dynamicconss;
   SCIP_Bool dynamicrows;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool dynamic;
   SCIP_Bool removable;

   assert( lpinput != NULL );
   assert( binvar != NULL );

   /* check that binvalue is 0 or 1 */
   if( !SCIPisFeasEQ(scip, binvalue, 0.0) && !SCIPisFeasEQ(scip, binvalue, 1.0) )
   {
      syntaxError(scip, lpinput, "value for binary variable must be '0' or '1'.");
      return SCIP_OKAY;
   }

   if( SCIPisFeasEQ(scip, binvalue, 0.0) )
   {
      SCIP_VAR* negbinvar;
      SCIP_Bool infeasible;

      /* At this point we force the variable binvar to be binary, since we need the negated variable. We have to check
       * later whether the type of the variable specified in the file agrees with this specification. 
       */
      /* check whether bounds are correct - might already been set if variable is used in another indicator constraint */
      if( SCIPvarGetLbGlobal(binvar) < 0.0 )
         SCIP_CALL( SCIPchgVarLb(scip, binvar, 0.0) );
      if( SCIPvarGetUbGlobal(binvar) > 1.0 )
         SCIP_CALL( SCIPchgVarUb(scip, binvar, 1.0) );
      SCIP_CALL( SCIPchgVarType(scip, binvar, SCIP_VARTYPE_BINARY, &infeasible) );
      /* don't assert feasibility here because the presolver will and should detect a infeasibility */

      SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &negbinvar) );
      binvar = negbinvar;
      assert( binvar != NULL );
   }

   /* read linear constraint */
   SCIP_CALL( readCoefficients(scip, lpinput, FALSE, name2, &linvars, &lincoefs, &nlincoefs, &quadvars1, &quadvars2, &quadcoefs, &nquadcoefs, &newsection) );
   if( hasError(lpinput) )
      goto TERMINATE;
   if( newsection )
   {
      syntaxError(scip, lpinput, "expected constraint.");
      goto TERMINATE;
   }
   if( nquadcoefs > 0 )
   {
      /* @todo could introduce auxiliary variable and move quadratic part into quadratic constraint? */
      syntaxError(scip, lpinput, "quadratic indicator constraints not supported.");
      goto TERMINATE;
   }
   if( name2[0] != '\0' )
   {
      syntaxError(scip, lpinput, "did not expect name for linear constraint.");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if( !getNextToken(scip, lpinput) || !isSense(lpinput, &linsense) )
   {
      syntaxError(scip, lpinput, "expected constraint sense '<=', '=', or '>='.");
      goto TERMINATE;
   }

   /* read the right hand side */
   linsidesign = +1;
   if( !getNextToken(scip, lpinput) )
   {
      syntaxError(scip, lpinput, "missing right hand side.");
      goto TERMINATE;
   }
   if( isSign(lpinput, &linsidesign) )
   {
      if( !getNextToken(scip, lpinput) )
      {
         syntaxError(scip, lpinput, "missing value of right hand side.");
         goto TERMINATE;
      }
   }
   if( !isValue(scip, lpinput, &linsidevalue) )
   {
      syntaxError(scip, lpinput, "expected value for right hand side.");
      goto TERMINATE;
   }
   linsidevalue *= linsidesign;

   /* assign the left and right hand side, depending on the constraint sense */
   linConsEQ = FALSE;
   switch( linsense )
   {
   case LP_SENSE_GE:
      linrhs = linsidevalue;
      for( j = 0; j < nlincoefs; ++j )
         lincoefs[j] *= -1;
      break;
   case LP_SENSE_LE:
      linrhs = linsidevalue;
      break;
   case LP_SENSE_EQ:
      linConsEQ = TRUE;
      linrhs = linsidevalue;
      break;
   case LP_SENSE_NOTHING:
   default:
      SCIPerrorMessage("invalid constraint sense <%d>\n", linsense);
      return SCIP_INVALIDDATA;
   }

   /* create and add the indicator constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamicrows", &dynamicrows) );
   initial = !dynamicrows && !lpinput->inlazyconstraints && !lpinput->inusercuts;
   separate = TRUE;
   enforce = !lpinput->inusercuts;
   check = !lpinput->inusercuts;
   propagate = TRUE;
   local = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows || lpinput->inusercuts;

   SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, binvar, nlincoefs, linvars, lincoefs, linrhs,
         initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugMessage("(line %d) created constraint%s: ", lpinput->linenumber,
      lpinput->inlazyconstraints ? " (lazy)" : (lpinput->inusercuts ? " (user cut)" : ""));
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create second constraint if it was an equation */
   if( linConsEQ )
   {
      for( j = 0; j < nlincoefs; ++j )
         lincoefs[j] *= -1;
      linrhs *= -1;
      SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, binvar, nlincoefs, linvars, lincoefs, linrhs,
            initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("(line %d) created constraint%s: ", lpinput->linenumber,
         lpinput->inlazyconstraints ? " (lazy)" : (lpinput->inusercuts ? " (user cut)" : ""));
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

 TERMINATE:
   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &linvars);
   SCIPfreeMemoryArrayNull(scip, &lincoefs);
   SCIPfreeMemoryArrayNull(scip, &quadvars1);
   SCIPfreeMemoryArrayNull(scip, &quadvars2);
   SCIPfreeMemoryArrayNull(scip, &quadcoefs);

   return SCIP_OKAY;
}

/** reads the constraints section 
 *
 *  Read linear and indicator constraints.
 *
 *  The CPLEX manual says that indicator constraints are of the following form:
 *
 *  [constraintname:]  binaryvariable = value  ->  linear constraint
 *
 *  We also accept "<->".
 */
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
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
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
   int nquadcoefs;
   int sidesign;
   SCIP_Bool isIndicatorCons;

   assert(lpinput != NULL);

   /* read coefficients */
   SCIP_CALL( readCoefficients(scip, lpinput, FALSE, name, &vars, &coefs, &ncoefs, &quadvars1, &quadvars2, &quadcoefs, &nquadcoefs, &newsection) );
   if( hasError(lpinput) )
      goto TERMINATE;
   if( newsection )
   {
      if( ncoefs > 0 || nquadcoefs > 0 )
         syntaxError(scip, lpinput, "expected constraint sense '<=', '=', or '>='.");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if( !getNextToken(scip, lpinput) || !isSense(lpinput, &sense) )
   {
      syntaxError(scip, lpinput, "expected constraint sense '<=', '=', or '>='.");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if( !getNextToken(scip, lpinput) )
   {
      syntaxError(scip, lpinput, "missing right hand side.");
      goto TERMINATE;
   }
   if( isSign(lpinput, &sidesign) )
   {
      if( !getNextToken(scip, lpinput) )
      {
         syntaxError(scip, lpinput, "missing value of right hand side.");
         goto TERMINATE;
      }
   }
   if( !isValue(scip, lpinput, &sidevalue) )
   {
      syntaxError(scip, lpinput, "expected value as right hand side.");
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
      SCIPerrorMessage("invalid constraint sense <%d>.\n", sense);
      return SCIP_INVALIDDATA;
   }

   /* check whether we read the first part of an indicator constraint */
   isIndicatorCons = FALSE;
   if ( getNextToken(scip, lpinput) && !isNewSection(scip, lpinput) )
   {
      /* check whether we have '<' from a "<->" string */
      if ( *lpinput->token == '<' )
      {
         SCIP_Bool haveequiv = FALSE;
         int linepos = lpinput->linepos-1;

         /* check next token - cannot be a new section */
         if ( getNextToken(scip, lpinput) )
         {
            /* check for "<-" */
            if ( *lpinput->token == '-' )
            {
               /* check next token - cannot be a new section */
               if ( getNextToken(scip, lpinput) )
               {
                  /* check for "<->" */
                  if ( *lpinput->token == '>' )
                  {
                     lpinput->linepos = linepos;
                     strcpy(lpinput->token, "<");
                     syntaxError(scip, lpinput,
                        "SCIP does not support equivalence (<->) indicator constraints; consider using the \"->\" form.");
                     haveequiv = TRUE;
                  }
               }
            }
         }
         if ( ! haveequiv )
         {
            lpinput->linepos = linepos;
            strcpy(lpinput->token, "<");
            syntaxError(scip, lpinput, "unexpected \"<\".");
         }
         goto TERMINATE;
      }

      /* check for "->" */
      if ( *lpinput->token == '-' )
      {
         /* remember '-' in token buffer */
         swapTokenBuffer(lpinput);

         /* check next token - cannot be a new section */
         if( getNextToken(scip, lpinput) )
         {
            /* check for "->" */
            if ( *lpinput->token == '>' )
               isIndicatorCons = TRUE;
            else
            {
               /* push back last token and '-' */
               pushToken(lpinput);
               pushBufferToken(lpinput);
            }
         }
         else
            pushBufferToken(lpinput);
      }
      else
         pushToken(lpinput);
   }

   if( !isIndicatorCons )
   {
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
      if( nquadcoefs == 0 )
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, lhs, rhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsQuadratic(scip, &cons, name, ncoefs, vars, coefs,
               nquadcoefs, quadvars1, quadvars2, quadcoefs, lhs, rhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
      }
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("(line %d) created constraint%s: ", lpinput->linenumber,
         lpinput->inlazyconstraints ? " (lazy)" : (lpinput->inusercuts ? " (user cut)" : ""));
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      /* now we should have an indicator constraint */
      if( ncoefs != 1 || nquadcoefs > 0 )
      {
         syntaxError(scip, lpinput, "Indicator part can only consist of one binary variable.");
         goto TERMINATE;
      }
      if( !SCIPisEQ(scip, coefs[0], 1.0) )
      {
         syntaxError(scip, lpinput, "There cannot be a coefficient before the binary indicator variable.");
         goto TERMINATE;
      }
      if( sense != LP_SENSE_EQ )
      {
         syntaxError(scip, lpinput, "Indicator part cannot handle equations.");
         goto TERMINATE;
      }

      SCIP_CALL( createIndicatorConstraint(scip, lpinput, name, vars[0], lhs) );
   }

 TERMINATE:
   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &vars);
   SCIPfreeMemoryArrayNull(scip, &coefs);
   SCIPfreeMemoryArrayNull(scip, &quadvars1);
   SCIPfreeMemoryArrayNull(scip, &quadvars2);
   SCIPfreeMemoryArrayNull(scip, &quadcoefs);

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

   while( getNextToken(scip, lpinput) )
   {
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Real lb;
      SCIP_Real ub;
      int sign;
      SCIP_Bool hassign;
      LPSENSE leftsense;

      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
         return SCIP_OKAY;

      /* default bounds are [0,+inf] */
      lb = 0.0;
      ub = SCIPinfinity(scip);
      leftsense = LP_SENSE_NOTHING;

      /* check if the first token is a sign */
      sign = +1;
      hassign = isSign(lpinput, &sign);
      if( hassign && !getNextToken(scip, lpinput) )
      {
         syntaxError(scip, lpinput, "expected value.");
         return SCIP_OKAY;
      }

      /* the first token must be either a value or a variable name */
      if( isValue(scip, lpinput, &value) )
      {
         /* first token is a value: the second token must be a sense */
         if( !getNextToken(scip, lpinput) || !isSense(lpinput, &leftsense) )
         {
            syntaxError(scip, lpinput, "expected bound sense '<=', '=', or '>='.");
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
         syntaxError(scip, lpinput, "expected value.");
         return SCIP_OKAY;
      }
      else
         pushToken(lpinput);

      /* the next token must be a variable name */
      if( !getNextToken(scip, lpinput) )
      {
         syntaxError(scip, lpinput, "expected variable name.");
         return SCIP_OKAY;
      }
      SCIP_CALL( getVariable(scip, lpinput->token, &var, NULL) );

      /* the next token might be another sense, or the word "free" */
      if( getNextToken(scip, lpinput) )
      {
         LPSENSE rightsense;

         if( isSense(lpinput, &rightsense) )
         {
            /* check, if the senses fit */
            if( leftsense == LP_SENSE_NOTHING
               || (leftsense == LP_SENSE_LE && rightsense == LP_SENSE_LE)
               || (leftsense == LP_SENSE_GE && rightsense == LP_SENSE_GE) )
            {
               if( !getNextToken(scip, lpinput) )
               {
                  syntaxError(scip, lpinput, "expected value or sign.");
                  return SCIP_OKAY;
               }

               /* check if the next token is a sign */
               sign = +1;
               hassign = isSign(lpinput, &sign);
               if( hassign && !getNextToken(scip, lpinput) )
               {
                  syntaxError(scip, lpinput, "expected value.");
                  return SCIP_OKAY;
               }

               /* the next token must be a value */
               if( !isValue(scip, lpinput, &value) )
               {
                  syntaxError(scip, lpinput, "expected value.");
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
               syntaxError(scip, lpinput, "the two bound senses do not fit.");
               return SCIP_OKAY;
            }
         }
         else if( strcasecmp(lpinput->token, "FREE") == 0 )
         {
            if( leftsense != LP_SENSE_NOTHING )
            {
               syntaxError(scip, lpinput, "variable with bound is marked as 'free'.");
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

      /* change the bounds of the variable if bounds have been given (do not destroy earlier specification of bounds) */
      if( lb != 0.0 )
         SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      /*lint --e{777}*/
      if( ub != SCIPinfinity(scip) )
         SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      SCIPdebugMessage("(line %d) new bounds: <%s>[%g,%g]\n", lpinput->linenumber, SCIPvarGetName(var),
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
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

   while( getNextToken(scip, lpinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;
      SCIP_Bool infeasible;

      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpinput, "unknown variable in generals section.");
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
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   while( getNextToken(scip, lpinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;
      SCIP_Bool infeasible;

      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpinput, "unknown variable in binaries section.");
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

/** reads the semi-continuous section */
static
SCIP_RETCODE readSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Real oldlb;
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_Bool created;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamiccols;

   SCIP_VAR* vars[2];
   SCIP_BOUNDTYPE boundtypes[2];
   SCIP_Real bounds[2];

   assert(lpinput != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/lpreader/dynamiccols", &dynamiccols) );

   /* if section is titles "semi-continuous", then the parser breaks this into parts */
   if( strcasecmp(lpinput->token, "SEMI") == 0 )
   {
      if( !getNextToken(scip, lpinput) )
      {
         syntaxError(scip, lpinput, "unexpected end.");
         return SCIP_OKAY;
      }

      if( strcasecmp(lpinput->token, "-") == 0 )
      {
         if( !getNextToken(scip, lpinput) || strcasecmp(lpinput->token, "CONTINUOUS") != 0 )
         {
            syntaxError(scip, lpinput, "expected 'CONTINUOUS' after 'SEMI-'.");
            return SCIP_OKAY;
         }
      }
      else
      {
         pushToken(lpinput);
      }
   }

   while( getNextToken(scip, lpinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpinput, "unknown variable in semi-continuous section.");
         return SCIP_OKAY;
      }

      if( SCIPvarGetLbGlobal(var) <= 0.0 )
      {
         SCIPdebugMessage("ignore semi-continuity of variable <%s> with negative lower bound %g\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
         continue;
      }

      oldlb = SCIPvarGetLbGlobal(var);

      /* change the lower bound to 0.0 */
      SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );

      /* add a bound disjunction constraint to say var <= 0.0 or var >= oldlb */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "semicont_%s", SCIPvarGetName(var));

      vars[0] = var;
      vars[1] = var;
      boundtypes[0] = SCIP_BOUNDTYPE_UPPER;
      boundtypes[1] = SCIP_BOUNDTYPE_LOWER;
      bounds[0] = 0.0;
      bounds[1] = oldlb;

      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, 2, vars, boundtypes, bounds,
            !dynamiccols, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamiccols, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );

      SCIPdebugMessage("add bound disjunction constraint for semi-continuity of <%s>:\n\t", SCIPvarGetName(var));
      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** reads the sos section
 *
 *  The format is as follows:
 *
 *  SOS
 *  \<constraint name\>: [S1|S2]:: {\<variable name\>:\<weight\>}
 *  ...
 *  \<constraint name\>: [S1|S2]:: {\<variable name\>:\<weight\>}
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

   while( getNextToken(scip, lpinput) )
   {
      int type = -1;
      SCIP_CONS* cons;

      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
         return SCIP_OKAY;

      /* check for an SOS constraint name */
      *name = '\0';

      /* remember the token in the token buffer */
      swapTokenBuffer(lpinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(scip, lpinput) )
      {
         if( strcmp(lpinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the constraint name */
	    (void)SCIPmemccpy(name, lpinput->tokenbuf, '\0', SCIP_MAXSTRLEN);

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
      if( !getNextToken(scip, lpinput) )
      {
         syntaxError(scip, lpinput, "expected SOS type: 'S1::' or 'S2::'.");
         return SCIP_OKAY;
      }
      /* check whether constraint name was left out */
      if( strcmp(lpinput->token, ":") == 0 )
      {
         /* we have to push twice ':' and once the type: */
         pushToken(lpinput);
         lpinput->token[0] = ':';
         lpinput->token[1] = '\0';
         pushToken(lpinput);
         swapTokenBuffer(lpinput);

         /* set artificial name */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SOS%d", ++cnt);
      }

      /* check whether it is type 1 or type 2 */
      if( strcmp(lpinput->token, "S1") == 0 )
      {
         type = 1;
         SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
               local, modifiable, dynamic, removable) );
      }
      else if( strcmp(lpinput->token, "S2") == 0 )
      {
         type = 2;
         SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
               local, modifiable, dynamic, removable) );
      }
      else
      {
         syntaxError(scip, lpinput, "SOS constraint type other than 1 or 2 appeared.");
         return SCIP_OKAY;
      }
      assert( type == 1 || type == 2 );

      SCIPdebugMessage("created SOS%d constraint <%s>\n", type, name);

      /* make sure that a colons follows */
      if( !getNextToken(scip, lpinput) || strcmp(lpinput->token, ":") != 0 )
      {
         syntaxError(scip, lpinput, "SOS constraint type has to be followed by two colons.");
         return SCIP_OKAY;
      }

      /* make sure that another colons follows */
      if( !getNextToken(scip, lpinput) || strcmp(lpinput->token, ":") != 0 )
      {
         syntaxError(scip, lpinput, "SOS constraint type has to be followed by two colons.");
         return SCIP_OKAY;
      }

      /* parse elements of SOS constraint */
      while( getNextToken(scip, lpinput) )
      {
         SCIP_VAR* var;
         SCIP_Real weight;

         /* check if we reached a new section */
         if( isNewSection(scip, lpinput) )
            break;

         /* remember the token in the token buffer */
         swapTokenBuffer(lpinput);

         /* get variable and colon */
         var = SCIPfindVar(scip, lpinput->tokenbuf);

         /* if token is a variable name */
         if( var == NULL )
         {
            pushBufferToken(lpinput);
            break;
         }
         else
         {
            SCIPdebugMessage("found variable <%s>\n", SCIPvarGetName(var));
            if( !getNextToken(scip, lpinput) || strcmp(lpinput->token, ":") != 0 )
            {
               syntaxError(scip, lpinput, "expected colon and weight.");
               return SCIP_OKAY;
            }
            /* check next token */
            if( !getNextToken(scip, lpinput) )
            {
               /* push back token, since it could be the name of a new constraint */
               pushToken(lpinput);
               pushBufferToken(lpinput);
               break;
            }
            else
            {
               /* get weight */
               if( !isValue(scip, lpinput, &weight) )
               {
                  /* push back token, since it could be the name of a new constraint */
                  pushToken(lpinput);
                  pushBufferToken(lpinput);
                  break;
               }
               else
               {
                  /* we now know that we have a variable/weight pair -> add variable*/
                  switch( type )
                  {
                  case 1: 
                     SCIP_CALL( SCIPaddVarSOS1(scip, cons, var, weight) ); 
                     break;
                  case 2: 
                     SCIP_CALL( SCIPaddVarSOS2(scip, cons, var, weight) );
                     break;
                  default: 
                     SCIPerrorMessage("unknown SOS type: <%d>\n", type); /* should not happen */
                     SCIPABORT();
                  }
                  SCIPdebugMessage("added variable <%s> with weight %g.\n", SCIPvarGetName(var), weight);
               }
            }
         }
      }

      /* add the SOS constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMessage("(line %d) added constraint <%s>: ", lpinput->linenumber, SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** reads an LP file
 *
 *  @todo check whether variables forced to be binary for the creation of indicator constraints are
 *  really specified to be binary (or general with 0/1 bounds) in the file.
 */
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
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* parse the file */
   lpinput->section = LP_START;
   while( lpinput->section != LP_END && !hasError(lpinput) )
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


#if 0
/* prints variable name LP format conform; always use this method to stay consistent
 *
 * 1) variable names should not start with a digit
 * 2) avoid variable name starting with an 'e' or 'E' since this notation is reserved for exponential entries
 */
static
void printVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Bool             genericnames        /**< use generic variable names? */
   )
{
   const char* name;

   assert( scip != NULL );
   assert( var != NULL );

   name = SCIPvarGetName(var);
   assert( name != NULL );

   if( genericnames || name[0] == '\0' )
      SCIPinfoMessage(scip, file, "x%d", SCIPvarGetProbindex(var) + 1);
   else
   {
      if( isdigit((unsigned char)name[0]) || name[0] == 'e' || name[0] == 'E' )
         SCIPinfoMessage(scip, file, "_%s", name);
      else
         SCIPinfoMessage(scip, file, "%s", name);
   }
}
#endif

/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(scalars != NULL);
   assert(*vars != NULL);
   assert(*scalars != NULL);
   assert(nvars != NULL);
   assert(constant != NULL);

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );
         
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&(*vars)[v], &(*scalars)[v], constant) );
      }
   }
   return SCIP_OKAY;
}

/** clears the given line buffer */
static
void clearLine(
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of characters in line */
   )
{
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   (*linecnt) = 0;
   linebuffer[0] = '\0';
}

/** ends the given line with '\\0' and prints it to the given file stream */
static
void endLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of characters in line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   if( (*linecnt) > 0 )
   {
      linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** appends extension to line and prints it to the give file stream if the
 *  line exceeded the length given in the define LP_PRINTLEN */
static
void appendLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extent the line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );
   assert( strlen(linebuffer) + strlen(extension) < LP_MAX_PRINTLEN );

   /* NOTE: avoid
    *   sprintf(linebuffer, "%s%s", linebuffer, extension); 
    * because of overlapping memory areas in memcpy used in sprintf.
    */
   strncat(linebuffer, extension, LP_MAX_PRINTLEN - strlen(linebuffer));

   (*linecnt) += (int) strlen(extension);

   SCIPdebugMessage("linebuffer <%s>, length = %lu\n", linebuffer, (unsigned long)strlen(linebuffer));
   
   if( (*linecnt) > LP_PRINTLEN )
      endLine(scip, file, linebuffer, linecnt);
}


/* print row in LP format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=", "<=", or ">=") */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficient values */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nquadvarterms,      /**< number of quadratic variable terms */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int v;
   char linebuffer[LP_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   SCIP_VAR* var;
   char varname[LP_MAX_NAMELEN];
   char varname2[LP_MAX_NAMELEN];
   char consname[LP_MAX_NAMELEN + 1]; /* an extra character for ':' */
   char buffer[LP_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strcmp(type, "=") == 0 || strcmp(type, "<=") == 0 || strcmp(type, ">=") == 0 );
   assert( nlinvars == 0 || (linvars != NULL && linvals != NULL) );
   assert( nquadvarterms == 0 || quadvarterms != NULL );

   /* if there is a bilinear term, then there need to be at least two quadratic variables */
   assert( nbilinterms == 0 || (bilinterms != NULL && nquadvarterms >= 2) );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if( strlen(rowname) > 0 || strlen(rownameextension) > 0 )
   {
      (void) SCIPsnprintf(consname, LP_MAX_NAMELEN + 1, "%s%s:", rowname, rownameextension);
      appendLine(scip, file, linebuffer, &linecnt, consname);
   }

   /* print coefficients */
   for( v = 0; v < nlinvars; ++v )
   {
      var = linvars[v];
      assert( var != NULL );

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");

      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s", linvals[v], varname);
      
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }
   
   /* print quadratic part */
   if( nquadvarterms > 0 )
   {
      /* print linear coefficients of quadratic variables */
      for( v = 0; v < nquadvarterms; ++v )
      {
         if( quadvarterms[v].lincoef == 0.0 )
            continue;

         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(quadvarterms[v].var));
         (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s", quadvarterms[v].lincoef, varname);
         
         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      
      /* start quadratic part */
      appendLine(scip, file, linebuffer, &linecnt, " + [");
      
      /* print square terms */
      for( v = 0; v < nquadvarterms; ++v )
      {
         if( quadvarterms[v].sqrcoef == 0.0 )
            continue;

         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(quadvarterms[v].var));
         (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s^2", quadvarterms[v].sqrcoef, varname);
         
         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      
      /* print bilinear terms */
      for( v = 0; v < nbilinterms; ++v )
      {
         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname,  LP_MAX_NAMELEN, "%s", SCIPvarGetName(bilinterms[v].var1));
         (void) SCIPsnprintf(varname2, LP_MAX_NAMELEN, "%s", SCIPvarGetName(bilinterms[v].var2));
         (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s %s", bilinterms[v].coef, varname, varname2);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      
      /* end quadratic part */
      appendLine(scip, file, linebuffer, &linecnt, " ]");
   }

   /* print left hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %s %+.15g", type, rhs);

   /* we start a new line; therefore we tab this line */
   if( linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, " ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);
}


/** prints given (linear or) quadratic constraint information in LP format to file stream */
static
SCIP_RETCODE printQuadraticCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficients values (or NULL if all linear coefficient values are 1) */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nquadvarterms,      /**< number of quadratic variable terms */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   SCIP_VAR** activevars = NULL;
   SCIP_Real* activevals = NULL;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( rowname != NULL );

   /* The LP format does not forbid that the variable array is empty */
   assert( nlinvars == 0 || linvars != NULL );
   assert( nquadvarterms == 0 || quadvarterms != NULL );
   assert( nbilinterms == 0 || bilinterms != NULL );

   assert( lhs <= rhs );
   
   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   nactivevars = nlinvars;
   if( nlinvars > 0 ) 
   {
      /* duplicate variable and value array */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, linvars, nactivevars ) );
      if( linvals != NULL )
      {
         SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, linvals, nactivevars ) );
      }
      else
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );
         
         for( v = 0; v < nactivevars; ++v )
            activevals[v] = 1.0;
      }

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activeconstant, transformed) );
   }
   
   /* print row(s) in LP format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equal constraint */
      printRow(scip, file, rowname, "", "=", activevars, activevals, nactivevars,
         quadvarterms, nquadvarterms, bilinterms, nbilinterms,
         rhs - activeconstant);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", ">=",
            activevars, activevals, nactivevars,
            quadvarterms, nquadvarterms, bilinterms, nbilinterms,
            lhs - activeconstant);
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         printRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "<=",
            activevars, activevals, nactivevars,
            quadvarterms, nquadvarterms, bilinterms, nbilinterms,
            rhs - activeconstant);
      }
   }

   if( nlinvars > 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevars);
      SCIPfreeBufferArray(scip, &activevals);
   }
   
   return SCIP_OKAY;
}


/** prints given SOS constraint information in LP format to file stream */
static
void printSosCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            weights,            /**< array of weight values (or NULL) */
   int                   nvars,              /**< number of variables */
   int                   type                /**< SOS type (SOS1 or SOS2) */
   )
{
   int v;

   char linebuffer[LP_MAX_PRINTLEN];
   int linecnt;
   char buffer[LP_MAX_PRINTLEN];
   char varname[LP_MAX_NAMELEN];

   assert( scip != NULL );
   assert( file != NULL );
   assert( type == 1 || type == 2 );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");
   assert( strlen(rowname) < LP_MAX_NAMELEN );

   if( strlen(rowname) > 0 )
   {
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, "%s:", rowname);
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   /* SOS type */
   (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " S%d::", type);
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   for( v = 0; v < nvars; ++v )
   {
      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[v]));

      if( weights != NULL )
         (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %s:%.15g", varname, weights[v]);
      else
         (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %s:%d", varname, v);

      if(linecnt == 0 )
      {
         /* we start a new line; therefore we tab this line */
         appendLine(scip, file, linebuffer, &linecnt, " ");
      }
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   endLine(scip, file, linebuffer, &linecnt);
}

/** prints given soc constraint in LP format to file stream */
static
SCIP_RETCODE printSOCCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   SCIP_CONS*            cons                /**< second order cone constraint */
   )
{
   int v;
   char linebuffer[LP_MAX_PRINTLEN] = { '\0' };
   int linecnt;
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Real offset;
   char varname[LP_MAX_NAMELEN];
   char consname[LP_MAX_NAMELEN + 1]; /* an extra character for ':' */
   char buffer[LP_MAX_PRINTLEN];
   
   SCIP_Real rhs;

   assert( scip != NULL );
   assert( rowname != NULL );
   assert( cons != NULL );

   /* print constraint in LP format
    * the SOC constraint is given as
    *   sqrt(constant + sum_i (lhscoef_i(lhsvar_i+lhsoffset_i))^2) <= rhscoef(rhsvar+rhsoffset)
    * and is printed as
    *    sum_i (2*lhscoef_i^2 lhs_offset_i) lhsvar_i - (2 * rhscoef^2 * rhsoffset) rhsvar 
    *    + [ sum_i lhscoef_i^2 lhsvar_i^2 - rhscoef^2 rhsvar^2 ]
    *  <=
    *    - sum_i lhscoef_i^2 lhs_offset_i^2 - constant + rhscoef^2 rhsoffset^2 
    */

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if( strlen(rowname) > 0 )
   {
      (void) SCIPsnprintf(consname, LP_MAX_NAMELEN + 1, "%s:", rowname);
      appendLine(scip, file, linebuffer, &linecnt, consname);
   }
   
   rhs = -SCIPgetLhsConstantSOC(scip, cons);

   /* print linear part of left hand side and add constant parts to rhs */
   for( v = 0; v < SCIPgetNLhsVarsSOC(scip, cons); ++v )
   {
      var = SCIPgetLhsVarsSOC(scip, cons)[v];
      assert( var != NULL );
      offset = SCIPgetLhsOffsetsSOC(scip, cons)[v];
      coef = SCIPgetLhsCoefsSOC(scip, cons)[v];

      rhs -= coef * coef * offset * offset;

      if( offset == 0.0 || coef == 0.0 )
         continue;

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");
      
      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s", 2*offset*coef*coef, varname);
      
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }
   
   /* print linear part from right hand side and add constant part to rhs */
   offset = SCIPgetRhsOffsetSOC(scip, cons);
   coef = SCIPgetRhsCoefSOC(scip, cons);
   if( offset != 0.0 && coef != 0.0 )
   {
      var = SCIPgetRhsVarSOC(scip, cons);
      assert( var != NULL );
      
      rhs += coef * coef * offset * offset; 
      
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");

      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s", -2*offset*coef*coef, varname);

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   /* start quadratic part */
   appendLine(scip, file, linebuffer, &linecnt, " + [");

   /* print quadratic part of left hand side  */
   for( v = 0; v < SCIPgetNLhsVarsSOC(scip, cons); ++v )
   {
      var = SCIPgetLhsVarsSOC(scip, cons)[v];
      assert( var != NULL );
      coef = SCIPgetLhsCoefsSOC(scip, cons)[v];

      if( coef == 0.0 )
         continue;

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");
      
      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s^2", coef*coef, varname);
      
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }
   
   /* print quadratic part of right hand side  */
   coef = SCIPgetRhsCoefSOC(scip, cons);
   if( coef != 0.0 )
   {
      var = SCIPgetRhsVarSOC(scip, cons);
      assert( var != NULL );

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");
      
      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s^2", -coef*coef, varname);
      
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }
   
   /* end quadratic part */
   appendLine(scip, file, linebuffer, &linecnt, " ]");

   /* print right hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " <= %+.15g", rhs);

   /* we start a new line; therefore we tab this line */
   if( linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, " ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/** check whether given variables are aggregated and put them into an array without duplication */
static
SCIP_RETCODE collectAggregatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of mutable variables in the problem */
   SCIP_VAR**            vars,               /**< variable array */
   int*                  nAggregatedVars,    /**< number of aggregated variables on output */
   SCIP_VAR***           aggregatedVars,     /**< array storing the aggregated variables on output */
   SCIP_HASHTABLE**      varAggregated       /**< hashtable for checking duplicates */
   )
{
   int j;

   /* check variables */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VARSTATUS status;
      SCIP_VAR* var;

      var = vars[j];
      status = SCIPvarGetStatus(var);

      /* collect aggregated variables in a list */
      if( status >= SCIP_VARSTATUS_AGGREGATED )
      {
         assert( status == SCIP_VARSTATUS_AGGREGATED || 
            status == SCIP_VARSTATUS_MULTAGGR ||
            status == SCIP_VARSTATUS_NEGATED );
         
         if( !SCIPhashtableExists(*varAggregated, (void*) var) )
         {
            (*aggregatedVars)[(*nAggregatedVars)++] = var;
            SCIP_CALL( SCIPhashtableInsert(*varAggregated, (void*) var) );
         }
      }
   }

   return SCIP_OKAY;
}


/** print aggregated variable-constraints */
static
SCIP_RETCODE printAggregatedCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   int                   nvars,              /**< number of mutable variables in the problem */
   int                   nAggregatedVars,    /**< number of aggregated variables */
   SCIP_VAR**            aggregatedVars      /**< array storing the aggregated variables */
   )
{
   int j;

   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;
   char consname[LP_MAX_NAMELEN];

   assert( scip != NULL );

   /* write aggregation constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nvars) );

   for( j = 0; j < nAggregatedVars; ++j )
   {
      /* set up list to obtain substitution variables */
      nactivevars = 1;

      activevars[0] = aggregatedVars[j];
      activevals[0] = 1.0;
      activeconstant = 0.0;

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activeconstant, transformed) );

      activevals[nactivevars] = -1.0;
      activevars[nactivevars] = aggregatedVars[j];
      ++nactivevars;

      /* output constraint */
      (void) SCIPsnprintf(consname, LP_MAX_NAMELEN, "aggr_%s", SCIPvarGetName(aggregatedVars[j]));
      printRow(scip, file, consname, "", "=", activevars, activevals, nactivevars, NULL, 0, NULL, 0, - activeconstant);
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/** method check if the variable names are not longer than LP_MAX_NAMELEN */
static
void checkVarnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars               /**< number of variables */
   )
{
   int v;

   assert(scip != NULL);
   assert(vars != NULL || nvars == 0);

   /* check if the variable names are not to long */
   for( v = 0; v < nvars; ++v )
   {
      if( strlen(SCIPvarGetName(vars[v])) > LP_MAX_NAMELEN )  /*lint !e613*/
      {
         SCIPwarningMessage(scip, "there is a variable name which has to be cut down to %d characters; LP might be corrupted\n", 
            LP_MAX_NAMELEN - 1);
         return;
      }
   }
}

/** method check if the constraint names are not longer than LP_MAX_NAMELEN */
static
void checkConsnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   int c;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   assert( scip != NULL );
   assert( conss != NULL );

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL );

      /* in case the transformed is written only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_Real lhs = SCIPgetLhsLinear(scip, cons);
         SCIP_Real rhs = SCIPgetLhsLinear(scip, cons);

         if( (SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > LP_MAX_NAMELEN)
            || ( !SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > LP_MAX_NAMELEN -  4) )
         {
            SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
               LP_MAX_NAMELEN  - 1);
            return;
         }
      }
      else if( strlen(SCIPconsGetName(conss[c])) > LP_MAX_NAMELEN )
      {
         SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
            LP_MAX_NAMELEN  - 1);
         return;
      }
   }
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyLp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderLp(scip) );
 
   return SCIP_OKAY;
}


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
   SCIP_CALL( SCIPwriteLp(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );

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
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyLp) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadLp) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteLp) );

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


/** reads problem from file */
SCIP_RETCODE SCIPreadLp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   LPINPUT lpinput;
   int i;

   /* initialize LP input data */
   lpinput.file = NULL;
   lpinput.linebuf[0] = '\0';
   lpinput.probname[0] = '\0';
   lpinput.objname[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpinput.token, LP_MAX_LINELEN) ); /*lint !e506*/
   lpinput.token[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpinput.tokenbuf, LP_MAX_LINELEN) ); /*lint !e506*/
   lpinput.tokenbuf[0] = '\0';
   for( i = 0; i < LP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(lpinput.pushedtokens[i]), LP_MAX_LINELEN) );  /*lint !e866 !e506*/
   }

   lpinput.npushedtokens = 0;
   lpinput.linenumber = 0;
   lpinput.linepos = 0;
   lpinput.section = LP_START;
   lpinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   lpinput.inlazyconstraints = FALSE;
   lpinput.inusercuts = FALSE;
   lpinput.haserror = FALSE;
   lpinput.comment = FALSE;
   lpinput.endline = FALSE;

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
      return SCIP_READERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, lpinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** writes problem to file */
SCIP_RETCODE SCIPwriteLp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of mutable variables in the problem */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of general integer variables */
   int                   nimplvars,          /**< number of implicit integer variables */
   int                   ncontvars,          /**< number of continuous variables */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   int c;
   int v;

   int linecnt;
   char linebuffer[LP_MAX_PRINTLEN];

   char varname[LP_MAX_NAMELEN];
   char buffer[LP_MAX_PRINTLEN];

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLR* conshdlrInd;
   const char* conshdlrname;
   SCIP_CONS* cons;
   SCIP_CONS** consSOS1;
   SCIP_CONS** consSOS2;
   SCIP_CONS** consQuadratic;
   SCIP_CONS** consSOC;
   int nConsSOS1 = 0;
   int nConsSOS2 = 0;
   int nConsQuadratic = 0;
   int nConsSOC = 0;
   char consname[LP_MAX_NAMELEN];

   SCIP_VAR** aggregatedVars;
   int nAggregatedVars = 0;
   SCIP_HASHTABLE* varAggregated;
   SCIP_HASHMAP* consHidden;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;

   assert( scip != NULL );

   /* find indicator constraint handler */
   conshdlrInd = SCIPfindConshdlr(scip, "indicator");
   consHidden = NULL;

   /* if indicator constraint handler is present */
   if( conshdlrInd != NULL )
   {
      /* create hashtable storing linear constraints that should not be output */
      SCIP_CALL( SCIPhashmapCreate(&consHidden, SCIPblkmem(scip), SCIPcalcHashtableSize(1000)) );

      /* loop through indicator constraints (works only in transformed problem) */
      if( transformed )
      {
         SCIP_CONS** consInd;
         int nConsInd;

         consInd = SCIPconshdlrGetConss(conshdlrInd);
         nConsInd = SCIPconshdlrGetNConss(conshdlrInd);
         SCIPdebugMessage("Number of indicator constraints: %d\n", nConsInd);

         for( c = 0; c < nConsInd; ++c )
         {
            assert( consInd[c] != NULL );
            cons = SCIPgetLinearConsIndicator(consInd[c]);

            assert( !SCIPhashmapExists(consHidden, (void*) cons) );
            SCIP_CALL( SCIPhashmapSetImage(consHidden, (void*) cons, (void*) TRUE) );
            SCIPdebugMessage("Marked linear constraint <%s> as hidden.\n", SCIPconsGetName(cons));
         }
      }
      else
      {
         /* otherwise we have to pass through all constraints */
         for( c = 0; c < nconss; ++c )
         {
            cons = conss[c];
            assert( cons != NULL);

            conshdlr = SCIPconsGetHdlr(cons);
            assert( conshdlr != NULL );
            conshdlrname = SCIPconshdlrGetName(conshdlr);

            if( strcmp(conshdlrname, "indicator") == 0 )
            {
               SCIP_CONS* lincons;

               lincons = SCIPgetLinearConsIndicator(cons);
               assert( lincons != NULL );

               assert( !SCIPhashmapExists(consHidden, (void*) lincons) );
               SCIP_CALL( SCIPhashmapSetImage(consHidden, (void*) lincons, (void*) TRUE) );
               SCIPdebugMessage("Marked linear constraint <%s> as hidden.\n", SCIPconsGetName(lincons));
            }
         }
      }
   }

   /* check if the variable names are not to long */
   checkVarnames(scip, vars, nvars);
   /* check if the constraint names are to long */
   checkConsnames(scip, conss, nconss, transformed);

   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "\\ SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "\\   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "\\   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "\\   Constraints      : %d\n", nconss);
   SCIPinfoMessage(scip, file, "\\   Obj. scale       : %.15g\n", objscale);
   SCIPinfoMessage(scip, file, "\\   Obj. offset      : %.15g\n", objoffset);

   /* print objective sense */
   SCIPinfoMessage(scip, file, "%s\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "Minimize" : "Maximize");

   clearLine(linebuffer, &linecnt);
   appendLine(scip, file, linebuffer, &linecnt, " Obj:");

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];

#ifndef NDEBUG
      /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
      if( !transformed )
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );
#endif
      
      if( SCIPisZero(scip, SCIPvarGetObj(var)) )
         continue;

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, "     ");

      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %+.15g %s", SCIPvarGetObj(var), varname );

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   endLine(scip, file, linebuffer, &linecnt);

   /* print "Subject to" section */
   SCIPinfoMessage(scip, file, "Subject to\n");

   /* collect SOS, quadratic, and SOC constraints in array for later output */
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS1, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS2, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consQuadratic, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOC, nconss) );

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed is written only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));
      
      /* skip marked constraints in connection with indicator constraints */
      if( conshdlrInd != NULL && SCIPhashmapExists(consHidden, (void*) cons) )
      {
         assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0 );
         continue;
      }

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      (void) SCIPsnprintf(consname, LP_MAX_NAMELEN, "%s", SCIPconsGetName(cons));
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
               NULL, 0, NULL, 0, SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 0, NULL, 0, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 0, NULL, 0, -SCIPinfinity(scip), 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 0, NULL, 0, 1.0, SCIPinfinity(scip), transformed) );
            break;
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
               NULL, 0, NULL, 0, 1.0, SCIPinfinity(scip), transformed) );
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = (SCIP_Real)weights[v];

         SCIP_CALL( printQuadraticCons(scip, file, consname, consvars, consvals, nconsvars,
               NULL, 0, NULL, 0, -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( printQuadraticCons(scip, file, consname, consvars, consvals, 2, NULL, 0, NULL, 0, 
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "SOS1") == 0 )
      {
         /* store constraint */
         consSOS1[nConsSOS1++] = cons;
      }
      else if( strcmp(conshdlrname, "SOS2") == 0 )
      {
         /* store constraint */
         consSOS2[nConsSOS2++] = cons;
      }
      else if( strcmp(conshdlrname, "indicator") == 0 )
      {
         SCIP_CONS* lincons;
         SCIP_VAR* binvar;
         SCIP_VAR* slackvar;
         SCIP_VAR** linvars;
         SCIP_Real* linvals;
         int nlinvars;
         int cnt;
         int rhs;

         assert( conshdlrInd != NULL );

         lincons = SCIPgetLinearConsIndicator(cons);
         binvar = SCIPgetBinaryVarIndicator(cons);
         slackvar = SCIPgetSlackVarIndicator(cons);

         assert( lincons != NULL );
         assert( binvar != NULL );
         assert( slackvar != NULL );

         rhs = 1;
         if( SCIPvarGetStatus(binvar) == SCIP_VARSTATUS_NEGATED )
         {
            rhs = 0;
            binvar = SCIPvarGetNegatedVar(binvar);
         }

         /* collect linear constraint information (remove slack variable) */
         linvars = SCIPgetVarsLinear(scip, lincons);
         linvals = SCIPgetValsLinear(scip, lincons);
         nlinvars = SCIPgetNVarsLinear(scip, lincons);
         assert( linvars != NULL );
         assert( linvals != NULL );

         /* linvars always contains slack variable, thus nlinvars >= 1 */
         if( nlinvars > 1 && !SCIPconsIsDeleted(lincons) )
         {
            (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(binvar) );
            if( strlen(consname) > 0 )
               SCIPinfoMessage(scip, file, " %s: %s = %d ->", consname, varname, rhs);
            else
               SCIPinfoMessage(scip, file, " %s = %d ->", varname, rhs);

            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nlinvars-1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nlinvars-1) );

            cnt = 0;
            for( v = 0; v < nlinvars; ++v )
            {
               var = linvars[v];
               if( var != slackvar )
               {
                  consvars[cnt] = var;
                  consvals[cnt++] = linvals[v];
               }
            }
            assert( nlinvars == 0 || cnt == nlinvars-1 );

            SCIP_CALL( printQuadraticCons(scip, file, "", consvars, consvals, cnt, NULL, 0, NULL, 0, 
                  SCIPgetLhsLinear(scip, lincons), SCIPgetRhsLinear(scip, lincons), transformed) );

            SCIPfreeBufferArray(scip, &consvars);
            SCIPfreeBufferArray(scip, &consvals);
         }
      }
      else if( strcmp(conshdlrname, "quadratic") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
               SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
               SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
               SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetLhsQuadratic(scip, cons),
               SCIPgetRhsQuadratic(scip, cons), transformed) );

         consQuadratic[nConsQuadratic++] = cons;
      }
      else if( strcmp(conshdlrname, "soc") == 0 )
      {
         SCIP_CALL( printSOCCons(scip, file, consname, cons) );

         consSOC[nConsSOC++] = cons;
      }
      else
      {
         SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "\\ ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
         SCIPinfoMessage(scip, file, ";\n");
      }
   }

   /* create hashtable for storing aggregated variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &aggregatedVars, nvars) );
   SCIP_CALL( SCIPhashtableCreate(&varAggregated, SCIPblkmem(scip), 1000, hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );

   /* check for aggregated variables in SOS1 constraints and output aggregations as linear constraints */
   for( c = 0; c < nConsSOS1; ++c )
   {
      cons = consSOS1[c];
      consvars = SCIPgetVarsSOS1(scip, cons);
      nconsvars = SCIPgetNVarsSOS1(scip, cons);

      SCIP_CALL( collectAggregatedVars(scip, nconsvars, consvars, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }

   /* check for aggregated variables in SOS2 constraints and output aggregations as linear constraints */
   for( c = 0; c < nConsSOS2; ++c )
   {
      cons = consSOS2[c];
      consvars = SCIPgetVarsSOS2(scip, cons);
      nconsvars = SCIPgetNVarsSOS2(scip, cons);

      SCIP_CALL( collectAggregatedVars(scip, nconsvars, consvars, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }
   
   /* check for aggregated variables in quadratic parts of quadratic constraints and output aggregations as linear constraints */
   for( c = 0; c < nConsQuadratic; ++c )
   {
      cons = consQuadratic[c];
      for( v = 0; v < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++v )
      {
         SCIP_CALL( collectAggregatedVars(scip, 1, &SCIPgetQuadVarTermsQuadratic(scip, cons)[v].var,
               &nAggregatedVars, &aggregatedVars, &varAggregated) );
      }         
   }

   /* check for aggregated variables in second order cone constraints and output aggregations as linear constraints */
   for( c = 0; c < nConsSOC; ++c )
   {
      cons = consSOC[c];
      
      SCIP_CALL( collectAggregatedVars(scip, SCIPgetNLhsVarsSOC(scip, cons), SCIPgetLhsVarsSOC(scip, cons),
            &nAggregatedVars, &aggregatedVars, &varAggregated) );
      var = SCIPgetRhsVarSOC(scip, cons);
      SCIP_CALL( collectAggregatedVars(scip, 1, &var, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }

   /* print aggregation constraints */
   SCIP_CALL( printAggregatedCons(scip, file, transformed, nvars, nAggregatedVars, aggregatedVars) );

   /* print "Bounds" section */
   SCIPinfoMessage(scip, file, "Bounds\n");
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );
      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted which are valid in the current node */
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }
      else
      {
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
      }

      if( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
         SCIPinfoMessage(scip, file, " %s free\n", varname);
      else
      {
         /* print lower bound */
         if( SCIPisInfinity(scip, -lb) )
            SCIPinfoMessage(scip, file, " -inf <= ");
         else
         {
            if( SCIPisZero(scip, lb) )
            {
               /* variables are nonnegative by default - so we skip these variables */
               if( SCIPisInfinity(scip, ub) )
                  continue;
               lb = 0.0;
            }

            SCIPinfoMessage(scip, file, " %.15g <= ", lb);
         }
         /* print variable name */
         SCIPinfoMessage(scip, file, "%s", varname);

         /* print upper bound as far this one is not infinity */
         if( !SCIPisInfinity(scip, ub) )
            SCIPinfoMessage(scip, file, " <= %.15g", ub);

         SCIPinfoMessage(scip, file, "\n");
      }
   }

   /* output aggregated variables as 'free' */
   for( v = 0; v < nAggregatedVars; ++v )
   {
      var = aggregatedVars[v];
      assert( var != NULL );
      (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      SCIPinfoMessage(scip, file, " %s free\n", varname);
   }

   /* free space */
   SCIPfreeBufferArray(scip, &aggregatedVars);
   SCIPhashtableFree(&varAggregated);
   if( conshdlrInd != NULL )
      SCIPhashmapFree(&consHidden);

   /* print binaries section */
   if( nbinvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Binaries\n");

      clearLine(linebuffer, &linecnt);

      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         assert( var != NULL );

         if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         {
            (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );
            (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %s", varname);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
      }

      endLine(scip, file, linebuffer, &linecnt);
   }

   /* print generals section */
   if( nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Generals\n");

      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         assert( var != NULL );

         if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
         {
            (void) SCIPsnprintf(varname, LP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );
            (void) SCIPsnprintf(buffer, LP_MAX_PRINTLEN, " %s", varname);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
      }
      endLine(scip, file, linebuffer, &linecnt);
   }

   /* print SOS section */
   if( nConsSOS1 > 0 || nConsSOS2 > 0 )
   {
      SCIP_Real* weights;
      SCIPinfoMessage(scip, file, "SOS\n");

      /* first output SOS1 constraints */
      for( c = 0; c < nConsSOS1; ++c )
      {
         cons = consSOS1[c];
         consvars = SCIPgetVarsSOS1(scip, cons);
         nconsvars = SCIPgetNVarsSOS1(scip, cons);
         weights = SCIPgetWeightsSOS1(scip, cons);

         (void) SCIPsnprintf(consname, LP_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
         printSosCons(scip, file, consname, consvars, weights, nconsvars, 1);
      }

      /* next output SOS2 constraints */
      for( c = 0; c < nConsSOS2; ++c )
      {
         cons = consSOS2[c];
         consvars = SCIPgetVarsSOS2(scip, cons);
         nconsvars = SCIPgetNVarsSOS2(scip, cons);
         weights = SCIPgetWeightsSOS2(scip, cons);

         (void) SCIPsnprintf(consname, LP_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
         printSosCons(scip, file, consname, consvars, weights, nconsvars, 2);
      }
   }

   /* free space */
   SCIPfreeBufferArray(scip, &consSOS1);
   SCIPfreeBufferArray(scip, &consSOS2);
   SCIPfreeBufferArray(scip, &consQuadratic);
   SCIPfreeBufferArray(scip, &consSOC);

   /* end of lp format */
   SCIPinfoMessage(scip, file, "%s\n", "End");

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
