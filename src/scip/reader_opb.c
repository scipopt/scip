/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_opb.c,v 1.66 2010/10/28 17:31:23 bzfwinkm Exp $"

/**@file   reader_opb.c
 * @ingroup FILEREADERS 
 * @brief  pseudo-Boolean file reader (opb format)
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/* http://www.cril.univ-artois.fr/PB07/solver_req.html
 * http://www.cril.univ-artois.fr/PB10/format.pdf
 *
 * The syntax of the input file format can be described by a simple Backus-Naur form. <formula> is the start symbol of
 * this grammar.
 *
 *  <formula>::= <sequence_of_comments> 
 *               [<objective>] | [<softheader>]
 *               <sequence_of_comments_or_constraints>
 *
 *  <sequence_of_comments>::= <comment> [<sequence_of_comments>]
 *  <comment>::= "*" <any_sequence_of_characters_other_than_EOL> <EOL>
 *  <sequence_of_comments_or_constraints>::=<comment_or_constraint> [<sequence_of_comments_or_constraints>]
 *  <comment_or_constraint>::=<comment>|<constraint> 
 *
 *  <objective>::= "min:" <zeroOrMoreSpace> <sum>  ";"
 *  <constraint>::= <sum> <relational_operator> <zeroOrMoreSpace> <integer> <zeroOrMoreSpace> ";"
 *  
 *  <sum>::= <weightedterm> | <weightedterm> <sum>
 *  <weightedterm>::= <integer> <oneOrMoreSpace> <term> <oneOrMoreSpace>
 *  
 *  <integer>::= <unsigned_integer> | "+" <unsigned_integer> | "-" <unsigned_integer>
 *  <unsigned_integer>::= <digit> | <digit><unsigned_integer>
 *  
 *  <relational_operator>::= ">=" | "="
 *  
 *  <variablename>::= "x" <unsigned_integer>
 *  
 *  <oneOrMoreSpace>::= " " [<oneOrMoreSpace>]
 *  <zeroOrMoreSpace>::= [" " <zeroOrMoreSpace>]
 *  
 *  For linear pseudo-Boolean instances, <term> is defined as
 *  
 *  <term>::=<variablename>
 *  
 *  For non-linear instances, <term> is defined as
 *  
 *  <term>::= <oneOrMoreLiterals>
 *  <oneOrMoreLiterals>::= <literal> | <literal> <oneOrMoreSpace> <oneOrMoreLiterals>
 *  <literal>::= <variablename> | "~"<variablename>
 *  
 * For wbo-files are the following additional/changed things possible.
 *  
 *  <softheader>::= "soft:" [<unsigned integer>] ";"
 *  
 *  <comment_or_constraint>::=<comment>|<constraint>|<softconstraint> 
 *
 *  <softconstraint>::= "[" <zeroOrMoreSpace> <unsigned integer> <zeroOrMoreSpace> "]" <contraint>
 *  
 */

/**@note Our parser should also be lax by handling variable names and it's possible to read doubles instead of integer
 *       and possible some more :). 
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

#include "scip/reader_opb.h"
#include "scip/cons_and.h"
#include "scip/cons_indicator.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/pub_misc.h"

#define READER_NAME             "opbreader"
#define READER_DESC             "file reader for pseudo-Boolean problem in opb format"
#define READER_EXTENSION        "opb"

#define USEINDICATOR            TRUE
#define BACKIMPLICATION         FALSE
#define GENCONSNAMES                      /* remove if no constraint names should be generated */


/*
 * Data structures
 */

#define HASHSIZE_OPBANDCONS   131101 /**< minimal size of hash table in and constraint tables */
#define OPB_MAX_LINELEN       65536  /**< size of the line buffer for reading or writing */
#define OPB_MAX_PUSHEDTOKENS  2
#define OPB_INIT_COEFSSIZE    8192


/** Section in OPB File */
enum OpbExpType 
{
   OPB_EXP_NONE,
   OPB_EXP_UNSIGNED, 
   OPB_EXP_SIGNED
};
typedef enum OpbExpType OPBEXPTYPE;

enum OpbSense 
{
   OPB_SENSE_NOTHING,
   OPB_SENSE_LE,
   OPB_SENSE_GE,
   OPB_SENSE_EQ
};
typedef enum OpbSense OPBSENSE;

/* struct for reading an opb files, used to find fast whether an and-constraint is new or already existing */
struct ConsAndData
{
   SCIP_VAR**            vars;
   SCIP_VAR*             resultant;
   int                   nvars;
};

typedef struct ConsAndData CONSANDDATA;

/** OPB reading data */
struct OpbInput
{
   SCIP_FILE*            file;
   char                  linebuf[OPB_MAX_LINELEN+1];
   char*                 token;
   char*                 tokenbuf;
   char*                 pushedtokens[OPB_MAX_PUSHEDTOKENS];
   int                   npushedtokens;
   int                   linenumber;
   int                   linepos;
   int                   bufpos;
   SCIP_OBJSENSE         objsense;
   SCIP_Bool             comment;
   SCIP_Bool             endline;
   SCIP_Bool             eof;
   SCIP_Bool             haserror;
   CONSANDDATA**         consanddata;
   int                   nconsanddata;
   int                   sconsanddata;
   int                   maxvarsperand;
   int                   nproblemcoeffs;
   SCIP_HASHTABLE*       hashtable;
   int                   hashtablesize;
   SCIP_Bool             wbo;
   SCIP_Real             topcost;
   int                   nindvars;
#ifdef GENCONSNAMES
   int                   consnumber;
#endif
};

typedef struct OpbInput OPBINPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=;[]";
static const char commentchars[] = "*";

/*
 * Local methods (for reading)
 */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyOpbAndcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqOpbAndcons)
{
   SCIP* scip;
   CONSANDDATA* consanddata1;
   CONSANDDATA* consanddata2;
   int i;

   consanddata1 = (CONSANDDATA*)key1;
   consanddata2 = (CONSANDDATA*)key2;
   scip = (SCIP*)userptr; 

   assert(scip != NULL);
   assert(consanddata1 != NULL);
   assert(consanddata2 != NULL);
   assert(consanddata1->vars != NULL);
   assert(consanddata1->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      int v;
      for( v = consanddata1->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(consanddata1->vars[v]) >= SCIPvarGetIndex(consanddata1->vars[v - 1]));
   }
#endif
   assert(consanddata2->vars != NULL);
   assert(consanddata2->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      int v;
      for( v = consanddata2->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(consanddata2->vars[v]) >= SCIPvarGetIndex(consanddata2->vars[v - 1]));
   }
#endif

   /* checks trivial case */
   if( consanddata1->nvars != consanddata2->nvars )
      return FALSE;

   for( i = 0; i < consanddata1->nvars ; ++i )
   {
      /* tests if variables are equal */
      if( consanddata1->vars[i] != consanddata2->vars[i] )
      {
         assert(SCIPvarCompare(consanddata1->vars[i], consanddata2->vars[i]) == 1 || 
            SCIPvarCompare(consanddata1->vars[i], consanddata2->vars[i]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompare(consanddata1->vars[i], consanddata2->vars[i]) == 0); 
   } 
   
   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValOpbAndcons)
{  /*lint --e{715}*/
   CONSANDDATA* consanddata;
   int minidx;
   int mididx;
   int maxidx;
   
   consanddata = (CONSANDDATA*)key;
   
   assert(consanddata != NULL);
   assert(consanddata->vars != NULL);
   assert(consanddata->nvars > 1);
#ifndef NDEBUG
   {
      /* check that these variables are sorted */
      int v;
      for( v = consanddata->nvars - 1; v > 0; --v )
         assert(SCIPvarGetIndex(consanddata->vars[v]) >= SCIPvarGetIndex(consanddata->vars[v - 1]));
   }
#endif

   minidx = SCIPvarGetIndex(consanddata->vars[0]);
   mididx = SCIPvarGetIndex(consanddata->vars[consanddata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consanddata->vars[consanddata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   return (consanddata->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx; /*lint !e701*/
}

/** issues an error message and marks the OPB data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           msg                 /**< error message */
   )
{
#if 0
   char formatstr[256];
#endif
   assert(opbinput != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error in line %d: %s found <%s>\n",
      opbinput->linenumber, msg, opbinput->token);
   if( opbinput->linebuf[strlen(opbinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", opbinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", opbinput->linebuf);
   }

#if 0
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", opbinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, formatstr, "^");
#endif

   opbinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   assert(opbinput != NULL);

   return opbinput->haserror;
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
   OPBEXPTYPE*           exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit(c) )
      return TRUE;
   else if( (*exptype == OPB_EXP_NONE) && !(*hasdot) && (c == '.') )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == OPB_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = OPB_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit(nextc) )
      {
         *exptype = OPB_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == OPB_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = OPB_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   int i;
   char* last;
  
   assert(opbinput != NULL);

   /* if we previously detected a comment we have to parse the remaining line away if there is something left */
   if( !opbinput->endline && opbinput->comment )
   { 
      SCIPdebugMessage("Throwing rest of comment away.\n");
   
      do
      {
         opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';
         (void)SCIPfgets(opbinput->linebuf, sizeof(opbinput->linebuf), opbinput->file);
      }
      while( opbinput->linebuf[OPB_MAX_LINELEN-2] != '\0' );
         
      opbinput->comment = FALSE;
      opbinput->endline = TRUE;
   }

   /* clear the line */
   BMSclearMemoryArray(opbinput->linebuf, OPB_MAX_LINELEN);
   opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';
   
   /* set line position */
   if( opbinput->endline )
   {
      opbinput->linepos = 0;
      opbinput->linenumber++;
   }
   else
      opbinput->linepos += OPB_MAX_LINELEN - 2;
   
   if( SCIPfgets(opbinput->linebuf, sizeof(opbinput->linebuf), opbinput->file) == NULL )
      return FALSE;
   
   opbinput->bufpos = 0;
      
   if( opbinput->linebuf[OPB_MAX_LINELEN-2] != '\0' )
   {
      /* overwrite the character to search the last blank from this position backwards */
      opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';

      /* buffer is full; erase last token since it might be incomplete */
      opbinput->endline = FALSE;
      last = strrchr(opbinput->linebuf, ' ');

      if( last == NULL )
      {
	 SCIPwarningMessage("we read %d character from the file; these might indicates an corrupted input file!", 
            OPB_MAX_LINELEN - 2);
	 opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';
	 SCIPdebugMessage("the buffer might be corrupted\n");
      }
      else 
      {
	 SCIPfseek(opbinput->file, -(long) strlen(last) - 1, SEEK_CUR);
	 SCIPdebugMessage("correct buffer, reread the last %ld characters\n", (long) strlen(last) + 1);
	 *last = '\0';
      }
   }
   else 
   {
      /* found end of line */
      opbinput->endline = TRUE;
   }
   
   opbinput->linebuf[OPB_MAX_LINELEN-1] = '\0';
   opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   opbinput->comment = FALSE;

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(opbinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
         opbinput->comment = TRUE;
         break;
      }
   }

   SCIPdebugMessage("%s\n", opbinput->linebuf);

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
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   SCIP_Bool hasdot;
   OPBEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(opbinput != NULL);
   assert(opbinput->bufpos < OPB_MAX_LINELEN);

   /* check the token stack */
   if( opbinput->npushedtokens > 0 )
   {
      swapPointers(&opbinput->token, &opbinput->pushedtokens[opbinput->npushedtokens-1]);
      opbinput->npushedtokens--;
      SCIPdebugMessage("(line %d) read token again: '%s'\n", opbinput->linenumber, opbinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = opbinput->linebuf;
   while( isDelimChar(buf[opbinput->bufpos]) )
   {
      if( buf[opbinput->bufpos] == '\0' )
      {
         if( !getNextLine(opbinput) )
         {
            SCIPdebugMessage("(line %d) end of file\n", opbinput->linenumber);
            return FALSE;
         }
         assert(opbinput->bufpos == 0);
      }
      else
      {
         opbinput->bufpos++;
         opbinput->linepos++;
      }
   }
   assert(opbinput->bufpos < OPB_MAX_LINELEN);
   assert(!isDelimChar(buf[opbinput->bufpos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = OPB_EXP_NONE;
   if( isValueChar(buf[opbinput->bufpos], buf[opbinput->bufpos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < OPB_MAX_LINELEN);
         assert(!isDelimChar(buf[opbinput->bufpos]));
         opbinput->token[tokenlen] = buf[opbinput->bufpos];
         tokenlen++;
         opbinput->bufpos++;
         opbinput->linepos++;
      }
      while( isValueChar(buf[opbinput->bufpos], buf[opbinput->bufpos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < OPB_MAX_LINELEN);
         opbinput->token[tokenlen] = buf[opbinput->bufpos];
         tokenlen++;
         opbinput->bufpos++;
         opbinput->linepos++;
         if( tokenlen == 1 && isTokenChar(opbinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[opbinput->bufpos]) && !isTokenChar(buf[opbinput->bufpos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', 
       * replace the token by the inequality sense
       */
      if( tokenlen >= 1
         && (opbinput->token[tokenlen-1] == '<' || opbinput->token[tokenlen-1] == '>' || opbinput->token[tokenlen-1] == '=')
         && buf[opbinput->bufpos] == '=' )
      {
         opbinput->bufpos++;
         opbinput->linepos++;
      }
      else if( opbinput->token[tokenlen-1] == '=' && (buf[opbinput->bufpos] == '<' || buf[opbinput->bufpos] == '>') )
      {
         opbinput->token[tokenlen-1] = buf[opbinput->bufpos];
         opbinput->bufpos++;
         opbinput->linepos++;
      }
   }
   assert(tokenlen < OPB_MAX_LINELEN);
   opbinput->token[tokenlen] = '\0';

   SCIPdebugMessage("(line %d) read token: '%s'\n", opbinput->linenumber, opbinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   assert(opbinput != NULL);
   assert(opbinput->npushedtokens < OPB_MAX_PUSHEDTOKENS);

   swapPointers(&opbinput->pushedtokens[opbinput->npushedtokens], &opbinput->token);
   opbinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   assert(opbinput != NULL);
   assert(opbinput->npushedtokens < OPB_MAX_PUSHEDTOKENS);

   swapPointers(&opbinput->pushedtokens[opbinput->npushedtokens], &opbinput->tokenbuf);
   opbinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   assert(opbinput != NULL);

   swapPointers(&opbinput->token, &opbinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isEndLine(
   OPBINPUT*              opbinput             /**< OPB reading data */
   )
{
   assert(opbinput != NULL);
   
   if( *(opbinput->token) ==  ';')
      return TRUE;
   
   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   OPBINPUT*             opbinput,           /**< OPB reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(opbinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( strlen(opbinput->token) == 1 && opbinput->token[1] == '\0' )
   {
      if( *opbinput->token == '+' )
         return TRUE;
      else if( *opbinput->token == '-' )
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
   OPBINPUT*             opbinput,           /**< OPB reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(opbinput != NULL);
   assert(value != NULL);

   if( strcasecmp(opbinput->token, "INFINITY") == 0 || strcasecmp(opbinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(opbinput->token, &endptr);
      if( endptr != opbinput->token && *endptr == '\0' )
      {
         *value = val;
         if( strlen(opbinput->token) > 18 )
            opbinput->nproblemcoeffs++;
         return TRUE;
      }
   }
   
   return FALSE;
}

/** returns whether the current token is an equation sense */
static
SCIP_Bool isSense(
   OPBINPUT*              opbinput,            /**< OPB reading data */
   OPBSENSE*              sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(opbinput != NULL);

   if( strcmp(opbinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = OPB_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(opbinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = OPB_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(opbinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = OPB_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isStartingSoftConstraintWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(scip != NULL);
   assert(opbinput != NULL);

   if( strcmp(opbinput->token, "[") == 0 )
      return TRUE;

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isEndingSoftConstraintWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(scip != NULL);
   assert(opbinput != NULL);

   if( strcmp(opbinput->token, "]") == 0 )
      return TRUE;

   return FALSE;
}

/** create binary variable with given name */
static
SCIP_RETCODE createVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   char*                 name                /**< name for the variable */
   )
{   
   SCIP_VAR* newvar;
   SCIP_Bool dynamiccols;
   SCIP_Bool initial;
   SCIP_Bool removable;
   
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/dynamiccols", &dynamiccols) );
   initial = !dynamiccols;
   removable = dynamiccols;
   
   /* create new variable of the given name */
   SCIPdebugMessage("creating new variable: <%s>\n", name);
   
   SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
         initial, removable, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, newvar) );
   *var = newvar;
   
   /* because the variable was added to the problem, it is captured by SCIP and we
    * can safely release it right now without making the returned *var invalid */
   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

   return SCIP_OKAY;
}

/** returns the variable with the given name, or creates a new variable if it does not exist */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Bool*            isAndResultant      /**< pointer to store if we found an AND */
   )
{
   SCIP_Bool negated;
   SCIP_Bool created = FALSE;
   char* name;

   CONSANDDATA* newdata;
   int svars;

   assert(var != NULL);
   
   SCIP_CALL( SCIPallocBlockMemory(scip, &newdata) );

   newdata->nvars = 0;
   svars = opbinput->maxvarsperand;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newdata->vars), svars) );
   
   name = opbinput->token; 
   assert(name != NULL);
   
   *isAndResultant = FALSE;

   /* parse AND terms */
   while(!isdigit( *name ) && !isTokenChar(*name) && !opbinput->haserror )
   {
      negated = FALSE;
      if( *name == '~' )
      {
         negated = TRUE;
         ++name;
      }

      *var = SCIPfindVar(scip, name);
      if( *var == NULL )
      {
         SCIP_CALL( createVariable(scip, var, name) );
         created = TRUE;
      }
      
      if( negated )
      {
         SCIP_VAR* negvar;
         SCIP_CALL( SCIPgetNegatedVar(scip, *var, &negvar) );
         
         *var = negvar;
      }
      
      /* reallocated memory */
      if( newdata->nvars == svars )
      {
         svars *= 2;
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(newdata->vars), svars) );
      }
      
      newdata->vars[newdata->nvars] = *var;
      ++(newdata->nvars);
      
      if( !getNextToken(opbinput) )
         opbinput->haserror = TRUE;

      name = opbinput->token;
   }
   
   /* check if we found at least on variable */
   if( newdata->nvars == 0 )
      syntaxError(scip, opbinput, "expected a variable name");

   pushToken(opbinput);
   
   if( newdata->nvars > 1 )
   {
      *isAndResultant = TRUE;

      /* sort vars to be able to check easily whether this andconstraint exists again */
      SCIPsortPtr((void**)(newdata->vars), SCIPvarComp, newdata->nvars);
      
      if (!created)
      {
	CONSANDDATA* tmpdata;

	/* get constraint from current hash table with same variables as cons0 */
	tmpdata = (CONSANDDATA*)(SCIPhashtableRetrieve(opbinput->hashtable, (void*)newdata));
	if( tmpdata != NULL )
	  *var = tmpdata->resultant;
	else 
	  *var = NULL;
      }
      
      if( created || *var == NULL )
      {
         SCIP_Bool separate;
         SCIP_Bool propagate;
         SCIP_Bool removable;
         
         SCIP_CONS* cons;
         char varname[SCIP_MAXSTRLEN];
         
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "andresultant%d", opbinput->nconsanddata);
         SCIP_CALL( createVariable(scip, var, varname) );
         assert( var != NULL );
        
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/nlcseparate", &separate) );
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/nlcpropagate", &propagate) );
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/nlcremovable", &removable) );
         
         SCIP_CALL( SCIPcreateConsAnd(scip, &cons, "", *var, newdata->nvars, newdata->vars,
               TRUE, separate, TRUE, TRUE, propagate, FALSE, FALSE, FALSE, removable, FALSE) );
         
	 newdata->resultant = *var;

         if(opbinput->nconsanddata == opbinput->sconsanddata)
         {
            opbinput->sconsanddata *= 2;
            SCIP_CALL( SCIPreallocMemoryArray(scip, &(opbinput->consanddata), opbinput->sconsanddata) );
         }
         
         opbinput->consanddata[opbinput->nconsanddata] = newdata;
         ++(opbinput->nconsanddata);
         /* no such constraint in current hash table: insert cons0 into hash table */  
         SCIP_CALL( SCIPhashtableInsert(opbinput->hashtable, (void*) newdata) );

         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
         
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

	 return SCIP_OKAY;
      }
   }
   
   SCIPfreeMemoryArray(scip, &(newdata->vars));
   SCIPfreeBlockMemory(scip, &newdata);

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   OPB_MAX_LINELEN */
   SCIP_VAR***           vars,               /**< pointer to store the array with variables (must be freed by caller) */
   SCIP_Real**           coefs,              /**< pointer to store the array with coefficients (must be freed by caller) */
   int*                  ncoefs,             /**< pointer to store the number of coefficients */
   SCIP_Bool*            newsection,         /**< pointer to store whether a new section was encountered */
   SCIP_Bool*            isNonlinear,        /**< pointer to store if we have an nonlinear constraint */
   SCIP_Bool*            issoftcons,         /**< pointer to store whether it is a soft constraint (for wbo files) */
   SCIP_Real*            weight              /**< pointer to store the weight of the softconstraint */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Bool haveweightstart;
   SCIP_Bool haveweightend;
   SCIP_Bool isAndResultant;
   SCIP_Real coef;
   int coefsign;
   int coefssize;

   assert(opbinput != NULL);
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
   *isNonlinear = FALSE;
   *issoftcons = FALSE;

   SCIPdebugMessage("read coefficients\n");

   /* read the first token, which may be the name of the line */
   if( getNextToken(opbinput) )
   {
      /* remember the token in the token buffer */
      swapTokenBuffer(opbinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(opbinput) )
      {
         if( strcmp(opbinput->token, ":") == 0 )
         {
            /* the second token was a colon ':' the first token is a constraint name */
            (void)strncpy(name, opbinput->tokenbuf, SCIP_MAXSTRLEN);
            name[SCIP_MAXSTRLEN-1] = '\0';
            SCIPdebugMessage("(line %d) read constraint name: '%s'\n", opbinput->linenumber, name);

            /* all but the first coefficient need a sign */
            if( strcmp(name, "soft") == 0 && (SCIPgetNVars(scip) > 0 || SCIPgetNConss(scip) > 0) )
            {
               syntaxError(scip, opbinput, "Soft top cost line needs to be the first non-comment line, and without any objective function.\n");
               return SCIP_OKAY;
            }
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            SCIPdebugMessage("(line %d) constraint has no name\n", opbinput->linenumber);
            pushToken(opbinput);
            pushBufferToken(opbinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(opbinput);
      }
   }
   else
   {
      assert(SCIPfeof( opbinput->file ) );
      opbinput->eof = TRUE;
      return SCIP_OKAY;
   }

   /* initialize buffers for storing the coefficients */
   coefssize = OPB_INIT_COEFSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, vars, coefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, coefs, coefssize) );

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   haveweightstart = FALSE;
   haveweightend = FALSE;
   *ncoefs = 0;
   while( getNextToken(opbinput) && !hasError(opbinput) )
   {
      SCIP_VAR* var;
      
      if( isEndLine(opbinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* check if we reached an equation sense */
      if( isSense(opbinput, NULL) )
      {
         /* put the sense back onto the token stack */
         pushToken(opbinput);
         break;
      }

      /* check if we read a sign */
      if( isSign(opbinput, &coefsign) )
      {
         SCIPdebugMessage("(line %d) read coefficient sign: %+d\n", opbinput->linenumber, coefsign);
         havesign = TRUE;
         continue;
      }

      /* check if we read a value */
      if( isValue(scip, opbinput, &coef) )
      {
         /* all but the first coefficient need a sign */
         if( *ncoefs > 0 && !havesign )
         {
            syntaxError(scip, opbinput, "expected sign ('+' or '-') or sense ('<' or '>')");
            return SCIP_OKAY;
         }

         SCIPdebugMessage("(line %d) read coefficient value: %g with sign %+d\n", opbinput->linenumber, coef, coefsign);
         if( havevalue )
         {
            syntaxError(scip, opbinput, "two consecutive values");
            return SCIP_OKAY;
         }
         havevalue = TRUE;

         /* if we read a wbo file */
         if( strcmp(name, "soft") == 0 )
         {
            assert(*ncoefs == 0);
            
            (*coefs)[*ncoefs] = coefsign * coef;
            ++(*ncoefs);
         }

         continue;
      }

      /* check if we are reading a soft constraint line */
      if( *ncoefs == 0 && !havesign && !havevalue && strcmp(name, "soft") != 0 && isStartingSoftConstraintWeight(scip, opbinput) )
      {
         if( !opbinput->wbo )
         {
            SCIPwarningMessage("Found in line %d a soft constraint, without having read a starting top-cost line.\n", opbinput->linenumber);
         }
         haveweightstart = TRUE;

         continue;
      }
      if( *ncoefs == 0 && havevalue && haveweightstart && isEndingSoftConstraintWeight(scip, opbinput) )
      {
         *weight = coefsign * coef;
#if (USEINDICATOR == TRUE && BACKIMPLICATION == FALSE)
         assert(*weight > 0);
#endif
         SCIPdebugMessage("(line %d) found soft constraint weight: %g\n", opbinput->linenumber, *weight);

         coefsign = +1;
         havesign = FALSE;
         havevalue = FALSE;
         haveweightend = TRUE;
         *issoftcons = TRUE;

         continue;
      }

      /* if we read a '[' we should already read a ']', which indicates that we read a soft constraint, 
       * we have a parsing error */
      if( haveweightstart != haveweightend )
      {
         syntaxError(scip, opbinput, "Wrong soft constraint.");
         return SCIP_OKAY;
      }

      /* if we read the first non-comment line of a wbo file we should never be here */
      if( strcmp(name, "soft") == 0 )
      {
         syntaxError(scip, opbinput, "Wrong soft top cost line.");
         return SCIP_OKAY;
      }

      /* the token is a variable name: get the corresponding variable (or create a new one) */
      SCIP_CALL( getVariable(scip, opbinput, &var, &isAndResultant) );
      
      if( isAndResultant )
         *isNonlinear = TRUE;

      /* insert the coefficient */
      SCIPdebugMessage("(line %d) found term: %+g<%s>\n", opbinput->linenumber, coefsign * coef, SCIPvarGetName(var));
      if( !SCIPisZero(scip, coef) )
      {
         /* resize the vars and coefs array if needed */
         if( *ncoefs >= coefssize )
         {
            coefssize *= 2;
            coefssize = MAX(coefssize, (*ncoefs)+1);
            SCIP_CALL( SCIPreallocBufferArray(scip, vars, coefssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, coefs, coefssize) );
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

/** set the objective section */
static
SCIP_RETCODE setObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           sense,              /**< objective sense */ 
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            coefs,              /**< array of objective values */      
   int                   ncoefs              /**< number of coefficients */ 
   )
{
   assert(opbinput != NULL);
   assert( isEndLine(opbinput) );

   if( strcmp(sense, "max" ) == 0 )
      opbinput->objsense = SCIP_OBJSENSE_MAXIMIZE;
  
   if( !hasError(opbinput) )
   {
      int i;
      
      /* set the objective values */
      for( i = 0; i < ncoefs; ++i )
      {
	SCIP_CALL( SCIPchgVarObj(scip, vars[i], SCIPvarGetObj(vars[i]) + coefs[i]) );
      }
   }

   return SCIP_OKAY;
}

/** reads the constraints section */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   int*                  nNonlinearConss     /**< pointer to store number of nonlinear constraints */
   )
{
   char name[OPB_MAX_LINELEN];
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Bool newsection;
   OPBSENSE sense;
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
   SCIP_Bool isNonlinear;
   int ncoefs;
   int sidesign;
   SCIP_Bool issoftcons;
   SCIP_Real weight;

   assert(scip != NULL);
   assert(opbinput != NULL);
   assert(nNonlinearConss != NULL);
   
   weight = -SCIPinfinity(scip);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, opbinput, name, &vars, &coefs, &ncoefs, &newsection, &isNonlinear, &issoftcons, &weight) );

   if( hasError(opbinput) || opbinput->eof )
      goto TERMINATE;
   if( newsection )
   {
      if( strcmp(name, "min") == 0 || strcmp(name, "max") == 0 )
      {
         if( opbinput->wbo )
         {
            syntaxError(scip, opbinput, "Cannot have an objective function when having soft constraints.\n");
            goto TERMINATE;
         }

         /* set objective function  */
         SCIP_CALL( setObjective(scip, opbinput, name, vars, coefs, ncoefs) );
      }
      else if( strcmp(name, "soft") == 0 )
      {
         /* we have a "weighted boolean optimization"-file(wbo) */
         opbinput->wbo = TRUE;
         if( ncoefs == 0 )
            opbinput->topcost = SCIPinfinity(scip);
         else
         {
            assert(ncoefs == 1);
            opbinput->topcost = coefs[0];
         } 
         SCIPdebugMessage("Weighted Boolean Optimization problem has topcost of %g\n", opbinput->topcost);
      }
      else if( ncoefs > 0 )
         syntaxError(scip, opbinput, "expected constraint sense '=' or '>='");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if( !getNextToken(opbinput) || !isSense(opbinput, &sense) )
   {
      syntaxError(scip, opbinput, "expected constraint sense '=' or '>='");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if( !getNextToken(opbinput) )
   {
      syntaxError(scip, opbinput, "missing right hand side");
      goto TERMINATE;
   }
   if( isSign(opbinput, &sidesign) )
   {
      if( !getNextToken(opbinput) )
      {
         syntaxError(scip, opbinput, "missing value of right hand side");
         goto TERMINATE;
      }
   }
   if( !isValue(scip, opbinput, &sidevalue) )
   {
      syntaxError(scip, opbinput, "expected value as right hand side");
      goto TERMINATE;
   }
   sidevalue *= sidesign;

   /* check if we reached the line end */
   if( !getNextToken(opbinput) || !isEndLine(opbinput) )
   {
      /*
       *(opbinput->token) = '\0';
       *(opbinput->tokenbuf) = '\0';
       */
      syntaxError(scip, opbinput, "expected endline character ';'");
      goto TERMINATE;
   }

   /* assign the left and right hand side, depending on the constraint sense */
   switch( sense )
   {
   case OPB_SENSE_GE:
      lhs = sidevalue;
      rhs = SCIPinfinity(scip);
      break;
   case OPB_SENSE_LE:
      lhs = -SCIPinfinity(scip);
      rhs = sidevalue;
      break;
   case OPB_SENSE_EQ:
      lhs = sidevalue;
      rhs = sidevalue;
      break;
   case OPB_SENSE_NOTHING:
   default:
      SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
      return SCIP_INVALIDDATA;
   }

   /* create and add the linear constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/dynamicrows", &dynamicrows) );
   initial = !dynamicrows;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows;
   
#ifdef GENCONSNAMES
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif

   /* create corresponding constraint */
   if( issoftcons )
   {
      SCIP_VAR* indvar;
      SCIP_VAR* negindvar;
      SCIP_Bool created;
      int v;
      char indname[SCIP_MAXSTRLEN];
#if USEINDICATOR == FALSE
      SCIP_Real maxact;
      SCIP_Real minact;
      SCIP_Real lb;
      SCIP_Real ub;
#endif
      (void) SCIPsnprintf(indname, SCIP_MAXSTRLEN, "indicatorvar%d", opbinput->nindvars);
      ++(opbinput->nindvars);
      SCIP_CALL( createVariable(scip, &indvar, indname) );

      assert(!SCIPisInfinity(scip, -weight));
      SCIP_CALL( SCIPchgVarObj(scip, indvar, weight) );
      SCIP_CALL( SCIPgetNegatedVar(scip, indvar, &negindvar) );

#if USEINDICATOR == FALSE
      /* @todo check whether it's better to set the initial flag to false */         
      initial = FALSE;
      created = FALSE;

      maxact = 0.0;
      minact = 0.0;
      for( v = ncoefs - 1; v >= 0; --v )
         if( coefs[v] > 0 )
            maxact += coefs[v];
         else
            minact += coefs[v];

      if( SCIPisInfinity(scip, maxact) )
      {
         SCIPwarningMessage("maxactivity = %g exceed infinity value.\n", maxact);
      }
      if( SCIPisInfinity(scip, -minact) )
      {
         SCIPwarningMessage("minactivity = %g exceed -infinity value.\n", minact);
      }

      /* first soft constraints for lhs */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* first we are modelling the feasibility of the soft contraint by adding a slack variable */
         /* we ensure that if indvar == 1 => (a^T*x + ub*indvar >= lhs) */
         ub = lhs - minact;

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, lhs, SCIPinfinity(scip),
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, indvar, ub) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         created = TRUE;

#ifdef GENCONSNAMES
         ++(opbinput->consnumber);
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif
         if( isNonlinear )
            (*nNonlinearConss)++;

         /* second we are modelling the implication that if the slack variable is on( negation is off), the constraint
          * is disabled, so only the cost arise if the slack variable is necessary */
         /* indvar == 1 => (a^T*x (+ ub * negindvar) <= lhs - 1) */
         ub = lhs - maxact - 1;
         
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, -SCIPinfinity(scip), lhs - 1,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, negindvar, ub) );
      }

      /* second soft constraints for rhs */
      if( !SCIPisInfinity(scip, rhs) )
      {
         if( created )
         {
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
#ifdef GENCONSNAMES
            ++(opbinput->consnumber);
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif
            if( isNonlinear )
               (*nNonlinearConss)++;
         }

         /* first we are modelling the feasibility of the soft-constraint by adding a slack variable */
         /* indvar == 1 => (a^T*x + lb * indvar <= rhs) */
         lb = rhs - maxact;

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, -SCIPinfinity(scip), rhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, indvar, lb) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

#ifdef GENCONSNAMES
         ++(opbinput->consnumber);
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif
         if( isNonlinear )
            (*nNonlinearConss)++;

         /* second we are modelling the implication that if the slack variable is on( negation is off), the constraint
          * is disabled, so only the cost arise if the slack variable is necessary */
         /* indvar == 1 => (a^T*x (+ lb * negindvar) >= rhs + 1) */
         lb = rhs - minact + 1;
         
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, rhs + 1, SCIPinfinity(scip),
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, negindvar, lb) );
      }
#else /* with indicator */
      /* @todo check whether it's better to set the initial flag to false */         
      created = FALSE;

      if( !SCIPisInfinity(scip, rhs) )
      {
         /* first we are modelling the implication that if the negation of the indicator variable is on, the constraint
          * is enabled */
         /* indvar == 0 => a^T*x <= rhs */
         SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, negindvar, ncoefs, vars, coefs, rhs,
               initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         created = TRUE;

#if BACKIMPLICATION == TRUE
         {
            SCIP_Real* tmpcoefs;

#ifdef GENCONSNAMES
            ++(opbinput->consnumber);
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif
            if( isNonlinear )
               (*nNonlinearConss)++;

            /* allocate temporary memory */
            SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpcoefs, coefs, ncoefs) );

            /* second we are modelling the implication that if the indicator variable is on, the constraint is disabled */
            /* indvar == 1 =>  a^T*x >= rhs + 1 */
            /* change the a^T*x >= rhs + 1 to -a^Tx<= -rhs -1, for indicator constraint */
            for( v = ncoefs - 1; v >= 0; --v )
               tmpcoefs[v] *= -1;

            SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, indvar, ncoefs, vars, tmpcoefs, -rhs - 1,
                  initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
            /* free temporary memory */
            SCIPfreeBufferArray(scip, & tmpcoefs);
         }
#endif /* BACKIMPLICATION == TRUE */
      }

      if( !SCIPisInfinity(scip, -lhs) )
      {
         if( created )
         {
            /* we have already created one indicator constraint with the same indicator variable */
            SCIP_CALL( SCIPaddCons(scip, cons) ); /*lint !e644*/
            SCIPdebugMessage("(line %d) created constraint: ", opbinput->linenumber);
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
#ifdef GENCONSNAMES
            ++(opbinput->consnumber);
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif
            if( isNonlinear )
               (*nNonlinearConss)++;
         }

#if BACKIMPLICATION == TRUE
         /* first we are modelling the implication that if the indicator variable is on, the constraint is disabled */
         /* indvar == 1 => a^T*x <= lhs - 1 */
         SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, indvar, ncoefs, vars, coefs, lhs - 1,
               initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

#ifdef GENCONSNAMES
         ++(opbinput->consnumber);
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", opbinput->consnumber);
#endif
         if( isNonlinear )
            (*nNonlinearConss)++;
#endif /* BACKIMPLICATION == TRUE */

         /* second we are modelling the implication that if the negation of the indicator variable is on, the constraint
          * is enabled */
         /* change the a^T*x >= lhs to -a^Tx<= -lhs, for indicator constraint */
         for( v = ncoefs - 1; v >= 0; --v )
            coefs[v] *= -1;
         SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name, negindvar, ncoefs, vars, coefs, -lhs,
               initial, separate, enforce, check, propagate, local, dynamic, removable, FALSE) );
      }
#endif
   }
   else
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
   }
#ifdef GENCONSNAMES
   ++(opbinput->consnumber);
#endif
   
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugMessage("(line %d) created constraint: ", opbinput->linenumber);
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   if( isNonlinear )
      (*nNonlinearConss)++;

 TERMINATE:
   /* free memory */
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** tries to read the first comment line which usually contains information about the max size of "and" products */
static
SCIP_RETCODE getMaxAndConsDim(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           filename            /**< name of the input file */
   )
{   
   SCIP_Bool stop;
   char* commentstart;
   char* nproducts;
   int i;

   assert(scip != NULL);
   assert(opbinput != NULL);

   stop = FALSE;
   commentstart = NULL;
   nproducts = NULL;

   do
   {   
      if( SCIPfgets(opbinput->linebuf, sizeof(opbinput->linebuf), opbinput->file) == NULL )
      {
         assert(SCIPfeof( opbinput->file ) );
         break;
      }
      
      /* read characters after comment symbol */
      for( i = 0; commentchars[i] != '\0'; ++i )
      {
         commentstart = strchr(opbinput->linebuf, commentchars[i]);
         
         /* found a commentline */
         if( commentstart != NULL )
         { 
	    /* search for "#product= xxx" in commentline, where xxx represents the number of and constraints */
            nproducts = strstr(opbinput->linebuf, "#product= ");
            if( nproducts != NULL )
            {
               char* pos;
               nproducts += strlen("#product= ");

               pos = strtok(nproducts, delimchars);
               
               if( pos != NULL )
               { 
                  opbinput->sconsanddata = atoi(pos);
                  SCIPdebugMessage("%d of and constraints supposed to be in file.\n", opbinput->sconsanddata);
               }
               if( opbinput->sconsanddata < 0 )
                  opbinput->sconsanddata = 10;
               
               pos = strtok (NULL, delimchars);
               
               if( pos != NULL && strcmp(pos, "sizeproduct=") == 0 )
               {
                  pos = strtok (NULL, delimchars);
                  if( pos != NULL )
                  {
                     opbinput->maxvarsperand = atoi(pos) / opbinput->sconsanddata + 1;
                     SCIPdebugMessage("%d max length of and constraints in file.\n", opbinput->maxvarsperand);
                  }
               }
                  
	       if( opbinput->maxvarsperand < 0 )
                  opbinput->maxvarsperand = 10;

               stop = TRUE;
            }
            break;
         }
      }
   }
   while(commentstart != NULL && !stop);

   opbinput->linebuf[0] = '\0';

#if 0 /* following lines should be correct, but gzseek seems to not reseting status beeing at the end of file */
   /* reset filereader pointer to the beginning */
   (void) SCIPfseek(opbinput->file, 0, SEEK_SET);
#else
   SCIPfclose(opbinput->file);
   opbinput->file = SCIPfopen(filename, "r");
#endif

   return SCIP_OKAY;
}

/** reads an OPB file */
static
SCIP_RETCODE readOPBFile(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           filename            /**< name of the input file */
   )
{
   int nNonlinearConss;
   int i;

   assert(scip != NULL);
   assert(opbinput != NULL);
   
   /* open file */
   opbinput->file = SCIPfopen(filename, "r");
   if( opbinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }
   
   /* tries to read the first comment line which usually contains information about the max size of "and" products */
   SCIP_CALL( getMaxAndConsDim(scip, opbinput, filename) );

   /* reading additional information about the number of and constraints in comments to avoid reallocating
      "opbinput.andconss" */
   BMSclearMemoryArray(opbinput->linebuf, OPB_MAX_LINELEN);
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &(opbinput->consanddata), opbinput->sconsanddata ) );

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   nNonlinearConss = 0;

   while( !SCIPfeof( opbinput->file ) && !hasError(opbinput) )
   {
      SCIP_CALL( readConstraints(scip, opbinput, &nNonlinearConss) );
   }

   /* if we read a wbo file we need to make sure thta the top cost won't be exceeded */
   if( opbinput->wbo )
   {
      SCIP_VAR** topcostvars;
      SCIP_Real* topcosts;
      SCIP_VAR** vars;
      int nvars;
      int ntopcostvars;
      SCIP_Longint topcostrhs;
      SCIP_CONS* topcostcons;

      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
      assert(nvars > 0 || vars != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &topcostvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &topcosts, nvars) );

      ntopcostvars = 0;
      for( i = nvars - 1; i >= 0; --i )
         if( !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
         {
            topcostvars[ntopcostvars] = vars[i];
            topcosts[ntopcostvars] = SCIPvarGetObj(vars[i]);
            ++ntopcostvars;
         }

      if( SCIPisIntegral(scip, opbinput->topcost) )
         topcostrhs = (SCIP_Longint) SCIPfloor(scip, opbinput->topcost - 1);
      else
         topcostrhs = (SCIP_Longint) SCIPfloor(scip, opbinput->topcost);

      SCIP_CALL( SCIPcreateConsLinear(scip, &topcostcons, "topcost", ntopcostvars, topcostvars, topcosts, -SCIPinfinity(scip),
            (SCIP_Real) topcostrhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, topcostcons) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, topcostcons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &topcostcons) );

      SCIPfreeBufferArray(scip, &topcosts);
      SCIPfreeBufferArray(scip, &topcostvars);
   }

   for( i = opbinput->nconsanddata - 1; i >= 0; --i )
   {
      SCIPfreeMemoryArray(scip, &((opbinput->consanddata)[i]->vars) );
      SCIPfreeBlockMemory(scip, &((opbinput->consanddata)[i]) );
   }

   /* free dynamically allocated memory */
   SCIPfreeMemoryArrayNull(scip, &(opbinput->consanddata) );

   /* close file */
   SCIPfclose(opbinput->file);

   return SCIP_OKAY;
}


/*
 * Local methods (for writing)
 */

/** transforms given and constraint variables to the corresponding active or negated variables */
static
SCIP_RETCODE getBinVarsRepresentatives(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< vars array to get active variables for */
   int const             nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Bool const       transformed         /**< transformed constraint? */
   )
{
   SCIP_Bool negated;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );

   if( transformed )
   {
      for( v = nvars - 1; v >= 0; --v )
      {
         /** gets a binary variable that is equal to the given binary variable, and that is either active, fixed, or
          *  multi-aggregated, or the negated variable of an active, fixed, or multi-aggregated variable
          */
         SCIP_CALL( SCIPgetBinvarRepresentative( scip, vars[v], &vars[v], &negated) );
      }
   }
   else
   {
      SCIP_Real scalar;
      SCIP_Real constant;
      
      for( v = nvars - 1; v >= 0; --v )
      {
         scalar = 1.0;
         constant = 0.0;

         /** retransforms given variable, scalar and constant to the corresponding original variable, scalar and constant,
          *  if possible;
          *  if the retransformation is impossible, NULL is returned as variable
          */
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalar, &constant) );

         if( vars[v] == NULL )
         {
            SCIPdebugMessage("A variable coundn't retransformed to an original variable.\n");
            return SCIP_INVALIDDATA;
         }
         if( SCIPisEQ(scip, scalar, -1.0) && SCIPisEQ(scip, constant, 1.0) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &vars[v]) );
         }
         else
         {
            if( !SCIPisEQ(scip, scalar, 1.0) || !SCIPisZero(scip, constant) )
            {
               SCIPdebugMessage("A variable coundn't retransformed to an original variable or a negated variable of an original variable (scalar = %g, constant = %g).\n", scalar, constant);
               return SCIP_INVALIDDATA;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
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
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, requiredsize) );
         
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );

         if( vars[v] == NULL )
            return SCIP_INVALIDDATA;
      }

   return SCIP_OKAY;
}

/* computes all and-resultants and their corresonding constraint variables */
static
SCIP_RETCODE computeAndConstraintInfos(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_Bool const       transformed,        /**< transformed problem? */
   SCIP_VAR***           resvars,            /**< pointer to store all resultant variables */
   int*                  nresvars,           /**< pointer to store the number of all resultant variables */ 
   SCIP_VAR****          andvars,            /**< pointer to store to all resultant variables their corresponding active( or negated) and-constraint variables */
   int**                 nandvars,           /**< pointer to store the number of all corresponding and-variables to their corresponding resultant variable */
   SCIP_Bool*const       existandconshdlr,   /**< pointer to store whether the and-constrainthandler exists*/
   SCIP_Bool*const       existands           /**< pointer to store if their exists some and-constraints */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(resvars != NULL);
   assert(nresvars != NULL);
   assert(andvars != NULL);
   assert(nandvars != NULL);
   assert(existandconshdlr != NULL);
   assert(existands != NULL);

   *resvars = NULL;
   *nandvars = NULL;
   *andvars = NULL;
   *nresvars = 0;

   /* detect all and-resultants */
   conshdlr = SCIPfindConshdlr(scip, "and");
   if( conshdlr != NULL )
   {
      SCIP_CONS** andconss;
      int nandconss;
      int* shouldnotbeinand;
      int a;
      int c;
      int r;
      int v;
      int pos;
      int ncontainedands;

      andconss = NULL;
      nandconss = 0;
      *existandconshdlr = TRUE;

      /* if we write the original problem we need to get the original and constraints */
      if( !transformed )
      {
         SCIP_CONS** origconss;
         int norigconss;

         origconss = SCIPgetOrigConss(scip);
         norigconss = SCIPgetNOrigConss(scip);

         /* allocate memory for all possible and-constraints */
         SCIP_CALL( SCIPallocBufferArray(scip, &andconss, norigconss) );

         /* collect all original and-constraints */
         for( c = norigconss - 1; c >= 0; --c )
         {
            conshdlr = SCIPconsGetHdlr(origconss[c]);
            assert( conshdlr != NULL );

            if( strcmp(SCIPconshdlrGetName(conshdlr), "and") == 0 )
            {
               andconss[nandconss] = origconss[c];
               ++nandconss;
            }
         }
      }
      else
      {
         nandconss = SCIPconshdlrGetNConss(conshdlr);
         andconss = SCIPconshdlrGetConss(conshdlr);
      }
      
      assert(andconss != NULL || nandconss == 0);

      *nresvars = nandconss;

      if( nandconss > 0 )
      {
         *existands = TRUE;

         SCIP_CALL( SCIPallocMemoryArray(scip, resvars, *nresvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, andvars, *nresvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, nandvars, *nresvars) );

         /* collect all and-constraint variables */
         for( c = nandconss - 1; c >= 0; --c )
         {
	    (*nandvars)[c] = SCIPgetNVarsAnd(scip, andconss[c]);
	    SCIP_CALL( SCIPduplicateMemoryArray(scip, &((*andvars)[c]), SCIPgetVarsAnd(scip, andconss[c]), (*nandvars)[c]) );
	    SCIP_CALL( getBinVarsRepresentatives(scip, (*andvars)[c], (*nandvars)[c], transformed) );

	    (*resvars)[c] = SCIPgetResultantAnd(scip, andconss[c]);

            assert((*andvars)[c] != NULL && (*nandvars)[c] > 0);
            assert((*resvars)[c] != NULL);
         }

         /* sorted the array */
         SCIPsortPtrPtrInt((void**)(*resvars), (void**)(*andvars), (*nandvars), SCIPvarComp, (*nresvars));
      }
      else
         *existands = FALSE;

      SCIP_CALL( SCIPallocBufferArray(scip, &shouldnotbeinand, *nresvars) );
      
      /* check that all and-constraints doesn't contain any and-resultants, if they do try to resolve this */
      /* attention: if resolving leads to x = x*y*... , we can't do anything here ( this only means (... >=x and) y >= x, so normally the and-constraint needs to be 
         deleted and the inequality from before needs to be added ) */
      for( r = *nresvars - 1; r >= 0; --r )
      {
         ncontainedands = 0;
         shouldnotbeinand[ncontainedands] = r;
         ++ncontainedands;
         
         for( v = 0; v < (*nandvars)[r]; ++v )
         {
            if( SCIPsortedvecFindPtr((void**)(*resvars), SCIPvarComp, (*andvars)[r][v], *nresvars, &pos) )
            {
               
               for( a = ncontainedands - 1; a >= 0; --a )
                  if( shouldnotbeinand[a] == pos )
                  {
                     SCIPwarningMessage("This should not happen here. The and-constraint with resultant variable: ");
                     SCIP_CALL( SCIPprintVar(scip, (*resvars)[r], NULL) );
                     SCIPwarningMessage("contains a loop with and resultant:");
                     SCIP_CALL( SCIPprintVar(scip, (*resvars)[pos], NULL) );
                     return SCIP_INVALIDDATA;
                  }
               SCIPdebugMessage("Another and-contraint contains and-resultant:");
               SCIPdebug( SCIP_CALL( SCIPprintVar(scip, (*resvars)[pos], NULL) ) );
               SCIPdebugMessage("Trying to resolve.\n");
               
               shouldnotbeinand[ncontainedands] = pos;
               ++ncontainedands;
               
               /* try to resolve containing ands */
               (*nandvars)[r] = (*nandvars)[r] + (*nandvars)[pos] - 1;
               SCIP_CALL( SCIPreallocMemoryArray(scip, &((*andvars)[r]), (*nandvars)[r]) ); /*lint !e866 */
               
               for( a = (*nandvars)[pos] - 1; a >= 0; --a )
                  (*andvars)[r][(*nandvars)[r] - a - 1] = (*andvars)[pos][a];
               /* check same position with new variable */ 
               --v;
            }
         }
      }
      SCIPfreeBufferArray(scip, &shouldnotbeinand);
      
      /* free memory iff necessary */
      if( !transformed )
      {
         SCIPfreeBufferArray(scip, &andconss);
      }
   }
   else
   {
      SCIPdebugMessage("found no and-constraint-handler\n");
      *existands = FALSE;
      *existandconshdlr = FALSE;
   }

   return SCIP_OKAY;
}

/** clears the given line buffer */
static
void clearBuffer(
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of charaters in line */
   )
{
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   (*linecnt) = 0;
   linebuffer[0] = '\0';
}


/** ends the given line with '\\0' and prints it to the given file stream */
static
void writeBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of charaters in line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   if( (*linecnt) > 0 )
   {
      linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s", linebuffer);
      clearBuffer(linebuffer, linecnt);
   }
}


/** appends extension to line and prints it to the give file stream if the line buffer get full */
static
void appendBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of charaters in line */
   const char*           extension           /**< string to extent the line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );
   
   if( (*linecnt) + strlen(extension) >= OPB_MAX_LINELEN )
      writeBuffer(scip, file, linebuffer, linecnt);
   
   /* append extension to linebuffer */
   strncat(linebuffer, extension, OPB_MAX_LINELEN - (unsigned int)(*linecnt));
   (*linecnt) += (int) strlen(extension);
}

/* write objective function */
static
SCIP_RETCODE writeOpbObjective(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**const       vars,               /**< array with active (binary) variables */
   int const             nvars,              /**< number of mutable variables in the problem */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_OBJSENSE const   objsense,           /**< objective sense */
   SCIP_Real const       objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real const       objoffset,          /**< objective offset from bound shifting and fixing */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       existands,          /**< does some and-constraints exist? */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   SCIP_VAR* var;
   char linebuffer[OPB_MAX_LINELEN];
   char buffer[OPB_MAX_LINELEN];
   SCIP_Longint mult;
   SCIP_Bool objective;
   int v;
   int linecnt;
   int pos;

   assert(scip != NULL);
   assert(file != NULL);
   assert(vars != NULL || nvars == 0);
   assert(resvars != NULL || nresvars == 0);
   assert(andvars != NULL || nandvars == 0);
   assert(multisymbol != NULL);

   mult = 1;
   objective = FALSE;
   
   /* check if a objective function exits and compute the multiplier to
    * shift the coefficients to integers */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      
#ifndef NDEBUG
      {
	 /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
         if( !transformed )
	    assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL ||
               SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );
      }
#endif
      
      if( !SCIPisZero(scip, SCIPvarGetObj(var)) )
      {
         objective = TRUE;
         while( !SCIPisIntegral(scip, SCIPvarGetObj(var) * mult) ) 
	 {
	    assert(mult * 10 > mult);
	    mult *= 10;
	 }
      }
   }
   
   if( objective )
   {
      /* there exist a objective function*/
      SCIPinfoMessage(scip, file, "*   Obj. scale       : %.15g\n", objscale * mult);
      SCIPinfoMessage(scip, file, "*   Obj. offset      : %.15g\n", objoffset);
      
      clearBuffer(linebuffer, &linecnt);
      
      /* opb format supports only minimization; therefore, a maximization problem has to be converted */
      if( objsense == SCIP_OBJSENSE_MAXIMIZE )
         mult *= -1;
      
      SCIPdebugMessage("print objective function multiplyed with %"SCIP_LONGINT_FORMAT"\n", mult);
      
      appendBuffer(scip, file, linebuffer, &linecnt, "min:");

#ifndef NDEBUG
      if( existands )
      {
         int c;
         /* check that these variables are sorted */
         for( c = nresvars - 1; c > 0; --c )
            assert(SCIPvarGetIndex(resvars[c]) >= SCIPvarGetIndex(resvars[c - 1]));
      }
#endif

      for( v = nvars - 1; v >= 0; --v )
      {
         SCIP_Bool negated;
         var = vars[v];
         
         assert(var != NULL);

         if( SCIPisZero(scip, SCIPvarGetObj(var)) )
            continue;
         
         negated = SCIPvarIsNegated(var);
         
         assert( linecnt != 0 );

         /* replace and-resultant with corresponding variables */
         if( existands && SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, var, nresvars, &pos) )
         {
	    int a;

            assert(pos >= 0 && nandvars[pos] > 0 && andvars[pos] != NULL);

	    negated = SCIPvarIsNegated(andvars[pos][nandvars[pos] - 1]);

            /* print and-vars */
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+"SCIP_LONGINT_FORMAT"%s%s%s", 
	       (SCIP_Longint) (SCIPvarGetObj(var) * mult), multisymbol, negated ? "~" : "",
	       strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][nandvars[pos] - 1]) : andvars[pos][nandvars[pos] - 1]), "x"));
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         
            for(a = nandvars[pos] - 2; a >= 0; --a )
            {
	       negated = SCIPvarIsNegated(andvars[pos][a]);

               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][a]) : andvars[pos][a]), "x"));
               appendBuffer(scip, file, linebuffer, &linecnt, buffer);
            }
            
            appendBuffer(scip, file, linebuffer, &linecnt, " ");
	 }
         else
         {
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " %+"SCIP_LONGINT_FORMAT"%s%s%s", 
	       (SCIP_Longint) (SCIPvarGetObj(var) * mult), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         }
      }
      
      /* and objective function line ends with a ';' */
      appendBuffer(scip, file, linebuffer, &linecnt, " ;\n");
      writeBuffer(scip, file, linebuffer, &linecnt);
   }

   return SCIP_OKAY;
}

/* print maybe non linear row in OPB format to file stream */
static
SCIP_RETCODE printNLRow(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file (or NULL for standard output) */
   char const*const      type,               /**< row type ("=" or ">=") */
   SCIP_VAR**const       vars,               /**< array of variables */
   SCIP_Real const*const vals,               /**< array of values */
   int const             nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_Longint*const    mult,               /**< multiplier for the coefficients */  
   char const*const      multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR* var;
   char buffer[OPB_MAX_LINELEN];
   char linebuffer[OPB_MAX_LINELEN + 1];
   int v;
   int pos;
   int linecnt;
   
   assert(scip != NULL);
   assert(strcmp(type, "=") == 0 || strcmp(type, ">=") == 0);
   assert(mult != NULL);
   assert(resvars != NULL); 
   assert(nresvars > 0); 
   assert(andvars != NULL && nandvars != NULL); 
   
   clearBuffer(linebuffer, &linecnt);

   /* check if all coefficients are internal; if not commentstart multiplier */
   for( v = 0; v < nvars; ++v )
   {
      while( !SCIPisIntegral(scip, vals[v] * (*mult)) )
         (*mult) *= 10;
   }

   while( !SCIPisIntegral(scip, lhs * (*mult)) )
      (*mult) *= 10;
   
   /* print comment line if we have to multiply the coefficients to get integrals */
   if( ABS(*mult) != 1 )
      SCIPinfoMessage(scip, file, "* the following constraint is multiplied by %"SCIP_LONGINT_FORMAT" to get integral coefficients\n", ABS(*mult) );

#ifndef NDEBUG
   /* check that these variables are sorted */
   for( v = nresvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(resvars[v]) >= SCIPvarGetIndex(resvars[v - 1]));
#endif
   
   /* print coefficients */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Bool negated;

      var = vars[v];
      assert( var != NULL );

      negated = SCIPvarIsNegated(var);

      /* replace and-resultant with corresponding variables */
      if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, var, nresvars, &pos) )
      {
	 int a;
	 
	 assert(pos >= 0 && nandvars[pos] > 0 && andvars[pos] != NULL);
	 
	 negated = SCIPvarIsNegated(andvars[pos][nandvars[pos] - 1]);

	 /* print and-vars */
	 (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+"SCIP_LONGINT_FORMAT"%s%s%s", 
	    (SCIP_Longint) (vals[v] * (*mult)), multisymbol, negated ? "~" : "",
	    strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][nandvars[pos] - 1]) : andvars[pos][nandvars[pos] - 1]), "x") );
	 appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         
	 for(a = nandvars[pos] - 2; a >= 0; --a )
	 {
	    negated = SCIPvarIsNegated(andvars[pos][a]);
	    
	    (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][a]) : andvars[pos][a]), "x"));
	    appendBuffer(scip, file, linebuffer, &linecnt, buffer);
	 }
            
	 appendBuffer(scip, file, linebuffer, &linecnt, " ");
      }
      else
      {
         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+"SCIP_LONGINT_FORMAT"%s%s%s ", 
            (SCIP_Longint) (vals[v] * (*mult)), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);
      }
   }
   
   /* print left hand side */
   if( SCIPisZero(scip, lhs) )
      lhs = 0.0;
   
   (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s %"SCIP_LONGINT_FORMAT" ;\n", type, (SCIP_Longint) (lhs * (*mult)) );
   appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   
   writeBuffer(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}


/** prints given maybe non-linear constraint information in OPB format to file stream */
static
SCIP_RETCODE printNonLinearCons(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**const       vars,               /**< array of variables */
   SCIP_Real*const       vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int const             nvars,              /**< number of variables */
   SCIP_Real const       lhs,                /**< left hand side */
   SCIP_Real const       rhs,                /**< right hand side */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int  const*const      nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_Bool const       transformed,        /**< transformed constraint? */
   char const*const      multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   SCIP_Real activeconstant;
   SCIP_Longint mult;
   int v;
   int nactivevars;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(lhs <= rhs);
   assert(resvars != NULL); 
   assert(nresvars > 0); 
   assert(andvars != NULL && nandvars != NULL); 

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   activeconstant = 0.0;
   nactivevars = nvars;
   
   /* duplicate variable and value array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );
      
      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }
   
   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );
   
   mult = 1;

   /* print row(s) in LP format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equal constrain */
      SCIP_CALL( printNLRow(scip, file, "=", activevars, activevals, nactivevars, rhs - activeconstant, resvars, nresvars,
            andvars, nandvars, &mult, multisymbol) );
   }
   else
   { 
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printNLRow(scip, file, ">=", activevars, activevals, nactivevars, lhs - activeconstant, resvars, nresvars,
               andvars, nandvars, &mult, multisymbol) );
      }

      
      if( !SCIPisInfinity(scip, rhs) )
      {
         mult *= -1;

         /* print inequality ">=" and multiplying all coefficients by -1 */
         SCIP_CALL( printNLRow(scip, file, ">=", activevars, activevals, nactivevars, rhs - activeconstant, resvars, nresvars,
               andvars, nandvars, &mult, multisymbol) );
      }
   }
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/* print row in OPB format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           type,               /**< row type ("=" or ">=") */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Longint*         mult,               /**< multiplier for the coefficients */  
   const char*           multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR* var;
   char buffer[OPB_MAX_LINELEN];
   char linebuffer[OPB_MAX_LINELEN + 1];
   int v;
   int linecnt;
   
   assert(scip != NULL);
   assert(strcmp(type, "=") == 0 || strcmp(type, ">=") == 0);
   assert(mult != NULL);
   
   clearBuffer(linebuffer, &linecnt);

   /* check if all coefficients are internal; if not commentstart multiplier */
   for( v = 0; v < nvars; ++v )
   {
      while( !SCIPisIntegral(scip, vals[v] * (*mult)) )
         (*mult) *= 10;
   }

   while ( !SCIPisIntegral(scip, lhs * (*mult)) )
      (*mult) *= 10;
   
   /* print comment line if we have to multiply the coefficients to get integrals */
   if( ABS(*mult) != 1 )
      SCIPinfoMessage(scip, file, "* the following constraint is multiplied by %"SCIP_LONGINT_FORMAT" to get integral coefficients\n", ABS(*mult) );
   
   /* print coefficients */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Bool negated;

      var = vars[v];
      assert( var != NULL );

      negated = SCIPvarIsNegated(var);

      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+"SCIP_LONGINT_FORMAT"%s%s%s ", 
         (SCIP_Longint) (vals[v] * (*mult)), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   }
   
   /* print left hand side */
   if( SCIPisZero(scip, lhs) )
      lhs = 0.0;
   
   (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s %"SCIP_LONGINT_FORMAT" ;\n", type, (SCIP_Longint) (lhs * (*mult)) );
   appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   
   writeBuffer(scip, file, linebuffer, &linecnt);
}


/** prints given linear constraint information in OPB format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   const char*           multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant;
   SCIP_Longint mult;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   activeconstant = 0.0;
   
   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );
      
      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }
   
   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );
   
   mult = 1;

   /* print row(s) in LP format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equal constrain */
      printRow(scip, file, "=", activevars, activevals, nactivevars, rhs - activeconstant, &mult, multisymbol);
   }
   else
   { 
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRow(scip, file, ">=", activevars, activevals, nactivevars, lhs - activeconstant, &mult, multisymbol);
      }

      if( !SCIPisInfinity(scip, rhs) )
      {
         mult *= -1;

         /* print inequality ">=" and multiplying all coefficients by -1 */
         printRow(scip, file, ">=", activevars, activevals, nactivevars, rhs - activeconstant, &mult, multisymbol);
      }
   }
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

static
SCIP_RETCODE writeOpbConstraints(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_CONS**const      conss,              /**< array with constraints of the problem */
   int const             nconss,             /**< number of constraints in the problem */
   SCIP_VAR**const       vars,               /**< array with active (binary) variables */
   int const             nvars,              /**< number of mutable variables in the problem */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       existandconshdlr,   /**< does and-constrainthandler exist? */
   SCIP_Bool const       existands,          /**< does some and-constraints exist? */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;
   int v, c;

   assert(scip != NULL);
   assert(file != NULL);
   assert(conss != NULL || nconss == 0);
   assert(vars != NULL || nvars == 0);
   assert(resvars != NULL || nresvars == 0);
   assert(andvars != NULL || nandvars == 0);
   assert(multisymbol != NULL);

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );
      
      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );
      
      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         if( existands )
         {
            SCIP_CALL( printNonLinearCons(scip, file,
                  SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), resvars, nresvars, andvars, nandvars, 
                  transformed, multisymbol) );
         }            
         else
         {
            SCIP_CALL( printLinearCons(scip, file,
                  SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed, multisymbol) );
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            if( existands )
            {
               SCIP_CALL( printNonLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, 1.0, resvars, nresvars, andvars, nandvars, transformed, multisymbol) );
            }            
            else
            {
               SCIP_CALL( printLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, 1.0, transformed, multisymbol) );
            }
            break;
         case SCIP_SETPPCTYPE_PACKING :
            if( existands )
            {
               SCIP_CALL( printNonLinearCons(scip, file,
                     consvars, NULL, nconsvars, -SCIPinfinity(scip), 1.0, resvars, nresvars, andvars, nandvars,
                     transformed, multisymbol) );
            }            
            else
            {
               SCIP_CALL( printLinearCons(scip, file,
                     consvars, NULL, nconsvars, -SCIPinfinity(scip), 1.0, transformed, multisymbol) );
            }
            break;
         case SCIP_SETPPCTYPE_COVERING :
            if( existands )
            {
               SCIP_CALL( printNonLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, SCIPinfinity(scip), resvars, nresvars, andvars, nandvars,
                     transformed, multisymbol) );
            }            
            else
            {
               SCIP_CALL( printLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, SCIPinfinity(scip), transformed, multisymbol) );
            }
            break;
         }
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         if( existands )
         {
            SCIP_CALL( printNonLinearCons(scip, file,
                  SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip), 
                  resvars, nresvars, andvars, nandvars, transformed, multisymbol) );
         }     
         else
         {       
            SCIP_CALL( printLinearCons(scip, file,
                  SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
                  1.0, SCIPinfinity(scip), transformed, multisymbol) );
         }
      }
      else if ( strcmp(conshdlrname, "knapsack") == 0 )
      {
	 SCIP_Longint* weights;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = weights[v];

         if( existands )
         {
            SCIP_CALL( printNonLinearCons(scip, file, consvars, consvals, nconsvars, -SCIPinfinity(scip), 
                  (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), resvars, nresvars, andvars, nandvars,
                  transformed, multisymbol) );
         }     
         else
         {       
            SCIP_CALL( printLinearCons(scip, file, consvars, consvals, nconsvars, -SCIPinfinity(scip), 
                  (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed, multisymbol) );
         }

         SCIPfreeBufferArray(scip, &consvals);
      }
      else if ( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         if( existands )
         {
            SCIP_CALL( printNonLinearCons(scip, file, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons), 
                  SCIPgetRhsVarbound(scip, cons), resvars, nresvars, andvars, nandvars, transformed, multisymbol) );
         }     
         else
         {       
            SCIP_CALL( printLinearCons(scip, file, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons), 
                  SCIPgetRhsVarbound(scip, cons), transformed, multisymbol) );
         }

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if ( strcmp(conshdlrname, "and") == 0 )
      {
         /* all resultants of the and constraint will be replaced by all corresponding variables of this constraint, 
          * so no and-constraint will be printed directly */
	 assert(existandconshdlr);
      }
      else
      {
         SCIPwarningMessage("constraint handler <%s> can not print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "* ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
      }
   }

   return SCIP_OKAY;
}

/* write and constraints of inactive but relevant and-resulants and and variables which are fixed to one */
static
SCIP_RETCODE writeOpbRelevantAnds(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**const       resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   SCIP_VAR* resvar;
   SCIP_Longint rhslhs;
   char linebuffer[OPB_MAX_LINELEN];
   char buffer[OPB_MAX_LINELEN];
   int linecnt;
   int r, v;

   assert(scip != NULL);
   assert(file != NULL);
   assert(resvars != NULL || nresvars == 0);
   assert(andvars != NULL || nandvars == 0);
   assert(multisymbol != NULL);

   clearBuffer(linebuffer, &linecnt);

   /* print and-variables which are fixed to 1 where the and-resultant is not fixed, maybe doesn't appear and should only be asserted*/
   for( r = nresvars - 1; r >= 0; --r )
   {
      SCIP_VAR* var;
      SCIP_Bool neg;
      
      resvar = resvars[r];
      if( SCIPvarGetLbLocal(resvar) > 0.5 || SCIPvarGetUbLocal(resvar) < 0.5 )
         continue;
      for( v = nandvars[r] - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, andvars[r][v], &var, &neg) );
         
         if( SCIPvarGetLbLocal(andvars[r][v]) > 0.5 )
         {
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s = 1;\n", neg ? "~" : "", strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(var) : var), "x"));
	    appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         }
      }   
   }

   /* print and-constraints with fixed andresultant to zero and all and-constraints with 
      aggregated resultant, otherwise we would loose this information */
   for( r = nresvars - 1; r >= 0; --r )
   {
      resvar = resvars[r];
      rhslhs = (SCIPvarGetUbLocal(resvar) < 0.5) ? 0 : ((SCIPvarGetLbLocal(resvar) > 0.5) ? 1 : -1);

      if( rhslhs == 0 )
      {
         SCIP_Bool cont;
      
         cont = FALSE;
         /* if resultant variable and one other and variable is already zero, so we did'nt need to print this and constraint because all other variables are free*/
         for( v = nandvars[r] - 1; v >= 0; --v )
         {
             if( SCIPvarGetUbLocal(andvars[r][v]) < 0.5 )
             {
                cont = TRUE;
                break;
             }
         }
         if( cont )
            continue;
      }

      /* rhslhs equals to 0 means the and constraint is relavant due to it's not clear on which values the and variables are
       * rhslhs equals to 1 means the and constraint is irrelavant cause all and variables have to be 1 too
       * rhslhs equals to -1 means the and constraint is relavant cause the variable is only aggregated */
      if( rhslhs != 1 && !SCIPvarIsActive(resvar) )
      {
	 SCIP_VAR* var;
         SCIP_Bool neg;
         SCIP_Bool firstprinted;

         firstprinted = FALSE;

	 for( v = nandvars[r] - 1; v >= 0; --v )
	 {
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, andvars[r][v], &var, &neg) );
            
	    /* all fixed variables doesn't need printing */ 
	    if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
	       continue;
            
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", (firstprinted) ? multisymbol : "", neg ? "~" : "", strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(var) : var), "x"));
	    appendBuffer(scip, file, linebuffer, &linecnt, buffer);
            
            firstprinted = TRUE;
	 }
         
         /* if all and variables were fixed the resultant has to be fixed to we cannot be here */
         assert(firstprinted);
 
	 /* if the resultant is aggregated we need to print his binary representation */
	 if( rhslhs == -1 )
	 {
	    SCIP_VAR* myvar;
	    int pos;

	    assert(transformed);

	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, resvar, &myvar, &neg) );

	    if(neg)
	       assert(SCIPvarIsActive(SCIPvarGetNegationVar(myvar)));
	    else
	       assert(SCIPvarIsActive(myvar));

	    resvar = myvar;

	    /* replace and-resultant with corresponding variables */
	    if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, neg ? SCIPvarGetNegationVar(resvar) : resvar, nresvars, &pos) )
	    {
	       SCIP_Bool negated;
	       int a;
	 
	       assert(pos >= 0 && nandvars[pos] > 0 && andvars[pos] != NULL);
	 
	       negated = SCIPvarIsNegated(andvars[pos][nandvars[pos] - 1]);

	       /* print and-vars */
	       (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, neg ? " +1%s%s%s" : " -1%s%s%s", multisymbol, negated ? "~" : "",
                  strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][nandvars[pos] - 1]) : andvars[pos][nandvars[pos] - 1]), "x"));
	       appendBuffer(scip, file, linebuffer, &linecnt, buffer);
	       
	       for(a = nandvars[pos] - 2; a >= 0; --a )
	       {
		  negated = SCIPvarIsNegated(andvars[pos][a]);
	    
		  (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][a]) : andvars[pos][a]), "x"));
		  appendBuffer(scip, file, linebuffer, &linecnt, buffer);
	       }
            
	       appendBuffer(scip, file, linebuffer, &linecnt, " ");
	       
	       if( neg )
		  rhslhs = 1;
	       else
		  rhslhs = 0;
	    }
	    else
	    {
	       (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " -1%s%s%s", multisymbol, neg ? "~" : "", 
		  strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(resvar) : resvar), "x"));
	       appendBuffer(scip, file, linebuffer, &linecnt, buffer);

	       rhslhs = 0;
	    }
	 }
	 
	 /* print rhslhs */
	 (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " = %"SCIP_LONGINT_FORMAT" ;\n", rhslhs);
	 appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         
	 writeBuffer(scip, file, linebuffer, &linecnt);
	 
      }
   }

   return SCIP_OKAY;
}

/* writes problem to file */
static
SCIP_RETCODE writeOpb(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
					          extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active (binary) variables */
   int                   nvars,              /**< number of mutable variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_Bool const       existandconshdlr,   /**< does and-constrainthandler exist? */
   SCIP_Bool const       existands,          /**< does some and-constraints exist? */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   char multisymbol[OPB_MAX_LINELEN];
   SCIP_Bool usesymbole;

   assert( scip != NULL );
   assert( vars != NULL || nvars == 0 ); 
   assert( nvars == SCIPgetNBinVars(scip) );
   assert( conss != NULL || nconss == 0 ); 
   assert( result != NULL );

   /* check if should use a multipliers symbol star '*' between coefficients and variables */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/"READER_NAME"/multisymbol", &usesymbole) );
   (void) SCIPsnprintf(multisymbol, OPB_MAX_LINELEN, "%s", usesymbole ? " * " : " ");
   
   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (all binary)\n", nvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n", nconss);

   /* write objective function */
   SCIP_CALL( writeOpbObjective(scip, file, vars, nvars, resvars, nresvars, andvars, nandvars, 
				objsense, objscale, objoffset, multisymbol, existands, transformed) );

   /* write constraints */
   SCIP_CALL( writeOpbConstraints(scip, file, conss, nconss, vars, nvars, resvars, nresvars, andvars, nandvars, 
				  multisymbol, existandconshdlr, existands, transformed) );

   if( existands )
   {
      /* write and constraints of inactive but relevant and-resulants and and variables which are fixed to one 
         with no fixed and resultant */
      SCIP_CALL( writeOpbRelevantAnds(scip, file, resvars, nresvars, andvars, nandvars, multisymbol, transformed) );
   }

   *result = SCIP_SUCCESS;
   return  SCIP_OKAY;
}


/*
 * extern methods
 */

/* reads problem from file */
SCIP_RETCODE SCIPreadOpb(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   OPBINPUT opbinput;
   int i;

   /* initialize OPB input data */
   opbinput.file = NULL;
   opbinput.linebuf[0] = '\0';
   SCIP_CALL( SCIPallocBufferArray(scip, &opbinput.token, OPB_MAX_LINELEN) );
   opbinput.token[0] = '\0';
   SCIP_CALL( SCIPallocBufferArray(scip, &opbinput.tokenbuf, OPB_MAX_LINELEN) );
   opbinput.tokenbuf[0] = '\0';
   for( i = 0; i < OPB_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &opbinput.pushedtokens[i], OPB_MAX_LINELEN) );
   }

   opbinput.npushedtokens = 0;
   opbinput.linenumber = 1;
   opbinput.bufpos = 0;
   opbinput.linepos = 0;
   opbinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   opbinput.comment = FALSE;
   opbinput.endline = FALSE;
   opbinput.eof = FALSE;
   opbinput.haserror = FALSE;
   opbinput.consanddata = NULL;
   opbinput.nconsanddata = 0;
   opbinput.sconsanddata = 10;
   opbinput.nproblemcoeffs = 0;
   opbinput.maxvarsperand = 10;
   opbinput.wbo = FALSE;
   opbinput.topcost = -SCIPinfinity(scip);
   opbinput.nindvars = 0;
#ifdef GENCONSNAMES
   opbinput.consnumber = 0;
#endif

   /* create a hash table for the constraint set */
   opbinput.hashtablesize = SCIPcalcHashtableSize(HASHSIZE_OPBANDCONS);
   SCIP_CALL( SCIPhashtableCreate(&(opbinput.hashtable), SCIPblkmem(scip), opbinput.hashtablesize,
         hashGetKeyOpbAndcons, hashKeyEqOpbAndcons, hashKeyValOpbAndcons, (void*) scip) );

   /* read the file */
   SCIP_CALL( readOPBFile(scip, &opbinput, filename) );

   /* free hash table */
   SCIPhashtableFree(&(opbinput.hashtable));
   
   /* free dynamically allocated memory */
   SCIPfreeBufferArrayNull(scip, &opbinput.token);
   SCIPfreeBufferArrayNull(scip, &opbinput.tokenbuf);
   for( i = 0; i < OPB_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeBufferArrayNull(scip, &opbinput.pushedtokens[i]);
   }

   if( opbinput.nproblemcoeffs > 0 )
   {
      SCIPwarningMessage("there might be <%d> coefficients or weight out of range!\n", opbinput.nproblemcoeffs); 
   }

   /* evaluate the result */
   if( opbinput.haserror )
      return SCIP_READERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, opbinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}

/* writes problem to file */
SCIP_RETCODE SCIPwriteOpb(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file, or NULL if standard output should be used */
   const char*        name,               /**< problem name */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE      objsense,           /**< objective sense */
   SCIP_Real          objscale,           /**< scalar applied to objective function; external objective value is
					     extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real          objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                nvars,              /**< number of mutable variables in the problem */
   int                nbinvars,           /**< number of binary variables */
   int                nintvars,           /**< number of general integer variables */
   int                nimplvars,          /**< number of implicit integer variables */
   int                ncontvars,          /**< number of continuous variables */
   int                nfixedvars,         /**< number of fixed and aggregated variables in the problem */
   SCIP_CONS**        conss,              /**< array with constraints of the problem */
   int                nconss,             /**< number of constraints in the problem */
   SCIP_Bool          genericnames,       /**< should generic variable and constraint names be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file writing call */
   )
{  /*lint --e{715}*/
   if( nvars != nbinvars )
   {
      SCIPwarningMessage("OPB format is only capable for binary problems.\n");
      *result = SCIP_DIDNOTRUN;
   }
   else
   {
      SCIP_VAR*** andvars;
      SCIP_VAR** resvars;
      int* nandvars;
      SCIP_Bool existands;
      SCIP_Bool existandconshdlr;
      int nresvars;

      /* computes all and-resultants and their corresonding constraint variables */
      SCIP_CALL( computeAndConstraintInfos(scip, transformed, &resvars, &nresvars, &andvars, &nandvars, &existandconshdlr, &existands) );

      if( genericnames )
      {
#ifndef NDEBUG
         /* check for correct names for opb-format */
         int v;
         int idx;
         int pos;

         for( v = nvars - 1; v >= 0; --v )
         {
            if( existands )
            {
               if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
                  continue;
            }
            
            assert(sscanf(SCIPvarGetName(vars[v]), "x%d", &idx) == 1);
         }
#endif
	 SCIP_CALL( writeOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
               nvars, conss, nconss, resvars, nresvars, andvars, nandvars, existandconshdlr, existands, result) );
      }
      else
      {
         SCIP_Bool printed;
         int v;
         int idx;
         int pos;

         printed = FALSE;

         /* check if there are already generic names for all (not fixed variables)*/
         for( v = nvars - 1; v >= 0; --v )
            if( !existands || !SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
            {
               if( sscanf(SCIPvarGetName(vars[v]), transformed ? "t_x%d" : "x%d", &idx) != 1 )
               {
                  SCIPwarningMessage("At least following variable name isn't allowed in opb format.\n");
                  SCIP_CALL( SCIPprintVar(scip, vars[v], NULL) );
                  SCIPwarningMessage("OPB format needs generic variable names!\n");
                  
                  if( transformed )
                  {
                     SCIPwarningMessage("write transformed problem with generic variable names.\n");
                     SCIP_CALL( SCIPprintTransProblem(scip, file, "opb", TRUE) );
                  }
                  else
                  {
                     SCIPwarningMessage("write original problem with generic variable names.\n");
                     SCIP_CALL( SCIPprintOrigProblem(scip, file, "opb", TRUE) );
                  }
                  printed = TRUE;
                  break;
               }
            }

         /* check if there are already generic names for all (fixed variables)*/
         for( v = nfixedvars - 1; v >= 0; --v )
            if( !existands || !SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
            {
               if( sscanf(SCIPvarGetName(vars[v]), transformed ? "t_x%d" : "x%d", &idx) != 1 )
               {
                  SCIPwarningMessage("At least following variable name isn't allowed in opb format.\n");
                  SCIP_CALL( SCIPprintVar(scip, vars[v], NULL) );
                  SCIPwarningMessage("OPB format needs generic variable names!\n");
                  
                  if( transformed )
                  {
                     SCIPwarningMessage("write transformed problem with generic variable names.\n");
                     SCIP_CALL( SCIPprintTransProblem(scip, file, "opb", TRUE) );
                  }
                  else
                  {
                     SCIPwarningMessage("write original problem with generic variable names.\n");
                     SCIP_CALL( SCIPprintOrigProblem(scip, file, "opb", TRUE) );
                  }
                  printed = TRUE;
                  break;
               }
            }

         if( !printed )
         {
#ifndef NDEBUG
            for( v = nvars - 1; v >= 0; --v )
            {
               if( existands )
               {
                  if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
                     continue;
               }
               
               assert(sscanf(SCIPvarGetName(vars[v]), transformed ? "t_x%d" : "x%d", &idx) == 1);
            }
#endif
            SCIP_CALL( writeOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
                  nvars, conss, nconss, resvars, nresvars, andvars, nandvars, existandconshdlr, existands, result) );
         }

         if( existands )
         {
            /* free temporary buffers */
            assert(resvars != NULL);
            assert(andvars != NULL);
            assert(nandvars != NULL);
            
            for( v = nresvars - 1; v >= 0; --v )
            {
               assert(andvars[v] != NULL);
               SCIPfreeMemoryArray(scip, &andvars[v]);
            }
            SCIPfreeMemoryArray(scip, &resvars);
            SCIPfreeMemoryArray(scip, &andvars);
            SCIPfreeMemoryArray(scip, &nandvars);
         }
      }

      *result = SCIP_SUCCESS;
   }
   
   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyOpb)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderOpb(scip) );
 
   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeOpb NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadOpb)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadOpb(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteOpb)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPwriteOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, nfixedvars, conss, nconss, genericnames, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the opb file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderOpb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create opb reader data */
   readerdata = NULL;

   /* include opb reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyOpb,
         readerFreeOpb, readerReadOpb, readerWriteOpb,
         readerdata) );

   /* add opb reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/dynamicrows", "should rows be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/nlcseparate", "should the nonlinear constraint be separated during LP processing?",
         NULL, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/nlcpropagate", "should the nonlinear constraint be propagated during node processing?",
         NULL, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/nlcremovable", "should the nonlinear constraints be removable?",
         NULL, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/"READER_NAME"/multisymbol", "use '*' between coefficients and variables by writing to problem?",
         NULL, TRUE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
