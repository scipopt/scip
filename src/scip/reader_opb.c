/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_opb.c,v 1.25 2009/02/27 14:12:47 bzfheinz Exp $"

/**@file   reader_opb.c
 * @ingroup FILEREADERS 
 * @brief  pseudo-Boolean file reader (opb format)
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/* http://www.cril.univ-artois.fr/PB07/solver_req.html
 *
 * The syntax of the input file format can be described by a simple Backus-Naur
 *  form. <formula> is the start symbol of this grammar.
 *
 *  <formula>::= <sequence_of_comments> 
 *               [<objective>]
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
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/pub_misc.h"

#define READER_NAME             "opbreader"
#define READER_DESC             "file reader for pseudo-Boolean problem in opb format"
#define READER_EXTENSION        "opb"


/*
 * Data structures
 */
#define OPB_MAX_LINELEN       65536     /**< size of the line buffer for reading or writing */
#define OPB_MAX_PUSHEDTOKENS  2
#define OPB_INIT_COEFSSIZE    8192

/** Section in OPB File */
enum OpbExpType {
   OPB_EXP_NONE, OPB_EXP_UNSIGNED, OPB_EXP_SIGNED
};
typedef enum OpbExpType OPBEXPTYPE;

enum OpbSense {
   OPB_SENSE_NOTHING, OPB_SENSE_LE, OPB_SENSE_GE, OPB_SENSE_EQ
};
typedef enum OpbSense OPBSENSE;

/** OPB reading data */
struct OpbInput
{
   SCIP_FILE*           file;
   char                 linebuf[OPB_MAX_LINELEN];
   char*                token;
   char*                tokenbuf;
   char*                pushedtokens[OPB_MAX_PUSHEDTOKENS];
   int                  npushedtokens;
   int                  linenumber;
   int                  linepos;
   int                  bufpos;
   SCIP_OBJSENSE        objsense;
   SCIP_Bool            endline;
   SCIP_Bool            eof;
   SCIP_Bool            haserror;
   SCIP_CONS**          andconss;
   int                  nandconss;
   int                  sandconss;
   int                  nproblemcoeffs;
};

typedef struct OpbInput OPBINPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=;";
static const char commentchars[] = "*";

/*
 * Local methods (for reading)
 */

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
      /* buffer is full; erase last token since it might be incomplete */
      opbinput->endline = FALSE;
      last = strrchr(opbinput->linebuf, ' ');

      if( last == NULL )
      {
         SCIPwarningMessage("we read %d character from the file; these might indicates an corrupted input file!", OPB_MAX_LINELEN - 2);
         opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';
         SCIPdebugMessage("the buffer might be currented\n");
      }
      else
      {
         SCIPfseek(opbinput->file, -(long) strlen(last), SEEK_CUR);
         *last = '\0';
         SCIPdebugMessage("correct buffer\n");
      }
   }
   else 
   {
      /* found end of line */
      opbinput->endline = TRUE;
   }
   
   opbinput->linebuf[OPB_MAX_LINELEN-1] = '\0';
   opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(opbinput->linebuf, commentchars[i]);
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

   if( opbinput->token[1] == '\0' )
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
         if (strlen(opbinput->token)>18)
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

/** check if an and constraint exist with these variables */
static
SCIP_VAR* exitsAndCons(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of variables */
   )
{
   int i;
   int nandconss;
   int j,k;
   SCIP_CONS* cons;
   SCIP_VAR** andvars;
   SCIP_Bool found;
   
   assert( opbinput != NULL );
   assert( vars != NULL );
   assert( nvars > 1 );
   
   nandconss = opbinput->nandconss;

   for( i = 0; i < nandconss; ++i )
   {
      cons = opbinput->andconss[i];

      if( SCIPgetNVarsAnd(scip, cons) == nvars )
      {
         andvars = SCIPgetVarsAnd(scip, cons);
         
         found = TRUE;
         for( k = 0; k < nvars && found; ++k )
         {
            found = FALSE;
            for( j = 0; j < nvars && !found; ++j )
               if( andvars[k] == vars[j] )
                  found=TRUE;
         }
         if (found)
            return SCIPgetResultantAnd(scip, cons);
      }
   }
   
   return NULL;
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
   
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/dynamiccols", &dynamiccols) );
   initial = !dynamiccols;
   removable = dynamiccols;
   
   /* create new variable of the given name */
   SCIPdebugMessage("creating new variable: <%s>\n", name);
   SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
         initial, removable, NULL, NULL, NULL, NULL) );
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
   OPBINPUT*             opbinput,            /**< OPB reading data */
   SCIP_VAR**            var                 /**< pointer to store the variable */
   )
{
   SCIP_Bool negated;
   SCIP_Bool created = FALSE;
   char* name;

   SCIP_VAR** vars;
   int nvars;
   int svars;

   assert(var != NULL);
   
   nvars = 0;
   svars = 10;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, svars) );
   
   name = opbinput->token; 
   assert(name != NULL);
   

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
      if( nvars == svars )
      {
         svars *= 2;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, svars) );
      }
      
      vars[nvars] = *var;
      nvars++;
      
      if( !getNextToken(opbinput) )
         opbinput->haserror = TRUE;

      name = opbinput->token;
   }
   
   /* check if we found at least on variable */
   if( nvars == 0 )
      syntaxError(scip, opbinput, "expected a variable name");

   pushToken(opbinput);
   
   if( nvars > 1 )
   {
      if (!created)
         *var = exitsAndCons(scip, opbinput, vars, nvars);
      
      if( created || *var == NULL )
      {
         SCIP_Bool separate;
         SCIP_Bool propagate;
         SCIP_Bool removable;
         SCIP_Bool enforced;
         
         SCIP_CONS* cons;
         char varname[128];
         
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "andresultant%d", opbinput->nandconss);
         SCIP_CALL( createVariable(scip, var, varname) );
         assert( var != NULL );
        
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/nlcseparate", &separate) );
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/nlcpropagate", &propagate) );
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/nlcremovable", &removable) );
         SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/nlcenforced", &enforced) );
         
         SCIP_CALL( SCIPcreateConsAnd(scip, &cons, "", *var, nvars, vars,
               TRUE, separate, enforced, TRUE, propagate, FALSE, FALSE, FALSE, removable, FALSE) );
         
         if(opbinput->nandconss == opbinput->sandconss)
         {
            opbinput->sandconss *= 2;
            SCIP_CALL( SCIPreallocBufferArray(scip, &(opbinput->andconss), opbinput->sandconss) );
         }
         
         opbinput->andconss[opbinput->nandconss] = cons;
         opbinput->nandconss++;

         SCIPdebug( SCIPprintCons(scip, cons, NULL) );
         
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      }
   }
   
   SCIPfreeBufferArray(scip, &vars);

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
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
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
            strncpy(name, opbinput->tokenbuf, SCIP_MAXSTRLEN);
            name[SCIP_MAXSTRLEN-1] = '\0';
            SCIPdebugMessage("(line %d) read constraint name: '%s'\n", opbinput->linenumber, name);
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
         continue;
      }

      /* the token is a variable name: get the corresponding variable (or create a new one) */
      SCIP_CALL( getVariable(scip, opbinput, &var) );
      
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
         SCIP_CALL( SCIPchgVarObj(scip, vars[i], coefs[i]) );
      }
   }

   return SCIP_OKAY;
}

/** reads the constraints section */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput            /**< OPB reading data */
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
   int ncoefs;
   int sidesign;

   assert(opbinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, opbinput, name, &vars, &coefs, &ncoefs, &newsection) );

   if( hasError(opbinput) || opbinput->eof )
      goto TERMINATE;
   if( newsection )
   {
      if ( strcmp(name, "min") == 0 || strcmp(name, "max") == 0 )
      {
         /* set objective function  */
         SCIP_CALL( setObjective(scip, opbinput, name, vars, coefs, ncoefs) );
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
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/opbreader/dynamicrows", &dynamicrows) );
   initial = !dynamicrows;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows;
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncoefs, vars, coefs, lhs, rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugMessage("(line %d) created constraint: ", opbinput->linenumber);
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

 TERMINATE:
   /* free memory */
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &coefs);

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
   assert(opbinput != NULL);

   /* open file */
   opbinput->file = SCIPfopen(filename, "r");
   if( opbinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL) );

   while( !SCIPfeof( opbinput->file ) && !hasError(opbinput) )
   {
      SCIP_CALL( readConstraints(scip, opbinput) );
   }

   /* close file */
   SCIPfclose(opbinput->file);

   return SCIP_OKAY;
}


/* reads problem from file */
static
SCIP_RETCODE readFile(
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
   opbinput.endline = FALSE;
   opbinput.eof = FALSE;
   opbinput.haserror = FALSE;
   opbinput.andconss = NULL;
   opbinput.nandconss = 0;
   opbinput.sandconss = 10;
   opbinput.nproblemcoeffs = 0;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &(opbinput.andconss), opbinput.sandconss ) );
   
   /* read the file */
   SCIP_CALL( readOPBFile(scip, &opbinput, filename) );
   
   /* free dynamically allocated memory */
   SCIPfreeBufferArrayNull(scip, &(opbinput.andconss) );
   SCIPfreeBufferArrayNull(scip, &opbinput.token);
   SCIPfreeBufferArrayNull(scip, &opbinput.tokenbuf);
   for( i = 0; i < OPB_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeBufferArrayNull(scip, &opbinput.pushedtokens[i]);
   }

   if( opbinput.nproblemcoeffs > 0 )
   {
      SCIPwarningMessage("there might be <%d> coefficients out of range!\n", opbinput.nproblemcoeffs); 
   }

   /* evaluate the result */
   if( opbinput.haserror )
      return SCIP_PARSEERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, opbinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/*
 * Local methods (for writing)
 */


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
         *nvars = requiredsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, *nvars ) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, *nvars ) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
      for( v = 0; v < *nvars; ++v )
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
   
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
   strncat(linebuffer, extension, OPB_MAX_LINELEN - (*linecnt));
   (*linecnt) += (int) strlen(extension);
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
   SCIP_Longint*         mult                /**< multiplier for the coefficients */  
   )
{
   int v;
   char linebuffer[OPB_MAX_LINELEN + 1];
   int linecnt;

   SCIP_VAR* var;
   char buffer[OPB_MAX_LINELEN];
   
   assert( scip != NULL );
   assert( strcmp(type, "=") == 0 || strcmp(type, ">=") == 0 );
   assert( mult != NULL );
   
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
      var = vars[v];
      assert( var != NULL );
      
      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+"SCIP_LONGINT_FORMAT" %s ", 
         (SCIP_Longint) (vals[v] * (*mult)), SCIPvarGetName(var) );
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
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;
   SCIP_Longint mult;


   assert( scip != NULL );
   assert( vars != NULL || nvars == 0 );
   assert( nvars > 0 );
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;
   
   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
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
      printRow(scip, file, "=", activevars, activevals, nactivevars, rhs - activeconstant, &mult);
   }
   else
   { 
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRow(scip, file, ">=", activevars, activevals, nactivevars, lhs - activeconstant, &mult);
      }

      
      if( !SCIPisInfinity(scip, rhs) )
      {
         mult *= -1;

         /* print inequality ">=" and multiplying all coefficients by -1 */
         printRow(scip, file, ">=", activevars, activevals, nactivevars, rhs - activeconstant, &mult);
      }
   }
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/* writes problem to file */
static
SCIP_RETCODE writeOpb(
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
   SCIP_CONS**        conss,              /**< array with constraints of the problem */
   int                nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*       result              /**< pointer to store the result of the file writing call */
   )
{
   int c,v;
   SCIP_Longint mult;
   SCIP_Bool objective;
   
   int linecnt;
   char linebuffer[OPB_MAX_LINELEN];
   char buffer[OPB_MAX_LINELEN];
   
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;
   
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;

   assert( scip != NULL );
   
   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n", nconss);

   mult = 1;
   objective = FALSE;
   
   /* check if a objective function exits and compute the multiplier to
    * shift the coefficients to integers */
   for (v = 0; v < nvars; ++v)
   {
      var = vars[v];
      
#ifndef NDEBUG
      {
	 /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
	 int idx;
	 
	 if ( !transformed )
	    assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL ||
		    SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );

	 /* the variable name have to be of the form x%d */
	 assert( sscanf(SCIPvarGetName(var), "x%d", &idx) == 1 );
      }
#endif
      
      if ( !SCIPisZero(scip, SCIPvarGetObj(var)) )
      {
         objective = TRUE;
         while( !SCIPisIntegral(scip, SCIPvarGetObj(var) * mult) ) 
            mult *= 10;
         
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
      
      for (v = 0; v < nvars; ++v)
      {
         var = vars[v];
         
         
         if (SCIPisZero(scip, SCIPvarGetObj(var)) )
            continue;
         
         assert( linecnt != 0 );
         
         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " %+"SCIP_LONGINT_FORMAT" %s", 
            (SCIP_Longint) (SCIPvarGetObj(var) * mult), SCIPvarGetName(var) );
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);
      }
      
      /* ane dobjective function line with an ';' */
      appendBuffer(scip, file, linebuffer, &linecnt, " ;\n");
      writeBuffer(scip, file, linebuffer, &linecnt);
   }
   
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL);
      
      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      if( transformed && !SCIPconsIsEnabled(cons) )
         continue;
      
      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file,
               SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
               SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printLinearCons(scip, file,
                  consvars, NULL, nconsvars, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printLinearCons(scip, file,
                  consvars, NULL, nconsvars, -SCIPinfinity(scip), 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printLinearCons(scip, file,
                  consvars, NULL, nconsvars, 1.0, SCIPinfinity(scip), transformed) );
            break;
         }
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file,
               SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
               1.0, SCIPinfinity(scip), transformed) );
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

         SCIP_CALL( printLinearCons(scip, file, consvars, consvals, nconsvars, -SCIPinfinity(scip), 
				    (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed) );

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

         SCIP_CALL( printLinearCons(scip, file,
               consvars, consvals, 2,
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else
      {
         SCIPwarningMessage("constraint handler <%s> can not print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "* ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
      }
   }

   *result = SCIP_SUCCESS;
   return  SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeOpb NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadOpb)
{  /*lint --e{715}*/

   SCIP_CALL( readFile(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteOpb)
{  /*lint --e{715}*/
   if( nvars != nbinvars )
   {
      SCIPwarningMessage("OPB format is only capable for binary problems.\n");
      *result = SCIP_DIDNOTRUN;
   }
   else
   {
      if( genericnames )
      {
         SCIP_CALL( writeOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
               nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );
      }
      else
      {
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
      }
      *result = SCIP_SUCCESS;
   }
   
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
         readerFreeOpb, readerReadOpb, readerWriteOpb, readerdata) );

   /* add opb reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/dynamiccols", "should columns be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/dynamicrows", "should rows be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/nlcseparate", "should the nonlinear constraint be separated during LP processing?",
         NULL, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/nlcpropagate", "should the nonlinear constraint be propagated during node processing?",
         NULL, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/nlcremovable", "should the nonlinear constraints be removable?",
         NULL, TRUE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/opbreader/nlcenforced", "should the nonlinear constraints be enforced?",
         NULL, TRUE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
