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

/**@file   reader_rtp.c
 * @brief  Rtp file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "reader_rtp.h"
#include "cons_bitvar.h"
#include "cons_bitarith.h"


#define READER_NAME             "rtpreader"
#define READER_DESC             "rtp file reader"
#define READER_EXTENSION        "rtp"




/*
 * Data structures
 */

/* TODO: (optional) fill in the necessary reader data */

/** data for rtp reader */
struct ReaderData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */

#define MAX_LINE_LEN   65536
#define MAX_FNAME_LEN  256
#define RTP_MAGIC      "012CEA02"

enum reading_states { RTP_ERROR = 0, RTP_START, RTP_SYMB, RTP_CONS, RTP_PROP, RTP_END };

typedef enum reading_states STATE;

static
void parseError(
   SCIP*            scip,               /**< SCIP data structure */
   int              lineno,             /**< current line number of input file */
   const char*      msg,                /**< error message to display */
   const char*      erritem,            /**< token where the error occured, or NULL */
   STATE*           state               /**< pointer to current reading state */
   )
{
   char s[MAXSTRLEN];

   assert(msg != NULL);
   assert(state != NULL);

   if( erritem != NULL )
      sprintf(s, "Line %d: %s <%s>", lineno, msg, erritem);
   else
      sprintf(s, "Line %d: %s", lineno, msg);
   SCIPmessage(scip, SCIP_VERBLEVEL_MINIMAL, s);
   *state = RTP_ERROR;
}
   
static
RETCODE getStart(
   SCIP*            scip,               /**< SCIP data structure */
   int              lineno,             /**< current line number of input file */
   char*            s,                  /**< current line */
   STATE*           state               /**< pointer to current reading state */
   )
{
   assert(s != NULL);
   assert(state != NULL);

   if (lineno == 1 && !strncmp(s, RTP_MAGIC, 8))
      *state = RTP_START;
   else if (lineno == 1)
      parseError(scip, lineno, "wrong magic number found", NULL, state);
   else if (!strncmp(s, "SYMBOL", 6))
      *state = RTP_SYMB;
   else
      parseError(scip, lineno, "unknown keyword", s, state);

   return SCIP_OKAY;
}

static
RETCODE getSymbol(
   SCIP*            scip,               /**< SCIP data structure */
   int              lineno,             /**< current line number of input file */
   char*            s,                  /**< current line */
   STATE*           state               /**< pointer to current reading state */
   )
{
   CONS* cons;
   int tag;
   int width;
   int i;
   int p;
   char* cstart;
   char key;
   char name[32];
   
   assert(s != NULL);
   assert(state != NULL);

   if (!strncmp(s, "CONS", 4))
   {
      *state = RTP_CONS;
      return SCIP_OKAY;
   }

   if (2 != sscanf(s, "@%d [%d]", &tag, &width))
   {
      parseError(scip, lineno, "invalid symbol", NULL, state);
      return SCIP_OKAY;
   }

   s = strchr(s, ']');
   s++;

   sprintf(name, "x%d", tag);
    
   if (width <= 0)
   {
      parseError(scip, lineno, "invalid symbol width", NULL, state);
      return SCIP_OKAY;
   }    
   while(isspace(*s))
      s++;

   /* Binary Constant 
    */
   key = tolower(*s);
    
   if (key == 'b' || key == 'o' || key == 'd' || key == 'h')
   {
      cstart = s;
      while(*s != '\0' && !isspace(*s))
         s++;

      *s = '\0';
        
      CHECK_OKAY( SCIPcreateConsBitconstString(scip, &cons, name, width, 0.0, cstart) );
      debugMessage("constant <%s> [%d] = %s\n", SCIPconsGetName(cons), width, cstart);
   }
   else
   {
      CHECK_OKAY( SCIPcreateConsBitvar(scip, &cons, name, width, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE) );
      debugMessage("variable <%s> [%d]\n", SCIPconsGetName(cons), width);
   }
   CHECK_OKAY( SCIPaddCons(scip, cons) );
   CHECK_OKAY( SCIPreleaseCons(scip, &cons) );

   *state = RTP_SYMB;

   return SCIP_OKAY;
}

static
RETCODE parseFunction(
   SCIP*            scip,               /**< SCIP data structure */
   int              lineno,             /**< current line number of input file */
   int              ftag,               /**< function number tag */
   char*            s,                  /**< current line */
   STATE*           state               /**< pointer to current reading state */
   )
{
   static const struct
   {
      const char*  fname;
      BITARITHTYPE opcode;
   } optab[] = {
      { "Add", SCIP_BITARITHTYPE_ADD },
      { "Sub", SCIP_BITARITHTYPE_SUB },
      { "Eq", SCIP_BITARITHTYPE_EQ },
      { NULL, -1 }
   };
   char fname[MAX_FNAME_LEN + 51];
   char* fargsizes;
   int tag1;
   int tag2;
   int i;
   int cnt;
   char name[32];
   CONS* cons;
   CONS* cons1 = NULL;
   CONS* cons2 = NULL;
   CONS* rcons = NULL;

   assert(s != NULL);
   assert(state != NULL);
   assert(*state == RTP_CONS || *state == RTP_PROP);

   cnt = sscanf(s, "%63s @%d @%d", fname, &tag1, &tag2);
   cnt--; /* don't count the function's name */
   if( cnt < 0 )
   {
      parseError(scip, lineno, "syntax error", NULL, state);
      return SCIP_OKAY;
   }

   fargsizes = strchr(fname, '_');
   if( fargsizes != NULL )
   {
      *fargsizes = '\0';
      fargsizes++;
   }

   for( i = 0; optab[i].fname != NULL; i++ )
      if( !strcmp(fname, optab[i].fname) )
         break;
   
   if( optab[i].fname == NULL )
   {
      parseError(scip, lineno, "unknown function", fname, state);
      return SCIP_OKAY;
   }
   
   if( cnt != SCIPgetArityBitarith(optab[i].opcode) )
   {
      parseError(scip, lineno, "wrong number of parameters", NULL, state);
      return SCIP_OKAY;
   }
   if( cnt > 0 )
   {
      sprintf(name, "x%d", tag1);
      cons1 = SCIPfindCons(scip, name);
      if( cons1 == NULL )
      {
         parseError(scip, lineno, "unknown variable", name, state);
         return SCIP_OKAY;
      }

      if (cnt > 1)
      {
         sprintf(name, "x%d", tag2);
         cons2 = SCIPfindCons(scip, name);
         if( cons2 == NULL )
         {
            parseError(scip, lineno, "unknown variable", name, state);
            return SCIP_OKAY;
         }
      }
   }
   /* Is this a property ? */
   if (ftag < 0)
   {
      /*int width = MAX(SCIPgetNBitsBitvar(cons1), SCIPgetNBitsBitvar(cons2));*/ /*???????????????????*/
      CHECK_OKAY( SCIPcreateConsBitconstString(scip, &rcons, "property", 1, 0.0, "b0") );
      debugMessage("property <%s> [1]\n", SCIPconsGetName(rcons));
   }
   else
   {
      /* Return value */
      sprintf(name, "x%d", ftag);
      rcons = SCIPfindCons(scip, name);
      if( rcons == NULL )
      {
         parseError(scip, lineno, "unknown variable", name, state);
         return SCIP_OKAY;
      }
   }
   sprintf(name, "l%d", lineno);
   
#ifdef DEBUG
   switch( cnt )
   {
   case 0:
      debugMessage("constraint <%s>: <%s> == %s()\n", name, SCIPconsGetName(rcons), optab[i].fname);
      break;
   case 1:
      debugMessage("constraint <%s>: <%s> == %s(<%s>)\n", name, SCIPconsGetName(rcons), optab[i].fname,
         SCIPconsGetName(cons1));
      break;
   case 2:
      debugMessage("constraint <%s>: <%s> == %s(<%s>,<%s>)\n", name, SCIPconsGetName(rcons), optab[i].fname,
         SCIPconsGetName(cons1), SCIPconsGetName(cons2));
      break;
   }
#endif

   CHECK_OKAY( SCIPcreateConsBitarith(scip, &cons, name, optab[i].opcode, cons1, cons2, rcons,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
   CHECK_OKAY( SCIPaddCons(scip, cons) );
   CHECK_OKAY( SCIPreleaseCons(scip, &cons) );

   if (ftag < 0)
   {
      CHECK_OKAY( SCIPreleaseCons(scip, &rcons) );
   }

   return SCIP_OKAY;
}

static
RETCODE getConstraint(
   SCIP*            scip,               /**< SCIP data structure */
   int              lineno,             /**< current line number of input file */
   char*            s,                  /**< current line */
   STATE*           state               /**< pointer to current reading state */
   )
{
   int tag;
   
   assert(s != NULL);
   assert(state != NULL);
   assert(*state = RTP_CONS);

   if (!strncmp(s, "PROP", 4))
   {
      *state = RTP_PROP;
      return SCIP_OKAY;
   }

   if (1 != sscanf(s, "@%d", &tag))
   {
      parseError(scip, lineno, "invalid constraint tag", s, state);
      return SCIP_OKAY;
   }
   while((*s != '\0') && !isspace(*s))
      s++;

   while(isspace(*s))
      s++;

   CHECK_OKAY( parseFunction(scip, lineno, tag, s, state) );

   return SCIP_OKAY;
}

static
RETCODE getProperty(
   SCIP*            scip,               /**< SCIP data structure */
   int              lineno,             /**< current line number of input file */
   char*            s,                  /**< current line */
   STATE*           state               /**< pointer to current reading state */
   )
{
   assert(s != NULL);
   assert(state != NULL);
   assert(*state = RTP_PROP);

   if (!strncmp(s, "END", 4))
   {
      *state = RTP_END;
      return SCIP_OKAY;
   }

   CHECK_OKAY( parseFunction(scip, lineno, -1, s, state) );
   
   return SCIP_OKAY;
}

/* register transfer property */
static
RETCODE readFile(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< name of input file */
   )
{
   FILE* fp;
   char  buf[MAX_LINE_LEN];
   int   lineno = 0;
   char* s;
   STATE state = RTP_START;

   assert(filename != NULL);

   state = RTP_START;

   if( NULL == (fp = fopen(filename, "r")) )
   {
      perror(filename);
      return SCIP_READERROR;
   }

   CHECK_OKAY( SCIPcreateProb(scip, "RTP", NULL, NULL, NULL) );

   while( state != RTP_END && state != RTP_ERROR && (NULL != fgets(buf, sizeof(buf), fp)) )
   {
      lineno++;

      if( NULL != (s = strpbrk(buf, "#\r\n")) )
         *s = '\0';
      else
      {
         parseError(scip, lineno, "line truncated", NULL, &state);
         break;
      }
      s = buf;

      while(isspace(*s))
         s++;

      if (*s == '\0')
         continue;

      switch( state )
      {
      case RTP_ERROR:
         break;
      case RTP_START:
         CHECK_OKAY( getStart(scip, lineno, s, &state) );
         break;
      case RTP_SYMB:
         CHECK_OKAY( getSymbol(scip, lineno, s, &state) );
         break;
      case RTP_CONS:
         CHECK_OKAY( getConstraint(scip, lineno, s, &state) );
         break;
      case RTP_PROP:
         CHECK_OKAY( getProperty(scip, lineno, s, &state) );
         break;
      default:
         errorMessage("invalid RTP state");
         abort();
      }	
   }
   fclose(fp);

   if( state != RTP_END && state != RTP_ERROR )
      parseError(scip, lineno, "unexpected EOF", NULL, &state);

   if( state == RTP_ERROR )
      return SCIP_PARSEERROR;
   else
      return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeRtp NULL


/** problem reading method of reader */
static
DECL_READERREAD(readerReadRtp)
{
   CHECK_OKAY( readFile(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * reader specific interface methods
 */

/** includes the rtp file reader in SCIP */
RETCODE SCIPincludeReaderRtp(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   READERDATA* readerdata;

   /* create rtp reader data */
   readerdata = NULL;
   
   /* include rtp reader */
   CHECK_OKAY( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
                  readerFreeRtp, readerReadRtp, readerdata) );

   return SCIP_OKAY;
}
