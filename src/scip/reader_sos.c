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

/**@file   reader_sos.c
 * @ingroup FILEREADERS 
 * @brief  SOS file reader
 * @author Marc Pfetsch
 *
 * This basically a modification of the MPS file reader. The code can
 * probably be streamlined a lot by removing unnecessary code left
 * over from the MPS reader.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/reader_sos.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/pub_misc.h"


#define READER_NAME             "sosreader"
#define READER_DESC             "file reader for specially ordered sets (SOS)"
#define READER_EXTENSION        "sos"



/*
 * sos reader internal methods
 */

#define SOS_MAX_LINELEN 256

#define PATCH_CHAR    '_'
#define BLANK         ' '

enum SosSection
{
   SOS_NAME,
   SOS_SOSSECTION,
   SOS_ENDATA
};
typedef enum SosSection SOSSECTION;

struct SosInput
{
   SOSSECTION           section;
   SCIP_FILE*           fp;
   int                  lineno;
   SCIP_Bool            haserror;
   char                 buf[SOS_MAX_LINELEN];
   const char*          f0;
   const char*          f1;
   const char*          f2;
   const char*          f3;
   const char*          f4;
   const char*          f5;
};
typedef struct SosInput SOSINPUT;



static
SCIP_RETCODE SOSinputCreate(
   SCIP*                 scip,
   SOSINPUT**            sosi,
   SCIP_FILE*            fp
   )
{
   assert(sosi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocMemory(scip, sosi) );

   (*sosi)->section     = SOS_NAME;
   (*sosi)->fp          = fp;
   (*sosi)->lineno      = 0;
   (*sosi)->haserror    = FALSE;
   (*sosi)->buf     [0] = '\0';
   (*sosi)->f0          = NULL;
   (*sosi)->f1          = NULL;
   (*sosi)->f2          = NULL;
   (*sosi)->f3          = NULL;
   (*sosi)->f4          = NULL;
   (*sosi)->f5          = NULL;

   return SCIP_OKAY;
}

static
void SOSinputFree(
   SCIP*                 scip,
   SOSINPUT**            sosi
   )
{
   SCIPfreeMemory(scip, sosi);
}

static
SOSSECTION SOSinputSection(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->section;
}

static
const char* SOSinputField0(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->f0;
}

static
const char* SOSinputField1(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->f1;
}

static
const char* SOSinputField2(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->f2;
}

static
const char* SOSinputField3(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->f3;
}

static
const char* SOSinputField4(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->f4;
}

static
const char* SOSinputField5(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->f5;
}

static
SCIP_Bool SOSinputHasError(
   const SOSINPUT*       sosi
   )
{
   assert(sosi != NULL);

   return sosi->haserror;
}

static
void SOSinputSetSection(
   SOSINPUT*             sosi,
   SOSSECTION            section
   )
{
   assert(sosi != NULL);

   sosi->section = section;
}


static
void SOSinputSyntaxerror(
   SOSINPUT*             sosi
   )
{
   assert(sosi != NULL);

   SCIPwarningMessage("Syntax error in line %d\n", sosi->lineno);
   sosi->section  = SOS_ENDATA;
   sosi->haserror = TRUE;
}

/* fill the line from \p pos up to column 80 with blanks. */
static
void clearFrom(
   char*                 buf,
   unsigned int          pos
   )
{
   unsigned int i;

   for(i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/* change all blanks inside a field to #PATCH_CHAR. */
static
void patchField(
   char*                 buf,
   int                   beg,
   int                   end
   )
{
   int i;

   while((beg <= end) && (buf[end] == BLANK))
      end--;

   while((beg <= end) && (buf[beg] == BLANK))
      beg++;

   for(i = beg; i <= end; i++)
      if (buf[i] == BLANK)
         buf[i] = PATCH_CHAR;
}

/* read a sos format data line and parse the fields. */
static
SCIP_Bool SOSinputReadLine(
   SOSINPUT*             sosi
   )
{
   unsigned int len;
   unsigned int i;
   int space;
   char* s;
   SCIP_Bool is_marker;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      sosi->f0  = sosi->f1 = sosi->f2 = sosi->f3 = sosi->f4 = sosi->f5 = 0;
      is_marker = FALSE;
      is_empty = FALSE;

      /* Read until we have a non comment line. */
      do
      {
         if (NULL == SCIPfgets(sosi->buf, sizeof(sosi->buf), sosi->fp))
            return FALSE;
         sosi->lineno++;
      }
      while (*sosi->buf == '*');

      /* Normalize line */
      len = strlen(sosi->buf);

      for(i = 0; i < len; i++)
      {
         if ((sosi->buf[i] == '\t') || (sosi->buf[i] == '\n') || (sosi->buf[i] == '\r'))
            sosi->buf[i] = BLANK;
      }

      if (len < 80)
         clearFrom(sosi->buf, len);

      assert(strlen(sosi->buf) >= 80);

      /* Look for new section */
      if (*sosi->buf != BLANK)
      {
         sosi->f0 = SCIPstrtok(&sosi->buf[0], " ", &nexttok);

         assert(sosi->f0 != 0);

         sosi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      /* Test for fixed format comments */
      if ((sosi->buf[14] == '$') && (sosi->buf[13] == ' '))
	 clearFrom(sosi->buf, 14);
      else if ((sosi->buf[39] == '$') && (sosi->buf[38] == ' '))
	 clearFrom(sosi->buf, 39);

      /* Test for fixed format */
      space = sosi->buf[12] | sosi->buf[13]
	 | sosi->buf[22] | sosi->buf[23]
	 | sosi->buf[36] | sosi->buf[37] | sosi->buf[38]
	 | sosi->buf[47] | sosi->buf[48]
	 | sosi->buf[61] | sosi->buf[62] | sosi->buf[63];

      if (space == BLANK)
      {
	 /* Now we have space at the right positions.
	  * But are there also the non space where they
	  * should be ?
	  */
	 SCIP_Bool number = isdigit(sosi->buf[24]) || isdigit(sosi->buf[25])
	    || isdigit(sosi->buf[26]) || isdigit(sosi->buf[27])
	    || isdigit(sosi->buf[28]) || isdigit(sosi->buf[29])
	    || isdigit(sosi->buf[30]) || isdigit(sosi->buf[31])
	    || isdigit(sosi->buf[32]) || isdigit(sosi->buf[33])
	    || isdigit(sosi->buf[34]) || isdigit(sosi->buf[35]);

	 /* len < 13 is handle ROW lines with embedded spaces
	  * in the names correctly
	  */
	 if (number || len < 13)
	 {
	    /* We assume fixed format, so we patch possible embedded spaces. */
	    patchField(sosi->buf,  4, 12);
	    patchField(sosi->buf, 14, 22);
	    patchField(sosi->buf, 39, 47);
	 }
      }
      s = &sosi->buf[1];

      /* At this point it is not clear if we have a indicator field.
       * If there is none (e.g. empty) f1 will be the first name field.
       * If there is one, f2 will be the first name field.
       *
       * Initially comment marks '$' ar only allowed in the beginning
       * of the 2nd and 3rd name field. We test all fields but the first.
       * This makes no difference, since if the $ is at the start of a value
       * field, the line will be errornous anyway.
       */
      do
      {
	 if (NULL == (sosi->f1 = SCIPstrtok(s, " ", &nexttok)))
	    break;

	 if ((NULL == (sosi->f2 = SCIPstrtok(NULL, " ", &nexttok))) || (*sosi->f2 == '$'))
	 {
	    sosi->f2 = 0;
	    break;
	 }
	 if (!strcmp(sosi->f2, "'MARKER'"))
	    is_marker = TRUE;

	 if ((NULL == (sosi->f3 = SCIPstrtok(NULL, " ", &nexttok))) || (*sosi->f3 == '$'))
	 {
	    sosi->f3 = 0;
	    break;
	 }
	 if (is_marker)
	 {
	    break; /* unknown marker */
	 }
	 if (!strcmp(sosi->f3, "'MARKER'"))
	    is_marker = TRUE;

	 if ((NULL == (sosi->f4 = SCIPstrtok(NULL, " ", &nexttok))) || (*sosi->f4 == '$'))
	 {
	    sosi->f4 = 0;
	    break;
	 }
	 if (is_marker)
	 {
	    break; /* unknown marker */
	 }
	 if ((NULL == (sosi->f5 = SCIPstrtok(NULL, " ", &nexttok))) || (*sosi->f5 == '$'))
	    sosi->f5 = 0;
      }
      while(FALSE);

      /* check for empty lines */
      is_empty = (sosi->f0 == NULL && sosi->f1 == NULL);
   }
   while (is_marker || is_empty);

   return TRUE;
}


/* Process NAME section.
 *
 * The result is currently ignored.
 */
static
SCIP_RETCODE readName(
      SOSINPUT*             sosi
      )
{
   assert( sosi != NULL );

   /* This has to be the Line with the NAME section. */
   if (!SOSinputReadLine(sosi) || SOSinputField0(sosi) == NULL || strcmp(SOSinputField0(sosi), "NAME"))
   {
      SOSinputSyntaxerror(sosi);
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Problem name   : %s\n", SOSinputField1(sosi));

   /* set SOS section */
   SOSinputSetSection(sosi, SOS_SOSSECTION);

   return SCIP_OKAY;
}



/* Process SOS section. */
static
SCIP_RETCODE readSOS(
   SOSINPUT*             sosi,
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool initial, separate, enforce, check, propagate;
   SCIP_Bool local, modifiable, dynamic, removable;
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS*  cons = NULL;
   int consType = -1;
   int cnt = 0;

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

   /* loop through section */
   while ( SOSinputReadLine(sosi) )
   {
      int type = -1;

      /* check if next is found */
      if ( SOSinputField0(sosi) != NULL )
      {
         if ( ! strcmp(SOSinputField0(sosi), "ENDATA") )
	 {
            SOSinputSetSection(sosi, SOS_ENDATA);
	    break;
	 }
         else
	 {
	    SCIPerrorMessage("next section unkown <%s>.\n", SOSinputField0(sosi));
	    SOSinputSyntaxerror(sosi);
	    return SCIP_PARSEERROR;
	 }
      }
      if ( SOSinputField1(sosi) == NULL && SOSinputField2(sosi) == NULL )
      {
	 SCIPerrorMessage("empty data in a non-comment line.\n");
	 SOSinputSyntaxerror(sosi);
	 return SCIP_PARSEERROR;
      }

      /* check for new SOS set */
      if ( strcmp(SOSinputField1(sosi), "S1") == 0 )
	 type = 1;
      if ( strcmp(SOSinputField1(sosi), "S2") == 0 )
	 type = 2;

      /* add last constraint and create a new one */
      if ( type > 0 )
      {
	 assert( type == 1 || type == 2 );
	 if ( cons != NULL )
	 {
	    /* add last constraint */
	    SCIP_CALL( SCIPaddCons(scip, cons) );
	    SCIPdebugMessage("(line %d) added constraint <%s>: ", sosi->lineno, SCIPconsGetName(cons));
	    SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL, NULL) ) );
	    SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	 }
	 /* create new name, since we do not get a name from the file */
	 (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "SOS%d", ++cnt);

	 /* create new SOS1 constraint */
	 if ( type == 1 )
	 {
	    /* we do not know the name of the constraint */
	    SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
					  local, modifiable, dynamic, removable) );
	 }
	 else
	 {
	    assert( type == 2 );
	    SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
					  local, modifiable, dynamic, removable) );
	 }
	 consType = type;
	 SCIPdebugMessage("created constraint <%s> of type %d.\n", name, type);
	 /* note: we ignore the priorities! */
      }
      else
      {
	 /* otherwise we are in the section given variables */
	 SCIP_VAR* var = NULL;
	 SCIP_Real weight = 0.0;
	 char* endptr;

	 if ( consType != 1 && consType != 2 )
	 {
	    SCIPerrorMessage("missing SOS type specification.\n");
	    SOSinputSyntaxerror(sosi);
	    return SCIP_PARSEERROR;
	 }

	 /* get variable */
	 var = SCIPfindVar(scip, SOSinputField1(sosi));
	 if ( var == NULL )
	 {
	    SCIPerrorMessage("variable <%s> unkown - is correct problem loaded?\n", SOSinputField1(sosi));
	    SOSinputSyntaxerror(sosi);
	    return SCIP_PARSEERROR;
	 }

	 /* get weight */
	 weight = strtod(SOSinputField2(sosi), &endptr);
	 if ( endptr == SOSinputField2(sosi) || *endptr != '\0' )
	 {
	    SCIPerrorMessage("weight for variable <%s> not specified.\n", SOSinputField1(sosi));
	    SOSinputSyntaxerror(sosi);
	    return SCIP_PARSEERROR;
	 }

	 /* add variable and weight */
	 assert( consType == 1 || consType == 2 );
	 switch (consType)
	 {
	 case 1: SCIP_CALL( SCIPaddVarSOS1(scip, cons, var, weight) ); break;
	 case 2: SCIP_CALL( SCIPaddVarSOS2(scip, cons, var, weight) ); break;
	 default: abort(); /* should not happen */
	 }
	 SCIPdebugMessage("added variable <%s> with weight %g.\n", SCIPvarGetName(var), weight);

	 /* check other fields */
	 if ( ( SOSinputField3(sosi) != NULL && SOSinputField3(sosi) != '\0' ) ||
	      ( SOSinputField4(sosi) != NULL && SOSinputField4(sosi) != '\0' ) ||
	      ( SOSinputField5(sosi) != NULL && SOSinputField5(sosi) != '\0' ) )
	 {
	    SCIPwarningMessage("ignoring data in fields 3-5 <%s> <%s> <%s>.\n",
			       SOSinputField3(sosi), SOSinputField4(sosi), SOSinputField5(sosi));
	 }
      }
   }

   return SCIP_OKAY;
}

/* Read LP in "SOS File Format" */
static
SCIP_RETCODE readSOSFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;
   SOSINPUT* sosi;
   SCIP_Bool error;

   assert(scip != NULL);
   assert(filename != NULL);

   fp = SCIPfopen(filename, "r");
   if (fp == NULL)
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   SCIP_CALL( SOSinputCreate(scip, &sosi, fp) );

   if (SOSinputSection(sosi) == SOS_NAME)
   {
      SCIP_CALL( readName(sosi) );
   }
   if (SOSinputSection(sosi) == SOS_SOSSECTION)
   {
      SCIP_CALL( readSOS(sosi, scip) );
   }
   if (SOSinputSection(sosi) != SOS_ENDATA)
      SOSinputSyntaxerror(sosi);

   SCIPfclose(fp);

   error = SOSinputHasError(sosi);

   SOSinputFree(scip, &sosi);

   if( error )
      return SCIP_PARSEERROR;
   else
      return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSOS NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSOS)
{  /*lint --e{715}*/
   assert( reader != NULL );
   assert( strcmp(SCIPreaderGetName(reader), READER_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   if ( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPwarningMessage("reading of solution file is only possible after a problem was created.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
   {
      SCIPwarningMessage("reading of solution file is only possible in problem creating stage.\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( readSOSFile(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteSOS NULL


/*
 * SOS file reader specific interface methods
 */

/** includes the sos file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSOS(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   assert( scip != NULL );

   /* create sos reader data */
   readerdata = NULL;

   /* include sos reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeSOS, readerReadSOS, readerWriteSOS, readerdata) );

   return SCIP_OKAY;
}
