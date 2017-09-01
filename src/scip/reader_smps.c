/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_smps.c
 * @brief  SMPS file reader (MPS format of the core problem for stochastic programs)
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/reader_smps.h"
#include "scip/reader_cor.h"
#include "scip/reader_tim.h"
#include "scip/reader_sto.h"

#define READER_NAME             "smpsreader"
#define READER_DESC             "file reader for core problem of stochastic programs in the SMPS file format"
#define READER_EXTENSION        "smps"

#define SMPS_MAX_LINELEN  1024
#define BLANK         ' '

/** smps input structure */
struct SmpsInput
{
   SCIP_FILE*            fp;
   int                   lineno;
   SCIP_Bool             haserror;
   char                  buf[SMPS_MAX_LINELEN];
   const char*           f0;
   const char*           f1;
};
typedef struct SmpsInput SMPSINPUT;


/** creates the smps input structure */
static
SCIP_RETCODE smpsinputCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SMPSINPUT**            smpsi,               /**< smps input structure */
   SCIP_FILE*            fp                  /**< file object for the input file */
   )
{
   assert(smpsi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, smpsi) );

   (*smpsi)->fp          = fp;
   (*smpsi)->lineno      = 0;
   (*smpsi)->haserror    = FALSE;
   (*smpsi)->buf     [0] = '\0';
   (*smpsi)->f0          = NULL;
   (*smpsi)->f1          = NULL;


   return SCIP_OKAY;
}

/** free the smps input structure */
static
void smpsinputFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SMPSINPUT**            smpsi                /**< smps input structure */
   )
{
   SCIPfreeBlockMemory(scip, smpsi);
}

/** return the current value of field 0 */
static
const char* smpsinputField0(
   const SMPSINPUT*       smpsi                /**< smps input structure */
   )
{
   assert(smpsi != NULL);

   return smpsi->f0;
}

/** fill the line from \p pos up to column 80 with blanks. */
static
void clearFrom(
   char*                 buf,                /**< buffer to clear */
   unsigned int          pos                 /**< position to start the clearing process */
   )
{
   unsigned int i;

   for(i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/** read a smps format data line and parse the fields. */
static
SCIP_Bool smpsinputReadLine(
   SMPSINPUT*             smpsi                /**< smps input structure */
   )
{
   unsigned int len;
   unsigned int i;
   SCIP_Bool is_marker;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      smpsi->f0 = smpsi->f1 = 0;
      is_marker = FALSE;

      /* Read until we have not a comment line. */
      do
      {
         smpsi->buf[SMPS_MAX_LINELEN-1] = '\0';
         if( NULL == SCIPfgets(smpsi->buf, (int) sizeof(smpsi->buf), smpsi->fp) )
            return FALSE;
         smpsi->lineno++;
      }
      while( *smpsi->buf == '*' );

      /* Normalize line */
      len = (unsigned int) strlen(smpsi->buf);

      for( i = 0; i < len; i++ )
         if( (smpsi->buf[i] == '\t') || (smpsi->buf[i] == '\n') || (smpsi->buf[i] == '\r') )
            smpsi->buf[i] = BLANK;

      if( len < 80 )
         clearFrom(smpsi->buf, len);

      SCIPdebugMessage("line %d: <%s>\n", smpsi->lineno, smpsi->buf);

      assert(strlen(smpsi->buf) >= 80);

      /* Look for new section */
      if( *smpsi->buf != BLANK )
      {
         smpsi->f0 = SCIPstrtok(&smpsi->buf[0], " ", &nexttok);

         assert(smpsi->f0 != 0);

         smpsi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      /* check for empty lines */
      is_empty = (smpsi->f0 == NULL && smpsi->f1 == NULL);
   }
   while( is_marker || is_empty );

   return TRUE;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySmps)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderSmps(scip) );

   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
/**! [SnippetReaderFreeSmps] */
static
SCIP_DECL_READERFREE(readerFreeSmps)
{
   return SCIP_OKAY;
}
/**! [SnippetReaderFreeSmps] */


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSmps)
{  /*lint --e{715}*/
   SCIP_FILE* fp;
   SMPSINPUT* smpsi;
   SCIP_RETCODE retcode = SCIP_OKAY;

   char newfilename[SCIP_MAXSTRLEN];
   char* fromlastslash;
   char* parent;

   assert(scip != NULL);
   assert(filename != NULL);

   fromlastslash = strrchr(filename, '/');

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &parent, filename, strlen(filename) - (strlen(fromlastslash) - 1)) );

   fp = SCIPfopen(filename, "r");
   if( fp == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   SCIP_CALL( smpsinputCreate(scip, &smpsi, fp) );

   while( smpsinputReadLine(smpsi) )
   {
      (void) SCIPsnprintf(newfilename, SCIP_MAXSTRLEN, "%s%s", parent, smpsinputField0(smpsi));
      SCIPinfoMessage(scip, NULL, "read problem <%s>\n", newfilename);
      SCIPinfoMessage(scip, NULL, "============\n");
      SCIP_CALL_TERMINATE( retcode, SCIPreadProb(scip, newfilename, NULL), TERMINATE );
      SCIPinfoMessage(scip, NULL, "\n\n");
   }

   SCIPfclose(fp);

 /* cppcheck-suppress unusedLabel */
 TERMINATE:
   smpsinputFree(scip, &smpsi);
   SCIPfreeBlockMemoryArray(scip, &parent, strlen(parent));

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
      return retcode;

   SCIP_CALL( retcode );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
#if 0
static
SCIP_DECL_READERWRITE(readerWriteSmps)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}
#else
#define readerWriteSmps NULL
#endif


/*
 * reader specific interface methods
 */

/** includes the smps file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSmps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySmps) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeSmps) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSmps) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteSmps) );

   return SCIP_OKAY;
}
