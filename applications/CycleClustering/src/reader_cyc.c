/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_cyc.c
 * @brief  file reader for cycle clustering instances
 * @author Leon Eifler
 *
 * This file implements the reader for the cycle clustering problem. The data is read from a matrix, entries separated
 * by whitespace. The first line in the file has to be of the form "# p nstates ncluster",
 * where nstates is the size of the matrix and ncluster is the number of clusters that should be used.
 * The file has to have the ending ".cyc" to be recognized by the reader.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "reader_cyc.h"

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "probdata_cyc.h"

#define READER_NAME             "cycreader"
#define READER_DESC             "file reader for a .cyc-file with a transition matrix for a cycle clustering problem"
#define READER_EXTENSION        "cyc"

#define COL_MAX_LINELEN 10000


/*
 * Local methods
 */

/** get next number from string s */
static
SCIP_Real getNextNumber(
   char**                s                   /**< pointer to the pointer of the current position in the string */
   )
{
   SCIP_Real tmp;

   /* skip whitespaces */
   while( isspace(**s) )
      ++(*s);
   /* read number */
   tmp = atof(*s);
   /* skip whitespaces */
   while( (**s != 0) && (!isspace(**s)) )
      ++(*s);

   return tmp;
}

/** read LP in Cyc File Format.
 * That means first line is "p edges nbins ncluster".
 * Then a matrix with whitespace-separated entries of size nbins x nbins
*/
static
SCIP_RETCODE readCyc(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;                            /* file-reader */
   char buf[COL_MAX_LINELEN];                /* maximal length of line */
   char* char_p;                             /* current char */
   SCIP_Real** cmatrix;                      /* transition matrix */
   int nbins;                                /* number of states */
   int ncluster;                             /* number of clusters */
   int i;
   int col;

   assert(scip != NULL);
   assert(filename != NULL);

   if( NULL == (fp = SCIPfopen(filename, "r")) )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }

   /* Get problem name from filename and save it */
   if( SCIPfgets(buf, (int) sizeof(buf), fp) == NULL )
      return SCIP_READERROR;

   while( !SCIPfeof(fp) && (buf[0] != '#' || buf[2] != 'p') )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
   }

   /* no graph information in file! */
   if( SCIPfeof(fp) )
   {
      SCIPerrorMessage("Error! Could not find line starting with 'p'.\n");
      return SCIP_READERROR;
   }
   char_p = &buf[3];

   /* read out number of nodes and edges, the pointer char_p will be changed */
   nbins = (int) getNextNumber(&char_p);
   ncluster = (int) getNextNumber(&char_p);

   if( nbins <= 0 )
   {
      SCIPerrorMessage("Number of bins must be positive!\n");
      return SCIP_READERROR;
   }

   if( ncluster <= 0 || nbins <= ncluster )
   {
         SCIPerrorMessage("Number of cluster must be positive and smaller than number of bins!\n");
         return SCIP_READERROR;
   }

   /* create cmatrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &cmatrix, nbins) );
   for( i = 0; i < nbins; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(cmatrix[i]), nbins) ); /*lint !e866*/
   }

   /* fill array the cmatrix */
   i = 0;
   while( !SCIPfeof(fp) && i < nbins )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
      char_p = &buf[0];
      for( col = 0; col < nbins; ++col )
      {
         cmatrix[i][col] = (SCIP_Real) getNextNumber(&char_p);
      }

      if( i >= nbins )
      {
         SCIPerrorMessage( "more lines than expected: expected %d many, but got already %d'th (non-duplicate) edge",
            nbins, i+1 );

         return SCIP_READERROR;
      }

      i++;
   }

   /* create problem data */
   SCIP_CALL( SCIPcreateProbCyc(scip, filename, nbins, ncluster, cmatrix) );

   SCIPinfoMessage(scip, NULL, "Original problem: \n");

   for( i = nbins - 1; i >= 0; i-- )
   {
      SCIPfreeMemoryArray(scip, &(cmatrix[i]));
   }

   SCIPfreeMemoryArray(scip, &cmatrix);

   SCIPfclose(fp);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCyc)
{
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp( SCIPreaderGetName(reader), READER_NAME) == 0);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCyc)
{
   assert(reader != NULL);
   assert(strcmp( SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_CALL( readCyc( scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * cyc file reader specific interface methods
 */

/** includes the cyc file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCyc(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create cyc reader data */
   readerdata = NULL;

   /* include cyc reader */
   SCIP_CALL( SCIPincludeReaderBasic( scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION,
      readerdata) );

   SCIP_CALL( SCIPsetReaderCopy( scip, reader, readerCopyCyc) );
   SCIP_CALL( SCIPsetReaderRead( scip, reader, readerReadCyc ) );

   SCIP_CALL( SCIPaddRealParam(scip,"cycleclustering/scale_coherence",
      "factor to scale the cohrence in the target function", NULL, FALSE, 0.001, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddCharParam(scip, "cycleclustering/model",
      "the model variant", NULL, FALSE, 's', "seqt", NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "cycleclustering/usecutselection",
      "true if cut selection should be used in cyc-separators", NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/goodscorefac", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.8, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/badscorefac", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.0, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/goodmaxparall", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.1, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/maxparall", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.5, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/dircutoffdist", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.5, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/efficacyweight", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.4, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/objparalweight", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.1, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "cycleclustering/intsuppweight", "used for cut-selection in cycle-clustering",
      NULL, FALSE, 0.3, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
