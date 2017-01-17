/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_col.c
 * @brief  file reader for vertex coloring instances
 * @author Gerald Gamrath
 *
 * This file implements the reader for vertex coloring problems in DIMACS standard format.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "reader_spa.h"
#include "probdata_spa.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_indicator.h"

#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>



#define READER_NAME             "spareader"
#define READER_DESC             "file reader for a .spa-file representing an adjacency matrix for a clustering problem"
#define READER_EXTENSION        "spa"

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
   while ( isspace(**s) )
      ++(*s);
   /* read number */
   tmp = atof(*s);
   /* skip whitespaces */
   while ( (**s != 0) && (!isspace(**s)) )
      ++(*s);
   return tmp;
}

/** read LP in Spa File Format.
 * That means first line is "p edges nbins nedges(in pmatrix) ncluster".
 * Then a list of all edges in pmatrix with <startnode> <endnode> <weight>
 * Then a line "sd" that specifies the stationairy distribution vector starts
 * Then all weights in the stationary distribution vector
*/
static
SCIP_RETCODE readSpa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
)
{
   SCIP_FILE* fp;               /* file-reader */
   char buf[COL_MAX_LINELEN];   /* maximal length of line */
   int nbins;
   int ncluster;
   char* char_p;
   SCIP_Real** cmatrix;
   int i;
   int col;

   assert(scip != NULL);
   assert(filename != NULL);

   if ( NULL == (fp = SCIPfopen(filename, "r")) )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }

   /* Get problem name from filename and save it */
   if( SCIPfgets(buf, (int) sizeof(buf), fp) == NULL)
      return SCIP_READERROR;

   while( !SCIPfeof(fp) && (buf[0] != '#' || buf[2] != 'p') )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
   }
   /* no graph information in file! */
   if ( SCIPfeof(fp) )
   {
      SCIPerrorMessage("Error! Could not find line starting with 'p'.\n");
      return SCIP_READERROR;
   }
   char_p = &buf[3];
   /* if line reads 'edges' (non-standard!), instead of 'edge'. */

   /* read out number of nodes and edges, the pointer char_p will be changed */

   nbins = (int) getNextNumber(&char_p);
   ncluster = (int) getNextNumber(&char_p);

   if ( nbins <= 0 )
   {
      SCIPerrorMessage("Number of bins must be positive!\n");
      return SCIP_READERROR;
   }

   if ( ncluster <= 0 || nbins <= ncluster )
   {
         SCIPerrorMessage("Number of cluster must be positive and smaller than number of bins!\n");
         return SCIP_READERROR;
   }
   /* create cmatrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &cmatrix, nbins) );
   for( i = 0; i < nbins; i++)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(cmatrix[i]), nbins) ); /*lint !e866*/
   }
   /* fill array for edges */
   i = 0;
   while ( !SCIPfeof(fp) && i < nbins )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
      char_p = &buf[0];
      for( col = 0; col < nbins; ++col )
      {
         cmatrix[i][col] = (SCIP_Real) getNextNumber(&char_p);
      }

      if( i >= nbins )
      {
         SCIPerrorMessage( "more lines than expected: expected %d many, but got already %d'th (non-duplicate) edge", nbins, i+1 );
         return SCIP_READERROR;
      }
      i++;
   }

   /* create problem data */
   SCIP_CALL( SCIPcreateProbSpa(scip, filename, nbins, ncluster, cmatrix) );
   SCIPinfoMessage(scip, NULL, "Original problem: \n");

   for ( i = nbins - 1; i >= 0; i-- )
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
SCIP_DECL_READERCOPY(readerCopySpa)
{  /*lint --e{715}*/
   assert( scip != NULL);
   assert(reader != NULL);
   assert(strcmp( SCIPreaderGetName(reader), READER_NAME) == 0);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSpa)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp( SCIPreaderGetName(reader), READER_NAME) == 0);
   assert( scip != NULL);
   assert(result != NULL);

   SCIP_CALL( readSpa( scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * col file reader specific interface methods
 */

/** includes the spa file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSpa(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create col reader data */
   readerdata = NULL;

   /* include col reader */
   SCIP_CALL( SCIPincludeReaderBasic( scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION,
      readerdata) );

   SCIP_CALL( SCIPsetReaderCopy( scip, reader, readerCopySpa) );
   SCIP_CALL( SCIPsetReaderRead( scip, reader, readerReadSpa ) );

   SCIP_CALL( SCIPaddRealParam(scip,"scale_coherence","factor to scale the cohrence in the target function", NULL, FALSE, 0.001, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "ncluster", "the amount of clusters allowed", NULL, FALSE, 3, 1, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "model", "the model variant", NULL, FALSE, 's', "se", NULL, NULL) );

   return SCIP_OKAY;
}
