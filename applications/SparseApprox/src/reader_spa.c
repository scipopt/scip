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

#define COL_MAX_LINELEN 1024


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
   int nedges;
   int nbins;
   char* char_p;
   SCIP_Real** edges;
   SCIP_Real* sd;
   SCIP_Real begin;
   SCIP_Real end;
   SCIP_Real weight;
   int i;

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

   while( !SCIPfeof(fp) && (buf[0] != 'p') )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
   }
   /* no graph information in file! */
   if ( SCIPfeof(fp) )
   {
      SCIPerrorMessage("Error! Could not find line starting with 'p'.\n");
      return SCIP_READERROR;
   }
   if ( buf[2] != 'e' || buf[3] != 'd' || buf[4] != 'g' || buf[5] != 'e' )
   {
      SCIPerrorMessage("Line starting with 'p' must continue with 'edge'!\n");
      return SCIP_READERROR;
   }
   char_p = &buf[6];
   /* if line reads 'edges' (non-standard!), instead of 'edge'. */
   if ( *char_p == 's' )
      ++(char_p);

   /* read out number of nodes and edges, the pointer char_p will be changed */

   nbins = (int) getNextNumber(&char_p);
   nedges = (int) getNextNumber(&char_p);

   if ( nbins <= 0 )
   {
      SCIPerrorMessage("Number of bins must be positive!\n");
      return SCIP_READERROR;
   }
   if ( nedges < 0 )
   {	  
      SCIPerrorMessage("Number of edges must be nonnegative!\n");
      return SCIP_READERROR;
   }
   /* create array for edges */
   SCIP_CALL( SCIPallocMemoryArray(scip, &edges, nedges) );
   for( i = 0; i < nedges; i++)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(edges[i]), 3) ); /*lint !e866*/
   }
   /* fill array for edges */
   i = 0;
   while ( !SCIPfeof(fp) && i < nedges )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/


      char_p = &buf[0];
      assert(buf[0] != 's');
      begin = (SCIP_Real) getNextNumber(&char_p);
      end = (SCIP_Real) getNextNumber(&char_p);
      if( (int) begin > nbins || (int) end > nbins )
      {
         SCIPerrorMessage("Invalid edge from %d to %d. Only %d bins were specified", begin, end, nbins);
         return SCIP_READERROR;
      }
      weight = (SCIP_Real) getNextNumber(&char_p);

      if( i >= nedges )
      {
         SCIPerrorMessage( "more edges than expected: expected %d many, but got already %d'th (non-duplicate) edge", nedges, i+1 );
         return SCIP_READERROR;
      }
      edges[i][0] = begin;
      edges[i][1] = end;
      edges[i][2] = weight;
      assert((edges[i][0] > 0) && (edges[i][0] <= nbins));
      assert((edges[i][1] > 0) && (edges[i][1] <= nbins));
      i++;
   }
   /* search for the section with stationairy distribution */
   while( !SCIPfeof(fp) && (buf[0] != 's') )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
   }
   /* no graph information in file! */
   if ( SCIPfeof(fp) )
   {
      SCIPerrorMessage("Error! Could not find line starting with 's'.\n");
      return SCIP_READERROR;
   }
   if ( buf[1] != 'd' )
   {
      SCIPerrorMessage("The section for the stationary distribution must be initialized by sd \n");
      return SCIP_READERROR;
   }
   /* create array for stationary distribution */
   SCIP_CALL( SCIPallocMemoryArray(scip, &sd, nbins) );

   i = 0;
   while( !SCIPfeof(fp) && i < nbins )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
      char_p = &buf[0];

      sd[i] = (SCIP_Real) getNextNumber(&char_p);
      i++;
   }
   if( i != nbins )
   {
      SCIPerrorMessage("error in the stationary distribution part. Expected %d values but found only %d many\n", nbins, i );
      return SCIP_READERROR;
   }

   /* create problem data */
   SCIP_CALL( SCIPcreateProbSpa(scip, filename, nbins, nedges, edges, sd) );
   SCIPinfoMessage(scip, NULL, "Original problem: \n");

   for ( i = nedges-1; i >= 0; i-- )
   {
      SCIPfreeMemoryArray(scip, &(edges[i]));
   }

   SCIPfreeMemoryArray(scip, &sd);
   SCIPfreeMemoryArray(scip, &edges);
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

   return SCIP_OKAY;
}
