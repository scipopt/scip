/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
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
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "reader_col.h"


#define READER_NAME             "colreader"
#define READER_DESC             "file reader for a .col-file representing a graph that should be colored"
#define READER_EXTENSION        "col"

#define COL_MAX_LINELEN 1024



/*
 * Local methods
 */

/** get next number from string s */
static
long getNextNumber(
   char**                s                   /**< pointer to the pointer of the current position in the string */
   )
{
  long tmp;
  /* skip whitespaces */
  while ( isspace((unsigned char)**s) )
    ++(*s);
  /* read number */
  tmp = atol(*s);
  /* skip whitespaces */
  while ( (**s != 0) && (!isspace((unsigned char)**s)) )
    ++(*s);
  return tmp;
}

/** read LP in "COL File Format" */
static
SCIP_RETCODE readCol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;               /* file-reader */
   char buf[COL_MAX_LINELEN];   /* maximal length of line */
   int nedges;
   int nnodes;
   char* char_p;
   char* probname;
   int** edges;
   int i;
   int j;
   int begin;
   int end;
   int nduplicateedges;
   SCIP_Bool duplicateedge;


   assert(scip != NULL);
   assert(filename != NULL);

   if (NULL == (fp = SCIPfopen(filename, "r")))
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }

   /* Get problem name from filename and save it */
   if( SCIPfgets(buf, (int) sizeof(buf), fp) == NULL)
      return SCIP_READERROR;

   i = 1;
   while ( (filename[i] != '/') && (filename[i] != '\0') )
   {
      i++;
   }
   if ( filename[i] != '/' )
   {
      j = i;
      i = -1;
   }
   else
   {
      j = i+1;
      while ( filename[i] == '/' && filename[j] != '\0' )
      {
         j = i+1;
         while ( filename[j] != '\0' )
         {
            j++;
            if ( filename[j] == '/' )
            {
               i = j;
               break;
            }
         }
      }
   }

   if( j-i-4 <= 0 )
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocBufferArray(scip, &probname, (j-i-4)) );
   (void) strncpy(probname, &filename[i+1], (j-i-5)); /*lint !e732 !e776*/
   probname[j-i-5]= '\0';

   /* Read until information about graph starts */
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

   /* wrong format of the line containig number of nodes and edges */
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
   nduplicateedges = 0;
   nnodes = (int) getNextNumber(&char_p);
   nedges = (int) getNextNumber(&char_p);

   if ( nnodes <= 0 )
   {
      SCIPerrorMessage("Number of vertices must be positive!\n");
      return SCIP_READERROR;
   }

   if ( nedges < 0 )
   {
      SCIPerrorMessage("Number of edges must be nonnegative!\n");
      return SCIP_READERROR;
   }

   /* create array for edges */
   SCIP_CALL( SCIPallocBufferArray(scip, &edges, nedges) );
   for( i = 0; i < nedges; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(edges[i]), 2) ); /*lint !e866*/
   }

   /* fill array for edges */
   i = 0;
   while ( !SCIPfeof(fp) )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
      if ( buf[0] == 'e')
      {
         duplicateedge = FALSE;
         char_p = &buf[2];

         begin = (int) getNextNumber(&char_p);
         end = (int) getNextNumber(&char_p);
         for ( j = 0; j < i; j++)
         {
            if ( ((edges[j][0] == begin) && (edges[j][1] == end))
               || ((edges[j][1] == begin) && (edges[j][0] == end)) )
            {
               duplicateedge = TRUE;
               nduplicateedges++;
               break;
            }
         }
         if ( !duplicateedge )
         {
            if( i >= nedges )
            {
               SCIPerrorMessage("more edges than expected: expected %d many, but got already %d'th (non-duplicate) edge", nedges, i+1);
               return SCIP_READERROR;
            }
            edges[i][0] = begin;
            edges[i][1] = end;
            assert((edges[i][0] > 0) && (edges[i][0] <= nnodes));
            assert((edges[i][1] > 0) && (edges[i][1] <= nnodes));
            i++;
         }
      }
   }
   if( i + nduplicateedges != nedges ) /*lint !e845*/
   {
      SCIPerrorMessage("incorrect number of edges: expected %d many, but got %d many\n", nedges, i + nduplicateedges); /*lint !e845*/
      return SCIP_ERROR;
   }

   printf("Read graph: %d nodes, %d edges (%d duplicates)\n", nnodes, nedges, nduplicateedges); /*lint !e845*/

   /* create problem data */
   SCIP_CALL( SCIPcreateProbColoring(scip, probname, nnodes, nedges-nduplicateedges, edges) );

   /* create LP */
   SCIPdebugMessage("Create LP...\n");
   SCIP_CALL( COLORprobSetUpArrayOfCons(scip) );

   /* activate the pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "coloring")) );
   SCIP_CALL( SCIPsetObjIntegral(scip) );
   for ( i = nedges-1; i >= 0; i--)
   {
      SCIPfreeBufferArray(scip, &(edges[i]));
   }
   SCIPfreeBufferArray(scip, &edges);
   SCIPfreeBufferArray(scip, &probname);
   SCIPfclose(fp);

   return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCol)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIP_CALL( readCol(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * col file reader specific interface methods
 */

/** includes the col file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
  SCIP_READERDATA* readerdata;
  SCIP_READER* reader;

  /* create col reader data */
  readerdata = NULL;

  /* include col reader */
  SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

  SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCol) );
  SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCol) );

  return SCIP_OKAY;
}
