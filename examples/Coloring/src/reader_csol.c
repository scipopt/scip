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

/**@file   reader_csol.c
 * @brief  file reader and writer for vertex coloring solutions
 * @author Gerald Gamrath
 *
 * This file implements the reader and writer for coloring solution files.
 *
 * These files have the following structure:@n The first line contains the name of the problem, the
 * number of colors used in the solution, and - optional - the name of the algorithm that computed
 * this solution.  The second line lists the colors of the nodes, separated by spaces. It is sorted
 * increasingly by the node indices. The numbers for the colors start with 0.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "reader_csol.h"
#include "reader_col.h"
#include "probdata_coloring.h"


#define READER_NAME             "csolreader"
#define READER_DESC             "file reader which reads and writes csol-files"
#define READER_EXTENSION        "csol"

#define COL_MAX_LINELEN 65535



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
  while ( isspace(**s) )
    ++(*s);
  /* read number */
  tmp = atol(*s);
  /* skip whitespaces */
  while ( (**s != 0) && (!isspace(**s)) )
    ++(*s);
  return tmp;
}


/* put your local methods here, and declare them static */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCsol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
 
   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCsol)
{   
   SCIP_FILE* fp;               /* file-reader */
   char buf[COL_MAX_LINELEN];   /* maximal length of line */
   char* char_p;
   char* solprobname;
   const char* probname;

   SCIP_Bool correctinstance;
   TCLIQUE_GRAPH* graph;
   SCIP_VAR* var;
   SCIP_CONS** constraints;     

   int** sets;
   int* setlengths;
   int nsets;
   int setindex;

   int i;
   int j;
   int k;
   int color;
   int node;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(filename != NULL);
   *result = SCIP_SUCCESS;

   if ( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPerrorMessage("Please read in problem before reading the solution!\n");
      return SCIP_OKAY;
   }

   if (NULL == (fp = SCIPfopen(filename, "r")))
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }
   
   /* Read out the name of the problem belonging to this solution*/
   SCIPfgets(buf, sizeof(buf), fp);
   i = 1;
   while ( !isspace(buf[i]) )
   {
      i++;
   }
   SCIPallocMemoryArray(scip, &solprobname, i+2);
   strncpy(solprobname, &buf[0], i);
   solprobname[i]= '\0';

   printf("Reading solution for %s...\n", solprobname);

   /* get the name of the current problem */
   probname = SCIPgetProbName(scip);

   /* check whether the solution belongs to the current problem */
   correctinstance = TRUE;
   for ( j = 0; j <= i; j++ )
   {
      if ( solprobname[j] != probname[j] )
      {
         correctinstance = FALSE;
      }
   }
   if ( !correctinstance )
   {
      SCIPerrorMessage("The selected solution file doesn't belong to the current problem!\n");
      return SCIP_OKAY;
   }

   /* get the graph of the current problem */
   graph = COLORprobGetGraph(scip);
   assert(graph != NULL);

   /* read out number of colors */
   char_p = &buf[i];
   nsets = getNextNumber(&char_p);
   assert(nsets > 0);

   /* allocate memory for the stable sets */
   SCIP_CALL( SCIPallocBufferArray(scip, &sets, nsets) );
   SCIP_CALL( SCIPallocBufferArray(scip, &setlengths, nsets) );
   for ( i = 0; i < nsets; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(sets[i]), COLORprobGetOriginalNNodes(scip)+1) );
      setlengths[i] = 0;
   }

   /* read out the colors for the nodes */
   SCIPfgets(buf, sizeof(buf), fp);
   char_p = &buf[0];
   for ( i = 0; i < COLORprobGetOriginalNNodes(scip); i++ )
   {
      color = getNextNumber(&char_p);
      sets[color][setlengths[color]] = i;
      sets[color][setlengths[color]+1] = -1;
      setlengths[color]++;
   }
   
   /* the given coloring is a coloring for the original graph, now transform it into a coloring for the transformed graph */
   for ( i = 0; i < nsets; i++ )
   {
      j = 0;
      k = 0;
      while ( sets[i][j] != -1 )
      {
         node = COLORprobGetNewNodeForOriginalNode(scip, sets[i][j]);
         if ( node == -1 )
         {
            j++;
         }
         else
         {
            sets[i][k] = node;
            setlengths[i] = k+1;
            k++;
            j++;
         }
      }
      while ( k < j )
      {
         sets[i][k] = -1;
         k++;
      }
   }

   printf("testing validity...\n");
   /* check solution */
   for ( i = 0; i < nsets; i++ )
   {
      for ( j = 0; j < setlengths[i]; j++ )
      {
         for ( k = j+1; k < setlengths[i]; k++ )
         {
            if ( tcliqueIsEdge(graph, sets[i][j], sets[i][k]) )
            {
               SCIPerrorMessage("The solution is not valid!\n");
               return SCIP_OKAY;
            }
         }
      }
   }
   printf("valid!\n");

   /* get the node-constraits */
   constraints = COLORprobGetConstraints(scip);
   assert(constraints != NULL);
   /* try to add nodes to the stable sets */
   for ( i = 0; i < nsets; i++ )
   {
      for ( node = 0; node < COLORprobGetNNodes(scip); node++ )
      {
         for ( j = 0; j < setlengths[i]; j++ )
         {
            if ( sets[i][j] == node )
            {
               break;
            }
            if ( tcliqueIsEdge(graph, sets[i][j], node) )
            {
               break;
            }
         }
         if ( j == setlengths[i] )
         {
            sets[i][setlengths[i]] = node;
            sets[i][setlengths[i]+1] = -1;
            setlengths[i]++;
         }
      }
   }

   /* sort the sets and add them to the problem, creating one variable for each set */
   for ( i = 0; i < nsets; i++ )
   {
      SCIPsortDownInt(sets[i], setlengths[i]);
      COLORprobAddNewStableSet(scip, sets[i], setlengths[i], &setindex);
      assert(setindex == i);

      SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0, 1, 1, SCIP_VARTYPE_BINARY, 
            TRUE, FALSE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setindex) );

      SCIP_CALL( COLORprobAddVarForStableSet(scip, setindex, var) );
      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

      /* add variable to node constraints of nodes in the set */
      for ( j = 0; j < setlengths[i]; j++ )
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, constraints[sets[i][j]], var) );
      }

   }
   
   
   /* free memory for the stable sets */
   for ( i = nsets-1; i >= 0; i-- )
   {
      SCIPfreeBufferArray(scip, &(sets[i]));
   }
   SCIPfreeBufferArray(scip, &setlengths);
   SCIPfreeBufferArray(scip, &sets);
   SCIPfreeMemoryArray(scip, &solprobname);

   return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCsol)
{  
   SCIP_SOL* sol;
   SCIP_Bool colorpossible;
   TCLIQUE_GRAPH* oldgraph;
   int** sets;
   int* nsetelements;
   int nsets;
   int nnodes;
   int i;
   int j;
   int actcolor;
   int node;
   int* originalnodes;
   int* deletednodes;
   int* firstedge;
   int* lastedge;
   int* colors;

   *result = SCIP_DIDNOTRUN;

   /* get the data of the original graph, the preprocessing information and the array stable sets in the preprocessed graph */
   nnodes = COLORprobGetOriginalNNodes(scip);
   originalnodes = COLORprobGetOriginalNodesForNewNodes(scip);
   assert(originalnodes != NULL);
   deletednodes = COLORprobGetDeletedNodes(scip);
   assert(deletednodes != NULL);
   oldgraph = COLORprobGetOriginalGraph(scip);
   COLORprobGetStableSets(scip, &sets, &nsetelements, &nsets);
   assert(sets != NULL && nsetelements != NULL);
   
   /* get the solution */
   sol = SCIPgetBestSol(scip);

   /* create array for the colors of the nodes and initialize it with -1 */
   SCIP_CALL( SCIPallocBufferArray(scip, &colors, nnodes) );
   for ( i = 0; i < nnodes; i++ )
   {
      colors[i] = -1;
   }

   /* for all stable sets in the solution, color all nodes, that are in the set and not yet colored with the same, new color */
   actcolor = 0;
   for ( i = 0; i < nsets; i++ )
   {
      if ( SCIPgetSolVal(scip, sol, COLORprobGetVarForStableSet(scip, i)) > 0 )
      {
         assert(SCIPgetSolVal( scip, sol, COLORprobGetVarForStableSet(scip, i)) == 1);
         for ( j = 0; j < nsetelements[i]; j++ )
         {
            if ( colors[originalnodes[sets[i][j]]] == -1 )
            {
               colors[originalnodes[sets[i][j]]] = actcolor;
            }
         }
         actcolor++;
      }
   }

   /* set i to the index of the last node deleted during preprocessing */
   i = COLORprobGetOriginalNNodes(scip)-1;
   while ( deletednodes[i] == -1 )
   {
      i--;
   }

   /*compute colors for nodes deleted during preprocessing */
   while ( i >= 0 )
   {
      node = deletednodes[i];
      j = 0;
      while ( colors[node] == -1 )
      {
         colorpossible = TRUE;
         firstedge = tcliqueGetFirstAdjedge(oldgraph, node);
         lastedge = tcliqueGetLastAdjedge(oldgraph, node);
         while ( firstedge <= lastedge )
         {
            if ( colors[*firstedge] == j )
            {
               colorpossible = FALSE;
               break;
            }
            firstedge++;
         }
         if ( colorpossible == TRUE )
         {
            colors[node] = j;
         }
         else
         {
            j++;
         }
      }
      i--;
   }
   
   SCIPinfoMessage(scip, file, "%s %d generated by ColumnGenerationColoring\n", name, actcolor);
   for ( i = 0; i < nnodes; i++ )
   {
      SCIPinfoMessage(scip, file, "%d ", colors[i]);
   }

   SCIPfreeBufferArray(scip, &colors);

   *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the csol file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create csol reader data */
   readerdata = NULL;
   
   /* include csol reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION,
         readerdata) );

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCsol) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCsol) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCsol) );

   return SCIP_OKAY;
}
