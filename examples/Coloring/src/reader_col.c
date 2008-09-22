/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_col.c,v 1.2 2008/09/22 16:21:32 bzfgamra Exp $"

/**@file   reader_col.c
 * @brief  COL file reader
 * @author Gerald Gamrath
 *
 * This file implements the reader for coloring files in
 * DIMACS standard format.
 *
 * Additionally, it provides two sorting functions and a method,
 * which ensures, that all nodes in the graph are covered by
 * at least one stable set.
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
 * Data structures
 */

/** data for col reader */
struct SCIP_ReaderData
{
};


/*
 * Local methods
 */

/** get next number from string s */
static
long getNextNumber(
   char**                s                   /**< pointer to the pointer of the current position in the string */
   )
{
  /* skip whitespaces */
  while ( isspace(**s) )
    ++(*s);
  /* read number */
  long tmp = atol(*s);
  /* skip whitespaces */
  while ( (**s != 0) && (!isspace(**s)) )
    ++(*s);
  return tmp;
}


/** returns a node with maximal weight which is not yet colored, 
 *  used to assure that each node is covered with stable sets after fixing variables to 0 
 */
static
SCIP_Bool hasUncoloredNode(
   TCLIQUE_GRAPH*        graph,              /**< the graph that should be colored */
   SCIP_Bool*            colored             /**< array of booleans, colored[i] == TRUE iff node i is colored */
   )
{
   int i;
   const TCLIQUE_WEIGHT* weights;

   assert(graph != NULL);
   assert(colored != NULL);
   weights = tcliqueGetWeights(graph);
   assert(weights != NULL);

   for ( i = 0; i < tcliqueGetNNodes(graph); i++)
   {
      /* node not yet colored */
      if (!colored[i])
      {
	return TRUE;
      }
   }
   return FALSE;
}


/** computes a stable set with a greedy-method.  attention: the weight of the maximum stable set is not computed! */
static 
SCIP_RETCODE greedyStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< pointer to graph data structure */
   int*                  maxstablesetnodes,  /**< pointer to store nodes of the maximum weight stableset */
   int*                  nmaxstablesetnodes  /**< pointer to store number of nodes in the maximum weight stableset */
   )
{
   const TCLIQUE_WEIGHT* weights; 
   SCIP_Bool indnode;
   int nnodes;
   int i;
   int j;
   int* degrees;
   int* sortednodes;
   int* values;    /* values for sorting the nodes: deg(v)+w(v)*nnodes  */

   assert(scip != NULL);
   assert(graph != NULL);
   assert(maxstablesetnodes != NULL);
   assert(nmaxstablesetnodes != NULL);

   /* get number of nodes */
   nnodes = tcliqueGetNNodes(graph);
   *nmaxstablesetnodes = 0;

   /* get the  degrees and weights for the nodes in the graph */
   degrees = tcliqueGetDegrees(graph);
   weights = tcliqueGetWeights(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &values, nnodes) );   
   SCIP_CALL( SCIPallocBufferArray(scip, &sortednodes, nnodes) );

   /* set values to the nodes which are used for sorting them */
   /* value = degree of the node + weight of the node * number of nodes, therefore the yet colored nodes
      (which have weight 0) have lower values than the not yet colored nodes which have weight 1 */
   for ( i = 0; i < nnodes; i++ )
   {
      sortednodes[i] = i;
      values[i] = degrees[i]+weights[i]*nnodes;
   }

   /* sort the nodes w.r.t. the computed values */
   COLORreaderBubbleSortIntInt(sortednodes, values, nnodes);

   /* insert first node */
   maxstablesetnodes[0] = sortednodes[0];
   (*nmaxstablesetnodes) = 1;
   for ( i = 1; i < nnodes; i++)
   {
      /* check whether node is independent to nodes in the set */
      indnode = TRUE;
      for ( j = 0; j < (*nmaxstablesetnodes); j++ )
      {
         if ( tcliqueIsEdge(graph, sortednodes[i], maxstablesetnodes[j]) )
         {
            indnode = FALSE;
            break;
         }
      }
      if ( indnode == TRUE )
      {
         /* node is independent, thus add it to the set */
         maxstablesetnodes[*nmaxstablesetnodes] = sortednodes[i];
         (*nmaxstablesetnodes) = (*nmaxstablesetnodes)+1;
      }

   }
   SCIPfreeBufferArray(scip, &sortednodes);
   SCIPfreeBufferArray(scip, &values);   
   
   return SCIP_OKAY;
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
   long nedges;
   long nnodes;
   int line_nr;
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
   SCIPfgets(buf, sizeof(buf), fp);
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
   
   SCIPallocMemoryArray(scip, &probname, j-i-4);
   strncpy(probname, &filename[i+1], j-i-5);
   probname[j-i-5]= '\0';

   /* Read until information about graph starts */
   while( !SCIPfeof(fp) && (buf[0] != 'p') )
   {
      SCIPfgets(buf, sizeof(buf), fp);
      line_nr++;
   } 
   /* no graph information in file! */
   if ( SCIPfeof(fp) )
   {
      SCIPerrorMessage("Error! Could not find line starting with 'p'.\n");
      return SCIP_PARSEERROR;
   }
   /* wrong format of the line containig number of nodes and edges */
   if ( buf[2] != 'e' || buf[3] != 'd' || buf[4] != 'g' || buf[5] != 'e' )
   {
      SCIPerrorMessage("Line starting with 'p' must continue with 'edge'!\n");
      return SCIP_PARSEERROR;
   }
   char_p = &buf[6];
   /* if line reads 'edges' (non-standard!), instead of 'edge'. */
   if ( *char_p == 's' )
      ++(char_p);

   /* read out number of nodes and edges, the pointer char_p will be changed */
   nduplicateedges = 0;
   nnodes = getNextNumber(&char_p);
   nedges = getNextNumber(&char_p);
   if ( nnodes <= 0 )
   {
      SCIPerrorMessage("Number of vertices must be positive!\n");
      return SCIP_PARSEERROR;
   }
   if ( nedges < 0 )
   {	  
      SCIPerrorMessage("Number of edges must be nonnegative!\n");
      return SCIP_PARSEERROR;
   }
   /* create array for edges */
   SCIP_CALL( SCIPallocMemoryArray(scip, &edges, nedges) );
   for( i = 0; i < nedges; i++)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(edges[i]), 2) );
   }
   /* fill array for edges */
   SCIPfgets(buf, sizeof(buf), fp);
   line_nr++;
   i = 0;
   while ( !SCIPfeof(fp) )
   {
      if ( buf[0] == 'e')
      {
         duplicateedge = FALSE;
         char_p = &buf[2];
         
         begin = getNextNumber(&char_p);
         end = getNextNumber(&char_p);
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
            edges[i][0] = begin;
            edges[i][1] = end;
            assert((edges[i][0] > 0) && (edges[i][0] <= nnodes));
            assert((edges[i][1] > 0) && (edges[i][1] <= nnodes));
            i++;
         }
      }
      SCIPfgets(buf, sizeof(buf), fp);
      line_nr++;
   }
   if ( nduplicateedges > 0 )
   {
      printf("%d duplicate edges!\n", nduplicateedges);
   }
   
   /* create problem data */
   SCIP_CALL( SCIPcreateProbColoring(scip, probname, nnodes, nedges-nduplicateedges, edges) );

   /* create LP */
   SCIPdebugMessage("Erstelle LP...\n");
   COLORprobSetUpArrayOfCons(scip);

   
   /* activate the pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "coloring")) );
   SCIP_CALL( SCIPsetObjIntegral(scip) );
   for ( i = nedges-1; i >= 0; i--)
   {
      SCIPfreeMemoryArray(scip, &(edges[i]));
   }
   SCIPfreeMemoryArray(scip, &edges);
   SCIPfreeMemoryArray(scip, &probname);
   SCIPfclose(fp);

   return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */
#define readerFreeCol NULL

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

  /* create col reader data */
  readerdata = NULL;

  /* include col reader */
  SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
        readerFreeCol, readerReadCol, NULL, readerdata) );


  return SCIP_OKAY;
}


/* ----------------------------------- external methods -------------------------- */

/** creates the initial LP */
SCIP_RETCODE COLORreaderCreateSetsForUncoveredNodes(
   SCIP*                 scip,               /**< SCIP data structure */   
   TCLIQUE_GRAPH*        graph               /**< pointer to graph data structure */
   )
{
   int nnodes;                  /* number of nodes */
   int* maxstablesetnodes;      /* array containig the nodes of the max stable set */
   int nmaxstablesetnodes;      /* number of nodes in stable set */
   int setnumber;               /* number of already found stable sets */
   SCIP_VAR* var;               /* var for the actual stable set */
   SCIP_CONS** constraints;     /* array of added constraints */
   SCIP_Bool* colored;          /* array for marking of yet colored nodes
                                   colored_i = true iff node i is already colored */
   int**              stablesets;
   int*               nstablesetelements;
   int                nstablesets;
   int i;
   int j;
   assert(scip != NULL);
   assert(graph != NULL);

   nnodes = COLORprobGetNNodes(scip);
   assert(nnodes > 0);
   /* get the node-constraits */
   constraints = COLORprobGetConstraints(scip);
   assert(constraints != NULL);
   /* get all yet computed stable sets */
   COLORprobGetStableSets(scip, &stablesets, &nstablesetelements, &nstablesets);
   assert(stablesets != NULL && nstablesetelements != NULL);
   assert(nstablesets >= 0);
   assert(nnodes == tcliqueGetNNodes(graph));

   /* set weight of all nodes to 1 */
   for ( i = 0; i < nnodes; i++)
   {
      tcliqueChangeWeight(graph, i, 1);
   }

   /* allocate memory for arrays */
   SCIP_CALL( SCIPallocBufferArray( scip, &colored, nnodes) );
   SCIP_CALL( SCIPallocBufferArray( scip, &maxstablesetnodes, nnodes) );
   nmaxstablesetnodes = 0;

   /* fill colored-array with FALSE */
   BMSclearMemoryArray(colored, nnodes);
   /* go through all stable sets and set colored to true for nodes in them, also set the weight to 0 for the covered nodes */
   for ( i = 0; i < nstablesets; i++ )
   {
      if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal( COLORprobGetVarForStableSet(scip, i))) 
         && (SCIPgetNNodes(scip) == 0 || SCIPvarIsInLP(COLORprobGetVarForStableSet(scip, i))
            || SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ) )
      {
         for ( j = 0; j < nstablesetelements[i]; j++ )
         {
            colored[stablesets[i][j]] = TRUE;
	    tcliqueChangeWeight(graph, stablesets[i][j], 0);
         }
      }
   }

   /* create maximal Stable Sets until all Nodes are covered */
   while ( hasUncoloredNode(graph, colored) )
   {
      greedyStableSet(scip, graph, maxstablesetnodes, &nmaxstablesetnodes);
      COLORreaderBubbleSortIntInt(maxstablesetnodes, maxstablesetnodes, nmaxstablesetnodes);
      SCIP_CALL( COLORprobAddNewStableSet(scip, maxstablesetnodes, nmaxstablesetnodes, &setnumber) );
      assert(setnumber != -1);
      
      /* create variable for the stable set and add it to SCIP*/
      SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0, SCIPinfinity(scip), 1, SCIP_VARTYPE_INTEGER, 
            TRUE, FALSE, NULL, NULL, NULL, (SCIP_VARDATA*)setnumber) );
      COLORprobAddVarForStableSet(scip, setnumber, var);
      SCIP_CALL( SCIPaddVar(scip, var) );

      for ( i = 0; i < nmaxstablesetnodes; i++ )
      {
         /* add variable to node constraints of nodes in the set */
         SCIP_CALL( SCIPaddCoefLinear(scip, constraints[maxstablesetnodes[i]], var, -1) );
         /* mark node as colored */
         colored[maxstablesetnodes[i]] = TRUE;
         /* set weight of already colored nodes to 0 */
         tcliqueChangeWeight(graph, maxstablesetnodes[i], 0);
      }

   }
   /* free memory */
   SCIPfreeBufferArray(scip, &maxstablesetnodes);
   SCIPfreeBufferArray(scip, &colored);
   return SCIP_OKAY;

}


/** bubble sort of two joint arrays of int, sorted s.t. the second array is in non-increasing order */
void COLORreaderBubbleSortIntInt(
   int*                  values,             /**< int array to be permuted in the same way */
   int*                  keys,               /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   int tmpValue;
   int tmpKey;

   assert(len == 0 || values != NULL);
   assert(len == 0 || keys != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && keys[pos] >= keys[pos+1] )
            pos++;
         if( pos >= lastpos )
            break;
         tmpValue = values[pos];
         tmpKey = keys[pos];
         do
         {
            values[pos] = values[pos+1];
            keys[pos] = keys[pos+1];
            pos++;
         }
         while( pos < lastpos && tmpKey < keys[pos+1] );
         values[pos] = tmpValue;
         keys[pos] = tmpKey;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && keys[pos-1] >= keys[pos] )
            pos--;
         if( pos <= firstpos )
            break;
         tmpValue = values[pos];
         tmpKey = keys[pos];
         do
         {
            values[pos] = values[pos-1];
            keys[pos] = keys[pos-1];
            pos--;
         }
         while( pos > firstpos && keys[pos-1] < tmpKey );
         values[pos] = tmpValue;
         keys[pos] = tmpKey;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

/** bubble sort of two joint arrays of int and reals, sorted s.t. the second array is in non-increasing order */
void COLORreaderBubbleSortIntReal(
   int*                  values,             /**< int array to be permuted in the same way */
   SCIP_Real*            keys,               /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   )
{
   int firstpos;
   int lastpos;
   int pos;
   int sortpos;
   int tmpValue;
   SCIP_Real tmpKey;

   assert(len == 0 || values != NULL);
   assert(len == 0 || keys != NULL);

   firstpos = 0;
   lastpos = len-1;
   while( firstpos < lastpos )
   {
      /* bubble from left to right */
      pos = firstpos;
      sortpos = firstpos;
      while( pos < lastpos )
      {
         while( pos < lastpos && keys[pos] >= keys[pos+1] )
            pos++;
         if( pos >= lastpos )
            break;
         tmpValue = values[pos];
         tmpKey = keys[pos];
         do
         {
            values[pos] = values[pos+1];
            keys[pos] = keys[pos+1];
            pos++;
         }
         while( pos < lastpos && tmpKey < keys[pos+1] );
         values[pos] = tmpValue;
         keys[pos] = tmpKey;
         sortpos = pos;
         pos++;
      }
      lastpos = sortpos-1;

      /* bubble from right to left */
      pos = lastpos;
      sortpos = lastpos;
      while( pos > firstpos )
      {
         while( pos > firstpos && keys[pos-1] >= keys[pos] )
            pos--;
         if( pos <= firstpos )
            break;
         tmpValue = values[pos];
         tmpKey = keys[pos];
         do
         {
            values[pos] = values[pos-1];
            keys[pos] = keys[pos-1];
            pos--;
         }
         while( pos > firstpos && keys[pos-1] < tmpKey );
         values[pos] = tmpValue;
         keys[pos] = tmpKey;
         sortpos = pos;
         pos--;
      }
      firstpos = sortpos+1;
   }
}

