/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   graph_save.c
 * @brief  Saving routines for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file includes several saving routines for Steiner problems
 *
 * A list of all interface methods can be found in graph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include "portab.h"
#include "graph.h"
#include "reduce.h"






/*---------------------------------------------------------------------------*/
/*--- Name     : STP Save Graph                                           ---*/
/*--- Function : Writes a graph to a file in STP format.                  ---*/
/*--- Parameter: Graph and Filepointer.                                   ---*/
/*--- Returns  : Nothing.                                                 ---*/
/*---------------------------------------------------------------------------*/
static void stp_save(
   const GRAPH* g,
   FILE*        fp
   )
{
   int   i;

   assert(g  != NULL);
   assert(fp != NULL);

   fprintf(fp, "%8x STP File, STP Format Version %2d.%02d\n",
      STP_FILE_MAGIC, STP_FILE_VERSION_MAJOR, STP_FILE_VERSION_MINOR);

   fprintf(fp, "Section Comment\n");
   fprintf(fp, "Name \"%s\"\n", "noname");
   fprintf(fp, "End\n\n");

   fprintf(fp, "Section Graph\n");
   fprintf(fp, "Nodes %d\n", g->knots);
   fprintf(fp, "Edges %d\n", g->edges / 2);

   for(i = 0; i < g->edges; i += 2)
   {
      if (g->ieat[i] != EAT_FREE)
      {
         assert(g->oeat[i] != EAT_FREE);

         fprintf(fp, "E %d %d %g\n",
            g->tail[i] + 1,
            g->head[i] + 1,
            g->cost[i]);
      }
   }
   fprintf(fp, "End\n\n");
   fprintf(fp, "Section Terminals\n");
   fprintf(fp, "Terminals %d\n", g->terms);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]))
         fprintf(fp, "T %d\n", i + 1);
   }
   fprintf(fp, "End\n\n");
/*
   if (g->flags & GRAPH_HAS_COORDINATES)
   {
      fprintf(fp, "Section Coordinates\n");

      for(i = 0; i < g->knots; i++)
         fprintf(fp, "DD %d %d %d\n", i + 1, g->xpos[i], g->ypos[i]);

      fprintf(fp, "End\n\n");
   }
*/
   fprintf(fp, "EOF\n");
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Beasley Save Graph                                       ---*/
/*--- Function : Writes a graph to a file in Beasleys format.             ---*/
/*--- Parameter: Graph and Filepointer.                                   ---*/
/*--- Returns  : Nothing.                                                 ---*/
/*---------------------------------------------------------------------------*/
static void bea_save(
   const GRAPH* g,
   FILE*        fp)
{
   int i;

   assert(g  != NULL);
   assert(fp != NULL);

   fprintf(fp, "%d %d\n",
      g->knots,
      g->edges / 2);

   for(i = 0; i < g->edges; i += 2)
   {
      if (g->ieat[i] != EAT_FREE)
      {
         assert(g->oeat[i] != EAT_FREE);

         fprintf(fp, "%d %d %g\n",
            g->tail[i] + 1,
            g->head[i] + 1,
            g->cost[i]);
      }
   }
   fprintf(fp, "%d\n", g->terms);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]))
         fprintf(fp, "%d\n", i + 1);
   }
}

/** gets node map */
static
void getOrgNodeToNodeMap(
   const GRAPH* g,
   int* orgToNewNode
   )
{
   const SCIP_Bool isPcMw = graph_pc_isPcMw(g);
   const int nnodes = graph_get_nNodes(g);
   int nodecount = 0;

   for( int k = 0; k < nnodes; ++k )
   {
      if( isPcMw && graph_pc_knotIsDummyTerm(g, k) )
      {
         orgToNewNode[k] = -1;
         continue;
      }

      if( g->grad[k] > 0 || (isPcMw && GT(g->prize[k], 0.0)) )
      {
         orgToNewNode[k] = nodecount;
         nodecount++;
      }
      else
      {
         orgToNewNode[k] = -1;
      }
   }
}



/** Writes node to gml file.
 *  Discriminates between root, terminals and Steiner nodes */
static inline
void gmlWriteNode(
   const GRAPH*          graph,              /**< Graph to be printed */
   int                   k,                  /**< the node */
   FILE*                 file                /**> the gml file */
)
{
   char label[SCIP_MAXSTRLEN];

   if( k == graph->source )
   {
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", k);
      SCIPgmlWriteNode(file, (unsigned int) k, label, "rectangle", "#666666", NULL);
   }
   else if( graph->term[k] == 0 )
   {
      (void) SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal", k);
      SCIPgmlWriteNode(file, (unsigned int) k, label, "circle", "#ff0000", NULL);
   }
   else
   {
      (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node", k);
      SCIPgmlWriteNode(file, (unsigned int)k, label, "circle", "#336699", NULL);
   }
}

/** Writes edge to gml file */
static inline
void gmlWriteEdge(
   const GRAPH*          graph,              /**< Graph to be printed */
   int                   edge,               /**< the edge */
   int                   tail_id,            /**< tail ID for gml file */
   int                   head_id,            /**< head ID for gml file */
   const SCIP_Bool*      edgemark,           /**< Array of (undirected) edges to highlight or NULL */
   FILE*                 file                /**> the gml file */
)
{
   char label[SCIP_MAXSTRLEN];

   assert(graph_edge_isInRange(graph, edge));
   assert(graph_knot_isInRange(graph, tail_id));
   assert(graph_knot_isInRange(graph, head_id));

   (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[edge]);

   if( edgemark != NULL && edgemark[edge / 2] )
      SCIPgmlWriteEdge(file, (unsigned int)tail_id, (unsigned int)head_id, label, "#ff0000");
   else
      SCIPgmlWriteEdge(file, (unsigned int)tail_id, (unsigned int)head_id, label, NULL);
}


/**  gets ratios of remaining nodes/edges */
static
void getReductionRatiosPcMw(
  const GRAPH*        graph,              /**< graph */
  SCIP_Real*          ratio_nodes,        /**< nodes ratio */
  SCIP_Real*          ratio_edges         /**< edges ratio */
)
{
   assert(graph_pc_isPcMw(graph));
   assert(graph->extended);

   graph_pc_getReductionRatios(graph, ratio_nodes, ratio_edges);
}


/**  Call before graph packing! */
static
void getReductionRatios(
   const GRAPH*        graph,              /**< graph */
   SCIP_Real*	        ratio_nodes,        /**< nodes ratio */
   SCIP_Real*	        ratio_edges         /**< edges ratio */
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   int nnodes_real;
   int nedges_real;

   assert(!graph_pc_isPcMw(graph));

   graph_get_nVET(graph, &nnodes_real, &nedges_real, NULL);

   assert(nnodes_real <= nnodes);
   assert(nedges_real <= nedges);
   assert(nedges % 2 == 0);
   assert(nedges_real % 2 == 0);
   assert(nnodes >= 1);
   assert(nedges >= 2);

   *ratio_nodes = (SCIP_Real) nnodes_real / (SCIP_Real) nnodes;
   *ratio_edges = (SCIP_Real) nedges_real / (SCIP_Real) nedges;
}

/*
 * Interface methods
 */



/** Write (append) reduction ratio statistics of current graph to file.
 *  Call before graph packing!*/
void graph_writeReductionRatioStats(
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           probname,           /**< Name of the problem */
   const char*           filename            /**< Name of the output file */
)
{
   FILE* file;
   SCIP_Real ratio_nodes;
   SCIP_Real ratio_edges;

   assert(filename && probname);

   if( graph_pc_isPcMw(graph) )
   {
	   getReductionRatiosPcMw(graph, &ratio_nodes, &ratio_edges);
   }
   else
   {
	   getReductionRatios(graph, &ratio_nodes, &ratio_edges);
   }

   file = fopen(filename, "a+");

   fprintf(file, "%s: %f %f   \n", probname, ratio_nodes, ratio_edges);
   fclose(file);
}


/** Write (append) reduction ratio statistics of current graph to file.
 *  NOTE: Only for testing! Graph will be killed */
// todo to properly, copy graph etc.
SCIP_RETCODE graph_writeReductionRatioStatsLive(
   SCIP*                 scip,               /**< SCIP data */
   GRAPH*                graph,              /**< Graph to be printed */
   const char*           probname            /**< Name of the problem */
)
{
   REDSOL* redsol;
   FILE* file;
   char* filename;
   SCIP_Real ratio_nodes;
   SCIP_Real ratio_edges;
   int minelims;
   const SCIP_Bool useNodeSol = (graph_pc_isPcMw(graph) || graph_typeIsSpgLike(graph));
   SCIP_Real time_start;
   SCIP_Real time_end;

   assert(graph && probname);

   SCIP_CALL( reduce_solInit(scip, graph, useNodeSol, &redsol) );
   SCIP_CALL( SCIPgetIntParam(scip, "stp/minelims", &(minelims)) );

   time_start = clock();
   SCIP_CALL( reduce_exec(scip, graph, redsol, STP_REDUCTION_ADVANCED, minelims, TRUE) );
   time_end = clock();

   SCIP_CALL( SCIPgetStringParam(scip, "stp/redstatsfile", &filename) );
   assert(filename);

   if( graph_pc_isPcMw(graph) )
   {
      getReductionRatiosPcMw(graph, &ratio_nodes, &ratio_edges);
   }
   else
   {
      getReductionRatios(graph, &ratio_nodes, &ratio_edges);
   }

   file = fopen(filename, "a+");

   fprintf(file, "%s: %f %f %f  \n", probname, ratio_nodes, ratio_edges, (time_end - time_start) / CLOCKS_PER_SEC);
   fclose(file);

   reduce_solFree(scip, &redsol);

   return SCIP_OKAY;
}


/*---------------------------------------------------------------------------*/
/*--- Name     : Graph Save                                               ---*/
/*--- Function : Write a graph to a file.                                 ---*/
/*--- Parameter: Graph, Filename, Filetype (Fileformat).                  ---*/
/*--- Returns  : Nothing.                                                 ---*/
/*---------------------------------------------------------------------------*/
void graph_save(
   SCIP* scip,
   const GRAPH* g,
   const char*  filename,
   FILETYPE     type)
{
   const char* msg_writing_s = "Writing Graph to File %s:\n";

   FILE* fp;

   assert(g && scip && filename);
   assert(graph_valid(scip, g));
   assert(strlen(filename) > 0);
   assert((type == FF_BEA) || (type == FF_STP));

   if ((fp = fopen(filename, "w")) == NULL)
      perror(filename);
   else
   {
      printf(msg_writing_s, filename);

      switch(type)
      {
      case FF_BEA :
         bea_save(g, fp);
         break;
      case FF_STP :
         stp_save(g, fp);
         break;
      case FF_PRB :
      case FF_GRD :
      default     :
         /* CONSTCOND */
         assert(FALSE);
      }
      fclose(fp);
   }
}


/** Write (append) reduction statistics of current graph to file.
 *  Call before graph packing!*/
void graph_writeReductionStats(
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           probname,           /**< Name of the problem */
   const char*           filename            /**< Name of the output file */
)
{
   FILE* file;
   int nnodes_real;
   int nedges_real;
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   assert(filename && probname);

   graph_get_nVET(graph, &nnodes_real, &nedges_real, NULL);

   assert(nnodes_real <= nnodes);
   assert(nedges_real <= nedges);
   assert(nedges % 2 == 0);
   assert(nedges_real % 2 == 0);

   file = fopen(filename, "a+");

   fprintf(file, "%s: %d %d     %d %d \n", probname, nnodes, nedges / 2, nnodes_real, nedges_real / 2);
   fclose(file);
}


/** print graph (in undirected form) in GML format */
SCIP_RETCODE graph_writeGml(
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   const SCIP_Bool*      edgemark            /**< Array of (undirected) edges to highlight */
   )
{
   FILE* file;

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "stpgraph.gml", "w");

#ifndef NDEBUG
   for( int e = 0; e < graph->edges; e += 2 )
   {
      assert(graph->tail[e] == graph->head[e + 1]);
      assert(graph->tail[e + 1] == graph->head[e]);
   }
#endif

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   for( int k = 0; k < graph->knots; k++ )
   {
      if( graph->grad[k] > 0 || graph->source == k )
      {
         gmlWriteNode(graph, k, file);
      }
   }

   /* write all edges (undirected) */
   for( int e = 0; e < graph->edges; e += 2 )
   {
      if( graph->oeat[e] != EAT_FREE )
      {
         gmlWriteEdge(graph, e, graph->tail[e], graph->head[e], edgemark, file);
      }
   }

   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   fclose(file);

   return SCIP_OKAY;
}


/** prints subgraph of given graph (in undirected form) in GML format */
SCIP_RETCODE graph_writeGmlSub(
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   const SCIP_Bool*      subnodesmark        /**< Array to mark the nodes of the subgraph*/
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "stpgraph.gml", "w");

#ifndef NDEBUG
   for( e = 0; e < graph->edges; e += 2 )
   {
      assert(graph->tail[e] == graph->head[e + 1]);
      assert(graph->tail[e + 1] == graph->head[e]);
   }
#endif

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   for( int k = 0; k < graph->knots; k++ )
   {
      if( graph->grad[k] > 0 || graph->source == k )
      {
         if( !subnodesmark[k] )
            continue;
         gmlWriteNode(graph, k, file);
      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e += 2 )
   {
      if( graph->oeat[e] != EAT_FREE )
      {
         const int head = graph->head[e];
         const int tail = graph->tail[e];

         if( !subnodesmark[tail] || !subnodesmark[head] )
            continue;

         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);

         SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, NULL);
      }
   }

   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}


/** writes graph in .stp format to file */
void graph_writeStpByName(
   SCIP* scip,
   const GRAPH* g,
   const char*  filename,
   SCIP_Real    offset
   )
{
   FILE *fp;
   assert(filename);

   fp = fopen(filename, "a+");
   graph_writeStp(scip, g, fp, offset);

   fclose(fp);
}


/** writes graph in .stp format to file */
void graph_writeStp(
   SCIP* scip,
   const GRAPH* g,
   FILE*        fp,
   SCIP_Real    offset
   )
{
   const int nnodes = graph_get_nNodes(g);
   int nnodes_curr;
   int nedges_curr;
   int hopfactor;
   int* orgToNewNode;

   assert(fp != NULL);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &orgToNewNode, nnodes) );
   getOrgNodeToNodeMap(g, orgToNewNode);

   fprintf(fp, "%8x STP File, STP Format Version %2d.%02d\n",
      STP_FILE_MAGIC, STP_FILE_VERSION_MAJOR, STP_FILE_VERSION_MINOR);

   fprintf(fp, "Section Comment\n");
   fprintf(fp, "Problem ");
   switch( g->stp_type )
   {
      case STP_SPG:
         fprintf(fp, "\"%s\"\n", "SPG");
         break;

      case STP_PCSPG:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      case STP_RPCSPG:
         fprintf(fp, "\"%s\"\n", "RPCST");
         break;

      case STP_NWSPG:
         fprintf(fp, "\"%s\"\n", "NWSPG");
         break;

      case STP_DCSTP:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      case STP_NWPTSPG:
         fprintf(fp, "\"%s\"\n", "STP_NWPTSPG");
         break;

      case STP_RSMT:
         fprintf(fp, "\"%s\"\n", "RSMT");
         break;

      case STP_OARSMT:
         fprintf(fp, "\"%s\"\n", "OARSMT");
         break;

      case STP_MWCSP:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      case STP_DHCSTP:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      default:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
   }
   fprintf(fp, "End\n\n");

   if( !SCIPisEQ(scip, offset, 0.0) )
   {
      fprintf(fp, "Section Presolve\n");
      fprintf(fp, "Fixed %f \n", offset);
      fprintf(fp, "End\n\n");
   }

   if( graph_pc_isPcMw(g) )
   {
      nnodes_curr = 0;
      nedges_curr = 0;
      assert(!g->extended);

      for( int i = 0; i < nnodes; i++ )
      {
         if( !graph_pc_knotIsDummyTerm(g, i) && (g->grad[i] > 0 || GT(g->prize[i], 0.0)) )
            nnodes_curr++;
      }

      for( int i = 0; i < g->edges; i += 2 )
      {
         if( g->ieat[i] != EAT_FREE )
         {
            const int tail = g->tail[i];
            const int head = g->head[i];

            if( graph_pc_knotIsDummyTerm(g, tail) || graph_pc_knotIsDummyTerm(g, head) )
               continue;

            assert(EQ(g->cost[i], g->cost[i + 1]));

            nedges_curr += 2;
         }
      }
   }
   else
   {
      graph_get_nVET(g, &nnodes_curr, &nedges_curr, NULL);
   }

   fprintf(fp, "Section Graph\n");
   fprintf(fp, "Nodes %d\n", nnodes_curr);
   fprintf(fp, "Edges %d\n", nedges_curr / 2);

   for( int i = 0; i < g->edges; i += 2 )
   {
      if (g->ieat[i] != EAT_FREE)
      {
         const int tail = g->tail[i];
         const int head = g->head[i];
         const int tail_curr = orgToNewNode[tail];
         const int head_curr = orgToNewNode[head];

         if( graph_pc_isPcMw(g) )
         {
            if( graph_pc_knotIsDummyTerm(g, tail) || graph_pc_knotIsDummyTerm(g, head) )
               continue;

            assert(EQ(g->cost[i], g->cost[i + 1]));
         }

         assert(0 <= tail_curr && tail_curr < nnodes_curr);
         assert(0 <= head_curr && head_curr < nnodes_curr);
         assert(g->oeat[i] != EAT_FREE);

         if( g->stp_type == STP_SPG || g->stp_type == STP_DCSTP || g->stp_type == STP_RSMT || g->stp_type == STP_OARSMT || graph_pc_isPcMw(g) )
            fprintf(fp, "E ");
         else
            fprintf(fp, "AA ");

         fprintf(fp, "%d %d ", tail_curr + 1, head_curr + 1);

         if( graph_typeIsSpgLike(g) || g->stp_type == STP_DCSTP || graph_pc_isPcMw(g) )
            fprintf(fp, "%f\n", g->cost[i]);
         else
            fprintf(fp, "%f %f\n", g->cost[i], g->cost[Edge_anti(i)]);
      }
   }

   fprintf(fp, "End\n\n");
   fprintf(fp, "Section Terminals\n");
   fprintf(fp, "Terminals %d\n", g->terms);

   if( g->stp_type == STP_RPCSPG )
   {
      assert(0 <= orgToNewNode[g->source] && orgToNewNode[g->source] < nnodes_curr);
      fprintf(fp, "RootP %d\n", orgToNewNode[g->source] + 1);
   }

   for( int i = 0; i < g->knots; i++ )
   {
      const int i_curr = orgToNewNode[i];

      if( g->grad[i] == 0 )
      {
         if( graph_pc_isPcMw(g) && GT(g->prize[i], 0.0) )
         {
            assert(0 <= i_curr && i_curr < nnodes_curr);
            assert(graph_pc_knotIsNonLeafTerm(g, i));
            fprintf(fp, "TP %d %f\n", i_curr + 1, g->prize[i]);
         }
         continue;
      }

      if( graph_pc_isPcMw(g) )
      {
         if( graph_pc_knotIsDummyTerm(g, i) )
            continue;

         if( !Is_anyTerm(g->term[i]) )
            continue;

         if( !graph_pc_knotIsFixedTerm(g, i) )
         {
            assert(0 <= i_curr && i_curr < nnodes_curr);
            fprintf(fp, "TP %d %f\n", i_curr + 1, g->prize[i]);
            continue;
         }

         assert(Is_term(g->term[i]) && g->stp_type == STP_RPCSPG);
         assert(0 <= i_curr && i_curr < nnodes_curr);

         fprintf(fp, "TF %d\n", i_curr + 1);
         continue;
      }

      assert(0 <= i_curr && i_curr < nnodes_curr);

      if (Is_term(g->term[i]) )
         fprintf(fp, "T %d\n", i_curr + 1);
   }
   fprintf(fp, "End\n\n");

   /* Hop-Constrained STP */
   if( g->stp_type == STP_DHCSTP )
   {
      fprintf(fp, "Section Hop Constraint\n");
      fprintf(fp, "limit %d\n", g->hoplimit);
      for( int e = 0; e < nedges_curr / 2; e++ )
      {
         hopfactor = 1;
         fprintf(fp, "HC %d %d\n", e + 1, hopfactor);
      }
      fprintf(fp, "End\n\n");
   }

   /* Degree-Constrained STP */
   /* It is just enough to know that the degree constraint is required */
   if( g->stp_type == STP_DCSTP )
   {
      fprintf(fp, "Section Degree Constraint\n");
      fprintf(fp, "End\n\n");
   }

   fprintf(fp, "EOF\n");

   SCIPfreeBufferArray(scip, &orgToNewNode);
}


/** writes SPG (no variant!) to a file; deprecated */
void graph_writeStpOrg(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const char*           filename            /**< file name */
   )
{
   FILE *fptr;

   assert(scip != NULL);
   assert(graph != NULL);

   fptr = fopen(filename, "w");
   assert(fptr != NULL);

   fprintf(fptr, "33D32945 STP File, STP Format Version 1.0\n");
   fprintf(fptr, "SECTION Comment\n");
   fprintf(fptr, "END\n\n");
   fprintf(fptr, "SECTION Graph\n");
   fprintf(fptr, "Nodes %d\n", graph->knots);
   fprintf(fptr, "Edges %d\n", graph->edges);

   for( int e = 0; e < graph->edges; e += 2 )
      fprintf(fptr, "E %d %d %f\n", graph->tail[e] + 1, graph->head[e] + 1, graph->cost[e]);
   fprintf(fptr, "END\n\n");

   fprintf(fptr, "SECTION Terminals\n");
   fprintf(fptr, "Terminals %d\n", graph->terms);

   for( int k = 0; k < graph->knots; k++ )
      if( Is_term(graph->term[k]) )
         fprintf(fptr, "T %d\n", k + 1);

   fprintf(fptr, "END\n\n");

   fprintf(fptr, "EOF\n");

   fclose(fptr);
}
