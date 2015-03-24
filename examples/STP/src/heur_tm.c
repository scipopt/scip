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

/**@file   heur_tm.c
 * @brief  TM primal heuristic
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "heur_tm.h"
#include "probdata_stp.h"
#include "portab.h"
#include "scip/misc.h"
#include <math.h>
#define HEUR_NAME             "TM"
#define HEUR_DESC             "takahashi matsuyama primal heuristic for steiner trees"
#define HEUR_DISPCHAR         '+'
#define HEUR_PRIORITY         1000 /* TODO 0 */
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_EVALRUNS 10 /*10*/
#define DEFAULT_INITRUNS 100 /*100*/
#define DEFAULT_LEAFRUNS 10 /*10*/
#define DEFAULT_ROOTRUNS 50
#define DEFAULT_DURINGLPFREQ 10
#define DEFAULT_TYPE  0
#define DEFAULT_RANDSEED 0
#define VEL 1
#define AUTO     0
#define TM       1
#define TMPOLZIN 2
#define TMX 1

#ifdef WITH_UG
extern
int getUgRank();
#endif

/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint ncalls;
   SCIP_Longint nexecs;
   SCIP_Real hopfactor;
   int stp_type;
   int evalruns;
   int initruns;
   int leafruns;
   int rootruns;
   int duringlpfreq;
   int beststartnode;
   int* rootedges_t;
   int*                  rootedges_z;
   unsigned int          randseed;           /**< seed value for random number generator                              */
   unsigned int          timing;
};

/*
 * Static methods
 */

/** information method for a parameter change of random seed */
static
SCIP_DECL_PARAMCHGD(paramChgdRandomseed)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int newrandseed;

   newrandseed = SCIPparamGetInt(param);

   heurdata = (SCIP_HEURDATA*)SCIPparamGetData(param);
   assert(heurdata != NULL);

   heurdata->randseed = (unsigned int)newrandseed;

   return SCIP_OKAY;
}


#if 1
/** for debug purposes only */
static
SCIP_RETCODE printGraph(
   SCIP* scip,
   const GRAPH*          graph,              /**< Graph to be printed */
   const char*           filename,           /**< Name of the output file */
   int*                  result
   )
{
   char label[SCIP_MAXSTRLEN];
   FILE* file;
   int e;
   int n;
   int m;
   char* stnodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &stnodes, graph->knots ) );

   assert(graph != NULL);
   file = fopen((filename != NULL) ? filename : "graphX.gml", "w");

   for( e = 0; e < graph->knots; e++ )
   {
      stnodes[e] = FALSE;
   }
   for( e = 0; e < graph->edges; e++ )
   {
      if( result[e] == CONNECT )
      {
	 stnodes[graph->tail[e]] = TRUE;
	 stnodes[graph->head[e]] = TRUE;
      }
   }

   /* write GML format opening, undirected */
   SCIPgmlWriteOpening(file, FALSE);

   /* write all nodes, discriminate between root, terminals and the other nodes */
   e = 0;
   m = 0;
   for( n = 0; n < graph->knots; ++n )
   {
      if( stnodes[n] )
      {
         if( n == graph->source[0] )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Root", n);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "rectangle", "#666666", NULL);
            m = 1;
         }
         else if( graph->term[n] == 0 )
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Terminal %d", n, e + 1);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#ff0000", NULL);
            e += 1;
         }
         else
         {
            (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "(%d) Node %d", n, n + 1 - e - m);
            SCIPgmlWriteNode(file, (unsigned int)n, label, "circle", "#336699", NULL);
         }

      }
   }

   /* write all edges (undirected) */
   for( e = 0; e < graph->edges; e ++ )
   {
      if( result[e] == CONNECT )
      {
         (void)SCIPsnprintf(label, SCIP_MAXSTRLEN, "%8.2f", graph->cost[e]);

	 SCIPgmlWriteEdge(file, (unsigned int)graph->tail[e], (unsigned int)graph->head[e], label, "#ff0000");
      }
   }
   SCIPfreeBufferArray(scip, &stnodes);
   /* write GML format closing */
   SCIPgmlWriteClosing(file);

   return SCIP_OKAY;
}

#endif

/* prune the (rooted) prize collecting Steiner tree in such a way that all leaves are terminals */
SCIP_RETCODE do_pcprune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,               /**< edge costs */
   int*                  result,             /**< ST edges */
   int                   root,
   char*                 connected           /**< ST nodes */
   )
{
   PATH*  mst;
   int i;
   int j;
   int e1;
   int e2;
   int k1;
   int k2;
   int count;
   int nnodes;
   nnodes = g->knots;
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );

   /* compute the MST, exclude all terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( connected[i] && !Is_term(g->term[i]) )
         g->mark[i] = TRUE;
      else
         g->mark[i] = FALSE;
   }

   if( g->stp_type == STP_ROOTED_PRIZE_COLLECTING )
   {
      root = g->source[0];
      g->mark[root] = TRUE;
   }
   else
   {
      assert(!Is_term(g->term[root]));
   }
   //printf("root: %d, \n", root);
   assert(root >= 0);
   assert(root < nnodes);

   graph_path_exec(g, MST_MODE, root, cost, mst);

   for( i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && (mst[i].edge != -1) )
      {
	 assert(g->path_state[i] == CONNECT);
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);
         result[mst[i].edge] = CONNECT;
      }
   }

   /* connect all terminals */
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != g->source[0] )
      {
         e1 = g->inpbeg[i];
         assert(e1 >= 0);
         e2 = g->ieat[e1];

         if( e2 == EAT_LAST )
         {
            result[e1] = CONNECT;
         }
         else
	 {
            assert(e2 >= 0);

            assert(g->ieat[e2] == EAT_LAST);
            k1 = g->tail[e1];
            k2 = g->tail[e2];
            assert(k1 == g->source[0] || k2 == g->source[0]);
            if( k1 != g->source[0] && g->path_state[k1] == CONNECT )
            {
               result[e1] = CONNECT;
            }
            else if( k2 != g->source[0] && g->path_state[k2] == CONNECT )
            {
               result[e2] = CONNECT;
            }
            else if( k1 == g->source[0] )
            {
               result[e1] = CONNECT;
            }
            else if( k2 == g->source[0] )
            {
               result[e2] = CONNECT;
            }
	 }
      }
      else if( i == root && g->stp_type != STP_ROOTED_PRIZE_COLLECTING )
      {
	 for( e1 = g->inpbeg[i]; e1 != EAT_LAST; e1 = g->ieat[e1] )
            if( g->tail[e1] == g->source[0] )
               break;
	 assert(e1 != EAT_LAST);
	 result[e1] = CONNECT;
      }
   }
   if( 0 )
   {
      char varname[SCIP_MAXSTRLEN];
      (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "AA%d.gml", 0);
      SCIP_CALL( printGraph(scip, g, varname, result) );
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = 0; i < nnodes; i++ )
      {
         if( !g->mark[i] )
            continue;
         if( g->path_state[i] != CONNECT )
            continue;

         if( Is_term(g->term[i]) )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == CONNECT )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == 0 )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
		  //printf("disconned %d (term %d)\n", i, g->term[i] );
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }

            if( j == EAT_LAST )
	    {
               /* (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "AA1%d.gml", rruns);
                  SCIP_CALL( printGraph(scip, g, varname, result) );*/
               //printf("in %d \n", i);
	    }
            assert(j != EAT_LAST);
         }
      }
   }
   while( count > 0 );
   //assert(graph_sol_valid(g, result));
   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}


/* prune the Steiner Tree in such a way that all leaves are terminals */
SCIP_RETCODE do_prune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,               /**< edge costs */
   int                   layer,
   int*                  result,             /**< ST edges */
   char*                 connected           /**< ST nodes */
   )
{
   PATH*  mst;
   int i;
   int j;
   int count;
   int nnodes;
   nnodes = g->knots;
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nnodes) );
   assert(layer == 0);
   /*
     printf("in \n");
     if(connected[1053])
     {
     printf("connected 1053!! term? : %d\n", g->term[1053]);
     for( j = g->outbeg[1053]; j != EAT_LAST; j = g->oeat[j] )
     if( connected[g->head[j]] )
     printf("out 1053: %d->%d \n", g->tail[j], g->head[j]);

     }
   */
   j = 0;
   /* compute the MST */
   for( i = 0; i < nnodes; i++ )
   {
      if( connected[i] )
         j++;
      g->mark[i] = connected[i];
   }
   assert(g->source[layer] >= 0);
   assert(g->source[layer] <  nnodes);

   graph_path_exec(g, MST_MODE, g->source[layer], cost, mst);
   for( i = 0; i < nnodes; i++ )
   {
      if( connected[i] && (mst[i].edge != -1) )
      {
         assert(g->head[mst[i].edge] == i);
         assert(result[mst[i].edge] == -1);
         result[mst[i].edge] = layer;
      }
   }

   /* prune */
   do
   {
      SCIPdebug(fputc('C', stdout));
      SCIPdebug(fflush(stdout));

      count = 0;

      for( i = 0; i < nnodes; i++ )
      {
         if( !g->mark[i] )
            continue;

         if( g->term[i] == layer )
            continue;

         for( j = g->outbeg[i]; j != EAT_LAST; j = g->oeat[j] )
            if( result[j] == layer )
               break;

         if( j == EAT_LAST )
         {
            /* there has to be exactly one incoming edge
             */
            for( j = g->inpbeg[i]; j != EAT_LAST; j = g->ieat[j] )
            {
               if( result[j] == layer )
               {
                  result[j]    = -1;
                  g->mark[i]   = FALSE;
                  connected[i] = FALSE;
                  count++;
                  break;
               }
            }
            if( g->stp_type != STP_ROOTED_PRIZE_COLLECTING )
	    {
	       if( j == EAT_LAST )
                  printf("in %d \n", i);
               assert(j != EAT_LAST);
	    }
         }
      }
   }
   while( count > 0 );

   SCIPfreeBufferArray(scip, &mst);

   return SCIP_OKAY;
}



/* pure TM heuristic */
static
SCIP_RETCODE do_tm(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH**                path,
   SCIP_Real*            cost,
   SCIP_Real*            costrev,
   int                   layer,
   int                   start,
   int*                  result,
   char*                 connected
   )
{
   int*   cluster;
   int    csize = 0;
   int    k;
   int    e;
   int    i;
   int    j;
   int    old;
   int    newval;
   int nnodes;
   SCIP_Real min;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(path      != NULL);
   assert(layer == 0);
   nnodes = g->knots;

   SCIPdebugMessage("Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );

   cluster[csize++] = start;

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }
   connected[start] = TRUE;

   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to the current ST
       */
      min = FARAWAY;
      old = -1;
      newval = -1;

      for( i = 0; i < nnodes; i++ )
      {
         if( g->grad[i] == 0 || g->term[i] != layer || connected[i] )
            continue;

         /*
          */
         if( path[i] == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(path[i]), nnodes) );

            assert(path[i] != NULL);
            if( g->source[0] == i )
               graph_path_exec(g, FSP_MODE, i, cost, path[i]);
            else
               graph_path_exec(g, FSP_MODE, i, costrev, path[i]);
         }
         for( k = 0; k < csize; k++ )
         {
            j = cluster[k];
            assert(i != j);
            assert(connected[j]);
            if (LT(path[i][j].dist, min))
            {
               min = path[i][j].dist;
               newval = i;
               old = j;
            }
         }
      }
      /* Nichts mehr gefunden, also fertig
       */
      if (newval == -1)
         break;

      /* Weg setzten
       */
      assert((old > -1) && (newval > -1));
      assert(path[newval] != NULL);
      assert(path[newval][old].dist < FARAWAY);
      assert(g->term[newval] == layer);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /*    printf("Connecting Knot %d-%d dist=%d\n", newval, old, path[newval][old].dist);
       */
      /* Gegen den Strom schwimmend alles markieren
       */
      k = old;

      while(k != newval)
      {
         e = path[newval][k].edge;
         k = g->tail[e];

         if (!connected[k])
         {
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }

   SCIPdebug(fputc('M', stdout));
   SCIPdebug(fflush(stdout));
   SCIPfreeBufferArray(scip, &cluster);

   SCIP_CALL( do_prune(scip, g, cost, layer, result, connected) );

   return SCIP_OKAY;
}


static
SCIP_RETCODE do_tmX(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,
   SCIP_Real*            costrev,
   SCIP_Real**           pathdist,
   int                   start,
   int                   artroot,
   int*                  result,
   int**                 pathedge,
   char*                 connected
   )
{
   SCIP_Real min;
   int    csize = 0;
   int    k;
   int    e;
   int    i;
   int    j;
   int    old;
   int    newval;
   int    nnodes;
   int*   cluster;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   nnodes = g->knots;
   assert(start < nnodes);
   SCIPdebugMessage("Heuristic: Start=%5d ", start);


   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );

   cluster[csize++] = start;

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }

   connected[start] = TRUE;

   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to the current ST
       */
      min = FARAWAY;
      old = -1;
      newval = -1;
      for( i = 0; i < nnodes; i++ )
      {
         if( !Is_term(g->term[i]) || connected[i] || g->grad[i] == 0 )
            continue;

         if( pathdist[i] == NULL )
         {
	    assert(pathedge[i] == NULL);
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathdist[i]), nnodes) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathedge[i]), nnodes) );
            assert(pathedge[i] != NULL);
	    assert(pathdist[i] != NULL);
            if( g->source[0] == i )
               graph_path_execX(scip, g, i, cost,  pathdist[i], pathedge[i]);
            else
               graph_path_execX(scip, g, i, costrev, pathdist[i], pathedge[i]);
         }
         for( k = 0; k < csize; k++ )
         {
            j = cluster[k];

            assert(i != j);
            assert(connected[j]);

            if( SCIPisLT(scip, pathdist[i][j], min) )
            {
               min = pathdist[i][j];
               newval = i;
               old = j;
            }
         }
      }
      /* Nichts mehr gefunden, also fertig
       */
      if( newval == -1 )
         break;

      /* Weg setzten
       */
      assert(old > -1);
      assert(newval > -1);
      assert(pathdist[newval] != NULL);
      assert(pathdist[newval][old] < FARAWAY);
      assert(g->term[newval] == 0);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /*    printf("Connecting Knot %d-%d dist=%d\n", newval, old, path[newval][old].dist);
       */
      /* Gegen den Strom schwimmend alles markieren
       */
      k = old;
      //printf("connect from %d\n", old);
      while(k != newval)
      {
         e = pathedge[newval][k];
         k = g->tail[e];
         if (!connected[k])
         {
            //printf("connected (%d->)%d  cost: %f (rev) %f\n", g->head[e], k, cost[k], costrev[k] );
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
   }

   SCIPdebug(fputc('M', stdout));
   SCIPdebug(fflush(stdout));
   SCIPfreeBufferArray(scip, &cluster);

   for( e = 0; e < g->edges; e++ )
   {
      assert(cost[e] >= 0);
      if( SCIPisLT(scip, cost[e], 1e+9) )
         cost[e] = g->cost[e];
   }
   if( g->stp_type == STP_UNDIRECTED || g->stp_type == STP_GRID || g->stp_type == STP_OBSTACLES_GRID )
      SCIP_CALL( do_prune(scip, g, cost, 0, result, connected) );
   else
      SCIP_CALL( do_pcprune(scip, g, cost, result, artroot, connected) );
   for( e = 0; e < g->edges; e++ )
      cost[e] = costrev[flipedge(e)];

   return SCIP_OKAY;
}


/* pure TM heuristic for degree constrained STPs */
static
SCIP_RETCODE do_tm_degcons(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_Real*            cost,
   SCIP_Real*            costrev,
   SCIP_Real**           pathdist,
   int                   start,
   int*                  result,
   int**                 pathedge,
   char*                 connected,
   char*                 solfound
   )
{
   SCIP_QUEUE* queue;
   SCIP_Real min;
   int    csize = 0;
   int    k;
   int    e;
   int    i;
   int    j;
   int    t;
   int    u;
   int    old;
   int    degcount;
   int    degmax;
   int    newval;
   int    nnodes;
   int    nterms;
   int*   cluster;
   int*   realterms;
   int*   degs;
   int*   maxdegs;
   assert(scip      != NULL);
   assert(g         != NULL);
   assert(g->maxdeg != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   assert(pathdist   != NULL);
   assert(pathedge   != NULL);
   nnodes = g->knots;

   realterms = SCIPprobdataGetRTerms(scip);
   nterms = SCIPprobdataGetRNTerms(scip) + 1;

   SCIPdebugMessage("Heuristic: Start=%5d ", start);

   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &degs, nnodes) );
   maxdegs = g->maxdeg;
   cluster[csize++] = start;

   for( i = 0; i < nnodes; i++ )
   {
      degs[i] = 0;
      g->mark[i]   = (g->grad[i] > 0);
      connected[i] = FALSE;
   }

   for( e = 0; e < g->edges; e++ )
      result[e] = -1;
   connected[start] = TRUE;

   //printf("root %d ,start: %d \n", g->source[0], start);
   /* CONSTCOND */
   for( ;; )
   {
      /* Find a terminal with minimal distance to the current ST
       */
      min = FARAWAY;
      degmax = -1;
      old = -1;
      newval = -1;
      for( t = 0; t < nterms; t++ )
      {
	 if( t != 0 )
            i = realterms[t - 1];
         else
	    i = g->source[0];
         if( connected[i] || g->grad[i] == 0 )
            continue;

         if( pathdist[i] == NULL )
         {
	    assert(pathedge[i] == NULL);
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathdist[i]), nnodes) );
            SCIP_CALL( SCIPallocBufferArray(scip, &(pathedge[i]), nnodes) );
            assert(pathedge[i] != NULL);
	    assert(pathdist[i] != NULL);
            if( g->source[0] == i )
               graph_path_execX(scip, g, i, cost,  pathdist[i], pathedge[i]);
            else
               graph_path_execX(scip, g, i, costrev, pathdist[i], pathedge[i]);
         }
         for( k = 0; k < csize; k++ )
         {
            j = cluster[k];

            assert(i != j);
            assert(connected[j]);

            if( SCIPisLT(scip, pathdist[i][j], min) && degs[j] < maxdegs[j])
            {
	       u = j;
	       degcount = 0;
	       while( u != i )
	       {
	          u = g->tail[pathedge[i][u]];
		  if( !connected[u] )
		     degcount += MIN(g->grad[u] - 1, maxdegs[u] - 1);
	       }
	       if( degcount >= degmax || degcount > 2 )
	       {
		  degmax = degcount;
                  min = pathdist[i][j];
                  newval = i;
                  old = j;
	       }
            }
         }
      }
      /* Nichts mehr gefunden, also fertig
       */
      if( newval == -1 )
         break;
      //printf("LONGCONNECT: %d %d \n", old, newval);
      /* Weg setzten
       */
      assert(old > -1);
      assert(newval > -1);
      assert(pathdist[newval] != NULL);
      assert(pathdist[newval][old] < FARAWAY);
      assert(g->term[newval] == 0);
      assert(!connected[newval]);
      assert(connected[old]);

      SCIPdebug(fputc('R', stdout));
      SCIPdebug(fflush(stdout));

      /*    printf("Connecting Knot %d-%d dist=%d\n", newval, old, path[newval][old].dist);
       */
      /* Gegen den Strom schwimmend alles markieren
       */
      k = old;

      while(k != newval)
      {
         e = pathedge[newval][k];
	 u = k;
         k = g->tail[e];

         //    printf("connect: %d->%d \n", g->tail[e], g->head[e]);
         if( !connected[k])
         {
	    result[flipedge(e)] = CONNECT;
            result[e] = CONNECT;
	    degs[u]++;
            connected[k] = TRUE;
            cluster[csize++] = k;
         }
      }
      assert(degs[newval] == 0);
      degs[newval] = 1;
   }

   SCIPdebug(fputc('M', stdout));
   SCIPdebug(fflush(stdout));
   SCIPfreeBufferArray(scip, &cluster);

   *solfound = TRUE;

   for( t = 0; t < nterms; t++ )
   {
      if( t != 0 )
         i = realterms[t - 1];
      else
         i = g->source[0];
      if( !connected[i] )
      {
         //printf("fail! \n");
         *solfound = FALSE;
         break;

      }
   }

   if( *solfound )
   {
      int* pnode;
      int termcount;
      /* BFS until all terminals are reached */
      SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2) );

      SCIP_CALL( SCIPqueueInsert(queue, &(g->source[0])) );
      termcount = 1;

      for( i = 0; i < nnodes; i++ )
         connected[i] = FALSE;

      connected[g->source[0]] = TRUE;
      while( !SCIPqueueIsEmpty(queue) )
      {
         pnode = (SCIPqueueRemove(queue));
         for( e = g->outbeg[*pnode]; e != EAT_LAST; e = g->oeat[e] )
         {
            if( result[e] == CONNECT && !(connected[g->head[e]]) )
            {
               //printf("%d->%d \n", g->tail[e], g->head[e]);
               i = g->head[e];
               result[flipedge(e)] = -1;
               connected[i] = TRUE;
               if( Is_term(g->term[i]) )
               {
                  termcount++;
               }
               SCIP_CALL( SCIPqueueInsert(queue, &(g->head[e])) );
            }
         }
      }

      SCIPqueueFree(&queue);
      // assert(graph_sol_valid(g, result));
      //assert(0);
      for( t = 0; t < nnodes; t++ )
         if( degs[t] > maxdegs[t] )
	 {
	    //printf("deg fail: %d (%d, %d)\n ", t, degs[t], maxdegs[t] );
            *solfound = FALSE;
	 }
   }

   //SCIP_CALL( do_prune(scip, g, cost, 0, result, connected) );
   SCIPfreeBufferArray(scip, &degs);
   return SCIP_OKAY;
}

static
SCIP_RETCODE do_tm_polzin(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   SCIP_PQUEUE*          pqueue,
   GNODE**               gnodearr,
   SCIP_Real*            cost,
   SCIP_Real*            costrev,
   int                   layer,
   SCIP_Real**           node_dist,
   int                   start,
   int*                  result,
   int*                  vcount,
   int*                  nodenterms,
   int**                 node_base,
   int**                 node_edge,
   char                  firstrun,
   char*                 connected
   )
{
   int    k;
   int    i;
   int    j;
   int    best;
   int    term;
   int    count;
   int   nnodes;
   int   nterms;

   assert(scip      != NULL);
   assert(g         != NULL);
   assert(result    != NULL);
   assert(connected != NULL);
   assert(cost      != NULL);
   assert(costrev   != NULL);
   nnodes = g->knots;
   nterms = g->terms;

   SCIPdebugMessage("TM_Polzin Heuristic: Start=%5d ", start);

   /* if the heuristic is called for the first time several data structures have to be set up */
   if( firstrun )
   {
      PATH* vnoi;
      SCIP_Real* vcost;
      int old;
      int oedge;
      int root = g->source[0];
      int   ntovisit;
      int   nneighbnodes;
      int   nneighbterms;
      int   nreachednodes;
      int*  state;
      int*  vbase;
      int*  terms;
      int*  tovisit;
      int*  reachednodes;
      char* termsmark;
      char* visited;
      int e;
      /* PHASE I: */
      for( i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);

      /* allocate memory needed in PHASE I */
      SCIP_CALL( SCIPallocBufferArray(scip, &terms, nterms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &termsmark, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &visited, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &reachednodes, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tovisit, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vcost, nnodes) );

      j = 0;
      for( i = 0; i < nnodes; i++ )
      {
         visited[i] = FALSE;
         if( Is_term(g->term[i]) )
         {
            termsmark[i] = TRUE;
            terms[j++] = i;
         }
         else
         {
            termsmark[i] = FALSE;
         }
      }

      for( e = 0; e < g->edges; e++)
      {
         assert(SCIPisGE(scip, cost[e], 0));
         assert(SCIPisGE(scip, costrev[e], 0));
      }

      assert(j == nterms);
      voronoi(scip, g, cost, costrev, termsmark, vbase, vnoi);
      state = g->path_state;

      for( i = 0; i < nnodes; i++ )
         if( Is_term(g->term[i]) )
            assert(vbase[i] == i);

      for( k = 0; k < nnodes; k++ )
	 assert(termsmark[vbase[k]]);

      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
	 vcount[k] = 0;
	 SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[k]) );
	 gnodearr[k]->number = k;
         if( !Is_term(g->term[k]) )
         {
	    node_dist[k][0] = vnoi[k].dist;
            node_edge[k][0] = vnoi[k].edge;
            node_base[k][0] = vbase[k];
            nodenterms[k] = 1;
         }
         else
         {
            nodenterms[k] = 0;
	    node_edge[k][0] = UNKNOWN;
            termsmark[k] = FALSE;
         }
         state[k] = UNKNOWN;
         vcost[k] = vnoi[k].dist;
         vnoi[k].dist = FARAWAY;
      }

      /* for each terminal: extend the voronoi regions until all neighbouring terminals have been visited */
      for( i = 0; i < nterms; i++ )
      {
	 term = terms[i];
         nneighbterms = 0;
         nneighbnodes = 0;
         nreachednodes = 0;
         for( k = 0; k < nnodes; k++ )
            assert(termsmark[k] == FALSE);
         /* DFS (starting from terminal i) until the entire voronoi region has been visited */
         tovisit[0] = term;
         ntovisit = 1;
         visited[term] = TRUE;
         state[term] = CONNECT;
         while( ntovisit > 0 )
         {
            /* iterate all incident edges */
            old = tovisit[--ntovisit];

	    for( oedge = g->outbeg[old]; oedge != EAT_LAST; oedge = g->oeat[oedge] )
            {
               k = g->head[oedge];

               /* is node k in the voronoi region of the i-th terminal ? */
               if( vbase[k] == term )
	       {
                  if( !visited[k] )
                  {
                     state[k] = CONNECT;
		     assert(nnodes - (nneighbnodes + 1) > ntovisit);
                     tovisit[ntovisit++] = k;
                     visited[k] = TRUE;
                     reachednodes[nreachednodes++] = k;
                  }
	       }
               else
	       {
                  if( !visited[k] )
                  {
                     visited[k] = TRUE;
                     vnoi[k].dist = vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge]);
                     vnoi[k].edge = oedge;

                     if( termsmark[vbase[k]] == FALSE )
                     {
                        termsmark[vbase[k]] = TRUE;
                        nneighbterms++;
                     }
                     assert(nnodes - (nneighbnodes + 1) > ntovisit - 1);
                     tovisit[nnodes - (++nneighbnodes)] = k;
                  }
                  else
                  {
                     /* if edge 'oedge' allows a shorter connection of node k, update */
                     if( SCIPisGT(scip, vnoi[k].dist, vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge])) )
                     {
                        vnoi[k].dist = vcost[old] + ((vbase[k] == root)? cost[oedge] : costrev[oedge]);
                        vnoi[k].edge = oedge;
                     }
                  }
	       }
            }
         }

         count = 0;
         for( j = 0; j < nneighbnodes; j++ )
         {
	    assert(termsmark[vbase[tovisit[nnodes - j - 1]]]);
            heap_add(g->path_heap, state, &count, tovisit[nnodes - j - 1], vnoi);
         }
         SCIP_CALL( voronoi_extend2(scip, g, ((term == root)? cost : costrev), vnoi, node_dist, node_base, node_edge, termsmark, reachednodes, &nreachednodes, nodenterms,
               nneighbterms, term, nneighbnodes) );

         reachednodes[nreachednodes++] = term;

         for( j = 0; j < nreachednodes; j++ )
         {
            vnoi[reachednodes[j]].dist = FARAWAY;
            state[reachednodes[j]] = UNKNOWN;
            visited[reachednodes[j]] = FALSE;
         }

         for( j = 0; j < nneighbnodes; j++ ) /* TODO AVOID DOUBLE WORK */
         {
            vnoi[tovisit[nnodes - j - 1]].dist = FARAWAY;
            state[tovisit[nnodes - j - 1]] = UNKNOWN;
            visited[tovisit[nnodes - j - 1]] = FALSE;
         }
      }

      /* for each node v: sort the terminal arrays according to their distance to v */
      for( i = 0; i < nnodes; i++ )
	 SCIPsortRealIntInt(node_dist[i], node_base[i], node_edge[i], nodenterms[i]);

      /* free memory */
      SCIPfreeBufferArray(scip, &vnoi);
      SCIPfreeBufferArray(scip, &vbase);
      SCIPfreeBufferArray(scip, &terms);
      SCIPfreeBufferArray(scip, &termsmark);
      SCIPfreeBufferArray(scip, &visited);
      SCIPfreeBufferArray(scip, &tovisit);
      SCIPfreeBufferArray(scip, &reachednodes);
      SCIPfreeBufferArray(scip, &vcost);
   }

   /** PHASE II **/
   else
   {
      for( k = 0; k < nnodes; k++ )
      {
         connected[k] = FALSE;
	 vcount[k] = 0;
      }
   }

   connected[start] = TRUE;
   gnodearr[start]->dist = node_dist[start][0];
   SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[start]) );

   while( SCIPpqueueNElems(pqueue) > 0 )
   {
      best = ((GNODE*) SCIPpqueueRemove(pqueue))->number;

      term = node_base[best][vcount[best]];
      assert( Is_term(g->term[term]) );
      /* has the terminal already been connected? */
      if( !connected[term] )
      {
	 /*printf("term: %d \n", term ); */
         /* connect the terminal */
         k = g->tail[node_edge[best][vcount[best]]];
         while( k != term )
         {
            j = 0;

	    while( node_base[k][vcount[k] + j] != term )
               j++;

            if(!(vcount[k] + j < nodenterms[k]))
	    {
	       int z;
               printf("  vcount[k] + j: %d, nodenterms: %d\n", vcount[k] + j,nodenterms[k]);
               for( z = 0; z < nodenterms[k]; z++)
		  if( node_base[k][z] == term )
		     printf("ok! for z: %d\n", z);
		  else
		     printf("  node_base[k][z]: %d ", node_base[k][z]);
               printf("\n term: %d \n", term );
               assert(0);
	    }

            if( !connected[k] )
            {
	       assert(vcount[k] == 0);

               connected[k] = TRUE;
	       while( vcount[k] < nodenterms[k] && connected[node_base[k][vcount[k]]] )
               {
	          vcount[k]++;
		  j--;
               }

               if( vcount[k] < nodenterms[k] )
	       {
		  gnodearr[k]->dist = node_dist[k][vcount[k]];
		  SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
	       }
	    }

            assert( vcount[k] + j < nodenterms[k] );
	    assert(node_base[k][vcount[k] + j] == term);
	    k = g->tail[node_edge[k][vcount[k] + j]];
         }
         /* finally, connected the terminal reached*/

         assert( k == term );
	 assert( !connected[k] );
         connected[k] = TRUE;

         assert( vcount[k] == 0 );
         while( vcount[k] < nodenterms[k] && connected[node_base[k][vcount[k]]] )
         {
            vcount[k]++;
         }
         if( vcount[k] < nodenterms[k] )
         {
            gnodearr[k]->dist = node_dist[k][vcount[k]];
            SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[k]) );
            /*SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[k][vcount[k]]) );*/
         }
      }

      while( vcount[best] + 1 < nodenterms[best] )
      {
         if( !connected[node_base[best][++vcount[best]]] )
         {
	    gnodearr[best]->dist = node_dist[best][vcount[best]];
	    SCIP_CALL( SCIPpqueueInsert(pqueue, gnodearr[best]) );
            /*SCIP_CALL( SCIPpqueueInsert(pqueue, nodeterms[best][vcount[best]]) );*/
            break;
         }
      }
   }

   /* prune the ST, so that all leaves are terminals */
   SCIP_CALL( do_prune(scip, g, cost, layer, result, connected) );
   /*
     for( k = 0; k < g->edges; k++ )
     if( SCIPisLT(scip, cost[k], 1e+9) )
     cost[k] = g->cost[k];
     SCIP_CALL( do_prune(scip, g, cost, 0, result, connected) );
     for( k = 0; k < g->edges; k++ )
     cost[k] = costrev[flipedge(k)];
   */

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPtmHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,                  /**< graph structure */
   PATH**                path,
   SCIP_Real*            cost,
   SCIP_Real*            costrev,
   //SCIP_Real**           pathdist,
   int*                  result
   //int**                 pathedge
   )
{
   int e;
   int nnodes;
   char* connected;
   assert(scip != NULL);
   assert(graph != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   nnodes = graph->knots;
   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );

   for( e = 0; e < graph->edges; e++ )
      result[e] = UNKNOWN;
   SCIP_CALL( do_tm(scip, graph, path, cost, costrev, 0, graph->source[0], result, connected) );
   SCIPfreeBufferArray(scip, &connected);
   return SCIP_OKAY;
}

static
void do_prizecoll_trivial(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*  graph,
   int*          best_result
   )
{
   double maxcost = -1;
   int e;
   int i = -1;
   int maxedge = UNKNOWN;
   int maxterm;
   int root;
   assert(graph != NULL);
   assert(best_result != NULL);
   root = graph->source[0];

   for( e = graph->outbeg[root]; e != EAT_LAST; e = graph->oeat[e] )
   {
      if( Is_term(graph->term[graph->head[e]]) )
      {
	 best_result[e] = CONNECT;
         if( SCIPisLT(scip, maxcost, graph->cost[e]) )
	 {
	    maxcost = graph->cost[e];
	    maxedge = e;
	 }
      }
   }

   assert(maxedge >= 0);
   maxterm = graph->head[maxedge];
   for( e = graph->inpbeg[maxterm]; e != EAT_LAST; e = graph->ieat[e] )
   {
      if( graph->tail[e] != root )
      {
	 best_result[e] = CONNECT;
	 i = graph->tail[e];
	 break;
      }
   }
   assert(i >= 0);
   for( e = graph->inpbeg[i]; e != EAT_LAST; e = graph->ieat[e] )
      if( graph->tail[e] == root )
	 break;
   assert(e != EAT_LAST);
   best_result[maxedge] = UNKNOWN;
   best_result[e] = CONNECT;
}



SCIP_RETCODE do_layer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,
   const GRAPH*          graph,
   int*          starts,
   int*          bestnewstart,
   int*          best_result,
   int           runs,
   int           bestincstart,
   SCIP_Real*    cost,
   SCIP_Real*    costrev,
   SCIP_Real      maxcost
   )
{
#if !TMX
   PATH** path;
#endif
   SCIP_PQUEUE* pqueue;
   SCIP_Longint nexecs;
   GNODE** gnodearr;
   SCIP_Real obj;
   SCIP_Real objt;
   SCIP_Real min = FARAWAY;
   SCIP_Real** pathdist;
   SCIP_Real** node_dist;
   int best;
   int k;
   int r;
   int e;
   int mode;
   int nedges;
   int nnodes;
   int nterms;
   int* start;
   int* vcount;
   int* result;
   int* nodenterms;
   int** node_base;

   int** node_edge;
   int** pathedge;

   char printfs = FALSE;
   char* connected;

   best = bestincstart;
   for( e = 0; e < graph->edges; e++)
   {
      assert(SCIPisGE(scip, cost[e], 0));
      assert(SCIPisGE(scip, costrev[e], 0));
   }
   assert(scip != NULL);
   assert(graph != NULL);
   assert(best_result != NULL);
   assert(cost != NULL);
   assert(costrev != NULL);
   assert(heurdata != NULL);
   assert(runs > 0);

   nexecs = heurdata->nexecs;
   nnodes = graph->knots;
   nedges = graph->edges;
   nterms = graph->terms;

   assert(nexecs >= 0);
   assert(nedges > 0);
   assert(nnodes > 0);
   assert(nterms > 0);
   SCIP_CALL( SCIPallocBufferArray(scip, &connected, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &start, MIN(runs, nnodes)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );

   /* get user parameter */
   SCIP_CALL( SCIPgetIntParam(scip, "heuristics/"HEUR_NAME"/type", &mode) );
   if( printfs )
      printf(" tmmode: %d -> ", mode);
   assert(mode == AUTO || mode == TM || mode == TMPOLZIN);

   if( mode == AUTO )
   {
      /* are there enough terminals for the TM Polzin variant to (expectably) be advantageous? */
      if( SCIPisGE(scip, ((double) nterms) / ((double) nnodes ), 0.1) )
         mode = TMPOLZIN;
      else
         mode = TM;
   }

   if( graph->stp_type == STP_ROOTED_PRIZE_COLLECTING )
      mode = TM;
   if( printfs )
      printf(" %d \n ", mode);

   if( graph->layers > 1 )
   {
      /*  SCIP_CALL( do_heuristic(scip, graph, layer, best_result, graph->source[layer], connected, cost, costrev, path) );*/
      assert(0);
   }
   else
   {
#if 0
      if( graph->stp_type == STP_ROOTED_PRIZE_COLLECTING && runs != heurdata->initruns )
      {
	 runs = graph->terms;
	 printf("NO INITRUNS\n");
      }
#endif
      if( graph->stp_type == STP_DEG_CONS )
         mode = TM;
      if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
      {
	 if( runs > (nterms - 1) )
	    runs = nterms - 1;
	 mode = TM;
      }
      else if( starts != NULL )
      {
	 for( k = 0; k < MIN(runs, nnodes); k++ )
	    start[k] = starts[k];
      }
      else if( runs < nnodes )
      {
	 int* perm = NULL; /* permutation array */
	 if( SCIPgetRandomInt(0, 2, &(heurdata->randseed)) == 2 || best == -1 )
	    best = graph->source[0];

         r = 0;

         /* allocate memory for permutation array */
         SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );
         for( k = 0; k < nnodes; k++ )
            perm[k] = k;
         SCIPpermuteIntArray(perm, 0, nnodes, &(heurdata->randseed));

         /* use terminals (randomly permutated) as starting points for TM heuristic */
         for( k = 0; k < nnodes; k++ )
         {
	    if( r >= runs || r >= nterms )
	       break;

            if( Is_term(graph->term[perm[k]]) )
            {
               //printf(" perm1[%d): %d \n", r, perm[k]);
               start[r++] = perm[k];
            }
	 }

	 /* still empty slots in start array? */
	 if( r < runs )
	 {
            /* fill empty slots randomly */
            for( k = 0; k < nnodes; k++ )
            {
               if( r >= runs )
                  break;
	       if( !Is_term(graph->term[perm[k]]) )
	       {
		  //printf(" perm2[%d): %d \n", r, perm[k]);
                  start[r++] = perm[k];
	       }
            }
	 }

         /* check whether we have a already selected the best starting node */
         for( r = 0; r < runs; r++ )
	 {
            if( start[r] == best )
               break;
	 }
         /* no, we still have to */
         if( r == runs )
            start[nexecs % runs] = best;

	 /* free memory */
         SCIPfreeBufferArray(scip, &perm);
      }
      else
      {
         runs = nnodes;
	 for( k = 0; k < nnodes; k++ )
            start[k] = k;
      }

      if( mode == TM )
      {
#if TMX
	 SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
	 SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
	 BMSclearMemoryArray(pathdist, nnodes);
	 BMSclearMemoryArray(pathedge, nnodes);
#else
         SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
         BMSclearMemoryArray(path, nnodes);
#endif
         /* for( k = 0; k < nnodes; k++ )
            path[k] = NULL;  TODO why???*/
      }
      else
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &nodenterms, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &node_base, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &node_dist, nnodes) );
         SCIP_CALL( SCIPallocBufferArray(scip, &node_edge, nnodes) );

         for( k = 0; k < nnodes; k++ )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &node_base[k], nterms) );
            SCIP_CALL( SCIPallocBufferArray(scip, &node_dist[k], nterms) );
            SCIP_CALL( SCIPallocBufferArray(scip, &node_edge[k], nterms) );
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &vcount, nnodes) );
         SCIP_CALL( SCIPpqueueCreate( &pqueue, nnodes, 2, GNODECmpByDist) );
      }

      /* incorrect if layers > 1 ! */
      assert(graph->layers == 1);
      if( graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT )
      {
	 SCIP_Bool firstrun;
         int root;
         int rootedge;
         int t;
         int z;
         int nzterms;
         int* perm;
         int* edges_tz;
         int* rootedges_t;
         int* rootedges_z;
         k = 0;
         r = 0;
         z = 0;
         nzterms = 0;
         root = graph->source[0];

         SCIP_CALL( SCIPallocBufferArray(scip, &(rootedges_t), graph->terms - 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(rootedges_z), graph->terms - 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(edges_tz), graph->terms - 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &perm, graph->terms - 1) );
         for( k = 0; k < graph->terms - 1; k++ )
         {
            rootedges_t[k] = UNKNOWN;
            rootedges_z[k] = UNKNOWN;
            edges_tz[k] = UNKNOWN;
            perm[k] = k;
         }

         for( k = 0; k < graph->knots; k++ )
         {
            if( Is_term(graph->term[k]) && k != root )
            {
               t = UNKNOWN;
               assert(graph->grad[k] > 0);
               for( e = graph->outbeg[k]; e != EAT_LAST; e = graph->oeat[e] )
               {
                  if( graph->head[e] == root )
                  {
                     rootedges_t[z] = e;
		  }
                  else
		  {

		     if( SCIPisEQ(scip, costrev[e], 0.0) )
		     {
                        assert(costrev[e] == 0);
                        assert(cost[flipedge(e)] == 0);
                        edges_tz[z] = e;
		     }
                     else
		     {
		        edges_tz[z] = UNKNOWN;
			//printf("costrevnorm %d: %f\n", e, costrev[e] );
		     }
                     t = graph->head[e];

                  }
               }
               if( t != UNKNOWN )
               {
                  for( e = graph->outbeg[t]; e != EAT_LAST; e = graph->oeat[e] )
                  {
                     if( graph->head[e] == root )
                     {
		        assert(SCIPisEQ(scip, graph->cost[flipedge(e)], 0.0));
                        cost[flipedge(e)] = FARAWAY;
                        costrev[e] = FARAWAY;
                        rootedges_z[z] = e;
                        nzterms++;
                     }

                  }
               }
               assert(rootedges_t[z] != UNKNOWN);
               z++;
            }
         }
         assert(z == nterms - 1);

         if( runs > nzterms )
            runs = nzterms;
         if( runs < nzterms )
            SCIPpermuteIntArray(perm, 0, nterms - 1, &(heurdata->randseed));
         z = 0;
	 firstrun = TRUE;
         for( r = 0; r < runs; r++ )
         {
            while( rootedges_z[perm[z]] == UNKNOWN && z < nterms - 2 )
               z++;
            assert(z <= nterms - 2);

	    /* if the edge has been fixed, continue*/
	    if( edges_tz[perm[z]] == UNKNOWN )
	       continue;

            rootedge = rootedges_z[perm[z]];

            costrev[rootedge] = 0.0;
            cost[flipedge(rootedge)] = 0.0;

            for( e = 0; e < nedges; e++ )
               result[e] = UNKNOWN;

            if( !firstrun )
            {
               for( k = 0; k < nterms - 1; k++ )
               {
                  e = rootedges_t[k];
                  assert( e != UNKNOWN );
#if TMX
                  pathdist[graph->tail[e]][root] = cost[flipedge(e)];
                  pathedge[graph->tail[e]][root] = e;
#else
                  path[graph->tail[e]][root].dist = cost[flipedge(e)];
                  path[graph->tail[e]][root].edge = e;
#endif
               }
               assert(edges_tz[perm[z]] != UNKNOWN);
               k = graph->tail[edges_tz[perm[z]]]; //k = graph->tail[rootedge + 2];
#if TMX
               pathdist[k][root] = 0.0;
               pathedge[k][root] = rootedge;
#else
               path[k][root].dist = 0.0;
               path[k][root].edge = rootedge;
#endif
               //  assert(graph->head[rootedge + 2] == root);
            }
#if TMX
            SCIP_CALL( do_tmX(scip, graph, cost, costrev, pathdist, root, graph->tail[rootedge], result, pathedge, connected) );
#else
            SCIP_CALL( do_tm(scip, graph, path, cost, costrev, 0, root, result, connected) );
#endif


            firstrun = FALSE;
            obj = 0.0;

            /* here another measure than in the do_(...) heuristics is being used*/
            for( e = 0; e < nedges; e++)
            {
               obj += (result[e] > -1) ? graph->cost[e] : 0.0;
               //  printf("obj: %f \n ", obj);
            }

            if( SCIPisLT(scip, obj, min) )
            {
               min = obj;
               if( printfs )
                  printf(" Obj(run: %d, ncall: %d)=%.12e\n", r, (int) nexecs, obj);
               for( e = 0; e < nedges; e++ )
               {
                  best_result[e] = result[e];
               }
            }
            cost[flipedge(rootedge)] = FARAWAY;
            costrev[rootedge] = FARAWAY;
            z++;
         }

         for( r = 0; r < nterms - 1; r++ )
         {

            rootedge = rootedges_z[r];
            if( rootedge != UNKNOWN )
            {
               assert(costrev[rootedge] == FARAWAY);
               costrev[rootedge] = 0.0;
               cost[flipedge(rootedge)] = 0.0;
            }
         }
         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &edges_tz);
         SCIPfreeBufferArray(scip, &rootedges_z);
         SCIPfreeBufferArray(scip, &rootedges_t);
      }
      else
      {
         SCIP_Real hopfactor = heurdata->hopfactor;
         double* orgcost = NULL;
         int edgecount;
         char solfound = FALSE;
         assert(SCIPisGT(scip, hopfactor, 0.0));
         if( SCIPisLE(scip, maxcost, 0.0) )
            maxcost = 1.0;

         if( graph->stp_type == STP_HOP_CONS )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &orgcost, nedges) );
            BMScopyMemoryArray(orgcost, cost, nedges);
         }
         for( r = 0; r < runs; r++ )
         {
            for( e = 0; e < nedges; e++ )
               result[e] = UNKNOWN;

            if( graph->stp_type == STP_HOP_CONS && (hopfactor != heurdata->hopfactor || r == 0)  )
            {
               heurdata->hopfactor = hopfactor;
               for( e = 0; e < nedges; e++ )
               {
                  if( Is_term(graph->term[graph->tail[e]]) && graph->tail[e] != graph->source[0] && (SCIPisLT(scip, graph->cost[e], 1e+9 )) )
                     assert(0);
                  if( (SCIPisLT(scip, cost[e], 1e+8 )) )
                     cost[e] = 1 + orgcost[e] / (hopfactor * maxcost);
               }
               for( e = 0; e < nedges; e++)
                  costrev[e] = cost[flipedge(e)];
            }
            assert(start[r] >= 0);
            assert(start[r] < nnodes);
            if( graph->stp_type == STP_DEG_CONS )
            {
#if TMX
               SCIP_CALL( do_tm_degcons(scip, graph, cost, costrev, pathdist, start[r], result, pathedge, connected, &solfound) );
               /*if( solfound )
                 printf("=) yeah");
                 else
                 printf("oh noo \n");*/
#endif
            }
            else if( mode == TM )
#if TMX
               SCIP_CALL( do_tmX(scip, graph, cost, costrev, pathdist, start[r], -1, result, pathedge, connected) );
#else
            SCIP_CALL( do_tm(scip, graph, path, cost, costrev, 0, start[r], result, connected) );
#endif
            else
               SCIP_CALL( do_tm_polzin(scip, graph, pqueue, gnodearr, cost, costrev, 0, node_dist, start[r], result, vcount,
                     nodenterms, node_base, node_edge, (r == 0), connected) );
            obj = 0.0;
            objt = 0.0;
            edgecount = 0;
            /* here another measure than in the do_(...) heuristics is being used*/
            for( e = 0; e < nedges; e++)
            {
               if(result[e] > -1)
               {
                  obj += graph->cost[e];
                  edgecount++;
               }
            }
            SCIPdebugMessage(" Obj=%.12e\n", obj);

            if( SCIPisLT(scip, obj, min) && (graph->stp_type != STP_DEG_CONS || solfound) )
            {
               if( graph->stp_type != STP_HOP_CONS || edgecount <= graph->hoplimit )
               {
                  min = obj;
                  *bestnewstart = start[r];
                  SCIPdebugMessage(" Objt=%.12e    ", objt);
                  if( 0 || printfs )
                     printf(" Obj(run: %d, ncall: %d)=%.12e, count: %d  limit: %d \n", r, (int) nexecs, obj, edgecount, graph->hoplimit);
                  if( printfs )
                     printf(" Objt: %.12e\n", objt);

                  for( e = 0; e < nedges; e++ )
                     best_result[e] = result[e];

               }

               /*    else
                     {
                     assert(edgecount > graph->hoplimit);
                     //printf(" before hopfactor: %f int: %d ", hopfactor, edgecount - graph->hoplimit);
                     hopfactor = hopfactor * (1.0 + MIN(1.0, (((double)(edgecount - graph->hoplimit))) / graph->hoplimit));
                     if( SCIPisLE(scip, hopfactor, 0.0) )
                     hopfactor = 1.0;
                     //printf(" aft hopfactor: %f int: %d ", hopfactor, edgecount - graph->hoplimit);
                     }*/

            }
            if( graph->stp_type == STP_HOP_CONS && (edgecount > graph->hoplimit || SCIPisGT(scip, fabs(edgecount - graph->hoplimit)/(double) graph->hoplimit, 0.2)) )
            {
               // printf("old hopfactor: %f fac: %f hoplimit: %d ", hopfactor, fabs((double) (edgecount - graph->hoplimit)), graph->hoplimit);
               hopfactor = hopfactor * (1.0 + fabs(edgecount - graph->hoplimit) / (double) graph->hoplimit);
               //|| SCIPisGT(scip, fabs(edgecount - graph->hoplimit)/(double) graph->hoplimit, 0.05))
               //    printf("new hopfactor: %f  ", hopfactor);
               if( SCIPisLE(scip, hopfactor, 0.0) )
                  hopfactor = 1.0;
               //  printf(" aft hopfactor: %f count: %d  limit: %d  \n", hopfactor, edgecount, graph->hoplimit);
               //	  printf("new1 hopfactor: %f  ", hopfactor);
            }
         }
         if( graph->stp_type == STP_HOP_CONS )
            SCIPfreeBufferArray(scip, &orgcost);
         //	 printf(" hopfactor: %f obj: %.12e\n", hopfactor, obj);
      }
   }


   /* free allocated memory */
   if( mode == TM )
   {
      for( k = nnodes - 1; k >= 0; k-- )
      {
#if TMX
         SCIPfreeBufferArrayNull(scip, &(pathdist[k]));
         SCIPfreeBufferArrayNull(scip, &(pathedge[k]));
#else
         assert(path[k] == NULL || graph->term[k] == 0);
         SCIPfreeBufferArrayNull(scip, &(path[k]));
#endif
      }
#if TMX
      SCIPfreeBufferArray(scip, &pathdist);
      SCIPfreeBufferArray(scip, &pathedge);
#else
      SCIPfreeBufferArray(scip, &path);
#endif
   }
   else if( mode == TMPOLZIN )
   {
      SCIPpqueueFree(&pqueue);
      for( k = nnodes - 1; k >= 0; k-- )
      {
         SCIPfreeBuffer(scip, &gnodearr[k]);
         SCIPfreeBufferArray(scip, &node_dist[k]);
         SCIPfreeBufferArray(scip, &node_edge[k]);
         SCIPfreeBufferArray(scip, &node_base[k]);
      }
      SCIPfreeBufferArray(scip, &node_dist);
      SCIPfreeBufferArray(scip, &node_edge);
      SCIPfreeBufferArray(scip, &node_base);
      SCIPfreeBufferArray(scip, &gnodearr);
      SCIPfreeBufferArray(scip, &vcount);
      SCIPfreeBufferArray(scip, &nodenterms);
   }

   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &start);
   SCIPfreeBufferArray(scip, &connected);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTM)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTM(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_PROBDATA* probdata;
   GRAPH* graph;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   if( graph == NULL )
   {
      heurdata->stp_type = STP_UNDIRECTED;
      return SCIP_OKAY;
   }
   heurdata->stp_type = graph->stp_type;

   if( 0 && (graph->stp_type == STP_PRIZE_COLLECTING || graph->stp_type == STP_MAX_NODE_WEIGHT) )
   {
      int e;
      int k = 0;
      int r = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(heurdata->rootedges_t), graph->terms - 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(heurdata->rootedges_z), graph->terms - 1) );

      for( e = graph->inpbeg[graph->source[0]]; e != EAT_LAST; e = graph->ieat[e] )
      {
         if(  Is_term(graph->term[graph->tail[e]]) )
         {
            heurdata->rootedges_t[k++] = e;
         }
         else
         {
            //  costrootedges_z[r] = -graph->cost[flipedge(e + 2)];
            if (graph->head[e+2] != graph->source[0])
               printf("PC root edge missing \n");
            /* if( e + 2 != heurdata->rootedges_t[k-1] )
               printf("NEQ! \n");
               else
               printf("cost: %f \n", costrootedges_z[r]);*/

            heurdata->rootedges_z[r++] = e;

         }
      }
   }
   else
   {
      heurdata->rootedges_t = NULL;
      heurdata->rootedges_z = NULL;
   }
   heurdata->beststartnode = -1;
   heurdata->ncalls = 0;
   heurdata->nexecs = 0;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#else
   heurdata->randseed = 0;
#endif

   SCIPheurSetTimingmask(heur, (SCIP_HEURTIMING) heurdata->timing);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitTM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->stp_type == STP_MAX_NODE_WEIGHT ||  heurdata->stp_type == STP_PRIZE_COLLECTING )
   {
      SCIPfreeMemoryArrayNull(scip, &(heurdata->rootedges_t));
      SCIPfreeMemoryArrayNull(scip, &(heurdata->rootedges_z));
   }

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolTM NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolTM)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of TM primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolTM NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTM)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   GRAPH* graph;
   SCIP_Real* cost;
   SCIP_Real* costrev;
   SCIP_Real* nval;
   SCIP_Real* xval;
   SCIP_Real maxcost = 0.0;
   int* results;
   SCIP_Real pobj;
   int best_start = -1;
   int nedges;
   int edgecount;
   int nvars;
   int layer;
   int runs;
   int e;
   int v;
   int pctrivialbound = 100000;
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;
   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   graph = SCIPprobdataGetGraph(probdata);
   assert(graph != NULL);
   nedges = graph->edges;
   assert(nedges >= 0);
   if( graph->stp_type == STP_PRIZE_COLLECTING && graph->knots > pctrivialbound && !(heurtiming & SCIP_HEURTIMING_BEFORENODE) )
      return SCIP_OKAY;

   runs = 0;

   /* set the runs, i.e. number of different starting points for the heuristic */
   if( heurtiming & SCIP_HEURTIMING_BEFORENODE )
   {
      if( SCIPgetDepth(scip) > 0 )
         return SCIP_OKAY;

      runs = heurdata->initruns;
   }
   else if( ((heurtiming & SCIP_HEURTIMING_DURINGLPLOOP) && (heurdata->ncalls % heurdata->duringlpfreq == 0)) || (heurtiming & SCIP_HEURTIMING_AFTERLPLOOP) )
      runs = heurdata->evalruns;
   else if( heurtiming & SCIP_HEURTIMING_AFTERNODE )
   {
      if( SCIPgetDepth(scip) == 0 )
         runs = heurdata->rootruns;
      else
         runs = heurdata->leafruns;
   }

   /* increase counter for number of (TM) calls */
   heurdata->ncalls++;

   if( runs == 0 )
      return SCIP_OKAY;
   heurdata->nexecs++;

   SCIPdebugMessage("Heuristic Start\n");

   /* get all variables (corresponding to the edges) */
   nvars = SCIPprobdataGetNVars(scip);
   vars = SCIPprobdataGetVars(scip);
   assert(vars != NULL);
   assert(vars[0] != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &results, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nval, nvars) );

   *result = SCIP_DIDNOTFIND;

   /* */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      sol = NULL;
      xval = NULL;
   }
   else
   {
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );

      xval = SCIPprobdataGetXval(scip, sol);

      SCIPfreeSol(scip, &sol);
   }

   /* set (edge) result array to default */
   for( e = 0; e < nedges; e++ )
      results[e] = UNKNOWN;

   /* prize collecting problem and too many edges for the heuristic to handle? */
   if( graph->stp_type == STP_PRIZE_COLLECTING && graph->edges > pctrivialbound )
   {
      printf("tm trivial\n");
      do_prizecoll_trivial(scip, graph, results);
   }
   else
   {
      SCIP_Real rand;
      SCIP_Real randupper;
      SCIP_Real randlower;
      randupper = SCIPgetRandomReal(1.1, 2.5, &(heurdata->randseed));
      randlower = SCIPgetRandomReal(1.1, randupper, &(heurdata->randseed));

      for( layer = 0; layer < 1; layer++ ) /*graph->layers */
      {
         if( xval == NULL )
         {
            int fixed = 0;
            BMScopyMemoryArray(cost, graph->cost, nedges);

            /* hop constraint problem? */
            if( graph->stp_type == STP_HOP_CONS )
            {
               for( e = 0; e < nedges; e++)
               {
                  if( SCIPvarGetUbGlobal(vars[e] ) < 0.5 )
                  {
                     cost[e] = 1e+10;
                  }
                  if( SCIPisLT(scip, graph->cost[e], 1e+8 ) && SCIPisGT(scip, graph->cost[e], maxcost) )
                     maxcost = graph->cost[e];

               }
               for( e = 0; e < nedges; e++)
                  costrev[e] = cost[flipedge(e)];
            }
            else
            {
               /* TODO chg. for asymmetric graphs */
               for( e = 0; e < nedges; e += 2)
               {
                  if( SCIPvarGetUbGlobal(vars[layer * nedges + e + 1]) < 0.5 )
                  {
                     costrev[e] = 1e+10;
                     cost[e + 1] = 1e+10;
                     fixed++;
                  }
                  else
                  {
                     costrev[e] = cost[e + 1];
                     costrev[e + 1] = cost[e];
                  }

                  if( SCIPvarGetUbGlobal(vars[layer * nedges + e]) < 0.5 )
                  {
                     fixed++;
                     costrev[e + 1] = 1e+10;
                     cost[e] = 1e+10;
                  }
                  else
                  {
                     costrev[e] = cost[e + 1];
                     costrev[e + 1] = cost[e];
                  }
               }
            }
            //printf("nvars: %d fixed: %d \n", nvars,  fixed);
         }
         else
         {
            /*
              if( heurdata->nexecs % 10 == 0 )
              printf("tm totally randomized: \n\n");
              if( heurdata->nexecs % 4 == 0 )
              printf("tm partly randomized: (%f  %f)\n\n", randlower, randupper);
	    */
            if( graph->stp_type == STP_HOP_CONS )
            {

               for( e = 0; e < nedges; e++)
               {
                  if( SCIPvarGetUbGlobal(vars[e] ) < 0.5 )
                  {
                     cost[e] = 1e+10;
                  }
                  else
                  {
                     if( heurdata->nexecs % 10 == 0 )
                     {
                        rand = SCIPgetRandomReal(randlower, randupper, &(heurdata->randseed));
                        cost[e] = graph->cost[e] * rand;
                     }
                     else
                     {
                        cost[e] = ((1.0 - xval[e]) * graph->cost[e]);
                     }
                  }
                  if( heurdata->nexecs % 4 == 0 )
                  {
                     rand = SCIPgetRandomReal(randlower, randupper, &(heurdata->randseed));
                     cost[e] = cost[e] * rand;
                  }
                  /* TODO graphcost to cost*/
                  if( SCIPisLT(scip, cost[e], 1e+8 ) && SCIPisGT(scip, cost[e], maxcost) )
                     maxcost = cost[e];
                  assert(SCIPisGE(scip, cost[e], 0));
               }
               for( e = 0; e < nedges; e++)
                  costrev[e] = cost[flipedge(e)];
            }
            else
            {
               /* swap costs; set a high cost if the variable is fixed to 0 */
               for( e = 0; e < nedges; e += 2)
               {
                  rand = SCIPgetRandomReal(randlower, randupper, &(heurdata->randseed));

                  if( SCIPvarGetUbLocal(vars[layer * nedges + e + 1]) < 0.5 )
                  {
                     costrev[e] = 1e+10;
                     cost[e + 1] = 1e+10;
                  }
                  else
                  {
                     if( heurdata->nexecs % 10 == 0 )
                        costrev[e] = graph->cost[e + 1] * rand;
                     else
                        costrev[e] = ((1.0 - xval[layer * nedges + e + 1]) * graph->cost[e + 1]);

                     if( heurdata->nexecs % 4 == 0 )
                     {
                        costrev[e] = costrev[e] * rand;
                     }
                     cost[e + 1] = costrev[e];
                  }

                  if( SCIPvarGetUbLocal(vars[layer * nedges + e]) < 0.5 )
                  {
                     costrev[e + 1] = 1e+10;
                     cost[e] = 1e+10;
                  }
                  else
                  {
                     if( heurdata->nexecs % 10 == 0 )
                        costrev[e + 1] = graph->cost[e] * rand;
                     else
                        costrev[e + 1] = ((1.0 - xval[layer * nedges + e]) * graph->cost[e]);

                     if( heurdata->nexecs % 4 == 0 )
                        costrev[e + 1] = costrev[e + 1]  * rand;
                     cost[e] = costrev[e + 1];
                  }
                  assert(SCIPisGE(scip, cost[e], 0));
                  assert(SCIPisGE(scip, costrev[e], 0));
               }
            }
         }

         /*
           if( graph->stp_type == STP_HOP_CONS )
           {
           for( e = 0; e < nedges; e++ )
           {
           cost[e] = 1 + cost[e] / maxcost;
           costrev[e] = 1 + costrev[e] / maxcost;
           }
           }
           hop constraint problem? */
         /* can we connect the network */
         SCIP_CALL( do_layer(scip, heurdata, graph, NULL, &best_start, results, runs, heurdata->beststartnode, cost, costrev, maxcost) );
#if 0
         /* take the path */
         if( graph->layers > 1 )
         {
            for( e = 0; e < nedges; e += 2)
            {
               if( (results[e] == layer) || (results[e + 1] == layer) )
                  graph_edge_hide(graph, e);
            }
         }
#endif
      }
#if 0
      if( graph->layers > 1 )
         graph_uncover(graph);
#endif
   }
   edgecount = 0;
   for( v = 0; v < nvars; v++ )
   {
      nval[v] = (results[v % nedges] == (v / nedges)) ? 1.0 : 0.0;
      if( SCIPisEQ(scip, nval[v], 1.0) )
         edgecount++;
   }

   if( graph->stp_type == STP_HOP_CONS )
   {
      //printf("edgecount: %d root: %d \n ", edgecount, graph->source[0]);
      for( v = graph->inpbeg[ graph->source[0]]; v != EAT_LAST; v = graph->ieat[v] )
         assert( SCIPisEQ(scip, nval[v], 0.0) );
      /*
        if( edgecount < graph->hoplimit )
        {
        heurdata->hopfactor = heurdata->hopfactor * (1.0 - MIN(0.5, (double)(graph->hoplimit - edgecount) / 10.0) );
        assert(heurdata->hopfactor > 0);
        printf("reduced hopfactor: %f \n ", heurdata->hopfactor );
        }*/
   }
   if( validate(graph, nval) )
   {
      pobj = 0.0;

      for( v = 0; v < nvars; v++ )
         pobj += graph->cost[v % nedges] * nval[v];

      //assert(graph_valid2(scip, graph, cost));
      printf("tm: %f \n", pobj + SCIPprobdataGetOffset(scip));
      if( SCIPisLE(scip, pobj, SCIPgetPrimalbound(scip)) )
      {
         heurdata->beststartnode = best_start;
      }
      if( 1 || SCIPisLE(scip, pobj, SCIPgetPrimalbound(scip)) ) /* best_start != graph->source[0] || */
      {
         SCIP_Bool success;

         SCIP_CALL( SCIPprobdataAddNewSol(scip, nval, sol, heur, &success) );

         if( success )
            *result = SCIP_FOUNDSOL;
      }
   }

   SCIPfreeBufferArray(scip, &nval);
   SCIPfreeBufferArray(scip, &results);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &costrev);
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the TM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
   char paramdesc[SCIP_MAXSTRLEN];

   /* create TM primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTM, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTM) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTM) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTM) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitTM) );
#if 0
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolTM) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolTM) );
#endif
   heurdata->ncalls = 0;
   heurdata->nexecs = 0;

#ifdef WITH_UG
   heurdata->randseed += getUgRank();
#else
   heurdata->randseed = 0;
#endif
   /* add TM primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/evalruns",
         "number of runs for eval",
         &heurdata->evalruns, FALSE, DEFAULT_EVALRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/randseed",
         "random seed for heuristic",
         NULL, FALSE, DEFAULT_RANDSEED, 0, INT_MAX, paramChgdRandomseed, (SCIP_PARAMDATA*)heurdata) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/initruns",
         "number of runs for init",
         &heurdata->initruns, FALSE, DEFAULT_INITRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/leafruns",
         "number of runs for leaf",
         &heurdata->leafruns, FALSE, DEFAULT_LEAFRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/rootruns",
         "number of runs for root",
         &heurdata->rootruns, FALSE, DEFAULT_ROOTRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/duringlpfreq",
         "frequency for calling heuristic during LP loop",
         &heurdata->duringlpfreq, FALSE, DEFAULT_DURINGLPFREQ, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/type",
         "Heuristic: 0 automatic, 1 TM, 2 TMPOLZIN",
         NULL, FALSE, DEFAULT_TYPE, 0, 2, NULL, NULL) );
   heurdata->hopfactor = 0.33;

   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when heuristc should be called (%u:BEFORENODE, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:AFTERNODE)", SCIP_HEURTIMING_BEFORENODE, SCIP_HEURTIMING_DURINGLPLOOP, SCIP_HEURTIMING_AFTERLPLOOP, SCIP_HEURTIMING_AFTERNODE);
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/timing", paramdesc,
         (int*) &heurdata->timing, TRUE, (int) HEUR_TIMING, (int) SCIP_HEURTIMING_BEFORENODE, 2 * (int) SCIP_HEURTIMING_AFTERNODE - 1, NULL, NULL) ); /*lint !e713*/

   return SCIP_OKAY;
}
