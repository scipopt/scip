/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   validate.c
 * @brief  Method to validate Steiner problem solutions
 * @author Thorsten Koch
 * @author Gerald Gamrath
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"

#if 0
static void show(
   GRAPH*  g,
   int     vars,
   STP_Bool*   state,
   double* xval)
{
   int i;
   for(i = 0; i < vars; i++)
   {
      if (xval[i] > 1e-6)
      {
	 fprintf(stderr, "%d-%d, xval[%d]=%g (%d)\n",
            g->tail[i % g->edges] + 1,
            g->head[i % g->edges] + 1,
            i,
            xval[i],
            state[i]);
      }
   }
}
#endif

static int nail(
   const GRAPH*  g,
   const double* xval)
{
   double  sum;
   int     i;
   int     l;
   int     ind;

   assert(g      != NULL);
   assert(xval   != NULL);

   if (g->layers == 1)
      return(TRUE);

   assert(0 && "currently not used!");

   for(i = 0; i < g->edges; i += 2)
   {
      sum   = 0.0;

      for(l = 0; l < g->layers; l++)
      {
         ind  = l * g->edges + i;
         sum += xval[ind] + xval[ind + 1];
      }
      assert(sum >= -EPSILON);

#ifdef WITH_CAPACITIES
      assert(g->capa[i] > 0);

      if (sum - EPSILON > (double)g->capa[i])
         return(FALSE);
#else
      if (sum - EPSILON > 1.0)
         return(FALSE);
#endif
   }

   return(TRUE);
}

#ifndef NDEBUG
/*---------------------------------------------------------------------------*/
/*--- Name     : Trail                                                    ---*/
/*--- Function : Durchlaeuft einen Graphen entsprechend einer Loesung und ---*/
/*---            stellt fest ob er bei allen Knoten vorbeikommt           ---*/
/*--- Parameter: Startknoten, Loesung, Herkunft, Schongewesenliste        ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
static void trail_old(
   const GRAPH*  g,
   int           i,
   const double* xval,
   int           tail,
   int*          visitcount
   )
{
   /* node visited at most one time so far? */
   if( visitcount[i] <= 1 )
   {
      assert(visitcount[i] == 0 || visitcount[i] == 1);

      ++(visitcount[i]);

      /* node visited the first time? */
      if( visitcount[i] <= 1 )
      {
         assert(visitcount[i] == 1);

         for( int k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k] )
         {
            if( (xval[k] + EPSILON > 1.0) ) // && g->head[k] != tail
            {
               trail_old(g, g->head[k], xval, i, visitcount);
            }
         }
      }
   }
}
#endif


/** traverses the graph from vertex 'start' and marks all reached nodes (counts up to to at most 2) */
static
SCIP_RETCODE trail(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the new graph */
   const double*         xval,               /**< (LP) solution */
   int                   start,              /**< node to start from */
   int*                  visitcount          /**< marks which node has been visited */
   )
{
   int* stackarr;
   int stacksize;
   const int nnodes = g->knots;

   for( int i = 0; i < nnodes; i++ )
      visitcount[i] = 0;

   visitcount[start] = 1;

   if( g->grad[start] == 0 )
      return SCIP_OKAY;

   stacksize = 0;

   SCIP_CALL(SCIPallocBufferArray(scip, &stackarr, nnodes));

   stackarr[stacksize++] = start;

   /* DFS loop */
   while( stacksize != 0 )
   {
      const int node = stackarr[--stacksize];

      /* traverse outgoing arcs */
      for( int a = g->outbeg[node]; a != EAT_LAST; a = g->oeat[a] )
      {
         if( (xval[a] + EPSILON > 1.0) )
         {
            const int head = g->head[a];

            /* not visited yet? */
            if( visitcount[head] == 0 )
              stackarr[stacksize++] = head;

            if( visitcount[head] <= 1 )
               visitcount[head]++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &stackarr);

#ifndef NDEBUG
   if( nnodes <= 10000 )
   {
      int* visitcount_dbg;
      SCIP_CALL( SCIPallocClearBufferArray(scip, &visitcount_dbg, nnodes) );

      trail_old(g, start, xval, -1, visitcount_dbg);

      for( int i = 0; i < nnodes; i++ )
      {
         assert(visitcount_dbg[i] == visitcount[i]);
      }

      SCIPfreeBufferArray(scip, &visitcount_dbg);
   }
#endif

   return SCIP_OKAY;
}





/*---------------------------------------------------------------------------*/
/*--- Name     : Validate Solution                                        ---*/
/*--- Function : Stellt fuer eine (Teil-)Loesung fest, ob sie zulaessig   ---*/
/*---            unzulaessig oder schon eine vollstaendige Loesung ist    ---*/
/*--- Parameter: Graph, Loesung.                                          ---*/
/*--- Returns  : TRUE / FALSE                                             ---*/
/*---------------------------------------------------------------------------*/
/** validates whether a (LP) solution is feasible */
SCIP_RETCODE SCIPStpValidateSol(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the new graph */
   const double*         xval,               /**< (LP) solution */
   SCIP_Bool             allow_cyles,        /**< allow cycles? */
   SCIP_Bool*            feasible            /**< is feasible? */
)
{
   int* visitcount;
   SCIP_Bool isValid = TRUE;
   const int nnodes = graph_get_nNodes(g);

   assert(scip && xval && feasible);
   assert(graph_valid(scip, g));

   if( nnodes == 1 )
   {
      *feasible = TRUE;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &visitcount, nnodes) );

   trail(scip, g, xval, g->source, visitcount);

   /* check whether solution is feasible */
   for( int i = 0; i < nnodes; i++ )
   {
      if( g->stp_type == STP_DCSTP )
      {
         int deg = 0;

         for( int e = g->outbeg[i]; e != EAT_LAST ; e = g->oeat[e] )
         {
            if( GE(xval[e], 1.0) || GE(xval[flipedge(e)], 1.0) )
               deg++;
         }

         if( deg > g->maxdeg[i] )
         {
            isValid = FALSE;
            SCIPdebugMessage("deg condition violated \n");

            break;
         }
      }

      /* with cycle? */
      if( !allow_cyles && visitcount[i] >= 2 )
      {
         isValid = FALSE;
         SCIPdebugMessage("cycle found \n");

         break;
      }

      /* terminal not reached? */
      if ((g->grad[i] > 0) && (g->term[i] == 0) && visitcount[i] == 0 )
      {
         isValid = FALSE;
         SCIPdebugMessage("terminal %d not reached \n", i);

         break;
      }
   }

   SCIPfreeBufferArray(scip, &visitcount);

   if( isValid )
      isValid = nail(g, xval);

   /* SCIP-Heuristiken können (z.B. durch Runden) Lösungen konstruieren, die einen Kreis aus Steiner-Knoten enthalten, der
    * nicht zum Baum gehört
    */
#ifndef NDEBUG
   /* Teste, ob alle Kanten nur in eine Richtung benutzt werden.
    */
   if( isValid && !allow_cyles )
   {
      for( int e = 0; e < g->edges; e += 2)
      {
         assert((xval[e] <= EPSILON) || (xval[e + 1] <= EPSILON));
      }
   }
#endif

   *feasible = isValid;

   return SCIP_OKAY;
}
