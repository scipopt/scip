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

/**@file   portab.h
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
#include "grph.h"
#include "portab.h"

#if 0
static void show(
   GRAPH*  g,
   int     vars,
   char*   state,
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
#if 1
/*---------------------------------------------------------------------------*/
/*--- Name     : Trail                                                    ---*/
/*--- Function : Durchlaeuft einen Graphen entsprechend einer Loesung und ---*/
/*---            stellt fest ob er bei allen Knoten vorbeikommt           ---*/
/*--- Parameter: Startknoten, Loesung, Herkunft, Schongewesenliste        ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
static void trail(
   const GRAPH*  g,
   int           i,
   const double* xval,
   int           tail,
   char*         connected,
   int           hop,
   int           max_hops)
{
   int k;

   assert(connected[i] >= 0);

   if( connected[i] < 2 )
   {
      ++(connected[i]);

      if ((connected[i] < 2) && (hop < max_hops))
      {
         for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
            if ((g->head[k] != tail) && (xval[k] + EPSILON > 1.0))
               trail(g, g->head[k], xval, i, connected, hop + 1, max_hops);
      }
   }
}
#endif
#if 0
/*---------------------------------------------------------------------------*/
/*--- Name     : Trail                                                    ---*/
/*--- Function : Durchlaeuft einen Graphen entsprechend einer Loesung und ---*/
/*---            stellt fest ob er bei allen Knoten vorbeikommt           ---*/
/*--- Parameter: Startknoten, Loesung, Herkunft, Schongewesenliste        ---*/
/*--- Returns  : Nichts                                                   ---*/
/*---------------------------------------------------------------------------*/
static void trail2(
   const GRAPH*  g,
   int           layer,
   const double* xval,
   char*         connected,
   int           max_hops)
{
   int* stackstart = malloc((size_t)g->knots * sizeof(int));
   int* stackedge = malloc((size_t)g->knots * sizeof(int));
   int* stacktail = malloc((size_t)g->knots * sizeof(int));
   int k;
   int i = g->source[layer];
   int stacksize = 1;

   stackstart[0] = i;
   stackedge[0] = g->outbeg[i];
   stacktail[0] = -1;

   while( stacksize > 0 )
   {
      k = stackedge[stacksize-1];

      printf("stacksize: %d, edge %d: %d --> %d\n", stacksize, k, g->tail[k], g->head[k]);

      if ( k == EAT_LAST )
      {
         --stacksize;
      }
      else if ((g->head[k] != stacktail[stacksize-1]) && (xval[k] + EPSILON > 1.0))
      {
         assert(connected[g->head[k]] >= 0);
         if( connected[g->head[k]] < 2 )
         {
            ++(connected[g->head[k]]);

            if ((connected[i] < 2) && (stacksize < max_hops + 1))
            {
               stackstart[stacksize] = g->head[k];
               stackedge[stacksize] = g->outbeg[g->head[k]];
               stacktail[stacksize] = stackstart[stacksize-1];

               printf("stack: size=%d start[%d]=%d, edge[%d]=%d, tail[%d]=%d\n",
                  stacksize+1, stacksize, stackstart[stacksize], stacksize, stackedge[stacksize],
                  stacksize, stacktail[stacksize]);

               ++stacksize;
               continue;
            }
         }
      }
      if (stacksize > 0)
         stackedge[stacksize-1] = g->oeat[k];
   }

   free(stacktail);
   free(stackedge);
   free(stackstart);

   printf("done\n");
}
#endif
/*---------------------------------------------------------------------------*/
/*--- Name     : Validate Solution                                        ---*/
/*--- Function : Stellt fuer eine (Teil-)Loesung fest, ob sie zulaessig   ---*/
/*---            unzulaessig oder schon eine vollstaendige Loesung ist    ---*/
/*--- Parameter: Graph, Loesung.                                          ---*/
/*--- Returns  : TRUE / FALSE                                             ---*/
/*---------------------------------------------------------------------------*/
/** validates whether a (LP) solution is feasible */
SCIP_RETCODE SCIPvalidateStpSol(
   SCIP* scip,
   const GRAPH*  g,
   const double* xval,
   SCIP_Bool*    feasible
		     )
{
   char* connected;
   int   ret       = TRUE;
   int   i;
   int   layer;
   int   deg;
   int   e;
   assert(g         != NULL);

   if( g->knots == 1 )
   {
      *feasible = TRUE;
      return SCIP_OKAY;
   }

   assert(xval      != NULL);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &connected, g->knots) );

   *feasible = FALSE;
   for(layer = 0; ret && (layer < g->layers); layer++)
   {
#if 0
      if (layer > 0)
         memset(connected, 0, (size_t)g->knots * sizeof(char));
#endif
#if 1
      trail(g, g->source[layer], xval + layer * g->edges, -1,
         connected,
         0, 1000000000);
#else
      trail2(g, layer, xval + layer * g->edges, connected, 1000000000);
#endif
      for(i = 0; i < g->knots; i++)
      {
         if( g->stp_type == STP_DEG_CONS )
         {
            deg = 0;
            for( e = g->outbeg[i]; e != EAT_LAST ; e = g->oeat[e] )
               if( EQ(xval[e], 1.0) || EQ(xval[flipedge(e)], 1.0) )
                  deg++;

            if( deg > g->maxdeg[i] )
            {
               ret = FALSE;
               SCIPdebugMessage("deg condition violated \n");
               break;
            }

         }
         /* circle? */
         if( connected[i] >= 2)
         {
            ret = FALSE;
            break;
         }

         /* vertex not reached? */
         if ((g->grad[i] > 0) && (g->term[i] == layer) && !connected[i])
         {
            ret = FALSE;
            break;
         }
      }
   }
   SCIPfreeBufferArray(scip, &connected);

   if( ret )
      ret = nail(g, xval);

   /* SCIP-Heuristiken können (z.B. durch Runden) Lösungen konstruieren, die einen Kreis aus Steiner-Knoten enthalten, der
    * nicht zum Baum gehört
    */
#if 0
#ifndef NDEBUG
   /* Test ob alle Kanten nur in eine Richtung benutzt werden.
    */
   if (ret)
   {
      int e;

      for(e = 0; e < g->edges; e += 2)
         if ((xval[e] > EPSILON) && (xval[e + 1] > EPSILON))
            /* CONSTCOND */
            assert(FALSE);
   }
#endif
#endif
   *feasible = (SCIP_Bool) ret;
   return SCIP_OKAY;
}
