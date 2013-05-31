/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: validate.c                                                    */
/*   Name....: STP Solution validation                                       */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

   if ((++connected[i] < 2) && (hop < max_hops))
      for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
         if ((g->head[k] != tail) && (xval[k] + EPSILON > 1.0))
            trail(g, g->head[k], xval, i, connected, hop + 1, max_hops);
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Validate Solution                                        ---*/
/*--- Function : Stellt fuer eine (Teil-)Loesung fest, ob sie zulaessig   ---*/
/*---            unzulaessig oder schon eine vollstaendige Loesung ist    ---*/
/*--- Parameter: Graph, Loesung.                                          ---*/
/*--- Returns  : TRUE / FALSE                                             ---*/
/*---------------------------------------------------------------------------*/
int validate(
   const GRAPH*  g,
   const double* xval)
{
   char* connected = malloc((size_t)g->knots * sizeof(char));
   int   ret       = TRUE;
   int   i;
   int   layer;

   assert(g         != NULL);
   assert(xval      != NULL);
   assert(connected != NULL);

   printf("source=%d\n", g->source[0]);

   for(layer = 0; ret && (layer < g->layers); layer++)
   {
#if 0
      trail(g, g->source[layer], xval + layer * g->edges, -1,
         memset(connected, 0, (size_t)g->knots * sizeof(char)),
         0, param_get("MAX_HOPS")->i);
#endif
      trail(g, g->source[layer], xval + layer * g->edges, -1,
         memset(connected, 0, (size_t)g->knots * sizeof(char)),
         0, 1000000000);


      for(i = 0; i < g->knots; i++)
      {
         /* Etwa ein Kreis ?
          */
         if (connected[i] >= 2)
            ret = FALSE;

         /* Jemand noch einsam
          */
         if ((g->grad[i] > 0) && (g->term[i] == layer) && !connected[i])
            ret = FALSE;
      }
   }
   free(connected);

   if (ret)
      ret = nail(g, xval);

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
   return(ret);
}
