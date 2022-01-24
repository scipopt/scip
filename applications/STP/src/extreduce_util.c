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

/**@file   extreduce_util.c
 * @brief  utility methods for Steiner tree extended reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem extended reduction techniques.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// #define SCIP_DEBUG
#include "extreduce.h"
#include "misc_stp.h"
#include "portab.h"


/**@name Local methods
 *
 * @{
 */




/**@} */

/**@name Interface methods
 *
 * @{
 */



/** reverts extension component (makes root the only extension leaf) */
void extreduce_extCompRevert(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   EXTCOMP*              extcomp             /**< component to be cleaned for */
)
{
   int* const extleaves = extcomp->extleaves;
   int* const compedges = extcomp->compedges;
   const int comproot_org = extcomp->comproot;
   const SCIP_Bool compIsSingleEdge = (extcomp->ncompedges == 1);

   assert(extcomp->ncompedges >= 1);
   assert(extcomp->comproot == graph->tail[extcomp->compedges[0]]);

   if( compIsSingleEdge )
   {
      assert(extcomp->nextleaves == 1);
   }
   else
   {
      const int firstedge = compedges[0];

      assert(extcomp->ncompedges >= 3);
      assert(extcomp->nextleaves >= 2);
      assert(graph->head[compedges[1]] == extleaves[0]);

      compedges[0] = compedges[1];
      compedges[1] = flipedge(firstedge);
   }

   compedges[0] = flipedge(compedges[0]);

   assert(extleaves[0] != extcomp->comproot);

   extcomp->comproot = extleaves[0];
   extleaves[0] = comproot_org;
   extcomp->nextleaves = 1;
}


/** is extension component promising candidate for extension? */
SCIP_Bool extreduce_extCompIsPromising(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   const EXTCOMP*        extcomp             /**< component to be cleaned for */)
{
   const int* const extleaves = extcomp->extleaves;
   const SCIP_Bool* const isterm = extperma->isterm;
   const int nextleaves = extcomp->nextleaves;

   assert(extcomp->ncompedges == 1 || extcomp->ncompedges >= 3);

   /* we want to at least check each star! */
   if( extcomp->ncompedges >= 3 )
   {
      return TRUE;
   }

   /* go over all possible extensions and see whether any of them are promising */
   for( int i = 0; i < nextleaves; ++i )
   {
      const int leaf = extleaves[i];

      if( extLeafIsExtendable(graph, isterm, leaf) )
         return TRUE;
   }

   return FALSE;
}


/** is extension component or its reversed version promising candidate for extension? */
SCIP_Bool extreduce_extCompFullIsPromising(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   const EXTCOMP*        extcomp             /**< component to be cleaned for */)
{
   assert(graph && extperma && extcomp);

   if( extreduce_extCompIsPromising(graph, extperma, extcomp) )
      return TRUE;

   if( extLeafIsExtendable(graph, extperma->isterm, extcomp->comproot) )
      return TRUE;

   return FALSE;
}

/**@} */
