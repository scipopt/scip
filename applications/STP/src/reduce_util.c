/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_util.c
 * @brief  utility methods for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem reduction techniques.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include "reduce.h"
#include "portab.h"


/** storage for edge on complete graph */
typedef struct complete_edge
{
   int                   tail;              /**< tail vertex */
   int                   head;              /**< head vertex */
   SCIP_Real             cost;              /**< edge cost */
} CEDGE;


/** lightweight minimum spanning tree structure that allows to add vertices to given MST on complete graph (in CSR format) */
struct dynamic_complete_minimum_spanning_tree
{
   CEDGE*                edgestore;         /**< storage for edges (of size maxnnodes) */
   SCIP_Real*            adjcost_buffer;    /**< distances buffer (of size maxnnodes) */
   SCIP_Bool*            nodemark;          /**< array for marking nodes (of size maxnnodes) */
   int                   maxnnodes;         /**< maximum number of nodes that can be handled */
};


/** recursive method for adding node to MST */
static
void dcmstInsert(
   const CSR*            org_mst,            /**< the base MST */
   const SCIP_Real       adjcosts[],         /**< (undirected) adjacency costs for new node */
   int                   root,               /**< the current root */
   CEDGE                 new_mst[],          /**< new MST */
   SCIP_Bool             new_nodemarked[],   /**< array */
   CEDGE*                max_path_edge,      /**< pointer to maximum edge on path to new node */
   int*                  new_nedges          /**< pointer to current number of edges */
)
{
   CEDGE root2new = { .tail = root, .head = org_mst->nnodes, .cost = adjcosts[root] };
   const int* const org_start = org_mst->start;
   const int* const org_head = org_mst->head;
   const SCIP_Real* const org_cost = org_mst->cost;

   assert(new_nodemarked[root]);

   /* visit all neighbors or root in the original MST */
   for( int i = org_start[root]; i != org_start[root + 1]; ++i )
   {
      const int w = org_head[i];

      /* node not visited yet? */
      if( !new_nodemarked[w] )
      {
         const SCIP_Real costroot2w = org_cost[i];

         new_nodemarked[w] = TRUE;
         dcmstInsert(org_mst, adjcosts, w, new_mst, new_nodemarked, max_path_edge, new_nedges);

         assert(max_path_edge->tail >= 0);
         assert(*new_nedges >= 0 && *new_nedges < org_mst->nnodes);

         if( max_path_edge->cost < costroot2w )
         {
            new_mst[(*new_nedges)++] = *max_path_edge;

            if( costroot2w < root2new.cost )
            {
               root2new.tail = root;
               root2new.head = w;
               root2new.cost = costroot2w;
            }
         }
         else
         {
            const int nedges = (*new_nedges);

            new_mst[nedges].tail = root;
            new_mst[nedges].head = w;
            new_mst[nedges].cost = costroot2w;

            (*new_nedges)++;

            if( max_path_edge->cost < root2new.cost )
            {
               root2new = *max_path_edge;
            }
         }
      }
   }

   *max_path_edge = root2new;
}


/** add node to MST */
static inline
void dcmstAddNode(
   const CSR*            mst_in,             /**< source */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst                /**< underlying structure */
)
{
   CEDGE max_path_edge = { .tail = -1, .head = -1, .cost = -1.0 };
   CEDGE* const edgestore = dmst->edgestore;
   SCIP_Bool* const nodemark = dmst->nodemark;
   int nedges_new = 0;
   const int nnodes_in = mst_in->nnodes;

   assert(nnodes_in >= 1);

   nodemark[0] = TRUE;

   for( int i = 1; i < nnodes_in; ++i )
      nodemark[i] = FALSE;

   dcmstInsert(mst_in, adjcosts, 0, edgestore, nodemark, &max_path_edge, &nedges_new);

   assert(nedges_new == nnodes_in - 1);

   edgestore[nedges_new] = max_path_edge;
}


/** transforms edge-store to CSR  */
static inline
void dcmstGetCSRfromStore(
   const DCMST*          dmst,               /**< underlying structure */
   CSR*                  mst_out             /**< target */
)
{
   const CEDGE* const edgestore = dmst->edgestore;
   int* const mst_start = mst_out->start;
   int* const mst_head = mst_out->head;
   SCIP_Real* const mst_cost = mst_out->cost;
   const int mst_nnodes = mst_out->nnodes;

   /* undirected edges */
   const int mst_nedges = mst_nnodes - 1;

   assert(mst_nnodes <= dmst->maxnnodes);
   assert(2 * mst_nedges == mst_out->nedges);

   BMSclearMemoryArray(mst_start, mst_nnodes + 1);

   for( int i = 0; i < mst_nedges; ++i )
   {
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;

      assert(v1 >= 0 && v1 < mst_nnodes);
      assert(v2 >= 0 && v2 < mst_nnodes);

      mst_start[v1]++;
      mst_start[v2]++;
   }

   assert(mst_start[mst_nnodes] == 0);

   for( int i = 1; i <= mst_nnodes; ++i )
   {
      mst_start[i] += mst_start[i - 1];
   }

   assert(mst_start[mst_nnodes] == mst_out->nedges);

   for( int i = 0; i < mst_nedges; ++i )
   {
      const int v1 = edgestore[i].tail;
      const int v2 = edgestore[i].head;
      const SCIP_Real cost = edgestore[i].cost;

      assert(mst_start[v1] >= 1);
      assert(mst_start[v2] >= 1);

      mst_head[--mst_start[v1]] = v2;
      mst_cost[mst_start[v1]] = cost;

      mst_head[--mst_start[v2]] = v1;
      mst_cost[mst_start[v2]] = cost;
   }
}


/** initializes dynamic MST structure */
SCIP_RETCODE reduce_dcmstInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnnodes,          /**< maximum number of nodes that can be handled */
   DCMST**               dcmst               /**< to be initialized */
)
{
   DCMST* mst;

   assert(scip && dcmst);
   assert(maxnnodes >= 1);

   SCIP_CALL( SCIPallocMemory(scip, dcmst) );

   mst = *dcmst;

   mst->maxnnodes = maxnnodes;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->edgestore), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->adjcost_buffer), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->nodemark), maxnnodes) );


   return SCIP_OKAY;
}


/** adds node to CSR "mst_in" and saves result in "mst_out" */
void reduce_dcmstAddNode(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst_in,             /**< source */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst,               /**< underlying structure */
   CSR*                  mst_out             /**< target */
)
{
   assert(mst_in && adjcosts && dmst && mst_out);

   assert(reduce_dcmstMstIsValid(scip, mst_in));

   assert(mst_out->nnodes == mst_in->nnodes + 1);
   assert(mst_out->nedges == mst_in->nedges + 2);
   assert(mst_in->nnodes < dmst->maxnnodes);

   dcmstAddNode(mst_in, adjcosts, dmst);

   dcmstGetCSRfromStore(dmst, mst_out);

   assert(mst_out->nnodes == mst_in->nnodes + 1);
   assert(reduce_dcmstMstIsValid(scip, mst_out));
}


/** Adds node to CSR "mst".
 *  NOTE: There needs to be enough space in CSR arrays for one more node! */
void reduce_dcmstAddNodeInplace(
   SCIP*                 scip,               /**< SCIP */
   const SCIP_Real*      adjcosts,           /**< (undirected) adjacency costs for new node */
   DCMST*                dmst,               /**< underlying structure */
   CSR*                  mst                 /**< source/target */
)
{
   assert(mst && adjcosts && dmst);

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(mst->nnodes < dmst->maxnnodes);

   dcmstAddNode(mst, adjcosts, dmst);

   mst->nnodes += 1;
   mst->nedges += 2;

   dcmstGetCSRfromStore(dmst, mst);

   assert(reduce_dcmstMstIsValid(scip, mst));
}

/** computes MST on 0 node */
void reduce_dcmstGet0NodeMst(
   SCIP*                 scip,               /**< SCIP */
   CSR*                  mst                 /**< MST */
)
{
   int* const start = mst->start;

   assert(mst->nnodes == 0);
   assert(mst->nedges == 0);

   start[0] = 0;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(0.0, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 1 node */
void reduce_dcmstGet1NodeMst(
   SCIP*                 scip,               /**< SCIP */
   CSR*                  mst                 /**< MST */
)
{
   int* const start = mst->start;

   assert(mst->nnodes == 1);
   assert(mst->nedges == 0);

   start[0] = 0;
   start[1] = 0;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(0.0, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 2 nodes */
void reduce_dcmstGet2NodeMst(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real             edgecost,           /**< edge cost */
   CSR*                  mst                 /**< MST */
)
{
   SCIP_Real* const cost = mst->cost;
   int* const start = mst->start;
   int* const head = mst->head;

   assert(edgecost > 0.0);
   assert(mst->nnodes == 2);
   assert(mst->nedges == 2);

   start[0] = 0;
   start[1] = 1;
   start[2] = 2;

   head[0] = 1;
   head[1] = 0;

   cost[0] = edgecost;
   cost[1] = edgecost;

   assert(reduce_dcmstMstIsValid(scip, mst));
   assert(EQ(edgecost, reduce_dcmstGetWeight(scip, mst)));
}


/** computes MST on 3 nodes */
void reduce_dcmstGet3NodeMst(
   SCIP*                 scip,               /**< SCIP */
   SCIP_Real             edgecost01,         /**< edge cost */
   SCIP_Real             edgecost02,         /**< edge cost */
   SCIP_Real             edgecost12,         /**< edge cost */
   CSR*                  mst                 /**< MST */
)
{
   assert(0 && "implement me");
}


/** gets weight of MST */
SCIP_Real reduce_dcmstGetWeight(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            mst_in              /**< source */
)
{
   SCIP_Real weight = 0.0;
   const int nedges = mst_in->nedges;
   const SCIP_Real* cost = mst_in->cost;

   assert(scip);
   assert(reduce_dcmstMstIsValid(scip, mst_in));

   for( int i = 0; i < nedges; i++ )
   {
      assert(cost[i] >= 0.0);

      weight += cost[i];
   }

   weight /= 2.0;

   assert(GE(weight, 0.0));

   if( GT(weight, FARAWAY) )
      weight = FARAWAY;

   return weight;
}


/** returns maximum number of nodes */
int reduce_dcmstGetMaxnnodes(
   const DCMST*          dmst                /**< underlying structure */
)
{
   assert(dmst);

   return dmst->maxnnodes;
}


/** Returns buffer of size 'reduce_dcmstGetMaxnnodes'.
  * NOTE: buffer is never used within any other function, apart from allocation and freeing.
  * NOTE: in debug mode the array is initialized to -1.0 */
SCIP_Real* reduce_dcmstGetAdjcostBuffer(
   const DCMST*          dmst                /**< underlying structure */
)
{
   assert(dmst);
   assert(dmst->adjcost_buffer);

#ifndef NDEBUG
   for( int i = 0; i < dmst->maxnnodes; i++ )
      dmst->adjcost_buffer[i] = -1.0;
#endif

   return dmst->adjcost_buffer;
}


/** frees dynamic MST structure */
void reduce_dcmstFree(
   SCIP*                 scip,               /**< SCIP */
   DCMST**               dcmst               /**< to be initialized */
)
{
   assert(scip && dcmst);

   SCIPfreeMemoryArray(scip, &((*dcmst)->nodemark));
   SCIPfreeMemoryArray(scip, &((*dcmst)->adjcost_buffer));
   SCIPfreeMemoryArray(scip, &((*dcmst)->edgestore));

   SCIPfreeMemory(scip, dcmst);
}


/** is the CSR a valid MST on any underlying graph (with number of nodes and edges of the CSR)? */
SCIP_Bool reduce_dcmstMstIsValid(
   SCIP*                 scip,               /**< SCIP */
   const CSR*            cmst                /**< the MST candidate */
)
{
   SCIP_Bool* visited;
   const int* const start_csr = cmst->start;
   const int* const head_csr = cmst->head;
   const int nnodes = cmst->nnodes;
   SCIP_Bool isValid = TRUE;

   if( nnodes == 0 )
   {
      assert(cmst->nedges == 0);
      assert(start_csr[0] == 0);

      return TRUE;
   }

   assert(nnodes >= 1);
   assert(cmst->nedges % 2 == 0);
   assert(start_csr[0] == 0);

   if( !graph_csr_isValid(cmst, FALSE) )
   {
      SCIPdebugMessage("CSR is broken! \n");
      return FALSE;
   }

   if( cmst->nnodes != (cmst->nedges / 2) + 1 )
   {
      SCIPdebugMessage("wrong nodes/edges ratio \n");
      return FALSE;
   }

   if( nnodes == 1 )
   {
      return TRUE;
   }

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &visited, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      visited[i] = FALSE;

   for( int i = 0; i < cmst->nedges; i++ )
   {
      const int head = head_csr[i];

      assert(head >= 0 && head < nnodes);

      visited[head] = TRUE;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !visited[i] )
      {
         SCIPdebugMessage("mst does not contain node %d \n", i);

         isValid = FALSE;
         break;
      }
   }

   SCIPfreeMemoryArray(scip, &visited);

   return isValid;
}
