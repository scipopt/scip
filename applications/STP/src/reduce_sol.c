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

/**@file   reduce_sol.c
 * @brief  Reduction solution storage methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes methods to save and retain solutions and sollocal bounds during the reduction
 * process
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h) */
//#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "heur_tm.h"
#include "solstp.h"
#include "scip/scip.h"
#include "portab.h"
#include "stpvector.h"


#define REDSOLVAL_UNSET (-FARAWAY)

struct reduction_local_solution_storage
{
   int*                  nodesol;
   SCIP_Real             nodesol_ub;
   SCIP_Real             offset;             /**< offset NOTE: offset is only used within this structure! */
   SCIP_Real             primalbound;        /**< best sollocal bound */
   int                   nnodes;
   SCIP_Bool             isPcMw;
   SCIP_Bool             nodesol_use;
};


/** single level representation (needed for decomposition of problem graph) */
typedef struct reduction_solution_level
{
   REDSOLLOCAL*          redsollocal;        /**< local storage */
   int*                  nodesol;            /**< solution nodes marker */
   SCIP_Real             solval_postred;     /**< value after reduction loop */
   SCIP_Real             solval_incomplete;  /**< incomplete value (to be summed up during decomposition ) */
   int                   nnodes;
   SCIP_Real             nodesol_ub;
   SCIP_Bool             nodesol_use;
} REDSOLLEVEL;


struct reduction_solution_storage
{
   STP_Vectype(REDSOLLEVEL*) levels;          /**< levels */
   SCIP_Real             offset;              /**< offset  */
   SCIP_Bool             isPcMw;
   SCIP_Bool             nodesol_use;
};



/*
 * Local methods
 */


/** updates node solution */
static
SCIP_RETCODE nodesolUpdate(
   SCIP*                scip,               /**< SCIP data structure */
   GRAPH*               g,
   SCIP_Real*           solval,             /**< FARAWAY if no valid solution build */
   PATH*                solpath,
   int*                 nodesol             /**< solution array to be filled */
   )
{
   assert(graph_typeIsSpgLike(g));

   /*
         if( pcmw )
         SCIP_CALL( SCIPStpHeurTMBuildTreePcMw(scip, prunegraph, path, prunegraph->cost, &objold, solnode) );
      else
         SCIPStpHeurTMBuildTree(scip, prunegraph, path, prunegraph->cost, &objold, solnode);

    */

   SCIPStpHeurTMBuildTree(scip, g, solpath, g->cost, solval, nodesol);

   return SCIP_OKAY;
}


/** gets edge solution after nodesolUpdate call */
static
void nodesolGetEdgeSol(
   const GRAPH*         graph,            /**< graph */
   const PATH*          solpath,          /**< stores solution */
   int*                 edgesol           /**< solution array to be filled */
   )
{
   const int nnodes = graph_get_nNodes(graph);
   const int nedges = graph_get_nEdges(graph);

   for( int e = 0; e < nedges; e++ )
      edgesol[e] = UNKNOWN;

   for( int k = 0; k < nnodes; k++ )
   {
      const int e = solpath[k].edge;

      if( e >= 0 )
         edgesol[e] = CONNECT;
   }
}


/** initializes */
static
SCIP_RETCODE redlevelInit(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nnodes,             /**< number of nodes */
   REDSOLLEVEL**         redlevel            /**< to be initialized */
   )
{
   REDSOLLEVEL* rl;

   assert(scip);
   assert(nnodes > 0);

   SCIP_CALL( SCIPallocMemory(scip, redlevel) );
   rl = *redlevel;

   rl->redsollocal = NULL;
   rl->solval_postred = REDSOLVAL_UNSET;
   rl->solval_incomplete = 0.0;
   rl->nnodes = nnodes;
   rl->nodesol_use = TRUE;
   rl->nodesol_ub = REDSOLVAL_UNSET;
   /* NOTE: will later be moved from sollocal*/
   rl->nodesol = NULL;

   return SCIP_OKAY;
}


/** frees */
static
void redlevelFree(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOLLEVEL**         redlevel            /**< to be freed */
   )
{
   REDSOLLEVEL* rl;

   assert(scip && redlevel);
   assert(*redlevel);

   rl = *redlevel;

   if( rl->redsollocal )
      reduce_sollocalFree(scip, &(rl->redsollocal));

   SCIPfreeMemoryArrayNull(scip, &(rl->nodesol));

   SCIPfreeMemory(scip, redlevel);
}


/** cleans */
static
void redlevelClean(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOLLEVEL*          redlevel            /**< to be cleaned */
   )
{
   assert(scip && redlevel);
   assert(!redlevel->redsollocal);

   redlevel->solval_postred = REDSOLVAL_UNSET;
   redlevel->solval_incomplete = 0.0;
}


/** initializes incomplete info */
static
SCIP_RETCODE redlevelInitIncomplete(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDSOL*         redsol,             /**< reduction solution */
   REDSOLLEVEL*          redlevel            /**< to be initialized */
   )
{
   assert(scip && redlevel);
   assert(redlevel->nnodes > 0);
   assert(EQ(redlevel->solval_incomplete, 0.0));

   redlevel->solval_incomplete = redsol->offset;

   SCIPdebugMessage("setting initial solval_incomplete=%f  \n", redsol->offset);

   return SCIP_OKAY;
}


/** adds local solution */
static
SCIP_RETCODE redlevelAddLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   REDSOLLEVEL*          redlevel,           /**< level */
   REDSOLLOCAL**         redsollocal_out     /**< pointer to newly initialized local */
   )
{
   assert(scip && redlevel);
   assert(!redlevel->redsollocal);
   assert(redlevel->nnodes > 0);

   SCIP_CALL( reduce_sollocalInit(scip, g, &(redlevel->redsollocal)) );

   if( redlevel->nodesol_use )
   {
      SCIP_CALL( reduce_sollocalInitNodesol(scip, redlevel->redsollocal) );
   }

   if( redsollocal_out )
   {
      *redsollocal_out = redlevel->redsollocal;
   }

   return SCIP_OKAY;
}


/** gets number of levels */
static
int redsolGetNlevels(
   const REDSOL*         redsol             /**< solution */
   )
{
   int levels;
   assert(redsol);

   levels = StpVecGetSize(redsol->levels);
   assert(levels > 0);

   return (levels);
}


/** gets top level */
static
REDSOLLEVEL* redsolGetTopLevel(
   const REDSOL*         redsol             /**< solution */
   )
{
   const int nlevels = redsolGetNlevels(redsol);

   return (redsol->levels[nlevels - 1]);
}





/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE reduce_sollocalInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   REDSOLLOCAL**         sollocal            /**< to initialize */
   )
{
   REDSOLLOCAL* rp;
   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, sollocal) );
   rp = *sollocal;

   rp->offset = 0.0;
   rp->primalbound = FARAWAY;
   rp->nnodes = graph_get_nNodes(g);
   rp->nodesol = NULL;
   rp->nodesol_use = FALSE;
   rp->nodesol_ub = REDSOLVAL_UNSET;
   rp->isPcMw = graph_pc_isPcMw(g);

   return SCIP_OKAY;
}


/** frees */
void reduce_sollocalFree(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOLLOCAL**         sollocal            /**< to free */
   )
{
   assert(scip && sollocal);
   assert(*sollocal);

   SCIPfreeMemoryArrayNull(scip, &((*sollocal)->nodesol));
   SCIPfreeMemory(scip, sollocal);
}


/** initializes node solution */
SCIP_RETCODE reduce_sollocalInitNodesol(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOLLOCAL*          sollocal            /**< solution */
   )
{
   int nnodes;
   assert(sollocal);
   assert(!sollocal->nodesol);
   assert(sollocal->nnodes > 0);
   assert(!sollocal->nodesol_use);

   nnodes = sollocal->nnodes;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sollocal->nodesol), nnodes) );
   sollocal->nodesol_use = TRUE;
   sollocal->nodesol_ub = REDSOLVAL_UNSET;

   for( int i = 0; i < nnodes; i++ )
   {
      sollocal->nodesol[i] = UNKNOWN;
   }

   return SCIP_OKAY;
}


/** sets offset
 *  NOTE: offset is only used within this structure! */
void reduce_sollocalSetOffset(
   SCIP_Real             offsetnew,          /**< new offset */
   REDSOLLOCAL*          sollocal            /**< sollocal */
   )
{
   assert(sollocal);
   assert(GE(offsetnew, 0.0));

   if( sollocal->isPcMw && !EQ(sollocal->offset, offsetnew) )
   {
      sollocal->primalbound = FARAWAY;
   }

   sollocal->offset = offsetnew;
}


/** updates */
SCIP_Bool reduce_sollocalUsesNodesol(
   const REDSOLLOCAL*          sollocal      /**< sollocal */
   )
{
   assert(sollocal);

   return sollocal->nodesol_use;
}


/** updates */
SCIP_RETCODE reduce_sollocalUpdateNodesol(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgesol,            /**< incoming solution */
   GRAPH*                g,                  /**< graph data structure */
   REDSOLLOCAL*          sollocal            /**< sollocal */
   )
{
   const int nnodes = graph_get_nNodes(g);

   assert(scip && edgesol && sollocal);
   assert(reduce_sollocalUsesNodesol(sollocal));
   assert(nnodes == sollocal->nnodes);
   assert(solstp_isValid(scip, g, edgesol));


 //  printf("updating nodesol... \n");

   {
      PATH* solpath;
      const SCIP_Real solval_new = solstp_getObj(g, edgesol, 0.0);
      SCIP_Real solval_old = FARAWAY;

      if( !EQ(sollocal->nodesol_ub, REDSOLVAL_UNSET) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &solpath, nnodes) );
         SCIP_CALL( nodesolUpdate(scip, g, &solval_old, solpath, sollocal->nodesol) );
         SCIPfreeBufferArray(scip, &solpath);
      }

   //   printf("old vs new: %f vs %f  \n", solval_old, solval_new);

      if( LT(solval_new, solval_old) )
      {
       //  printf("updating! \n");

         solstp_setVertexFromEdgeConn(g, edgesol, sollocal->nodesol);
         sollocal->nodesol_ub = solval_new;
      }
   }

   return SCIP_OKAY;
}


/** sets new sollocal bound if better than old one */
void reduce_sollocalUpdateUpperBound(
   SCIP_Real             ubnew,              /**< new upper bound, not including offset! */
   REDSOLLOCAL*          sollocal            /**< sollocal */
   )
{
   SCIP_Real ubnew_scaled;

   assert(sollocal);

   ubnew_scaled = ubnew + sollocal->offset;

   assert(GE(ubnew_scaled, 0.0));
   assert(LE(ubnew_scaled, FARAWAY));

   if( ubnew_scaled < sollocal->primalbound )
   {
      SCIPdebugMessage("updating upper bound %f->%f \n", sollocal->primalbound, ubnew_scaled);
      sollocal->primalbound = ubnew_scaled;
   }
}


/** gets upper bound; not including (last set) offset */
SCIP_Real reduce_sollocalGetUpperBound(
   const REDSOLLOCAL*    sollocal            /**< sollocal */
   )
{
   assert(sollocal);

   if( EQ(sollocal->primalbound, FARAWAY) )
      return FARAWAY;

   assert(GE(sollocal->primalbound - sollocal->offset, 0.0));

   SCIPdebugMessage("returning best bound: %f (%f-%f) \n", sollocal->primalbound - sollocal->offset, sollocal->primalbound, sollocal->offset);
   return (sollocal->primalbound - sollocal->offset);
}



/** gets array */
int* reduce_sollocalGetSolnode(
   REDSOLLOCAL*         sollocal            /**< sollocal */
   )
{
   assert(sollocal);
   assert(reduce_sollocalUsesNodesol(sollocal));

   return sollocal->nodesol;
}



/** gets upper bound; including (last set) offset */
SCIP_Real reduce_sollocalGetUpperBoundWithOffset(
   const REDSOLLOCAL*    sollocal            /**< sollocal */
   )
{
   assert(sollocal);

   return (sollocal->primalbound);
}


/** do we have a (non-trivial) primal bound? */
SCIP_Bool reduce_sollocalHasUpperBound(
   const REDSOLLOCAL*    sollocal            /**< sollocal */
   )
{
   assert(sollocal);
   assert(LE(sollocal->primalbound, FARAWAY));

   return !EQ(sollocal->primalbound, FARAWAY);
}


/** initializes */
SCIP_RETCODE reduce_solInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             useNodeSol,         /**< should solution be used? (additionally to solution value) */
   REDSOL**              redsol              /**< to be initialized */
   )
{
   REDSOLLEVEL* redlevel;
   REDSOL* rs;
   assert(scip && g);

   SCIP_CALL( SCIPallocMemory(scip, redsol) );
   rs = *redsol;

   rs->offset = 0.0;
   rs->isPcMw = graph_pc_isPcMw(g);
   rs->levels = NULL;
   rs->nodesol_use = useNodeSol;

   SCIP_CALL( redlevelInit(scip, g->knots, &redlevel) );
   StpVecPushBack(scip, rs->levels, redlevel);
   redlevel->nodesol_use = useNodeSol;

   assert(redsolGetNlevels(rs) == 1);

   return SCIP_OKAY;
}


/** adds local for given level; call before reduction loop */
SCIP_RETCODE reduce_solInitLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   REDSOL*               redsol,             /**< reduction solution */
   REDSOLLOCAL**         redsollocal_out     /**< pointer to newly initialized local */
   )
{
   const int nlevels = redsolGetNlevels(redsol);
   assert(redsol);
   assert(nlevels >= 1);

   SCIP_CALL( redlevelAddLocal(scip, g, redsol->levels[nlevels - 1], redsollocal_out) );

   return SCIP_OKAY;
}


/** finalizes local for given level; call after reduction loop */
void reduce_solFinalizeLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   REDSOL*               redsol              /**< reduction solution */
   )
{
   REDSOLLEVEL* toplevel = redsolGetTopLevel(redsol);
   REDSOLLOCAL* redsollocal = toplevel->redsollocal;
   const SCIP_Real offset = reduce_solGetOffset(redsol);

   assert(redsollocal);

   if( g->terms == 1 )
   {
      reduce_sollocalSetOffset(offset, redsollocal);
      reduce_sollocalUpdateUpperBound(0.0, redsollocal);
      // todo also set trivial node solution
   }

   assert(EQ(toplevel->solval_postred, REDSOLVAL_UNSET));

   if( reduce_sollocalHasUpperBound(redsollocal) )
      toplevel->solval_postred = reduce_sollocalGetUpperBoundWithOffset(redsollocal);
   else
      toplevel->solval_postred = FARAWAY;

   assert(!toplevel->nodesol);

   if( redsol->nodesol_use )
   {
      toplevel->nodesol_ub = redsollocal->nodesol_ub;
      toplevel->nodesol = redsollocal->nodesol;
      redsollocal->nodesol = NULL;
   }

   reduce_sollocalFree(scip, &(toplevel->redsollocal));
}


/** frees */
void reduce_solFree(
   SCIP*                 scip,               /**< SCIP data structure */
   REDSOL**              redsol              /**< to be freed */
   )
{
   REDSOL* rs;

   assert(scip && redsol);
   assert(*redsol);

   rs = *redsol;
   assert(!StpVecIsEmpty(rs->levels));

   while( !StpVecIsEmpty(rs->levels) )
   {
      const int size = StpVecGetSize(rs->levels);
      assert(size >= 1);

      redlevelFree(scip, &(rs->levels[size - 1]));
      StpVecPopBack(rs->levels);
   }

   StpVecFree(scip, rs->levels);
   SCIPfreeMemory(scip, redsol);
}


/** packs solution */
void reduce_solPack(
   const GRAPH*          g,                /**< graph data structure */
   const int*            nodes_old2packed, /**< map */
   int                   nnodes_packed,
   REDSOL*               redsol            /**< sollocal */
   )
{
   REDSOLLEVEL* redlevel;
   int* nodesol;
   const int nnodes = graph_get_nNodes(g);

   assert(g && redsol && nodes_old2packed);

   redlevel = redsol->levels[0];
   assert(redlevel);
   assert(nnodes == redlevel->nnodes);
   assert(1 <= nnodes_packed && nnodes_packed <= nnodes);

   redlevel->nnodes = nnodes_packed;

   if( !redsol->nodesol_use )
      return;

   nodesol = redlevel->nodesol;
   assert(nodesol);

   for( int i = 0; i < nnodes; i++ )
   {
      const int node_packed = nodes_old2packed[i];

      if( node_packed >= 0 )
      {
         assert(node_packed <= i);
         assert(node_packed < nnodes_packed);

         nodesol[node_packed] = nodesol[i] ;
      }
   }
}


/** adds level */
SCIP_RETCODE reduce_solLevelAdd(
   SCIP*                 scip,             /**< SCIP data structure */
   const GRAPH*          g,                /**< graph data structure */
   REDSOL*               redsol            /**< sollocal */
   )
{
   REDSOLLEVEL* redlevel_new;
   REDSOLLEVEL* redlevel_top = redsolGetTopLevel(redsol);
   int todo;

   redsol->nodesol_use = FALSE; // todo remove
   //printf("SETTING SOL FALSE (nodes=%d) \n", g->knots);


   SCIPdebugMessage("adding level %d \n", redsolGetNlevels(redsol));

   assert(redsol);
   assert(EQ(redsolGetTopLevel(redsol)->solval_incomplete, 0.0));

   SCIP_CALL( redlevelInitIncomplete(scip, redsol, redlevel_top) );

   SCIP_CALL( redlevelInit(scip, g->knots, &redlevel_new) );
   StpVecPushBack(scip, redsol->levels, redlevel_new);
   redlevel_new->nodesol_use = redsol->nodesol_use;

   return SCIP_OKAY;
}


/** removes top level */
void reduce_solLevelTopRemove(
   SCIP*                 scip,             /**< SCIP data structure */
   REDSOL*               redsol            /**< sollocal */
   )
{
   const int nlevels = redsolGetNlevels(redsol);
   assert(nlevels >= 2);

   redlevelFree(scip, &(redsol->levels[nlevels - 1]));
   StpVecPopBack(redsol->levels);

   assert(redsolGetNlevels(redsol) == nlevels - 1);
}


/** finalizes top level; also removes the level! */
void reduce_solLevelTopFinalize(
   SCIP*                 scip,             /**< SCIP data structure */
   REDSOL*               redsol            /**< sollocal */
   )
{
   REDSOLLOCAL* redsollocal_top;
   REDSOLLEVEL* redlevel_top;
   const int nlevels = redsolGetNlevels(redsol);
   assert(nlevels >= 2);

   redlevelFree(scip, &(redsol->levels[nlevels - 1]));
   StpVecPopBack(redsol->levels);

   redlevel_top = redsolGetTopLevel(redsol);
   redsollocal_top = redlevel_top->redsollocal;
   assert(redsollocal_top);

   reduce_sollocalSetOffset(0.0, redsollocal_top);
   reduce_sollocalUpdateUpperBound(redlevel_top->solval_incomplete, redsollocal_top);

   SCIPdebugMessage("finalizing top level %d; solval=%f offset=%f \n", nlevels - 1, redlevel_top->solval_incomplete, redsol->offset);

   assert(redsolGetNlevels(redsol) == nlevels - 1);
}



/** removes level */
void reduce_solLevelTopClean(
   SCIP*                 scip,             /**< SCIP data structure */
   REDSOL*               redsol            /**< sollocal */
   )
{
   const int nlevels = redsolGetNlevels(redsol);
   assert(nlevels >= 1);

   redlevelClean(scip, redsol->levels[nlevels - 1]);
}


/** merges level up */
void reduce_solLevelTopMergeUp(
   REDSOL*               redsol            /**< sollocal */
   )
{
   REDSOLLEVEL* levelTop;
   REDSOLLEVEL* levelParent;
   const int nlevels = redsolGetNlevels(redsol);
   assert(nlevels >= 2);

   levelTop = redsol->levels[nlevels - 1];
   levelParent = redsol->levels[nlevels - 2];

   assert(!levelTop->redsollocal);
   assert(levelParent->redsollocal);
   assert(GE(levelTop->solval_postred, 0.0));

   SCIPdebugMessage("merge levels %d->%d \n", nlevels - 1, nlevels - 2);
   SCIPdebugMessage("...update solval_incomplete: %f->%f \n", levelParent->solval_incomplete, levelParent->solval_incomplete + levelTop->solval_postred);

   levelParent->solval_incomplete += levelTop->solval_postred;

   if( GT(levelParent->solval_incomplete, FARAWAY) )
      levelParent->solval_incomplete = FARAWAY;
}


/** sets offset */
void reduce_solSetOffset(
   SCIP_Real             offsetnew,          /**< new offset */
   REDSOL*               redsol              /**< solution */
   )
{
   assert(redsol);
   assert(GE(offsetnew, 0.0));

   redsol->offset = offsetnew;
}


/** gets */
SCIP_Real reduce_solGetOffset(
   const REDSOL*         redsol            /**< solution */
   )
{
   assert(redsol);


   return (redsol->offset);
}



/** adds (and possibly overwrites) nodesol */

#ifdef XXXXX
SCIP_RETCODE reduce_solAddNodesol(
   const int*           nodesol,          /**< incoming solution */
   REDSOL*              redsol            /**< solution */
   )
{
   assert(redsol);


   return SCIP_OKAY;
}
#endif


/** initializes */
SCIP_RETCODE reduce_solGetEdgesol(
   SCIP*                scip,             /**< SCIP data structure */
   GRAPH*               g,                /**< graph data structure */
   REDSOL*              redsol,           /**< solution */
   SCIP_Real*           solval,           /**< FARAWAY if no solution available */
   int*                 edgesol           /**< solution array to be filled */
   )
{
   REDSOLLEVEL* level;

   assert(scip && g && redsol && solval && edgesol);

   *solval = FARAWAY;
   level = redsol->levels[0];
   assert(level);

   if( !redsol->nodesol_use )
      return SCIP_OKAY;

   if( !EQ(level->nodesol_ub, REDSOLVAL_UNSET) && LT(level->nodesol_ub, FARAWAY) )
   {
      PATH* solpath;
      const int nnodes = graph_get_nNodes(g);

      SCIP_CALL( SCIPallocBufferArray(scip, &solpath, nnodes) );
      SCIP_CALL( nodesolUpdate(scip, g, solval, solpath, level->nodesol) );

      if( LT(*solval, FARAWAY) )
         nodesolGetEdgeSol(g, solpath, edgesol);

      SCIPfreeBufferArray(scip, &solpath);
   }

   return SCIP_OKAY;
}


/** is node solution in use? */
SCIP_Bool reduce_solUsesNodesol(
   const REDSOL*         redsol            /**< solution */
   )
{

   assert(redsol);

   return (redsol->nodesol_use);
}


/** gets */
SCIP_Real reduce_solGetUpperBoundWithOffset(
   const REDSOL*         redsol            /**< solution */
   )
{
   REDSOLLEVEL* level;
   SCIP_Real ub;

   assert(redsol);
   assert(redsolGetNlevels(redsol) >= 1);
   assert(redsol->levels);

   level = redsol->levels[0];
   ub = level->solval_postred;

   if( EQ(ub, REDSOLVAL_UNSET) )
      return FARAWAY;

   return (ub);
}


/** gets */
SCIP_Real* reduce_solGetOffsetPointer(
  REDSOL*                redsol            /**< solution */
   )
{
   assert(redsol);


   return &(redsol->offset);
}
