/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepastore.c,v 1.19 2004/07/29 09:01:56 bzfpfend Exp $"

/**@file   sepastore.c
 * @brief  methods for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "set.h"
#include "lp.h"
#include "var.h"
#include "tree.h"
#include "sepastore.h"

#include "struct_sepastore.h"



/*
 * dynamic memory arrays
 */

/** resizes cuts and score arrays to be able to store at least num entries */
static
RETCODE sepastoreEnsureCutsMem(
   SEPASTORE*       sepastore,          /**< separation storage */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(sepastore != NULL);
   assert(set != NULL);

   if( num > sepastore->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&sepastore->cuts, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->score, newsize) );
      sepastore->cutssize = newsize;
   }
   assert(num <= sepastore->cutssize);

   return SCIP_OKAY;
}

/** resizes bdchgvars, bdchgvals, and bdchgtypes arrays to be able to store at least num entries */
static
RETCODE sepastoreEnsureBdchgsMem(
   SEPASTORE*       sepastore,          /**< separation storage */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(sepastore != NULL);
   assert(set != NULL);

   if( num > sepastore->bdchgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&sepastore->bdchgvars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->bdchgvals, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->bdchgtypes, newsize) );
      sepastore->bdchgssize = newsize;
   }
   assert(num <= sepastore->bdchgssize);

   return SCIP_OKAY;
}




/** creates separation storage */
RETCODE SCIPsepastoreCreate(
   SEPASTORE**           sepastore                /**< pointer to store separation storage */
   )
{
   assert(sepastore != NULL);
   
   ALLOC_OKAY( allocMemory(sepastore) );
   
   (*sepastore)->cuts = NULL;
   (*sepastore)->score = NULL;
   (*sepastore)->bdchgvars = NULL;
   (*sepastore)->bdchgvals = NULL;
   (*sepastore)->bdchgtypes = NULL;
   (*sepastore)->cutssize = 0;
   (*sepastore)->bdchgssize = 0;
   (*sepastore)->ncuts = 0;
   (*sepastore)->nbdchgs = 0;
   (*sepastore)->ncutsfound = 0;
   (*sepastore)->ncutsapplied = 0;
   (*sepastore)->initiallp = FALSE;

   return SCIP_OKAY;
}

/** frees separation storage */
RETCODE SCIPsepastoreFree(
   SEPASTORE**           sepastore                /**< pointer to store separation storage */
   )
{
   assert(sepastore != NULL);

   freeMemoryArrayNull(&(*sepastore)->cuts);
   freeMemoryArrayNull(&(*sepastore)->score);
   freeMemoryArrayNull(&(*sepastore)->bdchgvars);
   freeMemoryArrayNull(&(*sepastore)->bdchgvals);
   freeMemoryArrayNull(&(*sepastore)->bdchgtypes);
   freeMemory(sepastore);

   return SCIP_OKAY;
}

/** informs separation storage, that the setup of the initial LP starts now */
void SCIPsepastoreStartInitialLP(
   SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(!sepastore->initiallp);
   assert(sepastore->ncuts == 0);

   sepastore->initiallp = TRUE;
}

/** informs separation storage, that the setup of the initial LP is now finished */
void SCIPsepastoreEndInitialLP(
   SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->initiallp);
   assert(sepastore->ncuts == 0);

   sepastore->initiallp = FALSE;
}

/** adds cut stored as LP row to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
static
RETCODE sepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             root                /**< are we at the root node? */
   )
{
   int maxsepacuts;
   int c;
   int i;

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(cut != NULL);
   assert(!SCIProwIsInLP(cut));
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));

   /**@todo don't add redundant cuts to the separation store */

   /* get maximum of separated cuts at this node:
    *  - for root LP, use all cuts
    *  - use all cuts that have infinite score
    *  - the remaining cuts are only forwarded to the LP, if they fit into the maximal separation size
    */
   if( sepastore->initiallp || SCIPsetIsInfinity(set, score) )
      maxsepacuts = INT_MAX;
   else
      maxsepacuts = SCIPsetGetMaxsepacuts(set, root);
   assert(sepastore->ncuts <= maxsepacuts);
   if( maxsepacuts == 0 )
      return SCIP_OKAY;

   /* get enough memory to store the cut */
   CHECK_OKAY( sepastoreEnsureCutsMem(sepastore, set, sepastore->ncuts+1) );
   assert(sepastore->ncuts < sepastore->cutssize);

   /* check, if cut belongs to the best "maxsepacuts" separation cuts */
   if( sepastore->ncuts < maxsepacuts || score > sepastore->score[maxsepacuts-1] )
   {
      debugMessage("adding cut to separation storage of size %d/%d\n", sepastore->ncuts, maxsepacuts);
      debug(SCIProwPrint(cut, NULL));

      /* capture the cut */
      SCIProwCapture(cut);

      /* search the correct position of the cut in the cuts array */
      for( c = 0; c < sepastore->ncuts && score <= sepastore->score[c]; ++c )
      {
      }
      assert(c <= sepastore->ncuts);
      assert(c < maxsepacuts);

      /* if the array consists of "maxsepacuts" cuts, release the worst cut */
      if( sepastore->ncuts == maxsepacuts )
      {
         CHECK_OKAY( SCIProwRelease(&sepastore->cuts[sepastore->ncuts-1], memhdr, set, lp) );
         sepastore->ncuts--;
      }
      assert(sepastore->ncuts < maxsepacuts);

      /* insert cut in the sorted arrays */
      for( i = sepastore->ncuts; i > c; --i )
      {
         sepastore->cuts[i] = sepastore->cuts[i-1];
         sepastore->score[i] = sepastore->score[i-1];
      }
      sepastore->cuts[c] = cut;
      sepastore->score[c] = score;
      sepastore->ncuts++;
   }

   return SCIP_OKAY;
}

/** adds cut stored as bound change to separation storage */
static
RETCODE sepastoreAddBdchg(
   SEPASTORE*       sepastore,          /**< separation storage */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new bound value */
   BOUNDTYPE        boundtype           /**< type of bound to change */
   )
{
   debugMessage("adding bound change to separation storage: variable <%s>, new %s bound: %g\n",
      SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound);

   /* get enough memory to store the cut */
   CHECK_OKAY( sepastoreEnsureBdchgsMem(sepastore, set, sepastore->nbdchgs+1) );
   assert(sepastore->nbdchgs < sepastore->bdchgssize);

   /* insert bound change in last slot of bdchgs arrays */
   sepastore->bdchgvars[sepastore->nbdchgs] = var;
   sepastore->bdchgvals[sepastore->nbdchgs] = newbound;
   sepastore->bdchgtypes[sepastore->nbdchgs] = boundtype;
   sepastore->nbdchgs++;

   return SCIP_OKAY;
}

/** adds cut to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
RETCODE SCIPsepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             root                /**< are we at the root node? */
   )
{
   assert(sepastore != NULL);
   assert(cut != NULL);
   assert(!SCIProwIsInLP(cut));
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));

   /* update statistics of total number of found cuts */
   if( !sepastore->initiallp )
      sepastore->ncutsfound++;

   /* check, if the cut is a bound change (i.e. a row with only one variable) */
   if( !SCIProwIsModifiable(cut) && SCIProwGetNNonz(cut) == 1 )
   {
      COL** cols;
      VAR* var;
      Real* vals;
      Real lhs;
      Real rhs;
      Real bound;

      /* get the single variable and its coefficient of the cut */
      cols = SCIProwGetCols(cut);
      assert(cols != NULL);
      var = SCIPcolGetVar(cols[0]);
      vals = SCIProwGetVals(cut);
      assert(vals != NULL);
      assert(!SCIPsetIsZero(set, vals[0]));

      /* get the left hand side of the cut and convert it to a bound */
      lhs = SCIProwGetLhs(cut);
      if( !SCIPsetIsInfinity(set, -lhs) )
      {
         lhs -= SCIProwGetConstant(cut);
         if( vals[0] > 0.0 )
         {
            /* coefficient is positive -> lhs corresponds to lower bound */
            bound = lhs/vals[0];
            if( SCIPsetIsGT(set, bound, SCIPvarGetLbLocal(var)) )
            {
               CHECK_OKAY( sepastoreAddBdchg(sepastore, set, var, bound, SCIP_BOUNDTYPE_LOWER) );
            }
         }
         else
         {
            /* coefficient is negative -> lhs corresponds to upper bound */
            bound = lhs/vals[0];
            if( SCIPsetIsLT(set, bound, SCIPvarGetUbLocal(var)) )
            {
               CHECK_OKAY( sepastoreAddBdchg(sepastore, set, var, bound, SCIP_BOUNDTYPE_UPPER) );
            }
         }
      }

      /* get the right hand side of the cut and convert it to a bound */
      rhs = SCIProwGetRhs(cut);
      if( !SCIPsetIsInfinity(set, rhs) )
      {
         rhs -= SCIProwGetConstant(cut);
         if( vals[0] > 0.0 )
         {
            /* coefficient is positive -> rhs corresponds to upper bound */
            bound = rhs/vals[0];
            if( SCIPsetIsLT(set, bound, SCIPvarGetUbLocal(var)) )
            {
               CHECK_OKAY( sepastoreAddBdchg(sepastore, set, var, bound, SCIP_BOUNDTYPE_UPPER) );
            }
         }
         else
         {
            /* coefficient is negative -> rhs corresponds to lower bound */
            bound = rhs/vals[0];
            if( SCIPsetIsGT(set, bound, SCIPvarGetLbLocal(var)) )
            {
               CHECK_OKAY( sepastoreAddBdchg(sepastore, set, var, bound, SCIP_BOUNDTYPE_LOWER) );
            }
         }
      }
   }
   else
   {
      /* add LP row cut to separation storage */
      CHECK_OKAY( sepastoreAddCut(sepastore, memhdr, set, lp, cut, score, root) );
   }

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
RETCODE SCIPsepastoreApplyCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   VAR* var;
   Real val;
   int i;

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);

   /**@todo avoid applying redundant cuts (see sepastoreAddCut()) */
   debugMessage("applying %d bound changes and %d cuts\n", sepastore->nbdchgs, sepastore->ncuts);

   /* apply cuts stored as bound changes */
   *cutoff = FALSE;
   for( i = 0; i < sepastore->nbdchgs && !(*cutoff); ++i )
   {
      var = sepastore->bdchgvars[i];
      val = sepastore->bdchgvals[i];
      if( sepastore->bdchgtypes[i] == SCIP_BOUNDTYPE_LOWER )
      {
         if( SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var)) )
         {
            debugMessage("apply bound change: <%s>: [%g,%g] -> [%g,%g]\n", 
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), val, SCIPvarGetUbLocal(var));

            if( SCIPsetIsLE(set, val, SCIPvarGetUbLocal(var)) )
            {
               CHECK_OKAY( SCIPnodeAddBoundchg(tree->actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue,
                              var, val, SCIP_BOUNDTYPE_LOWER) );
            }
            else
               *cutoff = TRUE;

            if( !sepastore->initiallp )
               sepastore->ncutsapplied++;
         }
      }
      else
      {
         if( SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var)) )
         {
            debugMessage("apply bound change: <%s>: [%g,%g] -> [%g,%g]\n", 
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var), val);

            if( SCIPsetIsGE(set, val, SCIPvarGetLbLocal(var)) )
            {
               CHECK_OKAY( SCIPnodeAddBoundchg(tree->actnode, memhdr, set, stat, tree, lp, branchcand, eventqueue,
                              var, val, SCIP_BOUNDTYPE_UPPER) );
            }
            else
               *cutoff = TRUE;

            if( !sepastore->initiallp )
               sepastore->ncutsapplied++;
         }
      }
   }

   /* apply cuts stored as LP rows */
   for( i = 0; i < sepastore->ncuts; ++i )
   {
      /* a row could have been added twice to the separation store; add it only once! */
      if( !SCIProwIsInLP(sepastore->cuts[i]) )
      {
         debugMessage("apply cut: ");
         debug( SCIProwPrint(sepastore->cuts[i], NULL) );
         
         /* add cut to the LP and capture it */
         CHECK_OKAY( SCIPlpAddRow(lp, set, sepastore->cuts[i]) );
         
         /* release the row */
         CHECK_OKAY( SCIProwRelease(&sepastore->cuts[i], memhdr, set, lp) );
         
         if( !sepastore->initiallp )
            sepastore->ncutsapplied++;
      }
   }

   /* clear the separation storage */
   sepastore->ncuts = 0;
   sepastore->nbdchgs = 0;

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
RETCODE SCIPsepastoreClearCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepastore != NULL);

   debugMessage("clearing %d cuts\n", sepastore->ncuts);

   for( c = 0; c < sepastore->ncuts; ++c )
   {
      /* release the row */
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[c], memhdr, set, lp) );
   }

   /* clear the separation storage */
   sepastore->ncuts = 0;

   return SCIP_OKAY;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreGetNCuts(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncuts;
}

/** get total number of cuts found so far */
int SCIPsepastoreGetNCutsFound(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfound;
}

/** get total number of cuts applied to the LPs */
int SCIPsepastoreGetNCutsApplied(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsapplied;
}
