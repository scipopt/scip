/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepastore.c,v 1.37 2005/05/03 14:48:03 bzfpfend Exp $"

/**@file   sepastore.c
 * @brief  methods for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/sepastore.h"
#include "scip/debug.h"

#include "scip/struct_sepastore.h"



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
      ALLOC_OKAY( reallocMemoryArray(&sepastore->efficacies, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->objparallelisms, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->orthogonalities, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&sepastore->scores, newsize) );
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
   (*sepastore)->efficacies = NULL;
   (*sepastore)->objparallelisms = NULL;
   (*sepastore)->orthogonalities = NULL;
   (*sepastore)->scores = NULL;
   (*sepastore)->bdchgvars = NULL;
   (*sepastore)->bdchgvals = NULL;
   (*sepastore)->bdchgtypes = NULL;
   (*sepastore)->cutssize = 0;
   (*sepastore)->bdchgssize = 0;
   (*sepastore)->ncuts = 0;
   (*sepastore)->nbdchgs = 0;
   (*sepastore)->ncutsfound = 0;
   (*sepastore)->ncutsfoundround = 0;
   (*sepastore)->ncutsstored = 0;
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
   assert(*sepastore != NULL);
   assert((*sepastore)->ncuts == 0);
   assert((*sepastore)->nbdchgs == 0);

   freeMemoryArrayNull(&(*sepastore)->cuts);
   freeMemoryArrayNull(&(*sepastore)->efficacies);
   freeMemoryArrayNull(&(*sepastore)->objparallelisms);
   freeMemoryArrayNull(&(*sepastore)->orthogonalities);
   freeMemoryArrayNull(&(*sepastore)->scores);
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
   assert(sepastore->nbdchgs == 0);

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
   assert(sepastore->nbdchgs == 0);

   sepastore->initiallp = FALSE;
}

/** checks cut for redundancy due to activity bounds */
static
Bool sepastoreIsCutRedundant(
   SEPASTORE*       sepastore,          /**< separation storage */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   ROW*             cut                 /**< separated cut */
   )
{
   Real minactivity;
   Real maxactivity;
   Real lhs;
   Real rhs;

   assert(sepastore != NULL);
   assert(cut != NULL);

   /* modifiable cuts cannot be declared redundant, since we don't know all coefficients */
   if( SCIProwIsModifiable(cut) )
      return FALSE;

   /* check for activity redundancy */
   lhs = SCIProwGetLhs(cut);
   rhs = SCIProwGetRhs(cut);
   minactivity = SCIProwGetMinActivity(cut, set, stat);
   maxactivity = SCIProwGetMaxActivity(cut, set, stat);
   if( SCIPsetIsLE(set, lhs, minactivity) && SCIPsetIsLE(set, maxactivity, rhs) )
   {
      debugMessage("ignoring activity redundant cut <%s> (sides=[%g,%g], act=[%g,%g]\n",
         SCIProwGetName(cut), lhs, rhs, minactivity, maxactivity);
      debug(SCIProwPrint(cut, NULL));
      return TRUE;
   }

   return FALSE;
}

/** adds cut stored as LP row to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
static
RETCODE sepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   Bool             root                /**< are we at the root node? */
   )
{
   Real mincutorthogonality;
   Real cutefficacy;
   Real cutobjparallelism;
   Real cutorthogonality;
   Real cutscore;
   int maxsepacuts;
   int c;
   int i;
   int j;

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(cut != NULL);
   assert(!SCIProwIsInLP(cut));
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));

   /* check cut for redundancy
    * in each separation round, make sure that at least one (even redundant) cut enters the LP to avoid cycling
    */
   if( !forcecut && sepastore->ncuts + sepastore->nbdchgs > 0 && sepastoreIsCutRedundant(sepastore, set, stat, cut) )
      return SCIP_OKAY;

   /* if only one cut is currently present in the cut store, it could be redundant; in this case, it can now be removed
    * again, because now a non redundant cut enters the store
    */
   if( sepastore->ncuts == 1 && sepastore->nbdchgs == 0
      && sepastoreIsCutRedundant(sepastore, set, stat, sepastore->cuts[0]) )
   {
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[0], blkmem, set, lp) );
      sepastore->ncuts--;
   }

   /* get maximum of separated cuts at this node:
    *  - for initial LP, use all cuts
    *  - use all cuts that have infinite score factor
    *  - the remaining cuts are only forwarded to the LP, if they fit into the maximal separation size
    */
   forcecut = forcecut || sepastore->initiallp;
   if( forcecut )
      maxsepacuts = INT_MAX;
   else
      maxsepacuts = SCIPsetGetSepaMaxcuts(set, root);
   assert(sepastore->ncuts <= maxsepacuts);
   if( maxsepacuts == 0 )
      return SCIP_OKAY;

   /* calculate cut's efficacy and objective hyperplane parallelism;
    * orthogonality starts with 1.0 and is reduced when the cut is inserted
    */
   cutobjparallelism = SCIProwGetObjParallelism(cut, set, lp);
   cutorthogonality = 1.0;
   if( maxsepacuts == INT_MAX )
   {
      cutefficacy = SCIPsetInfinity(set);
      cutscore = SCIPsetInfinity(set);
   }
   else
   {
      cutefficacy = SCIProwGetEfficacy(cut, set, stat, lp);
      cutscore = cutefficacy + set->sepa_objparalfac * cutobjparallelism + set->sepa_orthofac * cutorthogonality;
   }

   /* check, if cut has potential to belong to the best "maxsepacuts" separation cuts */
   if( sepastore->ncuts >= maxsepacuts && cutscore <= sepastore->scores[maxsepacuts-1] )
      return SCIP_OKAY;

   /* calculate minimal cut orthogonality */
   mincutorthogonality = (root ? set->sepa_minorthoroot : set->sepa_minortho);

   /* get enough memory to store the cut */
   CHECK_OKAY( sepastoreEnsureCutsMem(sepastore, set, sepastore->ncuts+1) );
   assert(sepastore->ncuts < sepastore->cutssize);

   debugMessage("adding cut to separation storage of size %d/%d (forcecut=%d, efficacy=%g, objparallelism=%g, score=%g)\n", 
      sepastore->ncuts, maxsepacuts, forcecut, cutefficacy, cutobjparallelism, cutscore);
   debug(SCIProwPrint(cut, NULL));

   /* search the correct position of the cut in the cuts array */
   for( c = 0; c < sepastore->ncuts && cutscore <= sepastore->scores[c]; ++c )
   {
      /* update the minimal orthogonality of the cut and it's score */
      if( set->sepa_orthofac > 0.0 && !forcecut )
      {
         Real thisortho;
         
         assert(!SCIPsetIsInfinity(set, cutscore));
         thisortho = SCIProwGetOrthogonality(cut, sepastore->cuts[c]);
         if( thisortho < cutorthogonality )
         {
            cutorthogonality = thisortho;
            cutscore = cutefficacy + set->sepa_objparalfac * cutobjparallelism + set->sepa_orthofac * cutorthogonality;
            
            /* check, if the cut (after regarding orthogonality) is still good enough */
            if( (sepastore->ncuts >= maxsepacuts && cutscore <= sepastore->scores[maxsepacuts-1])
               || cutorthogonality < mincutorthogonality )
               return SCIP_OKAY;
         }
      }
   }
   assert(c <= sepastore->ncuts);
   assert(c < maxsepacuts);

   /* the cut is good enough to keep at this point: increase stored counter */
   sepastore->ncutsstored++;

   /* insert cut in the sorted arrays, and update the orthogonalities and scores of the inferior cuts;
    * the inferior cuts must be resorted, because their score could have been changed;
    *
    * Here, we make an error:
    * If a cut is moved down due to a reduced score, the cut that moved up has now a potentially better
    * orthogonality value because the orthogonality to the cut that was moved down leaves the minimum function.
    * However, we do not account for this, because otherwise we would have to recalculate the minimum, which is
    * too costly.
    */
   for( i = sepastore->ncuts; i > c; --i )
   {
      ROW* currentcut;
      Real currentefficacy;
      Real currentobjparallelism;
      Real currentorthogonality;
      Real currentscore;

      /* slot i is empty; slot i-1 has to fall down to a slot in {i,...,ncuts} */
      currentcut = sepastore->cuts[i-1];
      currentefficacy = sepastore->efficacies[i-1];
      currentobjparallelism = sepastore->objparallelisms[i-1];
      currentorthogonality = sepastore->orthogonalities[i-1];
      currentscore = sepastore->scores[i-1];
      assert(!SCIPsetIsInfinity(set, currentscore));

      /* update orthogonality and score of cut at slot i-1 due to new cut */
      if( set->sepa_orthofac > 0.0 )
      {
         Real thisortho;
               
         thisortho = SCIProwGetOrthogonality(currentcut, cut);
         if( thisortho < currentorthogonality )
         {
            currentorthogonality = thisortho;
            if( currentorthogonality < mincutorthogonality )
               currentscore = -SCIPsetInfinity(set);
            else
               currentscore = currentefficacy + set->sepa_objparalfac * currentobjparallelism
                  + set->sepa_orthofac * currentorthogonality;
         }
      }

      /* insert the current cut in the list of inferior cuts */
      for( j = i; j < sepastore->ncuts && sepastore->scores[j+1] > currentscore; ++j )
      {
         sepastore->cuts[j] = sepastore->cuts[j+1];
         sepastore->efficacies[j] = sepastore->efficacies[j+1];
         sepastore->objparallelisms[j] = sepastore->objparallelisms[j+1];
         sepastore->orthogonalities[j] = sepastore->orthogonalities[j+1];
         sepastore->scores[j] = sepastore->scores[j+1];

         /* update orthogonality and score of cut at slot i-1 */
         if( set->sepa_orthofac > 0.0 && !SCIPsetIsInfinity(set, -currentscore) )
         {
            Real thisortho;
                  
            thisortho = SCIProwGetOrthogonality(currentcut, sepastore->cuts[j]);
            if( thisortho < currentorthogonality )
            {
               currentorthogonality = thisortho;
               currentscore = currentefficacy + set->sepa_objparalfac * currentobjparallelism
                  + set->sepa_orthofac * currentorthogonality;
            }
         }
      }
      sepastore->cuts[j] = currentcut;
      sepastore->efficacies[j] = currentefficacy;
      sepastore->objparallelisms[j] = currentobjparallelism;
      sepastore->orthogonalities[j] = currentorthogonality;
      sepastore->scores[j] = currentscore;
   }
   sepastore->cuts[c] = cut;
   sepastore->efficacies[c] = cutefficacy;
   sepastore->objparallelisms[c] = cutobjparallelism;
   sepastore->orthogonalities[c] = cutorthogonality;
   sepastore->scores[c] = cutscore;
   sepastore->ncuts++;
   assert(c == 0 || sepastore->scores[c-1] >= sepastore->scores[c]);
   assert(c == sepastore->ncuts-1 || sepastore->scores[c] >= sepastore->scores[c+1]);

   /* capture the cut */
   SCIProwCapture(cut);

   /* if the array consists of more than "maxsepacuts" cuts, release the worst cut */
   if( sepastore->ncuts > maxsepacuts )
   {
      assert(sepastore->ncuts == maxsepacuts+1);
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[sepastore->ncuts-1], blkmem, set, lp) );
      sepastore->ncuts--;
   }
   assert(sepastore->ncuts <= maxsepacuts);

   /* remove cuts that are no longer orthogonal enough */
   while( sepastore->ncuts > 0 && sepastore->orthogonalities[sepastore->ncuts-1] < mincutorthogonality )
   {
      assert(SCIPsetIsInfinity(set, -sepastore->scores[sepastore->ncuts-1]));
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[sepastore->ncuts-1], blkmem, set, lp) );
      sepastore->ncuts--;
   }

   debugMessage(" -> cut was added at position %d/%d (efficacy=%g, objparallelism=%g, orthogonality=%g, score=%g)\n",
      c, sepastore->ncuts, cutefficacy, cutobjparallelism, cutorthogonality, cutscore);

   return SCIP_OKAY;
}

/** adds cut stored as bound change to separation storage */
static
RETCODE sepastoreAddBdchg(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< LP data */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new bound value */
   BOUNDTYPE        boundtype           /**< type of bound to change */
   )
{
   debugMessage("adding bound change to separation storage: variable <%s>, new %s bound: %g\n",
      SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", newbound);

   /* if only one cut is currently present in the cut store, it could be redundant; in this case, it can now be removed
    * again, because now a non redundant cut enters the store
    */
   if( sepastore->ncuts == 1 && sepastore->nbdchgs == 0
      && sepastoreIsCutRedundant(sepastore, set, stat, sepastore->cuts[0]) )
   {
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[0], blkmem, set, lp) );
      sepastore->ncuts--;
   }

   /* get enough memory to store the cut */
   CHECK_OKAY( sepastoreEnsureBdchgsMem(sepastore, set, sepastore->nbdchgs+1) );
   assert(sepastore->nbdchgs < sepastore->bdchgssize);

   /* insert bound change in last slot of bdchgs arrays */
   sepastore->bdchgvars[sepastore->nbdchgs] = var;
   sepastore->bdchgvals[sepastore->nbdchgs] = newbound;
   sepastore->bdchgtypes[sepastore->nbdchgs] = boundtype;
   sepastore->nbdchgs++;

   /* count the bound change as stored cut (bound changes are always good enough to keep) */
   sepastore->ncutsstored++;

   return SCIP_OKAY;
}

/** adds cut to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
RETCODE SCIPsepastoreAddCut(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   Bool             root                /**< are we at the root node? */
   )
{
   assert(sepastore != NULL);
   assert(cut != NULL);
   assert(!SCIProwIsInLP(cut));
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));

   /* debug: check cut for feasibility */
   CHECK_OKAY( SCIPdebugCheckRow(cut, set) ); /*lint !e506 !e774*/

   /* update statistics of total number of found cuts */
   if( !sepastore->initiallp )
   {
      sepastore->ncutsfound++;
      sepastore->ncutsfoundround++;
   }

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
               CHECK_OKAY( sepastoreAddBdchg(sepastore, blkmem, set, stat, lp, var, bound, SCIP_BOUNDTYPE_LOWER) );
            }
         }
         else
         {
            /* coefficient is negative -> lhs corresponds to upper bound */
            bound = lhs/vals[0];
            if( SCIPsetIsLT(set, bound, SCIPvarGetUbLocal(var)) )
            {
               CHECK_OKAY( sepastoreAddBdchg(sepastore, blkmem, set, stat, lp, var, bound, SCIP_BOUNDTYPE_UPPER) );
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
               CHECK_OKAY( sepastoreAddBdchg(sepastore, blkmem, set, stat, lp, var, bound, SCIP_BOUNDTYPE_UPPER) );
            }
         }
         else
         {
            /* coefficient is negative -> rhs corresponds to lower bound */
            bound = rhs/vals[0];
            if( SCIPsetIsGT(set, bound, SCIPvarGetLbLocal(var)) )
            {
               CHECK_OKAY( sepastoreAddBdchg(sepastore, blkmem, set, stat, lp, var, bound, SCIP_BOUNDTYPE_LOWER) );
            }
         }
      }
   }
   else
   {
      /* add LP row cut to separation storage */
      CHECK_OKAY( sepastoreAddCut(sepastore, blkmem, set, stat, lp, cut, forcecut, root) );
   }

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
RETCODE SCIPsepastoreApplyCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   NODE* node;
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

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL);

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
            debugMessage(" -> applying bound change: <%s>: [%g,%g] -> [%g,%g]\n", 
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), val, SCIPvarGetUbLocal(var));

            if( SCIPsetIsLE(set, val, SCIPvarGetUbLocal(var)) )
            {
               CHECK_OKAY( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue,
                     var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );
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
            debugMessage(" -> applying bound change: <%s>: [%g,%g] -> [%g,%g]\n", 
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var), val);

            if( SCIPsetIsGE(set, val, SCIPvarGetLbLocal(var)) )
            {
               CHECK_OKAY( SCIPnodeAddBoundchg(node, blkmem, set, stat, tree, lp, branchcand, eventqueue,
                     var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );
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
         debugMessage(" -> applying cut %d/%d (efficacy=%g, objparallelism=%g, orthogonality=%g, score=%g)\n",
            i, sepastore->ncuts, sepastore->efficacies[i], sepastore->objparallelisms[i], 
            sepastore->orthogonalities[i], sepastore->scores[i]);
         /*debug(SCIProwPrint(sepastore->cuts[i], NULL));*/
         assert(i == 0 || sepastore->scores[i] <= sepastore->scores[i-1]);

         /* add cut to the LP and capture it */
         CHECK_OKAY( SCIPlpAddRow(lp, set, sepastore->cuts[i], SCIPnodeGetDepth(node)) );
         
         /* release the row */
         CHECK_OKAY( SCIProwRelease(&sepastore->cuts[i], blkmem, set, lp) );
         
         if( !sepastore->initiallp )
            sepastore->ncutsapplied++;
      }
   }

   /* clear the separation storage */
   sepastore->ncuts = 0;
   sepastore->nbdchgs = 0;
   sepastore->ncutsfoundround = 0;
   sepastore->ncutsstored = 0;

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
RETCODE SCIPsepastoreClearCuts(
   SEPASTORE*       sepastore,          /**< separation storage */
   BLKMEM*          blkmem,             /**< block memory */
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
      CHECK_OKAY( SCIProwRelease(&sepastore->cuts[c], blkmem, set, lp) );
   }

   /* clear the separation storage */
   sepastore->ncuts = 0;
   sepastore->nbdchgs = 0;
   sepastore->ncutsfoundround = 0;
   sepastore->ncutsstored = 0;

   return SCIP_OKAY;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreGetNCuts(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncuts + sepastore->nbdchgs;
}

/** get total number of cuts found so far */
int SCIPsepastoreGetNCutsFound(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfound;
}

/** get number of cuts found so far in current separation round */
int SCIPsepastoreGetNCutsFoundRound(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfoundround;
}

/** get total number of cuts stored (and possibly removed again) in current separation round */
int SCIPsepastoreGetNCutsStored(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsstored;
}

/** get total number of cuts applied to the LPs */
int SCIPsepastoreGetNCutsApplied(
   SEPASTORE*            sepastore                /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsapplied;
}
