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
#pragma ident "@(#) $Id: var.c,v 1.99 2004/06/24 15:34:37 bzfpfend Exp $"

/**@file   var.c
 * @brief  methods for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "history.h"
#include "event.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "primal.h"
#include "scip.h"
#include "cons.h"




/*
 * hole, holelist, and domain methods
 */

/** creates a new holelist element */
static
RETCODE holelistCreate(
   HOLELIST**       holelist,           /**< pointer to holelist to create */
   MEMHDR*          memhdr,             /**< block memory for target holelist */
   SET*             set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(holelist != NULL);
   assert(memhdr != NULL);
   assert(SCIPsetIsLT(set, left, right));

   ALLOC_OKAY( allocBlockMemory(memhdr, holelist) );
   (*holelist)->hole.left = left;
   (*holelist)->hole.right = right;
   (*holelist)->next = NULL;

   return SCIP_OKAY;
}

/** frees all elements in the holelist */
static
void holelistFree(
   HOLELIST**       holelist,           /**< pointer to holelist to free */
   MEMHDR*          memhdr              /**< block memory for target holelist */
   )
{
   assert(holelist != NULL);
   assert(memhdr != NULL);

   while( *holelist != NULL )
   {
      HOLELIST* next;

      next = (*holelist)->next;
      freeBlockMemory(memhdr, holelist);
      *holelist = next;
   }
}

/** duplicates a list of holes */
static
RETCODE holelistDuplicate(
   HOLELIST**       target,             /**< pointer to target holelist */
   MEMHDR*          memhdr,             /**< block memory for target holelist */
   SET*             set,                /**< global SCIP settings */
   HOLELIST*        source              /**< holelist to duplicate */
   )
{
   assert(target != NULL);

   while( source != NULL )
   {
      assert(source->next == NULL || SCIPsetIsGE(set, source->next->hole.left, source->hole.right));
      CHECK_OKAY( holelistCreate(target, memhdr, set, source->hole.left, source->hole.right) );
      source = source->next;
      target = &(*target)->next;
   }

   return SCIP_OKAY;
}

/** adds a hole to the domain */
static
RETCODE domAddHole(
   DOM*             dom,                /**< domain to add hole to */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   HOLELIST** insertpos;
   HOLELIST* nextholelist;

   assert(dom != NULL);

   /* sort new hole in holelist of variable */
   insertpos = &dom->holelist;
   while( *insertpos != NULL && (*insertpos)->hole.left < left )
      insertpos = &(*insertpos)->next;

   nextholelist = *insertpos;

   CHECK_OKAY( holelistCreate(insertpos, memhdr, set, left, right) );
   (*insertpos)->next = nextholelist;
   
   return SCIP_OKAY;
}

#if 0 /* for future use */
/** merges overlapping holes into single holes, moves bounds respectively */
static
RETCODE domMerge(
   DOM*             dom,                /**< domain to merge */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   HOLELIST** holelistptr;
   Real* lastrightptr;

   assert(dom != NULL);
   assert(SCIPsetIsLE(set, dom->lb, dom->ub));

   lastrightptr = &dom->lb;  /* lower bound is the right bound of the hole (-infinity,lb) */
   holelistptr = &dom->holelist;

   while( *holelistptr != NULL )
   {
      assert(SCIPsetIsLT(set, (*holelistptr)->hole.left, (*holelistptr)->hole.right));

      if( SCIPsetIsGE(set, (*holelistptr)->hole.left, dom->ub) )
      {
         /* the remaining holes start behind the upper bound: kill them */
         holelistFree(holelistptr, memhdr);
         assert(*holelistptr == NULL);
      }
      else if( SCIPsetIsGT(set, (*holelistptr)->hole.right, dom->ub) )
      {
         /* the hole overlaps the upper bound: decrease upper bound, kill this and all remaining holes */
         dom->ub = (*holelistptr)->hole.left;
         holelistFree(holelistptr, memhdr);
         assert(*holelistptr == NULL);
      }
      else if( SCIPsetIsGT(set, *lastrightptr, (*holelistptr)->hole.left) )
      {
         /* the right bound of the last hole is greater than the left bound of this hole:
          * increase the right bound of the last hole, delete this hole
          */
         HOLELIST* next;

         *lastrightptr = MAX(*lastrightptr, (*holelistptr)->hole.right);
         next = (*holelistptr)->next;
         (*holelistptr)->next = NULL;
         holelistFree(holelistptr, memhdr);
         *holelistptr = next;
      }
      else
      {
         /* the holes do not overlap: update lastrightptr and holelistptr to the next hole */
         lastrightptr = &(*holelistptr)->hole.right;
         holelistptr = &(*holelistptr)->next;
      }
   }

   return SCIP_OKAY;
}
#endif




/*
 * domain change methods
 */

/** applies single bound change */
RETCODE SCIPboundchgApply(
   BOUNDCHG*        boundchg,           /**< bound change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              depth,              /**< depth in the tree, where the bound change takes place */
   int              index               /**< index of the bound change in its bound change array */
   )
{
   assert(boundchg != NULL);
   assert(boundchg->var != NULL);
   assert(SCIPvarGetStatus(boundchg->var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(boundchg->var) == SCIP_VARSTATUS_COLUMN);
   assert(depth >= 0);
   assert(index >= 0);

   /* apply bound change */
   switch( boundchg->boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      if( boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
      {
         debugMessage("branching: new lower bound of <%s> = %g\n", SCIPvarGetName(boundchg->var), boundchg->newbound);
         CHECK_OKAY( SCIPvarChgLbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->newbound,
                        NULL, NULL, 0, depth, index, boundchg->boundchgtype) );
         stat->lastbranchvar = boundchg->var;
         stat->lastbranchdir = SCIP_BRANCHDIR_UPWARDS;
      }
      else
      {
         assert(boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_INFERENCE);
         debugMessage("  -> inference: new lower bound of <%s> = %g\n", SCIPvarGetName(boundchg->var), boundchg->newbound);
         CHECK_OKAY( SCIPvarChgLbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->newbound,
                        boundchg->data.inferencedata.var, boundchg->data.inferencedata.cons,
                        boundchg->data.inferencedata.info, depth, index, boundchg->boundchgtype) );
      }
      break;

   case SCIP_BOUNDTYPE_UPPER:
      if( boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
      {
         debugMessage("branching: new upper bound of <%s> = %g\n", SCIPvarGetName(boundchg->var), boundchg->newbound);
         CHECK_OKAY( SCIPvarChgUbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->newbound, 
                        NULL, NULL, 0, depth, index, boundchg->boundchgtype) );
         stat->lastbranchvar = boundchg->var;
         stat->lastbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
      }
      else
      {
         assert(boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_INFERENCE);
         debugMessage("  -> inference: new upper bound of <%s> = %g\n", SCIPvarGetName(boundchg->var), boundchg->newbound);
         CHECK_OKAY( SCIPvarChgUbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->newbound, 
                        boundchg->data.inferencedata.var, boundchg->data.inferencedata.cons,
                        boundchg->data.inferencedata.info, depth, index, boundchg->boundchgtype) );
      }
      break;

   default:
      errorMessage("unknown bound type\n");
      abort();
   }

   return SCIP_OKAY;
}

/** undoes single bound change */
RETCODE SCIPboundchgUndo(
   BOUNDCHG*        boundchg,           /**< bound change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(boundchg != NULL);
   assert(boundchg->var != NULL);
   assert(SCIPvarGetStatus(boundchg->var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(boundchg->var) == SCIP_VARSTATUS_COLUMN);

   /* undo bound change */
   switch( boundchg->boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      CHECK_OKAY( SCIPvarChgLbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->oldbound, 
                     NULL, NULL, 0, -1, -1, SCIP_BOUNDCHGTYPE_BRANCHING) );
      break;
   case SCIP_BOUNDTYPE_UPPER:
      CHECK_OKAY( SCIPvarChgUbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->oldbound,
                     NULL, NULL, 0, -1, -1, SCIP_BOUNDCHGTYPE_BRANCHING) );
      break;
   default:
      errorMessage("unknown bound type\n");
      abort();
   }

   /* update last branching variable */
   if( boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
      stat->lastbranchvar = NULL;

   return SCIP_OKAY;
}

/** captures branching and inference data of bound change */
static
RETCODE boundchgCaptureData(
   BOUNDCHG*        boundchg            /**< bound change to remove */
   )
{
   assert(boundchg != NULL);

   switch( boundchg->boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      break;

   case SCIP_BOUNDCHGTYPE_INFERENCE:
      assert(boundchg->data.inferencedata.var != NULL);
      if( boundchg->data.inferencedata.cons != NULL )
      {
         SCIPconsCapture(boundchg->data.inferencedata.cons);
      }
      break;
      
   default:
      errorMessage("invalid bound change type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** releases branching and inference data of bound change */
static
RETCODE boundchgReleaseData(
   BOUNDCHG*        boundchg,           /**< bound change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(boundchg != NULL);

   switch( boundchg->boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      break;

   case SCIP_BOUNDCHGTYPE_INFERENCE:
      assert(boundchg->data.inferencedata.var != NULL);
      if( boundchg->data.inferencedata.cons != NULL )
      {
         CHECK_OKAY( SCIPconsRelease(&boundchg->data.inferencedata.cons, memhdr, set) );
      }
      break;
      
   default:
      errorMessage("invalid bound change type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates empty domain change data with dynamic arrays */
static
RETCODE domchgCreate(
   DOMCHG**         domchg,             /**< pointer to domain change data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemorySize(memhdr, domchg, sizeof(DOMCHGDYN)) );
   (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_DYNAMIC; /*lint !e641*/
   (*domchg)->domchgdyn.nboundchgs = 0;
   (*domchg)->domchgdyn.boundchgs = NULL;
   (*domchg)->domchgdyn.nholechgs = 0;
   (*domchg)->domchgdyn.holechgs = NULL;
   (*domchg)->domchgdyn.boundchgssize = 0;
   (*domchg)->domchgdyn.holechgssize = 0;

   return SCIP_OKAY;
}

/** frees domain change data */
RETCODE SCIPdomchgFree(
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   if( *domchg != NULL )
   {
      int i;

      /* release branching and inference data associated with the bound changes */
      for( i = 0; i < (int)(*domchg)->domchgbound.nboundchgs; ++i )
      {
         CHECK_OKAY( boundchgReleaseData(&(*domchg)->domchgbound.boundchgs[i], memhdr, set) );
      }

      /* free memory for bound and hole changes */
      switch( (*domchg)->domchgdyn.domchgtype )
      {
      case SCIP_DOMCHGTYPE_BOUND:
         freeBlockMemoryArrayNull(memhdr, &(*domchg)->domchgbound.boundchgs, (*domchg)->domchgbound.nboundchgs);
         freeBlockMemorySize(memhdr, domchg, sizeof(DOMCHGBOUND));
         break;
      case SCIP_DOMCHGTYPE_BOTH:
         freeBlockMemoryArrayNull(memhdr, &(*domchg)->domchgboth.boundchgs, (*domchg)->domchgboth.nboundchgs);
         freeBlockMemoryArrayNull(memhdr, &(*domchg)->domchgboth.holechgs, (*domchg)->domchgboth.nholechgs);
         freeBlockMemorySize(memhdr, domchg, sizeof(DOMCHGBOTH));
         break;
      case SCIP_DOMCHGTYPE_DYNAMIC:
         freeBlockMemoryArrayNull(memhdr, &(*domchg)->domchgdyn.boundchgs, (*domchg)->domchgdyn.boundchgssize);
         freeBlockMemoryArrayNull(memhdr, &(*domchg)->domchgdyn.holechgs, (*domchg)->domchgdyn.holechgssize);
         freeBlockMemorySize(memhdr, domchg, sizeof(DOMCHGDYN));
         break;
      default:
         errorMessage("invalid domain change type\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** converts a static domain change data into a dynamic one */
static
RETCODE domchgMakeDynamic(
   DOMCHG**         domchg,             /**< pointer to domain change data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   debugMessage("making domain change data %p pointing to %p dynamic\n", domchg, *domchg);

   if( *domchg == NULL )
   {
      CHECK_OKAY( domchgCreate(domchg, memhdr) );
   }
   else
   {
      switch( (*domchg)->domchgdyn.domchgtype )
      {
      case SCIP_DOMCHGTYPE_BOUND:
         ALLOC_OKAY( reallocBlockMemorySize(memhdr, domchg, sizeof(DOMCHGBOUND), sizeof(DOMCHGDYN)) );
         (*domchg)->domchgdyn.nholechgs = 0;
         (*domchg)->domchgdyn.holechgs = NULL;
         (*domchg)->domchgdyn.boundchgssize = (*domchg)->domchgdyn.nboundchgs;
         (*domchg)->domchgdyn.holechgssize = 0;
         (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_DYNAMIC; /*lint !e641*/
         break;
      case SCIP_DOMCHGTYPE_BOTH:
         ALLOC_OKAY( reallocBlockMemorySize(memhdr, domchg, sizeof(DOMCHGBOTH), sizeof(DOMCHGDYN)) );
         (*domchg)->domchgdyn.boundchgssize = (*domchg)->domchgdyn.nboundchgs;
         (*domchg)->domchgdyn.holechgssize = (*domchg)->domchgdyn.nholechgs;
         (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_DYNAMIC; /*lint !e641*/
         break;
      case SCIP_DOMCHGTYPE_DYNAMIC:
         break;
      default:
         errorMessage("invalid domain change type\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** converts a dynamic domain change data into a static one, using less memory than for a dynamic one */
RETCODE SCIPdomchgMakeStatic(
   DOMCHG**         domchg,             /**< pointer to domain change data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   debugMessage("making domain change data %p pointing to %p static\n", domchg, *domchg);

   if( *domchg != NULL )
   {
      switch( (*domchg)->domchgdyn.domchgtype )
      {
      case SCIP_DOMCHGTYPE_BOUND:
         if( (*domchg)->domchgbound.nboundchgs == 0 )
         {
            CHECK_OKAY( SCIPdomchgFree(domchg, memhdr, set) );
         }
         break;
      case SCIP_DOMCHGTYPE_BOTH:
         if( (*domchg)->domchgboth.nholechgs == 0 )
         {
            if( (*domchg)->domchgbound.nboundchgs == 0 )
            {
               CHECK_OKAY( SCIPdomchgFree(domchg, memhdr, set) );
            }
            else
            {
               ALLOC_OKAY( reallocBlockMemorySize(memhdr, domchg, sizeof(DOMCHGBOTH), sizeof(DOMCHGBOUND)) );
               (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_BOUND; /*lint !e641*/
            }
         }
         break;
      case SCIP_DOMCHGTYPE_DYNAMIC:
         if( (*domchg)->domchgboth.nholechgs == 0 )
         {
            if( (*domchg)->domchgbound.nboundchgs == 0 )
            {
               CHECK_OKAY( SCIPdomchgFree(domchg, memhdr, set) );
            }
            else
            {
               /* shrink dynamic size arrays to their minimal sizes */
               ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchg)->domchgdyn.boundchgs,
                              (*domchg)->domchgdyn.boundchgssize, (*domchg)->domchgdyn.nboundchgs) );
               freeBlockMemoryArrayNull(memhdr, &(*domchg)->domchgdyn.holechgs, (*domchg)->domchgdyn.holechgssize);
            
               /* convert into static domain change */
               ALLOC_OKAY( reallocBlockMemorySize(memhdr, domchg, sizeof(DOMCHGDYN), sizeof(DOMCHGBOUND)) );
               (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_BOUND; /*lint !e641*/
            }
         }
         else
         {
            /* shrink dynamic size arrays to their minimal sizes */
            ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchg)->domchgdyn.boundchgs,
                           (*domchg)->domchgdyn.boundchgssize, (*domchg)->domchgdyn.nboundchgs) );
            ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchg)->domchgdyn.holechgs,
                           (*domchg)->domchgdyn.holechgssize, (*domchg)->domchgdyn.nholechgs) );

            /* convert into static domain change */
            ALLOC_OKAY( reallocBlockMemorySize(memhdr, domchg, sizeof(DOMCHGDYN), sizeof(DOMCHGBOTH)) );
            (*domchg)->domchgdyn.domchgtype = SCIP_DOMCHGTYPE_BOTH; /*lint !e641*/
         }
         break;
      default:
         errorMessage("invalid domain change type\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** ensures, that boundchgs array can store at least num entries */
static
RETCODE domchgEnsureBoundchgsSize(
   DOMCHG*          domchg,             /**< domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(domchg != NULL);
   assert(domchg->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   if( num > domchg->domchgdyn.boundchgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &domchg->domchgdyn.boundchgs, domchg->domchgdyn.boundchgssize, newsize) );
      domchg->domchgdyn.boundchgssize = newsize;
   }
   assert(num <= domchg->domchgdyn.boundchgssize);

   return SCIP_OKAY;
}

/** ensures, that holechgs array can store at least num additional entries */
static
RETCODE domchgEnsureHolechgsSize(
   DOMCHG*          domchg,             /**< domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(domchg != NULL);
   assert(domchg->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   if( num > domchg->domchgdyn.holechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &domchg->domchgdyn.holechgs, domchg->domchgdyn.holechgssize, newsize) );
      domchg->domchgdyn.holechgssize = newsize;
   }
   assert(num <= domchg->domchgdyn.holechgssize);

   return SCIP_OKAY;
}

/** applies domain change */
RETCODE SCIPdomchgApply(
   DOMCHG*          domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              depth               /**< depth in the tree, where the domain change takes place */
   )
{
   int i;

   assert(lp != NULL);

   debugMessage("applying domain changes at %p\n", domchg);
   if( domchg == NULL )
      return SCIP_OKAY;

   /* apply bound changes */
   for( i = 0; i < (int)domchg->domchgbound.nboundchgs; ++i )
   {
      CHECK_OKAY( SCIPboundchgApply(&domchg->domchgbound.boundchgs[i], memhdr, set, stat, lp,
                     branchcand, eventqueue, depth, i) );
   }
   debugMessage(" -> %d bound changes\n", domchg->domchgbound.nboundchgs);
   
   /* apply holelist changes */
   if( domchg->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_BOUND ) /*lint !e641*/
   {
      for( i = 0; i < domchg->domchgboth.nholechgs; ++i )
         *(domchg->domchgboth.holechgs[i].ptr) = domchg->domchgboth.holechgs[i].newlist;
      debugMessage(" -> %d hole changes\n", domchg->domchgboth.nholechgs);
   }

   return SCIP_OKAY;
}
   
/** undoes domain change */
RETCODE SCIPdomchgUndo(
   DOMCHG*          domchg,             /**< domain change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   int i;

   assert(lp != NULL);

   debugMessage("undoing domain changes at %p\n", domchg);
   if( domchg == NULL )
      return SCIP_OKAY;

   /* undo holelist changes */
   if( domchg->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_BOUND ) /*lint !e641*/
   {
      for( i = domchg->domchgboth.nholechgs-1; i >= 0; --i )
         *(domchg->domchgboth.holechgs[i].ptr) = domchg->domchgboth.holechgs[i].oldlist;
      debugMessage(" -> %d hole changes\n", domchg->domchgboth.nholechgs);
   }

   /* undo bound changes */
   for( i = domchg->domchgbound.nboundchgs-1; i >= 0; --i )
   {
      CHECK_OKAY( SCIPboundchgUndo(&domchg->domchgbound.boundchgs[i], memhdr, set, stat, lp,
                     branchcand, eventqueue) );
   }
   debugMessage(" -> %d bound changes\n", domchg->domchgbound.nboundchgs);

   return SCIP_OKAY;
}

/** adds bound change to domain changes */
RETCODE SCIPdomchgAddBoundchg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   BOUNDCHGTYPE     boundchgtype,       /**< type of bound change: branching decision or inference */
   Real             lpsolval,           /**< solval of variable in last LP on path to node, or SCIP_INVALID if unknown */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo           /**< user information for inference to help resolving the conflict */
   )
{
   BOUNDCHG* boundchg;

   assert(domchg != NULL);
   assert(stat != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
   assert((int)SCIP_BOUNDTYPE_LOWER == 0 && (int)SCIP_BOUNDTYPE_UPPER == 1); /* must be one bit */
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || boundchgtype == SCIP_BOUNDCHGTYPE_INFERENCE);
   assert((int)SCIP_BOUNDCHGTYPE_BRANCHING == 0 && (int)SCIP_BOUNDCHGTYPE_INFERENCE == 1); /* must be one bit */
   assert(boundchgtype != SCIP_BOUNDCHGTYPE_INFERENCE || infervar != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || infercons == NULL);

   debugMessage("adding %s bound change <%s>: %g -> %g of variable <%s> to domain change at %p pointing to %p\n",
      boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference", 
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", oldbound, newbound, var->name, domchg, *domchg);

   /* if domain change data doesn't exist, create it;
    * if domain change is static, convert it into dynamic change
    */
   if( *domchg == NULL )
   {
      CHECK_OKAY( domchgCreate(domchg, memhdr) );
   }
   else if( (*domchg)->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_DYNAMIC ) /*lint !e641*/
   {
      CHECK_OKAY( domchgMakeDynamic(domchg, memhdr) );
   }
   assert(*domchg != NULL && (*domchg)->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   /* get memory for additional bound change */
   CHECK_OKAY( domchgEnsureBoundchgsSize(*domchg, memhdr, set, (*domchg)->domchgdyn.nboundchgs+1) );

   /* fill in the bound change data */
   boundchg = &(*domchg)->domchgdyn.boundchgs[(*domchg)->domchgdyn.nboundchgs];
   boundchg->var = var;
   switch( boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      boundchg->data.branchingdata.lpsolval = lpsolval;
      break;
   case SCIP_BOUNDCHGTYPE_INFERENCE:
      boundchg->data.inferencedata.var = infervar;
      boundchg->data.inferencedata.cons = infercons;
      boundchg->data.inferencedata.info = inferinfo; 
      break;
   default:
      errorMessage("invalid bound change type\n");
      return SCIP_INVALIDDATA;
   }

   boundchg->newbound = newbound;
   boundchg->oldbound = oldbound;
   boundchg->boundtype = boundtype; /*lint !e641*/
   boundchg->boundchgtype = boundchgtype; /*lint !e641*/
   (*domchg)->domchgdyn.nboundchgs++;

   /* capture branching and inference data associated with the bound changes */
   CHECK_OKAY( boundchgCaptureData(boundchg) );

   stat->nboundchgs++;

   return SCIP_OKAY;
}

/** adds hole change to domain changes */
RETCODE SCIPdomchgAddHolechg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   HOLECHG* holechg;

   assert(domchg != NULL);
   assert(ptr != NULL);

   /* if domain change data doesn't exist, create it;
    * if domain change is static, convert it into dynamic change
    */
   if( *domchg == NULL )
   {
      CHECK_OKAY( domchgCreate(domchg, memhdr) );
   }
   else if( (*domchg)->domchgdyn.domchgtype != SCIP_DOMCHGTYPE_DYNAMIC ) /*lint !e641*/
   {
      CHECK_OKAY( domchgMakeDynamic(domchg, memhdr) );
   }
   assert(*domchg != NULL && (*domchg)->domchgdyn.domchgtype == SCIP_DOMCHGTYPE_DYNAMIC); /*lint !e641*/

   /* get memory for additional hole change */
   CHECK_OKAY( domchgEnsureHolechgsSize(*domchg, memhdr, set, (*domchg)->domchgdyn.nholechgs+1) );

   /* fill in the hole change data */
   holechg = &(*domchg)->domchgdyn.holechgs[(*domchg)->domchgdyn.nholechgs];
   holechg->ptr = ptr;
   holechg->newlist = newlist;
   holechg->oldlist = oldlist;
   (*domchg)->domchgdyn.nholechgs++;

   stat->nholechgs++;

   return SCIP_OKAY;
}




/*
 * methods for variable bounds
 */

/** creates a variable bounds data structure */
static
RETCODE vboundsCreate(
   VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(vbounds != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, vbounds) );
   (*vbounds)->vars = NULL;
   (*vbounds)->coefs = NULL;
   (*vbounds)->constants = NULL;
   (*vbounds)->len = 0;
   (*vbounds)->size = 0;

   return SCIP_OKAY;
}

/** frees a variable bounds data structure */
static
void vboundsFree(
   VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(vbounds != NULL);

   if( *vbounds != NULL )
   {
      freeBlockMemoryArrayNull(memhdr, &(*vbounds)->vars, (*vbounds)->size);
      freeBlockMemoryArrayNull(memhdr, &(*vbounds)->coefs, (*vbounds)->size);
      freeBlockMemoryArrayNull(memhdr, &(*vbounds)->constants, (*vbounds)->size);
      freeBlockMemory(memhdr, vbounds);
   }
}

/** ensures, that variable bounds arrays can store at least num entries */
static
RETCODE vboundsEnsureSize(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(vbounds != NULL);
   
   /* create variable bounds data structure, if not yet existing */
   if( *vbounds == NULL )
   {
      CHECK_OKAY( vboundsCreate(vbounds, memhdr) );
   }
   assert(*vbounds != NULL);
   assert((*vbounds)->len <= (*vbounds)->size);

   if( num > (*vbounds)->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*vbounds)->vars, (*vbounds)->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*vbounds)->coefs, (*vbounds)->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*vbounds)->constants, (*vbounds)->size, newsize) );
      (*vbounds)->size = newsize;
   }
   assert(num <= (*vbounds)->size);

   return SCIP_OKAY;
}

/** adds a variable bound to the variable bounds data structure */
static
RETCODE vboundsAdd(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   Real             coef,               /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   Real             constant            /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   )
{
   assert(vbounds != NULL);
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE);
   assert(var->vartype != SCIP_VARTYPE_CONTINUOUS);

   CHECK_OKAY( vboundsEnsureSize(vbounds, memhdr, set, *vbounds != NULL ? (*vbounds)->len+1 : 1) );
   assert(*vbounds != NULL);

   (*vbounds)->vars[(*vbounds)->len] = var;
   (*vbounds)->coefs[(*vbounds)->len] = coef;
   (*vbounds)->constants[(*vbounds)->len] = constant;
   (*vbounds)->len++;

   return SCIP_OKAY;
}

/** replaces bounding variables in variable bounds by their active problem variable counterparts */
static
RETCODE vboundsUseActiveVars(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   int i;

   if( vbounds == NULL )
      return SCIP_OKAY;

   for( i = 0; i < vbounds->len; ++i )
   {
      /* transform linear sum  b*z + d  into  b'*z' + d'  with active problem variable z' */
      CHECK_OKAY( SCIPvarGetProbvarSum(&vbounds->vars[i], &vbounds->coefs[i], &vbounds->constants[i]) );
      
      /* if the bounding variable was reduced to a constant, remove the entry from the vbounds */
      if( vbounds->vars[i] == NULL )
      {
         vbounds->vars[i] = vbounds->vars[vbounds->len-1];
         vbounds->coefs[i] = vbounds->coefs[vbounds->len-1];
         vbounds->constants[i] = vbounds->constants[vbounds->len-1];
         vbounds->len--;
         i--;
      }
   }

   return SCIP_OKAY;
}




/*
 * methods for variables 
 */

/** returns adjusted lower bound value, which is rounded for integral variable types */
static
Real adjustedLb(
   SET*             set,                /**< global SCIP settings */
   VARTYPE          vartype,            /**< type of variable */
   Real             lb                  /**< lower bound to adjust */
   )
{
   if( SCIPsetIsInfinity(set, -lb) )
      return -set->infinity;
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      return SCIPsetCeil(set, lb);
   else if( SCIPsetIsZero(set, lb) )
      return 0.0;
   else
      return lb;
}

/** returns adjusted upper bound value, which is rounded for integral variable types */
static
Real adjustedUb(
   SET*             set,                /**< global SCIP settings */
   VARTYPE          vartype,            /**< type of variable */
   Real             ub                  /**< upper bound to adjust */
   )
{
   if( SCIPsetIsInfinity(set, ub) )
      return set->infinity;
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      return SCIPsetFloor(set, ub);
   else if( SCIPsetIsZero(set, ub) )
      return 0.0;
   else
      return ub;
}

/** creates variable; if variable is of integral type, fractional bounds are automatically rounded; an integer variable
 *  with bounds zero and one is automatically converted into a binary variable
 */
static
RETCODE varCreate(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

   /* adjust bounds of variable */
   lb = adjustedLb(set, vartype, lb);
   ub = adjustedUb(set, vartype, ub);

   /* convert [0,1]-integers into binary variables */
   if( vartype == SCIP_VARTYPE_INTEGER
      && (SCIPsetIsEQ(set, lb, 0.0) || SCIPsetIsEQ(set, lb, 1.0))
      && (SCIPsetIsEQ(set, ub, 0.0) || SCIPsetIsEQ(set, ub, 1.0)) )
      vartype = SCIP_VARTYPE_BINARY;

   ALLOC_OKAY( allocBlockMemory(memhdr, var) );

   if( name == NULL )
   {
      char s[MAXSTRLEN];
      sprintf(s, "_var%d_", stat->nvaridx);
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*var)->name, s, strlen(s)+1) );
   }
   else
   {
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*var)->name, name, strlen(name)+1) );
   }
   (*var)->parentvars = NULL;
   (*var)->negatedvar = NULL;
   (*var)->eventfilter = NULL;
   (*var)->infervar = NULL;
   (*var)->infercons = NULL;
   (*var)->inferinfo = 0;
   (*var)->glbdom.holelist = NULL;
   (*var)->glbdom.lb = lb;
   (*var)->glbdom.ub = ub;
   (*var)->locdom.holelist = NULL;
   (*var)->locdom.lb = lb;
   (*var)->locdom.ub = ub;
   (*var)->vlbs = NULL;
   (*var)->vubs = NULL;
   (*var)->obj = obj;
   (*var)->branchfactor = 1.0;
   (*var)->branchpriority = 0;
   (*var)->rootsol = SCIP_INVALID;
   (*var)->index = stat->nvaridx++;
   (*var)->probindex = -1;
   (*var)->pseudocandindex = -1;
   (*var)->eventqueueindexobj = -1;
   (*var)->eventqueueindexlb = -1;
   (*var)->eventqueueindexub = -1;
   (*var)->parentvarssize = 0;
   (*var)->nparentvars = 0;
   (*var)->nuses = 0;
   (*var)->nlocksdown = 0;
   (*var)->nlocksup = 0;
   (*var)->fixdepth = SCIPsetIsEQ(set, lb, ub) ? 0 : -1;
   (*var)->fixindex = -1;
   (*var)->conflictsetcount = 0;
   (*var)->boundchgtype = SCIP_BOUNDCHGTYPE_BRANCHING;
   (*var)->initial = initial;
   (*var)->removeable = removeable;
   (*var)->vartype = vartype; /*lint !e641*/
   (*var)->pseudocostflag = FALSE;
   (*var)->vardelorig = vardelorig;
   (*var)->vartrans = vartrans;
   (*var)->vardeltrans = vardeltrans;
   (*var)->vardata = vardata;

   /* adjust bounds, if variable is integral */
   SCIPvarAdjustLb(*var, set, &lb);
   SCIPvarAdjustUb(*var, set, &ub);
   assert(lb <= ub);

   /* create branching and inference history entries */
   CHECK_OKAY( SCIPhistoryCreate(&(*var)->history, memhdr) );

   return SCIP_OKAY;
}

/** creates and captures an original problem variable; an integer variable with bounds
 *  zero and one is automatically converted into a binary variable
 */
RETCODE SCIPvarCreateOriginal(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

   /* create variable */
   CHECK_OKAY( varCreate(var, memhdr, set, stat, name, lb, ub, obj, vartype, initial, removeable,
                  vardelorig, vartrans, vardeltrans, vardata) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_ORIGINAL; /*lint !e641*/
   (*var)->data.transvar = NULL;

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;
}

/** creates and captures a loose variable belonging to the transformed problem; an integer variable with bounds
 *  zero and one is automatically converted into a binary variable
 */
RETCODE SCIPvarCreateTransformed(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);

   /* create variable */
   CHECK_OKAY( varCreate(var, memhdr, set, stat, name, lb, ub, obj, vartype, initial, removeable,
                  vardelorig, vartrans, vardeltrans, vardata) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_LOOSE; /*lint !e641*/

   /* create event filter for transformed variable */
   CHECK_OKAY( SCIPeventfilterCreate(&(*var)->eventfilter, memhdr) );

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;   
}

/** ensures, that parentvars array of var can store at least num entries */
static
RETCODE varEnsureParentvarsSize(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(var->nparentvars <= var->parentvarssize);
   
   if( num > var->parentvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &var->parentvars, var->parentvarssize, newsize) );
      var->parentvarssize = newsize;
   }
   assert(num <= var->parentvarssize);

   return SCIP_OKAY;
}

/** adds variable to parent list of a variable and captures parent variable */
static
RETCODE varAddParent(
   VAR*             var,                /**< variable to add parent to */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   VAR*             parentvar           /**< parent variable to add */
   )
{
   assert(var != NULL);
   assert(parentvar != NULL);

   debugMessage("adding parent <%s>[%p] to variable <%s>[%p] in slot %d\n", 
      parentvar->name, parentvar, var->name, var, var->nparentvars);

   CHECK_OKAY( varEnsureParentvarsSize(var, memhdr, set, var->nparentvars+1) );

   var->parentvars[var->nparentvars] = parentvar;
   var->nparentvars++;

   SCIPvarCapture(parentvar);

   return SCIP_OKAY;
}

/** deletes and releases all variables from the parent list of a variable, frees the memory of parents array */
static
RETCODE varFreeParents(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data (or NULL, if it's an original variable) */
   )
{
   VAR* parentvar;
   int i;

   debugMessage("free parents of <%s>\n", (*var)->name);

   /* release the parent variables and remove the link from the parent variable to the child */
   for( i = 0; i < (*var)->nparentvars; ++i )
   {
      assert((*var)->parentvars != NULL);
      parentvar = (*var)->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         assert(parentvar->data.transvar == *var);
         assert(&parentvar->data.transvar != var);
         parentvar->data.transvar = NULL;
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         assert(parentvar->data.aggregate.var == *var);
         assert(&parentvar->data.aggregate.var != var);
         parentvar->data.aggregate.var = NULL;
         break;

#if 0
      case SCIP_VARSTATUS_MULTAGGR:
         assert(parentvar->data.multaggr.vars != NULL);
         for( v = 0; v < parentvar->data.multaggr.nvars && parentvar->data.multaggr.vars[v] != *var; ++v )
         {}
         assert(v < parentvar->data.multaggr.nvars && parentvar->data.multaggr.vars[v] == *var);
         if( v < parentvar->data.multaggr.nvars-1 )
         {
            parentvar->data.multaggr.vars[v] = parentvar->data.multaggr.vars[parentvar->data.multaggr.nvars-1];
            parentvar->data.multaggr.scalars[v] = parentvar->data.multaggr.scalars[parentvar->data.multaggr.nvars-1];
         }
         parentvar->data.multaggr.nvars--;
         break;
#endif

      case SCIP_VARSTATUS_NEGATED:
         assert(parentvar->negatedvar == *var);
         assert((*var)->negatedvar == parentvar);
         parentvar->negatedvar = NULL;
         (*var)->negatedvar = NULL;
         break;

      default:
         errorMessage("parent variable is neither ORIGINAL, AGGREGATED nor NEGATED\n");
         return SCIP_INVALIDDATA;
      }  /*lint !e788*/

      CHECK_OKAY( SCIPvarRelease(&(*var)->parentvars[i], memhdr, set, lp) );
   }

   /* free parentvars array */
   freeBlockMemoryArrayNull(memhdr, &(*var)->parentvars, (*var)->parentvarssize);

   return SCIP_OKAY;
}

/** frees a variable */
static
RETCODE varFree(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data (may be NULL, if it's not a column variable) */
   )
{
   assert(memhdr != NULL);
   assert(var != NULL);
   assert(*var != NULL);
   assert(SCIPvarGetStatus(*var) != SCIP_VARSTATUS_COLUMN || &(*var)->data.col->var != var);
   assert((*var)->nuses == 0);
   assert((*var)->probindex == -1);

   debugMessage("free variable <%s> with status=%d\n", (*var)->name, SCIPvarGetStatus(*var));
   switch( (*var)->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert((*var)->data.transvar == NULL);  /* cannot free variable, if transformed variable is still existing */
      break;
   case SCIP_VARSTATUS_LOOSE:
      break;
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( SCIPcolFree(&(*var)->data.col, memhdr, set, lp) );  /* free corresponding LP column */
      break;
   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
      break;
   case SCIP_VARSTATUS_MULTAGGR:
      freeBlockMemoryArray(memhdr, &(*var)->data.multaggr.vars, (*var)->data.multaggr.varssize);
      freeBlockMemoryArray(memhdr, &(*var)->data.multaggr.scalars, (*var)->data.multaggr.varssize);
      break;
   case SCIP_VARSTATUS_NEGATED:
      break;
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   /* release all parent variables and free the parentvars array */
   CHECK_OKAY( varFreeParents(var, memhdr, set, lp) );

   /* free user data */
   if( (*var)->varstatus == SCIP_VARSTATUS_ORIGINAL )
   {
      if( (*var)->vardelorig != NULL )
      {
         CHECK_OKAY( (*var)->vardelorig(set->scip, *var, &(*var)->vardata) );
      }
   }
   else
   {
      if( (*var)->vardeltrans != NULL )
      {
         CHECK_OKAY( (*var)->vardeltrans(set->scip, *var, &(*var)->vardata) );
      }
   }

   /* free event filter */
   if( (*var)->eventfilter != NULL )
   {
      assert(SCIPvarGetStatus(*var) != SCIP_VARSTATUS_ORIGINAL);
      CHECK_OKAY( SCIPeventfilterFree(&(*var)->eventfilter, memhdr, set) );
   }
   assert((*var)->eventfilter == NULL);

   /* free variable bounds data structures */
   vboundsFree(&(*var)->vlbs, memhdr);
   vboundsFree(&(*var)->vubs, memhdr);

   /* free branching and inference history entries */
   SCIPhistoryFree(&(*var)->history, memhdr);

   /* free variable data structure */
   freeBlockMemoryArray(memhdr, &(*var)->name, strlen((*var)->name)+1);
   freeBlockMemory(memhdr, var);

   return SCIP_OKAY;
}

/** increases usage counter of variable */
void SCIPvarCapture(
   VAR*             var                 /**< variable */
   )
{
   assert(var != NULL);
   assert(var->nuses >= 0);

   debugMessage("capture variable <%s> with nuses=%d\n", var->name, var->nuses);
   var->nuses++;
}

/** decreases usage counter of variable, and frees memory if necessary */
RETCODE SCIPvarRelease(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data (or NULL, if it's an original variable) */
   )
{
   assert(memhdr != NULL);
   assert(var != NULL);
   assert(*var != NULL);
   assert((*var)->nuses >= 1);

   debugMessage("release variable <%s> with nuses=%d\n", (*var)->name, (*var)->nuses);
   (*var)->nuses--;
   if( (*var)->nuses == 0 )
   {
      CHECK_OKAY( varFree(var, memhdr, set, lp) );
   }

   *var = NULL;

   return SCIP_OKAY;
}

/** outputs variable information into file stream */
void SCIPvarPrint(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(var != NULL);

   if( file == NULL )
      file = stdout;

   switch( var->vartype )
   {
   case SCIP_VARTYPE_BINARY:
      fprintf(file, "  [binary]");
      break;
   case SCIP_VARTYPE_INTEGER:
      fprintf(file, "  [integer]");
      break;
   case SCIP_VARTYPE_IMPLINT:
      fprintf(file, "  [implicit]");
      break;
   case SCIP_VARTYPE_CONTINUOUS:
      fprintf(file, "  [continuous]");
      break;
   default:
      errorMessage("unknown variable type\n");
      abort();
   }

   fprintf(file, " <%s>: ", var->name);

   if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
      fprintf(file, " [-inf,");
   else
      fprintf(file, " [%g,", var->glbdom.lb);
   if( SCIPsetIsInfinity(set, var->glbdom.ub) )
      fprintf(file, "+inf]");
   else
      fprintf(file, "%g]", var->glbdom.ub);
   
   /**@todo print holes */

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      break;

   case SCIP_VARSTATUS_FIXED:
      fprintf(file, " == %g", var->glbdom.lb);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      fprintf(file, " ==");
      if( !SCIPsetIsZero(set, var->data.aggregate.constant) )
         fprintf(file, " %g", var->data.aggregate.constant);
      fprintf(file, " %+g<%s>", var->data.aggregate.scalar, SCIPvarGetName(var->data.aggregate.var));
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      fprintf(file, " ==");
      if( !SCIPsetIsZero(set, var->data.multaggr.constant) )
         fprintf(file, " %g", var->data.multaggr.constant);
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         fprintf(file, " %+g<%s>", var->data.multaggr.scalars[i], SCIPvarGetName(var->data.multaggr.vars[i]));
      break;

   case SCIP_VARSTATUS_NEGATED:
      fprintf(file, " == %g - <%s>", var->data.negate.constant, SCIPvarGetName(var->negatedvar));
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   fprintf(file, "\n");
}

/** increases lock numbers for rounding */
static
void varAddRoundLocks(
   VAR*             var,                /**< problem variable */
   int              addnlocksdown,      /**< increase in number of rounding down locks */
   int              addnlocksup         /**< increase in number of rounding up locks */
   )
{
   int i;

   assert(var != NULL);
   assert(var->nlocksup >= 0);
   assert(var->nlocksdown >= 0);

   debugMessage("add rounding locks %d/%d to variable <%s> (locks=%d/%d)\n",
      addnlocksdown, addnlocksup, var->name, var->nlocksdown, var->nlocksup);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         varAddRoundLocks(var->data.transvar, addnlocksdown, addnlocksup);
      else
      {
         var->nlocksdown += addnlocksdown;
         var->nlocksup += addnlocksup;
      }
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      var->nlocksdown += addnlocksdown;
      var->nlocksup += addnlocksup;
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         varAddRoundLocks(var->data.aggregate.var, addnlocksdown, addnlocksup);
      else
         varAddRoundLocks(var->data.aggregate.var, addnlocksup, addnlocksdown);
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            varAddRoundLocks(var->data.multaggr.vars[i], addnlocksdown, addnlocksup);
         else
            varAddRoundLocks(var->data.multaggr.vars[i], addnlocksup, addnlocksdown);
      }
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      varAddRoundLocks(var->negatedvar, addnlocksup, addnlocksdown);
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   assert(var->nlocksdown >= 0);
   assert(var->nlocksup >= 0);
}

/** increases lock number for rounding down by one; tells variable, that rounding its value down will make the
 *  solution infeasible
 */
void SCIPvarLockDown(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("forbid rounding down of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 1, 0);
}

/** increases lock number for rounding up by one; tells variable, that rounding its value up will make the
 *  solution infeasible
 */
void SCIPvarLockUp(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("forbid rounding up of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 0, 1);
}

/** increases lock number for rounding down and up by one; tells variable, that rounding value in either direction will
 *  make the solution infeasible
 */
void SCIPvarLockBoth(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("forbid both roundings of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 1, 1);
}

/** declares that rounding down the given variable would destroy the feasibility of the given constraint;
 *  locks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
void SCIPvarLockDownCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   varAddRoundLocks(var, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
}

/** declares that rounding up the given variable would destroy the feasibility of the given constraint;
 *  locks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
void SCIPvarLockUpCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   varAddRoundLocks(var, (int)SCIPconsIsLockedNeg(cons), (int)SCIPconsIsLockedPos(cons));
}

/** declares that rounding the given variable in any direction would destroy the feasibility of the given constraint;
 *  locks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
void SCIPvarLockBothCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   )
{
   int nlocks;

   assert(cons != NULL);

   nlocks = (int)SCIPconsIsLockedPos(cons) + (int)SCIPconsIsLockedNeg(cons);
   varAddRoundLocks(var, nlocks, nlocks);
}

/** increases lock number for roundings of variable; tells variable, that rounding value in a direction set to
 *  a positive value will make the solution infeasible
 */
void SCIPvarLock(
   VAR*             var,                /**< problem variable */
   int              nlocksdown,         /**< increase in number of rounding down locks */
   int              nlocksup            /**< increase in number of rounding up locks */
   )
{
   debugMessage("forbid rounding (%d/%d) of <%s> (locks=%d/%d)\n", 
      nlocksdown, nlocksup, var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, nlocksdown, nlocksup);
}

/** decreases lock number for rounding down by one; cancels a prior SCIPvarLockDown() */
void SCIPvarUnlockDown(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("allow rounding down of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, -1, 0);
}

/** decreases lock number for rounding up by one; cancels a prior SCIPvarLockUp() */
void SCIPvarUnlockUp(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("allow rounding up of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 0, -1);
}

/** decreases lock number for rounding down and up by one; cancels a prior SCIPvarLockBoth() */
void SCIPvarUnlockBoth(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("allow both roundings of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, -1, -1);
}

/** declares that rounding down the given variable would no longer destroy the feasibility of the given constraint;
 *  unlocks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
void SCIPvarUnlockDownCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   varAddRoundLocks(var, -(int)SCIPconsIsLockedPos(cons), -(int)SCIPconsIsLockedNeg(cons));
}

/** declares that rounding up the given variable would no longer destroy the feasibility of the given constraint;
 *  unlocks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
void SCIPvarUnlockUpCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   varAddRoundLocks(var, -(int)SCIPconsIsLockedNeg(cons), -(int)SCIPconsIsLockedPos(cons));
}

/** declares that rounding the given variable in any direction would no longer destroy the feasibility of the given
 *  constraint; unlocks the roundings of the variable corresponding to the lock status of the constraint and its negation
 */
void SCIPvarUnlockBothCons(
   VAR*             var,                /**< problem variable */
   CONS*            cons                /**< constraint */
   )
{
   int nlocks;

   assert(cons != NULL);

   nlocks = (int)SCIPconsIsLockedPos(cons) + (int)SCIPconsIsLockedNeg(cons);
   varAddRoundLocks(var, -nlocks, -nlocks);
}

/** decreases lock number for roundings of variable; cancels a prior call to SCIPvarLock() */
void SCIPvarUnlock(
   VAR*             var,                /**< problem variable */
   int              nunlocksdown,       /**< decrease in number of rounding down locks */
   int              nunlocksup          /**< decrease in number of rounding up locks */
   )
{
   debugMessage("allow rounding (%d/%d) of <%s> (locks=%d/%d)\n", 
      nunlocksdown, nunlocksup, var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, -nunlocksdown, -nunlocksup);
}

/** gets number of locks for rounding down */
int SCIPvarGetNLocksDown(
   VAR*             var                 /**< problem variable */
   )
{
   int nlocks;
   int i;

   assert(var != NULL);
   assert(var->nlocksdown >= 0);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         return SCIPvarGetNLocksDown(var->data.transvar);
      else
         return var->nlocksdown;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      return var->nlocksdown;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNLocksDown(var->data.aggregate.var);
      else
         return SCIPvarGetNLocksUp(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      nlocks = 0;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            nlocks += SCIPvarGetNLocksDown(var->data.multaggr.vars[i]);
         else
            nlocks += SCIPvarGetNLocksUp(var->data.multaggr.vars[i]);
      }
      return nlocks;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return SCIPvarGetNLocksUp(var->negatedvar);

   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets number of locks for rounding up */
int SCIPvarGetNLocksUp(
   VAR*             var                 /**< problem variable */
   )
{
   int nlocks;
   int i;

   assert(var != NULL);
   assert(var->nlocksup >= 0);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         return SCIPvarGetNLocksUp(var->data.transvar);
      else
         return var->nlocksup;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      return var->nlocksup;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNLocksUp(var->data.aggregate.var);
      else
         return SCIPvarGetNLocksDown(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      nlocks = 0;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            nlocks += SCIPvarGetNLocksUp(var->data.multaggr.vars[i]);
         else
            nlocks += SCIPvarGetNLocksDown(var->data.multaggr.vars[i]);
      }
      return nlocks;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return SCIPvarGetNLocksDown(var->negatedvar);

   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** is it possible, to round variable down and stay feasible? */
Bool SCIPvarMayRoundDown(
   VAR*             var                 /**< problem variable */
   )
{
   return (SCIPvarGetNLocksDown(var) == 0);
}

/** is it possible, to round variable up and stay feasible? */
Bool SCIPvarMayRoundUp(
   VAR*             var                 /**< problem variable */
   )
{
   return (SCIPvarGetNLocksUp(var) == 0);
}

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
RETCODE SCIPvarTransform(
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   char name[MAXSTRLEN];

   assert(origvar != NULL);
   assert(SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_ORIGINAL);
   assert(origvar->glbdom.lb == origvar->locdom.lb); /*lint !e777*/
   assert(origvar->glbdom.ub == origvar->locdom.ub); /*lint !e777*/
   assert(origvar->vlbs == NULL);
   assert(origvar->vubs == NULL);
   assert(transvar != NULL);

   /* check if variable is already transformed */
   if( origvar->data.transvar != NULL )
   {
      *transvar = origvar->data.transvar;
      SCIPvarCapture(*transvar);
   }
   else
   {
      /* create transformed variable */
      sprintf(name, "t_%s", origvar->name);
      CHECK_OKAY( SCIPvarCreateTransformed(transvar, memhdr, set, stat, name,
                     origvar->glbdom.lb, origvar->glbdom.ub, (Real)objsense * origvar->obj,
                     origvar->vartype, origvar->initial, origvar->removeable,
                     NULL, NULL, origvar->vardeltrans, origvar->vardata) );
      
      /* copy the branch factor and priority */
      (*transvar)->branchfactor = origvar->branchfactor;
      (*transvar)->branchpriority = origvar->branchpriority;

      /* duplicate hole lists */
      CHECK_OKAY( holelistDuplicate(&(*transvar)->glbdom.holelist, memhdr, set, origvar->glbdom.holelist) );
      CHECK_OKAY( holelistDuplicate(&(*transvar)->locdom.holelist, memhdr, set, origvar->locdom.holelist) );
      
      /* link original and transformed variable */
      origvar->data.transvar = *transvar;
      CHECK_OKAY( varAddParent(*transvar, memhdr, set, origvar) );
      
      /* copy rounding locks */
      (*transvar)->nlocksdown = origvar->nlocksdown;
      (*transvar)->nlocksup = origvar->nlocksup;
      assert((*transvar)->nlocksdown >= 0);
      assert((*transvar)->nlocksup >= 0);

      /* transform user data */
      if( origvar->vartrans != NULL )
      {
         CHECK_OKAY( origvar->vartrans(set->scip, origvar, origvar->vardata, *transvar, &(*transvar)->vardata) );
      }
   }

   debugMessage("transformed variable: <%s>[%p] -> <%s>[%p]\n", origvar->name, origvar, (*transvar)->name, *transvar);

   return SCIP_OKAY;
}

/** gets corresponding transformed variable of an original or negated original variable */
RETCODE SCIPvarGetTransformed(
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR**            transvar            /**< pointer to store the transformed variable, or NULL if not existing yet */
   )
{
   assert(origvar != NULL);
   assert(SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_NEGATED);

   if( SCIPvarGetStatus(origvar) == SCIP_VARSTATUS_NEGATED )
   {
      assert(origvar->negatedvar != NULL);
      assert(SCIPvarGetStatus(origvar->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

      if( origvar->negatedvar->data.transvar == NULL )
         *transvar = NULL;
      else
      {
         CHECK_OKAY( SCIPvarNegate(origvar->negatedvar->data.transvar, memhdr, set, stat, transvar) );
      }
   }
   else
      *transvar = origvar->data.transvar;

   return SCIP_OKAY;
}

/** converts loose transformed variable into column variable, creates LP column */
RETCODE SCIPvarColumn(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< current LP data */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   debugMessage("creating column for variable <%s>\n", var->name);

   /* switch variable status */
   var->varstatus = SCIP_VARSTATUS_COLUMN; /*lint !e641*/

   /* create column of variable */
   CHECK_OKAY( SCIPcolCreate(&var->data.col, memhdr, set, stat, var, 0, NULL, NULL, var->removeable) );

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, NULL, var) );
      
      /* inform LP, that problem variable is now a column variable and no longer loose */
      CHECK_OKAY( SCIPlpUpdateVarColumn(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** converts column transformed variable back into loose variable, frees LP column */
RETCODE SCIPvarLoose(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< current LP data */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(var->data.col != NULL);
   assert(var->data.col->len == var->data.col->nunlinked); /* no rows can include the column */
   assert(var->data.col->lppos == -1);
   assert(var->data.col->lpipos == -1);

   debugMessage("deleting column for variable <%s>\n", var->name);

   /* free column of variable */
   CHECK_OKAY( SCIPcolFree(&var->data.col, memhdr, set, lp) );

   /* switch variable status */
   var->varstatus = SCIP_VARSTATUS_LOOSE; /*lint !e641*/

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, NULL, var) );
      
      /* inform LP, that problem variable is now a loose variable and no longer a column */
      CHECK_OKAY( SCIPlpUpdateVarLoose(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** issues a VARFIXED event on the given variable and all its parents (except ORIGINAL parents);
 *  the event issuing on the parents is necessary, because unlike with bound changes, the parent variables
 *  are not informed about a fixing of an active variable they are pointing to
 */
static
RETCODE varEventVarFixed(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   EVENT* event;
   int i;

   assert(var != NULL);

   /* issue VARFIXED event on variable */
   CHECK_OKAY( SCIPeventCreateVarFixed(&event, memhdr, var) );
   CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, NULL, NULL, NULL, NULL, &event) );

   /* process all parents */
   for( i = 0; i < var->nparentvars; ++i )
   {
      if( var->parentvars[i]->varstatus != SCIP_VARSTATUS_ORIGINAL )
      {
         CHECK_OKAY( varEventVarFixed(var->parentvars[i], memhdr, set, eventqueue) );
      }
   }

   return SCIP_OKAY;
}

/** converts variable into fixed variable */
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             fixedval,           /**< value to fix variable at */
   Bool*            infeasible          /**< pointer to store whether the fixing is infeasible */
   )
{
   Real obj;

   assert(var != NULL);
   assert(var->glbdom.lb == var->locdom.lb); /*lint !e777*/
   assert(var->glbdom.ub == var->locdom.ub); /*lint !e777*/
   assert(SCIPsetIsLE(set, var->glbdom.lb, fixedval));
   assert(SCIPsetIsLE(set, fixedval, var->glbdom.ub));
   assert(infeasible != NULL);
   
   debugMessage("fix variable <%s>[%g,%g] to %g\n", var->name, var->glbdom.lb, var->glbdom.ub, fixedval);

   if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPsetIsIntegral(set, fixedval))
      || SCIPsetIsFeasLT(set, fixedval, var->locdom.lb)
      || SCIPsetIsFeasGT(set, fixedval, var->locdom.ub) )
   {
      debugMessage(" -> fixing infeasible: locdom=[%g,%g], fixedval=%g\n", var->locdom.lb, var->locdom.ub, fixedval);
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
   {
      *infeasible = !SCIPsetIsFeasEQ(set, fixedval, var->locdom.lb);
      debugMessage(" -> variable already fixed to %g (fixedval=%g): infeasible=%d\n", 
         var->locdom.lb, fixedval, *infeasible);
      return SCIP_OKAY;
   }

   *infeasible = FALSE;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("Cannot fix an untransformed original variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarFix(var->data.transvar, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, 
                     fixedval, infeasible) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* set the fixed variable's objective value to 0.0 */
      obj = var->obj;
      CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, primal, lp, eventqueue, 0.0) );

      /* change variable's bounds to fixed value */
      holelistFree(&var->glbdom.holelist, memhdr);
      CHECK_OKAY( SCIPvarSetLbGlobal(var, set, fixedval) );
      CHECK_OKAY( SCIPvarSetUbGlobal(var, set, fixedval) );
      holelistFree(&var->locdom.holelist, memhdr);
      CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, fixedval,
                     NULL, NULL, 0, 0, -1, SCIP_BOUNDCHGTYPE_INFERENCE) );
      CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, fixedval,
                     NULL, NULL, 0, 0, -1, SCIP_BOUNDCHGTYPE_INFERENCE) );

      /* explicitly set variable's bounds, even if the fixed value is in epsilon range of the old bound */
      var->glbdom.lb = fixedval;
      var->glbdom.ub = fixedval;
      var->locdom.lb = fixedval;
      var->locdom.ub = fixedval;

      /* delete variable bounds information */
      vboundsFree(&var->vlbs, memhdr);
      vboundsFree(&var->vubs, memhdr);

      /* clear the history of the variable */
      SCIPhistoryReset(var->history);

      /* convert variable into fixed variable */
      var->varstatus = SCIP_VARSTATUS_FIXED; /*lint !e641*/

      /* inform problem about the variable's status change */
      if( var->probindex != -1 )
      {
         CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, branchcand, var) );
      }

      /* issue VARFIXED event */
      CHECK_OKAY( varEventVarFixed(var, memhdr, set, eventqueue) );

      /* reset the objective value of the fixed variable, thus adjusting the problem's objective offset */
      CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, stat, prob, primal, lp, eventqueue, obj) );
      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("Cannot fix a column variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot fix a fixed variable again\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      /* fix aggregation variable y in x = a*y + c, instead of fixing x directly */
      assert(SCIPsetIsZero(set, var->obj));
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      CHECK_OKAY( SCIPvarFix(var->data.aggregate.var, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue,
                     (fixedval - var->data.aggregate.constant)/var->data.aggregate.scalar, infeasible) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot fix a multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      /* fix negation variable x in x' = offset - x, instead of fixing x' directly */
      assert(SCIPsetIsZero(set, var->obj));
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarFix(var->negatedvar, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue,
                     var->data.negate.constant - fixedval, infeasible) );
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }
   
   return SCIP_OKAY;
}

/** tightens the bounds of both variables in aggregation x = a*y + c */
static
RETCODE varUpdateAggregationBounds(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PRIMAL*          primal,             /**< primal data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             aggvar,             /**< variable y in aggregation x = a*y + c */
   Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   Real varlb;
   Real varub;
   Real aggvarlb;
   Real aggvarub;
   Bool aggvarbdschanged;

   assert(var != NULL);
   assert(aggvar != NULL);
   assert(!SCIPsetIsZero(set, scalar));
   assert(infeasible != NULL);

   debugMessage("updating bounds of variables in aggregation <%s> == %g*<%s> %+g\n",
      var->name, scalar, aggvar->name, constant);
   debugMessage("  old bounds: <%s> [%g,%g]   <%s> [%g,%g]\n",
      var->name, var->glbdom.lb, var->glbdom.ub, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub);

   /* loop as long additional changes may be found */
   do
   {
      aggvarbdschanged = FALSE;

      /* update the bounds of the aggregated variable x in x = a*y + c */
      if( scalar > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -aggvar->glbdom.lb) )
            varlb = -set->infinity;
         else
            varlb = aggvar->glbdom.lb * scalar + constant;
         if( SCIPsetIsInfinity(set, aggvar->glbdom.ub) )
            varub = set->infinity;
         else
            varub = aggvar->glbdom.ub * scalar + constant;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -aggvar->glbdom.lb) )
            varub = set->infinity;
         else
            varub = aggvar->glbdom.lb * scalar + constant;
         if( SCIPsetIsInfinity(set, aggvar->glbdom.ub) )
            varlb = -set->infinity;
         else
            varlb = aggvar->glbdom.ub * scalar + constant;
      }
      varlb = MAX(varlb, var->glbdom.lb);
      varub = MIN(varub, var->glbdom.ub);
      SCIPvarAdjustLb(var, set, &varlb);
      SCIPvarAdjustUb(var, set, &varub);

      /* check the new bounds */
      if( SCIPsetIsGT(set, varlb, varub) )
      {
         /* the aggregation is infeasible */
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPsetIsEQ(set, varlb, varub) )
      {
         /* the aggregated variable is fixed -> fix both variables */
         CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, varlb, infeasible) );
         if( !(*infeasible) )
         {
            CHECK_OKAY( SCIPvarFix(aggvar, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, 
                           (varlb-constant)/scalar, infeasible) );
         }
         return SCIP_OKAY;
      }
      else
      {
         if( SCIPsetIsGT(set, varlb, var->glbdom.lb) )
         {
            CHECK_OKAY( SCIPvarSetLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, varlb) );
            CHECK_OKAY( SCIPvarSetLbGlobal(var, set, varlb) );
         }
         if( SCIPsetIsLT(set, varub, var->glbdom.ub) )
         {
            CHECK_OKAY( SCIPvarSetUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, varub) );
            CHECK_OKAY( SCIPvarSetUbGlobal(var, set, varub) );
         }

         /* update the hole list of the aggregation variable */
         /**@todo update hole list of aggregation variable */
      }

      /* update the bounds of the aggregation variable y in x = a*y + c  ->  y = (x-c)/a */
      if( scalar > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarlb = -set->infinity;
         else
            aggvarlb = (var->glbdom.lb - constant) / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarub = set->infinity;
         else
            aggvarub = (var->glbdom.ub - constant) / scalar;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarub = set->infinity;
         else
            aggvarub = (var->glbdom.lb - constant) / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarlb = -set->infinity;
         else
            aggvarlb = (var->glbdom.ub - constant) / scalar;
      }
      aggvarlb = MAX(aggvarlb, aggvar->glbdom.lb);
      aggvarub = MIN(aggvarub, aggvar->glbdom.ub);
      SCIPvarAdjustLb(aggvar, set, &aggvarlb);
      SCIPvarAdjustUb(aggvar, set, &aggvarub);

      /* check the new bounds */
      if( SCIPsetIsGT(set, aggvarlb, aggvarub) )
      {
         /* the aggregation is infeasible */
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPsetIsEQ(set, aggvarlb, aggvarub) )
      {
         /* the aggregation variable is fixed -> fix both variables */
         CHECK_OKAY( SCIPvarFix(aggvar, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, aggvarlb, 
                        infeasible) );
         if( !(*infeasible) )
         {
            CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, 
                           aggvarlb * scalar + constant, infeasible) );
         }
         return SCIP_OKAY;
      }
      else
      {
         if( SCIPsetIsGT(set, aggvarlb, aggvar->glbdom.lb) )
         {
            CHECK_OKAY( SCIPvarSetLbLocal(aggvar, memhdr, set, stat, lp, branchcand, eventqueue, aggvarlb) );
            CHECK_OKAY( SCIPvarSetLbGlobal(aggvar, set, aggvarlb) );
            aggvarbdschanged = TRUE;
         }
         if( SCIPsetIsLT(set, aggvarub, aggvar->glbdom.ub) )
         {
            CHECK_OKAY( SCIPvarSetUbLocal(aggvar, memhdr, set, stat, lp, branchcand, eventqueue, aggvarub) );
            CHECK_OKAY( SCIPvarSetUbGlobal(aggvar, set, aggvarub) );
            aggvarbdschanged = TRUE;
         }

         /* update the hole list of the aggregation variable */
         /**@todo update hole list of aggregation variable */
      }
   }
   while( aggvarbdschanged );

   debugMessage("  new bounds: <%s> [%g,%g]   <%s> [%g,%g]\n",
      var->name, var->glbdom.lb, var->glbdom.ub, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub);

   return SCIP_OKAY;
}

/** converts loose variable into aggregated variable */
RETCODE SCIPvarAggregate(
   VAR*             var,                /**< loose problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             aggvar,             /**< loose variable y in aggregation x = a*y + c */
   Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   EVENT* event;
   VAR** vars;
   Real* coefs;
   Real* constants;
   Real obj;
   Real branchfactor;
   int branchpriority;
   int nlocksdown;
   int nlocksup;
   int nvbds;
   int i;

   assert(var != NULL);
   assert(var->glbdom.lb == var->locdom.lb); /*lint !e777*/
   assert(var->glbdom.ub == var->locdom.ub); /*lint !e777*/
   assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */
   assert(infeasible != NULL);

   /* aggregation is a fixing, if the scalar is zero */
   if( SCIPsetIsZero(set, scalar) )
      return SCIPvarFix(var, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, constant, infeasible);

   assert(aggvar != NULL);
   assert(aggvar->glbdom.lb == aggvar->locdom.lb); /*lint !e777*/
   assert(aggvar->glbdom.ub == aggvar->locdom.ub); /*lint !e777*/

   debugMessage("aggregate variable <%s>[%g,%g] == %g*<%s>[%g,%g] %+g\n", var->name, var->glbdom.lb, var->glbdom.ub,
      scalar, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub, constant);

   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetStatus(aggvar) == SCIP_VARSTATUS_LOOSE);

   *infeasible = FALSE;

   /* tighten the bounds of aggregated and aggregation variable */
   CHECK_OKAY( varUpdateAggregationBounds(var, memhdr, set, stat, prob, lp, branchcand, primal, eventqueue,
                  aggvar, scalar, constant, infeasible) );
   if( *infeasible )
      return SCIP_OKAY;

   /* set the aggregated variable's objective value to 0.0 */
   obj = var->obj;
   CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, primal, lp, eventqueue, 0.0) );

   /* unlock all rounding locks */
   nlocksdown = var->nlocksdown;
   nlocksup = var->nlocksup;
   var->nlocksdown = 0;
   var->nlocksup = 0;

   /* convert variable into aggregated variable */
   var->varstatus = SCIP_VARSTATUS_AGGREGATED; /*lint !e641*/
   var->data.aggregate.var = aggvar;
   var->data.aggregate.scalar = scalar;
   var->data.aggregate.constant = constant;

   /* make aggregated variable a parent of the aggregation variable */
   CHECK_OKAY( varAddParent(aggvar, memhdr, set, var) );

   /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variable */
   varAddRoundLocks(var, nlocksdown, nlocksup);

   /* move the variable bounds to the aggregation variable:
    *  - add all variable bounds again to the variable, thus adding it to the aggregated variable
    *  - free the variable bounds data structures
    */
   if( var->vlbs != NULL )
   {
      nvbds = var->vlbs->len;
      vars = var->vlbs->vars;
      coefs = var->vlbs->coefs;
      constants = var->vlbs->constants;
      for( i = 0; i < nvbds; ++i )
      {
         CHECK_OKAY( SCIPvarAddVlb(var, memhdr, set, vars[i], coefs[i], constants[i]) );
      }
   }
   if( var->vubs != NULL )
   {
      nvbds = var->vubs->len;
      vars = var->vubs->vars;
      coefs = var->vubs->coefs;
      constants = var->vubs->constants;
      for( i = 0; i < nvbds; ++i )
      {
         CHECK_OKAY( SCIPvarAddVub(var, memhdr, set, vars[i], coefs[i], constants[i]) );
      }
   }
   vboundsFree(&var->vlbs, memhdr);
   vboundsFree(&var->vubs, memhdr);

   /* add the history entries to the aggregation variable and clear the history of the aggregated variable */
   SCIPhistoryUnite(aggvar->history, var->history, scalar < 0.0);
   SCIPhistoryReset(var->history);

   /* update flags of aggregation variable */
   aggvar->removeable &= var->removeable;

   /* update branching factors and priorities of both variables to be the maximum of both variables */
   branchfactor = MAX(aggvar->branchfactor, var->branchfactor);
   branchpriority = MAX(aggvar->branchpriority, var->branchpriority);
   SCIPvarChgBranchFactor(aggvar, set, branchfactor);
   SCIPvarChgBranchPriority(aggvar, set, branchpriority);
   SCIPvarChgBranchFactor(var, set, branchfactor);
   SCIPvarChgBranchPriority(var, set, branchpriority);

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, branchcand, var) );
   }

   /* issue VARFIXED event */
   CHECK_OKAY( SCIPeventCreateVarFixed(&event, memhdr, var) );
   CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, NULL, NULL, NULL, NULL, &event) );

   /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
    * variable and the problem's objective offset
    */
   CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, stat, prob, primal, lp, eventqueue, obj) );

   return SCIP_OKAY;
}

/** converts variable into multi-aggregated variable */
RETCODE SCIPvarMultiaggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   EVENT* event;
   Real obj;
   Real branchfactor;
   int branchpriority;
   int nlocksdown;
   int nlocksup;
   int v;

   assert(var != NULL);
   assert(var->glbdom.lb == var->locdom.lb); /*lint !e777*/
   assert(var->glbdom.ub == var->locdom.ub); /*lint !e777*/
   assert(naggvars == 0 || aggvars != NULL);
   assert(naggvars == 0 || scalars != NULL);
   assert(infeasible != NULL);

   debugMessage("multi-aggregate variable <%s> == ...%d vars... %+g\n", var->name, naggvars, constant);

   *infeasible = FALSE;

   /* check, if we are in one of the simple cases */
   if( naggvars == 0 )
      return SCIPvarFix(var, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, constant, infeasible);
   else if( naggvars == 1)
      return SCIPvarAggregate(var, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue, 
         aggvars[0], scalars[0], constant, infeasible);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("Cannot multi-aggregate an untransformed original variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarMultiaggregate(var->data.transvar, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue,
                     naggvars, aggvars, scalars, constant, infeasible) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* set the aggregated variable's objective value to 0.0 */
      obj = var->obj;
      CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, primal, lp, eventqueue, 0.0) );

      /* unlock all rounding locks */
      nlocksdown = var->nlocksdown;
      nlocksup = var->nlocksup;
      var->nlocksdown = 0;
      var->nlocksup = 0;

      /* convert variable into multi-aggregated variable */
      var->varstatus = SCIP_VARSTATUS_MULTAGGR; /*lint !e641*/
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &var->data.multaggr.vars, aggvars, naggvars) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &var->data.multaggr.scalars, scalars, naggvars) );
      var->data.multaggr.constant = constant;
      var->data.multaggr.nvars = naggvars;
      var->data.multaggr.varssize = naggvars;

      /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variables */
      varAddRoundLocks(var, nlocksdown, nlocksup);

      /* free the variable bounds data structures */
      vboundsFree(&var->vlbs, memhdr);
      vboundsFree(&var->vubs, memhdr);

      /* update flags and branching factors and priorities of aggregation variables */
      branchfactor = var->branchfactor;
      branchpriority = var->branchpriority;
      for( v = 0; v < naggvars; ++v )
      {
         aggvars[v]->removeable &= var->removeable;
         branchfactor = MAX(aggvars[v]->branchfactor, branchfactor);
         branchpriority = MAX(aggvars[v]->branchpriority, branchpriority);
      }
      for( v = 0; v < naggvars; ++v )
      {
         SCIPvarChgBranchFactor(aggvars[v], set, branchfactor);
         SCIPvarChgBranchPriority(aggvars[v], set, branchpriority);
      }
      SCIPvarChgBranchFactor(var, set, branchfactor);
      SCIPvarChgBranchPriority(var, set, branchpriority);

      if( var->probindex != -1 )
      {
         /* inform problem about the variable's status change */
         CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, branchcand, var) );
      }

      /* issue VARFIXED event */
      CHECK_OKAY( SCIPeventCreateVarFixed(&event, memhdr, var) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, NULL, NULL, NULL, NULL, &event) );

      /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
       * variables and the problem's objective offset
       */
      CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, stat, prob, primal, lp, eventqueue, obj) );

      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("Cannot multi-aggregate a column variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot multi-aggregate a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      errorMessage("Cannot multi-aggregate an aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot multi-aggregate a multiple aggregated variable again\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      /* aggregate negation variable x in x' = offset - x, instead of aggregating x' directly:
       *   x' = a_1*y_1 + ... + a_n*y_n + c  ->  x = offset - x' = offset - a_1*y_1 - ... - a_n*y_n - c
       */
      assert(SCIPsetIsZero(set, var->obj));
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);

      /* switch the signs of the aggregation scalars */
      for( v = 0; v < naggvars; ++v )
         scalars[v] *= -1.0;
         
      /* perform the multi aggregation on the negation variable */
      CHECK_OKAY( SCIPvarMultiaggregate(var->negatedvar, memhdr, set, stat, prob, primal, lp, branchcand, eventqueue,
                     naggvars, aggvars, scalars, var->data.negate.constant - constant, infeasible) );

      /* switch the signs of the aggregation scalars again, to reset them to their original values */
      for( v = 0; v < naggvars; ++v )
         scalars[v] *= -1.0;
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }
   
   return SCIP_OKAY;
}

/** gets negated variable x' = offset - x of problem variable x; the negated variable is created if not yet existing;
 *  the negation offset of binary variables is always 1, the offset of other variables is fixed to lb + ub when the
 *  negated variable is created
 */
RETCODE SCIPvarNegate(
   VAR*             var,                /**< problem variable to negate */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR**            negvar              /**< pointer to store the negated variable */
   )
{
   assert(var != NULL);
   assert(negvar != NULL);

   /* check, if we already created the negated variable */
   if( var->negatedvar == NULL )
   {
      char negvarname[MAXSTRLEN];

      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED);

      debugMessage("creating negated variable of <%s>\n", var->name);
      
      /* negation is only possible for bounded variables */
      if( SCIPsetIsInfinity(set, -var->glbdom.lb) || SCIPsetIsInfinity(set, var->glbdom.ub) )
      {
         errorMessage("cannot negate unbounded variable\n");
         return SCIP_INVALIDDATA;
      }

      sprintf(negvarname, "%s_neg", var->name);

      /* create negated variable */
      CHECK_OKAY( varCreate(negvar, memhdr, set, stat, negvarname, var->glbdom.lb, var->glbdom.ub, 0.0,
                     SCIPvarGetType(var), var->initial, var->removeable, NULL, NULL, NULL, NULL) );
      (*negvar)->varstatus = SCIP_VARSTATUS_NEGATED; /*lint !e641*/
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         (*negvar)->data.negate.constant = 1.0;
      else
         (*negvar)->data.negate.constant = var->glbdom.lb + var->glbdom.ub;

      /* set the bounds corresponding to the negation variable */
      (*negvar)->glbdom.lb = (*negvar)->data.negate.constant - var->glbdom.ub;
      (*negvar)->glbdom.ub = (*negvar)->data.negate.constant - var->glbdom.lb;
      (*negvar)->locdom.lb = (*negvar)->data.negate.constant - var->locdom.ub;
      (*negvar)->locdom.ub = (*negvar)->data.negate.constant - var->locdom.lb;
      /**@todo create holes in the negated variable corresponding to the holes of the negation variable */
      
      /* copy the fixing information to the negated variable */
      (*negvar)->fixdepth = var->fixdepth;
      (*negvar)->fixindex = var->fixindex;

      /* create event filter, if the negated variable belongs to the transformed problem */
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         CHECK_OKAY( SCIPeventfilterCreate(&(*negvar)->eventfilter, memhdr) );
      }

      /* link the variables together */
      var->negatedvar = *negvar;
      (*negvar)->negatedvar = var;

      /* copy the branch factor and priority */
      (*negvar)->branchfactor = var->branchfactor;
      (*negvar)->branchpriority = var->branchpriority;

      /* make negated variable a parent of the negation variable (negated variable is captured as a parent) */
      CHECK_OKAY( varAddParent(var, memhdr, set, *negvar) );
      assert((*negvar)->nuses == 1);
   }
   assert(var->negatedvar != NULL);

   /* return the negated variable */
   *negvar = var->negatedvar;

   /* exactly one variable of the negation pair has to be marked as negated variable */
   assert((SCIPvarGetStatus(*negvar) == SCIP_VARSTATUS_NEGATED) != (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED));

   return SCIP_OKAY;
}

/** changes type of variable; cannot be called, if var belongs to a problem */
RETCODE SCIPvarChgType(
   VAR*             var,                /**< problem variable to change */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(var != NULL);

   debugMessage("change type of <%s> from %d to %d\n", var->name, SCIPvarGetType(var), vartype);

   if( var->probindex >= 0 )
   {
      errorMessage("cannot change type of variable already in the problem\n");
      return SCIP_INVALIDDATA;
   }
   
   var->vartype = vartype; /*lint !e641*/
   if( var->negatedvar != NULL )
      var->negatedvar->vartype = vartype; /*lint !e641*/

   return SCIP_OKAY;
}

/** appends OBJCHANGED event to the event queue */
static
RETCODE varEventObjChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             oldobj,             /**< old objective value for variable */
   Real             newobj              /**< new objective value for variable */
   )
{
   EVENT* event;

   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert(!SCIPsetIsEQ(set, oldobj, newobj));

   CHECK_OKAY( SCIPeventCreateObjChanged(&event, memhdr, var, oldobj, newobj) );
   CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, primal, lp, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** changes objective value of variable */
RETCODE SCIPvarChgObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newobj              /**< new objective value for variable */
   )
{
   Real oldobj;

   assert(var != NULL);

   debugMessage("changing objective value of <%s> from %g to %g\n", var->name, var->obj, newobj);

   if( !SCIPsetIsEQ(set, var->obj, newobj) )
   {
      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarChgObj(var->data.transvar, memhdr, set, primal, lp, eventqueue, newobj) );
         }
         else
         {
            assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
            var->obj = newobj;
         }
         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         oldobj = var->obj;
         var->obj = newobj;
         CHECK_OKAY( varEventObjChanged(var, memhdr, set, primal, lp, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_MULTAGGR:
      case SCIP_VARSTATUS_NEGATED:
         errorMessage("cannot change objective value of a fixed, aggregated, multi-aggregated, or negated variable\n");
         return SCIP_INVALIDDATA;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** adds value to objective value of variable */
RETCODE SCIPvarAddObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             addobj              /**< additional objective value for variable */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(SCIPstage(set->scip) < SCIP_STAGE_INITSOLVE);

   debugMessage("adding %g to objective value %g of <%s>\n", addobj, var->obj, var->name);

   if( !SCIPsetIsZero(set, addobj) )
   {
      Real oldobj;
      int i;

      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarAddObj(var->data.transvar, memhdr, set, stat, prob, primal, lp, eventqueue, addobj) );
         }
         else
         {
            assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
            var->obj += addobj;
         }
         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         oldobj = var->obj;
         var->obj += addobj;
         CHECK_OKAY( varEventObjChanged(var, memhdr, set, primal, lp, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
         assert(SCIPsetIsEQ(set, var->locdom.lb, var->locdom.ub));
         SCIPprobAddObjoffset(prob, var->locdom.lb * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, NULL, lp) );
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         /* x = a*y + c  ->  add a*addobj to obj. val. of y, and c*addobj to obj. offset of problem */
         SCIPprobAddObjoffset(prob, var->data.aggregate.constant * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, NULL, lp) );
         CHECK_OKAY( SCIPvarAddObj(var->data.aggregate.var, memhdr, set, stat, prob, primal, lp, eventqueue,
                        var->data.aggregate.scalar * addobj) );
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         /* x = a_1*y_1 + ... + a_n*y_n  + c  ->  add a_i*addobj to obj. val. of y_i, and c*addobj to obj. offset */
         SCIPprobAddObjoffset(prob, var->data.multaggr.constant * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, NULL, lp) );
         for( i = 0; i < var->data.multaggr.nvars; ++i )
         {
            CHECK_OKAY( SCIPvarAddObj(var->data.multaggr.vars[i], memhdr, set, stat, prob, primal, lp, 
                           eventqueue, var->data.multaggr.scalars[i] * addobj) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED:
         /* x' = offset - x  ->  add -addobj to obj. val. of x and offset*addobj to obj. offset of problem */
         assert(var->negatedvar != NULL);
         assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(var->negatedvar->negatedvar == var);
         SCIPprobAddObjoffset(prob, var->data.negate.constant * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, NULL, lp) );
         CHECK_OKAY( SCIPvarAddObj(var->negatedvar, memhdr, set, stat, prob, primal, lp, eventqueue,
                        -addobj) );
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes objective value of variable in current dive */
RETCODE SCIPvarChgObjDive(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newobj              /**< new objective value for variable */
   )
{
   assert(var != NULL);
   assert(lp != NULL);

   debugMessage("changing objective of <%s> to %g in current dive\n", var->name, newobj);

   if( SCIPsetIsZero(set, newobj) )
      newobj = 0.0;

   /* mark the LP's objective function invalid */
   lp->divingobjchg = TRUE;

   /* change objective value of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      CHECK_OKAY( SCIPvarChgObjDive(var->data.transvar, set, lp, newobj) );
      break;
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      CHECK_OKAY( SCIPcolChgObj(var->data.col, set, lp, newobj) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      /* nothing to do here: only the constant shift in objective function would change */
      break;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      CHECK_OKAY( SCIPvarChgObjDive(var->data.aggregate.var, set, lp, newobj / var->data.aggregate.scalar) );
      /* the constant can be ignored, because it would only affect the objective shift */
      break;
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot change diving objective value of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgObjDive(var->negatedvar, set, lp, -newobj) );
      /* the offset can be ignored, because it would only affect the objective shift */
      break;
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** adjust lower bound to integral value, if variable is integral */
void SCIPvarAdjustLb(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real*            lb                  /**< pointer to lower bound to adjust */
   )
{
   assert(var != NULL);
   assert(lb != NULL);

   debugMessage("adjust lower bound %g of <%s>\n", *lb, var->name);

   *lb = adjustedLb(set, var->vartype, *lb);
}

/** adjust upper bound to integral value, if variable is integral */
void SCIPvarAdjustUb(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real*            ub                  /**< pointer to upper bound to adjust */
   )
{
   assert(var != NULL);
   assert(ub != NULL);

   debugMessage("adjust upper bound %g of <%s>\n", *ub, var->name);

   *ub = adjustedUb(set, var->vartype, *ub);
}


/* forward declaration, because both methods call each other recursively */

/** performs the current change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   );

/** performs the current change in lower bound, changes all parents accordingly */
static            
RETCODE varProcessChgLbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(SCIPsetIsEQ(set, newbound, adjustedLb(set, var->vartype, newbound)));

   debugMessage("process changing global lower bound of <%s> from %f to %f\n", var->name, var->glbdom.lb, newbound);

   if( SCIPsetIsEQ(set, newbound, var->glbdom.lb) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->glbdom.lb;
   var->glbdom.lb = newbound;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change bounds across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            assert((SCIPsetIsInfinity(set, -parentvar->glbdom.lb) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsEQ(set, parentvar->glbdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else
         {
            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, parentvar->glbdom.ub) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsEQ(set, parentvar->glbdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         CHECK_OKAY( varProcessChgUbGlobal(parentvar, set, parentvar->data.negate.constant - newbound) );
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** performs the current change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(SCIPsetIsEQ(set, newbound, adjustedUb(set, var->vartype, newbound)));

   debugMessage("process changing global upper bound of <%s> from %f to %f\n", var->name, var->glbdom.ub, newbound);

   if( SCIPsetIsEQ(set, newbound, var->glbdom.ub) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->glbdom.ub;
   var->glbdom.ub = newbound;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of y */
            assert((SCIPsetIsInfinity(set, parentvar->glbdom.ub) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsEQ(set, parentvar->glbdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else 
         {
            /* a < 0 -> change lower bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, -parentvar->glbdom.lb) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsEQ(set, parentvar->glbdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         CHECK_OKAY( varProcessChgLbGlobal(parentvar, set, parentvar->data.negate.constant - newbound) );
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** sets global lower bound of variable; if possible, adjusts bound to integral value */
RETCODE SCIPvarSetLbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   /* adjust bound for integral variables */
   SCIPvarAdjustLb(var, set, &newbound);

   debugMessage("changing global lower bound of <%s> from %g to %g\n", var->name, var->glbdom.lb, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->glbdom.lb, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarSetLbGlobal(var->data.transvar, set, newbound) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         CHECK_OKAY( varProcessChgLbGlobal(var, set, newbound) );
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      CHECK_OKAY( varProcessChgLbGlobal(var, set, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, -var->glbdom.lb) && SCIPsetIsInfinity(set, -var->data.aggregate.var->glbdom.lb))
            || SCIPsetIsEQ(set, var->glbdom.lb,
               var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarSetLbGlobal(var->data.aggregate.var, set,
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, -var->glbdom.lb) && SCIPsetIsInfinity(set, var->data.aggregate.var->glbdom.ub))
            || SCIPsetIsEQ(set, var->glbdom.lb,
               var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarSetUbGlobal(var->data.aggregate.var, set,
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo change the sides of the corresponding linear constraint */
      errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarSetUbGlobal(var->negatedvar, set,  var->data.negate.constant - newbound) );
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** sets global upper bound of variable; if possible, adjusts bound to integral value */
RETCODE SCIPvarSetUbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   /* adjust bound for integral variables */
   SCIPvarAdjustUb(var, set, &newbound);

   debugMessage("changing global upper bound of <%s> from %g to %g\n", var->name, var->glbdom.ub, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->glbdom.ub, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarSetUbGlobal(var->data.transvar, set, newbound) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         CHECK_OKAY( varProcessChgUbGlobal(var, set, newbound) );
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      CHECK_OKAY( varProcessChgUbGlobal(var, set, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, var->glbdom.ub) && SCIPsetIsInfinity(set, var->data.aggregate.var->glbdom.ub))
            || SCIPsetIsEQ(set, var->glbdom.ub,
               var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarSetUbGlobal(var->data.aggregate.var, set,
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, var->glbdom.ub) && SCIPsetIsInfinity(set, -var->data.aggregate.var->glbdom.lb))
            || SCIPsetIsEQ(set, var->glbdom.ub,
               var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarSetLbGlobal(var->data.aggregate.var, set,
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo change the sides of the corresponding linear constraint */
      errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarSetLbGlobal(var->negatedvar, set,  var->data.negate.constant - newbound) );
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** sets global bound of variable; if possible, adjusts bound to integral value */
RETCODE SCIPvarSetBdGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarSetLbGlobal(var, set, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarSetUbGlobal(var, set, newbound);
   default:
      errorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }
}

/** appends LBTIGHTENED or LBRELAXED event to the event queue */
static
RETCODE varEventLbChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             oldbound,           /**< old lower bound for variable */
   Real             newbound            /**< new lower bound for variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert(!SCIPsetIsEQ(set, oldbound, newbound));

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because row activities and LP changes have to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_LBCHANGED) != 0)
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;

      debugMessage("issue LBCHANGED event for variable <%s>: %g -> %g\n", var->name, oldbound, newbound);

      CHECK_OKAY( SCIPeventCreateLbChanged(&event, memhdr, var, oldbound, newbound) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, NULL, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/** appends UBTIGHTENED or UBRELAXED event to the event queue */
static
RETCODE varEventUbChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             oldbound,           /**< old upper bound for variable */
   Real             newbound            /**< new upper bound for variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert(!SCIPsetIsEQ(set, oldbound, newbound));

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because row activities and LP changes have to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_UBCHANGED) != 0)
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;
   
      debugMessage("issue UBCHANGED event for variable <%s>: %g -> %g\n", var->name, oldbound, newbound);

      CHECK_OKAY( SCIPeventCreateUbChanged(&event, memhdr, var, oldbound, newbound) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, NULL, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/* forward declaration, because both methods call each other recursively */

/** performs the current change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place */
   int              fixindex,           /**< bound change index for each node representing the order of changes */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   );

/** performs the current change in lower bound, changes all parents accordingly */
static            
RETCODE varProcessChgLbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place */
   int              fixindex,           /**< bound change index for each node representing the order of changes */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(SCIPsetIsEQ(set, newbound, adjustedLb(set, var->vartype, newbound)));

   debugMessage("process changing lower bound of <%s> from %g to %g (depth:%d, index:%d)\n",
      var->name, var->locdom.lb, newbound, fixdepth, fixindex);

   if( SCIPsetIsEQ(set, newbound, var->locdom.lb) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->locdom.lb;
   var->locdom.lb = newbound;
   var->infervar = infervar;
   var->infercons = infercons;
   var->inferinfo = inferinfo;
   var->fixdepth = fixdepth;
   var->fixindex = fixindex;
   var->boundchgtype = boundchgtype;

   /* issue bound change event */
   assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL) == (var->eventfilter == NULL));
   if( var->eventfilter != NULL )
   {
      CHECK_OKAY( varEventLbChanged(var, memhdr, set, lp, branchcand, eventqueue, oldbound, newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change bounds across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            assert((SCIPsetIsInfinity(set, -parentvar->locdom.lb) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant, 
                           infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
         }
         else 
         {
            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, parentvar->locdom.ub) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant, 
                           infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue,
                        parentvar->data.negate.constant - newbound, 
                        infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** performs the current change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place */
   int              fixindex,           /**< bound change index for each node representing the order of changes */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(SCIPsetIsEQ(set, newbound, adjustedUb(set, var->vartype, newbound)));

   debugMessage("process changing upper bound of <%s> from %g to %g (depth:%d, index:%d)\n",
      var->name, var->locdom.ub, newbound, fixdepth, fixindex);

   if( SCIPsetIsEQ(set, newbound, var->locdom.ub) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->locdom.ub;
   var->locdom.ub = newbound;
   var->infervar = infervar;
   var->infercons = infercons;
   var->inferinfo = inferinfo;
   var->fixdepth = fixdepth;
   var->fixindex = fixindex;
   var->boundchgtype = boundchgtype;

   /* issue bound change event */
   assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL) == (var->eventfilter == NULL));
   if( var->eventfilter != NULL )
   {
      CHECK_OKAY( varEventUbChanged(var, memhdr, set, lp, branchcand, eventqueue, oldbound, newbound) );
   }

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change bounds across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of x */
            assert((SCIPsetIsInfinity(set, parentvar->locdom.ub) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant, 
                           infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
         }
         else
         {
            /* a < 0 -> change lower bound of x */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, -parentvar->locdom.lb) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant, 
                           infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue,
                        parentvar->data.negate.constant - newbound, 
                        infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes current local lower bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
RETCODE SCIPvarChgLbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place, or -1 */
   int              fixindex,           /**< bound change index for each node representing the order of changes, or -1 */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   )
{
   assert(var != NULL);

   /* adjust bound for integral variables */
   SCIPvarAdjustLb(var, set, &newbound);

   debugMessage("changing lower bound of <%s> from %g to %g (depth:%d, index:%d)\n",
      var->name, var->locdom.lb, newbound, fixdepth, fixindex);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   assert(SCIPstage(set->scip) != SCIP_STAGE_SOLVING || SCIPgetDepth(set->scip) == 0
      || var->vartype != SCIP_VARTYPE_BINARY || newbound < 0.5 || fixdepth >= 0);

   if( SCIPsetIsEQ(set, var->locdom.lb, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.transvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                        newbound, infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         stat->nboundchanges++;
         var->locdom.lb = newbound;
         var->infervar = infervar;
         var->infercons = infercons;
         var->inferinfo = inferinfo;
         var->fixdepth = fixdepth;
         var->fixindex = fixindex;
         var->boundchgtype = boundchgtype;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      stat->nboundchanges++;
      CHECK_OKAY( varProcessChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, 
                     newbound, infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, -var->locdom.lb) && SCIPsetIsInfinity(set, -var->data.aggregate.var->locdom.lb))
            || SCIPsetIsEQ(set, var->locdom.lb,
               var->data.aggregate.var->locdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.aggregate.var, memhdr, set, stat, lp, branchcand, eventqueue, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar, 
                        infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, -var->locdom.lb) && SCIPsetIsInfinity(set, var->data.aggregate.var->locdom.ub))
            || SCIPsetIsEQ(set, var->locdom.lb,
               var->data.aggregate.var->locdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.aggregate.var, memhdr, set, stat, lp, branchcand, eventqueue, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar, 
                        infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo change the sides of the corresponding linear constraint */
      errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgUbLocal(var->negatedvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                     var->data.negate.constant - newbound, 
                     infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      break;
         
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
   
   return SCIP_OKAY;
}

/** changes current local upper bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
RETCODE SCIPvarChgUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place, or -1 */
   int              fixindex,           /**< bound change index for each node representing the order of changes, or -1 */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   )
{
   assert(var != NULL);

   /* adjust bound for integral variables */
   SCIPvarAdjustUb(var, set, &newbound);

   debugMessage("changing upper bound of <%s> from %g to %g (depth:%d, index:%d)\n",
      var->name, var->locdom.ub, newbound, fixdepth, fixindex);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;
   
   assert(SCIPstage(set->scip) != SCIP_STAGE_SOLVING || SCIPgetDepth(set->scip) == 0
      || var->vartype != SCIP_VARTYPE_BINARY || newbound > 0.5 || fixdepth >= 0);

   if( SCIPsetIsEQ(set, var->locdom.ub, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.transvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                        newbound, infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         stat->nboundchanges++;
         var->locdom.ub = newbound;
         var->infervar = infervar;
         var->infercons = infercons;
         var->inferinfo = inferinfo;
         var->fixdepth = fixdepth;
         var->fixindex = fixindex;
         var->boundchgtype = boundchgtype;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      stat->nboundchanges++;
      CHECK_OKAY( varProcessChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, 
                     newbound, infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, var->locdom.ub) && SCIPsetIsInfinity(set, var->data.aggregate.var->locdom.ub))
            || SCIPsetIsEQ(set, var->locdom.ub,
               var->data.aggregate.var->locdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.aggregate.var, memhdr, set, stat, lp, branchcand, eventqueue, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar,
                        infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, var->locdom.ub) && SCIPsetIsInfinity(set, -var->data.aggregate.var->locdom.lb))
            || SCIPsetIsEQ(set, var->locdom.ub,
               var->data.aggregate.var->locdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.aggregate.var, memhdr, set, stat, lp, branchcand, eventqueue, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar, 
                        infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo change the sides of the corresponding linear constraint */
      errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgLbLocal(var->negatedvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                     var->data.negate.constant - newbound, 
                     infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype) );
      break;
         
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** changes current local bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
RETCODE SCIPvarChgBdLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place, or -1 */
   int              fixindex,           /**< bound change index for each node representing the order of changes, or -1 */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound, 
         infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound, 
         infervar, infercons, inferinfo, fixdepth, fixindex, boundchgtype);
   default:
      errorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }
}

/** sets current local lower bound of variable; if possible, adjusts bound to integral value; doesn't change the
 *  inference information stored in the variable
 */
RETCODE SCIPvarSetLbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   )
{
   CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound,
                  SCIPvarGetInferVar(var), SCIPvarGetInferCons(var), SCIPvarGetInferInfo(var),
                  SCIPvarGetFixDepth(var), SCIPvarGetFixIndex(var), SCIPvarGetBoundchgType(var)) );

   return SCIP_OKAY;
}

/** sets current local upper bound of variable; if possible, adjusts bound to integral value; doesn't change the
 *  inference information stored in the variable
 */
RETCODE SCIPvarSetUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   )
{
   CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound,
                  SCIPvarGetInferVar(var), SCIPvarGetInferCons(var), SCIPvarGetInferInfo(var),
                  SCIPvarGetFixDepth(var), SCIPvarGetFixIndex(var), SCIPvarGetBoundchgType(var)) );

   return SCIP_OKAY;
}

/** sets current local bound of variable; if possible, adjusts bound to integral value; doesn't change the
 *  inference information stored in the variable
 */
RETCODE SCIPvarSetBdLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarSetLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarSetUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound);
   default:
      errorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
   }
}

/** changes lower bound of variable in current dive; if possible, adjusts bound to integral value */
RETCODE SCIPvarChgLbDive(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(lp != NULL);
   assert(lp->diving);

   /* adjust bound for integral variables */
   SCIPvarAdjustLb(var, set, &newbound);

   debugMessage("changing lower bound of <%s> to %g in current dive\n", var->name, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      CHECK_OKAY( SCIPvarChgLbDive(var->data.transvar, set, lp, newbound) );
      break;
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      CHECK_OKAY( SCIPcolChgLb(var->data.col, set, lp, newbound) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      errorMessage("cannot change variable's bounds in dive for LOOSE variables\n");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         CHECK_OKAY( SCIPvarChgLbDive(var->data.aggregate.var, set, lp, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         CHECK_OKAY( SCIPvarChgUbDive(var->data.aggregate.var, set, lp, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
      
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo change the sides of the corresponding linear constraint */
      errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgUbDive(var->negatedvar, set, lp, var->data.negate.constant - newbound) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** changes upper bound of variable in current dive; if possible, adjusts bound to integral value */
RETCODE SCIPvarChgUbDive(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);
   assert(lp != NULL);
   assert(lp->diving);

   /* adjust bound for integral variables */
   SCIPvarAdjustUb(var, set, &newbound);

   debugMessage("changing upper bound of <%s> to %g in current dive\n", var->name, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      CHECK_OKAY( SCIPvarChgUbDive(var->data.transvar, set, lp, newbound) );
      break;
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      CHECK_OKAY( SCIPcolChgUb(var->data.col, set, lp, newbound) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      errorMessage("cannot change variable's bounds in dive for LOOSE variables\n");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable\n");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change upper bound of y */
         CHECK_OKAY( SCIPvarChgUbDive(var->data.aggregate.var, set, lp, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change lower bound of y */
         CHECK_OKAY( SCIPvarChgLbDive(var->data.aggregate.var, set, lp, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
      
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo change the sides of the corresponding linear constraint */
      errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgLbDive(var->negatedvar, set, lp, var->data.negate.constant - newbound) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** adds a hole to the variable's global domain and to its current local domain */
RETCODE SCIPvarAddHoleGlobal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(var != NULL);

   debugMessage("adding global hole (%g,%g) to <%s>\n", left, right, var->name);

   CHECK_OKAY( domAddHole(&var->glbdom, memhdr, set, left, right) );
   CHECK_OKAY( domAddHole(&var->locdom, memhdr, set, left, right) );

   /**@todo add hole in parent and child variables (just like with bound changes) */

   return SCIP_OKAY;
}

/** adds a hole to the variable's current local domain */
RETCODE SCIPvarAddHoleLocal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(var != NULL);

   debugMessage("adding local hole (%g,%g) to <%s>\n", left, right, var->name);

   CHECK_OKAY( domAddHole(&var->locdom, memhdr, set, left, right) );

   /**@todo add hole in parent and child variables (just like with bound changes) */

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z */
RETCODE SCIPvarAddVlb(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   Real             vlbconstant         /**< constant d    in x >= b*z + d */
   )
{
   assert(var != NULL);
   assert(vlbvar->vartype != SCIP_VARTYPE_CONTINUOUS);

   debugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
      SCIPvarGetName(var), vlbcoef, SCIPvarGetName(vlbvar), vlbconstant);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      CHECK_OKAY( SCIPvarAddVlb(var->data.transvar, memhdr, set, vlbvar, vlbcoef, vlbconstant) );
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      /* transform b*z + d into the corresponding sum after transforming z to an active problem variable */
      CHECK_OKAY( SCIPvarGetProbvarSum(&vlbvar, &vlbcoef, &vlbconstant) );
      debugMessage(" -> transformed to variable lower bound <%s> >= %g<%s> + %g\n", 
         SCIPvarGetName(var), vlbcoef, SCIPvarGetName(vlbvar), vlbconstant);

      if( vlbvar != NULL )
      {
         /* add variable bound to the variable bounds list */
         CHECK_OKAY( vboundsAdd(&var->vlbs, memhdr, set, vlbvar, vlbcoef, vlbconstant) );
      }
      break;
      
   case SCIP_VARSTATUS_FIXED:
      /* nothing to do here */
      break;
      
   case SCIP_VARSTATUS_AGGREGATED:
      /* x = a*y + c:  x >= b*z + d  <=>  a*y + c >= b*z + d  <=>  y >= b/a * z + (d-c)/a, if a > 0
       *                                                           y <= b/a * z + (d-c)/a, if a < 0
       */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> add variable lower bound */
         CHECK_OKAY( SCIPvarAddVlb(var->data.aggregate.var, memhdr, set, vlbvar,
                        vlbcoef/var->data.aggregate.scalar,
                        (vlbconstant - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> add variable upper bound */
         CHECK_OKAY( SCIPvarAddVub(var->data.aggregate.var, memhdr, set, vlbvar,
                        vlbcoef/var->data.aggregate.scalar,
                        (vlbconstant - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
      
   case SCIP_VARSTATUS_MULTAGGR:
      /* nothing to do here */
      break;
      
   case SCIP_VARSTATUS_NEGATED:
      /* x = offset - x':  x >= b*z + d  <=>  offset - x' >= b*z + d  <=>  x' <= -b*z + (offset-d) */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarAddVub(var->negatedvar, memhdr, set, vlbvar,
                     -vlbcoef, var->data.negate.constant - vlbconstant) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z */
RETCODE SCIPvarAddVub(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   Real             vubconstant         /**< constant d    in x <= b*z + d */
   )
{
   assert(var != NULL);
   assert(vubvar->vartype != SCIP_VARTYPE_CONTINUOUS);

   debugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
      SCIPvarGetName(var), vubcoef, SCIPvarGetName(vubvar), vubconstant);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      CHECK_OKAY( SCIPvarAddVub(var->data.transvar, memhdr, set, vubvar, vubcoef, vubconstant) );
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      /* transform b*z + d into the corresponding sum after transforming z to an active problem variable */
      CHECK_OKAY( SCIPvarGetProbvarSum(&vubvar, &vubcoef, &vubconstant) );
      debugMessage(" -> transformed to variable upper bound <%s> <= %g<%s> + %g\n", 
         SCIPvarGetName(var), vubcoef, SCIPvarGetName(vubvar), vubconstant);

      if( vubvar != NULL )
      {
         /* add variable bound to the variable bounds list */
         CHECK_OKAY( vboundsAdd(&var->vubs, memhdr, set, vubvar, vubcoef, vubconstant) );
      }
      break;
      
   case SCIP_VARSTATUS_FIXED:
      /* nothing to do here */
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      /* x = a*y + c:  x <= b*z + d  <=>  a*y + c <= b*z + d  <=>  y <= b/a * z + (d-c)/a, if a > 0
       *                                                           y >= b/a * z + (d-c)/a, if a < 0
       */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> add variable upper bound */
         CHECK_OKAY( SCIPvarAddVub(var->data.aggregate.var, memhdr, set, vubvar,
                        vubcoef/var->data.aggregate.scalar,
                        (vubconstant - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> add variable lower bound */
         CHECK_OKAY( SCIPvarAddVlb(var->data.aggregate.var, memhdr, set, vubvar,
                        vubcoef/var->data.aggregate.scalar,
                        (vubconstant - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         return SCIP_INVALIDDATA;
      }
      break;
      
   case SCIP_VARSTATUS_MULTAGGR:
      /* nothing to do here */
      break;
      
   case SCIP_VARSTATUS_NEGATED:
      /* x = offset - x':  x <= b*z + d  <=>  offset - x' <= b*z + d  <=>  x' >= -b*z + (offset-d) */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarAddVlb(var->negatedvar, memhdr, set, vubvar,
                     -vubcoef, var->data.negate.constant - vubconstant) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** replaces bounding variables in variable bounds of variable by their active problem variable counterparts */
RETCODE SCIPvarUseActiveVbds(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   CHECK_OKAY( vboundsUseActiveVars(var->vlbs) );
   CHECK_OKAY( vboundsUseActiveVars(var->vubs) );

   return SCIP_OKAY;
}

/** actually changes the branch factor of the variable and of all parent variables */
static
void varProcessChgBranchFactor(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   VAR* parentvar;
   int i;

   assert(var != NULL);
   assert(set != NULL);

   /* only use positive values */
   branchfactor = MAX(branchfactor, set->epsilon);

   debugMessage("process changing branch factor of <%s> from %f to %f\n", var->name, var->branchfactor, branchfactor);

   if( SCIPsetIsEQ(set, branchfactor, var->branchfactor) )
      return;

   /* change the branch factor */
   var->branchfactor = branchfactor;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change priorities across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_NEGATED:
         varProcessChgBranchFactor(parentvar, set, branchfactor);
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }
}

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
void SCIPvarChgBranchFactor(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   int v;

   assert(var != NULL);
   assert(branchfactor >= 0.0);

   debugMessage("changing branch factor of <%s> from %g to %g\n", var->name, var->branchfactor, branchfactor);

   if( SCIPsetIsEQ(set, var->branchfactor, branchfactor) )
      return;

   /* change priorities of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarChgBranchFactor(var->data.transvar, set, branchfactor);
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         var->branchfactor = branchfactor;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      varProcessChgBranchFactor(var, set, branchfactor);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      SCIPvarChgBranchFactor(var->data.aggregate.var, set, branchfactor);
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      for( v = 0; v < var->data.multaggr.nvars; ++v )
         SCIPvarChgBranchFactor(var->data.multaggr.vars[v], set, branchfactor);
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIPvarChgBranchFactor(var->negatedvar, set, branchfactor);
      break;
         
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** actually changes the branch priority of the variable and of all parent variables */
static
void varProcessChgBranchPriority(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   int              branchpriority      /**< branching priority of the variable */
   )
{
   VAR* parentvar;
   int i;

   assert(var != NULL);

   debugMessage("process changing branch priority of <%s> from %d to %d\n", 
      var->name, var->branchpriority, branchpriority);

   if( branchpriority == var->branchpriority )
      return;

   /* change the branch priority */
   var->branchpriority = branchpriority;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change priorities across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_NEGATED:
         varProcessChgBranchPriority(parentvar, set, branchpriority);
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }
}

/** sets the branch priority of the variable; variables with higher branch priority are always prefered to variables
 *  with lower priority in selection of branching variable
 */
void SCIPvarChgBranchPriority(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   int              branchpriority      /**< branching priority of the variable */
   )
{
   int v;

   assert(var != NULL);

   debugMessage("changing branch priority of <%s> from %d to %d\n", var->name, var->branchpriority, branchpriority);

   if( var->branchpriority == branchpriority )
      return;

   /* change priorities of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarChgBranchPriority(var->data.transvar, set, branchpriority);
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         var->branchpriority = branchpriority;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      varProcessChgBranchPriority(var, set, branchpriority);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      SCIPvarChgBranchPriority(var->data.aggregate.var, set, branchpriority);
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      for( v = 0; v < var->data.multaggr.nvars; ++v )
         SCIPvarChgBranchPriority(var->data.multaggr.vars[v], set, branchpriority);
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIPvarChgBranchPriority(var->negatedvar, set, branchpriority);
      break;
         
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** compares the index of two variables, returns -1 if first is smaller than, and +1 if first is greater than second
 *  variable index; returns 0 if both indices are equal, which means both variables are equal
 */
int SCIPvarCmp(
   VAR*             var1,               /**< first problem variable */
   VAR*             var2                /**< second problem variable */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);

   if( var1->index < var2->index )
      return -1;
   else if( var1->index > var2->index )
      return +1;
   else
   {
      assert(var1 == var2);
      return 0;
   }
}

/** gets corresponding active problem variable of a variable; returns NULL for fixed variables */
VAR* SCIPvarGetProbvar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   debugMessage("get problem variable of <%s>\n", var->name);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("original variable has no transformed variable attached\n");
         return NULL;
      }
      return SCIPvarGetProbvar(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var;

   case SCIP_VARSTATUS_FIXED:
      return NULL;

   case SCIP_VARSTATUS_AGGREGATED:
      return SCIPvarGetProbvar(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("multi-aggregated variable has no single active problem variable\n");
      abort();

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetProbvar(var->negatedvar);

   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** transforms given variable, boundtype and bound to the corresponding active variable values */
RETCODE SCIPvarGetProbvarBound(
   VAR**            var,                /**< pointer to problem variable */
   Real*            bound,              /**< pointer to bound value to transform */
   BOUNDTYPE*       boundtype           /**< pointer to type of bound: lower or upper bound */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert(bound != NULL);
   assert(boundtype != NULL);

   debugMessage("get probvar bound %g of type %d of variable <%s>\n", *bound, *boundtype, (*var)->name);

   switch( SCIPvarGetStatus(*var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( (*var)->data.transvar == NULL )
      {
         errorMessage("original variable has no transformed variable attached\n");
         return SCIP_INVALIDDATA;
      }
      *var = (*var)->data.transvar;
      CHECK_OKAY( SCIPvarGetProbvarBound(var, bound, boundtype) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("fixed variable has no corresponding active problem variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:  /* x = a*y + c  ->  y = x/a - c/a */
      assert((*var)->data.aggregate.var != NULL);
      assert((*var)->data.aggregate.scalar != 0.0);
      
      (*bound) /= (*var)->data.aggregate.scalar;
      (*bound) -= (*var)->data.aggregate.constant/(*var)->data.aggregate.scalar;
      if( (*var)->data.aggregate.scalar < 0.0 )
      {
         if( *boundtype == SCIP_BOUNDTYPE_LOWER )
            *boundtype = SCIP_BOUNDTYPE_UPPER;
         else
            *boundtype = SCIP_BOUNDTYPE_LOWER;
      }
      *var = (*var)->data.aggregate.var;
      CHECK_OKAY( SCIPvarGetProbvarBound(var, bound, boundtype) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("multiple aggregated variable has no single corresponding active problem variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert((*var)->negatedvar != NULL);
      assert(SCIPvarGetStatus((*var)->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert((*var)->negatedvar->negatedvar == *var);
      (*bound) = (*var)->data.negate.constant - *bound;
      if( *boundtype == SCIP_BOUNDTYPE_LOWER )
         *boundtype = SCIP_BOUNDTYPE_UPPER;
      else
         *boundtype = SCIP_BOUNDTYPE_LOWER;
      *var = (*var)->negatedvar;
      CHECK_OKAY( SCIPvarGetProbvarBound(var, bound, boundtype) );
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   return SCIP_OKAY;
}

/** transforms given variable, scalar and constant to the corresponding active variable, scalar and constant;
 *  if the variable resolves to a fixed variable, the returned variable will be NULL, "scalar" will be 0.0 and
 *  the value of the sum will be stored in "constant"
 */
RETCODE SCIPvarGetProbvarSum(
   VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   Real*            constant            /**< pointer to constant c in sum a*x + c */
   )
{
   assert(var != NULL);
   assert(scalar != NULL);
   assert(constant != NULL);

   while( *var != NULL )
   {
      switch( SCIPvarGetStatus(*var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( (*var)->data.transvar == NULL )
         {
            errorMessage("original variable has no transformed variable attached\n");
            return SCIP_INVALIDDATA;
         }
         *var = (*var)->data.transvar;
         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         return SCIP_OKAY;

      case SCIP_VARSTATUS_FIXED:       /* x = c'          =>  a*x + c ==             (a*c' + c) */
         (*constant) += *scalar * (*var)->glbdom.lb;
         *scalar = 0.0;
         *var = NULL;
         return SCIP_OKAY;

      case SCIP_VARSTATUS_AGGREGATED:  /* x = a'*x' + c'  =>  a*x + c == (a*a')*x' + (a*c' + c) */
         assert((*var)->data.aggregate.var != NULL);
         (*constant) += *scalar * (*var)->data.aggregate.constant;
         (*scalar) *= (*var)->data.aggregate.scalar;
         *var = (*var)->data.aggregate.var;
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("multiple aggregated variable has no single corresponding active problem variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_NEGATED:     /* x =  - x' + c'  =>  a*x + c ==   (-a)*x' + (a*c' + c) */
         assert((*var)->negatedvar != NULL);
         assert(SCIPvarGetStatus((*var)->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert((*var)->negatedvar->negatedvar == *var);
         (*constant) += *scalar * (*var)->data.negate.constant;
         (*scalar) *= -1.0;
         *var = (*var)->negatedvar;
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }
   *scalar = 0.0;

   return SCIP_OKAY;
}

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** get name of variable */
const char* SCIPvarGetName(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->name;
}

/** returns the user data of the variable */
VARDATA* SCIPvarGetData(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vardata;
}

/** gets status of variable */
VARSTATUS SCIPvarGetStatus(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (VARSTATUS)(var->varstatus);
}

/** returns whether the variable belongs to the original problem */
Bool SCIPvarIsOriginal(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED || var->negatedvar != NULL);

   return (SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED && SCIPvarGetStatus(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL));
}

/** returns whether the variable belongs to the transformed problem */
Bool SCIPvarIsTransformed(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED || var->negatedvar != NULL);

   return (SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL
      && (SCIPvarGetStatus(var) != SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_ORIGINAL));
}

/** returns whether the variable was created by negation of a different variable */
Bool SCIPvarIsNegated(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
}

/** gets type of variable */
VARTYPE SCIPvarGetType(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (VARTYPE)(var->vartype);
}

/** returns whether variable is of integral type (binary, integer, or implicit integer) */
Bool SCIPvarIsIntegral(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (var->vartype != SCIP_VARTYPE_CONTINUOUS);
}

/** returns whether variable's column should be present in the initial root LP */
Bool SCIPvarIsInitial(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->initial;
}

/** returns whether variable's column is removeable from the LP (due to aging or cleanup) */
Bool SCIPvarIsRemoveable(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->removeable;
}

/** returns whether variable is an active (neither fixed nor aggregated) variable */
Bool SCIPvarIsActive(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (var->probindex >= 0);
}

/** gets unique index of variable */
int SCIPvarGetIndex(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->index;
}

/** gets position of variable in problem, or -1 if variable is not active */
int SCIPvarGetProbindex(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->probindex;
}

/** gets transformed variable of ORIGINAL variable */
VAR* SCIPvarGetTransVar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   return var->data.transvar;
}

/** gets column of COLUMN variable */
COL* SCIPvarGetCol(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   return var->data.col;
}

/** gets aggregation variable y of an aggregated variable x = a*y + c */
VAR* SCIPvarGetAggrVar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.var;
}

/** gets aggregation scalar a of an aggregated variable x = a*y + c */
Real SCIPvarGetAggrScalar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.scalar;
}

/** gets aggregation constant c of an aggregated variable x = a*y + c */
Real SCIPvarGetAggrConstant(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.constant;
}

/** gets number n of aggregation variables of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
int SCIPvarGetMultaggrNVars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.nvars;
}

/** gets vector of aggregation variables y of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
VAR** SCIPvarGetMultaggrVars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.vars;
}

/** gets vector of aggregation scalars a of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
Real* SCIPvarGetMultaggrScalars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.scalars;
}

/** gets aggregation constant c of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
Real SCIPvarGetMultaggrConstant(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.constant;
}

/** gets the negation of the given variable; may return NULL, if no negation is existing yet */
VAR* SCIPvarGetNegatedVar(
   VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);

   return var->negatedvar;
}

/** gets the negation variable x of a negated variable x' = offset - x */
VAR* SCIPvarGetNegationVar(
   VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

   return var->negatedvar;
}

/** gets the negation offset of a negated variable x' = offset - x */
Real SCIPvarGetNegationConstant(
   VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

   return var->data.negate.constant;
}

/** gets objective function value of variable */
Real SCIPvarGetObj(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->obj;
}
   
/** gets global lower bound of variable */
Real SCIPvarGetLbGlobal(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->glbdom.lb;
}
   
/** gets global upper bound of variable */
Real SCIPvarGetUbGlobal(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->glbdom.ub;
}

/** gets current lower bound of variable */
Real SCIPvarGetLbLocal(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->locdom.lb;
}
   
/** gets current upper bound of variable */
Real SCIPvarGetUbLocal(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->locdom.ub;
}

/** gets inference variable of variable (variable that was assigned: parent of var, or var itself), or NULL */
VAR* SCIPvarGetInferVar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->infervar;
}

/** gets inference constraint of variable (constraint that deduced the current assignment), or NULL */
CONS* SCIPvarGetInferCons(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->infercons;
}

/** gets user information for inference to help resolving the conflict */
int SCIPvarGetInferInfo(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->inferinfo;
}

/** gets depth level, where the binary variable was fixed, or -1 if unfixed */
int SCIPvarGetFixDepth(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_FIXED || (var->fixdepth == 0 && var->fixindex == -1));

   return var->fixdepth;
}

/** gets fixing index of the variable in the depth level, where the binary variable was fixed, or -1 if unfixed or
 *  fixed during preprocessing
 */
int SCIPvarGetFixIndex(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_FIXED || (var->fixdepth == 0 && var->fixindex == -1));

   return var->fixindex;
}

/** returns TRUE iff first variable was fixed earlier than second variable */
Bool SCIPvarWasFixedEarlier(
   VAR*             var1,               /**< first problem variable */
   VAR*             var2                /**< second problem variable */
   )
{
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(var1->varstatus != SCIP_VARSTATUS_FIXED || (var1->fixdepth == 0 && var1->fixindex == -1));
   assert(var2->varstatus != SCIP_VARSTATUS_FIXED || (var2->fixdepth == 0 && var2->fixindex == -1));

   return (var1->fixdepth >= 0
      && (var2->fixdepth == -1
         || var1->fixdepth < var2->fixdepth
         || (var1->fixdepth == var2->fixdepth && var1->fixindex < var2->fixindex)));
}

/** gets type of bound change of fixed binary variable (fixed due to branching or due to inference) */
BOUNDCHGTYPE SCIPvarGetBoundchgType(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->boundchgtype;
}

/** gets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
Real SCIPvarGetBranchFactor(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->branchfactor;
}

/** gets the branch priority of the variable; variables with higher priority should always be prefered to variables
 *  with lower priority
 */
int SCIPvarGetBranchPriority(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->branchpriority;
}

/** gets number of variable lower bounds x >= b_i*z_i + d_i of given variable x */
int SCIPvarGetNVlbs(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vlbs != NULL ? var->vlbs->len : 0;
}

/** gets array with bounding variables z_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
VAR** SCIPvarGetVlbVars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vlbs != NULL ? var->vlbs->vars : NULL;
}

/** gets array with bounding coefficients b_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
Real* SCIPvarGetVlbCoefs(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vlbs != NULL ? var->vlbs->coefs : NULL;
}

/** gets array with bounding constants d_i in variable lower bounds x >= b_i*z_i + d_i of given variable x */
Real* SCIPvarGetVlbConstants(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vlbs != NULL ? var->vlbs->constants : NULL;
}

/** gets number of variable upper bounds x <= b_i*z_i + d_i of given variable x */
int SCIPvarGetNVubs(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vubs != NULL ? var->vubs->len : 0;
}

/** gets array with bounding variables z_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
VAR** SCIPvarGetVubVars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vubs != NULL ? var->vubs->vars : NULL;
}

/** gets array with bounding coefficients b_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
Real* SCIPvarGetVubCoefs(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vubs != NULL ? var->vubs->coefs : NULL;
}

/** gets array with bounding constants d_i in variable upper bounds x <= b_i*z_i + d_i of given variable x */
Real* SCIPvarGetVubConstants(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vubs != NULL ? var->vubs->constants : NULL;
}

#endif

/** gets objective value of variable in current dive */
Real SCIPvarGetObjDive(
   VAR*             var,                /**< problem variable */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetObjDive(var->data.transvar, set);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetObj(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->obj;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetObjDive(var->data.aggregate.var, set);
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot get the objective value of a multiple aggregated variable\n");
      abort();
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return -SCIPvarGetObjDive(var->negatedvar, set);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets lower bound of variable in current dive */
Real SCIPvarGetLbDive(
   VAR*             var,                /**< problem variable */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetLbDive(var->data.transvar, set);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetLb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->locdom.lb;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbDive(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbDive(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         abort();
      }
      
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo get the sides of the corresponding linear constraint */
      errorMessage("getting the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetUbDive(var->negatedvar, set);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets upper bound of variable in current dive */
Real SCIPvarGetUbDive(
   VAR*             var,                /**< problem variable */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetUbDive(var->data.transvar, set);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetUb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->locdom.ub;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbDive(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbDive(var->data.aggregate.var, set) + var->data.aggregate.constant;
      }
      else
      {
         errorMessage("scalar is zero in aggregation\n");
         abort();
      }
      
   case SCIP_VARSTATUS_MULTAGGR:
      /**@todo get the sides of the corresponding linear constraint */
      errorMessage("getting the bounds of a multiple aggregated variable is not implemented yet\n");
      abort();
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetLbDive(var->negatedvar, set);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets best local bound of variable with respect to the objective function */
Real SCIPvarGetBestBound(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->locdom.lb;
   else
      return var->locdom.ub;
}

/** gets worst local bound of variable with respect to the objective function */
Real SCIPvarGetWorstBound(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->locdom.ub;
   else
      return var->locdom.lb;
}

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
BOUNDTYPE SCIPvarGetBestBoundType(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return SCIP_BOUNDTYPE_LOWER;
   else
      return SCIP_BOUNDTYPE_UPPER;
}

/** gets type (lower or upper) of worst bound of variable with respect to the objective function */
BOUNDTYPE SCIPvarGetWorstBoundType(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return SCIP_BOUNDTYPE_UPPER;
   else
      return SCIP_BOUNDTYPE_LOWER;
}

/** gets primal LP solution value of variable */
Real SCIPvarGetLPSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real primsol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetLPSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
      return SCIPvarGetBestBound(var);

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetPrimsol(var->data.col);

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetLPSol(var->data.aggregate.var) + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      primsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         primsol += var->data.multaggr.scalars[i] * SCIPvarGetLPSol(var->data.multaggr.vars[i]);
      return primsol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetLPSol(var->negatedvar);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets pseudo solution value of variable at current node */
Real SCIPvarGetPseudoSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real pseudosol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetPseudoSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPvarGetBestBound(var);

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetPseudoSol(var->data.aggregate.var) + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      pseudosol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         pseudosol += var->data.multaggr.scalars[i] * SCIPvarGetPseudoSol(var->data.multaggr.vars[i]);
      return pseudosol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetPseudoSol(var->negatedvar);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets current LP or pseudo solution value of variable */
Real SCIPvarGetSol(
   VAR*             var,                /**< problem variable */
   Bool             getlpval            /**< should the LP solution value be returned? */
   )
{
   if( getlpval )
      return SCIPvarGetLPSol(var);
   else
      return SCIPvarGetPseudoSol(var);
}

/** remembers the current solution as root solution in the problem variables */
void SCIPvarStoreRootSol(
   VAR*             var,                /**< problem variable */
   Bool             roothaslp           /**< is the root solution from LP? */
   )
{
   assert(var != NULL);

   var->rootsol = SCIPvarGetSol(var, roothaslp);
}

/** returns the solution of the variable in the root node's relaxation, returns SCIP_INVALID if the root relaxation
 *  is not yet completely solved
 */
Real SCIPvarGetRootSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real rootsol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetRootSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var->rootsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetRootSol(var->data.aggregate.var) + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      rootsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         rootsol += var->data.multaggr.scalars[i] * SCIPvarGetRootSol(var->data.multaggr.vars[i]);
      return rootsol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetRootSol(var->negatedvar);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** resolves variable to columns and adds them with the coefficient to the row */
RETCODE SCIPvarAddToRow(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   int i;

   assert(var != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsInfinity(set, ABS(val)));

   debugMessage("adding coefficient %g<%s> to row <%s>\n", val, var->name, row->name);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("Cannot add untransformed original variable <%s> to LP row <%s>\n", var->name, row->name);
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarAddToRow(var->data.transvar, memhdr, set, stat, prob, lp, row, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
      /* convert loose variable into column */
      CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat, prob, lp) );
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      /*lint -fallthrough*/

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      assert(var->data.col->var == var);
      CHECK_OKAY( SCIProwIncCoef(row, memhdr, set, lp, var->data.col, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(var->glbdom.lb == var->glbdom.ub); /*lint !e777*/
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      assert(var->locdom.lb == var->glbdom.lb); /*lint !e777*/
      assert(!SCIPsetIsInfinity(set, ABS(var->locdom.lb)));
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, val * var->locdom.lb) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      CHECK_OKAY( SCIPvarAddToRow(var->data.aggregate.var, memhdr, set, stat, prob, lp,
                     row, var->data.aggregate.scalar * val) );
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, var->data.aggregate.constant * val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         CHECK_OKAY( SCIPvarAddToRow(var->data.multaggr.vars[i], memhdr, set, stat, prob, lp,
                        row, var->data.multaggr.scalars[i] * val) );
      }
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, var->data.multaggr.constant * val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarAddToRow(var->negatedvar, memhdr, set, stat, prob, lp, row, -val) );
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, var->data.negate.constant * val) );
      return SCIP_OKAY;

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** includes event handler with given data in variable's event filter */
RETCODE SCIPvarCatchEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert((eventtype & ~SCIP_EVENTTYPE_VARCHANGED) == 0);
   assert((eventtype & SCIP_EVENTTYPE_VARCHANGED) != 0);

   debugMessage("catch event of type 0x%x of variable <%s> with handler %p and data %p\n", 
      eventtype, var->name, eventhdlr, eventdata);

   CHECK_OKAY( SCIPeventfilterAdd(var->eventfilter, memhdr, set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}

/** deletes event handler with given data from variable's event filter */
RETCODE SCIPvarDropEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);

   debugMessage("drop event of variable <%s> with handler %p and data %p\n", var->name, eventhdlr, eventdata);

   CHECK_OKAY( SCIPeventfilterDel(var->eventfilter, memhdr, set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of
 *  "solvaldelta" in the variable's solution value and resulting change of "objdelta" in the in the LP's objective value
 */
RETCODE SCIPvarUpdatePseudocost(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   )
{
   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("cannot update pseudo costs of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarUpdatePseudocost(var->data.transvar, set, stat, solvaldelta, objdelta, weight) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      SCIPhistoryUpdatePseudocost(var->history, set, solvaldelta, objdelta, weight);
      SCIPhistoryUpdatePseudocost(stat->glbhistory, set, solvaldelta, objdelta, weight);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot update pseudo cost values of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      CHECK_OKAY( SCIPvarUpdatePseudocost(var->data.aggregate.var, set, stat,
                     solvaldelta/var->data.aggregate.scalar, objdelta, weight) );
      return SCIP_OKAY;
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot update pseudo cost values of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      CHECK_OKAY( SCIPvarUpdatePseudocost(var->negatedvar, set, stat, -solvaldelta, objdelta, weight) );
      return SCIP_OKAY;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value */
Real SCIPvarGetPseudocost(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   BRANCHDIR dir;

   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIPhistoryGetPseudocost(stat->glbhistory, solvaldelta);
      else
         return SCIPvarGetPseudocost(var->data.transvar, stat, solvaldelta);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      dir = (solvaldelta >= 0.0 ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);
      
      return SCIPhistoryGetPseudocostCount(var->history, dir) > 0.0
         ? SCIPhistoryGetPseudocost(var->history, solvaldelta)
         : SCIPhistoryGetPseudocost(stat->glbhistory, solvaldelta);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      return SCIPvarGetPseudocost(var->data.aggregate.var, stat, var->data.aggregate.scalar * solvaldelta);
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocost(var->negatedvar, stat, -solvaldelta);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
Real SCIPvarGetPseudocostCount(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction: 0 (down), or 1 (up) */
   )
{
   assert(var != NULL);
   assert(dir == 0 || dir == 1);
   
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetPseudocostCount(var->data.transvar, dir);
      
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetPseudocostCount(var->history, dir);
      
   case SCIP_VARSTATUS_FIXED:
      return 0.0;
      
   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetPseudocostCount(var->data.aggregate.var, dir);
      else
         return SCIPvarGetPseudocostCount(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostCount(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** increases the number of branchings counter of the variable */
RETCODE SCIPvarIncNBranchings(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   int              depth,              /**< depth at which the bound change took place */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("cannot update branching counter of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarIncNBranchings(var->data.transvar, stat, depth, dir) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      SCIPhistoryIncNBranchings(var->history, depth, dir);
      SCIPhistoryIncNBranchings(stat->glbhistory, depth, dir);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot update pseudo cost values of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
      {
         CHECK_OKAY( SCIPvarIncNBranchings(var->data.aggregate.var, stat, depth, dir) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         CHECK_OKAY( SCIPvarIncNBranchings(var->data.aggregate.var, stat, depth, SCIPbranchdirOpposite(dir)) );
      }
      return SCIP_OKAY;
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot update pseudo cost values of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      CHECK_OKAY( SCIPvarIncNBranchings(var->negatedvar, stat, depth, SCIPbranchdirOpposite(dir)) );
      return SCIP_OKAY;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases the number of inferences counter of the variable */
RETCODE SCIPvarIncNInferences(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("cannot update branching counter of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarIncNInferences(var->data.transvar, stat, dir) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      SCIPhistoryIncNInferences(var->history, dir);
      SCIPhistoryIncNInferences(stat->glbhistory, dir);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot update pseudo cost values of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
      {
         CHECK_OKAY( SCIPvarIncNInferences(var->data.aggregate.var, stat, dir) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         CHECK_OKAY( SCIPvarIncNInferences(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir)) );
      }
      return SCIP_OKAY;
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot update pseudo cost values of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      CHECK_OKAY( SCIPvarIncNInferences(var->negatedvar, stat, SCIPbranchdirOpposite(dir)) );
      return SCIP_OKAY;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases the number of cutoffs counter of the variable */
RETCODE SCIPvarIncNCutoffs(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("cannot update branching counter of original untransformed variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarIncNCutoffs(var->data.transvar, stat, dir) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      SCIPhistoryIncNCutoffs(var->history, dir);
      SCIPhistoryIncNCutoffs(stat->glbhistory, dir);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot update pseudo cost values of a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
      {
         CHECK_OKAY( SCIPvarIncNCutoffs(var->data.aggregate.var, stat, dir) );
      }
      else
      {
         assert(var->data.aggregate.scalar < 0.0);
         CHECK_OKAY( SCIPvarIncNCutoffs(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir)) );
      }
      return SCIP_OKAY;
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot update pseudo cost values of a multi-aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      CHECK_OKAY( SCIPvarIncNCutoffs(var->negatedvar, stat, SCIPbranchdirOpposite(dir)) );
      return SCIP_OKAY;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** returns the number of times, a bound of the variable was changed in given direction due to branching */
Longint SCIPvarGetNBranchings(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNBranchings(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNBranchings(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNBranchings(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNBranchings(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the average depth of bound changes in given direction due to branching on the variable */
Real SCIPvarGetAvgBranchdepth(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetAvgBranchdepth(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return  SCIPhistoryGetAvgBranchdepth(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgBranchdepth(var->data.aggregate.var, dir);
      else
         return SCIPvarGetAvgBranchdepth(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgBranchdepth(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the number of inferences branching on this variable in given direction triggered */
Longint SCIPvarGetNInferences(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNInferences(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNInferences(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNInferences(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNInferences(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNInferences(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the average number of inferences found after branching on the variable in given direction */
Real SCIPvarGetAvgInferences(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIPhistoryGetAvgInferences(stat->glbhistory, dir);
      else
         return SCIPvarGetAvgInferences(var->data.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->history, dir) > 0
         ? SCIPhistoryGetAvgInferences(var->history, dir)
         : SCIPhistoryGetAvgInferences(stat->glbhistory, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgInferences(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgInferences(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgInferences(var->negatedvar, stat, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the number of cutoffs branching on this variable in given direction produced */
Longint SCIPvarGetNCutoffs(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNCutoffs(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNCutoffs(var->history, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNCutoffs(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNCutoffs(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNCutoffs(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the average number of cutoffs found after branching on the variable in given direction */
Real SCIPvarGetAvgCutoffs(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction */
   )
{
   assert(var != NULL);
   assert(stat != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIPhistoryGetAvgCutoffs(stat->glbhistory, dir);
      else
         return SCIPvarGetAvgCutoffs(var->data.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->history, dir) > 0
         ? SCIPhistoryGetAvgCutoffs(var->history, dir)
         : SCIPhistoryGetAvgCutoffs(stat->glbhistory, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgCutoffs(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgCutoffs(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgCutoffs(var->negatedvar, stat, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
DECL_HASHGETKEY(SCIPhashGetKeyVar)
{  /*lint --e{715}*/
   VAR* var = (VAR*)elem;

   assert(var != NULL);
   return var->name;
}

