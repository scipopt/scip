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
#pragma ident "@(#) $Id: var.c,v 1.134 2005/01/25 09:59:32 bzfpfend Exp $"

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

/** ensures, that bound change info array for lower bound changes can store at least num entries */
static
RETCODE varEnsureLbchginfosSize(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(var != NULL);
   assert(var->nlbchginfos <= var->lbchginfossize);

   if( num > var->lbchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &var->lbchginfos, var->lbchginfossize, newsize) );
      var->lbchginfossize = newsize;
   }
   assert(num <= var->lbchginfossize);

   return SCIP_OKAY;
}

/** ensures, that bound change info array for upper bound changes can store at least num entries */
static
RETCODE varEnsureUbchginfosSize(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(var != NULL);
   assert(var->nubchginfos <= var->ubchginfossize);

   if( num > var->ubchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &var->ubchginfos, var->ubchginfossize, newsize) );
      var->ubchginfossize = newsize;
   }
   assert(num <= var->ubchginfossize);

   return SCIP_OKAY;
}

/** adds domain change info to the variable's lower bound change info array */
static
RETCODE varAddLbchginfo(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             oldbound,           /**< old value for bound */
   Real             newbound,           /**< new value for bound */
   int              depth,              /**< depth in the tree, where the bound change takes place */
   int              pos,                /**< position of the bound change in its bound change array */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   CONS*            infercons,          /**< constraint that infered this bound change, or NULL */
   PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   BOUNDTYPE        inferboundtype,     /**< type of bound for inference var: lower or upper bound */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type: branching decision or infered bound change */
   )
{
   assert(var != NULL);
   assert(SCIPsetIsLT(set, oldbound, newbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, oldbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, oldbound, 0.0));
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, newbound, 1.0));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || infervar != NULL);
   assert((boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER) == (infercons != NULL));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER || inferprop == NULL);

   debugMessage("adding lower bound change info to var <%s>[%g,%g]: depth=%d, pos=%d, %g -> %g\n",
      SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, depth, pos, oldbound, newbound);

   CHECK_OKAY( varEnsureLbchginfosSize(var, memhdr, set, var->nlbchginfos+1) );
   var->lbchginfos[var->nlbchginfos].oldbound = oldbound;
   var->lbchginfos[var->nlbchginfos].newbound = newbound;
   var->lbchginfos[var->nlbchginfos].var = var;
   var->lbchginfos[var->nlbchginfos].bdchgidx.depth = depth;
   var->lbchginfos[var->nlbchginfos].bdchgidx.pos = pos;
   var->lbchginfos[var->nlbchginfos].boundchgtype = boundchgtype;
   var->lbchginfos[var->nlbchginfos].boundtype = SCIP_BOUNDTYPE_LOWER;
   var->lbchginfos[var->nlbchginfos].inferboundtype = inferboundtype;
   var->lbchginfos[var->nlbchginfos].inferencedata.var = infervar;
   var->lbchginfos[var->nlbchginfos].inferencedata.info = inferinfo;

   switch( boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(infercons != NULL);
      var->lbchginfos[var->nlbchginfos].inferencedata.reason.cons = infercons;
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      var->lbchginfos[var->nlbchginfos].inferencedata.reason.prop = inferprop;
      break;
   default:
      errorMessage("invalid bound change type %d\n", boundchgtype);
      return SCIP_INVALIDDATA;
   }

   var->nlbchginfos++;

   assert(var->nlbchginfos < 2
      || SCIPbdchgidxIsEarlier(&var->lbchginfos[var->nlbchginfos-2].bdchgidx,
         &var->lbchginfos[var->nlbchginfos-1].bdchgidx));

   return SCIP_OKAY;
}

/** adds domain change info to the variable's upper bound change info array */
static
RETCODE varAddUbchginfo(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             oldbound,           /**< old value for bound */
   Real             newbound,           /**< new value for bound */
   int              depth,              /**< depth in the tree, where the bound change takes place */
   int              pos,                /**< position of the bound change in its bound change array */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   CONS*            infercons,          /**< constraint that infered this bound change, or NULL */
   PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   BOUNDTYPE        inferboundtype,     /**< type of bound for inference var: lower or upper bound */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type: branching decision or infered bound change */
   )
{
   assert(var != NULL);
   assert(SCIPsetIsGT(set, oldbound, newbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, oldbound));
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, oldbound, 1.0));
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, newbound, 0.0));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || infervar != NULL);
   assert((boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER) == (infercons != NULL));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER || inferprop == NULL);

   debugMessage("adding upper bound change info to var <%s>[%g,%g]: depth=%d, pos=%d, %g -> %g\n",
      SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, depth, pos, oldbound, newbound);

   CHECK_OKAY( varEnsureUbchginfosSize(var, memhdr, set, var->nubchginfos+1) );
   var->ubchginfos[var->nubchginfos].oldbound = oldbound;
   var->ubchginfos[var->nubchginfos].newbound = newbound;
   var->ubchginfos[var->nubchginfos].var = var;
   var->ubchginfos[var->nubchginfos].bdchgidx.depth = depth;
   var->ubchginfos[var->nubchginfos].bdchgidx.pos = pos;
   var->ubchginfos[var->nubchginfos].boundchgtype = boundchgtype;
   var->ubchginfos[var->nubchginfos].boundtype = SCIP_BOUNDTYPE_UPPER;
   var->ubchginfos[var->nubchginfos].inferboundtype = inferboundtype;
   var->ubchginfos[var->nubchginfos].inferencedata.var = infervar;
   var->ubchginfos[var->nubchginfos].inferencedata.info = inferinfo;

   switch( boundchgtype )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(infercons != NULL);
      var->ubchginfos[var->nubchginfos].inferencedata.reason.cons = infercons;
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      var->ubchginfos[var->nubchginfos].inferencedata.reason.prop = inferprop;
      break;
   default:
      errorMessage("invalid bound change type %d\n", boundchgtype);
      return SCIP_INVALIDDATA;
   }

   var->nubchginfos++;

   assert(var->nubchginfos < 2
      || SCIPbdchgidxIsEarlier(&var->ubchginfos[var->nubchginfos-2].bdchgidx,
         &var->ubchginfos[var->nubchginfos-1].bdchgidx));

   return SCIP_OKAY;
}

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
   int              pos,                /**< position of the bound change in its bound change array */
   Bool*            cutoff              /**< pointer to store whether an infeasible bound change was detected */
   )
{
   VAR* var;

   assert(boundchg != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(pos >= 0);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   var = boundchg->var;
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* apply bound change */
   switch( boundchg->boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      /* check, if the bound change is still active (could be replaced by inference due to repropagation of higher node) */
      if( SCIPsetIsGT(set, boundchg->newbound, var->locdom.lb) )
      {
         if( SCIPsetIsLE(set, boundchg->newbound, var->locdom.ub) )
         {
            /* add the bound change info to the variable's bound change info array */
            switch( boundchg->boundchgtype )
            {
            case SCIP_BOUNDCHGTYPE_BRANCHING:
               debugMessage(" -> branching: new lower bound of <%s>[%g,%g]: %g\n",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               CHECK_OKAY( varAddLbchginfo(var, memhdr, set, var->locdom.lb, boundchg->newbound, depth, pos,
                     NULL, NULL, NULL, 0, SCIP_BOUNDTYPE_LOWER, SCIP_BOUNDCHGTYPE_BRANCHING) );
               stat->lastbranchvar = var;
               stat->lastbranchdir = SCIP_BRANCHDIR_UPWARDS;
               break;

            case SCIP_BOUNDCHGTYPE_CONSINFER:
               assert(boundchg->data.inferencedata.reason.cons != NULL);
               debugMessage(" -> constraint <%s> inference: new lower bound of <%s>[%g,%g]: %g\n",
                  SCIPconsGetName(boundchg->data.inferencedata.reason.cons),
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               CHECK_OKAY( varAddLbchginfo(var, memhdr, set, var->locdom.lb, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, boundchg->data.inferencedata.reason.cons, NULL,
                     boundchg->data.inferencedata.info,
                     (BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_CONSINFER) );
               break;

            case SCIP_BOUNDCHGTYPE_PROPINFER:
               debugMessage(" -> propagator <%s> inference: new lower bound of <%s>[%g,%g]: %g\n",
                  boundchg->data.inferencedata.reason.prop != NULL
                  ? SCIPpropGetName(boundchg->data.inferencedata.reason.prop) : "-",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               CHECK_OKAY( varAddLbchginfo(var, memhdr, set, var->locdom.lb, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, NULL, boundchg->data.inferencedata.reason.prop,
                     boundchg->data.inferencedata.info,
                     (BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_PROPINFER) );
               break;

            default:
               errorMessage("invalid bound change type %d\n", boundchg->boundchgtype);
               return SCIP_INVALIDDATA;
            }
            
            /* change local bound of variable */
            CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->newbound) );
         }
         else
         {
            debugMessage(" -> cutoff: new lower bound of <%s>[%g,%g]: %g\n",
               SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
            boundchg->var = NULL;
            *cutoff = TRUE;
         }
      }
      else
      {
         /* mark bound change to be inactive */
         debugMessage(" -> inactive %s: new lower bound of <%s>[%g,%g]: %g\n",
            boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference",
            SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
         boundchg->var = NULL;
      }
      break;

   case SCIP_BOUNDTYPE_UPPER:
      /* check, if the bound change is still active (could be replaced by inference due to repropagation of higher node) */
      if( SCIPsetIsLT(set, boundchg->newbound, var->locdom.ub) )
      {
         if( SCIPsetIsGE(set, boundchg->newbound, var->locdom.lb) )
         {
            /* add the bound change info to the variable's bound change info array */
            switch( boundchg->boundchgtype )
            {
            case SCIP_BOUNDCHGTYPE_BRANCHING:
               debugMessage(" -> branching: new upper bound of <%s>[%g,%g]: %g\n",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               CHECK_OKAY( varAddUbchginfo(var, memhdr, set, var->locdom.ub, boundchg->newbound, depth, pos,
                     NULL, NULL, NULL, 0, SCIP_BOUNDTYPE_UPPER, SCIP_BOUNDCHGTYPE_BRANCHING) );
               stat->lastbranchvar = var;
               stat->lastbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
               break;

            case SCIP_BOUNDCHGTYPE_CONSINFER:
               assert(boundchg->data.inferencedata.reason.cons != NULL);
               debugMessage(" -> constraint <%s> inference: new upper bound of <%s>[%g,%g]: %g\n",
                  SCIPconsGetName(boundchg->data.inferencedata.reason.cons),
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               CHECK_OKAY( varAddUbchginfo(var, memhdr, set, var->locdom.ub, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, boundchg->data.inferencedata.reason.cons, NULL,
                     boundchg->data.inferencedata.info,
                     (BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_CONSINFER) );
               break;

            case SCIP_BOUNDCHGTYPE_PROPINFER:
               debugMessage(" -> propagator <%s> inference: new upper bound of <%s>[%g,%g]: %g\n",
                  boundchg->data.inferencedata.reason.prop != NULL
                  ? SCIPpropGetName(boundchg->data.inferencedata.reason.prop) : "-",
                  SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
               CHECK_OKAY( varAddUbchginfo(var, memhdr, set, var->locdom.ub, boundchg->newbound, depth, pos,
                     boundchg->data.inferencedata.var, NULL, boundchg->data.inferencedata.reason.prop,
                     boundchg->data.inferencedata.info,
                     (BOUNDTYPE)(boundchg->inferboundtype), SCIP_BOUNDCHGTYPE_PROPINFER) );
               break;

            default:
               errorMessage("invalid bound change type %d\n", boundchg->boundchgtype);
               return SCIP_INVALIDDATA;
            }

            /* change local bound of variable */
            CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, boundchg->newbound) );
         }
         else
         {
            debugMessage(" -> cutoff: new upper bound of <%s>[%g,%g]: %g\n",
               SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
            boundchg->var = NULL;
            *cutoff = TRUE;
         }
      }
      else
      {
         /* mark bound change to be inactive */
         debugMessage(" -> inactive %s: new upper bound of <%s>[%g,%g]: %g\n",
            boundchg->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference",
            SCIPvarGetName(var), var->locdom.lb, var->locdom.ub, boundchg->newbound);
         boundchg->var = NULL;
      }
      break;

   default:
      errorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
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
   VAR* var;

   assert(boundchg != NULL);
   assert(stat != NULL);

   var = boundchg->var;
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* undo bound change: apply the previous bound change of variable */
   switch( boundchg->boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      var->nlbchginfos--;
      assert(var->nlbchginfos >= 0);
      assert(var->lbchginfos != NULL);
      assert(var->lbchginfos[var->nlbchginfos].newbound == var->locdom.lb);
      assert(boundchg->newbound == var->locdom.lb);

      debugMessage("removed lower bound change info of var <%s>[%g,%g]: depth=%d, pos=%d, %g -> %g\n",
         SCIPvarGetName(var), var->locdom.lb, var->locdom.ub,
         var->lbchginfos[var->nlbchginfos].bdchgidx.depth, var->lbchginfos[var->nlbchginfos].bdchgidx.pos, 
         var->lbchginfos[var->nlbchginfos].oldbound, var->lbchginfos[var->nlbchginfos].newbound);

      /* reinstall the previous local bound */
      CHECK_OKAY( SCIPvarChgLbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue,
            var->lbchginfos[var->nlbchginfos].oldbound) );
      break;

   case SCIP_BOUNDTYPE_UPPER:
      var->nubchginfos--;
      assert(var->nubchginfos >= 0);
      assert(var->ubchginfos != NULL);
      assert(var->ubchginfos[var->nubchginfos].newbound == var->locdom.ub);
      assert(boundchg->newbound == var->locdom.ub);

      debugMessage("removed upper bound change info of var <%s>[%g,%g]: depth=%d, pos=%d, %g -> %g\n",
         SCIPvarGetName(var), var->locdom.lb, var->locdom.ub,
         var->ubchginfos[var->nubchginfos].bdchgidx.depth, var->ubchginfos[var->nubchginfos].bdchgidx.pos, 
         var->ubchginfos[var->nubchginfos].oldbound, var->ubchginfos[var->nubchginfos].newbound);

      /* reinstall the previous local bound */
      CHECK_OKAY( SCIPvarChgUbLocal(boundchg->var, memhdr, set, stat, lp, branchcand, eventqueue,
            var->ubchginfos[var->nubchginfos].oldbound) );
      break;

   default:
      errorMessage("unknown bound type\n");
      return SCIP_INVALIDDATA;
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
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      break;

   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(boundchg->data.inferencedata.var != NULL);
      assert(boundchg->data.inferencedata.reason.cons != NULL);
      SCIPconsCapture(boundchg->data.inferencedata.reason.cons);
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
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      break;

   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(boundchg->data.inferencedata.var != NULL);
      assert(boundchg->data.inferencedata.reason.cons != NULL);
      CHECK_OKAY( SCIPconsRelease(&boundchg->data.inferencedata.reason.cons, memhdr, set) );
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
         return SCIP_INVALIDDATA;
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
         return SCIP_INVALIDDATA;
      }
   }
#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < (*domchg)->domchgbound.nboundchgs; ++i )
         assert(SCIPvarGetType((*domchg)->domchgbound.boundchgs[i].var) == SCIP_VARTYPE_CONTINUOUS
            || EPSISINT((*domchg)->domchgbound.boundchgs[i].newbound, 1e-06));
   }
#endif

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
         return SCIP_INVALIDDATA;
      }
#ifndef NDEBUG
      if( *domchg != NULL )
      {
         int i;
         for( i = 0; i < (*domchg)->domchgbound.nboundchgs; ++i )
            assert(SCIPvarGetType((*domchg)->domchgbound.boundchgs[i].var) == SCIP_VARTYPE_CONTINUOUS
               || SCIPsetIsFeasIntegral(set, (*domchg)->domchgbound.boundchgs[i].newbound));
      }
#endif
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

#if 0
static
RETCODE domchgCleanupBdchgs(
   DOMCHG*          domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              firstremoved        /**< array position of first deactivated bound change */
   )
{
   VAR* var;
   int i;
   int j;
   int k;

   assert(domchg != NULL);
   assert(0 <= firstremoved && firstremoved < domchg->domchgbound.nboundchgs);
   assert(domchg->domchgbound.boundchgs[firstremoved].var == NULL);

#ifndef NDEBUG
   for( i = 0; i < firstremoved; ++i )
      assert(domchg->domchgbound.boundchgs[i].var != NULL);
#endif

   /* remove empty slots in bound change arrays by pushing the active bound changes to the front */
   i = firstremoved;
   for( j = i+1; j < domchg->domchgbound.nboundchgs; ++j )
   {
      var = domchg->domchgbound.boundchgs[j].var;
      if( var != NULL )
      {
         domchg->domchgbound.boundchgs[i] = domchg->domchgbound.boundchgs[j];
      }            
   }
}
#endif

/** shrinks bound change array of domain change to the given number of elements */
static
RETCODE domchgShrinkBdchgs(
   DOMCHG*          domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   int              newnbdchgs          /**< new number of bound changes */
   )
{
   assert(domchg != NULL);
   assert(newnbdchgs < (int)domchg->domchgbound.nboundchgs);

   if( domchg->domchgbound.domchgtype != SCIP_DOMCHGTYPE_DYNAMIC )
   {
      if( newnbdchgs == 0 )
      {
         freeBlockMemoryArray(memhdr, &domchg->domchgbound.boundchgs, domchg->domchgbound.nboundchgs);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &domchg->domchgbound.boundchgs, domchg->domchgbound.nboundchgs, 
               newnbdchgs) );
      }
   }
   domchg->domchgbound.nboundchgs = newnbdchgs;

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
   int              depth,              /**< depth in the tree, where the domain change takes place */
   Bool*            cutoff              /**< pointer to store whether an infeasible domain change was detected */
   )
{
   int i;
   int j;

   assert(cutoff != NULL);

   *cutoff = FALSE;

   debugMessage("applying domain changes at %p in depth %d\n", domchg, depth);
   if( domchg == NULL )
      return SCIP_OKAY;

   /* apply bound changes, removing inactive bound changes immediately */
   j = 0;
   for( i = 0; i < (int)domchg->domchgbound.nboundchgs; ++i )
   {
      assert(j <= i);
      CHECK_OKAY( SCIPboundchgApply(&domchg->domchgbound.boundchgs[i], memhdr, set, stat, lp,
            branchcand, eventqueue, depth, j, cutoff) );
      if( *cutoff )
         break;
      if( domchg->domchgbound.boundchgs[i].var != NULL )
      {
         if( j < i )
            domchg->domchgbound.boundchgs[j] = domchg->domchgbound.boundchgs[i];
         j++;
      }
      else
      {
         /* release inference data of inactive bound change */
         CHECK_OKAY( boundchgReleaseData(&domchg->domchgbound.boundchgs[i], memhdr, set) );
      }
   }
   debugMessage(" -> %d bound changes (%d active)\n", domchg->domchgbound.nboundchgs, j);

   /* release inference data of all bound changes from the point on where a cutoff was detected */
   for( ; i < (int)domchg->domchgbound.nboundchgs; ++i )
   {
      assert(*cutoff);
      CHECK_OKAY( boundchgReleaseData(&domchg->domchgbound.boundchgs[i], memhdr, set) );
   }

   /* clean up deactivated bound changes and all bound changes from the point on where a cutoff was detected */
   if( j < (int)domchg->domchgbound.nboundchgs )
   {
      CHECK_OKAY( domchgShrinkBdchgs(domchg, memhdr, j) );
   }
   assert((int)domchg->domchgbound.nboundchgs == j);

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
      CHECK_OKAY( SCIPboundchgUndo(&domchg->domchgbound.boundchgs[i], memhdr, set, stat, lp, branchcand, eventqueue) );
   }
   debugMessage(" -> %d bound changes\n", domchg->domchgbound.nboundchgs);

   return SCIP_OKAY;
}

/** adds bound change to domain changes */
RETCODE SCIPdomchgAddBoundchg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   BOUNDCHGTYPE     boundchgtype,       /**< type of bound change: branching decision or inference */
   Real             lpsolval,           /**< solval of variable in last LP on path to node, or SCIP_INVALID if unknown */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   BOUNDTYPE        inferboundtype      /**< type of bound for inference var: lower or upper bound */
   )
{
   BOUNDCHG* boundchg;

   assert(domchg != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsFeasIntegral(set, newbound));
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY
      || SCIPsetIsEQ(set, newbound, boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 : 0.0));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING || infervar != NULL);
   assert((boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER) == (infercons != NULL));
   assert(boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER || inferprop == NULL);

   debugMessage("adding %s bound change <%s: %g> of variable <%s> to domain change at %p pointing to %p\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", 
      boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING ? "branching" : "inference", 
      newbound, var->name, domchg, *domchg);

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
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      assert(infercons != NULL);
      boundchg->data.inferencedata.var = infervar;
      boundchg->data.inferencedata.reason.cons = infercons;
      boundchg->data.inferencedata.info = inferinfo; 
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      boundchg->data.inferencedata.var = infervar;
      boundchg->data.inferencedata.reason.prop = inferprop;
      boundchg->data.inferencedata.info = inferinfo; 
      break;
   default:
      errorMessage("invalid bound change type %d\n", boundchgtype);
      return SCIP_INVALIDDATA;
   }

   boundchg->newbound = newbound;
   boundchg->boundchgtype = boundchgtype; /*lint !e641*/
   boundchg->boundtype = boundtype; /*lint !e641*/
   boundchg->inferboundtype = inferboundtype;
   (*domchg)->domchgdyn.nboundchgs++;

   /* capture branching and inference data associated with the bound changes */
   CHECK_OKAY( boundchgCaptureData(boundchg) );

#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < (*domchg)->domchgbound.nboundchgs; ++i )
         assert(SCIPvarGetType((*domchg)->domchgbound.boundchgs[i].var) == SCIP_VARTYPE_CONTINUOUS
            || SCIPsetIsFeasIntegral(set, (*domchg)->domchgbound.boundchgs[i].newbound));
   }
#endif

   return SCIP_OKAY;
}

/** adds hole change to domain changes */
RETCODE SCIPdomchgAddHolechg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
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
   /**@todo don't add redundant variable bounds */
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
 * methods for implications
 */

/** creates an implications data structure */
static
RETCODE implicsCreate(
   IMPLICS**        implics,            /**< pointer to store implications data structure */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(implics != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, implics) );

   (*implics)->implvars[0] = NULL;
   (*implics)->impltypes[0] = NULL;
   (*implics)->implbounds[0] = NULL;
   (*implics)->implsize[0] = 0;
   (*implics)->nimpls[0] = 0;
   (*implics)->nbinimpls[0] = 0;

   (*implics)->implvars[1] = NULL;
   (*implics)->impltypes[1] = NULL;
   (*implics)->implbounds[1] = NULL;
   (*implics)->implsize[1] = 0;
   (*implics)->nimpls[1] = 0;
   (*implics)->nbinimpls[1] = 0;

   return SCIP_OKAY;
}

/** frees an implications data structure */
static
void implicsFree(
   IMPLICS**        implics,            /**< pointer of implications data structure to free */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(implics != NULL);

   if( *implics != NULL )
   {
      freeBlockMemoryArrayNull(memhdr, &(*implics)->implbounds[1], (*implics)->implsize[1]);
      freeBlockMemoryArrayNull(memhdr, &(*implics)->impltypes[1], (*implics)->implsize[1]);
      freeBlockMemoryArrayNull(memhdr, &(*implics)->implvars[1], (*implics)->implsize[1]);
      freeBlockMemoryArrayNull(memhdr, &(*implics)->implbounds[0], (*implics)->implsize[0]);
      freeBlockMemoryArrayNull(memhdr, &(*implics)->impltypes[0], (*implics)->implsize[0]);
      freeBlockMemoryArrayNull(memhdr, &(*implics)->implvars[0], (*implics)->implsize[0]);
      freeBlockMemory(memhdr, implics);
   }
}

/** ensures, that arrays for x <= 0 or x >= 1 in implications data structure can store at least num entries */
static
RETCODE implicsEnsureSize(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             i,                  /**< FALSE if size of arrays for x <= 0 has to be ensured, TRUE for x >= 1 */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(implics != NULL);
   
   /* create implications data structure, if not yet existing */
   if( *implics == NULL )
   {
      CHECK_OKAY( implicsCreate(implics, memhdr) );
   }
   assert(*implics != NULL);
   assert((*implics)->nimpls[i] <= (*implics)->implsize[i]);

   if( num > (*implics)->implsize[i] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*implics)->implvars[i], (*implics)->implsize[i], newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*implics)->impltypes[i], (*implics)->implsize[i], newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*implics)->implbounds[i], (*implics)->implsize[i], newsize) );
      (*implics)->implsize[i] = newsize;
   }
   assert(num <= (*implics)->implsize[i]);

   return SCIP_OKAY;
}

/** searches if variable y is allready contained in implications for x <= 0 or x >= 1
 *  y can be contained in structure with y >= b (y_lower) and y <= b (y_upper) 
 */
static
RETCODE implicsSearchVar(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   VAR*             implvar,            /**< variable y to search for */
   BOUNDTYPE        impltype,           /**< type of implication y <=/>= b to search for */
   Bool             i,                  /**< FALSE if y is searched in implications for x <= 0, TRUE for x >= 1 */
   int*             poslower,           /**< pointer to store position of y_lower (inf if not found) */
   int*             posupper,           /**< pointer to store position of y_upper (inf if not found) */
   int*             posadd              /**< pointer to store correct position (with respect to impltype) to add y */
   )
{
   int left;
   int right;
   int middle;

   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);

   /* set left and right pointer */
   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      /* no implications with binary variable y */
      if( (*implics)->nbinimpls[i] == 0 )
      {
         *posadd = 0;
         *poslower = INT_MAX;
         *posupper = INT_MAX;
          return SCIP_OKAY;
      }      
      left = 0;
      right = (*implics)->nbinimpls[i] - 1;
      assert(left <= right);

   }
   else
   {
      /* no implications with nonbinary variable y */
      if( (*implics)->nimpls[i] == (*implics)->nbinimpls[i] )
      {
         *posadd = (*implics)->nbinimpls[i];
         *poslower = INT_MAX;
         *posupper = INT_MAX;
         return SCIP_OKAY;
      }
      left = (*implics)->nbinimpls[i];
      right = (*implics)->nimpls[i] - 1;
      assert(left <= right);
   }

   /* searches for y */
   middle = (left + right) / 2;

   while( left <= right && (*implics)->implvars[i][middle] != implvar)
   {
      if( implvar < (*implics)->implvars[i][middle] ) 
         right = middle - 1;
      else
         left = middle + 1;
      middle = (left + right) / 2;
   }

   /* y was not found */
   if( left == right + 1 )
   {
      assert(right == -1 || (*implics)->implvars[i][right] < implvar);
      *poslower = INT_MAX;
      *posupper = INT_MAX;
      *posadd = left;
   }
   /* y was found */
   else
   {
      assert(implvar == (*implics)->implvars[i][middle]);

      /* sets poslower and posupper */
      /* y was found as y_lower (on position middle) */
      if( (*implics)->impltypes[i][middle] == SCIP_BOUNDTYPE_LOWER )
      {
         *poslower = middle;
         if( (*implics)->implvars[i][middle+1] == implvar )
         {  
            assert((*implics)->impltypes[i][middle+1] == SCIP_BOUNDTYPE_UPPER);
            *posupper = middle + 1;
         }
         else
            *posupper = INT_MAX;
      }
      /* y was found as y_upper (on position middle) */
      else
      {
         *posupper = middle;
         if( (*implics)->implvars[i][middle-1] == implvar )
         {  
            assert((*implics)->impltypes[i][middle-1] == SCIP_BOUNDTYPE_LOWER);
            *poslower = middle - 1;
         }
         else
            *poslower = INT_MAX;
      }
      /* sets posadd */
      if( impltype == SCIP_BOUNDTYPE_LOWER )
      {
         if( *poslower < INT_MAX )
            *posadd = *poslower;
         else
            *posadd = *posupper;
      }     
      else
      {
         if( *posupper < INT_MAX )
            *posadd = *posupper;
         else
            *posadd = *poslower;
      }
      assert(*posadd < INT_MAX);
   }

   return SCIP_OKAY;
}

/** adds an implication y <=/>= b for x <= 0 or x >= 1 (x binary) to the implications data structure */
static
RETCODE implicsAdd(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             i,                  /**< FALSE if implication for x <= 0 has to be added, TRUE for x >= 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   Bool*            conflict            /**< pointer to store whether implication causes a conflict for variable x */
   )
{
   int poslower;
   int posupper;
   int posadd;
   int k;

   assert(implics != NULL);
   assert(implvar != NULL);
   assert(conflict != NULL);
   assert(implvar->varstatus == SCIP_VARSTATUS_COLUMN || implvar->varstatus == SCIP_VARSTATUS_LOOSE); 
   assert(*implics == NULL || (*implics)->nbinimpls[i] <= (*implics)->nimpls[i]);
 
   *conflict = FALSE;

   /* don't add redundant implications with respect to global bounds:
    *    y >= b     if b <= lb (global lower bound of y) and
    *    y <= b     if b >= ub (global upper bound of y) */
   if( ( impltype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLE(set, implbound, SCIPvarGetLbGlobal(implvar)) ) 
      || ( impltype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGE(set, implbound, SCIPvarGetUbGlobal(implvar)) ) )
   {
      return SCIP_OKAY;
   }
   /* searches if variable is allready contained in implications data structure */
   if( *implics != NULL )
   {
      CHECK_OKAY( implicsSearchVar(implics, implvar, impltype, i, &poslower, &posupper, &posadd) );
      assert(poslower >= 0);
      assert(posupper >= 0);
      assert(posadd >= 0 && posadd <= (*implics)->nimpls[i]);
   }
   else
   {
      poslower = INT_MAX;
      posupper = INT_MAX;
      posadd = 0;
   }

   /* adds y >= b if not redundant and if it does not cause a conflict in x */
   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      /* y >= b is redundant */
      if( poslower < INT_MAX && SCIPsetIsLE(set, implbound,(*implics)->implbounds[i][poslower]) )
      {
         return SCIP_OKAY;
      }      

      /* y >= b causes conflict for x (i.e. y <= a (with a < b) is also valid) */
      if( posupper < INT_MAX && SCIPsetIsGT(set, implbound, (*implics)->implbounds[i][posupper]) )
      {      
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      /* adds y >= b by changing old entry on poslower */
      if( posadd == poslower )
      {
         assert(SCIPsetIsGT(set, implbound, (*implics)->implbounds[i][poslower]));
         (*implics)->implbounds[i][poslower] = implbound;

         return SCIP_OKAY;
      }
      
      /* adds y >= b by creating a new entry on posadd */
      assert(poslower == INT_MAX);

      CHECK_OKAY( implicsEnsureSize(implics, memhdr, set, i, *implics != NULL ? (*implics)->nimpls[i]+1 : 1) );
      assert(*implics != NULL);
      
      for( k = (*implics)->nimpls[i]; k > posadd; k-- )
      {
         assert((*implics)->implvars[i][k-1] >= implvar);
         (*implics)->implvars[i][k] = (*implics)->implvars[i][k-1];
         (*implics)->impltypes[i][k] = (*implics)->impltypes[i][k-1];
         (*implics)->implbounds[i][k] = (*implics)->implbounds[i][k-1];
      }
      assert(posadd == k);
      (*implics)->implvars[i][posadd] = implvar;
      (*implics)->impltypes[i][posadd] = impltype;
      (*implics)->implbounds[i][posadd] = implbound;
      if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
         (*implics)->nbinimpls[i]++;
      (*implics)->nimpls[i]++;
   }
   /* adds y <= b if not redundant and if it does not cause a conflict in x */
   else
   {
      /* y <= b is redundant */
      if( posupper < INT_MAX && SCIPsetIsGE(set, implbound,(*implics)->implbounds[i][posupper]) )
      {
         return SCIP_OKAY;
      }      

      /* y <= b causes conflict for x (i.e. y >= a (with a > b) is also valid) */
      if( poslower < INT_MAX && SCIPsetIsLT(set, implbound, (*implics)->implbounds[i][poslower]) )
      {      
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      /* adds y <= b by changing old entry on posupper */
      if( posadd == posupper )
      {
         assert(SCIPsetIsLT(set, implbound,(*implics)->implbounds[i][posupper]));
         (*implics)->implbounds[i][posupper] = implbound;

         return SCIP_OKAY;
      }
      
      /* adds y <= b by creating a new entry on posadd */
      assert(posupper == INT_MAX);

      CHECK_OKAY( implicsEnsureSize(implics, memhdr, set, i, *implics != NULL ? (*implics)->nimpls[i]+1 : 1) );
      assert(*implics != NULL);
      
      for( k = (*implics)->nimpls[i]; k > posadd; k-- )
      {
         assert((*implics)->implvars[i][k-1] >= implvar);
         (*implics)->implvars[i][k] = (*implics)->implvars[i][k-1];
         (*implics)->impltypes[i][k] = (*implics)->impltypes[i][k-1];
         (*implics)->implbounds[i][k] = (*implics)->implbounds[i][k-1];
      }
      assert(posadd == k);
      (*implics)->implvars[i][posadd] = implvar;
      (*implics)->impltypes[i][posadd] = impltype;
      (*implics)->implbounds[i][posadd] = implbound;
      if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
         (*implics)->nbinimpls[i]++;
      (*implics)->nimpls[i]++;
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
      return -SCIPsetInfinity(set);
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      return SCIPsetFeasCeil(set, lb);
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
      return SCIPsetInfinity(set);
   else if( vartype != SCIP_VARTYPE_CONTINUOUS )
      return SCIPsetFeasFloor(set, ub);
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

   assert(vartype != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, lb, 0.0) || SCIPsetIsEQ(set, lb, 1.0));
   assert(vartype != SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, ub, 0.0) || SCIPsetIsEQ(set, ub, 1.0));

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
   (*var)->obj = obj;
   (*var)->branchfactor = 1.0;
   (*var)->rootsol = 0.0;
   (*var)->primsolavg = 0.5 * (lb + ub);
   (*var)->glbdom.holelist = NULL;
   (*var)->glbdom.lb = lb;
   (*var)->glbdom.ub = ub;
   (*var)->locdom.holelist = NULL;
   (*var)->locdom.lb = lb;
   (*var)->locdom.ub = ub;
   (*var)->vardelorig = vardelorig;
   (*var)->vartrans = vartrans;
   (*var)->vardeltrans = vardeltrans;
   (*var)->vardata = vardata;
   (*var)->parentvars = NULL;
   (*var)->negatedvar = NULL;
   (*var)->vlbs = NULL;
   (*var)->vubs = NULL;
   (*var)->implics = NULL;
   (*var)->eventfilter = NULL;
   (*var)->lbchginfos = NULL;
   (*var)->ubchginfos = NULL;
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
   (*var)->branchpriority = 0;
   (*var)->branchdirection = SCIP_BRANCHDIR_AUTO;
   (*var)->lbchginfossize = 0;
   (*var)->nlbchginfos = 0;
   (*var)->ubchginfossize = 0;
   (*var)->nubchginfos = 0;
   (*var)->conflictsetcount = 0;
   (*var)->initial = initial;
   (*var)->removeable = removeable;
   (*var)->vartype = vartype; /*lint !e641*/
   (*var)->pseudocostflag = FALSE;

   /* create branching and inference history entries */
   CHECK_OKAY( SCIPhistoryCreate(&(*var)->history, memhdr) );
   CHECK_OKAY( SCIPhistoryCreate(&(*var)->historycrun, memhdr) );

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
      return SCIP_INVALIDDATA;
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

   /* free implications data structures */
   implicsFree(&(*var)->implics, memhdr);

   /* free bound change information arrays */
   freeBlockMemoryArrayNull(memhdr, &(*var)->lbchginfos, (*var)->lbchginfossize);
   freeBlockMemoryArrayNull(memhdr, &(*var)->ubchginfos, (*var)->ubchginfossize);

   /* free branching and inference history entries */
   SCIPhistoryFree(&(*var)->history, memhdr);
   SCIPhistoryFree(&(*var)->historycrun, memhdr);

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

/** initializes variable data structure for solving */
void SCIPvarInitSolve(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   SCIPhistoryReset(var->historycrun);
   var->conflictsetcount = 0;

   /* the negations of active problem variables can also be member of a conflict set */
   if( var->negatedvar != NULL )
      var->negatedvar->conflictsetcount = 0;
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

   /* type of variable */
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

   /* name */
   fprintf(file, " <%s>:", var->name);

   /* objective value */
   fprintf(file, " obj=%g", var->obj);

   /* bounds */
   fprintf(file, ", bounds=");
   if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
      fprintf(file, "[-inf,");
   else
      fprintf(file, "[%g,", var->glbdom.lb);
   if( SCIPsetIsInfinity(set, var->glbdom.ub) )
      fprintf(file, "+inf]");
   else
      fprintf(file, "%g]", var->glbdom.ub);

   /* holes */
   /**@todo print holes */

   /* fixings and aggregations */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      break;

   case SCIP_VARSTATUS_FIXED:
      fprintf(file, ", fixed: %g", var->glbdom.lb);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      fprintf(file, ", aggregated:");
      if( !SCIPsetIsZero(set, var->data.aggregate.constant) )
         fprintf(file, " %g", var->data.aggregate.constant);
      fprintf(file, " %+g<%s>", var->data.aggregate.scalar, SCIPvarGetName(var->data.aggregate.var));
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      fprintf(file, ", aggregated:");
      if( !SCIPsetIsZero(set, var->data.multaggr.constant) )
         fprintf(file, " %g", var->data.multaggr.constant);
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         fprintf(file, " %+g<%s>", var->data.multaggr.scalars[i], SCIPvarGetName(var->data.multaggr.vars[i]));
      break;

   case SCIP_VARSTATUS_NEGATED:
      fprintf(file, ", negated: %g - <%s>", var->data.negate.constant, SCIPvarGetName(var->negatedvar));
      break;

   default:
      errorMessage("unknown variable status\n");
      abort();
   }

   fprintf(file, "\n");
}

/** issues a LOCKSCHANGED event on the given variable */
static
RETCODE varEventLocksChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   EVENT* event;

   assert(var != NULL);

   /* issue LOCKSCHANGED event on variable */
   CHECK_OKAY( SCIPeventCreateLocksChanged(&event, memhdr, var) );
   CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, NULL, NULL, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** modifies lock numbers for rounding */
RETCODE SCIPvarAddLocks(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              addnlocksdown,      /**< increase in number of rounding down locks */
   int              addnlocksup         /**< increase in number of rounding up locks */
   )
{
   int i;

   assert(var != NULL);
   assert(var->nlocksup >= 0);
   assert(var->nlocksdown >= 0);

   if( addnlocksdown == 0 && addnlocksup == 0 )
      return SCIP_OKAY;

   debugMessage("add rounding locks %d/%d to variable <%s> (locks=%d/%d)\n",
      addnlocksdown, addnlocksup, var->name, var->nlocksdown, var->nlocksup);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarAddLocks(var->data.transvar, memhdr, set, eventqueue, addnlocksdown, addnlocksup) );
      }
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
      CHECK_OKAY( varEventLocksChanged(var, memhdr, set, eventqueue) );
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
      {
         CHECK_OKAY( SCIPvarAddLocks(var->data.aggregate.var, memhdr, set, eventqueue, addnlocksdown, addnlocksup) );
      }
      else
      {
         CHECK_OKAY( SCIPvarAddLocks(var->data.aggregate.var, memhdr, set, eventqueue, addnlocksup, addnlocksdown) );
      }
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
         {
            CHECK_OKAY( SCIPvarAddLocks(var->data.multaggr.vars[i], memhdr, set, eventqueue, 
                  addnlocksdown, addnlocksup) );
         }
         else
         {
            CHECK_OKAY( SCIPvarAddLocks(var->data.multaggr.vars[i], memhdr, set, eventqueue, 
                  addnlocksup, addnlocksdown) );
         }
      }
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarAddLocks(var->negatedvar, memhdr, set, eventqueue, addnlocksup, addnlocksdown) );
      break;

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   assert(var->nlocksdown >= 0);
   assert(var->nlocksup >= 0);

   return SCIP_OKAY;
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
            (VARTYPE)(origvar->vartype), origvar->initial, origvar->removeable,
            NULL, NULL, origvar->vardeltrans, origvar->vardata) );
      
      /* copy the branch factor and priority */
      (*transvar)->branchfactor = origvar->branchfactor;
      (*transvar)->branchpriority = origvar->branchpriority;
      (*transvar)->branchdirection = origvar->branchdirection;

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
   TREE*            tree,               /**< branch and bound tree */
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

   if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPsetIsFeasIntegral(set, fixedval))
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
         errorMessage("cannot fix an untransformed original variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarFix(var->data.transvar, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 
                     fixedval, infeasible) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* set the fixed variable's objective value to 0.0 */
      obj = var->obj;
      CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, primal, lp, eventqueue, 0.0) );

      /* change variable's bounds to fixed value */
      holelistFree(&var->glbdom.holelist, memhdr);
      CHECK_OKAY( SCIPvarChgLbGlobal(var, set, fixedval) );
      CHECK_OKAY( SCIPvarChgUbGlobal(var, set, fixedval) );
      holelistFree(&var->locdom.holelist, memhdr);
      CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, fixedval) );
      CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, fixedval) );

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
      SCIPhistoryReset(var->historycrun);

      /* convert variable into fixed variable */
      var->varstatus = SCIP_VARSTATUS_FIXED; /*lint !e641*/

      /* inform problem about the variable's status change */
      if( var->probindex != -1 )
      {
         CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, branchcand, var) );
      }

      /* reset the objective value of the fixed variable, thus adjusting the problem's objective offset */
      CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, stat, prob, primal, tree, lp, eventqueue, obj) );

      /* issue VARFIXED event */
      CHECK_OKAY( varEventVarFixed(var, memhdr, set, eventqueue) );
      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("cannot fix a column variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot fix a fixed variable again\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      /* fix aggregation variable y in x = a*y + c, instead of fixing x directly */
      assert(SCIPsetIsZero(set, var->obj));
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      CHECK_OKAY( SCIPvarFix(var->data.aggregate.var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue,
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
      CHECK_OKAY( SCIPvarFix(var->negatedvar, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue,
                     var->data.negate.constant - fixedval, infeasible) );
      break;

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
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
            varlb = -SCIPsetInfinity(set);
         else
            varlb = aggvar->glbdom.lb * scalar + constant;
         if( SCIPsetIsInfinity(set, aggvar->glbdom.ub) )
            varub = SCIPsetInfinity(set);
         else
            varub = aggvar->glbdom.ub * scalar + constant;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -aggvar->glbdom.lb) )
            varub = SCIPsetInfinity(set);
         else
            varub = aggvar->glbdom.lb * scalar + constant;
         if( SCIPsetIsInfinity(set, aggvar->glbdom.ub) )
            varlb = -SCIPsetInfinity(set);
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
         CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 
               varlb, infeasible) );
         if( !(*infeasible) )
         {
            CHECK_OKAY( SCIPvarFix(aggvar, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 
                  (varlb-constant)/scalar, infeasible) );
         }
         return SCIP_OKAY;
      }
      else
      {
         if( SCIPsetIsGT(set, varlb, var->glbdom.lb) )
         {
            CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, varlb) );
            CHECK_OKAY( SCIPvarChgLbGlobal(var, set, varlb) );
         }
         if( SCIPsetIsLT(set, varub, var->glbdom.ub) )
         {
            CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, varub) );
            CHECK_OKAY( SCIPvarChgUbGlobal(var, set, varub) );
         }

         /* update the hole list of the aggregation variable */
         /**@todo update hole list of aggregation variable */
      }

      /* update the bounds of the aggregation variable y in x = a*y + c  ->  y = (x-c)/a */
      if( scalar > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarlb = -SCIPsetInfinity(set);
         else
            aggvarlb = (var->glbdom.lb - constant) / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarub = SCIPsetInfinity(set);
         else
            aggvarub = (var->glbdom.ub - constant) / scalar;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarub = SCIPsetInfinity(set);
         else
            aggvarub = (var->glbdom.lb - constant) / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarlb = -SCIPsetInfinity(set);
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
         CHECK_OKAY( SCIPvarFix(aggvar, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, aggvarlb, 
               infeasible) );
         if( !(*infeasible) )
         {
            CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 
                  aggvarlb * scalar + constant, infeasible) );
         }
         return SCIP_OKAY;
      }
      else
      {
         if( SCIPsetIsGT(set, aggvarlb, aggvar->glbdom.lb) )
         {
            CHECK_OKAY( SCIPvarChgLbLocal(aggvar, memhdr, set, stat, lp, branchcand, eventqueue, aggvarlb) );
            CHECK_OKAY( SCIPvarChgLbGlobal(aggvar, set, aggvarlb) );
            aggvarbdschanged = TRUE;
         }
         if( SCIPsetIsLT(set, aggvarub, aggvar->glbdom.ub) )
         {
            CHECK_OKAY( SCIPvarChgUbLocal(aggvar, memhdr, set, stat, lp, branchcand, eventqueue, aggvarub) );
            CHECK_OKAY( SCIPvarChgUbGlobal(aggvar, set, aggvarub) );
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
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             aggvar,             /**< loose variable y in aggregation x = a*y + c */
   Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   BOUNDTYPE** impltypes;
   VAR*** implvars;
   VAR** vars;
   Real** implbounds;
   Real* coefs;
   Real* constants;
   Real obj;
   Real branchfactor;
   Bool implconflict[2];
   Bool conflict;
   int* nimpls;
   int branchpriority;
   int nlocksdown;
   int nlocksup;
   int nvbds;
   int i;
   int j;

   assert(var != NULL);
   assert(var->glbdom.lb == var->locdom.lb); /*lint !e777*/
   assert(var->glbdom.ub == var->locdom.ub); /*lint !e777*/
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */
   assert(infeasible != NULL);

   /* get active problem variable of aggregation variable */
   CHECK_OKAY( SCIPvarGetProbvarSum(&aggvar, &scalar, &constant) );
   assert((aggvar == NULL) == SCIPsetIsZero(set, scalar));

   /* aggregation is a fixing, if the scalar is zero */
   if( SCIPsetIsZero(set, scalar) )
      return SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, constant, infeasible);

   assert(aggvar != NULL);
   assert(aggvar->glbdom.lb == aggvar->locdom.lb); /*lint !e777*/
   assert(aggvar->glbdom.ub == aggvar->locdom.ub); /*lint !e777*/
   assert(SCIPvarGetStatus(aggvar) == SCIP_VARSTATUS_LOOSE);

   debugMessage("aggregate variable <%s>[%g,%g] == %g*<%s>[%g,%g] %+g\n", var->name, var->glbdom.lb, var->glbdom.ub,
      scalar, aggvar->name, aggvar->glbdom.lb, aggvar->glbdom.ub, constant);

   *infeasible = FALSE;

   /* if variable and aggregation variable are equal, the variable can be fixed: x == a*x + c  =>  x == c/(1-a) */
   if( var == aggvar )
   {
      if( SCIPsetIsEQ(set, scalar, 1.0) )
         *infeasible = !SCIPsetIsZero(set, constant);
      else
      {
         CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue,
               constant/(1.0-scalar), infeasible) );
      }
      return SCIP_OKAY;
   }

   /* tighten the bounds of aggregated and aggregation variable */
   CHECK_OKAY( varUpdateAggregationBounds(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue,
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

   /* check, if variable should be used as NEGATED variable of the aggregation variable */
   if( var->vartype == SCIP_VARTYPE_BINARY && aggvar->vartype == SCIP_VARTYPE_BINARY
      && var->negatedvar == NULL && aggvar->negatedvar == NULL
      && SCIPsetIsEQ(set, scalar, -1.0) && SCIPsetIsEQ(set, constant, 1.0) )
   {
      /* link both variables as negation pair */
      var->varstatus = SCIP_VARSTATUS_NEGATED;
      var->data.negate.constant = 1.0;
      var->negatedvar = aggvar;
      aggvar->negatedvar = var;
   }
   else
   {
      /* convert variable into aggregated variable */
      var->varstatus = SCIP_VARSTATUS_AGGREGATED; /*lint !e641*/
      var->data.aggregate.var = aggvar;
      var->data.aggregate.scalar = scalar;
      var->data.aggregate.constant = constant;
   }

   /* make aggregated variable a parent of the aggregation variable */
   CHECK_OKAY( varAddParent(aggvar, memhdr, set, var) );

   /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variable */
   CHECK_OKAY( SCIPvarAddLocks(var, memhdr, set, eventqueue, nlocksdown, nlocksup) );

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


   /* move the implications to the aggregation variable:
    *  - add all implications again to the variable, thus adding it to the aggregated variable
    *  - free the implications data structures
    */
   implconflict[0] = FALSE;
   implconflict[1] = FALSE;

   if( var->implics != NULL  && SCIPvarGetType(aggvar) == SCIP_VARTYPE_BINARY )
   {
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      nimpls = var->implics->nimpls;
      implvars = var->implics->implvars;
      impltypes = var->implics->impltypes;
      implbounds = var->implics->implbounds;
      for( i = 0; i < 2; ++i )
      {
         for( j = 0; j < nimpls[i]; ++j )
         {
            CHECK_OKAY( SCIPvarAddImplic(var, memhdr, set, i, implvars[i][j], impltypes[i][j], implbounds[i][j],
                  &conflict) );
            implconflict[i] = implconflict[i] || conflict;
         }
      }
   }
   implicsFree(&var->implics, memhdr);

   /* add the history entries to the aggregation variable and clear the history of the aggregated variable */
   SCIPhistoryUnite(aggvar->history, var->history, scalar < 0.0);
   SCIPhistoryUnite(aggvar->historycrun, var->historycrun, scalar < 0.0);
   SCIPhistoryReset(var->history);
   SCIPhistoryReset(var->historycrun);

   /* update flags of aggregation variable */
   aggvar->removeable &= var->removeable;

   /* update branching factors and priorities of both variables to be the maximum of both variables */
   branchfactor = MAX(aggvar->branchfactor, var->branchfactor);
   branchpriority = MAX(aggvar->branchpriority, var->branchpriority);
   SCIPvarChgBranchFactor(aggvar, set, branchfactor);
   SCIPvarChgBranchPriority(aggvar, branchpriority);
   SCIPvarChgBranchFactor(var, set, branchfactor);
   SCIPvarChgBranchPriority(var, branchpriority);

   /* update branching direction of both variables to agree to a single direction */
   if( scalar >= 0.0 )
   {
      if( var->branchdirection == SCIP_BRANCHDIR_AUTO )
         SCIPvarChgBranchDirection(var, aggvar->branchdirection);
      else if( aggvar->branchdirection == SCIP_BRANCHDIR_AUTO )
         SCIPvarChgBranchDirection(aggvar, var->branchdirection);
      else if( var->branchdirection != aggvar->branchdirection )
         SCIPvarChgBranchDirection(var, SCIP_BRANCHDIR_AUTO);
   }
   else
   {
      if( var->branchdirection == SCIP_BRANCHDIR_AUTO )
         SCIPvarChgBranchDirection(var, SCIPbranchdirOpposite(aggvar->branchdirection));
      else if( aggvar->branchdirection == SCIP_BRANCHDIR_AUTO )
         SCIPvarChgBranchDirection(aggvar, SCIPbranchdirOpposite(var->branchdirection));
      else if( var->branchdirection != aggvar->branchdirection )
         SCIPvarChgBranchDirection(var, SCIP_BRANCHDIR_AUTO);
   }

   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      CHECK_OKAY( SCIPprobVarChangedStatus(prob, set, branchcand, var) );
   }

   /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
    * variable and the problem's objective offset
    */
   CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, stat, prob, primal, tree, lp, eventqueue, obj) );

   /* issue VARFIXED event */
   CHECK_OKAY( varEventVarFixed(var, memhdr, set, eventqueue) );

   /* fix variable if the implications detected a conflict */
   if( !(*infeasible) )
   {
      if( implconflict[0] && implconflict[1] )
         *infeasible = TRUE;
      else if( implconflict[0] ) 
      {
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 1.0, infeasible) );
      }
      else if( implconflict[1] ) 
      {
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         CHECK_OKAY( SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 0.0, infeasible) );
      }
   }

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
   TREE*            tree,               /**< branch and bound tree */
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
   BRANCHDIR branchdirection;
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
      return SCIPvarFix(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, constant, infeasible);
   else if( naggvars == 1)
      return SCIPvarAggregate(var, memhdr, set, stat, prob, primal, tree, lp, branchcand, eventqueue, 
         aggvars[0], scalars[0], constant, infeasible);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("cannot multi-aggregate an untransformed original variable\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarMultiaggregate(var->data.transvar, memhdr, set, stat, prob, primal, tree, lp,
            branchcand, eventqueue, naggvars, aggvars, scalars, constant, infeasible) );
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
      CHECK_OKAY( SCIPvarAddLocks(var, memhdr, set, eventqueue, nlocksdown, nlocksup) );

      /* free the variable bounds data structures */
      vboundsFree(&var->vlbs, memhdr);
      vboundsFree(&var->vubs, memhdr);

      /* free the implications data structures */
      implicsFree(&var->implics, memhdr);

      /* update flags and branching factors and priorities of aggregation variables;
       * update preferred branching direction of all aggregation variables that don't have a preferred direction yet
       */
      branchfactor = var->branchfactor;
      branchpriority = var->branchpriority;
      branchdirection = var->branchdirection;
      for( v = 0; v < naggvars; ++v )
      {
         assert(aggvars[v] != NULL);
         aggvars[v]->removeable &= var->removeable;
         branchfactor = MAX(aggvars[v]->branchfactor, branchfactor);
         branchpriority = MAX(aggvars[v]->branchpriority, branchpriority);
      }
      for( v = 0; v < naggvars; ++v )
      {
         SCIPvarChgBranchFactor(aggvars[v], set, branchfactor);
         SCIPvarChgBranchPriority(aggvars[v], branchpriority);
         if( aggvars[v]->branchdirection == SCIP_BRANCHDIR_AUTO )
         {
            if( scalars[v] >= 0.0 )
               SCIPvarChgBranchDirection(aggvars[v], branchdirection);
            else
               SCIPvarChgBranchDirection(aggvars[v], SCIPbranchdirOpposite(branchdirection));
         }
      }
      SCIPvarChgBranchFactor(var, set, branchfactor);
      SCIPvarChgBranchPriority(var, branchpriority);

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
      CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, stat, prob, primal, tree, lp, eventqueue, obj) );

      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("cannot multi-aggregate a column variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot multi-aggregate a fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      errorMessage("cannot multi-aggregate an aggregated variable\n");
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
      CHECK_OKAY( SCIPvarMultiaggregate(var->negatedvar, memhdr, set, stat, prob, primal, tree, lp,
            branchcand, eventqueue, naggvars, aggvars, scalars, var->data.negate.constant - constant, infeasible) );

      /* switch the signs of the aggregation scalars again, to reset them to their original values */
      for( v = 0; v < naggvars; ++v )
         scalars[v] *= -1.0;
      break;

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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
      
      /* create event filter, if the negated variable belongs to the transformed problem */
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         CHECK_OKAY( SCIPeventfilterCreate(&(*negvar)->eventfilter, memhdr) );
      }

      /* link the variables together */
      var->negatedvar = *negvar;
      (*negvar)->negatedvar = var;

      /* copy the branch factor and priority, and use the negative preferred branching direction */
      (*negvar)->branchfactor = var->branchfactor;
      (*negvar)->branchpriority = var->branchpriority;
      (*negvar)->branchdirection = SCIPbranchdirOpposite(var->branchdirection);

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

/** informs variable that its position in problem's vars array changed */
void SCIPvarSetProbindex(
   VAR*             var,                /**< problem variable */
   int              probindex           /**< new problem index of variable */
   )
{
   assert(var != NULL);

   var->probindex = probindex;
   if( var->varstatus == SCIP_VARSTATUS_COLUMN )
   {
      assert(var->data.col != NULL);
      var->data.col->var_probindex = probindex;
   }
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
            assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
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
         return SCIP_INVALIDDATA;
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
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             addobj              /**< additional objective value for variable */
   )
{
   assert(var != NULL);
   assert(set != NULL);
   assert(SCIPgetStage(set->scip) < SCIP_STAGE_INITSOLVE);

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
            CHECK_OKAY( SCIPvarAddObj(var->data.transvar, memhdr, set, stat, prob, primal, tree, lp, eventqueue, addobj) );
         }
         else
         {
            assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
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
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, tree, lp) );
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         /* x = a*y + c  ->  add a*addobj to obj. val. of y, and c*addobj to obj. offset of problem */
         SCIPprobAddObjoffset(prob, var->data.aggregate.constant * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, tree, lp) );
         CHECK_OKAY( SCIPvarAddObj(var->data.aggregate.var, memhdr, set, stat, prob, primal, tree, lp, eventqueue,
                        var->data.aggregate.scalar * addobj) );
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         /* x = a_1*y_1 + ... + a_n*y_n  + c  ->  add a_i*addobj to obj. val. of y_i, and c*addobj to obj. offset */
         SCIPprobAddObjoffset(prob, var->data.multaggr.constant * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, tree, lp) );
         for( i = 0; i < var->data.multaggr.nvars; ++i )
         {
            CHECK_OKAY( SCIPvarAddObj(var->data.multaggr.vars[i], memhdr, set, stat, prob, primal, tree, lp, 
                           eventqueue, var->data.multaggr.scalars[i] * addobj) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED:
         /* x' = offset - x  ->  add -addobj to obj. val. of x and offset*addobj to obj. offset of problem */
         assert(var->negatedvar != NULL);
         assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(var->negatedvar->negatedvar == var);
         SCIPprobAddObjoffset(prob, var->data.negate.constant * addobj);
         CHECK_OKAY( SCIPprimalUpdateUpperbound(primal, memhdr, set, stat, prob, tree, lp) );
         CHECK_OKAY( SCIPvarAddObj(var->negatedvar, memhdr, set, stat, prob, primal, tree, lp, eventqueue, -addobj) );
         break;

      default:
         errorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
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
   SCIPlpMarkDivingObjChanged(lp);

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
      return SCIP_INVALIDDATA;
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

   *lb = adjustedLb(set, (VARTYPE)(var->vartype), *lb);
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

   *ub = adjustedUb(set, (VARTYPE)(var->vartype), *ub);
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
   assert(SCIPsetIsEQ(set, newbound, adjustedLb(set, (VARTYPE)(var->vartype), newbound)));

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
         return SCIP_INVALIDDATA;
      
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
         return SCIP_INVALIDDATA;
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
   assert(SCIPsetIsEQ(set, newbound, adjustedUb(set, (VARTYPE)(var->vartype), newbound)));

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
         /* do not change bounds across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;
      
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
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** changes global lower bound of variable; if possible, adjusts bound to integral value */
RETCODE SCIPvarChgLbGlobal(
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
         CHECK_OKAY( SCIPvarChgLbGlobal(var->data.transvar, set, newbound) );
      }
      else
      {
         assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
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
         CHECK_OKAY( SCIPvarChgLbGlobal(var->data.aggregate.var, set,
               (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, -var->glbdom.lb) && SCIPsetIsInfinity(set, var->data.aggregate.var->glbdom.ub))
            || SCIPsetIsEQ(set, var->glbdom.lb,
               var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbGlobal(var->data.aggregate.var, set,
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
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgUbGlobal(var->negatedvar, set,  var->data.negate.constant - newbound) );
      break;

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes global upper bound of variable; if possible, adjusts bound to integral value */
RETCODE SCIPvarChgUbGlobal(
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
         CHECK_OKAY( SCIPvarChgUbGlobal(var->data.transvar, set, newbound) );
      }
      else
      {
         assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
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
         CHECK_OKAY( SCIPvarChgUbGlobal(var->data.aggregate.var, set,
               (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, var->glbdom.ub) && SCIPsetIsInfinity(set, -var->data.aggregate.var->glbdom.lb))
            || SCIPsetIsEQ(set, var->glbdom.ub,
               var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbGlobal(var->data.aggregate.var, set,
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
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgLbGlobal(var->negatedvar, set,  var->data.negate.constant - newbound) );
      break;

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** changes global bound of variable; if possible, adjusts bound to integral value */
RETCODE SCIPvarChgBdGlobal(
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
      return SCIPvarChgLbGlobal(var, set, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUbGlobal(var, set, newbound);
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
   Real             newbound            /**< new bound for variable */
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
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(SCIPsetIsEQ(set, newbound, adjustedLb(set, (VARTYPE)(var->vartype), newbound)));

   debugMessage("process changing lower bound of <%s> from %g to %g\n", var->name, var->locdom.lb, newbound);

   if( SCIPsetIsEQ(set, newbound, var->locdom.lb) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->locdom.lb;
   var->locdom.lb = newbound;

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
         return SCIP_INVALIDDATA;
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            assert((SCIPsetIsInfinity(set, -parentvar->locdom.lb) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                  parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else 
         {
            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, parentvar->locdom.ub) && SCIPsetIsInfinity(set, -oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                  parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue,
               parentvar->data.negate.constant - newbound) );
         break;

      default:
         errorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
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
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(SCIPsetIsEQ(set, newbound, adjustedUb(set, (VARTYPE)(var->vartype), newbound)));

   debugMessage("process changing upper bound of <%s> from %g to %g\n", var->name, var->locdom.ub, newbound);

   if( SCIPsetIsEQ(set, newbound, var->locdom.ub) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->locdom.ub;
   var->locdom.ub = newbound;

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
         return SCIP_INVALIDDATA;
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of x */
            assert((SCIPsetIsInfinity(set, parentvar->locdom.ub) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.ub,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                  parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else
         {
            /* a < 0 -> change lower bound of x */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert((SCIPsetIsInfinity(set, -parentvar->locdom.lb) && SCIPsetIsInfinity(set, oldbound))
               || SCIPsetIsEQ(set, parentvar->locdom.lb,
                  oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue, 
                  parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = offset - x'  ->  x' = offset - x */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, lp, branchcand, eventqueue,
               parentvar->data.negate.constant - newbound) );
         break;

      default:
         errorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
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
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   /* adjust bound for integral variables */
   SCIPvarAdjustLb(var, set, &newbound);
   assert(SCIPsetIsLE(set, newbound, var->locdom.ub));

   debugMessage("changing lower bound of <%s>[%g,%g] to %g\n", var->name, var->locdom.lb, var->locdom.ub, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->locdom.lb, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.transvar, memhdr, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      else
      {
         assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
         stat->domchgcount++;
         var->locdom.lb = newbound;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      stat->domchgcount++;
      CHECK_OKAY( varProcessChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound) );
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
               (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert((SCIPsetIsInfinity(set, -var->locdom.lb) && SCIPsetIsInfinity(set, var->data.aggregate.var->locdom.ub))
            || SCIPsetIsEQ(set, var->locdom.lb,
               var->data.aggregate.var->locdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.aggregate.var, memhdr, set, stat, lp, branchcand, eventqueue, 
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
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgUbLocal(var->negatedvar, memhdr, set, stat, lp, branchcand, eventqueue, 
            var->data.negate.constant - newbound) );
      break;
         
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   /* adjust bound for integral variables */
   SCIPvarAdjustUb(var, set, &newbound);
   assert(SCIPsetIsGE(set, newbound, var->locdom.lb));

   debugMessage("changing upper bound of <%s>[%g,%g] to %g\n", var->name, var->locdom.lb, var->locdom.ub, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;
   
   if( SCIPsetIsEQ(set, var->locdom.ub, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.transvar, memhdr, set, stat, lp, branchcand, eventqueue, newbound) );
      }
      else
      {
         assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
         stat->domchgcount++;
         var->locdom.ub = newbound;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      stat->domchgcount++;
      CHECK_OKAY( varProcessChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound) );
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
               (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change lower bound of y */
         assert((SCIPsetIsInfinity(set, var->locdom.ub) && SCIPsetIsInfinity(set, -var->data.aggregate.var->locdom.lb))
            || SCIPsetIsEQ(set, var->locdom.ub,
               var->data.aggregate.var->locdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.aggregate.var, memhdr, set, stat, lp, branchcand, eventqueue, 
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
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgLbLocal(var->negatedvar, memhdr, set, stat, lp, branchcand, eventqueue, 
            var->data.negate.constant - newbound) );
      break;
         
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, newbound);
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
   assert(SCIPlpDiving(lp));

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
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgUbDive(var->negatedvar, set, lp, var->data.negate.constant - newbound) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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
   assert(SCIPlpDiving(lp));

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
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      CHECK_OKAY( SCIPvarChgLbDive(var->negatedvar, set, lp, var->data.negate.constant - newbound) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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
      return SCIP_INVALIDDATA;
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
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** informs binary variable x about a globally valid implication:  x <= 0 or x >= 1  ==>  y <= b  or  y >= b */
RETCODE SCIPvarAddImplic(
   VAR*             var,                /**< problem variable  */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             i,                  /**< FALSE if y should be added in implications for x <= 0, TRUE for x >= 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   Bool*            conflict            /**< pointer to store whether fixing variable to given value results in conflict */ 
   )
{
   assert(var != NULL);
   assert(conflict != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY); 
   
   *conflict = FALSE;

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      CHECK_OKAY( SCIPvarAddImplic(var->data.transvar, memhdr, set, i, implvar, impltype, implbound, conflict) );
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      CHECK_OKAY( SCIPvarGetProbvarBound(&implvar, &implbound, &impltype) );
         
      /* add implication to the implications list */
      if( (implbound == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLT(set, implvar->glbdom.ub, implbound))
         || (implbound == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGT(set, implvar->glbdom.lb, implbound)) )
         *conflict = TRUE;
      else
      {
         assert(var != implvar);
         CHECK_OKAY( implicsAdd(&var->implics, memhdr, set, i, implvar, impltype, implbound, conflict) );
      }

      debugMessage("implication added: "); 
      debugMessage("<%s> %s  ==>  <%s> %s %g\n", SCIPvarGetName(var), i == 0 ? "<= 0" : ">= 1", SCIPvarGetName(implvar), 
         impltype == SCIP_BOUNDTYPE_UPPER ? "<=" : ">=", implbound);
      break;
      
   case SCIP_VARSTATUS_FIXED:
      /* nothing to do here */
      if( i == (var->glbdom.lb > 0.5) ) /**@todo implement implication */ 
         printf("implication can be applied: <%s> %s %g\n", SCIPvarGetName(implvar), 
            impltype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", implbound);
      break;
      
   case SCIP_VARSTATUS_AGGREGATED:
      /* implication added for x >= 1:
       *   x >= 1 && x =  1*z + 0  ==>  y <= b or y >= b    <==>    z >= 1  ==>  y <= b or y >= b 
       *   x >= 1 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z <= 0  ==>  y <= b or y >= b
       * implication added for x <= 0:
       *   x <= 0 && x =  1*z + 0  ==>  y <= b or y >= b    <==>    z <= 0  ==>  y <= b or y >= b
       *   x <= 0 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z >= 1  ==>  y <= b or y >= b
       *
       * use only binary variables z 
       */
      assert(var->data.aggregate.var != NULL);
      if( SCIPvarGetType(var->data.aggregate.var) == SCIP_VARTYPE_BINARY )
      {   
         assert((var->data.aggregate.scalar == 1 && var->data.aggregate.constant == 0) || 
            (var->data.aggregate.scalar == -1 && var->data.aggregate.constant == 1));
       
         if( var->data.aggregate.scalar > 0 )
         {
            CHECK_OKAY( SCIPvarAddImplic(var->data.aggregate.var, memhdr, set, i, implvar, impltype, implbound, 
                  conflict) );
         }
         else
         {
            CHECK_OKAY( SCIPvarAddImplic(var->data.aggregate.var, memhdr, set, !i, implvar, impltype, implbound,
                  conflict) );
         }
      }
      break;
      
   case SCIP_VARSTATUS_MULTAGGR:
      /* nothing to do here */
      break;
      
   case SCIP_VARSTATUS_NEGATED:
      /* implication added for x >= 1:
       *   x >= 1 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z <= 0  ==>  y <= b or y >= b
       * implication added for x <= 0:
       *   x <= 0 && x = -1*z + 1  ==>  y <= b or y >= b    <==>    z >= 1  ==>  y <= b or y >= b
       */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      assert(SCIPvarGetType(var->negatedvar) == SCIP_VARTYPE_BINARY);

      CHECK_OKAY( SCIPvarAddImplic(var->negatedvar, memhdr, set, !i, implvar, impltype, implbound, conflict) );
      break;
      
   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
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

/** replaces variables in implications of binary variable by their active problem variable counterparts */
RETCODE SCIPvarUseActiveImplics(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            infeasible          /**< pointer to store TRUE, if an infeasibility was detected */
   )
{
   IMPLICS* oldimplics;
   int j;
   Bool conflict;

   assert(var != NULL);
   assert(infeasible != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   if( var->implics == NULL )
      return SCIP_OKAY;

   oldimplics = var->implics;
   var->implics = NULL;

   for( j = 0; j < oldimplics->nimpls[0]; ++j )
   {
      CHECK_OKAY( SCIPvarAddImplic(var, memhdr, set, 0, oldimplics->implvars[0][j], oldimplics->impltypes[0][j],
            oldimplics->implbounds[0][j], &conflict) );

      if( conflict )
      {
         if( SCIPvarGetUbGlobal(var) > 0.5 )
         {
            CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, 1.0) );
            CHECK_OKAY( SCIPvarChgLbGlobal(var, set, 1.0) );
         }
         else
            *infeasible = TRUE;
      }
   }
     
   for( j = 0; j < oldimplics->nimpls[1]; ++j )
   {
      CHECK_OKAY( SCIPvarAddImplic(var, memhdr, set, 1, oldimplics->implvars[1][j], oldimplics->impltypes[1][j],
            oldimplics->implbounds[1][j], &conflict) );

      if( conflict )
      {
         if( SCIPvarGetLbGlobal(var) < 0.5 )
         {
            CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, 0.0) );
            CHECK_OKAY( SCIPvarChgUbGlobal(var, set, 0.0) );
         }
         else
            *infeasible = TRUE;
      }
   }

   implicsFree(&oldimplics, memhdr);
      
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
   branchfactor = MAX(branchfactor, SCIPsetEpsilon(set));

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
         assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);
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
         varProcessChgBranchPriority(parentvar, branchpriority);
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }
}

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 */
void SCIPvarChgBranchPriority(
   VAR*             var,                /**< problem variable */
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
         SCIPvarChgBranchPriority(var->data.transvar, branchpriority);
      else
         var->branchpriority = branchpriority;
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      varProcessChgBranchPriority(var, branchpriority);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      SCIPvarChgBranchPriority(var->data.aggregate.var, branchpriority);
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      for( v = 0; v < var->data.multaggr.nvars; ++v )
         SCIPvarChgBranchPriority(var->data.multaggr.vars[v], branchpriority);
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIPvarChgBranchPriority(var->negatedvar, branchpriority);
      break;
         
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** actually changes the branch direction of the variable and of all parent variables */
static
void varProcessChgBranchDirection(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   VAR* parentvar;
   int i;

   assert(var != NULL);

   debugMessage("process changing branch direction of <%s> from %d to %d\n", 
      var->name, var->branchdirection, branchdirection);

   if( branchdirection == var->branchdirection )
      return;

   /* change the branch direction */
   var->branchdirection = branchdirection;

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         /* do not change directions across the border between transformed and original problem */
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED:
         if( parentvar->data.aggregate.scalar > 0.0 )
            varProcessChgBranchDirection(parentvar, branchdirection);
         else
            varProcessChgBranchDirection(parentvar, SCIPbranchdirOpposite(branchdirection));
         break;

      case SCIP_VARSTATUS_NEGATED:
         varProcessChgBranchDirection(parentvar, SCIPbranchdirOpposite(branchdirection));
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }
}

/** sets the branch direction of the variable; variables with higher branch direction are always preferred to variables
 *  with lower direction in selection of branching variable
 */
void SCIPvarChgBranchDirection(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   int v;

   assert(var != NULL);

   debugMessage("changing branch direction of <%s> from %d to %d\n", var->name, var->branchdirection, branchdirection);

   if( var->branchdirection == branchdirection )
      return;

   /* change directions of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarChgBranchDirection(var->data.transvar, branchdirection);
      else
         var->branchdirection = branchdirection;
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      varProcessChgBranchDirection(var, branchdirection);
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
         SCIPvarChgBranchDirection(var->data.aggregate.var, branchdirection);
      else
         SCIPvarChgBranchDirection(var->data.aggregate.var, SCIPbranchdirOpposite(branchdirection));
      break;
         
   case SCIP_VARSTATUS_MULTAGGR:
      for( v = 0; v < var->data.multaggr.nvars; ++v )
      {
         /* only update branching direction of aggregation variables, if they don't have a preferred direction yet */
         assert(var->data.multaggr.vars[v] != NULL);
         if( var->data.multaggr.vars[v]->branchdirection == SCIP_BRANCHDIR_AUTO )
         {
            if( var->data.multaggr.scalars[v] > 0.0 )
               SCIPvarChgBranchDirection(var->data.multaggr.vars[v], branchdirection);
            else
               SCIPvarChgBranchDirection(var->data.multaggr.vars[v], SCIPbranchdirOpposite(branchdirection));
         }
      }
      break;

   case SCIP_VARSTATUS_NEGATED:
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      SCIPvarChgBranchDirection(var->negatedvar, SCIPbranchdirOpposite(branchdirection));
      break;
         
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** compares the index of two variables, returns -1 if first is smaller than, and +1 if first is greater than second
 *  variable index; returns 0 if both indices are equal, which means both variables are equal
 */
int SCIPvarCompare(
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

/** gets corresponding active problem variable of a binary variable and updates the given negation status;
 *  for fixed variables, NULL is returned and the negation status is switched iff the variable is fixed to one
 */
RETCODE SCIPvarGetProbvarBinary(
   VAR**            var,                /**< pointer to binary problem variable */
   Bool*            negated             /**< pointer to update the negation status */
   )
{
   assert(var != NULL);
   assert(negated != NULL);

   while( *var != NULL )
   {
      assert(SCIPvarGetType(*var) == SCIP_VARTYPE_BINARY);

      switch( SCIPvarGetStatus(*var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( (*var)->data.transvar == NULL )
            return SCIP_OKAY;
         *var = (*var)->data.transvar;
         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         return SCIP_OKAY;

      case SCIP_VARSTATUS_FIXED:
         assert(SCIPvarGetLbGlobal(*var) > 0.5 || SCIPvarGetUbGlobal(*var) < 0.5);
         *negated = *negated ^ (SCIPvarGetLbGlobal(*var) > 0.5);
         *var = NULL;
         break;

      case SCIP_VARSTATUS_AGGREGATED:  /* x = a'*x' + c'  =>  a*x + c == (a*a')*x' + (a*c' + c) */
         assert((*var)->data.aggregate.var != NULL);
         assert(SCIPvarGetType((*var)->data.aggregate.var) == SCIP_VARTYPE_BINARY);
         assert(EPSEQ((*var)->data.aggregate.constant, 0.0, 1e-06) || EPSEQ((*var)->data.aggregate.constant, 1.0, 1e-06));
         assert(EPSEQ((*var)->data.aggregate.scalar, 1.0, 1e-06) || EPSEQ((*var)->data.aggregate.scalar, -1.0, 1e-06));
         assert(EPSEQ((*var)->data.aggregate.constant, 0.0, 1e-06) == EPSEQ((*var)->data.aggregate.scalar, 1.0, 1e-06));
         *negated = *negated ^ ((*var)->data.aggregate.scalar < 0.0);
         *var = (*var)->data.aggregate.var;
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("multiple aggregated variable has no single corresponding active problem variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_NEGATED:     /* x =  - x' + c'  =>  a*x + c ==   (-a)*x' + (a*c' + c) */
         assert((*var)->negatedvar != NULL);
         *negated = !(*negated);
         *var = (*var)->negatedvar;
         break;

      default:
         errorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** transforms given variable, boundtype and bound to the corresponding active or fixed variable values */
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
   case SCIP_VARSTATUS_FIXED:
      break;
      
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
      return SCIP_INVALIDDATA;
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
         return SCIP_INVALIDDATA;
      }
   }
   *scalar = 0.0;

   return SCIP_OKAY;
}

/** gets objective value of variable in current LP; the value can be different from the bound stored in the variable's own
 *  data due to diving, that operate only on the LP without updating the variables
 */
Real SCIPvarGetObjLP(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetObjLP(var->data.transvar);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetObj(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->obj;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetObjLP(var->data.aggregate.var);
      
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot get the objective value of a multiple aggregated variable\n");
      abort();
      
   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return -SCIPvarGetObjLP(var->negatedvar);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets lower bound of variable in current LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
Real SCIPvarGetLbLP(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetLbLP(var->data.transvar);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetLb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->locdom.lb;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbLP(var->data.aggregate.var) + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbLP(var->data.aggregate.var) + var->data.aggregate.constant;
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
      return var->data.negate.constant - SCIPvarGetUbLP(var->negatedvar);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets upper bound of variable in current LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
Real SCIPvarGetUbLP(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetUbLP(var->data.transvar);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetUb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->locdom.ub;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbLP(var->data.aggregate.var) + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbLP(var->data.aggregate.var) + var->data.aggregate.constant;
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
      return var->data.negate.constant - SCIPvarGetLbLP(var->negatedvar);
      
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

/** returns the solution of the variable in the root node's relaxation, if the root relaxation is not yet completely
 *  solved, zero is returned
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
         return 0.0;
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

/** returns a weighted average solution value of the variable in all feasible primal solutions found so far */
Real SCIPvarGetAvgSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real avgsol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0.0;
      return SCIPvarGetAvgSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      avgsol = var->primsolavg;
      avgsol = MIN(avgsol, var->glbdom.lb);
      avgsol = MAX(avgsol, var->glbdom.ub);
      return avgsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->locdom.lb;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      return var->data.aggregate.scalar * SCIPvarGetAvgSol(var->data.aggregate.var)
         + var->data.aggregate.constant;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      avgsol = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
         avgsol += var->data.multaggr.scalars[i] * SCIPvarGetAvgSol(var->data.multaggr.vars[i]);
      return avgsol;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPvarGetAvgSol(var->negatedvar);
      
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
   assert(!SCIPsetIsInfinity(set, REALABS(val)));

   debugMessage("adding coefficient %g<%s> to row <%s>\n", val, var->name, row->name);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("cannot add untransformed original variable <%s> to LP row <%s>\n", var->name, row->name);
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
      assert(!SCIPsetIsInfinity(set, REALABS(var->locdom.lb)));
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
      SCIPhistoryUpdatePseudocost(var->historycrun, set, solvaldelta, objdelta, weight);
      SCIPhistoryUpdatePseudocost(stat->glbhistory, set, solvaldelta, objdelta, weight);
      SCIPhistoryUpdatePseudocost(stat->glbhistorycrun, set, solvaldelta, objdelta, weight);
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

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value,
 *  only using the pseudo cost information of the current run
 */
Real SCIPvarGetPseudocostCurrentRun(
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
         return SCIPhistoryGetPseudocost(stat->glbhistorycrun, solvaldelta);
      else
         return SCIPvarGetPseudocostCurrentRun(var->data.transvar, stat, solvaldelta);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      dir = (solvaldelta >= 0.0 ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);
      
      return SCIPhistoryGetPseudocostCount(var->historycrun, dir) > 0.0
         ? SCIPhistoryGetPseudocost(var->historycrun, solvaldelta)
         : SCIPhistoryGetPseudocost(stat->glbhistorycrun, solvaldelta);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      return SCIPvarGetPseudocostCurrentRun(var->data.aggregate.var, stat, var->data.aggregate.scalar * solvaldelta);
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostCurrentRun(var->negatedvar, stat, -solvaldelta);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
Real SCIPvarGetPseudocostCount(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   
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

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
Real SCIPvarGetPseudocostCountCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetPseudocostCountCurrentRun(var->data.transvar, dir);
      
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetPseudocostCount(var->historycrun, dir);
      
   case SCIP_VARSTATUS_FIXED:
      return 0.0;
      
   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetPseudocostCountCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetPseudocostCountCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetPseudocostCountCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));
      
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
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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
      SCIPhistoryIncNBranchings(var->historycrun, depth, dir);
      SCIPhistoryIncNBranchings(stat->glbhistory, depth, dir);
      SCIPhistoryIncNBranchings(stat->glbhistorycrun, depth, dir);
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
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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
      SCIPhistoryIncNInferences(var->historycrun, dir);
      SCIPhistoryIncNInferences(stat->glbhistory, dir);
      SCIPhistoryIncNInferences(stat->glbhistorycrun, dir);
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
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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
      SCIPhistoryIncNCutoffs(var->historycrun, dir);
      SCIPhistoryIncNCutoffs(stat->glbhistory, dir);
      SCIPhistoryIncNCutoffs(stat->glbhistorycrun, dir);
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
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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

/** returns the number of times, a bound of the variable was changed in given direction due to branching 
 *  in the current run
 */
Longint SCIPvarGetNBranchingsCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNBranchingsCurrentRun(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNBranchingsCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNBranchingsCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNBranchingsCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the average depth of bound changes in given direction due to branching on the variable */
Real SCIPvarGetAvgBranchdepth(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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

/** returns the average depth of bound changes in given direction due to branching on the variable
 *  in the current run
 */
Real SCIPvarGetAvgBranchdepthCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0.0;
      else
         return SCIPvarGetAvgBranchdepthCurrentRun(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return  SCIPhistoryGetAvgBranchdepth(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgBranchdepthCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetAvgBranchdepthCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgBranchdepthCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the number of inferences branching on this variable in given direction triggered */
Longint SCIPvarGetNInferences(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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

/** returns the number of inferences branching on this variable in given direction triggered
 *  in the current run
 */
Longint SCIPvarGetNInferencesCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNInferencesCurrentRun(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNInferences(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNInferencesCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNInferencesCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNInferencesCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the average number of inferences found after branching on the variable in given direction */
Real SCIPvarGetAvgInferences(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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

/** returns the average number of inferences found after branching on the variable in given direction
 *  in the current run
 */
Real SCIPvarGetAvgInferencesCurrentRun(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIPhistoryGetAvgInferences(stat->glbhistorycrun, dir);
      else
         return SCIPvarGetAvgInferencesCurrentRun(var->data.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->historycrun, dir) > 0
         ? SCIPhistoryGetAvgInferences(var->historycrun, dir)
         : SCIPhistoryGetAvgInferences(stat->glbhistorycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgInferencesCurrentRun(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgInferencesCurrentRun(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgInferencesCurrentRun(var->negatedvar, stat, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the number of cutoffs branching on this variable in given direction produced */
Longint SCIPvarGetNCutoffs(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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

/** returns the number of cutoffs branching on this variable in given direction produced in the current run */
Longint SCIPvarGetNCutoffsCurrentRun(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return 0;
      else
         return SCIPvarGetNCutoffsCurrentRun(var->data.transvar, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNCutoffs(var->historycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetNCutoffsCurrentRun(var->data.aggregate.var, dir);
      else
         return SCIPvarGetNCutoffsCurrentRun(var->data.aggregate.var, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNCutoffsCurrentRun(var->negatedvar, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns the average number of cutoffs found after branching on the variable in given direction */
Real SCIPvarGetAvgCutoffs(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

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

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run */
Real SCIPvarGetAvgCutoffsCurrentRun(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   assert(var != NULL);
   assert(stat != NULL);
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIPhistoryGetAvgCutoffs(stat->glbhistorycrun, dir);
      else
         return SCIPvarGetAvgCutoffsCurrentRun(var->data.transvar, stat, dir);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPhistoryGetNBranchings(var->historycrun, dir) > 0
         ? SCIPhistoryGetAvgCutoffs(var->historycrun, dir)
         : SCIPhistoryGetAvgCutoffs(stat->glbhistorycrun, dir);

   case SCIP_VARSTATUS_FIXED:
      return 0.0;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         return SCIPvarGetAvgCutoffsCurrentRun(var->data.aggregate.var, stat, dir);
      else
         return SCIPvarGetAvgCutoffsCurrentRun(var->data.aggregate.var, stat, SCIPbranchdirOpposite(dir));
      
   case SCIP_VARSTATUS_MULTAGGR:
      return 0.0;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetAvgCutoffsCurrentRun(var->negatedvar, stat, SCIPbranchdirOpposite(dir));
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}




/*
 * information methods for bound changes
 */

/** returns the bound change information for the last lower bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower bound was applied up to this point of time
 */
BDCHGINFO* SCIPvarGetLbchgInfo(
   VAR*             var,                /**< active problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   int i;

   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   /* search the correct bound change information for the given bound change index */
   if( after )
   {
      for( i = var->nlbchginfos-1; i >= 0; --i )
      {
         assert(var->lbchginfos[i].var == var);
         assert(var->lbchginfos[i].boundtype == SCIP_BOUNDTYPE_LOWER);
         
         if( !SCIPbdchgidxIsEarlier(bdchgidx, &var->lbchginfos[i].bdchgidx) )
            return &var->lbchginfos[i];
      }
   }
   else
   {
      for( i = var->nlbchginfos-1; i >= 0; --i )
      {
         assert(var->lbchginfos[i].var == var);
         assert(var->lbchginfos[i].boundtype == SCIP_BOUNDTYPE_LOWER);
         
         if( SCIPbdchgidxIsEarlier(&var->lbchginfos[i].bdchgidx, bdchgidx) )
            return &var->lbchginfos[i];
      }
   }

   return NULL;
}

/** returns the bound change information for the last upper bound change on given active problem variable before or
 *  after the bound change with the given index was applied;
 *  returns NULL, if no change to the upper bound was applied up to this point of time
 */
BDCHGINFO* SCIPvarGetUbchgInfo(
   VAR*             var,                /**< active problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   int i;

   assert(var != NULL);
   assert(SCIPvarIsActive(var));

   /* search the correct bound change information for the given bound change index */
   if( after )
   {
      for( i = var->nubchginfos-1; i >= 0; --i )
      {
         assert(var->ubchginfos[i].var == var);
         assert(var->ubchginfos[i].boundtype == SCIP_BOUNDTYPE_UPPER);
         
         if( !SCIPbdchgidxIsEarlier(bdchgidx, &var->ubchginfos[i].bdchgidx) )
            return &var->ubchginfos[i];
      }
   }
   else
   {
      for( i = var->nubchginfos-1; i >= 0; --i )
      {
         assert(var->ubchginfos[i].var == var);
         assert(var->ubchginfos[i].boundtype == SCIP_BOUNDTYPE_UPPER);
         
         if( SCIPbdchgidxIsEarlier(&var->ubchginfos[i].bdchgidx, bdchgidx) )
            return &var->ubchginfos[i];
      }
   }

   return NULL;
}

/** returns the bound change information for the last lower or upper bound change on given active problem variable
 *  before or after the bound change with the given index was applied;
 *  returns NULL, if no change to the lower/upper bound was applied up to this point of time
 */
BDCHGINFO* SCIPvarGetBdchgInfo(
   VAR*             var,                /**< active problem variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      return SCIPvarGetLbchgInfo(var, bdchgidx, after);
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      return SCIPvarGetUbchgInfo(var, bdchgidx, after);
   }
}

/** returns lower bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
Real SCIPvarGetLbAtIndex(
   VAR*             var,                /**< problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetLbAtIndex(var->data.transvar, bdchgidx, after);
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      if( bdchgidx == NULL )
         return SCIPvarGetLbLocal(var);
      else
      {
         BDCHGINFO* bdchginfo;
         
         bdchginfo = SCIPvarGetLbchgInfo(var, bdchgidx, after);
         if( bdchginfo != NULL )
            return SCIPbdchginfoGetNewbound(bdchginfo);
         else
            return var->glbdom.lb;
      }

   case SCIP_VARSTATUS_FIXED:
      return var->glbdom.lb;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
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
      return var->data.negate.constant - SCIPvarGetUbAtIndex(var->negatedvar, bdchgidx, after);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
Real SCIPvarGetUbAtIndex(
   VAR*             var,                /**< problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetUbAtIndex(var->data.transvar, bdchgidx, after);
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      if( bdchgidx == NULL )
         return SCIPvarGetUbLocal(var);
      else
      {
         BDCHGINFO* bdchginfo;
         
         bdchginfo = SCIPvarGetUbchgInfo(var, bdchgidx, after);
         if( bdchginfo != NULL )
            return SCIPbdchginfoGetNewbound(bdchginfo);
         else
            return var->glbdom.ub;
      }

   case SCIP_VARSTATUS_FIXED:
      return var->glbdom.ub;
      
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         /* a > 0 -> get lower bound of y */
         return var->data.aggregate.scalar * SCIPvarGetUbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         /* a < 0 -> get upper bound of y */
         return var->data.aggregate.scalar * SCIPvarGetLbAtIndex(var->data.aggregate.var, bdchgidx, after)
            + var->data.aggregate.constant;
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
      return var->data.negate.constant - SCIPvarGetLbAtIndex(var->negatedvar, bdchgidx, after);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** returns lower or upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
Real SCIPvarGetBdAtIndex(
   VAR*             var,                /**< problem variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      return SCIPvarGetLbAtIndex(var, bdchgidx, after);
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      return SCIPvarGetUbAtIndex(var, bdchgidx, after);
   }
}

/** returns whether the binary variable was fixed at the time given by the bound change index */
Bool SCIPvarWasFixedAtIndex(
   VAR*             var,                /**< problem variable */
   BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   Bool             after               /**< should the bound change with given index be included? */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   /* check the current bounds first in order to decide at which bound change information we have to look
    * (which is expensive because we have to follow the aggregation tree to the active variable)
    */
   return ((SCIPvarGetLbLocal(var) > 0.5 && SCIPvarGetLbAtIndex(var, bdchgidx, after) > 0.5)
      || (SCIPvarGetUbLocal(var) < 0.5 && SCIPvarGetUbAtIndex(var, bdchgidx, after) < 0.5));
}

/** bound change index representing the initial time before any bound changes took place */
static BDCHGIDX initbdchgidx = {-2, 0};

/** bound change index representing the presolving stage */
static BDCHGIDX presolvebdchgidx = {-1, 0};

/** returns the last bound change index, at which the bounds of the given variable were tightened */
BDCHGIDX* SCIPvarGetLastBdchgIndex(
   VAR*             var                 /**< problem variable */
   )
{
   BDCHGIDX* lbchgidx;
   BDCHGIDX* ubchgidx;

   assert(var != NULL);

   var = SCIPvarGetProbvar(var);

   /* check, if variable was fixed in presolving */
   if( var == NULL )
      return &presolvebdchgidx;

   /* get depths of last bound change infos for the lower and upper bound */
   lbchgidx = (var->nlbchginfos > 0 ? &var->lbchginfos[var->nlbchginfos-1].bdchgidx : &initbdchgidx);
   ubchgidx = (var->nubchginfos > 0 ? &var->ubchginfos[var->nubchginfos-1].bdchgidx : &initbdchgidx);

   if( SCIPbdchgidxIsEarlierNonNull(lbchgidx, ubchgidx) )
      return ubchgidx;
   else
      return lbchgidx;
}

/** returns the last depth level, at which the bounds of the given variable were tightened;
 *  returns -2, if the variable's bounds are still the global bounds
 *  returns -1, if the variable was fixed in presolving
 */
int SCIPvarGetLastBdchgDepth(
   VAR*             var                 /**< problem variable */
   )
{
   BDCHGIDX* bdchgidx;

   bdchgidx = SCIPvarGetLastBdchgIndex(var);
   assert(bdchgidx != NULL);

   return bdchgidx->depth;
}

/** returns whether the first binary variable was fixed earlier than the second one;
 *  returns FALSE, if the first variable is not fixed, and returns TRUE, if the first variable is fixed, but the
 *  second one is not fixed
 */
Bool SCIPvarWasFixedEarlier(
   VAR*             var1,               /**< first binary variable */
   VAR*             var2                /**< second binary variable */
   )
{
   BDCHGIDX* bdchgidx1;
   BDCHGIDX* bdchgidx2;

   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarGetType(var1) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarGetType(var2) == SCIP_VARTYPE_BINARY);

   var1 = SCIPvarGetProbvar(var1);
   var2 = SCIPvarGetProbvar(var2);

   /* check, if variables are globally fixed */
   if( var2 == NULL || var2->glbdom.lb > 0.5 || var2->glbdom.ub < 0.5 )
      return FALSE;
   if( var1 == NULL || var1->glbdom.lb > 0.5 || var1->glbdom.ub < 0.5 )
      return TRUE;

   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarGetType(var1) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarGetType(var2) == SCIP_VARTYPE_BINARY);
   assert(var1->nlbchginfos + var1->nubchginfos <= 1);
   assert(var2->nlbchginfos + var2->nubchginfos <= 1);

   if( var1->nlbchginfos == 1 )
      bdchgidx1 = &var1->lbchginfos[0].bdchgidx;
   else if( var1->nubchginfos == 1 )
      bdchgidx1 = &var1->ubchginfos[0].bdchgidx;
   else
      bdchgidx1 = NULL;

   if( var2->nlbchginfos == 1 )
      bdchgidx2 = &var2->lbchginfos[0].bdchgidx;
   else if( var2->nubchginfos == 1 )
      bdchgidx2 = &var2->ubchginfos[0].bdchgidx;
   else
      bdchgidx2 = NULL;

   return SCIPbdchgidxIsEarlier(bdchgidx1, bdchgidx2);
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




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPvarGetName
#undef SCIPvarGetData
#undef SCIPvarGetStatus
#undef SCIPvarIsOriginal
#undef SCIPvarIsTransformed
#undef SCIPvarIsNegated
#undef SCIPvarGetType
#undef SCIPvarIsIntegral
#undef SCIPvarIsInitial
#undef SCIPvarIsRemoveable
#undef SCIPvarIsActive
#undef SCIPvarGetIndex
#undef SCIPvarGetProbindex
#undef SCIPvarGetTransVar
#undef SCIPvarGetCol
#undef SCIPvarIsInLP
#undef SCIPvarGetAggrVar
#undef SCIPvarGetAggrScalar
#undef SCIPvarGetAggrConstant
#undef SCIPvarGetMultaggrNVars
#undef SCIPvarGetMultaggrVars
#undef SCIPvarGetMultaggrScalars
#undef SCIPvarGetMultaggrConstant
#undef SCIPvarGetNegatedVar
#undef SCIPvarGetNegationVar
#undef SCIPvarGetNegationConstant
#undef SCIPvarGetObj
#undef SCIPvarGetLbGlobal
#undef SCIPvarGetUbGlobal
#undef SCIPvarGetLbLocal
#undef SCIPvarGetUbLocal
#undef SCIPvarGetBranchFactor
#undef SCIPvarGetBranchPriority
#undef SCIPvarGetBranchDirection
#undef SCIPvarGetNVlbs
#undef SCIPvarGetVlbVars
#undef SCIPvarGetVlbCoefs
#undef SCIPvarGetVlbConstants
#undef SCIPvarGetNVubs
#undef SCIPvarGetVubVars
#undef SCIPvarGetVubCoefs
#undef SCIPvarGetVubConstants
#undef SCIPvarGetNimpl
#undef SCIPvarGetImplvars
#undef SCIPvarGetImpltypes
#undef SCIPvarGetImplbounds
#undef SCIPvarCatchEvent
#undef SCIPvarDropEvent
#undef SCIPbdchgidxIsEarlierNonNull
#undef SCIPbdchgidxIsEarlier
#undef SCIPbdchginfoGetOldbound
#undef SCIPbdchginfoGetNewbound
#undef SCIPbdchginfoGetVar
#undef SCIPbdchginfoGetChgtype
#undef SCIPbdchginfoGetBoundtype
#undef SCIPbdchginfoGetDepth
#undef SCIPbdchginfoGetPos
#undef SCIPbdchginfoGetIdx
#undef SCIPbdchginfoGetInferVar
#undef SCIPbdchginfoGetInferCons
#undef SCIPbdchginfoGetInferProp
#undef SCIPbdchginfoGetInferInfo
#undef SCIPbdchginfoGetInferBoundtype
#undef SCIPbdchginfoHasInferenceReason

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
   assert(var->varstatus != SCIP_VARSTATUS_NEGATED || var->negatedvar != NULL);

   return (var->varstatus == SCIP_VARSTATUS_ORIGINAL
      || (var->varstatus == SCIP_VARSTATUS_NEGATED && var->negatedvar->varstatus == SCIP_VARSTATUS_ORIGINAL));
}

/** returns whether the variable belongs to the transformed problem */
Bool SCIPvarIsTransformed(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_NEGATED || var->negatedvar != NULL);

   return (var->varstatus != SCIP_VARSTATUS_ORIGINAL
      && (var->varstatus != SCIP_VARSTATUS_NEGATED || var->negatedvar->varstatus != SCIP_VARSTATUS_ORIGINAL));
}

/** returns whether the variable was created by negation of a different variable */
Bool SCIPvarIsNegated(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (var->varstatus == SCIP_VARSTATUS_NEGATED);
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
   assert(var->varstatus == SCIP_VARSTATUS_ORIGINAL);

   return var->data.transvar;
}

/** gets column of COLUMN variable */
COL* SCIPvarGetCol(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_COLUMN);

   return var->data.col;
}

/** returns whether the variable is a COLUMN variable that is member of the current LP */
Bool SCIPvarIsInLP(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (var->varstatus == SCIP_VARSTATUS_COLUMN && SCIPcolIsInLP(var->data.col));
}

/** gets aggregation variable y of an aggregated variable x = a*y + c */
VAR* SCIPvarGetAggrVar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.var;
}

/** gets aggregation scalar a of an aggregated variable x = a*y + c */
Real SCIPvarGetAggrScalar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.scalar;
}

/** gets aggregation constant c of an aggregated variable x = a*y + c */
Real SCIPvarGetAggrConstant(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.constant;
}

/** gets number n of aggregation variables of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
int SCIPvarGetMultaggrNVars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.nvars;
}

/** gets vector of aggregation variables y of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
VAR** SCIPvarGetMultaggrVars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.vars;
}

/** gets vector of aggregation scalars a of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
Real* SCIPvarGetMultaggrScalars(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_MULTAGGR);

   return var->data.multaggr.scalars;
}

/** gets aggregation constant c of a multi aggregated variable x = a0*y0 + ... + a(n-1)*y(n-1) + c */
Real SCIPvarGetMultaggrConstant(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_MULTAGGR);

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
   assert(var->varstatus == SCIP_VARSTATUS_NEGATED);

   return var->negatedvar;
}

/** gets the negation offset of a negated variable x' = offset - x */
Real SCIPvarGetNegationConstant(
   VAR*             var                 /**< negated problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_NEGATED);

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

/** gets the branch priority of the variable; variables with higher priority should always be preferred to variables
 *  with lower priority
 */
int SCIPvarGetBranchPriority(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->branchpriority;
}

/** gets the preferred branch direction of the variable (downwards, upwards, or auto) */
BRANCHDIR SCIPvarGetBranchDirection(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->branchdirection;
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

/** gets number of implications  y <= b or y >= b for x <= 0 or x >= 1 of given variable x, 
 *  there are no implications for nonbinary variable x
 */
int SCIPvarGetNimpl(
   VAR*             var,                /**< problem variable */
   Bool             i                   /**< FALSE for implications for x <= 0, TRUE for x >= 1 */
   )
{
   assert(var != NULL);

   return var->implics != NULL ? var->implics->nimpls[i] : 0;
}

/** gets array with implication variables y of implications  y <= b or y >= b for x <= 0 or x >= 1 of given variable x,  
 *  there are no implications for nonbinary variable x
 */
VAR** SCIPvarGetImplvars(
   VAR*             var,                /**< problem variable */
   Bool             i                   /**< FALSE for implications for x <= 0, TRUE for x >= 1 */
   )
{
   assert(var != NULL);

   return var->implics != NULL ? var->implics->implvars[i] : NULL;
}

/** gets array with implication types of implications  y <= b or y >= b for x <= 0 or x >= 1 of given variable x
 *  (SCIP_BOUNDTYPE_UPPER if y <= b, SCIP_BOUNDTYPE_LOWER if y >= b), 
 *  there are no implications for nonbinary variable x
 */
BOUNDTYPE* SCIPvarGetImpltypes(
   VAR*             var,                /**< problem variable */
   Bool             i                   /**< FALSE for implications for x <= 0, TRUE for x >= 1 */
   )
{
   assert(var != NULL);

   return var->implics != NULL ? var->implics->impltypes[i] : NULL;
}

/** gets array with implication bounds b of implications  y <= b or y >= b for x <= 0 or x >= 1 of given variable x,  
 *  there are no implications for nonbinary variable x
 */
Real* SCIPvarGetImplbounds(
   VAR*             var,                /**< problem variable */
   Bool             i                   /**< FALSE for implications for x <= 0, TRUE for x >= 1 */
   )
{
   assert(var != NULL);

   return var->implics != NULL ? var->implics->implbounds[i] : NULL;
}

/** includes event handler with given data in variable's event filter */
RETCODE SCIPvarCatchEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*             filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert((eventtype & ~SCIP_EVENTTYPE_VARCHANGED) == 0);
   assert((eventtype & SCIP_EVENTTYPE_VARCHANGED) != 0);

   debugMessage("catch event of type 0x%x of variable <%s> with handler %p and data %p\n", 
      eventtype, var->name, eventhdlr, eventdata);

   CHECK_OKAY( SCIPeventfilterAdd(var->eventfilter, memhdr, set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** deletes event handler with given data from variable's event filter */
RETCODE SCIPvarDropEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int              filterpos           /**< position of event filter entry returned by SCIPvarCatchEvent(), or -1 */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);

   debugMessage("drop event of variable <%s> with handler %p and data %p\n", var->name, eventhdlr, eventdata);

   CHECK_OKAY( SCIPeventfilterDel(var->eventfilter, memhdr, set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** returns whether first bound change index belongs to an earlier applied bound change than second one */
Bool SCIPbdchgidxIsEarlierNonNull(
   BDCHGIDX*        bdchgidx1,          /**< first bound change index */
   BDCHGIDX*        bdchgidx2           /**< second bound change index */
   )
{
   assert(bdchgidx1 != NULL);
   assert(bdchgidx1->depth >= -2);
   assert(bdchgidx1->pos >= 0);
   assert(bdchgidx2 != NULL);
   assert(bdchgidx2->depth >= -2);
   assert(bdchgidx2->pos >= 0);

   return (bdchgidx1->depth < bdchgidx2->depth)
      || (bdchgidx1->depth == bdchgidx2->depth && (bdchgidx1->pos < bdchgidx2->pos));
}

/** returns whether first bound change index belongs to an earlier applied bound change than second one;
 *  if a bound change index is NULL, the bound change index represents the current time, i.e. the time after the
 *  last bound change was applied to the current node
 */
Bool SCIPbdchgidxIsEarlier(
   BDCHGIDX*        bdchgidx1,          /**< first bound change index, or NULL */
   BDCHGIDX*        bdchgidx2           /**< second bound change index, or NULL */
   )
{
   assert(bdchgidx1 == NULL || bdchgidx1->depth >= -2);
   assert(bdchgidx1 == NULL || bdchgidx1->pos >= 0);
   assert(bdchgidx2 == NULL || bdchgidx2->depth >= -2);
   assert(bdchgidx2 == NULL || bdchgidx2->pos >= 0);

   if( bdchgidx1 == NULL )
      return FALSE;
   else if( bdchgidx2 == NULL )
      return TRUE;
   else
      return (bdchgidx1->depth < bdchgidx2->depth)
         || (bdchgidx1->depth == bdchgidx2->depth && (bdchgidx1->pos < bdchgidx2->pos));
}

/** returns old bound that was overwritten for given bound change information */
Real SCIPbdchginfoGetOldbound(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->oldbound;
}

/** returns new bound installed for given bound change information */
Real SCIPbdchginfoGetNewbound(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->newbound;
}

/** returns variable that belongs to the given bound change information */
VAR* SCIPbdchginfoGetVar(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->var;
}

/** returns whether the bound change information belongs to a branching decision or a deduction */
BOUNDCHGTYPE SCIPbdchginfoGetChgtype(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (BOUNDCHGTYPE)(bdchginfo->boundchgtype);
}

/** returns whether the bound change information belongs to a lower or upper bound change */
BOUNDTYPE SCIPbdchginfoGetBoundtype(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (BOUNDTYPE)(bdchginfo->boundtype);
}

/** returs depth level of given bound change information */
int SCIPbdchginfoGetDepth(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->bdchgidx.depth;
}

/** returs bound change position in its depth level of given bound change information */
int SCIPbdchginfoGetPos(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return bdchginfo->bdchgidx.pos;
}

/** returs bound change index of given bound change information */
BDCHGIDX* SCIPbdchginfoGetIdx(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return &bdchginfo->bdchgidx;
}

/** returs inference variable of given bound change information */
VAR* SCIPbdchginfoGetInferVar(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
      || bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return bdchginfo->inferencedata.var;
}

/** returs inference constraint of given bound change information */
CONS* SCIPbdchginfoGetInferCons(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER);
   assert(bdchginfo->inferencedata.reason.cons != NULL);

   return bdchginfo->inferencedata.reason.cons;
}

/** returs inference propagator of given bound change information, or NULL if no propagator was responsible */
PROP* SCIPbdchginfoGetInferProp(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return bdchginfo->inferencedata.reason.prop;
}

/** returs inference user information of given bound change information */
int SCIPbdchginfoGetInferInfo(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
      || bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return bdchginfo->inferencedata.info;
}

/** returs inference bound of inference variable of given bound change information */
BOUNDTYPE SCIPbdchginfoGetInferBoundtype(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER
      || bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER);

   return (BOUNDTYPE)(bdchginfo->inferboundtype);
}

/** returns whether the bound change has an inference reason (constraint or propagator), that can be resolved */
Bool SCIPbdchginfoHasInferenceReason(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_CONSINFER)
      || (bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_PROPINFER && bdchginfo->inferencedata.reason.prop != NULL);
}
