/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   var.c
 * @brief  Methods and datastructures for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "var.h"



/** hole in a domain */
struct Hole
{
   Real             left;               /**< left bound of open interval defining the hole $(left,right)$ */
   Real             right;              /**< right bound of open interval defining the hole $(left,right)$ */
};

/** list of domain holes */
struct Holelist
{
   HOLE             hole;               /**< this hole */
   HOLELIST*        next;               /**< next hole in list */
};

/** change in a hole list */
struct HoleChg
{
   HOLELIST**       ptr;                /**< changed list pointer */
   HOLELIST*        newlist;            /**< new value of list pointer */
   HOLELIST*        oldlist;            /**< old value of list pointer */
};

/** change in one bound of a variable */
struct BoundChg
{
   VAR*             var;                /**< variable to change the bounds for */
   Real             newbound;           /**< new value for bound */
   Real             oldbound;           /**< old value for bound */
   BOUNDTYPE        boundtype;          /**< type of bound: lower or upper bound */
};

/** dynamic size attachment for domain change data */
struct DomChgDyn
{
   DOMCHG**         domchg;             /**< pointer to domain change data */
   int              boundchgsize;       /**< size of bound change array */
   int              holechgsize;        /**< size of hole change array */
};




/*
 * memory growing methods for dynamically allocated arrays
 */

/** creates empty fixed size domain change data */
static
RETCODE domchgCreate(
   DOMCHG**         domchg,             /**< pointer to fixed size domain change data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, domchg) );
   (*domchg)->boundchg = NULL;
   (*domchg)->holechg = NULL;
   (*domchg)->nboundchg = 0;
   (*domchg)->nholechg = 0;

   return SCIP_OKAY;
}

/** ensures, that boundchg array can store at least num entries */
static
RETCODE ensureBoundchgSize(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   DOMCHG** domchg;

   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg != NULL);

   domchg = domchgdyn->domchg;
   assert(*domchg != NULL || domchgdyn->boundchgsize == 0);
   assert(*domchg == NULL || (*domchg)->nboundchg <= domchgdyn->boundchgsize);

   if( num > domchgdyn->boundchgsize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      if( *domchg == NULL )
      {
         CHECK_OKAY( domchgCreate(domchg, memhdr) );
      }
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchg)->boundchg, domchgdyn->boundchgsize, newsize) );
      domchgdyn->boundchgsize = newsize;
   }
   assert(num <= domchgdyn->boundchgsize);

   return SCIP_OKAY;
}

/** ensures, that holechg array can store at least num additional entries */
static
RETCODE ensureHolechgSize(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   DOMCHG** domchg;

   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg != NULL);

   domchg = domchgdyn->domchg;
   assert(*domchg != NULL || domchgdyn->holechgsize == 0);
   assert(*domchg == NULL || (*domchg)->nholechg <= domchgdyn->holechgsize);

   if( num > domchgdyn->holechgsize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      if( *domchg == NULL )
      {
         CHECK_OKAY( domchgCreate(domchg, memhdr) );
      }
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchg)->holechg, domchgdyn->holechgsize, newsize) );
      domchgdyn->holechgsize = newsize;
   }
   assert(num <= domchgdyn->holechgsize);

   return SCIP_OKAY;
}



/*
 * hole, holelist, and domain methods
 */

/** creates a new holelist element */
static
RETCODE holelistCreate(
   HOLELIST**       holelist,           /**< pointer to holelist to create */
   MEMHDR*          memhdr,             /**< block memory for target holelist */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(holelist != NULL);
   assert(memhdr != NULL);
   assert(SCIPsetIsL(set, left, right));

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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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

/** merges overlapping holes into single holes, moves bounds respectively */
static
RETCODE domMerge(
   DOM*             dom,                /**< domain to merge */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
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
      assert(SCIPsetIsL(set, (*holelistptr)->hole.left, (*holelistptr)->hole.right));

      if( SCIPsetIsGE(set, (*holelistptr)->hole.left, dom->ub) )
      {
         /* the remaining holes start behind the upper bound: kill them */
         holelistFree(holelistptr, memhdr);
         assert(*holelistptr == NULL);
      }
      else if( SCIPsetIsG(set, (*holelistptr)->hole.right, dom->ub) )
      {
         /* the hole overlaps the upper bound: decrease upper bound, kill this and all remaining holes */
         dom->ub = (*holelistptr)->hole.left;
         holelistFree(holelistptr, memhdr);
         assert(*holelistptr == NULL);
      }
      else if( SCIPsetIsG(set, *lastrightptr, (*holelistptr)->hole.left) )
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


/*
 * domain change methods
 */

/** frees fixed size domain change data */
void SCIPdomchgFree(
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   if( *domchg != NULL )
   {
      freeBlockMemoryArrayNull(memhdr, &(*domchg)->boundchg, (*domchg)->nboundchg);
      freeBlockMemoryArrayNull(memhdr, &(*domchg)->holechg, (*domchg)->nholechg);
      freeBlockMemory(memhdr, domchg);
   }
}

/** applies domain change */
RETCODE SCIPdomchgApply(
   const DOMCHG*    domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   int i;

   assert(lp != NULL);

   debugMessage("applying domain changes at %p\n", domchg);
   if( domchg == NULL )
      return SCIP_OKAY;
   debugMessage(" -> %d bound changes, %d hole changes\n", domchg->nboundchg, domchg->nholechg);

   /* apply bound changes */
   for( i = 0; i < domchg->nboundchg; ++i )
   {
      var = domchg->boundchg[i].var;
      assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

      /* apply bound change to the LP data */
      switch( domchg->boundchg[i].boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:
         CHECK_OKAY( SCIPvarChgLb(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, domchg->boundchg[i].newbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUb(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, domchg->boundchg[i].newbound) );
         break;
      default:
         errorMessage("Unknown bound type");
         abort();
      }
   }

   /* apply holelist changes */
   for( i = 0; i < domchg->nholechg; ++i )
      *(domchg->holechg[i].ptr) = domchg->holechg[i].newlist;

   return SCIP_OKAY;
}
   
/** undoes domain change */
RETCODE SCIPdomchgUndo(
   const DOMCHG*    domchg,             /**< domain change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   int i;

   assert(lp != NULL);

   debugMessage("undoing domain changes at %p\n", domchg);
   if( domchg == NULL )
      return SCIP_OKAY;
   debugMessage(" -> %d bound changes, %d hole changes\n", domchg->nboundchg, domchg->nholechg);

   /* undo bound changes */
   for( i = domchg->nboundchg-1; i >= 0; --i )
   {
      var = domchg->boundchg[i].var;
      assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

      /* apply bound change to the LP data */
      switch( domchg->boundchg[i].boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:
         CHECK_OKAY( SCIPvarChgLb(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, domchg->boundchg[i].oldbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUb(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, domchg->boundchg[i].oldbound) );
         break;
      default:
         errorMessage("Unknown bound type");
         abort();
      }
   }

   /* undo holelist changes */
   for( i = domchg->nholechg-1; i >= 0; --i )
      *(domchg->holechg[i].ptr) = domchg->holechg[i].oldlist;

   return SCIP_OKAY;
}


/*
 * dynamic size attachment methods
 */

/** creates a dynamic size attachment for a domain change data structure */
RETCODE SCIPdomchgdynCreate(
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchgdyn != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, domchgdyn) );

   (*domchgdyn)->domchg = NULL;
   (*domchgdyn)->boundchgsize = 0;
   (*domchgdyn)->holechgsize = 0;

   return SCIP_OKAY;
}

/** frees a dynamic size attachment for a domain change data structure */
void SCIPdomchgdynFree(
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchgdyn != NULL);
   assert(*domchgdyn != NULL);
   assert((*domchgdyn)->domchg == NULL);
   assert((*domchgdyn)->boundchgsize == 0);
   assert((*domchgdyn)->holechgsize == 0);

   freeBlockMemory(memhdr, domchgdyn);
}

/** attaches dynamic size information to domain change data */
void SCIPdomchgdynAttach(
   DOMCHGDYN*       domchgdyn,          /**< dynamic size information */
   DOMCHG**         domchg              /**< pointer to static domain change */
   )
{
   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg == NULL);
   assert(domchgdyn->boundchgsize == 0);
   assert(domchgdyn->holechgsize == 0);
   assert(domchg != NULL);

   debugMessage("attaching dynamic size information at %p to domain change data at %p\n", domchgdyn, domchg);
   domchgdyn->domchg = domchg;
   if( *domchg != NULL )
   {
      domchgdyn->boundchgsize = (*domchg)->nboundchg;
      domchgdyn->holechgsize = (*domchg)->nholechg;
   }
   else
   {
      domchgdyn->boundchgsize = 0;
      domchgdyn->holechgsize = 0;
   }
}

/** detaches dynamic size information and shrinks domain change data */
RETCODE SCIPdomchgdynDetach(
   DOMCHGDYN*       domchgdyn,          /**< dynamic size information */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchgdyn != NULL);
   debugMessage("detaching dynamic size information at %p from domain change data at %p\n", domchgdyn, domchgdyn->domchg);
   assert(domchgdyn->domchg != NULL);
   assert(domchgdyn->boundchgsize == 0 || *domchgdyn->domchg != NULL);
   assert(domchgdyn->boundchgsize == 0 || (*domchgdyn->domchg)->nboundchg <= domchgdyn->boundchgsize);
   assert(domchgdyn->holechgsize == 0 || *domchgdyn->domchg != NULL);
   assert(domchgdyn->holechgsize == 0 || (*domchgdyn->domchg)->nholechg <= domchgdyn->holechgsize);

   /* shrink static domain change data to the size of the used elements */
   if( *domchgdyn->domchg != NULL )
   {
      if( (*domchgdyn->domchg)->nboundchg == 0 )
      {
         freeBlockMemoryArrayNull(memhdr, &(*domchgdyn->domchg)->boundchg, domchgdyn->boundchgsize);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchgdyn->domchg)->boundchg,
                        domchgdyn->boundchgsize, (*domchgdyn->domchg)->nboundchg) );
      }
      if( (*domchgdyn->domchg)->nholechg == 0 )
      {
         freeBlockMemoryArrayNull(memhdr, &(*domchgdyn->domchg)->holechg, domchgdyn->holechgsize);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*domchgdyn->domchg)->holechg,
                        domchgdyn->holechgsize, (*domchgdyn->domchg)->nholechg) );
      }
      if( (*domchgdyn->domchg)->nboundchg == 0 && (*domchgdyn->domchg)->nholechg == 0 )
         SCIPdomchgFree(domchgdyn->domchg, memhdr);
   }
   else
   {
      assert(domchgdyn->boundchgsize == 0);
      assert(domchgdyn->holechgsize == 0);
   }

   /* detach domain change data */
   domchgdyn->domchg = NULL;
   domchgdyn->holechgsize = 0;
   domchgdyn->boundchgsize = 0;

   return SCIP_OKAY;
}

/** frees attached domain change data and detaches dynamic size attachment */
void SCIPdomchgdynDiscard(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchgdyn != NULL);
   debugMessage("discarding dynamic size information at %p from domain change data at %p\n", domchgdyn, domchgdyn->domchg);
   assert(domchgdyn->domchg != NULL);
   assert(domchgdyn->boundchgsize == 0 || *domchgdyn->domchg != NULL);
   assert(domchgdyn->boundchgsize == 0 || (*domchgdyn->domchg)->nboundchg <= domchgdyn->boundchgsize);
   assert(domchgdyn->holechgsize == 0 || *domchgdyn->domchg != NULL);
   assert(domchgdyn->holechgsize == 0 || (*domchgdyn->domchg)->nholechg <= domchgdyn->holechgsize);

   /* free static domain change data */
   if( *domchgdyn->domchg != NULL )
   {      
      freeBlockMemoryArrayNull(memhdr, &(*domchgdyn->domchg)->boundchg, domchgdyn->boundchgsize);
      (*domchgdyn->domchg)->nboundchg = 0;
      freeBlockMemoryArrayNull(memhdr, &(*domchgdyn->domchg)->holechg, domchgdyn->holechgsize);
      (*domchgdyn->domchg)->nholechg = 0;
      SCIPdomchgFree(domchgdyn->domchg, memhdr);
      domchgdyn->boundchgsize = 0;
      domchgdyn->holechgsize = 0;
   }

   /* detach domain change data */
   domchgdyn->domchg = NULL;
}

/** adds bound change to domain changes */
RETCODE SCIPdomchgdynAddBoundchg(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   DOMCHG* domchg;
   int nboundchg;

   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg != NULL);
   assert(var != NULL);

   debugMessage("(1) adding bound change <%s>: %g -> %g of variable <%s> to domain change at %p pointing to %p\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", oldbound, newbound, var->name, domchgdyn->domchg,
      *domchgdyn->domchg);
   nboundchg = (*domchgdyn->domchg == NULL ? 0 : (*domchgdyn->domchg)->nboundchg);
   CHECK_OKAY( ensureBoundchgSize(domchgdyn, memhdr, set, nboundchg+1) );
   debugMessage("(2) adding bound change <%s>: %g -> %g of variable <%s> to domain change at %p pointing to %p\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", oldbound, newbound, var->name, domchgdyn->domchg,
      *domchgdyn->domchg);

   domchg = *domchgdyn->domchg;
   assert(domchg != NULL);

   domchg->boundchg[domchg->nboundchg].var = var;
   domchg->boundchg[domchg->nboundchg].newbound = newbound;
   domchg->boundchg[domchg->nboundchg].oldbound = oldbound;
   domchg->boundchg[domchg->nboundchg].boundtype = boundtype;
   domchg->nboundchg++;

   return SCIP_OKAY;
}

/** adds hole change to domain changes */
RETCODE SCIPdomchgdynAddHolechg(
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   DOMCHG* domchg;
   int nholechg;

   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg != NULL);
   assert(ptr != NULL);

   nholechg = (*domchgdyn->domchg == NULL ? 0 : (*domchgdyn->domchg)->nholechg);
   CHECK_OKAY( ensureHolechgSize(domchgdyn, memhdr, set, nholechg+1) );

   domchg = *domchgdyn->domchg;
   assert(domchg != NULL);

   domchg->holechg[domchg->nholechg].ptr = ptr;
   domchg->holechg[domchg->nholechg].newlist = newlist;
   domchg->holechg[domchg->nholechg].oldlist = oldlist;
   domchg->nholechg++;

   return SCIP_OKAY;
}

/** gets pointer to domain change data the dynamic size information references */
DOMCHG** SCIPdomchgdynGetDomchgPtr(
   DOMCHGDYN*       domchgdyn           /**< dynamically sized domain change data structure */
   )
{
   assert(domchgdyn != NULL);

   return domchgdyn->domchg;
}


/*
 * methods for variables 
 */

/** creates a variable */
static
RETCODE varCreate(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, var) );

   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*var)->name, name, strlen(name)+1) );
   (*var)->parentvars = NULL;
   (*var)->eventfilter = NULL;
   (*var)->dom.holelist = NULL;
   (*var)->dom.lb = lb;
   (*var)->dom.ub = ub;
   (*var)->obj = obj;
   (*var)->index = stat->nvaridx++;
   (*var)->probindex = -1;
   (*var)->pseudocandindex = -1;
   (*var)->eventqueueindexlb = -1;
   (*var)->eventqueueindexub = -1;
   (*var)->parentvarssize = 0;
   (*var)->nparentvars = 0;
   (*var)->nuses = 0;
   (*var)->nlocksdown = 0;
   (*var)->nlocksup = 0;
   (*var)->vartype = vartype;

   return SCIP_OKAY;
}

/** creates and captures an original problem variable */
RETCODE SCIPvarCreate(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

   /* create variable */
   CHECK_OKAY( varCreate(var, memhdr, set, stat, name, lb, ub, obj, vartype) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_ORIGINAL;
   (*var)->data.transvar = NULL;

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;
}

/** creates and captures a loose variable belonging to the transformed problem */
RETCODE SCIPvarCreateTransformed(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);

   /* create variable */
   CHECK_OKAY( varCreate(var, memhdr, set, stat, name, lb, ub, obj, vartype) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_LOOSE;

   /* create event filter for transformed variable */
   CHECK_OKAY( SCIPeventfilterCreate(&(*var)->eventfilter, memhdr) );

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;   
}

/** frees a variable */
static
RETCODE varFree(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (may be NULL, if it's not a column variable) */
   )
{
   VAR* parentvar;
   int i;

   assert(memhdr != NULL);
   assert(var != NULL);
   assert(*var != NULL);
   assert((*var)->varstatus != SCIP_VARSTATUS_COLUMN || (VAR**)(&(*var)->data.col->var) != var);
   assert((*var)->nuses == 0);
   assert((*var)->probindex == -1);

   debugMessage("free variable <%s> with status=%d\n", (*var)->name, (*var)->varstatus);
   switch( (*var)->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert((*var)->nparentvars == 0); /* original variable cannot have a parent variable */
      assert((*var)->data.transvar == NULL);  /* cannot free variable, if transformed variable is still existing */
      break;
   case SCIP_VARSTATUS_LOOSE:
      break;
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( SCIPcolFree(&(*var)->data.col, memhdr, set, lp) );  /* free corresponding LP column */
      break;
   case SCIP_VARSTATUS_FIXED:
      break;
   case SCIP_VARSTATUS_AGGREGATED:
      break;
   case SCIP_VARSTATUS_MULTAGGR:
      assert((*var)->data.multaggr.vars != NULL);
      assert((*var)->data.multaggr.scalars != NULL);
      assert((*var)->data.multaggr.nvars >= 2);
      freeBlockMemoryArray(memhdr, &(*var)->data.multaggr.vars, (*var)->data.multaggr.nvars);
      freeBlockMemoryArray(memhdr, &(*var)->data.multaggr.scalars, (*var)->data.multaggr.nvars);
      break;
   default:
      errorMessage("Unknown variable status");
      abort();
   }

   /* remove the possibly existing links in the aggregation tree */
   for( i = 0; i < (*var)->nparentvars; ++i )
   {
      assert((*var)->varstatus != SCIP_VARSTATUS_ORIGINAL);
      assert((*var)->parentvars != NULL);
      parentvar = (*var)->parentvars[i];
      assert(parentvar != NULL);

      switch( parentvar->varstatus )
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
      case SCIP_VARSTATUS_MULTAGGR:
         todoMessage("remove variable from multiple aggregation");
         errorMessage("remove variable from multiple aggregation not implemented yet");
         break;
      default:
         errorMessage("parent variable is neither ORIGINAL, AGGREGATED or MULTAGGR");
         return SCIP_INVALIDDATA;
      }
   }

   /* free parentvars array */
   freeBlockMemoryArrayNull(memhdr, &(*var)->parentvars, (*var)->parentvarssize);

   /* free event filter */
   if( (*var)->varstatus != SCIP_VARSTATUS_ORIGINAL )
   {
      CHECK_OKAY( SCIPeventfilterFree(&(*var)->eventfilter, memhdr, set) );
   }

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
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
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

/** adds a hole to the variables domain */
RETCODE SCIPvarAddHole(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(var != NULL);

   CHECK_OKAY( domAddHole(&var->dom, memhdr, set, left, right) );

   return SCIP_OKAY;
}

/** ensures, that parentvars array of var can store at least num entries */
static
RETCODE varEnsureParentvarsSize(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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

/** adds variable to parent list of a variable */
static
RETCODE varAddParent(
   VAR*             var,                /**< variable to add parent to */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
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

   return SCIP_OKAY;
}

/** copies original variable into loose transformed variable, that is captured */
RETCODE SCIPvarTransform(
   VAR**            transvar,           /**< pointer to store the transformed variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR*             origvar             /**< original problem variable */
   )
{
   char name[255];

   assert(origvar != NULL);
   assert(origvar->varstatus == SCIP_VARSTATUS_ORIGINAL);
   assert(origvar->data.transvar == NULL);
   assert(transvar != NULL);

   sprintf(name, "t_%s", origvar->name);
   CHECK_OKAY( SCIPvarCreateTransformed(transvar, memhdr, set, stat,
                  name, origvar->dom.lb, origvar->dom.ub, objsense * origvar->obj, origvar->vartype) );

   CHECK_OKAY( holelistDuplicate(&(*transvar)->dom.holelist, memhdr, set, origvar->dom.holelist) );

   /* link original and transformed variable */
   origvar->data.transvar = *transvar;
   CHECK_OKAY( varAddParent(*transvar, memhdr, set, origvar) );

   /* copy rounding locks */
   (*transvar)->nlocksdown = origvar->nlocksdown;
   (*transvar)->nlocksup = origvar->nlocksup;

   debugMessage("transformed variable: <%s>[%p] -> <%s>[%p]\n", origvar->name, origvar, (*transvar)->name, *transvar);

   return SCIP_OKAY;
}

/** converts loose transformed variable into column variable, creates LP column */
RETCODE SCIPvarColumn(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_LOOSE);

   var->varstatus = SCIP_VARSTATUS_COLUMN;

   CHECK_OKAY( SCIPcolCreate(&var->data.col, memhdr, set, stat, lp, var, 0, NULL, NULL) );

   return SCIP_OKAY;
}

/** converts variable into fixed variable, updates LP respectively */
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data (needed, if variable has a column attached) */
   STAT*            stat,               /**< problem statistics */
   Real             fixedval            /**< value to fix variable at */
   )
{
   assert(var != NULL);
   assert(SCIPsetIsLE(set, var->dom.lb, fixedval));
   assert(SCIPsetIsLE(set, fixedval, var->dom.ub));

   errorMessage("fixing of variables not implemented yet");
   abort();

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      errorMessage("Cannot fix an original variable");
      return SCIP_INVALIDDATA;
   case SCIP_VARSTATUS_LOOSE:
      todoMessage("apply fixings to problem: objoffset, move var from vars to fixedvars array");
      abort();
      break;
   case SCIP_VARSTATUS_COLUMN:
      /* fix column */
      todoMessage("apply fixings to LP (change rows, offset obj)");
      todoMessage("apply fixings to problem: objoffset, move var from vars to fixedvars array");
      abort();
      break;
   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot fix a fixed variable again");
      return SCIP_INVALIDDATA;
   case SCIP_VARSTATUS_AGGREGATED:
      todoMessage("fix aggregation variable of aggregated variable, transform aggregated to fixed");
      abort();
   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot fix a multiple aggregated variable");
      return SCIP_INVALIDDATA;
   default:
      errorMessage("Unknown variable status");
      abort();
   }
   
   var->varstatus = SCIP_VARSTATUS_FIXED;
   var->dom.lb = fixedval;
   var->dom.ub = fixedval;
   holelistFree(&var->dom.holelist, memhdr);

   return SCIP_OKAY;
}

/** converts variable into aggregated variable, updates LP respectively */
RETCODE SCIPvarAggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data (needed, if variable has a column attached) */
   STAT*            stat,               /**< problem statistics */
   VAR*             aggvar,             /**< variable $y$ in aggregation $x = a*y + c$ */
   Real             scalar,             /**< multiplier $a$ in aggregation $x = a*y + c$ */
   Real             constant            /**< constant shift $c$ in aggregation $x = a*y + c$ */
   )
{
   errorMessage("aggregation of variables not yet implemented");
   abort();

   todoMessage("don't forget to move the rounding locks to the aggregation variable");

   return SCIP_OKAY;
}

/** changes type of variable; cannot be called, if var belongs to a problem */
RETCODE SCIPvarChgType(
   VAR*             var,                /**< problem variable to change */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(var != NULL);

   if( var->probindex >= 0 )
   {
      errorMessage("cannot change type of variable already in the problem");
      return SCIP_INVALIDDATA;
   }
   
   var->vartype = vartype;

   return SCIP_OKAY;
}

/** increases lock number for rounding down; tells variable, that rounding its value down will make the
 *  solution infeasible
 */
void SCIPvarForbidRoundDown(
   VAR*             var                 /**< problem variable */
   )
{
   int i;

   assert(var != NULL);
   assert(var->nlocksdown >= 0);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarForbidRoundDown(var->data.transvar);
      else
         var->nlocksdown++;
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      var->nlocksdown++;
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         SCIPvarForbidRoundDown(var->data.aggregate.var);
      else
         SCIPvarForbidRoundUp(var->data.aggregate.var);
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            SCIPvarForbidRoundDown(var->data.multaggr.vars[i]);
         else
            SCIPvarForbidRoundUp(var->data.multaggr.vars[i]);
      }
      break;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** increases lock number for rounding up; tells variable, that rounding its value up will make the solution infeasible */
void SCIPvarForbidRoundUp(
   VAR*             var                 /**< problem variable */
   )
{
   int i;

   assert(var != NULL);
   assert(var->nlocksup >= 0);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarForbidRoundUp(var->data.transvar);
      else
         var->nlocksup++;
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      var->nlocksup++;
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         SCIPvarForbidRoundUp(var->data.aggregate.var);
      else
         SCIPvarForbidRoundDown(var->data.aggregate.var);
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            SCIPvarForbidRoundUp(var->data.multaggr.vars[i]);
         else
            SCIPvarForbidRoundDown(var->data.multaggr.vars[i]);
      }
      break;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/**< increases lock number for rounding down and up; tells variable, that rounding value in either direction will
 *   make the solution infeasible
 */
void SCIPvarForbidRound(
   VAR*             var                 /**< problem variable */
   )
{
   SCIPvarForbidRoundDown(var);
   SCIPvarForbidRoundUp(var);
}

/** decreases lock number for rounding down; cancels a prior forbidRoundDown() */
void SCIPvarAllowRoundDown(
   VAR*             var                 /**< problem variable */
   )
{
   int i;

   assert(var != NULL);
   assert(var->nlocksdown >= 0);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarAllowRoundDown(var->data.transvar);
      else
         var->nlocksdown--;
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      var->nlocksdown--;
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         SCIPvarAllowRoundDown(var->data.aggregate.var);
      else
         SCIPvarAllowRoundUp(var->data.aggregate.var);
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            SCIPvarAllowRoundDown(var->data.multaggr.vars[i]);
         else
            SCIPvarAllowRoundUp(var->data.multaggr.vars[i]);
      }
      break;

   default:
      errorMessage("unknown variable status");
      abort();
   }

   assert(var->nlocksdown >= 0);
}

/** decreases lock number for rounding up; cancels a prior forbidRoundUp() */
void SCIPvarAllowRoundUp(
   VAR*             var                 /**< problem variable */
   )
{
   int i;

   assert(var != NULL);
   assert(var->nlocksup >= 0);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
         SCIPvarAllowRoundUp(var->data.transvar);
      else
         var->nlocksup--;
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_FIXED:
      var->nlocksup--;
      break;

   case SCIP_VARSTATUS_AGGREGATED:
      if( var->data.aggregate.scalar > 0.0 )
         SCIPvarAllowRoundUp(var->data.aggregate.var);
      else
         SCIPvarAllowRoundDown(var->data.aggregate.var);
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         if( var->data.multaggr.scalars[i] > 0.0 )
            SCIPvarAllowRoundUp(var->data.multaggr.vars[i]);
         else
            SCIPvarAllowRoundDown(var->data.multaggr.vars[i]);
      }
      break;

   default:
      errorMessage("unknown variable status");
      abort();
   }

   assert(var->nlocksup >= 0);
}

/** decreases lock number for rounding down & up; cancels a prior forbidRound() */
void SCIPvarAllowRound(
   VAR*             var                 /**< problem variable */
   )
{
   SCIPvarAllowRoundDown(var);
   SCIPvarAllowRoundUp(var);
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

   switch( var->varstatus )
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
      break;

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

   default:
      errorMessage("unknown variable status");
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

   switch( var->varstatus )
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
      break;

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

   default:
      errorMessage("unknown variable status");
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

/** appends LBTIGHTENED or LBRELAXED event to the event queue */
static
RETCODE varEventLbChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             oldbound,           /**< old lower bound for variable */
   Real             newbound            /**< new lower bound for variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert(!SCIPsetIsEQ(set, oldbound, newbound));

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because row activities and LP changes have to be updated
    */
   if( var->eventfilter->len > 0 || var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;
   
      CHECK_OKAY( SCIPeventCreateLbChanged(&event, memhdr, var, oldbound, newbound) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, tree, lp, branchcand, &event) );
   }

   return SCIP_OKAY;
}

/** appends UBTIGHTENED or UBRELAXED event to the event queue */
static
RETCODE varEventUbChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             oldbound,           /**< old upper bound for variable */
   Real             newbound            /**< new upper bound for variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert(!SCIPsetIsEQ(set, oldbound, newbound));

   /* check, if the variable is being tracked for bound changes
    * COLUMN and LOOSE variables are tracked always, because row activities and LP changes have to be updated
    */
   if( var->eventfilter->len > 0 || var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;
   
      CHECK_OKAY( SCIPeventCreateUbChanged(&event, memhdr, var, oldbound, newbound) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, tree, lp, branchcand, &event) );
   }

   return SCIP_OKAY;
}

/* forward declaration, because both methods call each other recursively */

/** performs the actual change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   );

/** performs the actual change in lower bound, changes all parents accordingly */
static            
RETCODE varProcessChgLb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);

   debugMessage("process changing lower bound of <%s> from %g to %g\n", var->name, var->dom.lb, newbound);

   /* change the bound */
   oldbound = var->dom.lb;
   var->dom.lb = newbound;

#if 0
   /* inform LP and tree about bound change */
   if( var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      if( var->varstatus == SCIP_VARSTATUS_COLUMN )
      {
         CHECK_OKAY( SCIPcolBoundChanged(var->data.col, memhdr, set, lp, SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      }
      CHECK_OKAY( SCIPtreeBoundChanged(tree, memhdr, set, var, SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
   }
#endif

   /* issue bound change event */
   CHECK_OKAY( varEventLbChanged(var, memhdr, set, tree, lp, branchcand, eventqueue, oldbound, newbound) );

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( parentvar->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multaggr variable cannot be the parent of a variable");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPos(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            CHECK_OKAY( varProcessChgLb(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else if( SCIPsetIsNeg(set, var->data.aggregate.scalar) )
         {
            /* a < 0 -> change upper bound of y */
            CHECK_OKAY( varProcessChgUb(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else
         {
            errorMessage("scalar is zero in aggregation");
            return SCIP_INVALIDDATA;
         }
         break;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** performs the actual change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);

   debugMessage("process changing upper bound of <%s> from %g to %g\n", var->name, var->dom.ub, newbound);

   /* change the bound */
   oldbound = var->dom.ub;
   var->dom.ub = newbound;

#if 0
   /* inform LP and tree about bound change */
   if( var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      if( var->varstatus == SCIP_VARSTATUS_COLUMN )
      {
         CHECK_OKAY( SCIPcolBoundChanged(var->data.col, memhdr, set, lp, SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      }
      CHECK_OKAY( SCIPtreeBoundChanged(tree, memhdr, set, var, SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
   }
#endif

   /* issue bound change event */
   CHECK_OKAY( varEventUbChanged(var, memhdr, set, tree, lp, branchcand, eventqueue, oldbound, newbound) );

   /* process parent variables */
   for( i = 0; i < var->nparentvars; ++i )
   {
      parentvar = var->parentvars[i];
      assert(parentvar != NULL);

      switch( parentvar->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multaggr variable cannot be the parent of a variable");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPos(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of y */
            CHECK_OKAY( varProcessChgUb(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else if( SCIPsetIsNeg(set, var->data.aggregate.scalar) )
         {
            /* a < 0 -> change lower bound of y */
            CHECK_OKAY( varProcessChgLb(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else
         {
            errorMessage("scalar is zero in aggregation");
            return SCIP_INVALIDDATA;
         }
         break;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes lower bound of variable */
RETCODE SCIPvarChgLb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing lower bound of <%s> from %g to %g\n", var->name, var->dom.lb, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( !SCIPsetIsEQ(set, var->dom.lb, newbound) )
   {
      /* change bounds of attached variables */
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarChgLb(var->data.transvar, memhdr, set, stat, lp, tree, branchcand, eventqueue, newbound) );
         }
         else
         {
            assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
            stat->nboundchanges++;
            var->dom.lb = newbound;
         }
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
         stat->nboundchanges++;
         CHECK_OKAY( varProcessChgLb(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, newbound) );
         break;

      case SCIP_VARSTATUS_FIXED:
         errorMessage("cannot change the bounds of a fixed variable");
         return SCIP_INVALIDDATA;
         
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(var->data.aggregate.var != NULL);
         if( SCIPsetIsPos(set, var->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            CHECK_OKAY( SCIPvarChgLb(var->data.aggregate.var, memhdr, set, stat, lp, tree, branchcand, eventqueue, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else if( SCIPsetIsNeg(set, var->data.aggregate.scalar) )
         {
            /* a < 0 -> change upper bound of y */
            CHECK_OKAY( SCIPvarChgUb(var->data.aggregate.var, memhdr, set, stat, lp, tree, branchcand, eventqueue, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else
         {
            errorMessage("scalar is zero in aggregation");
            return SCIP_INVALIDDATA;
         }
         break;
         
      case SCIP_VARSTATUS_MULTAGGR:
         todoMessage("change the sides of the corresponding linear constraint");
         errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet");
         abort();

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes upper bound of variable */
RETCODE SCIPvarChgUb(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing upper bound of <%s> from %g to %g\n", var->name, var->dom.ub, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( !SCIPsetIsEQ(set, var->dom.ub, newbound) )
   {
      /* change bounds of attached variables */
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarChgUb(var->data.transvar, memhdr, set, stat, lp, tree, branchcand, eventqueue, newbound) );
         }
         else
         {
            assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
            stat->nboundchanges++;
            var->dom.ub = newbound;
         }
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
         stat->nboundchanges++;
         CHECK_OKAY( varProcessChgUb(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, newbound) );
         break;

      case SCIP_VARSTATUS_FIXED:
         errorMessage("cannot change the bounds of a fixed variable");
         return SCIP_INVALIDDATA;
         
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(var->data.aggregate.var != NULL);
         if( SCIPsetIsPos(set, var->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of y */
            CHECK_OKAY( SCIPvarChgUb(var->data.aggregate.var, memhdr, set, stat, lp, tree, branchcand, eventqueue, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else if( SCIPsetIsNeg(set, var->data.aggregate.scalar) )
         {
            /* a < 0 -> change lower bound of y */
            CHECK_OKAY( SCIPvarChgLb(var->data.aggregate.var, memhdr, set, stat, lp, tree, branchcand, eventqueue, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else
         {
            errorMessage("scalar is zero in aggregation");
            return SCIP_INVALIDDATA;
         }
         break;
         
      case SCIP_VARSTATUS_MULTAGGR:
         todoMessage("change the sides of the corresponding linear constraint");
         errorMessage("changing the bounds of a multiple aggregated variable is not implemented yet");
         abort();

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes bound of variable */
RETCODE SCIPvarChgBd(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
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
      return SCIPvarChgLb(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUb(var, memhdr, set, stat, lp, tree, branchcand, eventqueue, newbound);
   default:
      errorMessage("Unknown bound type");
      return SCIP_INVALIDDATA;
   }
}

/** adjust lower bound to integral value, if variable is integral */
void SCIPvarAdjustLb(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            lb                  /**< pointer to lower bound to adjust */
   )
{
   assert(var != NULL);
   assert(lb != NULL);

   if( var->vartype != SCIP_VARTYPE_CONTINOUS )
   {
      /* adjust new bound to integral value */
      *lb = SCIPsetCeil(set, *lb);
   }
}

/** adjust upper bound to integral value, if variable is integral */
void SCIPvarAdjustUb(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            ub                  /**< pointer to upper bound to adjust */
   )
{
   assert(var != NULL);
   assert(ub != NULL);

   if( var->vartype != SCIP_VARTYPE_CONTINOUS )
   {
      /* adjust new bound to integral value */
      *ub = SCIPsetFloor(set, *ub);
   }
}

/** changes objective value of variable */
RETCODE SCIPvarChgObj(
   VAR*             var,                /**< variable to change, must not be member of the problem */
   Real             newobj              /**< new objective value for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing objective value of <%s> from %g to %g\n", var->name, var->obj, newobj);

   if( var->probindex >= 0 )
   {
      errorMessage("cannot change the objective value of a variable already in the problem");
      return SCIP_INVALIDDATA;
   }
   
   var->obj = newobj;
      
   return SCIP_OKAY;
}

/** get name of variable */
const char* SCIPvarGetName(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->name;
}

/** gets status of variable */
VARSTATUS SCIPvarGetStatus(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->varstatus;
}

/** gets type of variable */
VARTYPE SCIPvarGetType(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->vartype;
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
int SCIPvarGetProbIndex(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->probindex;
}

/** gets corresponding transformed variable of an original variable */
VAR* SCIPvarGetTransformed(
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

/** gets objective function value of variable */
Real SCIPvarGetObj(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->obj;
}
   
/** gets lower bound of variable */
Real SCIPvarGetLb(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->dom.lb;
}
   
/** gets upper bound of variable */
Real SCIPvarGetUb(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->dom.ub;
}

/** gets best bound of variable with respect to the objective function */
static
Real varGetBestBound(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->dom.lb;
   else
      return var->dom.ub;
}

/** get primal LP solution value of variable */
Real SCIPvarGetLPSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real primsol;
   int i;

   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetLPSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
      assert(var->dom.lb <= 0.0 && var->dom.ub >= 0.0);
      return 0.0;

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return var->data.col->primsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->dom.lb == var->dom.ub);
      return var->dom.lb;

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

   default:
      errorMessage("Unknown variable status");
      return SCIP_INVALIDDATA;
   }
}

/** get pseudo solution value of variable at actual node */
Real SCIPvarGetPseudoSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real pseudosol;
   int i;

   assert(var != NULL);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetPseudoSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return varGetBestBound(var);

   case SCIP_VARSTATUS_FIXED:
      assert(var->dom.lb == var->dom.ub);
      return var->dom.lb;

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

   default:
      errorMessage("Unknown variable status");
      return SCIP_INVALIDDATA;
   }
}

/**< get solution value of variable at actual node: if LP was solved at the node, the method returns the
 *   LP primal solution value, otherwise the pseudo solution
 */
Real SCIPvarGetSol(
   VAR*             var,                /**< problem variable */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   if( tree->actnodehaslp )
      return SCIPvarGetLPSol(var);
   else
      return SCIPvarGetPseudoSol(var);
}

/** resolves variable to columns and adds them with the coefficient to the row */
RETCODE SCIPvarAddToRow(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   STAT*            stat,               /**< problem statistics */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   int i;

   assert(var != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsInfinity(set, ABS(val)));

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         char s[255];
         sprintf(s, "Cannot add untransformed original variable <%s> to LP row <%s>", var->name, row->name);
         errorMessage(s);
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarAddToRow(var->data.transvar, memhdr, set, lp, stat, row, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
      /* convert loose variable into column */
      CHECK_OKAY( SCIPvarColumn(var, memhdr, set, lp, stat) );
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
      /* fallthrough */

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      assert(var->data.col->var == var);
      CHECK_OKAY( SCIProwIncCoeff(row, memhdr, set, lp, var->data.col, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(var->dom.lb == var->dom.ub);
      assert(!SCIPsetIsInfinity(set, ABS(var->dom.lb)));
      CHECK_OKAY( SCIProwAddConst(row, set, lp, val * var->dom.lb) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      CHECK_OKAY( SCIPvarAddToRow(var->data.aggregate.var, memhdr, set, lp, stat, row, var->data.aggregate.scalar * val) );
      CHECK_OKAY( SCIProwAddConst(row, set, lp, var->data.aggregate.constant * val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         CHECK_OKAY( SCIPvarAddToRow(var->data.multaggr.vars[i], memhdr, set, lp, stat, row, 
                        var->data.multaggr.scalars[i] * val) );
      }
      CHECK_OKAY( SCIProwAddConst(row, set, lp, var->data.multaggr.constant * val) );
      return SCIP_OKAY;

   default:
      errorMessage("Unknown variable status");
      return SCIP_INVALIDDATA;
   }
}

/** includes event handler in variable's event filter */
RETCODE SCIPvarCatchEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert((eventtype & !SCIP_EVENTTYPE_DOMCHANGED) == 0);
   assert((eventtype & SCIP_EVENTTYPE_DOMCHANGED) != 0);

   CHECK_OKAY( SCIPeventfilterAdd(var->eventfilter, memhdr, set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}


/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
DECL_HASHGETKEY(SCIPhashGetKeyVar)
{
   VAR* var = (VAR*)elem;

   assert(var != NULL);
   return var->name;
}

