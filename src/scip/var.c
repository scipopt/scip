/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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

static
RETCODE domchgCreate(                   /**< creates empty fixed size domain change data */
   DOMCHG**         domchg,             /**< pointer to fixed size domain change data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, *domchg) );
   (*domchg)->boundchg = NULL;
   (*domchg)->holechg = NULL;
   (*domchg)->nboundchg = 0;
   (*domchg)->nholechg = 0;

   return SCIP_OKAY;
}

static
RETCODE ensureBoundchgSize(             /**< ensures, that boundchg array can store at least num entries */
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
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, (*domchg)->boundchg, domchgdyn->boundchgsize, newsize) );
      domchgdyn->boundchgsize = newsize;
   }
   assert(num <= domchgdyn->boundchgsize);

   return SCIP_OKAY;
}

static
RETCODE ensureHolechgSize(              /**< ensures, that holechg array can store at least num additional entries */
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
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, (*domchg)->holechg, domchgdyn->holechgsize, newsize) );
      domchgdyn->holechgsize = newsize;
   }
   assert(num <= domchgdyn->holechgsize);

   return SCIP_OKAY;
}



/*
 * hole, holelist, and domain methods
 */

static
RETCODE holelistCreate(                 /**< creates a new holelist element */
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

   ALLOC_OKAY( allocBlockMemory(memhdr, *holelist) );
   (*holelist)->hole.left = left;
   (*holelist)->hole.right = right;
   (*holelist)->next = NULL;

   return SCIP_OKAY;
}

static
void holelistFree(                      /**< frees all elements in the holelist */
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
      freeBlockMemory(memhdr, *holelist);
      *holelist = next;
   }
}

static
RETCODE holelistDuplicate(              /**< duplicates a list of holes */
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

static
RETCODE domAddHole(                     /**< adds a hole to the domain */
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

static
RETCODE domMerge(                       /**< merges overlapping holes into single holes, moves bounds respectively */
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

void SCIPdomchgFree(                    /**< frees fixed size domain change data */
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);

   if( *domchg != NULL )
   {
      freeBlockMemoryArrayNull(memhdr, (*domchg)->boundchg, (*domchg)->nboundchg);
      freeBlockMemoryArrayNull(memhdr, (*domchg)->holechg, (*domchg)->nholechg);
      freeBlockMemory(memhdr, *domchg);
   }
}

RETCODE SCIPdomchgApply(                /**< applies domain change */
   const DOMCHG*    domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
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
         CHECK_OKAY( SCIPvarChgLb(var, memhdr, set, stat, lp, tree, domchg->boundchg[i].newbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUb(var, memhdr, set, stat, lp, tree, domchg->boundchg[i].newbound) );
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
   
RETCODE SCIPdomchgUndo(                 /**< undoes domain change */
   const DOMCHG*    domchg,             /**< domain change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree                /**< branch-and-bound tree */
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
         CHECK_OKAY( SCIPvarChgLb(var, memhdr, set, stat, lp, tree, domchg->boundchg[i].oldbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUb(var, memhdr, set, stat, lp, tree, domchg->boundchg[i].oldbound) );
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

RETCODE SCIPdomchgdynCreate(            /**< creates a dynamic size attachment for a domain change data structure */
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchgdyn != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, *domchgdyn) );

   (*domchgdyn)->domchg = NULL;
   (*domchgdyn)->boundchgsize = 0;
   (*domchgdyn)->holechgsize = 0;

   return SCIP_OKAY;
}

void SCIPdomchgdynFree(                 /**< frees a dynamic size attachment for a domain change data structure */
   DOMCHGDYN**      domchgdyn,          /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(domchgdyn != NULL);
   assert(*domchgdyn != NULL);
   assert((*domchgdyn)->domchg == NULL);
   assert((*domchgdyn)->boundchgsize == 0);
   assert((*domchgdyn)->holechgsize == 0);

   freeBlockMemory(memhdr, *domchgdyn);
}

void SCIPdomchgdynAttach(               /**< attaches dynamic size information to domain change data */
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

RETCODE SCIPdomchgdynDetach(            /**< detaches dynamic size information and shrinks domain change data */
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
         freeBlockMemoryArrayNull(memhdr, (*domchgdyn->domchg)->boundchg, domchgdyn->boundchgsize);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, (*domchgdyn->domchg)->boundchg,
                        domchgdyn->boundchgsize, (*domchgdyn->domchg)->nboundchg) );
      }
      if( (*domchgdyn->domchg)->nholechg == 0 )
      {
         freeBlockMemoryArrayNull(memhdr, (*domchgdyn->domchg)->holechg, domchgdyn->holechgsize);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, (*domchgdyn->domchg)->holechg,
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

void SCIPdomchgdynDiscard(              /**< frees attached domain change data and detaches dynamic size attachment */
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
      freeBlockMemoryArrayNull(memhdr, (*domchgdyn->domchg)->boundchg, domchgdyn->boundchgsize);
      (*domchgdyn->domchg)->nboundchg = 0;
      freeBlockMemoryArrayNull(memhdr, (*domchgdyn->domchg)->holechg, domchgdyn->holechgsize);
      (*domchgdyn->domchg)->nholechg = 0;
      SCIPdomchgFree(domchgdyn->domchg, memhdr);
      domchgdyn->boundchgsize = 0;
      domchgdyn->holechgsize = 0;
   }

   /* detach domain change data */
   domchgdyn->domchg = NULL;
}

RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
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

   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg != NULL);
   assert(var != NULL);

   debugMessage("(1) adding bound change <%s>: %g -> %g of variable <%s> to domain change at %p pointing to %p\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", oldbound, newbound, var->name, domchgdyn->domchg,
      *domchgdyn->domchg);
   CHECK_OKAY( ensureBoundchgSize(domchgdyn, memhdr, set, domchgdyn->boundchgsize+1) );
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

RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   DOMCHG* domchg;

   assert(domchgdyn != NULL);
   assert(domchgdyn->domchg != NULL);
   assert(ptr != NULL);

   CHECK_OKAY( ensureHolechgSize(domchgdyn, memhdr, set, domchgdyn->holechgsize+1) );

   domchg = *domchgdyn->domchg;
   assert(domchg != NULL);

   domchg->holechg[domchg->nholechg].ptr = ptr;
   domchg->holechg[domchg->nholechg].newlist = newlist;
   domchg->holechg[domchg->nholechg].oldlist = oldlist;
   domchg->nholechg++;

   return SCIP_OKAY;
}

DOMCHG** SCIPdomchgdynGetDomchgPtr(     /**< gets pointer to domain change data the dynamic size information references */
   DOMCHGDYN*       domchgdyn           /**< dynamically sized domain change data structure */
   )
{
   assert(domchgdyn != NULL);

   return domchgdyn->domchg;
}


/*
 * methods for variables 
 */

RETCODE SCIPvarCreate(                  /**< creates and captures an original problem variable */
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

   ALLOC_OKAY( allocBlockMemory(memhdr, *var) );

   (*var)->data.transvar = NULL;
   (*var)->origvar = NULL;
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*var)->name, name, strlen(name)+1) );
   (*var)->dom.holelist = NULL;
   (*var)->dom.lb = lb;
   (*var)->dom.ub = ub;
   (*var)->obj = obj;
   (*var)->index = stat->nvaridx++;
   (*var)->probindex = -1;
   (*var)->nuses = 0;
   (*var)->vartype = vartype;
   (*var)->varstatus = SCIP_VARSTATUS_ORIGINAL;

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;
}

RETCODE SCIPvarCreateTransformed(       /**< creates and captures a loose variable belonging to the transformed problem */
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

   ALLOC_OKAY( allocBlockMemory(memhdr, *var) );

   (*var)->origvar = NULL;
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*var)->name, name, strlen(name)+1) );
   (*var)->dom.holelist = NULL;
   (*var)->dom.lb = lb;
   (*var)->dom.ub = ub;
   (*var)->obj = obj;
   (*var)->index = stat->nvaridx++;
   (*var)->probindex = -1;
   (*var)->nuses = 0;
   (*var)->vartype = vartype;
   (*var)->varstatus = SCIP_VARSTATUS_LOOSE;

   /* capture variable */
   SCIPvarCapture(*var);

   return SCIP_OKAY;   
}

RETCODE SCIPvarFree(                    /**< frees a variable */
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (may be NULL, if it's not a column variable) */
   )
{
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
      assert((*var)->origvar == NULL); /* original variable cannot have link to other original variable */
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
      freeBlockMemoryArray(memhdr, (*var)->data.multaggr.vars, (*var)->data.multaggr.nvars);
      freeBlockMemoryArray(memhdr, (*var)->data.multaggr.scalars, (*var)->data.multaggr.nvars);
      break;
   default:
      errorMessage("Unknown variable status");
      abort();
   }

   /* remove the possibly existing link between transformed and original var */
   if( (*var)->origvar != NULL )
   {
      assert((*var)->varstatus != SCIP_VARSTATUS_ORIGINAL);
      assert((*var)->origvar->varstatus == SCIP_VARSTATUS_ORIGINAL);
      assert((*var)->origvar->data.transvar == *var);
      assert(&(*var)->origvar->data.transvar != var);
      (*var)->origvar->data.transvar = NULL;
   }

   freeBlockMemoryArray(memhdr, (*var)->name, strlen((*var)->name)+1);
   freeBlockMemory(memhdr, *var);

   return SCIP_OKAY;
}

void SCIPvarCapture(                    /**< increases usage counter of variable */
   VAR*             var                 /**< variable */
   )
{
   assert(var != NULL);
   assert(var->nuses >= 0);

   debugMessage("capture variable <%s> with nuses=%d\n", var->name, var->nuses);
   var->nuses++;
}

RETCODE SCIPvarRelease(                 /**< decreases usage counter of variable, and frees memory if necessary */
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
      CHECK_OKAY( SCIPvarFree(var, memhdr, set, lp) );
   }

   *var = NULL;

   return SCIP_OKAY;
}

RETCODE SCIPvarAddHole(                 /**< adds a hole to the variables domain */
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

RETCODE SCIPvarTransform(               /**< copies original variable into loose transformed variable, that is captured */
   VAR**            transvar,           /**< pointer to store the transformed variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR*             origvar             /**< original problem variable */
   )
{
   assert(origvar != NULL);
   assert(origvar->varstatus == SCIP_VARSTATUS_ORIGINAL);
   assert(origvar->data.transvar == NULL);
   assert(transvar != NULL);

   CHECK_OKAY( SCIPvarCreateTransformed(transvar, memhdr, set, stat,
                  origvar->name, origvar->dom.lb, origvar->dom.ub, objsense * origvar->obj, origvar->vartype) );

   CHECK_OKAY( holelistDuplicate(&(*transvar)->dom.holelist, memhdr, set, origvar->dom.holelist) );

   /* link original and transformed variable */
   origvar->data.transvar = *transvar;
   (*transvar)->origvar = origvar;

   debugMessage("transformed variable: %s[%p] -> %s[%p]\n", origvar->name, origvar, (*transvar)->name, *transvar);

   return SCIP_OKAY;
}

RETCODE SCIPvarColumn(                  /**< converts loose transformed variable into column variable, creates LP column */
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

RETCODE SCIPvarFix(                     /**< converts variable into fixed variable, updates LP respectively */
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

RETCODE SCIPvarAggregate(               /**< converts variable into aggregated variable, updates LP respectively */
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
   errorMessage("Method not yet implemented");
   abort();

   return SCIP_OKAY;
}

RETCODE SCIPvarChgType(                 /**< changes type of variable; cannot be called, if var belongs to a problem */
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

RETCODE SCIPvarChgLb(                   /**< changes lower bound of variable */
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing lower bound of <%s> from %g to %g\n", var->name, var->dom.lb, newbound);

   if( !SCIPsetIsEQ(set, var->dom.lb, newbound) )
   {
      Real oldbound;
   
      stat->nboundchanges++;

      oldbound = var->dom.lb;
      var->dom.lb = newbound;
      
      /* change bounds of attached variables */
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarChgLb(var->data.transvar, memhdr, set, stat, lp, tree, newbound) );
         }
         break;
         
      case SCIP_VARSTATUS_COLUMN:
         /* notifiy LP of the bound change */
         CHECK_OKAY( SCIPcolBoundChanged(var->data.col, memhdr, set, lp, SCIP_BOUNDTYPE_LOWER) );
         /* fallthrough */
         
      case SCIP_VARSTATUS_LOOSE:
         CHECK_OKAY( SCIPtreeBoundChanged(tree, memhdr, set, var, SCIP_BOUNDTYPE_LOWER, oldbound) );
         break;

      case SCIP_VARSTATUS_FIXED:
         errorMessage("cannot change the bounds of a fixed variable");
         return SCIP_INVALIDDATA;
         
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(var->data.aggregate.var != NULL);
         if( SCIPsetIsPos(set, var->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            CHECK_OKAY( SCIPvarChgLb(var->data.aggregate.var, memhdr, set, stat, lp, tree, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else if( SCIPsetIsNeg(set, var->data.aggregate.scalar) )
         {
            /* a < 0 -> change upper bound of y */
            CHECK_OKAY( SCIPvarChgUb(var->data.aggregate.var, memhdr, set, stat, lp, tree, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else
         {
            errorMessage("scalar is zero in aggregation");
            return SCIP_INVALIDDATA;
         }
         break;
         
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("cannot change the bounds of a multiple aggregated variable");
         return SCIP_INVALIDDATA;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

RETCODE SCIPvarChgUb(                   /**< changes upper bound of variable */
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing upper bound of <%s> from %g to %g\n", var->name, var->dom.ub, newbound);

   if( !SCIPsetIsEQ(set, var->dom.ub, newbound) )
   {
      Real oldbound;

      stat->nboundchanges++;

      oldbound = var->dom.ub;
      var->dom.ub = newbound;
      
      /* change bounds of attached variables */
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarChgUb(var->data.transvar, memhdr, set, stat, lp, tree, newbound) );
         }
         break;
         
      case SCIP_VARSTATUS_COLUMN:
         /* notifiy LP of the bound change */
         CHECK_OKAY( SCIPcolBoundChanged(var->data.col, memhdr, set, lp, SCIP_BOUNDTYPE_UPPER) );
         /* fallthrough */
         
      case SCIP_VARSTATUS_LOOSE:
         CHECK_OKAY( SCIPtreeBoundChanged(tree, memhdr, set, var, SCIP_BOUNDTYPE_UPPER, oldbound) );
         break;

      case SCIP_VARSTATUS_FIXED:
         errorMessage("cannot change the bounds of a fixed variable");
         return SCIP_INVALIDDATA;
         
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(var->data.aggregate.var != NULL);
         if( SCIPsetIsPos(set, var->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of y */
            CHECK_OKAY( SCIPvarChgUb(var->data.aggregate.var, memhdr, set, stat, lp, tree, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else if( SCIPsetIsNeg(set, var->data.aggregate.scalar) )
         {
            /* a < 0 -> change lower bound of y */
            CHECK_OKAY( SCIPvarChgLb(var->data.aggregate.var, memhdr, set, stat, lp, tree, 
                           (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
         }
         else
         {
            errorMessage("scalar is zero in aggregation");
            return SCIP_INVALIDDATA;
         }
         break;
         
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("cannot change the bounds of a multiple aggregated variable");
         return SCIP_INVALIDDATA;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

RETCODE SCIPvarChgBd(                   /**< changes bound of variable */
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   /* apply bound change to the LP data */
   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      return SCIPvarChgLb(var, memhdr, set, stat, lp, tree, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUb(var, memhdr, set, stat, lp, tree, newbound);
   default:
      errorMessage("Unknown bound type");
      return SCIP_INVALIDDATA;
   }
}

RETCODE SCIPvarChgObj(                  /**< changes objective value of variable */
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

const char* SCIPvarGetName(             /**< get name of variable */
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->name;
}

Real SCIPvarGetLb(                      /**< gets lower bound of variable */
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->dom.lb;
}

Real SCIPvarGetUb(                      /**< gets upper bound of variable */
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->dom.ub;
}

static
Real varGetBestBound(                   /**< gets best bound of variable with respect to the objective function */
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( var->obj >= 0.0 )
      return var->dom.lb;
   else
      return var->dom.ub;
}

Real SCIPvarGetLPSol(                   /**< get primal LP solution value of variable */
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

Real SCIPvarGetPseudoSol(               /**< get pseudo solution value of variable at actual node */
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

Real SCIPvarGetSol(                     /**< get solution value of variable at actual node: if LP was solved at the node,
                                           the method returns the LP primal solution value, otherwise the pseudo solution */
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

RETCODE SCIPvarAddToRow(                /**< resolves variable to columns and adds them with the coefficient to the row */
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



/*
 * Hash functions
 */

DECL_HASHGETKEY(SCIPhashGetKeyVar)      /**< gets the key (i.e. the name) of the given variable */
{
   VAR* var = (VAR*)elem;

   assert(var != NULL);
   return var->name;
}

