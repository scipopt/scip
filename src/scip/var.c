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


/*
 * domain change methods
 */

/** frees fixed size domain change data */
RETCODE SCIPdomchgFree(
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

   return SCIP_OKAY;
}

/** applies domain change */
RETCODE SCIPdomchgApply(
   const DOMCHG*    domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
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
         CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        domchg->boundchg[i].newbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        domchg->boundchg[i].newbound) );
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
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
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
         CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        domchg->boundchg[i].oldbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        domchg->boundchg[i].oldbound) );
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
      {
         CHECK_OKAY( SCIPdomchgFree(domchgdyn->domchg, memhdr) );
      }
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
RETCODE SCIPdomchgdynDiscard(
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
      CHECK_OKAY( SCIPdomchgFree(domchgdyn->domchg, memhdr) );
      domchgdyn->boundchgsize = 0;
      domchgdyn->holechgsize = 0;
   }

   /* detach domain change data */
   domchgdyn->domchg = NULL;

   return SCIP_OKAY;
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
   assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

   debugMessage("adding bound change <%s>: %g -> %g of variable <%s> to domain change at %p pointing to %p\n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", oldbound, newbound, var->name, domchgdyn->domchg,
      *domchgdyn->domchg);

   nboundchg = (*domchgdyn->domchg == NULL ? 0 : (*domchgdyn->domchg)->nboundchg);
   CHECK_OKAY( ensureBoundchgSize(domchgdyn, memhdr, set, nboundchg+1) );

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
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

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
   (*var)->eventfilter = NULL;
   (*var)->glbdom.holelist = NULL;
   (*var)->glbdom.lb = lb;
   (*var)->glbdom.ub = ub;
   (*var)->actdom.holelist = NULL;
   (*var)->actdom.lb = lb;
   (*var)->actdom.ub = ub;
   (*var)->obj = obj;
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
   (*var)->vartype = vartype;
   (*var)->removeable = removeable;

   return SCIP_OKAY;
}

/** creates and captures an original problem variable */
RETCODE SCIPvarCreate(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

   /* create variable */
   CHECK_OKAY( varCreate(var, memhdr, set, stat, name, lb, ub, obj, vartype, removeable) );

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
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
   )
{
   assert(var != NULL);
   assert(memhdr != NULL);

   /* create variable */
   CHECK_OKAY( varCreate(var, memhdr, set, stat, name, lb, ub, obj, vartype, removeable) );

   /* set variable status and data */
   (*var)->varstatus = SCIP_VARSTATUS_LOOSE;

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

/** adds variable to parent list of a variable and captures parent variable */
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

   SCIPvarCapture(parentvar);

   return SCIP_OKAY;
}

/** deletes and releases all variables from the parent list of a variable, frees the memory of parents array */
static
RETCODE varFreeParents(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
   )
{
   VAR* parentvar;
   int i;
   int v;

   debugMessage("free parents of <%s>\n", (*var)->name);

   /* release the parent variables and remove the link from the parent variable to the child */
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
      default:
         errorMessage("parent variable is neither ORIGINAL nor AGGREGATED");
         return SCIP_INVALIDDATA;
      }

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
      freeBlockMemoryArray(memhdr, &(*var)->data.multaggr.vars, (*var)->data.multaggr.varssize);
      freeBlockMemoryArray(memhdr, &(*var)->data.multaggr.scalars, (*var)->data.multaggr.varssize);
      break;
   default:
      errorMessage("Unknown variable status");
      abort();
   }

   /* release all parent variables and free the parentvars array */
   CHECK_OKAY( varFreeParents(var, memhdr, set, lp) );

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
   char name[MAXSTRLEN];
   VARTYPE vartype;

   assert(origvar != NULL);
   assert(origvar->varstatus == SCIP_VARSTATUS_ORIGINAL);
   assert(origvar->glbdom.lb == origvar->actdom.lb);
   assert(origvar->glbdom.ub == origvar->actdom.ub);
   assert(origvar->data.transvar == NULL);
   assert(transvar != NULL);

   /* convert 0/1 integer variables into binary variables */
   vartype = (VARTYPE)(origvar->vartype);
   if( vartype == SCIP_VARTYPE_INTEGER
      && SCIPsetIsEQ(set, origvar->glbdom.lb, 0.0) && SCIPsetIsEQ(set, origvar->glbdom.ub, 1.0) )
      vartype = SCIP_VARTYPE_BINARY;

   /* create transformed variable */
   sprintf(name, "t_%s", origvar->name);
   CHECK_OKAY( SCIPvarCreateTransformed(transvar, memhdr, set, stat,
                  name, origvar->glbdom.lb, origvar->glbdom.ub, objsense * origvar->obj, vartype, origvar->removeable) );

   /* duplicate hole lists */
   CHECK_OKAY( holelistDuplicate(&(*transvar)->glbdom.holelist, memhdr, set, origvar->glbdom.holelist) );
   CHECK_OKAY( holelistDuplicate(&(*transvar)->actdom.holelist, memhdr, set, origvar->actdom.holelist) );

   /* link original and transformed variable */
   origvar->data.transvar = *transvar;
   CHECK_OKAY( varAddParent(*transvar, memhdr, set, origvar) );

   /* copy rounding locks */
   (*transvar)->nlocksdown = origvar->nlocksdown;
   (*transvar)->nlocksup = origvar->nlocksup;

   debugMessage("transformed variable: <%s>[%p] -> <%s>[%p]\n", origvar->name, origvar, (*transvar)->name, *transvar);

   return SCIP_OKAY;
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

   switch( var->varstatus )
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

   default:
      errorMessage("unknown variable status");
      abort();
   }

   assert(var->nlocksdown >= 0);
   assert(var->nlocksup >= 0);
}

/** increases lock number for rounding down; tells variable, that rounding its value down will make the
 *  solution infeasible
 */
void SCIPvarForbidRoundDown(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("forbid rounding down of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 1, 0);
}

/** increases lock number for rounding up; tells variable, that rounding its value up will make the solution infeasible */
void SCIPvarForbidRoundUp(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("forbid rounding up of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 0, 1);
}

/**< increases lock number for rounding down and up; tells variable, that rounding value in either direction will
 *   make the solution infeasible
 */
void SCIPvarForbidRound(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("forbid rounding of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 1, 1);
}

/** decreases lock number for rounding down; cancels a prior forbidRoundDown() */
void SCIPvarAllowRoundDown(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("allow rounding down of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, -1, 0);
}

/** decreases lock number for rounding up; cancels a prior forbidRoundUp() */
void SCIPvarAllowRoundUp(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("allow rounding up of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, 0, -1);
}

/** decreases lock number for rounding down & up; cancels a prior forbidRound() */
void SCIPvarAllowRound(
   VAR*             var                 /**< problem variable */
   )
{
   debugMessage("allow rounding of <%s> (locks=%d/%d)\n", var->name, var->nlocksdown, var->nlocksup);

   varAddRoundLocks(var, -1, -1);
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

   debugMessage("get down locks of <%s>\n", var->name);

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

   debugMessage("get up locks of <%s>\n", var->name);

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

/** converts loose transformed variable into column variable, creates LP column */
RETCODE SCIPvarColumn(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat                /**< problem statistics */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_LOOSE);

   debugMessage("creating column for variable <%s>\n", var->name);

   var->varstatus = SCIP_VARSTATUS_COLUMN;

   CHECK_OKAY( SCIPcolCreate(&var->data.col, memhdr, set, stat, var, 0, NULL, NULL, var->removeable) );

   return SCIP_OKAY;
}

/** converts variable into fixed variable */
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             fixedval            /**< value to fix variable at */
   )
{
   assert(var != NULL);
   assert(var->glbdom.lb == var->actdom.lb);
   assert(var->glbdom.ub == var->actdom.ub);
   assert(SCIPsetIsLE(set, var->glbdom.lb, fixedval));
   assert(SCIPsetIsLE(set, fixedval, var->glbdom.ub));
   assert(var->vartype == SCIP_VARTYPE_CONTINOUS || SCIPsetIsIntegral(set, fixedval));

   debugMessage("fix variable <%s> to %g\n", var->name, fixedval);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("Cannot fix an untransformed original variable");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarFix(var->data.transvar, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue, fixedval) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* move the objective value to the problem's objective offset */
      SCIPprobIncObjoffset(prob, set, var->obj * fixedval);
      CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, tree, lp, branchcand, eventqueue, 0.0) );

      /* change variable's bounds to fixed value */
      holelistFree(&var->glbdom.holelist, memhdr);
      CHECK_OKAY( SCIPvarChgLbGlobal(var, set, fixedval) );
      CHECK_OKAY( SCIPvarChgUbGlobal(var, set, fixedval) );
      holelistFree(&var->actdom.holelist, memhdr);
      CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, fixedval) );
      CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, fixedval) );

      /* convert variable into fixed variable */
      var->varstatus = SCIP_VARSTATUS_FIXED;
      
      /* update the problem's vars array */
      CHECK_OKAY( SCIPprobVarFixed(prob, set, branchcand, var) );
      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("Cannot fix a column variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot fix a fixed variable again");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      /* fix aggregation variable y in x = a*y + c, instead of fixing x directly */
      assert(SCIPsetIsZero(set, var->obj));
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      CHECK_OKAY( SCIPvarFix(var->data.aggregate.var, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue,
                     (fixedval - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot fix a multiple aggregated variable");
      return SCIP_INVALIDDATA;

   default:
      errorMessage("Unknown variable status");
      abort();
   }
   
   return SCIP_OKAY;
}

/** converts variable into aggregated variable */
RETCODE SCIPvarAggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             aggvar,             /**< variable $y$ in aggregation $x = a*y + c$ */
   Real             scalar,             /**< multiplier $a$ in aggregation $x = a*y + c$ */
   Real             constant,           /**< constant shift $c$ in aggregation $x = a*y + c$ */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   Real varlb;
   Real varub;
   Real aggvarlb;
   Real aggvarub;
   Real obj;
   int nlocksdown;
   int nlocksup;

   assert(var != NULL);
   assert(var->glbdom.lb == var->actdom.lb);
   assert(var->glbdom.ub == var->actdom.ub);
   assert(aggvar->glbdom.lb == aggvar->actdom.lb);
   assert(aggvar->glbdom.ub == aggvar->actdom.ub);
   assert(!SCIPsetIsZero(set, scalar));
   assert(infeasible != NULL);

   debugMessage("aggregate variable <%s> == %g*<%s> %+g\n", var->name, scalar, aggvar->name, constant);

   *infeasible = FALSE;

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("Cannot aggregate an untransformed original variable");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarAggregate(var->data.transvar, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue,
                     aggvar, scalar, constant, infeasible) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* set the aggregated variable's objective value to 0.0 */
      obj = var->obj;
      CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, tree, lp, branchcand, eventqueue, 0.0) );

      /* unlock all rounding locks */
      nlocksdown = var->nlocksdown;
      nlocksup = var->nlocksup;
      var->nlocksdown = 0;
      var->nlocksup = 0;

      /* convert variable into aggregated variable */
      var->varstatus = SCIP_VARSTATUS_AGGREGATED;
      var->data.aggregate.var = aggvar;
      var->data.aggregate.scalar = scalar;
      var->data.aggregate.constant = constant;

      /* make aggregated variable a parent of the aggregation variable */
      CHECK_OKAY( varAddParent(aggvar, memhdr, set, var) );

      /* update the problem's vars array */
      CHECK_OKAY( SCIPprobVarFixed(prob, set, branchcand, var) );

      /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
       * variable and the problem's objective offset
       */
      CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, prob, tree, lp, branchcand, eventqueue, obj) );

      /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variable */
      varAddRoundLocks(var, nlocksdown, nlocksup);

      /* update the bounds of the aggregation variable y in x = a*y + c  ->  y = x/a - c/a */
      if( scalar > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarlb = -set->infinity;
         else
            aggvarlb = var->glbdom.lb / scalar - constant / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarub = set->infinity;
         else
            aggvarub = var->glbdom.ub / scalar - constant / scalar;
      }
      else
      {
         if( SCIPsetIsInfinity(set, -var->glbdom.lb) )
            aggvarub = set->infinity;
         else
            aggvarub = var->glbdom.lb / scalar - constant / scalar;
         if( SCIPsetIsInfinity(set, var->glbdom.ub) )
            aggvarlb = -set->infinity;
         else
            aggvarlb = var->glbdom.ub / scalar - constant / scalar;
      }
      SCIPvarAdjustLb(aggvar, set, &aggvarlb);
      SCIPvarAdjustUb(aggvar, set, &aggvarub);
      aggvarlb = MAX(aggvarlb, aggvar->glbdom.lb);
      aggvarub = MIN(aggvarub, aggvar->glbdom.ub);

      /* check the new bounds */
      if( SCIPsetIsGT(set, aggvarlb, aggvarub) )
      {
         /* the aggregation is infeasible */
         *infeasible = TRUE;
      }
      else if( SCIPsetIsEQ(set, aggvarlb, aggvarub) )
      {
         /* the aggregation variable is fixed */
         CHECK_OKAY( SCIPvarFix(aggvar, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue, aggvarlb) );
      }
      else
      {
         if( SCIPsetIsGT(set, aggvarlb, aggvar->glbdom.lb) )
         {
            CHECK_OKAY( SCIPvarChgLbLocal(aggvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, aggvarlb) );
            CHECK_OKAY( SCIPvarChgLbGlobal(aggvar, set, aggvarlb) );
         }
         if( SCIPsetIsLT(set, aggvarub, aggvar->glbdom.ub) )
         {
            CHECK_OKAY( SCIPvarChgUbLocal(aggvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, aggvarub) );
            CHECK_OKAY( SCIPvarChgUbGlobal(aggvar, set, aggvarub) );
         }

         /* update the hole list of the aggregation variable */
         todoMessage("update hole list of aggregation variable");
      }

      /* update flags of aggregation variable */
      aggvar->removeable &= var->removeable;
      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("Cannot aggregate a column variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot aggregate a fixed variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      errorMessage("Cannot aggregate an aggregated variable again");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot aggregate a multiple aggregated variable");
      return SCIP_INVALIDDATA;

   default:
      errorMessage("Unknown variable status");
      abort();
   }

   return SCIP_OKAY;
}

/** converts variable into multi-aggregated variable */
RETCODE SCIPvarMultiaggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              naggvars,           /**< number $n$ of variables in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   VAR**            aggvars,            /**< variables $y_i$ in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   Real*            scalars,            /**< multipliers $a_i$ in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   Real             constant,           /**< constant shift $c$ in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   Real obj;
   int nlocksdown;
   int nlocksup;
   int v;

   assert(var != NULL);
   assert(var->glbdom.lb == var->actdom.lb);
   assert(var->glbdom.ub == var->actdom.ub);
   assert(naggvars == 0 || aggvars != NULL);
   assert(naggvars == 0 || scalars != NULL);
   assert(infeasible != NULL);

   debugMessage("multi-aggregate variable <%s> == ...%d vars... %+g\n", var->name, naggvars, constant);

   *infeasible = FALSE;

   /* check, if we are in one of the simple cases */
   if( naggvars == 0 )
      return SCIPvarFix(var, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue, constant);
   else if( naggvars == 1)
      return SCIPvarAggregate(var, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue, 
         aggvars[0], scalars[0], constant, infeasible);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("Cannot multi-aggregate an untransformed original variable");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarMultiaggregate(var->data.transvar, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue,
                     naggvars, aggvars, scalars, constant, infeasible) );
      break;

   case SCIP_VARSTATUS_LOOSE:
      assert(!SCIPeventqueueIsDelayed(eventqueue)); /* otherwise, the pseudo objective value update gets confused */

      /* set the aggregated variable's objective value to 0.0 */
      obj = var->obj;
      CHECK_OKAY( SCIPvarChgObj(var, memhdr, set, tree, lp, branchcand, eventqueue, 0.0) );

      /* unlock all rounding locks */
      nlocksdown = var->nlocksdown;
      nlocksup = var->nlocksup;
      var->nlocksdown = 0;
      var->nlocksup = 0;

      /* convert variable into multi-aggregated variable */
      var->varstatus = SCIP_VARSTATUS_MULTAGGR;
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &var->data.multaggr.vars, aggvars, naggvars) );
      ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &var->data.multaggr.scalars, scalars, naggvars) );
      var->data.multaggr.constant = constant;
      var->data.multaggr.nvars = naggvars;
      var->data.multaggr.varssize = naggvars;

      /* update the problem's vars array */
      CHECK_OKAY( SCIPprobVarFixed(prob, set, branchcand, var) );

      /* reset the objective value of the aggregated variable, thus adjusting the objective value of the aggregation
       * variables and the problem's objective offset
       */
      CHECK_OKAY( SCIPvarAddObj(var, memhdr, set, prob, tree, lp, branchcand, eventqueue, obj) );

      /* relock the rounding locks of the variable, thus increasing the locks of the aggregation variables */
      varAddRoundLocks(var, nlocksdown, nlocksup);

      /* update flags of aggregation variables */
      for( v = 0; v < naggvars; ++v )
         aggvars[v]->removeable &= var->removeable;
      break;

   case SCIP_VARSTATUS_COLUMN:
      errorMessage("Cannot multi-aggregate a column variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot multi-aggregate a fixed variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED:
      errorMessage("Cannot multi-aggregate an aggregated variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot multi-aggregate a multiple aggregated variable again");
      return SCIP_INVALIDDATA;

   default:
      errorMessage("Unknown variable status");
      abort();
   }
   
   return SCIP_OKAY;
}

/** changes type of variable; cannot be called, if var belongs to a problem */
RETCODE SCIPvarChgType(
   VAR*             var,                /**< problem variable to change */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(var != NULL);

   debugMessage("change type of <%s> from %d to %d\n", var->name, var->vartype, vartype);

   if( var->probindex >= 0 )
   {
      errorMessage("cannot change type of variable already in the problem");
      return SCIP_INVALIDDATA;
   }
   
   var->vartype = vartype;

   return SCIP_OKAY;
}

/** appends OBJCHANGED event to the event queue */
static
RETCODE varEventObjChanged(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             oldobj,             /**< old objective value for variable */
   Real             newobj              /**< new objective value for variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);
   assert(var->eventfilter != NULL);
   assert(!SCIPsetIsEQ(set, oldobj, newobj));

   /* check, if the variable is being tracked for objective changes
    * COLUMN and LOOSE variables are tracked always, because the pseudo objective value has to be updated
    */
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_OBJCHANGED) != 0)
      || var->varstatus == SCIP_VARSTATUS_COLUMN
      || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;

      CHECK_OKAY( SCIPeventCreateObjChanged(&event, memhdr, var, oldobj, newobj) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, tree, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/** changes objective value of variable */
RETCODE SCIPvarChgObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newobj              /**< new objective value for variable */
   )
{
   Real oldobj;

   assert(var != NULL);

   debugMessage("changing objective value of <%s> from %g to %g\n", var->name, var->obj, newobj);

   if( !SCIPsetIsEQ(set, var->obj, newobj) )
   {
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarChgObj(var->data.transvar, memhdr, set, tree, lp, branchcand, eventqueue, newobj) );
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
         CHECK_OKAY( varEventObjChanged(var, memhdr, set, tree, lp, branchcand, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("cannot change objective value of a fixed, aggregated, or multi aggregated variable");
         return SCIP_INVALIDDATA;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** adds value to objective value of variable */
RETCODE SCIPvarAddObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             addobj              /**< additional objective value for variable */
   )
{
   Real oldobj;
   int i;

   assert(var != NULL);

   debugMessage("adding %g to objective value %g of <%s>\n", addobj, var->obj, var->name);

   if( !SCIPsetIsZero(set, addobj) )
   {
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.transvar != NULL )
         {
            CHECK_OKAY( SCIPvarAddObj(var->data.transvar, memhdr, set, prob, tree, lp, branchcand, eventqueue, addobj) );
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
         CHECK_OKAY( varEventObjChanged(var, memhdr, set, tree, lp, branchcand, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
         assert(SCIPsetIsEQ(set, var->actdom.lb, var->actdom.ub));
         SCIPprobIncObjoffset(prob, set, var->actdom.lb * addobj);
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         /* x = a*y + c  ->  add a*addobj to obj. val. of y, and c*addobj to obj. offset of problem */
         SCIPprobIncObjoffset(prob, set, var->data.aggregate.constant * addobj);
         CHECK_OKAY( SCIPvarAddObj(var->data.aggregate.var, memhdr, set, prob, tree, lp, branchcand, eventqueue,
                        var->data.aggregate.scalar * addobj) );
         break;

      case SCIP_VARSTATUS_MULTAGGR:
         /* x = a_1*y_1 + ... + a_n*y_n  + c  ->  add a_i*addobj to obj. val. of y_i, and c*addobj to obj. offset */
         SCIPprobIncObjoffset(prob, set, var->data.multaggr.constant * addobj);
         for( i = 0; i < var->data.multaggr.nvars; ++i )
         {
            CHECK_OKAY( SCIPvarAddObj(var->data.multaggr.vars[i], memhdr, set, prob, tree, lp, branchcand, eventqueue,
                           var->data.multaggr.scalars[i] * addobj) );
         }
         break;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/* forward declaration, because both methods call each other recursively */

/** performs the actual change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   );

/** performs the actual change in lower bound, changes all parents accordingly */
static            
RETCODE varProcessChgLbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);

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

      switch( parentvar->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            assert(SCIPsetIsEQ(set, parentvar->glbdom.lb,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else
         {
            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert(SCIPsetIsEQ(set, parentvar->glbdom.ub,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
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
RETCODE varProcessChgUbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   VAR* parentvar;
   Real oldbound;
   int i;

   assert(var != NULL);
   assert(var->varstatus != SCIP_VARSTATUS_ORIGINAL);

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

      switch( parentvar->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;
         
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of y */
            assert(SCIPsetIsEQ(set, parentvar->glbdom.ub,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else 
         {
            /* a < 0 -> change lower bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert(SCIPsetIsEQ(set, parentvar->glbdom.lb,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbGlobal(parentvar, set,
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         break;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes global lower bound of variable */
RETCODE SCIPvarChgLbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing global lower bound of <%s> from %g to %g\n", var->name, var->glbdom.lb, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->glbdom.lb, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgLbGlobal(var->data.transvar, set, newbound) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         var->glbdom.lb = newbound;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      CHECK_OKAY( varProcessChgLbGlobal(var, set, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         assert(SCIPsetIsEQ(set, var->glbdom.lb,
                   var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbGlobal(var->data.aggregate.var, set,
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert(SCIPsetIsEQ(set, var->glbdom.lb,
                   var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbGlobal(var->data.aggregate.var, set,
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

   return SCIP_OKAY;
}

/** changes global upper bound of variable */
RETCODE SCIPvarChgUbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing global upper bound of <%s> from %g to %g\n", var->name, var->glbdom.ub, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->glbdom.ub, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgUbGlobal(var->data.transvar, set, newbound) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         var->glbdom.ub = newbound;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      CHECK_OKAY( varProcessChgUbGlobal(var, set, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         assert(SCIPsetIsEQ(set, var->glbdom.ub,
                   var->data.aggregate.var->glbdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbGlobal(var->data.aggregate.var, set,
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert(SCIPsetIsEQ(set, var->glbdom.ub,
                   var->data.aggregate.var->glbdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbGlobal(var->data.aggregate.var, set,
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

   return SCIP_OKAY;
}

/** changes global bound of variable */
RETCODE SCIPvarChgBdGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
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
      errorMessage("Unknown bound type");
      return SCIP_INVALIDDATA;
   }
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
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_LBCHANGED) != 0)
      || var->varstatus == SCIP_VARSTATUS_COLUMN
      || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;

      CHECK_OKAY( SCIPeventCreateLbChanged(&event, memhdr, var, oldbound, newbound) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, tree, lp, branchcand, NULL, &event) );
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
   if( (var->eventfilter->len > 0 && (var->eventfilter->eventmask & SCIP_EVENTTYPE_UBCHANGED) != 0)
      || var->varstatus == SCIP_VARSTATUS_COLUMN
      || var->varstatus == SCIP_VARSTATUS_LOOSE )
   {
      EVENT* event;
   
      CHECK_OKAY( SCIPeventCreateUbChanged(&event, memhdr, var, oldbound, newbound) );
      CHECK_OKAY( SCIPeventqueueAdd(eventqueue, memhdr, set, tree, lp, branchcand, NULL, &event) );
   }

   return SCIP_OKAY;
}

/* forward declaration, because both methods call each other recursively */

/** performs the actual change in upper bound, changes all parents accordingly */
static            
RETCODE varProcessChgUbLocal(
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
RETCODE varProcessChgLbLocal(
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

   debugMessage("process changing lower bound of <%s> from %f to %f\n", var->name, var->actdom.lb, newbound);

   if( SCIPsetIsEQ(set, newbound, var->actdom.lb) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->actdom.lb;
   var->actdom.lb = newbound;

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
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change lower bound of y */
            assert(SCIPsetIsEQ(set, parentvar->actdom.lb,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else 
         {
            /* a < 0 -> change upper bound of y */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert(SCIPsetIsEQ(set, parentvar->actdom.ub,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
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
RETCODE varProcessChgUbLocal(
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

   debugMessage("process changing upper bound of <%s> from %f to %f\n", var->name, var->actdom.ub, newbound);

   if( SCIPsetIsEQ(set, newbound, var->actdom.ub) )
      return SCIP_OKAY;

   /* change the bound */
   oldbound = var->actdom.ub;
   var->actdom.ub = newbound;

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
         errorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable");
         abort();
      
      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c */
         assert(parentvar->data.aggregate.var == var);
         if( SCIPsetIsPositive(set, parentvar->data.aggregate.scalar) )
         {
            /* a > 0 -> change upper bound of x */
            assert(SCIPsetIsEQ(set, parentvar->actdom.ub,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgUbLocal(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         else
         {
            /* a < 0 -> change lower bound of x */
            assert(SCIPsetIsNegative(set, parentvar->data.aggregate.scalar));
            assert(SCIPsetIsEQ(set, parentvar->actdom.lb,
                      oldbound * parentvar->data.aggregate.scalar + parentvar->data.aggregate.constant));
            CHECK_OKAY( varProcessChgLbLocal(parentvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                           parentvar->data.aggregate.scalar * newbound + parentvar->data.aggregate.constant) );
         }
         break;

      default:
         errorMessage("unknown variable status");
         abort();
      }
   }

   return SCIP_OKAY;
}

/** changes current local lower bound of variable */
RETCODE SCIPvarChgLbLocal(
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
   assert(var != NULL);

   debugMessage("changing lower bound of <%s> from %g to %g\n", var->name, var->actdom.lb, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->actdom.lb, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.transvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        newbound) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         stat->nboundchanges++;
         var->actdom.lb = newbound;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      stat->nboundchanges++;
      CHECK_OKAY( varProcessChgLbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change lower bound of y */
         assert(SCIPsetIsEQ(set, var->actdom.lb,
                   var->data.aggregate.var->actdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.aggregate.var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change upper bound of y */
         assert(SCIPsetIsEQ(set, var->actdom.lb,
                   var->data.aggregate.var->actdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.aggregate.var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
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
   
   return SCIP_OKAY;
}

/** changes current local upper bound of variable */
RETCODE SCIPvarChgUbLocal(
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
   assert(var != NULL);

   debugMessage("changing upper bound of <%s> from %g to %g\n", var->name, var->actdom.ub, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   if( SCIPsetIsEQ(set, var->actdom.ub, newbound) )
      return SCIP_OKAY;

   /* change bounds of attached variables */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar != NULL )
      {
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.transvar, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        newbound) );
      }
      else
      {
         assert(SCIPstage(set->scip) == SCIP_STAGE_PROBLEM);
         stat->nboundchanges++;
         var->actdom.ub = newbound;
      }
      break;
         
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      stat->nboundchanges++;
      CHECK_OKAY( varProcessChgUbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, newbound) );
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable");
      return SCIP_INVALIDDATA;
         
   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( SCIPsetIsPositive(set, var->data.aggregate.scalar) )
      {
         /* a > 0 -> change upper bound of y */
         assert(SCIPsetIsEQ(set, var->actdom.ub,
                   var->data.aggregate.var->actdom.ub * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgUbLocal(var->data.aggregate.var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
                        (newbound - var->data.aggregate.constant)/var->data.aggregate.scalar) );
      }
      else if( SCIPsetIsNegative(set, var->data.aggregate.scalar) )
      {
         /* a < 0 -> change lower bound of y */
         assert(SCIPsetIsEQ(set, var->actdom.ub,
                   var->data.aggregate.var->actdom.lb * var->data.aggregate.scalar + var->data.aggregate.constant));
         CHECK_OKAY( SCIPvarChgLbLocal(var->data.aggregate.var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 
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

   return SCIP_OKAY;
}

/** changes current local bound of variable */
RETCODE SCIPvarChgBdLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
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
      return SCIPvarChgLbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, newbound);
   case SCIP_BOUNDTYPE_UPPER:
      return SCIPvarChgUbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, newbound);
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

   debugMessage("adjust lower bound %g of <%s>\n", *lb, var->name);

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

   debugMessage("adjust upper bound %g of <%s>\n", *ub, var->name);

   if( var->vartype != SCIP_VARTYPE_CONTINOUS )
   {
      /* adjust new bound to integral value */
      *ub = SCIPsetFloor(set, *ub);
   }
}

/** changes lower bound of variable in current dive */
RETCODE SCIPvarChgLbDive(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing lower bound of <%s> to %g in current dive\n", var->name, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   /* change bounds of attached variables */
   switch( var->varstatus )
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
      errorMessage("cannot change variable's bounds in dive for LOOSE variables");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable");
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

   return SCIP_OKAY;
}

/** changes upper bound of variable in current dive */
RETCODE SCIPvarChgUbDive(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   )
{
   assert(var != NULL);

   debugMessage("changing upper bound of <%s> to %g in current dive\n", var->name, newbound);

   if( SCIPsetIsZero(set, newbound) )
      newbound = 0.0;

   /* change bounds of attached variables */
   switch( var->varstatus )
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
      errorMessage("cannot change variable's bounds in dive for LOOSE variables");
      return SCIP_INVALIDDATA;
      
   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot change the bounds of a fixed variable");
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

   return SCIP_OKAY;
}

/** adds a hole to the variable's global domain and to its current local domain */
RETCODE SCIPvarAddHoleGlobal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(var != NULL);

   debugMessage("adding global hole (%g,%g) to <%s>\n", left, right, var->name);

   CHECK_OKAY( domAddHole(&var->glbdom, memhdr, set, left, right) );
   CHECK_OKAY( domAddHole(&var->actdom, memhdr, set, left, right) );

   return SCIP_OKAY;
}

/** adds a hole to the variable's current local domain */
RETCODE SCIPvarAddHoleLocal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(var != NULL);

   debugMessage("adding local hole (%g,%g) to <%s>\n", left, right, var->name);

   CHECK_OKAY( domAddHole(&var->actdom, memhdr, set, left, right) );

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

   return (VARSTATUS)(var->varstatus);
}

/** gets type of variable */
VARTYPE SCIPvarGetType(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return (VARTYPE)(var->vartype);
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

/** gets corresponding transformed variable of an original variable */
VAR* SCIPvarGetTransformed(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_ORIGINAL);

   return var->data.transvar;
}

/** gets corresponding active problem variable of a variable */
VAR* SCIPvarGetProbvar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   debugMessage("get problem variable of <%s>\n", var->name);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         errorMessage("original variable has no transformed variable attached");
         return NULL;
      }
      return SCIPvarGetProbvar(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return var;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("fixed variable has no corresponding active problem variable");
      return NULL;

   case SCIP_VARSTATUS_AGGREGATED:
      return SCIPvarGetProbvar(var->data.aggregate.var);

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("multiple aggregated variable has no single corresponding active problem variable");
      return NULL;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** transforms given variable, boundtype and bound to the corresponding active variable values */
RETCODE SCIPvarTransformBound(
   VAR**            var,                /**< pointer to problem variable */
   Real*            bound,              /**< pointer to bound value to transform */
   BOUNDTYPE*       boundtype           /**< pointer to type of bound: lower or upper bound */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert(bound != NULL);
   assert(boundtype != NULL);

   debugMessage("transform bound %g of type %d of variable <%s>\n", *bound, *boundtype, (*var)->name);

   switch( (*var)->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( (*var)->data.transvar == NULL )
      {
         errorMessage("original variable has no transformed variable attached");
         return SCIP_INVALIDDATA;
      }
      *var = (*var)->data.transvar;
      CHECK_OKAY( SCIPvarTransformBound(var, bound, boundtype) );
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      break;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("fixed variable has no corresponding active problem variable");
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
      CHECK_OKAY( SCIPvarTransformBound(var, bound, boundtype) );
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("multiple aggregated variable has no single corresponding active problem variable");
      return SCIP_INVALIDDATA;

   default:
      errorMessage("unknown variable status");
      abort();
   }

   return SCIP_OKAY;
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

/** gets aggregation variable $y$ of an aggregated variable $x = a*y + c$ */
VAR* SCIPvarGetAggrVar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.var;
}

/** gets aggregation scalar $a$ of an aggregated variable $x = a*y + c$ */
Real SCIPvarGetAggrScalar(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.scalar;
}

/** gets aggregation constant $c$ of an aggregated variable $x = a*y + c$ */
Real SCIPvarGetAggrConstant(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_AGGREGATED);

   return var->data.aggregate.constant;
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

   return var->actdom.lb;
}
   
/** gets current upper bound of variable */
Real SCIPvarGetUbLocal(
   VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   return var->actdom.ub;
}

/** gets lower bound of variable in current dive */
Real SCIPvarGetLbDive(
   VAR*             var,                /**< problem variable */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetLbDive(var->data.transvar, set);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetLb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->actdom.lb;
      
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
         errorMessage("scalar is zero in aggregation");
         abort();
      }
      
   case SCIP_VARSTATUS_MULTAGGR:
      todoMessage("get the sides of the corresponding linear constraint");
      errorMessage("getting the bounds of a multiple aggregated variable is not implemented yet");
      abort();
      
   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** gets upper bound of variable in current dive */
Real SCIPvarGetUbDive(
   VAR*             var,                /**< problem variable */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(var != NULL);

   /* get bounds of attached variables */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.transvar != NULL);
      return SCIPvarGetUbDive(var->data.transvar, set);
         
   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolGetUb(var->data.col);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      return var->actdom.ub;
      
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
         errorMessage("scalar is zero in aggregation");
         abort();
      }
      
   case SCIP_VARSTATUS_MULTAGGR:
      todoMessage("get the sides of the corresponding linear constraint");
      errorMessage("getting the bounds of a multiple aggregated variable is not implemented yet");
      abort();
      
   default:
      errorMessage("unknown variable status");
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
      return var->actdom.lb;
   else
      return var->actdom.ub;
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

/** gets primal LP solution value of variable */
Real SCIPvarGetLPSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real primsol;
   int i;

   assert(var != NULL);

   debugMessage("get LP solution of <%s>\n", var->name);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetLPSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
      assert(var->actdom.lb <= 0.0 && var->actdom.ub >= 0.0);
      return 0.0;

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return var->data.col->primsol;

   case SCIP_VARSTATUS_FIXED:
      assert(var->actdom.lb == var->actdom.ub);
      return var->actdom.lb;

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
      abort();
   }
}

/** gets pseudo solution value of variable at actual node */
Real SCIPvarGetPseudoSol(
   VAR*             var                 /**< problem variable */
   )
{
   Real pseudosol;
   int i;

   assert(var != NULL);

   debugMessage("get pseudo solution of <%s>\n", var->name);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
         return SCIP_INVALID;
      return SCIPvarGetPseudoSol(var->data.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPvarGetBestBound(var);

   case SCIP_VARSTATUS_FIXED:
      assert(var->actdom.lb == var->actdom.ub);
      return var->actdom.lb;

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
      abort();
   }
}

/** gets solution value of variable at actual node: if LP was solved at the node, the method returns the
 *  LP primal solution value, otherwise the pseudo solution
 */
Real SCIPvarGetSol(
   VAR*             var,                /**< problem variable */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(tree != NULL);

   debugMessage("get solution of <%s>\n", var->name);

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
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   )
{
   int i;

   assert(var != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsInfinity(set, ABS(val)));

   debugMessage("adding coefficient %g<%s> to row <%s>\n", val, var->name, row->name);

   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.transvar == NULL )
      {
         char s[MAXSTRLEN];
         sprintf(s, "Cannot add untransformed original variable <%s> to LP row <%s>", var->name, row->name);
         errorMessage(s);
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPvarAddToRow(var->data.transvar, memhdr, set, stat, lp, row, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
      /* convert loose variable into column */
      CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat) );
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
      /* fallthrough */

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      assert(var->data.col->var == var);
      CHECK_OKAY( SCIProwIncCoeff(row, memhdr, set, lp, var->data.col, val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(var->glbdom.lb == var->glbdom.ub);
      assert(var->actdom.lb == var->actdom.ub);
      assert(var->actdom.lb == var->glbdom.lb);
      assert(!SCIPsetIsInfinity(set, ABS(var->actdom.lb)));
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, val * var->actdom.lb) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED:
      assert(var->data.aggregate.var != NULL);
      CHECK_OKAY( SCIPvarAddToRow(var->data.aggregate.var, memhdr, set, stat, lp, row, var->data.aggregate.scalar * val) );
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, var->data.aggregate.constant * val) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      assert(var->data.multaggr.vars != NULL);
      assert(var->data.multaggr.scalars != NULL);
      assert(var->data.multaggr.nvars >= 2);
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         CHECK_OKAY( SCIPvarAddToRow(var->data.multaggr.vars[i], memhdr, set, stat, lp, row, 
                        var->data.multaggr.scalars[i] * val) );
      }
      CHECK_OKAY( SCIProwAddConstant(row, set, stat, lp, var->data.multaggr.constant * val) );
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
   assert((eventtype & !SCIP_EVENTTYPE_VARCHANGED) == 0);
   assert((eventtype & SCIP_EVENTTYPE_VARCHANGED) != 0);

   debugMessage("catch event %d of variable <%s>\n", eventtype, var->name);

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

