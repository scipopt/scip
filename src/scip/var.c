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

/** tracks changes of the variable's domains (fixed sized arrays) */
struct DomChg
{
   BOUNDCHG*        boundchg;           /**< array with changes in bounds of variables */
   HOLECHG*         holechg;            /**< array with changes in hole lists */
   int              nboundchg;          /**< number of bound changes */
   int              nholechg;           /**< number of hole list changes */
};

/** tracks changes of the variable's domains (dynamically sized arrays) */
struct DomChgDyn
{
   DOMCHG           domchg;             /**< domain changes */
   int              boundchgsize;       /**< size of bound change array */
   int              holechgsize;        /**< size of hole change array */
};




/*
 * memory growing methods for dynamically allocated arrays
 */

static
RETCODE ensureBoundchgSize(             /**< ensures, that boundchg array can store at least num entries */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(domchgdyn->domchg.nboundchg <= domchgdyn->boundchgsize);
   
   if( num > domchgdyn->boundchgsize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(domchgdyn->domchg.boundchg, newsize) );
      domchgdyn->boundchgsize = newsize;
   }
   assert(num <= domchgdyn->boundchgsize);

   return SCIP_OKAY;
}

static
RETCODE ensureHolechgSize(              /**< ensures, that holechg array can store at least num additional entries */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of additional entries to store */
   )
{
   assert(domchgdyn->domchg.nholechg <= domchgdyn->holechgsize);
   
   if( num > domchgdyn->holechgsize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(domchgdyn->domchg.holechg, newsize) );
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

RETCODE SCIPdomchgdynCreate(            /**< creates a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn           /**< pointer to dynamically sized domain change data structure */
   )
{
   assert(domchgdyn != NULL);

   ALLOC_OKAY( allocMemory(*domchgdyn) );

   (*domchgdyn)->domchg.boundchg = NULL;
   (*domchgdyn)->domchg.holechg = NULL;
   (*domchgdyn)->domchg.nboundchg = 0;
   (*domchgdyn)->domchg.nholechg = 0;
   (*domchgdyn)->boundchgsize = 0;
   (*domchgdyn)->holechgsize = 0;

   return SCIP_OKAY;
};

void SCIPdomchgdynFree(                 /**< frees a dynamically sized domain change data structure */
   DOMCHGDYN**      domchgdyn           /**< pointer to dynamically sized domain change data structure */
   )
{
   assert(domchgdyn != NULL);
   assert(*domchgdyn != NULL);

   freeMemoryArrayNull((*domchgdyn)->domchg.boundchg);
   freeMemoryArrayNull((*domchgdyn)->domchg.holechg);
   freeMemory(*domchgdyn);
}

RETCODE SCIPdomchgdynCopy(              /**< copies data from fixed size domain change into dynamically sized one */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   DOMCHG*          domchg              /**< static domain change */
   )
{
   assert(domchgdyn != NULL);
   assert(set != NULL);

   if( domchg == NULL )
   {
      domchgdyn->domchg.nboundchg = 0;
      domchgdyn->domchg.nholechg = 0;
   }
   else
   {
      CHECK_OKAY( ensureBoundchgSize(domchgdyn, set, domchg->nboundchg) );
      CHECK_OKAY( ensureHolechgSize(domchgdyn, set, domchg->nholechg) );
      
      copyMemoryArray(domchgdyn->domchg.boundchg, domchg->boundchg, domchg->nboundchg);
      copyMemoryArray(domchgdyn->domchg.holechg, domchg->holechg, domchg->nholechg);
      domchgdyn->domchg.nboundchg = domchg->nboundchg;
      domchgdyn->domchg.nholechg = domchg->nholechg;
   }

   return SCIP_OKAY;
}

RETCODE SCIPdomchgdynAddBoundchg(       /**< adds bound change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   assert(domchgdyn != NULL);
   assert(var != NULL);

   CHECK_OKAY( ensureBoundchgSize(domchgdyn, set, domchgdyn->boundchgsize+1) );
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].var = var;
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].newbound = newbound;
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].oldbound = oldbound;
   domchgdyn->domchg.boundchg[domchgdyn->domchg.nboundchg].boundtype = boundtype;

   return SCIP_OKAY;
}

RETCODE SCIPdomchgdynAddHolechg(        /**< adds hole change to domain changes */
   DOMCHGDYN*       domchgdyn,          /**< dynamically sized domain change data structure */
   const SET*       set,                /**< global SCIP settings */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   )
{
   assert(domchgdyn != NULL);
   assert(ptr != NULL);

   CHECK_OKAY( ensureHolechgSize(domchgdyn, set, domchgdyn->holechgsize+1) );
   domchgdyn->domchg.holechg[domchgdyn->domchg.nholechg].ptr = ptr;
   domchgdyn->domchg.holechg[domchgdyn->domchg.nholechg].newlist = newlist;
   domchgdyn->domchg.holechg[domchgdyn->domchg.nholechg].oldlist = oldlist;

   return SCIP_OKAY;
}

RETCODE SCIPdomchgCreate(               /**< creates domain change data (fixed size) from dynamically sized data */
   DOMCHG**         domchg,             /**< pointer to fixed size domain change data */
   MEMHDR*          memhdr,             /**< block memory */
   const DOMCHGDYN* domchgdyn           /**< dynamically sized domain change data structure */
   )
{
   assert(domchg != NULL);
   assert(memhdr != NULL);
   assert(domchgdyn != NULL);

   if( domchgdyn->domchg.nboundchg > 0 || domchgdyn->domchg.nholechg > 0 )
   {
      ALLOC_OKAY( allocBlockMemory(memhdr, *domchg) );
      
      if( domchgdyn->domchg.nboundchg > 0 )
      {
         ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*domchg)->boundchg, domchgdyn->domchg.boundchg,
                        domchgdyn->domchg.nboundchg) );
      }
      else
         (*domchg)->boundchg = NULL;
      
      if( domchgdyn->domchg.nholechg > 0 )
      {
         ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*domchg)->holechg, domchgdyn->domchg.holechg,
                        domchgdyn->domchg.nholechg) );
      }
      else
         (*domchg)->holechg = NULL;
      
      (*domchg)->nboundchg = domchgdyn->domchg.nboundchg;
      (*domchg)->nholechg = domchgdyn->domchg.nholechg;
   }
   else
   {
      assert(domchgdyn->domchg.nboundchg == 0);
      assert(domchgdyn->domchg.nholechg == 0);
      *domchg = NULL;
   }

   return SCIP_OKAY;
}

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
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   VAR* var;
   int i;

   assert(lp != NULL);

   if( domchg == NULL )
      return SCIP_OKAY;

   /* apply bound changes */
   for( i = 0; i < domchg->nboundchg; ++i )
   {
      var = domchg->boundchg[i].var;
      assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

      /* apply bound change to the LP data */
      switch( domchg->boundchg[i].boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:
         CHECK_OKAY( SCIPvarChgLb(var, set, lp, domchg->boundchg[i].newbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUb(var, set, lp, domchg->boundchg[i].newbound) );
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
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   VAR* var;
   int i;

   assert(lp != NULL);

   if( domchg == NULL )
      return SCIP_OKAY;

   /* undo bound changes */
   for( i = domchg->nboundchg-1; i >= 0; --i )
   {
      var = domchg->boundchg[i].var;
      assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

      /* apply bound change to the LP data */
      switch( domchg->boundchg[i].boundtype )
      {
      case SCIP_BOUNDTYPE_LOWER:
         CHECK_OKAY( SCIPvarChgLb(var, set, lp, domchg->boundchg[i].oldbound) );
         break;
      case SCIP_BOUNDTYPE_UPPER:
         CHECK_OKAY( SCIPvarChgUb(var, set, lp, domchg->boundchg[i].oldbound) );
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
 * methods for variables 
 */

RETCODE SCIPvarCreate(                  /**< creates an original problem variable */
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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

   (*var)->data.transvar = NULL;
   (*var)->origvar = NULL;
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*var)->name, name, strlen(name)+1) );
   (*var)->dom.holelist = NULL;
   (*var)->dom.lb = lb;
   (*var)->dom.ub = ub;
   (*var)->obj = obj;
   (*var)->vartype = vartype;
   (*var)->varstatus = SCIP_VARSTATUS_ORIGINAL;

   return SCIP_OKAY;
}

RETCODE SCIPvarCreateTransformed(       /**< creates a loose variable belonging only to the transformed problem */
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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
   (*var)->vartype = vartype;
   (*var)->varstatus = SCIP_VARSTATUS_LOOSE;

   return SCIP_OKAY;   
}

void SCIPvarFree(                       /**< frees a variable */
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
   )
{
   int i;

   assert(memhdr != NULL);
   assert(var != NULL);
   assert(*var != NULL);
   assert((*var)->varstatus == SCIP_VARSTATUS_ORIGINAL || lp != NULL);

   switch( (*var)->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert((*var)->origvar == NULL); /* original variable cannot have link to other original variable */
      assert((*var)->data.transvar == NULL);  /* cannot free variable, if transformed variable is still existing */
      break;
   case SCIP_VARSTATUS_LOOSE:
      break;
   case SCIP_VARSTATUS_COLUMN:
      SCIPcolFree(&(*var)->data.col, memhdr, set, lp);  /* free corresponding LP column */
      break;
   case SCIP_VARSTATUS_FIXED:
      break;
   case SCIP_VARSTATUS_AGGREGATED:
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

RETCODE SCIPvarTransform(               /**< copies original variable into loose transformed variable */
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   VAR**            transvar            /**< pointer to transformed variable */
   )
{
   assert(origvar != NULL);
   assert(origvar->varstatus == SCIP_VARSTATUS_ORIGINAL);
   assert(origvar->data.transvar == NULL);
   assert(transvar != NULL);

   CHECK_OKAY( SCIPvarCreateTransformed(transvar, memhdr, set, 
                  origvar->name, origvar->dom.lb, origvar->dom.ub, origvar->obj, origvar->vartype) );

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

   CHECK_OKAY( SCIPcolCreate(&var->data.col, memhdr, set, lp, stat, var, 0, NULL, NULL) );

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
      abort();
   case SCIP_VARSTATUS_LOOSE:
      break;
   case SCIP_VARSTATUS_COLUMN:
      /* fix column */
      todoMessage("apply fixings to LP (change rows, offset obj)");
      abort();
      break;
   case SCIP_VARSTATUS_FIXED:
      errorMessage("Cannot fix a fixed variable again");
      abort();
   case SCIP_VARSTATUS_AGGREGATED:
      errorMessage("Cannot fix an aggregated variable");
      abort();
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
}

RETCODE SCIPvarChgLb(                   /**< changes lower bound of variable */
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   )
{
   Bool changed;

   assert(var != NULL);

   changed = !SCIPsetIsEQ(set, var->dom.lb, newbound);
   var->dom.lb = newbound;

   /* if variable is a column, inform LP that bounds of column changed */
   if( changed && var->varstatus == SCIP_VARSTATUS_COLUMN )
   {
      CHECK_OKAY( SCIPcolBoundChanged(var->data.col, set, lp, SCIP_BOUNDTYPE_LOWER) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPvarChgUb(                   /**< changes upper bound of variable */
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   )
{
   Bool changed;

   assert(var != NULL);

   changed = !SCIPsetIsEQ(set, var->dom.ub, newbound);
   var->dom.ub = newbound;

   /* if variable is a column, inform LP that bounds of column changed */
   if( changed && var->varstatus == SCIP_VARSTATUS_COLUMN )
   {
      CHECK_OKAY( SCIPcolBoundChanged(var->data.col, set, lp, SCIP_BOUNDTYPE_UPPER) );
   }

   return SCIP_OKAY;
}

