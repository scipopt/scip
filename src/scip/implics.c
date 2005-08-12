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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: implics.c,v 1.5 2005/08/12 12:36:22 bzfpfend Exp $"

/**@file   implics.c
 * @brief  methods for implications, variable bounds, and clique tables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/var.h"
#include "scip/implics.h"

#ifndef NDEBUG
#include "scip/struct_implics.h"
#endif




/*
 * methods for variable bounds
 */

/** creates a variable bounds data structure */
static
RETCODE vboundsCreate(
   VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(vbounds != NULL);

   ALLOC_OKAY( allocBlockMemory(blkmem, vbounds) );
   (*vbounds)->vars = NULL;
   (*vbounds)->coefs = NULL;
   (*vbounds)->constants = NULL;
   (*vbounds)->len = 0;
   (*vbounds)->size = 0;

   return SCIP_OKAY;
}

/** frees a variable bounds data structure */
void SCIPvboundsFree(
   VBOUNDS**        vbounds,            /**< pointer to store variable bounds data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(vbounds != NULL);

   if( *vbounds != NULL )
   {
      freeBlockMemoryArrayNull(blkmem, &(*vbounds)->vars, (*vbounds)->size);
      freeBlockMemoryArrayNull(blkmem, &(*vbounds)->coefs, (*vbounds)->size);
      freeBlockMemoryArrayNull(blkmem, &(*vbounds)->constants, (*vbounds)->size);
      freeBlockMemory(blkmem, vbounds);
   }
}

/** ensures, that variable bounds arrays can store at least num entries */
static
RETCODE vboundsEnsureSize(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(vbounds != NULL);
   
   /* create variable bounds data structure, if not yet existing */
   if( *vbounds == NULL )
   {
      CHECK_OKAY( vboundsCreate(vbounds, blkmem) );
   }
   assert(*vbounds != NULL);
   assert((*vbounds)->len <= (*vbounds)->size);

   if( num > (*vbounds)->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*vbounds)->vars, (*vbounds)->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*vbounds)->coefs, (*vbounds)->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*vbounds)->constants, (*vbounds)->size, newsize) );
      (*vbounds)->size = newsize;
   }
   assert(num <= (*vbounds)->size);

   return SCIP_OKAY;
}

/** binary searches the insertion position of the given variable in the vbounds data structure */
static
RETCODE vboundsSearchPos(
   VBOUNDS*         vbounds,            /**< variable bounds data structure, or NULL */
   VAR*             var,                /**< variable to search in vbounds data structure */
   int*             insertpos,          /**< pointer to store position where to insert new entry */
   Bool*            found               /**< pointer to store whether the same variable was found at the returned pos */
   )
{
   int varidx;
   int left;
   int right;

   assert(insertpos != NULL);
   assert(found != NULL);

   /* check for empty vbounds data */
   if( vbounds == NULL )
   {
      *insertpos = 0;
      *found = FALSE;
      return SCIP_OKAY;
   }
   assert(vbounds->len >= 0);

   /* binary search for the given variable */
   varidx = SCIPvarGetIndex(var);
   left = -1;
   right = vbounds->len;
   while( left < right-1 )
   {
      int middle;
      int idx;

      middle = (left+right)/2;
      assert(0 <= middle && middle < vbounds->len);
      idx = SCIPvarGetIndex(vbounds->vars[middle]);

      if( varidx < idx )
         right = middle;
      else if( varidx > idx )
         left = middle;
      else
      {
         assert(var == vbounds->vars[middle]);
         *insertpos = middle;
         *found = TRUE;
         return SCIP_OKAY;
      }
   }

   *insertpos = right;
   *found = FALSE;

   return SCIP_OKAY;
}

/** adds a variable bound to the variable bounds data structure */
RETCODE SCIPvboundsAdd(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   BOUNDTYPE        vboundtype,         /**< type of variable bound (LOWER or UPPER) */
   VAR*             var,                /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   Real             coef,               /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   Real             constant            /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   )
{
   int insertpos;
   Bool found;

   assert(vbounds != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);

   /* identify insertion position of variable */
   CHECK_OKAY( vboundsSearchPos(*vbounds, var, &insertpos, &found) );
   if( found )
   {
      /* the same variable already exists in the vbounds data structure: use the better vbound */
      assert(*vbounds != NULL);
      assert(0 <= insertpos && insertpos < (*vbounds)->len);
      assert((*vbounds)->vars[insertpos] == var);

      if( vboundtype == SCIP_BOUNDTYPE_UPPER )
      {
         if( constant + MIN(coef, 0.0) < (*vbounds)->constants[insertpos] + MIN((*vbounds)->coefs[insertpos], 0.0) )
         {
            (*vbounds)->coefs[insertpos] = coef;
            (*vbounds)->constants[insertpos] = constant;
         }
      }
      else
      {
         if( constant + MAX(coef, 0.0) > (*vbounds)->constants[insertpos] + MAX((*vbounds)->coefs[insertpos], 0.0) )
         {
            (*vbounds)->coefs[insertpos] = coef;
            (*vbounds)->constants[insertpos] = constant;
         }
      }
   }
   else
   {
      int i;

      /* the given variable does not yet exist in the vbounds */
      CHECK_OKAY( vboundsEnsureSize(vbounds, blkmem, set, *vbounds != NULL ? (*vbounds)->len+1 : 1) );
      assert(*vbounds != NULL);
      assert(0 <= insertpos && insertpos <= (*vbounds)->len);
      assert(0 <= insertpos && insertpos < (*vbounds)->size);

      /* insert variable at the correct position */
      for( i = (*vbounds)->len; i > insertpos; --i )
      {
         (*vbounds)->vars[i] = (*vbounds)->vars[i-1];
         (*vbounds)->coefs[i] = (*vbounds)->coefs[i-1];
         (*vbounds)->constants[i] = (*vbounds)->constants[i-1];
      }
      (*vbounds)->vars[insertpos] = var;
      (*vbounds)->coefs[insertpos] = coef;
      (*vbounds)->constants[insertpos] = constant;
      (*vbounds)->len++;
   }

   return SCIP_OKAY;
}

/** removes from variable x a variable bound x >=/<= b*z + d with binary or integer z */
RETCODE SCIPvboundsDel(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   VAR*             vbdvar              /**< variable z    in x >=/<= b*z + d */
   )
{
   Bool found;
   int pos;
   int i;

   assert(vbounds != NULL);
   assert(*vbounds != NULL);

   /* searches for variable z in variable bounds of x */
   CHECK_OKAY( vboundsSearchPos(*vbounds, vbdvar, &pos, &found) );
   if( !found )
      return SCIP_OKAY;

   assert(0 <= pos && pos < (*vbounds)->len);
   assert((*vbounds)->vars[pos] == vbdvar);

   /* removes z from variable bounds of x */
   for( i = pos; i < (*vbounds)->len - 1; i++ )
   {
      (*vbounds)->vars[i] = (*vbounds)->vars[i+1];
      (*vbounds)->coefs[i] = (*vbounds)->coefs[i+1];
      (*vbounds)->constants[i] = (*vbounds)->constants[i+1];
   }
   (*vbounds)->len--;

#ifndef NDEBUG
   CHECK_OKAY( vboundsSearchPos(*vbounds, vbdvar, &pos, &found) );
   assert(!found);
#endif

   /* free vbounds data structure, if it is empty */
   if( (*vbounds)->len == 0 )
      SCIPvboundsFree(vbounds, blkmem);

   return SCIP_OKAY;
}

/** reduces the number of variable bounds stored in the given variable bounds data structure */
void SCIPvboundsShrink(
   VBOUNDS**        vbounds,            /**< pointer to variable bounds data structure */
   BLKMEM*          blkmem,             /**< block memory */
   int              newnvbds            /**< new number of variable bounds */
   )
{
   assert(vbounds != NULL);
   assert(*vbounds != NULL);
   assert(newnvbds <= (*vbounds)->len);

   if( newnvbds == 0 )
      SCIPvboundsFree(vbounds, blkmem);
   else
      (*vbounds)->len = newnvbds;
}




/*
 * methods for implications
 */

#ifndef NDEBUG
/** comparator function for implication variables in the implication data structure */
static
DECL_SORTPTRCOMP(compVars)
{  /*lint --e{715}*/
   VAR* var1;
   VAR* var2;
   VARTYPE var1type;
   VARTYPE var2type;
   int var1idx;
   int var2idx;

   var1 = (VAR*)elem1;
   var2 = (VAR*)elem2;
   assert(var1 != NULL);
   assert(var2 != NULL);
   var1type = SCIPvarGetType(var1);
   var2type = SCIPvarGetType(var2);
   var1idx = SCIPvarGetIndex(var1);
   var2idx = SCIPvarGetIndex(var2);

   if( var1type == var2type )
   {
      if( var1idx < var2idx )
         return -1;
      else if( var1idx > var2idx )
         return +1;
      else
         return 0;
   }
   else
   {
      if( var1type == SCIP_VARTYPE_BINARY && var2type != SCIP_VARTYPE_BINARY )
         return -1;
      if( var1type != SCIP_VARTYPE_BINARY && var2type == SCIP_VARTYPE_BINARY )
         return +1;
      else if( var1idx < var2idx )
         return -1;
      else if( var1idx > var2idx )
         return +1;
      else
      {
         assert(var1 == var2);
         return 0;
      }
   }
}

/** performs integrity check on implications data structure */
static
void checkImplics(
   IMPLICS*         implics,            /**< implications data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   Bool varfixing;

   if( implics == NULL )
      return;

   varfixing = FALSE;
   do
   {
      VAR** vars;
      BOUNDTYPE* types;
      Real* bounds;
      int nimpls;
      int nbinimpls;
      int i;
      
      vars = implics->vars[varfixing];
      types = implics->types[varfixing];
      bounds = implics->bounds[varfixing];
      nimpls = implics->nimpls[varfixing];
      nbinimpls = implics->nbinimpls[varfixing];

      assert(0 <= nbinimpls && nbinimpls <= nimpls && nimpls <= implics->size[varfixing]);
      assert(nimpls == 0 || vars != NULL);
      assert(nimpls == 0 || types != NULL);
      assert(nimpls == 0 || bounds != NULL);

      for( i = 0; i < nbinimpls; ++i )
      {
         int cmp;

         assert(SCIPvarGetType(implics->vars[varfixing][i]) == SCIP_VARTYPE_BINARY);
         assert((types[i] == SCIP_BOUNDTYPE_LOWER) == (bounds[i] > 0.5));
         assert(SCIPsetIsFeasEQ(set, bounds[i], 0.0) || SCIPsetIsFeasEQ(set, bounds[i], 1.0));

         if( i == 0 )
            continue;

         cmp = compVars(vars[i-1], vars[i]);
         assert(cmp <= 0);
         assert((cmp == 0) == (vars[i-1] == vars[i]));
         assert(cmp < 0 || (types[i-1] == SCIP_BOUNDTYPE_LOWER && types[i] == SCIP_BOUNDTYPE_UPPER));
      }

      for( i = nbinimpls; i < nimpls; ++i )
      {
         int cmp;
         
         assert(SCIPvarGetType(implics->vars[varfixing][i]) != SCIP_VARTYPE_BINARY);

         if( i == 0 )
            continue;

         cmp = compVars(vars[i-1], vars[i]);
         assert(cmp <= 0);
         assert((cmp == 0) == (vars[i-1] == vars[i]));
         assert(cmp < 0 || (types[i-1] == SCIP_BOUNDTYPE_LOWER && types[i] == SCIP_BOUNDTYPE_UPPER));
      }

      varfixing = !varfixing;
   }
   while( varfixing == TRUE );
}
#else
#define checkImplics(implics,set) /**/
#endif

/** creates an implications data structure */
static
RETCODE implicsCreate(
   IMPLICS**        implics,            /**< pointer to store implications data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(implics != NULL);

   ALLOC_OKAY( allocBlockMemory(blkmem, implics) );

   (*implics)->vars[0] = NULL;
   (*implics)->types[0] = NULL;
   (*implics)->bounds[0] = NULL;
   (*implics)->ids[0] = NULL;
   (*implics)->size[0] = 0;
   (*implics)->nimpls[0] = 0;
   (*implics)->nbinimpls[0] = 0;

   (*implics)->vars[1] = NULL;
   (*implics)->types[1] = NULL;
   (*implics)->bounds[1] = NULL;
   (*implics)->ids[1] = NULL;
   (*implics)->size[1] = 0;
   (*implics)->nimpls[1] = 0;
   (*implics)->nbinimpls[1] = 0;

   return SCIP_OKAY;
}

/** frees an implications data structure */
void SCIPimplicsFree(
   IMPLICS**        implics,            /**< pointer of implications data structure to free */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(implics != NULL);

   if( *implics != NULL )
   {
      freeBlockMemoryArrayNull(blkmem, &(*implics)->vars[0], (*implics)->size[0]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->types[0], (*implics)->size[0]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->bounds[0], (*implics)->size[0]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->ids[0], (*implics)->size[0]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->vars[1], (*implics)->size[1]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->types[1], (*implics)->size[1]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->bounds[1], (*implics)->size[1]);
      freeBlockMemoryArrayNull(blkmem, &(*implics)->ids[1], (*implics)->size[1]);
      freeBlockMemory(blkmem, implics);
   }
}

/** ensures, that arrays for x == 0 or x == 1 in implications data structure can store at least num entries */
static
RETCODE implicsEnsureSize(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             varfixing,          /**< FALSE if size of arrays for x == 0 has to be ensured, TRUE for x == 1 */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(implics != NULL);
   
   /* create implications data structure, if not yet existing */
   if( *implics == NULL )
   {
      CHECK_OKAY( implicsCreate(implics, blkmem) );
   }
   assert(*implics != NULL);
   assert((*implics)->nimpls[varfixing] <= (*implics)->size[varfixing]);

   if( num > (*implics)->size[varfixing] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*implics)->vars[varfixing], (*implics)->size[varfixing],
            newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*implics)->types[varfixing], (*implics)->size[varfixing], 
            newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*implics)->bounds[varfixing], (*implics)->size[varfixing],
            newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*implics)->ids[varfixing], (*implics)->size[varfixing],
            newsize) );
      (*implics)->size[varfixing] = newsize;
   }
   assert(num <= (*implics)->size[varfixing]);

   return SCIP_OKAY;
}

/** returns whether variable y is already contained in implications for x == 0 or x == 1 with the given impltype
 *  y can be contained in structure with y >= b (y_lower) and y <= b (y_upper) 
 */
static
Bool implicsSearchVar(
   IMPLICS*         implics,            /**< implications data structure */
   Bool             varfixing,          /**< FALSE if y is searched in implications for x == 0, TRUE for x == 1 */
   VAR*             implvar,            /**< variable y to search for */
   BOUNDTYPE        impltype,           /**< type of implication y <=/>= b to search for */
   int*             poslower,           /**< pointer to store position of y_lower (inf if not found) */
   int*             posupper,           /**< pointer to store position of y_upper (inf if not found) */
   int*             posadd              /**< pointer to store correct position (with respect to impltype) to add y */
   )
{
   int implvaridx;
   int left;
   int right;
   int middle;
   Bool found;

   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);

   /* set left and right pointer */
   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      if( implics->nbinimpls[varfixing] == 0 )
      {
         /* there are no implications with binary variable y */
         *posadd = 0;
         *poslower = INT_MAX;
         *posupper = INT_MAX;
          return FALSE;
      }      
      left = 0;
      right = implics->nbinimpls[varfixing] - 1;
   }
   else
   {
      if( implics->nimpls[varfixing] == implics->nbinimpls[varfixing] )
      {
         /* there are no implications with nonbinary variable y */
         *posadd = implics->nbinimpls[varfixing];
         *poslower = INT_MAX;
         *posupper = INT_MAX;
         return FALSE;
      }
      left = implics->nbinimpls[varfixing];
      right = implics->nimpls[varfixing] - 1;
   }
   assert(left <= right);

   /* search for y */
   implvaridx = SCIPvarGetIndex(implvar);
   do
   {
      int idx;

      middle = (left + right) / 2;
      idx = SCIPvarGetIndex(implics->vars[varfixing][middle]);
      if( implvaridx < idx )
         right = middle - 1;
      else if( implvaridx > idx )
         left = middle + 1;
      else
      {
         assert(implvar == implics->vars[varfixing][middle]);
         break;
      }
   }
   while( left <= right );
   assert(left <= right+1);

   found = FALSE;
   if( left > right )
   {
      /* y was not found */
      assert(right == -1 || compVars((void*)implics->vars[varfixing][right], (void*)implvar) < 0);
      assert(left >= implics->nimpls[varfixing] || implics->vars[varfixing][left] != implvar);
      *poslower = INT_MAX;
      *posupper = INT_MAX;
      *posadd = left;
      found = FALSE;
   }
   else
   {
      /* y was found, but do we have the correct impltype? */
      assert(implvar == implics->vars[varfixing][middle]);

      /* set poslower and posupper */
      if( implics->types[varfixing][middle] == SCIP_BOUNDTYPE_LOWER )
      {
         /* y was found as y_lower (on position middle) */
         *poslower = middle;
         if( middle + 1 < implics->nimpls[varfixing] && implics->vars[varfixing][middle+1] == implvar )
         {  
            assert(implics->types[varfixing][middle+1] == SCIP_BOUNDTYPE_UPPER);
            *posupper = middle + 1;
         }
         else
            *posupper = INT_MAX;
      }
      else
      {
         /* y was found as y_upper (on position middle) */
         *posupper = middle;
         if( middle - 1 >= 0 && implics->vars[varfixing][middle-1] == implvar )
         {  
            assert(implics->types[varfixing][middle-1] == SCIP_BOUNDTYPE_LOWER);
            *poslower = middle - 1;
         }
         else
            *poslower = INT_MAX;
      }

      /* set posadd */
      if( impltype == SCIP_BOUNDTYPE_LOWER )
      {
         if( *poslower < INT_MAX )
         {
            *posadd = *poslower;
            found = TRUE;
         }
         else
         {
            *posadd = *posupper;
            found = FALSE;
         }
      }     
      else
      {
         if( *posupper < INT_MAX )
         {
            *posadd = *posupper;
            found = TRUE;
         }
         else
         {
            *posadd = (*poslower)+1;
            found = FALSE;
         }
      }
      assert(*posadd < INT_MAX);
   }

   return found;
}

/** adds an implication x == 0/1 -> y <= b or y >= b to the implications data structure;
 *  the implication must be non-redundant
 */
RETCODE SCIPimplicsAdd(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Bool             varfixing,          /**< FALSE if implication for x == 0 has to be added, TRUE for x == 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   Bool*            conflict            /**< pointer to store whether implication causes a conflict for variable x */
   )
{
   int poslower;
   int posupper;
   int posadd;
   Bool found;
   int k;

   assert(implics != NULL);
   assert(*implics == NULL || (*implics)->nbinimpls[varfixing] <= (*implics)->nimpls[varfixing]);
   assert(stat != NULL);
   assert(SCIPvarIsActive(implvar));
   assert(SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(implvar) == SCIP_VARSTATUS_LOOSE); 
   assert((impltype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGT(set, implbound, SCIPvarGetLbGlobal(implvar)))
      || (impltype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLT(set, implbound, SCIPvarGetUbGlobal(implvar))));
   assert(conflict != NULL);

   checkImplics(*implics, set);

   *conflict = FALSE;

   /* check if variable is already contained in implications data structure */
   if( *implics != NULL )
   {
      found = implicsSearchVar(*implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
      assert(poslower >= 0);
      assert(posupper >= 0);
      assert(posadd >= 0 && posadd <= (*implics)->nimpls[varfixing]);
   }
   else
   {
      poslower = INT_MAX;
      posupper = INT_MAX;
      posadd = 0;
   }

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      /* check if y >= b is redundant */
      if( poslower < INT_MAX && SCIPsetIsFeasLE(set, implbound, (*implics)->bounds[varfixing][poslower]) )
         return SCIP_OKAY;

      /* check if y >= b causes conflict for x (i.e. y <= a (with a < b) is also valid) */
      if( posupper < INT_MAX && SCIPsetIsFeasGT(set, implbound, (*implics)->bounds[varfixing][posupper]) )
      {      
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      /* check if entry of the same type already exists */
      if( posadd == poslower )
      {
         /* add y >= b by changing old entry on poslower */
         assert((*implics)->vars[varfixing][poslower] == implvar);
         assert(SCIPsetIsFeasGT(set, implbound, (*implics)->bounds[varfixing][poslower]));
         (*implics)->bounds[varfixing][poslower] = implbound;

         return SCIP_OKAY;
      }
      
      /* add y >= b by creating a new entry on posadd */
      assert(poslower == INT_MAX);

      CHECK_OKAY( implicsEnsureSize(implics, blkmem, set, varfixing,
            *implics != NULL ? (*implics)->nimpls[varfixing]+1 : 1) );
      assert(*implics != NULL);
      
      for( k = (*implics)->nimpls[varfixing]; k > posadd; k-- )
      {
         assert(compVars((void*)(*implics)->vars[varfixing][k-1], (void*)implvar) >= 0);
         (*implics)->vars[varfixing][k] = (*implics)->vars[varfixing][k-1];
         (*implics)->types[varfixing][k] = (*implics)->types[varfixing][k-1];
         (*implics)->bounds[varfixing][k] = (*implics)->bounds[varfixing][k-1];
         (*implics)->ids[varfixing][k] = (*implics)->ids[varfixing][k-1];
      }
      assert(posadd == k);
      (*implics)->vars[varfixing][posadd] = implvar;
      (*implics)->types[varfixing][posadd] = impltype;
      (*implics)->bounds[varfixing][posadd] = implbound;
      (*implics)->ids[varfixing][posadd] = stat->nimplications;
      if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
         (*implics)->nbinimpls[varfixing]++;
      (*implics)->nimpls[varfixing]++;
#ifndef NDEBUG
      for( k = posadd-1; k >= 0; k-- )
         assert(compVars((void*)(*implics)->vars[varfixing][k], (void*)implvar) <= 0);
#endif
      stat->nimplications++;
   }
   else
   {
      /* check if y <= b is redundant */
      if( posupper < INT_MAX && SCIPsetIsFeasGE(set, implbound, (*implics)->bounds[varfixing][posupper]) )
         return SCIP_OKAY;

      /* check if y <= b causes conflict for x (i.e. y >= a (with a > b) is also valid) */
      if( poslower < INT_MAX && SCIPsetIsFeasLT(set, implbound, (*implics)->bounds[varfixing][poslower]) )
      {      
         *conflict = TRUE;
         return SCIP_OKAY;
      }

      /* check if entry of the same type already exists */
      if( posadd == posupper )
      {
         /* add y <= b by changing old entry on posupper */
         assert((*implics)->vars[varfixing][posupper] == implvar);
         assert(SCIPsetIsFeasLT(set, implbound,(*implics)->bounds[varfixing][posupper]));
         (*implics)->bounds[varfixing][posupper] = implbound;

         return SCIP_OKAY;
      }
      
      /* add y <= b by creating a new entry on posadd */
      assert(posupper == INT_MAX);

      CHECK_OKAY( implicsEnsureSize(implics, blkmem, set, varfixing,
            *implics != NULL ? (*implics)->nimpls[varfixing]+1 : 1) );
      assert(*implics != NULL);
      
      for( k = (*implics)->nimpls[varfixing]; k > posadd; k-- )
      {
         assert(compVars((void*)(*implics)->vars[varfixing][k-1], (void*)implvar) >= 0);
         (*implics)->vars[varfixing][k] = (*implics)->vars[varfixing][k-1];
         (*implics)->types[varfixing][k] = (*implics)->types[varfixing][k-1];
         (*implics)->bounds[varfixing][k] = (*implics)->bounds[varfixing][k-1];
         (*implics)->ids[varfixing][k] = (*implics)->ids[varfixing][k-1];
      }
      assert(posadd == k);
      (*implics)->vars[varfixing][posadd] = implvar;
      (*implics)->types[varfixing][posadd] = impltype;
      (*implics)->bounds[varfixing][posadd] = implbound;
      (*implics)->ids[varfixing][posadd] = stat->nimplications;
      if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
         (*implics)->nbinimpls[varfixing]++;
      (*implics)->nimpls[varfixing]++;
#ifndef NDEBUG
      for( k = posadd-1; k >= 0; k-- )
         assert(compVars((void*)(*implics)->vars[varfixing][k], (void*)implvar) <= 0);
#endif
      stat->nimplications++;
   }
    
   checkImplics(*implics, set);

   return SCIP_OKAY;
}

/** removes the implication  x <= 0 or x >= 1  ==>  y <= b  or  y >= b  from the implications data structure */
RETCODE SCIPimplicsDel(
   IMPLICS**        implics,            /**< pointer to implications data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             varfixing,          /**< FALSE if y should be removed from implications for x <= 0, TRUE for x >= 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype            /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   )
{
   int i;
   int poslower;
   int posupper; 
   int posadd;
   Bool found;

   assert(implics != NULL);
   assert(*implics != NULL);
   assert(implvar != NULL);

   /* searches for y in implications of x */
   found = implicsSearchVar(*implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
   if( !found )
      return SCIP_OKAY;

   assert((impltype == SCIP_BOUNDTYPE_LOWER && poslower < INT_MAX && posadd == poslower) 
      || (impltype == SCIP_BOUNDTYPE_UPPER && posupper < INT_MAX && posadd == posupper));
   assert(0 <= posadd && posadd < (*implics)->nimpls[varfixing]);
   assert((SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY) == (posadd < (*implics)->nbinimpls[varfixing]));
   assert((*implics)->vars[varfixing][posadd] == implvar);
   assert((*implics)->types[varfixing][posadd] == impltype);

   /* removes y from implications of x */
   for( i = posadd; i < (*implics)->nimpls[varfixing] - 1; i++ )
   {
      (*implics)->vars[varfixing][i] = (*implics)->vars[varfixing][i+1];
      (*implics)->types[varfixing][i] = (*implics)->types[varfixing][i+1];
      (*implics)->bounds[varfixing][i] = (*implics)->bounds[varfixing][i+1];
   }
   (*implics)->nimpls[varfixing]--;
   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      assert(posadd < (*implics)->nbinimpls[varfixing]);
      (*implics)->nbinimpls[varfixing]--;
   }

   /* free implics data structure, if it is empty */
   if( (*implics)->nimpls[0] == 0 && (*implics)->nimpls[1] == 0 )
      SCIPimplicsFree(implics, blkmem);

   return SCIP_OKAY;
}

/** returns whether an implication y <= b or y >= b is contained in implications for x == 0 or x == 1 */
Bool SCIPimplicsContainsImpl(
   IMPLICS*         implics,            /**< implications data structure */
   Bool             varfixing,          /**< FALSE if y should be searched in implications for x == 0, TRUE for x == 1 */
   VAR*             implvar,            /**< variable y to search for */
   BOUNDTYPE        impltype            /**< type of implication y <=/>= b to search for */
   )
{
   int poslower;
   int posupper;
   int posadd;

   return implicsSearchVar(implics, varfixing, implvar, impltype, &poslower, &posupper, &posadd);
}




/*
 * methods for cliques
 */

/** creates a clique data structure */
static
RETCODE cliqueCreate(
   CLIQUE**         clique,             /**< pointer to store clique data structure */
   BLKMEM*          blkmem,             /**< block memory */
   int              size,               /**< initial size of clique */
   int              id                  /**< unique identifier of the clique */
   )
{
   assert(clique != NULL);

   ALLOC_OKAY( allocBlockMemory(blkmem, clique) );
   if( size > 0 )
   {
      ALLOC_OKAY( allocBlockMemoryArray(blkmem, &(*clique)->vars, size) );
      ALLOC_OKAY( allocBlockMemoryArray(blkmem, &(*clique)->values, size) );
   }
   else
   {
      (*clique)->vars = NULL;
      (*clique)->values = NULL;
   }
   (*clique)->nvars = 0;
   (*clique)->size = size;
   (*clique)->id = id;

   return SCIP_OKAY;
}

/** frees a clique data structure */
static
void cliqueFree(
   CLIQUE**         clique,             /**< pointer to store clique data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(clique != NULL);

   if( *clique != NULL )
   {
      freeBlockMemoryArrayNull(blkmem, &(*clique)->vars, (*clique)->size);
      freeBlockMemoryArrayNull(blkmem, &(*clique)->values, (*clique)->size);
      freeBlockMemory(blkmem, clique);
   }
}

/** ensures, that clique arrays can store at least num entries */
static
RETCODE cliqueEnsureSize(
   CLIQUE*          clique,             /**< clique data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(clique != NULL);
   
   if( num > clique->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &clique->vars, clique->size, newsize) );
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &clique->values, clique->size, newsize) );
      clique->size = newsize;
   }
   assert(num <= clique->size);

   return SCIP_OKAY;
}

/** adds a single variable to the given clique */
RETCODE SCIPcliqueAddVar(
   CLIQUE*          clique,             /**< clique data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to the clique */
   Bool             value,              /**< value of the variable in the clique */
   Bool*            doubleentry,        /**< pointer to store whether the variable and value occurs twice in the clique */
   Bool*            oppositeentry       /**< pointer to store whether the variable with opposite value is in the clique */
   )
{
   int varidx;
   int i;

   assert(clique != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(doubleentry != NULL);
   assert(oppositeentry != NULL);

   debugMessage("adding variable <%s> == %d to clique %d\n", SCIPvarGetName(var), value, clique->id);

   *doubleentry = FALSE;
   *oppositeentry = FALSE;

   /* allocate memory */
   CHECK_OKAY( cliqueEnsureSize(clique, blkmem, set, clique->nvars+1) );

   /* store variable in clique, sorted by index and values */
   varidx = SCIPvarGetIndex(var);
   assert(varidx >= 0);
   for( i = clique->nvars; i > 0 && SCIPvarGetIndex(clique->vars[i-1]) > varidx; --i )
   {
      clique->vars[i] = clique->vars[i-1];
      clique->values[i] = clique->values[i-1];
   }
   clique->vars[i] = var;
   for( ; i > 0 && clique->vars[i-1] == var && clique->values[i-1] > value; --i )
      clique->values[i] = clique->values[i-1];
   clique->values[i] = value;
   clique->nvars++;

   /* check whether the variable is contained twice in the clique */
   for( i--; i >= 0 && clique->vars[i] == var; --i )
   {
      *doubleentry = *doubleentry || (clique->values[i] == value);
      *oppositeentry = *oppositeentry || (clique->values[i] != value);
   }

   return SCIP_OKAY;
}

/** gets the position of the given variable in the clique; returns -1 if variable is not member of clique */
static
int cliqueSearchVar(
   CLIQUE*          clique,             /**< clique data structure */
   VAR*             var,                /**< variable to search for */
   Bool             value               /**< value of the variable in the clique */
   )
{
   int varidx;
   int left;
   int right;

   assert(clique != NULL);

   varidx = SCIPvarGetIndex(var);
   left = -1;
   right = clique->nvars;
   while( left < right-1 )
   {
      int middle;
      int idx;

      middle = (left+right)/2;
      idx = SCIPvarGetIndex(clique->vars[middle]);
      assert(idx >= 0);
      if( varidx < idx )
         right = middle;
      else if( varidx > idx )
         left = middle;
      else
      {
         assert(var == clique->vars[middle]);

         /* now watch out for the correct value */
         if( clique->values[middle] < value )
         {
            int i;
            for( i = middle+1; i < clique->nvars && clique->vars[i] == var; ++i )
            {
               if( clique->values[i] == value )
                  return i;
            }
            return -1;
         }
         if( clique->values[middle] > value )
         {
            int i;
            for( i = middle-1; i >= 0 && clique->vars[i] == var; --i )
            {
               if( clique->values[i] == value )
                  return i;
            }
         }
         return middle;
      }
   }

   return -1;
}

/** removes a single variable from the given clique */
RETCODE SCIPcliqueDelVar(
   CLIQUE*          clique,             /**< clique data structure */
   VAR*             var,                /**< variable to remove from the clique */
   Bool             value               /**< value of the variable in the clique */
   )
{
   int pos;

   assert(clique != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   debugMessage("deleting variable <%s> == %d from clique %d\n", SCIPvarGetName(var), value, clique->id);

   /* find variable in clique */
   pos = cliqueSearchVar(clique, var, value);
   assert(0 <= pos && pos < clique->nvars);
   assert(clique->vars[pos] == var);
   assert(clique->values[pos] == value);

   /* remove entry from clique */
   for( ; pos < clique->nvars-1; ++pos )
   {
      clique->vars[pos] = clique->vars[pos+1];
      clique->values[pos] = clique->values[pos+1];
   }
   clique->nvars--;

   return SCIP_OKAY;
}

/** gets the position of the given clique in the cliques array; returns -1 if clique is not member of cliques array */
static
int cliquesSearchClique(
   CLIQUE**         cliques,            /**< array of cliques */
   int              ncliques,           /**< number of cliques in the cliques array */
   CLIQUE*          clique              /**< clique to search for */
   )
{
   int cliqueid;
   int left;
   int right;

   assert(cliques != NULL || ncliques == 0);
   assert(clique != NULL);

   cliqueid = clique->id;
   assert(cliqueid >= 0);
   left = -1;
   right = ncliques;
   while( left < right-1 )
   {
      int middle;
      int id;

      middle = (left+right)/2;
      id = cliques[middle]->id;
      assert(id >= 0);
      if( cliqueid < id )
         right = middle;
      else if( cliqueid > id )
         left = middle;
      else
      {
         assert(clique == cliques[middle]);
         return middle;
      }
   }

   return -1;
}

#ifndef NDEBUG
/** checks whether clique appears in all clique lists of the involved variables */
static
void cliqueCheck(
   CLIQUE*          clique              /**< clique data structure */
   )
{
   int i;

   assert(clique != NULL);

   for( i = 0; i < clique->nvars; ++i )
   {
      CLIQUE** cliques;
      int ncliques;
      int pos;

      assert(i == 0 || SCIPvarGetIndex(clique->vars[i-1]) <= SCIPvarGetIndex(clique->vars[i]));
      assert(i == 0 || clique->vars[i-1] != clique->vars[i] || clique->values[i-1] <= clique->values[i]);
      ncliques = SCIPvarGetNCliques(clique->vars[i], clique->values[i]);
      cliques = SCIPvarGetCliques(clique->vars[i], clique->values[i]);
      pos = cliquesSearchClique(cliques, ncliques, clique);
      assert(0 <= pos && pos < ncliques);
      assert(cliques[pos] == clique);
   }
}
#else
#define cliqueCheck(clique) /**/
#endif

/** creates a clique list data structure */
static
RETCODE cliquelistCreate(
   CLIQUELIST**     cliquelist,         /**< pointer to store clique list data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(cliquelist != NULL);

   ALLOC_OKAY( allocBlockMemory(blkmem, cliquelist) );
   (*cliquelist)->cliques[0] = NULL;
   (*cliquelist)->cliques[1] = NULL;
   (*cliquelist)->ncliques[0] = 0;
   (*cliquelist)->ncliques[1] = 0;
   (*cliquelist)->size[0] = 0;
   (*cliquelist)->size[1] = 0;

   return SCIP_OKAY;
}

/** frees a clique list data structure */
void SCIPcliquelistFree(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   assert(cliquelist != NULL);

   if( *cliquelist != NULL )
   {
      freeBlockMemoryArrayNull(blkmem, &(*cliquelist)->cliques[0], (*cliquelist)->size[0]);
      freeBlockMemoryArrayNull(blkmem, &(*cliquelist)->cliques[1], (*cliquelist)->size[1]);
      freeBlockMemory(blkmem, cliquelist);
   }
}

/** ensures, that clique list arrays can store at least num entries */
static
RETCODE cliquelistEnsureSize(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             value,              /**< value of the variable for which the clique list should be extended */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(cliquelist != NULL);

   if( *cliquelist == NULL )
   {
      CHECK_OKAY( cliquelistCreate(cliquelist, blkmem) );
   }
   assert(*cliquelist != NULL);

   if( num > (*cliquelist)->size[value] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(blkmem, &(*cliquelist)->cliques[value], (*cliquelist)->size[value],
            newsize) );
      (*cliquelist)->size[value] = newsize;
   }
   assert(num <= (*cliquelist)->size[value]);

   return SCIP_OKAY;
}

/** adds a clique to the clique list */
RETCODE SCIPcliquelistAdd(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Bool             value,              /**< value of the variable for which the clique list should be extended */
   CLIQUE*          clique              /**< clique that should be added to the clique list */
   )
{
   int id;
   int i;

   assert(cliquelist != NULL);

   /* allocate memory */
   CHECK_OKAY( cliquelistEnsureSize(cliquelist, blkmem, set, value, SCIPcliquelistGetNCliques(*cliquelist, value)+1) );
   assert(*cliquelist != NULL);
   assert((*cliquelist)->cliques[value] != NULL);

   debugMessage("adding clique %d to cliquelist %p value %d (length: %d)\n", 
      clique->id, *cliquelist, value, (*cliquelist)->ncliques[value]);
   
   /* insert clique into list, sorted by clique id */
   id = clique->id;
   assert(id >= 0);
   for( i = (*cliquelist)->ncliques[value]; i > 0 && (*cliquelist)->cliques[value][i-1]->id > id; --i )
      (*cliquelist)->cliques[value][i] = (*cliquelist)->cliques[value][i-1];
   (*cliquelist)->cliques[value][i] = clique;
   (*cliquelist)->ncliques[value]++;

   return SCIP_OKAY;
}

/** removes a clique from the clique list */
RETCODE SCIPcliquelistDel(
   CLIQUELIST**     cliquelist,         /**< pointer to the clique list data structure */
   BLKMEM*          blkmem,             /**< block memory */
   Bool             value,              /**< value of the variable for which the clique list should be reduced */
   CLIQUE*          clique              /**< clique that should be deleted from the clique list */
   )
{
   int pos;

   assert(cliquelist != NULL);
   assert(*cliquelist != NULL);

   debugMessage("deleting clique %d from cliquelist %p value %d (length: %d)\n", 
      clique->id, *cliquelist, value, (*cliquelist)->ncliques[value]);
   
   pos = cliquesSearchClique((*cliquelist)->cliques[value], (*cliquelist)->ncliques[value], clique);
   assert(0 <= pos && pos < (*cliquelist)->ncliques[value]);
   assert((*cliquelist)->cliques[value][pos] == clique);

   /* remove clique from list */
   for( ; pos < (*cliquelist)->ncliques[value] - 1; ++pos )
      (*cliquelist)->cliques[value][pos] = (*cliquelist)->cliques[value][pos+1];
   (*cliquelist)->ncliques[value]--;

   /* free cliquelist if it is empty */
   if( (*cliquelist)->ncliques[0] == 0 && (*cliquelist)->ncliques[1] == 0 )
      SCIPcliquelistFree(cliquelist, blkmem);

   return SCIP_OKAY;
}

/** returns whether the given clique lists have a non-empty intersection, i.e. whether there is a clique that appears
 *  in both lists
 */
Bool SCIPcliquelistsHaveCommonClique(
   CLIQUELIST*      cliquelist1,        /**< first clique list data structure */
   Bool             value1,             /**< value of first variable */
   CLIQUELIST*      cliquelist2,        /**< second clique list data structure */
   Bool             value2              /**< value of second variable */
   )
{
   CLIQUE** cliques1;
   CLIQUE** cliques2;
   int ncliques1;
   int ncliques2;
   int i1;
   int i2;

   if( cliquelist1 == NULL || cliquelist2 == NULL )
      return FALSE;

   ncliques1 = cliquelist1->ncliques[value1];
   cliques1 = cliquelist1->cliques[value1];
   ncliques2 = cliquelist2->ncliques[value2];
   cliques2 = cliquelist2->cliques[value2];
   i1 = 0;
   i2 = 0;
   while( i1 < ncliques1 && i2 < ncliques2 )
   {
      int cliqueid;

      cliqueid = SCIPcliqueGetId(cliques2[i2]);
      while( i1 < ncliques1 && SCIPcliqueGetId(cliques1[i1]) < cliqueid )
         i1++;
      if( i1 == ncliques1 )
         break;

      cliqueid = SCIPcliqueGetId(cliques1[i1]);
      while( i2 < ncliques2 && SCIPcliqueGetId(cliques2[i2]) < cliqueid )
         i2++;
      if( i2 == ncliques2 )
         break;

      if( SCIPcliqueGetId(cliques2[i2]) == cliqueid )
         return TRUE;
   }

   return FALSE;
}

/** removes all listed entries from the cliques */
void SCIPcliquelistRemoveFromCliques(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   VAR*             var                 /**< active problem variable the clique list belongs to */
   )
{
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   if( cliquelist != NULL )
   {
      int value;

      debugMessage("removing variable <%s> from cliques (%d with value 0, %d with value 1)\n",
         SCIPvarGetName(var), cliquelist->ncliques[0], cliquelist->ncliques[1]);

      for( value = 0; value < 2; ++value )
      {
         int i;

         assert(SCIPvarGetCliques(var, (Bool)value) == cliquelist->cliques[value]);
         assert(SCIPvarGetNCliques(var, (Bool)value) == cliquelist->ncliques[value]);
         for( i = 0; i < cliquelist->ncliques[value]; ++i )
         {
            CLIQUE* clique;
            int pos;

            clique = cliquelist->cliques[value][i];
            assert(clique != NULL);

            debugMessage(" -> removing variable <%s> == %d from clique %d (size %d)\n",
               SCIPvarGetName(var), value, clique->id, clique->nvars);

            /* binary search the position of the variable in the clique */
            pos = cliqueSearchVar(clique, var, (Bool)value);
            assert(0 <= pos && pos < clique->nvars);
            assert(clique->vars[pos] == var);
            assert(clique->values[pos] == (Bool)value);

            /* remove the entry from the clique */
            for( ; pos < clique->nvars-1; ++pos )
            {
               clique->vars[pos] = clique->vars[pos+1];
               clique->values[pos] = clique->values[pos+1];
            }
            clique->nvars--;

            cliqueCheck(clique);
         }
      }
   }
}

/** creates a clique table data structure */
RETCODE SCIPcliquetableCreate(
   CLIQUETABLE**    cliquetable         /**< pointer to store clique table data structure */
   )
{
   assert(cliquetable != NULL);

   ALLOC_OKAY( allocMemory(cliquetable) );
   (*cliquetable)->cliques = NULL;
   (*cliquetable)->ncliques = 0;
   (*cliquetable)->size = 0;
   (*cliquetable)->ncreatedcliques = 0;

   return SCIP_OKAY;
}

/** frees a clique table data structure */
RETCODE SCIPcliquetableFree(
   CLIQUETABLE**    cliquetable,        /**< pointer to store clique table data structure */
   BLKMEM*          blkmem              /**< block memory */
   )
{
   int i;

   assert(cliquetable != NULL);
   assert(*cliquetable != NULL);

   /* free all cliques */
   for( i = 0; i < (*cliquetable)->ncliques; ++i )
   {
      cliqueFree(&(*cliquetable)->cliques[i], blkmem);
   }

   /* free clique table data */
   freeMemoryArrayNull(&(*cliquetable)->cliques);
   freeMemory(cliquetable);

   return SCIP_OKAY;
}

/** ensures, that clique table arrays can store at least num entries */
static
RETCODE cliquetableEnsureSize(
   CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(cliquetable != NULL);

   if( num > cliquetable->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&cliquetable->cliques, newsize) );
      cliquetable->size = newsize;
   }
   assert(num <= cliquetable->size);

   return SCIP_OKAY;
}

/** adds a clique to the clique table; performs implications if the clique contains the same variable twice */
RETCODE SCIPcliquetableAdd(
   CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int              nvars,              /**< number of variables in the clique */
   Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*             nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   )
{
   CLIQUE* clique;
   int i;

   assert(cliquetable != NULL);
   assert(vars != NULL);

   debugMessage("adding clique %d with %d vars to clique table\n", cliquetable->ncliques, nvars);

   /* create the clique data structure */
   CHECK_OKAY( cliqueCreate(&clique, blkmem, nvars, cliquetable->ncreatedcliques) );
   cliquetable->ncreatedcliques++;

   /* add clique to clique table */
   CHECK_OKAY( cliquetableEnsureSize(cliquetable, set, cliquetable->ncliques+1) );
   cliquetable->cliques[cliquetable->ncliques] = clique;
   cliquetable->ncliques++;

   /* add the corresponding active problem variables to the clique */
   for( i = 0; i < nvars; ++i )
   {
      /* put the clique into the sorted clique table of the variable */
      CHECK_OKAY( SCIPvarAddClique(vars[i], blkmem, set, stat, lp, branchcand, eventqueue, TRUE, clique,
            infeasible, nbdchgs) );
   }

   cliqueCheck(clique);

   return SCIP_OKAY;
}

/** removes all empty and single variable cliques from the clique table, and converts all two variable cliques
 *  into implications
 */
RETCODE SCIPcliquetableCleanup(
   CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Bool*            infeasible          /**< pointer to store whether an infeasibility was detected */
   )
{
   int i;

   assert(cliquetable != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;
   i = 0;
   while( i < cliquetable->ncliques && !(*infeasible) )
   {
      CLIQUE* clique;

      clique = cliquetable->cliques[i];
      if( clique->nvars == 2 )
      {
         /* add the 2-clique as implication */
         CHECK_OKAY( SCIPvarAddImplic(clique->vars[0], blkmem, set, stat, lp, branchcand, eventqueue, clique->values[0],
               clique->vars[1], clique->values[1] ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
               (Real)(!clique->values[1]), infeasible, NULL) );

         /* delete the clique: remove one variable - the rest is done below */
         CHECK_OKAY( SCIPvarDelClique(clique->vars[0], blkmem, clique->values[0], clique) );
         assert(clique->nvars == 1);
      }
      if( clique->nvars == 1 )
      {
         /* a clique with only one variable is redundant */
         CHECK_OKAY( SCIPvarDelClique(clique->vars[0], blkmem, clique->values[0], clique) );
         assert(clique->nvars == 0);
      }
      if( clique->nvars == 0 )
      {
         /* remove empty cliques from clique table */
         cliqueFree(&cliquetable->cliques[i], blkmem);
         cliquetable->cliques[i] = cliquetable->cliques[cliquetable->ncliques-1];
         cliquetable->ncliques--;
      }
      else
         i++;
   }

   return SCIP_OKAY;
}


/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPvboundsGetNVbds
#undef SCIPvboundsGetVars
#undef SCIPvboundsGetCoefs
#undef SCIPvboundsGetConstants
#undef SCIPimplicsGetNImpls
#undef SCIPimplicsGetNBinImpls
#undef SCIPimplicsGetVars
#undef SCIPimplicsGetTypes
#undef SCIPimplicsGetBounds
#undef SCIPimplicsGetIds
#undef SCIPcliqueGetNVars
#undef SCIPcliqueGetVars
#undef SCIPcliqueGetValues
#undef SCIPcliqueGetId
#undef SCIPcliquelistGetNCliques
#undef SCIPcliquelistGetCliques
#undef SCIPcliquetableGetNCliques
#undef SCIPcliquelistCheck
#undef SCIPcliquetableGetNCliques

/** gets number of variable bounds contained in given variable bounds data structure */
int SCIPvboundsGetNVbds(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   assert(vbounds != NULL);

   return vbounds->len;
}

/** gets array of variables contained in given variable bounds data structure */
VAR** SCIPvboundsGetVars(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   assert(vbounds != NULL);

   return vbounds->vars;
}

/** gets array of coefficients contained in given variable bounds data structure */
Real* SCIPvboundsGetCoefs(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   assert(vbounds != NULL);

   return vbounds->coefs;
}

/** gets array of constants contained in given variable bounds data structure */
Real* SCIPvboundsGetConstants(
   VBOUNDS*         vbounds             /**< variable bounds data structure */
   )
{
   assert(vbounds != NULL);

   return vbounds->constants;
}

/** gets number of implications for a given binary variable fixing */
int SCIPimplicsGetNImpls(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   assert(implics != NULL);

   return implics->nimpls[varfixing];
}

/** gets number of implications on binary variables for a given binary variable fixing */
int SCIPimplicsGetNBinImpls(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   assert(implics != NULL);

   return implics->nbinimpls[varfixing];
}

/** gets array with implied variables for a given binary variable fixing */
VAR** SCIPimplicsGetVars(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   assert(implics != NULL);

   return implics->vars[varfixing];
}

/** gets array with implication types for a given binary variable fixing */
BOUNDTYPE* SCIPimplicsGetTypes(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   assert(implics != NULL);

   return implics->types[varfixing];
}

/** gets array with implication bounds for a given binary variable fixing */
Real* SCIPimplicsGetBounds(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   assert(implics != NULL);

   return implics->bounds[varfixing];
}

/** gets array with unique implication identifiers for a given binary variable fixing */
int* SCIPimplicsGetIds(
   IMPLICS*         implics,            /**< implication data */
   Bool             varfixing           /**< should the implications on var == FALSE or var == TRUE be returned? */
   )
{
   assert(implics != NULL);

   return implics->ids[varfixing];
}

/** gets number of variables in the cliques */
int SCIPcliqueGetNVars(
   CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->nvars;
}

/** gets array of active problem variables in the cliques */
VAR** SCIPcliqueGetVars(
   CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->vars;
}

/** gets array of values of active problem variables in the cliques, i.e. whether the variable is fixed to FALSE or
 *  to TRUE in the clique
 */
Bool* SCIPcliqueGetValues(
   CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->values;
}

/** gets unique identifier of the clique */
int SCIPcliqueGetId(
   CLIQUE*          clique              /**< clique data structure */
   )
{
   assert(clique != NULL);

   return clique->id;
}
   
/** returns the number of cliques stored in the clique list */
int SCIPcliquelistGetNCliques(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   Bool             value               /**< value of the variable for which the cliques should be returned */
   )
{
   return cliquelist != NULL ? cliquelist->ncliques[value] : 0;
}

/** returns the cliques stored in the clique list, or NULL if the clique list is empty */
CLIQUE** SCIPcliquelistGetCliques(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   Bool             value               /**< value of the variable for which the cliques should be returned */
   )
{
   return cliquelist != NULL ? cliquelist->cliques[value] : NULL;
}

/** checks whether variable is contained in all cliques of the cliquelist */
void SCIPcliquelistCheck(
   CLIQUELIST*      cliquelist,         /**< clique list data structure */
   VAR*             var                 /**< variable, the clique list belongs to */
   )
{
   int value;

   assert(cliquelist != NULL);
   assert(SCIPvarGetNCliques(var, FALSE) == cliquelist->ncliques[0]);
   assert(SCIPvarGetCliques(var, FALSE) == cliquelist->cliques[0]);
   assert(SCIPvarGetNCliques(var, TRUE) == cliquelist->ncliques[1]);
   assert(SCIPvarGetCliques(var, TRUE) == cliquelist->cliques[1]);

   for( value = 0; value < 2; ++value )
   {
      int i;
      
      for( i = 0; i < cliquelist->ncliques[value]; ++i )
      {
         CLIQUE* clique;
         int pos;

         clique = cliquelist->cliques[value][i];
         pos = cliqueSearchVar(clique, var, (Bool)value);
         assert(0 <= pos && pos < clique->nvars);
         assert(clique->vars[pos] == var);
         assert(clique->values[pos] == (Bool)value);
      }
   }
}

/** gets the number of cliques stored in the clique table */
int SCIPcliquetableGetNCliques(
   CLIQUETABLE*     cliquetable         /**< clique table data structure */
   )
{
   assert(cliquetable != NULL);

   return cliquetable->ncliques;
}
