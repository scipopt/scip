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
#pragma ident "@(#) $Id: implics.c,v 1.2 2005/08/09 16:27:06 bzfpfend Exp $"

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

/** searches if variable y is already contained in implications for x == 0 or x == 1
 *  y can be contained in structure with y >= b (y_lower) and y <= b (y_upper) 
 */
static
RETCODE implicsSearchVar(
   IMPLICS*         implics,            /**< implications data structure */
   VAR*             implvar,            /**< variable y to search for */
   BOUNDTYPE        impltype,           /**< type of implication y <=/>= b to search for */
   Bool             varfixing,          /**< FALSE if y is searched in implications for x == 0, TRUE for x == 1 */
   int*             poslower,           /**< pointer to store position of y_lower (inf if not found) */
   int*             posupper,           /**< pointer to store position of y_upper (inf if not found) */
   int*             posadd,             /**< pointer to store correct position (with respect to impltype) to add y */
   Bool*            found               /**< pointer to store whether an implication on the same bound exists */
   )
{
   int implvaridx;
   int left;
   int right;
   int middle;

   assert(implics != NULL);
   assert(poslower != NULL);
   assert(posupper != NULL);
   assert(posadd != NULL);
   assert(found != NULL);

   /* set left and right pointer */
   if( SCIPvarGetType(implvar) == SCIP_VARTYPE_BINARY )
   {
      if( implics->nbinimpls[varfixing] == 0 )
      {
         /* there are no implications with binary variable y */
         *posadd = 0;
         *poslower = INT_MAX;
         *posupper = INT_MAX;
         *found = FALSE;
          return SCIP_OKAY;
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
         *found = FALSE;
         return SCIP_OKAY;
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

   if( left > right )
   {
      /* y was not found */
      assert(right == -1 || compVars((void*)implics->vars[varfixing][right], (void*)implvar) < 0);
      assert(left >= implics->nimpls[varfixing] || implics->vars[varfixing][left] != implvar);
      *poslower = INT_MAX;
      *posupper = INT_MAX;
      *posadd = left;
      *found = FALSE;
   }
   else
   {
      /* y was found */
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
            *found = TRUE;
         }
         else
         {
            *posadd = *posupper;
            *found = FALSE;
         }
      }     
      else
      {
         if( *posupper < INT_MAX )
         {
            *posadd = *posupper;
            *found = TRUE;
         }
         else
         {
            *posadd = (*poslower)+1;
            *found = FALSE;
         }
      }
      assert(*posadd < INT_MAX);
   }

   return SCIP_OKAY;
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
      CHECK_OKAY( implicsSearchVar(*implics, implvar, impltype, varfixing, &poslower, &posupper, &posadd, &found) );
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
   CHECK_OKAY( implicsSearchVar(*implics, implvar, impltype, varfixing, &poslower, &posupper, &posadd, &found) );
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




/*
 * methods for cliques
 */

/** creates a clique data structure */
static
RETCODE cliqueCreate(
   CLIQUE**         clique,             /**< pointer to store clique data structure */
   BLKMEM*          blkmem,             /**< block memory */
   int              size                /**< initial size of clique */
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
   (*clique)->tablepos = -1;

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
   Bool             value               /**< value of the variable in the clique */
   )
{
   int probindex;
   int i;

   assert(clique != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

   /* allocate memory */
   CHECK_OKAY( cliqueEnsureSize(clique, blkmem, set, clique->nvars+1) );

   /* store variable in clique, sorted by probindex */
   probindex = SCIPvarGetProbindex(var);
   assert(probindex >= 0);
   for( i = clique->nvars; i > 0 && SCIPvarGetProbindex(clique->vars[i-1]) > probindex; --i )
   {
      clique->vars[i] = clique->vars[i-1];
      clique->values[i] = clique->values[i-1];
   }
   clique->vars[i] = var;
   clique->values[i] = value;
   clique->nvars++;

   return SCIP_OKAY;
}

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
   int tablepos;
   int i;

   assert(cliquelist != NULL);

   /* allocate memory */
   CHECK_OKAY( cliquelistEnsureSize(cliquelist, blkmem, set, value, SCIPcliquelistGetNCliques(*cliquelist, value)+1) );
   assert(*cliquelist != NULL);
   assert((*cliquelist)->cliques[value] != NULL);

   /* insert clique into list, sorted by pointer */
   tablepos = clique->tablepos;
   assert(tablepos >= 0);
   for( i = (*cliquelist)->ncliques[value]; i > 0 && (*cliquelist)->cliques[value][i-1]->tablepos > tablepos; --i )
      (*cliquelist)->cliques[value][i] = (*cliquelist)->cliques[value][i-1];
   (*cliquelist)->cliques[value][i] = clique;
   (*cliquelist)->ncliques[value]++;

   return SCIP_OKAY;
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
      int probindex;
      int value;

      probindex = SCIPvarGetProbindex(var);
      assert(probindex >= 0);

      for( value = 0; value < 2; ++value )
      {
         int i;

         assert(SCIPvarGetCliques(var, (Bool)value) == cliquelist->cliques[value]);
         assert(SCIPvarGetNCliques(var, (Bool)value) == cliquelist->ncliques[value]);
         for( i = 0; i < cliquelist->ncliques[value]; ++i )
         {
            CLIQUE* clique;
            int left;
            int right;

            clique = cliquelist->cliques[value][i];
            assert(clique != NULL);

            /* binary search the position of the variable in the clique */
            left = 0;
            right = clique->nvars-1;
            while( left < right )
            {
               int middle;
               int idx;

               middle = (left+right)/2;
               idx = SCIPvarGetProbindex(clique->vars[middle]);
               assert(idx >= 0);
               if( probindex < idx )
                  right = middle-1;
               else if( probindex > idx )
                  left = middle+1;
               else
                  break;
            }
            assert(var == clique->vars[left]);
            
            /* remove the entry from the clique (the variable might be contained in the clique with both values;
             * however, it does not matter which entry is deleted here, because we delete both entries anyways)
             */
            for( ; left < clique->nvars-1; ++left )
            {
               clique->vars[left] = clique->vars[left+1];
               clique->values[left] = clique->values[left+1];
            }
            clique->nvars--;
         }
      }
   }
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

/** adds a clique to the clique table */
RETCODE SCIPcliquetableAdd(
   CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int              nvars               /**< number of variables in the clique */
   )
{
   CLIQUE* clique;
   int i;

   assert(cliquetable != NULL);
   assert(vars != NULL);

   debugMessage("adding clique %d with %d vars to clique table:", cliquetable->ncliques, nvars);

   /* create the clique data structure */
   CHECK_OKAY( cliqueCreate(&clique, blkmem, nvars) );

   /* add clique to clique table */
   CHECK_OKAY( cliquetableEnsureSize(cliquetable, set, cliquetable->ncliques+1) );
   cliquetable->cliques[cliquetable->ncliques] = clique;
   clique->tablepos = cliquetable->ncliques;
   cliquetable->ncliques++;

   /* add the corresponding active problem variables to the clique */
   for( i = 0; i < nvars; ++i )
   {
      /* put the clique into the sorted clique table of the variable */
      CHECK_OKAY( SCIPvarAddClique(vars[i], blkmem, set, TRUE, clique) );
      debugPrintf(" <%s>", SCIPvarGetName(vars[i]));
   }
   debugPrintf("\n");

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
